using Test
using Random
using TorchPID

@testset "TorchPID Phase P0 Verification Suite" begin

    @testset "1. Geometry & Physics Constants" begin
        @test DETECTOR.active == 53.0
        @test DETECTOR.n_ypixels == 128
        @test DETECTOR.y_size == 53.0 / 128
        @test DETECTOR.x_size == 53.0 / 8
        @test CONSTANTS.N_AIR == 1.000273
        @test CONSTANTS.CLIGHT == 299.792458
    end

    @testset "2. Particle Construction & Rotation" begin
        p_kw = Particle(; pid = 211, pMag = 10.0, xDir = 0.0, yDir = 0.0, zDir = 1.0)
        @test p_kw.pid == 211
        @test p_kw.pMag == 10.0
        initial_rotation!(p_kw)
        rx, ry, rz = rotate(p_kw, 1.0, 2.0, 3.0)
        @test (rx, ry, rz) == (1.0, 2.0, 3.0)

        # Tilted particle rotation check (orthogonal transformation preserves norm)
        theta = 0.1
        phi = 0.3
        xd = sin(theta) * cos(phi)
        yd = sin(theta) * sin(phi)
        zd = cos(theta)
        p_tilt = Particle(211, 0.0, 0.0, 10.0, xd, yd, zd, 10xd, 10yd, 10zd, 10000.0, 0.0)
        initial_rotation!(p_tilt)
        # Rotating (0, 0, 1) in particle frame should yield particle direction (xd, yd, zd)
        ux, uy, uz = rotate(p_tilt, 0.0, 0.0, 1.0)
        @test isapprox(ux, xd; atol = 1e-14)
        @test isapprox(uy, yd; atol = 1e-14)
        @test isapprox(uz, zd; atol = 1e-14)
    end

    @testset "3. Photon Creation, Surface Roughness & Tracing" begin
        ctx = PhotonContext(RADIATOR, CONSTANTS)
        p = Particle(211, 0.0, 0.0, 10.0, 0.0, 0.01, sqrt(1.0 - 0.01^2), 0.0, 0.1, 10.0 * sqrt(1.0 - 0.01^2), 10000.0, 0.0)
        initial_rotation!(p)
        beta = particle_beta(p, get_particle_mass(211))

        energy = 3.5
        np = TorchPID.nphase_Corning(energy)
        ng = TorchPID.ngroup_Corning(energy)

        ph_exp = create_explicit_photon(p, beta, 0.5 * π, energy, np, ng, ctx)
        @test ph_exp.energy == energy
        @test isapprox(ph_exp.xdir^2 + ph_exp.ydir^2 + ph_exp.zdir^2, 1.0; atol = 1e-13)

        ph_approx1 = create_approximate_photon(p, beta, np, ng, energy, 5.0, 0.4, ctx)
        ph_approx2 = create_approximate_photon(p, beta, np, ng, energy, 5.0, cos(0.4), sin(0.4), ctx)
        @test isapprox(ph_approx1.xdir, ph_approx2.xdir; atol = 1e-15)
        @test isapprox(ph_approx1.ydir, ph_approx2.ydir; atol = 1e-15)
        @test isapprox(ph_approx1.zdir, ph_approx2.zdir; atol = 1e-15)

        # Surface roughness: 0 reflections always passes (r^0 == 1.0 > rand())
        @test test_z_surface_roughness(ph_exp, 0) == true
        @test test_z_surface_roughness(ph_exp, 0.5, 0) == true

        # PhotonMapper with surface_roughness = true must not throw MethodError
        mapper_sr = PhotonMapper(; surface_roughness = true)
        hit = trace_photon(mapper_sr, ph_exp, 0.0)
        @test hit isa Union{HitCoordinate, Nothing}
    end

    @testset "4. PhotonSpectrum Normalization & Thresholds" begin
        s = PhotonSpectrum()
        d_above = PhotonSpectrumDistribution(s, 0.999)
        @test TorchPID.spectrum_above_threshold(d_above)

        # Probability density integral over energy bins must equal 1.0
        integral = sum(spectrum_probability(s, d_above, e) * s.dE for e in s.energy)
        @test isapprox(integral, 1.0; atol = 1e-12)

        # Below Cherenkov threshold (beta = 0.5 < 1/n_max)
        d_below = PhotonSpectrumDistribution(s, 0.5)
        @test !TorchPID.spectrum_above_threshold(d_below)
        @test spectrum_probability(s, d_below, 3.5) == 0.0
        @test spectrum_random_sampling(s, d_below, 0.5) === nothing
    end

    @testset "5. PixelMapper Local MCP Coordinates" begin
        # Across all 11 MCPs (-330 to +330 mm), local MCP x-position must be within [-30, +30] mm
        for mcp_idx in 0:(DETECTOR.n_detectors - 1)
            mcp_center = -0.5 * RADIATOR.width + (mcp_idx + 0.5) * DETECTOR.width
            @test isapprox(TorchPID.get_mcp_xposition(mcp_center), 0.0; atol = 1e-12)
            @test isapprox(TorchPID.get_mcp_xposition(mcp_center + 12.3), 12.3; atol = 1e-12)
            @test isapprox(TorchPID.get_mcp_xposition(mcp_center - 15.0), -15.0; atol = 1e-12)
        end
    end

    @testset "6. ChargeDepositTester Lookup & Time Smearing" begin
        cdt = ChargeDepositTester()
        @test length(cdt.time_cache) == 2 * cdt.n_time + 1

        # Verify 1-based (i - 0.5) matches C++ 0-based (0 + 0.5) for first bin i = 1
        expected_xf1 = TorchPID._gaussian_integral(0.5 * cdt.space_step, cdt.point_spread, 0.0, DETECTOR.x_size)
        expected_yf1 = TorchPID._gaussian_integral(0.5 * cdt.space_step, cdt.point_spread, 0.0, DETECTOR.y_size)
        @test isapprox(cdt.xf[1], expected_xf1; atol = 1e-15)
        @test isapprox(cdt.yf[1], expected_yf1; atol = 1e-15)

        tc = TorchPID._set_resolution(cdt.time_resolution, cdt.n_time)
        @test tc == cdt.time_cache

        t_smeared = smear_time(cdt, 50.0)
        @test t_smeared isa Float64
        hit_smeared = smear_time(cdt, HitCoordinate(10.0, -5.0, 50.0))
        @test hit_smeared.x == 10.0
        @test hit_smeared.y == -5.0
    end

    @testset "7. TorchImage & TorchImageAccumulator" begin
        acc = TorchImageAccumulator(DETECTOR)
        cdt = ChargeDepositTester()
        @test fill!(acc, cdt, PixelHit(4, 64, 500, 1.0)) === true
        @test fill!(acc, cdt, PixelHit(-1, 64, 500, 1.0)) === false
        @test acc.filled == 1
        @test acc.integral > 0.0
    end

    @testset "8. FrontEnd Digitization, Charge Sharing & Deadtime" begin
        # Charge sharing uses physical 64-y-pixel anode pitch (as in PID_Algorithm/data/geometry.yaml)
        # or higher gain for 128-y-pixel pitch so single-pixel charge exceeds the 30 fC threshold
        sig_share = TorchPID.SignalParameters(gain = 1.5e6)
        cdt = ChargeDepositTester(DETECTOR, sig_share; charge_sharing = true)
        fe_share = FrontEnd(DETECTOR; sharing_mode = true, smear_time = false)
        mcps = create_mcp_images(fe_share)

        # Hit on MCP 5 (inside pixel ix=4, iy=64 at x = 3.0 mm, y = 0.2 mm, t = 50.0 ns)
        add_photon!(mcps, fe_share, cdt, HitCoordinate(3.0, 0.2, 50.0))
        hits = get_hits(mcps, fe_share)
        @test hits !== nothing
        @test !isempty(hits)
        @test all(h -> 0 <= h.x < DETECTOR.n_xtotal, hits)
        @test all(h -> 0 <= h.y < DETECTOR.n_ypixels, hits)

        # Deadtime test on a single pixel (no sharing)
        fe_noshare = FrontEnd(DETECTOR; sharing_mode = false, smear_time = false)
        mcps_dt = create_mcp_images(fe_noshare)
        # Add two hits separated by 1.0 ns (< 10.0 ns deadtime) and a third separated by 12.0 ns (> 10.0 ns deadtime)
        add_photon!(mcps_dt, fe_noshare, cdt, HitCoordinate(0.0, 0.0, 42.0))
        add_photon!(mcps_dt, fe_noshare, cdt, HitCoordinate(0.0, 0.0, 43.0))
        add_photon!(mcps_dt, fe_noshare, cdt, HitCoordinate(0.0, 0.0, 56.0))
        dt_hits = get_hits(mcps_dt, fe_noshare)
        @test length(dt_hits) == 2
        @test dt_hits[1].t < dt_hits[2].t
    end

    @testset "9. PatternMatcher Uniform Sampling & End-to-End" begin
        for _ in 1:100
            rp = TorchPID.random_point()
            @test 0.0 <= rp.e < 1.0
            @test 0.0 <= rp.emissionZ < RADIATOR.depth
            @test 0.0 <= rp.phi < 2 * π
            @test isapprox(rp.cosPhi, cos(rp.phi); atol = 1e-15)
            @test isapprox(rp.sinPhi, sin(rp.phi); atol = 1e-15)
        end

        ctx = PhotonContext(RADIATOR, CONSTANTS)
        s = PhotonSpectrum()
        mapper = PhotonMapper()
        acc = TorchImageAccumulator(DETECTOR)
        cdt = ChargeDepositTester()
        p = Particle(211, 0.0, 0.0, 10.0, 0.0, 0.05, sqrt(1.0 - 0.05^2), 0.0, 0.5, 10.0 * sqrt(1.0 - 0.05^2), 10000.0, 0.0)
        initial_rotation!(p)
        beta = particle_beta(p, get_particle_mass(211))

        coords = project_pattern(p, beta, mapper, s, ctx)
        @test coords isa Vector{HitCoordinate}

        yield_est = make_pattern(p, beta, mapper, s, acc, cdt, ctx, 1000)
        @test yield_est >= 0.0
    end

    @testset "10. TestBeam Config & Photon Generation" begin
        cfg_path = joinpath(@__DIR__, "..", "..", "PID_Algorithm", "config", "testBeamConfig.yaml")
        if isfile(cfg_path)
            tb_cfg = TestBeamSimulator(cfg_path)
            @test tb_cfg.entryX == 330.0
            @test tb_cfg.entryY == 233.6
            @test tb_cfg.sigmaEntryX == 2.9
            @test tb_cfg.momentumMean == 5.0
            @test tb_cfg.thetaMean == 0.087
            @test tb_cfg.pid == 2212
        end

        Random.seed!(42)
        tb = TestBeamSimulator()
        ctx = PhotonContext(RADIATOR, CONSTANTS)
        s = PhotonSpectrum()
        mapper = PhotonMapper()
        fe = FrontEnd(DETECTOR)
        cdt = ChargeDepositTester()

        particles = generate_particles(tb, RADIATOR, CONSTANTS, 5)
        @test length(particles) == 5
        batches = generate_photons_batch(particles, s, mapper, ctx, fe, cdt)
        @test length(batches) == 5
        for b in batches
            if b !== nothing
                @test length(b.xpixels) == b.npixels
                @test all(x -> 0.0 <= x < DETECTOR.n_xtotal, b.xpixels)
                @test all(y -> 0.0 <= y < DETECTOR.n_ypixels, b.ypixels)
                @test all(t -> 0.0 <= t < DETECTOR.n_tpixels, b.tpixels)
            end
        end
    end

    @testset "11. PIDAlgorithm Photon Selection & Coordinates" begin
        pidalg = PIDAlg()
        dht = DetectorHitTester()
        in_time = 0.5 * (DETECTOR.t_window_min + DETECTOR.t_window_max)
        out_time = DETECTOR.t_window_max + 10.0

        ph_in = PhotonHit(0.0, 0.0, 0.0, in_time, 3.5, 0.0, 0.0, 5.0, 0.0, in_time, Int32(1), Int32(1), false, false, false)
        ph_out = PhotonHit(0.0, 0.0, 0.0, out_time, 3.5, 0.0, 0.0, 5.0, 0.0, out_time, Int32(1), Int32(1), false, false, false)

        @test TorchPID.select_photon_hit(pidalg, RADIATOR, WEDGE, DETECTOR, MASK, ph_in) == true
        @test TorchPID.select_photon_hit(pidalg, RADIATOR, WEDGE, DETECTOR, MASK, ph_out) == false
    end

    @testset "12. Phase P1: Type Stability, isbits & Zero Allocations" begin
        @test isbitstype(Particle)
        @test isbitstype(PhotonHit)
        @test isbitstype(PixelHit)
        @test isbitstype(HitCoordinate)

        acc = TorchImageAccumulator()
        @test size(acc.image.data) == (DETECTOR.n_tpixels, DETECTOR.n_ypixels, DETECTOR.n_xtotal)
        @test stride(acc.image.data, 1) == 1

        # Warmup and check zero allocations & type inference on hot-path functions
        @test (@inferred TorchPID.epotek_305_interpolator(3.5)) isa Float64
        @test (@inferred TorchPID.mirror_reflect(3.5)) isa Float64
        @test (@inferred TorchPID.QE_interpolator(3.5)) isa Float64
        @test (@inferred get_particle_mass(211)) isa Float64

        @test (@allocated TorchPID.epotek_305_interpolator(3.5)) == 0
        @test (@allocated TorchPID.mirror_reflect(3.5)) == 0
        @test (@allocated TorchPID.QE_interpolator(3.5)) == 0
        @test (@allocated get_particle_mass(211)) == 0

        p = (@inferred Particle(211, 0.0, 0.0, 10.0, 0.1, 0.2, sqrt(1.0 - 0.05), 1.0, 2.0, 9.746794, 10000.0, 0.0))
        @test (@inferred rotate(p, 1.0, 2.0, 3.0)) isa Tuple{Float64,Float64,Float64}
        @test (@allocated rotate(p, 1.0, 2.0, 3.0)) == 0

        cdt = ChargeDepositTester()
        TorchPID.fill_smeared!(acc, cdt, 5, 10, 50)
        @test (@allocated TorchPID.fill_smeared!(acc, cdt, 5, 10, 50)) == 0

        fe = FrontEnd(DETECTOR)
        mcps = create_mcp_images(fe)
        add_photon!(mcps, fe, cdt, 0.0, 0.0, 25.0)
        reset!(mcps, fe)
        @test (@allocated reset!(mcps, fe)) == 0
    end

end
