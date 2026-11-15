"""
    TestBeamSimulator(
        entryX::Float64=100.0,
        entryY::Float64=100.0,
        sigmaEntryX::Float64=1.0,
        sigmaEntryY::Float64=1.0,
        momentumMean::Float64=5.0,
        momentumSigma::Float64=0.005,
        thetaMean::Float64=-0.001,
        thetaSigma::Float64=0.0002,
        phiMean::Float64=-0.001,
        phiSigma::Float64=0.0002,
        timingDistance::Float64=8000.0,
        timingResolution::Float64=0.015,
        pid::Int64=2212,
    )

Type epresenting a test beam simulator configuration for particle generation.

# Fields
- `entryX::Float64`: X-coordinate of the beam entry point.
- `entryY::Float64`: Y-coordinate of the beam entry point.
- `sigmaEntryX::Float64`: Standard deviation of the beam entry point in X.
- `sigmaEntryY::Float64`: Standard deviation of the beam entry point in Y.
- `momentumMean::Float64`: Mean momentum of the beam particles.
- `momentumSigma::Float64`: Standard deviation of the beam momentum.
- `thetaMean::Float64`: Mean theta angle of the beam direction.
- `thetaSigma::Float64`: Standard deviation of the theta angle.
- `phiMean::Float64`: Mean phi angle of the beam direction.
- `phiSigma::Float64`: Standard deviation of the phi angle.
- `timingDistance::Float64`: Distance from timing reference.
- `timingResolution::Float64`: Resolution of the timing measurement.
- `pid::Int64`: Particle ID code.

# Constructors

    TestBeamSimulator(
        entryX::Float64=100.0,
        entryY::Float64=100.0,
        sigmaEntryX::Float64=1.0,
        sigmaEntryY::Float64=1.0,
        momentumMean::Float64=5.0,
        momentumSigma::Float64=0.005,
        thetaMean::Float64=-0.001,
        thetaSigma::Float64=0.0002,
        phiMean::Float64=-0.001,
        phiSigma::Float64=0.0002,
        timingDistance::Float64=8000.0,
        timingResolution::Float64=0.015,
        pid::Int64=2212,
    )

Constructs a TestBeamSimulator with default or specified parameters.

## Keywords
- `entryX::Float64`: X-coordinate of the beam entry point. Default: 100.0.
- `entryY::Float64`: Y-coordinate of the beam entry point. Default: 100.0.
- `sigmaEntryX::Float64`: Standard deviation of the beam entry point in X. Default: 1.0.
- `sigmaEntryY::Float64`: Standard deviation of the beam entry point in Y. Default: 1.0.
- `momentumMean::Float64`: Mean momentum of the beam particles. Default: 5.0.
- `momentumSigma::Float64`: Standard deviation of the beam momentum. Default: 0.005.
- `thetaMean::Float64`: Mean theta angle of the beam direction. Default: -0.001.
- `thetaSigma::Float64`: Standard deviation of the theta angle. Default: 0.0002.
- `phiMean::Float64`: Mean phi angle of the beam direction. Default: -0.001.
- `phiSigma::Float64`: Standard deviation of the phi angle. Default: 0.0002.
- `timingDistance::Float64`: Distance from timing reference. Default: 8000.0.
- `timingResolution::Float64`: Resolution of the timing measurement. Default: 0.015.
- `pid::Int64`: Particle ID code. Default: 2212 (proton).

    TestBeamSimulator(config_file::String)
    
Constructs a TestBeamSimulator from a YAML configuration file.
The configuration file should contain sections for entry, beam, timingRef, and particle parameters.
Default values are used for any parameters not specified in the configuration file.

## Arguments
- `config_file::String`: Path to the YAML configuration file.
"""
struct TestBeamSimulator
    entryX::Float64
    entryY::Float64
    sigmaEntryX::Float64
    sigmaEntryY::Float64

    # Beam momentum properties
    momentumMean::Float64
    momentumSigma::Float64

    # Beam angle properties
    thetaMean::Float64
    thetaSigma::Float64
    phiMean::Float64
    phiSigma::Float64

    # Distance from timing reference and its resolution
    timingDistance::Float64
    timingResolution::Float64

    pid::Int64
end

function TestBeamSimulator(;
    entryX::Float64 = 100.0,
    entryY::Float64 = 100.0,
    sigmaEntryX::Float64 = 1.0,
    sigmaEntryY::Float64 = 1.0,
    momentumMean::Float64 = 5.0,
    momentumSigma::Float64 = 0.005,
    thetaMean::Float64 = -0.001,
    thetaSigma::Float64 = 0.0002,
    phiMean::Float64 = -0.001,
    phiSigma::Float64 = 0.0002,
    timingDistance::Float64 = 8000.0,
    timingResolution::Float64 = 0.015,
    pid::Int64 = 2212,
)
    TestBeamSimulator(
        entryX,
        entryY,
        sigmaEntryX,
        sigmaEntryY,
        momentumMean,
        momentumSigma,
        thetaMean,
        thetaSigma,
        phiMean,
        phiSigma,
        timingDistance,
        timingResolution,
        pid,
    )
end

function TestBeamSimulator(config_file::String)
    config = YAML.load_file(config_file)
    entry = get(config, "entry", Dict{String,Any}())
    beam = get(config, "beam", Dict{String,Any}())
    angle = get(beam, "angle", Dict{String,Any}())
    timing = get(config, "timingRef", Dict{String,Any}())
    particle = get(config, "particle", Dict{String,Any}())
    TestBeamSimulator(
        entryX = Float64(get(entry, "x", get(entry, "entryX", 100.0))),
        entryY = Float64(get(entry, "y", get(entry, "entryY", 100.0))),
        sigmaEntryX = Float64(get(entry, "sx", get(entry, "sigmaEntryX", 1.0))),
        sigmaEntryY = Float64(get(entry, "sy", get(entry, "sigmaEntryY", 1.0))),
        momentumMean = Float64(get(beam, "momentum", get(beam, "momentumMean", 5.0))),
        momentumSigma = Float64(get(beam, "smomentum", get(beam, "momentumSigma", 0.005))),
        thetaMean = Float64(get(angle, "theta", get(beam, "thetaMean", -0.001))),
        thetaSigma = Float64(get(angle, "stheta", get(beam, "thetaSigma", 0.0002))),
        phiMean = Float64(get(angle, "phi", get(beam, "phiMean", -0.001))),
        phiSigma = Float64(get(angle, "sphi", get(beam, "phiSigma", 0.0002))),
        timingDistance = Float64(get(timing, "distance", get(timing, "timingDistance", 8000.0))),
        timingResolution = Float64(get(timing, "resolution", get(timing, "timingResolution", 0.015))),
        pid = Int64(get(beam, "id", get(particle, "pid", 211))),
    )
end

"""
    TestBeamParticle(
        particle::Particle,
        energy::Float64,
        beta::Float64,
        timingPosition::Tuple{Float64,Float64},
        flightTime::Float64,
    )

Type storing the particle information from a simulated test beam event.

# Fields
- `particle::Particle`: The base particle object.
- `energy::Float64`: Total energy of the particle.
- `beta::Float64`: Relativistic beta (velocity/c) of the particle.
- `timingPosition::Tuple{Float64,Float64}`: Position (X,Y) where timing reference was measured.
- `flightTime::Float64`: Time of flight from timing reference to entry point.
"""
struct TestBeamParticle
    particle::Particle
    energy::Float64
    beta::Float64
    timingPosition::Tuple{Float64,Float64}
    flightTime::Float64
end

"""
    TestBeamPhotons(
        photon_yield::Int,
        npixels::Int,
        xpixels::Vector{Float64},
        ypixels::Vector{Float64},
        tpixels::Vector{Float64},
        mcps::Vector{Int},
        mcpcolums::Vector{Int},
        charge::Vector{Float64},
    )

Type storing the photon detection results from a simulated test beam event.

# Fields
- `photon_yield::Int`: Total number of photons generated.
- `npixels::Int`: Number of pixels with photon hits.
- `xpixels::Vector{Float64}`: X-positions of hit pixels.
- `ypixels::Vector{Float64}`: Y-positions of hit pixels.
- `tpixels::Vector{Float64}`: Time values of hit pixels.
- `mcps::Vector{Int}`: MCP numbers for each hit.
- `mcpcolums::Vector{Int}`: Column within MCP for each hit.
- `charge::Vector{Float64}`: Charge deposited for each hit.
"""
struct TestBeamPhotons
    photon_yield::Int
    npixels::Int
    xpixels::Vector{Float64}
    ypixels::Vector{Float64}
    tpixels::Vector{Float64}
    mcps::Vector{Int}
    mcpcolums::Vector{Int}
    charge::Vector{Float64}
end

"""
    TestBeamWorkspace(spectrum::PhotonSpectrum, fe::FrontEnd)

Per-thread reusable workspace for zero-temporary-allocation TestBeam photon simulation and digitization.
"""
struct TestBeamWorkspace
    mcp_images::MCPImageArray
    dist_buf::Vector{Float64}
    cum_buf::Vector{Float64}
    photons_buf::Vector{HitCoordinate}
    pixels_buf::Vector{PixelHit}
end

function TestBeamWorkspace(spectrum::PhotonSpectrum, fe::FrontEnd)::TestBeamWorkspace
    mcp_images = create_mcp_images(fe)
    dist_buf = Vector{Float64}(undef, spectrum.nbins)
    cum_buf = Vector{Float64}(undef, spectrum.nbins)
    photons_buf = sizehint!(HitCoordinate[], 128)
    pixels_buf = sizehint!(PixelHit[], 128)
    return TestBeamWorkspace(mcp_images, dist_buf, cum_buf, photons_buf, pixels_buf)
end

@inline function _generate_particles_range!(
    particles::Vector{TestBeamParticle},
    lo::Int,
    hi::Int,
    tb::TestBeamSimulator,
    mass::Float64,
    half_width::Float64,
    half_height::Float64,
    clight::Float64,
    has_timing_res::Bool,
)
    @inbounds for i = lo:hi
        ex = (tb.entryX + tb.sigmaEntryX * randn()) - half_width
        ey = (tb.entryY + tb.sigmaEntryY * randn()) - half_height

        mom = tb.momentumMean + tb.momentumSigma * randn()
        en = hypot(mom, mass)
        beta = mom / en

        th = tb.thetaMean + tb.thetaSigma * randn()
        ph = tb.phiMean + tb.phiSigma * randn()

        sin_th, cos_th = sincos(th)
        sin_ph, cos_ph = sincos(ph)

        x_dir = sin_th * cos_ph
        y_dir = sin_th * sin_ph
        z_dir = cos_th

        px = mom * x_dir
        py = mom * y_dir
        pz = mom * z_dir

        plen = tb.timingDistance / cos_th
        ftime = plen / (beta * clight)

        tpos_x = ex - plen * x_dir
        tpos_y = ey - plen * y_dir

        t0 = has_timing_res ? tb.timingResolution * randn() : 0.0

        particle = Particle(
            tb.pid,
            ex,
            ey,
            mom,
            x_dir,
            y_dir,
            z_dir,
            px,
            py,
            pz,
            plen,
            t0,
        )
        particles[i] = TestBeamParticle(
            particle,
            en,
            beta,
            (tpos_x, tpos_y),
            ftime,
        )
    end
    return nothing
end

"""
    generate_particles(
        tb::TestBeamSimulator,
        radiator::Radiator,
        constants::Constants,
        N::Int;
        nthreads::Int = Threads.nthreads(),
    )::Vector{TestBeamParticle}

Simulates `N` test beam particles based on the beam configuration, using `nthreads` threads when `nthreads > 1`.
"""
function generate_particles(
    tb::TestBeamSimulator,
    radiator::Radiator,
    constants::Constants,
    N::Int;
    nthreads::Int = Threads.nthreads(),
)::Vector{TestBeamParticle}
    mass = get_particle_mass(tb.pid)
    half_width = radiator.half_width
    half_height = 0.5 * radiator.height
    clight = constants.CLIGHT
    has_timing_res = tb.timingResolution != 0.0

    particles = Vector{TestBeamParticle}(undef, N)
    nchunks = min(max(nthreads, 1), N)
    if nchunks <= 1
        _generate_particles_range!(
            particles, 1, N, tb, mass, half_width, half_height, clight, has_timing_res
        )
    else
        Threads.@threads for chunk = 1:nchunks
            lo = div((chunk - 1) * N, nchunks) + 1
            hi = div(chunk * N, nchunks)
            _generate_particles_range!(
                particles, lo, hi, tb, mass, half_width, half_height, clight, has_timing_res
            )
        end
    end
    return particles
end

"""
    generate_photons(
        particle::TestBeamParticle,
        spectrum::PhotonSpectrum,
        mapper::PhotonMapper,
        context::PhotonContext,
        fe::FrontEnd,
        cdt::ChargeDepositTester,
        ws::TestBeamWorkspace,
    )::Union{Nothing,TestBeamPhotons}

Simulates and digitizes photons from a test beam particle using a pre-allocated `TestBeamWorkspace` (or `MCPImageArray`).
"""
function generate_photons(
    particle::TestBeamParticle,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    ws::TestBeamWorkspace,
)::Union{Nothing,TestBeamPhotons}
    photon_distribution = update_distribution!(ws.dist_buf, ws.cum_buf, spectrum, particle.beta)
    photons = project_pattern!(
        ws.photons_buf,
        particle.particle,
        particle.beta,
        mapper,
        spectrum,
        photon_distribution,
        context,
    )
    photon_yield = length(photons)

    reset!(ws.mcp_images, fe)
    for photon in photons
        add_photon!(ws.mcp_images, fe, cdt, photon)
    end

    pixels = get_hits!(ws.pixels_buf, ws.mcp_images, fe)
    if isempty(pixels)
        return nothing
    end

    npixels = length(pixels)

    xpixels = Vector{Float64}(undef, npixels)
    ypixels = Vector{Float64}(undef, npixels)
    tpixels = Vector{Float64}(undef, npixels)
    mcp = Vector{Int}(undef, npixels)
    mcpcolum = Vector{Int}(undef, npixels)
    charge = Vector{Float64}(undef, npixels)

    @inbounds for i = 1:npixels
        p = pixels[i]
        xpixels[i] = Float64(p.x)
        ypixels[i] = Float64(p.y)
        tpixels[i] = Float64(p.t)
        mcp[i] = div(p.x, fe.n_xpixels)
        mcpcolum[i] = rem(p.x, fe.n_xpixels)
        charge[i] = p.ch
    end

    return TestBeamPhotons(photon_yield, npixels, xpixels, ypixels, tpixels, mcp, mcpcolum, charge)
end

function generate_photons(
    particle::TestBeamParticle,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    mcp_images::MCPImageArray = create_mcp_images(fe),
)::Union{Nothing,TestBeamPhotons}
    ws = TestBeamWorkspace(
        mcp_images,
        Vector{Float64}(undef, spectrum.nbins),
        Vector{Float64}(undef, spectrum.nbins),
        sizehint!(HitCoordinate[], 128),
        sizehint!(PixelHit[], 128),
    )
    return generate_photons(particle, spectrum, mapper, context, fe, cdt, ws)
end

function _generate_photons_range!(
    results::Vector{Union{Nothing,TestBeamPhotons}},
    particles::Vector{TestBeamParticle},
    lo::Int,
    hi::Int,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    ws::TestBeamWorkspace,
)
    @inbounds for i = lo:hi
        results[i] = generate_photons(particles[i], spectrum, mapper, context, fe, cdt, ws)
    end
    return nothing
end

function generate_photons_batch(
    particles::Vector{TestBeamParticle},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester;
    nthreads::Int = Threads.nthreads(),
)::Vector{Union{Nothing,TestBeamPhotons}}
    N = length(particles)
    results = Vector{Union{Nothing,TestBeamPhotons}}(undef, N)
    nchunks = min(max(nthreads, 1), N)
    if nchunks <= 1
        ws = TestBeamWorkspace(spectrum, fe)
        _generate_photons_range!(results, particles, 1, N, spectrum, mapper, context, fe, cdt, ws)
    else
        Threads.@threads for chunk = 1:nchunks
            ws = TestBeamWorkspace(spectrum, fe)
            lo = div((chunk - 1) * N, nchunks) + 1
            hi = div(chunk * N, nchunks)
            _generate_photons_range!(results, particles, lo, hi, spectrum, mapper, context, fe, cdt, ws)
        end
    end
    return results
end
