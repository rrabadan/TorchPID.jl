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
    generate_particle(tb::TestBeamSimulator)::TestBeamParticle

`generate_particle` simulates a particle based on the beam configuration. 
The entry position is sampled from a normal distribution with means `entryX` and `entryY` and standard deviations `sigmaEntryX` and `sigmaEntryY`. 
The momentum is sampled from a normal distribution with mean `momentumMean` and standard deviation `momentumSigma`. 
Direction angles (`theta` and `phi`) are also sampled from normal distributions with means `thetaMean` and `phiMean` 
and standard deviations `thetaSigma` and `phiSigma`. 
The timing reference position and flight time are calculated using the `timingDistance` parameter, accounting for the particle's velocity and path length.

# Arguments
- `tb::TestBeamSimulator`: Test beam simulator configuration.

# Returns
- `TestBeamParticle`: A simulated test beam particle.
"""
function generate_particles(
    tb::TestBeamSimulator,
    radiator::Radiator,
    constants::Constants,
    N::Int,
)::Vector{TestBeamParticle}
    mass = get_particle_mass(tb.pid)
    half_width = radiator.half_width
    half_height = 0.5 * radiator.height
    clight = constants.CLIGHT
    has_timing_res = tb.timingResolution != 0.0

    particles = Vector{TestBeamParticle}(undef, N)
    @inbounds for i = 1:N
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
        mcp_images::MCPImageArray = create_mcp_images(fe),
    )::Union{Nothing,TestBeamPhotons}

`generate_photons` simulates and processes photons resulting from a particle interaction. 
It first calculates the photon yield based on the particle's beta and the provided spectrum. 
The photon pattern is then projected using the mapper, and photon hits are added to MCP images. 
Hits are extracted from the MCP images based on the front-end configuration. 
If no valid pixel hits are detected, the function returns `nothing`.

# Arguments
- `particle::TestBeamParticle`: The particle generating the photons.
- `spectrum::PhotonSpectrum`: Spectrum parameters for photon generation.
- `mapper::PhotonMapper`: Mapper for photon projection.
- `context::PhotonContext`: Photon tracing context.
- `fe::FrontEnd`: Front-end electronics configuration.
- `cdt::ChargeDepositTester`: Configuration for charge deposit testing.
- `mcp_images::MCPImageArray`: Optional pre-allocated MCP image array to reuse across events.

# Returns
- `Union{Nothing,TestBeamPhotons}`: Photon detection results, or `nothing` if no valid hits were detected.
"""
function generate_photons(
    particle::TestBeamParticle,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    mcp_images::MCPImageArray = create_mcp_images(fe),
)::Union{Nothing,TestBeamPhotons}

    # Generate the number of photons
    photon_distribution = PhotonSpectrumDistribution(spectrum, particle.beta)
    photons = project_pattern(
        particle.particle,
        particle.beta,
        mapper,
        spectrum,
        photon_distribution,
        context,
    )
    photon_yield = length(photons)

    reset!(mcp_images, fe)
    for photon in photons
        add_photon!(mcp_images, fe, cdt, photon)
    end

    pixels = get_hits(mcp_images, fe)
    if isnothing(pixels)
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

    TestBeamPhotons(photon_yield, npixels, xpixels, ypixels, tpixels, mcp, mcpcolum, charge)
end

function generate_photons_batch(
    particles::Vector{TestBeamParticle},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
)::Vector{Union{Nothing,TestBeamPhotons}}
    N = length(particles)
    results = Vector{Union{Nothing,TestBeamPhotons}}(undef, N)
    mcp_images = create_mcp_images(fe)
    for i = 1:N
        results[i] = generate_photons(particles[i], spectrum, mapper, context, fe, cdt, mcp_images)
    end
    return results
end
