"""
Represents a test beam simulator configuration for particle generation.

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

"""
Constructs a TestBeamSimulator with default or specified parameters.

# Arguments
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

# Returns
- `TestBeamSimulator`: A configured test beam simulator instance.
"""
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

"""
Constructs a TestBeamSimulator from a YAML configuration file.

The configuration file should contain sections for entry, beam, timingRef, and particle parameters.
Default values are used for any parameters not specified in the configuration file.

# Arguments
- `config_file::String`: Path to the YAML configuration file.

# Returns
- `TestBeamSimulator`: A test beam simulator instance configured according to the file.
"""
function TestBeamSimulator(config_file::String)
    config = YAML.load_file(config_file)
    TestBeamSimulator(
        entryX = get(config["entry"], "entryX", 100.0),
        entryY = get(config["entry"], "entryY", 100.0),
        sigmaEntryX = get(config["entry"], "sigmaEntryX", 1.0),
        sigmaEntryY = get(config["entry"], "sigmaEntryY", 1.0),
        momentumMean = get(config["beam"], "momentumMean", 5.0),
        momentumSigma = get(config["beam"], "momentumSigma", 0.005),
        thetaMean = get(config["beam"], "thetaMean", -0.001),
        thetaSigma = get(config["beam"], "thetaSigma", 0.0002),
        phiMean = get(config["beam"], "phiMean", -0.001),
        phiSigma = get(config["beam"], "phiSigma", 0.0002),
        timingDistance = get(config["timingRef"], "timingDistance", 8000.0),
        timingResolution = get(config["timingRef"], "timingResolution", 0.015),
        pid = get(config["particle"], "pid", 211),
    )
end

"""
Represents a particle in a test beam experiment.

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
Stores the photon detection results from a test beam event.

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
Contains all data for a single test beam event.

# Fields
- `eventNumber::Int`: Sequential identifier for the event.
- `particle::TestBeamParticle`: The particle that generated the event.
- `photons::TestBeamPhotons`: Photon detection results for the event.
"""
struct TestBeamData
    eventNumber::Int
    particle::TestBeamParticle
    photons::TestBeamPhotons
end

"""
Generates a simulated particle according to the beam configuration.

Entry position is randomized according to a Normal distribution with means `entryX` and `entryY`
and standard deviations `sigmaEntryX` and `sigmaEntryY`.
Momentum is randomized according to a Normal distribution with mean `momentumMean`
and standard deviation `momentumSigma`.
Direction angles are randomized according to Normal distributions with means `thetaMean` and `phiMean`
and standard deviations `thetaSigma` and `phiSigma`.
The timing reference and flight time are calculated based on the `timingDistance` parameter.

# Arguments
- `tb::TestBeamSimulator`: Test beam simulator configuration.

# Returns
- `TestBeamParticle`: A simulated particle with randomized properties according to the beam configuration.
"""
function generate_particle(tb::TestBeamSimulator)::TestBeamParticle
    # Generate the entry point using Normal distribution
    entryX = tb.entryX + tb.sigmaEntryX * randn()
    entryY = tb.entryY + tb.sigmaEntryY * randn()

    entryX -= RADIATOR.half_width
    entryY -= 0.5 * RADIATOR.height / 2.0

    # Generate the momentum
    momentum = tb.momentumMean + tb.momentumSigma * randn()

    mass = get_particle_mass(tb.pid)
    energy = sqrt(momentum^2 + mass^2)
    beta = momentum / energy

    # Generate the angles
    theta = tb.thetaMean + tb.thetaSigma * randn()
    phi = tb.phiMean + tb.phiSigma * randn()

    # Generate the direction
    px = momentum * sin(theta) * cos(phi)
    py = momentum * sin(theta) * sin(phi)
    pz = momentum * cos(theta)

    pathlength = tb.timingDistance / cos(theta)
    flightTime = pathlength / (beta * CLIGHT)

    timingPositionX = entryX - pathlength * cos(phi) * sin(theta)
    timingPositionY = entryY - pathlength * sin(phi) * sin(theta)

    t0 = tb.timingResolution == 0 ? 0 : tb.timingResolution * randn()

    particle = Particle(
        pid = tb.pid,
        xCoord = entryX,
        yCoord = entryY,
        pMag = momentum,
        xDir = px / momentum,
        yDir = py / momentum,
        zDir = pz / momentum,
        recoPX = px,
        recoPY = py,
        recoPZ = pz,
        pathlength = pathlength,
        t0 = t0,
    )
    initRotation(particle)

    TestBeamParticle(particle, energy, beta, (timingPositionX, timingPositionY), flightTime)
end

"""
Generates and processes photons from a particle interaction.

First generates photons according to the spectrum and the particle's beta.
Projects the photon pattern using the provided mapper.
Creates MCP images and adds photon hits to them.
Extracts hits from the MCP images according to the front-end configuration.
Returns `nothing` if no valid pixel hits were found.

# Arguments
- `particle::TestBeamParticle`: The particle generating the photons.
- `spectrum::PhotonSpectrum`: Spectrum parameters for photon generation.
- `mapper::PhotonMapper`: Mapper for photon projection.
- `fe::FrontEnd`: Front-end electronics configuration.
- `cdt::ChargeDepositTester`: Configuration for charge deposit testing.

# Returns
- `Union{Nothing,TestBeamPhotons}`: Photon detection results, or `nothing` if no valid hits were detected.
"""
function generate_photons(
    particle::TestBeamParticle,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
)::Union{Nothing,TestBeamPhotons}

    # Generate the number of photons
    photon_distribution = PhotonSpectrumDistribution(spectrum, particle.beta)
    photons = project_pattern(
        particle.particle,
        particle.beta,
        spectrum,
        photon_distribution,
        mapper,
    )
    photon_yield = length(photons)
    println("Photon yield: ", photon_yield)
    #npixels = length(photons)
    #xpixels = Vector{Float64}(undef, npixels)
    #ypixels = Vector{Float64}(undef, npixels)
    #for i = 1:npixels
    #    xpixels[i] = photons[i].x
    #    ypixels[i] = photons[i].y
    #end

    mcp_images = create_mcp_images(fe)
    #println("MCP images created ", size(mcp_images[begin].pixelmap[begin]))
    for photon in photons
        add_photon!(mcp_images, fe, cdt, photon)
    end

    pixels = get_hits(mcp_images, fe)
    if isnothing(pixels)
        return nothing
    end

    npixels = length(pixels)

    xpixels = Vector{Float64}(undef, 0)
    ypixels = Vector{Float64}(undef, 0)
    tpixels = Vector{Float64}(undef, 0)

    mcp = Vector{Int}(undef, 0)
    mcpcolum = Vector{Int}(undef, 0)
    charge = Vector{Float64}(undef, 0)

    for pixel in pixels
        push!(xpixels, pixel.x)
        push!(ypixels, pixel.y)
        push!(tpixels, pixel.t)

        push!(mcp, pixel.x / fe.n_xpixels)
        push!(mcpcolum, pixel.x % fe.n_xpixels)
        push!(charge, pixel.charge)
    end

    TestBeamPhotons(
        photon_yield,
        npixels,
        xpixels,
        ypixels,
        tpixels,
        mcps,
        mcpcolums,
        charge,
    )
end
