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

struct TestBeamParticle
    particle::Particle
    energy::Float64
    beta::Float64
    timingPosition::Tuple{Float64,Float64}
    flightTime::Float64
end

struct TestBeamPhotons
    photon_yield::Int
    npixels::Int
    xpixels::Vector{Float64}
    ypixels::Vector{Float64}
    #mcps::Vector{Int}
    #mcpColums::Vector{Int}
    #charge::Vector{Float64}
end

struct TestBeamData
    eventNumber::Int
    particle::TestBeamParticle
    photons::TestBeamPhotons
end

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

function generate_photons(
    particle::TestBeamParticle,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
)::TestBeamPhotons
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
    npixels = length(photons)
    xpixels = Vector{Float64}(undef, npixels)
    ypixels = Vector{Float64}(undef, npixels)
    for i = 1:npixels
        xpixels[i] = photons[i].x
        ypixels[i] = photons[i].y
    end

    TestBeamPhotons(photon_yield, npixels, xpixels, ypixels)
end
