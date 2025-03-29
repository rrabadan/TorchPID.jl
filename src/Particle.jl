""" Represents a particle with kinematic and identification fields. """
struct Particle
    pid::Int
    eventId::Int
    trackId::Int
    moduleId::Int
    xCoord::Float64
    yCoord::Float64
    pMag::Float64
    pMagTrue::Float64
    xDir::Float64
    yDir::Float64
    zDir::Float64
    truePX::Float64
    truePY::Float64
    truePZ::Float64
    recoPX::Float64
    recoPY::Float64
    recoPZ::Float64
    pathlength::Float64
    truePathlength::Float64
    t0::Float64
    rotMatrix::Array{Float64,2}
end

""" Outer constructor for Particle with default parameters.
    Returns a new Particle initialized with specified or default values.
"""
function Particle(;
    pid::Int = 211,
    eventId::Int = 0,
    trackId::Int = 0,
    moduleId::Int = 0,
    xCoord::Float64 = 0.0,
    yCoord::Float64 = 0.0,
    pMag::Float64 = 0.0,
    pMagTrue::Float64 = 0.0,
    xDir::Float64 = 0.0,
    yDir::Float64 = 0.0,
    zDir::Float64 = 0.0,
    truePX::Float64 = 0.0,
    truePY::Float64 = 0.0,
    truePZ::Float64 = 0.0,
    recoPX::Float64 = 0.0,
    recoPY::Float64 = 0.0,
    recoPZ::Float64 = 0.0,
    pathlength::Float64 = 0.0,
    truePathlength::Float64 = 0.0,
    t0::Float64 = 0.0,
    rotMatrix::Array{Float64,2} = Array{Float64,2}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]),
)
    Particle(
        pid,
        eventId,
        trackId,
        moduleId,
        xCoord,
        yCoord,
        pMag,
        pMagTrue,
        xDir,
        yDir,
        zDir,
        truePX,
        truePY,
        truePZ,
        recoPX,
        recoPY,
        recoPZ,
        pathlength,
        truePathlength,
        t0,
        rotMatrix,
    )
end

""" Calculates the beta factor for the particle.
    beta = pMag / √(pMag² + m²)
"""
function beta(p::Particle, m::Float64)::Float64
    return p.pMag / sqrt(p.pMag^2 + m^2)
end

""" Calculates the gamma factor for the particle.
    gamma = 1 / √(1 - beta²)
"""
function gamma(p::Particle, m::Float64)::Float64
    return 1 / sqrt(1 - beta(p, m)^2)
end

""" Initializes the rotation matrix of the particle based on its direction components.
    Updates the rotMatrix field in-place.
"""
function initRotation(p::Particle)
    uperp::Float64 = sqrt(p.xDir^2 + p.yDir^2)

    if uperp != 0
        p.rotMatrix[1, 1] = p.xDir * p.zDir / uperp
        p.rotMatrix[1, 2] = -p.yDir / uperp
        p.rotMatrix[1, 3] = p.xDir
        p.rotMatrix[2, 1] = p.yDir * p.zDir / uperp
        p.rotMatrix[2, 2] = p.xDir / uperp
        p.rotMatrix[2, 3] = p.yDir
        p.rotMatrix[3, 1] = -uperp
        p.rotMatrix[3, 2] = 0
        p.rotMatrix[3, 3] = p.zDir
    end
end

""" Rotates a 3D vector (x, y, z) using the particle's rotation matrix.
    Returns the rotated (x, y, z) tuple.
"""
function rotate(
    p::Particle,
    x::Float64,
    y::Float64,
    z::Float64,
)::Tuple{Float64,Float64,Float64}
    return (
        p.rotMatrix[1, 1] * x + p.rotMatrix[1, 2] * y + p.rotMatrix[1, 3] * z,
        p.rotMatrix[2, 1] * x + p.rotMatrix[2, 2] * y + p.rotMatrix[2, 3] * z,
        p.rotMatrix[3, 1] * x + p.rotMatrix[3, 2] * y + p.rotMatrix[3, 3] * z,
    )
end

@enum ParticleType begin
    UNKNOWN
    PION
    KAON
    PROTON
    ELECTRON
    MUON
end

struct ParticleProperty
    name::String
    mass::Float64
    pdgid::Int
end

const PARTICLE_PROPERTIES = Dict(
    ELECTRON => ParticleProperty("electron", ELECTRON_MASS, 11),
    MUON => ParticleProperty("muon", MUON_MASS, 13),
    PION => ParticleProperty("pion", PION_MASS, 211),
    KAON => ParticleProperty("kaon", KAON_MASS, 321),
    PROTON => ParticleProperty("proton", PROTON_MASS, 2212),
)

""" Returns the particle type and its properties based on the particle's pid.
    Iterates through PARTICLE_PROPERTIES to match the absolute pid value.
"""
function get_type_and_properties(p::Particle)
    for (k, v) in PARTICLE_PROPERTIES
        if v.pdgid == abs(p.pid)
            return k, v
        end
    end
    return UNKNOWN, ParticleProperties(0.0, 0)
end

""" Looks up and returns the mass of a particle given its pid.
    Returns 0.0 if the pid is not found.
"""
function get_particle_mass(pid::Int64)::Float64
    for (k, v) in PARTICLE_PROPERTIES
        if v.pdgid == abs(pid)
            return v.mass
        end
    end
    return 0.0
end
