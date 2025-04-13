"""
Struct representing a particle with kinematic and identification fields.

"True" values refer to the original, undistorted values, while "reco" values are the reconstructed, measured values.
The units of position, momentum, and time depend on the specific detector system.
The rotation matrix is used to transform coordinates from the particle's local reference frame to the global 
reference frame.

# Fields
- `pid::Int`: Particle identification code (e.g., PDG ID).
- `eventId::Int`: Identifier for the event the particle belongs to.
- `trackId::Int`: Identifier for the particle's reconstructed track.
- `moduleId::Int`: Identifier for the detector module where the particle was detected.
- `xCoord::Float64`: x-coordinate of the particle's position.
- `yCoord::Float64`: y-coordinate of the particle's position.
- `pMag::Float64`: Magnitude of the particle's momentum.
- `pMagTrue::Float64`: True magnitude of the particle's momentum.
- `xDir::Float64`: x-component of the particle's direction vector.
- `yDir::Float64`: y-component of the particle's direction vector.
- `zDir::Float64`: z-component of the particle's direction vector.
- `truePX::Float64`: True x-component of the particle's momentum.
- `truePY::Float64`: True y-component of the particle's momentum.
- `truePZ::Float64`: True z-component of the particle's momentum.
- `recoPX::Float64`: Reconstructed x-component of the particle's momentum.
- `recoPY::Float64`: Reconstructed y-component of the particle's momentum.
- `recoPZ::Float64`: Reconstructed z-component of the particle's momentum.
- `pathlength::Float64`: Reconstructed path length of the particle's track.
- `truePathlength::Float64`: True path length of the particle's track.
- `t0::Float64`: Time of the particle's creation or detection.
- `rotMatrix::Array{Float64,2}`: Rotation matrix associated with the particle.
"""
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

"""
Constructs a `Particle` object with specified or default values.

Default values are provided for all arguments, allowing for partial initialization.
The default rotation matrix is the identity matrix, representing no rotation.

# Arguments
- `pid::Int=211`: Particle identification code (default: 211, pion).
- `eventId::Int=0`: Identifier for the event the particle belongs to.
- `trackId::Int=0`: Identifier for the particle's reconstructed track.
- `moduleId::Int=0`: Identifier for the detector module.
- `xCoord::Float64=0.0`: X-coordinate of the particle's position.
- `yCoord::Float64=0.0`: Y-coordinate of the particle's position.
- `pMag::Float64=0.0`: Magnitude of the particle's momentum.
- `pMagTrue::Float64=0.0`: True magnitude of the particle's momentum.
- `xDir::Float64=0.0`: X-component of the particle's direction vector.
- `yDir::Float64=0.0`: Y-component of the particle's direction vector.
- `zDir::Float64=0.0`: Z-component of the particle's direction vector.
- `truePX::Float64=0.0`: True x-component of the particle's momentum.
- `truePY::Float64=0.0`: True y-component of the particle's momentum.
- `truePZ::Float64=0.0`: True z-component of the particle's momentum.
- `recoPX::Float64=0.0`: Reconstructed x-component of the particle's momentum.
- `recoPY::Float64=0.0`: Reconstructed y-component of the particle's momentum.
- `recoPZ::Float64=0.0`: Reconstructed z-component of the particle's momentum.
- `pathlength::Float64=0.0`: Reconstructed path length of the particle's track.
- `truePathlength::Float64=0.0`: True path length of the particle's track.
- `t0::Float64=0.0`: Time of the particle's creation or detection.
- `rotMatrix::Array{Float64,2}=Array{Float64,2}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])`: Rotation matrix.

# Returns
- `Particle`: A new Particle instance with the specified properties.
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

""" 
Calculates the beta factor for the particle with a given mass hypothesis.

Beta is calculated as β = p / √(p² + m²), where p is momentum magnitude.
This factor represents the ratio of the particle's velocity to the speed of light.

# Arguments
- `p::Particle`: The particle instance.
- `m::Float64`: The mass hypothesis in appropriate units.

# Returns
- `Float64`: The relativistic beta factor (v/c) of the particle.
"""
function beta(p::Particle, m::Float64)::Float64
    return p.pMag / sqrt(p.pMag^2 + m^2)
end

"""
Calculates the gamma factor for the particle with a given mass hypothesis.

Gamma is calculated as γ = 1 / √(1 - β²), where β is the beta factor.
Internally calls the beta() function to calculate the beta factor.

# Arguments
- `p::Particle`: The particle instance.
- `m::Float64`: The mass hypothesis in appropriate units.

# Returns
- `Float64`: The relativistic gamma factor of the particle.
"""
function gamma(p::Particle, m::Float64)::Float64
    return 1 / sqrt(1 - beta(p, m)^2)
end

""" 
Initializes the rotation matrix of the particle based on its direction components.

The rotation matrix transforms coordinates from the particle's local reference frame to the global reference frame.
The calculation depends on the particle's direction components (xDir, yDir, zDir).
If the perpendicular component (uperp) is zero, the rotation matrix remains unchanged.

# Arguments
- `p::Particle`: The particle instance whose rotation matrix will be initialized.
"""
function initRotation!(p::Particle)
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

"""
Rotates a 3D vector from the particle's local reference frame to the global reference frame.

Uses the particle's rotation matrix to perform the coordinate transformation.
The rotation preserves the vector's magnitude.

# Arguments
- `p::Particle`: The particle instance containing the rotation matrix.
- `x::Float64`: The x-component of the vector in the local reference frame.
- `y::Float64`: The y-component of the vector in the local reference frame.
- `z::Float64`: The z-component of the vector in the local reference frame.

# Returns
- `Tuple{Float64, Float64, Float64}`: The rotated vector components in the global reference frame.
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

"""
An enumeration representing different particle types.

Used to categorize particles for identification and physics calculations.
Each enum value corresponds to a specific particle with well-defined properties.

# Values
- `UNKNOWN`: An unknown or undefined particle type.
- `PION`: A pion particle.
- `KAON`: A kaon particle.
- `PROTON`: A proton particle.
- `ELECTRON`: An electron particle.
- `MUON`: A muon particle.
"""
@enum ParticleType begin
    UNKNOWN
    PION
    KAON
    PROTON
    ELECTRON
    MUON
end

"""
Struct representing the properties of a particle type.

The PDG ID is a standardized identifier used in particle physics.
Mass values are constants defined elsewhere in the code.

# Fields
- `name::String`: The name of the particle (e.g., "electron", "pion").
- `mass::Float64`: The mass of the particle in appropriate units.
- `pdgid::Int`: The Particle Data Group ID (PDG ID) for the particle.
"""
struct ParticleProperty
    name::String
    mass::Float64
    pdgid::Int
end

"""
Dictionary mapping `ParticleType` values to their corresponding `ParticleProperty`.

Provides a convenient way to access particle properties by type.
Contains entries for standard particles: electron, muon, pion, kaon, and proton.
Each entry maps a ParticleType enum value to its corresponding ParticleProperty instance.
"""
const PARTICLE_PROPERTIES = Dict(
    ELECTRON => ParticleProperty("electron", ELECTRON_MASS, 11),
    MUON => ParticleProperty("muon", MUON_MASS, 13),
    PION => ParticleProperty("pion", PION_MASS, 211),
    KAON => ParticleProperty("kaon", KAON_MASS, 321),
    PROTON => ParticleProperty("proton", PROTON_MASS, 2212),
)

""" 
Determines the particle type and properties based on the particle's PDG ID.

Searches the PARTICLE_PROPERTIES dictionary for a matching PDG ID.
Returns UNKNOWN type with null properties if no match is found.
The absolute value of the PDG ID is used for the comparison to handle antiparticles.

# Arguments
- `p::Particle`: The particle instance whose type and properties are to be determined.

# Returns
- `Tuple{ParticleType, ParticleProperty}`: A tuple containing the particle type and its properties.
"""
function get_type_and_properties(p::Particle)
    for (k, v) in PARTICLE_PROPERTIES
        if v.pdgid == abs(p.pid)
            return k, v
        end
    end
    return UNKNOWN, ParticleProperties(0.0, 0)
end

""" 
Looks up the mass of a particle based on its PDG ID.

Searches the PARTICLE_PROPERTIES dictionary for a matching PDG ID.
The absolute value of the PDG ID is used for the comparison to handle antiparticles.
Returns 0.0 if no matching particle type is found.

# Arguments
- `pid::Int64`: The Particle Data Group ID (PDG ID).

# Returns
- `Float64`: The mass of the particle in appropriate units, or 0.0 if not found.
"""
function get_particle_mass(pid::Int64)::Float64
    for (k, v) in PARTICLE_PROPERTIES
        if v.pdgid == abs(pid)
            return v.mass
        end
    end
    return 0.0
end
