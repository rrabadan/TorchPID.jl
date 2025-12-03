mutable struct PhotonHit
    localX::Float64
    localY::Float64
    localZ::Float64
    arrivalTime::Float64
    energy::Float64
    emissionX::Union{Float64,Nothing}
    emissionY::Union{Float64,Nothing}
    emissionZ::Union{Float64,Nothing}
    emissionTime::Union{Float64,Nothing}
    timeOfFlight::Union{Float64,Nothing}

    moduleId::Int32
    trackId::Int32

    isDownward::Bool
    isSurfaceScattered::Bool
    isRayleighScattered::Bool

    cherenkovTheta::Union{Float64,Nothing}
    cherenkovPhi::Union{Float64,Nothing}
    pathlengthInRadiator::Union{Float64,Nothing}
    GlobalX::Union{Float64,Nothing}
    GlobalY::Union{Float64,Nothing}
    GlobalZ::Union{Float64,Nothing}

    function PhotonHit(
        localX,
        localY,
        localZ,
        arrivalTime,
        energy,
        emissionX,
        emissionY,
        emissionZ,
        emissionTime,
        timeOfFlight,
        moduleId,
        trackId,
        isDownward,
        isSurfaceScattered,
        isRayleighScattered,
    )
        return new(
            localX,
            localY,
            localZ,
            arrivalTime,
            energy,
            emissionX,
            emissionY,
            emissionZ,
            emissionTime,
            timeOfFlight,
            moduleId,
            trackId,
            isDownward,
            isSurfaceScattered,
            isRayleighScattered,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
end


"""
    update_photon_hit!(hit::PhotonHit, hit_data; Debug::Bool=false)

Update the properties of a `PhotonHit` object `hit` using the provided `hit_data`.
An optional `Debug` flag can be set to `true` for additional debugging information.

# Arguments
- `hit::PhotonHit`: The photon hit object to be updated.
- `hit_data`: Data used to update the photon hit.
- `Debug::Bool=false`: Optional flag to enable debugging output.

# Returns
Nothing. The function modifies the `hit` object in place.
"""
function update_photon_hit!(hit::PhotonHit, hit_data; Debug::Bool = false)
    hit.localX = hit_data.LocalX
    hit.localY = hit_data.LocalY
    hit.localZ = hit_data.LocalZ
    hit.arrivalTime = hit_data.ArrivalTime
    hit.energy = hit_data.Energy
    hit.emissionX = hit_data.EmissionX
    hit.emissionY = hit_data.EmissionY
    hit.emissionZ = hit_data.EmissionZ
    hit.emissionTime = hit_data.EmissionTime
    hit.timeOfFlight = hit_data.TimeOfFlight
    hit.moduleId = hit_data.Module
    hit.trackId = hit_data.Track
    hit.isDownward = hit_data.IsDownward
    hit.isSurfaceScattered = hit_data.IsSurfaceScattered
    hit.isRayleighScattered = hit_data.IsRayleighScattered

    if Debug
        hit.cherenkovTheta = hit_data.CherenkovTheta
        hit.cherenkovPhi = hit_data.CherenkovPhi
        hit.pathlengthInRadiator = hit_data.PathlengthInRadiator
        hit.GlobalX = hit_data.GlobalX
        hit.GlobalY = hit_data.GlobalY
        hit.GlobalZ = hit_data.GlobalZ
    end

    return hit
end

"""
    photon_from_radiator(xemission::Float64, yemission::Float64, zemission::Float64)::Bool

`photon_from_focus` determines if a photon originates from the focus region.

# Arguments
- `xemission::Float64`: The x-coordinate of the photon's emission.
- `yemission::Float64`: The y-coordinate of the photon's emission.
- `zemission::Float64`: The z-coordinate of the photon's emission.

# Returns
- `Bool`: `true` if the photon originates from the radiator region, `false` otherwise.
"""
function photon_from_radiator(xemission, yemission, zemission, radiator::Radiator)::Bool
    return (
        abs(xemission) < radiator.half_width &&
        abs(yemission) < radiator.half_height &&
        zemission > 0 &&
        zemission < radiator.depth
    )
end

function photon_from_radiator(hit::PhotonHit, radiator::Radiator)::Bool
    return photon_from_radiator(hit.emissionX, hit.emissionY, hit.emissionZ, radiator)
end

"""
    photon_from_focus(yemission::Float64)::Bool

`photon_from_focus` determines if a photon originates from the focus region.

# Arguments
- `yemission::Float64`: The y-coordinate of the photon's emission.

# Returns
- `Bool`: `true` if the photon originates from the focus region, `false` otherwise.
"""
function photon_from_focus(yemission::Float64, radiator::Radiator, wedge::Wedge)::Bool
    radiator_top = 0.5 * radiator.height - wedge.offset
    return yemission > radiator_top
end

function photon_from_focus(hit::PhotonHit, radiator::Radiator, wedge::Wedge)::Bool
    return photon_from_focus(hit.emissionY, radiator, wedge)
end

"""
    photon_on_detector(x::Float64, y::Float64)::Bool

`photon_on_detector` checks if a photon hits the detector.

# Arguments
- `x::Float64`: The x-coordinate of the photon.
- `y::Float64`: The y-coordinate of the photon.

# Returns
- `Bool`: `true` if the photon hits the detector, `false` otherwise.
"""
function photon_on_detector(x::Float64, y::Float64, detector::Detector)::Bool
    return (
        y > detector.y_min &&
        y < detector.y_max &&
        x > detector.x_min &&
        x < detector.x_max
    )
end

function photon_on_detector(hit::PhotonHit, detector::Detector)::Bool
    return photon_on_detector(hit.localX, hit.localY, detector)
end


"""
    is_photon_masked(y::Float64)::Bool

`is_photon_masked` checks if a photon is masked based on its y-coordinate.

# Arguments
- `y::Float64`: The y-coordinate of the photon.

# Returns
- `Bool`: `true` if the photon is masked, `false` otherwise.
"""
function is_photon_masked(y::Float64, mask::Mask)::Bool
    return mask.y_min < y < mask.y_max
end

function is_photon_masked(hit::PhotonHit, mask::Mask)::Bool
    return is_photon_masked(hit.localY, mask::Mask)
end
