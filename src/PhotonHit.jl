@inline _to_f64(::Nothing) = NaN
@inline _to_f64(x::Real) = Float64(x)

struct PhotonHit
    localX::Float64
    localY::Float64
    localZ::Float64
    arrivalTime::Float64
    energy::Float64
    emissionX::Float64
    emissionY::Float64
    emissionZ::Float64
    emissionTime::Float64
    timeOfFlight::Float64

    moduleId::Int32
    trackId::Int32

    isDownward::Bool
    isSurfaceScattered::Bool
    isRayleighScattered::Bool

    cherenkovTheta::Float64
    cherenkovPhi::Float64
    pathlengthInRadiator::Float64
    GlobalX::Float64
    GlobalY::Float64
    GlobalZ::Float64

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
        cherenkovTheta = NaN,
        cherenkovPhi = NaN,
        pathlengthInRadiator = NaN,
        GlobalX = NaN,
        GlobalY = NaN,
        GlobalZ = NaN,
    )
        return new(
            Float64(localX),
            Float64(localY),
            Float64(localZ),
            Float64(arrivalTime),
            Float64(energy),
            _to_f64(emissionX),
            _to_f64(emissionY),
            _to_f64(emissionZ),
            _to_f64(emissionTime),
            _to_f64(timeOfFlight),
            Int32(moduleId),
            Int32(trackId),
            Bool(isDownward),
            Bool(isSurfaceScattered),
            Bool(isRayleighScattered),
            _to_f64(cherenkovTheta),
            _to_f64(cherenkovPhi),
            _to_f64(pathlengthInRadiator),
            _to_f64(GlobalX),
            _to_f64(GlobalY),
            _to_f64(GlobalZ),
        )
    end
end


"""
    update_photon_hit!(hit::PhotonHit, hit_data; Debug::Bool=false)

Create an updated `PhotonHit` object using the provided `hit_data`.
An optional `Debug` flag can be set to `true` for additional debugging information.
"""
function update_photon_hit!(::PhotonHit, hit_data; Debug::Bool = false)
    return PhotonHit(
        hit_data.LocalX,
        hit_data.LocalY,
        hit_data.LocalZ,
        hit_data.ArrivalTime,
        hit_data.Energy,
        hit_data.EmissionX,
        hit_data.EmissionY,
        hit_data.EmissionZ,
        hit_data.EmissionTime,
        hit_data.TimeOfFlight,
        hit_data.Module,
        hit_data.Track,
        hit_data.IsDownward,
        hit_data.IsSurfaceScattered,
        hit_data.IsRayleighScattered,
        Debug ? hit_data.CherenkovTheta : NaN,
        Debug ? hit_data.CherenkovPhi : NaN,
        Debug ? hit_data.PathlengthInRadiator : NaN,
        Debug ? hit_data.GlobalX : NaN,
        Debug ? hit_data.GlobalY : NaN,
        Debug ? hit_data.GlobalZ : NaN,
    )
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
