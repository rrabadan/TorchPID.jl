"""
    PhotonMapper(blackened_sides blackened_bottom,
                 blackened, surface_roughness,
                 max_x_reflections)

Type representing the photon mapping configuration for photon propagation.

# Fields
- `blackened_sides::Bool`: If true, photons reflecting off the radiator's sides are absorbed.
- `blackened_bottom::Bool`: If true, photons reflecting off the radiator's bottom are ignored.
- `blackened_focus::Bool`: If true, reflections in the focus region are handled differently.
- `surface_roughness::Bool`: If true, surface imperfections of the radiator are considered.
- `max_x_reflections::INt`: Specifies the maximum allowable x-direction reflections before the photon is discarded.

# Constructors

    PhotonMapper(; 
        blackened_sides::Bool=false,
        blackened_bottom::Bool=false,
        blackened_focus::Bool=false,
        surface_roughness::Bool=false,
        max_x_reflections::Int=10
    )

## Keywords
- `blackened_sides::Bool=false`: Enable absorption of side reflections.
- `blackened_bottom::Bool=false`: Ignore bottom reflections if true.
- `blackened_focus::Bool=false`: Handle focus region reflections differently.
- `surface_roughness::Bool=false`: Account for surface imperfections.
- `max_x_reflections::Int=10`: Maximum number of x-direction reflections allowed.

## Returns
- `PhotonMapper`: A new PhotonMapper instance with the specified configuration.
"""
struct PhotonMapper
    blackened_sides::Bool
    blackened_bottom::Bool
    blackened_focus::Bool
    surface_roughness::Bool
    max_x_reflections::Int
end

function PhotonMapper(;
    blackened_sides::Bool = false,
    blackened_bottom::Bool = false,
    blackened_focus::Bool = false,
    surface_roughness::Bool = false,
    max_x_reflections::Int = 10,
)
    PhotonMapper(
        blackened_sides,
        blackened_bottom,
        blackened_focus,
        surface_roughness,
        max_x_reflections,
    )
end

"""
    trace_photon(mapper::PhotonMapper, photon::Photon, t0::Float64)

`trace_photon` simulates the photon's journey through the optical system, from emission to detection. 
It propagates the photon to the top of the plate, accounting for y-direction and z-direction reflections, surface roughness, and pathlength accumulation.
It then traces the photon to the focus block, simulating its interaction with the focusing optics, and finally to the detector plane, where it checks for acceptance criteria.
The function returns the detected hit coordinate if the photon is successfully traced, or `nothing` if any tracing step fails.

# Arguments
- `mapper::PhotonMapper`: The photon mapping configuration.
- `photon::Photon`: The photon to trace.
- `t0::Float64`: The production time of the photon.

# Returns
- `Union{HitCoordinate,Nothing}`: The detected hit coordinate if the photon is successfully traced,
  or `nothing` if any tracing step fails.
"""
function trace_photon(
    mapper::PhotonMapper,
    photon::Photon,
    t0::Float64,
)::Union{HitCoordinate,Nothing}

    # check if photon reflects in z
    if !photon.is_z_reflected || abs(photon.ydir) < 1e-5
        return nothing
    end

    if in_focus_acceptance(photon)
        return nothing
    end

    tot_pathlength = 0.0

    # trace to top of the plate
    result = _trace_to_top_of_plate(mapper, photon)
    if isnothing(result)
        return nothing
    end

    zp, typ, tzp, pathlength = result
    tot_pathlength += pathlength

    # trace to the focus block
    result = _trace_to_focus(zp, typ, tzp)
    if isnothing(result)
        return nothing
    end

    yf, tyf, tzf, pathlength = result
    tot_pathlength += pathlength

    # trace to the detector
    result = _trace_to_detector(yf, tyf, tzf)
    if isnothing(result)
        return nothing
    end

    ydetected, pathlength = result
    tot_pathlength += pathlength

    # project x-detector position
    xdetected = _project_x_detector_position(mapper, photon, tot_pathlength)
    if isnothing(xdetected)
        return nothing
    end

    # detected time is production time + time of propagation
    tdetected = t0 + (pathlength * photon.ng / CLIGHT)

    HitCoordinate(xdetected, ydetected, tdetected)

end

"""
    _trace_to_top_of_plate(mapper::PhotonMapper, photon::Photon)

`_trace_to_top_of_plate` propagates the photon to the top of the plate, accounting for y-direction and z-direction reflections, surface roughness, and pathlength accumulation.

# Arguments
- `mapper::PhotonMapper`: The mapping configuration.
- `photon::Photon`: The photon to be traced.

# Returns
- `Union{Tuple{Float64,Float64,Float64,Float64},Nothing}`: A tuple containing:
  - `zp::Float64`: Position in z
  - `typ::Float64`: Y-component of photon momentum
  - `tzp::Float64`: Z-component of photon momentum
  - `pathlength::Float64`: Accumulated pathlength
  
  Returns `nothing` if the photon's path is terminated.
"""
function _trace_to_top_of_plate(
    mapper::PhotonMapper,
    photon::Photon,
)::Union{Tuple{Float64,Float64,Float64,Float64},Nothing}
    ydistance = 0.0
    if photon.ydir > 0.0
        # photon is moving up
        ydistance = RADIATOR.height - WEDGE.offset - photon.ypos
        typ = photon.ydir
    else
        if mapper.blackened_bottom || !photon.is_y_reflected
            return nothing
        end
        ydistance = RADIATOR.height - WEDGE.offset + photon.ypos
        typ = -photon.ydir
    end

    # travelled diztance in z
    zdistance = photon.zpos + ydistance * photon.slope

    # number of reflections in zDir
    nz = Int(floor(zdistance / RADIATOR.depth))
    if mapper.surface_roughness && !test_z_surface_roughness(photon, nz)
        return nothing
    end

    zp = zdistance - nz * RADIATOR.depth

    nz_odd = nz % 2 == 1
    if nz_odd
        zp = RADIATOR.depth - zp
    end

    if zp < 0.0 || zp > RADIATOR.depth
        return nothing
    end

    tzp = (!nz_odd) ? photon.zdir : -photon.zdir

    pathlength = ydistance / typ

    return zp, typ, tzp, pathlength
end

"""
    _trace_to_focus(zp::Float64, typ::Float64, tzp::Float64)

`_trace_to_focus` propagates the photon from the top of the flat section to the entrance of the focus block.

# Arguments
- `zp::Float64`: Position at top of the plate
- `typ::Float64`: Y-component of photon momentum at top of the plate
- `tzp::Float64`: Z-component of photon momentum at top of the plate

# Returns
- `Union{Tuple{Float64,Float64,Float64,Float64},Nothing}`: A tuple containing:
  - `yf::Float64`: Position in focus block coordinates
  - `tyf::Float64`: Y-component of photon momentum in focus block coordinates
  - `tzf::Float64`: Z-component of photon momentum in focus block coordinates
  - `pathlength::Float64`: Accumulated pathlength
  
  Returns `nothing` if the photon's path is terminated.
"""
function _trace_to_focus(
    zp::Float64,
    typ::Float64,
    tzp::Float64,
)::Union{Tuple{Float64,Float64,Float64,Float64},Nothing}
    tanthetap = abs(typ / tzp)
    alpha = (tzp < 0.0) ? (RADIATOR.depth + zp) : (RADIATOR.depth - zp)

    # Use sine rule to calculate position on focus block
    yf =
        alpha * tanthetap * WEDGE.tan_theta /
        (WEDGE.sin_theta * (1.0 + tanthetap * WEDGE.tan_theta))

    # Use distance travelled in y to compute the path
    pathlength = yf * WEDGE.cos_theta / typ

    # Calculate position relative to axis
    yf = FOCUS.axis_to_wedge - yf

    # Rotate the photon to the focus block
    tyf = abs(tzp) * FOCUS.cos_theta - typ * FOCUS.sin_theta
    tzf = abs(tzp) * FOCUS.sin_theta + typ * FOCUS.cos_theta

    return yf, tyf, tzf, pathlength
end

"""
    _trace_to_detector(yf::Float64, tyf::Float64, tzf::Float64)

`_trace_to_detector` simulates the photon's propagation from the entrance of the focusing block to the detector plane, accounting for mirror reflections and acceptance criteria.


# Arguments
- `yf::Float64`: Position at entrance of the focus block
- `tyf::Float64`: Y-component of photon momentum at entrance of the focus block
- `tzf::Float64`: Z-component of photon momentum at entrance of the focus block

# Returns
- `Union{Tuple{Float64,Float64},Nothing}`: A tuple containing:
  - `ydetected::Float64`: Position of photon on detector plane
  - `pathlength::Float64`: Accumulated pathlength
  
  Returns `nothing` if the photon's path is terminated.
"""
function _trace_to_detector(
    yf::Float64,
    tyf::Float64,
    tzf::Float64,
)::Union{Tuple{Float64,Float64},Nothing}
    zf::Float64 = FOCUS.entrance

    # Compute mirror position
    # Solve yf + k tyf = r_c \cos\theta_m
    #       zf + k tzf = r_c \sin\theta_m
    # where r_c is the radius of curvature of the mirror
    # and \theta_m is the mirror angle

    # Normalize input direction components
    norm = sqrt(tyf^2 + tzf^2)
    tyfp = tyf / norm
    tzfp = tzf / norm

    # Compute coefficients for quadratic in k
    b = 2.0 * (tzfp * zf + tyfp * yf)
    c = zf^2 + yf^2 - FOCUS.rsquared

    # Solve for k (assuming the positive solution is desired)
    k = 0.5 * (-b + sqrt(b^2 - 4.0 * c))

    # Calculate the position on the mirror surface
    ys = yf + k * tyfp
    zs = zf + k * tzfp

    # Reject if the mirror hit lies outside the allowed axis limits
    if ys > FOCUS.axis_to_rhs || ys < -FOCUS.axis_to_lhs
        return nothing
    end

    # Compute the mirror angle components
    costhetam = ys / FOCUS.radius
    sinthetam = zs / FOCUS.radius

    # Compute the outgoing photon direction using reflection formula
    # (see https://en.wikipedia.org/wiki/Specular_reflection)
    # \vec{p}_{\text{out}} = \vec{p}_{\text{in}} - 2 (\vec{p}_{\text{in}} \cdot \hat{n}) \hat{n}
    # the paralle component of the photon momentum is conserved, the perpendicular component is reversed
    # \hat{n} = (-\cos(\theta_m), \sin(\theta_m))
    typ_out = tyf - 2.0 * costhetam * (costhetam * tyf + sinthetam * tzf)
    tzp_out = tzf - 2.0 * sinthetam * (costhetam * tyf + sinthetam * tzf)

    # Compute the z-distance from mirror to detector plane
    deltazf = zf - zs

    # Compute the y-position on the detector plane
    ydetected = ys + deltazf * typ_out / tzp_out

    # Shift to detector coordinates
    ydetected = -(ydetected + FOCUS.axis_to_detector)

    # Reject photons outside the acceptance region
    if abs(ydetected) > DETECTOR.y_max
        return nothing
    end

    # Reject photons in the image mask range
    if ydetected > MASK.y_min && ydetected < MASK.y_max
        return nothing
    end

    # Update the pathlength (accumulating the additional distance)
    pathlength = abs(deltazf / tzp_out) + abs(deltazf / tzf)

    return ydetected, pathlength
end

"""
    _project_x_detector_position(mapper::PhotonMapper, photon::Photon, pathlength::Float64)

`_project_x_detector_position` calculates the photon's x-coordinate on the detector plane by projecting its trajectory.
It considers x-direction reflections based on the mapper's configuration.
Returns `nothing` if side blackening is enabled and reflection is required, the photon cannot reflect in x, or the maximum allowed reflections is exceeded.

# Arguments
- `mapper::PhotonMapper`: The mapping configuration.
- `photon::Photon`: The photon instance.
- `pathlength::Float64`: Total travelled pathlength.

# Returns
- `Union{Float64,Nothing}`: The projected x position on the detector, or `nothing` if the
  photon does not reach the detector due to reflection conditions.
"""
function _project_x_detector_position(
    mapper::PhotonMapper,
    photon::Photon,
    pathlength::Float64,
)::Union{Float64,Nothing}
    # Compute distance travelled in x
    xdetected = photon.xpos + pathlength * photon.xdir

    # Counter for reflections in x
    nx = 0

    # Fold the pattern in x until the condition is met
    while abs(xdetected) > RADIATOR.half_width
        # If blackening is enabled, the photon does not reflect further,
        # or photon was not reflected in x, or we already exceeded max allowed reflections:
        if mapper.blackened_sides || !photon.is_x_reflected || nx > mapper.max_x_reflections
            return nothing
        end

        # Derive the x-position at the detector plane:
        if xdetected > DETECTOR.x_max
            xdetected = RADIATOR.width - xdetected
        else
            xdetected = -RADIATOR.width - xdetected
        end

        nx += 1
    end
    return xdetected
end
