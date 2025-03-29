""" 
Represents the photon mapping configuration for photon propagation.
Fields:
- blackened_sides: Enable blackening of side reflections.
- blackened_bottom: Ignore bottom reflections if true.
- blackened_focus: Apply focus reflection handling.
- surface_roughness: Consider surface roughness effects.
- max_x_reflections: Maximum number of x-direction reflections.
"""
struct PhotonMapper
    blackened_sides::Bool
    blackened_bottom::Bool
    blackened_focus::Bool
    surface_roughness::Bool
    max_x_reflections::Int
end

""" 
Construct a PhotonMapper with configurable options.
Keyword Arguments:
- blackened_sides::Bool (default: false)
- blackened_bottom::Bool (default: false)
- blackened_focus::Bool (default: false)
- surface_roughness::Bool (default: false)
- max_x_reflections::Int (default: 10)
Returns a new PhotonMapper instance.
"""
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
Traces the photon through the entire optical system from emission to detection.
Arguments:
- mapper::PhotonMapper: The photon mapping configuration.
- photon::Photon: The photon to trace.
- t0::Float64: The production time of the photon.
Returns:
- HitCoordinate: The detected hit coordinate if the photon is successfully traced.
- Nothing: if any tracing step fails.
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
Propagates the photon to the top of the plate.
Arguments:
- mapper::PhotonMapper: The mapping configuration.
- photon::Photon: The photon to be traced.
Returns:
- Tuple{Float64,Float64,Float64,Float64}: (zp, typ, tzp, pathlength) if the photon reaches the top.
    zp -- position in z
    typ -- y-componnet of photon momentum
    tzp -- z-componnet of photon momentum
    pathlength -- accumulated pathlength
- Nothing: if the photon's path is terminated.
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
Propagates the photon from the top of the flat section to the entrance of the focus block.
Arguments:
- zp -- position at top of the plate
- typ -- y slope at top of the plate
- tzp -- z slope at top of the plate
Returns:
- Tuple{Float64,Float64,Float64,Float64}: (yf, tyf, tzf, pathlength) if the photon reaches the top.
    yf -- position in focus block coordinates
    tyf -- y-componnet of photon momentum in focus block coordinates
    tzf -- z-componnet of photon momentum in focus block coordinates
    pathlength -- accumulated pathlength
- Nothing: if the photon's path is terminated.
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
Propagates the photon from the entrance of the focussing block to the MCPs
Arguments:
- yf -- position at entrance of the focus block
- tyf -- y-slope at entrance of the focus block
- tzp -- z-slope at entrance of the focus block
Returns:
- Tuple{Float64,Float64}: (ydetected,pathlength) if the photon reaches the MCPs.
    ydetected -- position of photon on detector plane
    pathlength -- accumulated pathlength
- Nothing: if the photon's path is terminated.
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
Calculates the x position on the detector by projecting the photon’s x coordinate
using the provided pathlength and accounting for reflections.
Arguments:
- mapper::PhotonMapper: The mapping configuration.
- photon::Photon: The photon instance.
- pathlength::Float64: Total travelled pathlength.
Returns:
- Float64: The projected x position or nothing if conditions are not met.
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
