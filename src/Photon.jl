"""
Struct representing a photon emitted by TORCH's radiator and propagated to the detector.

Position coordinates are in millimeters (mm).

# Fields
- `energy::Float64`: Photon energy.
- `xdir::Float64`: X-component of the photon direction in the lab frame.
- `ydir::Float64`: Y-component of the photon direction in the lab frame.
- `zdir::Float64`: Z-component of the photon direction in the lab frame.
- `xpos::Float64`: X-coordinate of the emission point.
- `ypos::Float64`: Y-coordinate of the emission point.
- `zpos::Float64`: Z-coordinate of the emission point.
- `slope::Float64`: Slope of the trajectory.
- `np::Float64`: Phase refractive index.
- `ng::Float64`: Group refractive index.
- `lambda::Float64`: Wavelength computed as LAMBDA/energy.
- `emissionTime::Float64`: Time at emission.
- `is_x_reflected::Bool`: Flag indicating reflection status in x-direction.
- `is_y_reflected::Bool`: Flag indicating reflection status in y-direction.
- `is_z_reflected::Bool`: Flag indicating reflection status in z-direction.
"""
struct Photon
    energy::Float64
    xdir::Float64
    ydir::Float64
    zdir::Float64
    xpos::Float64
    ypos::Float64
    zpos::Float64
    slope::Float64
    np::Float64
    ng::Float64
    lambda::Float64
    emissionTime::Float64
    is_x_reflected::Bool
    is_y_reflected::Bool
    is_z_reflected::Bool
end

""" 
Computes the direction of the photon emitted by the particle in the lab frame.

The photon direction is determined from the Cherenkov angle and azimuthal angle and 
rotated to the lab frame using the particle's direction.

# Arguments
- `p::Particle`: The particle instance.
- `costhetac::Float64`: Cosine of the Cherenkov angle.
- `sinthetac::Float64`: Sine of the Cherenkov angle.
- `phic::Float64`: Azimuthal angle around the particle direction.

# Returns
- `Tuple{Float64,Float64,Float64,Float64}`: A tuple containing (xdir, ydir, zdir, slope) of the photon direction.
"""
function _photon_direction(
    p::Particle,
    costhetac::Float64,
    sinthetac::Float64,
    phic::Float64,
)::Tuple{Float64,Float64,Float64,Float64}
    # photon direction
    xprime = sinthetac * cos(phic)
    yprime = sinthetac * sin(phic)
    zprime = costhetac

    # rotate to the lab frame
    xdir, ydir, zdir = rotate(p, xprime, yprime, zprime)
    slope = abs(zdir / xdir)

    return xdir, ydir, zdir, slope
end

""" 
Calculates the photon emission coordinates at a specified z depth in the radiator.

The emission point is determined by the particle's entry point and direction.
The given z-coordinate is fixed, while the x and y coordinates are adjusted based on the
particle's direction.

# Arguments
- `p::Particle`: The particle instance.
- `zemission::Float64`: Emission depth along the z-axis (mm).

# Returns
- `Tuple{Float64,Float64,Float64}`: A tuple containing (xpos, ypos, zpos) of the emission point.
"""
function _photon_emission(p::Particle, zemission::Float64)::Tuple{Float64,Float64,Float64}
    xpos = p.xCoord + zemission * p.xDir / p.zDir
    ypos = p.yCoord + zemission * p.yDir / p.zDir
    ypos += 0.5 * RADIATOR.height
    return xpos, ypos, zemission
end

""" 
Calculates the photon emission coordinates using the default radiator depth.

The z-coordinate is set to half the radiator depth.

# Arguments
- `p::Particle`: The particle instance.

# Returns
- `Tuple{Float64,Float64,Float64}`: A tuple containing (xpos, ypos, zpos) of the emission point.
"""
function _photon_emission(p::Particle)::Tuple{Float64,Float64,Float64}
    # emission point in z
    zemission = RADIATOR.depth * 0.5
    return _photon_emission(p, zemission)
end

"""
Generates a Photon instance with computed direction, emission, and reflection flags from
the particle and medium properties.

Cherenkov angle is calculated based on the particle's velocity and the medium's refractive index.
Phi angle is randomly generated (uniform distribution from 0 to 2π) to determine the azimuthal emission angle.
The emission point is randomly distributed along the radiator depth to simulate realistic photon production.
Reflection flags are determined by comparing angles with the critical angle for total internal reflection.
Surface roughness effects are not considered in this constructor but are handled separately.

# Arguments
- `p::Particle`: The particle instance emitting the photon.
- `beta::Float64`: The velocity factor of the particle (v/c).
- `nphase::Float64`: The phase refractive index of the medium for the specific photon energy.
- `ngroup::Float64`: The group refractive index of the medium for the specific photon energy.
- `energy::Float64`: The photon energy in appropriate units.

# Returns
- `Photon`: A new Photon object populated with calculated properties.
"""
function Photon(
    p::Particle,
    beta::Float64,
    nphase::Float64,
    ngroup::Float64,
    energy::Float64,
)
    # phi angle
    phi = 2.0 * π * rand()

    # cos and sinsqthetac
    costhetac = 1.0 / (beta * nphase)
    sinthetac = sqrt(1.0 - costhetac^2)

    # photon direction
    xdir, ydir, zdir, slope = _photon_direction(p, costhetac, sinthetac, phi)

    # emission point in z
    xpos, ypos, zpos = _photon_emission(p, RADIATOR.depth * rand())

    # critical angle
    sinsqthetacrit = (N_AIR / nphase)^2

    # check reflection
    is_x_reflected = ((1.0 - xdir * xdir) > sinsqthetacrit)
    is_y_reflected = ((1.0 - ydir * ydir) > sinsqthetacrit)
    is_z_reflected = ((1.0 - zdir * zdir) > sinsqthetacrit)

    Photon(
        energy,
        xdir,
        ydir,
        zdir,
        xpos,
        ypos,
        zpos,
        slope,
        nphase,
        ngroup,
        LAMBDA / energy,
        0.0,
        is_x_reflected,
        is_y_reflected,
        is_z_reflected,
    )
end

"""
Test the photon surface roughness in the z-direction.

Surface roughness is modeled as a Gaussian effect dependent on wavelength.
The probability of passing decreases with number of reflections.

# Arguments
- `photon::Photon`: The photon instance.
- `nreflec::Int`: The number of reflections.

# Returns
- `Bool`: Whether the photon passes the surface roughness condition (true) or not (false).
"""
function test_z_surface_roughness(photon::Photon, nreflec::Int)::Bool
    c = π * photon.np * ROUGHNESS / photon.lambda
    r = exp(-c * c)
    return rand() < r^nreflec
end

""" 
Checks whether the photon is within the focus acceptance region.

Acceptance is based on the slope of the photon trajectory.
Checks against minimum and maximum slope angles defined in FOCUS constants.

# Arguments
- `photon::Photon`: The photon instance.

# Returns
- `Bool`: Whether the photon is within the acceptance region (true) or not (false).

"""
function in_focus_acceptance(photon::Photon)::Bool
    # Check if the photon is within the acceptance region
    if photon.slope < FOCUS.tan_theta_min || photon.slope > FOCUS.tan_theta_max
        return false
    end
    return true
end
