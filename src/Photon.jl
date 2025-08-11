"""
    Photon(energy, xdir, ydir, zdir, xpos, ypos, zpos, slope,
           np, ng, lambda, emissionTime, is_x_reflected,
           is_y_reflected, is_z_reflected)

Type representing a photon emitted by TORCH's radiator and propagated to the detector.
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
    PhotonContext(radiator, constants)

Type holding constants for creating photons with common properties.

# Fields
- `radiator::Radiator`: The radiator material properties and dimensions.
- `constants::Constants`: The physical constants used in photon calculations.
"""
struct PhotonContext
    radiator::Radiator
    constants::Constants
end

"""
    create_random_photon(particle, beta, nphase, ngroup, energy, context)

Creates a Photon instance by determining its direction, emission point, and reflection flags based on the particle's properties, photon energy, and refractive indices. 
The Cherenkov angle is calculated using the particle's velocity and the phase refractive index. 
A random azimuthal angle (phi) is uniformly sampled between 0 and 2π to simulate isotropic emission. 
The emission point is randomly chosen along the radiator depth to reflect realistic photon production. 
Reflection flags are set by comparing the photon's direction with the critical angle for total internal reflection. 
Surface roughness effects are excluded from this constructor and are handled separately.

## Arguments
- `particle::Particle`: The particle instance emitting the photon.
- `beta::Float64`: The velocity factor of the particle (v/c).
- `nphase::Float64`: The phase refractive index for the specific photon energy.
- `ngroup::Float64`: The group refractive index for the specific photon energy.
- `energy::Float64`: The photon energy in appropriate units.
- `context::PhotonContext`: The constants object containing radiator and refractive index properties.

## Returns
- `Photon`: A new Photon object populated with calculated properties.
"""
function create_random_photon(
    particle::Particle,
    beta::Float64,
    nphase::Float64,
    ngroup::Float64,
    energy::Float64,
    context::PhotonContext,
)
    # phi angle
    phi = 2.0 * π * rand()

    # cos and sinsqthetac
    costhetac = 1.0 / (beta * nphase)
    sinthetac = sqrt(1.0 - costhetac^2)

    # photon direction
    xdir, ydir, zdir, slope = _photon_direction(particle, costhetac, sinthetac, phi)

    # emission point in z
    xpos, ypos, zpos = _photon_emission(particle, context.radiator.depth * rand(), context)

    # critical angle
    sinsqthetacrit = (context.constants.N_AIR / nphase)^2

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
        context.constants.LAMBDA / energy,
        0.0,
        is_x_reflected,
        is_y_reflected,
        is_z_reflected,
    )
end


"""
    create_explicit_photon(particle, beta, phic, energy, nphase, ngroup, context)

Creates a Photon instance with a specific azimuthal angle instead of random sampling.
The Cherenkov angle is calculated from beta and the phase refractive index.
The emission point is determined at the default location in the radiator.

# Arguments
- `particle::Particle`: The particle instance emitting the photon.
- `beta::Float64`: The velocity factor of the particle (v/c).
- `phic::Float64`: The specified azimuthal angle for the photon emission.
- `energy::Float64`: The photon energy in appropriate units.
- `nphase::Float64`: The phase refractive index for the specific photon energy.
- `ngroup::Float64`: The group refractive index for the specific photon energy.
- `context::PhotonContext`: The constants object containing radiator and refractive index properties.

# Returns
- `Photon`: A new Photon object populated with calculated properties.
"""
function create_explicit_photon(
    particle::Particle,
    beta::Float64,
    phic::Float64,
    energy::Float64,
    nphase::Float64,
    ngroup::Float64,
    context::PhotonContext,
)::Photon

    # emission point in z
    xpos, ypos, zpos = _photon_emission(particle, context)

    # cos and sinsqthetac
    costhetac = 1.0 / (beta * nphase)
    sinthetac = sqrt(1.0 - costhetac^2)

    # photon direction
    xdir, ydir, zdir, slope = _photon_direction(p, costhetac, sinthetac, phic)

    # critical angle
    sinsqthetacrit = (context.constants.N_AIR / nphase)^2

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
        context.constants.LAMBDA / energy,
        0.0,
        is_x_reflected,
        is_y_reflected,
        is_z_reflected,
    )
end

"""
    create_explicit_photon(particle, thetac, phic, energy, nphase, ngroup, xemission, yemission, zemission, context)

Creates a Photon instance with precisely specified Cherenkov angle, azimuthal angle, and emission point coordinates.
This function provides maximum control over the photon properties for simulation purposes.

# Arguments
- `particle::Particle`: The particle instance emitting the photon.
- `thetac::Float64`: The specified Cherenkov angle for the photon emission.
- `phic::Float64`: The specified azimuthal angle for the photon emission.
- `energy::Float64`: The photon energy in appropriate units.
- `nphase::Float64`: The phase refractive index for the specific photon energy.
- `ngroup::Float64`: The group refractive index for the specific photon energy.
- `xemission::Float64`: The x-coordinate of the emission point.
- `yemission::Float64`: The y-coordinate of the emission point.
- `zemission::Float64`: The z-coordinate of the emission point (depth in radiator).
- `context::PhotonContext`: The constants object containing radiator and refractive index properties.

# Returns
- `Photon`: A new Photon object populated with calculated properties.
"""
function create_explicit_photon(
    particle::Particle,
    thetac::Float64,
    phic::Float64,
    energy::Float64,
    nphase::Float64,
    ngroup::Float64,
    xemission::Float64,
    yemission::Float64,
    zemission::Float64,
    context::PhotonContext,
)::Photon

    # photon direction
    xdir, ydir, zdir, slope = _photon_direction(particle, cos(thetac), sin(thetac), phic)

    # critical angle
    sinsqthetacrit = (context.constants.N_AIR / nphase)^2

    # check reflection
    is_x_reflected = ((1.0 - xdir * xdir) > sinsqthetacrit)
    is_y_reflected = ((1.0 - ydir * ydir) > sinsqthetacrit)
    is_z_reflected = ((1.0 - zdir * zdir) > sinsqthetacrit)

    Photon(
        energy,
        xdir,
        ydir,
        zdir,
        xemission,
        yemission + 0.5 * context.radiator.height,
        zemission,
        slope,
        nphase,
        ngroup,
        context.constants.LAMBDA / energy,
        0.0,
        is_x_reflected,
        is_y_reflected,
        is_z_reflected,
    )
end


function create_approximate_photon(
    particle::Particle,
    beta::Float64,
    nphase::Float64,
    ngroup::Float64,
    energy::Float64,
    zemission::Float64,
    phi::Float64,
    context::PhotonContext,
)::Photon

    # cos and sinsqthetac
    costhetac = 1.0 / (beta * nphase)
    sinthetac = sqrt(1.0 - costhetac^2)

    # photon direction
    xdir, ydir, zdir, slope = _photon_direction(particle, costhetac, sinthetac, phi)

    # emission point in z
    xpos, ypos, zpos = _photon_emission(particle, zemission, context)

    # critical angle
    sinsqthetacrit = (context.constants.N_AIR / nphase)^2

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
        context.constants.LAMBDA / energy,
        0.0,
        is_x_reflected,
        is_y_reflected,
        is_z_reflected,
    )
end

"""
    test_z_surface_roughness(photon::Photon, roughness::Float64, nreflec::Int)::Bool

Tests the photon surface roughness in the z-direction. 
Surface roughness is modeled as a Gaussian effect dependent on wavelength.
The probability of passing decreases with number of reflections.

# Arguments
- `photon::Photon`: The photon instance.
- `roughness::Float64`: The surface roughness parameter in appropriate units.
- `nreflec::Int`: The number of reflections.

# Returns
- `Bool`: Whether the photon passes the surface roughness condition (true) or not (false).
"""
function test_z_surface_roughness(photon::Photon, roughness::Float64, nreflec::Int)::Bool
    c = π * photon.np * roughness / photon.lambda
    r = exp(-c * c)
    return rand() < r^nreflec
end

""" 
    in_focus_acceptance(photon::Photon)::Bool

`in_focus_acceptance` checks whether the photon is within the focus acceptance region. 
    Acceptance is based on the slope of the photon trajectory. 
    Checks against minimum and maximum slope angles defined in FOCUS constants.

# Arguments
- `photon::Photon`: The photon instance.
- `focus::Focus`: The focus instance.

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

""" 
    _photon_direction(p::Particle, costhetac::Float64, sinthetac::Float64, phic::Float64)

`_photon_direction` computes the direction of the photon emitted by the particle in the lab frame. 
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
    #    uperp = sqrt(p.xDir^2 + p.yDir^2)
    #    if uperp == 0.0
    #        xdir = xprime
    #        ydir = yprime
    #        zdir = zprime
    #    else
    #        xdir = (p.xDir * p.zDir / uperp) * xprime -
    #               (p.yDir / uperp) * yprime + p.xDir * zprime
    #        ydir = (p.yDir * p.zDir / uperp) * xprime +
    #               (p.xDir / uperp) * yprime + p.yDir * zprime
    #        zdir = -uperp * xprime + p.zDir * zprime
    #    end

    xdir, ydir, zdir = rotate(p, xprime, yprime, zprime)
    slope = abs(zdir / ydir)

    return xdir, ydir, zdir, slope
end

""" 
    _photon_emission(p::Particle, zemission::Float64, context::PhotonContext)::Tuple{Float64,Float64,Float64}
    _photon_emission(p::Particle, context::PhotonContext)::Tuple{Float64,Float64,Float64}

`_photon_emission` calculates the photon emission coordinates at a specified z depth in the radiator.
    The emission point is determined by the particle's entry point and direction.
    For x and y coordinates, the values are adjusted based on the particle's direction.
    A single-argument method uses half the radiator depth as the emission point.

# Arguments
- `p::Particle`: The particle instance.
- `zemission::Float64`: Optional. Emission depth along the z-axis (mm). 
   If not provided, defaults to half the radiator depth.
- `context::PhotonContext`: The context containing radiator properties.

# Returns
- `Tuple{Float64,Float64,Float64}`: A tuple containing (xpos, ypos, zpos) of the emission point.

"""
function _photon_emission(
    p::Particle,
    zemission::Float64,
    context::PhotonContext,
)::Tuple{Float64,Float64,Float64}
    xpos = p.xCoord + zemission * (p.xDir / p.zDir)
    ypos = p.yCoord + zemission * (p.yDir / p.zDir)
    ypos += 0.5 * context.radiator.height
    return xpos, ypos, zemission
end

function _photon_emission(
    p::Particle,
    context::PhotonContext,
)::Tuple{Float64,Float64,Float64}
    # emission point in z
    zemission = context.radiator.depth * 0.5
    return _photon_emission(p, zemission, context)
end
