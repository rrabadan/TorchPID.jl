"""
Represents a photon with its energy, direction, emission point, and reflection flags.

Fields:
- energy: Photon energy.
- xdir, ydir, zdir: Direction components in the lab frame.
- xpos, ypos, zpos: Emission point coordinates.
- slope: Slope of the trajectory.
- np: Phase refractive index.
- ng: Group refractive index.
- lambda: Wavelength computed as LAMBDA/energy.
- emissionTime: Time at emission.
- is_x_reflected, is_y_reflected, is_z_reflected: Flags indicating reflection status.
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
Compute the photon direction in the lab frame.
# Arguments
- p::Particle: Input particle.
- costhetac::Float64: Cosine of the critical angle.
- sinthetac::Float64: Sine of the critical angle.
- phi::Float64: Azimuthal angle.
# Returns
Tuple{Float64,Float64,Float64,Float64} representing (xdir, ydir, zdir, slope).
"""
function photon_direction(
    p::Particle,
    costhetac::Float64,
    sinthetac::Float64,
    phi::Float64,
)::Tuple{Float64,Float64,Float64,Float64}
    # photon direction
    xprime = sinthetac * cos(phi)
    yprime = sinthetac * sin(phi)
    zprime = costhetac

    # rotate to the lab frame
    xdir, ydir, zdir = rotate(p, xprime, yprime, zprime)
    slope = abs(xdir / zdir)

    return xdir, ydir, zdir, slope
end

""" 
Calculate the photon emission coordinates at a specified z depth.
# Arguments
- p::Particle: Input particle.
- zemission::Float64: Emission depth along the z-axis.
# Returns
Tuple{Float64,Float64,Float64} representing (xpos, ypos, zpos).
"""
function photon_emission(p::Particle, zemission::Float64)::Tuple{Float64,Float64,Float64}
    xpos = p.xCoord + zemission * p.xDir / p.zDir
    ypos = p.yCoord + zemission * p.yDir / p.zDir
    return xpos, ypos, zemission
end

""" 
Calculate the photon emission coordinates using the default radiator depth.
# Arguments
- p::Particle: Input particle.
# Returns
Tuple{Float64,Float64,Float64} representing (xpos, ypos, zpos).
"""
function photon_emission(p::Particle)::Tuple{Float64,Float64,Float64}
    # emission point in z
    zemission = RADIATOR.depth * 0.5
    return photon_emission(p, zemission)
end

"""
Test the photon surface roughness in the z-direction.

# Arguments
- photon::Photon: The Photon instance.
- nreflec::Int: The number of reflections.

# Returns
Bool indicating whether the photon passes the roughness condition.
"""
function photon_test_z_surface_roughness(photon::Photon, nreflec::Int)::Bool
    c = π * photon.np * ROUGHNESS / photon.lambda
    r = exp(-c * c)
    return rand() < r^nreflec
end

""" 
Generate a Photon instance with computed direction, emission, and reflection flags.
# Arguments
- p::Particle: Input particle.
- beta::Float64: Velocity factor.
- nphase::Float64: Phase refractive index.
- ngroup::Float64: Group refractive index.
- energy::Float64: Photon energy.
# Returns
A Photon struct populated with calculated properties.
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
    xdir, ydir, zdir, slope = photon_direction(p, costhetac, sinthetac, phi)

    # emission point in z
    xpos, ypos, zpos = photon_emission(p, RADIATOR.depth * rand())

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
