"""
    project_pattern(particle, beta, mapper, spectrum, distribution)
    project_pattern(particle, beta, mapper, spectrum)

`project_pattern` calculates the hit coordinate pattern for a particle defined by its velocity factor (`beta`) and trajectory. 
    It generates photons with energies sampled from the `PhotonSpectrum` and `PhotonSpectrumDistribution`, 
    then traces each photon to the detector plane using the specified `PhotonMapper`. 
    A convenience method is available to construct the photon distribution internally, 
    delegating the computation to the primary `project_pattern` function.

# Arguments
- `particle::Particle`: Particle struct representing the particle's trajectory and timing.
- `beta::Float64`: The particle's velocity factor relative to the speed of light.
- `mapper::PhotonMapper`: Used to map photon properties to detector hit coordinates.
- `spectrum::PhotonSpectrum`: Energy spectrum of photons.
- `distribution::PhotonSpectrumDistribution`: Details of the photon distribution.

# Returns
- `Vector{HitCoordinate}`: A vector of hit coordinates compatible with the provided spectrum and particle parameters.
"""
function project_pattern(
    particle::Particle,
    beta::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    distribution::PhotonSpectrumDistribution,
    factory::PhotonFactory,
)::Vector{HitCoordinate}
    hits = HitCoordinate[]

    if !spectrum_above_threshold(distribution)
        return hits
    end

    pathlength = factory.radiator.depth / particle.zDir
    photonYield = spectrum_yield(distribution, pathlength)

    # Poisson distribution
    nphotons = rand(Poisson(photonYield))

    for _ = 1:nphotons
        energy = spectrum_random_energy(spectrum, distribution)
        n_phase = nphase_Corning(energy)
        n_group = ngroup_Corning(energy)

        photon = create_random_photon(factory, particle, beta, n_phase, n_group, energy)

        # Calculate the time offset
        time_offset = _get_time_offset(particle, beta, photon.zpos)

        hit = trace_photon(mapper, photon, time_offset)

        if isnothing(hit)
            continue
        end

        push!(hits, hit)
    end
    return hits
end

function project_pattern(
    particle::Particle,
    beta::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    factory::PhotonFactory,
)::Vector{HitCoordinate}
    photonDistribution = PhotonSpectrumDistribution(spectrum, particle, beta)
    return project_pattern(particle, beta, mapper, spectrum, photonDistribution, factory)
end

"""
    _get_time_offset(particle, beta, depth, time)

`_get_time_offset` computes the time offset for a particle based on its velocity and traveled distance. 
    The calculation incorporates the particle's path length and its depth within the detector. 
    A simplified method is also available, which uses the particle's predefined `t0` value as the initial time.

# Arguments
- `particle::Particle`: Particle struct containing path length, initial time (t0) and z directional component.
- `beta::Float64`: The particle's velocity factor relative to the speed of light.
- `depth::Float64`: The measured depth in the detector in millimeters (mm).
- `time::Float64`: The initial time offset for the particle.

# Returns
- `Float64`: The computed time offset adjusted by the particle's path length and depth.
"""
function _get_time_offset(
    particle::Particle,
    beta::Float64,
    depth::Float64,
    time::Float64,
)::Float64
    distance = particle.pathlength + (depth / particle.zDir)
    return time + distance / (beta * CLIGHT)
end

function _get_time_offset(particle::Particle, beta::Float64, depth::Float64)::Float64
    return _get_time_offset(particle, beta, depth, particle.t0)
end


struct RandomPoint
    e::Float64
    emissionZ::Float64
    phi::Float64
    cosPhi::Float64
    sinPhi::Float64
end

function random_point()
    point = (randn(), randn(), randn())
    e = point[3]
    emissionZ = point[2] * 10 # Convert to mm
    phi = point[1] * 2 * π
    cosPhi = cos(phi)
    sinPhi = sin(phi)
    return RandomPoint(e, emissionZ, phi, cosPhi, sinPhi)
end

function make_pattern(
    particle::Particle,
    beta::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    distribution::PhotonSpectrumDistribution,
    factory::PhotonFactory,
    n_generations::Int=1000000,
)::Float64
    for i = 1:n_generations
        point = random_point()
        spectra = spectrum_random_sampling(spectrum, distribution, point.e)
        if isnothing(spectra)
            continue
        end
        energy, nphase, ngroup = spectra
        photon = create_approximate_photon(
            factory,
            particle,
            beta,
            nphase,
            ngroup,
            energy,
            point.emissionZ,
            point.phi,
        )
        time_offset = _get_time_offset(particle, beta, photon.zpos, particle.t0)

        hit = trace_photon(mapper, photon, time_offset)
        if isnothing(hit)
            continue
        end
    end
end
