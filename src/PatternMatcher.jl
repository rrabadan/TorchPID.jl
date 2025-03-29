""" 
Compute the hit coordinates pattern with the given spectrum, distribution, mapper, particle and beta.

Parameters:
  - s: PhotonSpectrum representing the energy spectrum of photons.
  - d: PhotonSpectrumDistribution containing details of the photon distribution.
  - mapper: PhotonMapper used to map photon properties to detector hit coordinates.
  - particle: Particle struct representing the particle's trajectory and timing.
  - beta: Float64 representing the particle's velocity factor.

Returns:
  - Vector{HitCoordinate}: A vector of hit coordinates compatible with the provided spectrum and particle parameters.
"""
function project_pattern(
    particle::Particle,
    beta::Float64,
    spectrum::PhotonSpectrum,
    distribution::PhotonSpectrumDistribution,
    mapper::PhotonMapper,
)::Vector{HitCoordinate}
    hits = HitCoordinate[]

    if !spectrum_above_threshold(distribution)
        return hits
    end

    pathlength = RADIATOR.depth / particle.zDir
    photonYield = spectrum_yield(distribution, pathlength)

    # Poisson distribution
    nphotons = rand(Poisson(photonYield))

    for _ = 1:nphotons
        energy = spectrum_random_energy(spectrum, distribution)
        n_phase = nphase_Corning(energy)
        n_group = ngroup_Corning(energy)

        photon = Photon(particle, beta, n_phase, n_group, energy)

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

""" 
Constructs the photon spectrum distribution based on the provided spectrum and particle parameters,
and computes the corresponding hit coordinates compatible with the provided spectrum and particle parameters.

Parameters:
  - s: PhotonSpectrum representing the energy spectrum.
  - mapper: PhotonMapper used to map photon properties to detector hit coordinates.
  - particle: Particle struct containing the particle's trajectory and timing information.
  - beta: Float64 representing the particle's velocity factor.

Returns:
  - Vector{HitCoordinate}: A vector of hit coordinates computed by this process.
"""
function project_pattern(
    particle::Particle,
    beta::Float64,
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
)::Vector{HitCoordinate}
    photonDistribution = PhotonSpectrumDistribution(spectrum, particle, beta)
    return projec_pattern(particle, beta, spectrum, photonDistribution, mapper)
end

""" 
Calculate the time offset for a given particle.

Parameters:
  - particle: Particle struct containing path length, initial time (t0) and z directional component.
  - beta: Float64 representing the particle's velocity factor relative to the speed of light.
  - depth: Float64 representing the measured depth in the detector.
  - time: Float64 representing the initial time offset for the particle.

Returns:
  - Float64: The computed time offset adjusted by the particle's path length and depth.
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

""" 
Get the time offset using the particle's default t0.

Parameters:
  - particle: Particle struct with default t0.
  - beta: Float64 representing the particle's velocity factor.
  - depth: Float64 representing the measured depth in the detector.

Returns:
  - Float64: The computed time offset using particle.t0.
"""
function _get_time_offset(particle::Particle, beta::Float64, depth::Float64)::Float64
    return _get_time_offset(particle, beta, depth, particle.t0)
end
