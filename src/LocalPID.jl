function run_local_pid(
    particles::Vector{Particle},
    hits::Vector{PixelHit},
    hypotheses::Vector{Int},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    geometry::Geometry = GEOMETRY,
    constants::Constants = CONSTANTS,
)
    # Check if particles and hits are non-empty
    if isempty(particles) || isempty(hits)
        error("Particles and hits must not be empty")
    end

    image = TorchImageAccumulator(geometry.detector)
    charge_tester = ChargeDepositTester()
    photon_context = PhotonContext(geometry.radiator, constants)

    for particle in particles
        for hypothesis in hypotheses
            mass = get_particle_mass(hypothesis)
            beta = particle_beta(particle, mass)
            spectrum_distribution = PhotonSpectrumDistribution(spectrum, particle, beta)
            reset!(image)
            expected_yield = make_pattern(
                particle,
                beta,
                mapper,
                spectrum,
                spectrum_distribution,
                image,
                charge_tester,
                photon_context,
            )
        end
    end

end
