module TorchPID

using Random
using YAML
using UnROOT
using Distributions
using Base.Threads
using StaticArrays
import Base: fill!

export CONSTANTS
export Focus, Detector, Mask, Signal
export GEOMETRY, RADIATOR, WEDGE, FOCUS, DETECTOR, MASK, SIGNAL
export HitCoordinate, PixelHit
export DetectorHitTester, photon_efficiency, test_photon
export PhotonSpectrum, PhotonSpectrumDistribution
export spectrum_nphase, spectrum_ngroup
export spectrum_yield,
    spectrum_random_energy,
    spectrum_probability,
    spectrum_probability_nocheck,
    propagation_properties,
    propagation_properties_nocheck,
    spectrum_random_sampling
export Particle, particle_beta, particle_gamma, initial_rotation!, rotate, get_particle_mass
export PARTICLE_PROPERTIES
export Photon,
    PhotonContext,
    create_approximate_photon,
    create_random_photon,
    create_explicit_photon,
    test_z_surface_roughness,
    in_focus_acceptance
export PhotonMapper, trace_photon
export ChargeDepositTester,
    charge_over_threshold, get_charge, get_smeared_time, smear_time, cached_time_weights
export TorchImage, TorchImageAccumulator, fill!, reset!
export project_pattern, make_pattern
export get_xpixel_position,
    get_ypixel_position, get_tpixel_position, get_hit_coordinate, is_in_pixel
export FrontEnd, create_mcp_images, add_photon!, reset!, get_hits
export TestBeamSimulator,
    TestBeamParticle, TestBeamPhotons, generate_particles, generate_photons_batch
export AbstractBackgroundModel, FlatBackground, get_bkg_probability, fill_bkg_probabilities!
export ReverseMappingConfig,
    MappedPhoton,
    focus_slope,
    focus_ypath,
    photon_jacobian,
    hit_signal_probability,
    accumulate_hit_signal_probs!,
    map_hit
export EventProbability,
    log_sum,
    pattern_efficiency,
    compute_event_probability!,
    solve_unbinned_likelihoods!,
    run_unbinned_pid!,
    run_unbinned_pid,
    run_unbinned_local_pid!
export run_local_pid
export PhotonHit,
    update_photon_hit!,
    photon_from_radiator,
    photon_from_focus,
    photon_on_detector,
    is_photon_masked
export create_event_reader,
    get_event_data, event_iterator, get_particles_in_event, get_hits_in_event
export PIDAlg, run_algorithm

include("PhysicsConstants.jl")
include("GeometryConstants.jl")
include("TorchFunctions.jl")
include("HitModels.jl")
include("DetectorHitTester.jl")
include("PhotonSpectrum.jl")
include("Particle.jl")
include("Photon.jl")
include("PhotonMapper.jl")
include("ChargeDepositTester.jl")
include("TorchImage.jl")
include("PatternMatcher.jl")
include("PixelMapper.jl")
include("FrontEnd.jl")
include("TestBeam.jl")
include("Background.jl")
include("ReverseMapping.jl")
include("SignalPDF.jl")
include("EventProbability.jl")
include("UnbinnedPID.jl")
include("LocalPID.jl")
include("PhotonHit.jl")
include("EventReader.jl")
include("PIDAlgorithm.jl")

end # module TorchPID
