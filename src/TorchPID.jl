module TorchPID

using Random
using YAML
using UnROOT
using Distributions
using Base.Threads
using StaticArrays

export CONSTANTS
export GEOMETRY, RADIATOR, WEDGE, FOCUS, DETECTOR, MASK, SIGNAL
export HitCoordinate, PixelHit
export DetectorHitTester, photon_efficiency, test_photon
export PhotonSpectrum, PhotonSpectrumDistribution
export spectrum_nphase, spectrum_ngroup
export spectrum_yield, spectrum_random_energy, spectrum_probability
export Particle, particle_beta, particle_gamma, initial_rotation!, rotate, get_particle_mass
export PARTICLE_PROPERTIES
export Photon,
    PhotonFactory,
    create_random_photon,
    create_explicit_photon,
    test_z_surface_roughness,
    in_focus_acceptance
export PhotonMapper, trace_photon
export project_pattern
export ChargeDepositTester, charge_over_threshold, get_charge, get_smeared_time, smear_time
export FrontEnd, create_mcp_images, add_photon!, reset!, get_hits
export TestBeamSimulator,
    TestBeamParticle, TestBeamPhotons, generate_particle, generate_photons
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
include("PatternMatcher.jl")
include("PixelMapper.jl")
include("ChargeDepositTester.jl")
include("FrontEnd.jl")
include("TestBeam.jl")
include("PhotonHit.jl")
include("EventReader.jl")
include("PIDAlgorithm.jl")

end # module TorchPID
