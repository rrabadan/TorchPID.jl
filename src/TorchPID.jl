module TorchPID

using Random
using YAML
using UnROOT
using Distributions

export RADIATOR, WEDGE, FOCUS, DETECTOR, MASK, SIGNAL
export HitCoordinate, PixelHit
export DetectorHitTester, photon_efficiency, test_photon
export PhotonSpectrum, PhotonSpectrumDistribution
export spectrum_nphase, spectrum_ngroup
export spectrum_yield, spectrum_random_energy, spectrum_probability
export Particle, particle_beta, particle_gamma, initial_rotation!, rotate, get_particle_mass
export PARTICLE_PROPERTIES
export Photon, test_z_surface_roughness, in_focus_acceptance
export PhotonMapper, trace_photon
export project_pattern
export ChargeDepositTester, charge_over_threshold, get_charge, get_smeared_time, smear_time
export FrontEnd, create_mcp_images, add_photon!, reset!, get_hits
export TestBeamSimulator, TestBeamData, generate_particle, generate_photons
# export EventReader, get_particle, photon_columns, track_columns
# export PIDAlgorithm, findHitCoordinates, runAlgorithm

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
#include("EventReader.jl")
#include("PIDAlgorithm.jl")

end # module TorchPID
