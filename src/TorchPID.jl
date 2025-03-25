module TorchPID

using Random
using YAML
using UnROOT
using Distributions

export RADIATOR, WEDGE, FOCUS, DETECTOR
export HitCoordinate, PixelHit
export DetectorHitTester, efficiency, testHit
export PhotonSpectrum, PhotonSpectrumDistribution
export spectrum_nphase, spectrum_ngroup
export spectrum_yield, spectrum_random_energy, spectrum_probability
export PARTICLE_PROPERTIES, Particle, beta, gamma, initRotation, rotate, get_particle_mass
export Photon, photon_test_z_surface_roughness
export TestBeamSimulator, TestBeamData, generate_particle
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
#include("PatternMatcher.jl")
#include("TestBeam.jl")
#include("EventReader.jl")
#include("PIDAlgorithm.jl")

end # module TorchPID
