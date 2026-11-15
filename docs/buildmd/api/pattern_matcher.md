
<a id='Pattern-Matcher-API'></a>

<a id='Pattern-Matcher-API-1'></a>

# Pattern Matcher API




Pattern matching algorithms and models used in TorchPID for particle identification.


<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.project_pattern!-Tuple{Vector{HitCoordinate}, Particle, Float64, PhotonMapper, PhotonSpectrum, PhotonSpectrumDistribution, PhotonContext}' href='#TorchPID.project_pattern!-Tuple{Vector{HitCoordinate}, Particle, Float64, PhotonMapper, PhotonSpectrum, PhotonSpectrumDistribution, PhotonContext}'>#</a>
**`TorchPID.project_pattern!`** &mdash; *Method*.



```julia
project_pattern(particle, beta, mapper, spectrum, distribution, context)
project_pattern(particle, beta, mapper, spectrum, context)
```

`project_pattern` calculates the hit coordinate pattern for a particle defined by its velocity factor (`beta`) and trajectory.      It generates photons with energies sampled from the `PhotonSpectrum` and `PhotonSpectrumDistribution`,      then traces each photon to the detector plane using the specified `PhotonMapper`.      A convenience method is available to construct the photon distribution internally,      delegating the computation to the primary `project_pattern` function.

**Arguments**

  * `particle::Particle`: Particle struct representing the particle's trajectory and timing.
  * `beta::Float64`: The particle's velocity factor relative to the speed of light.
  * `mapper::PhotonMapper`: Used to map photon properties to detector hit coordinates.
  * `spectrum::PhotonSpectrum`: Energy spectrum of photons.
  * `distribution::PhotonSpectrumDistribution`: Details of the photon distribution.
  * `context::PhotonContext`: Context containing the radiator and constants used for calculations.

**Returns**

  * `Vector{HitCoordinate}`: A vector of hit coordinates compatible with the provided spectrum and particle parameters.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/PatternMatcher.jl#L1-L21' class='documenter-source'>source</a><br>


<a id='Implementation-Details'></a>

<a id='Implementation-Details-1'></a>

### Implementation Details

<a id='TorchPID._get_time_offset-Tuple{Particle, Float64, Float64, Float64, PhotonContext}' href='#TorchPID._get_time_offset-Tuple{Particle, Float64, Float64, Float64, PhotonContext}'>#</a>
**`TorchPID._get_time_offset`** &mdash; *Method*.



```julia
_get_time_offset(particle, beta, depth, time, context)
```

`_get_time_offset` computes the time offset for a particle based on its velocity and traveled distance.      The calculation incorporates the particle's path length and its depth within the detector.      A simplified method is also available, which uses the particle's predefined `t0` value as the initial time.

**Arguments**

  * `particle::Particle`: Particle struct containing path length, initial time (t0) and z directional component.
  * `beta::Float64`: The particle's velocity factor relative to the speed of light.
  * `depth::Float64`: The measured depth in the detector in millimeters (mm).
  * `time::Float64`: The initial time offset for the particle.
  * `context::PhotonContext`: Context containing the radiator and constants used for calculations.

**Returns**

  * `Float64`: The computed time offset adjusted by the particle's path length and depth.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/PatternMatcher.jl#L87-L103' class='documenter-source'>source</a><br>

