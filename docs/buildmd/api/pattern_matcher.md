
<a id='Pattern-Matcher-API'></a>

<a id='Pattern-Matcher-API-1'></a>

# Pattern Matcher API




Pattern matching algorithms and models used in TorchPID for particle identification.


<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.project_pattern-Tuple{Particle, Float64, PhotonMapper, PhotonSpectrum, PhotonSpectrumDistribution, PhotonFactory}' href='#TorchPID.project_pattern-Tuple{Particle, Float64, PhotonMapper, PhotonSpectrum, PhotonSpectrumDistribution, PhotonFactory}'>#</a>
**`TorchPID.project_pattern`** &mdash; *Method*.



```julia
project_pattern(particle, beta, mapper, spectrum, distribution)
project_pattern(particle, beta, mapper, spectrum)
```

`project_pattern` calculates the hit coordinate pattern for a particle defined by its velocity factor (`beta`) and trajectory.      It generates photons with energies sampled from the `PhotonSpectrum` and `PhotonSpectrumDistribution`,      then traces each photon to the detector plane using the specified `PhotonMapper`.      A convenience method is available to construct the photon distribution internally,      delegating the computation to the primary `project_pattern` function.

**Arguments**

  * `particle::Particle`: Particle struct representing the particle's trajectory and timing.
  * `beta::Float64`: The particle's velocity factor relative to the speed of light.
  * `mapper::PhotonMapper`: Used to map photon properties to detector hit coordinates.
  * `spectrum::PhotonSpectrum`: Energy spectrum of photons.
  * `distribution::PhotonSpectrumDistribution`: Details of the photon distribution.

**Returns**

  * `Vector{HitCoordinate}`: A vector of hit coordinates compatible with the provided spectrum and particle parameters.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PatternMatcher.jl#L1-L20' class='documenter-source'>source</a><br>


<a id='Implementation-Details'></a>

<a id='Implementation-Details-1'></a>

### Implementation Details

<a id='TorchPID._get_time_offset-Tuple{Particle, Float64, Float64, Float64}' href='#TorchPID._get_time_offset-Tuple{Particle, Float64, Float64, Float64}'>#</a>
**`TorchPID._get_time_offset`** &mdash; *Method*.



```julia
_get_time_offset(particle, beta, depth, time)
```

`_get_time_offset` computes the time offset for a particle based on its velocity and traveled distance.      The calculation incorporates the particle's path length and its depth within the detector.      A simplified method is also available, which uses the particle's predefined `t0` value as the initial time.

**Arguments**

  * `particle::Particle`: Particle struct containing path length, initial time (t0) and z directional component.
  * `beta::Float64`: The particle's velocity factor relative to the speed of light.
  * `depth::Float64`: The measured depth in the detector in millimeters (mm).
  * `time::Float64`: The initial time offset for the particle.

**Returns**

  * `Float64`: The computed time offset adjusted by the particle's path length and depth.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PatternMatcher.jl#L73-L88' class='documenter-source'>source</a><br>

