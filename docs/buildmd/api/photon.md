
<a id='Photon-API'></a>

<a id='Photon-API-1'></a>

# Photon API




Core definitions and utilities for photon properties and behavior.


<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.Photon' href='#TorchPID.Photon'>#</a>
**`TorchPID.Photon`** &mdash; *Type*.



```julia
Photon(energy, xdir, ydir, zdir, xpos, ypos, zpos, slope,
       np, ng, lambda, emissionTime, is_x_reflected,
       is_y_reflected, is_z_reflected)
```

Type representing a photon emitted by TORCH's radiator and propagated to the detector. Position coordinates are in millimeters (mm).

**Fields**

  * `energy::Float64`: Photon energy.
  * `xdir::Float64`: X-component of the photon direction in the lab frame.
  * `ydir::Float64`: Y-component of the photon direction in the lab frame.
  * `zdir::Float64`: Z-component of the photon direction in the lab frame.
  * `xpos::Float64`: X-coordinate of the emission point.
  * `ypos::Float64`: Y-coordinate of the emission point.
  * `zpos::Float64`: Z-coordinate of the emission point.
  * `slope::Float64`: Slope of the trajectory.
  * `np::Float64`: Phase refractive index.
  * `ng::Float64`: Group refractive index.
  * `lambda::Float64`: Wavelength computed as LAMBDA/energy.
  * `emissionTime::Float64`: Time at emission.
  * `is_x_reflected::Bool`: Flag indicating reflection status in x-direction.
  * `is_y_reflected::Bool`: Flag indicating reflection status in y-direction.
  * `is_z_reflected::Bool`: Flag indicating reflection status in z-direction.

**Constructors**

```
Photon(p, beta, nphase, ngroup, energy)
```

Creates a Photon instance by determining its direction, emission point, and reflection flags based on the particle's properties, photon energy, and refractive indices.  The Cherenkov angle is calculated using the particle's velocity and the phase refractive index.  A random azimuthal angle (phi) is uniformly sampled between 0 and 2π to simulate isotropic emission.  The emission point is randomly chosen along the radiator depth to reflect realistic photon production.  Reflection flags are set by comparing the photon's direction with the critical angle for total internal reflection.  Surface roughness effects are excluded from this constructor and are handled separately.

**Arguments**

  * `p::Particle`: The particle instance emitting the photon.
  * `beta::Float64`: The velocity factor of the particle (v/c).
  * `nphase::Float64`: The phase refractive index for the specific photon energy.
  * `ngroup::Float64`: The group refractive index for the specific photon energy.
  * `energy::Float64`: The photon energy in appropriate units.

**Returns**

  * `Photon`: A new Photon object populated with calculated properties.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/8d8443ddd5309103f4b5f28bca51018042c83983/src/Photon.jl#L1-L45' class='documenter-source'>source</a><br>

<a id='TorchPID.in_focus_acceptance-Tuple{Photon}' href='#TorchPID.in_focus_acceptance-Tuple{Photon}'>#</a>
**`TorchPID.in_focus_acceptance`** &mdash; *Method*.



```julia
in_focus_acceptance(photon::Photon)::Bool
```

`in_focus_acceptance` checks whether the photon is within the focus acceptance region.      Acceptance is based on the slope of the photon trajectory.      Checks against minimum and maximum slope angles defined in FOCUS constants.

**Arguments**

  * `photon::Photon`: The photon instance.

**Returns**

  * `Bool`: Whether the photon is within the acceptance region (true) or not (false).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/8d8443ddd5309103f4b5f28bca51018042c83983/src/Photon.jl#L131-L145' class='documenter-source'>source</a><br>

<a id='TorchPID.test_z_surface_roughness-Tuple{Photon, Int64}' href='#TorchPID.test_z_surface_roughness-Tuple{Photon, Int64}'>#</a>
**`TorchPID.test_z_surface_roughness`** &mdash; *Method*.



```julia
test_z_surface_roughness(photon::Photon, nreflec::Int)::Bool
```

`test_z_surface_roughness` tests the photon surface roughness in the z-direction.      Surface roughness is modeled as a Gaussian effect dependent on wavelength.      The probability of passing decreases with number of reflections.

**Arguments**

  * `photon::Photon`: The photon instance.
  * `nreflec::Int`: The number of reflections.

**Returns**

  * `Bool`: Whether the photon passes the surface roughness condition (true) or not (false).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/8d8443ddd5309103f4b5f28bca51018042c83983/src/Photon.jl#L111-L124' class='documenter-source'>source</a><br>


<a id='Implementation-Details'></a>

<a id='Implementation-Details-1'></a>

### Implementation Details


The Photon constructor internally uses helper functions for:


  * Calculating photon directions (`_photon_direction`)
  * Determining emission positions (`_photon_emission`)


These internal functions are not part of the public API and may change without notice.

<a id='TorchPID._photon_direction-Tuple{Particle, Float64, Float64, Float64}' href='#TorchPID._photon_direction-Tuple{Particle, Float64, Float64, Float64}'>#</a>
**`TorchPID._photon_direction`** &mdash; *Method*.



```julia
_photon_direction(p::Particle, costhetac::Float64, sinthetac::Float64, phic::Float64)
```

`_photon_direction` computes the direction of the photon emitted by the particle in the lab frame.      The photon direction is determined from the Cherenkov angle and azimuthal angle and      rotated to the lab frame using the particle's direction.

**Arguments**

  * `p::Particle`: The particle instance.
  * `costhetac::Float64`: Cosine of the Cherenkov angle.
  * `sinthetac::Float64`: Sine of the Cherenkov angle.
  * `phic::Float64`: Azimuthal angle around the particle direction.

**Returns**

  * `Tuple{Float64,Float64,Float64,Float64}`: A tuple containing (xdir, ydir, zdir, slope) of the photon direction.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/8d8443ddd5309103f4b5f28bca51018042c83983/src/Photon.jl#L153-L169' class='documenter-source'>source</a><br>

<a id='TorchPID._photon_emission-Tuple{Particle, Float64}' href='#TorchPID._photon_emission-Tuple{Particle, Float64}'>#</a>
**`TorchPID._photon_emission`** &mdash; *Method*.



```julia
_photon_emission(p::Particle, zemission::Float64)::Tuple{Float64,Float64,Float64}
_photon_emission(p::Particle)::Tuple{Float64,Float64,Float64}
```

`_photon_emission` calculates the photon emission coordinates at a specified z depth in the radiator.     The emission point is determined by the particle's entry point and direction.     For x and y coordinates, the values are adjusted based on the particle's direction.     A single-argument method uses half the radiator depth as the emission point.

**Arguments**

  * `p::Particle`: The particle instance.
  * `zemission::Float64`: Optional. Emission depth along the z-axis (mm).   If not provided, defaults to half the radiator depth.

**Returns**

  * `Tuple{Float64,Float64,Float64}`: A tuple containing (xpos, ypos, zpos) of the emission point.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/8d8443ddd5309103f4b5f28bca51018042c83983/src/Photon.jl#L200-L218' class='documenter-source'>source</a><br>

