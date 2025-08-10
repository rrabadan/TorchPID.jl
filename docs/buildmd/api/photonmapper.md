
<a id='PhotonMapper-API'></a>

<a id='PhotonMapper-API-1'></a>

# PhotonMapper API




Tools for simulating the photon propagation through TORCH's optical system.


<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.PhotonMapper' href='#TorchPID.PhotonMapper'>#</a>
**`TorchPID.PhotonMapper`** &mdash; *Type*.



```julia
PhotonMapper(blackened_sides blackened_bottom,
             blackened, surface_roughness,
             max_x_reflections)
```

Type representing the photon mapping configuration for photon propagation.

**Fields**

  * `blackened_sides::Bool`: If true, photons reflecting off the radiator's sides are absorbed.
  * `blackened_bottom::Bool`: If true, photons reflecting off the radiator's bottom are ignored.
  * `blackened_focus::Bool`: If true, reflections in the focus region are handled differently.
  * `surface_roughness::Bool`: If true, surface imperfections of the radiator are considered.
  * `max_x_reflections::INt`: Specifies the maximum allowable x-direction reflections before the photon is discarded.

**Constructors**

```
PhotonMapper(; 
    blackened_sides::Bool=false,
    blackened_bottom::Bool=false,
    blackened_focus::Bool=false,
    surface_roughness::Bool=false,
    max_x_reflections::Int=10
)
```

**Keywords**

  * `blackened_sides::Bool=false`: Enable absorption of side reflections.
  * `blackened_bottom::Bool=false`: Ignore bottom reflections if true.
  * `blackened_focus::Bool=false`: Handle focus region reflections differently.
  * `surface_roughness::Bool=false`: Account for surface imperfections.
  * `max_x_reflections::Int=10`: Maximum number of x-direction reflections allowed.

**Returns**

  * `PhotonMapper`: A new PhotonMapper instance with the specified configuration.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L1-L34' class='documenter-source'>source</a><br>

<a id='TorchPID.trace_photon-Tuple{PhotonMapper, Photon, Float64}' href='#TorchPID.trace_photon-Tuple{PhotonMapper, Photon, Float64}'>#</a>
**`TorchPID.trace_photon`** &mdash; *Method*.



```julia
trace_photon(mapper::PhotonMapper, photon::Photon, t0::Float64)
```

`trace_photon` simulates the photon's journey through the optical system, from emission to detection.  It propagates the photon to the top of the plate, accounting for y-direction and z-direction reflections, surface roughness, and pathlength accumulation. It then traces the photon to the focus block, simulating its interaction with the focusing optics, and finally to the detector plane, where it checks for acceptance criteria. The function returns the detected hit coordinate if the photon is successfully traced, or `nothing` if any tracing step fails.

**Arguments**

  * `mapper::PhotonMapper`: The photon mapping configuration.
  * `photon::Photon`: The photon to trace.
  * `t0::Float64`: The production time of the photon.

**Returns**

  * `Union{HitCoordinate,Nothing}`: The detected hit coordinate if the photon is successfully traced, or `nothing` if any tracing step fails.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L59-L75' class='documenter-source'>source</a><br>


<a id='Implementation-Details'></a>

<a id='Implementation-Details-1'></a>

### Implementation Details


The `trace_photon` function internally uses helper functions to propagate the cherenkov emitted photon.


These internal functions are not part of the public API and may change without notice.

<a id='TorchPID._project_x_detector_position-Tuple{PhotonMapper, Photon, Float64}' href='#TorchPID._project_x_detector_position-Tuple{PhotonMapper, Photon, Float64}'>#</a>
**`TorchPID._project_x_detector_position`** &mdash; *Method*.



```julia
_project_x_detector_position(mapper::PhotonMapper, photon::Photon, pathlength::Float64)
```

`_project_x_detector_position` calculates the photon's x-coordinate on the detector plane by projecting its trajectory. It considers x-direction reflections based on the mapper's configuration. Returns `nothing` if side blackening is enabled and reflection is required, the photon cannot reflect in x, or the maximum allowed reflections is exceeded.

**Arguments**

  * `mapper::PhotonMapper`: The mapping configuration.
  * `photon::Photon`: The photon instance.
  * `pathlength::Float64`: Total travelled pathlength.

**Returns**

  * `Union{Float64,Nothing}`: The projected x position on the detector, or `nothing` if the photon does not reach the detector due to reflection conditions.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L329-L344' class='documenter-source'>source</a><br>

<a id='TorchPID._trace_to_detector-Tuple{Float64, Float64, Float64}' href='#TorchPID._trace_to_detector-Tuple{Float64, Float64, Float64}'>#</a>
**`TorchPID._trace_to_detector`** &mdash; *Method*.



```julia
_trace_to_detector(yf::Float64, tyf::Float64, tzf::Float64)
```

`_trace_to_detector` simulates the photon's propagation from the entrance of the focusing block to the detector plane, accounting for mirror reflections and acceptance criteria.

**Arguments**

  * `yf::Float64`: Position at entrance of the focus block
  * `tyf::Float64`: Y-component of photon momentum at entrance of the focus block
  * `tzf::Float64`: Z-component of photon momentum at entrance of the focus block

**Returns**

  * `Union{Tuple{Float64,Float64},Nothing}`: A tuple containing:

      * `ydetected::Float64`: Position of photon on detector plane
      * `pathlength::Float64`: Accumulated pathlength

    Returns `nothing` if the photon's path is terminated.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L240-L257' class='documenter-source'>source</a><br>

<a id='TorchPID._trace_to_focus-Tuple{Float64, Float64, Float64}' href='#TorchPID._trace_to_focus-Tuple{Float64, Float64, Float64}'>#</a>
**`TorchPID._trace_to_focus`** &mdash; *Method*.



```julia
_trace_to_focus(zp::Float64, typ::Float64, tzp::Float64)
```

`_trace_to_focus` propagates the photon from the top of the flat section to the entrance of the focus block.

**Arguments**

  * `zp::Float64`: Position at top of the plate
  * `typ::Float64`: Y-component of photon momentum at top of the plate
  * `tzp::Float64`: Z-component of photon momentum at top of the plate

**Returns**

  * `Union{Tuple{Float64,Float64,Float64,Float64},Nothing}`: A tuple containing:

      * `yf::Float64`: Position in focus block coordinates
      * `tyf::Float64`: Y-component of photon momentum in focus block coordinates
      * `tzf::Float64`: Z-component of photon momentum in focus block coordinates
      * `pathlength::Float64`: Accumulated pathlength

    Returns `nothing` if the photon's path is terminated.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L195-L213' class='documenter-source'>source</a><br>

<a id='TorchPID._trace_to_top_of_plate-Tuple{PhotonMapper, Photon}' href='#TorchPID._trace_to_top_of_plate-Tuple{PhotonMapper, Photon}'>#</a>
**`TorchPID._trace_to_top_of_plate`** &mdash; *Method*.



```julia
_trace_to_top_of_plate(mapper::PhotonMapper, photon::Photon)
```

`_trace_to_top_of_plate` propagates the photon to the top of the plate, accounting for y-direction and z-direction reflections, surface roughness, and pathlength accumulation.

**Arguments**

  * `mapper::PhotonMapper`: The mapping configuration.
  * `photon::Photon`: The photon to be traced.

**Returns**

  * `Union{Tuple{Float64,Float64,Float64,Float64},Nothing}`: A tuple containing:

      * `zp::Float64`: Position in z
      * `typ::Float64`: Y-component of photon momentum
      * `tzp::Float64`: Z-component of photon momentum
      * `pathlength::Float64`: Accumulated pathlength

    Returns `nothing` if the photon's path is terminated.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonMapper.jl#L133-L150' class='documenter-source'>source</a><br>

