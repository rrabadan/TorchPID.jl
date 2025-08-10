
<a id='Hit-Models-API'></a>

<a id='Hit-Models-API-1'></a>

# Hit Models API


Types and methods for modeling hit detection.




<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.HitCoordinate' href='#TorchPID.HitCoordinate'>#</a>
**`TorchPID.HitCoordinate`** &mdash; *Type*.



```julia
HitCoordinate(x, y, t; track=0)
```

Type representing a hit coordinate.

**Fields**

  * `x::Float64`: The x position of the hit.
  * `y::Float64`: The y position of the hit.
  * `t::Float64`: The time of the hit.
  * `track::Int=0`: An optional track identifier for the hit (default value is 0).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/HitModels.jl#L1-L11' class='documenter-source'>source</a><br>

<a id='TorchPID.PixelHit' href='#TorchPID.PixelHit'>#</a>
**`TorchPID.PixelHit`** &mdash; *Type*.



```julia
PixelHit(x, y, t; ch=-1)
```

Type representing a pixel hit.

**Fields**

  * `x::Int`: The x position of the pixel hit.
  * `y::Int`: The y position of the pixel hit.
  * `t::Int`: The time of the pixel hit.
  * `ch::Int=-1`: (Optional) The charge associated with the pixel hit (default value is -1).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/HitModels.jl#L23-L33' class='documenter-source'>source</a><br>

