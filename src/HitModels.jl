"""
Struct representing a hit coordinate with:
  - `x::Float64`: x position
  - `y::Float64`: y position
  - `t::Float64`: time
  - `track::Int`: track identifier
"""
struct HitCoordinate
    x::Float64
    y::Float64
    t::Float64
    track::Int
end

"""
HitCoordinate(x, y, t; track=0)

Construct a `HitCoordinate` object given x, y, and t values with an optional `track` identifier.

# Arguments
- `x::Float64`: The x position of the hit.
- `y::Float64`: The y position of the hit.
- `t::Float64`: The time of the hit.
- `track::Int=0`: An optional track identifier for the hit (default value is 0).

# Returns
- `HitCoordinate` struct.
"""
function HitCoordinate(x::Float64, y::Float64, t::Float64; track::Int = 0)
    HitCoordinate(x, y, t, track)
end

"""
Struct representing a pixel hit with:
  - `x::Int`: x position.
  - `y::Int`: y position.
  - `t::Int`: time.
  - `ch::Float64`: charge.
"""
struct PixelHit
    x::Int
    y::Int
    t::Int
    ch::Float64
end

"""
PixelHit(x, y, t; ch=-1)

Construct a `PixelHit` object given x, y, and t values and optional charge value.

# Arguments
- `x::Int`: The x position of the pixel hit.
- `y::Int`: The y position of the pixel hit.
- `t::Int`: The time of the pixel hit.
- `ch::Int=-1`: (Optional) The charge associated with the pixel hit (default value is -1).

# Returns
- `PixelHit` struct.
"""
function PixelHit(x::Int, y::Int, t::Int; ch::Int = -1)
    PixelHit(x, y, t, ch)
end
