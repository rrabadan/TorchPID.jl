"""
    HitCoordinate(x, y, t; track=0)

Type representing a hit coordinate.

# Fields
- `x::Float64`: The x position of the hit.
- `y::Float64`: The y position of the hit.
- `t::Float64`: The time of the hit.
- `track::Int=0`: An optional track identifier for the hit (default value is 0).
"""
struct HitCoordinate
    x::Float64
    y::Float64
    t::Float64
    track::Int
end

function HitCoordinate(x::Float64, y::Float64, t::Float64; track::Int = 0)
    HitCoordinate(x, y, t, track)
end

"""
    PixelHit(x, y, t; ch=-1)

Type representing a pixel hit.

# Fields
- `x::Int`: The x position of the pixel hit.
- `y::Int`: The y position of the pixel hit.
- `t::Int`: The time of the pixel hit.
- `ch::Int=-1`: (Optional) The charge associated with the pixel hit (default value is -1).
"""
struct PixelHit
    x::Int
    y::Int
    t::Int
    ch::Float64
end

function PixelHit(x::Int, y::Int, t::Int; ch::Int = -1)
    PixelHit(x, y, t, ch)
end
