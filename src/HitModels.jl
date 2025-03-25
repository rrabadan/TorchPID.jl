"""
    HitCoordinate

Struct representing a hit coordinate with:
      x :: Float64 - x position
      y :: Float64 - y position
      t :: Float64 - time
      track :: Int - track identifier
"""
struct HitCoordinate
    x::Float64
    y::Float64
    t::Float64
    track::Int
end

"""
    HitCoordinate(x, y, t; track=0)

Construct a HitCoordinate given x, y, and t values with an optional track (default=0).
"""
function HitCoordinate(x::Float64, y::Float64, t::Float64; track::Int = 0)
    HitCoordinate(x, y, t, track)
end

""" 
    PixelHit

Struct representing a pixel hit with:
      x :: Int - x position
      y :: Int - y position
      t :: Int - time
      ch :: Int - channel identifier
"""
struct PixelHit
    x::Int
    y::Int
    t::Int
    ch::Int
end

""" 
    PixelHit(x, y, t; ch=-1)

Construct a PixelHit given x, y, and t values with an optional channel (default=-1).
"""
function PixelHit(x::Int, y::Int, t::Int; ch::Int = -1)
    PixelHit(x, y, t, ch)
end
