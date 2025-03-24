struct HitCoordinate
    x::Float64
    y::Float64
    t::Float64
    track::Int
end

function HitCoordinate(x::Float64, y::Float64, t::Float64; track::Int = 0)
    HitCoordinate(x, y, t, track)
end

struct PixelHit
    x::Int
    y::Int
    t::Int
    ch::Int
end

function PixelHit(x::Int, y::Int, t::Int; ch::Int = -1)
    PixelHit(x, y, t, ch)
end
