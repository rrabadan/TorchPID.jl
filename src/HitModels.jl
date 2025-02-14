struct HitCoordinate
    x::Float64
    y::Float64
    t::Float64
    track::Int
end

function HitCoordinate(x::Float64, y::Float64, t::Float64; track::Int=0)
    HitCoordinate(x, y, t, track)
end

struct PixelHit
    x::Int
    y::Int
    t::Int
    ch::Int
end

function PixelHit(x::Int, y::Int, t::Int; ch::Int=-1)
    PixelHit(x, y, t, ch)
end

struct DetectorHitTester
    scaleFactor::Float64
    emin::Float64
    emax::Float64
end

function DetectorHitTester(; scaleFactor::Float64=0.2, emin::Float64=1.75, emax::Float64=7.00)
    DetectorHitTester(scaleFactor, emin, emax)
end

function efficiency(dht::DetectorHitTester, energy::Float64; glueLayers::Int=1)::Float64
    if energy < dht.emin || energy > dht.emax
        return 0.0
    end
    return 1.0
end

function testHit(dht::DetectorHitTester, energy::Float64; glueLayers::Int=1)::Bool
    return efficiency(dht, energy, glueLayers=glueLayers) > rand()
end