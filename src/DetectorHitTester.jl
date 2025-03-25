"""
    DetectorHitTester

Structure to hold configuration for performing detector hit tests.

# Fields
- `scaleFactor::Float64`: Factor applied to random value.
- `emin::Float64`: Minimum energy threshold.
- `emax::Float64`: Maximum energy threshold.
- `implement_epotek_cutoff::Bool`: Enable Epotek cutoff.
- `implement_imperfect_mirror::Bool`: Enable imperfect mirror behavior.
- `implement_QE::Bool`: Enable quantum efficiency.
- `implement_CE::Bool`: Enable collection efficiency.
"""
struct DetectorHitTester
    scaleFactor::Float64
    emin::Float64
    emax::Float64
    implement_epotek_cutoff::Bool
    implement_imperfect_mirror::Bool
    implement_QE::Bool
    implement_CE::Bool
end

"""
    DetectorHitTester(; scaleFactor, emin, emax, implement_epotek_cutoff,
    implement_imperfect_mirror, implement_QE, implement_CE)

Constructor for DetectorHitTester with default values.

# Keyword Arguments
- `scaleFactor::Float64=1.0`
- `emin::Float64=1.75`
- `emax::Float64=7.00`
- `implement_epotek_cutoff::Bool=true`
- `implement_imperfect_mirror::Bool=true`
- `implement_QE::Bool=true`
- `implement_CE::Bool=true`
"""
function DetectorHitTester(;
    scaleFactor::Float64 = 1.0,
    emin::Float64 = 1.75,
    emax::Float64 = 7.00,
    implement_epotek_cutoff::Bool = true,
    implement_imperfect_mirror::Bool = true,
    implement_QE::Bool = true,
    implement_CE::Bool = true,
)
    DetectorHitTester(
        scaleFactor,
        emin,
        emax,
        implement_epotek_cutoff,
        implement_imperfect_mirror,
        implement_QE,
        implement_CE,
    )
end

"""
    efficiency(hitTester, energy; glueLayers=1) -> Float64

Calculates the efficiency for a given detector hit based on energy and configuration options.

# Arguments
- `hitTester::DetectorHitTester`: Instance holding configuration parameters.
- `energy::Float64`: The energy value to test.
- `glueLayers::Int=1`: Number of glue layers affecting the efficiency.

# Returns
- `Float64`: Calculated efficiency (0.0 when energy is out of threshold bounds).
"""
function efficiency(
    hitTester::DetectorHitTester,
    energy::Float64;
    glueLayers::Int = 1,
)::Float64
    if energy < dht.emin || energy > dht.emax
        return 0.0
    end
    result = 1.0
    if hitTester.implement_QE
        result *= QE_interpolator(energy)
    end
    if hitTester.implement_CE
        result *= CE()
    end
    if hitTester.implement_epotek_cutoff
        result *= epotek_305_interpolator(energy)^glueLayers
    end
    if hitTester.implement_imperfect_mirror
        result *= mirror_reflect(energy)
    end
    return result
end

"""
    testHit(hitTester, energy; glueLayers=1) -> Bool

Determines if a hit is detected by comparing a random value to the calculated efficiency.

# Arguments
- `hitTester::DetectorHitTester`: Instance with detector configuration.
- `energy::Float64`: Energy value of the hit.
- `glueLayers::Int=1`: Number of glue layers used in the efficiency calculation.

# Returns
- `Bool`: `true` if a hit is registered; otherwise `false`.
"""
function testHit(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int = 1)::Bool
    random_val = rand() * hitTester.scaleFactor
    return random_val < efficiency(hitTester, energy, glueLayers = glueLayers)
end
