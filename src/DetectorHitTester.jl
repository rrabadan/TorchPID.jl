"""
Struct representing a detector hit testing configuration with various efficiency parameters.

Efficiency calculations combine multiple effects based on the enabled parameters
Energy thresholds define the valid range for photon detection

# Fields
- `scaleFactor::Float64`: Factor applied to random value in hit testing.
- `emin::Float64`: Minimum energy threshold for photon detection.
- `emax::Float64`: Maximum energy threshold for photon detection.
- `implement_epotek_cutoff::Bool`: Enable Epotek cutoff in efficiency calculation.
- `implement_imperfect_mirror::Bool`: Enable imperfect mirror behavior in efficiency calculation.
- `implement_QE::Bool`: Enable quantum efficiency in efficiency calculation.
- `implement_CE::Bool`: Enable collection efficiency in efficiency calculation.
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
Constructs a `DetectorHitTester` object with specified or default values.

# Keywords
- `scaleFactor::Float64=1.0`: Factor applied to random value in hit testing (default: 1.0).
- `emin::Float64=1.75`: Minimum energy threshold for photon detection (default: 1.75).
- `emax::Float64=7.00`: Maximum energy threshold for photon detection (default: 7.00).
- `implement_epotek_cutoff::Bool=true`: Enable Epotek cutoff in efficiency calculation (default: true).
- `implement_imperfect_mirror::Bool=true`: Enable imperfect mirror behavior in efficiency calculation (default: true).
- `implement_QE::Bool=true`: Enable quantum efficiency in efficiency calculation (default: true).
- `implement_CE::Bool=true`: Enable collection efficiency in efficiency calculation (default: true).

# Returns
A new `DetectorHitTester` object.
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
Calculates the detection efficiency for a photon of specified energy.

The efficiency calculation combines various detector effects based on the enabled options
in the detector configuration, including quantum efficiency, collection efficiency, 
Epotek cutoff, and mirror reflectivity.

# Arguments
- `hitTester::DetectorHitTester`: Instance holding configuration parameters.
- `energy::Float64`: The energy value of the photon in appropriate units.
- `glueLayers::Int=1`: Number of glue layers affecting the transmission efficiency.

# Returns
- `Float64`: Calculated efficiency value between 0.0 and 1.0. Returns 0.0 when energy is outside the threshold bounds.
"""
function efficiency(
    hitTester::DetectorHitTester,
    energy::Float64;
    glueLayers::Int = 1,
)::Float64
    if energy < hitTester.emin || energy > hitTester.emax
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
Determines if a photon hit is detected based on its energy and detector efficiency.

The function generates a random value and scales it by the detector's scale factor.
This scaled random value is compared against the calculated efficiency to simulate the probabilistic nature of photon detection.
Efficiency calculation takes into account various detector parameters including quantum efficiency and mirror reflectivity.
Multiple glue layers affect the transmission efficiency proportionally.

# Arguments
- `hitTester::DetectorHitTester`: Instance with detector configuration parameters.
- `energy::Float64`: Energy value of the photon in appropriate units.
- `glueLayers::Int=1`: Number of glue layers used in the efficiency calculation.

# Returns
- `Bool`: `true` if a hit is registered (random value < efficiency); otherwise `false`.
"""
function testHit(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int = 1)::Bool
    random_val = rand() * hitTester.scaleFactor
    return random_val < efficiency(hitTester, energy, glueLayers = glueLayers)
end
