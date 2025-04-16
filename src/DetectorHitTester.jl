"""
    DetectorHitTester

Type representing the configuration for photon detection efficiency testing with adjustable parameters.
Combines multiple effects such as quantum efficiency, collection efficiency, and mirror reflectivity.
to calculate detection probabilities. Energy thresholds define the valid photon detection range.

# Fields
- `scaleFactor::Float64`: Factor applied to random value in hit testing.
- `emin::Float64`: Minimum energy threshold for photon detection.
- `emax::Float64`: Maximum energy threshold for photon detection.
- `implement_epotek_cutoff::Bool`: Enable Epotek cutoff in efficiency calculation.
- `implement_imperfect_mirror::Bool`: Enable imperfect mirror behavior in efficiency calculation.
- `implement_QE::Bool`: Enable quantum efficiency in efficiency calculation.
- `implement_CE::Bool`: Enable collection efficiency in efficiency calculation.

# Constructors

    DetectorHitTester(;
        scaleFactor::Float64 = 1.0,
        emin::Float64 = 1.75,
        emax::Float64 = 7.00,
        implement_epotek_cutoff::Bool = true,
        implement_imperfect_mirror::Bool = true,
        implement_QE::Bool = true,
        implement_CE::Bool = true,
    )

Constructs a `DetectorHitTester` object with specified or default values.

## Keywords
- `scaleFactor::Float64=1.0`: Factor applied to random value in hit testing (default: 1.0).
- `emin::Float64=1.75`: Minimum energy threshold for photon detection (default: 1.75).
- `emax::Float64=7.00`: Maximum energy threshold for photon detection (default: 7.00).
- `implement_epotek_cutoff::Bool=true`: Enable Epotek cutoff in efficiency calculation (default: true).
- `implement_imperfect_mirror::Bool=true`: Enable imperfect mirror behavior in efficiency calculation (default: true).
- `implement_QE::Bool=true`: Enable quantum efficiency in efficiency calculation (default: true).
- `implement_CE::Bool=true`: Enable collection efficiency in efficiency calculation (default: true).

# Examples
```julia
# Create a DetectorHitTester instance with default parameters
detector = DetectorHitTester()
# Create a DetectorHitTester instance with custom parameters
detector_custom = DetectorHitTester(
    scaleFactor = 1.5,
    emin = 2.0,
    emax = 6.0,
    implement_epotek_cutoff = false,
    implement_imperfect_mirror = true,
    implement_QE = true,
    implement_CE = false,
)
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
    photon_efficiency(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int=1)

`photon_efficiency` computes the photon detection efficiency for a given energy value. 
The calculation integrates multiple detector effects, such as quantum efficiency, 
collection efficiency, Epotek cutoff, and mirror reflectivity, based on the 
configuration options enabled in the `DetectorHitTester` instance.

# Arguments
- `hitTester::DetectorHitTester`: Instance holding configuration parameters.
- `energy::Float64`: The energy value of the photon in appropriate units.
- `glueLayers::Int=1`: Number of glue layers affecting the transmission efficiency.

# Returns
- `Float64`: Calculated efficiency value between 0.0 and 1.0. Returns 0.0 when energy is outside the threshold bounds.
"""
function photon_efficiency(
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
    test_photon(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int=1)

`test_photon` simulates the detection of a Cherenkov photon based on its energy and the configuration of `DetectorHitTester`. 
The function generates a random value, scales it by the detector's scale factor, and compares it to the calculated efficiency. 
The efficiency calculation incorporates various detector parameters, such as quantum efficiency, mirror reflectivity, and the impact of glue layers on transmission efficiency. 
This probabilistic approach models the likelihood of photon detection.

# Arguments
- `hitTester::DetectorHitTester`: Instance with detector configuration parameters.
- `energy::Float64`: Energy value of the photon in appropriate units.
- `glueLayers::Int=1`: Number of glue layers used in the efficiency calculation.

# Returns
- `Bool`: `true` if a hit is registered (random value < efficiency); otherwise `false`.
"""
function test_photon(
    hitTester::DetectorHitTester,
    energy::Float64;
    glueLayers::Int = 1,
)::Bool
    random_val = rand() * hitTester.scaleFactor
    return random_val < photon_efficiency(hitTester, energy, glueLayers = glueLayers)
end
