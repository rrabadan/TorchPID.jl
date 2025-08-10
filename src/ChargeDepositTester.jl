import SpecialFunctions: erfc

"""
    gaussian_integral(v, s, min, max)

`_gaussian_integral` is a helper function that calculates the integral of a Gaussian distribution between specified bounds.
It leverages the complementary error function (erfc) for numerical stability. The result represents the area under the Gaussian curve within the given range.

# Arguments
- `v::Float64`: Mean value of the Gaussian distribution.
- `s::Float64`: Standard deviation of the Gaussian distribution.
- `min::Float64`: Lower bound of the integration range.
- `max::Float64`: Upper bound of the integration range.

# Returns
- `Float64`: The probability of a value from the Gaussian distribution falling between min and max.
"""
function _gaussian_integral(v::Float64, s::Float64, min::Float64, max::Float64)::Float64
    result = 0.5 * (erfc((min - v) / (sqrt(2.0) * s)) - erfc((max - v) / (sqrt(2.0) * s)))
    return result
end

"""
    _set_resolution(resolution)

`_set_resolution` is a helper (internal) function that generates a cache of Gaussian integral values for time binning.
Creates a cache for -n_time to n_time bins centered around zero.
Each value represents the probability of a hit falling within that time bin.

# Arguments
- `resolution::Float64`: Time resolution parameter for the Gaussian smearing.

# Returns
- `Vector{Float64}`: Vector containing probability values for time bins.
"""
function _set_resolution(resolution::Float64)::Vector{Float64}
    # time_cache = Vector{Float64}(undef, 2 * n_time + 1)
    time_cache = Vector{Float64}(undef, 0)

    for i = -n_time:n_time
        p = gaussian_integral(
            0.0,
            resolution,
            (i - 0.5) * DETECTOR.t_bin,
            (i + 0.5) * DETECTOR.t_bin,
        )
        push!(time_cache, p)
    end
    return time_cache
end

"""
    _get_max_bins(a, b, n, point_spread, gain, threshold)

`_get_max_bins` calculates the maximum number of bins that can contain a signal above threshold.
Iteratively calculates signal in neighboring bins until it falls below threshold.
Used to determine how many neighboring pixels might be activated by a single hit.

# Arguments
- `a::Float64`: First dimension size (typically detector bin size).
- `b::Float64`: Second dimension size (typically detector bin size).
- `n::Int`: Maximum number of bins to check.
- `point_spread::Float64`: Standard deviation of the point spread function.
- `gain::Float64`: Signal gain factor.
- `threshold::Float64`: Detection threshold value.

# Returns
- `Int`: Maximum number of bins that can have a signal above threshold.
"""
function _get_max_bins(
    a::Float64,
    b::Float64,
    n::Int,
    point_spread::Float64,
    gain::Float64,
    threshold::Float64,
)::Int
    for i = 0:(n-1)
        p =
            _gaussian_integral(0.0, point_spread, i * a, (i + 1) * a) *
            _gaussian_integral(0.5 * b, point_spread, 0.0, b)
        if p * gain * QQe < threshold
            return i
        end
    end
    return n
end

"""
    ChargeDepositTester(gain, threshold, point_spread, time_resolution,
                        space_step, n_time, n_space, n_xbins, n_ybins,
                        charge_sharing, clustering)

Represents a charge deposit tester with parameters for simulating detector response. 
This struct encapsulates all necessary configuration details, including point spread, 
time resolution, signal thresholds, and spatial-temporal binning, for realistic detector modeling.

# Fields
- `gain::Float64`: Signal amplification factor.
- `threshold::Float64`: Minimum signal threshold for detection.
- `point_spread::Float64`: Standard deviation of the point spread function.
- `time_resolution::Float64`: Standard deviation of the time resolution in ns.
- `space_step::Float64`: Spatial step size for simulation in mm.
- `n_time::Int`: Number of time bins in each direction.
- `n_space::Int`: Number of spatial bins in each direction.
- `n_xbins::Int`: Maximum number of x bins that can have signal above threshold.
- `n_ybins::Int`: Maximum number of y bins that can have signal above threshold.
- `charge_sharing::Bool`: Whether charge sharing between pixels is enabled.
- `clustering::Bool`: Whether clustering of signals is enabled.

# Constructors

    ChargeDepositTester(;n_time=40, n_space=100, charge_sharing=false, clustering=false)

Creates a `ChargeDepositTester` instance using global constants from `SIGNAL` and `DETECTOR` for configuration.
Determines the maximum number of bins affected in the x and y directions based on the point spread and threshold.
Sets the spatial step size to 10% of the smaller detector dimension.

# Keywords
- `n_time::Int = 40`: Number of time bins in each direction.
- `n_space::Int = 100`: Number of spatial bins for simulation.
- `charge_sharing::Bool = false`: Enable charge sharing between pixels.
- `clustering::Bool = false`: Enable clustering of signals.

# Returns
- `ChargeDepositTester`: A new ChargeDepositTester object with calculated detector parameters.
"""
struct ChargeDepositTester
    gain::Float64
    threshold::Float64
    point_spread::Float64
    time_resolution::Float64
    space_step::Float64
    n_time::Int
    n_space::Int
    n_xbins::Int
    n_ybins::Int
    charge_sharing::Bool
    clustering::Bool
    time_cache::Vector{Float64}
    xf::Vector{Float64}
    yf::Vector{Float64}
end

function ChargeDepositTester(;
    n_time::Int = 40,
    n_space::Int = 100,
    charge_sharing::Bool = false,
    clustering::Bool = false,
)
    gain = SIGNAL.gain
    threshold = SIGNAL.threshold
    point_spread = SIGNAL.point_spread
    resolution = SIGNAL.time_resolution

    space_step = 0.1 * min(DETECTOR.x_size, DETECTOR.y_size)

    xf = zeros(Float64, n_space)
    yf = zeros(Float64, n_space)

    for i = 1:n_space
        xf[i] =
            _gaussian_integral(space_step * (i + 0.5), point_spread, 0.0, DETECTOR.x_size)
        yf[i] =
            _gaussian_integral(space_step * (i + 0.5), point_spread, 0.0, DETECTOR.y_size)
    end

    n_xbins = _get_max_bins(
        DETECTOR.x_size,
        DETECTOR.y_size,
        DETECTOR.n_xpixels,
        point_spread,
        gain,
        threshold,
    )
    n_ybins = _get_max_bins(
        DETECTOR.y_size,
        DETECTOR.x_size,
        DETECTOR.n_ypixels,
        point_spread,
        gain,
        threshold,
    )

    # time template
    time_cache = zeros(Float64, 2 * n_time + 1)
    for i = -n_time:n_time
        p = _gaussian_integral(
            0.0,
            resolution,
            (i - 0.5) * DETECTOR.t_bin,
            (i + 0.5) * DETECTOR.t_bin,
        )
        time_cache[i+n_time+1] = p
    end

    ChargeDepositTester(
        gain,
        threshold,
        point_spread,
        resolution,
        space_step,
        n_time,
        n_space,
        n_xbins,
        n_ybins,
        charge_sharing,
        clustering,
        time_cache,
        xf,
        yf,
    )
end

# Functions operating on ChargeDepositTester
"""
    charge_over_threshold(cdt, charge, min, max)

`charge_over_threshold` checks if a charge exceeds the detection threshold.
It calculates the fraction of charge collected in the specified range using Gaussian spreading.
The result is scaled by the gain and quantum efficiency to determine if the signal surpasses the threshold.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `charge::Float64`: The charge value to test.
- `min::Float64`: Lower bound of the integration range.
- `max::Float64`: Upper bound of the integration range.

# Returns
- `Bool`: `true` if the charge produces a signal above threshold, otherwise `false`.
"""
function charge_over_threshold(
    cdt::ChargeDepositTester,
    charge::Float64,
    min::Float64,
    max::Float64,
)::Bool
    f = _gaussian_integral(charge, cdt.point_spread, min, max)
    signal = f * cdt.gain * QQe
    return signal > cdt.threshold
end

"""
    charge_over_threshold(cdt, x, xmin, xmax, y, ymin, ymax)

`charge_over_threshold` determines whether a charge at a specific position surpasses the detection threshold.
It computes the 2D Gaussian integral to evaluate the charge fraction within the defined region, considering both x and y dimensions independently using separable Gaussian distributions.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `x::Float64`: X coordinate of the charge.
- `xmin::Float64`: Minimum x coordinate of the detection region.
- `xmax::Float64`: Maximum x coordinate of the detection region.
- `y::Float64`: Y coordinate of the charge.
- `ymin::Float64`: Minimum y coordinate of the detection region.
- `ymax::Float64`: Maximum y coordinate of the detection region.

# Returns
- `Bool`: `true` if the charge produces a signal above threshold, otherwise `false`.
"""
function charge_over_threshold(
    cdt::ChargeDepositTester,
    x::Float64,
    xmin::Float64,
    xmax::Float64,
    y::Float64,
    ymin::Float64,
    ymax::Float64,
)::Bool
    xf = _gaussian_integral(x, cdt.point_spread, xmin, xmax)
    yf = _gaussian_integral(y, cdt.point_spread, ymin, ymax)
    signal = xf * yf * cdt.gain * QQe
    return signal > cdt.threshold
end

"""
    get_charge(cdt, charge, min, max)

`get_charge` computes the signal generated by a charge within a specified region.
Unlike `charge_over_threshold`, it returns the actual signal value instead of a boolean.
The calculation uses 2D Gaussian spreading to determine the signal contribution in the defined range.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `x::Float64`: X coordinate of the charge.
- `xmin::Float64`: Minimum x coordinate of the detection region.
- `xmax::Float64`: Maximum x coordinate of the detection region.
- `y::Float64`: Y coordinate of the charge.
- `ymin::Float64`: Minimum y coordinate of the detection region.
- `ymax::Float64`: Maximum y coordinate of the detection region.

# Returns
- `Float64`: The signal value produced by the charge.
"""
function get_charge(
    cdt::ChargeDepositTester,
    x::Float64,
    xmin::Float64,
    xmax::Float64,
    y::Float64,
    ymin::Float64,
    ymax::Float64,
)::Float64
    xf = _gaussian_integral(x, cdt.point_spread, xmin, xmax)
    yf = _gaussian_integral(y, cdt.point_spread, ymin, ymax)
    signal = xf * yf * cdt.gain * QQe
    return signal
end

"""
    smear_time(cdt::ChargeDepositTester, time::Float64)
    smear_time(cdt::ChargeDepositTester, hit::HitCoordinate)

`smear_time` introduces Gaussian smearing to a time value, simulating the detector's finite time resolution. 
It adds random noise drawn from a Gaussian distribution with a standard deviation equal to the detector's time resolution.

# Methods
1. `smear_time(cdt::ChargeDepositTester, time::Float64)`: 
   Takes a time value and returns the smeared time value.

2. `smear_time(cdt::ChargeDepositTester, hit::HitCoordinate)`:
   Takes a `HitCoordinate` object and returns a new `HitCoordinate` with the time value smeared.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `time::Float64`: The original time value to be smeared (for the first method).
- `hit::HitCoordinate`: The original hit coordinate with time to be smeared (for the second method).

# Returns
- `Float64`: The smeared time value (for the first method).
- `HitCoordinate`: A new `HitCoordinate` with smeared time (for the second method).
"""
function smear_time(cdt::ChargeDepositTester, time::Float64)::Float64
    # Apply Gaussian smearing to the time
    smeared_time = time + randn() * cdt.resolution
    return smeared_time
end

function smear_time(cdt::ChargeDepositTester, hit::HitCoordinate)::HitCoordinate
    # Apply Gaussian smearing to the time
    smeared_time = hit.time + randn() * cdt.resolution
    return HitCoordinate(hit.x, hit.y, smeared_time)
end

"""
    cached_time_weights(cdt::ChargeDepositTester, start::Integer)

Returns a view of the `time_cache` from the `start` index to the end.
This function is an efficient way to get a sub-array of time weights without allocating new memory.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `start::Integer`: The starting index for the view (1-based).

# Returns
- `SubArray`: A view of the `time_cache` vector.
"""
function cached_time_weights(cdt::ChargeDepositTester, start::Integer)::SubArray{Float64,1}
    return @view cdt.time_cache[start:end]
end
