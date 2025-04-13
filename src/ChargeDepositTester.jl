import SpecialFunctions: erfc

# Private helper functions
"""
Calculates the integral of a Gaussian distribution between specified bounds.

Uses the complementary error function (erfc) for numerical stability.
The result represents the area under the Gaussian curve between the specified bounds.

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
Generates a cache of Gaussian integral values for time binning.

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
Calculates the maximum number of bins that can contain a signal above threshold.

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
Struct representing a charge deposit tester with detector response parameters.

This struct contains all parameters necessary for realistic detector response simulation.
Parameters account for point spread, time resolution, and signal thresholds.

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
    #time_cache::Vector{Float64}
    #xf::Vector{Float64}
    #yf::Vector{Float64}
end

"""
Constructs a `ChargeDepositTester` object with specified or default values.

Uses global constants from SIGNAL and DETECTOR for configuration.
Calculates maximum affected bins in x and y directions.
Space step is set to 10% of the minimum detector dimension.

# Arguments
- `n_time::Int = 40`: Number of time bins in each direction.
- `n_space::Int = 100`: Number of spatial bins for simulation.
- `charge_sharing::Bool = false`: Enable charge sharing between pixels.
- `clustering::Bool = false`: Enable clustering of signals.

# Returns
- `ChargeDepositTester`: A new ChargeDepositTester object with calculated detector parameters.
"""
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
    )
end

# Functions operating on ChargeDepositTester

"""
Determines if a charge exceeds the detection threshold.

Uses Gaussian spreading to calculate the fraction of charge collected in the given range.
Multiplies by gain and quantum efficiency to determine if signal exceeds threshold.

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
Determines if a charge at a given position exceeds the detection threshold.

Calculates 2D Gaussian integral to determine charge fraction in the specified region.
Accounts for both x and y dimensions using separable Gaussian distributions.

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
Calculates the signal produced by a charge at a given position.

Similar to charge_over_threshold but returns the actual signal value rather than a boolean.
Uses 2D Gaussian spreading to calculate signal contribution in the specified region.

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
Applies Gaussian smearing to a time value based on the detector's time resolution.

Adds random Gaussian noise using the detector's time resolution parameter.
Models the detector's finite time resolution capabilities.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `time::Float64`: The original time value to be smeared.

# Returns
- `Float64`: The smeared time value.
"""
function smear_time(cdt::ChargeDepositTester, time::Float64)::Float64
    # Apply Gaussian smearing to the time
    smeared_time = time + randn() * cdt.resolution
    return smeared_time
end

"""
Applies Gaussian smearing to the time component of a hit coordinate.

Preserves the spatial coordinates while applying time smearing.
Creates a new HitCoordinate object rather than modifying the original.
Time smearing follows a Gaussian distribution with width determined by the detector's time resolution.

# Arguments
- `cdt::ChargeDepositTester`: The charge deposit tester configuration.
- `hit::HitCoordinate`: The original hit coordinate with time to be smeared.

# Returns
- `HitCoordinate`: A new hit coordinate with the same position but smeared time.
"""
function smear_time(cdt::ChargeDepositTester, hit::HitCoordinate)::HitCoordinate
    # Apply Gaussian smearing to the time
    smeared_time = hit.time + randn() * cdt.resolution
    return HitCoordinate(hit.x, hit.y, smeared_time)
end
