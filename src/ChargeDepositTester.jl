using SpecialFunctions

# Private helper functions
function _gaussian_integral(v::Float64, s::Float64, min::Float64, max::Float64)::Float64
    result = 0.5 * (erfc((min - v) / (sqrt(2.0) * s)) - erfc((max - v) / (sqrt(2.0) * s)))
    return result
end

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

# Type definition and constructor
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

function get_smeared_time(cdt::ChargeDepositTester, time::Float64)::Float64
    # Apply Gaussian smearing to the time
    smeared_time = time + randn() * cdt.resolution
    return smeared_time
end

function smear_time(cdt::ChargeDepositTester, hit::HitCoordinate)::HitCoordinate
    # Apply Gaussian smearing to the time
    smeared_time = hit.time + randn() * cdt.resolution
    return HitCoordinate(hit.x, hit.y, smeared_time)
end
