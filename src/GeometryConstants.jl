# Radiator struct
struct Radiator
    height::Float64
    width::Float64
    depth::Float64
    half_width::Float64
end

# constructor computing derived field half_width
function Radiator(; height::Float64 = 2500.0, width::Float64 = 660.0, depth::Float64 = 10.0)
    Radiator(height, width, depth, width / 2)
end

struct Wedge
    theta::Float64
    sin_theta::Float64
    cos_theta::Float64
    tan_theta::Float64
    offset::Float64
end

# Outer constructor with default value for theta
function Wedge(radiator::Radiator; theta::Float64 = 0.62831853)
    Wedge(theta, sin(theta), cos(theta), tan(theta), radiator.depth / tan(theta))
end

# Focus: computed fields are derived from the Radiator and Wedge structs
struct Focus
    theta_min::Float64
    theta_max::Float64
    tan_theta_min::Float64
    tan_theta_max::Float64
    theta::Float64
    sin_theta::Float64
    cos_theta::Float64
    radius::Float64
    rsquared::Float64
    wedge::Float64
    wedge_gap::Float64
    width::Float64
    entrance::Float64
    axis_to_rhs::Float64
    axis_to_lhs::Float64
    axis_to_wedge::Float64
    axis_to_detector::Float64
    air_gap::Float64
end

# Outer constructor that computes derived fields from the Radiator and Wedge structs
function Focus(
    radiator::Radiator,
    wedge::Wedge,
    theta_min::Float64 = 0.45,
    theta_max::Float64 = 0.85,
    radius::Float64 = 260.0,
    wedge_gap::Float64 = 3.0,
    width::Float64 = 91.0,
    entrance::Float64 = 133.7,
    axis_to_rhs::Float64 = 26.704,
    air_gap::Float64 = 0.0,
)
    tan_theta_min = tan(theta_min)
    tan_theta_max = tan(theta_max)
    theta = pi / 2.0 - wedge.theta
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    rsquared = radius^2
    axis_to_lhs = width - axis_to_rhs
    axis_to_wedge = axis_to_rhs - wedge_gap
    axis_to_detector = 64.0 - axis_to_rhs
    Focus(
        theta_min,
        theta_max,
        tan_theta_min,
        tan_theta_max,
        theta,
        sin_theta,
        cos_theta,
        radius,
        rsquared,
        radiator.depth / wedge.sin_theta,
        wedge_gap,
        width,
        entrance,
        axis_to_rhs,
        axis_to_lhs,
        axis_to_wedge,
        axis_to_detector,
        air_gap,
    )
end

struct Detector
    theta::Float64
    middle::Float64
    width::Float64
    active::Float64
    gap::Float64
    n_xpixels::Int
    n_ypixels::Int
    n_detectors::Int
    n_xtotal::Int
    x_size::Float64
    y_size::Float64
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    t_min::Float64
    t_max::Float64
    t_bin::Float64
    n_tpixels::Int
    t_deadtime::Float64
    n_deadtime::Int
    t_window_min::Float64
    t_window_max::Float64
    n_total_pixels::Int
end

function Detector(
    radiator::Radiator,
    wedge::Wedge,
    middle::Float64 = 34.405,
    width::Float64 = 60.0,
    active::Float64 = 53.0,
    n_xpixels::Int = 8,
    n_ypixels::Int = 128,
    n_detectors::Int = 11,
    t_min::Float64 = 40.0,
    t_max::Float64 = 65.0,
    t_bin::Float64 = 0.01,
    t_deadtime::Float64 = 10.0,
)
    theta = wedge.theta
    gap = 0.5 * (width - active)
    n_xtotal = n_xpixels * n_detectors
    x_size = active / n_xpixels
    y_size = width / n_ypixels
    x_min = -radiator.half_width
    x_max = +radiator.half_width
    y_min = -0.5 * active
    y_max = +0.5 * active
    n_tpixels = Int((t_max - t_min) / t_bin)
    n_deadtime = Int(t_deadtime / t_bin)
    t_window_min = t_min - t_deadtime - 1.0
    t_window_max = t_max + 1.0
    n_total_pixels = n_xtotal * n_ypixels * n_tpixels
    Detector(
        theta,
        middle,
        width,
        active,
        gap,
        n_xpixels,
        n_ypixels,
        n_detectors,
        n_xtotal,
        x_size,
        y_size,
        x_min,
        x_max,
        y_min,
        y_max,
        t_min,
        t_max,
        t_bin,
        n_tpixels,
        t_deadtime,
        n_deadtime,
        t_window_min,
        t_window_max,
        n_total_pixels,
    )
end

struct Mask
    y_min::Float64
    y_max::Float64
end

function Mask(y_min::Float64 = -999.0, y_max::Float64 = -999.0)
    Mask(y_min, y_max)
end

# Create instances of the structs and export them
RADIATOR = Radiator(height = 2500.0, width = 660.0, depth = 10.0)
WEDGE = Wedge(RADIATOR)
FOCUS = Focus(RADIATOR, WEDGE)
DETECTOR = Detector(RADIATOR, WEDGE)
MASK = Mask()


function photonFromFocus(yemission)
    radiator_top = 0.5 * RADIATOR.height - WEDGE.height
    return yemission > radiator_top
end

function photonOnDetector(x, y)
    return (
        y > DETECTOR.y_min &&
        y < DETECTOR.y_max &&
        x > DETECTOR.x_min &&
        x < DETECTOR.x_max
    )
end
