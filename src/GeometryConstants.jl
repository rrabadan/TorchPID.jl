# Radiator struct
struct Radiator
    height::Float64
    width::Float64
    depth::Float64
    half_width::Float64
end

# constructor computing derived field half_width
function Radiator(; height::Float64=2500.0, width::Float64=660.0, depth::Float64=10.0)
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
function Wedge(radiator::Radiator; theta::Float64=0.62831853)
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
function Focus(radiator::Radiator, wedge::Wedge,
    theta_min::Float64=0.45,
    theta_max::Float64=0.85,
    radius::Float64=260.0,
    wedge_gap::Float64=3.0,
    width::Float64=91.0,
    entrance::Float64=133.7,
    axis_to_rhs::Float64=26.704,
    air_gap::Float64=0.0,
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

# Create instances of the structs and export them
RADIATOR = Radiator(height=2500.0, width=660.0, depth=10.0)
WEDGE = Wedge(RADIATOR)
FOCUS = Focus(RADIATOR, WEDGE)
