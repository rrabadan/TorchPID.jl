"""
# GeometryConstants.jl - TORCH Detector Geometry

This file defines the geometric components and physical parameters of the TORCH 
detector system used in LHCb.

This file 
includes definitions for the major physical components:

- `Radiator`: The quartz plate where Cherenkov photons are generated
- `Wedge`: The focusing block that directs photons
- `Focus`: The focusing optics configuration 
- `Detector`: The MCP-PMT photon detector array
- `Mask`: Optional masking geometry
- `SignalParameters`: Electronic signal detection parameters

These components form a complete geometric model of the TORCH detector system,
with default parameters matching the LHCb experimental setup. The file also provides
utility functions for photon propagation and detection.

## Physical Units

Unless otherwise specified:
- Distances are in millimeters (mm)
- Times are in nanoseconds (ns)
- Angles are in radians (rad)
"""

"""
Struct representing the TORCH's Radiator geometry.

# Fields
- `height::Float64`: The height of the radiator.
- `width::Float64`: The width of the radiator.
- `depth::Float64`: The depth of the radiator.
- `half_width::Float64`: Half of the width of the radiator.
"""
struct Radiator
    height::Float64
    width::Float64
    depth::Float64
    half_width::Float64
end

"""
Constructs a `Radiator` instance.

# Arguments
- `height::Float64`: The height of the radiator (default: 2500.0).
- `width::Float64`: The width of the radiator (default: 660.0).
- `depth::Float64`: The depth of the radiator (default: 10.0).

# Returns
- `Radiator`: An instance of the `Radiator` struct.
"""
function Radiator(; height::Float64 = 2500.0, width::Float64 = 660.0, depth::Float64 = 10.0)
    Radiator(height, width, depth, width / 2)
end

"""
Struct representing the TORCH's Wedge geometry.

# Fields
- `theta::Float64`: The angle of the wedge in radians.
- `sin_theta::Float64`: The sine of the angle `theta` for optimization.
- `cos_theta::Float64`: The cosine of the angle `theta` for optimization.
- `tan_theta::Float64`: The tangent of the angle `theta` for optimization.
- `offset::Float64`: Offset value calculated as radiator depth divided by tangent of `theta`.

The wedge is a critical component directing photons from the radiator to the focusing optics.
Trigonometric values are pre-calculated for computational efficiency.
"""
struct Wedge
    theta::Float64
    sin_theta::Float64
    cos_theta::Float64
    tan_theta::Float64
    offset::Float64
end

"""
Constructs a `Wedge` instance with parameters derived from `Radiator` geometry.

Pre-calculates and stores trigonometric values to improve computational efficiency.
The offset is calculated based on the radiator depth and wedge angle.

# Arguments
- `radiator::Radiator`: The radiator geometry providing dimensional constraints.
- `theta::Float64=0.62831853`: The angle of the wedge in radians (default: 0.62831853, approximately π/5).

# Returns
- `Wedge`: An instance of the `Wedge` struct with calculated trigonometric values.
"""
function Wedge(radiator::Radiator; theta::Float64 = 0.62831853)
    Wedge(theta, sin(theta), cos(theta), tan(theta), radiator.depth / tan(theta))
end

"""
Struct representing the TORCH's Focus block geometry.

# Fields
- `theta_min::Float64`: Minimum angle for the focus.
- `theta_max::Float64`: Maximum angle for the focus.
- `tan_theta_min::Float64`: Tangent of the minimum angle.
- `tan_theta_max::Float64`: Tangent of the maximum angle.
- `theta::Float64`: Derived angle for the focus.
- `sin_theta::Float64`: Sine of the derived angle.
- `cos_theta::Float64`: Cosine of the derived angle.
- `radius::Float64`: Radius of the focus.
- `rsquared::Float64`: Square of the radius.
- `wedge::Float64`: Wedge parameter derived from the radiator and wedge.
- `wedge_gap::Float64`: Gap between wedges.
- `width::Float64`: Width of the focus.
- `entrance::Float64`: Entrance parameter for the focus.
- `axis_to_rhs::Float64`: Distance from axis to the right-hand side.
- `axis_to_lhs::Float64`: Distance from axis to the left-hand side.
- `axis_to_wedge::Float64`: Distance from axis to the wedge.
- `axis_to_detector::Float64`: Distance from axis to the detector.
- `air_gap::Float64`: Air gap parameter.
"""
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

"""
Constructs a `Focus` instance with parameters derived from `Radiator` and `Wedge` geometries.

# Arguments
- `radiator::Radiator`: The radiator geometry.
- `wedge::Wedge`: The wedge geometry.
- `theta_min::Float64`: Minimum angle for the focus (default: 0.45).
- `theta_max::Float64`: Maximum angle for the focus (default: 0.85).
- `radius::Float64`: Radius of the focus (default: 260.0).
- `wedge_gap::Float64`: Gap between wedges (default: 3.0).
- `width::Float64`: Width of the focus (default: 91.0).
- `entrance::Float64`: Entrance parameter (default: 133.7).
- `axis_to_rhs::Float64`: Distance from axis to the right-hand side (default: 26.704).
- `air_gap::Float64`: Air gap parameter (default: 0.0).

# Returns
- `Focus`: An instance of the `Focus` struct.
"""
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

"""
A struct represents the TORCH's Micro-channel Plate (MCP) Photomultiplier (PMT) Detector geometry and parameters.

# Fields
- `theta::Float64`: Angle of the detector.
- `middle::Float64`: Middle position of the detector.
- `width::Float64`: Width of the detector.
- `active::Float64`: Active area of the detector.
- `gap::Float64`: Gap in the detector.
- `n_xpixels::Int`: Number of pixels in the x-direction.
- `n_ypixels::Int`: Number of pixels in the y-direction.
- `n_detectors::Int`: Number of detectors.
- `n_xtotal::Int`: Total number of pixels in the x-direction.
- `x_size::Float64`: Size of a pixel in the x-direction.
- `y_size::Float64`: Size of a pixel in the y-direction.
- `x_min::Float64`: Minimum x-coordinate.
- `x_max::Float64`: Maximum x-coordinate.
- `y_min::Float64`: Minimum y-coordinate.
- `y_max::Float64`: Maximum y-coordinate.
- `t_min::Float64`: Minimum time window.
- `t_max::Float64`: Maximum time window.
- `t_bin::Float64`: Time bin size.
- `n_tpixels::Int`: Number of time pixels.
- `t_deadtime::Float64`: Dead time for the detector.
- `n_deadtime::Int`: Number of dead time bins.
- `t_window_min::Float64`: Minimum time window including dead time.
- `t_window_max::Float64`: Maximum time window.
- `n_total_pixels::Int`: Total number of pixels.
"""
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

"""
Constructs a `Detector` instance with parameters derived from `Radiator` and `Wedge` geometries.

# Arguments
- `radiator::Radiator`: The radiator geometry.
- `wedge::Wedge`: The wedge geometry.
- `middle::Float64`: Middle position of the detector (default: 34.405).
- `width::Float64`: Width of the detector (default: 60.0).
- `active::Float64`: Active area of the detector (default: 53.0).
- `n_xpixels::Int`: Number of pixels in the x-direction (default: 8).
- `n_ypixels::Int`: Number of pixels in the y-direction (default: 128).
- `n_detectors::Int`: Number of detectors (default: 11).
- `t_min::Float64`: Minimum time window (default: 40.0).
- `t_max::Float64`: Maximum time window (default: 65.0).
- `t_bin::Float64`: Time bin size (default: 0.01).
- `t_deadtime::Float64`: Dead time for the detector (default: 10.0).

# Returns
- `Detector`: An instance of the `Detector` struct.
"""
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

"""
A struct representing the TORCH's Mask.

# Fields
- `y_min::Float64`: Minimum y-coordinate for the mask.
- `y_max::Float64`: Maximum y-coordinate for the mask.
"""
struct Mask
    y_min::Float64
    y_max::Float64
end

"""
Constructs a `Mask` instance.

# Arguments
- `y_min::Float64`: Minimum y-coordinate for the mask (default: -999.0).
- `y_max::Float64`: Maximum y-coordinate for the mask (default: -999.0).

# Returns
- `Mask`: A new instance of the `Mask` struct.
"""
function Mask(; y_min::Float64 = -999.0, y_max::Float64 = -999.0)
    Mask(y_min, y_max)
end

"""
A struct representing TORCH's signal detection parameters.

# Fields
- `gain::Float64`: Gain of the signal.
- `threshold::Float64`: Threshold for the signal.
- `point_spread::Float64`: Point spread of the signal.
- `time_resolution::Float64`: Time resolution of the signal.
- `effective_constant::Float64`: Effective constant parameter.
- `effective_linear::Float64`: Effective linear parameter.
"""
struct SignalParameters
    gain::Float64
    threshold::Float64
    point_spread::Float64
    time_resolution::Float64
    effective_constant::Float64
    effective_linear::Float64
end

"""
Constructs a `SignalParameters` instance.

# Arguments
- `gain::Float64`: Gain of the signal (default: 6.5e5).
- `threshold::Float64`: Threshold for the signal (default: 3.0e-14).
- `point_spread::Float64`: Point spread of the signal (default: 0.8).
- `time_resolution::Float64`: Time resolution of the signal (default: 0.055).
- `effective_constant::Float64`: Effective constant parameter (default: 0.060).
- `effective_linear::Float64`: Effective linear parameter (default: 1.4e-5).

# Returns
- `SignalParameters`: An instance of the `SignalParameters` struct.
"""
function SignalParameters(;
    gain::Float64 = 6.5e5,
    threshold::Float64 = 3.0e-14,     # fC
    point_spread::Float64 = 0.8,      # mm
    time_resolution::Float64 = 0.055, # 55 ps
    effective_constant::Float64 = 0.060,
    effective_linear::Float64 = 1.4e-5,
)
    SignalParameters(
        gain,
        threshold,
        point_spread,
        time_resolution,
        effective_constant,
        effective_linear,
    )
end

# Create instances of the structs and export them
RADIATOR = Radiator(height = 2500.0, width = 660.0, depth = 10.0)
WEDGE = Wedge(RADIATOR)
FOCUS = Focus(RADIATOR, WEDGE)
DETECTOR = Detector(RADIATOR, WEDGE)
MASK = Mask()
SIGNAL = SignalParameters()

"""
Determines if a photon originates from the focus region.

# Arguments
- `yemission::Float64`: The y-coordinate of the photon's emission.

# Returns
- `Bool`: `true` if the photon originates from the focus region, `false` otherwise.
"""
function photon_from_focus(yemission)
    radiator_top = 0.5 * RADIATOR.height - WEDGE.height
    return yemission > radiator_top
end

"""
Checks if a photon hits the detector.

# Arguments
- `x::Float64`: The x-coordinate of the photon.
- `y::Float64`: The y-coordinate of the photon.

# Returns
- `Bool`: `true` if the photon hits the detector, `false` otherwise.
"""
function photon_on_detector(x, y)
    return (
        y > DETECTOR.y_min &&
        y < DETECTOR.y_max &&
        x > DETECTOR.x_min &&
        x < DETECTOR.x_max
    )
end
