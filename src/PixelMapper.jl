"""
Determines the MCP (Micro-Channel Plate) number based on the x-coordinate.

The function applies an offset of half the radiator width to the x-coordinate.
The MCP number is calculated by dividing the position by the detector width.

# Arguments
- `x::Float64`: The x-coordinate in the detector coordinate system.

# Returns
- `Int`: The MCP number.
"""
function get_mcp_number(x::Float64)::Int
    xpos = x + 0.5 * RADIATOR.width
    return floor(Int, xpos / DETECTOR.width)
end

"""
Calculates the x-pixel number within a specific MCP from the x-coordinate.

The function first calculates which MCP the coordinate belongs to.
It then calculates the local position within that MCP, accounting for gaps.
The pixel number is determined by dividing by the x-pixel size.

# Arguments
- `x::Float64`: The x-coordinate in the detector coordinate system.

# Returns
- `Int`: The x-pixel number within the MCP.
"""
function get_mcp_xpixel(x::Float64)::Int
    xpos = x + 0.5 * RADIATOR.width
    mcp = floor(Int, xpos / DETECTOR.width)

    xpos -= mcp * DETECTOR.width
    xpos -= DETECTOR.gap

    return floor(Int, xpos / DETECTOR.x_size)
end

"""
Calculates the y-pixel number from the y-coordinate.

The function calculates the pixel number based on the y-coordinate's position relative to the minimum y-value.
Returns -1 if the calculated pixel number is outside the detector's y-pixel range.

# Arguments
- `y::Float64`: The y-coordinate in the detector coordinate system.

# Returns
- `Int`: The y-pixel number, or -1 if the coordinate is outside the valid range.
"""
function get_ypixel(y::Float64)::Int
    pixel_num = floor(Int, (y - DETECTOR.y_min) / DETECTOR.y_size)
    if (pixel_num < 0) || (pixel_num >= DETECTOR.n_ypixels)
        return -1
    end
    return pixel_num
end

"""
Calculates the global x-pixel number from the x-coordinate.

First determines which detector the coordinate belongs to.
Returns -1 if the detector number is invalid.
Calculates the local position on the detector, accounting for gaps.
Returns -1 if the local pixel number is outside the valid range.
Converts the local pixel number to a global pixel number based on detector number.

# Arguments
- `x::Float64`: The x-coordinate in the detector coordinate system.

# Returns
- `Int`: The global x-pixel number, or -1 if the coordinate is outside the valid range.
"""
function get_xpixel(x::Float64)::Int
    # Get the detector number
    ndetector = floor(Int, (x + 0.5 * RADIATOR.width) / DETECTOR.width)

    if (ndetector < 0) || (ndetector >= DETECTOR.n_detectors)
        return -1
    end

    # Position on the detector based on the detector number
    xdetector = (x + 0.5 * RADIATOR.width - ndetector * DETECTOR.width)

    # Local pixel number within detector
    pixel_num = floor(Int, (xdetector - DETECTOR.gap) / DETECTOR.x_size)

    # Check if pixel is in the gap
    if (pixel_num < 0) || (pixel_num >= DETECTOR.n_xpixels)
        return -1
    end

    # Convert to a global number
    pixel_num += DETECTOR.n_xpixels * ndetector
    return pixel_num
end

"""
Calculates the x-position relative to the MCP center.

The function applies an offset of half the radiator width to the x-coordinate.
The MCP is determined, and the position is adjusted relative to that MCP's center.

# Arguments
- `x::Float64`: The x-coordinate in the detector coordinate system.

# Returns
- `Float64`: The x-position relative to the MCP center.
"""
function get_mcp_xposition(x::Float64)::Float64
    xpos = x + 0.5 * RADIATOR.width
    mcp = floor(Int, DETECTOR.width)
    xpos -= (0.5 + mcp) * DETECTOR.width
    return xpos
end

"""
Calculates the time pixel number from a time value without checking bounds.

The function calculates the pixel number based on the time value's position relative to the minimum time value.
No bounds checking is performed in this function.

# Arguments
- `t::Float64`: The time value.

# Returns
- `Int`: The time pixel number.
"""
function get_tpixel_nocheck(t::Float64)::Int
    # Get the time pixel number without checking bounds
    tpixel = floor(Int, t - DETECTOR.t_min) / DETECTOR.t_bin
    return tpixel
end

"""
Calculates the time pixel number from a time value with bounds checking.

Uses get_tpixel_nocheck to calculate the initial pixel number.
Returns -1 if the calculated pixel number is outside the detector's time pixel range.

# Arguments
- `t::Float64`: The time value.

# Returns
- `Int`: The time pixel number, or -1 if the time value is outside the valid range.
"""
function get_tpixel(t::Float64)::Int
    pixelnum = get_tpixel_nocheck(t)
    # Check if pixel is in the gap
    if (pixelnum < 0) || (pixelnum >= DETECTOR.n_tpixels)
        return -1
    end
    return pixelnum
end
