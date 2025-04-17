"""
    get_mcp_number(x)

`get_mcp_number` calculates the MCP (Micro-Channel Plate) number corresponding to a given x-coordinate. 
It offsets the x-coordinate by half the radiator width and determines the MCP number by dividing the adjusted position by the detector width.

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
    get_mcp_xpixel(x)
    
`get_mcp_xpixel` determines the x-pixel number within a specific MCP based on the x-coordinate. 
The function identifies the MCP containing the coordinate, calculates the local position within that MCP while accounting for gaps, 
and computes the pixel number by dividing the local position by the x-pixel size.

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
    get_mcp_ypixel(y)

`get_mcp_ypixel` calculates the y-pixel number within a specific MCP from the y-coordinate. 
It determines the pixel number relative to the minimum y-value in the detector coordinate system. 
If the computed pixel number is outside the valid y-pixel range, the function returns -1.

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
    get_xpixel(x)

`get_xpixel` computes the global x-pixel number from the x-coordinate. 
It identifies the detector, calculates the local pixel position, and converts it to a global pixel number. 
Returns -1 if the detector or pixel is invalid.

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
    get_mcp_xposition(x)

`get_mcp_xposition` calculates the x-position relative to the center of the MCP. 
The function offsets the x-coordinate by half the radiator width, determines the MCP, 
and adjusts the position relative to the center of that MCP.

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
    get_tpixel_nocheck(t)

`get_tpixel_nocheck` computes the time pixel number from a given time value without performing bounds checking. 
It determines the pixel number based on the time value's offset from the minimum time value and the time bin size.

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
    get_tpixel(t)

`get_tpixel` computes the time pixel number from a given time value, ensuring it falls within valid bounds. 
It leverages `get_tpixel_nocheck` for the initial calculation and returns -1 if the result is outside the detector's time pixel range.

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
