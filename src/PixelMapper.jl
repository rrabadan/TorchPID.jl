function get_mcp_number(x::Float64)::Int
    xpos = x + 0.5 * RADIATOR.width
    return floor(Int, xpos / DETECTOR.width)
end

function get_mcp_xpixel(x::Float64)::Int
    xpos = x + 0.5 * RADIATOR.width
    mcp = floor(Int, xpos / DETECTOR.width)

    xpos -= mcp * DETECTOR.width
    xpos -= DETECTOR.gap

    return floor(Int, xpos / DETECTOR.x_size)
end

function get_ypixel(y::Float64)::Int
    pixel_num = floor(Int, (y - DETECTOR.y_min) / DETECTOR.y_size)
    if (pixel_num < 0) || (pixel_num >= DETECTOR.n_ypixels)
        return -1
    end
    return pixel_num
end

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

function get_mcp_xposition(x::Float64)::Float64
    xpos = x + 0.5 * RADIATOR.width
    mcp = floor(Int, DETECTOR.width)
    xpos -= (0.5 + mcp) * DETECTOR.width
    return xpos
end

function get_tpixel_nocheck(t::Float64)::Int
    # Get the time pixel number without checking bounds
    tpixel = floor(Int, t - DETECTOR.t_min) / DETECTOR.t_bin
    return tpixel
end

function get_tpixel(t::Float64)::Int
    pixelnum = get_tpixel_nocheck(t)
    # Check if pixel is in the gap
    if (pixelnum < 0) || (pixelnum >= DETECTOR.n_tpixels)
        return -1
    end
    return pixelnum
end
