# Include dependencies
include("PixelMapper.jl")

# ===== Type Definitions =====
const MCPPixelMap = Vector{Vector{Vector{Tuple{Int,Float64}}}}

struct MCPImage
    n_xpixels::Int
    n_ypixels::Int
    n_tpixels::Int
    n_deadtime::Int
    pixelmap::MCPPixelMap
end

function MCPImage(n_xpixels::Int, n_ypixels::Int, n_tpixels::Int, n_deadtime::Int)::MCPImage
    pixelmap = [[Tuple{Int,Float64}[] for _ = 1:n_ypixels] for _ = 1:n_xpixels]
    MCPImage(n_xpixels, n_ypixels, n_tpixels, n_deadtime, pixelmap)
end

function _fill_hit!(
    mcpimage::MCPImage,
    xpixel::Int,
    ypixel::Int,
    tpixel::Int;
    charge::Float64 = -1.0,
)
    if (tpixel < -mcpimage.n_deadtime) ||
       (tpixel >= mcpimage.n_tpixels) ||
       (xpixel < 0) ||
       (xpixel >= mcpimage.n_xpixels) ||
       (ypixel < 0) ||
       (ypixel >= mcpimage.n_ypixels)
        return
    end
    pixelmap = mcpimage.pixelmap
    # Add the hit to the pixel map
    if xpixel + 1 > 0 &&
       xpixel + 1 <= length(pixelmap) &&
       ypixel + 1 > 0 &&
       ypixel + 1 <= length(pixelmap[xpixel+1])
        push!(pixelmap[xpixel+1][ypixel+1], (tpixel, charge))
    end
end

function _fill_hit!(mcpimage::MCPImage, hit::PixelHit)
    _fill_hit!(mcpimage, hit.xpixel, hit.ypixel, hit.tpixel; charge = hit.charge)
end

function _is_occupied(pixelmap::MCPPixelMap, xpixel::Int, ypixel::Int)

    if (xpixel < 0) ||
       (xpixel >= length(pixelmap)) ||
       (ypixel < 0) ||
       (ypixel >= length(pixelmap[1]))
        return false
    end

    # Check if the pixel is occupied
    return !isempty(pixelmap[xpixel+1][ypixel+1]) # && any(x -> x[1] == tpixel, pixelMap[xpixel][ypixel])
end

function _reset!(pixelmap::MCPPixelMap)
    # Reset the pixel map to empty
    for col in pixelmap
        for row in col
            empty!(row)
        end
    end
end

function _get_hits(
    mcpimage::MCPImage,
    xpixel::Int,
    ypixel::Int,
    mcp::Int,
)::Union{Nothing,Vector{PixelHit}}

    pixelmap = mcpimage.pixelmap

    # Account for 1-based indexing in Julia
    x = xpixel + 1
    y = ypixel + 1

    if isempty(pixelmap[x][y])
        # If the pixel is empty, return nothing
        return nothing
    end

    hitcollection = Vector{PixelHit}()

    offset = mcp * length(pixelmap)

    # Sort the hits in descending order by the first element of the tuple (assumes hit tuple is (tpixel, charge))
    sorted_hits = sort(pixelmap[x][y], by = first, rev = true)

    # If the first hit's t value is > 0, add it to the collection
    if sorted_tpixels[1][1] > 0
        push!(
            hitcollection,
            PixelHit(xpixel + offset, ypixel, sorted_hits[1][1], sorted_hits[1][2]),
        )
    end

    # Loop over the rest of the hits and add them if the gap to the previous hit is sufficient
    for t in Iterators.drop(eachindex(sorted_hits), 1)
        tPixel = sorted_hits[t][1]
        tPrev = sorted_hits[prevind(sorted_hits, t)][1]
        if (tPixel > 0) && (tPixel - tPrev > mcpimage.n_deadtime)
            push!(hitcollection, PixelHit(x + offset, y, tPixel, sorted_hits[t][2]))
        end
    end
    return hitcollection
end

"""
    get_hits(pixelmap::MCPPixelMap, mcpimage::MCPImage, mcp::Int)::Union{Nothing,Vector{PixelHit}}

Retrieve all hits from the pixel map for a given MCP (Microchannel Plate).

# Arguments
    hitcollection = PixelHit[]
- `mcpimage::MCPImage`: The MCP image metadata.
- `mcp::Int`: The index of the MCP to retrieve hits for.

# Returns
- A vector of `PixelHit` objects or `nothing` if no hits are found.
"""
function _get_hits(mcpimage::MCPImage, mcp::Int)::Union{Nothing,Vector{PixelHit}}
    pixelmap = mcpimage.pixelmap
    hitcollection = Vector{PixelHit}(undef, 0)
    for x = 1:mcpimage.n_xpixels
        for y = 1:mcpimage.n_ypixels
            if _is_occupied(pixelmap, x - 1, y - 1)
                pixelhits = _get_hits(mcpimage, x - 1, y - 1, mcp)
                if pixelhits !== nothing
                    append!(hitcollection, pixelhits)
                end
            end
        end
    end
    if isempty(hitcollection)
        return nothing
    else
        return hitcollection
    end
end

function _get_hits(
    mcpimage::MCPImage,
    xpixel::Int,
    mcp::Int,
)::Union{Nothing,Vector{PixelHit}}

    hitcollection = Vector{PixelHit}(undef, 0)
    if xpixel < 0 || xpixel >= mcpimage.n_xpixels
        return nothing
    end
    for y = 1:mcpimage.n_ypixels
        if _is_occupied(pixelmap, xpixel, y - 1)
            pixelhits = _get_hits(mcpimage, xpixel, y - 1, mcp)
            if pixelhits !== nothing
                append!(hitcollection, pixelhits)
            end
        end
    end
    if isempty(hitcollection)
        return nothing
    else
        return hitcollection
    end
end

"""
    get_number_hits(pixelmap::MCPPixelMap)::Int

Calculate the total number of hits in the given pixel map.

# Arguments
- `pixelmap::MCPPixelMap`: The pixel map containing hit data.

# Returns
- The total number of hits as an integer.
"""
function _get_number_hits(pixelmap::MCPPixelMap)::Int
    total_hits = 0
    # Iterate over all pixels in the pixelmap to count the total number of hits
    for x in pixelmap
        for y in x
            total_hits += length(y)
        end
    end
    return total_hits
end

import .PixelMapper:
    get_mcp_number,
    get_mcp_xpixel,
    get_ypixel,
    get_xpixel,
    get_mcp_xposition,
    get_tpixel_nocheck

# ===== FrontEnd Definition and Constructor =====
struct FrontEnd
    n_mcps::Int
    n_xpixels::Int
    n_ypixels::Int
    sharing_mode::Bool
    smear_time::Bool
    pixelmap::Vector{MCPImage}
end

const MCPImageArray = Vector{MCPImage}

function FrontEnd(;
    n_mcps::Int = DETECTOR.n_detectors,
    n_xpixels::Int = DETECTOR.n_xpixels,
    n_ypixels::Int = DETECTOR.n_ypixels,
    n_tpixels::Int = DETECTOR.n_tpixels,
    n_deadtime::Int = DETECTOR.n_deadtime,
    sharing_mode::Bool = false,
    smear_time::Bool = false,
)::FrontEnd
    pixelmap = Vector{MCPImage}(undef, n_mcps)
    for i = 1:n_mcps
        pixelmap[i] = MCPImage(n_xpixels, n_ypixels, n_tpixels, n_deadtime)
    end
    FrontEnd(n_mcps, n_xpixels, n_ypixels, sharing_mode, smear_time, pixelmap)
end

function find_pixels(
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    xpos::Float64,
    ypos::Float64,
)::Union{Nothing,Vector{PixelHit}}
    mcp = get_mcp_number(xpos)
    xcentre = get_mcp_xpixel(xpos)
    ycentre = get_ypixel(ypos)
    if (mcp < 0) ||
       (mcp >= fe.n_mcps) ||
       (xcentre < 0) ||
       (xcentre >= fe.n_xpixels) ||
       (ycentre < 0) ||
       (ycentre >= fe.n_ypixels)
        return nothing
    end

    xval = get_mcp_xposition(xpos)
    nx = cdt.n_xbins
    ny = cdt.n_ybins

    # Define the coordinate limits for looping over pixels.
    ix_low = max(xcentre - nx, 0)
    ix_upp = min(xcentre + nx, m_nx - 1)
    iy_low = max(ycentre - ny, 0)
    iy_upp = min(ycentre + ny, m_ny - 1)

    # Initialize the pixel collection.
    pixels = Vector{PixelHit}(undef, 0)

    # Loop over pixels in the defined ranges.
    for ix = ix_low:ix_upp
        for iy = iy_low:iy_upp
            # Compute the lower boundaries for this pixel.
            xlow = ix * DETECTOR.x_size
            ylow = iy * DETECTOR.y_size

            # Shift the boundaries by subtracting half the active area.
            xlow -= 0.5 * DETECTOR.active
            ylow -= 0.5 * DETECTOR.active

            # Compute the upper boundaries.
            xupp = xlow + DETECTOR.x_size
            yupp = ylow + DETECTOR.y_size

            # Check if the charge is over threshold.
            if charge_over_threshold(cdt, xval, xlow, xupp, ypos, ylow, yupp)
                # Each hit is smeared independently.
                tprime = fe.smear_time ? smear_time(cdt, time) : time

                # Get the time pixel (without bounds checking).
                tpixel = get_tpixel_noheck(tprime)

                # Compute the x-pixel in the module pixel coordinate system.
                # (Assuming 'MCP' and Geometry.Detector.n_xpixels are defined.)
                xpixel = mcp * DETECTOR.n_xpixels + ix

                # Get the charge value.
                charge = get_charge(xval, xlow, xupp, ypos, ylow, yupp)

                # Append the new hit to the pixels collection.
                push!(pixels, PixelHit(xpixel, iy, tpixel, charge))
            end
        end
    end
    return pixels
end

function add_without_sharing!(
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    xpos::Float64,
    ypos::Float64,
    time::Float64,
)
    mcp = get_mcp_number(xpos)
    xcentre = get_mcp_xpixel(xpos)
    ycentre = get_ypixel(ypos)
    if (mcp < 0) ||
       (mcp >= fe.n_mcps) ||
       (xcentre < 0) ||
       (xcentre >= fe.n_xpixels) ||
       (ycentre < 0) ||
       (ycentre >= fe.n_ypixels)
        return
    end
    # smear_time
    tprime = fe.smear_time ? smear_time(cdt, time) : time
    tpixel = get_tpixel_nocheck(tprime)
    _fill_hit!(fe.pixelmap[mcp], xcentre, ycentre, tpixel)
    return
end

function add_without_sharing!(fe::FrontEnd, cdt::ChargeDepositTester, hit::HitCoordinate)
    add_without_sharing!(fe, cdt, hit.x, hit.y, hit.t)
end

function add_photon!(
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    xpos::Float64,
    ypos::Float64,
    time::Float64,
)
    if fe.sharing_mode
        mcp = get_mcp_number(xpos)
        pixels = find_pixels(fe, cdt, xpos, ypos)
        if pixels !== nothing
            for pixel in pixels
                _fill_hit!(fe.pixelmap[mcp], pixel)
            end
        end
    else
        add_without_sharing!(fe, cdt, xpos, ypos, time)
    end
end

function add_photon!(fe::FrontEnd, cdt::ChargeDepositTester, hit::HitCoordinate)
    add_photon!(fe, cdt, hit.x, hit.y, hit.t)
end

function reset(fe::FrontEnd)
    for i = 1:fe.n_mcps
        _reset!(fe.pixelmap[i].pixelmap)
    end
end

function get_hits(fe::FrontEnd)::Union{Nothing,Vector{PixelHit}}
    hitcollection = Vector{PixelHit}()
    for i = 1:fe.n_mcps
        pixelhits = _get_hits(fe.pixelmap[i], i - 1)
        if pixelhits !== nothing
            append!(hitcollection, pixelhits)
        end
    end
    return isempty(hitcollection) ? nothing : hitcollection
end
