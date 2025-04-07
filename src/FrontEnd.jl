struct MCPImage
    n_xpixels::Int
    n_ypixels::Int
    n_tpixels::Int
    n_deadtime::Int
end

const MCPPixelMap = Vector{Vector{Vector{Tuple{Int,Float64}}}}

function create_pixel_map(mcp_image::MCPImage)::MCPPixelMap
    return [
        [Tuple{Int,Float64}[] for _ = 1:mcp_image.n_ypixels] for _ = 1:mcp_image.n_xpixels
    ]
end

function fill_hit!(
    pixelmap::MCPPixelMap,
    mcpimage::MCPImage,
    xpixel::Int,
    ypixel::Int,
    tpixel::Int,
    charge::Float64,
)

    if (tpixel < -mcpimage.n_deadtime) ||
       (tpixel >= mcpimage::n_tpixels) ||
       (xpixel < 0) ||
       (xpixel >= mcpimage.n_xpixels) ||
       (ypixel < 0) ||
       (ypixel >= mcpimage.n_ypixels)
        return
    end

    # Add the hit to the pixel map
    push!(pixelmap[xpixel+1][ypixel+1], (tpixel, charge))
end

function fill_hit!(pixelmap::MCPPixelMap, mcpimage::MCPImage, hit::PixelHit)
    fill_hit!(pixelmap, mcpimage, hit.xpixel, hit.ypixel, hit.tpixel, hit.charge)
end


function is_occupied(pixelmap::MCPPixelMap, xpixel::Int, ypixel::Int)

    if (xpixel < 0) ||
       (xpixel >= length(pixelmap)) ||
       (ypixel < 0) ||
       (ypixel >= length(pixelmap[1]))
        return false
    end

    # Check if the pixel is occupied
    return !isempty(pixelmap[xpixel+1][ypixel+1]) # && any(x -> x[1] == tpixel, pixelMap[xpixel][ypixel])
end

function get_hits(
    pixelmap::MCPPixelMap,
    mcpimage::MCPImage,
    xpixel::Int,
    ypixel::Int,
    mcp::Int,
)::Union{Nothing,Vector{PixelHit}}

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
    sorted_hits = sort(pixelmap[x][y], by=first, rev=true)

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
function get_hits(
    pixelmap::MCPPixelMap,
    mcpimage::MCPImage,
    mcp::Int,
)::Union{Nothing,Vector{PixelHit}}

    hitcollection = Vector{PixelHit}(undef, 0)
    for x = 1:mcpimage.n_xpixels
        for y = 1:mcpimage.n_ypixels
            if is_occupied(pixelmap, x - 1, y - 1)
                pixelhits = get_hits(pixelmap, mcpimage, x - 1, y - 1, mcp)
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

function get_hits(
    pixelmap::MCPPixelMap,
    mcpimage::MCPImage,
    xpixel::Int,
    mcp::Int,
)::Union{Nothing,Vector{PixelHit}}

    hitcollection = Vector{PixelHit}(undef, 0)
    if xpixel < 0 || xpixel >= mcpimage.n_xpixels
        return nothing
    end
    for y = 1:mcpimage.n_ypixels
        if is_occupied(pixelmap, xpixel, y - 1)
            pixelhits = get_hits(pixelmap, mcpimage, xpixel, y - 1, mcp)
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
function get_number_hits(pixelmap::MCPPixelMap)::Int
    total_hits = 0
    # Iterate over all pixels in the pixelmap to count the total number of hits
    for x in pixelmap
        for y in x
            total_hits += length(y)
        end
    end
    return total_hits
end
