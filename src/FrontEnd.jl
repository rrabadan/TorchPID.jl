"""
    MCPPixelMap

Type representing a Microchannel Plate (MCP) pixel map.

# Fields
- `data::Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}`: Dictionary mapping pixel coordinates to a vector of hit tuples.
- `n_xpixels::Int`: Number of pixels in the x-direction.
- `n_ypixels::Int`: Number of pixels in the y-direction.

# Constructor
    `MCPPixelMap(n_xpixels::Int, n_ypixels::Int)`:

Initializes an empty pixel map with the specified dimensions.

## Arguments
- `n_xpixels::Int`: Number of pixels in the x-direction.
- `n_ypixels::Int`: Number of pixels in the y-direction.

## Returns
- `MCPPixelMap`: An initialized MCPPixelMap object with empty data.
"""
mutable struct MCPPixelMap
    data::Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}
    n_xpixels::Int
    n_ypixels::Int
end

function MCPPixelMap(n_xpixels::Int, n_ypixels::Int)::MCPPixelMap
    MCPPixelMap(Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}(), n_xpixels, n_ypixels)
end

"""
    MCPImage

Type representing a Microchannel Plate (MCP) image with pixel data.

# Fields
- `n_tpixels::Int`: Number of pixels in the time dimension.
- `n_deadtime::Int`: Dead time in time pixels.
- `pixelmap::MCPPixelMap`: 3D structure containing hit data for each pixel.

# Constructor
    `MCPImage(n_xpixels::Int, n_ypixels::Int, n_tpixels::Int, n_deadtime::Int)`: 
  
Initializes an MCP image with the specified dimensions and an empty pixel map.
The pixel map is a 3D array where each element is a vector of tuples,
representing the time pixel and charge associated with a hit.

## Arguments
- `n_xpixels::Int`: Number of pixels in the x-direction.
- `n_ypixels::Int`: Number of pixels in the y-direction.
- `n_tpixels::Int`: Number of pixels in the time dimension.
- `n_deadtime::Int`: Dead time in time pixels.

## Returns
- `MCPImage`: An initialized MCPImage object with an empty pixel map.
"""
struct MCPImage
    n_tpixels::Int
    n_deadtime::Int
    pixelmap::MCPPixelMap
end

function MCPImage(n_xpixels::Int, n_ypixels::Int, n_tpixels::Int, n_deadtime::Int)::MCPImage
    pixelmap = MCPPixelMap(n_xpixels, n_ypixels)
    MCPImage(n_tpixels, n_deadtime, pixelmap)
end

"""
    _is_pixelmap_occupied(pixelmap::MCPPixelMap, xpixel::Int, ypixel::Int)

`_is_pixelmap_occupied` determines whether a pixel at the specified coordinates contains any hits.  
It returns `false` if the coordinates are out of bounds.  
A pixel is deemed occupied if its hit list is non-empty.

# Arguments
- `pixelmap::MCPPixelMap`: The pixel map to check.
- `xpixel::Int`: X-coordinate of the pixel.
- `ypixel::Int`: Y-coordinate of the pixel.

# Returns
- `Bool`: `true` if the pixel has at least one hit recorded; otherwise `false`.
"""
function _is_occupied(pixelmap::MCPPixelMap, xpixel::Int, ypixel::Int)
    data = pixelmap.data
    n_xpixels = pixelmap.n_xpixels
    n_ypixels = pixelmap.n_ypixels
    if (xpixel < 0) || (xpixel >= n_xpixels) || (ypixel < 0) || (ypixel >= n_ypixels)
        return false
    end

    # Check if the pixel is occupied, accounting for missing keys
    return haskey(data, (xpixel, ypixel)) && !isempty(data[(xpixel, ypixel)])
end

"""
    _reset_pixelmap!(pixelmap::MCPPixelMap)

`reset_pixelmap` clears all hit data in the pixel map without altering its structure.  
It iterates through all pixels and empties their hit lists, leaving the pixel map intact.

# Arguments
- `pixelmap::MCPPixelMap`: The pixel map to reset.

# Returns
- `Nothing`: The function modifies the input pixel map in place.
"""
function _reset_pixelmap!(pixelmap::MCPPixelMap)
    pixelmap.data = Dict{Tuple{Int,Int},Vector{Tuple{Int,Float64}}}()
end

"""
    _get_hits(mcpimage::MCPImage, mcp::Int, n_deadtime::Int, xpixel::Int, ypixel::Int)
    _get_hits(mcpimage::MCPImage, mcp::Int, n_deadtime::Int, xpixel::Int)
    _get_hits(mcpimage::MCPImage, mcp::Int, n_deadtime::Int)

`_get_hits` retrieves hits from a specific pixel or set of pixels in the MCP image.  
It extracts hits from the specified pixel location(s) and returns them as a collection of `PixelHit` objects.  
If no hits are found, it returns `nothing`.  
Hits are sorted in descending order of time, and only those with positive time values are included.  
The function respects the detector's dead time, excluding hits that occur within the dead time of a previous hit.  
Each returned `PixelHit` includes global coordinates and associated charge information.

# Arguments
- `mcpimage::MCPImage`: The MCP image containing pixel data.
- `mcp::Int`: MCP detector index.
- `n_deadtime::Int`: Dead time in time pixels.
- `xpixel::Int`: X-coordinate of the pixel (optional, depending on the method).
- `ypixel::Int`: Y-coordinate of the pixel (optional, depending on the method).

# Returns
- `Union{Nothing,Vector{PixelHit}}`: A collection of pixel hits or `nothing` if no hits are found.
"""
function _get_hits(
    pixelmap::MCPPixelMap,
    mcp::Int,
    n_deadtime::Int,
    xpixel::Int,
    ypixel::Int,
)::Union{Nothing,Vector{PixelHit}}

    data = pixelmap.data

    # Check if the key exists in the pixelmap
    if !haskey(data, (xpixel, ypixel)) || isempty(data[(xpixel, ypixel)])
        # If the key does not exist or the pixel is empty, return nothing
        return nothing
    end

    hitcollection = Vector{PixelHit}()

    offset = mcp * pixelmap.n_xpixels

    # Sort the hits in descending order by the first element of the tuple (assumes hit tuple is (tpixel, charge))
    sorted_hits = sort(data[(xpixel, ypixel)], by = first, rev = true)

    # If the first hit's t value is > 0, add it to the collection
    if sorted_hits[1][1] > 0
        push!(
            hitcollection,
            PixelHit(xpixel + offset, ypixel, sorted_hits[1][1], sorted_hits[1][2]),
        )
    end

    # Loop over the rest of the hits and add them if the gap to the previous hit is sufficient
    for t in Iterators.drop(eachindex(sorted_hits), 1)
        tPixel = sorted_hits[t][1]
        tPrev = sorted_hits[prevind(sorted_hits, t)][1]
        if (tPixel > 0) && (tPixel - tPrev > n_deadtime)
            push!(
                hitcollection,
                PixelHit(xpixel + offset, ypixel, tPixel, sorted_hits[t][2]),
            )
        end
    end
    return hitcollection
end

function _get_hits(
    pixelmap::MCPPixelMap,
    mcp::Int,
    n_deadtime::Int,
)::Union{Nothing,Vector{PixelHit}}

    hitcollection = Vector{PixelHit}(undef, 0)

    data = pixelmap.data
    for keys in keys(data)
        x = keys[1]
        y = keys[2]
        pixelhits = _get_hits(pixelmap, mcp, n_deadtime, x, y)
        if pixelhits !== nothing
            append!(hitcollection, pixelhits)
        end
    end

    return isempty(hitcollection) ? nothing : hitcollection
end

function _get_hits(
    pixelmap::MCPPixelMap,
    mcp::Int,
    n_deadtime::Int,
    xpixel::Int,
)::Union{Nothing,Vector{PixelHit}}

    hitcollection = Vector{PixelHit}(undef, 0)
    if xpixel < 0 || xpixel >= pixelmap.n_xpixels
        return nothing
    end
    for y = 0:pixelmap.n_ypixels-1
        if _is_pixelmap_occupied(mcpimage, xpixel, y)
            pixelhits = _get_hits(mcpimage, mcp, n_deadtime, xpixel, y)
            if pixelhits !== nothing
                append!(hitcollection, pixelhits)
            end
        end
    end
    return isempty(hitcollection) ? nothing : hitcollection
end

"""
    _get_number_hits(pixelmap::MCPPixelMap)

`_get_number_hits` computes the total number of hits in a pixel map by iterating through all pixels and summing the lengths of their hit lists. 
Each entry in a pixel's hit list represents a single hit.

# Arguments
- `pixelmap::MCPPixelMap`: The pixel map containing hit data.

# Returns
- `Int`: The total number of hits across all pixels.
"""
function _get_number_hits(pixelmap::MCPPixelMap)::Int
    total_hits = 0
    # Iterate over all entries in the pixelmap dictionary to count the total number of hits
    for (_, hits) in pixelmap.data
        total_hits += length(hits)
    end
    return total_hits
end

"""
    _fill_pixelmap!(mcpimage::MCPImage, xpixel::Int, ypixel::Int, tpixel::Int; charge::Float64=-1.0)
    _fill_pixelmap!(mcpimage::MCPImage, hit::PixelHit)

`_fill_pixelmap!` adds a hit to the MCPImage at the specified pixel position.  
It returns early if the pixel coordinates are outside the valid range of the MCP image.  
Handles negative time values up to `-mcpimage.n_deadtime`.  
Updates the `pixelmap` field of the `MCPImage` by appending a tuple `(tpixel, charge)` to the specified pixel.

# Arguments Method 1
- `mcpimage::MCPImage`: The MCP image to update.
- `xpixel::Int`: X-coordinate of the pixel.
- `ypixel::Int`: Y-coordinate of the pixel.
- `tpixel::Int`: Time coordinate of the pixel.
- `charge::Float64`: Charge value associated with the hit (default: -1.0).

# Arguments Method 2
- `mcpimage::MCPImage`: The MCP image to update.
- `hit::PixelHit`: A `PixelHit` object containing the pixel coordinates and charge. (for the second method)

# Returns
- `Nothing`: The function modifies the input `MCPImage` in place.
"""
function _fill_pixelmap!(
    mcpimage::MCPImage,
    xpixel::Int,
    ypixel::Int,
    tpixel::Int;
    charge::Float64 = -1.0,
)
    if (tpixel < -mcpimage.n_deadtime) ||
       (tpixel >= mcpimage.n_tpixels) ||
       (xpixel < 0) ||
       (xpixel >= mcpimage.pixelmap.n_xpixels) ||
       (ypixel < 0) ||
       (ypixel >= mcpimage.pixelmap.n_ypixels)
        return
    end
    key = (xpixel, ypixel)
    pixelmap = mcpimage.pixelmap.data
    if !haskey(pixelmap, key)
        pixelmap[key] = Vector{Tuple{Int,Float64}}()
    end
    # Append the hit to the pixelmap
    push!(pixelmap[key], (tpixel, charge))
end

function _fill_pixelmap!(mcpimage::MCPImage, hit::PixelHit)
    _fill_pixelmap!(mcpimage, hit.xpixel, hit.ypixel, hit.tpixel; charge = hit.charge)
end

"""
    FrontEnd(
        n_mcps::Int, n_xpixels::Int, n_ypixels::Int, n_tpixels::Int, n_deadtime::Int,
        sharing_mode::Bool, smear_time::Bool
    )

Type representing the front-end electronics configuration for the TORCH detector.

# Fields
- `n_mcps::Int`: Number of Microchannel Plate (MCP) detectors.
- `n_xpixels::Int`: Number of pixels in the x-direction per MCP.
- `n_ypixels::Int`: Number of pixels in the y-direction per MCP.
- `n_tpixels::Int`: Number of pixels in the time dimension.
- `n_deadtime::Int`: Dead time in time pixels.
- `sharing_mode::Bool`: Whether charge sharing between adjacent pixels is enabled.
- `smear_time::Bool`: Whether time smearing is applied to hits.

# Constructor
    `FrontEnd(; n_mcps::Int, n_xpixels::Int, n_ypixels::Int, n_tpixels::Int, n_deadtime::Int,
              sharing_mode::Bool, smear_time::Bool)`

Constructs a FrontEnd object with the specified parameters or defaults to values from the global DETECTOR configuration.
The `sharing_mode` parameter determines whether charge sharing between adjacent pixels is enabled, and the `smear_time` parameter enables time smearing to simulate timing imprecision.

# Keywords
- `n_mcps::Int`: Number of MCP detectors (default: DETECTOR.n_detectors).
- `n_xpixels::Int`: Number of pixels in the x-direction (default: DETECTOR.n_xpixels).
- `n_ypixels::Int`: Number of pixels in the y-direction (default: DETECTOR.n_ypixels).
- `n_tpixels::Int`: Number of pixels in the time dimension (default: DETECTOR.n_tpixels).
- `n_deadtime::Int`: Dead time in time pixels (default: DETECTOR.n_deadtime).
- `sharing_mode::Bool`: Whether to enable charge sharing between pixels (default: false).
- `smear_time::Bool`: Whether to apply time smearing to hits (default: false).

# Returns
- `FrontEnd`: A configured front-end electronics object.
"""
struct FrontEnd
    n_mcps::Int
    n_xpixels::Int
    n_ypixels::Int
    n_tpixels::Int
    n_deadtime::Int
    sharing_mode::Bool
    smear_time::Bool
end

function FrontEnd(;
    n_mcps::Int = DETECTOR.n_detectors,
    n_xpixels::Int = DETECTOR.n_xpixels,
    n_ypixels::Int = DETECTOR.n_ypixels,
    n_tpixels::Int = DETECTOR.n_tpixels,
    n_deadtime::Int = DETECTOR.n_deadtime,
    sharing_mode::Bool = false,
    smear_time::Bool = false,
)::FrontEnd
    FrontEnd(n_mcps, n_xpixels, n_ypixels, n_tpixels, n_deadtime, sharing_mode, smear_time)
end

"""
    MCPImageArray = Vector{MCPImage}
"""
const MCPImageArray = Vector{MCPImage}

"""
    create_mcp_images(fe::FrontEnd)::MCPImageArray

`create_mcp_images` initializes an array of MCP image objects based on the front-end configuration.  
For each MCP detector specified in the front-end configuration, it creates an MCP image initialized with the specified pixel dimensions and dead time.

# Arguments
- `fe::FrontEnd`: Front-end electronics configuration.

# Returns
- `MCPImageArray`: An array of initialized MCP image objects.
"""
function create_mcp_images(fe::FrontEnd)::MCPImageArray
    mcpimages = Vector{MCPImage}(undef, fe.n_mcps)
    for i = 1:fe.n_mcps
        mcpimages[i] = MCPImage(fe.n_xpixels, fe.n_ypixels, fe.n_tpixels, fe.n_deadtime)
    end
    return mcpimages
end

"""
    find_pixels(fe::FrontEnd, cdt::ChargeDepositTester, xpos::Float64, ypos::Float64)

`find_pixels` determines the pixels activated by a hit at a given position.  
It validates the MCP number and pixel coordinates to ensure they are within bounds.  
The function iterates over nearby pixels, checking if the charge exceeds the threshold.  
For each activated pixel, it calculates the time pixel and charge, creating a `PixelHit` object.  
This process identifies all pixels affected by the hit.

# Arguments
- `fe::FrontEnd`: Front-end electronics configuration.
- `cdt::ChargeDepositTester`: Charge deposit test configuration.
- `xpos::Float64`: X-position of the hit.
- `ypos::Float64`: Y-position of the hit.

# Returns
- `Union{Nothing,Vector{PixelHit}}`: Collection of activated pixels or nothing if out of bounds.
"""
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

    # Maximum possible number of pixels.
    max_pixels = (ix_upp - ix_low + 1) * (iy_upp - iy_low + 1)

    # Pre-allocate the maximum number of pixels.
    pixels = Vector{PixelHit}(undef, max_pixels)
    pixel_count = 0

    # Loop over pixels in the defined ranges.
    @inbounds for ix = ix_low:ix_upp
        @inbounds for iy = iy_low:iy_upp
            # Compute the lower boundaries for this pixel.
            xlow = ix * DETECTOR.x_size
            ylow = iy * DETECTOR.y_size

            @fastmath begin
                # Shift the boundaries by subtracting half the active area.
                xlow -= 0.5 * DETECTOR.active
                ylow -= 0.5 * DETECTOR.active

                # Compute the upper boundaries.
                xupp = xlow + DETECTOR.x_size
                yupp = ylow + DETECTOR.y_size
            end

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
                charge = get_charge(cdt, xval, xlow, xupp, ypos, ylow, yupp)

                # Directly insert into the pre-allocated array.
                pixel_count += 1
                pixels[pixel_count] = PixelHit(xpixel, iy, tpixel, charge)
            end
        end
    end
    return pixel_count == 0 ? nothing : view(pixels, 1:pixel_count)
end

"""
    add_without_sharing!(
        mcps::MCPImageArray, fe::FrontEnd, cdt::ChargeDepositTester,
        xpos::Float64, ypos::Float64, time::Float64
    )
    add_without_sharing!(
        mcps::MCPImageArray, fe::FrontEnd, cdt::ChargeDepositTester,
        hit::HitCoordinate
    )

`add_without_sharing!` adds a photon hit to the MCP image array without charge sharing.

# Arguments
- `mcps::MCPImageArray`: Array of MCP image objects to update.
- `fe::FrontEnd`: Front-end electronics configuration.
- `cdt::ChargeDepositTester`: Charge deposit test configuration.
- `xpos::Float64`: X-position of the hit.
- `ypos::Float64`: Y-position of the hit.
- `time::Float64`: Time of the hit.

# Returns
Nothing.
"""
function add_without_sharing!(
    mcps::MCPImageArray,
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
    _fill_pixelmap!(mcps[mcp+1], xcentre, ycentre, tpixel) # index should be mcp+1
    return
end

function add_without_sharing!(
    mcps::MCPImageArray,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    hit::HitCoordinate,
)
    add_without_sharing!(mcps, fe, cdt, hit.x, hit.y, hit.t)
end

"""
    add_photon!(
        mcps::MCPImageArray, fe::FrontEnd, cdt::ChargeDepositTester,
        xpos::Float64, ypos::Float64, time::Float64
    )
    add_photon!(
        mcps::MCPImageArray, fe::FrontEnd, cdt::ChargeDepositTester,
        hit::HitCoordinate
    )

`add_photon!` updates the MCP image array with a photon hit, considering the front-end configuration.  
If charge sharing is enabled, it identifies affected pixels and updates their pixel maps.  
Otherwise, it directly updates the pixel map of the corresponding MCP image.

# Arguments
- `mcps::MCPImageArray`: Array of MCP image objects to update.
- `fe::FrontEnd`: Front-end electronics configuration.
- `cdt::ChargeDepositTester`: Charge deposit test configuration.
- `xpos::Float64`: X-position of the hit.
- `ypos::Float64`: Y-position of the hit.
- `time::Float64`: Time of the hit.

# Returns
Nothing.
"""
function add_photon!(
    mcps::MCPImageArray,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    xpos::Float64,
    ypos::Float64,
    time::Float64,
)
    if fe.sharing_mode
        mcp = get_mcp_number(xpos)

        if (mcp < 0) || mcp >= fe.n_mcps
            return
        end

        # Only look up pixels if mcp is valid
        pixels = find_pixels(fe, cdt, xpos, ypos)
        if pixels !== nothing
            mcpimg = mcps[mcp+1] # should it be mcp+1 the index?
            for pixel in pixels
                _fill_pixelmap!(mcpimg, pixel)
            end
        end
    else
        add_without_sharing!(mcps, fe, cdt, xpos, ypos, time)
    end
end

function add_photon!(
    mcps::MCPImageArray,
    fe::FrontEnd,
    cdt::ChargeDepositTester,
    hit::HitCoordinate,
)
    add_photon!(mcps, fe, cdt, hit.x, hit.y, hit.t)
end

"""
    reset!(mcps::MCPImageArray, fe::FrontEnd)

`reset!` clears all hits from the MCP image array by resetting the pixel map of each MCP image.  
This function is useful for reinitializing the detector state before starting a new measurement or simulation.

# Arguments
- `mcps::MCPImageArray`: Array of MCP image objects to reset.
- `fe::FrontEnd`: Front-end electronics configuration.

# Returns
Nothing.
"""
function reset!(mcps::MCPImageArray, fe::FrontEnd)
    for i = 1:fe.n_mcps
        _reset_pixelmap!(mcps[i].pixelmap)
    end
end

"""
    get_hits(mcps::MCPImageArray, fe::FrontEnd)

`get_hits` retrieves all hits from the MCP image array by iterating through each MCP image and collecting hits from all pixels.  
If no hits are found, it returns `nothing`. This function aggregates hits from all MCPs in the detector.

# Arguments
- `mcps::MCPImageArray`: Array of MCP image objects to get hits from.
- `fe::FrontEnd`: Front-end electronics configuration.

# Returns
- `Union{Nothing,Vector{PixelHit}}`: Collection of all hits or nothing if no hits are found.
"""
function get_hits(mcps::MCPImageArray, fe::FrontEnd)::Union{Nothing,Vector{PixelHit}}
    hitcollection = Vector{PixelHit}()
    for i = 1:fe.n_mcps
        pixelmap = mcps[i].pixelmap
        n_deadtime = mcps[i].n_deadtime
        pixelhits = _get_hits(pixelmap, i - 1, n_deadtime)
        if pixelhits !== nothing
            append!(hitcollection, pixelhits)
        end
    end
    return isempty(hitcollection) ? nothing : hitcollection
end
