struct TorchImage
    nx::Int
    ny::Int
    nt::Int
    tmin::Float64
    tmax::Float64
    data::Vector{Float64}
end

mutable struct TorchImageAccumulator
    image::TorchImage
    integral::Float64
    filled::Int
end

function TorchImage(detector::Detector)::TorchImage
    nx = detector.n_xtotal
    ny = detector.n_ypixels
    nt = detector.n_tpixels
    tmin = detector.t_min
    tmax = detector.t_max

    n_bins = nx * ny * nt
    data = zeros(Float64, n_bins)

    return TorchImage(nx, ny, nt, tmin, tmax, vec(data))
end

function reset!(acc::TorchImageAccumulator)::Nothing
    acc.integral = 0.0
    acc.filled = 0
    fill!(acc.image.data, 0.0)
    return nothing
end

function is_validpixel(image::TorchImage, x::Int, y::Int, t::Int)::Bool
    return (0 <= x < image.nx) && (0 <= y < image.ny) && (0 <= t < image.nt)
end

function fill!(acc::TorchImageAccumulator, x::Int, y::Int, t::Int)::Bool
    image = acc.image
    if is_validpixel(image, x, y, t)
        pixel = (image.nt * image.ny) * x + image.nt * y + t
        image.data[pixel+1] += 1.0
        return true
    end
    return false
end

function fill!(acc::TorchImageAccumulator, hit::PixelHit)::Bool
    return fill!(acc, hit.x, hit.y, hit.t)
end

function fill_smeared!(
    acc::TorchImageAccumulator,
    cdt::ChargeDepositTester,
    x::Int,
    y::Int,
    t::Int,
)::Union{Float64,Nothing}
    image = acc.image
    if !is_validpixel(image, x, y, t)
        return false
    end

    nwin = cdt.n_time
    ntot = length(cdt.time_cache)

    tmin = max(0, t - nwin)
    tmax = min(image.nt - 1, t + nwin)  # image.nt is exclusive upper bound

    offset = (image.nt * image.ny) * x + image.nt * y

    # Calculate start and end indices for time weights (1-based indexing)
    tw_start = nwin + tmin - t + 1
    tw_end = nwin + tmax - t + 2  # +1 for inclusive range in Julia

    # Bounds check
    if tw_end - 1 > ntot || tw_start < 1
        error("Time weights indices out of bounds")
    end

    # Get views into the relevant slices
    time_weights = @view cdt.time_cache[tw_start:tw_end-1]
    image_slice = @view image.data[(offset+tmin+1):(offset+tmax+1)]

    # Add time weights to image slice (in-place)
    image_slice .+= time_weights

    acc.integral += sum(time_weights)
    acc.filled += 1

    return true
end

function fill_smeared!(
    acc::TorchImageAccumulator,
    cdt::ChargeDepositTester,
    hit::PixelHit,
)::Bool
    return fill_smeared!(acc, cdt, hit.x, hit.y, hit.t)
end

function fill_smeared!(
    acc::TorchImageAccumulator,
    cdt::ChargeDepositTester,
    hit::HitCoordinate,
)::Bool
    xpixel = get_xpixel(hit.x)
    ypixel = get_ypixel(hit.y)
    tpixel = get_tpixel(hit.t)
    return fill_smeared!(acc, cdt, xpixel, ypixel, tpixel)
end
