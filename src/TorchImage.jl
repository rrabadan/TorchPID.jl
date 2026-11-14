struct TorchImage{T<:AbstractFloat,A<:AbstractArray{T,3}}
    nx::Int
    ny::Int
    nt::Int
    tmin::T
    tmax::T
    data::A
end

mutable struct TorchImageAccumulator{T<:AbstractFloat,A<:AbstractArray{T,3}}
    image::TorchImage{T,A}
    integral::T
    filled::Int
end

function TorchImage(detector::Detector = DETECTOR)::TorchImage{Float64,Array{Float64,3}}
    nx = detector.n_xtotal
    ny = detector.n_ypixels
    nt = detector.n_tpixels
    tmin = detector.t_min
    tmax = detector.t_max

    # Column-major 3D layout (nt, ny, nx): time dimension is stride-1 contiguous
    data = zeros(Float64, nt, ny, nx)

    return TorchImage(nx, ny, nt, tmin, tmax, data)
end

function TorchImageAccumulator(
    detector::Detector = DETECTOR,
)::TorchImageAccumulator{Float64,Array{Float64,3}}
    image = TorchImage(detector)
    return TorchImageAccumulator(image, 0.0, 0)
end

function reset!(acc::TorchImageAccumulator)::Nothing
    acc.integral = 0.0
    acc.filled = 0
    fill!(acc.image.data, zero(eltype(acc.image.data)))
    return nothing
end

@inline function is_validpixel(image::TorchImage, x::Int, y::Int, t::Int)::Bool
    return (0 <= x < image.nx) && (0 <= y < image.ny) && (0 <= t < image.nt)
end

function fill!(acc::TorchImageAccumulator, x::Int, y::Int, t::Int)::Bool
    image = acc.image
    if is_validpixel(image, x, y, t)
        @inbounds image.data[t+1, y+1, x+1] += 1.0
        acc.integral += 1.0
        acc.filled += 1
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
)::Bool
    image = acc.image
    if !is_validpixel(image, x, y, t)
        return false
    end

    nwin = cdt.n_time
    ntot = length(cdt.time_cache)

    tmin = max(0, t - nwin)
    tmax = min(image.nt - 1, t + nwin)  # image.nt is exclusive upper bound

    # Calculate start and end indices for time weights (1-based indexing)
    tw_start = nwin + tmin - t + 1
    tw_end = nwin + tmax - t + 1

    # Bounds check
    if tw_end > ntot || tw_start < 1
        error("Time weights indices out of bounds")
    end

    # Stride-1 contiguous view along dimension 1 (time) of column-major (nt, ny, nx) array
    @inbounds begin
        time_weights = @view cdt.time_cache[tw_start:tw_end]
        image_slice = @view image.data[(tmin+1):(tmax+1), y+1, x+1]
        image_slice .+= time_weights
        acc.integral += sum(time_weights)
    end
    acc.filled += 1

    return true
end

function fill!(acc::TorchImageAccumulator, cdt::ChargeDepositTester, hit::PixelHit)::Bool
    return fill_smeared!(acc, cdt, hit.x, hit.y, hit.t)
end

function fill!(
    acc::TorchImageAccumulator,
    cdt::ChargeDepositTester,
    hit::HitCoordinate,
)::Bool
    xpixel = get_xpixel(hit.x)
    ypixel = get_ypixel(hit.y)
    tpixel = get_tpixel(hit.t)
    return fill_smeared!(acc, cdt, xpixel, ypixel, tpixel)
end
