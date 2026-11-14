"""
    EventProbability{T, M<:AbstractMatrix{T}, A<:AbstractArray{T,3}, V<:AbstractVector{T}}

Contiguous column-major storage for per-event unbinned PID probabilities:
- `yields::M`: Expected photon yields of size `(ntracks, npids)`
- `sig_prob::A`: Normalized signal PDF values of size `(nhits, ntracks, npids)` (`hit` is stride-1)
- `bkg_prob::V`: Background PDF values of size `(nhits,)`
"""
struct EventProbability{
    T<:AbstractFloat,
    M<:AbstractMatrix{T},
    A<:AbstractArray{T,3},
    V<:AbstractVector{T},
}
    ntracks::Int
    nhits::Int
    npids::Int
    yields::M
    sig_prob::A
    bkg_prob::V
end

function EventProbability(
    ntracks::Int,
    nhits::Int,
    npids::Int,
    ::Type{T} = Float64,
) where {T<:AbstractFloat}
    yields = zeros(T, ntracks, npids)
    sig_prob = zeros(T, nhits, ntracks, npids)
    bkg_prob = zeros(T, nhits)
    return EventProbability(ntracks, nhits, npids, yields, sig_prob, bkg_prob)
end

"""
    log_sum(prob::EventProbability{T}, ids::AbstractVector{<:Integer})::T

Zero-allocation evaluation of the unbinned extended log-likelihood for a given vector of
1-based hypothesis indices `ids` (one per track).
"""
@inline function log_sum(
    prob::EventProbability{T},
    ids::AbstractVector{<:Integer},
)::T where {T<:AbstractFloat}
    ntracks = prob.ntracks
    nhits = prob.nhits
    if nhits == 0
        return zero(T)
    end

    inv_nhits = one(T) / T(nhits)
    nsum = zero(T)
    @inbounds for t = 1:ntracks
        nsum += prob.yields[t, ids[t]]
    end
    bg_frac = max(T(nhits) - nsum, one(T)) * inv_nhits

    total_ll = zero(T)
    @inbounds for i = 1:nhits
        p = prob.bkg_prob[i] * bg_frac
        for t = 1:ntracks
            h = ids[t]
            fr = prob.yields[t, h] * inv_nhits
            p = muladd(fr, prob.sig_prob[i, t, h], p)
        end
        total_ll += log(max(p, floatmin(T)))
    end

    return total_ll
end

"""
    pattern_efficiency(
        particle, beta, mapper, spectrum, distribution, context, detector;
        n_photons = 1000, random_points = nothing,
    )::Float64

Estimates geometric and acceptance efficiency for a `(particle, beta)` hypothesis by tracing
`n_photons` Cherenkov photons. Accepts an optional pre-generated `random_points` vector for
deterministic, thread-safe, RNG-free execution.
"""
function pattern_efficiency(
    particle::Particle,
    beta::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    distribution::PhotonSpectrumDistribution,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR;
    n_photons::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
)::Float64
    if !spectrum_above_threshold(distribution) || n_photons <= 0
        return 0.0
    end

    radiator = context.radiator
    accepted = 0

    @inbounds for i = 1:n_photons
        pt = random_points === nothing ? random_point(radiator.depth) : random_points[i]
        sampled = spectrum_random_sampling(spectrum, distribution, pt.e)
        if sampled === nothing
            continue
        end
        energy, nphase, ngroup = sampled
        photon = create_approximate_photon(
            particle,
            beta,
            nphase,
            ngroup,
            energy,
            pt.emissionZ,
            pt.cosPhi,
            pt.sinPhi,
            context,
        )
        toffset = _get_time_offset(particle, beta, photon.zpos, particle.t0, context)
        hit = trace_photon(mapper, photon, toffset)
        if hit !== nothing && is_in_pixel(hit, detector, radiator)
            accepted += 1
        end
    end

    return accepted / n_photons
end

const PIXEL_QUAD_WEIGHTS = (0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 1.0)

@inline function _accumulate_pixel_signal_probs!(
    sig_prob::AbstractArray{T,3},
    pixel_idx::Int,
    track_idx::Int,
    pixel::PixelHit,
    particle::Particle,
    betas::AbstractVector{Float64},
    distributions::AbstractVector{PhotonSpectrumDistribution},
    inv_norms::AbstractVector{Float64},
    cfg::ReverseMappingConfig,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters,
    context::PhotonContext,
    detector::Detector,
    wedge::Wedge,
) where {T<:AbstractFloat}
    center = get_hit_coordinate(pixel, detector, context.radiator)
    dx = 0.5 * detector.x_size
    dy = 0.5 * detector.y_size
    xmin = center.x - dx
    xmax = center.x + dx
    ymin = center.y - dy
    ymax = center.y + dy

    xs = (xmin, xmin, xmax, xmax, center.x, center.x, xmin, xmax, center.x)
    ys = (ymin, ymax, ymin, ymax, ymin, ymax, center.y, center.y, center.y)

    area_factor = 0.25 * detector.x_size * detector.y_size * detector.t_bin

    @inbounds for q = 1:9
        sub_hit = HitCoordinate(xs[q], ys[q], center.t)
        slope = focus_slope(cfg, sub_hit.y)
        ypath = focus_ypath(cfg, sub_hit.y, context.radiator, wedge)
        w = PIXEL_QUAD_WEIGHTS[q] * area_factor

        for h in eachindex(betas)
            inv_n = inv_norms[h]
            if inv_n > 0.0
                p = hit_signal_probability(
                    particle,
                    betas[h],
                    distributions[h],
                    sub_hit,
                    cfg,
                    mapper,
                    spectrum,
                    signal,
                    context,
                    wedge;
                    slope = slope,
                    ypath = ypath,
                )
                sig_prob[pixel_idx, track_idx, h] += T(w * p * inv_n)
            end
        end
    end
    return nothing
end

function compute_event_probability!(
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    hits::AbstractVector{HitCoordinate},
    hypotheses::AbstractVector{<:Integer},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    background::AbstractBackgroundModel = FlatBackground(DETECTOR);
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_gen::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    threaded::Bool = true,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)

    fill!(event_prob.yields, zero(T))
    fill!(event_prob.sig_prob, zero(T))
    fill_bkg_probabilities!(event_prob.bkg_prob, background, hits)

    if ntracks == 0 || nhits == 0 || npids == 0
        return event_prob
    end

    pts =
        random_points === nothing ?
        [random_point(context.radiator.depth) for _ = 1:n_gen] : random_points

    # Precompute focus slopes and ypaths once per hit
    slopes = Vector{Float64}(undef, nhits)
    ypaths = Vector{Float64}(undef, nhits)
    @inbounds for i = 1:nhits
        y = hits[i].y
        slopes[i] = focus_slope(cfg, y)
        ypaths[i] = focus_ypath(cfg, y, context.radiator, wedge)
    end

    # Precompute distributions, efficiencies, and yields per (track, hypothesis)
    betas = Matrix{Float64}(undef, npids, ntracks)
    distributions = Matrix{PhotonSpectrumDistribution}(undef, npids, ntracks)
    inv_norms = zeros(Float64, npids, ntracks)

    @inbounds for t = 1:ntracks
        p = particles[t]
        plen = context.radiator.depth / p.zDir
        for h = 1:npids
            mass = get_particle_mass(Int(hypotheses[h]))
            beta = particle_beta(p, mass)
            dist = PhotonSpectrumDistribution(spectrum, beta)
            betas[h, t] = beta
            distributions[h, t] = dist

            eff = pattern_efficiency(
                p,
                beta,
                mapper,
                spectrum,
                dist,
                context,
                detector;
                n_photons = n_gen,
                random_points = pts,
            )
            if eff > 0.0
                inv_norms[h, t] = 1.0 / eff
                event_prob.yields[t, h] = T(spectrum_yield(dist, plen) * eff)
            end
        end
    end

    # Parallel evaluation over hits (stride-1 dimension) and tracks
    if threaded && nhits * ntracks > 16 && Threads.nthreads() > 1
        Threads.@threads for i = 1:nhits
            hit = @inbounds hits[i]
            sl = @inbounds slopes[i]
            yp = @inbounds ypaths[i]
            for t = 1:ntracks
                accumulate_hit_signal_probs!(
                    event_prob.sig_prob,
                    i,
                    t,
                    hit,
                    sl,
                    yp,
                    @inbounds(particles[t]),
                    @view(betas[:, t]),
                    @view(distributions[:, t]),
                    @view(inv_norms[:, t]),
                    cfg,
                    mapper,
                    spectrum,
                    signal,
                    context,
                )
            end
        end
    else
        @inbounds for t = 1:ntracks
            p = particles[t]
            b_view = @view betas[:, t]
            d_view = @view distributions[:, t]
            n_view = @view inv_norms[:, t]
            for i = 1:nhits
                accumulate_hit_signal_probs!(
                    event_prob.sig_prob,
                    i,
                    t,
                    hits[i],
                    slopes[i],
                    ypaths[i],
                    p,
                    b_view,
                    d_view,
                    n_view,
                    cfg,
                    mapper,
                    spectrum,
                    signal,
                    context,
                )
            end
        end
    end

    return event_prob
end

function compute_event_probability!(
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    pixels::AbstractVector{PixelHit},
    hypotheses::AbstractVector{<:Integer},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    background::AbstractBackgroundModel = FlatBackground(DETECTOR);
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_gen::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    threaded::Bool = true,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(pixels)
    npids = length(hypotheses)

    fill!(event_prob.yields, zero(T))
    fill!(event_prob.sig_prob, zero(T))
    fill_bkg_probabilities!(event_prob.bkg_prob, background, pixels)

    if ntracks == 0 || nhits == 0 || npids == 0
        return event_prob
    end

    pts =
        random_points === nothing ?
        [random_point(context.radiator.depth) for _ = 1:n_gen] : random_points

    betas = Matrix{Float64}(undef, npids, ntracks)
    distributions = Matrix{PhotonSpectrumDistribution}(undef, npids, ntracks)
    inv_norms = zeros(Float64, npids, ntracks)

    @inbounds for t = 1:ntracks
        p = particles[t]
        plen = context.radiator.depth / p.zDir
        for h = 1:npids
            mass = get_particle_mass(Int(hypotheses[h]))
            beta = particle_beta(p, mass)
            dist = PhotonSpectrumDistribution(spectrum, beta)
            betas[h, t] = beta
            distributions[h, t] = dist

            eff = pattern_efficiency(
                p,
                beta,
                mapper,
                spectrum,
                dist,
                context,
                detector;
                n_photons = n_gen,
                random_points = pts,
            )
            if eff > 0.0
                inv_norms[h, t] = 1.0 / eff
                event_prob.yields[t, h] = T(spectrum_yield(dist, plen) * eff)
            end
        end
    end

    if threaded && nhits * ntracks > 4 && Threads.nthreads() > 1
        Threads.@threads for i = 1:nhits
            pix = @inbounds pixels[i]
            for t = 1:ntracks
                _accumulate_pixel_signal_probs!(
                    event_prob.sig_prob,
                    i,
                    t,
                    pix,
                    @inbounds(particles[t]),
                    @view(betas[:, t]),
                    @view(distributions[:, t]),
                    @view(inv_norms[:, t]),
                    cfg,
                    mapper,
                    spectrum,
                    signal,
                    context,
                    detector,
                    wedge,
                )
            end
        end
    else
        @inbounds for t = 1:ntracks
            p = particles[t]
            b_view = @view betas[:, t]
            d_view = @view distributions[:, t]
            n_view = @view inv_norms[:, t]
            for i = 1:nhits
                _accumulate_pixel_signal_probs!(
                    event_prob.sig_prob,
                    i,
                    t,
                    pixels[i],
                    p,
                    b_view,
                    d_view,
                    n_view,
                    cfg,
                    mapper,
                    spectrum,
                    signal,
                    context,
                    detector,
                    wedge,
                )
            end
        end
    end

    return event_prob
end

"""
    solve_unbinned_likelihoods!(
        log_probs::AbstractMatrix{T},
        event_prob::EventProbability{T},
        particles::AbstractVector{Particle},
        hypotheses::AbstractVector{<:Integer};
        n_iter::Int = 3,
        cheated::Bool = false,
    )

Runs coordinate-ascent (or cheated MC truth) multi-track log-likelihood evaluation on
`event_prob`, writing the per-track, per-hypothesis log-likelihoods into `log_probs[t, h]`.
"""
function solve_unbinned_likelihoods!(
    log_probs::AbstractMatrix{T},
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    hypotheses::AbstractVector{<:Integer};
    n_iter::Int = 3,
    cheated::Bool = false,
) where {T<:AbstractFloat}
    ntracks = event_prob.ntracks
    npids = event_prob.npids
    fill!(log_probs, zero(T))
    if ntracks == 0 || npids == 0
        return log_probs
    end

    best_id = ones(Int, ntracks)

    if cheated
        @inbounds for t = 1:ntracks
            true_pid = abs(particles[t].pid)
            idx = findfirst(==(true_pid), hypotheses)
            best_id[t] = idx === nothing ? 1 : idx
        end
        @inbounds for t = 1:ntracks
            last_id = best_id[t]
            for h = 1:npids
                best_id[t] = h
                log_probs[t, h] = log_sum(event_prob, best_id)
            end
            best_id[t] = last_id
        end
    else
        @inbounds for iter = 1:n_iter
            if iter > 1
                for t = 1:ntracks
                    best_h = 1
                    best_val = log_probs[t, 1]
                    for h = 2:npids
                        if log_probs[t, h] > best_val
                            best_val = log_probs[t, h]
                            best_h = h
                        end
                    end
                    best_id[t] = best_h
                end
            end
            for t = 1:ntracks
                last_id = best_id[t]
                for h = 1:npids
                    best_id[t] = h
                    log_probs[t, h] = log_sum(event_prob, best_id)
                end
                best_id[t] = last_id
            end
        end
    end

    return log_probs
end

"""
    run_unbinned_pid!(
        log_probs::AbstractMatrix{T},
        particles::AbstractVector{Particle},
        hits::Union{AbstractVector{HitCoordinate},AbstractVector{PixelHit}};
        hypotheses = [211, 321, 2212],
        spectrum::PhotonSpectrum = PhotonSpectrum(),
        mapper::PhotonMapper = PhotonMapper(),
        background::AbstractBackgroundModel = FlatBackground(DETECTOR),
        cfg::ReverseMappingConfig = ReverseMappingConfig(),
        signal::SignalParameters = SIGNAL,
        context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
        detector::Detector = DETECTOR,
        wedge::Wedge = WEDGE,
        n_gen::Int = 1000,
        n_iter::Int = 3,
        cheated::Bool = false,
        random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
        threaded::Bool = true,
    )

Full unbinned PID reconstruction pipeline for a TORCH module. Populates `log_probs[t, h]`
for each track `t` and hypothesis `hypotheses[h]`.
"""
function run_unbinned_pid!(
    log_probs::AbstractMatrix{T},
    particles::AbstractVector{Particle},
    hits::Union{AbstractVector{HitCoordinate},AbstractVector{PixelHit}};
    hypotheses::AbstractVector{<:Integer} = [211, 321, 2212],
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    mapper::PhotonMapper = PhotonMapper(),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_gen::Int = 1000,
    n_iter::Int = 3,
    cheated::Bool = false,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    threaded::Bool = true,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)

    event_prob = EventProbability(ntracks, nhits, npids, T)
    compute_event_probability!(
        event_prob,
        particles,
        hits,
        hypotheses,
        spectrum,
        mapper,
        background;
        cfg = cfg,
        signal = signal,
        context = context,
        detector = detector,
        wedge = wedge,
        n_gen = n_gen,
        random_points = random_points,
        threaded = threaded,
    )
    solve_unbinned_likelihoods!(
        log_probs,
        event_prob,
        particles,
        hypotheses;
        n_iter = n_iter,
        cheated = cheated,
    )
    return log_probs
end

function run_unbinned_pid(
    particles::AbstractVector{Particle},
    hits::Union{AbstractVector{HitCoordinate},AbstractVector{PixelHit}};
    hypotheses::AbstractVector{<:Integer} = [211, 321, 2212],
    kwargs...,
)::Matrix{Float64}
    log_probs = zeros(Float64, length(particles), length(hypotheses))
    run_unbinned_pid!(log_probs, particles, hits; hypotheses = hypotheses, kwargs...)
    return log_probs
end

"""
    run_unbinned_local_pid!(
        log_probs::AbstractMatrix{T},
        particles::AbstractVector{Particle},
        hits::AbstractVector{HitCoordinate};
        hypotheses = [211, 321, 2212],
        ...
    )

Single-track (local) unbinned PID algorithm where each track's log-likelihood is evaluated
independently against the background model.
"""
function run_unbinned_local_pid!(
    log_probs::AbstractMatrix{T},
    particles::AbstractVector{Particle},
    hits::AbstractVector{HitCoordinate};
    hypotheses::AbstractVector{<:Integer} = [211, 321, 2212],
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    mapper::PhotonMapper = PhotonMapper(),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_gen::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    threaded::Bool = true,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)

    event_prob = EventProbability(ntracks, nhits, npids, T)
    compute_event_probability!(
        event_prob,
        particles,
        hits,
        hypotheses,
        spectrum,
        mapper,
        background;
        cfg = cfg,
        signal = signal,
        context = context,
        detector = detector,
        wedge = wedge,
        n_gen = n_gen,
        random_points = random_points,
        threaded = threaded,
    )

    if nhits == 0
        fill!(log_probs, zero(T))
        return log_probs
    end

    inv_nhits = one(T) / T(nhits)
    @inbounds for t = 1:ntracks
        for h = 1:npids
            yield = event_prob.yields[t, h]
            sig_frac = yield * inv_nhits
            bg_frac = max(T(nhits) - yield, one(T)) * inv_nhits
            ll = zero(T)
            for i = 1:nhits
                p = muladd(sig_frac, event_prob.sig_prob[i, t, h], bg_frac * event_prob.bkg_prob[i])
                ll += log(max(p, floatmin(T)))
            end
            log_probs[t, h] = ll
        end
    end

    return log_probs
end
