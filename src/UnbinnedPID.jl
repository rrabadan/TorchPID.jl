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
            particle, beta, nphase, ngroup, energy, pt.emissionZ, pt.cosPhi, pt.sinPhi, context,
        )
        toffset = _get_time_offset(particle, beta, photon.zpos, particle.t0, context)
        hit = trace_photon(mapper, photon, toffset)
        if hit !== nothing && is_in_pixel(hit, detector, radiator)
            accepted += 1
        end
    end

    return accepted / n_photons
end

function _prepare_event_probability!(
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    hits::AbstractVector,
    hypotheses::AbstractVector{<:Integer},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    background::AbstractBackgroundModel,
    context::PhotonContext,
    detector::Detector,
    n_gen::Int,
    random_points::Union{AbstractVector{RandomPoint},Nothing},
) where {T<:AbstractFloat}
    ntracks = length(particles)
    npids = length(hypotheses)

    fill!(event_prob.yields, zero(T))
    fill!(event_prob.sig_prob, zero(T))
    fill_bkg_probabilities!(event_prob.bkg_prob, background, hits)

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
                p, beta, mapper, spectrum, dist, context, detector;
                n_photons = n_gen, random_points = pts,
            )
            if eff > 0.0
                inv_norms[h, t] = 1.0 / eff
                event_prob.yields[t, h] = T(spectrum_yield(dist, plen) * eff)
            end
        end
    end

    return (betas, distributions, inv_norms)
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
    if ntracks == 0 || nhits == 0 || npids == 0
        fill!(event_prob.yields, zero(T))
        fill!(event_prob.sig_prob, zero(T))
        fill_bkg_probabilities!(event_prob.bkg_prob, background, hits)
        return event_prob
    end

    # Precompute focus slopes and ypaths once per hit
    slopes = Vector{Float64}(undef, nhits)
    ypaths = Vector{Float64}(undef, nhits)
    @inbounds for i = 1:nhits
        y = hits[i].y
        slopes[i] = focus_slope(cfg, y)
        ypaths[i] = focus_ypath(cfg, y, context.radiator, wedge)
    end

    betas, distributions, inv_norms = _prepare_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper,
        background, context, detector, n_gen, random_points,
    )

    # Parallel evaluation over hits (stride-1 dimension) and tracks
    if threaded && nhits * ntracks > 16 && Threads.nthreads() > 1
        Threads.@threads for i = 1:nhits
            hit = @inbounds hits[i]
            sl = @inbounds slopes[i]
            yp = @inbounds ypaths[i]
            for t = 1:ntracks
                accumulate_hit_signal_probs!(
                    event_prob.sig_prob, i, t, hit, sl, yp,
                    @inbounds(particles[t]),
                    @view(betas[:, t]), @view(distributions[:, t]), @view(inv_norms[:, t]),
                    cfg, mapper, spectrum, signal, context,
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
                    event_prob.sig_prob, i, t, hits[i], slopes[i], ypaths[i],
                    p, b_view, d_view, n_view,
                    cfg, mapper, spectrum, signal, context,
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
    if ntracks == 0 || nhits == 0 || npids == 0
        fill!(event_prob.yields, zero(T))
        fill!(event_prob.sig_prob, zero(T))
        fill_bkg_probabilities!(event_prob.bkg_prob, background, pixels)
        return event_prob
    end

    betas, distributions, inv_norms = _prepare_event_probability!(
        event_prob, particles, pixels, hypotheses, spectrum, mapper,
        background, context, detector, n_gen, random_points,
    )

    if threaded && nhits * ntracks > 4 && Threads.nthreads() > 1
        Threads.@threads for i = 1:nhits
            pix = @inbounds pixels[i]
            for t = 1:ntracks
                _accumulate_pixel_signal_probs!(
                    event_prob.sig_prob, i, t, pix, @inbounds(particles[t]),
                    @view(betas[:, t]), @view(distributions[:, t]), @view(inv_norms[:, t]),
                    cfg, mapper, spectrum, signal, context, detector, wedge,
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
                    event_prob.sig_prob, i, t, pixels[i], p,
                    b_view, d_view, n_view,
                    cfg, mapper, spectrum, signal, context, detector, wedge,
                )
            end
        end
    end

    return event_prob
end

"""
    run_unbinned_pid!(
        log_probs::AbstractMatrix{T},
        particles::AbstractVector{Particle},
        hits::Union{AbstractVector{HitCoordinate},AbstractVector{PixelHit}};
        hypotheses = [211, 321, 2212],
        spectrum = PhotonSpectrum(),
        mapper = PhotonMapper(),
        background = FlatBackground(DETECTOR),
        n_iter = 3,
        cheated = false,
        kwargs...,
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
    n_iter::Int = 3,
    cheated::Bool = false,
    kwargs...,
) where {T<:AbstractFloat}
    event_prob = EventProbability(length(particles), length(hits), length(hypotheses), T)
    compute_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper, background; kwargs...,
    )
    solve_unbinned_likelihoods!(
        log_probs, event_prob, particles, hypotheses; n_iter = n_iter, cheated = cheated,
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
        kwargs...,
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
    kwargs...,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)

    event_prob = EventProbability(ntracks, nhits, npids, T)
    compute_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper, background; kwargs...,
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
