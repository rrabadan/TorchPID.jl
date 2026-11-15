"""
    UnbinnedPIDWorkspace{T<:AbstractFloat}

Reusable per-thread scratch workspace for zero-allocation execution of
`compute_event_probability!`, `run_unbinned_pid!`, `run_unbinned_local_pid!`,
and `run_unbinned_pid_batch!`.
"""
struct UnbinnedPIDWorkspace{T<:AbstractFloat}
    yields_buf::Vector{T}
    sig_prob_buf::Vector{T}
    bkg_prob_buf::Vector{T}
    betas_buf::Vector{Float64}
    inv_norms_buf::Vector{Float64}
    distributions_buf::Vector{PhotonSpectrumDistribution}
    slopes_buf::Vector{Float64}
    ypaths_buf::Vector{Float64}
    best_id_buf::Vector{Int}
    random_points_buf::Vector{RandomPoint}
end

function UnbinnedPIDWorkspace(
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    ::Type{T} = Float64;
    max_tracks::Int = 16,
    max_hits::Int = 750,
    max_pids::Int = 3,
    n_gen::Int = 1000,
) where {T<:AbstractFloat}
    nt_np = max_tracks * max_pids
    nh_nt_np = max_hits * max_tracks * max_pids
    yields_buf = zeros(T, nt_np)
    sig_prob_buf = zeros(T, nh_nt_np)
    bkg_prob_buf = zeros(T, max_hits)
    betas_buf = zeros(Float64, nt_np)
    inv_norms_buf = zeros(Float64, nt_np)
    distributions_buf = [PhotonSpectrumDistribution(spectrum, 1.0) for _ = 1:nt_np]
    slopes_buf = Vector{Float64}(undef, max_hits)
    ypaths_buf = Vector{Float64}(undef, max_hits)
    best_id_buf = ones(Int, max_tracks)
    random_points_buf = Vector{RandomPoint}(undef, n_gen)
    return UnbinnedPIDWorkspace{T}(
        yields_buf,
        sig_prob_buf,
        bkg_prob_buf,
        betas_buf,
        inv_norms_buf,
        distributions_buf,
        slopes_buf,
        ypaths_buf,
        best_id_buf,
        random_points_buf,
    )
end

@inline function _workspace_event_probability!(
    ws::UnbinnedPIDWorkspace{T},
    ntracks::Int,
    nhits::Int,
    npids::Int,
) where {T<:AbstractFloat}
    nt_np = ntracks * npids
    nh_nt_np = nhits * ntracks * npids
    length(ws.yields_buf) < nt_np && resize!(ws.yields_buf, nt_np)
    length(ws.sig_prob_buf) < nh_nt_np && resize!(ws.sig_prob_buf, nh_nt_np)
    length(ws.bkg_prob_buf) < nhits && resize!(ws.bkg_prob_buf, nhits)
    yields = reshape(@view(ws.yields_buf[1:nt_np]), ntracks, npids)
    sig_prob = reshape(@view(ws.sig_prob_buf[1:nh_nt_np]), nhits, ntracks, npids)
    bkg_prob = @view(ws.bkg_prob_buf[1:nhits])
    return EventProbability(ntracks, nhits, npids, yields, sig_prob, bkg_prob)
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

function _prepare_hypotheses_threaded!(
    yields::AbstractMatrix{T},
    betas::AbstractMatrix{Float64},
    distributions::AbstractMatrix{PhotonSpectrumDistribution},
    inv_norms::AbstractMatrix{Float64},
    particles::AbstractVector{Particle},
    hypotheses::AbstractVector{<:Integer},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    context::PhotonContext,
    detector::Detector,
    n_gen::Int,
    pts::AbstractVector{RandomPoint},
    in_place_dists::Bool,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    npids = length(hypotheses)
    nt_np = npids * ntracks
    Threads.@threads for idx = 1:nt_np
        t = div(idx - 1, npids) + 1
        h = rem(idx - 1, npids) + 1
        p = @inbounds particles[t]
        plen = context.radiator.depth / p.zDir
        mass = get_particle_mass(Int(@inbounds hypotheses[h]))
        beta = particle_beta(p, mass)
        dist = if in_place_dists
            d = update_distribution!(@inbounds(distributions[h, t]), spectrum, beta)
            @inbounds distributions[h, t] = d
            d
        else
            d = PhotonSpectrumDistribution(spectrum, beta)
            @inbounds distributions[h, t] = d
            d
        end
        @inbounds betas[h, t] = beta
        eff = pattern_efficiency(
            p, beta, mapper, spectrum, dist, context, detector;
            n_photons = n_gen, random_points = pts,
        )
        if eff > 0.0
            @inbounds inv_norms[h, t] = 1.0 / eff
            @inbounds yields[t, h] = T(spectrum_yield(dist, plen) * eff)
        end
    end
    return nothing
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
    workspace::Union{UnbinnedPIDWorkspace,Nothing} = nothing,
    threaded::Bool = false,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    npids = length(hypotheses)
    nt_np = npids * ntracks

    fill!(event_prob.yields, zero(T))
    fill!(event_prob.sig_prob, zero(T))
    fill_bkg_probabilities!(event_prob.bkg_prob, background, hits)

    if random_points !== nothing
        pts = random_points
    elseif workspace !== nothing
        length(workspace.random_points_buf) < n_gen && resize!(workspace.random_points_buf, n_gen)
        depth = context.radiator.depth
        @inbounds for i = 1:n_gen
            workspace.random_points_buf[i] = random_point(depth)
        end
        pts = @view(workspace.random_points_buf[1:n_gen])
    else
        pts = [random_point(context.radiator.depth) for _ = 1:n_gen]
    end

    if workspace === nothing
        betas = Matrix{Float64}(undef, npids, ntracks)
        distributions = Matrix{PhotonSpectrumDistribution}(undef, npids, ntracks)
        inv_norms = zeros(Float64, npids, ntracks)
        if threaded && nt_np >= 6 && Threads.nthreads() > 1
            _prepare_hypotheses_threaded!(
                event_prob.yields, betas, distributions, inv_norms,
                particles, hypotheses, spectrum, mapper, context, detector, n_gen, pts, false,
            )
        else
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
        end
        return (betas, distributions, inv_norms)
    else
        length(workspace.betas_buf) < nt_np && resize!(workspace.betas_buf, nt_np)
        length(workspace.inv_norms_buf) < nt_np && resize!(workspace.inv_norms_buf, nt_np)
        while length(workspace.distributions_buf) < nt_np
            push!(workspace.distributions_buf, PhotonSpectrumDistribution(spectrum, 1.0))
        end

        betas_v = reshape(@view(workspace.betas_buf[1:nt_np]), npids, ntracks)
        dists_v = reshape(@view(workspace.distributions_buf[1:nt_np]), npids, ntracks)
        inv_norms_v = reshape(@view(workspace.inv_norms_buf[1:nt_np]), npids, ntracks)
        fill!(inv_norms_v, 0.0)

        if threaded && nt_np >= 6 && Threads.nthreads() > 1
            _prepare_hypotheses_threaded!(
                event_prob.yields, betas_v, dists_v, inv_norms_v,
                particles, hypotheses, spectrum, mapper, context, detector, n_gen, pts, true,
            )
        else
            @inbounds for t = 1:ntracks
                p = particles[t]
                plen = context.radiator.depth / p.zDir
                for h = 1:npids
                    mass = get_particle_mass(Int(hypotheses[h]))
                    beta = particle_beta(p, mass)
                    dist = update_distribution!(dists_v[h, t], spectrum, beta)
                    dists_v[h, t] = dist
                    betas_v[h, t] = beta

                    eff = pattern_efficiency(
                        p, beta, mapper, spectrum, dist, context, detector;
                        n_photons = n_gen, random_points = pts,
                    )
                    if eff > 0.0
                        inv_norms_v[h, t] = 1.0 / eff
                        event_prob.yields[t, h] = T(spectrum_yield(dist, plen) * eff)
                    end
                end
            end
        end
        return (betas_v, dists_v, inv_norms_v)
    end
end

function _accumulate_hits_threaded!(
    sig_prob::AbstractArray{T,3},
    particles::AbstractVector{Particle},
    hits::AbstractVector{HitCoordinate},
    slopes::Vector{Float64},
    ypaths::Vector{Float64},
    betas::AbstractMatrix{Float64},
    distributions::AbstractMatrix{PhotonSpectrumDistribution},
    inv_norms::AbstractMatrix{Float64},
    cfg::ReverseMappingConfig,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters,
    context::PhotonContext,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    Threads.@threads for i = 1:nhits
        hit = @inbounds hits[i]
        sl = @inbounds slopes[i]
        yp = @inbounds ypaths[i]
        for t = 1:ntracks
            accumulate_hit_signal_probs!(
                sig_prob, i, t, hit, sl, yp,
                @inbounds(particles[t]),
                @view(betas[:, t]), @view(distributions[:, t]), @view(inv_norms[:, t]),
                cfg, mapper, spectrum, signal, context,
            )
        end
    end
    return nothing
end

function _accumulate_pixels_threaded!(
    sig_prob::AbstractArray{T,3},
    particles::AbstractVector{Particle},
    pixels::AbstractVector{PixelHit},
    betas::AbstractMatrix{Float64},
    distributions::AbstractMatrix{PhotonSpectrumDistribution},
    inv_norms::AbstractMatrix{Float64},
    cfg::ReverseMappingConfig,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters,
    context::PhotonContext,
    detector::Detector,
    wedge::Wedge,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(pixels)
    Threads.@threads for i = 1:nhits
        pix = @inbounds pixels[i]
        for t = 1:ntracks
            _accumulate_pixel_signal_probs!(
                sig_prob, i, t, pix, @inbounds(particles[t]),
                @view(betas[:, t]), @view(distributions[:, t]), @view(inv_norms[:, t]),
                cfg, mapper, spectrum, signal, context, detector, wedge,
            )
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
    workspace::Union{UnbinnedPIDWorkspace,Nothing} = nothing,
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
    slopes = if workspace === nothing
        Vector{Float64}(undef, nhits)
    else
        length(workspace.slopes_buf) < nhits && resize!(workspace.slopes_buf, nhits)
        workspace.slopes_buf
    end
    ypaths = if workspace === nothing
        Vector{Float64}(undef, nhits)
    else
        length(workspace.ypaths_buf) < nhits && resize!(workspace.ypaths_buf, nhits)
        workspace.ypaths_buf
    end
    @inbounds for i = 1:nhits
        y = hits[i].y
        slopes[i] = focus_slope(cfg, y)
        ypaths[i] = focus_ypath(cfg, y, context.radiator, wedge)
    end

    betas, distributions, inv_norms = _prepare_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper,
        background, context, detector, n_gen, random_points, workspace, threaded,
    )

    # Parallel evaluation over hits (stride-1 dimension) and tracks
    if threaded && nhits * ntracks > 16 && Threads.nthreads() > 1
        _accumulate_hits_threaded!(
            event_prob.sig_prob, particles, hits, slopes, ypaths,
            betas, distributions, inv_norms, cfg, mapper, spectrum, signal, context,
        )
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
    workspace::Union{UnbinnedPIDWorkspace,Nothing} = nothing,
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
        background, context, detector, n_gen, random_points, workspace, threaded,
    )

    if threaded && nhits * ntracks > 4 && Threads.nthreads() > 1
        _accumulate_pixels_threaded!(
            event_prob.sig_prob, particles, pixels,
            betas, distributions, inv_norms, cfg, mapper, spectrum, signal, context, detector, wedge,
        )
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
        workspace = nothing,
        kwargs...,
    )

Full unbinned PID reconstruction pipeline for a TORCH module. Populates `log_probs[t, h]`
for each track `t` and hypothesis `hypotheses[h]`. Passing a reusable `workspace::UnbinnedPIDWorkspace{T}`
makes per-event evaluation completely zero-allocation.
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
    workspace::Union{UnbinnedPIDWorkspace{T},Nothing} = nothing,
    kwargs...,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)
    if workspace === nothing
        event_prob = EventProbability(ntracks, nhits, npids, T)
        best_id_buf = nothing
    else
        event_prob = _workspace_event_probability!(workspace, ntracks, nhits, npids)
        best_id_buf = workspace.best_id_buf
    end
    compute_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper, background;
        workspace = workspace, kwargs...,
    )
    solve_unbinned_likelihoods!(
        log_probs, event_prob, particles, hypotheses;
        n_iter = n_iter, cheated = cheated, best_id_buf = best_id_buf,
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
        workspace = nothing,
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
    workspace::Union{UnbinnedPIDWorkspace{T},Nothing} = nothing,
    kwargs...,
) where {T<:AbstractFloat}
    ntracks = length(particles)
    nhits = length(hits)
    npids = length(hypotheses)

    if workspace === nothing
        event_prob = EventProbability(ntracks, nhits, npids, T)
    else
        event_prob = _workspace_event_probability!(workspace, ntracks, nhits, npids)
    end
    compute_event_probability!(
        event_prob, particles, hits, hypotheses, spectrum, mapper, background;
        workspace = workspace, kwargs...,
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

function _run_unbinned_pid_range!(
    log_probs_batch::AbstractVector{<:AbstractMatrix{T}},
    particles_batch::AbstractVector{<:AbstractVector{Particle}},
    hits_batch::AbstractVector{<:AbstractVector},
    r::UnitRange{Int},
    ws::UnbinnedPIDWorkspace{T},
    hypotheses::AbstractVector{<:Integer},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    background::AbstractBackgroundModel,
    n_iter::Int,
    cheated::Bool,
    cfg::ReverseMappingConfig,
    signal::SignalParameters,
    context::PhotonContext,
    detector::Detector,
    wedge::Wedge,
    n_gen::Int,
    random_points::Union{AbstractVector{RandomPoint},Nothing},
) where {T<:AbstractFloat}
    @inbounds for ev in r
        run_unbinned_pid!(
            log_probs_batch[ev],
            particles_batch[ev],
            hits_batch[ev];
            hypotheses = hypotheses,
            spectrum = spectrum,
            mapper = mapper,
            background = background,
            n_iter = n_iter,
            cheated = cheated,
            workspace = ws,
            cfg = cfg,
            signal = signal,
            context = context,
            detector = detector,
            wedge = wedge,
            n_gen = n_gen,
            random_points = random_points,
            threaded = false,
        )
    end
    return nothing
end

"""
    run_unbinned_pid_batch!(
        log_probs_batch,
        particles_batch,
        hits_batch;
        hypotheses = [211, 321, 2212],
        spectrum = PhotonSpectrum(),
        mapper = PhotonMapper(),
        background = FlatBackground(DETECTOR),
        n_iter = 3,
        cheated = false,
        nthreads = Threads.nthreads(),
        kwargs...,
    )

Multithreaded batch evaluation of `run_unbinned_pid!` across `length(particles_batch)` events,
using one zero-allocation `UnbinnedPIDWorkspace` per thread chunk.
"""
function run_unbinned_pid_batch!(
    log_probs_batch::AbstractVector{<:AbstractMatrix{T}},
    particles_batch::AbstractVector{<:AbstractVector{Particle}},
    hits_batch::AbstractVector{<:AbstractVector};
    hypotheses::AbstractVector{<:Integer} = [211, 321, 2212],
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    mapper::PhotonMapper = PhotonMapper(),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    n_iter::Int = 3,
    cheated::Bool = false,
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_gen::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    nthreads::Int = Threads.nthreads(),
) where {T<:AbstractFloat}
    nevents = length(particles_batch)
    nevents == 0 && return log_probs_batch

    npids = length(hypotheses)
    nchunks = min(max(1, nthreads), nevents)
    if nchunks <= 1
        ws_seq = UnbinnedPIDWorkspace(spectrum, T; max_pids = npids, n_gen = n_gen)
        _run_unbinned_pid_range!(
            log_probs_batch, particles_batch, hits_batch, 1:nevents, ws_seq,
            hypotheses, spectrum, mapper, background, n_iter, cheated,
            cfg, signal, context, detector, wedge, n_gen, random_points,
        )
    else
        chunk_size = cld(nevents, nchunks)
        Threads.@threads for c = 1:nchunks
            i_start = (c - 1) * chunk_size + 1
            i_end = min(nevents, c * chunk_size)
            if i_start <= i_end
                ws_thr = UnbinnedPIDWorkspace(spectrum, T; max_pids = npids, n_gen = n_gen)
                _run_unbinned_pid_range!(
                    log_probs_batch, particles_batch, hits_batch, i_start:i_end, ws_thr,
                    hypotheses, spectrum, mapper, background, n_iter, cheated,
                    cfg, signal, context, detector, wedge, n_gen, random_points,
                )
            end
        end
    end

    return log_probs_batch
end

function run_unbinned_pid_batch(
    particles_batch::AbstractVector{<:AbstractVector{Particle}},
    hits_batch::AbstractVector{<:AbstractVector},
    ::Type{T} = Float64;
    hypotheses::AbstractVector{<:Integer} = [211, 321, 2212],
    kwargs...,
)::Vector{Matrix{T}} where {T<:AbstractFloat}
    npids = length(hypotheses)
    log_probs_batch = [zeros(T, length(ps), npids) for ps in particles_batch]
    run_unbinned_pid_batch!(
        log_probs_batch, particles_batch, hits_batch;
        hypotheses = hypotheses, kwargs...,
    )
    return log_probs_batch
end
