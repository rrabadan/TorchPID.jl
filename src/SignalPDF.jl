@inline function _trace_center_photon(
    particle::Particle, beta::Float64, energy::Float64,
    sinphi::Float64, cosphi::Float64, mapper::PhotonMapper,
    spectrum::PhotonSpectrum, context::PhotonContext,
)::Union{HitCoordinate,Nothing}
    nphase = spectrum_nphase(spectrum, energy)
    ngroup = spectrum_ngroup(spectrum, energy)
    if nphase <= 0.0 || ngroup <= 0.0 || beta * nphase <= 1.0
        return nothing
    end
    emission_z = 0.5 * context.radiator.depth
    photon = create_approximate_photon(
        particle, beta, nphase, ngroup, energy, emission_z, cosphi, sinphi, context,
    )
    toffset = _get_time_offset(particle, beta, photon.zpos, particle.t0, context)
    return trace_photon(mapper, photon, toffset)
end

@inline function _finite_diff_xy(
    mid::HitCoordinate,
    low::Union{HitCoordinate,Nothing},
    upp::Union{HitCoordinate,Nothing},
    step::Float64,
)::Tuple{Float64,Float64}
    low_ok = (low !== nothing) && (low.t > 0.0)
    upp_ok = (upp !== nothing) && (upp.t > 0.0)
    if low_ok && upp_ok
        inv_2s = 0.5 / step
        return ((upp.x - low.x) * inv_2s, (upp.y - low.y) * inv_2s)
    elseif low_ok
        inv_s = 1.0 / step
        return ((mid.x - low.x) * inv_s, (mid.y - low.y) * inv_s)
    elseif upp_ok
        inv_s = 1.0 / step
        return ((upp.x - mid.x) * inv_s, (upp.y - mid.y) * inv_s)
    else
        return (0.0, 0.0)
    end
end

@inline function photon_jacobian(
    particle::Particle, beta::Float64, hit::HitCoordinate,
    energy::Float64, sinphi::Float64, cosphi::Float64,
    mapper::PhotonMapper, spectrum::PhotonSpectrum, context::PhotonContext,
)::Float64
    dE = spectrum.dE
    dphi = 1.0e-6
    sindphi, cosdphi = sincos(dphi)

    h1 = _trace_center_photon(particle, beta, energy - dE, sinphi, cosphi, mapper, spectrum, context)
    h2 = _trace_center_photon(particle, beta, energy + dE, sinphi, cosphi, mapper, spectrum, context)
    h3 = _trace_center_photon(
        particle, beta, energy,
        sinphi * cosdphi - cosphi * sindphi, cosphi * cosdphi + sinphi * sindphi,
        mapper, spectrum, context,
    )
    h4 = _trace_center_photon(
        particle, beta, energy,
        sinphi * cosdphi + cosphi * sindphi, cosphi * cosdphi - sinphi * sindphi,
        mapper, spectrum, context,
    )

    dx_dE, dy_dE = _finite_diff_xy(hit, h1, h2, dE)
    dx_dphi, dy_dphi = _finite_diff_xy(hit, h3, h4, dphi)

    return abs(dx_dE * dy_dphi - dx_dphi * dy_dE)
end

@inline function _mapped_photon_pdf(
    particle::Particle, beta::Float64, hit::HitCoordinate,
    energy::Float64, rectime::Float64, sinphic::Float64, cosphic::Float64,
    pathlength::Float64, prob_energy::Float64,
    mapper::PhotonMapper, spectrum::PhotonSpectrum,
    signal::SignalParameters, context::PhotonContext,
)::Float64
    if energy < spectrum.emin || energy > spectrum.emax || prob_energy <= 0.0
        return 0.0
    end
    det_jac = photon_jacobian(particle, beta, hit, energy, sinphic, cosphic, mapper, spectrum, context)
    if det_jac <= 0.0
        return 0.0
    end
    sigma_t = signal.effective_constant + signal.effective_linear * pathlength
    z = (rectime - hit.t) / sigma_t
    prob_time = exp(-0.5 * z * z) / (sqrt(2.0 * π) * sigma_t)
    return (prob_time * prob_energy) / (2.0 * π * det_jac)
end

"""
    hit_signal_probability(
        particle, beta, distribution, hit, cfg, mapper, spectrum, signal, context, wedge;
        slope = focus_slope(cfg, hit.y),
        ypath = focus_ypath(cfg, hit.y, context.radiator, wedge),
    )::Float64

Stateless, zero-allocation evaluation of the unbinned signal PDF for a single `(particle, beta, hit)`
pair across all 26 reflection branches in registers.
"""
@inline function hit_signal_probability(
    particle::Particle,
    beta::Float64,
    distribution::PhotonSpectrumDistribution,
    hit::HitCoordinate,
    cfg::ReverseMappingConfig,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    wedge::Wedge = WEDGE;
    slope::Float64 = focus_slope(cfg, hit.y),
    ypath::Float64 = focus_ypath(cfg, hit.y, context.radiator, wedge),
)::Float64
    if !spectrum_above_threshold(distribution)
        return 0.0
    end

    radiator = context.radiator
    n_air = context.constants.N_AIR
    clight = context.constants.CLIGHT

    x0, y0, d0 = _track_emission_midpoint(particle, radiator)
    max_costhetac = 1.0 / (beta * spectrum.nphase_min)
    min_costhetac = 1.0 / (beta * spectrum.nphase_max)

    inv_slope_sq_plus_1 = 1.0 + 1.0 / (slope * slope)
    ny_branches = (cfg.is_blackened_bottom || mapper.blackened_bottom) ? 1 : 2
    nx = cfg.nx_reflections
    total_prob = 0.0

    @inbounds for ky = 1:ny_branches
        iy = (ky == 1) ? (ypath - y0) : (-ypath - y0)
        y_reflected = (ky == 2)
        iy_sq_factor = (iy * iy) * inv_slope_sq_plus_1
        abs_iy_over_slope = abs(iy / slope)

        for m = -nx:nx
            pathlength, dx, dy, dz, costhetac, x_reflected = _reconstruct_branch_path(
                m, hit.x, x0, iy, iy_sq_factor, abs_iy_over_slope, particle, radiator,
            )
            if costhetac < min_costhetac || costhetac > max_costhetac
                continue
            end

            valid, energy, rectime, prob_en = _evaluate_branch_hypothesis(
                beta, distribution, costhetac, dx, dy, dz, pathlength,
                x_reflected, y_reflected, particle.t0, d0, hit.t,
                cfg.time_window, n_air, clight, spectrum,
            )
            if valid
                sinthetac = sqrt(max(0.0, 1.0 - costhetac * costhetac))
                sinphic, cosphic = _sincos_phi_angle(dx, dy, dz, sinthetac, particle)
                total_prob += _mapped_photon_pdf(
                    particle, beta, hit, energy, rectime, sinphic, cosphic,
                    pathlength, prob_en, mapper, spectrum, signal, context,
                )
            end
        end
    end

    return total_prob
end

"""
    accumulate_hit_signal_probs!(
        sig_prob, hit_idx, track_idx, hit, slope, ypath,
        particle, betas, distributions, inv_norms,
        cfg, mapper, spectrum, signal, context, weight = 1.0,
    )

Fused multi-hypothesis kernel: computes the β-independent photon path `(dx, dy, dz, pathlength, costhetac)`
once per reflection branch in registers and immediately evaluates all mass hypotheses `h in eachindex(betas)`.
"""
@inline function accumulate_hit_signal_probs!(
    sig_prob::AbstractArray{T,3},
    hit_idx::Int,
    track_idx::Int,
    hit::HitCoordinate,
    slope::Float64,
    ypath::Float64,
    particle::Particle,
    betas::AbstractVector{Float64},
    distributions::AbstractVector{PhotonSpectrumDistribution},
    inv_norms::AbstractVector{Float64},
    cfg::ReverseMappingConfig,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters,
    context::PhotonContext,
    weight::Float64 = 1.0,
) where {T<:AbstractFloat}
    radiator = context.radiator
    n_air = context.constants.N_AIR
    clight = context.constants.CLIGHT
    nhypos = length(betas)

    x0, y0, d0 = _track_emission_midpoint(particle, radiator)
    inv_slope_sq_plus_1 = 1.0 + 1.0 / (slope * slope)
    ny_branches = (cfg.is_blackened_bottom || mapper.blackened_bottom) ? 1 : 2
    nx = cfg.nx_reflections

    @inbounds for ky = 1:ny_branches
        iy = (ky == 1) ? (ypath - y0) : (-ypath - y0)
        y_reflected = (ky == 2)
        iy_sq_factor = (iy * iy) * inv_slope_sq_plus_1
        abs_iy_over_slope = abs(iy / slope)

        for m = -nx:nx
            pathlength, dx, dy, dz, costhetac, x_reflected = _reconstruct_branch_path(
                m, hit.x, x0, iy, iy_sq_factor, abs_iy_over_slope, particle, radiator,
            )
            if costhetac <= 0.0 || costhetac >= 1.0
                continue
            end

            sinthetac = -1.0
            sinphic = 0.0
            cosphic = 1.0

            for h = 1:nhypos
                inv_norm = inv_norms[h]
                if inv_norm <= 0.0
                    continue
                end
                beta = betas[h]
                valid, energy, rectime, prob_en = _evaluate_branch_hypothesis(
                    beta, distributions[h], costhetac, dx, dy, dz, pathlength,
                    x_reflected, y_reflected, particle.t0, d0, hit.t,
                    cfg.time_window, n_air, clight, spectrum,
                )
                if valid
                    if sinthetac < 0.0
                        sinthetac = sqrt(max(0.0, 1.0 - costhetac * costhetac))
                        sinphic, cosphic = _sincos_phi_angle(dx, dy, dz, sinthetac, particle)
                    end
                    pdf_val = _mapped_photon_pdf(
                        particle, beta, hit, energy, rectime, sinphic, cosphic,
                        pathlength, prob_en, mapper, spectrum, signal, context,
                    )
                    if pdf_val > 0.0
                        sig_prob[hit_idx, track_idx, h] += T(weight * pdf_val * inv_norm)
                    end
                end
            end
        end
    end
    return nothing
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

    radiator = context.radiator
    sl_min = focus_slope(cfg, ymin)
    yp_min = focus_ypath(cfg, ymin, radiator, wedge)
    sl_max = focus_slope(cfg, ymax)
    yp_max = focus_ypath(cfg, ymax, radiator, wedge)
    sl_ctr = focus_slope(cfg, center.y)
    yp_ctr = focus_ypath(cfg, center.y, radiator, wedge)

    xs = (xmin, xmin, xmax, xmax, center.x, center.x, xmin, xmax, center.x)
    ys = (ymin, ymax, ymin, ymax, ymin, ymax, center.y, center.y, center.y)
    sls = (sl_min, sl_max, sl_min, sl_max, sl_min, sl_max, sl_ctr, sl_ctr, sl_ctr)
    yps = (yp_min, yp_max, yp_min, yp_max, yp_min, yp_max, yp_ctr, yp_ctr, yp_ctr)

    area_factor = 0.25 * detector.x_size * detector.y_size * detector.t_bin

    @inbounds for q = 1:9
        sub_hit = HitCoordinate(xs[q], ys[q], center.t)
        w = PIXEL_QUAD_WEIGHTS[q] * area_factor
        accumulate_hit_signal_probs!(
            sig_prob, pixel_idx, track_idx, sub_hit, sls[q], yps[q],
            particle, betas, distributions, inv_norms,
            cfg, mapper, spectrum, signal, context, w,
        )
    end
    return nothing
end
