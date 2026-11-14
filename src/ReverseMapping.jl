const DEFAULT_PARAM_SLOPE = (1.31035e+0, 2.11668e-2, 2.26502e-4, 2.88041e-6, 2.52544e-8)
const DEFAULT_PARAM_YPATH = (2.11098e+2, 1.36248e+0, -4.77368e-3, -8.17431e-5, -3.15632e-6)

struct ReverseMappingConfig
    nx_reflections::Int
    time_window::Float64
    is_blackened_bottom::Bool
    param_slope::NTuple{5,Float64}
    param_ypath::NTuple{5,Float64}
end

function ReverseMappingConfig(;
    nx_reflections::Int = 6,
    time_window::Float64 = 0.3,
    is_blackened_bottom::Bool = false,
    param_slope::NTuple{5,Float64} = DEFAULT_PARAM_SLOPE,
    param_ypath::NTuple{5,Float64} = DEFAULT_PARAM_YPATH,
)
    return ReverseMappingConfig(
        nx_reflections,
        time_window,
        is_blackened_bottom,
        param_slope,
        param_ypath,
    )
end

struct MappedPhoton
    time::Float64
    costhetac::Float64
    sinphic::Float64
    cosphic::Float64
    energy::Float64
    xdir::Float64
    ydir::Float64
    zdir::Float64
    pathlength::Float64
    energyprob::Float64
end

@inline function focus_slope(cfg::ReverseMappingConfig, yhit::Float64)::Float64
    return evalpoly(yhit, cfg.param_slope)
end

@inline function focus_ypath(
    cfg::ReverseMappingConfig,
    yhit::Float64,
    radiator::Radiator = RADIATOR,
    wedge::Wedge = WEDGE,
)::Float64
    return evalpoly(yhit, cfg.param_ypath) + radiator.height - wedge.offset
end

@inline function _x_offset(m::Int, xhit::Float64, x0::Float64, width::Float64)::Float64
    return isodd(m) ? (m * width - xhit - x0) : (m * width + xhit - x0)
end

@inline function _sincos_phi_angle(
    dx::Float64,
    dy::Float64,
    dz::Float64,
    sintheta::Float64,
    particle::Particle,
)::Tuple{Float64,Float64}
    if sintheta <= 0.0
        return (0.0, 1.0)
    end
    # Inverse rotation from lab frame to track Cherenkov frame using orthogonal rotMatrix^T
    R = particle.rotMatrix
    xprime = R[1, 1] * dx + R[2, 1] * dy + R[3, 1] * dz
    yprime = R[1, 2] * dx + R[2, 2] * dy + R[3, 2] * dz
    inv_sin = 1.0 / sintheta
    return (yprime * inv_sin, xprime * inv_sin)
end

@inline function _trace_center_photon(
    particle::Particle,
    beta::Float64,
    energy::Float64,
    sinphi::Float64,
    cosphi::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    context::PhotonContext,
)::Union{HitCoordinate,Nothing}
    nphase = spectrum_nphase(spectrum, energy)
    ngroup = spectrum_ngroup(spectrum, energy)
    if nphase <= 0.0 || ngroup <= 0.0 || beta * nphase <= 1.0
        return nothing
    end
    emission_z = 0.5 * context.radiator.depth
    photon = create_approximate_photon(
        particle,
        beta,
        nphase,
        ngroup,
        energy,
        emission_z,
        cosphi,
        sinphi,
        context,
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
    particle::Particle,
    beta::Float64,
    hit::HitCoordinate,
    energy::Float64,
    sinphi::Float64,
    cosphi::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    context::PhotonContext,
)::Float64
    dE = spectrum.dE
    dphi = 1.0e-6
    sindphi, cosdphi = sincos(dphi)

    h1 = _trace_center_photon(particle, beta, energy - dE, sinphi, cosphi, mapper, spectrum, context)
    h2 = _trace_center_photon(particle, beta, energy + dE, sinphi, cosphi, mapper, spectrum, context)
    h3 = _trace_center_photon(
        particle,
        beta,
        energy,
        sinphi * cosdphi - cosphi * sindphi,
        cosphi * cosdphi + sinphi * sindphi,
        mapper,
        spectrum,
        context,
    )
    h4 = _trace_center_photon(
        particle,
        beta,
        energy,
        sinphi * cosdphi + cosphi * sindphi,
        cosphi * cosdphi - sinphi * sindphi,
        mapper,
        spectrum,
        context,
    )

    dx_dE, dy_dE = _finite_diff_xy(hit, h1, h2, dE)
    dx_dphi, dy_dphi = _finite_diff_xy(hit, h3, h4, dphi)

    return abs(dx_dE * dy_dphi - dx_dphi * dy_dE)
end

@inline function _mapped_photon_pdf(
    particle::Particle,
    beta::Float64,
    hit::HitCoordinate,
    energy::Float64,
    rectime::Float64,
    sinphic::Float64,
    cosphic::Float64,
    pathlength::Float64,
    prob_energy::Float64,
    mapper::PhotonMapper,
    spectrum::PhotonSpectrum,
    signal::SignalParameters,
    context::PhotonContext,
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

    z0 = 0.5 * radiator.depth
    x0 = particle.xCoord + z0 * (particle.xDir / particle.zDir)
    y0 = particle.yCoord + z0 * (particle.yDir / particle.zDir) + 0.5 * radiator.height
    d0 = particle.pathlength + (z0 / particle.zDir)
    t0 = particle.t0 + d0 / (beta * clight)

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
            ix = _x_offset(m, hit.x, x0, radiator.width)
            pathlength = sqrt(ix * ix + iy_sq_factor)
            inv_path = 1.0 / pathlength

            dx = ix * inv_path
            dy = iy * inv_path
            dz = abs_iy_over_slope * inv_path

            costhetac = dx * particle.xDir + dy * particle.yDir + dz * particle.zDir
            if costhetac < min_costhetac || costhetac > max_costhetac
                continue
            end

            nphase = 1.0 / (beta * costhetac)
            if nphase <= spectrum.nphase_min || nphase >= spectrum.nphase_max
                continue
            end

            sinsqc = (n_air / nphase)^2
            if (1.0 - dz * dz) <= sinsqc
                continue
            end
            if y_reflected && ((1.0 - dy * dy) <= sinsqc)
                continue
            end
            x_reflected = abs(x0 + ix) > radiator.half_width
            if x_reflected && ((1.0 - dx * dx) <= sinsqc)
                continue
            end

            energy, vgroup, prob_en = propagation_properties_nocheck(spectrum, distribution, nphase)
            rectime = t0 + pathlength / vgroup
            if abs(rectime - hit.t) < cfg.time_window
                sinthetac = sqrt(max(0.0, 1.0 - costhetac * costhetac))
                sinphic, cosphic = _sincos_phi_angle(dx, dy, dz, sinthetac, particle)
                total_prob += _mapped_photon_pdf(
                    particle,
                    beta,
                    hit,
                    energy,
                    rectime,
                    sinphic,
                    cosphic,
                    pathlength,
                    prob_en,
                    mapper,
                    spectrum,
                    signal,
                    context,
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
        cfg, mapper, spectrum, signal, context,
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
) where {T<:AbstractFloat}
    radiator = context.radiator
    n_air = context.constants.N_AIR
    clight = context.constants.CLIGHT
    nhypos = length(betas)

    z0 = 0.5 * radiator.depth
    x0 = particle.xCoord + z0 * (particle.xDir / particle.zDir)
    y0 = particle.yCoord + z0 * (particle.yDir / particle.zDir) + 0.5 * radiator.height
    d0 = particle.pathlength + (z0 / particle.zDir)

    inv_slope_sq_plus_1 = 1.0 + 1.0 / (slope * slope)
    ny_branches = (cfg.is_blackened_bottom || mapper.blackened_bottom) ? 1 : 2
    nx = cfg.nx_reflections

    @inbounds for ky = 1:ny_branches
        iy = (ky == 1) ? (ypath - y0) : (-ypath - y0)
        y_reflected = (ky == 2)
        iy_sq_factor = (iy * iy) * inv_slope_sq_plus_1
        abs_iy_over_slope = abs(iy / slope)

        for m = -nx:nx
            ix = _x_offset(m, hit.x, x0, radiator.width)
            pathlength = sqrt(ix * ix + iy_sq_factor)
            inv_path = 1.0 / pathlength

            dx = ix * inv_path
            dy = iy * inv_path
            dz = abs_iy_over_slope * inv_path

            costhetac = dx * particle.xDir + dy * particle.yDir + dz * particle.zDir
            if costhetac <= 0.0 || costhetac >= 1.0
                continue
            end

            x_reflected = abs(x0 + ix) > radiator.half_width
            sinthetac = -1.0
            sinphic = 0.0
            cosphic = 1.0

            for h = 1:nhypos
                inv_norm = inv_norms[h]
                if inv_norm <= 0.0
                    continue
                end
                beta = betas[h]
                nphase = 1.0 / (beta * costhetac)
                if nphase <= spectrum.nphase_min || nphase >= spectrum.nphase_max
                    continue
                end

                sinsqc = (n_air / nphase)^2
                if (1.0 - dz * dz) <= sinsqc
                    continue
                end
                if y_reflected && ((1.0 - dy * dy) <= sinsqc)
                    continue
                end
                if x_reflected && ((1.0 - dx * dx) <= sinsqc)
                    continue
                end

                dist = distributions[h]
                energy, vgroup, prob_en = propagation_properties_nocheck(spectrum, dist, nphase)
                t0 = particle.t0 + d0 / (beta * clight)
                rectime = t0 + pathlength / vgroup

                if abs(rectime - hit.t) < cfg.time_window
                    if sinthetac < 0.0
                        sinthetac = sqrt(max(0.0, 1.0 - costhetac * costhetac))
                        sinphic, cosphic = _sincos_phi_angle(dx, dy, dz, sinthetac, particle)
                    end
                    pdf_val = _mapped_photon_pdf(
                        particle,
                        beta,
                        hit,
                        energy,
                        rectime,
                        sinphic,
                        cosphic,
                        pathlength,
                        prob_en,
                        mapper,
                        spectrum,
                        signal,
                        context,
                    )
                    if pdf_val > 0.0
                        sig_prob[hit_idx, track_idx, h] += T(pdf_val * inv_norm)
                    end
                end
            end
        end
    end
    return nothing
end

"""
    map_hit(particle, beta, distribution, hit, cfg, mapper, spectrum, context, wedge)::Vector{MappedPhoton}

Diagnostic helper returning all reconstructed `MappedPhoton` solutions for a single `HitCoordinate`.
"""
function map_hit(
    particle::Particle,
    beta::Float64,
    distribution::PhotonSpectrumDistribution,
    hit::HitCoordinate,
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    mapper::PhotonMapper = PhotonMapper(),
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    wedge::Wedge = WEDGE,
)::Vector{MappedPhoton}
    photons = MappedPhoton[]
    if !spectrum_above_threshold(distribution)
        return photons
    end

    radiator = context.radiator
    n_air = context.constants.N_AIR
    clight = context.constants.CLIGHT

    z0 = 0.5 * radiator.depth
    x0 = particle.xCoord + z0 * (particle.xDir / particle.zDir)
    y0 = particle.yCoord + z0 * (particle.yDir / particle.zDir) + 0.5 * radiator.height
    d0 = particle.pathlength + (z0 / particle.zDir)
    t0 = particle.t0 + d0 / (beta * clight)

    max_costhetac = 1.0 / (beta * spectrum.nphase_min)
    min_costhetac = 1.0 / (beta * spectrum.nphase_max)

    slope = focus_slope(cfg, hit.y)
    ypath = focus_ypath(cfg, hit.y, radiator, wedge)
    inv_slope_sq_plus_1 = 1.0 + 1.0 / (slope * slope)
    ny_branches = (cfg.is_blackened_bottom || mapper.blackened_bottom) ? 1 : 2
    nx = cfg.nx_reflections

    for ky = 1:ny_branches
        iy = (ky == 1) ? (ypath - y0) : (-ypath - y0)
        y_reflected = (ky == 2)
        for m = -nx:nx
            ix = _x_offset(m, hit.x, x0, radiator.width)
            pathlength = sqrt(ix * ix + (iy * iy) * inv_slope_sq_plus_1)
            dx = ix / pathlength
            dy = iy / pathlength
            dz = abs(iy / slope) / pathlength

            costhetac = dx * particle.xDir + dy * particle.yDir + dz * particle.zDir
            if costhetac < min_costhetac || costhetac > max_costhetac
                continue
            end
            nphase = 1.0 / (beta * costhetac)
            if nphase <= spectrum.nphase_min || nphase >= spectrum.nphase_max
                continue
            end
            sinsqc = (n_air / nphase)^2
            if (1.0 - dz * dz) <= sinsqc ||
               (y_reflected && (1.0 - dy * dy) <= sinsqc) ||
               (abs(x0 + ix) > radiator.half_width && (1.0 - dx * dx) <= sinsqc)
                continue
            end
            energy, vgroup, prob_en = propagation_properties_nocheck(spectrum, distribution, nphase)
            rectime = t0 + pathlength / vgroup
            if abs(rectime - hit.t) < cfg.time_window
                sinthetac = sqrt(max(0.0, 1.0 - costhetac * costhetac))
                sinphic, cosphic = _sincos_phi_angle(dx, dy, dz, sinthetac, particle)
                push!(
                    photons,
                    MappedPhoton(
                        rectime,
                        costhetac,
                        sinphic,
                        cosphic,
                        energy,
                        dx,
                        dy,
                        dz,
                        pathlength,
                        prob_en,
                    ),
                )
            end
        end
    end
    return photons
end
