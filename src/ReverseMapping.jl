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
    dx::Float64, dy::Float64, dz::Float64, sintheta::Float64, particle::Particle,
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

@inline function _track_emission_midpoint(
    particle::Particle, radiator::Radiator,
)::Tuple{Float64,Float64,Float64}
    z0 = 0.5 * radiator.depth
    inv_zdir = 1.0 / particle.zDir
    x0 = particle.xCoord + z0 * (particle.xDir * inv_zdir)
    y0 = particle.yCoord + z0 * (particle.yDir * inv_zdir) + 0.5 * radiator.height
    d0 = particle.pathlength + z0 * inv_zdir
    return (x0, y0, d0)
end

@inline function _reconstruct_branch_path(
    m::Int, xhit::Float64, x0::Float64, iy::Float64,
    iy_sq_factor::Float64, abs_iy_over_slope::Float64,
    particle::Particle, radiator::Radiator,
)::Tuple{Float64,Float64,Float64,Float64,Float64,Bool}
    ix = _x_offset(m, xhit, x0, radiator.width)
    pathlength = sqrt(ix * ix + iy_sq_factor)
    inv_path = 1.0 / pathlength
    dx = ix * inv_path
    dy = iy * inv_path
    dz = abs_iy_over_slope * inv_path
    costhetac = dx * particle.xDir + dy * particle.yDir + dz * particle.zDir
    x_reflected = abs(x0 + ix) > radiator.half_width
    return (pathlength, dx, dy, dz, costhetac, x_reflected)
end

@inline function _evaluate_branch_hypothesis(
    beta::Float64, dist::PhotonSpectrumDistribution, costhetac::Float64,
    dx::Float64, dy::Float64, dz::Float64, pathlength::Float64,
    x_reflected::Bool, y_reflected::Bool, t0_particle::Float64, d0::Float64,
    thit::Float64, time_window::Float64, n_air::Float64, clight::Float64,
    spectrum::PhotonSpectrum,
)::Tuple{Bool,Float64,Float64,Float64}
    nphase = 1.0 / (beta * costhetac)
    if nphase <= spectrum.nphase_min || nphase >= spectrum.nphase_max
        return (false, 0.0, 0.0, 0.0)
    end
    sinsqc = (n_air / nphase)^2
    if (1.0 - dz * dz) <= sinsqc ||
       (y_reflected && (1.0 - dy * dy) <= sinsqc) ||
       (x_reflected && (1.0 - dx * dx) <= sinsqc)
        return (false, 0.0, 0.0, 0.0)
    end
    energy, vgroup, prob_en = propagation_properties_nocheck(spectrum, dist, nphase)
    t0 = t0_particle + d0 / (beta * clight)
    rectime = t0 + pathlength / vgroup
    if abs(rectime - thit) >= time_window
        return (false, 0.0, 0.0, 0.0)
    end
    return (true, energy, rectime, prob_en)
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

    x0, y0, d0 = _track_emission_midpoint(particle, radiator)
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
                push!(
                    photons,
                    MappedPhoton(
                        rectime, costhetac, sinphic, cosphic, energy,
                        dx, dy, dz, pathlength, prob_en,
                    ),
                )
            end
        end
    end
    return photons
end
