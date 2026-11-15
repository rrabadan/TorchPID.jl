struct PIDAlg
    nModules::Int

    pids::Vector{Int}
    momentumRange::Tuple{Float64,Float64}
    checkTkAcceptance::Bool
    checkTkQuality::Bool
    smearT0::Bool

    blackenedBottom::Bool
    removeFocusPhotons::Bool
    photonTrackMatching::Bool
end

function PIDAlg(;
    nModules::Int = 18,
    pids::Vector{Int} = [211, 321, 2212],
    momentumRange::Tuple{Float64,Float64} = (1.0, 20.0),
    checkTkAcceptance::Bool = true,
    checkTkQuality::Bool = false,
    smearT0::Bool = false,
    blackenedBottom::Bool = false,
    removeFocusPhotons::Bool = true,
    photonTrackMatching::Bool = false,
)::PIDAlg
    return PIDAlg(
        nModules,
        pids,
        momentumRange,
        checkTkAcceptance,
        checkTkQuality,
        smearT0,
        blackenedBottom,
        removeFocusPhotons,
        photonTrackMatching,
    )
end

"""
    ModulePIDWorkspace{T<:AbstractFloat}

Pre-allocated, zero-allocation workspace for 18-module TORCH event partitioning and
load-balanced module-parallel Unbinned PID reconstruction on CPU.
Supports both 0-based (`0:(nModules-1)`) and 1-based (`1:nModules`) module IDs via
`nModules + 1` pre-allocated module slots (`slot = Int(mod_id) + 1`).
"""
struct ModulePIDWorkspace{T<:AbstractFloat}
    n_modules::Int
    module_particles::Vector{Vector{Particle}}
    module_hits::Vector{Vector{HitCoordinate}}
    module_log_probs_buf::Vector{Vector{T}}
    active_slots::Vector{Int}
    pid_workspaces::Vector{UnbinnedPIDWorkspace{T}}
end

function ModulePIDWorkspace(
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    ::Type{T} = Float64;
    n_modules::Int = 18,
    max_tracks::Int = 16,
    max_hits::Int = 750,
    max_pids::Int = 3,
    n_gen::Int = 1000,
) where {T<:AbstractFloat}
    n_slots = n_modules + 1
    module_particles = [sizehint!(Particle[], max_tracks) for _ = 1:n_slots]
    module_hits = [sizehint!(HitCoordinate[], max_hits) for _ = 1:n_slots]
    module_log_probs_buf = [zeros(T, max_tracks * max_pids) for _ = 1:n_slots]
    active_slots = sizehint!(Int[], n_slots)
    pid_workspaces = [
        UnbinnedPIDWorkspace(
            spectrum, T;
            max_tracks = max_tracks,
            max_hits = max_hits,
            max_pids = max_pids,
            n_gen = n_gen,
        ) for _ = 1:n_slots
    ]
    return ModulePIDWorkspace{T}(
        n_modules,
        module_particles,
        module_hits,
        module_log_probs_buf,
        active_slots,
        pid_workspaces,
    )
end

@inline function _module_slot(mod_id::Integer, n_modules::Int)::Int
    m = Int(mod_id)
    return (0 <= m <= n_modules) ? (m + 1) : 0
end

@inline function _module_id_from_slot(slot::Int)::Int
    return slot - 1
end

@inline function _module_log_probs_view!(
    ws::ModulePIDWorkspace{T},
    slot::Int,
    ntracks::Int,
    npids::Int,
) where {T<:AbstractFloat}
    buf = @inbounds ws.module_log_probs_buf[slot]
    needed = ntracks * npids
    length(buf) < needed && resize!(buf, needed)
    return reshape(@view(buf[1:needed]), ntracks, npids)
end

@inline function _sort_active_slots_desc!(ws::ModulePIDWorkspace)
    slots = ws.active_slots
    n = length(slots)
    @inbounds for i = 2:n
        key = slots[i]
        key_work = length(ws.module_particles[key]) * length(ws.module_hits[key])
        j = i - 1
        while j >= 1
             sj = slots[j]
            if length(ws.module_particles[sj]) * length(ws.module_hits[sj]) < key_work
                slots[j + 1] = sj
                j -= 1
            else
                break
            end
        end
        slots[j + 1] = key
    end
    return nothing
end

function select_photon_hit(
    pidalg::PIDAlg,
    radiator::Radiator,
    wedge::Wedge,
    detector::Detector,
    mask::Mask,
    photon::PhotonHit,
)::Bool
    return photon_on_detector(photon, detector) &&
           !(pidalg.removeFocusPhotons && photon_from_focus(photon, radiator, wedge)) &&
           !is_photon_masked(photon, mask) &&
           !(pidalg.blackenedBottom && photon.isDownward) &&
           (detector.t_window_min < photon.arrivalTime < detector.t_window_max)
end

"""
    partition_by_module!(
        ws::ModulePIDWorkspace,
        pidalg::PIDAlg,
        dht::DetectorHitTester,
        particles::AbstractVector{Particle},
        photons::AbstractVector{PhotonHit};
        radiator = RADIATOR,
        wedge = WEDGE,
        detector = DETECTOR,
        mask = MASK,
    )

Zero-allocation partitioning of event `particles` and `photons` into the 18 module buffers
of `ws`. Populates `ws.active_slots` with all modules having both tracks and hits, sorted in
descending order of work (`ntracks * nhits`) for LPT (Longest Processing Time first) scheduling.
"""
function partition_by_module!(
    ws::ModulePIDWorkspace,
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    particles::AbstractVector{Particle},
    photons::AbstractVector{PhotonHit};
    radiator::Radiator = RADIATOR,
    wedge::Wedge = WEDGE,
    detector::Detector = DETECTOR,
    mask::Mask = MASK,
)
    n_modules = min(pidalg.nModules, ws.n_modules)
    clear_modules!(ws)

    pmin, pmax = pidalg.momentumRange
    pids = pidalg.pids
    check_acc = pidalg.checkTkAcceptance

    @inbounds for i in eachindex(particles)
        p = particles[i]
        slot = _module_slot(p.moduleId, n_modules)
        if slot > 0 && (!check_acc || (pmin <= p.pMag <= pmax)) && (abs(p.pid) in pids)
            push!(ws.module_particles[slot], p)
        end
    end

    @inbounds for i in eachindex(photons)
        ph = photons[i]
        slot = _module_slot(ph.moduleId, n_modules)
        if slot > 0 &&
           !isempty(ws.module_particles[slot]) &&
           select_photon_hit(pidalg, radiator, wedge, detector, mask, ph) &&
           test_photon(dht, ph.energy)
            push!(ws.module_hits[slot], HitCoordinate(ph.localX, ph.localY, ph.arrivalTime))
        end
    end

    return update_active_modules!(ws)
end

"""
    clear_modules!(ws::ModulePIDWorkspace)

Empties all module track/hit buffers (keeping their capacity) and the active-module list.
"""
function clear_modules!(ws::ModulePIDWorkspace)
    @inbounds for s in eachindex(ws.module_particles)
        empty!(ws.module_particles[s])
        empty!(ws.module_hits[s])
    end
    empty!(ws.active_slots)
    return ws
end

"""
    set_module!(ws::ModulePIDWorkspace, mod_id, particles, hits)

Copies pre-selected `particles` and `hits` of module `mod_id` into the workspace buffers
(zero-allocation once buffer capacity is reached). Call `update_active_modules!(ws)` after
filling all modules.
"""
function set_module!(
    ws::ModulePIDWorkspace,
    mod_id::Integer,
    particles::AbstractVector{Particle},
    hits::AbstractVector{HitCoordinate},
)
    slot = _module_slot(mod_id, ws.n_modules)
    slot == 0 && throw(ArgumentError("module id $mod_id outside 0:$(ws.n_modules)"))
    mp = ws.module_particles[slot]
    mh = ws.module_hits[slot]
    resize!(mp, length(particles))
    resize!(mh, length(hits))
    copyto!(mp, particles)
    copyto!(mh, hits)
    return ws
end

"""
    update_active_modules!(ws::ModulePIDWorkspace)

Rebuilds `ws.active_slots` (modules with ≥1 track and ≥1 hit) sorted largest-work-first
(`ntracks * nhits`, LPT scheduling order).
"""
function update_active_modules!(ws::ModulePIDWorkspace)
    empty!(ws.active_slots)
    @inbounds for s in eachindex(ws.module_particles)
        if !isempty(ws.module_particles[s]) && !isempty(ws.module_hits[s])
            push!(ws.active_slots, s)
        end
    end
    _sort_active_slots_desc!(ws)
    return ws.active_slots
end

"""
    module_log_probs(ws::ModulePIDWorkspace, mod_id; npids = 3)

Returns a view `(ntracks, npids)` of the log-likelihoods computed for module `mod_id`.
"""
function module_log_probs(ws::ModulePIDWorkspace, mod_id::Integer; npids::Int = 3)
    slot = _module_slot(mod_id, ws.n_modules)
    return _module_log_probs_view!(ws, slot, length(ws.module_particles[slot]), npids)
end

function find_particles(
    pidalg::PIDAlg,
    particles::Vector{Particle},
)::Dict{Int,Vector{Particle}}
    module_particles = Dict{Int,Vector{Particle}}()
    pmin, pmax = pidalg.momentumRange

    for p in particles
        mod_id = Int(p.moduleId)
        if (0 <= mod_id <= pidalg.nModules) &&
           (abs(p.pid) in pidalg.pids) &&
           (!pidalg.checkTkAcceptance || (pmin <= p.pMag <= pmax))
            vec = get!(() -> Particle[], module_particles, mod_id)
            push!(vec, p)
        end
    end
    return module_particles
end

function find_coordinates(
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    photons::Vector{PhotonHit},
)::Dict{Int,Vector{HitCoordinate}}
    module_hits = Dict{Int,Vector{HitCoordinate}}()

    for p in photons
        mod_id = Int(p.moduleId)
        if (0 <= mod_id <= pidalg.nModules) &&
           select_photon_hit(pidalg, RADIATOR, WEDGE, DETECTOR, MASK, p) &&
           test_photon(dht, p.energy)
            vec = get!(() -> HitCoordinate[], module_hits, mod_id)
            push!(vec, HitCoordinate(p.localX, p.localY, p.arrivalTime))
        end
    end
    return module_hits
end

function _run_modules_dynamic!(
    ws::ModulePIDWorkspace{T},
    pids::Vector{Int},
    spectrum::PhotonSpectrum,
    mapper::PhotonMapper,
    background::AbstractBackgroundModel,
    cfg::ReverseMappingConfig,
    signal::SignalParameters,
    context::PhotonContext,
    detector::Detector,
    wedge::Wedge,
    n_iter::Int,
    cheated::Bool,
    n_gen::Int,
    random_points::Union{AbstractVector{RandomPoint},Nothing},
    nthreads::Int,
) where {T<:AbstractFloat}
    n_active = length(ws.active_slots)
    npids = length(pids)
    nworkers = min(nthreads, n_active)
    next_idx = Threads.Atomic{Int}(1)
    allow_nested = (nthreads >= Threads.nthreads())

    Threads.@threads for _ = 1:nworkers
        while true
            idx = Threads.atomic_add!(next_idx, 1)
            idx > n_active && break
            slot = @inbounds ws.active_slots[idx]
            mod_tracks = @inbounds ws.module_particles[slot]
            mod_hits = @inbounds ws.module_hits[slot]
            lp_view = _module_log_probs_view!(ws, slot, length(mod_tracks), npids)
            work = length(mod_tracks) * length(mod_hits)
            intra_threaded = allow_nested && ((work >= 1200) || (n_active < nthreads && work >= 400))
            run_unbinned_pid!(
                lp_view,
                mod_tracks,
                mod_hits;
                hypotheses = pids,
                spectrum = spectrum,
                mapper = mapper,
                background = background,
                n_iter = n_iter,
                cheated = cheated,
                workspace = @inbounds(ws.pid_workspaces[slot]),
                cfg = cfg,
                signal = signal,
                context = context,
                detector = detector,
                wedge = wedge,
                n_gen = n_gen,
                random_points = random_points,
                threaded = intra_threaded,
            )
        end
    end
    return nothing
end

"""
    run_module_pid!(
        ws::ModulePIDWorkspace{T},
        pids::AbstractVector{<:Integer};
        spectrum, mapper, background, cfg, signal, context, detector, wedge,
        n_iter = 3, cheated = false, n_gen = 1000, random_points = nothing,
        nthreads = Threads.nthreads(),
    )

CPU module dispatcher: runs Unbinned PID on every module listed in `ws.active_slots`
(already partitioned, sorted largest-work-first). Scheduling policy:
- `nthreads <= 1`: sequential loop over modules.
- one active module: intra-module threading (hits + (track, hypothesis) efficiencies).
- otherwise: dynamic work-stealing over modules (atomic counter, LPT order), with nested
  intra-module threading for heavy modules (`ntracks * nhits >= 1200`).
Results are read with `module_log_probs(ws, mod_id)`.
"""
function run_module_pid!(
    ws::ModulePIDWorkspace{T},
    pids::Vector{Int};
    spectrum::PhotonSpectrum = PhotonSpectrum(),
    mapper::PhotonMapper = PhotonMapper(),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    cfg::ReverseMappingConfig = ReverseMappingConfig(),
    signal::SignalParameters = SIGNAL,
    context::PhotonContext = PhotonContext(RADIATOR, CONSTANTS),
    detector::Detector = DETECTOR,
    wedge::Wedge = WEDGE,
    n_iter::Int = 3,
    cheated::Bool = false,
    n_gen::Int = 1000,
    random_points::Union{AbstractVector{RandomPoint},Nothing} = nothing,
    nthreads::Int = Threads.nthreads(),
) where {T<:AbstractFloat}
    n_active = length(ws.active_slots)
    n_active == 0 && return ws

    npids = length(pids)
    if nthreads <= 1 || n_active == 1
        intra = nthreads > 1
        @inbounds for idx = 1:n_active
            slot = ws.active_slots[idx]
            mod_tracks = ws.module_particles[slot]
            mod_hits = ws.module_hits[slot]
            lp_view = _module_log_probs_view!(ws, slot, length(mod_tracks), npids)
            run_unbinned_pid!(
                lp_view,
                mod_tracks,
                mod_hits;
                hypotheses = pids,
                spectrum = spectrum,
                mapper = mapper,
                background = background,
                n_iter = n_iter,
                cheated = cheated,
                workspace = ws.pid_workspaces[slot],
                cfg = cfg,
                signal = signal,
                context = context,
                detector = detector,
                wedge = wedge,
                n_gen = n_gen,
                random_points = random_points,
                threaded = intra,
            )
        end
    else
        _run_modules_dynamic!(
            ws, pids, spectrum, mapper, background,
            cfg, signal, context, detector, wedge,
            n_iter, cheated, n_gen, random_points, nthreads,
        )
    end
    return ws
end

"""
    run_algorithm!(
        ws::ModulePIDWorkspace{T},
        pidalg::PIDAlg,
        dht::DetectorHitTester,
        particles::AbstractVector{Particle},
        photons::AbstractVector{PhotonHit};
        spectrum = PhotonSpectrum(dht),
        mapper = PhotonMapper(blackened_bottom = pidalg.blackenedBottom),
        background = FlatBackground(DETECTOR),
        cfg = ReverseMappingConfig(is_blackened_bottom = pidalg.blackenedBottom),
        nthreads = Threads.nthreads(),
        kwargs...,
    )

Zero-allocation 18-module Unbinned PID reconstruction for a single TORCH event:
`partition_by_module!` followed by the CPU module dispatcher `run_module_pid!`.
"""
function run_algorithm!(
    ws::ModulePIDWorkspace{T},
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    particles::AbstractVector{Particle},
    photons::AbstractVector{PhotonHit};
    spectrum::PhotonSpectrum = PhotonSpectrum(dht),
    mapper::PhotonMapper = PhotonMapper(blackened_bottom = pidalg.blackenedBottom),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    cfg::ReverseMappingConfig = ReverseMappingConfig(is_blackened_bottom = pidalg.blackenedBottom),
    kwargs...,
) where {T<:AbstractFloat}
    partition_by_module!(ws, pidalg, dht, particles, photons)
    return run_module_pid!(
        ws, pidalg.pids;
        spectrum = spectrum, mapper = mapper, background = background, cfg = cfg, kwargs...,
    )
end

function run_algorithm(
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    particles::Vector{Particle},
    photons::Vector{PhotonHit};
    spectrum::PhotonSpectrum = PhotonSpectrum(dht),
    mapper::PhotonMapper = PhotonMapper(blackened_bottom = pidalg.blackenedBottom),
    background::AbstractBackgroundModel = FlatBackground(DETECTOR),
    cfg::ReverseMappingConfig = ReverseMappingConfig(is_blackened_bottom = pidalg.blackenedBottom),
    workspace::Union{ModulePIDWorkspace{Float64},Nothing} = nothing,
    kwargs...,
)::Dict{Int,Tuple{Vector{Particle},Matrix{Float64}}}
    ws = workspace === nothing ?
        ModulePIDWorkspace(spectrum, Float64; n_modules = pidalg.nModules, max_pids = length(pidalg.pids)) :
        workspace

    run_algorithm!(
        ws, pidalg, dht, particles, photons;
        spectrum = spectrum,
        mapper = mapper,
        background = background,
        cfg = cfg,
        kwargs...,
    )

    npids = length(pidalg.pids)
    results = Dict{Int,Tuple{Vector{Particle},Matrix{Float64}}}()
    for slot in ws.active_slots
        mod_id = _module_id_from_slot(slot)
        mod_tracks = copy(ws.module_particles[slot])
        lp_mat = Matrix{Float64}(_module_log_probs_view!(ws, slot, length(mod_tracks), npids))
        results[mod_id] = (mod_tracks, lp_mat)
    end
    return results
end
