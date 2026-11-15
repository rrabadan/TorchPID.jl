#
# 18-module TORCH event-stream benchmark for the CPU module dispatcher (`run_module_pid!`).
#
# Synthetic LHCb-like occupancy: per event, module m receives Poisson(λ_m) tracks
# (central modules up to 15) and its Cherenkov signal hits (`project_pattern`) plus uniform
# background hits so that central modules reach ~700 hits.
#
# Two execution modes are measured:
#   A. Latency mode    — events processed strictly one after another; the 18 modules of each
#                        event are parallelised (dynamic LPT dispatch + nested heavy-module threading).
#   B. Throughput mode — a window of consecutive events is processed concurrently, one event per
#                        task, modules run sequentially inside each task (no nested threading).
#
# Usage:
#   julia --project=. --threads=auto alg/benchmark_module_pid.jl [--events N] [--trials T]
#         [--thread-list 1,2,4,8,16] [--output benchmark_module_pid_report.md]
#
using Dates
using Printf
using Random
using Statistics
using TorchPID

const N_MODULES = 18
# Mean tracks per module (central modules near the beam pipe are most occupied)
const TRACK_PROFILE = (1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 8.5, 12.0, 15.0, 15.0, 12.0, 8.5, 6.0, 4.5, 3.0, 2.0, 1.5, 1.0)
const MAX_TRACKS = 15
const MAX_HITS = 700
const BKG_HITS_PER_TRACK = 15.0

struct ModuleEvent
    particles::Vector{Vector{Particle}}       # per module (index = module id 1:18)
    hits::Vector{Vector{HitCoordinate}}
end

function parse_commandline(args::Vector{String} = ARGS)
    opts = Dict{String,Any}(
        "events" => 200,
        "trials" => 3,
        "thread-list" => "1,2,4,8,16",
        "output" => "benchmark_module_pid_report.md",
        "seed" => 42,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if (arg == "--events" || arg == "-n") && i < length(args)
            opts["events"] = parse(Int, args[i + 1]); i += 2
        elseif (arg == "--trials" || arg == "-t") && i < length(args)
            opts["trials"] = parse(Int, args[i + 1]); i += 2
        elseif arg == "--thread-list" && i < length(args)
            opts["thread-list"] = args[i + 1]; i += 2
        elseif (arg == "--output" || arg == "-o") && i < length(args)
            opts["output"] = args[i + 1]; i += 2
        elseif arg == "--seed" && i < length(args)
            opts["seed"] = parse(Int, args[i + 1]); i += 2
        else
            error("Unknown or incomplete argument: $arg")
        end
    end
    return opts
end

function parse_thread_list(spec::String, max_threads::Int)::Vector{Int}
    req = [parse(Int, strip(s)) for s in split(spec, ",") if !isempty(strip(s))]
    out = unique(filter(t -> 1 <= t <= max_threads, req))
    1 in out || pushfirst!(out, 1)
    max_threads in out || push!(out, max_threads)
    return sort!(out)
end

function poisson(rng::AbstractRNG, λ::Float64)::Int
    L = exp(-λ); k = 0; p = 1.0
    while true
        p *= rand(rng)
        p <= L && return k
        k += 1
    end
end

function generate_event_stream(nevents::Int, rng::AbstractRNG, spectrum, mapper, context)
    tb = TestBeamSimulator()
    pool = generate_particles(tb, context.radiator, context.constants, 20_000)
    cursor = 0
    events = Vector{ModuleEvent}(undef, nevents)
    for ev = 1:nevents
        mparts = [Particle[] for _ = 1:N_MODULES]
        mhits = [HitCoordinate[] for _ = 1:N_MODULES]
        for m = 1:N_MODULES
            ntr = clamp(poisson(rng, TRACK_PROFILE[m]), 0, MAX_TRACKS)
            for k = 1:ntr
                cursor = mod1(cursor + 1, length(pool))
                bp = pool[cursor].particle
                p = Particle(
                    bp.pid, Int32(k), Int32(m), bp.xCoord + 40.0 * randn(rng), bp.yCoord,
                    bp.pMag, bp.pMagTrue, bp.xDir, bp.yDir, bp.zDir,
                    bp.truePX, bp.truePY, bp.truePZ, bp.recoPX, bp.recoPY, bp.recoPZ,
                    bp.pathlength, bp.truePathlength, bp.t0, bp.rotMatrix,
                )
                push!(mparts[m], p)
                beta = particle_beta(p, get_particle_mass(abs(p.pid)))
                append!(mhits[m], project_pattern(p, beta, mapper, spectrum, context))
            end
            isempty(mparts[m]) && continue
            sig = mhits[m]
            if !isempty(sig)
                xmin, xmax = extrema(h.x for h in sig)
                ymin, ymax = extrema(h.y for h in sig)
                tmin, tmax = extrema(h.t for h in sig)
                nbkg = poisson(rng, BKG_HITS_PER_TRACK * length(mparts[m]))
                for _ = 1:nbkg
                    push!(sig, HitCoordinate(
                        xmin + (xmax - xmin) * rand(rng),
                        ymin + (ymax - ymin) * rand(rng),
                        tmin + (tmax - tmin) * rand(rng),
                    ))
                end
                if length(sig) > MAX_HITS
                    shuffle!(rng, sig)
                    resize!(sig, MAX_HITS)
                end
            end
        end
        events[ev] = ModuleEvent(mparts, mhits)
    end
    return events
end

@inline function load_event!(ws::ModulePIDWorkspace, ev::ModuleEvent)
    clear_modules!(ws)
    @inbounds for m = 1:N_MODULES
        set_module!(ws, m, ev.particles[m], ev.hits[m])
    end
    update_active_modules!(ws)
    return ws
end

# ---------------------------------------------------------------- Mode A: latency
function run_latency!(ws, events, pids, kw, nthreads, lat_ns::Vector{Float64})
    @inbounds for i in eachindex(events)
        t0 = time_ns()
        load_event!(ws, events[i])
        run_module_pid!(ws, pids; kw..., nthreads = nthreads)
        lat_ns[i] = time_ns() - t0
    end
    return nothing
end

# ---------------------------------------------------------------- Mode B: throughput
function run_throughput!(wss, events, pids, kw, nthreads)
    nev = length(events)
    next = Threads.Atomic{Int}(1)
    Threads.@threads for w = 1:nthreads
        ws = wss[w]
        while true
            i = Threads.atomic_add!(next, 1)
            i > nev && break
            load_event!(ws, events[i])
            run_module_pid!(ws, pids; kw..., nthreads = 1)
        end
    end
    return nothing
end

function main()
    args = parse_commandline()
    nevents = args["events"]; ntrials = args["trials"]
    max_threads = Threads.nthreads()
    thread_counts = parse_thread_list(args["thread-list"], max_threads)

    rng = Xoshiro(args["seed"])
    Random.seed!(args["seed"])
    dht = DetectorHitTester()
    spectrum = PhotonSpectrum(dht)
    mapper = PhotonMapper()
    context = PhotonContext(RADIATOR, CONSTANTS)
    pids = [211, 321, 2212]
    pts = [TorchPID.random_point(context.radiator.depth) for _ = 1:1000]
    kw = (spectrum = spectrum, mapper = mapper, background = FlatBackground(DETECTOR),
          cfg = ReverseMappingConfig(), context = context, random_points = pts)

    println("Generating $nevents synthetic 18-module events ...")
    events = generate_event_stream(nevents, rng, spectrum, mapper, context)

    trk = [sum(length, e.particles) for e in events]
    hts = [sum(length, e.hits) for e in events]
    mod_t = [length(e.particles[m]) for e in events for m = 1:N_MODULES]
    mod_h = [length(e.hits[m]) for e in events for m = 1:N_MODULES]
    act = [count(m -> !isempty(e.particles[m]) && !isempty(e.hits[m]), 1:N_MODULES) for e in events]
    work = [length(e.particles[m]) * length(e.hits[m]) for e in events for m = 1:N_MODULES]
    heavy_frac = mean(begin
        w = [length(e.particles[m]) * length(e.hits[m]) for m = 1:N_MODULES]
        maximum(w) / max(sum(w), 1)
    end for e in events)
    @printf("Tracks/event %.1f, hits/event %.0f, active modules %.1f, max module (%d tracks, %d hits), heaviest-module work share %.1f%%\n",
        mean(trk), mean(hts), mean(act), maximum(mod_t), maximum(mod_h), 100 * heavy_frac)

    # Correctness: multi-threaded dispatcher reproduces sequential results bit-for-bit-ish
    ws_ref = ModulePIDWorkspace(spectrum, Float64; max_tracks = 16, max_hits = 750)
    ws_mt = ModulePIDWorkspace(spectrum, Float64; max_tracks = 16, max_hits = 750)
    for i = 1:min(5, nevents)
        load_event!(ws_ref, events[i]); run_module_pid!(ws_ref, pids; kw..., nthreads = 1)
        load_event!(ws_mt, events[i]); run_module_pid!(ws_mt, pids; kw..., nthreads = max_threads)
        for s in ws_ref.active_slots
            m = s - 1
            @assert module_log_probs(ws_ref, m) ≈ module_log_probs(ws_mt, m) "Mismatch in event $i module $m"
        end
    end
    alloc_seq = @allocated (load_event!(ws_ref, events[1]); run_module_pid!(ws_ref, pids; kw..., nthreads = 1))
    println("Consistency check passed; sequential per-event allocations = $alloc_seq bytes")

    # ---- Mode A
    lat = Vector{Float64}(undef, nevents)
    latency_rows = NamedTuple[]
    t1 = 0.0
    println("\n=== A. Latency mode (sequential events, module-parallel) ===")
    @printf("%-8s | %10s | %10s | %10s | %10s | %8s | %10s\n", "Threads", "Mean (ms)", "Median", "p95 (ms)", "Events/s", "Speedup", "Efficiency")
    for nt in thread_counts
        run_latency!(ws_mt, events, pids, kw, nt, lat)        # warmup
        best = Inf; med = 0.0; p95 = 0.0
        for _ = 1:ntrials
            GC.gc()
            run_latency!(ws_mt, events, pids, kw, nt, lat)
            tot = sum(lat) * 1e-9
            if tot < best
                best = tot
                s = sort(lat)
                med = s[cld(nevents, 2)] * 1e-6
                p95 = s[clamp(ceil(Int, 0.95 * nevents), 1, nevents)] * 1e-6
            end
        end
        nt == 1 && (t1 = best)
        row = (nt = nt, mean_ms = 1e3 * best / nevents, med_ms = med, p95_ms = p95,
               ev_s = nevents / best, speedup = t1 / best, eff = 100 * (t1 / best) / nt)
        push!(latency_rows, row)
        @printf("%-8d | %10.2f | %10.2f | %10.2f | %10.1f | %7.2fx | %9.1f%%\n",
            nt, row.mean_ms, row.med_ms, row.p95_ms, row.ev_s, row.speedup, row.eff)
    end

    # ---- Mode B
    wss = [ModulePIDWorkspace(spectrum, Float64; max_tracks = 16, max_hits = 750) for _ = 1:max_threads]
    thr_rows = NamedTuple[]
    t1 = 0.0
    println("\n=== B. Throughput mode (concurrent events, modules sequential per event) ===")
    @printf("%-8s | %10s | %10s | %8s | %10s\n", "Threads", "Total (s)", "Events/s", "Speedup", "Efficiency")
    for nt in thread_counts
        run_throughput!(wss, events, pids, kw, nt)            # warmup
        best = Inf
        for _ = 1:ntrials
            GC.gc()
            best = min(best, @elapsed run_throughput!(wss, events, pids, kw, nt))
        end
        nt == 1 && (t1 = best)
        row = (nt = nt, total_s = best, ev_s = nevents / best, speedup = t1 / best, eff = 100 * (t1 / best) / nt)
        push!(thr_rows, row)
        @printf("%-8d | %10.3f | %10.1f | %7.2fx | %9.1f%%\n", nt, row.total_s, row.ev_s, row.speedup, row.eff)
    end

    # ---- Report
    out = args["output"]
    cpu = isempty(Sys.cpu_info()) ? "Unknown CPU" : Sys.cpu_info()[1].model
    open(out, "w") do io
        println(io, "# 18-Module TORCH Event-Stream Benchmark (CPU module dispatcher)\n")
        println(io, "## Configuration")
        println(io, "- **Date**: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- **Julia**: `$(VERSION)`, **CPU**: $(cpu) (`$(Sys.CPU_THREADS)` logical, `$(max_threads)` Julia threads)")
        println(io, "- **Events**: `$(nevents)` synthetic 18-module events, `$(ntrials)` trials (best reported), hypotheses `[211, 321, 2212]`, `n_gen = 1000`")
        @printf(io, "- **Occupancy**: %.1f tracks/event, %.0f hits/event, %.1f active modules/event; largest module %d tracks / %d hits\n",
            mean(trk), mean(hts), mean(act), maximum(mod_t), maximum(mod_h))
        @printf(io, "- **Load skew**: heaviest module carries on average **%.1f%%** of the event work (`ntracks × nhits`) → Amdahl bound for pure module parallelism ≈ **%.1fx**\n",
            100 * heavy_frac, 1 / heavy_frac)
        println(io, "- **Allocations**: `$(alloc_seq)` bytes per event (sequential dispatcher, pre-sized `ModulePIDWorkspace`)\n")
        println(io, "## A. Latency mode — sequential events, modules in parallel (`run_module_pid!`)\n")
        println(io, "| Threads | Mean latency (ms) | Median (ms) | p95 (ms) | Events/s | Speedup | Efficiency |")
        println(io, "|--:|--:|--:|--:|--:|--:|--:|")
        for r in latency_rows
            @printf(io, "| %d | %.2f | %.2f | %.2f | %.1f | **%.2fx** | %.1f%% |\n", r.nt, r.mean_ms, r.med_ms, r.p95_ms, r.ev_s, r.speedup, r.eff)
        end
        println(io, "\n## B. Throughput mode — concurrent events, modules sequential per event\n")
        println(io, "| Threads | Total (s) | Events/s | Speedup | Efficiency |")
        println(io, "|--:|--:|--:|--:|--:|")
        for r in thr_rows
            @printf(io, "| %d | %.3f | %.1f | **%.2fx** | %.1f%% |\n", r.nt, r.total_s, r.ev_s, r.speedup, r.eff)
        end
    end
    println("\nMarkdown report written to: $out")
end

main()
