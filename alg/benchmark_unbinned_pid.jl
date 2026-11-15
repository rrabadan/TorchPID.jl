using Dates
using Printf
using Random
using Statistics
using TorchPID

struct PIDBenchmarkRow
    nthreads::Int
    median_ms::Float64
    min_ms::Float64
    events_per_sec::Float64
    tracks_per_sec::Float64
    alloc_mib::Float64
    speedup::Float64
    efficiency_pct::Float64
end

function parse_commandline(args::Vector{String} = ARGS)
    opts = Dict{String,Any}(
        "events" => 5_000,
        "pixel-events" => 1_000,
        "tracks-per-event" => 5,
        "trials" => 3,
        "thread-list" => "1,2,4,8,16",
        "output" => "benchmark_unbinned_pid_report.md",
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if (arg == "--events" || arg == "-n") && i < length(args)
            opts["events"] = parse(Int, args[i + 1])
            i += 2
        elseif arg == "--pixel-events" && i < length(args)
            opts["pixel-events"] = parse(Int, args[i + 1])
            i += 2
        elseif arg == "--tracks-per-event" && i < length(args)
            opts["tracks-per-event"] = parse(Int, args[i + 1])
            i += 2
        elseif (arg == "--trials" || arg == "-t") && i < length(args)
            opts["trials"] = parse(Int, args[i + 1])
            i += 2
        elseif arg == "--thread-list" && i < length(args)
            opts["thread-list"] = args[i + 1]
            i += 2
        elseif (arg == "--output" || arg == "-o") && i < length(args)
            opts["output"] = args[i + 1]
            i += 2
        else
            error("Unknown or incomplete argument: $arg")
        end
    end
    return opts
end

function parse_thread_list(spec::String, max_threads::Int)::Vector{Int}
    requested = [parse(Int, strip(s)) for s in split(spec, ",") if !isempty(strip(s))]
    filtered = unique(filter(t -> 1 <= t <= max_threads, requested))
    if isempty(filtered) || first(filtered) != 1
        pushfirst!(filtered, 1)
    end
    if !(max_threads in filtered)
        push!(filtered, max_threads)
    end
    return sort!(unique!(filtered))
end

function run_timed_trials(f::Function, ntrials::Int, nevents::Int, ntracks_total::Int)
    f()
    times_ms = Vector{Float64}(undef, ntrials)
    alloc_bytes = 0
    for t = 1:ntrials
        GC.gc()
        stats = @timed f()
        times_ms[t] = stats.time * 1e3
        alloc_bytes = stats.bytes
    end
    med_ms = median(times_ms)
    mn_ms = minimum(times_ms)
    ev_s = nevents / (med_ms * 1e-3)
    tr_s = ntracks_total / (med_ms * 1e-3)
    alloc_mib = alloc_bytes / (1024.0^2)
    return (med_ms, mn_ms, ev_s, tr_s, alloc_mib)
end

function benchmark_suite(
    label::String,
    thread_counts::Vector{Int},
    ntrials::Int,
    nevents::Int,
    ntracks_total::Int,
    runner::Function,
)::Vector{PIDBenchmarkRow}
    println("\n=== $label (N = $nevents events, $ntracks_total tracks, $ntrials trials) ===")
    @printf(
        "%-8s | %11s | %11s | %14s | %13s | %10s | %8s | %10s\n",
        "Threads",
        "Median (ms)",
        "Min (ms)",
        "Events / sec",
        "Tracks / sec",
        "Alloc (MiB)",
        "Speedup",
        "Efficiency"
    )
    println(repeat("-", 102))

    rows = Vector{PIDBenchmarkRow}(undef, length(thread_counts))
    t1_ms = 0.0

    for (idx, nt) in enumerate(thread_counts)
        med_ms, mn_ms, ev_s, tr_s, alloc_mib =
            run_timed_trials(() -> runner(nt), ntrials, nevents, ntracks_total)
        if idx == 1
            t1_ms = med_ms
        end
        speedup = t1_ms / med_ms
        eff = 100.0 * speedup / nt
        rows[idx] = PIDBenchmarkRow(nt, med_ms, mn_ms, ev_s, tr_s, alloc_mib, speedup, eff)
        @printf(
            "%-8d | %11.2f | %11.2f | %14.1f | %13.1f | %10.2f | %7.2fx | %9.1f%%\n",
            nt,
            med_ms,
            mn_ms,
            ev_s,
            tr_s,
            alloc_mib,
            speedup,
            eff
        )
    end
    return rows
end

function format_markdown_table(rows::Vector{PIDBenchmarkRow})::String
    io = IOBuffer()
    println(
        io,
        "| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Track Rate (tracks/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |",
    )
    println(
        io,
        "|--------:|-----------------:|--------------:|----------------------:|----------------------:|-------------:|--------------:|--------------------:|",
    )
    for r in rows
        @printf(
            io,
            "| %d | %.2f | %.2f | %.1f | %.1f | %.2f | **%.2fx** | %.1f%% |\n",
            r.nthreads,
            r.median_ms,
            r.min_ms,
            r.events_per_sec,
            r.tracks_per_sec,
            r.alloc_mib,
            r.speedup,
            r.efficiency_pct
        )
    end
    return String(take!(io))
end

function write_report(
    output_path::String,
    nevents::Int,
    npixel_events::Int,
    ntracks_per_ev::Int,
    ntrials::Int,
    mean_hits::Float64,
    alloc_no_ws::Int,
    alloc_ws::Int,
    hc_rows::Vector{PIDBenchmarkRow},
    px_rows::Vector{PIDBenchmarkRow},
)
    cpu_info = Sys.cpu_info()
    cpu_model = isempty(cpu_info) ? "Unknown CPU" : cpu_info[1].model

    open(output_path, "w") do io
        println(io, "# Unbinned PID Reconstruction (`run_unbinned_pid_batch!`) Benchmark Report")
        println(io)
        println(io, "## Benchmark Configuration")
        println(io, "- **Date**: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- **Julia Version**: `$(VERSION)`")
        println(io, "- **CPU**: $(cpu_model) (`$(Sys.CPU_THREADS)` logical threads, `$(Threads.nthreads())` active Julia threads)")
        println(io, "- **Tracks per Event**: `$(ntracks_per_ev)` (`3` mass hypotheses `[211, 321, 2212]`, `n_gen = 1000` MC efficiency photons per track-hypothesis)")
        @printf(io, "- **Mean Digitized Hits per Event**: `%.1f`\n", mean_hits)
        @printf(
            io,
            "- **Per-Event Heap Allocations**: `%d bytes` (with `UnbinnedPIDWorkspace`) vs. `%d bytes` (unbuffered)\n",
            alloc_ws,
            alloc_no_ws
        )
        println(io, "- **Trials per Thread Count**: `$(ntrials)` (`GC.gc()` before each trial)")
        println(io)
        println(io, "---")
        println(io)
        println(io, "## 1. Multi-Track Unbinned PID — `HitCoordinate` (`$(nevents)` events / `$(nevents * ntracks_per_ev)` tracks)")
        println(io)
        print(io, format_markdown_table(hc_rows))
        println(io)
        println(io, "## 2. Multi-Track Unbinned PID — `PixelHit` Fused 9-Point Quadrature (`$(npixel_events)` events / `$(npixel_events * ntracks_per_ev)` tracks)")
        println(io)
        print(io, format_markdown_table(px_rows))
    end
    println("\nMarkdown benchmark report written to: $output_path")
end

function main()
    args = parse_commandline()
    nevents = args["events"]
    npixel_events = min(args["pixel-events"], nevents)
    ntracks_per_ev = args["tracks-per-event"]
    ntrials = args["trials"]
    output_path = args["output"]
    max_threads = Threads.nthreads()
    thread_counts = parse_thread_list(args["thread-list"], max_threads)

    Random.seed!(42)
    tb = TestBeamSimulator()
    spectrum = PhotonSpectrum()
    mapper = PhotonMapper()
    context = PhotonContext(RADIATOR, CONSTANTS)
    fe = FrontEnd(DETECTOR)
    cdt = ChargeDepositTester()
    bkg = FlatBackground(DETECTOR)
    hypos = [211, 321, 2212]
    pts = [TorchPID.random_point(context.radiator.depth) for _ = 1:1000]

    n_raw = nevents * ntracks_per_ev
    tb_parts = generate_particles(tb, context.radiator, context.constants, n_raw; nthreads = max_threads)
    tb_phots = generate_photons_batch(tb_parts, spectrum, mapper, context, fe, cdt; nthreads = max_threads)

    ev_particles = Vector{Vector{Particle}}(undef, nevents)
    ev_hitcoords = Vector{Vector{HitCoordinate}}(undef, nevents)
    ev_pixelhits = Vector{Vector{PixelHit}}(undef, nevents)

    for ev = 1:nevents
        ps = Vector{Particle}(undef, ntracks_per_ev)
        hcs = HitCoordinate[]
        pxs = PixelHit[]
        for k = 1:ntracks_per_ev
            idx = (ev - 1) * ntracks_per_ev + k
            ps[k] = tb_parts[idx].particle
            ph = tb_phots[idx]
            if ph !== nothing
                for j = 1:ph.npixels
                    px = PixelHit(Int(ph.xpixels[j]), Int(ph.ypixels[j]), Int(ph.tpixels[j]), ph.charge[j])
                    push!(pxs, px)
                    push!(hcs, get_hit_coordinate(px, DETECTOR, context.radiator))
                end
            end
        end
        ev_particles[ev] = ps
        ev_hitcoords[ev] = hcs
        ev_pixelhits[ev] = pxs
    end

    mean_hits = mean(length.(ev_hitcoords))
    lp_single = zeros(Float64, ntracks_per_ev, length(hypos))
    ws = UnbinnedPIDWorkspace(spectrum, Float64; max_tracks = ntracks_per_ev + 2, max_hits = 300, max_pids = length(hypos), n_gen = 1000)

    # Warmup and check single-event allocations
    run_unbinned_pid!(lp_single, ev_particles[1], ev_hitcoords[1]; hypotheses = hypos, spectrum = spectrum, mapper = mapper, background = bkg, random_points = pts, threaded = false)
    run_unbinned_pid!(lp_single, ev_particles[1], ev_hitcoords[1]; hypotheses = hypos, spectrum = spectrum, mapper = mapper, background = bkg, random_points = pts, threaded = false, workspace = ws)
    alloc_no_ws = @allocated run_unbinned_pid!(lp_single, ev_particles[1], ev_hitcoords[1]; hypotheses = hypos, spectrum = spectrum, mapper = mapper, background = bkg, random_points = pts, threaded = false)
    alloc_ws = @allocated run_unbinned_pid!(lp_single, ev_particles[1], ev_hitcoords[1]; hypotheses = hypos, spectrum = spectrum, mapper = mapper, background = bkg, random_points = pts, threaded = false, workspace = ws)

    log_probs_hc = [zeros(Float64, ntracks_per_ev, length(hypos)) for _ = 1:nevents]
    log_probs_px = [zeros(Float64, ntracks_per_ev, length(hypos)) for _ = 1:npixel_events]
    ev_particles_px = @view(ev_particles[1:npixel_events])
    ev_pixelhits_px = @view(ev_pixelhits[1:npixel_events])

    hc_rows = benchmark_suite(
        "1. Multi-Track Unbinned PID — HitCoordinate (run_unbinned_pid_batch!)",
        thread_counts,
        ntrials,
        nevents,
        nevents * ntracks_per_ev,
        nt -> run_unbinned_pid_batch!(
            log_probs_hc,
            ev_particles,
            ev_hitcoords;
            hypotheses = hypos,
            spectrum = spectrum,
            mapper = mapper,
            background = bkg,
            random_points = pts,
            nthreads = nt,
        ),
    )

    px_rows = benchmark_suite(
        "2. Multi-Track Unbinned PID — PixelHit 9-Point Quadrature (run_unbinned_pid_batch!)",
        thread_counts,
        ntrials,
        npixel_events,
        npixel_events * ntracks_per_ev,
        nt -> run_unbinned_pid_batch!(
            log_probs_px,
            ev_particles_px,
            ev_pixelhits_px;
            hypotheses = hypos,
            spectrum = spectrum,
            mapper = mapper,
            background = bkg,
            random_points = pts,
            nthreads = nt,
        ),
    )

    write_report(
        output_path,
        nevents,
        npixel_events,
        ntracks_per_ev,
        ntrials,
        mean_hits,
        alloc_no_ws,
        alloc_ws,
        hc_rows,
        px_rows,
    )
end

main()
