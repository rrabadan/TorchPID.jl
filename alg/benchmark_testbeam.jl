using Dates
using Printf
using Random
using Statistics
using TorchPID

struct BenchmarkRow
    nthreads::Int
    median_ms::Float64
    min_ms::Float64
    events_per_sec::Float64
    mphotons_per_sec::Float64
    alloc_mib::Float64
    speedup::Float64
    efficiency_pct::Float64
end

function parse_commandline(args::Vector{String} = ARGS)
    opts = Dict{String,Any}(
        "events" => 50_000,
        "trials" => 5,
        "thread-list" => "1,2,4,8,16",
        "output" => "benchmark_testbeam_report.md",
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if (arg == "--events" || arg == "-n") && i < length(args)
            opts["events"] = parse(Int, args[i + 1])
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

function run_timed_trials(f::Function, ntrials::Int, nevents::Int, total_photons::Int)
    # Warmup
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
    mph_s = (total_photons * 1e-6) / (med_ms * 1e-3)
    alloc_mib = alloc_bytes / (1024.0^2)
    return (med_ms, mn_ms, ev_s, mph_s, alloc_mib)
end

function benchmark_suite(
    label::String,
    thread_counts::Vector{Int},
    ntrials::Int,
    nevents::Int,
    total_photons::Int,
    runner::Function,
)::Vector{BenchmarkRow}
    println("\n=== $label (N = $nevents events, $ntrials trials) ===")
    @printf(
        "%-8s | %11s | %11s | %14s | %13s | %10s | %8s | %10s\n",
        "Threads",
        "Median (ms)",
        "Min (ms)",
        "Events / sec",
        "MPhotons / s",
        "Alloc (MiB)",
        "Speedup",
        "Efficiency"
    )
    println(repeat("-", 102))

    rows = Vector{BenchmarkRow}(undef, length(thread_counts))
    t1_ms = 0.0

    for (idx, nt) in enumerate(thread_counts)
        med_ms, mn_ms, ev_s, mph_s, alloc_mib =
            run_timed_trials(() -> runner(nt), ntrials, nevents, total_photons)
        if idx == 1
            t1_ms = med_ms
        end
        speedup = t1_ms / med_ms
        eff = 100.0 * speedup / nt
        rows[idx] = BenchmarkRow(nt, med_ms, mn_ms, ev_s, mph_s, alloc_mib, speedup, eff)
        @printf(
            "%-8d | %11.2f | %11.2f | %14.0f | %13.2f | %10.2f | %7.2fx | %9.1f%%\n",
            nt,
            med_ms,
            mn_ms,
            ev_s,
            mph_s,
            alloc_mib,
            speedup,
            eff
        )
    end
    return rows
end

function format_markdown_table(rows::Vector{BenchmarkRow}; show_photons::Bool = true)::String
    io = IOBuffer()
    if show_photons
        println(
            io,
            "| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Photon Rate (Mphotons/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |",
        )
        println(
            io,
            "|--------:|-----------------:|--------------:|----------------------:|-------------------------:|-------------:|--------------:|--------------------:|",
        )
        for r in rows
            @printf(
                io,
                "| %d | %.2f | %.2f | %.0f | %.2f | %.2f | **%.2fx** | %.1f%% |\n",
                r.nthreads,
                r.median_ms,
                r.min_ms,
                r.events_per_sec,
                r.mphotons_per_sec,
                r.alloc_mib,
                r.speedup,
                r.efficiency_pct
            )
        end
    else
        println(
            io,
            "| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |",
        )
        println(
            io,
            "|--------:|-----------------:|--------------:|----------------------:|-------------:|--------------:|--------------------:|",
        )
        for r in rows
            @printf(
                io,
                "| %d | %.2f | %.2f | %.0f | %.2f | **%.2fx** | %.1f%% |\n",
                r.nthreads,
                r.median_ms,
                r.min_ms,
                r.events_per_sec,
                r.alloc_mib,
                r.speedup,
                r.efficiency_pct
            )
        end
    end
    return String(take!(io))
end

function write_report(
    output_path::String,
    nevents::Int,
    ntrials::Int,
    total_photons::Int,
    total_pixels_noshare::Int,
    total_pixels_share::Int,
    part_rows::Vector{BenchmarkRow},
    phot_noshare_rows::Vector{BenchmarkRow},
    phot_share_rows::Vector{BenchmarkRow},
    e2e_rows::Vector{BenchmarkRow},
)
    cpu_info = Sys.cpu_info()
    cpu_model = isempty(cpu_info) ? "Unknown CPU" : cpu_info[1].model
    mean_photons = total_photons / nevents
    mean_pixels_noshare = total_pixels_noshare / nevents
    mean_pixels_share = total_pixels_share / nevents

    open(output_path, "w") do io
        println(io, "# TestBeam Single-Thread vs. Multithreading Benchmark Report")
        println(io)
        println(io, "## Benchmark Configuration")
        println(io, "- **Date**: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "- **Julia Version**: `$(VERSION)`")
        println(io, "- **CPU**: $(cpu_model) (`$(Sys.CPU_THREADS)` logical threads, `$(Threads.nthreads())` active Julia threads)")
        println(io, "- **Events per Batch**: `$(nevents)` (`$(ntrials)` trials per thread count, `GC.gc()` before each trial)")
        @printf(
            io,
            "- **Mean Cherenkov Photons Traced per Event**: `%.2f` (`%d` total photons per batch)\n",
            mean_photons,
            total_photons
        )
        @printf(
            io,
            "- **Mean Digitized Pixel Hits per Event**: `%.2f` (no sharing) / `%.2f` (charge sharing)\n",
            mean_pixels_noshare,
            mean_pixels_share
        )
        println(io)
        println(io, "---")
        println(io)
        println(io, "## 1. End-to-End TestBeam Simulation (`generate_particles` + `generate_photons_batch`)")
        println(io)
        print(io, format_markdown_table(e2e_rows; show_photons = true))
        println(io)
        println(io, "## 2. Photon Ray-Tracing & FrontEnd Digitization (`generate_photons_batch`, `sharing_mode = false`)")
        println(io)
        print(io, format_markdown_table(phot_noshare_rows; show_photons = true))
        println(io)
        println(io, "## 3. Photon Ray-Tracing & Charge-Sharing Digitization (`generate_photons_batch`, `sharing_mode = true`)")
        println(io)
        print(io, format_markdown_table(phot_share_rows; show_photons = true))
        println(io)
        println(io, "## 4. Beam Particle Generation (`generate_particles`)")
        println(io)
        print(io, format_markdown_table(part_rows; show_photons = false))
    end
    println("\nMarkdown benchmark report written to: $output_path")
end

function main()
    args = parse_commandline()
    nevents = args["events"]
    ntrials = args["trials"]
    output_path = args["output"]
    max_threads = Threads.nthreads()
    thread_counts = parse_thread_list(args["thread-list"], max_threads)

    Random.seed!(42)
    tb = TestBeamSimulator()
    spectrum = PhotonSpectrum()
    mapper = PhotonMapper()
    context = PhotonContext(RADIATOR, CONSTANTS)
    fe_noshare = FrontEnd(DETECTOR; sharing_mode = false, smear_time = false)
    cdt_noshare = ChargeDepositTester()

    sig_share = TorchPID.SignalParameters(gain = 1.5e6)
    fe_share = FrontEnd(DETECTOR; sharing_mode = true, smear_time = false)
    cdt_share = ChargeDepositTester(DETECTOR, sig_share; charge_sharing = true)

    particles = generate_particles(tb, context.radiator, context.constants, nevents; nthreads = 1)
    sample_noshare = generate_photons_batch(
        particles, spectrum, mapper, context, fe_noshare, cdt_noshare; nthreads = 1
    )
    sample_share = generate_photons_batch(
        particles, spectrum, mapper, context, fe_share, cdt_share; nthreads = 1
    )

    total_photons = sum(b === nothing ? 0 : b.photon_yield for b in sample_noshare)
    total_pixels_noshare = sum(b === nothing ? 0 : b.npixels for b in sample_noshare)
    total_pixels_share = sum(b === nothing ? 0 : b.npixels for b in sample_share)

    part_rows = benchmark_suite(
        "1. Particle Generation (generate_particles)",
        thread_counts,
        ntrials,
        nevents,
        total_photons,
        nt -> generate_particles(tb, context.radiator, context.constants, nevents; nthreads = nt),
    )

    phot_noshare_rows = benchmark_suite(
        "2. Photon Generation & Digitization (sharing_mode = false)",
        thread_counts,
        ntrials,
        nevents,
        total_photons,
        nt -> generate_photons_batch(
            particles, spectrum, mapper, context, fe_noshare, cdt_noshare; nthreads = nt
        ),
    )

    phot_share_rows = benchmark_suite(
        "3. Photon Generation & Digitization (sharing_mode = true)",
        thread_counts,
        ntrials,
        nevents,
        total_photons,
        nt -> generate_photons_batch(
            particles, spectrum, mapper, context, fe_share, cdt_share; nthreads = nt
        ),
    )

    e2e_rows = benchmark_suite(
        "4. End-to-End TestBeam Pipeline (generate_particles + generate_photons_batch)",
        thread_counts,
        ntrials,
        nevents,
        total_photons,
        nt -> begin
            p = generate_particles(tb, context.radiator, context.constants, nevents; nthreads = nt)
            generate_photons_batch(
                p, spectrum, mapper, context, fe_noshare, cdt_noshare; nthreads = nt
            )
        end,
    )

    write_report(
        output_path,
        nevents,
        ntrials,
        total_photons,
        total_pixels_noshare,
        total_pixels_share,
        part_rows,
        phot_noshare_rows,
        phot_share_rows,
        e2e_rows,
    )
end

main()
