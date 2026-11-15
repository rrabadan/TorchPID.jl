using Plots
using Random

using TorchPID


function parse_commandline(raw_args::Vector{String} = ARGS)
    opts = Dict{String,Any}(
        "events" => nothing,
        "threads" => Threads.nthreads(),
        "plot-interval" => 1000,
        "no-plot" => false,
    )
    i = 1
    while i <= length(raw_args)
        arg = raw_args[i]
        if arg == "--threads" && i < length(raw_args)
            opts["threads"] = parse(Int, raw_args[i + 1])
            i += 2
        elseif arg == "--plot-interval" && i < length(raw_args)
            opts["plot-interval"] = parse(Int, raw_args[i + 1])
            i += 2
        elseif arg == "--no-plot"
            opts["no-plot"] = true
            i += 1
        elseif !startswith(arg, "-") && opts["events"] === nothing
            opts["events"] = parse(Int, arg)
            i += 1
        else
            error(
                "Usage: julia --project=. --threads=auto alg/testBeam.jl <events> [--threads N] [--plot-interval N] [--no-plot]",
            )
        end
    end
    if opts["events"] === nothing
        error(
            "Missing required positional argument <events>.\nUsage: julia --project=. --threads=auto alg/testBeam.jl <events> [--threads N] [--plot-interval N] [--no-plot]",
        )
    end
    return opts
end

function update_histograms(pMag_data, pixels_x, pixels_y)
    isempty(pMag_data) && return nothing
    try
        # Create histograms
        p1 = histogram(
            pMag_data,
            bins=50,
            title="Particle Momentum Magnitude",
            xlabel="pMag",
            ylabel="Count",
            legend=false,
        )
        p2 = scatter(
            pixels_x,
            pixels_y,
            title="Photon Hits",
            xlabel="X pixel",
            ylabel="Y pixel",
            markersize=2,
            markerstrokewidth=0,
            alpha=0.6,
            legend=false,
        )

        # Display the plots side by side with proper layout management
        plot(
            p1,
            p2,
            layout=(1, 2),
            size=(1200, 500),
            margin=5Plots.mm,     # Add margins to avoid cutting off labels
            dpi=300,              # Higher resolution
            titlefontsize=10,     # Adjust title size
            guidefontsize=8,      # Adjust axis label size
            tickfontsize=6,
        )       # Adjust tick label size

        display(current())  # Use display instead of gui() for better stability
    catch e
        println("Error in update_histograms: ", e)
    end
end

function main()
    args = parse_commandline()
    events = args["events"]
    nthreads = args["threads"]
    plot_interval = args["plot-interval"]
    no_plot = args["no-plot"]
    batch_size = 1000  # Choose a reasonable batch size

    tb_simulator = TestBeamSimulator()
    photon_spectrum = PhotonSpectrum()
    photon_mapper = PhotonMapper()
    photon_context = PhotonContext(RADIATOR, CONSTANTS)
    frontend = FrontEnd(DETECTOR)
    charge_tester = ChargeDepositTester()

    progress = 0
    pMag = Float64[]
    pixels_x = Float64[]
    pixels_y = Float64[]

    for batch_start in 1:batch_size:events
        this_batch = min(batch_size, events - batch_start + 1)
        tb_particles = generate_particles(
            tb_simulator,
            photon_context.radiator,
            photon_context.constants,
            this_batch;
            nthreads = nthreads,
        )
        tb_photons = generate_photons_batch(
            tb_particles,
            photon_spectrum,
            photon_mapper,
            photon_context,
            frontend,
            charge_tester;
            nthreads = nthreads,
        )

        for i in 1:this_batch
            push!(pMag, tb_particles[i].particle.pMag)
            photons = tb_photons[i]
            if photons !== nothing
                for j in 1:photons.npixels
                    push!(pixels_x, photons.xpixels[j])
                    push!(pixels_y, photons.ypixels[j])
                end
            end
        end

        progress += this_batch
        if !no_plot && (progress % plot_interval == 0)
            update_histograms(pMag, pixels_x, pixels_y)
        end
    end

    if !no_plot
        update_histograms(pMag, pixels_x, pixels_y)
    end
    println(
        "\nSimulation complete: $(length(pMag)) particles, $(length(pixels_x)) digitized pixel hits.",
    )
    if !no_plot && stdin isa Base.TTY
        println("Press Enter to exit.")
        readline()
    end
end

main()
