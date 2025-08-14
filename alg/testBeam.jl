using ArgParse
using BenchmarkTools
using Plots
using Random

using TorchPID


function parse_commandline()
    s = ArgParseSettings()
    s.description = "Test Beam Simulation"
    @add_arg_table s begin
        "events"
        help = "number of events to simulate"
        required = true
        arg_type = Int
    end
    @add_arg_table s begin
        "--threads"
        help = "number of threads to use (default: auto)"
        arg_type = Int
        default = Sys.CPU_THREADS
    end
    @add_arg_table s begin
        "--plot-interval"
        help = "Update plots every N events"
        arg_type = Int
        default = 100
    end
    return parse_args(s)
end

function update_histograms(pMag_data, pixels_x, pixels_y)
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
    plot_interval = args["plot-interval"]
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

    update_histograms(pMag, pixels_x, pixels_y)

    for batch_start in 1:batch_size:events
        this_batch = min(batch_size, events - batch_start + 1)
        tb_particles = generate_particles(
            tb_simulator,
            photon_context.radiator,
            photon_context.constants,
            this_batch,
        )
        tb_photons = generate_photons_batch(
            tb_particles,
            photon_spectrum,
            photon_mapper,
            photon_context,
            frontend,
            charge_tester,
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
        if progress % plot_interval == 0
            update_histograms(pMag, pixels_x, pixels_y)
        end
    end

    update_histograms(pMag, pixels_x, pixels_y)
    println("\nSimulation complete. Press Enter to exit.")
    readline()
end

main()
