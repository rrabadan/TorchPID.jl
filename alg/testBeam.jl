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
            bins = 50,
            title = "Particle Momentum Magnitude",
            xlabel = "pMag",
            ylabel = "Count",
            legend = false,
        )
        p2 = scatter(
            pixels_x,
            pixels_y,
            title = "Photon Hits",
            xlabel = "X pixel",
            ylabel = "Y pixel",
            markersize = 2,
            markerstrokewidth = 0,
            alpha = 0.6,
            legend = false,
        )

        # Display the plots side by side with proper layout management
        plot(
            p1,
            p2,
            layout = (1, 2),
            size = (1200, 500),
            margin = 5Plots.mm,     # Add margins to avoid cutting off labels
            dpi = 300,              # Higher resolution
            titlefontsize = 10,     # Adjust title size
            guidefontsize = 8,      # Adjust axis label size
            tickfontsize = 6,
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

    println("Running in single-threaded mode")

    # Create configuration
    tb_simulator = TestBeamSimulator()
    photon_spectrum = PhotonSpectrum()
    photon_mapper = PhotonMapper()
    photon_factory = PhotonFactory(RADIATOR, CONSTANTS)
    frontend = FrontEnd(DETECTOR)
    charge_tester = ChargeDepositTester()

    # Initialize progress counter
    progress = 0

    # Data containers for histograms
    pMag = Vector{Float64}()
    pixels_x = Vector{Float64}()
    pixels_y = Vector{Float64}()

    # Create initial empty plot
    update_histograms(pMag, pixels_x, pixels_y)

    # Process events sequentially
    for event = 1:events

        tb_particle = generate_particle(
            tb_simulator,
            photon_factory.radiator,
            photon_factory.constants,
        )
        tb_photons = generate_photons(
            tb_particle,
            photon_spectrum,
            photon_mapper,
            photon_factory,
            frontend,
            charge_tester,
        )

        # Collect histogram data
        push!(pMag, tb_particle.particle.pMag)
        for i = 1:tb_photons.npixels
            push!(pixels_x, tb_photons.xpixels[i])
            push!(pixels_y, tb_photons.ypixels[i])
        end

        # Progress update
        progress += 1
        if progress % 100000 == 0
            println("Processed $progress events")
        end

        # Update plots periodically
        if progress % plot_interval == 0
            update_histograms(pMag, pixels_x, pixels_y)
        end
    end

    # Final histogram update
    update_histograms(pMag, pixels_x, pixels_y)
    println("\nSimulation complete. Press Enter to exit.")
    readline()
end

main()
