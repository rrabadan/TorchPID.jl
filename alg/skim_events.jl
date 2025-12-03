#!/usr/bin/env julia

using Base.Threads
using UnROOT
using TorchPID
using DataFrames
using Arrow

"""
    skim_events(root_file::String, particles_output::String, hits_output::String; detector_kwargs...)

Loops over events in a ROOT file, filters particles and hits based on criteria,
and writes them to Parquet files.

For each event:
- Gets all particles
- Gets all hits that pass DetectorHitTester.test_photon AND have moduleId matching at least one Module from particles

# Arguments
- `root_file::String`: Path to the input ROOT file
- `particles_output::String`: Path prefix for particles Parquet file
- `hits_output::String`: Path prefix for filtered hits Parquet file
- `detector_kwargs...`: Keyword arguments for DetectorHitTester constructor

# Example
```julia
skim_events("input.root", "particles", "hits")
```
"""
function skim_events(root_file::String, particles_output::String, hits_output::String; detector_kwargs...)
    # Create detector tester
    detector = DetectorHitTester(; detector_kwargs...)

    # Create event reader
    reader = create_event_reader(root_file, useTruth=false, debugPhotons=true)

    # Prepare output data
    all_particles = Vector{Any}()
    all_filtered_hits = Vector{Any}()
    lock_obj = ReentrantLock()

    # Loop over events
    for event_data in event_iterator(reader)
        event = event_data.event
        particles = event_data.particles
        hits = event_data.hits

        # Get unique modules from particles
        particle_modules = Set(particle.moduleId for particle in particles)
        # particles_trackid = Set(particle.trackId for particle in particles)

        # Filter hits: must pass test_photon and have moduleId in particle_modules
        filtered_hits = []
        for hit in hits
            if test_photon(detector, hit.energy * 1.0e6) && (hit.moduleId in particle_modules)
                push!(filtered_hits, hit)
            end
        end

        # Collect particles and filtered hits (thread-safe if needed)
        lock(lock_obj) do
            append!(all_particles, particles)
            append!(all_filtered_hits, filtered_hits)
        end
    end

    println("Total particles: $(length(all_particles))")
    println("Total filtered hits: $(length(all_filtered_hits))")

    # Write particles to Parquet
    if !isempty(all_particles)
        particles_df = DataFrame(all_particles)
        Arrow.write(particles_output * ".parquet", particles_df)
        println("Particles written to $(particles_output).parquet")
    end

    # Write filtered hits to Parquet
    if !isempty(all_filtered_hits)
        hits_df = DataFrame(all_filtered_hits)
        Arrow.write(hits_output * ".parquet", hits_df)
        println("Filtered hits written to $(hits_output).parquet")
    end
end

# Example usage
if length(ARGS) >= 3
    root_file = ARGS[1]
    particles_output = ARGS[2]
    hits_output = ARGS[3]
    skim_events(root_file, particles_output, hits_output)
else
    println("Usage: julia alg/skim_events.jl <input.root> <particles_prefix> <hits_prefix>")
end