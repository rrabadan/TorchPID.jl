#!/usr/bin/env julia

using Base.Threads
using UnROOT
using TorchPID
using DataFrames
using Arrow

"""
    skim_hits(root_file::String, output_file::String; detector_config...)

Reads hits from a ROOT file, applies DetectorHitTester.test_photon to each hit,
and writes the passing hits to a Parquet file.

# Arguments
- `root_file::String`: Path to the input ROOT file
- `output_file::String`: Path to the output Parquet file for skimmed hits (without extension)
- `detector_config...`: Keyword arguments for DetectorHitTester constructor

# Example
```julia
skim_hits("input.root", "output")
```
"""
function skim_hits(root_file::String, output_file::String; detector_kwargs...)
    # Create detector tester
    detector = DetectorHitTester(; detector_kwargs...)

    # Open ROOT file
    f = ROOTFile(root_file)

    # Read hits tree
    hits_tree = LazyTree(f, "TorchTupleAlg/TorchHits")

    # Prepare output data with lock for thread safety
    passing_hits = Vector{Any}()
    lock_obj = ReentrantLock()

    # Process each hit
    @threads for hit in hits_tree
        if test_photon(detector, hit.Energy*1.0e6)  # Convert MeV to eV
            lock(lock_obj) do
                push!(passing_hits, hit)
            end
        end
    end

    println("Total hits: $(length(hits_tree))")
    println("Passing hits: $(length(passing_hits))")
    println("Passing fraction: $(length(passing_hits) / length(hits_tree))")

    # Convert to DataFrame and write to Parquet
    if !isempty(passing_hits)
        df = DataFrame(passing_hits)
        Arrow.write(output_file * ".parquet", df)
        println("Skimmed hits written to $(output_file).parquet")
    else
        println("No passing hits to write.")
    end

    close(f)
end

# Example usage
if length(ARGS) >= 2
    root_file = ARGS[1]
    output_file = ARGS[2]
    skim_hits(root_file, output_file)
else
    println("Usage: julia alg/skim_hits.jl <input.root> <output_prefix>")
end