"""
    ObjPool{T}

A simple object pool to reduce allocations by reusing objects.

# Fields
- `data::Vector{T}`: Storage for pooled objects
- `size::Ref{Int}`: Reference to the pool size
"""
struct ObjPool{T}
    data::Vector{T}
    size::Ref{Int}

    function ObjPool{T}(size::Int) where {T}
        return new(Vector{T}(undef, size), Ref(size))
    end
end

"""
    ThreadSafeObjPool{T}

Thread-safe version of ObjPool that maintains a separate pool for each thread.

# Fields
- `pools::Vector{ObjPool{T}}`: Vector of per-thread object pools
"""
struct ThreadSafeObjPool{T}
    pools::Vector{ObjPool{T}}

    function ThreadSafeObjPool{T}(per_thread_size::Int) where {T}
        pools = [ObjPool{T}(per_thread_size) for _ = 1:Threads.nthreads()]
        return new(pools)
    end
end

"""
    get_thread_pool(pool::ThreadSafeObjPool{T}) where {T}

Get the object pool assigned to the current thread.

# Arguments
- `pool::ThreadSafeObjPool{T}`: The thread-safe object pool

# Returns
- The object pool for the current thread
"""
function get_thread_pool(pool::ThreadSafeObjPool{T}) where {T}
    return pool.pools[Threads.threadid()]
end

"""
    hits_branches(; debug=false)

Get the list of branches needed for photon hits.

# Arguments
- `debug::Bool=false`: Whether to include additional debug branches

# Returns
- Vector{String}: List of branch names to read
"""
function hits_branches(; debug = false)
    br = [
        "LocalX",
        "LocalY",
        "LocalZ",
        "ArrivalTime",
        "Energy",
        "EmissionX",
        "EmissionY",
        "EmissionZ",
        "EmissionTime",
        "Module",
        "Track",
        "IsDownward",
        "IsSurfaceScattered",
        "IsRayleighScattered",
    ]
    return debug ? vcat(br, "CherenkovTheta", "CherenkovPhi") : br
end

"""
    tracks_branches(; useTruth=false)

Get the list of branches needed for tracks.

# Arguments
- `useTruth::Bool=false`: Whether to use truth information branches

# Returns
- Vector{String}: List of branch names to read
"""
function tracks_branches(; useTruth = false)
    br = [
        "True_ID",
        "Track",
        "True_Momentum",
        "True_PX",
        "True_PY",
        "True_PZ",
        "Reco_Momentum",
        "Reco_PX",
        "Reco_PY",
        "Reco_PZ",
        "Reco_Pathlength",
        "True_Pathlength",
        "True_Origin_T",
    ]
    if useTruth
        append!(
            br,
            [
                "True_Module",
                "True_Local_X",
                "True_Local_Y",
                "True_Local_Dir_X",
                "True_Local_Dir_Y",
                "True_Local_Dir_Z",
            ],
        )
    else
        append!(
            br,
            [
                "Reco_Module",
                "Reco_Local_X",
                "Reco_Local_Y",
                "Reco_Local_Dir_X",
                "Reco_Local_Dir_Y",
                "Reco_Local_Dir_Z",
            ],
        )
    end
    return br
end

"""
    EventReader

Main reader for TORCH events from ROOT files.

# Fields
- `file::ROOTFile`: The underlying ROOT file
- `events::LazyTree`: Tree containing event data
- `tracks::LazyTree`: Tree containing track data
- `hits::LazyTree`: Tree containing photon hit data
- `run_numbers::Vector{Int}`: Cache of run numbers
- `event_numbers::Vector{Int}`: Cache of event numbers
- `track_counts::Vector{Int}`: Cache of track counts per event
- `hit_counts::Vector{Int}`: Cache of hit counts per event
- `track_indices::Vector{Int}`: Pre-calculated indices for tracks
- `hit_indices::Vector{Int}`: Pre-calculated indices for hits
- `useTruth::Bool`: Whether to use truth information
- `debugPhotons::Bool`: Whether debug photon info is enabled
- `trackquality::Int`: Track quality threshold
- `hitPool::ObjPool{PhotonHit}`: Object pool for photon hits
"""
struct EventReader
    file::ROOTFile
    events::LazyTree
    tracks::LazyTree
    hits::LazyTree

    # Cache frequently used data
    run_numbers::Vector{Int}
    event_numbers::Vector{Int}
    track_counts::Vector{Int}
    hit_counts::Vector{Int}

    # Pre-calculated indices
    track_indices::Vector{Int}
    hit_indices::Vector{Int}

    useTruth::Bool
    debugPhotons::Bool
    trackquality::Int

    hitPool::ObjPool{PhotonHit}
end

"""
    EventReader(file::String, useTruth::Bool, debugPhotons::Bool, trackQuality::Int)::EventReader

Internal constructor for EventReader that initializes all caches and pre-calculated indices.

# Arguments
- `file::String`: Path to the ROOT file
- `useTruth::Bool`: Whether to use truth information
- `debugPhotons::Bool`: Enable photon debugging
- `trackQuality::Int`: Track quality threshold

# Returns
- EventReader: Fully initialized event reader
"""
function EventReader(
    file::String,
    useTruth::Bool,
    debugPhotons::Bool,
    trackQuality::Int,
)::EventReader
    f = ROOTFile(file)

    events = LazyTree(f, "TorchTupleAlg/TorchEvent")
    tracks = LazyTree(f, "TorchTupleAlg/TorchTracks", tracks_branches(useTruth = useTruth))
    hits = LazyTree(f, "TorchTupleAlg/TorchHits", hits_branches(debug = debugPhotons))

    # Pre-calculate the number of events
    n_events = length(events)

    # Pre-allocate all caches
    run_numbers = Vector{Int64}(undef, n_events)
    event_numbers = Vector{Int64}(undef, n_events)
    track_counts = Vector{Int64}(undef, n_events)
    hit_counts = Vector{Int64}(undef, n_events)

    track_indices = ones(Int64, n_events + 1)
    hit_indices = ones(Int64, n_events + 1)

    # Initialize the first indices

    if n_events > 1000
        Threads.@threads for i = 1:n_events
            event = events[i]
            run_numbers[i] = event.Run
            event_numbers[i] = event.Event
            track_counts[i] = event.Tracks
            hit_counts[i] = event.Photons
        end
    else
        for i = 1:n_events
            event = events[i]
            run_numbers[i] = event.Run
            event_numbers[i] = event.Event
            track_counts[i] = event.Tracks
            hit_counts[i] = event.Photons
        end
    end

    track_indices[2:end] = cumsum(track_counts) .+ 1
    hit_indices[2:end] = cumsum(hit_counts) .+ 1

    HitPool = ObjPool{PhotonHit}(maximum(hit_counts))
    for i = 1:maximum(hit_counts)
        HitPool.data[i] = PhotonHit(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            0,
            false,
            false,
            false,
        )
    end

    return EventReader(
        f,
        events,
        tracks,
        hits,
        run_numbers,
        event_numbers,
        track_counts,
        hit_counts,
        track_indices,
        hit_indices,
        useTruth,
        debugPhotons,
        trackQuality,
        HitPool,
    )
end

"""
    create_event_reader(file::String; useTruth::Bool=false, debugPhotons::Bool=false, trackQuality::Int=0)

Public constructor for EventReader with named parameters.

# Arguments
- `file::String`: Path to the ROOT file
- `useTruth::Bool=false`: Whether to use truth information
- `debugPhotons::Bool=false`: Enable photon debugging
- `trackQuality::Int=0`: Track quality threshold

# Returns
- EventReader: Fully initialized event reader
"""
function create_event_reader(
    file::String;
    useTruth::Bool = false,
    debugPhotons::Bool = false,
    trackQuality::Int = 0,
)
    # Use type parameter to select proper PhotonHit variant
    return EventReader(file, useTruth, debugPhotons, trackQuality)
end

"""
    get_particle(data; useTruth::Bool)::Particle

Converts a track entry from TorchTracks ROOT TTree to a Particle object.

# Arguments
- `data`: LazyEvent containing track information
- `useTruth::Bool`: Whether to use truth information

# Returns
- Particle: A new Particle object
"""
function get_particle(data; useTruth::Bool)::Particle
    return useTruth ? get_particle_true(data) : get_particle_reco(data)
end

"""
    get_particle_true(data)::Particle

Creates a Particle object from truth information.

# Arguments
- `data`: LazyEvent containing track information

# Returns
- Particle: A new Particle object with truth information
"""
@inline function get_particle_true(data)::Particle
    return Particle(
        data.True_ID,
        data.Track,
        data.True_Module,
        data.True_Local_X,
        data.True_Local_Y,
        data.True_Momentum,
        data.True_Momentum,
        data.True_Local_Dir_X,
        data.True_Local_Dir_Y,
        data.True_Local_Dir_Z,
        data.True_PX,
        data.True_PY,
        data.True_PZ,
        data.True_PX,
        data.True_PY,
        data.True_PZ,
        data.True_Pathlength,
        data.True_Pathlength,
        data.True_Origin_T,
        [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    )
end

"""
    get_particle_reco(data)::Particle

Creates a Particle object from reconstructed information.

# Arguments
- `data`: LazyEvent containing track information

# Returns
- Particle: A new Particle object with reconstructed information
"""
@inline function get_particle_reco(data)::Particle
    return Particle(
        data.True_ID,
        data.Track,
        data.Reco_Module,
        data.Reco_Local_X,
        data.Reco_Local_Y,
        data.Reco_Momentum,
        data.True_Momentum, # Still using true momentum
        data.Reco_Local_Dir_X,
        data.Reco_Local_Dir_Y,
        data.Reco_Local_Dir_Z,
        data.True_PX,      # Still using true values
        data.True_PY,
        data.True_PZ,
        data.Reco_PX,
        data.Reco_PY,
        data.Reco_PZ,
        data.Reco_Pathlength,
        data.True_Pathlength,
        data.True_Origin_T,
        [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    )
end

"""
    get_event_data(reader::EventReader, event_idx::Int)

Reads all tracks and hits for a specific event by index.

# Arguments
- `reader::EventReader`: The event reader instance
- `event_idx::Int`: The index of the event in the events tree

# Returns
- A named tuple with `event`, `particles`, and `hits` fields
"""
function get_event_data(reader::EventReader, event_idx::Int)
    # Read the event entry
    event = reader.events[event_idx]

    # Retrieve the tracks and hits
    particles = get_particles_in_event(reader, event_idx)
    hits = get_hits_in_event(reader, event_idx)

    return (event = event, particles = particles, hits = hits)
end

"""
    get_particles_in_event(reader::EventReader, event_idx::Int)::Vector{Particle}

Retrieves tracks for an event starting at index and with specified count.

# Arguments
- `reader::EventReader`: The event reader instance
- `event_idx::Int`: The index of the event in the events tree

# Returns
- Vector{Particle}: Array of particle objects
"""
function get_particles_in_event(reader::EventReader, event_idx::Int)::Vector{Particle}
    # Get track count and starting index from event
    count = reader.track_counts[event_idx]
    istart = reader.track_indices[event_idx]

    # Pre-allocate the tracks array
    tracks = Vector{Particle}(undef, count)

    for i = 1:count
        track_idx = istart + i - 1

        # Get track data and create particle
        track_data = reader.tracks[track_idx]
        tracks[i] = get_particle(track_data, useTruth = reader.useTruth)
    end

    return tracks
end

"""
    get_hits_in_event(reader::EventReader, event_idx::Int)::Vector{PhotonHit}

Retrieves photon hits for an event starting at index and with specified count.

# Arguments
- `reader::EventReader`: The event reader instance
- `event_idx::Int`: The index of the event in the events tree

# Returns
- Vector{PhotonHit}: Array of photon hits
"""
function get_hits_in_event(reader::EventReader, event_idx::Int)::Vector{PhotonHit}
    # Get hit count and starting index from event
    count = reader.hit_counts[event_idx]
    istart = reader.hit_indices[event_idx]

    # Pre-allocate the hits array
    hits_view = view(reader.hitPool.data, 1:count)

    for i = 1:count
        hit_idx = istart + i - 1
        hit_data = reader.hits[hit_idx]
        update_photon_hit!(hits_view[i], hit_data)
    end

    return hits_view
end

#
# Event Processing Functions
#

"""
    event_iterator(reader::EventReader)

Creates an iterator to process events one by one with their associated tracks and hits.

# Arguments
- `reader::EventReader`: The event reader instance

# Returns
- An iterator yielding (event, tracks, hits) for each event
"""
function event_iterator(reader::EventReader)
    return (get_event_data(reader, i) for i = 1:length(reader.events))
end

"""
    process_all_events(reader::EventReader, processor_fn::Function)

Processes all events with a provided function.

# Arguments
- `reader::EventReader`: The event reader instance
- `processor_fn::Function`: Function taking (event, tracks, hits) and returning a result

# Returns
- Vector of results from the processor function
"""
function process_all_events(reader::EventReader, processor_fn::Function)
    results = Vector{Any}(undef, length(reader.events))

    for (i, event_data) in enumerate(event_iterator(reader))
        results[i] = processor_fn(event_data.event, event_data.particles, event_data.hits)
    end

    return results
end

"""
    process_events_batch(reader::EventReader, start_idx::Int, batch_size::Int, processor_fn::Function)

Processes a batch of events for efficiency.

# Arguments
- `reader::EventReader`: The event reader instance
- `start_idx::Int`: The starting event index
- `batch_size::Int`: Number of events to process
- `processor_fn::Function`: Function taking (event, tracks, hits) and returning a result

# Returns
- Vector of results from the processor function
"""
function process_events_batch(
    reader::EventReader,
    start_idx::Int,
    batch_size::Int,
    processor_fn::Function,
)
    end_idx = min(start_idx + batch_size - 1, length(reader.events))
    actual_size = end_idx - start_idx + 1

    results = Vector{Any}(undef, actual_size)

    for i = 1:actual_size
        event_idx = start_idx + i - 1
        event_data = get_event_data(reader, event_idx)

        results[i] = processor_fn(event_data.event, event_data.particles, event_data.hits)
    end

    return results
end

"""
    process_all_events_parallel(reader::EventReader, processor_fn::Function)

Processes all events in parallel with multithreading.

# Arguments
- `reader::EventReader`: The event reader instance
- `processor_fn::Function`: Function taking (event, tracks, hits) and returning a result

# Returns
- Vector of results from the processor function
"""
function process_all_events_parallel(reader::EventReader, processor_fn::Function)
    n_events = length(reader.events)
    results = Vector{Any}(undef, n_events)

    Threads.@threads for i = 1:n_events
        event_data = get_event_data(reader, i)
        results[i] = processor_fn(event_data.event, event_data.particles, event_data.hits)
    end

    return results
end

"""
    process_events_batch_parallel(reader::EventReader, start_idx::Int, batch_size::Int, processor_fn::Function)

Processes a batch of events in parallel with multithreading.

# Arguments
- `reader::EventReader`: The event reader instance
- `start_idx::Int`: The starting event index
- `batch_size::Int`: Number of events to process
- `processor_fn::Function`: Function taking (event, tracks, hits) and returning a result

# Returns
- Vector of results from the processor function
"""
function process_events_batch_parallel(
    reader::EventReader,
    start_idx::Int,
    batch_size::Int,
    processor_fn::Function,
)
    end_idx = min(start_idx + batch_size - 1, length(reader.events))
    actual_size = end_idx - start_idx + 1

    results = Vector{Any}(undef, actual_size)

    Threads.@threads for i = 1:actual_size
        event_idx = start_idx + i - 1
        event_data = get_event_data(reader, event_idx)

        results[i] = processor_fn(event_data.event, event_data.particles, event_data.hits)
    end

    return results
end
