using UnROOT
include("Particle.jl")

function photon_columns(; debugPhotons::Bool=false)
    columns = ["Run", "Event", "Module", "Track", "Energy", "ArrivalTime", "LocalX", "LocalY"]
    if debugPhotons
        debug_columns = ["EmissionTime", "CherenkovTheta", "CherenkovPhi", "EmissionX", "EmissionY", "EmissionZ",
            "LocalTrackDirectionX", "LocalTrackDirectionY", "LocalTrackDirectionZ"]
        push!(columns, debug_columns...)
    end
    return columns
end

function track_columns(; useTruth::Bool=false)
    columns = ["Run", "Event", "Track", "True_ID", "True_PX", "True_PY", "True_PZ",
        "Reco_Momentum", "True_Momentum", "Reco_PathLength", "True_PathLength"]
    if useTruth
        truth_columns = ["True_Module", "True_Local_DirX", "True_Local_DirY", "True_Local_DirZ"]
        push!(columns, truth_columns...)
    else
        reco_columns = ["Reco_Module", "Reco_Local_DirX", "Reco_Local_DirY", "Reco_Local_DirZ"]
        push!(columns, reco_columns...)
    end
    return columns
end

struct EventReader
    events::LazyTree
    tracks::LazyTree
    hits::LazyTree

    useTruth::Bool
    debugPhotons::Bool

    #Inner constructor
    function EventReader(file::String; useTruth::Bool=false, debugPhotons::Bool=false)
        f = ROOTFile(file)
        new(LazyTree(f, "TorchTuple/TorchEvent"),
            LazyTree(f, "TorchTuple/TorchTracks"),
            LazyTree(f, "TorchTuple/TorchHits"),
            useTruth, debugPhotons)
    end
end

function get_particle(reader::EventReader, entry::UnROOT.LazyEvent)::Particle
    pid::Int = entry.True_ID
    eventId::UInt = entry.Event
    trackId::UInt = entry.Track
    moduleId::UInt = reader.useTruth ? entry.True_Module : entry.Reco_Module
    pMag::Float64 = entry.Reco_Momentum
    pMagTrue::Float64 = entry.True_Momentum
    xDir::Float64 = reader.useTruth ? entry.True_Local_Dir_X : entry.Reco_Local_Dir_X
    yDir::Float64 = reader.useTruth ? entry.True_Local_Dir_Y : entry.Reco_Local_Dir_Y
    zDir::Float64 = reader.useTruth ? entry.True_Local_Dir_Z : entry.Reco_Local_Dir_Z
    truePX::Float64 = entry.True_PX
    truePY::Float64 = entry.True_PY
    truePZ::Float64 = entry.True_PZ
    recoPX::Float64 = reader.useTruth ? 0.0 : entry.Reco_PX
    recoPY::Float64 = reader.useTruth ? 0.0 : entry.Reco_PY
    recoPZ::Float64 = reader.useTruth ? 0.0 : entry.Reco_PZ
    truePathlength::Float64 = entry.True_PathLength
    recoPathlength::Float64 = reader.useTruth ? 0.0 : entry.Reco_Pathlength
    Particle(pid=pid, eventId=eventId, trackId=trackId, moduleId=moduleId, pMag=pMag, pMagTrue=pMagTrue,
        xDir=xDir, yDir=yDir, zDir=zDir, truePX=truePX, truePY=truePY, truePZ=truePZ,
        recoPX=recoPX, recoPY=recoPY, recoPZ=recoPZ, recoPathlength=recoPathlength, truePathlength=truePathlength)
end