mutable struct PhotonHit
    localX::Float64
    localY::Float64
    localZ::Float64
    arrivalTime::Float64
    energy::Float64

    moduleId::Int32
    trackId::Int32

    isDownward::Bool
    isSurfaceScattered::Bool
    isRayleighScattered::Bool

    cherenkovTheta::Union{Float64,Nothing}
    cherenkovPhi::Union{Float64,Nothing}
    emissionX::Union{Float64,Nothing}
    emissionY::Union{Float64,Nothing}
    emissionZ::Union{Float64,Nothing}
    emissionTime::Union{Float64,Nothing}

    function PhotonHit(
        localX,
        localY,
        localZ,
        arrivalTime,
        energy,
        moduleId,
        trackId,
        isDownward,
        isSurfaceScattered,
        isRayleighScattered,
    )
        return new(
            localX,
            localY,
            localZ,
            arrivalTime,
            energy,
            moduleId,
            trackId,
            isDownward,
            isSurfaceScattered,
            isRayleighScattered,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
        )
    end
end

function update_photon_hit!(hit::PhotonHit, hit_data; Debug::Bool = false)
    hit.localX = hit_data.LocalX
    hit.localY = hit_data.LocalY
    hit.localZ = hit_data.LocalZ
    hit.arrivalTime = hit_data.ArrivalTime
    hit.energy = hit_data.Energy
    hit.moduleId = hit_data.Module
    hit.trackId = hit_data.Track
    hit.isDownward = hit_data.IsDownward
    hit.isSurfaceScattered = hit_data.IsSurfaceScattered
    hit.isRayleighScattered = hit_data.IsRayleighScattered

    if Debug
        hit.cherenkovTheta = hit_data.CherenkovTheta
        hit.cherenkovPhi = hit_data.CherenkovPhi
        hit.emissionX = hit_data.EmissionX
        hit.emissionY = hit_data.EmissionY
        hit.emissionZ = hit_data.EmissionZ
        hit.emissionTime = hit_data.EmissionTime
    end

    return hit
end
