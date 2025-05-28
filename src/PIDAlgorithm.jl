struct PIDAlg
    nModules::Int

    pids::Vector{Int}
    momentumRange::Tuple{Float64,Float64}
    checkTkAcceptance::Bool
    checkTkQuality::Bool
    smearT0::Bool

    blackenedBottom::Bool
    removeFocusPhotons::Bool
    photonTrackMatching::Bool
end

function PIDAlg(;
    nModules::Int = 18,
    pids::Vector{Int} = [211, 321, 2212],
    momentumRange::Tuple{Float64,Float64} = (1.0, 20.0),
    checkTkAcceptance::Bool = true,
    checkTkQuality::Bool = false,
    smearT0::Bool = false,
    blackenedBottom::Bool = false,
    removeFocusPhotons::Bool = true,
    photonTrackMatching::Bool = false,
)::PIDAlg
    PIDAlg(
        nModules,
        pids,
        momentumRange,
        checkTkAcceptance,
        checkTkQuality,
        smearT0,
        blackenedBottom,
        removeFocusPhotons,
        photonTrackMatching,
    )
end

"""
    find_particles(pidalg::PIDAlg, particles::Vector{Particle})::Dict{Int8, Vector{Particle}}

Filters and groups particles based on the provided `PIDAlg`.

# Arguments
- `pidalg::PIDAlg`: A pidalguration object that specifies filtering criteria, such as momentum range and a list of particle IDs (PIDs).
- `particles::Vector{Particle}`: A vector containing the particles to be filtered and grouped.

# Returns
- `Dict{Int8, Vector{Particle}}`: A dictionary where the keys represent module IDs (as `Int8`), and the values are vectors of particles (`Vector{Particle}`) that belong to the corresponding modules.
"""
function find_particles(
    pidalg::PIDAlg,
    particles::Vector{Particle},
)::Dict{Int,Vector{Particle}}
    module_particles = Dict{Int,Vector{Particle}}()
    for i = 1:pidalg.nModules
        module_particles[i] = Vector{Particle}()
    end
    pmin, pmax = pidalg.momentumRange

    for p in particles
        # Apply all filters in one go
        if abs(p.pid) in pidalg.pids &&
           (!pidalg.checkTkAcceptance || (pmin <= p.pMag <= pmax))

            moduleId = p.moduleId
            push!(module_particles[moduleId], p)
        end
    end
    filter!(pair -> !isempty(pair.second), module_particles)
    return module_particles
end

function find_coordinates(
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    photons::Vector{PhotonHit},
)::Dict{Int,Vector{HitCoordinate}}
    # Pre-allocate dictionary
    module_hits = Dict{Int,Vector{HitCoordinate}}()
    for i = 1:pidalg.nModules
        #module_hits[i] = Vector{HitCoordinate}(undef, 0)
        #sizehint!(module_hits[i], expected_capacity)
        module_hits[i] = Vector{HitCoordinate}()
    end

    # Combine filtering and processing in a single loop for better performance
    for p in photons
        # Only process if it passes selection criteria
        if select_photon_hit(pidalg, RADIATOR, WEDGE, DETECTOR, MASK, p) &&
           test_photon(dht, p.energy)
            moduleId = p.Module
            #    hit = HitCoordinate(p.LocalX, p.LocalY, p.arrivalTime)
            #    push!(module_hits[moduleId], hit)
        end
    end

    # Remove empty modules
    filter!(pair -> !isempty(pair.second), module_hits)
    return module_hits
end

function select_photon_hit(
    pidalg::PIDAlg,
    radiator::Radiator,
    wedge::Wedge,
    detector::Detector,
    mask::Mask,
    photon::PhotonHit,
)::Bool
    # Combine checks in a logical AND expression to benefit from short-circuit evaluation
    return photon_on_detector(photon, detector) &&
           !(pidalg.removeFocusPhotons && photon_from_focus(photon, radiator, wedge)) &&
           !is_photon_masked(photon, mask) &&
           !(pidalg.blackenedBottom && photon.isDownward) &&
           !(detector.t_window_min < photon.arrivalTime < detector.t_window_max)
end

function run_algorithm(
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    particles::Vector{Particle},
    photons::Vector{PhotonHit},
)
    # Find particles and photons in the modules
    #module_particles = find_particles(pidalg, particles)
    module_hits = find_coordinates(pidalg, dht, photons)

end
