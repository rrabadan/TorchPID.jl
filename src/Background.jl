abstract type AbstractBackgroundModel end

struct FlatBackground <: AbstractBackgroundModel
    detector::Detector
    prob_pixel::Float64
    prob_hit::Float64
end

function FlatBackground(detector::Detector = DETECTOR)
    v = (
        detector.n_detectors *
        (detector.active * detector.active) *
        (detector.t_max - detector.t_min)
    )

    # Uniform probability for each pixel, hit
    prob_hit = 1.0 / v
    prob_pixel = 1.0 / (detector.n_total_pixels)
    return FlatBackground(detector, prob_pixel, prob_hit)
end

@inline function get_bkg_probability(background::FlatBackground, pixel::PixelHit)::Float64
    return background.prob_pixel
end

@inline function get_bkg_probability(background::FlatBackground, hit::HitCoordinate)::Float64
    return background.prob_hit
end

function fill_bkg_probabilities!(
    bkg_prob::AbstractVector{T},
    background::AbstractBackgroundModel,
    hits::AbstractVector,
) where {T<:AbstractFloat}
    @inbounds for i in eachindex(hits, bkg_prob)
        bkg_prob[i] = T(get_bkg_probability(background, hits[i]))
    end
    return bkg_prob
end
