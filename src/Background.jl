struct FlatBackground
    detector::Detector
    prob_pixel::Float64
    prob_hit::Float64
end

function FlatBackground(detector::Detector)
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

function get_bkg_probability(background::FlatBackground, pixel::PixelHit)
    return background.prob_pixel
end

function get_bkg_probability(background::FlatBackground, hit::HitCoordinate)
    return background.prob_hit
end
