""" 
    PhotonSpectrum

A structure representing the photon spectrum parameters and arrays.
"""
struct PhotonSpectrum
    nbins::Int
    emin::Float64
    emax::Float64
    dE::Float64

    energy::Vector{Float64}
    efficiency::Vector{Float64}
    nphase::Vector{Float64}
    ngroup::Vector{Float64}
    vgroup::Vector{Float64}
    nphase_edge::Vector{Float64}
end

""" 
    PhotonSpectrum(dht::DetectorHitTester; nbins, emin, emax)

Construct a PhotonSpectrum from a DetectorHitTester using specified energy bin parameters.
"""
function PhotonSpectrum(
    dht::DetectorHitTester;
    nbins::Int = 525,
    emin::Float64 = 1.75,
    emax::Float64 = 7.00,
)

    dE = (emax - emin) / nbins
    energy = Vector{Float64}(undef, nbins)
    efficiency = Vector{Float64}(undef, nbins)
    nphase = Vector{Float64}(undef, nbins)
    ngroup = Vector{Float64}(undef, nbins)
    nphase_edge = Vector{Float64}(undef, nbins)

    for i = 1:nbins
        # In C++: E = m_emin + ( i + 0.5 ) * m_dE, where i starts at 0.
        # In Julia, adjust the index:
        E = emin + ((i - 1) + 0.5) * dE

        energy[i] = E
        efficiency[i] = efficiency(dht, E)             # call the efficiency function on dht
        nphase[i] = nphase_Corning(E)                  # using TorchFunctions.jl function
        ngroup[i] = ngroup_Corning(E, nphase[i])       # note: using nphase[i] computed above
        nphase_edge[i] = nphase_Corning(E - 0.5 * dE)
    end

    vgroup = CLIGHT ./ ngroup

    PhotonSpectrum(
        nbins,
        emin,
        emax,
        dE,
        energy,
        efficiency,
        nphase,
        ngroup,
        vgroup,
        nphase_edge,
    )
end

""" 
    PhotonSpectrum(; nbins, emin, emax)

Construct a PhotonSpectrum using a default DetectorHitTester with specified energy parameters.
"""
function PhotonSpectrum(; nbins::Int = 525, emin::Float64 = 1.75, emax::Float64 = 7.00)
    dht = DetectorHitTester(emin = emin, emax = emax)
    PhotonSpectrum(dht, nbins = nbins, emin = emin, emax = emax)
end

""" 
    PhotonSpectrumDistribution

A structure representing the photon spectrum distribution details.
"""
struct PhotonSpectrumDistribution
    beta::Float64
    yield_per_mm::Float64
    distribution::Vector{Float64}
    cumulative::Vector{Float64}
end

""" 
    PhotonSpectrumDistribution(s::PhotonSpectrum, beta::Float64)

Generate a photon spectrum distribution from a PhotonSpectrum using associated particle's beta parameter.
"""
function PhotonSpectrumDistribution(s::PhotonSpectrum, beta::Float64)
    #distribution = Vector{Float64}(undef, s.nbins)
    cumulative = Vector{Float64}(undef, s.nbins)

    distribution = map(
        (v1, v2) -> begin
            bn = beta * v1
            sinsqthetac = bn > 1.0 ? 1.0 - (1.0 / bn)^2 : 0.0
            37.0 * s.dE * sinsqthetac * v2
        end,
        s.nphase,
        s.efficiency,
    )

    cumsum!(cumulative, distribution)

    yield_per_mm = cumulative[end]
    cumulative ./= yield_per_mm

    PhotonSpectrumDistribution(beta, yield_per_mm, distribution, cumulative)
end

""" 
    spectrum_nphase(s::PhotonSpectrum, energy::Float64)::Float64

Return the nphase value corresponding to the given energy in the photon spectrum.
"""
function spectrum_nphase(s::PhotonSpectrum, energy::Float64)::Float64
    n = floor(Int, (energy - s.emin) / s.dE) + 1
    return (n ≥ 1 && n ≤ s.nbins) ? s.nphase[n] : 0.0
end

""" 
    spectrum_ngroup(s::PhotonSpectrum, energy::Float64)::Float64

Return the ngroup value corresponding to the given energy in the photon spectrum.
"""
function spectrum_ngroup(s::PhotonSpectrum, energy::Float64)::Float64
    n = floor(Int, (energy - s.emin) / s.dE) + 1
    return (n ≥ 1 && n ≤ s.nbins) ? s.ngroup[n] : 0.0
end

""" 
    spectrum_yield(d::PhotonSpectrumDistribution, pathlength::Float64)::Float64

Calculate the total yield over the specified pathlength.
"""
function spectrum_yield(d::PhotonSpectrumDistribution, pathlength::Float64)::Float64
    return pathlength * d.yield_per_mm
end

""" 
    spectrum_random_energy(s::PhotonSpectrum, d::PhotonSpectrumDistribution)::Float64

Generate a random energy value sampled from the photon spectrum distribution.
"""
function spectrum_random_energy(s::PhotonSpectrum, d::PhotonSpectrumDistribution)::Float64
    if d.yield_per_mm == 0.0
        return 0.0
    end
    v = rand()
    i = searchsortedfirst(d.cumulative, v)
    if i <= length(d.cumulative)
        n = i - 1
        l = i > 1 ? d.cumulative[i-1] : 0.0
        r = d.cumulative[i]
        return s.energy[n] + s.dE * (v - l) / (r - l)
    end
end

""" 
    spectrum_probability(s::PhotonSpectrum, d::PhotonSpectrumDistribution, energy::Float64)::Float64

Return the probability density for a given energy based on the spectrum distribution.
"""
function spectrum_probability(
    s::PhotonSpectrum,
    d::PhotonSpectrumDistribution,
    energy::Float64,
)::Float64
    n = floor(Int, (energy - s.emin) / s.dE) + 1
    return (n ≥ 1 && n ≤ s.nbins) ? d.distribution[n] / d.yield_per_mm : 0.0
end

"""
    spectrum_above_threshold(d::PhotonSpectrumDistribution) -> Bool

Check if the photon spectrum distribution has a yield per millimeter greater than zero.

# Arguments
- `d::PhotonSpectrumDistribution`: The photon spectrum distribution to check.

# Returns
- `Bool`: `true` if the yield per millimeter is greater than zero, `false` otherwise.
"""
function spectrum_above_threshold(d::PhotonSpectrumDistribution)::Bool
    return d.yield_per_mm > 0.0
end
