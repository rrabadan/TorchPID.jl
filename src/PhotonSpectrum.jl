""" 
Struct representing the photon spectrum parameters and arrays.

# Fields
- `nbins::Int`: Number of energy bins in the spectrum.
- `emin::Float64`: Minimum energy value in the spectrum.
- `emax::Float64`: Maximum energy value in the spectrum.
- `dE::Float64`: Energy bin width.
- `energy::Vector{Float64}`: Array of energy values for each bin.
- `efficiency::Vector{Float64}`: Array of detection efficiency values corresponding to each energy bin.
- `nphase::Vector{Float64}`: Array of phase refractive index values for each energy bin.
- `ngroup::Vector{Float64}`: Array of group refractive index values for each energy bin.
- `vgroup::Vector{Float64}`: Array of group velocities for each energy bin.
- `nphase_edge::Vector{Float64}`: Array of phase refractive index values at the left edge of each energy bin.
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
Constructs a PhotonSpectrum from a DetectorHitTester using specified energy bin parameters.

# Arguments
- `dht::DetectorHitTester`: Detector hit tester object containing efficiency parameters.

# Keywords
- `nbins::Int=525`: Number of energy bins in the spectrum (default: 525).
- `emin::Float64=1.75`: Minimum energy value in the spectrum (default: 1.75).
- `emax::Float64=7.00`: Maximum energy value in the spectrum (default: 7.00).

# Returns
A new `PhotonSpectrum` object with calculated spectrum properties.
"""
function PhotonSpectrum(
    dht::DetectorHitTester;
    nbins::Int = 525,
    emin::Float64 = 1.75,
    emax::Float64 = 7.00,
)

    dE = (emax - emin) / nbins
    energy = Vector{Float64}(undef, nbins)
    eff = Vector{Float64}(undef, nbins)
    nphase = Vector{Float64}(undef, nbins)
    ngroup = Vector{Float64}(undef, nbins)
    nphase_edge = Vector{Float64}(undef, nbins)

    for i = 1:nbins
        # In C++: E = m_emin + ( i + 0.5 ) * m_dE, where i starts at 0.
        # In Julia, adjust the index:
        E = emin + ((i - 1) + 0.5) * dE

        energy[i] = E
        eff[i] = efficiency(dht, E)             # call the efficiency function on dht
        nphase[i] = nphase_Corning(E)                  # using TorchFunctions.jl function
        ngroup[i] = ngroup_Corning(E, nphase[i])       # note: using nphase[i] computed above
        nphase_edge[i] = nphase_Corning(E - 0.5 * dE)
    end

    vgroup = CLIGHT ./ ngroup

    PhotonSpectrum(nbins, emin, emax, dE, energy, eff, nphase, ngroup, vgroup, nphase_edge)
end

""" 
Constructs a PhotonSpectrum using a default DetectorHitTester with specified energy parameters.

# Keywords
- `nbins::Int=525`: Number of energy bins in the spectrum (default: 525).
- `emin::Float64=1.75`: Minimum energy value in the spectrum (default: 1.75).
- `emax::Float64=7.00`: Maximum energy value in the spectrum (default: 7.00).

# Returns
A new `PhotonSpectrum` object with calculated spectrum properties using default detector parameters.
"""
function PhotonSpectrum(; nbins::Int = 525, emin::Float64 = 1.75, emax::Float64 = 7.00)
    dht = DetectorHitTester(emin = emin, emax = emax)
    PhotonSpectrum(dht, nbins = nbins, emin = emin, emax = emax)
end

""" 
Struct representing the photon spectrum distribution details.

# Fields
- `beta::Float64`: Relativistic beta factor of the particle.
- `yield_per_mm::Float64`: Expected photon yield per millimeter of radiator.
- `distribution::Vector{Float64}`: Probability distribution of photon energies.
- `cumulative::Vector{Float64}`: Cumulative distribution function for sampling photon energies.
"""
struct PhotonSpectrumDistribution
    beta::Float64
    yield_per_mm::Float64
    distribution::Vector{Float64}
    cumulative::Vector{Float64}
end

""" 
Generates a photon spectrum distribution from a PhotonSpectrum using associated particle's beta parameter.

The distribution is weighted by the Cherenkov emission formula and detector efficiency.
The yield calculation uses a factor of 37.0 photons per eV per mm as a base rate.

# Arguments
- `s::PhotonSpectrum`: The photon spectrum containing energy bin information.
- `beta::Float64`: Relativistic beta factor (v/c) of the particle.

# Returns
A new `PhotonSpectrumDistribution` object with calculated distribution properties.
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
Returns the phase refractive index value corresponding to the given energy in the photon spectrum.

# Arguments
- `s::PhotonSpectrum`: The photon spectrum object.
- `energy::Float64`: The energy value for which to find the phase refractive index.

# Returns
- `Float64`: The phase refractive index value for the given energy, or 0.0 if energy is outside the spectrum range.
"""
function spectrum_nphase(s::PhotonSpectrum, energy::Float64)::Float64
    n = floor(Int, (energy - s.emin) / s.dE) + 1
    return (n ≥ 1 && n ≤ s.nbins) ? s.nphase[n] : 0.0
end

""" 
Returns the group refractive index value corresponding to the given energy in the photon spectrum.

# Arguments
- `s::PhotonSpectrum`: The photon spectrum object.
- `energy::Float64`: The energy value for which to find the group refractive index.

# Returns
- `Float64`: The group refractive index value for the given energy, or 0.0 if energy is outside the spectrum range.
"""
function spectrum_ngroup(s::PhotonSpectrum, energy::Float64)::Float64
    n = floor(Int, (energy - s.emin) / s.dE) + 1
    return (n ≥ 1 && n ≤ s.nbins) ? s.ngroup[n] : 0.0
end

""" 
Calculates the total expected photon yield over the specified pathlength.

# Arguments
- `d::PhotonSpectrumDistribution`: The photon spectrum distribution.
- `pathlength::Float64`: Path length in millimeters traversed by the particle.

# Returns
- `Float64`: Expected total number of photons produced over the pathlength.
"""
function spectrum_yield(d::PhotonSpectrumDistribution, pathlength::Float64)::Float64
    return pathlength * d.yield_per_mm
end

""" 
Generates a random energy value sampled from the photon spectrum distribution.

Uses the cumulative distribution function for proper weighted sampling.
Performs linear interpolation within energy bins for accurate sampling.

# Arguments
- `s::PhotonSpectrum`: The photon spectrum object containing energy bin information.
- `d::PhotonSpectrumDistribution`: The photon spectrum distribution for sampling.

# Returns
- `Float64`: A randomly sampled energy value from the distribution, or 0.0 if yield is below threshold.
"""
function spectrum_random_energy(s::PhotonSpectrum, d::PhotonSpectrumDistribution)::Float64
    energy = 0.0
    if spectrum_above_threshold(d)
        v = rand()
        i = searchsortedfirst(d.cumulative, v)
        if i <= length(d.cumulative)
            n = i - 1
            l = i > 1 ? d.cumulative[i-1] : 0.0
            r = d.cumulative[i]
            energy = s.emin + s.dE * (n + (v - l) / (r - l))
        end
    end
    return energy
end

""" 
Returns the probability density for a given energy based on the spectrum distribution.

The function determines which energy bin the provided energy falls into.
Returns the corresponding normalized distribution value from that bin.
Returns 0.0 if the energy value is outside the spectrum's defined energy range.

# Arguments
- `s::PhotonSpectrum`: The photon spectrum object containing energy bin information.
- `d::PhotonSpectrumDistribution`: The photon spectrum distribution with normalized probabilities.
- `energy::Float64`: The energy value for which to calculate the probability.

# Returns
- `Float64`: The probability density at the given energy, normalized by the yield per millimeter,
  or 0.0 if energy is outside the spectrum range.
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
Checks if the photon spectrum distribution has a yield per millimeter greater than zero.

This function is used to determine if a distribution will produce any photons.
A distribution with zero yield indicates the particle is below the Cherenkov threshold.

# Arguments
- `d::PhotonSpectrumDistribution`: The photon spectrum distribution to check.

# Returns
- `Bool`: `true` if the yield per millimeter is greater than zero, `false` otherwise.
"""
function spectrum_above_threshold(d::PhotonSpectrumDistribution)::Bool
    return d.yield_per_mm > 0.0
end
