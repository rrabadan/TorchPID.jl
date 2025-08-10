
<a id='PhotonSpectrum-API'></a>

<a id='PhotonSpectrum-API-1'></a>

# PhotonSpectrum API




Tools to model the Cherenkov radiation photon spectra within the TORCH detector


<a id='Public-API'></a>

<a id='Public-API-1'></a>

### Public API

<a id='TorchPID.PhotonSpectrum' href='#TorchPID.PhotonSpectrum'>#</a>
**`TorchPID.PhotonSpectrum`** &mdash; *Type*.



```julia
PhotonSpectrum(nbins, emin, emax, dE, energy, efficiency, nphase, ngroup, vgroup, nphase_edge)
```

Type representing a discretized energy spectrum of Cherenkov photons, including detection efficiency and refractive indices derived from TORCH radiator materials.

**Fields**

  * `nbins::Int`: Number of energy bins in the spectrum.
  * `emin::Float64`: Minimum energy value in the spectrum.
  * `emax::Float64`: Maximum energy value in the spectrum.
  * `dE::Float64`: Energy bin width.
  * `energy::Vector{Float64}`: Array of energy values for each bin.
  * `efficiency::Vector{Float64}`: Array of detection efficiency values corresponding to each energy bin.
  * `nphase::Vector{Float64}`: Array of phase refractive index values for each energy bin.
  * `ngroup::Vector{Float64}`: Array of group refractive index values for each energy bin.
  * `vgroup::Vector{Float64}`: Array of group velocities for each energy bin.
  * `nphase_edge::Vector{Float64}`: Array of phase refractive index values at the left edge of each energy bin.

**Constructors**

```
PhotonSpectrum(dht::DetectorHitTester; nbins::Int = 525, emin::Float64 = 1.75, emax::Float64 = 7.00)
PhotonSpectrum(; nbins::Int=525, emin::Float64=1.75, emax::Float64=7.00)
```

Constructs a Cherenkov photon spectrum over a specified energy range, with properties computed for each energy bin. The refractive indices are derived from TORCH radiator materials, and detection efficiencies are determined using the provided `DetectorHitTester` object. A default `DetectorHitTester` is created if not provided.

**Arguments**

  * `dht::DetectorHitTester`: Detector hit tester object containing efficiency parameters.

**Keywords**

  * `nbins::Int=525`: Number of energy bins in the spectrum (default: 525).
  * `emin::Float64=1.75`: Minimum energy value in the spectrum (default: 1.75).
  * `emax::Float64=7.00`: Maximum energy value in the spectrum (default: 7.00).

**Examples**

```julia
# Create a PhotonSpectrum with default parameters
spectrum = PhotonSpectrum()

# Using convenience method is to construct a PhotonSpectrum
# with a default DetectorHitTester and specified energy parameters.
spectrum = PhotonSpectrum(nbins=1000, emin=2.0, emax=6.0)
```


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L1-L46' class='documenter-source'>source</a><br>

<a id='TorchPID.PhotonSpectrumDistribution' href='#TorchPID.PhotonSpectrumDistribution'>#</a>
**`TorchPID.PhotonSpectrumDistribution`** &mdash; *Type*.



```julia
PhotonSpectrumDistribution(beta, yield_per_mm, distribution, cumulative)
```

Type encapsulating the photon energy emission probability distribution for a charged particle traversing a TORCH radiator.

**Fields**

  * `beta::Float64`: Relativistic beta factor of the particle.
  * `yield_per_mm::Float64`: Expected photon yield per millimeter of radiator.
  * `distribution::Vector{Float64}`: Probability distribution of photon energies.
  * `cumulative::Vector{Float64}`: Cumulative distribution function for sampling photon energies.

**Constructors**

```
PhotonSpectrumDistribution(s::PhotonSpectrum, beta::Float64)
```

Constructs a distribution from a `PhotonSpectrum` and the particle's relativistic beta factor. The distribution is calculated using the Cherenkov emission formula, weighted by the detector efficiency. The photon yield is based on a factor of 37.0 photons per eV per mm as the base rate.

**Arguments**

  * `s::PhotonSpectrum`: The photon spectrum containing energy bin information.
  * `beta::Float64`: Relativistic beta factor (v/c) of the particle.

**Returns**

A new `PhotonSpectrumDistribution` object with photon energy emission probabilities calculated.

**Examples**

```julia
# Create a spectrum
spectrum = PhotonSpectrum()

# Create a photon energies probability distribution for a particle with beta = 0.99
distribution = PhotonSpectrumDistribution(spectrum, 0.99)

# Check if the distribution will produce photons
if spectrum_above_threshold(distribution)
    # Calculate expected yield over 10mm path
    yield = spectrum_yield(distribution, 10.0)
    println("Expected photon yield: $yield")
    
    # Sample a random energy from the distribution
    energy = spectrum_random_energy(spectrum, distribution)
end
```


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L96-L142' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_above_threshold-Tuple{PhotonSpectrumDistribution}' href='#TorchPID.spectrum_above_threshold-Tuple{PhotonSpectrumDistribution}'>#</a>
**`TorchPID.spectrum_above_threshold`** &mdash; *Method*.



```julia
spectrum_above_threshold(d::PhotonSpectrumDistribution)
```

`spectrum_above_threshold` determines whether the photon spectrum distribution has a non-zero yield per millimeter.      This indicates if the particle's velocity exceeds the Cherenkov threshold, allowing photon production.

**Arguments**

  * `d::PhotonSpectrumDistribution`: The photon spectrum distribution to check.

**Returns**

  * `Bool`: `true` if the yield per millimeter is greater than zero, `false` otherwise.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L274-L285' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_ngroup-Tuple{PhotonSpectrum, Float64}' href='#TorchPID.spectrum_ngroup-Tuple{PhotonSpectrum, Float64}'>#</a>
**`TorchPID.spectrum_ngroup`** &mdash; *Method*.



```julia
spectrum_ngroup(s::PhotonSpectrum, energy::Float64)
```

`spectrum_ngroup` returns the group refractive index value corresponding to the given energy in the photon spectrum.

**Arguments**

  * `s::PhotonSpectrum`: The photon spectrum object.
  * `energy::Float64`: The energy value for which to find the group refractive index.

**Returns**

  * `Float64`: The group refractive index value for the given energy, or 0.0 if energy is outside the spectrum range.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L188-L200' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_nphase-Tuple{PhotonSpectrum, Float64}' href='#TorchPID.spectrum_nphase-Tuple{PhotonSpectrum, Float64}'>#</a>
**`TorchPID.spectrum_nphase`** &mdash; *Method*.



```julia
spectrum_nphase(s::PhotonSpectrum, energy::Float64)
```

`spectrum_nphase` returns the phase refractive index value corresponding to the given energy in the photon spectrum.

**Arguments**

  * `s::PhotonSpectrum`: The photon spectrum object.
  * `energy::Float64`: The energy value for which to find the phase refractive index.

**Returns**

  * `Float64`: The phase refractive index value for the given energy, or 0.0 if energy is outside the spectrum range.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L171-L183' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_probability-Tuple{PhotonSpectrum, PhotonSpectrumDistribution, Float64}' href='#TorchPID.spectrum_probability-Tuple{PhotonSpectrum, PhotonSpectrumDistribution, Float64}'>#</a>
**`TorchPID.spectrum_probability`** &mdash; *Method*.



```julia
spectrum_probability(s::PhotonSpectrum, d::PhotonSpectrumDistribution, energy::Float64)
```

`spectrum_probability` computes the probability density for a given energy based on the photon spectrum distribution.     It identifies the energy bin corresponding to the input energy and retrieves the normalized probability density for that bin.      If the energy is outside the defined spectrum range, the function returns 0.0.

**Arguments**

  * `s::PhotonSpectrum`: The photon spectrum object containing energy bin information.
  * `d::PhotonSpectrumDistribution`: The photon spectrum distribution with normalized probabilities.
  * `energy::Float64`: The energy value for which to calculate the probability.

**Returns**

  * `Float64`: The probability density at the given energy, normalized by the yield per millimeter, or 0.0 if energy is outside the spectrum range.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L249-L265' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_random_energy-Tuple{PhotonSpectrum, PhotonSpectrumDistribution}' href='#TorchPID.spectrum_random_energy-Tuple{PhotonSpectrum, PhotonSpectrumDistribution}'>#</a>
**`TorchPID.spectrum_random_energy`** &mdash; *Method*.



```julia
spectrum_random_energy(s::PhotonSpectrum, d::PhotonSpectrumDistribution)
```

`spectrum_random_energy` samples a random energy value from the photon spectrum distribution using the cumulative distribution function.      It performs linear interpolation within energy bins to ensure accurate and weighted sampling based on the distribution.

**Arguments**

  * `s::PhotonSpectrum`: The photon spectrum object containing energy bin information.
  * `d::PhotonSpectrumDistribution`: The photon spectrum distribution for sampling.

**Returns**

  * `Float64`: A randomly sampled energy value from the distribution, or 0.0 if yield is below threshold.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L221-L234' class='documenter-source'>source</a><br>

<a id='TorchPID.spectrum_yield-Tuple{PhotonSpectrumDistribution, Float64}' href='#TorchPID.spectrum_yield-Tuple{PhotonSpectrumDistribution, Float64}'>#</a>
**`TorchPID.spectrum_yield`** &mdash; *Method*.



```julia
spectrum_yield(s::PhotonSpectrum, energy::Float64)
```

`spectrum_yield` calculates the total expected photon yield over the specified pathlength.

**Arguments**

  * `d::PhotonSpectrumDistribution`: The photon spectrum distribution.
  * `pathlength::Float64`: Path length in millimeters traversed by the particle.

**Returns**

  * `Float64`: Expected total number of photons produced over the pathlength.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/5f20ea0e22a6826c96b767257621959b44c97e4a/src/PhotonSpectrum.jl#L205-L217' class='documenter-source'>source</a><br>

