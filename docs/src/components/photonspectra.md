# Cherenkov Radiation Spectra

TorchPID.jl provides tools to model the Cherenkov radiation photon spectra within the TORCH detector.

These tools facilitate the generation and sampling of photon energy distributions, incorporating detection efficiency and refractive index variations caused by charged particles traversing the TORCH radiator.

## Key components

* PhotonSpectrum
* DetectorHitTester
* TorchFunctions

### Creating a Discretized Energy Spectrum

To model the Cherenkov photon energy spectrum, the `PhotonSpectrum` struct provides a discretized representation of photon energies. This includes detector efficiency values and refractive indices for each energy point. Below is an example of how to create and configure a `PhotonSpectrum`:

```julia
spectrum = PhotonSpectrum(;
    nbins = 525,  # Number of energy points
    emin = 1.75,  # Minimum photon energy (eV)
    emax = 7.00,  # Maximum photon energy (eV)
)
```
This discretized spectrum forms the foundation for further calculations, such as photon yield estimation, energy sampling, and threshold detection.

Generate the photon energy emission probability distribution for a particle (characterized by $\beta = 0.99$)
```julia
beta = 0.99
distribution = PhotonSpectrumDistribution(
    spectrum
    beta
)
```

Sample a random energy value from the photon spectrum distribution
```julia
photon_energy = spectrum_random_energy(
    spectrum,
    distribution
)
```

## Photon Detection Efficiency

Photon emission probabilities are calculated accounting for detection efficiency for a given photon energy. This considers quantum detector efficiency, collection efficiency, and mirror reflectivity.
`DetectorHitTester` encapsulates the configuration for photon detection efficiency.
```julia
dht = DetectorHitTester(
    emin::Float64 = 1.75,
    emax::Float64 = 7.00,
    implement_epotek_cutoff::Bool = true,
    implement_imperfect_mirror::Bool = true,
    implement_QE::Bool = true,
    implement_CE::Bool = true,
)
```

It can be used to construct a custom `PhotonSpectrum`
```julia
spectrum = PhotonSpectrum(
    dht;
    nbins = 525,  # Number of energy points
    emin = 1.75,  # Minimum photon energy (eV)
    emax = 7.00,  # Maximum photon energy (eV)
)

beta = 0.99
distribution = PhotonSpectrumDistribution(
    spectrum
    beta
)

photon_energy = spectrum_random_energy(
    spectrum,
    distribution
)
```

## TorchFunctions

Collection of functions to derive refraction indices derived from the radiatior's material properties.
It uncludes functionality to compute quatum efficiency and transmissions probabilities as function of the cherenkov photon's energy. 