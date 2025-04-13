"""
Physical constants used throughout the TORCH detector simulation.

Constants are defined with specific units relevant to the simulation:
- Distances in mm
- Times in ns
- Energies in GeV or eV as specified
"""

"""Speed of light in vacuum (mm/ns)"""
const CLIGHT = 299.792458

"""Elementary charge (C)"""
const QQe = 1.60217662e-19

"""Pion mass (GeV/c^2)"""
const PION_MASS = 0.139570

"""Kaon mass (GeV/c^2)"""
const KAON_MASS = 0.493677

"""Proton mass (GeV/c^2)"""
const PROTON_MASS = 0.938272

"""Electron mass (GeV/c^2)"""
const ELECTRON_MASS = 0.000511

"""Muon mass (GeV/c^2)"""
const MUON_MASS = 0.105658

"""Energy to wavelength conversion factor (nm * eV)"""
const LAMBDA = 1239.84193

"""Refractive index of air"""
const N_AIR = 1.00029

"""Surface roughness (nm)"""
const ROUGHNESS = 0.5
