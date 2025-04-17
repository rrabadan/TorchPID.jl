"""
Physical constants used throughout the TORCH detector simulation.

Constants are defined with specific units relevant to the simulation:
- Distances in mm
- Times in ns
- Energies in GeV or eV as specified
"""

"""Speed of light in vacuum (mm/ns): 299.792458"""
const CLIGHT = 299.792458

"""Elementary charge (C): 1.60217662e-19"""
const QQe = 1.60217662e-19

"""Pion mass (GeV/c^2): 0.139570"""
const PION_MASS = 0.139570

"""Kaon mass (GeV/c^2): 0.493677"""
const KAON_MASS = 0.493677

"""Proton mass (GeV/c^2): 0.938272"""
const PROTON_MASS = 0.938272

"""Electron mass (GeV/c^2): 0.000511"""
const ELECTRON_MASS = 0.000511

"""Muon mass (GeV/c^2): 0.105658"""
const MUON_MASS = 0.105658

"""Energy to wavelength conversion factor (nm * eV): 1239.84193"""
const LAMBDA = 1239.84193

"""Refractive index of air: 1.00029"""
const N_AIR = 1.00029

"""Surface roughness (nm): 0.5"""
const ROUGHNESS = 0.5
