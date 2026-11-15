# Electron-positron plasmas

Status: designed, 2026-11-15. Parts already run: Bhabha scattering, pair
annihilation and the thermal annihilation rate of
[QED and radiation](qed-radiation.md).

## Why

Electron-positron plasmas have two species of equal mass and opposite charge.
The APEX collaboration at IPP Garching and TU Munich prepares magnetically
confined pair plasmas at low temperature, first in a levitated dipole and then
in a stellarator (Stoneking et al. 2020), after injecting positrons from the
NEPOMUC source into a dipole trap without loss (Stenson et al. 2018). The mass
symmetry gives such plasmas gyrokinetic stability properties that
electron-ion plasmas lack, broken by curvature drifts in closed field-line
geometry (Helander 2014; Helander and Connor 2016; Kennedy et al. 2020).
Relativistic pair plasmas near compact objects belong to the
[QED branch](qed-radiation.md).

## What MILCA needs

- **Species.** Positrons as a species of charge +e and electron mass. No model
  may assume an ion-electron mass ratio; drift-kinetic and gyrokinetic models
  keep finite-Larmor-radius effects for both species alike.
- **Collisions.** Coulomb collisions of equal masses through the existing
  [collision hierarchy](collisions.md); at relativistic energies the Bhabha and
  Møller kernels of Q1.
- **Losses.** Direct annihilation (running for Juttner distributions) and,
  dominant at the electronvolt temperatures of the experiments, positronium
  formation by radiative capture and by charge exchange with neutrals. Radiative
  positronium formation is hydrogen-like with reduced mass m/2 and is checked
  against the scaled hydrogen rates.
- **Geometry.** Dipole and stellarator vacuum fields from the existing
  equilibrium lanes.
- **Stability.** Linear drift-kinetic and gyrokinetic stability in slab,
  dipole and stellarator geometry, checked against the analytic results cited
  above once the gyrokinetic branch exists.
- **Diagnostics.** The annihilation gamma-ray signal as a synthetic diagnostic
  from the annihilation and positronium rates.

## Overlap checks

- Charge conjugation: a symmetric pair plasma is invariant under exchanging
  the species; every pair-plasma result is checked against that symmetry.
- Equal-mass limit: the electron-ion models with the ion mass set to the
  electron mass and the charge reversed reproduce the pair results.
- Low energy: Bhabha scattering reduces to Coulomb scattering with exchange,
  and Dirac annihilation to π r_e² c (running, `qed_tree_level`, `qed_thermal`).

## Order

After P1g: positronium formation rates and the positron species in the
kinetic models, then linear stability in dipole geometry against Helander
(2014) and Kennedy et al. (2020).

## Sources

- P. Helander, Phys. Rev. Lett. 113, 135003 (2014).
- P. Helander and J. W. Connor, J. Plasma Phys. 82 (2016),
  doi:10.1017/S0022377816000490.
- E. V. Stenson et al., Phys. Rev. Lett. 121, 235005 (2018).
- D. Kennedy et al., J. Plasma Phys. 86 (2020), doi:10.1017/S0022377820000276.
- M. R. Stoneking et al., J. Plasma Phys. 86 (2020),
  doi:10.1017/S0022377820001385.

Chris&AI
