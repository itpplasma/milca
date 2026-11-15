# Nuclear reactions

Status: designed, 2026-11-15; the [proton-proton](pp-fusion.md) and
[pionless EFT](pp-fusion-eft.md) calculations run, and N3 runs: Maxwellian
reactivities from S factors (`test_fusion_reactivity`) for D-T (Bosch and Hale)
and p-¹¹B (Wang, Li, Wu and Cui, arXiv:2601.00241, 2026); N4 runs for beams on
Maxwellian targets (`test_fusion_beam_target`). Fusion reactions enter MILCA
at every level from few-body nuclear theory to reaction-rate kernels for the
plasma distributions MILCA computes.

## Hierarchy

**N0. Few-body nuclear theory.** Bound and scattering states of the compound
systems of light-ion fusion (D-T, D-³He, D-D) from nucleon-nucleon and
three-nucleon interactions of chiral or pionless effective field theory, by
few-body methods (Faddeev-Yakubovsky, hyperspherical harmonics, the no-core
shell model with continuum). QCD enters only through the low-energy constants
of these theories; a lattice-QCD calculation of these cross sections is out of
reach and not planned. p-¹¹B, with twelve nucleons and a three-alpha final
state, is beyond few-body exactness and enters at N1 and N2.

**N1. Potential and cluster models.** Two-body potentials with Coulomb
distortion, as in the pp calculation, and cluster models such as three-alpha
descriptions of ¹²C states for p-¹¹B.

**N2. R-matrix.** Resonance parameters fitted to measured cross sections
(Descouvemont and Baye 2010), giving the astrophysical S factor S(E) with its
uncertainty.

**N3. S factor and Gamow factor.** σ(E) = S(E)/E · exp(-(E_G/E)^½), with E
the centre-of-mass energy and E_G the Gamow energy, and the Maxwellian
reactivity, checked against Bosch and Hale (1992) for D-T, D-³He
and D-D, and against Nevins and Swain (2000) and Sikora and Weller (2016) for
p-¹¹B.

**N4. Rate kernels for non-Maxwellian plasmas.** Reaction rates and product
spectra for the distributions MILCA computes: beam-target, anisotropic and
fast-ion tails; birth spectra of alphas and of the three alphas of p-¹¹B as
kinetic sources. Plasma screening of the reaction (Salpeter 1954) connects to
the [quantum plasma](quantum-plasma.md) branch.

**N5. Nuclear elastic scattering.** Large-angle nuclear scattering of MeV ions
as a further collision operator next to the Coulomb operator. Its weight in
fast-ion slowing down depends on fuel and temperature and is computed, not
assumed; it grows with the ion temperatures of advanced fuels.

## Order

After the first radiation slice (ROADMAP P1g): N3 and N4 for D-T, D-³He and
p-¹¹B with published S factors, so that the p-¹¹B power balance of
[QED and radiation](qed-radiation.md) uses MILCA's own rate kernels; then N2,
N5, the extension of N1 from pp to D-T, and N0.

## Sources

- E. E. Salpeter, Aust. J. Phys. 7, 373 (1954).
- H.-S. Bosch and G. M. Hale, Nucl. Fusion 32, 611 (1992).
- W. M. Nevins and R. Swain, Nucl. Fusion 40, 865 (2000).
- P. Descouvemont and D. Baye, Rep. Prog. Phys. 73, 036301 (2010).
- M. H. Sikora and H. R. Weller, J. Fusion Energy 35, 538 (2016).

Chris&AI
