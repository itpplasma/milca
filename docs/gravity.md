# Gravity, general relativity and dark matter

Status: designed, 2026-11-15. Gravity and general relativity are inside MILCA's
scope from this date; the capability catalogue had left them to external
codes.

## Why

Self-gravitating systems are many-body systems with the structure of a
Coulomb system whose interaction has the opposite sign: like masses attract.
MILCA's many-body references, Vlasov and collision operators, and
fast-multipole and particle-mesh summation carry over with the sign of the
coupling reversed, while the physics differs: there is no screening, the
specific heat of a bound system is negative, and there is no true thermal
equilibrium. A separate inverse code, in the way TIAGO (our reconstruction
and inference code for fusion devices) consumes MILCA, will later
use MILCA as the differentiable forward model for Bayesian comparison of
dark-matter and gravity hypotheses with all available observations. That code
is not part of MILCA and is not started; MILCA provides the forward models with
the residual contract of [DESIGN.md](../DESIGN.md) and its Jacobian-vector and
vector-Jacobian products (forward and adjoint derivatives).

## Plasma and gravity side by side

| Plasma | Gravity |
|---|---|
| Coulomb interaction, both signs | Newtonian interaction, attractive only |
| Debye screening | no screening; Jeans instability at the Jeans length |
| Landau and Rosenbluth operators | dynamical friction (Chandrasekhar 1943) |
| Lenard-Balescu operator | its gravitational counterpart (Heyvaerts 2010) |
| Darwin Lagrangian, order (v/c)² | Einstein-Infeld-Hoffmann equations (1938) |
| Maxwell equations | linearized gravity, gravitomagnetism |
| collisionless relaxation | violent relaxation (Lynden-Bell 1967) |

## Hierarchy

**G0. General relativity with matter.** The 3+1 Einstein equations
(Baumgarte and Shapiro 2010) coupled to Vlasov, MHD (GRMHD) and
electromagnetic fields, and the cosmological background with its
perturbations.

**G1. Post-Newtonian many-body dynamics.** The first post-Newtonian
Einstein-Infeld-Hoffmann equations and higher orders with radiation reaction.

**G2. Newtonian many-body dynamics.** Point masses summed by the fast
multipole method (Greengard and Rokhlin 1987) and particle-mesh methods shared
with the Coulomb many-body code; softening is a declared model choice, as for
point charges. Point masses run on the shared N-body Hamiltonian with the
coupling -G and the masses as charges (`milca_self_gravity`); a Plummer
sphere sampled from its distribution function reproduces its kinetic and
potential energies 3π/64 and -3π/32 (G = M = a = 1) and stays in virial
equilibrium (`test_plummer`).

**G3. Vlasov-Poisson and kinetic models.** Collisionless dynamics in physical
and comoving coordinates, and collision operators for two-body relaxation in
star clusters. The linear kinetic response of a homogeneous self-gravitating
Maxwellian is the plasma response with the coupling reversed; its purely
growing Jeans root solves k² = 1 − √π y e^{y²} erfc(y), y = γ/(√2 k), in units
of ω_J = √(4πGρ) and the Jeans length v_th/ω_J, and runs in
`jeans_root` and `vlasov_density_response(..., coupling=-1)`
(`test_jeans`).

**G4. Self-gravitating fluids.** Hydrodynamics and MHD with gravity, with
cooling from the [radiation](radiation.md) and [QED](qed-radiation.md)
branches; star formation and feedback enter as declared subgrid models, never
as first principles.

## Dark-matter and gravity hypotheses

- Cold dark matter: collisionless Vlasov-Poisson and N-body (G2, G3).
- Warm dark matter: thermal velocities and free streaming in the initial
  state.
- Self-interacting dark matter: a collision operator from the scattering
  cross section (Tulin and Yu 2018), on the plasma collision machinery.
- Fuzzy dark matter: the Schrödinger-Poisson equations of an ultralight scalar
  field (Hu, Barkana and Gruzinov 2000; Schive, Chiueh and Broadhurst 2014),
  which reuse the [quantum branch](quantum-plasma.md) and reduce to
  Vlasov-Poisson through the Wigner representation (Widrow and Kaiser 1993).
  The solitonic core, the Schrödinger-Poisson ground state, runs
  (`milca_fuzzy_dark_matter`): it satisfies the virial theorem and the
  eigenvalue relation to 1e-12, its scaling symmetry, and r_1/2 M = 3.925
  ħ²/(G m²) of the literature.
- Modified gravity: MOND (Milgrom 1983) as an alternative to dark matter, not
  a reduction of general relativity. In its quasi-linear formulation QUMOND
  (Milgrom 2010) it is two linear Poisson problems with an algebraic step
  between them; in AQUAL (Bekenstein and Milgrom 1984) one nonlinear Poisson
  problem, both on the finite-element machinery.

## Speed

The target is computation per core and time to a named observable at matched
accuracy, against PKDGRAV3 (Potter, Stadel and Teyssier 2017), SWIFT (Schaller
et al. 2024), GADGET-4 (Springel et al. 2021) and AREPO (Springel 2010) for
N-body and hydrodynamics, against the Einstein Toolkit (Löffler et al. 2012)
for dynamical spacetimes, and against Athena++ (Stone et al. 2020) for GRMHD on
a fixed metric.

## Order

1. G2 and G3: a self-gravitating many-body slice on the shared fast-multipole
   code with the Plummer sphere as equilibrium, the Jeans growth rate, cold
   collapse and core collapse of a star cluster as checks.
2. The dark-matter hypotheses on that slice: fuzzy-dark-matter soliton cores,
   self-interaction, warm dark matter, QUMOND.
3. G1, then the comoving cosmological box with G4, then G0.

## Sources

- S. Chandrasekhar, Astrophys. J. 97, 255 (1943).
- A. Einstein, L. Infeld and B. Hoffmann, Ann. Math. 39, 65 (1938).
- D. Lynden-Bell, Mon. Not. R. Astron. Soc. 136, 101 (1967).
- M. Milgrom, Astrophys. J. 270, 365 (1983); Mon. Not. R. Astron. Soc. 403,
  886 (2010).
- J. D. Bekenstein and M. Milgrom, Astrophys. J. 286, 7 (1984).
- L. Greengard and V. Rokhlin, J. Comput. Phys. 73, 325 (1987).
- L. M. Widrow and N. Kaiser, Astrophys. J. 416, L71 (1993).
- W. Hu, R. Barkana and A. Gruzinov, Phys. Rev. Lett. 85, 1158 (2000).
- J. Heyvaerts, Mon. Not. R. Astron. Soc. 407, 355 (2010).
- T. W. Baumgarte and S. L. Shapiro, *Numerical Relativity* (Cambridge, 2010).
- V. Springel, Mon. Not. R. Astron. Soc. 401, 791 (2010).
- F. Löffler et al., Class. Quantum Grav. 29, 115001 (2012).
- H.-Y. Schive, T. Chiueh and T. Broadhurst, Nature Phys. 10, 496 (2014).
- D. Potter, J. Stadel and R. Teyssier, Comput. Astrophys. Cosmol. 4, 2 (2017).
- S. Tulin and H.-B. Yu, Phys. Rep. 730, 1 (2018).
- J. M. Stone et al., Astrophys. J. Suppl. 249, 4 (2020).
- V. Springel et al., Mon. Not. R. Astron. Soc. 506, 2871 (2021).
- M. Schaller et al., Mon. Not. R. Astron. Soc. 530, 2378 (2024).

Chris&AI
