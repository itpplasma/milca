# Equilibrium roadmap

One ideal-MHD model with two discretisations that are inverse to each
other, and separate physical models where flux surfaces do not exist. This
document is the plan; results live in the lane documents it links.

| Layer | Content |
|---|---|
| Shared infrastructure | spline de Rham spaces, Winslow meshes, coil fields, vacuum (BEM, exterior map), Newton with exact derivatives, gauge constraints, continuation, verification |
| Ideal model (surfaces exact) | the Cartesian flux map u: space → flux disc ([flux-map.md](flux-map.md)) and the Lagrangian map X: flux disc → space (`nested_equilibrium`, [ideal-map-3d.md](ideal-map-3d.md)); which one is the core is decided in P3 |
| Models without imposed surfaces | Grad–Shafranov (axisymmetric, X-points), NN1 (vacuum), Beltrami volumes, NN2 with transport closure ([NN2](equilibrium-nonnested-nn2.md)); used at separatrices, X-lines, open field lines and in chaotic volumes ([flux-map.md](flux-map.md), 5a) |
| Interface | one physical state and diagnostics: reconstruction, ideal stability, linear response, heat loads |

NN2 is not an alternative solver for the ideal equations: the physical
models stay distinct.

## Standing rules

- From the [audit](equilibrium-lane-audit.md) of 2026-11-22: a cause is
  written into a plan only with the run that isolates it, otherwise as a
  hypothesis; gates include a realistic B_φ/B_pol (≥ 10), a low-shear case
  and degree ≥ 3; null directions are removed by constraints from the
  first version.
- Mesh coordinates are never physical unknowns and are never added to cure
  conditioning. A mesh moves by its own equation (Winslow, alignment,
  refinement) between physics solves, and only where the physics has an
  interface (free boundary, region interface).
- No new representation is designed before P3 has its result.
- Record of what was tried and why it was left: the audit; the (χ, g)
  label solvers with axis unknowns in
  [eulerian-label-plan.md](eulerian-label-plan.md) (superseded).

## Way forward (state of 2026-11-22)

The flux map u: space → flux disc and the Lagrangian map X: flux disc →
space discretise one variational principle (the energy with the frozen flux
form) in opposite directions; exact solutions are inverse to each other up
to the gauge θ* → θ* + k(s). They are one ideal model with two
discretisations that check each other, not two development lines. No new
representation is designed before P3 has its result: a defect found on
the way is repaired where it is measured, or the existing Lagrangian map
is used.

Amended after an external review of the first version (same day): gauge
constraints and continuation moved before the comparison (P0c); checks of
the unprojected force balance and of the sign of the toroidal field added
(P0a); basins added to P1; the two kinds of interface distinguished in P4;
the triangular figure-eight as the final stress test (P6).

**P0. The flux map in 3D: state and what is still owed.**

Done ([flux-map.md](flux-map.md), section 4; #645): `milca_flux_map_3d`,
gates `flux_map_3d_gate_f2a` and `_f2b`, QA driver `test_flux_map_qa`.
- F2a: the F1 Solov'ev numbers through the 3D path, to three digits.
- F2b (sheared Landreman family, degree 2, start without axis
  information): B 3.49e-2, 7.43e-3 on 2, 4 cells (rate 2.23; the label
  lane had 4.80e-2, 1.159e-2), 11 Newton steps each.
- QA (Landreman–Paul, one period, pointwise field, splines): converges at
  every resolution tried; at 3,536 unknowns energy +0.40 %, |B| on the
  s = 0.25 and 0.5 surfaces 0.6–0.8 %, axis 1.5 mm, 22 s; at 14,416
  unknowns energy +0.06 %, |B| 0.2–0.4 %, 99 s. VMEC++: 9 s.

Not at parity, three findings, none explained by a run yet:
1. Beyond 4 cells (and at degree 3, and on 21 planes) the energy
   minimiser leaves the exact Landreman state for one of lower energy with
   the axis displaced; Newton on stationarity reaches the exact state at
   8 cells (B 2.48e-3). Hypothesis: the exact state is an unstable
   stationary point; the Lagrangian map and VMEC++ show an unstable mode on
   this family ([audit](equilibrium-lane-audit.md), section 1). No Hessian
   eigenvalue has been computed.
2. On QA the axis (1.3–1.8 mm) and the s = 0.25 surface (0.6–0.9 mm) do
   not improve with resolution while energy and |B| do. Candidates, none
   measured: the near-shearless ι, the section meshes (worst Jacobian
   ratio 0.04–0.06), the toroidal resolution of the wall.
3. No gain from the degree on QA (energy error 0.40, 0.30, 0.20 % at
   degrees 2, 3, 4 on 4 cells against 0.06 % at degree 2 on 8 cells).
   Suspects: the spline fit of the wall, the C0 block corners.

State of 2026-11-22 after the harmonic section meshes
([flux-map.md](flux-map.md), "Section meshes by harmonic maps"): finding 1
is a pressure-driven saddle (none of the standard decks with pressure has a
negative eigenvalue); finding 2 was an error of the flux map, not of VMEC
(DESC agrees with VMEC's limit), and the QA axis converges to that limit at
degree 3, the QH axis once 29 planes are used; the energy error of finding
3 was the area of the wall fit and is removed. W7-X standard and li383 mesh
and converge. Open from these: grading of the meshes towards tips, QH
against a third code, the decision on P0b.

P0a, owed before the flux map is compared (in this order):
- QA floor: plane scan at fixed section (13, 17, 21, 25), wall volume
  against VMEC, then section meshes by Winslow.
- F2b: lowest Hessian eigenvalues at the exact and at the lower-energy
  state; then either a stable exact 3D case for the minimiser or the gate
  in stationary mode on three meshes.
- Convergence on three mesh levels at degrees 2, 3, 4 of (i) B against
  the exact field and (ii) the residual of the unprojected force balance
  J×B − ∇p, evaluated independently of the energy that is minimised.
- At convergence: the sign of ∇ζ·(∇x×∇y) at all quadrature points and
  the winding of the wall angle, reported per run. Not a step-acceptance
  rule: with a projected energy a pointwise guard stalled the prototype.
- The order 1.8–1.9 of F1 (twice the projected exact label on 8 cells).

P0b, projection against guide-field locking: not on the critical path. In
F1 the pointwise solver at B_φ/B_pol ≈ 15 is at 2.2e-3 against 1.9e-3 at
1.5; the damage is at 150. No stable 3D projection exists (all per-cell
variants collapse in the prototype) and a projected field is neither
tangent nor divergence-free without a commuting construction. It is taken
up only if P1 shows amplification at realistic ratios, and then with the
force-balance check of P0a as its acceptance.

P0c, shared numerics first (was P5), because a comparison of solvers with
penalised gauges measures the penalties:
- S1 gauge by constraints: one module for moment constraints orthogonal
  to a continuous symmetry (k(s) of the flux map, the slide of the
  Lagrangian map, the closed-line gauge of Φ in NN2). Gate: F1 and F2b on
  meshes without symmetry and without the symmetry rows.
- S2 stiff parameters: one continuation driver and one step acceptance
  that tolerates a curved stiff constraint. Gates: F1 from a cold start at
  B_φ/B_pol = 150 without hand-set stages; the NN2 ratio 1e6 case from the
  seed.
- One solver loop for F1 and F2 (duplicated at present).

**P1. Head-to-head, flux map against the Lagrangian map**, VMEC++ as the
outside reference, on identical cases and inputs:

| Case | Reference | Why |
|---|---|---|
| H1 a 3D case with a stable exact state (or the sheared Landreman family in stationary mode) | exact field | absolute error |
| H2 Landreman–Paul QA, prescribed ι (`benchmarks/label-qa`) | VMEC++ at ns 101 | low shear, B_φ/B_pol ≈ 14 |
| H3 one strongly shaped deck (QI nfp 2 of nested-equilibrium-input.md, or the CTH-like deck) | VMEC++ | where the map lane folds or stalls |

Measured for each code at degrees 3 and 4 and at least three resolutions,
compared at equal physical error, not equal grid: unknowns, B error,
surface distance, axis, energy, Newton steps, wall time at 1 and 8
threads, peak memory, the convergence basin (perturbed starts: axis
offset, boundary-angle perturbation), the sign of the Jacobian, and the
failure mode where one fails; and the implicit sensitivities (boundary
and profile parameters) of each code against re-solved finite differences.
One table in `benchmarks/`.

**P2. Inverse check (gate V1).** With X_h from the Lagrangian map and u_h
from the flux map on the same case: s(u_h(X_h(q))) − s(q), free of the
gauge, and the angle part after removing its mean per surface. First on
the axisymmetric Solov'ev case, then on H1. Pass: both at the order of the
coarser discretisation.

**P3. Decision on the ideal core**, recorded here with the tables of P1
and P2: the map with the better robustness and time at matched error on
H1–H3 becomes the core; the other is kept as cross-check with parity
fixes. Until then neither is "reference only".

**P4. Region coupling where the answer is known (gates R1).** Two kinds
of interface, to be kept apart:
- numerical: two descriptions of the same ideal equilibrium meet at a good
  surface. Conditions: B·n = 0 from both sides, continuous pressure,
  continuous tangential B (no sheet current), matching fluxes;
- physical: two closures meet (ideal core against NN2, Beltrami or
  vacuum). Conditions from the variation of the shared interface: balance
  of p + B²/2, a sheet current where the model admits one, a
  transport-flux law towards NN2.
Axisymmetric first, numerical interface, ideal core map inside s = s_i and
ψ = R A_φ with f prescribed outside, so that the exact answer is
Grad–Shafranov:
- R1a: Solov'ev split at an interior flux surface; pass: the one-region
  F1 errors, interface at the exact surface.
- R1b: the exact double-null Solov'ev of the Grad–Shafranov lane, the
  separatrix and X-points in the outer region; pass: the errors of the
  one-region Grad–Shafranov solve.
- R1c: s_i = 0.8, 0.9, 0.95: no dependence beyond the discretisation
  error.
Then 3D with NN1 (vacuum) outside, then the physical interface to NN2
([flux-map.md](flux-map.md), 5a).

**P5. Difficult geometry, only after the core has passed P1–P3.** F3:
mesh swept along a guiding curve with a regular frame and a check that
the swept domain does not overlap itself. Strongly shaped QI and rounded
triangular sections need, separately from F3, local resolution at the
corners; the 3D spaces are global in the toroidal direction at present
and a ridge that rotates with the section may need a local basis there.

**P6. Stress test and beyond.** A triangular figure-eight: scan of the
corner rounding against GVEC and DESC on the same smooth equilibrium (an
exact cusp is a separate limit problem). Free boundary on the core chosen
in P3 (the exterior operator and shape derivatives exist); regions in 3D
(island interiors as their own maps inside good surfaces, the separatrix
layer in A); NN2 on the shared mesh.

**D1. Downstream: differentiable free-boundary equilibrium for
reconstruction and optimisation.** The gain (implicit tangent and adjoint
from the factorised Newton matrix instead of re-solved equilibria per
parameter) does not depend on the direction of the map. It is attached to
whichever map has a working free boundary: today the Lagrangian map
(axisymmetric free boundary at VMEC accuracy, exterior operator with
analytic derivatives); the flux map has none. First gate, independent of
P0–P3: the exact axisymmetric free-boundary Solov'ev case with coil
currents and profile parameters, sensitivities against re-solved finite
differences; then a provider for TIAGO behind the boundary of
[equilibria-inference.md](equilibria-inference.md). No optimiser is built
here. Fixed-boundary shape optimisation favours the Lagrangian map (the
boundary is its parameter; for the flux map it changes the domain), and
quasi-symmetry measures need Boozer angles on top of θ* in either map.
Whether TIAGO's present free-boundary derivatives are finite differences
was stated in a review and has not been checked here.

**D2. Downstream exploratory: exact-force-balance solution families and
Boozer/Hamada consistency.** See the scoped
[research plan](exact-force-balance-coordinate-plan.md). This asks a different
question from P0--P3: can an analytic/nonlinear representation satisfy
interior ideal-MHD force balance *by identity* while the fit error lies in
boundary, profile, coil or diagnostic data? Separate (E0) exact restricted
Grad--Shafranov/Trefftz families, (E1) published exact finite-pressure 3D
families, (E2) a high-risk divergence-free stress potential with **all**
nonlinear MHD constitutive constraints, and (E3) a nonlinear manifold
corrector which provides only **discrete** equilibrium accuracy unless
independently certified. No new representation displaces either map at P3.
Use [libneo's independent Boozer and Hamada routes](https://github.com/itpplasma/libneo/blob/47b1988cddf5e6544d92eb7418c87f83be3d7ea4/src/boozer_hybrid.f90)
as coordinate-compatibility oracles (C0--C2), always together with direct
force-balance checks. Do not claim universal machine precision near rational
surfaces or from current-based Hamada coordinates at p' = 0. The smallest
experimental slice is exact axisymmetric Solov'ev boundary/profile fitting
and coordinate-route agreement; the experimental track is **not active**
until it is explicitly admitted under [ROADMAP.md](../ROADMAP.md).

**Dropped:** further work on the (χ, g) solvers with axis unknowns (the
gauge work of branch `label/qa-recipe` stays unlanded), the frozen-map
relabelling (B1–B6, R1), meshes that follow the magnetic axis, gate C as
formulated.

**Not run so far in the flux map:** any case with pressure in 3D, QH and
QI decks, a non-symmetric case, free boundary, any region coupling. Speed:
the Hessian couples all planes of a section cell and is factorised
directly; a sparse toroidal basis or a Krylov solver is owed before a
time comparison means anything.

Effort split (equilibrium development), unchanged: common infrastructure
35 %, ideal model and parity with VMEC 35 %, NN2 and island benchmarks
20 %, stability and reconstruction interfaces 10 %. Until P3 the ideal
share goes to P0–P2 and the infrastructure share to P0c.

## Optional, not an immediate priority

- Harmonic boundary-conforming maps with proven invertibility for
  nonconvex sections: Babin, Hindenlang, Maj and Köberl, *Construction of an
  invertible mapping to boundary conforming coordinates for arbitrarily
  shaped toroidal domains*, arXiv:2411.04683, code `map2disc`
  (gitlab.mpcdf.mpg.de/gvec-group/map2disc). They solve Δξ = Δη = 0 (the
  Winslow map) by boundary integrals rather than the nonlinear inverse
  equations, then invert it. The construction is implemented from the
  paper in FortFEM (`fortfem_harmonic_disc`) and is the section mesh of the
  flux map; their code was not read or used.

## State of the Lagrangian map lane

- Axisymmetric free boundary at VMEC speed on the exact Solov'ev benchmark
  ([table](../benchmarks/free-boundary-exact/README.md)).
- 3D free boundary on a curve chart (`fb3d_synthetic`, R2 exterior): 44 s
  against 67 s for VMEC++. (An earlier "under-resolved at ntor 2" here was
  wrong: ntor 3 and 4 change the error by 0.25 mm; the isolated cause is
  the coil ripple of the oracle, free-boundary-nested-map-3d.md, 9.3.)
- CTH-like decks fold the shell of the cylindrical chart (1.5 cm margin
  against a 4.6 cm helical swing). `milca_deck_chart` builds the
  axis-following chart from any deck; wiring it into `nested_equilibrium`
  is not done and is needed for H3.
