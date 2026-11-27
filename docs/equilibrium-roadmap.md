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
and converge. Open from these: QH against a third code (the grading of the
meshes towards tips is measured under P0a). Decided 2026-11-22: P0b is dropped; the F2b gate in
stationary mode (level 1 from the polar start, no axis check) is accepted.

P0a, owed before the flux map is compared: done 2026-11-23, numbers in
[flux-map.md](flux-map.md), "Force-balance residual" and "Mesh grading
towards tips".
- QA floor: done. The plane scan is recorded, and the harmonic section
  meshes replaced the Winslow meshes that were planned here.
- F2b: done; the gate runs in stationary mode (accepted 2026-11-22).
- Convergence on three mesh levels at degrees 2, 3, 4: done on the
  Solov'ev case, B at order p (1.98, 2.93, 4.14) and the unprojected force
  balance J×B − ∇p, J = ∇×B_h from the second derivatives of the map
  (`milca_flux_map_force`), at order p − 1 (0.98, 2.03, 2.76). On the
  sheared Landreman family both stop beyond 4 cells at every degree (B
  2.0e-3 to 2.5e-3, residual 0.07 to 0.10 of |∇p|); planes or the saddle,
  not isolated.
- The same residual on the standard decks and for the VMEC++ states: the
  degree-2 states have 0.5 to 0.8 of |∇(B²/2)| where VMEC++ has 2e-3 to
  3.5e-2 (outside s = 0.1). It grows with the guide field at fixed error
  of B (isolated on Solov'ev), and on QA the layer s > 0.9 does not
  converge up to degree 3, 8 cells, 21 planes. Open, and it decides the
  standing benchmark below.
- Sign of ∇ζ·(∇x×∇y) at all quadrature points and winding of the wall
  angle: reported per run and one sign, winding one in every converged
  run; the gates check them.
- The order 1.8–1.9 of F1: pre-asymptotic. 1.97, 1.98 on 16, 32 cells;
  the factor against the projected exact label levels off at 2.3 to 2.4.
- Mesh grading towards tips (FortFEM `harmonic_disc_build`, `grading`):
  the W7-X wall distance at 4 cells falls 3.4 to 4.3 times, the axis gets
  3 to 4 times worse on 13 planes and is no better on 29. Off by default.

P0b, projection against guide-field locking: dropped (2026-11-22). The QA
axis converges to the VMEC/DESC limit with harmonic meshes at degree 3 and
the QH axis approaches VMEC once enough planes are used, without a
projection. No stable 3D projection exists (all per-cell variants collapse
in the prototype). Reopened only if P1 shows amplification at realistic
ratios B_φ/B_pol, then with the force-balance check of P0a as acceptance.

P0c, shared numerics first (was P5), because a comparison of solvers with
penalised gauges measures the penalties: done 2026-11-23 for the flux map,
numbers and the hypotheses with their outcomes in
[flux-map.md](flux-map.md), "Gauge constraint, one solver loop and
continuation"; generic parts in FortFEM (`fortfem_constrained_newton`,
`fortfem_parameter_continuation`, doc/CONSTRAINED_NEWTON.md there).
- S1 gauge by constraints: `milca_gauge_moments`, one module for moment
  constraints orthogonal to a continuous symmetry, with dependent
  generators dropped. The flux map uses it for the rotation of the disc,
  the one exact symmetry of its discrete energy; the penalty and symmetry
  rows are gone. Gate passed: F1 and F2b on meshes without symmetry and
  without symmetry rows, same errors as before. Ruled out: constraining
  k(s) beyond the rotation (it moves the solution by 0.5 to 33 %). Not
  done, each is more than a small change and is the next step of its lane:
  the slide of the Lagrangian map (removed there by reducing the axis
  basis, `nested_axis_basis`), the closed-line gauge of Φ in NN2 (proposed
  in the NN2 document, no code yet).
- S2 stiff parameters: one continuation driver, used by the flux map
  (guide field) and by NN2 (`nc_continue`, anisotropy). Gates passed: F1
  from a cold start at B_φ/B_pol = 150, stages chosen by the driver; NN2
  P2 at ratio 1e6 from the shifted seed (nr 4 and 8, N = 2) with the
  stages, steps and path differences it had. Step acceptance: three rules
  meant to tolerate the stiff constraint were tested and ruled out for the
  flux map (the Hessian is indefinite far from the solution; that, not the
  curvature of the constraint, refuses the steps); relaxed steps exist in
  FortFEM and are off here (fewer steps on F1 without symmetry, more on
  the decks). The NN2 Newton iteration is not on the shared loop.
- One solver loop for F1 and F2: done (`flux_map_newton`), Newton on the
  KKT system by LDLᵀ with the inertia on the constraints.
- Found on the way: a mesh without symmetry costs the minimiser up to
  three times the Newton steps (soft k(s) directions). Open.

**P1. Head-to-head, flux map against the Lagrangian map**, VMEC++ as the
outside reference, on identical cases and inputs:

| Case | Reference | Why |
|---|---|---|
| H1 a 3D case with a stable exact state (or the sheared Landreman family in stationary mode) | exact field | absolute error |
| H2 Landreman–Paul QA, prescribed ι (`benchmarks/label-qa`) | VMEC++ at ns 101 | low shear, B_φ/B_pol ≈ 14 |
| H3 one strongly shaped deck (QI nfp 2 of nested-equilibrium-input.md, or the CTH-like deck) | VMEC++ | where the map lane folds or stalls |

The integer-transform Landreman equilibrium is a fixed-boundary saddle even
within stellarator symmetry. Its flux and pressure labels are weakly
constrained by magnetic-field error. Use it as a stationary-solver test and
track normalized flux and angle as well as the field.
[Exact saddle and discrete curvature](flux-map.md#integer-transform-landreman-saddle).

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

**Later, not scheduled: both directions as modes of one stack** (recorded
2026-11-23 so it is not lost). Direction (VMEC mode X: flux disc → space,
the Lagrangian map; CEMV mode u: space → flux disc, the flux map) and
section discretisation (splines on the harmonic disc O-grid, curved
triangles P_k) are independent axes; today two of the four combinations
exist, in separate code. Target layering with duplication only where the
physics differs:
- one FortSym energy leaf per direction (the pointwise density of the
  map's jet with its generated derivatives): w(u, ∇u, R) for CEMV,
  |DX b|²/(2J) − pJ for the VMEC mode;
- basis jets per space from FortFEM (spline O-grid or triangles, times a
  trigonometric nodal or Fourier toroidal basis), which exist;
- one assembly over cells and quadrature points, one Newton with LDLᵀ and
  inertia, gauge by constraints and continuation (P0c);
- shared geometry and checks: charts (cylindrical, curve chart), harmonic
  meshes with their certificate, wall condition, free boundary through the
  exterior operator, VMEC import and state check, the inverse check of P2.
Features of either lane are wanted in both: curve chart, 3D free boundary
and matrix-free preconditioners of the Lagrangian map; harmonic meshes,
LDLᵀ inertia, exact wall condition and VMEC state check of CEMV. Not
measured, to be settled before any code: whether CEMV on C0 triangles
locks at a strong guide field (as C0 Lagrange quadrilaterals did, section
3 of flux-map.md), whether the VMEC mode keeps its orders on the spline
O-grid, how much of the free boundary depends on triangles, and whether
one assembly generic over spaces stays as fast as the specialised paths.
Working rule until then: when one lane hits a problem, look at how the
other solves it and unify that piece instead of writing a second
solution. A standing benchmark CEMV vs VMEC (both converged below their
difference, same independent force-balance residual, time at matched
accuracy including setup) belongs with P1.

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

**Not run so far in the flux map:** QI decks, an equilibrium without
symmetry (the meshes of gates F1 and F2b have none since P0c), free
boundary, any region coupling. Speed:
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
  axis-following chart from any deck; `nested_equilibrium` takes it with
  `chart = 'deck'` in fixed boundary. The free-boundary shell on it is open.
