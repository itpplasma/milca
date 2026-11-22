# Equilibrium: the way forward (plan of 2026-11-22)

One Eulerian computational infrastructure, two physical equilibrium models:
ideal MHD with flux labels (the label lane) and NN2 without imposed flux
surfaces. The Lagrangian map, Eulerian Clebsch, ALE–Winslow and NN2 are not
four competing solvers for one problem; they serve different purposes.

| Layer | Content |
|---|---|
| Shared infrastructure | spline de Rham spaces, ALE–Winslow mesh, coil fields, vacuum (BEM, exterior map), Newton–Krylov, verification |
| Ideal lane (labels) | Eulerian Clebsch, A = g∇χ + Σ F_j(χ) h_j + G(χ)∇φ ([plan](eulerian-label-plan.md)); nested surfaces, islands of prescribed topology, current sheets |
| NN2 lane | general A, p, Φ, U ([NN2](equilibrium-nonnested-nn2.md)); islands that form, stochastic regions, resistive steady states, transport |
| Interface | one physical state and diagnostics: reconstruction, ideal stability, linear response, heat loads |

NN2 is not an alternative nonlinear solver for the ideal equations: the
physical models stay distinct.

## Decisions

| Component | Decision | Reason |
|---|---|---|
| Lagrangian map (`nested_equilibrium`, [ideal-map-3d.md](ideal-map-3d.md)) | reference, no new development beyond parity fixes | verified, exact energy Hessian, ideal stability; its relabelling null needs gauges that fail in strong shaping |
| Eulerian Clebsch labels | develop: the ideal lane | exact tangency B·∇χ = 0, no multiplier, no locking (E0–E6c) |
| Constrained A ∈ V1 with B·∇χ = 0 weakly ([eulerian-3d.md](eulerian-3d.md)) | alternative, decided by gate C below | more flexible, but the weak constraint locks |
| ALE–Winslow | shared infrastructure | mesh quality independent of the physical unknowns |
| NN2 | continue on the shared mesh | coupled 3D Newton and island/transport tests exist |
| Beltrami / MRxMHD | specialised regions | relaxed volumes and current sheets |
| Vacuum | one coil and exterior operator for both lanes | |
| Newton–Krylov | primary scalable solver | exact derivatives, indefinite systems, reused for stability |

This resolves a conflict in the documents: section 0 of
[free-boundary-ale.md](free-boundary-ale.md) chose the constrained (A, χ)
and its section 4 still develops the frozen-map relabelling; the label plan
uses direct Clebsch, which removes the locking documented in
[eulerian-3d.md](eulerian-3d.md). The frozen-map relabelling (gates B1–B3,
R1) is closed; direct Clebsch is the working ideal lane, and gate C compares
it with the constrained A numerically, not by another redesign.

**Architectural rule.** Mesh coordinates are never physical unknowns, and
never added to cure conditioning. The mesh moves by its own equation
(Winslow, alignment, refinement) between physics solves; the transfer
preserves the represented fields and constraints, and mesh quality is
measured independently of mesh motion. Physical points such as magnetic
axes and X-points are physical unknowns of the field, not mesh vertices; the
mesh follows them afterwards.

## Sequence

Each step on cases that run in seconds, against analytic solutions first.

1. **ALE label solver, one axis (A1).** Module (not a test) for the
   axisymmetric label lane on a spline O-grid: physics Newton on (χ, g, X)
   on a frozen mesh, with X inside a cell (singular quadrature at X), then
   Winslow update of the mesh with its centre attracted to X and the wall
   fixed, then transfer of χ and g by L2 projection; loop until X and the
   mesh are stationary. Passed when the Solov'ev error equals the error of
   E1 (mesh tied to X), Newton converges quadratically on each frozen mesh,
   and the mesh loop converges in a few passes from an axis offset of 10 % of
   the minor radius. Then the same on the per-plane 3D O-grids (E4 case).
   *State:* passed (`milca_label_ale`, `milca_spline_ale`, gate
   `label_ale_axisym_gate_a1`, 10 s). X is a fixed parametric point of the
   core block; in Newton the core translates with X and the rest follows a
   fixed harmonic blend (1 on the core, 0 on the wall), so the residual is
   smooth in X and Newton is quadratic (4–8 steps, one mesh pass). B error
   2.23e-2, 4.52e-3, 9.47e-4 on 2, 4, 8 cells (E1: 2.15e-2, 4.44e-3,
   9.25e-4), offsets up to 30 % of the minor radius in all directions (the
   centring backtracks where the coarse ring would fold). Two findings:
   pinning one point of the Winslow map to X (a point constraint on a
   harmonic map, singular in the continuum) made the axis drift and B
   converge at order 1; the core must be held as a region. Winslow between
   core and wall converges at order p but with 1.6–2× the B error of the
   blend on this case, so the blend is the default and Winslow is optional
   until a case needs it (large distortion, A2).
2. **ALE label solver, several axes (A2).** Doublet and two-axis
   configurations: one region per lobe and one outer region, the separatrix a
   region boundary of the mesh, each O-point a physical unknown followed by
   its own mesh centre, the X-point followed by the separatrix vertex.
   First on the exact Solov'ev doublet with fixed labels (E2b), then Newton
   (E2, open: the log-singular flux at the separatrix shrinks the basin to
   the graded layer; the ALE separatrix alignment is the planned remedy).
   *State:* A2a passed (`milca_label_ale_multi`, gate
   `label_ale_doublet_gate_a2a`, 8 s): the exact Solov'ev doublet on
   FortFEM's doublet layout in the Grad–Shafranov form of the label
   (ψ = ψ_O + sχ, f prescribed, no winding form), unknowns χ, both
   O-points and their multipliers. Each core follows its O-point by its own
   blend (0 on the wall, the other core and the separatrix sides), the X
   vertex follows the saddle of χ_h between passes by its own blend, and the
   pins are renewed each pass in the original mesh so a core never moves by
   more than half a cell. From O-points offset by 10 % of the lobe half
   widths (up to 30 % tested) and the X vertex off by half that:

   | cells | B_pol | projected exact χ | B_pol in the lobes | O error | X vertex error | passes, Newton steps |
   |---|---|---|---|---|---|---|
   | 2 | 8.83e-2 | 9.32e-2 | 3.34e-2 | 3.4e-3 | 1.5e-4 | 5, 8 |
   | 4 | 2.00e-2 | 2.03e-2 | 8.26e-3 | 1.7e-4 | 3.1e-5 | 5, 7 |
   | 8 | 8.75e-3 | 8.80e-3 | 1.94e-3 | 1.5e-5 | 7.7e-6 | 5, 8 |

   Newton is quadratic in every pass. B_pol equals the projected exact label
   on the same mesh; in the lobes it converges at order 2, the outer ring to
   the wall is pre-asymptotic for both (order 1.2 from 4 to 8 cells, 2.5 from
   8 to 16). The O-points converge at order 3.5–4, the X vertex at order 2.
   Findings: χ_h is only C0, so at the X vertex its saddle is a kink
   (Newton on ∇χ_h circles the vertex); the saddle comes from a
   kernel-weighted quadratic fit of χ_h (smooth in χ_h, exact for a
   quadratic label), which also makes the X follow converge. A fixed pin far
   from the O-point folds the thin lobe rings (the lobes have aspect 2.3,
   the cores are 0.4 of them); renewing the pins removes that. E2 in the
   direct-Clebsch form (gate `label_ale_doublet_clebsch_gate_e2`, O and X
   fixed): Newton converges quadratically at B_φ/B_pol ≈ 1 (toroidal field
   scaled by κ = 1e-3; B at the fixed-label fit on 4 and 8 cells) and fails
   at the doublet's B_φ/B_pol ≈ 1e3, equally with the log terms of the
   profiles removed. The obstruction is the guide field, not the separatrix:
   after eliminating g the label Hessian has O(B_φ²) negative grid-scale
   modes that swamp the O(B_pol²) poloidal stiffness. E1 degrades the same
   way (f0 = 30: labels at order 0.6; f0 = 300: Newton diverges), and every
   label gate so far ran at B_φ/B_pol ≲ 1 or ι ≈ 5. This blocks realistic
   tokamaks and stellarators (B_φ/B_pol ≈ 10–20, the QA deck) and comes
   before E2 with moving O-points, E5 and VMEC parity
   ([plan, E2 diagnosis](eulerian-label-plan.md)).
3. **Gate C: Clebsch against constrained A.** On the same ALE mesh:
   Solov'ev, a sheared 3D nested case and a simple island. Compare
   conditioning, Newton convergence, locking, cost and topology handling.
   *Input from E5u* ([plan](eulerian-label-plan.md), `label_e5_unaligned`):
   on a mesh not aligned with a separatrix, with only O- and X-points
   tracked, the label lane's B converges at h^(1/2) in the separatrix band
   even with g broken along the level set and the log slope factored out
   (2.2e-2 at 32 cells against 1.9e-3 aligned, the gap growing like
   h^(-1.3)). Field-line labels are log-singular at separatrices; B and A are
   not. Rule: only codimension-2 lines (O and X lines) may be tracked, so
   islands and separatrices go to the A representation (constrained A, NN2),
   and the label lane stays in regions bounded by tracked surfaces. Gate C
   decides whether the constrained A converges at order p across an
   unaligned separatrix without locking.
4. **VMEC parity in the ideal lane.** Fixed boundary, then free boundary
   (E6 with the exterior map and coils from mgrid), prescribed ι or I,
   realistic pressure, stellarator geometry. Against VMEC2000, VMEC++, GVEC
   and DESC at matched physical error ([exact benchmark](../benchmarks/free-boundary-exact/README.md)).
   Correctness and robust convergence first, speed after.
5. **NN2 on the shared ALE mesh**, with Winslow alignment against
   perpendicular pollution in anisotropic transport; existing 3D NN2 tests,
   then island formation and stochastic regions.
6. **Coupled regions**: labels where surfaces survive, NN2 or Beltrami
   elsewhere, with physical jump and conservation conditions at the
   interfaces (island divertor, mixed topology).
7. **Stability and reconstruction** on the same linearisations, keeping the
   ideal energy Hessian distinct from the dynamical linearisation of NN2.

Effort split (equilibrium development): ALE and common FE infrastructure
35 %, ideal Eulerian lane and VMEC parity 35 %, NN2 and island benchmarks
20 %, stability and reconstruction interfaces 10 %.

Milestone: a robust Eulerian nested equilibrium solver on an ALE–Winslow
mesh that reproduces VMEC fixed- and free-boundary results without the
slide gauges of the Lagrangian map. NN2 then inherits mesh and exterior.

## Optional, not an immediate priority

- Harmonic boundary-conforming maps with proven invertibility for
  nonconvex sections: Babin, Hindenlang, Maj and Köberl, *Construction of an
  invertible mapping to boundary conforming coordinates for arbitrarily
  shaped toroidal domains*, arXiv:2411.04683, code `map2disc`
  (gitlab.mpcdf.mpg.de/gvec-group/map2disc). They solve Δξ = Δη = 0 (the
  Winslow map) by boundary integrals rather than the nonlinear inverse
  equations, then invert it. Candidate for the initial mesh and for
  remeshing after the spline Winslow of step 1 works; check licence and
  record the commit before reusing code (AGENTS.md provenance rule). It
  supplies geometry only: the coupling to the Clebsch Newton, field transfer
  and adaptation stay MILCA's work.

## State of the map lane (reference)

- Axisymmetric free boundary at VMEC speed on the exact Solov'ev benchmark
  ([table](../benchmarks/free-boundary-exact/README.md)).
- 3D free boundary on a curve chart (`fb3d_synthetic`, R2 exterior): 44 s
  against 67 s for VMEC++, under-resolved at ntor 2.
- CTH-like decks fold the shell of the cylindrical chart (1.5 cm margin
  against a 4.6 cm helical swing). `milca_deck_chart` builds the
  axis-following chart from any deck; wiring it into `nested_equilibrium`
  is parked behind steps 1–2.
