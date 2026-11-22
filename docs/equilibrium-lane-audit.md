# Equilibrium lanes: audit of 2026-11-22

What was tried, what is a dead end, what was a mistake of method, and what
is specific to a representation. Facts are taken from the lane documents
named in each row; the verdicts are a judgement and are marked as such. No
runs were made for this audit.

## 1. Verdict per lane

| Lane | Verdict | Basis in the documents |
|---|---|---|
| (χ, g) labels with axis unknowns ([plan](eulerian-label-plan.md)) | superseded discretisation of the right field | polar form of the [flux map](flux-map.md); every gate ran at B_φ/B_pol ≲ 1 or ι ≈ 5 and degree 2; no realistic case converged |
| One label field across a separatrix on an unaligned mesh | dead end | order 1/2 in three variants (E2a, E5u); the label singularity is physical ([flux-map.md](flux-map.md), 5a) |
| Newton on shearless, all-rational cases (Landreman integer family, ι = 2) | limit of the physics, not a defect | fails from the exact state; the cause (kernel of B·∇ beyond flux functions) is argued, not isolated; the map lane and VMEC++ show an unstable mode on the same family ([benchmarks/nested](../benchmarks/nested/README.md)) |
| Constrained A with a label and a multiplier ([eulerian-3d.md](eulerian-3d.md)) | abandoned, not shown dead | "no stable pairing" is inferred from two failures; the enrichment V1 ⊕ span{f_k(χ_h) h} announced in the label plan was never run |
| Frozen-map relabelling ([free-boundary-ale.md](free-boundary-ale.md), B1–B6, R1) | abandoned untested | no result exists for any of its gates |
| Mesh that follows the magnetic axis (A1) | not needed | required only by the polar form; moving meshes remain for free boundaries and interfaces |
| Lagrangian map (`nested_equilibrium`, [ideal-map-3d.md](ideal-map-3d.md)) | working; shelved on weak grounds | section 3 |
| NN1, one Beltrami volume, Grad–Shafranov incl. X-points and two axes | working | orders and comparisons (SPEC, FreeGS) on record |
| NN2 ([equilibrium-nonnested-nn2.md](equilibrium-nonnested-nn2.md)) | working, heavy and fragile | degree 3 needed at anisotropy ≥ 1e6; 12–53 GB and up to 43 min on the Landreman family; Newton basin shrinking like 1/ratio |
| Coupled regions, multi-region relaxed states | not started | interface conditions derived by variation only (island-equilibrium.md C, current-sheet-equilibrium.md) |

## 2. Mistakes of method that repeat

- Null directions damped, not removed. A ridge or penalty first and the
  constraint later: ω = f(χ) in the map lane (penalty 1e-4, then moments:
  9/4/5 steps against 27/40/118), g → g + k(χ) in the label lane (ridge
  1e-12 and 1e-8; the constraint only in E6c), the closed-line gauge of Φ in
  NN2 (open). Each time Newton crept along the direction until it was
  removed.
- A stiff parameter shrinks the Newton basin: the guide field in the label
  lane and the flux map, the anisotropy ratio in NN2. Solved three times by
  ad hoc continuation; one shared treatment (continuation driver, merit
  function that tolerates the curved constraint) is owed.
- Gates only in easy regimes: weak guide field, large ι, degree 2. The main
  obstacle of the label lane stayed hidden until the QA deck.
- Causes recorded without an isolating run, later revised or still
  unsupported:
  - E2 doublet: first the separatrix, then the guide field (κ scan);
  - QA: the guide-field cause was by analogy with E1/E2, no run on QA;
  - "13 % at every mesh" is attributed to G2a in the label plan and
    reported under G2 in eulerian-3d.md, and no table shows it;
  - "gauges fail in strong shaping" and "under-resolved at ntor 2" for the
    map lane (section 3).

Rule from now on: a cause enters a plan document only with the run that
isolates it, otherwise it is written as a hypothesis; every new gate
includes a realistic ratio B_φ/B_pol (≥ 10), a low-shear case and degree
≥ 3; null directions are removed by constraints from the first version.

## 3. The Lagrangian map lane

- The reason given for stopping it (the roadmap before this audit; free-boundary-ale.md
  section 0: its gauges "fail in strong shaping") has no case and no number
  in either document. Measured items nearby are a weight-dependent end
  point of the shaped P3 seed ladder, the N1 gauge drifting on the gs_free
  fixture, and one unconverged ITER elongation step; the slide moments
  converge all six fixture cases to 4.4e-11.
- "Under-resolved at ntor 2" for the 3D synthetic free boundary is
  contradicted by free-boundary-nested-map-3d.md, section 9.3: ntor 3 and 4
  change the error by 0.25 mm, and the isolated cause is the n = 24 coil
  ripple of the oracle.
- It has: orders 3–4 on analytic cases (6.8 at degree 6), axisymmetric free
  boundary at the speed and accuracy of VMEC2000 and VMEC++, the exact
  energy Hessian, a curve chart, 3D free boundary with an exterior operator.
- Measured weaknesses: 70–150 times the time of VMEC++ at equal accuracy in
  3D fixed boundary (4.17 s against 0.028 s at 2,300 coefficients), memory
  like ntor², shell folds on CTH-like and DIII-D decks (the axis-following
  chart `milca_deck_chart` exists and is not wired in), 3D stalls traced to
  globalisation of the shell map.

## 4. What is specific to a representation

Not specific to the flux map (the map lane has them, and any inverse map
can): Newton on the exact Hessian and derivatives for optimisation; the
axis as an ordinary mesh point (disc mesh, not a polar grid); finite
elements or splines in the section; a toroidal coordinate along a guiding
curve; gauge by moment constraints; stability from the equilibrium
Hessian; an exterior or neighbouring model through an interface operator;
one map per region with interfaces.

Specific to the flux map: a mesh fixed in physical space (no chart or
shell folds; shared with A-based regions); a gauge reduced to
θ* → θ* + k(s), where an inverse map can relabel the poloidal angle by an
arbitrary function; a polynomial energy density (the inverse map divides
by the Jacobian). Its own costs: guide-field locking, the exact wall
condition, 3D accuracy not yet shown.

The comparison that decides the ideal lane is therefore flux map against
the Lagrangian map of this code, not against VMEC.
