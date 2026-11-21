# Label-field equilibrium without locking: plan

Question: which discrete spaces make the label formulation of
[eulerian-3d.md](eulerian-3d.md) stable and free of locking on a mesh that is
not aligned with the labels, in 3D and for non-nested surfaces?

## 1. Defect

Symbols as in [eulerian-3d.md](eulerian-3d.md): A the vector potential, χ the
label, λ̃ the multiplier of B·∇χ = 0, p(χ) and C(χ) the profiles, h_θ the
winding form about an O-line, h_φ the toroidal harmonic form.

The kernel of B·∇ is the flux functions f(χ). On a mesh not aligned with χ_h,
f∘χ_h is not a spline for nonlinear f, so:

- a field tangent to the surfaces with a nonlinear flux profile is not in the
  trial space of A, and the weak constraint with an equal-degree multiplier
  locks (G2a: B error 13 % at every mesh),
- a smaller multiplier space cannot absorb the pressure equation, and χ
  takes the error at O(h) (G2: label error 2e-1, 4e-2, 2e-2).

Measurements: [eulerian-3d.md](eulerian-3d.md), section 5.

## 2. Tangent fields in closed form

Every divergence-free field tangent to the level sets of χ is, locally,
B = ∇α×∇χ. Around O-lines and the torus, α is multivalued, and its
multivalued part is a flux function times an angle:

    α = g + Σ_j a_j(χ) θ_j + b(χ) φ,   g single valued.                 (1)

With ∇θ_j = h_θ,j (one winding form per O-line) and ∇φ = h_φ:

    A = g ∇χ + Σ_j F_j(χ) h_θ,j + G(χ) h_φ,   F_j' = −a_j,  G' = −b,     (2)

    B = ∇×A = ∇g×∇χ + Σ_j F_j'(χ) ∇χ×h_θ,j + G'(χ) ∇χ×h_φ,               (3)

    B·∇χ = 0  pointwise, for every g, F_j, G.                            (4)

(3) uses ∇×(g∇χ) = ∇g×∇χ and ∇×(F(χ)h) = F'∇χ×h, because ∇×h = 0 away from
the O-lines. The line current 2π F_j(χ(X_j)) δ_X at an O-line X_j vanishes
for F_j(χ(X_j)) = 0. F_j and G are the toroidal and poloidal flux profiles.
In axisymmetry the poloidal part of (2) is ψ = G(χ), the Grad–Shafranov
form, and F with g carries the toroidal field.

(2) is a direct (Eulerian) Clebsch representation. χ and g are fields in
physical space on a mesh that does not follow them. VMEC, GVEC and DESC are
inverse: they solve for the map x(s, θ, ζ) from flux coordinates. In (2)
the topology of the level sets of χ is free:

- several magnetic axes and island chains: one h_θ,j per O-line, profiles
  per region,
- separatrices: χ has a saddle, g and χ are smooth there, and h_θ,j is
  smooth away from its own O-line, so the q → ∞ of flux coordinates does
  not appear in the fields.

This corrects the earlier claim that a Clebsch field does not generalise to
non-nested surfaces: that holds for Clebsch coordinates, not for Clebsch
fields.

## 3. Discrete spaces

- χ, g ∈ V0 (degree p), profiles F_j, G as 1D splines in χ with
  F_j(χ(X_j)) = 0, O-lines X_j as unknown curves (two numbers per toroidal
  node each).
- A from (2), evaluated pointwise at the quadrature points. No multiplier,
  because (4) holds exactly: the locking of section 1 cannot occur.
- Outside the plasma (scrape-off layer, vacuum) A ∈ V1 as in G1, matched to
  (2) on Γ_p by the tangential trace of A.

Energy and stationarity (μ0 = 1):

    W[χ, g, F, G, X] = ∫ |B|²/2 dV − ∫ p(χ) dV.                          (5)

- δg: ∫ (∇δg×∇χ)·B dV = ∫ δg ∇χ·∇×B dV (no wall term, since ∇χ ∥ n on
  the wall), so J·∇χ = 0 weakly.
- δχ: the force balance J×B = ∇p across the surfaces.
- δF_j, δG: one equation per profile coefficient, the flux-surface averaged
  current balance. In current mode (VMEC's I(s)) the profiles of the toroidal
  current are constraints with scalar multipliers per coefficient.

Redundancies of (2), to be removed, not regularised:

- g → g + k(χ) leaves B unchanged. Remedy: g orthogonal to the 1D splines in
  χ_h (K scalar constraints, the same K as the profile splines).
- relabelling χ → m(χ) with p, F, G carried along. Remedy: the label is fixed
  by a profile, e.g. the toroidal flux F_1(χ) = χ Φ_edge.

Pointwise evaluation of f(χ_h) needs a quadrature that resolves the profile
splines. Use profiles one degree smoother than V0 and check the quadrature
error per gate.

## 4. Gates (each runs in seconds)

| Gate | Case | Pass criterion |
|---|---|---|
| E0 | G2a geometry, fixed χ_h, nonlinear label, unknowns g and profiles | B at the interpolation error, no locking at any mesh |
| E1 | G2: Solov'ev, χ, g, profiles and axis unknown, both labels | B and χ at order p, quadratic Newton, profile error at the 1D spline error |
| E2 | Solov'ev doublet (two magnetic axes, one separatrix) | as E1, both axes found, separatrix flux to the interpolation error |
| E3 | E1 through the 3D complex (M ≥ 2, n = 0 solution) | identical to E1 |
| E4 | low-β stellarator, fixed boundary | VMEC++ and GVEC within their convergence |
| E5 | tokamak with a resonant boundary perturbation, island chain | island width against a published SPEC or HINT case |
| E6 | free boundary: plasma (2) coupled to V1 vacuum and the wall operator | the 2D spline free-boundary results, then Landreman–Paul QA (L5) |

E0 and E1 also run the multiplier formulation with the same enrichment of A
(V1 ⊕ span{f_k(χ_h) h}) for comparison. The representation that passes E1
with fewer unknowns and better Newton convergence is kept.

## 5. Order of work

1. Literature: Zotero and web for Eulerian (direct) Clebsch equilibrium
   solvers, partition-of-unity or level-set enrichment (GFEM, XFEM) used for
   flux functions, asymptotic-preserving treatments of B·∇. Note what exists
   and where (2) differs. Novelty is stated only after this check.
2. FortSym leaves for the pointwise terms of (3) and (5) and their first and
   second derivatives in (χ, g, profile coefficients), with identities.
3. E0, E1 in axisymmetry (MILCA tests next to G2a, G2).
4. O-line handling for several axes, E2.
5. 3D: h_θ about an axis curve, θ_X = atan2(Z − Z_X(φ), R − R_X(φ)) per
   toroidal angle, closed away from X. E3, then E4.
6. Islands and rational surfaces, E5. Open whether Newton from nested
   initial labels finds the island topology, or whether the topology has to
   be seeded.
7. Free boundary, E6, and the solver: Newton–Krylov with AMS for the V1 part
   and AMG for χ and g.

## 6. Open

- Stochastic regions have no label. Candidate: regions with p' = 0, where the
  plasma part of (5) has no constraint, as in SPEC's relaxed regions.
- Inf–sup of the δχ equation near O-points and saddles, where ∇χ = 0.
- Uniqueness of the representation (2) with several O-lines in one region.
