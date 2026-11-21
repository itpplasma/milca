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

Clebsch coordinates do not generalise to non-nested surfaces, Clebsch
fields do.

## 3. Place among the lanes and prior work

In the lane table of [equilibrium-core.md](equilibrium-core.md) (D2) this is
a fourth lane, the label lane:

| Lane | Unknown U | B = P(U) | Constraint set | Gauge | Seed |
|---|---|---|---|---|---|
| label | χ, g ∈ V0, F_j, G splines in χ, O-lines X_j | (3) | B·∇χ = 0 exactly, p = p(χ), topology of the level sets of χ | g → g + k(χ), relabelling | ψ lane or map lane |

It is ideal and keeps the topology of χ, like the map lane, but on a fixed
or ALE mesh in physical space, like the free-B lane. It shares the D3
Maxwell-stress residual and the D8 cross-lane checks. Against the free-B
lane's NN2 closure (parallel pressure diffusion) it is the ideal limit with
labelled pressure.

Prior work (reading state: R method sections read, A abstract only):

- Kruskal and Kulsrud, Phys. Fluids 1, 265 (1958) [R], section H: B =
  ∇ψ×∇ν with ν = λ + η χ'(ψ) − θ, single-valued λ and multivalued angles,
  and W stationary over (p, λ, ψ). This is (2) for nested surfaces and one
  axis, as an analytic variational principle. The representation is not
  new.
- Terradas and Neukirch, A&A 671, A31 (2023) [R]: Euler potentials (α, β)
  as finite-element fields in a fixed 3D box, solved in physical space, with
  open, line-tied field lines and single-valued potentials. No toroidal
  topology, winding parts or net current.
- Cheng, GRL 22, 2401 (1995) [A] and Zaharia et al., Ann. Geophys. 22, 251
  (2004) [R]: B = ∇ψ×∇α for the magnetosphere, solved in inverse form.
- Chodura and Schlüter, J. Comput. Phys. 41, 68 (1981) [A], HINT, PIES,
  SIESTA [A]: Eulerian or mixed equilibrium codes without a label field.
  SIESTA's V3FIT coupling uses normalised pressure as a diagnostic label.
- Degond, Lozinski, Narski and Negulescu, J. Comput. Phys. 231, 2724 (2012)
  [R], and Deluzet and Narski, Multiscale Model. Simul. (2019) [R]: the
  kernel of B·∇ is reached by penalty, multiplier or iteration for a given
  b. The latter name the failure without functions constant along b in the
  discrete space locking (citing Babuška and Suri 1992), the counterpart of
  section 1.
- MRX (arXiv:2510.26986, cited in [equilibrium-nonnested.md](equilibrium-nonnested.md)):
  FEEC relaxation in physical space with δB = ∇×(v×B), no label.

Not found in this search (searches of questions on enrichment were thin):
a discretisation of (2) on a mesh in physical space with Newton, one winding
form per O-line for several axes and island chains, the O-lines as unknowns,
and the use of (2) against parallel-constraint locking. Separatrices and
changes of topology are not treated in any work read.

## 4. Discrete spaces

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

## 5. Gates (each runs in seconds)

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

## 6. Order of work

1. Literature: done, section 3. Repeat the enrichment search with a
   working web search before any publication.
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

## 7. Open

- Stochastic regions have no label. Candidate: regions with p' = 0, where the
  plasma part of (5) has no constraint, as in SPEC's relaxed regions.
- Inf–sup of the δχ equation near O-points and saddles, where ∇χ = 0.
- Uniqueness of the representation (2) with several O-lines in one region.
