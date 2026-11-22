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
- separatrices: χ has a saddle and is smooth there, h_θ,j is smooth away
  from its own O-line. The q → ∞ of flux coordinates moves into g: with θ*
  the straight-field-line angle, g = Φ'(χ)(θ* − θ) + k(χ), and Φ' grows like
  log|χ − χ_s|. g jumps across the separatrix by a function of this
  amplitude. B does not see the jump (its gradient is normal to the
  surfaces), but a continuous g on an unaligned mesh does (E2a).

Clebsch coordinates do not generalise to non-nested surfaces, Clebsch
fields do.

Relation to flux coordinates. VMEC writes B = ∇ζ×∇Ψ + ∇Φ×∇(θ + λ) with Ψ,
Φ the poloidal and toroidal flux. With α = g − F'(χ)θ − G'(χ)φ the field-line
label of (2), F ↔ −Φ/(2π), G ↔ Ψ and g = Φ'(χ) λ up to k(χ):

| | VMEC | label lane |
|---|---|---|
| surfaces | R(s, θ, ζ), Z(s, θ, ζ): positions of fixed labels | χ(R, φ, Z): the label at fixed positions |
| field-line angle | λ(s, θ, ζ) | g(R, φ, Z) |
| profiles | Φ(s), ι(s) or current, p(s) | F(χ), G(χ), p(χ) |

In axisymmetry the lane is Grad–Shafranov with ψ = G(χ). Straight-field-line
angles follow from θ* = θ + g/F', Boozer angles from them by the usual
transform. SPEC keeps a few nested interfaces and relaxed (Beltrami) volumes
with flat pressure, the lane keeps labels and a continuous p(χ) everywhere.
Its field is integrable by construction (B·∇χ_h = 0 pointwise): a thin
chaotic layer at a separatrix is replaced by the separatrix, wide chaos needs
an unlabelled region (section 7).

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
- Chodura and Schlüter, J. Comput. Phys. 41, 68 (1981) [R]: the same
  argument against flux coordinates (p. 69: nested surfaces as Lagrangian
  coordinates allow one axis, a space-fixed Eulerian grid several, at the
  price of flux not attached to the plasma motion). W = ∫(B²/2 + p/(γ − 1))
  is minimised under δρ = −∇·(ρξ), δB = ∇×(ξ×B) by friction, gradient and
  conjugate-gradient iterations on a staggered grid (mass at cell centres,
  magnetic fluxes on faces, E on edges, force from the discrete energy). The
  Solov'ev test converges at order 2 in the displacement (Fig. 7), ι drifts
  by flux loss until a 20³ grid (Fig. 8), and a double-star equilibrium with
  two axes and a separatrix is computed (Fig. 13, quadrupole on an l = 3
  field, A = 16, β = 0.03). No label: tangency and ι are not preserved
  discretely, and topology is whatever the relaxation reaches.
- Hender, Carreras, Garcia, Rome and Lynch, J. Comput. Phys. 60, 76 (1985)
  [R], NEAR: the Chodura–Schlüter relaxation on vacuum flux coordinates used
  as a fixed (Eulerian) grid, Fourier series in both angles and finite
  differences in radius, chosen "because it allows the study of equilibria
  whose flux surface topology changes with pressure" and shown to be far
  more economical than the cylindrical grid; it agrees with Chodura–Schlüter
  at low β and departs where their resolution is short. Panici et al.,
  J. Plasma Phys. 89, 955890303 (2023) [R, sections 1–2] place both in the
  move to spectral angle representations.
- Reiman, Pomphrey and Boozer, Phys. Fluids B 1, 555 (1989) [R]: the
  limiting form of a general 3D equilibrium near a separatrix with nonzero
  toroidal field. In hyperbolic coordinates ξ = aρ cosh α, η = cρ sinh α about
  the x-line (a, c from the second derivatives of the helical flux), the
  magnetic angle is Θ = −a c B_0 α/(2 q_h) + Θ_1 and q_h ∝ ln(1/√Ψ_h): the
  straight-field-line angle carries the hyperbolic angle α, which is the
  singular part of g near X (section 4a, E2b). ∇p vanishes on the x-line and
  the pressure-driven current stays finite at the separatrix, unlike on a
  rational surface. They propose using this local solution in equilibrium
  codes, as axis expansions are used.
- HINT, PIES, SIESTA [A]: Eulerian or mixed equilibrium codes without a label
  field. SIESTA's V3FIT coupling uses normalised pressure as a diagnostic
  label.
- Degond, Lozinski, Narski and Negulescu, J. Comput. Phys. 231, 2724 (2012)
  [R], and Deluzet and Narski, Multiscale Model. Simul. (2019) [R]: the
  kernel of B·∇ is reached by penalty, multiplier or iteration for a given
  b. The latter name the failure without functions constant along b in the
  discrete space locking (citing Babuška and Suri 1992), the counterpart of
  section 1.
- Cianciosa, Hirshman, Seal and Shafer, Plasma Phys. Control. Fusion 60,
  044017 (2018) [R]: V3FIT with VMEC + SIESTA reconstructs DIII-D 154921
  with islands from 377 signals. The island size is the amplitude of an
  artificial helical seed perturbation (with artificial resistivity) fitted
  to the data, the phase follows from fitted n = 1 coil currents, the
  Jacobian is finite differences, and stellarator island divertors are named
  out of reach. In (2) the island size is the island toroidal flux, a
  physical invariant, the nested solution is its zero limit, and Newton gives
  the derivatives (section 4a).
- MRX (arXiv:2510.26986, cited in [equilibrium-nonnested.md](equilibrium-nonnested.md)):
  FEEC relaxation in physical space with δB = ∇×(v×B), no label.

Eulerian equilibria with several axes and a separatrix are not new
(Chodura and Schlüter 1981). Not found in this search (searches of questions
on enrichment were thin): a discretisation of (2) on a mesh in physical space
with Newton, which keeps the integrable topology exactly (B·∇χ_h = 0), one
winding form per O-line for several axes and island chains, the O-lines as
unknowns, the island flux as a parameter, convergence at order p, and the use
of (2) against parallel-constraint locking. Each claim stands only with its
gate (E2b, E5).

Measures beyond the force residual, after Chodura and Schlüter's flux-loss
experience: enclosed toroidal and poloidal flux per region, ι on the
surfaces, and the topology (number of O- and X-lines) against the reference.
Their double star is a candidate nonlinear multi-axis case once E2 runs, and
the same two-axis equilibrium in the free-field lane would compare the two
formulations.

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

## 4a. Regions, islands and the island size

With several O-lines, A = g∇χ + Σ_j F_j^r(χ) h_j + G^r(χ)∇φ in region r
(each lobe or island, the region between separatrix and wall):

- in region r, F_j^r is constant for every O-line j that the surfaces of r
  do not enclose (h_j is closed there, so it adds no B), and F_j^r(χ_s) is
  continuous across the separatrix, so A is continuous,
- outside a doublet's separatrix the split of the toroidal flux between h_1
  and h_2 is a gauge choice absorbed by g; up-down symmetry fixes it,
- region membership of a point is topological (the side of the
  separatrix), not a value of χ: inside an island χ takes values of the
  outer region,
- at an elliptic O-point θ_j is the angle in the frame where the Hessian of
  χ_h is isotropic, θ_j = atan2 of H^(1/2)(x − X_j). The circular angle
  leaves θ* − θ = O(1) at the axis, so g is discontinuous there (E2a).

Ideal MHD conserves the flux of each region, so the island size is an input:
the island toroidal flux ΔΦ_I = −2π F_I(χ_s), with the label inside the island
the normalised island flux. It is prescribed (from data or from a resistive
or two-fluid response), zero (ideal response: no island, a shielding current
sheet at the rational surface), or fixed by ∂W/∂ΔΦ_I = 0 at fixed helicity
(one scalar equation per chain). ΔΦ_I → 0 joins the nested solution, and the
Newton factorisation gives the derivatives with respect to ΔΦ_I, the
profiles and the boundary: island width and O-line phase become parameters
of a reconstruction from measured profiles and magnetics.

## 5. Gates (each runs in seconds)

| Gate | Case | Pass criterion |
|---|---|---|
| E0 | G2a geometry, fixed χ_h, nonlinear label, unknowns g and profiles | B at the interpolation error, no locking at any mesh |
| E1 | G2: Solov'ev, χ, g, profiles and axis unknown, both labels | B and χ at order p, quadratic Newton, profile error at the 1D spline error |
| E2a | doublet with fixed labels, g fitted, unaligned mesh | measurement: order p in the lobes and outside, the separatrix band at h^(1/2) |
| E2b | E2a on a separatrix-aligned mesh, g broken across the separatrix | acceptance below; before any Newton with several O-lines |
| E2 | Solov'ev doublet (two magnetic axes, one separatrix) | as E1, both axes found, separatrix flux to the interpolation error |
| E3 | E1 through the 3D complex (M ≥ 2, n = 0 solution) | identical to E1 |
| E3b | E3 with one axis point per toroidal plane, non-axisymmetric start | identical to E3, axis curve a circle |
| E4a | exact 3D Landreman equilibria, fixed labels, g fitted | B at order p against the exact field |
| E4 | E4a with χ, g and the axis curve by Newton (sheared family), then the registered integer case against VMEC++, DESC and GVEC | B and χ at order p; time to accuracy |
| E5 | tokamak with a resonant boundary perturbation, island chain, ΔΦ_I as input | continuous to the nested solution at ΔΦ_I → 0, derivatives against finite differences, minimum-energy width against a published SPEC or HINT case |
| E6 | free boundary: plasma (2) coupled to V1 vacuum and the wall operator | E6a: axisymmetric, exact manufactured solution; E6b: 3D, exact (oracle field as source); then Landreman–Paul QA (L5) |

E0 and E1 also run the multiplier formulation with the same enrichment of A
(V1 ⊕ span{f_k(χ_h) h}) for comparison. The representation that passes E1
with fewer unknowns and better Newton convergence is kept.

E0 result (degree 2, cubic profile splines on 8 + 4n intervals, relative
ridge 1e-12 for g → g + k(χ_h), relative L2 error of B):

| cells per block side | 2 | 4 | 8 | 16 |
|---|---|---|---|---|
| χ = ψ_N | 2.17e-2 | 4.53e-3 | 9.36e-4 | 2.20e-4 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 2.45e-2 | 5.62e-3 | 1.31e-3 | 3.28e-4 |
| projection of the exact A | 2.20e-2 | 4.37e-3 | 9.02e-4 | 2.12e-4 |

The nonlinear label converges at order p (rates 2.12, 2.10, 2.00) instead
of locking at 13 %.

E1 (Newton on (χ, g, X, μ) with given F, G, p, FortSym leaf
`clebsch_energy_axisym`):

- ∇×(F(χ) h_θ) contains 2π F(χ(X)) δ_X, which the quadrature does not see.
  Without the constraint χ_h(X) = 0 Newton converges to a state with B off
  by 10 % at every mesh.
- With X fixed while argmin χ_h moves, F' ∇χ×h_θ ~ 1/ρ at X. X is therefore
  the label minimum, ∇χ_h(X) = 0, as a Newton unknown. An outer update of X
  oscillates.
- The gradient terms F' ∇N×h_θ are ~ 1/ρ. Gauss quadrature gets the cells
  at X wrong at O(1). A polar (Duffy) rule in the cells whose closure
  contains X fixes them.
- With X free in a fixed mesh, the residual jumps whenever X crosses a cell
  edge, and Newton stalls from perturbed labels. ALE removes this: the core
  corners of the O-grid move with X, so X stays the central vertex of the
  core block and the polar cells never change. The X columns of the Jacobian
  are central differences of the residual with the mesh moving.

E1 result (degree 2, 10 Gauss points per direction, labels perturbed by 5 %,
FortFEM's polar rule at the vertex X):

| cells per block side | 2 | 4 | 8 | 16 |
|---|---|---|---|---|
| L2 error of B | 2.15e-2 | 4.44e-3 | 9.25e-4 | 2.17e-4 |
| interpolant | 2.20e-2 | 4.37e-3 | 9.02e-4 | 2.12e-4 |
| max \|χ_h − Π0 ψ_N\| | 3.6e-3 | 3.8e-4 | 6.6e-5 | 9.5e-6 |
| Newton steps | 6 | 4 | 3 | 5 |

The labels converge at about order p + 1 (rates 3.3, 2.5, 2.8). With 6
Gauss points per direction the polar rule is too coarse and the label error
at 16 cells rises to 1.5e-4.

E2 case: ψ = R⁴/8 + 1.2373 R² − 3.3090 ψ3 − 2.0931 ψ4 + 0.4209 ψ5 +
0.1550 ψ6 + 0.00694 ψ7 (R0 = 1, ψ_k the up-down symmetric homogeneous
solutions of Cerfon and Freidberg, Phys. Plasmas 17, 032502 (2010), each
checked numerically): magnetic axes at (1, ±0.170), X-point at (1, 0), wall
level ψ_X + 9.3e-4 closed within 0.27 of the X-point.

E3 and E3b (3D leaf `clebsch_energy_3d`, nt = 3 planes):

- E3 reproduces E1 (B error 2.148e-2, 4.443e-3, 9.251e-4 at 2, 4, 8 cells),
  labels axisymmetric to 1e-11.
- E3b gives each plane its own axis X_k. The section map is s(ξ, t) =
  Σ_k L_k(t) s_k(ξ) with V0 node values s_k, the core following X_k, so the
  axis conditions stay linear in χ. Exact checks: W of g = R, χ = Z is
  independent of the axis curve to 3e-15, and of a section translated with φ
  to 1e-14.
- From labels differing by ±50 % of a 5 % perturbation between planes and
  axes offset by 2e-3, Newton reaches the E3 errors in 9 (2 cells) and 8
  (4 cells) steps with a non-monotone line search; the planes agree to
  1e-11. Axis offsets of 1e-2 leave the basin.
- A piecewise-cubic profile table (F' jumping between pieces) leaves a
  Jacobian mismatch of 1e-5 and costs quadratic convergence. A Chebyshev
  series fixes it (mismatch 1e-9).
- A wall translated helically by 0.1(cos φ, sin φ), Solov'ev profiles: Newton
  converges in 6 steps (2 cells), the axis follows with ±5.7e-3 in Z. No
  reference yet (E4).

E2a (doublet, fixed χ_h, exact profiles per region from level-curve tracing,
relative L2 error of B_φ, degree 2, O-grid centred at X):

| cells per block side | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| total | 2.20e-1 | 1.53e-1 | 1.11e-1 | 8.47e-2 |
| lobes (χ < 0.3) | 7.70e-2 | 2.94e-2 | 8.71e-3 | 2.30e-3 |
| outside (1.7 < χ < 1.95) | 6.10e-2 | 2.47e-2 | 5.27e-3 | 1.25e-3 |
| lobes, circular winding angle | | | 2.96e-2 | 1.79e-2 |

The total converges at about h^(1/2) (rates 0.52, 0.47, 0.39) from the band
at the separatrix, where g jumps. Remedy: a mesh with the separatrix as block
edges and the X-point as a vertex, V0 discontinuous across it and graded
towards it, moving with χ (requested from FortFEM).

E2b acceptance (fixed labels, B against the exact Solov'ev doublet):

| Check | Required |
|---|---|
| B away from the separatrix and X | order p, as in E2a |
| B in a band along the smooth separatrix branches | order p or close to it |
| B in a disc around the X-point | order p or close to it, B_φ(X) = f0/R_X |
| jump of B_φ across the separatrix | at the approximation order: no spurious current sheet |
| labels exact χ, projected χ_h, χ_h constrained to the exact separatrix | separates the error of a misplaced separatrix from that of g |

The exact g is log-singular near X, g ≈ (B_0/2) log|(x + y)/(x − y)| for
χ − χ_X = (x² − y²)/2. B sees only |∇χ| times its tangential derivative, which
is bounded, so the aligned, graded, broken space should reach order p. If
the disc around X falls short, the fallback is a local enrichment
g = g_sing + g_h with this form, evaluated together with the winding term;
Reiman, Pomphrey and Boozer (1989) give its general 3D version through the
hyperbolic angle α about the x-line, with a, c and γ from the Hessian of χ_h
there.
A smooth carrier B_reg of the guide field is a diagnostic in axisymmetry
only: in 3D a divergence-free field tangent to χ is again of Clebsch form and
carries the same singularity.

E2b result (`test_eulerian_clebsch_doublet_aligned`, FortFEM `layout_doublet`
on the traced exact separatrix, χ in continuous V0, g in V0 cut along the
separatrix, degree 2, graded layers of power 4, relative L2 error of B_φ;
band |χ − 1| < 0.2 outside the X disc of radius 0.04; layer: relative RMS of
R B_φ/f0 − 1 in 0 < |χ − 1| < 1/(2n) outside the disc):

| cells per block side | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| total | 6.06e-2 | 3.06e-2 | 6.84e-3 | 1.91e-3 |
| lobes (χ < 0.3) | 8.85e-3 | 3.32e-3 | 8.16e-4 | 2.02e-4 |
| outside (1.7 < χ < 1.95) | 1.76e-2 | 1.21e-2 | 2.32e-3 | 4.86e-4 |
| band | 3.63e-2 | 1.71e-2 | 3.99e-3 | 1.28e-3 |
| X disc | 2.65e-2 | 1.05e-2 | 3.52e-3 | 1.07e-3 |
| layer | 9.02e-2 | 4.89e-2 | 1.66e-2 | 8.53e-3 |
| mesh separatrix, max \|χ − 1\| | 8.0e-2 | 1.2e-2 | 1.3e-3 | 8.3e-5 |

Projected label; the label constrained to 1 on the separatrix gives the same
numbers to three digits. Rates 16 → 32: total 1.84, lobes 2.02, outside 2.25,
band 1.64, X disc 1.71, layer 0.96. Against E2a (total 1.11e-1 at 16 cells,
order 1/2) the aligned broken space restores order p, at X included, and the
falling layer error excludes a current sheet; the singular enrichment is not
needed. The registered gate runs 4, 8, 16 in 42 s and checks the 8 → 16 rates.

- Grading. FortFEM's power grading (i/n)^β along a whole block chain keeps
  cells of size about β/n away from the separatrix, so the smooth regions
  lose order before the singular ones gain it (β = 8: lobes at rate 0.8 from
  4 to 8 cells; β = 4: outside at 0.4 from 8 to 16). A graded layer fixes it:
  a fraction θ = 1/2 of the chain's cells at τ (i/m)^β, the rest uniform, with
  τ = θ/(β(1 − θ) + θ) so the cell size is continuous. β = 4 suffices at
  p = 2 (β = 6 is similar, β = 8 is no better at 32 cells).
- Conditioning. The graded cells spread the diagonal of the normal equations
  over many decades. A ridge relative to the largest entry swamps the rows
  of the smallest cells, and the error then grows under refinement (β = 8,
  32 cells: 33 % in the lobes). Symmetric Jacobi scaling with the ridge
  relative to the scaled diagonal, plus iterative refinement against the
  unscaled matrix, removes it: at β = 8, 32 cells the factor leaves a relative
  residual of 2e-2 (condition about 1e14), two corrections 1e-5.
- Profiles at X. Graded cells at X reach |χ − 1| < 1e-10. dΦ/dχ follows
  A log|χ − 1| + B there (1.12e-2 per two decades in each lobe from 1e-2 to
  1e-8), so the reference profile is extrapolated along it below the
  tabulated 1e-10 instead of clamped.
- The exact label χ on the mesh, whose separatrix is the degree-p
  isoparametric image of the traced one (misplacement at order p + 1, 2.7 to
  4.0 above), puts points of the band on the wrong side of the
  log-singular profiles: band rate 1.2, against 2.1 for the discrete labels
  that the mesh separatrix carries. In Newton (E2) χ_h is the unknown, so the
  discrete case is the relevant one.
- The one-sided traces of B_φ are no sheet measure: each Clebsch term is
  log-singular on the separatrix and only their sum is finite.

E2 attempt (Newton on the doublet; open, not converged). Formulation tried:
the axisymmetric system with regions (χ in continuous V0, g cut along the
separatrix, per quadrature point the region's profile and winding form:
lobe j about O_j in its H^(1/2) frame, outside the mean of both windings with
the full outer flux slope, so the generated leaf `clebsch_energy_axisym` is
unchanged), χ = 2 on the wall, multipliers for χ and grad χ at both O-points,
staggered ALE (Newton on a fixed mesh, then O, X and the separatrix of the
discrete label for the next mesh), Jacobi-scaled KKT solves. Profile
derivatives are those of the tabulated interpolant. Findings at 4, 8 and 16
cells per block side, started from the projected exact label on the exact
mesh:

- χ free on the separatrix: the residual is differentiable only on the scale
  1e-10 (J v against central differences agrees to 1e-6 at h = 1e-10, 140 %
  at h = 1e-8): Gauss points of the graded layer sit at |χ − 1| ≈ 1e-10 and
  cross the log singularity of F' under any step.
- χ = 1 pinned on the mesh separatrix (each region a Dirichlet problem; the
  separatrix rows become the sheet-current measure of a staggered shape
  update) and |χ − 1| smoothed to ((χ − 1)² + ε²)^(1/2) in the profiles: J
  consistent to 1e-5 at ε = 1e-4, yet Newton stalls. Along the Newton
  direction the second-order term is 1e2 to 1e5 times the linear one,
  concentrated in the separatrix layer and on the O–X line of each lobe.
- Spectrum of the scaled Newton matrix (4 cells): six exact null modes of g
  (a constant and χ per region; removed by g = 0 on one transversal node line
  per region), then weak coupled χ–g modes (scaled λ 1e-5 to 1e-4, about
  95 % g) in the graded layer next to the separatrix, which carry most of
  the step. Exact elimination of g per trial step (W is quadratic in g) and
  a best-point line search do not restore convergence; ε = 0.1 does not
  converge in 16 steps either.
- Cause: in the layer g ≈ F'(χ)(θ* − θ) with F' ~ log|χ − 1|, so the residual
  responds to δχ like δχ/(χ − 1): the Newton basin shrinks to the smallest
  |χ − 1| of the layer nodes, and refinement (stronger grading) makes it
  smaller. Candidates: unknowns adapted to the layer (χ − 1 in logarithmic
  form, or the separatrix as an interface with χ − 1 proportional to the
  distance and its position an unknown), continuation in ε with the layer
  released from the outside in, or an axisymmetric doublet in the
  Grad–Shafranov form (f = R B_φ prescribed) to separate the separatrix
  topology from the label-lane representation.

E4a (fixed χ_h = Π0 ψ/δ per plane, per-plane O-grid on the exact axis fitted to
the exact boundary, winding angle in the Hessian frame of each plane,
relative L2 error of B against the exact field, degree 2):

| case | planes | 2 cells | 4 cells | 8 cells |
|---|---|---|---|---|
| integer family, ε = 0.1 | 13 | 1.07e-2 | 2.58e-3 | 5.59e-4 |
| integer family, ε = 0.1 | 21 | | | 4.96e-4 |
| sheared family, ε = 1.2 | 21 | 2.81e-2 | 1.03e-2 | 2.77e-3 |
| sheared family, ε = 1.2 | 31 | 2.81e-2 | 1.02e-2 | |

- The sign pair F = −Φ/(2π), G = −ιΦ/(2π) is selected by the data: the
  other three pairs miss B by 35–200 % at 2 cells.
- The sheared case needs 21 planes (axis R from 1.02 to 1.44 over φ); with
  13 the planes limit the error.
- Newton on the integer family does not converge, even from the exact state
  and at ε = 0: ι = 2 everywhere, so every surface is rational and shearless,
  B·∇ has a kernel far beyond functions of χ, the g block is nearly singular
  and the label equation loses its stiffness. The integer family stays a
  fixed-label (E4a) and cross-code case; Newton moves to the sheared family.
- The multiplier of χ(X) = 0 is nonzero and mesh independent (0.103 for
  ε = 0): the winding form leaves a point term 4πκF'² δχ(X) in the first
  variation of W, which the constraint absorbs.
- The O-grid of each plane is centred on its axis (FortFEM
  `layout_set_curve(..., centre_core=.true.)`, #277): before, the core
  corners were scaled about the axis with the wall, so the core block's
  centre left the axis wherever the boundary is not centred on it, and on
  the sheared sections the axis rows stayed at 0.5 from the exact state.
  The integer row above is with the centred core (before: 1.09e-2, 2.47e-3,
  4.99e-4; still order p). The Solov'ev gates keep the default core (scaled
  about r0): centring it there crowds the outboard cells and costs E1 its
  quadratic convergence at 4 cells.
- The 3D assembly is threaded over cells (OpenMP): the E3b gate takes 1.4 s
  and E4a 2.4 s on 16 threads, against 11 s and 17 s serial. The finite
  difference axis columns (2 nt gradient assemblies per Newton step) remain
  the cost at 21 planes.

E4 result (`test_eulerian_clebsch_landreman_newton`, sheared family ε = 0.6,
S = 2.4, δ = 0.2: ι from 4.73 to 5.11, the axis from R = 1.456 to 1.646 over
a field period; 13 planes, degree 2; χ = ψ/δ, F = −Φ/(2π), G' = −ιΦ'/(2π),
p = −ψ/λ², with Φ(ψ) and ι(ψ) Chebyshev series of 24 terms fitted to the
oracle (last coefficients 1e-16); start: labels perturbed by 2 % per plane,
the axis offset helically by 2e-3, g fitted):

| cells | 2 | 4 |
|---|---|---|
| L2 error of B, Newton | 4.86e-2 | 1.16e-2 |
| fixed-label fit (E4a) on the same mesh | 4.99e-2 | 1.16e-2 |
| max \|χ − Π0 ψ/δ\| | 2.3e-2 | 3.1e-3 |
| max \|X_k − X_exact\| | 8.1e-3 | 2.3e-3 |
| Newton steps (\|r\| to 1e-8 of the start) | 4 | 4 |

Rates: B 2.06, χ 2.89, axis 1.81. \|r\| at 4 cells: 7.5e-2, 1.0e-3, 8.9e-6,
1.0e-7, 2.8e-11. The gate runs in 17 s on 28 cores (the finite-difference
axis columns, 4 nt assemblies per step, dominate).

- The Newton matrix agrees with central differences of the residual to
  1e-12 in χ, g and μ and 1e-9 in X (argument `jacobian`), with the wall
  rows excluded (the residual is zero there and the solve fixes those
  labels).
- The discrete g → g + k(χ_h) is a near-null direction: with the ridge of
  the g solve (1e-12 relative) Newton steps drift along it (\|Δg\| 3e-2
  every other step at \|r\| = 1e-7). A ridge of 1e-8 on the g block of the
  Newton matrix only (not the residual, so the solution is unchanged)
  removes the drift. The proper remedy is the gauge constraint of section 4.
- Unresolved: the sheared cases with ε = 1.2 do not converge. With δ = 1/64
  (the E4a case, ι from 2.7997 to 2.8007) Newton diverges at 2 cells even
  from the exact state; with S = 2.4, δ = 0.2 (ι from 4.55 to 4.78) it
  converges at 2 cells in 4 steps and stagnates, then diverges, at 4 cells.
  The Newton matrix is consistent there too, so the linearised problem is
  nearly singular. Weak shear explains the first case only; candidates for
  the second are a near-marginal pressure-driven mode of the label Hessian
  (which is not the ideal δW, section 7) and the coarse planes. To be
  measured by the smallest eigenvalues of the reduced Hessian against ε.
- The ι of the sheared family (scan over ε, S, λ, δ): shear grows with S and
  δ and is independent of λ; at S = 1.4 it is below 2e-2 for every ε.

E6a (`test_eulerian_clebsch_free_boundary`, axisymmetric free boundary).
Formulation:

- Plasma: E1 inside Γ_p (χ = 1 on Γ_p, χ, g, X, μ), on an O-grid whose wall
  is Γ_p and whose core is centred on X (FortFEM `centre_core`), so X is
  always the central vertex. With the default core X lies on a cell edge
  (up-down symmetry), a boundary deformation moves it across, the polar
  rule switches cells and the residual is not differentiable in Γ_p.
- Shell: ψ in V0 of an O-grid with an annulus between Γ_p (interface) and
  the fixed Γ_c (wall); the functions inside Γ_p, interface included, are
  tied to one unknown ψ_b, so ψ is constant there and adds no energy, and
  continuity across Γ_p is exact. Energy ∫ |∇ψ|²/(2R) − fψ plus the
  exterior terms of `gs_spline_vacuum_setup` on Γ_c (Steklov–Poincaré map
  and coil data). B_φ = F_vac/R in the shell.
- Γ_p: R = (r0, 0) + r(θ)(cos θ, sin θ), r a Fourier series of 2M + 1
  coefficients. Rows: the pressure balance ∮ (p + |B_in|²/2) w_m =
  ∮ |B_out|²/2 w_m tested with the same modes; ψ_b: the current through
  Γ_p from both sides. Each side on its own representation of Γ_p.
- Newton: the trace rows differentiated by the chain rule, the X and Γ_p
  columns by central differences with both meshes rebuilt (ALE).
  Globalisation: Γ_p frozen first (E1 with the axis free), then all free;
  a monotone backtracking search in a trust region; mesh continuation (the
  previous level's solution projected onto the next mesh), because the
  axis basin shrinks with h (E3b).

Manufactured exact case: the Solov'ev field of E1 everywhere inside Γ_c,
Γ_p = {ψ = ψ_b}, Γ_c = {ψ = 2ψ_b}; the shell carries the known toroidal
current of a distributed coil (f = −(8α + 2β) R, Δ*ψ = (8α + 2β) R²), the
exterior sees ψ and ∂_nψ on Γ_c as coil data, p = 0 and B continuous on
Γ_p. Start: Γ_p radius, its m = 1 modes and ψ_b by 1 %, the axis by
(6e-3, 3e-3), the labels by 1.5 % χ(1 − χ). M = 10 (truncation of the
exact Γ_p 5e-5). Degree 2:

| cells | 6 | 12 |
|---|---|---|
| L2 error of B in the plasma | 2.59e-3 | 4.78e-4 |
| L2 error of ψ in the shell | 4.98e-4 | 1.06e-4 |
| max radial error of Γ_p | 6.88e-4 | 2.14e-4 |
| axis | 1.20e-2 | 1.41e-3 |
| Newton steps (stage 2, \|r\| < 5e-8) | 4 | 4 (from 6 cells) |

Rates: B 2.44, ψ 2.23, Γ_p 1.68, axis 3.09 (4 → 8 cells with M = 8: 3.13, 2.25,
1.98, 3.35; at 8 cells B is 1.2 times the E1 interpolation error). \|r\|
at 12 cells: 8.3e-1, 2.8e-2, 2.6e-4, 1.5e-7, 3.1e-8. The gate takes 2 min
(slow): the 2M + 3 difference columns rebuild both meshes twice each.

- The Newton matrix agrees with differences of the residual to 1e-8 or
  better in every block (argument `jacobian`).
- 2 and 4 cells are pre-asymptotic: the discrete equilibrium lies far from
  the exact one along the weakly restored rigid displacements (the
  smallest singular values at 2 cells belong to the vertical and the
  radial shift of plasma and axis), so Newton from the exact state needs
  4 cells and the gate starts at 6.
- FortFEM's `spline_interface_quadrature` finds no points on the interface
  next to the periodic annulus; the shell traces use the side quadrature of
  the interface tag, restricted to the annulus block.
- To do: analytic shape-derivative columns (fast gate), the plasma current
  as an input (I_p instead of G'), and a coil-only Γ_c (filament coils as
  in `gs_spline_free`).

E6b result (`test_eulerian_clebsch_free_boundary_3d`, the design below;
sheared family eps 0.6, S 2.4, delta 0.2, B_s the oracle field, Gamma_c the
oracle surface chi = 1.5, p = (delta - psi)/lambda^2; 4 cells, 9 planes,
M = 6 per plane, 3456 unknowns; start: Gamma_p radius by 1 % and its
m = 1 modes, the axis by 2e-3 helically, the labels by 1 %):

| | exact state on the mesh | start | Newton |
|---|---|---|---|
| L2 error of B in the plasma | 2.09e-2 (g fitted) | 3.31e-2 | 1.67e-2 |
| max radial error of Gamma_p | 5.87e-3 (truncation at M = 6) | 1.15e-2 | 5.65e-3 |
| axis | 0 | 2.0e-3 | 4.2e-3 |
| \|grad Phi\| relative to \|B_s\| (exact 0) | 0 | 0 | 9.2e-3 |

Stage 1 (Gamma_p frozen) 6 steps, stage 2 (free) 3 steps: \|r\| 3.1e-2,
2.1e-3, 1.7e-5, 8.2e-9. The Newton matrix agrees with differences of the
residual to 1e-8 or better in every block (argument `jacobian`, 2 cells).
Slow (4.5 min on 24 threads): the 2 nt + nt(2M + 1) difference columns each
evaluate the whole 3D residual.

- The Fourier series about the plane's axis converges slowly for these
  shaped sections (truncation 4e-2 at M = 3, 6e-3 at M = 6), and modes the
  boundary cells cannot carry (M above about a third of the 4n boundary
  cells) make the Newton matrix nearly singular (2 cells, M = 6: stage 2
  stalls). A refinement ladder needs M with n, or a boundary parametrised
  by the mesh itself.
- To do: analytic shape-derivative columns; the poloidal circulation of
  Phi (plasma current) and filament coils as B_s; the per-mode exterior map
  on Gamma_c instead of the conducting wall; the solver moved from the test
  into a module and an application.

E6b (3D, design). The plasma is E4 on per-plane O-grids whose wall is
Γ_p(θ, φ) and whose core follows the axis curve. The shell field is split,
B = B_s + ∇Φ: B_s a given divergence-free source field (the coil field by
Biot–Savart, or in the manufactured case the oracle field itself) and Φ
single valued in V0 of the shell section times the trigonometric nodal
basis in φ, with ΔΦ = 0, ∂_nΦ = −B_s·n on Γ_p (B·n = 0) and, first, a
conducting wall ∂_nΦ = 0 on Γ_c (the per-mode exterior map later). Γ_p
is a double Fourier series of normal offsets about the reference axis
curve, its rows the pressure balance tested with the same modes. With the
oracle field as B_s, Γ_c an oracle flux surface and the edge pressure
zero at ψ = δ, the exact solution has Φ = 0; a perturbed Γ_p sees a
nonzero vacuum response and must return.

FortFEM supplies the separatrix-aligned layouts for E2b and E5:
`layout_doublet` (18 blocks, X of valence 8) and `layout_island_chain` (any
m ≥ 2), V0 cut across `LAYOUT_SEPARATRIX`, one-sided traces and grading
towards the separatrix and X. For a log singularity in g, power grading with
exponent above 2(p + 1) is needed for order p + 1 in L2 of g; geometric
grading with a fixed ratio stalls. For B (E2b), power 4 in a graded layer
reaches order p at p = 2.

## 6. Order of work

1. Literature: done, section 3. Repeat the enrichment search with a
   working web search before any publication.
2. FortSym leaves for the pointwise terms of (3) and (5) and their first and
   second derivatives in (χ, g, profile coefficients), with identities.
3. E0, E1 in axisymmetry (done, `test_eulerian_clebsch_fit`,
   `test_eulerian_clebsch_newton`).
4. O-line handling for several axes: E2a and E2b done (graded layers of
   power 4 on FortFEM's doublet layout); E2 with Newton open (section 5,
   E2 attempt: the log-singular flux profile at the separatrix leaves
   Newton a basin of the size of the smallest |χ − 1| in the graded
   layer). The graded layer belongs in FortFEM's `spline_grade_knots` once
   E2 needs it too.
5. 3D: E3, E3b, E4a and E4 (Newton on the sheared family, ε = 0.6) done;
   ε = 1.2 open (above); then the registered integer case against
   VMEC++, DESC and GVEC.
6. Islands and rational surfaces, E5. Open whether Newton from nested
   initial labels finds the island topology, or whether the topology has to
   be seeded.
7. Free boundary, E6: E6a (axisymmetric) and E6b (3D, scalar-potential
   shell with a conducting wall, exact oracle case) done; next the plasma
   current and filament coils, the exterior map, analytic shape columns. Then the solver: Newton–Krylov with AMS
   for the V1 part and AMG for χ and g.

## 7. Open

- Regions and closures. Representation (labels (2), free A ∈ V1, vacuum)
  and physical closure (prescribed profiles, ideal force balance, Taylor
  relaxation, anisotropic transport, kinetic closures) are separate choices.
  Ideal force balance gives B·∇p = 0: where ∇p ≠ 0 the level sets of p are
  local magnetic surfaces, so ideal pressure gradients sit where labels are
  exact; on a fully ergodic volume p is constant and current continuity makes
  α in ∇×B = αB constant (Beltrami). The long-term model is one residual
  system with region blocks: labels for the core and island interiors, the
  Taylor closure (`nonnested_beltrami`) for fully relaxed volumes, NN2 for
  stochastic and open edge regions with finite transport, vacuum and coils
  through the free-boundary coupling. Interfaces are good magnetic surfaces
  with continuous tangential A, matched harmonic fluxes and total-pressure
  balance from the variation of the shared interface, plus a transport-flux
  law where NN2 meets them; at material targets sheath conditions replace
  B·n = 0. Regions are chosen explicitly or from a verified topology analysis
  between solves, never inside Newton. Order after E5: label + Taylor region
  on one movable interface, then label + NN2 (staggered, then monolithic),
  aimed at the W7-X island divertor (E9). Island size is an input in
  reconstruction and needs a Taylor, NN2 or resistive closure in prediction.
  NN2 with vanishing pressure gradient does not reduce to the helicity
  constrained Beltrami state by itself (different invariants); the same
  closed stochastic region solved with both, on matched boundary fields,
  calibrates when the Taylor closure suffices.
- Stability: the Newton Hessian of W in (χ, g) at fixed F, G and p(χ) is
  not the ideal δW (no compressional γp(∇·ξ)² term, not the space of
  displacements ξ); global stability of a label equilibrium goes through the
  ξ-based stability machinery until the relation is derived.
- Kinetics: the regional labels and θ* = θ + g/F' give the flux coordinates
  that a nested-surface drift-kinetic solver (K2) needs inside island
  interiors; orbits that cross separatrices or stochastic regions remain
  open. With anisotropic or flowing kinetic pressure, ∇p ≠ 0 no longer
  implies surfaces.
- Globalisation: the Newton basin for the axis is about 2e-3 of the minor
  radius; continuation in β, boundary and ΔΦ_I.
- Inf–sup of the δχ equation near O-points and saddles, where ∇χ = 0.
- Uniqueness of the representation (2) with several O-lines in one region.
