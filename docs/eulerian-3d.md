# Eulerian 3D equilibrium: vector potential, labels and current potential

Question: which discrete formulation solves 3D ideal-MHD equilibria with
nested or several regions in Eulerian unknowns on FortFEM's spline de Rham
complex? It must have no gauge in the physics and no material map
(decision in [free-boundary-ale.md](free-boundary-ale.md), section 0, ladder
steps L5 and L6).

## 1. Model

Symbols:
- Ω: the domain inside the wall Γ_c.
- A: the vector potential, a 1-form, with B = ∇×A.
- χ: a label per region (a 0-form) whose level sets are the flux surfaces.
- p(χ): the pressure as a function of the label.
- λ: the multiplier of the constraint B·∇χ = 0.

The potential energy and the Lagrangian:

    W[A, χ] = ∫_Ω |∇×A|²/(2μ0) dV − ∫_Ω p(χ) dV + W_wall[n×A on Γ_c]   (1)

    L[A, χ, λ] = W + ∫_Ω λ B·∇χ dV.                                     (2)

## 2. Stationarity

The variation in A uses ∫ λ (∇×A)·∇χ = ∫ A·∇×(λ∇χ) + boundary terms, and
∇×(λ∇χ) = ∇λ×∇χ:

    ∇×B/μ0 = J = ∇χ×∇λ          (variation in A, inside Ω)             (3)

The variation in χ uses ∫ λ B·∇δχ = −∫ δχ ∇·(λB) and ∇·B = 0:

    B·∇λ = −p'(χ)                (variation in χ)                      (4)

    B·∇χ = 0                     (variation in λ)                      (5)

From (3) to (5):

    J×B = (∇χ×∇λ)×B = ∇λ (B·∇χ) − ∇χ (B·∇λ) = p'(χ) ∇χ = ∇p.           (6)

So the stationary points are ideal-MHD equilibria with nested surfaces
χ = const. λ is the current potential: (3) is the Clebsch form of J, and (4)
is the magnetic differential equation of the Pfirsch–Schlüter current.

## 3. Net currents

With single-valued λ, J = ∇×(χ∇λ), so the toroidal current inside the
surface χ vanishes:

    I(χ) = ∮_{∂S(χ)} χ ∇λ·dl = χ ∮ ∇λ·dl = 0.                           (7)

The net currents come from multivalued parts of λ, one per angle:

    J = ∇χ×∇λ̃ + ∇×(C(χ) h_θ) + ∇×(G(χ) h_φ),   λ̃ single valued.       (8)

- h_φ = ∇φ is closed on the whole torus. It is FortFEM's toroidal harmonic
  class (`spline_derham_harmonic`), and G(χ) gives the poloidal current.
- h_θ = ∇θ is the winding form about the magnetic axis X. Poloidal loops
  are contractible in the solid torus, so the complex has no such class:
  ∇×h_θ = 2π δ_X. Then

      ∇×(C h_θ) = C'(χ) ∇χ×h_θ + 2π C(χ(X)) δ_X,                       (9)

  and the line current at X vanishes only for C(χ_axis) = 0, so X must be
  the magnetic axis. X is an unknown per toroidal plane, fixed by Newton
  together with the axis, as in the 2D solver. The toroidal current is
  I(χ) = 2π C(χ).

Both terms enter weakly, ∫ A·∇×(C h) dV = ∫ B·C(χ) h dV (plus a wall
term), with h evaluated pointwise. h_θ ~ 1/ρ at the axis is integrable.
A gradient added to h changes the load by −∫ f ∇·(C B) dV =
−∫ f C' B·∇χ dV, which vanishes under (5), so only the winding matters.

Current mode (VMEC's I(s)): C(χ) is the given profile. ι mode: C follows
from ι(χ) as an unknown profile with the poloidal flux as its constraint.

## 4. Discrete spaces and the weak constraint

- A ∈ V1, ungauged. B = ∇×A ∈ V2 is divergence free exactly (d∘d = 0).
- χ ∈ V0 and λ̃ ∈ V0 of one degree lower (section 5), one pair per region. C(χ) and G(χ) are splines in χ,
  and X is a curve sampled at the toroidal nodes.
- The constraint holds weakly: ∫ ν B·∇χ dV = 0 for all ν ∈ V0.
- Label normalisation: (1)–(5) are invariant under χ → g(χ) with p and C
  carried along, so a profile constraint fixes the label, e.g. the poloidal
  flux ψ_p(χ) = χ ψ_p,edge. One scalar per region is not enough.
- Discrete flux functions of χ_h ∈ V0 are only approximately in V0. With
  an equal-degree multiplier ν ∈ V0 the constraint locks for a profile that
  is not affine in χ_h. The multiplier space is therefore one degree lower
  than χ (section 5, G2a).
- Fixed values or constraints on the harmonic coefficients of A: the
  toroidal flux and, by mode, the linked poloidal flux.

## 5. Inf–sup and rational surfaces

The constraint operator χ ↦ (ν ↦ ∫ ν B·∇χ) has the flux functions as its
kernel. On a rational surface it also has the field-line labels, so the
magnetic differential operator B·∇ has no bounded inverse there: (4) has
the Pfirsch–Schlüter singularity ~ 1/x at a resonance. The discrete inf–sup
constant therefore falls with resolution near rational surfaces. This is
the physics, not a discretisation defect. Options, to be measured, not
assumed:

1. **Axisymmetry and quasi-symmetry:** in axisymmetry the resonant parts
   vanish. Gate G2 below checks that the constant stays bounded there.
2. **Regularisation:** ε ∇⊥·∇⊥λ added to (4), with the limit ε → 0
   controlled by the mesh. This is a resistive-layer model, and its convergence
   must be shown.
3. **Stepped pressure:** p' = 0 inside regions with pressure jumps at
   interfaces (SPEC). Then (4) is homogeneous, and the jumps carry sheet
   currents with [[p + B²/2μ0]] = 0, as in the edge-model decision.
4. **Islands:** opening resonant surfaces into islands (several regions,
   ladder L6), where the topology allows it.

VMEC hides this by truncating Fourier modes. MILCA has to state it.

Measured in axisymmetry (G2a, labels fixed, A and λ̃ unknown, −ε I on the
λ̃ block relative to max|K|, A and χ of degree 2, error of B relative to
the interpolant's):

| label | multiplier degree | ε | 2 cells | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| χ = ψ_N | 2 | 1e-14 | 1.00 | 1.00 | 1.00 | 1.00 |
| χ = ψ_N | 2 | 1e-10 | 1.00 | 1.00 | 1.00 | 2.10 |
| χ = ψ_N | 2 | 1e-8 | 1.00 | 1.00 | 1.91 | 38.6 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 2 | 1e-16 | 6.06 | 30.6 | 88.4 | 11.1 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 2 | 1e-14 | 6.06 | 30.3 | 10.4 | 1.88 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 2 | 1e-12 | 6.06 | 16.7 | 1.93 | 1.24 |
| χ = ψ_N | 1 | 0 | 1.04 | 1.01 | 1.00 | 1.00 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 1 | 0 | 1.18 | 1.11 | 1.07 | 1.04 |
| χ = (ψ_N + ψ_N²/2)/1.5 | 1 | 1e-8 | 1.17 | 1.09 | 1.37 | 31.9 |

- Equal degrees: D^T has a near kernel of discrete flux functions, so
  ε > 0 is needed, and a fixed ε pollutes B more on finer meshes. For the
  nonlinear label the exact field violates the discrete constraint by
  |D a_I|/(|D||a_I|) = 3.3e-4, 1.6e-5, 4.9e-7, 1.7e-8, and the solution
  locks until that falls below ε.
- Multiplier one degree lower: stable without regularisation (only the
  constants remain in the kernel and one node is pinned), no locking for
  either label, and B within 18 % of the interpolation error from the
  coarsest mesh on. A penalty still pollutes and is not used.
- This matters for 3D, where both Φ(χ) and ψ_p(χ) are constrained and one
  of them is nonlinear in any label. Whether the lower-degree multiplier
  stays stable without axisymmetry and near rational surfaces is the open
  part of the inf–sup question (G2, G3).

With the labels as unknowns neither multiplier degree is stable and free of
locking. The plan that replaces the multiplier by a direct Clebsch
representation is [eulerian-label-plan.md](eulerian-label-plan.md).

## 6. Wall and exterior

W_wall depends on A only through B·n = ∇_s×(n×A) on Γ_c (the surface curl
of FortFEM's `spline_derham_surface`), so it is gauge invariant. The
exterior vacuum enters through a 3D Steklov–Poincaré or BEM operator on
B·n. The coil split enters through the coil fluxes (section 2 of the plan),
and the fast far field (FMM) is needed there. The surface has two loop
classes (toroidal, poloidal), and their fluxes are the two numbers of the wall.

## 7. Solver

- **Production:** Newton–Krylov on the saddle-point system in (A, χ, λ̃, X)
  and the scalar multipliers.
  - A stays ungauged: the energy depends on A only through B, so every Newton
    right-hand side lies in range(∇×ᵀ) and Krylov iterations converge there.
  - Preconditioner: AMS (Hiptmair–Xu) for the curl-curl block, AMG for χ and
    λ̃, and a Schur complement approximation for the constraint block.
  - The harmonic directions are kept out of the Krylov space and enter as
    constraints.
- **Small verification cases:** a direct solve with the tree gauge
  (`spline_derham_tree`).

## 8. Gates (each runs in seconds)

| Gate | Case | Reference | Checks |
|---|---|---|---|
| G1 | vacuum torus, n×A fixed on Γ_c, prescribed harmonic fluxes | B = B0 R0/R φ̂ and a Dommaschk field | B to the projection error, AMS iterations independent of the mesh |
| G2a | axisymmetric Solov'ev, labels fixed to χ_h = Π0 ψ_N, A and λ̃ unknown | exact field | B to the interpolation error, (4) satisfied weakly without being imposed |
| G2 | axisymmetric Solov'ev, n = 0 only, through (1)–(5) | the spline Grad–Shafranov solver | identical ψ and q, inf–sup estimate bounded |
| G3 | 3D low-β stellarator with fixed boundary | VMEC++ and GVEC | E4 observables |
| L5 | Landreman–Paul QA, free boundary | VMEC++ | as in the plan |

G1 result (direct solve, tree gauge grown from the wall nodes, n×A of the
exact field on Γ_c, A = (−B0 R0 ln R + 2ε x y) e_z, R0 = 3, a = 1, ε = 0.1,
M = 2):

| degree | cells per block side | rel. L2 error of B, Galerkin | interpolant | rate |
|---|---|---|---|---|
| 2 | 2 / 4 / 8 | 1.39e-2 / 3.38e-3 / 8.16e-4 | 1.40e-2 / 3.39e-3 / 8.16e-4 | 2.04, 2.05 |
| 3 | 2 / 4 / 8 | 3.12e-3 / 3.55e-4 / 4.33e-5 | 3.15e-3 / 3.78e-4 / 4.48e-5 | 3.14, 3.03 |

The AMS part of G1 waits for the fortsparse preconditioner.

G2a result (ψ = α(R² − R0²)² + βR²Z², α = 1/4, β = 1, R0 = 3, wall
ψ = 6.25, μ0 = 1, F = 3, A and χ of degree 2, λ̃ and ν of degree 1, no
regularisation, the toroidal current through C(χ) h_θ about the exact
axis):

| cells per block side | 2 | 4 | 8 | 16 |
|---|---|---|---|---|
| rel. L2 error of B, χ = ψ_N | 2.28e-2 | 4.41e-3 | 9.03e-4 | 2.12e-4 |
| rel. L2 error of B, χ = (ψ_N + ψ_N²/2)/1.5 | 2.58e-2 | 4.85e-3 | 9.64e-4 | 2.20e-4 |
| interpolant | 2.20e-2 | 4.37e-3 | 9.02e-4 | 2.12e-4 |
| weak residual of (4), χ = ψ_N | 1.3e-2 | 1.8e-3 | 3.0e-4 | 6.9e-5 |
| weak residual of (4), nonlinear label | 3.0e-2 | 6.1e-3 | 1.9e-3 | 6.1e-4 |

The residual of (4) is tested with interior ν of degree 2 and taken
relative to the p' term.

## 9. Open

- The size of the inf–sup degradation for smooth p on a QA case, and which
  option of section 5 to choose. G2 and G3 answer this.
- Moving surfaces in 3D: Γ_p as a mesh surface needs FortFEM's 3D chart
  with a moving interface. Until then the plasma boundary is the level
  χ = 1, not a mesh surface, and the interface integrals need a level-set
  quadrature.
- Several regions: χ jumps at interfaces, with a multiplier per region
  boundary.
