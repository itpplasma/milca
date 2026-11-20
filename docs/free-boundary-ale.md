# Free boundary on a moving spline mesh with Eulerian fields

Plan of 2026-11-20 for the free-boundary gates E4, E5 and E8 of roadmap item
P1b (3D equilibrium, [procedure](../ROADMAP.md#p1b-procedure)). It makes the
non-nested field of [NN2](equilibrium-nonnested-nn2.md) the free-boundary
path for nested and island cases. NN2 solves for a vector potential, the
pressure, an electric potential and the loop voltage on a mesh that is fixed
in space, with pressure transport and Ohm's law as closure.

The frozen-flux map ([ideal-map-3d.md](ideal-map-3d.md)) stays the path for
fixed-boundary nested cases and for ideal stability (P1e). The free-boundary
version of that map ([ideal-free-boundary.md](ideal-free-boundary.md),
[free-boundary-nested-map-3d.md](free-boundary-nested-map-3d.md)) is no
longer developed.

## 1. Reason

The frozen-flux map takes the mesh node positions as unknowns and minimises
the plasma energy W. Each cell carries a fixed magnetic flux (the closed
2-form b), and so do the nodes of the plasma boundary Γ_p. A slide of these
nodes along Γ_p moves flux and is therefore physical, but it changes W only
at second order.

- DIII-D, coarsest mesh of the resolution ladder (48 boundary nodes): the
  smallest eigenvalues σ of the Newton matrix (the Hessian of W with the
  vacuum coupling) belong to tangential motions of the boundary nodes, with
  σ/σ_max = 4.5·10⁻⁶ to 6·10⁻⁵.
- From seed to solution the boundary nodes slide up to 0.52 m along Γ_p and
  move 0.15 m across it.
- Of 99 Newton iterations, about 60 refuse the full step, because W rises by
  about 3·10⁻⁷ along it. A watchdog acceptance (Nocedal and Wright §15.6)
  reduces the total to 65 iterations but does not remove the effect.

For a functional F of a domain Ω with boundary ∂Ω, outward normal n and a
boundary velocity field w, the first variation depends only on the normal
velocity (structure theorem of Hadamard and Zolésio),

$$
dF(\Omega)[\mathbf w] = \int_{\partial\Omega} h\,(\mathbf w\cdot\mathbf n)\,dS , \tag{1}
$$

with h a scalar density on ∂Ω. A tangential motion of boundary nodes is
therefore free unless the mesh carries physical data, as the frozen-flux map
does. CEDRES++ and FEEQS.m avoid this by storing no physics in the mesh
motion: the flux ψ lives on a fixed mesh, the boundary is a level set of ψ,
and the exterior couples on a fixed curve. The plan below follows the same
principle with NN2's fields, which also admit islands.

## 2. Formulation

Every field is Eulerian, i.e. a function of position. The mesh may move, but
by its own equation, and a tangential mesh motion changes no physical
quantity.

**Field split.** With B_coil the Biot–Savart field of the coils (including
the toroidal-field coils) and A_p the vector potential of the plasma current
J_p,

$$
\mathbf B = \mathbf B_\text{coil} + \nabla\times\mathbf A_p . \tag{2}
$$

**Interior.** Inside a fixed toroidal surface Γ_c that encloses J_p, the
unknowns A_p, the pressure p, the electric potential Φ and the loop voltage
U satisfy the NN2 equations (Ampère's law, anisotropic pressure transport,
current continuity and the loop condition) with B from (2). Coils and vessel
may cross Γ_c and the mesh, because B_coil is never discretised.

**Exterior.** Outside Γ_c the field of J_p is current-free and decays, so on
Γ_c

$$
\mathbf n\times\nabla\times\mathbf A_p = \mathcal S\,[\mathbf A_p] , \tag{3}
$$

with 𝒮 the exterior Dirichlet-to-Neumann map (DtN) of Γ_c, computed by a
boundary-element method (BEM). Γ_c is fixed, so 𝒮 is assembled and factored
once per geometry and reused in every Newton step and every parameter scan.
For an axisymmetric Γ_c, 𝒮 splits into independent blocks S_n, one per
toroidal mode number n. Γ_p does not enter the exterior problem.

**Mesh.** In each plane φ = const the mesh is a multipatch spline map
X(ξ, φ) = Σ_i c_i(φ) N_i(ξ) of a fixed 2D parameter domain, with spline
basis functions N_i and control points c_i that are Fourier series in φ. The
patch layout has no polar point, so no coordinate is singular at the
magnetic axis. The c_i minimise

$$
M[\mathbf c] = \int \Big( \frac{\operatorname{tr} g}{\sqrt{\det g}}
+ \kappa\, a\,(\boldsymbol\tau\cdot\nabla p)^2 \Big)\, d^2\xi\, d\varphi , \tag{4}
$$

with g the 2×2 metric of X in the plane, τ the unit tangent of the patch
lines that should follow pressure surfaces, a ≥ 0 a weight that vanishes
where p has no closed level sets, and κ ≥ 0. The first term is the Winslow
smoothing functional. Where the plasma edge (the last closed surface or the
limiter contour) exists, one patch interface is constrained to lie on it.
Without this constraint, the jumps of J_p and ∇p at Γ_p fall inside
elements and are integrated exactly on the cells cut by Γ_p.

**Coupling.** At first the field Newton with a fixed mesh and the mesh
update alternate, with the fields transferred by L² projection between the
spline spaces of the old and the new mesh. A monolithic Newton that includes
the mesh sensitivity follows only once the alternating loop passes gate A3.

**Profiles.** Mode B of [islands](island-equilibrium.md): sources and
conductivities are fitted on the seed so that target profiles of p and of
μ₀J·B/B² are reproduced, then frozen, and the target and fit defects are
reported.

## 3. Gates

A1 to A3 run locally in seconds. A4 and A5 are single cluster jobs. G7 and G8
are the island and free-boundary gates of
[islands](island-equilibrium.md#verification-ladder-before-machines).

| Gate | Case | Compared | Passed when |
|---|---|---|---|
| A1 discretisation | 2D free-boundary Solov'ev equilibrium with exact solution | Lagrange P3 triangles of the existing Grad–Shafranov solver `gs_free` against C² cubic and degree-5 splines | splines reach the P3 error with fewer unknowns and less wall time, Newton converges quadratically in at most 8 steps |
| A2 anisotropy | anisotropic diffusion in a prescribed field with islands, χ_∥/χ_⊥ from 10⁶ to 10⁹ | P3 against degree-4 splines, each with and without mesh alignment | spurious perpendicular diffusion (pollution) below 10⁻³ of the exact temperature range at χ_∥/χ_⊥ = 10⁹ |
| A3 moving mesh | A1 with the mesh moved by (4) in the alternating loop | interface on Γ_p against cut cells | convergence order of the spline degree at Γ_p with the interface, at most 8 Newton steps per mesh update, no refused full step |
| A4 coupled 2D | nested tokamak with NN2 on the spline mesh and the exterior (3) | `gs_free`, FreeGSNKE | agreement with `gs_free` within the estimated discretisation error, wall time at most FreeGSNKE's at the same accuracy |
| A5 3D exterior | field of a known current inside a stellarator coil set, with axisymmetric and shaped Γ_c | Biot–Savart, virtual casing (the field outside a surface from B on it) | exterior error converges at the rate of the spline trace space |

E4 and E5 follow A4 and A5. E8 also needs G7 on the island cases of E6.

## 4. Work by repository

| Item | Repository | Exists | To do |
|---|---|---|---|
| F1 spline curl complex (spaces for A, B and p that commute with grad, curl, div) on a mapped multipatch meridional domain, Fourier in φ | FortFEM | per-mode complex on one rectangular patch with exact per-mode exterior (`fortfem_meridian_bspline_n`), multipatch spline maps in 2D and 3D, spline extraction at a polar point | mapped geometry and several patches in the per-mode complex, patch layout without polar point |
| F2 matrix-free operator actions by sum factorisation (one 1D sweep per direction, cost O(p) per unknown at degree p) | FortFEM | not checked yet | actions on tensor patches with a low-order preconditioner |
| F3 DtN (3) in the spline trace space | FortFEM | per-mode maps on an axisymmetric curve (`toroidal_mode_boundary_operators_2d`), 3D Laplace and Maxwell BEM on curved tori, curved Maxwell DtN, virtual casing | spline traces on a curved axisymmetric Γ_c, then on a non-axisymmetric Γ_c |
| F4 minimisation of (4) with a patch interface on a level set | FortFEM | none | functional, constraint and Newton on c |
| F5 cut cells on spline elements | FortFEM | level-set interfaces on triangles and tetrahedra | only if A3 fails with the interface constraint |
| K pointwise kernels of (4) with gradient and Hessian | fortsym, generated into FortFEM | NN2 kernels `coupled_current`, `anisotropic_flux` | Winslow and alignment densities |
| C coil field B_coil and its trace on Γ_c | MILCA | `coils_3d`, `coil_quadrature_field`, reader of tabulated coil fields (mgrid) | evaluation in (2) and on Γ_c |
| S per-mode block preconditioner on spline spaces | fortsparse, fortnum | sparse direct LU, FGMRES, Newton-MR | blocks per toroidal mode |
| P parallel operator by field-line tracing between planes (flux-coordinate independent, FCI) | fortnum, MILCA | ODE integrators | only if A2 fails with splines |
| M1 NN2 on F1 spaces with (2) and (3) | MILCA | NN2 on Lagrange triangles times Fourier | residuals on the new spaces, field split, exterior coupling |
| M2 alternating mesh loop | MILCA | boundary refit loop of `gs_free` | loop for spline maps and NN2 fields |
| M3 profile mode B | MILCA | design | fit, freeze, defect report |

## 5. Order

1. A1 and A2 on the existing rectangular-patch complex and `gs_free`. P only
   if A2 fails.
2. F1, K and F4, then M2 and A3. F5 only if A3 fails.
3. C, F2, S and M1, with F3 for the axisymmetric Γ_c, then A4.
4. F3 for a shaped Γ_c, then A5.
5. M3, then E4, E5 and, after G7, E8.
