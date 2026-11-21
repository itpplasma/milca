# Free boundary on a moving mesh for both equilibrium lanes

Plan of 2026-11-20 for the free-boundary gates E4, E5 and E8 of roadmap item
P1b (3D equilibrium, [procedure](../ROADMAP.md#p1b-procedure)). MILCA has two
equilibrium lanes, and both get the same free boundary: a fixed outer
surface Γ_c with an exact exterior, the coil field split off, and a mesh
whose tangential motion is a gauge.

| Lane | Model | Topology | Comparators |
|---|---|---|---|
| Frozen lane ([ideal-map-3d.md](ideal-map-3d.md)) | ideal MHD, unknown the map of a reference mesh, fluxes given as functions of a label per region | nested, or islands whose topology the reference fixes | VMEC family, GVEC, DESC (nested), SIESTA ideal steps and SPEC with prescribed current (fixed-topology islands) |
| NN lane ([NN2](equilibrium-nonnested-nn2.md)) | vector potential, pressure, electric potential and loop voltage on a mesh fixed in space, closed by pressure transport and Ohm's law (or Taylor regions) | any, including islands that form and stochastic regions | SPEC and SPECTRE (Taylor regions), HINT and PIES (published results), SIESTA resistive steps |

VMEC is a comparator only for nested states of the frozen lane, because only
there both solve the same model. In 3D, nested surfaces with a smooth
pressure are not a generic solution of ideal MHD, and VMEC imposes them
through its coordinates. NN2 at finite resistivity opens islands there.

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
and the exterior couples on a fixed curve. Sections 2 and 3 apply this to the
NN lane, section 4 to the frozen lane.

## 2. NN lane: formulation

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

## 3. NN lane: gates

A1 to A3 run locally in seconds. A4 and A5 are single cluster jobs. G7 and G8
are the island and free-boundary gates of
[islands](island-equilibrium.md#verification-ladder-before-machines).

| Gate | Case | Compared | Passed when |
|---|---|---|---|
| A1 discretisation | 2D free-boundary Solov'ev equilibrium with exact solution | Lagrange P3 triangles of the existing Grad–Shafranov solver `gs_free` against C² cubic and degree-5 splines | splines reach the P3 error with fewer unknowns and less wall time, Newton converges quadratically in at most 8 steps |
| A2 anisotropy | anisotropic diffusion in a prescribed field with islands, χ_∥/χ_⊥ from 10⁶ to 10⁹ | P3 against degree-4 splines, each with and without mesh alignment | spurious perpendicular diffusion (pollution) below 10⁻³ of the exact temperature range at χ_∥/χ_⊥ = 10⁹ |
| A3 moving mesh | A1 with the mesh moved by (4) in the alternating loop | interface on Γ_p against cut cells | convergence order of the spline degree at Γ_p with the interface, at most 8 Newton steps per mesh update, no refused full step |
| A4 coupled 2D | nested tokamak with NN2 on the spline mesh and the exterior (3) | `gs_free`, FreeGSNKE | agreement with `gs_free` within the estimated discretisation error, wall time at most FreeGSNKE's at the same accuracy |
| A5 3D exterior | field of a known current inside a W7-X coil set, with axisymmetric and shaped Γ_c | Biot–Savart, virtual casing (the field outside a surface from B on it) | exterior error converges at the rate of the spline trace space |

E4 and E5 follow A4 and A5. E8 also needs G7 on the island cases of E6.

## 4. Frozen lane

The frozen lane stores b as a closed 2-form on the reference mesh,
b = (−∂₂χ, ∂₁χ, b₃) in master coordinates, with χ the reference label
(poloidal flux) and b₃ the toroidal-flux density. Today b₃ is sampled from
the seed point by point, so it is frozen per material element. Ideal motion
conserves only the toroidal flux between two χ surfaces,

$$
\Phi'(\chi) = \frac{d}{d\chi}\int_{\chi' < \chi} b_3 \, d^2s , \tag{5}
$$

for nested surfaces and irrational transform (Kruskal and Kulsrud). Every
other feature of b₃ along a surface is a relabelling. The equilibrium must
therefore slide material along the surfaces until b₃ is consistent with
the new geometry. That slide is the near-null of section 1.

The change is to keep only (5) frozen. The relabelling

$$
b \;\to\; b + \nabla_s f \times \nabla_s \chi \tag{6}
$$

with any scalar f(s₁, s₂, s₃) keeps χ and every flux (5), because
∇f × ∇χ is divergence-free and tangent to the χ surfaces. In axisymmetry
(6) changes only b₃, by −b_pol·∇f. With f a further unknown of the map,
a tangential node slide combined with the matching f leaves B unchanged.
The slide becomes an exact null, and a gauge fixes it: the transversal
gauge of [ideal-map-3d.md](ideal-map-3d.md) (each node moves along a fixed
direction across the surfaces), or the mesh equation (4) once the
transversals become too skewed, after which the map is frozen again on the
solved geometry. In 3D, f is the counterpart of VMEC's stream function λ
and of the toroidal displacement ω of the transversal gauge.

With islands the same holds per region, with χ the label of that region
and (5) per region. The separatrices are region boundaries of the mesh.
Regions without surfaces (stochastic) have no label and belong to the NN
lane.

The frozen lane has two vacuum variants, which must agree within the
discretisation error (gate B2b):

| Variant | Unknowns of the vacuum | Dense part | Shared with |
|---|---|---|---|
| V1 shell | A_p on shell cells between Γ_p and Γ_c with a harmonic mesh, the map 𝒮 of (3) on the fixed Γ_c | 𝒮, assembled and factored once per Γ_c | NN lane |
| V2 sheet | the potential φ of the surface current K = n × ∇_s Φ on Γ_p ([derivation](../theory/derivations/free-boundary-plasma-boundary-sheet-3d.md)), eliminated per state | the sheet self-energy on the moving Γ_p, by fast summation | none |

The sheet cannot serve the NN lane: it needs a closed Γ_p with B·n = 0
outside all plasma current, which islands and stochastic regions at the
edge do not provide.

In V1 the coil field never enters a shell cell. With J = J_p + J_c and A
the vector potential of all currents, the field energy is

$$
\frac12\int \mathbf J\cdot\mathbf A\,dV
= \frac12\int \mathbf J_p\cdot\mathbf A_p\,dV
+ \int \mathbf J_p\cdot\mathbf A_c\,dV + \text{const} . \tag{7}
$$

The first term needs only A_p (mesh inside Γ_c, 𝒮 outside), the coupling
term only A_c where J_p ≠ 0: in D and on Γ_p as the sheet current of the
pressure jump. Both are free of the coils, so coils and vessel may cross
the shell and Γ_c. As in the sheet derivation, (7) is evaluated through
coil fluxes across surfaces spanned inside D, without A_c.

Both variants evaluate B_coil only in D and on Γ_p, directly from the coil
filaments by Biot–Savart (no mgrid interpolation), and reject every Newton
step that brings Γ_p closer to a coil than a set distance.

The sheet self-energy is a pair sum over the N_b boundary points, O(N_b²)
per Hessian-vector product. Above a few thousand points it is summed in
O(N_b log N_b) by a fast multipole method (FMM) for the Laplace kernel with
charges and dipoles, which returns potential, gradient and Hessian, so the
same two passes give the gradient and the Hessian-vector product. The near
field keeps the corrected trapezoidal rule.

Gates of the frozen lane:

| Gate | Case | Compared | Passed when |
|---|---|---|---|
| B1 relabelling null | 2D free-boundary fixture `gs_free` (seconds) | b₃ frozen per element against (5) with f free | smallest eigenvalue of the gauged Newton matrix at least 10⁻³ of the largest at fixed mesh, Newton converges quadratically with no refused full step |
| B2 coil split | B1 with Γ_c crossing a coil | Γ_c inside the coils | same solution within the discretisation error |
| B3 DIII-D coarse | 48 boundary nodes, one level | 65 iterations with the watchdog today | at most 10 Newton steps, same energy |
| B4 3D nested | Landreman–Paul QA free boundary, coarse | VMEC++ on the same coils and profiles | E4 observables within the error estimate |
| B2b two vacua | DIII-D (2D) and Landreman–Paul QA (3D), coarse | V1 shell against V2 sheet | same equilibrium within the discretisation error; wall time per variant at matched accuracy |
| B6 fast sheet | sheet self-energy on the W7-X boundary, 10³ to 10⁵ points | FMM against direct pair sum | energy, gradient and Hessian-vector product agree to 10⁻¹⁰ relative; cost O(N log N); faster than direct above 5·10³ points |
| B5 fixed-topology islands | classical stellarator with 1/2 islands | SIESTA ideal steps, SPEC with prescribed current | E7 observables |

## 5. Work by repository

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
| R1 flux-function form (5) with relabelling (6) in the frozen map, 2D then 3D | MILCA, kernels by fortsym | b per element, toroidal displacement ω, slide gauge | f as unknown, gauge, refreezing on the solved geometry |
| R2 coil split and fixed Γ_c in the frozen-lane free boundary (V1) | MILCA | shell and per-mode exterior on Γ_c | coupling (7), B_coil from C, Γ_c free of the coils |
| R3 sheet vacuum (V2) | MILCA, kernels by fortsym | 2D on Γ_p (axisymmetric single layer), 3D by the corrected trapezoidal rule | B_coil from the filaments, coil-distance guard, FMM far field |
| F6 Laplace FMM with charges and dipoles, potential to Hessian, nfp symmetry | FortFEM | none | far field of the sheet self-energy, gate B6 |
| D dense BLAS and LAPACK from fortnum only, at OpenBLAS speed | fortnum, FortFEM | fortnum_lapack77 under the standard names; GEMM and SYRK at OpenBLAS speed, TRSM and Cholesky within 1.5× | FortFEM without system BLAS; threaded LU and LDLᵀ panels, two-stage symmetric eigensolver |

## 6. Order

Every step iterates on cases that run in seconds and at most two minutes.
Larger cases follow only once the step is complete, in rising size, up to
the decks of the [equilibrium benchmark](https://github.com/itpplasma/benchmark_vmec).

1. B1 with R1 in 2D. Then B2 with R2, B2b and B3, all on the DIII-D case
   against VMEC2000. D in parallel.
2. A1 and A2 on the existing rectangular-patch complex and `gs_free`. P only
   if A2 fails.
3. F1, K and F4, then M2 and A3. F5 only if A3 fails.
4. R1 and R3 in 3D, F6 and B6, then B4 with both vacua. C, F2, S and M1, with F3 for the axisymmetric Γ_c,
   then A4.
5. F3 for a shaped Γ_c, then A5.
6. B5 and M3, then E4, E5 and, after G7, E8.
