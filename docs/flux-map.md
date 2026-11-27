# Cartesian flux map: the ideal model as a map into the flux disc

One of the two discretisations of the ideal model in
[equilibrium-roadmap.md](equilibrium-roadmap.md); the other is its inverse,
the Lagrangian map. The unknown is the map from physical space to a model
torus in which the field is fixed by the profiles alone. Its polar form,
the pair (χ, g) with axis unknowns of
[eulerian-label-plan.md](eulerian-label-plan.md), is superseded. What is
specific to this direction of the map and what is not:
[audit](equilibrium-lane-audit.md), section 4.

**Name: CEMV**, *Code d'Équilibre par Moments Variationnels*, VMEC à
l'envers. VMEC (Variational Moments Equilibrium Code) solves for the map
from flux coordinates to physical space; CEMV solves for the map from
physical space to flux coordinates. The French word order reverses the
acronym as the method reverses the map, and since VMEC's direction is
traditionally called the inverse coordinate representation, CEMV is the
inverse of the inverse. "Moments" is kept for the pun: the section is
discretised by splines on a disc O-grid, only the toroidal direction by a
trigonometric basis. Comparisons against VMEC are written "CEMV vs VMEC".

## 1. Formulation

Model torus: the unit disc times a circle, with Cartesian disc coordinates
(x, y), s = x² + y² the normalised flux label and ζ the toroidal angle. On it
the profiles define a closed 2-form,

    β = (Φ'(s)/π) dx∧dy + ψ'(s) dζ∧ds,

Φ the toroidal flux inside s and 2πψ the poloidal flux (ψ' = ι Φ'/(2π)). The
unknown is the map u = (x, y) from the plasma volume Ω to the disc; ζ is a
given closed form of the toroidal class, not an unknown. The field is the
pullback of β,

    B = (Φ'(s)/π) ∇x×∇y + ψ'(s) ∇ζ×∇s,   s = x² + y²,                  (1)

    W[u] = ∫_Ω |B|²/2 dV − ∫_Ω p(s) dV.                                  (2)

By construction, for every u:

- ∇·B = 0 (the pullback of a closed form is closed);
- B·∇s = 0, so the level sets of s are magnetic surfaces;
- the fluxes and ι(s) are those of β, provided u maps the wall onto the unit
  circle with winding number one;
- the magnetic axis is the zero of u. It is a result, not an unknown, and the
  mesh need not know where it is.

(1) is polynomial in u and ∇u apart from the profiles: no winding form, no
1/ρ terms and no special quadrature. VMEC, GVEC and DESC solve for the
inverse map (flux coordinates → space); (1) is the direct map on a mesh in
physical space.

Relation to the label lane. With polar coordinates of the disc, x = √s cos θ*,
y = √s sin θ*, (1) is B = (Φ'/2π) ∇s×∇θ* + ψ' ∇ζ×∇s: the Clebsch form of the
label plan with χ = s and θ* = θ − g/F' the straight-field-line angle
(Kruskal and Kulsrud 1958, as cited there). (χ, g) are the polar coordinates
of the same map, and the constructions of the label lane are the artefacts
of the polar singularity: the axis X as an unknown with multipliers for
χ(X) = 0 and ∇χ(X) = 0, the winding form about X and its polar quadrature,
the frame H^(1/2) (the linearisation of u at the axis), the adapted forms
θ ← θ − g/F' (θ* is the unknown here), and the mesh tied to X.

Gauge. A rotation of each flux circle, θ* → θ* + k(s), leaves β and B
unchanged (the gauge g → g + k(χ) of the label lane). A further freedom is
the choice of ζ, which is fixed by giving ζ. Of this family the discrete
energy keeps one exact symmetry, the rigid rotation k = const; it is removed
by a moment constraint (section 4, "Gauge constraint, one solver loop and
continuation").

Variations. δu = −ξ·∇u for a displacement ξ; the component of ξ along B
drops out, so δu and ξ_⊥ correspond one to one, and stable equilibria are
minima of W at fixed profiles (pressure as a function of the label). The
Hessian of W in u is therefore expected to be the ideal δW at γ = 0. This is
not derived or tested yet.

Axisymmetry (gate F1). With ζ = φ and (x, y) functions of (R, Z):

    B_φ = τ(s) (x_R y_Z − x_Z y_R),  τ = Φ'/π,
    B_pol = ψ'(s) ∇s×∇φ,
    W = 2π ∫ [B_φ²/2 + |B_pol|²/2 − p(s)] R dR dZ.

## 2. Wall condition

On the wall |u| = 1. Imposed at nodes only, the image of the wall cuts
chords of the unit circle, the toroidal flux inside the image is smaller
than Φ, and it shrinks further whenever that lowers the energy: in the
prototype the map collapsed to a point of the circle (B = 0, W = 0) along a
path of monotonically decreasing energy.

The wall therefore carries an angle Θ, and u = (cos Θ_h, sin Θ_h) holds
along the whole wall. With N_b the V0 functions with a trace on the wall and
σ = Σ_b N_b,

    u_h = Σ_interior N_i u_i + σ (cos Θ_h, sin Θ_h),  Θ_h = Σ_b N_b Θ_b/σ.

On the wall σ = 1 and the interior functions vanish. The toroidal flux is
then ∫ τ(s) dx dy over the disc for every state: in gate F1 it changes by
1e-10 to 1e-15 between the start and the solution.

## 3. Guide field: locking and its cure

B_φ is τ times the Jacobian determinant of u, so the energy has the
structure of a nearly incompressible material: a bulk modulus of order B_φ²
on det ∇u against a stiffness of order B_pol² on ∇s. The angle part of the
map (θ*) is controlled by B_φ alone. With B_φ evaluated pointwise, the
discrete space cannot satisfy R τ(s) det ∇u_h = f(s) at every quadrature
point, and the defect is paid by the label: the error of B_pol grows with
B_φ/B_pol (the guide-field obstruction of the label plan, section 4b).

Cure: the energy uses a projection of B_φ onto a smaller space that keeps
the flux. Per cell, with the area measure dA (not R dA: the flux ∫ B_φ dA
must be reproduced, so constants must be in the test space of that measure),

    ⟨b̄, q⟩_dA = ⟨τ det ∇u_h, q⟩_dA  for all q in Q,   W_tor = ∫ R b̄²/2 dA.

Which Q works depends on the spaces (measured, section 4):

- C0 Lagrange Q_p elements (prototype): Q = Q_(p−1) per element removes the
  amplification completely. P_(p−1) with the R dA measure loses flux and
  oscillates wildly.
- FortFEM splines, C^(p−1) inside a block (gate F1): Q_(p−1) per cell still
  locks at the strongest field; the cell mean (Q_0) removes the
  amplification at p = 2 and p = 3. The count explains it: smooth splines
  have about one function per cell and component, Q_m has (m + 1)²
  constraints per cell.

The cell mean makes b̄ itself first order. The pointwise τ det ∇u_h of the
same solution stays accurate and is the field that is reported. The
order-consistent choice for splines is the global projection onto the
splines of degree p − 1 (V2 of the complex, one function per cell), which
needs a mixed system; it is the candidate for 3D (open, section 6).

A weaker projection space frees spurious modes of det ∇u_h (folded maps of
low energy). Far from the solution at a strong guide field the minimiser can
fall into them; continuation in the guide field from B_φ/B_pol of order one
avoids this (section 4).

## 4. Results

### Gate F1 (`test_flux_map_axisym`, `flux_map_axisym_gate_f1`, 15 s)

Module `milca_flux_map`, kernels `flux_map_toroidal` and `flux_map_poloidal`
(FortSym, `gen/gen_flux_map.f90`, with proofs of tangency, of the rotation
invariance of both kernels and of Hessian symmetry). Solov'ev case and
profiles of gate E1 (ψ(s) = ψ_b s, p linear, τ = −2F' from the toroidal
flux), degree 2, FortFEM O-grid centred on the middle of the wall (15 % of
the minor radius from the axis) and fixed. Start: the polar map about the
mesh centre from the control points, u = (x − c)/ρ_wall(angle). Newton with
an identity shift, accepted on energy decrease. The numbers of this
subsection were measured on the mesh centred on the midplane with the gauge
by up-down symmetry (penalty rows y_h = 0 on Z = 0 and Θ_h = 0 at the
outboard wall node, removed 2026-11-23). The gate now runs on a mesh without
symmetry (centre 10 % of the minor radius above the midplane), with the
moment constraint of the rotation and from a cold start at R B_φ = 300
through the continuation driver; its numbers are in "Gauge constraint, one
solver loop and continuation" below.

Pointwise B_φ, R B_φ = 3 (relative L2 errors, s in the maximum norm, axis
in metres):

| cells | Newton steps | B | B_pol | projected exact ψ_N on this mesh | s | axis |
|---|---|---|---|---|---|---|
| 2 | 14 | 2.42e-2 | 2.35e-2 | 2.25e-2 | 3.5e-2 | 5.0e-3 |
| 4 | 18 | 7.11e-3 | 6.61e-3 | 4.44e-3 | 3.3e-3 | 1.2e-4 |
| 8 | 23 | 1.93e-3 | 1.89e-3 | 9.08e-4 | 4.0e-4 | 1.3e-5 |

Rates: B 1.77, 1.88; s 3.4, 3.1; axis 5.4, 3.2. B_pol is 1.04, 1.49 and 2.08
times the B_pol of the projected exact label on the same mesh (open: the
constant at 8 cells; more quadrature points change nothing).

Guide field, B_pol error at 2, 4, 8 cells (R B_φ = 3, 30, 300 is
B_φ/B_pol ≈ 1.5, 15, 150; continuation 3 → 30 → 300 takes 17–41 steps for
the first stage and 3–6 for each further one):

| B_φ in the energy | R B_φ = 3 | 30 | 300 |
|---|---|---|---|
| pointwise | 2.35e-2, 6.61e-3, 1.89e-3 | 3.49e-2, 1.05e-2, 2.21e-3 | 1.95e-1, 6.44e-2, 2.05e-2 |
| projected on Q_1 per cell | 2.35e-2, 6.63e-3, 1.89e-3 | 2.99e-2, 8.96e-3, 2.05e-3 | 4.65e-2, 3.02e-2, 6.87e-3 |
| projected on the cell mean (gate) | 2.20e-2, 6.17e-3, 1.83e-3 | 2.28e-2, 7.18e-3, 1.95e-3 | 2.28e-2, 7.21e-3, 1.95e-3 |

With the cell mean the error at R B_φ = 300 equals that at 30 to 0.4 % and is
within 17 % of that at 3; the axis error at 8 cells is 6.8e-6 against 4.7e-3
pointwise. Degree 3, 8 cells: pointwise 2.54e-4, 2.77e-4, 1.34e-3; cell mean
2.40e-4, 2.54e-4, 2.55e-4. Even pointwise, R B_φ = 30 (the range of real
devices) converges at order p here, where the label lane reached order 0.6
(label plan, section 4b).

### Gates F2a, F2b and the QA driver (3D, `milca_flux_map_3d`)

Discretisation: the spline section of F1 times the trigonometric nodal
basis of nt planes over one field period (the complex of
`clebsch_3d_system`), the section of each plane isoparametric in t, vector
components (R, φ, Z), ζ the cylindrical angle. Kernel: the generated leaf
`flux_map_energy_3d` (FortSym proves tangency, rotation invariance and
Hessian symmetry). The wall condition is exact at every point and every t:
Θ_h = Σ_k L_k(t) Σ_b N_b Θ_bk/σ. B enters pointwise. Meshes
(`test/support/flux_map_sections.f90`, since 2026-11-22): every plane is the
pullback of one O-grid of the unit disc by the harmonic map G of its section
onto the disc (FortFEM `fortfem_harmonic_disc`), the mesh wall encloses the
area of the wall curve, and the start is u = G ("Section meshes by harmonic
maps" below). No axis unknown, multiplier, winding form or singular
quadrature. The toroidal flux through every plane is conserved to
round-off in all runs below (defect 1e-14 to 6e-10) and B_φ keeps its sign.

The numbers of this subsection were measured on the former meshes (ruled
O-grid blocks of `layout_set_curve` on the area centroid of each wall curve,
polar start after an affine normalisation by the second moments), unless
marked "harmonic". They are kept as the record; the gates run on the
harmonic meshes.

F2a (`flux_map_3d_gate_f2a`, 5 s): the F1 Solov'ev case through the 3D path
on 3 identical planes reproduces F1: B 2.420e-2, 7.110e-3, 1.932e-3 on 2, 4,
8 cells (F1: 2.42e-2, 7.11e-3, 1.93e-3), 13, 17, 26 Newton steps, axis
5.0e-3, 1.1e-4, 1.3e-5. Harmonic (the gate, start u = G): B 2.394e-2,
7.067e-3, 1.912e-3, 11, 17, 25 steps, axis 4.6e-3, 8.7e-5, 1.1e-5.

F2b (`flux_map_3d_gate_f2b`, 9 s): the exact sheared Landreman family of
E4 and A1-3D (ε 0.6, 13 planes over the torus, degree 2), wall curves as the
only geometric input (mesh centres 6 % of the minor radius off the axis),
polar start at each level:

| cells | unknowns | Newton steps | B error | fixed-label fit of E4a | label lane (A1-3D) | axis | time |
|---|---|---|---|---|---|---|---|
| 2 | 1352 | 11 | 3.49e-2 | 4.98e-2 | 4.80e-2 | 2.0e-3 | 1.9 s |
| 4 | 3536 | 11 | 7.43e-3 | 1.161e-2 | 1.159e-2 | 1.7e-3 | 6.6 s |

Rate of B 2.23. A1-3D takes 21 s from an axis guess within 5 % of the minor
radius, E4 133 s. The axis error falls only from 2.0e-3 to 1.7e-3.

The exact state is not the state of lowest energy, and beyond this
resolution the minimisation leaves it (same case, B error and axis error
against the exact state; energies compare only on the same mesh, because
the spline fit of the wall changes the volume):

| run | B error | axis | W |
|---|---|---|---|
| 4 cells, polar start | 7.43e-3 | 1.7e-3 | 9.310776e-2 |
| 4 cells, from the 2-cell solution | 3.36e-2 | 3.8e-2 | 9.306395e-2 |
| 8 cells, polar start | 1.89e-2 | 3.3e-2 | 9.296213e-2 |
| 8 cells, from the 4-cell solution, minimisation | 2.55e-2 | 2.9e-2 | 9.296084e-2 |
| 8 cells, from the 4-cell solution, Newton on stationarity | 2.48e-3 | 5.6e-3 | 9.298362e-2 |
| degree 3, 2 cells, from the degree-2 solution | 1.60e-2 | 1.6e-3 | 9.331707e-2 |
| degree 3, 4 cells, polar start or nested | 3.2e-2 | 3.3e-2 to 3.7e-2 | 9.3049e-2 |
| 21 planes, 4 cells, polar start | 1.18e-1 | 1.8e-1 | 9.145135e-2 |
| 21 planes, 4 cells, from the 13-plane solution | 7.59e-3 | 5.7e-3 | 9.310760e-2 |

So on the 4-cell and 8-cell meshes there is a second stationary state with
the same fluxes and lower energy (by 4.4e-5 and 2.3e-5), its axis about
3e-2 from the exact one. With `stationary` (Newton with a line search on
the gradient norm, any Hessian signature) the 8-cell solve started from the
4-cell solution converges quadratically to the state next to the exact one
(B 2.48e-3, rate 1.6 from 4 cells) while its axis error grows to 5.6e-3:
the axis displacement is a very soft direction. The label lane met the same
case as "Newton stagnates, then diverges" (label plan, E4, ε = 1.2). The gate
therefore ran 2 and 4 cells from the polar start; options of the test:
`nested`, `stationary`, degree, planes.

Harmonic (the gate since 2026-11-22). From u = G the minimisation goes to the
state of lower energy at 2 cells already (B 4.2e-2, axis 3.1e-2, 4 negative
eigenvalues), and from the polar start at 4 cells too (B 3.6e-2, axis
3.9e-2, none negative). The gate therefore takes the first level from the
polar start (`fm3_start_polar`) and every further level from the previous
one by Newton on stationarity: B 3.658e-2, 7.064e-3, 2.534e-3 on 2, 4, 8
cells (rates 2.37, 1.48; the gate runs 2 and 4), 9, 16, 12 steps, six
negative eigenvalues at every level, axis 1.6e-3, 3.3e-3, 6.2e-3 (not
gated: it lies along the unstable directions).

Hessian spectra (2026-11-22; `MILCA_FM3_HESSIAN=prefix` makes the driver
write the Hessian, the state and a per-plane L² mass matrix of δu at every
level; `benchmarks/flux-map-spike/hessian_spectrum.py` gives the lowest
eigenvalues μ of H v = μ M v over the free slots, with the gauge penalty
included). μ compares between meshes; the volume of the metric converges
(2.835, 2.869, 2.871 on 2, 4, 8 cells). For scale, W/V ≈ 3.2e-2.

| state (13 planes, degree 2) | negative μ | lowest μ |
|---|---|---|
| next to the exact state, 2 cells, polar start | 6 | −1.13e-2 … −3.0e-3 |
| next to the exact state, 4 cells, polar start | 6 | −1.35e-2 … −5.0e-3 |
| next to the exact state, 8 cells, stationary from 4 cells | 6 | −1.23e-2 … −3.9e-3 |
| lower state, 4 cells, from 2 cells | 0 | +7.0e-4 |
| lower state, 8 cells, from 4 cells (nested) / polar | 1 / 1 | −2.5e-4 / −9.3e-5 |

Inertia without the dump (2026-11-22): `flux_map3d_system` assembles only the
upper triangle of the Hessian (local blocks la <= lb, blocks (1,1), (1,2),
(2,2) of the section-pair sums, the (2,1) block being the transpose of (1,2)),
so it is symmetric by construction, and `flux_map3d_solve` factorises that
triangle with FortSparse's LDL^T (`options%triangle`), which reports the
Sylvester inertia. At convergence the solver factorises the unshifted
Hessian once more and sets `flux_map3d_t%negative_eigenvalues`, printed with
`verbose`. It reproduces the dense spectra: 6 on the 4-cell, 13-plane
Landreman state, 0 on the QA vacuum stages (9, 13, 29 planes), 2, 2, 4 on the
Solov'ev states of 2, 4, 8 cells (also saddles of W; the gate checks these).
`flux_map3d_expand_symmetric` mirrors the triangle where a full matrix is
needed (the dump).

- The state next to the exact one is a saddle with six unstable directions
  whose eigenvalues converge in the section mesh. It is not marginal, and
  the polar-start minimiser stops there only because its shift (1e-9)
  makes it Newton, which converges to any stationary point. The earlier
  reading "one soft axis-displacing mode" is withdrawn.
- The lower states are minima up to one small negative μ at 8 cells
  (2 % of the unstable ones), half of it along the rigid rotation of u;
  presumably the θ* gauge against the symmetric penalty rows, not checked.
- The count depends on the toroidal resolution, not on the section mesh:
  8 negative μ on 17 planes at 2 cells (−1.27e-2 … −2.9e-3). On 21 planes
  the polar start goes directly to a state with axis 0.14 from the exact
  one and lower energy (2 cells: W 9.44446e-2 against 9.49015e-2, one
  μ = −6e-5; 4 cells, stationary from it: W 9.2164e-2, gradient 2.6e-5,
  not converged).
- Hypothesis tested and ruled out as stated: the unstable modes are the
  ideal internal kink at ι = 5 (ι runs from 4.73 on the axis to 5.11 at the
  wall) and the lower state is its saturated helical core. Then the axis
  displacement of the lower state would be dominated by |n| = 5. Measured
  (FFT of δR + iδZ over the planes): 13 planes, n = −6 (2.5e-2), +6
  (7.3e-3), −4 (5.6e-3), ±5 below 1e-3; 21 planes, n = −4 (8.7e-2),
  +4, −6. On 13 planes n = 6 is the highest harmonic of the planes.
- Toroidal resolution (4 cells, Newton on stationarity from the 13-plane
  state, `MILCA_FM3_START_PLANES=13`; B error 7.3e-3 to 7.7e-3 throughout,
  so the state stays next to the exact one): negative μ 6, 12, 16, 24 on
  13, 17, 21, 29 planes, lowest μ −1.35e-2, −1.60e-2, −2.54e-2, −3.04e-2.
  Neither count nor lowest μ converges. 4 instead of 2 toroidal quadrature
  points per plane give the same spectra to 4–5 digits (13 and 21 planes):
  not under-integration in t.
- QA (vacuum, below) has no negative μ on 9 to 29 planes, its lowest μ
  converging (2.29e-2 on 2 cells; 1.48e-2, 1.45e-2, 1.43e-2, 1.41e-2 on 13,
  17, 21, 29 planes, 4 cells). A vacuum field has δW = ∫|δB|² ≥ 0, so the
  3D discretisation does not create unstable modes there; the growing
  spectrum belongs to this case (pressure, ι crossing 5 with weak shear).
- Open: a physical pressure-driven instability with an unstable spectrum
  unbounded in n, or a defect of the pressure term. Deciding case: a
  standard configuration with pressure.

QA driver (`test_flux_map_qa`, opt-in, not a registered test): the
Landreman–Paul QA deck of [benchmarks/label-qa](../benchmarks/label-qa/README.md),
one field period, p = 0, ι prescribed, gauge by stellarator symmetry
(symmetry rows until 2026-11-23, since then the moment constraint),
stages "degree:cells:planes" solved in turn from the previous one
(`fm3_prolong`). Compared by `compare.py` with VMEC++ (9 s, energy
0.281794): axis and surfaces on the 64 reference planes, surfaces traced as
preimages of circles. Wall times on 5 threads with other jobs running:

| last stage (start) | unknowns | Newton steps | time | energy | axis max | s = 0.25: distance rms, \|B\| rms | s = 0.50: distance rms, \|B\| rms |
|---|---|---|---|---|---|---|---|
| degree 2, 2 cells, 9 planes (polar) | 936 | 16 | 2 s | +2.06 % | 11.4 mm | 4.8 mm, 2.7 % | 4.2 mm, 2.8 % |
| degree 2, 4 cells, 13 planes (polar) | 3536 | 23 | 22 s | +0.40 % | 1.5 mm | 0.81 mm, 0.58 % | 1.27 mm, 0.75 % |
| degree 2, 4 cells, 13 planes (from 2 cells) | 3536 | 16 + 13 | 13 s | +0.40 % (same W to 11 digits) | | | |
| degree 3, 4 cells, 13 planes (from degree 2, 2 cells) | 5018 | 58 | 68 s | +0.30 % | 1.3 mm | 0.65 mm, 0.52 % | 0.99 mm, 0.65 % |
| degree 2, 8 cells, 17 planes (from 4 cells, 9 planes) | 14416 | 34 | 99 s | +0.06 % | 1.8 mm | 0.77 mm, 0.23 % | 0.43 mm, 0.37 % |
| degree 3, 4 cells, 21 planes (polar; the ladder 2:2:9, 2:4:13, 3:4:13, 3:4:21 gives the same W to 11 digits in 57 steps, 116 s) | 8106 | 27 | 175 s | +0.30 % | 1.5 mm | 0.59 mm, 0.45 % | 0.96 mm, 0.47 % |
| degree 4, 4 cells, 13 planes (polar; the ladder 2:2:9, 4:2:9, 4:2:13, 4:4:13 ends at W lower by 6e-6 relative in 93 steps, 193 s) | 6760 | 43 | 215 s | +0.20 % | 1.5 mm | 0.91 mm, 0.44 % | 0.92 mm, 0.62 % |

Every stage converges (scaled gradient 1e-10 or below), and the polar start
and the nested start reach the same state at degrees 2 and 3 (at degree 4
the energies differ by 6e-6 relative, not followed up): no second state
was met on this vacuum case, where the label lane did not converge. It is not yet competitive with VMEC++:

- Energy: explained by the wall. The volume inside the spline wall
  (the per-plane L² metric of the Hessian dump sums to it) against VMEC's
  0.564712: −1.867, −0.364, −0.055 % on 2, 4, 8 cells (degree 2), −0.276 %
  at degree 3, 4 cells; the energy error is +2.056, +0.395, +0.065, +0.302 %,
  i.e. ΔW/W ≈ −1.1 ΔV/V in every run. On 13 to 29 planes (4 cells) W
  jitters in the 5th digit with the volume (21 planes: ΔV/V −1.0e-4, ΔW/W
  +1.2e-4 against 17). The wall fit, not the solver, sets the energy error;
  it converges with the cells but hardly with the degree.
- Axis and inner surfaces: not explained by the wall. The axis lies
  inward of VMEC++'s at every φ (mean dR −1.29, −1.35, −1.29, −1.02 mm for
  4 cells/13 planes, 8 cells/13 planes, 4 cells/29 planes, degree 3), while
  the volume error falls 6.6× from 4 to 8 cells. Ruled out: a φ offset (fit
  1e-7 rad), the toroidal resolution (13 to 29 planes), the VMEC reference's
  radial resolution (0.06 mm at the axis, ns 101 against 151) and
  VMEC2000 against VMEC++ (1 µm), a uniform ι error (the axis moves 26 mm per
  unit relative change of ι, so 5 % would be needed; the fit error is 1e-4).
- VMEC's state in the flux map's functional (`MILCA_FM3_VMEC`, 
  `milca_vmec_state`, `flux_map_vmec`): W(VMEC) − W(flux map) = 2.3e-3 and
  5.2e-4 on 2 and 4 cells; its scaled gradient by bands s < 0.05, < 0.25,
  < 0.5, < 0.8, wall is 3e-4, 1.6e-3, 3.5e-3, 1.1e-2, 3.7e-2 (4 cells),
  largest at the wall; Newton from it returns the flux-map state (|Δz|
  7e-11). One discrete minimum, VMEC's state relaxes into it. Open: which of
  the two axes is nearer the continuum solution (VMEC's angular truncation
  MPOL 16/NTOR 12 is shared by both VMEC codes); deciding checks: VMEC++ at
  higher MPOL/NTOR, an independent code (GVEC), and j = ∇×B = 0 of both
  states.
- Speed. The Hessian couples all planes of a section cell (dense blocks of
  (2 (p + 1)² nt)²) and is factorised directly: 22 s at 3536 unknowns, 99 s
  at 14416. A toroidal basis with local support or a Krylov solver is
  needed before a time comparison means anything.
- Not run: degree 2, 8 cells, 25 planes from the polar start (stopped
  unfinished), QH and QI decks.

The percent-level \|B\| deviations of the prototype below (4 to 13 % rms
at 10272 unknowns) are not reproduced with splines (0.6 to 0.8 % at 3536
unknowns), so they were not an amplification by the guide field; which of
the prototype's meshes and C0 Lagrange elements caused them was not
separated. A scan of ι → κι in the prototype has its energy minimum at
κ ≈ 0.95 to 0.98, which fixes the orientation (u counter-clockwise in
(R, Z), ι > 0, ζ = φ). Every per-cell projection tried there in 3D (all
reference components of ∇x×∇y or only the plane-normal one, degrees p − 1
and p) is unstable, so the 3D solver has none.

### Integer-transform Landreman saddle

The exact integer-transform equilibrium is a saddle of the fixed-boundary
energy in (2), including for perturbations that preserve stellarator symmetry.
Its field alone is insensitive to a class of perturbations that change the
flux label and pressure surfaces.

Let $D$ be the unit disc, $s=\psi_L/\delta$ the normalized toroidal flux,
and $\varphi$ the geometric toroidal angle. Here $\psi_L$ is Landreman's
positive flux label with boundary value $\delta$. Write the inverse map
$u=(u_1,u_2)$ in the counterclockwise straight-field-line angle, where
the rotational transform is $-2$. Equation (24) of Landreman uses a
clockwise angle and gives $+2$. Rotating the disc by $2\varphi$ gives
the normalized material coordinates $w$ of Landreman's closed field lines.
Here $a=\sqrt{1+\epsilon}$, $b=\sqrt{1-\epsilon}$, and $\epsilon$ controls
the horizontal stretch. The constant toroidal flux derivative gives
$B=-\delta ab\,\nabla w_1\times\nabla w_2$. Equations (14) and (27) give
the volume element $\delta ab\,d^2w\,d\zeta$ and the mean magnetic energy
density $E(w)=(1-\epsilon^2/2)/2+\delta|w|^2$, where $\zeta$ runs once
around a closed field line. [Landreman, sections 2.2–2.5](https://arxiv.org/html/2609.26742v1#S2).

A diffeomorphism $g:D\to D$ equal to the identity at the wall replaces $s$
by $|g(w)|^2$ and multiplies $B$ by $\det Dg$. For
$p(s)=2\delta(1-s)$ and $C=2\pi\delta ab$, its energy is

$$
W[g]=C\int_D \bigl[E(w)(\det Dg)^2-p(|g(w)|^2)\bigr]\,d^2w . \tag{A1}
$$

Let $H(v,h)$ be the mixed derivative of $W[\mathrm{Id}+\alpha v+\beta h]$
at $\alpha=\beta=0$. Take $v$ divergence-free and both increments zero at
the wall. Differentiating the determinant and pressure gives

$$
H(v,h)=C\int_D\bigl[-2E\,\mathrm{tr}(Dv\,Dh)+4\delta\,v\cdot h\bigr]\,d^2w . \tag{A2}
$$

Integrating the first term by parts and using $\nabla E=2\delta w$ gives

$$
H(v,h)=4\delta C\int_D \nabla(w\cdot v)\cdot h\,d^2w . \tag{A3}
$$

Set $q=1-|w|^2$, $v=(\partial_{w_2},-\partial_{w_1})(w_2q^3)$, and
$h=q^2\nabla(w\cdot v)$. Under the stellarator reflection
$(w_1,w_2)\mapsto(w_1,-w_2)$, the first component of each increment is even
and the second is odd.
Equation (A3) and integration by parts give

$$
H(v,v)=0,\qquad H(v,h)=4\delta C\int_D q^2|\nabla(w\cdot v)|^2\,d^2w
=\frac{50\pi^2\delta^2 ab}{63}>0 . \tag{A4}
$$

The determinant of the Hessian on these two directions is
$-H(v,h)^2$. The pressure surfaces change along $v$ while the magnetic
field is unchanged to first order.

For $\epsilon=0.5$, $\delta=1/64$, direct volume quadrature gives
$H(v,h)=1.65615\times10^{-3}$ and $H(h,h)=0.468829$ in the oracle's
energy normalization. At the projected analytic inverse map, the CEMV
Hessian adds positive curvature along $v$. This curvature drops from
$1.03\times10^{-3}$ with four cells and 13 planes at degree two to
$2.41\times10^{-6}$ with 16 cells and 29 planes. The determinant on the
two projected directions changes from $+4.74\times10^{-4}$ to
$-1.62\times10^{-6}$, compared with the exact value
$-2.74\times10^{-6}$. On three interior surfaces the root-mean-square error
in normalized toroidal flux is 0.0273 at the four-cell energy minimizer
and 0.00219 at the analytic projection. A magnetic-field error check alone
can miss departure from the target branch.

### Prototype (Python, C0 Lagrange elements, [benchmarks/flux-map-spike](../benchmarks/flux-map-spike/README.md))

Solov'ev (R0 = 1, own parameters), degree 2, 8 elements per block side,
B_φ/B_pol ≈ 1.6, 16, 160, 1600: pointwise 3.6e-3, 1.4e-2, 8.6e-2, 3.0e-1;
projected on Q_1 per element with the area measure 2.7e-3 at all four, by
continuation in the guide field. Cold starts at B_φ/B_pol ≈ 160 failed with
and without a stabilising term, and continuation failed on the two coarsest
meshes at the highest ratios.

Landreman–Paul QA ([benchmarks/label-qa](../benchmarks/label-qa/README.md)),
3D form of (1) with ζ the cylindrical angle, hexahedral Q_2 elements over
one field period, pointwise B, gauge by stellarator symmetry, cold start from
the mesh's disc coordinates at the prescribed ι:

| mesh | unknowns | Newton steps | W (VMEC++ 0.281794) | axis error max / mean |
|---|---|---|---|---|
| 2 cells, 8 planes | 1296 | 18 | 0.279974 | 5.9 mm / 3.3 mm |
| 4 cells, 16 planes | 10272 | 31 | 0.281640 | 3.8 mm / 3.5 mm |

The gradient reaches round-off in both; the label lane did not converge on
this case. The mean axis error does not fall with refinement, which is read
as the amplification at B_φ/B_pol ≈ 14 without projection (not verified). A
first 3D projection (all three reference components of ∇x×∇y on Q_1 per
element) is unstable: the energy collapses.

### Section meshes by harmonic maps (default since 2026-11-22)

`fm3_build` of `test/support/flux_map_sections.f90` on FortFEM
`fortfem_harmonic_disc` (FortFEM #284, doc/HARMONIC_DISC_MAP.md there).
Check scripts in `benchmarks/flux-map-standard/checks/`. Implementation and
gate numbers changed together; the former numbers are above.

Construction. Per plane, G = (ξ, η) is harmonic on the section with the wall
mapped to the unit circle by normalised arc length. The target is convex, so
G is a diffeomorphism for any simple closed wall curve (Radó–Kneser–Choquet),
nonconvex ones included. G is computed without a mesh (double layer as the
real part of a Cauchy integral, trapezoidal rule on the wall samples,
barycentric evaluation inside). The mesh of a plane is the pullback of one
O-grid of the unit disc (tensor Greville interpolation of G⁻¹ ∘ disc block
map), and the start is u = G. The wall samples start at the outboard wall
point at the height of the area centroid (bisection), a smooth function of
the plane; nothing co-rotates with the section.

Reported per run (the error control of the mesh): the boundary residual of
the harmonic maps (the wall samples are doubled from 256 until it is below
1e-9; 256 sufficed for every deck here, residual 1e-11), the largest
distance of the mesh wall from the wall curve, the error of the section
area, and the certificate that the Bernstein coefficients of every block
Jacobian are positive (a mesh that is not certified is refused). Not
certified: the mesh between the planes (trigonometric interpolation of the
control points; the Jacobian ratio at the quadrature points is printed).

Wall area. With the wall interpolating the curve the energy error is the
error of the section area (QA, 4 and 8 cells: +3.96e-3, +5.9e-4). With the
wall data moved along the normal by one constant per plane so that the area
is exact (`preserve_area`, the default; MILCA_FM3_FREE_AREA switches it
off) it is +5.3e-5 and +2.9e-5, and the axis does not move (4 cells: +0.031
against +0.035 mm; 8 cells: −0.782 against −0.780 mm). The wall distance
stays large on coarse meshes where a tip is not resolved (4 cells: QA
17 mm, li383 33 mm, CTH 2 mm): a harmonic map spreads cells
away from convex tips, and no grading is implemented.

Standard decks, degree 2, stages 2:2:9 then 2:4:13 from u = G, against
VMEC++ (`benchmarks/flux-map-standard`; the former meshes fail to build on
li383):

| deck | Newton steps | energy | axis max dR | s = 0.25 rms, \|B\| rms | s = 0.50 rms, \|B\| rms | negative eigenvalues |
|---|---|---|---|---|---|---|
| QA (vacuum) | 15 + 16 | +5.3e-5 | 0.11 mm | 0.76 mm, 0.18 % | 1.02 mm, 0.39 % | 0 |
| QH (vacuum) | 12 + 10 | +3.6e-4 | 1.66 mm | 1.01 mm, 0.44 % | 1.37 mm, 0.91 % | 0 |
| CTH (pressure, shear) | 11 + 10 | −3.0e-5 | 0.76 mm | 1.2 mm, 0.22 % | 1.4 mm, 0.37 % | 0 |
| li383 (pressure, current) | 13 + 11 | −1.7e-4 | 1.8 mm | 2.4 mm, 0.42 % | 1.8 mm, 0.68 % | 0 |

QA axis (mean dR against VMEC++ ns 101, mm). References: VMEC++ MPOL/NTOR
16/12 → 24/16 → 32/20 changes the axis by 0.004 mm; ns 101 → 151 → 201 by
−0.045, −0.064 (first order, limit about −0.13). DESC 0.17.4 (zero current,
L = M = N = 6, 8, 10): −1.47, −0.139, −0.121. The limit is −0.12 to −0.13.

| mesh, 13 planes | 4 cells | 6 | 8 | 12 |
|---|---|---|---|---|
| former ruled O-grid, core 0.5, degree 2 | −1.36 | −1.67 | −1.34 | |
| former, core 0.35 / 0.65 | −2.52 / −0.46 | / −1.38 | / −1.15 | |
| former, degree 3 | −1.02 | −0.73 | | |
| harmonic, degree 2 | +0.03 | −0.88 | −0.78 | −0.47 |
| harmonic, degree 2, core 0.35 / 0.65 | −0.32 / +0.29 | | | |
| harmonic, degree 3 | −0.36 | −0.22 | −0.18 | |

Harmonic, degree 3, 8 cells, 21 planes: −0.160 (max 0.24), energy +5.7e-6
(13 planes: +2.7e-5). QH (NFP 4, against VMEC++ MPOL 12/NTOR 10, ns 101;
its own radial error 0.08 mm, no third code yet): 13 planes, degree 2, 4, 8,
12 cells +1.05, +1.45, +1.62 and degree 3, 4, 6, 8 cells +1.07, +1.32,
+1.50, with the energy rising to +9e-4; 4 cells on 13, 21, 29 planes +1.05,
+0.72, +0.65 (energy +3.6e-4, −2.1e-4, −2.4e-4); 29 planes: degree 2, 8
cells +0.33 (energy −1.4e-5), degree 3, 6 cells +0.25 (max 0.38).

- The flux map was the inaccurate side on QA, and there is no difference of
  model: with ι scaled by κ = 2 and 3 the axes agreed with VMEC++ to 0.03
  and 0.08 mm on 4 cells of the former meshes, at κ = 0.5 they differed by
  +3.4 mm.
- With harmonic meshes at degree 3 the QA axis converges to the limit of
  VMEC++ and DESC (distance 0.23, 0.10, 0.05 mm at 4, 6, 8 cells on 13
  planes, 0.03 mm on 21). At degree 2 it is not monotone up to 12 cells.
- The floor on QH was the toroidal resolution: 13 planes carry 6 harmonics
  per period. On 29 planes the axis error falls with cells and degree.
- The axis is weakly held: W(κ) of both codes has curvature 1.9e-3 per unit
  κ² (0.7 % of W), VMEC's axis moves 24 mm from κ = 0.5 to 1 and 13 mm from
  1 to 2. The minimum of W_h(κ) lay at κ = 0.993 (VMEC++: 1.0003).
- Enclosed toroidal current (`current.txt` of the driver: the circulation
  of B_h along a surface s_h = const in a plane, mean of 4 planes, over the
  current 2 ι s Φ_edge/R that would carry the transform; zero for a vacuum
  field). QA, degree 2, 13 planes: 0.3 to 0.8 % for s ≤ 0.5 and −1.7 to
  −1.9 % at s = 0.75, 0.9 on 4 cells; 0.1 to 0.6 % and −0.5 % on 8 cells;
  degree 3, 8 cells: up to 0.5 % on 13 planes and below 0.1 % (0.26 % at
  s = 0.9) on 21 planes. QH, 4 cells, 13 planes: up to 2 % inside, 5 % at
  s = 0.9; degree 3, 6 cells, 29 planes: below 0.2 % (0.6 % at s = 0.9).
  With the sensitivity of 26 mm per unit relative ι these sizes match the
  axis errors; the current falls with the axis error and is not a sharper
  predictor of it than the axis itself.
- Not the cause: the sample jitter of the former wall origin (0.06 mm), the
  section quadrature (degree + 8 points: unchanged), the wall area.
- Not measured: the guide-field amplification of section 3 as the mechanism
  behind the remaining error. A 3D projection (roadmap P0b) was not needed
  for the QA/QH axis and is dropped (2026-11-22).

Pressure-driven saddles (Hessian inertia of #653, modes by
`checks/hessian_toroidal_content.py`; pressure scaled by α at fixed ι,
MILCA_FM3_PSCALE).

- Landreman family (F2b, former mesh, 2 cells, 13 planes, polar start): 6
  negative μ at α = 1 (pairs dominated by n = ±6, ±5, ±4; −1.13e-2 …
  −3.0e-3), none at α = 0.5, 0 and −1. W is 0.095 at α = 1 and 0.193 at
  α = 0: the pressure term is as large as the magnetic energy, at ι ≈ 5.
- Solov'ev of F1/F2a (harmonic mesh, 3 planes): W = −237.8 at α = 1 and
  +166 at α = 0, the pressure term is 2.4 times the magnetic energy.
  α = 1: μ = −0.120, −0.117, −0.116 (n = ±1) and −0.042, −0.040, −0.040
  (n = 0 mixed with ±1) on 2, 4, 8 cells, and on 8 cells also −1.4e-2
  (n = ±1) and −3.9e-3 (n = 0 mixed), whose counterparts on 4 cells are
  +1.7e-2 and above. α = 0.5, 4 cells: −1.5e-2 (n = ±1), −4.7e-3. α = 0,
  4 cells: none (lowest +3.9e-3); 2 cells: 2 negative; 8 cells and α = −1:
  Newton did not converge from u = G. The gauge rows act in plane 1 only
  and mix n = 0 with ±1, so the n = 0 label is not clean.
- CTH (p₀ 432 Pa, strong shear) and li383 at α = 1: no
  negative eigenvalue on 2:2:9 and 2:4:13.
- Read: the negative eigenvalues of both analytic cases are driven by the
  pressure, in cases whose pressure energy equals or exceeds the magnetic
  energy; the standard decks with pressure have none. The Mercier criterion
  of the two analytic states has not been evaluated.

### Force-balance residual, independent of the energy (2026-11-23)

Module `milca_flux_map_force`, leaf `flux_map_force_3d` (FortSym, in
`gen/gen_flux_map.f90`; proved there: ∇·B = 0 and F·B = 0), test
`flux_map_force_balance`. Nothing of it enters the solver.

Definition. From the 2-jet of u = (x, y) in the cylindrical coordinates
(R, φ, Z) the leaf forms B of (1), J = ∇×B, F = J×B − p′(s)∇s and
∇(B²/2), all pointwise; J comes from the second derivatives of the map and
not from the gradient of W. For the discrete state the jet is that of u_h
inside each cell (the section functions with their second derivatives, the
trigonometric basis in t, the isoparametric planes and the wall part
σ (cos Θ_h, sin Θ_h), changed to (R, φ, Z) by the chain rule). Reported:

- η_B = ‖F‖/‖∇(B²/2)‖ and η_p = ‖F‖/‖∇p‖, ratios of L² norms over the
  volume (η_p only with pressure; a vacuum state has η_B = ‖J×B‖/‖∇B²/2‖),
  the same over s > 0.1, and per band of s (ten bands) the root mean
  squares;
- the parts of F along ∇s (radial balance) and along B×∇s (equal to
  |B| |J·∇s|/|∇s|: current crossing the surfaces); F·B = 0 identically;
- separately, because the volume norm does not see it: the jump of B²/2
  across the C0 block interfaces of the section (B·n is continuous there,
  the tangential field need not be), relative rms.

The same measure for a VMEC state (`fm3_vmec_force`, `vmec_state_jet`): no
mesh; u = ρ (cos θ*, sin θ*), θ* = θ + λ, as a function of (R, φ, Z)
through the 2-jets of R, Z, λ in (ρ, θ, ζ), on the midpoints of 200 equal
intervals of s. Radially a coefficient is ρ^parity h(s), parity = mod(m, 2),
h the C2 cubic spline through the grid values (VMEC's own variables); the
profiles are those of the deck. The first version divided by ρ^m as
`vmec_state_eval` does and gave residuals of order 100 inside s = 0.2 on
QA: the high modes of a wout file are not ∝ ρ^m near the axis.

Oracles (`test_flux_map_force`): the leaf on B = e_φ/R (J = 0) and on a
poloidal field with J_φ = −ψ′ Δ*s/R; the change of variables of the 2-jet
on an analytic function; J and F of a discrete 3D state (Landreman
sections, 13 planes, degree 3, 15 points with wall cells) against the curl
of B_h of `flux_map3d_point` by central differences (7e-10 relative); and
the VMEC route on an exact state, the vacuum field B₀R₀/R in a torus of
circular section, whose flux coordinates are closed forms (2πB₀R₀ (R₀ −
√(R₀² − r²)) = sΦ, θ* = 2 atan(√((R₀ − r)/(R₀ + r)) tan(θ/2))): η_B =
1.8e-9 and 1.5e-10 on 26 and 51 surfaces, the axis region included.

Convergence on the exact cases (relative L² error of B; η_p). Solov'ev of
F2a, 3 planes, start u = G:

| degree | cells 2 | 4 | 8 | 16 | 32 | rates |
|---|---|---|---|---|---|---|
| 2, B | 2.394e-2 | 7.067e-3 | 1.912e-3 | 4.879e-4 | 1.233e-4 | 1.76, 1.89, 1.97, 1.98 |
| 2, η_p | 1.03e-1 | 6.43e-2 | 3.73e-2 | 1.96e-2 | 9.94e-3 | 0.68, 0.79, 0.93, 0.98 |
| 3, B | 7.098e-3 | 1.801e-3 | 2.364e-4 | | | 1.98, 2.93 |
| 3, η_p | 3.68e-2 | 1.87e-2 | 4.58e-3 | | | 0.98, 2.03 |
| 4, B | 2.338e-3 | 5.424e-4 | 3.084e-5 | | | 2.11, 4.14 |
| 4, η_p | 1.17e-2 | 4.16e-3 | 6.15e-4 | | | 1.49, 2.76 |

B converges at order p and the residual at order p − 1 (J takes one more
derivative), from 4 to 8 cells on; 2 cells are below the asymptotic range
at degrees 3 and 4. Jump of B²/2 at the block interfaces, degree 2: 5.9e-3,
2.1e-3 on 2, 4 cells.

The order 1.8–1.9 of F1 is not an order loss (prediction: if it is
pre-asymptotic, the rate from 8 to 16 cells is at least 1.95; confirmed).
On the F1 mesh itself (`test_flux_map_axisym`, pointwise, cells 8, 16, 32):
B 1.93e-3, 4.90e-4, 1.23e-4 (rates 1.98, 1.99), B_pol 1.89e-3, 4.89e-4,
1.23e-4, the projected exact label 9.08e-4, 2.13e-4, 5.23e-5 (2.09, 2.03).
B_pol over the projected label is 1.04, 1.49, 2.08, 2.30, 2.35 on 2 to 32
cells: it levels off, a constant of 2.3 to 2.4 against the L² projection
of the exact label. The projection converges faster than its order on the
coarse meshes (2.34, 2.29), which made the ratio grow there.

Guide field (isolating run, Solov'ev, degree 2, `MILCA_FM3_GUIDE` scales
the toroidal field; the state stays exact). R B_φ = 3, 9, 30:

| R B_φ | B error, 4 cells | η_B on 2, 4, 8 cells | η_p |
|---|---|---|---|
| 3 | 7.07e-3 | 0.19, 0.12, 0.070 | 0.103, 0.064, 0.037 |
| 9 | 7.56e-3 | 0.31, 0.21, 0.11 | 0.180, 0.116, 0.061 |
| 30 | 6.89e-3 | 0.62, 0.45, 0.27 | 1.05, 0.66, 0.37 |

At the same error of B the residual grows with the guide field (3.8 times
in η_B and 10 times in η_p from R B_φ = 3 to 30 on 4 cells; prediction: at
least 4 times in η_B), the part along B×∇s most (0.034 → 0.37 of
‖∇B²/2‖). The order stays p − 1. A first-order current times a strong
toroidal field: at the field ratio of devices a degree-2 state has no
usable pointwise force balance on these meshes.

Sheared Landreman family of F2b (13 planes, stationary mode, first level
from the polar start, degrees 3 and 4 from the degree-2 state on 2 cells):

| degree | B on 2, 4, 8 cells | η_p | negative eigenvalues |
|---|---|---|---|
| 2 | 3.658e-2, 7.064e-3, 2.534e-3 | 0.221, 0.106, 0.071 | 6, 6, 6 |
| 3 | 1.534e-2, 6.017e-3, 2.008e-3 | 0.217, 0.066, 0.081 | 6, 7, 6 |
| 4 | 7.263e-3, 3.381e-3, 2.351e-3 | 0.111, 0.089, 0.103 | 6, 6, 6 |

No order beyond 4 cells: B stops at 2.0e-3 to 2.5e-3 and η_p at 0.07 to
0.10 for every degree, the part of F along B×∇s still falling (0.085,
0.042 at degree 2; 0.0086 at degree 4, 8 cells) and the radial part not.
A floor that does not depend on the section resolution; hypothesis, not
isolated: the 13 planes, or the position along the six unstable directions
(axis error 5e-3 to 7e-3 at 8 cells for all degrees). The runs meant to
isolate it left the state: from the 13-plane state on 2 cells, 17 and 21
planes give B 2.3e-2 and 1.0e-2 at 4 cells (axis 2.6e-2, 2.8e-2); from the
polar start on 4 cells the solve goes to the lower state as before.

Standard decks, degree 2, stages 2:2:9 then 2:4:13 (the table of "Section
meshes by harmonic maps"; axis and energy reproduce it), against the
VMEC++ states (ns 101) on the same measure, over s > 0.1:

| deck | CEMV η_B | radial, B×∇s parts | VMEC η_B | CEMV η_p | VMEC η_p | interface jump |
|---|---|---|---|---|---|---|
| QA | 0.73 | 0.57, 0.43 | 1.8e-3 | | | 8.9e-3 |
| QH | 0.82 | 0.71, 0.38 | 1.7e-3 | | | 1.5e-2 |
| CTH | 0.48 | 0.36, 0.30 | 4.3e-3 | 31 | 0.24 | 7.7e-3 |
| li383 | 0.66 | 0.57, 0.33 | 3.5e-2 | 3.0 | 0.11 | 1.8e-2 |

(The parts are of the whole volume.) VMEC++'s own residuals of these
states: FSQR, FSQZ, FSQL at 1e-14 to 1e-15 (the discrete equations are
solved); `equif`, the surface-averaged radial balance, rms over s > 0.1:
CTH 4.2e-4, li383 1.1e-4 (undefined for the vacuum decks).
Inside s = 0.1 the VMEC states have η_B of order one, confined to the six
grid surfaces next to the axis (QA: s < 0.06 at ns 101, s < 0.04 at
ns 151; the whole volume gives 0.33 and 0.23); the exact vacuum state
through the same route has none, so this belongs to the states' radial
grid data as differentiated twice, not to the evaluation.

- On the pointwise measure the degree-2 flux-map states are 20 to 500
  times worse than VMEC++ at ns 101, although axis, surfaces and |B| agree
  to millimetres and below one per cent: the error of B_h is small but
  varies on the cell scale, J_h is of order p − 1, and the strong toroidal
  field amplifies it (the Solov'ev scan above).
- Prediction ruled out: at least 3 times less at degree 3 on 4 cells. QA,
  13 planes, η_B over s > 0.1: 0.84 (2:2:9), 0.73 (2:4:13), 0.78 (2:8:13),
  0.64 (3:4:13), 0.72 (4:4:13), 0.81 (3:8:13).
- Toroidal resolution is one part (prediction: at least 2 times less from
  13 to 21 planes at degree 3; measured 1.4): 4 cells on 13, 21, 29 planes
  0.73, 0.67, 0.67 at degree 2 and 0.64, 0.46, 0.43 at degree 3.
- On 21 planes: 0.46 (3:4), 0.35 (3:8), 0.43 (4:4), 0.39 (4:6). Per band
  the interior converges and the wall layer does not: from 3:4:29 to
  3:8:21 the bands of s < 0.5 fall from 0.07–0.18 to 0.015–0.07 (order 2
  and above), the band s > 0.9 stays at 0.89 and 0.93. The volume number
  is that layer. Cause not isolated (hypotheses: the unresolved tips of the
  wall, the toroidal content of the wall; the graded meshes below lower
  the band to 0.62).
- Parity with VMEC on this measure is therefore open, and it is the
  measure of the standing benchmark of the roadmap.

Sign of the map and winding of the wall, at convergence (reported per
stage by both drivers; not a step-acceptance rule). ∇φ·(∇x×∇y) has one
sign at all quadrature points (section points times toroidal points) of
every run above: Solov'ev −0.77 to −0.22; Landreman −8.86 to −4.13; QA
−60.3 to −24.8; QH −68.7 to −28.1; CTH −99.4 to −20.3; li383 −10.6 to
−4.0 (2:4:13). The wall angle winds once on every
plane (1.000000) and is monotone: the smallest step between neighbouring
wall points is 0.39 to 0.58 of the mean step on the decks. The gates F2a
and F2b check one sign and winding one.

### Mesh grading towards tips (2026-11-23)

FortFEM `harmonic_disc_build(..., grading=c)` (doc/HARMONIC_DISC_MAP.md
there): the wall is mapped to the circle by the measure
(1 + (c κL/2π)²)^(1/3) dl instead of the arc length l, κ the curvature and
L the length of the wall curve. More wall cells where the wall bends; the
correspondence is still a homeomorphism onto the circle, so the map is
still a diffeomorphism and the certificate is unchanged. Opt-in here
(`MILCA_FM3_GRADING=c`, argument `grading` of `fm3_build`); the default
stays the arc length.

On QA, 3:4:21: wall distance 15 → 4.4 mm (c = 0.5), 2.9 mm (c = 1); η_B
over s > 0.1 0.46 → 0.26, 0.39; band s > 0.9 1.02 → 0.62, 1.05; axis max
dR 0.63 → 0.45, 0.34 mm; energy +5.2e-5 → −1.4e-5, −8.9e-6. Here the
grading helps, c = 0.5 most in the residual; prediction (band s > 0.9
at least halved) not met.

### Fourier wall sample angles

For a boundary supplied as a Fourier series, uniformly spaced poloidal
samples give a periodic boundary correspondence directly. The benchmark
driver maps them to uniform angles on the unit disc, starting at the outboard
anchor of each plane. `MILCA_FM3_SAMPLE_ANGLE=0` restores arc-length
placement and `MILCA_FM3_GRADING=c` with c > 0 curvature grading. Both
reparameterize the wall, which can introduce higher harmonics before it is
interpolated on the mesh.

Arc length against sample angles on the four standard decks, degree four,
stages `2:2:9,3:4:21,4:6:21`, one sCluster core each. The distance is the
RMS distance of the s = 0.5 surface from the VMEC++ reference, and η_B is
the force ratio over the whole volume.

| Deck | Distance, arc length [mm] | Distance, sample angle [mm] | η_B, arc length | η_B, sample angle | Wall time, sample angle [s] |
|---|---|---|---|---|---|
| QA | 0.28 | 0.055 | 0.51 | 3.6e-3 | 40 |
| QH | 0.73 | 0.16 | 0.73 | 3.2e-2 | 55 |
| CTH | 0.41 | 0.17 | 6.4e-2 | 3.3e-2 | 45 |
| li383 | 1.2 | 0.39 | 0.36 | 0.11 | 60 |

On the quasi-axisymmetric vacuum deck with degree four, six cells and
21 planes in the final stage, changing from curvature grading with $c=1$
to the input sample angles reduced the maximum wall distance from
710 to 0.86 µm and the volume-integrated force ratio fell from 0.145 to
0.00351. The half-flux surface distance from stock VMEC fell from 342 to
100 µm RMS. This run used the prescribed wall and profiles in both fits.

### Gauge constraint, one solver loop and continuation (P0c, 2026-11-23)

What changed. `flux_map_solve` and `flux_map3d_solve` were two copies of one
loop. Both now call `flux_map_newton`, which is FortFEM's
`fortfem_constrained_newton` (FortFEM doc/CONSTRAINED_NEWTON.md): Newton on
the KKT system of W with linear constraints, factorised by LDLᵀ, with the
inertia of the Hessian on the constraints. Both Hessians are assembled as
upper triangles; the axisymmetric solve gained the inertia and the
stationary mode. The penalty rows, their weight and the symmetry rows of the
tests are deleted. Implementation and gate conditions changed together: F1
and F2b moved to meshes without symmetry, and the F1 guide-field part starts
cold at the strongest field. The earlier numbers stay above as the record.

Gauge (S1). The rotation of the disc by one angle maps the discrete space
onto itself (interior values rotate, wall angles shift), so it is an exact
symmetry of W_h and the Hessian is singular along it at every stationary
point. A rotation by k(s) with k not constant is a symmetry of W but not of
W_h: the rotated u_h is not in the spline space. The rotation is removed by
one linear constraint from `milca_gauge_moments` (moment constraints
orthogonal to a continuous symmetry: rows ⟨g_j(z_ref), z − z_ref⟩_M = 0 for
generators g_j, generators that are dependent in the Gram matrix dropped;
`test_gauge_moments` checks it on the six rigid motions of a point cloud).
Here g is the generator (−y, x), whose pointwise moment with a variation is
the generated leaf `flux_map_gauge_moment` (the generator of the invariance
proofs of `gen_flux_map`), M is restricted to the wall functions (of plane 1
in 3D), and the reference follows the iterate: every Newton step keeps the
mean wall angle of that plane, weighted by the moment of the generator with
each wall function. The multiplier vanishes at convergence (1e-16 to 1e-18
in the runs below).

| Hypothesis and prediction | Result |
|---|---|
| The KKT matrix with the gauge row is factorised without shift or penalty; prediction: gates unchanged, inertia reliable, linear solves at most 1.5 times as long | confirmed for the row on the wall functions of one plane (0.56 s against 0.48 s without a row, F2b at 3536 unknowns); ruled out for a row over the whole volume (4.14 s) and over the section of one plane (0.89 s) |
| On symmetric meshes the constraint reproduces the numbers of the symmetry rows | confirmed: F1 B 2.42e-2, 7.11e-3, 1.93e-3 and 14, 17, 23 steps (rows: 14, 18, 23); F2a identical with inertia 2, 2, 4; F2b on the symmetric planes 3.658e-2, 7.064e-3, 9 + 16 steps, 6 negative eigenvalues; the decks below |
| The same accuracy on meshes without symmetry, without symmetry rows (the gate) | confirmed: F1 with the centre 10 % above the midplane B 2.42e-2, 7.11e-3, 1.93e-3, axis 4.8e-3, 1.2e-4, 1.3e-5, rates 1.77, 1.88; F2b on planes shifted by 0.37 of their spacing B 3.658e-2, 7.066e-3 (rate 2.37), 10 + 19 steps, 6 negative eigenvalues, axis 1.6e-3, 3.4e-3 |
| Constraining several k_j(s) is harmless below the radial resolution and removes the slow steps | ruled out as stated (table below) |

In 2D the Hessian on the constraint has no negative eigenvalue on any of the
three meshes; the 2, 2, 4 of F2a belong to the three planes (n = ±1).

Moments of several gauge functions k_j(s) (F1, pointwise B_φ, R B_φ = 3,
mesh without symmetry, energy acceptance; measured during development with
moments over the whole section, code that is not landed). "Start": the rows
are fixed at the start state; "iterate": they follow the state.

| gauge functions | rows at | B error on 2, 4, 8 cells | Newton steps |
|---|---|---|---|
| none (the shift holds the rotation) | | 2.42e-2, 7.11e-3, 1.93e-3 | 14, 22, 33 |
| rotation, moment over the section | iterate | 2.42e-2, 7.11e-3, 1.93e-3 | 15, 28, 45 |
| rotation, wall moment (landed) | iterate | 2.42e-2, 7.11e-3, 1.93e-3 | 15, 27, 68 |
| cells + 1 hat functions of √s | start | 2.46e-2, 7.43e-3, 2.65e-3 | 15, 21, 29 |
| cells + 1 hat functions of √s | iterate | 2.44e-2, 7.42e-3, 2.06e-3 | 12, 21, 28 |
| 2 Legendre polynomials in s | iterate | 2.42e-2, 7.11e-3, 1.93e-3 | 13, 21, 36 |
| 3 | iterate | 2.43e-2, 7.14e-3, 1.94e-3 | 14, 22, 33 |
| 4 | iterate | 2.47e-2, 7.16e-3, 1.94e-3 | not converged, 20, 30 |
| 6 | iterate | 3.21e-2, 8.17e-3, 2.01e-3 | 17, not converged, 25 |
| 4, then the rotation alone from that state (relaxed steps) | iterate | 2.42e-2, 7.11e-3, 1.93e-3 | 13 + 4, 20 + 4, 30 + 4 |
| 6, then the rotation alone (relaxed steps) | iterate | 2.42e-2, 7.11e-3, 1.93e-3 | 17 + 7, 22 + 8, 25 + 7 |

The last two rows were run with relaxed steps (below), with which the
rotation alone over the section takes 13, 28, 44 steps. On the symmetric
mesh the rows with 0, 1, cells + 1 and 16 functions were run and give the
same numbers (14, 17, 23 steps). Read:

- The k(s) directions are soft, not null. A mesh without symmetry excites
  them and the minimiser needs up to three times the steps at 8 cells (68
  against 23); the count moves by a third with the choice of the single row
  and is lowest without one (33).
- A constraint on them is a constraint on the solution: 0.5 % in B with
  three smooth functions, 4 to 33 % with six, and the axis by more (3.0e-3
  against 1.2e-4 at 4 cells with hat functions fixed at the start).
- Constraining and then releasing them returns the unconstrained state
  and saves steps at 8 cells (34 and 32 against 44). Not landed: rows over
  the section cost 1.9 to 8.6 times in the 3D factorisation. Open.

Step acceptance (S2). One loop, two modes as before: minimisation (energy
test with an identity shift) and the stationary point next to the start
(line search on the scaled gradient norm).

| Hypothesis and prediction | Result |
|---|---|
| The natural monotonicity test (the simplified Newton correction from the factors of the step must contract; affine covariant, so blind to the stiffness) converges from the cold start at R B_φ = 300 | ruled out: as the only test it finds no admissible damping after 1 to 26 steps from the polar start already at R B_φ = 3, on both meshes, and after 3, 6, 10 and 20 energy steps at R B_φ = 300; as a fallback of the energy test it accepted none of the refused steps |
| A second-order correction (the simplified correction added to the refused trial point) rescues the refused steps | ruled out: none accepted at R B_φ = 300 |
| The steps refused at a strong guide field are refused for the curvature of the constraint R τ det ∇u = f(s) | ruled out: the matrix is indefinite there. Negative eigenvalues on the gauge constraint along the cold path at R B_φ = 300, 4 cells, at the first trial shift of steps 2, 6, 10, 14, 30: 72, 17, 5, 1, 51; the accepted shift stays at 1e-3 of the smallest diagonal, which grows with B_φ², for 25 steps. At R B_φ = 3 (2 cells) it falls to 1e-5 within 14 steps |
| The gradient-norm line search of the stationary mode can be replaced by the natural monotonicity test | ruled out: F2b stops after 3 steps on the second level |
| Relaxed steps (watchdog: a refused step on a positive matrix is taken on trial and judged with the next one; FortFEM `NEWTON_ACCEPT_ENERGY_RELAXED`) reduce the steps | confirmed on a curved stiff valley (FortFEM test: 568 → 109 steps at stiffness 1e6) and on F1 without symmetry (pointwise, 8 cells: 68 → 38); ruled out on the decks (QA second stage 16 → 23 to 27 steps, QA wall time +20 %). Off by default (`relaxed_steps`) |

So no acceptance rule tested makes the cold start at a strong guide field
direct; the obstacle is the indefinite Hessian far from the solution
(section 3: the spurious modes freed by the projection), which only a state
near the solution avoids. That is what the continuation supplies.

Continuation (S2). `fortfem_parameter_continuation`: the caller solves one
stage P(s), the driver chooses s: P(0) from the start, then steps that
double after a converged stage and halve after a failed one. F1 uses it
with R B_φ = 3 · 100^s from the polar start, at most 150 Newton steps per
stage, no stage set by hand (B_φ projected on the cell mean; B_pol error at
the three guide fields, R B_φ = 30 solved from the state at 300):

| mesh | cells | stages (failed) | Newton steps, s = 0 + step to 1 | B_pol at R B_φ = 3, 30, 300 |
|---|---|---|---|---|
| without symmetry (gate) | 2 | 2 (0) | 16 + 20 | 2.20e-2, 2.28e-2, 2.28e-2 |
| | 4 | 2 (0) | 24 + 77 | 6.18e-3, 7.17e-3, 7.21e-3 |
| | 8 | 2 (0) | 70 + 58 | 1.83e-3, 1.95e-3, 1.95e-3 |
| symmetric | 2 | 2 (0) | 17 + 18 | 2.20e-2, 2.28e-2, 2.28e-2 |
| | 4 | 2 (0) | 21 + 25 | 6.14e-3, 7.20e-3, 7.24e-3 |
| | 8 | 2 (0) | 37 + 10 | 1.83e-3, 1.95e-3, 1.95e-3 |

In B_pol the states are those of the hand-set stages 3 → 30 → 300 (2.28e-2,
7.21e-3, 1.95e-3; 17 + 6 + 4, 23 + 5 + 4 and 41 + 4 + 3 steps). The error of
B with the pointwise B_φ differs between paths and meshes at 4 cells (3.05e-2
by the hand-set stages, 2.26e-2 and 2.83e-2 in the two rows above) and at
R B_φ = 3 (1.03e-2 and 1.14e-2): with the cell mean in the energy the
converged states differ in the pointwise B_φ. Not followed up. The step from 3 to 300
converges in one stage but is not cheap; a stage limit of a few steps would
make the driver insert R B_φ = 30 as the hand-set path did. Without the
continuation the cold start at R B_φ = 300 converges in 44 steps on 2 cells
and is not converged after 60 on 4 and 8 cells (B_pol 7.8e-3 and 1.6e-2
there). Rates of B_pol at the three guide fields on the gate mesh: 1.83,
1.67, 1.66 and 1.75, 1.88, 1.88. The gate takes 33 s (10 s before).

NN2. `nc_continue` (anisotropy ratio) runs on the same driver. P2, N = 2,
nr = 4 and 8 at ratio 1e6 from the shifted seed: 8 stages and 16 Newton
steps on each mesh, last backward error 8.6e-16 and 1.0e-15, path check
1.08e-9, 9.79e-8, 4.76e-9 in A, p, V: every printed digit as before the
change. The NN2 Newton iteration itself (a non-symmetric residual with its
own line search) is not on the shared loop.

Standard decks after the change (stages 2:2:9 then 2:4:13, 8 threads; the
table of "Section meshes by harmonic maps"):

| deck | Newton steps | energy | axis max dR | time (before) |
|---|---|---|---|---|
| QA | 15 + 16 | +5.33e-5 | 0.109 mm | 7.9 s (7.4) |
| QH | 12 + 10 | +3.64e-4 | 1.66 mm | 6.5 s (5.9) |
| CTH | 11 + 10 | −3.0e-5 | 0.76 mm | 6.3 s (5.8) |
| li383 | 13 + 11 | −1.67e-4 | 1.77 mm | 6.6 s (6.2) |

Steps, energies (QA to 11 digits) and axes are unchanged, no negative
eigenvalue; the gauge row and its rebuild at every step cost 3 to 10 % of
the time.

## 5. Sequence

The order of work, the gates still owed for F2 and the comparison with the
Lagrangian map are in the [roadmap](equilibrium-roadmap.md), "Way forward"
(P0–P6). F1 and F2 are done as far as section 4 records.

## 5a. Separatrices, X-lines and the model per region

Decision of 2026-11-22 (argued here; nothing of it is implemented or
measured in the flux-map form).

The singularity. Any field-line label is singular at a separatrix: a field
line needs infinitely many toroidal turns to pass the X-line, so q → ∞
like log|ψ − ψ_X|. For the flux map this means

- θ* is constant along each separatrix branch away from the X-line and its
  whole variation is squeezed into the neighbourhood of the X-line,
- det ∇u → 0 like 1/log along the whole separatrix, while B (a diverging
  profile times a vanishing determinant) and A stay smooth.

It does not depend on the direction of the map: the inverse map of VMEC has
the same singularity. Measured in the (χ, g) form
([plan](eulerian-label-plan.md), E2a, E2b, E5u): B at order 1/2 near a
separatrix on an unaligned mesh, order p on an aligned mesh graded towards
it with fixed labels; Newton on that mesh only at a weak guide field.

Why Grad–Shafranov has no such problem. Symmetry supplies (i) a surface
label that is a component of A (ψ = R A_φ, smooth at a saddle) and (ii) the
toroidal field as f(ψ)/R, without a field-line label. Without a continuous
symmetry neither exists, and B·∇ has to be inverted somewhere:

| Route | Where B·∇ is inverted | Cost |
|---|---|---|
| Ampère iterated with J from force balance (the Picard loop of GS in 3D; PIES, HINT) | magnetic differential equation for J_∥ | 1/(m ι − n) currents; not developed here |
| A ∈ V1 with a label and B·∇χ = 0 as a constraint ([eulerian-3d.md](eulerian-3d.md)) | the constraint | locks with an equal-degree multiplier; not shown dead ([audit](equilibrium-lane-audit.md)) |
| Clebsch field, flux map | the representation | labels singular at separatrices |

Region rule that follows.

- Flux map: volumes of closed surfaces bounded by good surfaces. The core
  map ends at a chosen surface s = s_i inside the separatrix, not on it;
  an island interior likewise ends inside its own separatrix.
- A-based model (NN2, Beltrami or vacuum): the separatrix layer, X-lines,
  open field lines (scrape-off layer, private flux region) and chaotic
  volumes, where B and A are smooth and labels are not. In 3D that layer is
  chaotic anyway, so nested ideal surfaces up to the separatrix are not the
  physics there.
- Interface at s = s_i: B·n = 0 from both sides, continuity of the
  tangential trace of A (fluxes), balance of p + B²/2, and a transport-flux
  law where NN2 meets the map. Its position is an unknown; the mesh follows
  it.
- In axisymmetry the combination reduces to Grad–Shafranov. An
  axisymmetric diverted tokamak is solved in the Grad–Shafranov lane (gate
  A2a has two axes and an X-point); the flux map adds something there only
  when q or ι is prescribed instead of f.
- Kept as alternatives, not planned first: the aligned, graded mesh up to
  the separatrix (E2b, fixed labels only); the analytic form of the map at
  an X-line built into the discretisation (Reiman, Pomphrey and Boozer 1989,
  cited in the label plan; untested).

Why not A everywhere. That is the NN2 lane, a different physical model and
not the same model in other variables:

- B·∇p = 0 needs a closure: a constraint (locks), flat pressure per relaxed
  region (the SPEC model) or strongly anisotropic transport (NN2);
- profiles are outcomes of sources and transport, not inputs, so p(s) and
  ι(s) of a VMEC deck cannot be posed one to one;
- topology is not preserved: islands open where the model allows;
- anisotropic transport on an unaligned mesh pollutes the perpendicular
  balance, so alignment returns in another form (NN2 on an aligned mesh, roadmap P6);
- no flux coordinates for kinetic and neoclassical consumers and no ideal δW
  from the solve.

The flux map gives exact surfaces, prescribed profiles, the ideal limit and
flux coordinates where surfaces are good; NN2 is strong exactly where the
map is singular. Which lane leads depends on the application: nested
stellarator equilibria at VMEC quality need the map in 3D (F2); islands,
X-lines and edges need NN2 with the map as its ideal core.

Open: the coupling does not exist; the interface conditions above are to be
derived from the variation of the shared interface (label plan, section 7);
the choice of s_i and its effect on the solution; the state of NN2 against
these needs has not been reviewed.

## 6. Open

Each with the run that would settle it in the roadmap (P0a–P0c).

- Equilibria that are not minima (F2b): the default solver minimises W and
  leaves a state with a soft or unstable direction; `stationary` finds it
  from a good start. Owed: the lowest eigenvalues of the Hessian at both
  states, and a rule for which solver a run uses. Distinct stationary
  states are distinct equilibria.
- QA: the floor of the axis and inner surfaces, the missing gain from the
  degree, the cost of the plane-dense Hessian (section 4).
- Pointwise force balance: the degree-2 states have J×B − ∇p of the order
  of ∇(B²/2) on the standard decks (VMEC++: 2e-3 to 3.5e-2), and on QA the
  wall layer s > 0.9 does not converge up to degree 3, 8 cells, 21 planes
  ("Force-balance residual"). Owed: the cause in the wall layer, then
  parity at degree 3 or 4. The floor of F2b beyond 4 cells (planes or
  the saddle) is not isolated.
- Gauge: done for the rotation, the one exact symmetry of W_h, on meshes
  without symmetry ("Gauge constraint, one solver loop and continuation").
  Open: the soft k(s) directions cost up to three times the Newton steps on
  a mesh without symmetry; constraining and releasing them saves a quarter
  of the steps in the one case measured and needs rows that keep the 3D
  factorisation sparse.
- Cold starts at a strong guide field: done by the continuation driver. No
  acceptance rule tested makes the direct start work; the Hessian is
  indefinite far from the solution there.
- A state without symmetry (the meshes of F1 and F2b have none now, the
  equilibria still do).
- Projection in 3D: none is stable so far; dropped (roadmap P0b,
  2026-11-22). A candidate is B_h = Π_2(u_h*β) with a commuting
  projection, in mixed form.
- Section meshes: the curvature grading brings the wall 3 to 5 times
  closer but does not improve the equilibria and needs more planes ("Mesh
  grading towards tips"; off by default); a correspondence that does not
  slide with φ is not tried. A certificate of the mesh between the planes,
  the origin of the wall samples for sections without a unique outboard
  point at centroid height.
- QH against a third code; the QA and QH axis at degree 3 on more planes.
- The relation of the Hessian to the ideal δW (section 1) is argued, not
  verified.
- Not run: QI decks, a non-symmetric case, free boundary, any region
  coupling.
- Search for prior work on the direct map before any publication.

## Sources

- Kruskal and Kulsrud, Phys. Fluids 1, 265 (1958): the Clebsch form with a
  single-valued part and multivalued angles and the variational principle
  (read for the label plan, section 3).
- Elguedj, Bazilevs, Calo and Hughes, B-bar and F-bar projection methods for
  nearly incompressible elasticity with higher-order NURBS elements,
  Comput. Methods Appl. Mech. Engrg. (2008): projection of the volumetric
  part onto a lower-degree spline space. Cited from memory, to be confirmed
  (volume, pages, DOI) and read before it is relied on.
- Babin, Hindenlang, Maj and Köberl, Construction of an invertible mapping
  to boundary conforming coordinates for arbitrarily shaped toroidal
  domains, arXiv:2411.04683: section meshes as inverted harmonic maps onto
  the disc (the idea; implementation and its further sources in FortFEM,
  doc/HARMONIC_DISC_MAP.md, commit 5e881d50).
- The exact vacuum state of `test_flux_map_force` (field B₀R₀/R in a torus
  of circular section, its flux label and straight-field-line angle) and
  curl(∇s×∇φ) = −Δ*s ∇φ are derived here, not taken from a source.
- Step acceptance: Deuflhard, Newton Methods for Nonlinear Problems
  (Springer 2004), for the natural monotonicity test; Chamberlain, Powell,
  Lemaréchal and Pedersen, The watchdog technique for forcing convergence in
  algorithms for constrained optimization, Math. Programming Study 16 (1982),
  for relaxed steps; the inertia of a KKT matrix from Nocedal and Wright,
  Numerical Optimization (2nd ed., Springer 2006). All three cited from
  memory, to be confirmed before they are relied on. Implementation in
  FortFEM (doc/CONSTRAINED_NEWTON.md).
- No code of others was used. Reference data: VMEC++ 0.9.3 results of
  benchmarks/label-qa and benchmarks/flux-map-standard; DESC 0.17.4
  (github.com/PlasmaControl/DESC, MIT) for the QA axis
  (benchmarks/flux-map-standard/checks/desc_axis.py).
