# Cartesian flux map: the ideal model as a map into the flux disc

One of the two discretisations of the ideal model in
[equilibrium-roadmap.md](equilibrium-roadmap.md); the other is its inverse,
the Lagrangian map. The unknown is the map from physical space to a model
torus in which the field is fixed by the profiles alone. Its polar form,
the pair (χ, g) with axis unknowns of
[eulerian-label-plan.md](eulerian-label-plan.md), is superseded. What is
specific to this direction of the map and what is not:
[audit](equilibrium-lane-audit.md), section 4.

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
the choice of ζ, which is fixed by giving ζ.

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
an identity shift, accepted on energy decrease. Gauge by up-down symmetry:
penalty rows y_h = 0 on Z = 0 and Θ_h = 0 at the outboard wall node.

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
(`test/support/flux_map_sections.f90`): `layout_ogrid` and
`layout_set_curve(centre_core)` on the area centroid of each wall curve; the
start is the polar map of the section after an affine normalisation by its
second moments. No axis unknown, multiplier, winding form or singular
quadrature. The toroidal flux through every plane is conserved to
round-off in all runs below (defect 1e-14 to 6e-10) and B_φ keeps its sign.

F2a (`flux_map_3d_gate_f2a`, 5 s): the F1 Solov'ev case through the 3D path
on 3 identical planes reproduces F1: B 2.420e-2, 7.110e-3, 1.932e-3 on 2, 4,
8 cells (F1: 2.42e-2, 7.11e-3, 1.93e-3), 13, 17, 26 Newton steps, axis
5.0e-3, 1.1e-4, 1.3e-5.

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
therefore runs 2 and 4 cells from the polar start; options of the test:
`nested`, `stationary`, degree, planes.

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
  standard configuration with pressure (W7-X, which does not mesh yet,
  below).

QA driver (`test_flux_map_qa`, opt-in, not a registered test): the
Landreman–Paul QA deck of [benchmarks/label-qa](../benchmarks/label-qa/README.md),
one field period, p = 0, ι prescribed, gauge by stellarator symmetry,
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
- W7-X standard (`benchmarks/flux-map-standard/w7x`, pressure, ι 0.878 to
  0.943 prescribed from VMEC++) does not mesh: `spline_quadrature_build`
  reports a singular block map on the bean plane φ = 0 (the O-grid corner
  quadrilaterals are valid; the interior folds). A smaller core does not
  help. This contradicts the flux map's premise of robustness on shaped
  sections and is the first open item of the section discretisation.
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
- B_pol at twice the projected exact label on 8 cells, order 1.8–1.9 (F1).
- Gauge without a symmetry: constraints orthogonal to θ* → θ* + k(s)
  instead of the symmetry penalty rows (untested in 2D and 3D).
- Cold starts at a strong guide field: continuation works; a step
  acceptance that tolerates the curved constraint R τ det ∇u = f(s) is
  owed.
- Projection in 3D: none is stable so far and it is off the critical path
  (roadmap P0b). A candidate is B_h = Π_2(u_h*β) with a commuting
  projection, in mixed form.
- Section meshes on strongly elongated or concave sections: Winslow
  (`milca_spline_ale`) or map2disc.
- The relation of the Hessian to the ideal δW (section 1) is argued, not
  verified.
- Not run: pressure in 3D, QH and QI decks, free boundary, any region
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
- No code of others was used. Reference data: VMEC++ results of
  benchmarks/label-qa.
