# Cartesian flux map: the ideal lane as a map into the flux disc

The ideal lane of [equilibrium-roadmap.md](equilibrium-roadmap.md) in the
form that is developed from here on. The unknown is the map from physical
space to a model torus in which the field is fixed by the profiles alone.
It replaces the pair (χ, g) with axis unknowns of
[eulerian-label-plan.md](eulerian-label-plan.md), which stays as reference.

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

### Prototype (Python, C0 Lagrange elements, not in the repository)

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

- F1 (done): axisymmetric, this document.
- F2: 3D on the per-plane meshes of the label lane (ζ cylindrical), the 3D
  projection space, QA, QH and QI decks against VMEC++.
- F3: mesh swept along a guiding curve with ζ its parameter, for axes that
  are not single valued in the cylindrical angle (figure-eight) and other
  exotic shapes; section meshes by Winslow.
- Then free boundary (the wall becomes Γ_p, exterior and shape derivatives
  of gates E6) and regions (one map per island or lobe, separatrices as
  aligned interfaces, label plan section 4a).

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
  balance, so alignment returns in another form (roadmap step 5);
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

- Projection space in 3D. Candidate: B_h = Π_2(u_h*β), the 2-form space of
  the complex, in mixed form; the per-cell mean of F1 and its first-order
  b̄ are the fallback.
- Globalisation at a strong guide field from a poor start: continuation in
  the guide field works; a merit function that tolerates the curved
  constraint R τ det ∇u = f(s) (mixed Newton) is the candidate for cold
  starts.
- Gauge without a symmetry: constraints orthogonal to θ* → θ* + k(s)
  instead of the symmetry rows.
- B_pol at 2 times the projected exact label on 8 cells (F1).
- Section meshes: harmonic extension of the boundary folds on a strongly
  elongated bean (QA); Winslow (`milca_spline_ale`) or map2disc
  (roadmap, optional).
- The relation of the Hessian to the ideal δW (section 1).
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
