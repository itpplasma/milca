# Free boundary: vacuum on the plasma boundary in 3D (sheet current)

- Question: extend the vacuum on the plasma boundary
  ([axisymmetric note](free-boundary-plasma-boundary-bem.md)) to non-axisymmetric Γ_p, with
  the exact gradient and Hessian of the discrete energy, so that no vacuum mesh is needed in 3D.
- Scope: static ideal MHD, fixed coil currents, fixed plasma current, closed toroidal Γ_p with
  nfp-fold symmetry, cylindrical chart. Sources are published work only; the open-source codes
  consulted are listed with repository, commit and licence in section 7.
- Status labels as in [first variation](free-boundary-first-variation.md).

## 1. Model

As in the axisymmetric note, (1) there:

\[
L = W_p[X] - \frac{1}{\mu_0}\, e(\Gamma_p), \qquad e = \frac12\int_E |B|^2\,dV + \text{const}, \tag{1}
\]

with B the vacuum field in the exterior E: curl- and divergence-free away from the coils,
B·n = 0 on Γ_p, toroidal-current circulation μ0 I_p, and the coils at fixed currents.

## 2. Sheet-current representation (proved)

Extend B by zero into the plasma region D. The extended field is the field of the coils and of
the surface current K = n × [B]/μ0 = −n × B_+/μ0 on Γ_p (B_+ the exterior trace). K is
tangential and surface-divergence-free, because B_+·n = 0 and B_+ is curl-free at Γ_p.

- Conversely, a divergence-free tangential K for which coils plus sheet give zero field in D and
  B·n = 0 from outside is the sheet of the vacuum problem. Inside D the field of coils plus sheet
  is harmonic, tangent to Γ_p from inside and has zero circulation on both cycles exactly when the
  net sheet currents cancel the linked coil currents, so it vanishes.
- On the torus, divergence-free tangential fields are K = n × ∇_s Φ, Φ = φ + I_t u + I_p' v/2π
  with φ single valued, (u, v) ∈ [0, 1) × [0, 2π) the reference poloidal and toroidal parameters,
  and I_t, I_p' the net currents through a poloidal and a toroidal cut. The two net currents are
  fixed: I_t = I_p (the plasma current) and I_p' cancels the toroidal-field coil current linked by
  Γ_p, I_p' = −2πF/μ0, F = R B_φ of the coils (signs fixed by the orientation of (u, v), checked
  against the axisymmetric limit, section 6).
- In the parameters, with N = X_u × X_v outward and dA = |N| du dv,

\[
K\, dA = (\Phi_u X_v - \Phi_v X_u)\, du\, dv , \tag{2}
\]

  polynomial in the geometry derivatives; no normal or area element is formed.

## 3. Energy (proved)

The field energy of coils and sheet is ½∫A·J dV. With κ = μ0 K (tesla) and A_c the coil vector
potential (B_c = ∇ × A_c), dropping the coil self-energy,

\[
e = \frac{1}{8\pi}\oint\!\!\oint \frac{\kappa(x)\cdot\kappa(y)}{|x-y|}\,dA_x\,dA_y
  + \oint \kappa\cdot A_c\,dA . \tag{3}
\]

Because the field vanishes in D, (3) equals ½∫_E|B|² + const; the toroidal-field term that the
axisymmetric note adds by hand, −½F²∫_D|∇φ|², is contained in (3). Its first variation is
de[ξ] = −½∮|B_+|² ξ·n dA as before (same function of Γ_p).

The sheet minimises (3) over φ at fixed net currents (the self term is a positive inductance), so
L is concave in φ and no multiplier is needed.

Coil term without A_c. In the parameters κ dA = (Φ_u X_v − Φ_v X_u) du dv by (2), and for single-
valued φ, with (X_u·∇A·X_v − X_v·∇A·X_u) = (∇×A)·(X_u × X_v),

\[
\iint A\cdot(\phi_u X_v - \phi_v X_u)\,du\,dv = -\iint \phi\, B\cdot N\,du\,dv . \tag{4a}
\]

The ramps of Φ = φ + I_t u + I_p' v/2π are not integrated by parts: Φ_u = φ_u + I_t and
Φ_v = φ_v + I_p'/2π are periodic, and directly

\[
\oint \kappa\cdot A_c\,dA = -\iint \phi\, B_c\cdot N\,du\,dv
  + I_t \int_0^1 \Psi_t(u)\,du - \frac{I_p'}{2\pi}\int_0^{2\pi} \Psi_p(v)\,dv , \tag{4}
\]

with Ψ_t(u) = ∮ A_c·X_v dv the flux through the toroidal loop at u and Ψ_p(v) = ∮ A_c·X_u du the
flux through the section at v. Integrating by parts instead moves the ramps to a cut, where they
jump; on a spectral grid that jump spreads into oscillating node forces. Both fluxes are coil fluxes
through spanning surfaces in D, which is coil-free:

- Ψ_p(v): the fan P(s, u) = L(v) + s(X(u, v) − L(v)) in the meridional plane at v,
  Ψ_p = ∫∫ s B_φ (ΔZ R_u − ΔR Z_u) ds du;
- Ψ_t(u): the flux through a fixed loop L inside D (the initial magnetic axis) plus the ribbon
  P(s, v) = L(v) + s(X(u, v) − L(v)), ∫∫ B_c·(P_s × P_v) ds dv. The flux through L is constant
  and drops out with the fixed I_t.

The net currents follow from Ampère's law in D: I_t = μ0 I_φ along +v (κ units) and, for a
counter-clockwise chain in (R, Z), I_p' = −2πF, which cancels the coil circulation 2πF of B_c on a
toroidal loop in D (o = −1 for a clockwise chain flips the sign).

So the energy needs B_c and its jets only inside D and on Γ_p, which the mgrid interpolant
supplies.

## 4. Discrete energy

- Grid: the boundary is sampled at n_u = 2 n_b + 1 equidistant points of the chain parameter (odd,
  so no Nyquist mode) and n_v per field period, from the P_k chain and the map's Fourier modes, and
  taken as its trigonometric interpolant: u derivatives spectral, v derivatives from the modes. All
  grid values are linear in the boundary coefficients. The cells are at most twice as long
  toroidally as poloidally.
- Self term: the corrected trapezoidal rule (FortFEM `fortfem_toroidal_sheet_trapezoid`): the
  punctured double sum minus Z_i|k_i|² per point, Z_i the Epstein zeta function at s = 1/2 of the
  local lattice form |m h_u X_u + n h_v X_v|²; O(h³) for the singular inner sum, spectral outer sum.
- Coil term: (4) by the trapezoidal rule on the same grid, Gauss in s along fans and ribbons; the
  integrands and their derivatives are generated with FortSym (gen/gen_sheet_flux.f90).
- φ: Fourier modes in u up to n_u/4 and the map's toroidal modes in v; the constant is the gauge.
  The rule is accurate only for densities the grid resolves: with φ nodal on the grid its block is
  indefinite (eigenvalues down to −0.3 of the largest on the NCSX boundary). e is quadratic in φ; its
  φ block is assembled in this basis and φ is eliminated exactly with its Cholesky factor at every
  state.

## 5. Newton system

Unknowns: the plasma map (position gauge). The reduced vacuum Hessian H_xx − H_xφ M⁻¹ H_φx enters
as a product (two Hessian-vector products of the trapezoidal energy and one Cholesky solve). The
Newton step solves the bordered system by GMRES, right-preconditioned by the factor of the plasma
block, whose inertia control sets the shift; target inertia with k unstable directions: positive
n_p − k, negative n_slide + k.

## 6. Checks

Test `ideal_fb_sheet_vacuum` (P2 chain around an ellipse):

| check | result |
|---|---|
| F alone: L_v = πF²∫∫dR dZ/R/μ0, the toroidal-field term of the axisymmetric form | 6.2e-7 |
| axisymmetric limit (PF rings, F = 0): L_v − I_t Ψ(L)/μ0 against the single layer | 8.3e-6 (24 nodes), 8.7e-7 (48) |
| shape derivative along smooth displacements against the single layer | 1.3e-4, 1.2e-5 |
| 3D surface (nfp = 3) in the window-coil mgrid field: gradient, Hessian-vector product, reduced product vs central differences | 1.8e-9, 3.9e-10, 3.9e-10 |

On a toroidal-field coil the rule converges at third order (4.9e-3, 6.1e-4, 7.6e-5, 9.5e-6,
1.2e-6 for 8 … 128 poloidal points); on cells of aspect 4 the error is 1.3e-3 against 3e-5 at
aspect 2.

## 7. Sources

Published methods: the relative-coordinate pair rule follows S. A. Sauter and C. Schwab, Numer.
Math. 78, 211 (1997), doi:10.1007/s002110050311; surface-current potentials on toroidal surfaces as
in NESCOIL and REGCOIL. No code was copied. Open-source codes read for the survey of vacuum
treatments (repository, commit, licence):

| code | repository | commit | licence | read for |
|---|---|---|---|---|
| BIEST | github.com/dmalhotra/BIEST | 2aeea0d | Apache-2.0 | partition-of-unity singular quadrature |
| DESC | github.com/PlasmaControl/DESC | d644e62 | MIT | free-boundary objective and quadrature |
| VMEC++ | github.com/proximafusion/vmecpp | 95ec135 | MIT | NESTOR vacuum |
| VMEC2000 | github.com/hiddenSymmetries/VMEC2000 | 70c637a | MIT | NESTOR vacuum |
| SIMSOPT | github.com/hiddenSymmetries/simsopt | 912bdcc | MIT | virtual casing |
| virtual-casing | github.com/hiddenSymmetries/virtual-casing | 81cf137 | Apache-2.0 | virtual casing |
| REGCOIL | github.com/landreman/regcoil | 528c5d0 | BSD-2-Clause | current potential on a torus |
| FreeGS | github.com/freegs-plasma/freegs | 5b41fe5 | LGPL-3.0 | study only |
| FMM3D | github.com/flatironinstitute/FMM3D | 0e13b3d | Apache-2.0 | licence only |
| ZetaTrap3D | github.com/bobbielf2/ZetaTrap3D | 26052f4 | MIT | licence only |
