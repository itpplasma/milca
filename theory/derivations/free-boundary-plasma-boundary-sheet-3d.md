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

The ramps of Φ = φ + I_t u + I_p' v/2π on the torus cut at u = 0 and v = 0 give in addition the
boundary terms of the same integration by parts, so that

\[
\oint \kappa\cdot A_c\,dA = -\iint_{[0,1)\times[0,2\pi)} \Phi\, B_c\cdot N\,du\,dv
  + I_t\,\Psi_t - I_p'\,\Psi_p , \tag{4}
\]

with Ψ_t = ∮ A_c·dl along the loop u = 0 in the direction of increasing v and Ψ_p = ∮ A_c·dl
along the loop v = 0 in the direction of increasing u. Both are coil fluxes through spanning
surfaces, which can lie in D, which is coil-free:

- Ψ_p: the fan P(s, u) = c + s(X(u, 0) − c) from a point c of D in the plane v = 0, with
  Ψ_p = ∫∫ B_c·(P_s × P_u) ds du = ∫∫ s B_φ (ΔZ R_u − ΔR Z_u) ds du;
- Ψ_t: the flux through a fixed loop L inside D (the initial magnetic axis) plus the ribbon
  P(s, v) = L(v) + s(X(0, v) − L(v)), ∫∫ B_c·(P_s × P_v) ds dv. The flux through L is constant
  and drops out with the fixed I_t.

The net currents follow from Ampère's law in D: I_t = μ0 I_φ along +v (κ units) and, for a
counter-clockwise chain in (R, Z), I_p' = −2πF, which cancels the coil circulation 2πF of B_c on a
toroidal loop in D (o = −1 for a clockwise chain flips the sign).

So the energy needs B_c and its jets only inside D and on Γ_p, which the mgrid interpolant
supplies.

## 4. Discrete energy

- Geometry: X(u, v) = H(q(u, v)) with q the plasma map on its boundary chain (P_k panels in u,
  real Fourier modes in v); in the position gauge q3 = v, so X is linear in the boundary
  coefficients for the cylindrical chart.
- φ: continuous P_k in u on the same chain, real Fourier modes in v (n ∈ nfp ℤ, |n| ≤ n_φ).
- Self term: Galerkin double integral of (2)·(2)/|X − Y| over pairs of patches (u-panel ×
  v-interval). At fixed quadrature points in the patch parameters every input is linear in the
  coefficients; near coincidence X − Y = ρ Q with Q a divided difference, so the integrand is
  κ_x·κ_y/|Q| with ρ absorbed in the singular rule. The value, gradient and Hessian of this leaf
  are generated with FortSym; the derivatives of the discrete energy are exact.
- Coil term: (4) with Φ B_c·N on the same surface quadrature (Gauss per patch), Ψ_p on the fan
  and Ψ_t on the ribbon (Gauss in s and along u or v), all with B_c jets; the integrands and their
  derivatives in the boundary inputs are generated with FortSym (gen/gen_sheet_flux.f90).
- Gauge: K vanishes for constant φ, so φ of mode 1 at chain node 1 is fixed.

## 5. Newton system

Unknowns: plasma map (position gauge) and φ. The vacuum block in φ is −(1/μ0) times the Galerkin
inductance matrix, negative definite. Target inertia with k unstable directions: positive
n_p − k, negative n_φ + n_slide + k.

## 6. Checks

Test `ideal_fb_sheet_vacuum` (P2 chain around an ellipse):

| check | result |
|---|---|
| F alone: L_v = πF²∫∫dR dZ/R/μ0, the toroidal-field term of the axisymmetric form | 1.8e-6 |
| axisymmetric limit (PF rings, F = 0): L_v − I_t Ψ(L)/μ0 against the single layer | 4.3e-7 (24 nodes), 2.3e-7 (48) |
| node gradient against the single layer (different discretisations, O(h²)) | 2.7e-3, 6.9e-4 |
| 3D surface (nfp = 3, one toroidal mode) in a fitted window-coil field: gradient, Hessian vs central differences | 4.9e-10, 7.4e-10 |

The toroidal intervals are at most three poloidal panel lengths long; the energy error decays
exponentially with nv (2.4e-2, 2.4e-3, 1.6e-4, 2.0e-6, 2.4e-9 for nv = 4 … 64 on the ellipse).

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
