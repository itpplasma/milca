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

Coil term without A_c. With Φ = φ + I_t u + I_p' v/2π, integrating ∮∇_sΦ·(A_c × n) by parts on
the torus cut along u = 0 and v = 0, and using ∇_s·(A_c × n) = B_c·n,

\[
\oint \kappa\cdot A_c\,dA = \mu_0\Big[\oint \phi\, B_c\cdot n\, dA + I_t\,\Psi_t + I_p'\,\Psi_p\Big], \tag{4}
\]

with Ψ_t the coil flux through a surface spanning the toroidal loop u = 0 and Ψ_p the coil flux
through the poloidal cross-section v = 0 (orientation by the cut). Both spanning surfaces can lie
in D, which is coil-free:

- Ψ_p: the cross-section of D at v = 0, the image of the reference disc under the plasma map;
- Ψ_t: the ribbon from a fixed loop L_0 inside D (the reference axis) to the loop u = 0, plus the
  flux through L_0, which is constant and drops out with the fixed I_t.

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
- Coil term: (4) with φ B_c·n on the same surface quadrature, Ψ_p on the plasma cross-section
  quadrature and Ψ_t on a ribbon quadrature, all with B_c jets.

## 5. Newton system

Unknowns: plasma map (position gauge) and φ. The vacuum block in φ is −(1/μ0) times the Galerkin
inductance matrix, negative definite. Target inertia with k unstable directions: positive
n_p − k, negative n_φ + n_slide + k.

## 6. Checks

- Axisymmetric limit: for an axisymmetric Γ_p and coils, (3)–(4) must give the energy, gradient
  and solution of the axisymmetric single-layer formulation (same equilibrium, the tokamak and
  D-shape decks).
- Field inside D: the coils plus the solved sheet give |B| at interior points small against |B_+|.
- Central differences of the gradient and Hessian for every block.
- Vacuum oracle: Γ_p a flux surface of a known vacuum field (coils plus a current filament) gives
  zero pressure-balance defect.

## 7. Sources

To be completed with the survey of published methods and open-source codes (repository, commit,
licence) before the first landing.
