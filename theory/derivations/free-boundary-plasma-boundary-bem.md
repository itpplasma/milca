# Free boundary: vacuum on the plasma boundary (axisymmetric)

- Question: replace the vacuum shell between Γ_p and Γ_c by a boundary integral on the moving plasma boundary Γ_p, so that no vacuum mesh exists, and keep the exact gradient and Hessian that the inertia-controlled Newton solve needs.
- Scope: axisymmetric equilibria (n = 0, cylindrical chart), static ideal MHD, fixed coil currents, fixed plasma current. The general 3D case needs the full-surface operator with a scalar potential (open, section 7).
- Status labels as in [first variation](free-boundary-first-variation.md): **proved**, **cited**, **hypothesis**, **open**.

## 1. Model

The Lagrangian at fixed coil currents I_k is ([first variation](free-boundary-first-variation.md) §2)

\[
L = W_p[X] - \frac{1}{\mu_0}\, e(\Gamma_p), \qquad e(\Gamma_p) = \frac12\int_{E} |B|^2\, dV + \text{const}, \tag{1}
\]

with E = ℝ³∖D the exterior of the plasma region D, and B the vacuum field in E. B is curl-free and divergence-free in E away from the coils, and B·n = 0 on Γ_p. Its circulation around the plasma is μ0 I_p, and its sources are the coil currents. Constants independent of Γ_p, such as coil self-energies, are dropped.

Symbols. Cylindrical (R, φ, Z). ψ is the poloidal flux per radian, B = ∇ψ × ∇φ + F ∇φ. F = R B_φ is the vacuum toroidal function (`f_vac`). G is the free-space Green function of div((1/R)∇·) (FortFEM `fortfem_axisymmetric_gs_kernel`), so a filament current I gives ψ = μ0 I G. The surface density is τ = μ0 K, with K the toroidal sheet current per unit arc length (A/m). s is arc length on Γ_p, and ψ_c = μ0 Σ_k I_k ψ_k is the coil flux.

## 2. Exterior energy as a boundary functional (proved)

The exterior field is the field of the coils plus a sheet current on Γ_p:

\[
\psi(x) = \psi_c(x) + \oint_{\Gamma_p} G(x,y)\,\tau(y)\, ds_y, \qquad \oint_{\Gamma_p}\tau\, ds = \mu_0 I_p, \qquad \psi|_{\Gamma_p} = \psi_b . \tag{2}
\]

- (2) gives Δ*ψ = 0 in E away from the coils, the circulation μ0 I_p and B·n = 0 (ψ constant on Γ_p). The decaying exterior solution with these data is unique, so (2) reproduces the vacuum field in E.
- Inside D, ψ is GS-harmonic (no current in the equivalent problem) with ψ = ψ_b on Γ_p, hence ψ ≡ ψ_b and the poloidal field vanishes in D.

The poloidal field energy over all space equals the energy of the currents, ½∫A·J dV. A_φ = ψ/R, and a toroidal current element covers 2πR. The coil self-terms are constant:

\[
\frac12\int_E |B_{pol}|^2 dV = \frac12\int_{\mathbb R^3} |B_{pol}|^2 dV
= 2\pi\Big[\tfrac12 \oint\oint \tau G \tau + \oint \tau\,\psi_c\Big] + \text{const}. \tag{3}
\]

The step uses B_pol = 0 in D and the reciprocity Σ_k I_k ψ_τ(x_k) = μ0^{-1}∮τψ_c.

The toroidal part has a Γ_p-independent total:

\[
\frac12\int_E \frac{F^2}{R^2} dV = \text{const} - \frac12\int_D \frac{F^2}{R^2}\, dV . \tag{4}
\]

The sheet current minimises (3) at fixed ∮τ = μ0 I_p, whose multiplier enforces ψ = ψ_b. V is symmetric positive definite (field energy of a current distribution), so the minimum is unique:

\[
e(\Gamma_p) = \min_{\oint\tau = \mu_0 I_p} 2\pi\Big[\tfrac12\langle \tau, \mathcal V\tau\rangle + \langle\tau,\psi_c\rangle\Big] - \frac{F^2}{2}\int_D |\nabla\varphi|^2 dV . \tag{5}
\]

## 3. First variation (proved)

(5) is the same function of Γ_p as the shell value e = E_v[X, Φ̃*] of [first variation](free-boundary-first-variation.md) §5: both equal ½∫_E|B|² up to constants for the same exterior field. Hence

\[
de[\xi] = -\oint_{\Gamma_p} \frac{|B_v|^2}{2}\, \xi\cdot n\, dS, \tag{6}
\]

and stationarity of L gives [[p + |B|²/2μ0]] = 0 on Γ_p as before (§6 there).

Check (circular wire, no coils, minor radius a ≪ R): e = μ0² I_p² R(ln(8R/a) − 2)/2. Then de/da = −μ0² I_p² R/(2a), and (6) with B_v = μ0 I_p/(2πa) on the area 4π²Ra gives the same value.

## 4. Discrete Lagrangian

Unknowns: the plasma map X (position gauge), τ (discontinuous P_d per panel of the plasma boundary chain, curved P_k panels from the boundary nodes of the plasma mesh) and the multiplier ψ̂ (= 2πψ_b):

\[
L_h = W_p - \frac{2\pi}{\mu_0}\Big[\tfrac12\,\tau^T V(X)\tau + \tau^T g(X)\Big] + \frac{\hat\psi}{\mu_0}\big(s(X)^T\tau - \mu_0 I_p\big) + \frac{F^2}{2\mu_0}\int_{D} |\nabla\varphi|^2 dV , \tag{7}
\]

- V_ij = ∮∮ G φ_i φ_j ds ds, g_i = ∮ ψ_c φ_i ds, and s_i = ∮ φ_i ds.
- V, g and s come from FortFEM `fortfem_gs_single_layer_shape_2d`, with a quadrature plan frozen per solve. Their first and second derivatives in the boundary node coordinates are exact derivatives of the discrete operator, from generated integrand jets.
- The last term is the shell leaf `gen/gen_vacuum_shell.f90` evaluated on plasma cells with the potential F q3, so it shares the shell's generated derivatives.

Stationarity in τ gives V τ + g = (ψ̂/2π) s, the Galerkin form of ψ = ψ_b. In ψ̂ it gives s^T τ = μ0 I_p.

## 5. KKT structure and inertia (proved for the vacuum block)

The vacuum block of the Newton matrix in (τ, ψ̂) is

\[
\frac{1}{\mu_0}\begin{pmatrix} -2\pi V & s \\ s^T & 0 \end{pmatrix}, \tag{8}
\]

with inertia (1 positive, n_τ negative): −2πV contributes n_τ negative eigenvalues, and the Schur complement s^T V^{-1} s/(2π) > 0 contributes one positive. The positivity of V in the discrete space is a per-state hypothesis, checked by the factorisation (FB-A2 analogue). The target inertia of the full matrix with k unstable directions of the reduced plasma Hessian is therefore

- positive: n_p + 1 − k,
- negative: n_τ + n_slide + k,

with n_p the plasma unknowns and n_slide the slide-gauge rows, as in the shell lane.

## 6. Coil flux from mgrid

The mgrid file has B per coil group on an (R, Z) grid and no vector potential. Per group k, ψ_k is fitted as a sum of ring sources outside the fit region,

\[
\psi_k(x) = \sum_j a_{kj}\, G(x, y_j), \tag{9}
\]

by least squares to (R B_Z, −R B_R) at the grid nodes inside the coil-free region around the plasma (method of fundamental solutions).

- (9) is exactly GS-harmonic in the fit region, and its jets come from the kernel.
- The fit residual is reported per group as part of the input error budget, as is the mgrid interpolation error in the VMEC codes.
- An additive constant in ψ_k changes L_h by a constant times μ0 I_p and drops out.

## 7. Open

- 3D: the full-surface operator (R2) on Γ_p with the multivalued scalar potential, and its shape Hessian.
- The discrete positivity of V per state (hypothesis above).
- The fit (9) and its error: number and placement of ring sources against the grid resolution.
