# Coil field split in the fixed-surface free boundary

Question: with the vacuum coupled through a fixed surface Γ_c, can the coils
lie anywhere, including inside Γ_c and inside the shell between the plasma
boundary Γ_p and Γ_c?

## Model

Ω_p is the plasma, Ω_o = ℝ³ \ Ω_p the vacuum, n the normal of Γ_p pointing
out of Ω_p. The vacuum field is split into the known field B_c of all coils
(including the toroidal-field coils) and the field of the plasma current,

$$
\mathbf B = \mathbf B_c + \nabla\Phi_p ,\qquad
\Phi_p = c\,\frac{\theta}{2\pi} + \tilde\Phi , \tag{1}
$$

with c the poloidal circulation (μ₀ I_p), θ a poloidal angle about a curve
inside Ω_p and Φ̃ single-valued, harmonic in Ω_o and decaying at infinity.
∇·B_c = 0 everywhere, and ∇×B_c = μ₀J_c vanishes outside the coils.

The vacuum co-energy at fixed coil currents is

$$
E_v = \frac12\int_{\Omega_o} |\mathbf B|^2 \, dV . \tag{2}
$$

## Derivation

Expanding (1) in (2):

$$
E_v = \frac12\int_{\Omega_o} |\nabla\Phi_p|^2
+ \int_{\Omega_o} \mathbf B_c\cdot\nabla\Phi_p
+ \frac12\int_{\Omega_o} |\mathbf B_c|^2 . \tag{3}
$$

First term, with Ω_s the shell between Γ_p and Γ_c and 𝒮 the exterior
Steklov–Poincaré map of Γ_c:

$$
\frac12\int_{\Omega_o} |\nabla\Phi_p|^2
= \frac12\int_{\Omega_s} |\nabla\Phi_p|^2 + \frac12\langle\mathcal S u, u\rangle_{\Gamma_c},
\qquad u = \tilde\Phi|_{\Gamma_c} . \tag{4}
$$

This is the present shell energy without the coil traces.

Second term, single-valued part. ∇·B_c = 0 in Ω_o, Φ̃ = O(1/r) and
B_c = O(1/r³) at infinity:

$$
\int_{\Omega_o} \mathbf B_c\cdot\nabla\tilde\Phi
= -\oint_{\Gamma_p} \tilde\Phi\, \mathbf B_c\cdot\mathbf n \, dS . \tag{5}
$$

Second term, multivalued part. Cut Ω_o along a surface Σ from Γ_p to
infinity, across which θ jumps by 2π:

$$
\int_{\Omega_o} \mathbf B_c\cdot\nabla\frac{c\,\theta}{2\pi}
= -\oint_{\Gamma_p} \frac{c\,\theta}{2\pi}\, \mathbf B_c\cdot\mathbf n \, dS
+ c\int_\Sigma \mathbf B_c\cdot d\mathbf S . \tag{6}
$$

The flux of B_c through Σ is the coil flux linked between Γ_p and infinity.
In axisymmetry it is −2π ψ_c(x_Σ), with ψ_c the coil poloidal flux per
radian at the point x_Σ where Σ meets Γ_p.

Third term, with Ω_c the interior of Γ_c:

$$
\frac12\int_{\Omega_o} |\mathbf B_c|^2
= \frac12\int_{\mathbb R^3\setminus\Omega_c} |\mathbf B_c|^2
+ \frac12\int_{\Omega_c} |\mathbf B_c|^2
- \frac12\int_{\Omega_p} |\mathbf B_c|^2 . \tag{7}
$$

The first two terms of (7) do not depend on the plasma, because Γ_c and the
coil currents are fixed. They are infinite for filament coils but drop out
of every variation.

## Result

Up to a plasma-independent constant,

$$
E_v = \frac12\int_{\Omega_s} |\nabla\Phi_p|^2 + \frac12\langle\mathcal S u, u\rangle_{\Gamma_c}
- \oint_{\Gamma_p} \Phi_p\, \mathbf B_c\cdot\mathbf n \, dS
+ c\int_\Sigma \mathbf B_c\cdot d\mathbf S
- \frac12\int_{\Omega_p} |\mathbf B_c|^2 . \tag{8}
$$

- B_c enters only on Γ_p, on Σ and in Ω_p, where it is smooth. Γ_c carries
  no coil data, so coils and vessel may lie inside Γ_c and cross the shell.
- Γ_c can be chosen once per device with a wide margin, and 𝒮 is assembled
  and factored once.
- The shape derivative of (8) in Γ_p needs B_c and ∇B_c on Γ_p and in Ω_p,
  which the coil-field spline provides.

## Open

- Discrete check of (8) against the present trace formulation on a case
  whose coils lie outside Γ_c (the two must agree to quadrature accuracy).
- The cut Σ in 3D, where the linked coil flux is a surface integral rather
  than a flux-function value.
