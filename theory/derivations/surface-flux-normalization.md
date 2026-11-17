# SI fluxes and toroidal torque from the surface matrix

The matrix returned by `surface_transport` needs geometric and velocity-moment
factors before it gives particle flux or torque. `surface_fluxes` applies
those factors to each species and keeps the parallel electric force distinct
from the potential gradient.

Let \(\psi_t\) be toroidal magnetic flux divided by \(2\pi\), with units
of webers, and \(\theta,\zeta\) Boozer angles. Use

\[
\mathbf B=\nabla\psi_t\times\nabla\theta
+\iota\nabla\zeta\times\nabla\psi_t,\qquad
\iota=1/q,\qquad B=B_0b. \tag{1}
\]

Here \(q\) is safety factor, \(B_0\) a field scale, and \(G,I\) the
covariant toroidal and poloidal magnetic-field components. Define

\[
H=G+\iota I,\qquad c=qG/H=\frac{q}{1+(I/G)/q}. \tag{2}
\]

Species \(a\) has density \(n_a\), temperature \(T_a\) in joules, mass
\(m_a\), charge \(Z_ae\), and thermal speed
\(v_{ta}=\sqrt{2T_a/m_a}\), with elementary charge \(e>0\).
Primes below mean \(d/d\psi_t\). The potential is \(\Phi\).
Let \(E_0\) specify the parallel field as \(E_\parallel=E_0b\).
The three forces are

\[
A_{a1}=\frac{n_a'}{n_a}+\frac{Z_ae}{T_a}\Phi'
-\frac32\frac{T_a'}{T_a},\quad
A_{a2}=\frac{T_a'}{T_a},\quad A_{a3}=\frac{Z_ae}{T_a}E_0. \tag{3}
\]

The first two forces have units of inverse webers and the third inverse
metres. A normalized radial coordinate needs its flux Jacobian before
its gradients can be passed here. The transport matrix must be calculated
at the physical electric rate \(\Omega_E=\Phi'G/H\) for an actual plasma
response. Independent matrix-basis tests can prescribe the forces separately.

## Distribution and moment factors

Write normalized speed \(x=v/v_{ta}\), pitch \(\xi=v_\parallel/v\), and
\(r=[\partial_\theta b-(I/G)\partial_\zeta b]/(qb)\). The radial magnetic
drift in the local model is

\[
\mathbf v_{ma}\cdot\nabla\psi_t
=-\frac{T_a}{Z_ae}\,c\,x^2(1+\xi^2)r. \tag{4}
\]

Let \(h_a^{bj}\) be the response to drive \(j\) of species \(b\) used in
[the matrix derivation](multispecies-electric-drift.md). The physical
Maxwellian-weighted perturbation is

\[
\frac{f_{1a}}{F_a}
=\sum_b\left[c h_a^{b1}A_{b1}+c h_a^{b2}A_{b2}
+h_a^{b3}A_{b3}\right],\qquad
F_a=\frac{n_a e^{-x^2}}{\pi^{3/2}v_{ta}^3}. \tag{5}
\]

The factor \(c\) in the first two drives follows from (4). There is no
such factor for the parallel drive. For Legendre polynomials \(P_l\),

\[
\int_{-1}^1(1+\xi^2)\sum_lh_lP_l(\xi)\,d\xi
=\frac83h_0+\frac4{15}h_2. \tag{6}
\]
\[
\int_{-1}^1\xi\sum_lh_lP_l(\xi)\,d\xi=\frac23h_1. \tag{7}
\]

The matrix's moment normalization is half the corresponding physical
velocity integral. Combining (4)-(7) gives, for the matrix \(D_{ab,ij}\),

\[
\Gamma_a\equiv\langle\boldsymbol\Gamma_a\cdot\nabla\psi_t\rangle
=-2n_ac\sum_b[c(D_{ab,11}A_{b1}+D_{ab,12}A_{b2})+D_{ab,13}A_{b3}]. \tag{8}
\]
\[
Q_a\equiv\langle\mathbf Q_a\cdot\nabla\psi_t\rangle
=-2n_aT_ac\sum_b[c(D_{ab,21}A_{b1}+D_{ab,22}A_{b2})+D_{ab,23}A_{b3}]. \tag{9}
\]
\[
U_a\equiv\langle b u_{\parallel a}\rangle
=2\sum_b[c(D_{ab,31}A_{b1}+D_{ab,32}A_{b2})+D_{ab,33}A_{b3}]. \tag{10}
\]

Brackets denote the flux-surface average with weight \(b^{-2}\).
The units of \(\Gamma_a,Q_a,U_a\) are
\(\mathrm{Wb\,m^{-3}s^{-1}}\),
\(\mathrm{J\,Wb\,m^{-3}s^{-1}}\), and \(\mathrm{m\,s^{-1}}\).
Equation (9) is kinetic-energy flux, before any enthalpy subtraction.
The receiving density \(n_a\), rather than the drive-species density
\(n_b\), multiplies the response.

## Current and field torque

Let \(\mathbf e_\zeta=\partial\mathbf x/\partial\zeta\) and
\(\mathbf j\) be electric current density. Equation (1) gives

\[
\mathbf B\times\mathbf e_\zeta=\iota\nabla\psi_t. \tag{11}
\]
\[
J_\psi\equiv\langle\mathbf j\cdot\nabla\psi_t\rangle
=\sum_a Z_ae\Gamma_a. \tag{12}
\]
\[
\tau_{\zeta,a}=Z_ae\Gamma_a/q,\qquad
\tau_\zeta=\langle(\mathbf j\times\mathbf B)\cdot\mathbf e_\zeta\rangle
=J_\psi/q. \tag{13}
\]

Torque is positive along increasing \(\zeta\) and has units of
newton metres per cubic metre. It is field torque on the plasma.
An ambipolar zero requires the total current in (12) to vanish, not each
species flux. Wall-carried momentum, external sources and momentum transport
need their own balance. They are not additional copies of (13).

## Verification

FortSym generates (3), (8)-(10), (12)-(13), and their forward and reverse
derivatives. It verifies the three pitch integrals in (6)-(7) by exact
differentiation and endpoint evaluation before generating the kernels.
The native checks compare all three SI moments for each matrix
drive against SFINCS, test a two-species Boltzmann equilibrium and the
receiving-density factor, and check generated derivatives by centered
differences and tangent-adjoint duality. The fixture's source, normalizations
and refinement evidence are recorded in
[its reference note](../../test/data/surface_fluxes_sfincs.md).

The conversion does not add electric velocity-space terms, tangential
magnetic drift, global orbit losses or a rotation-evolution equation to the
incompressible surface model.
