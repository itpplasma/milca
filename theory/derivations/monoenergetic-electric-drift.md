# Incompressible electric drift on a three-dimensional surface

## Question

The local monoenergetic transport operator needs electric precession at finite
magnetic ripple amplitude. The first reduction is the incompressible model of
Landreman et al., *Physics of Plasmas* **21**, 042503 (2014),
[equations (16), (18) and (20)](https://arxiv.org/html/1312.6058v2).
Its electric velocity-space terms and tangential magnetic drift are omitted.
It cannot establish physical loss torque or the defect of that reduction.

## Model and conventions

Let \(\psi_t\) be toroidal flux divided by \(2\pi\), \(\theta,\zeta\)
right-handed Boozer angles, \(\iota=1/q\) rotational transform, \(G,I\)
covariant field components, \(B=B_0 b>0\), \(v\) fixed speed,
\(\xi=v_\parallel/v\) pitch and \(\Phi(\psi_t)\) electric potential.
Primes denote \(d/d\psi_t\). SI units are used. The response \(g\) solves
\(Ag=s\), with the differential operator \(A\) defined by the left-hand side:

\[
\left[\xi\kappa b D_\parallel
-\frac{\kappa D_\parallel b}{2}(1-\xi^2)\partial_\xi
+\omega_E\frac{b^2}{\langle b^2\rangle}D_\perp
-\frac{\nu}{2}\partial_\xi(1-\xi^2)\partial_\xi\right]g=s. \tag{1}
\]

Here \(D_\parallel=\partial_\theta+q\partial_\zeta\),
\(D_\perp=\partial_\theta-(I/G)\partial_\zeta\),
\(\kappa=\iota B_0/(G+\iota I)\),
\(\nu=\nu_D/v>0\) is deflection frequency per speed and \(s\) a prescribed
drive divided by speed. Both \(\kappa\) and \(\nu\) have units of inverse
length. The signed electric parameter is
\(\omega_E=\Phi' G/[v(G+\iota I)]\), also in inverse metres.
The surface average is weighted by \(b^{-2}\).

## Derivation

The Boozer Jacobian \(\mathcal J\) and the covariant field give

\[
\mathcal J=(G+\iota I)/B^2. \tag{2}
\]
\[
(\mathbf B\times\nabla\psi_t)\cdot\nabla
=\mathcal J^{-1}(G\partial_\theta-I\partial_\zeta). \tag{3}
\]
\[
\mathbf v_E=\Phi'\mathbf B\times\nabla\psi_t/B^2. \tag{4}
\]
\[
\mathbf v_E\longmapsto\mathbf v_E^{\rm inc}
=\Phi'\mathbf B\times\nabla\psi_t/\langle B^2\rangle.
\quad\text{Incompressible approximation.} \tag{5}
\]
\[
v^{-1}\mathbf v_E^{\rm inc}\cdot\nabla
=\omega_E b^2\langle b^2\rangle^{-1}D_\perp. \tag{6}
\]
\[
\dot v_E\partial_v g\longmapsto0.
\quad\text{Monoenergetic approximation.} \tag{7}
\]
\[
\dot\xi_E\partial_\xi g\longmapsto0.
\quad\text{Electric pitch transport omitted.} \tag{8}
\]
\[
\mathbf v_m\cdot\nabla g\longmapsto0.
\quad\text{Magnetic drift retained only in the radial drive.} \tag{9}
\]

The radial derivative of the unknown is absent. Finite ripple amplitude enters
the streaming, mirror, drive and surface measure through the full supplied
\(b(\theta,\zeta)\). It is not expanded in the ripple amplitude.

For periodic functions \(u,w\), the electric operator \(A_E\) obeys

\[
\langle u A_E w\rangle
=\frac{\omega_E}{\langle b^2\rangle\int b^{-2}d\theta d\zeta}
\int u D_\perp w\,d\theta d\zeta. \tag{10}
\]
\[
\langle u A_E w\rangle=-\langle (A_Eu)w\rangle,
\qquad A_E1=0. \tag{11}
\]

The density constraint remains \(\langle g_0\rangle=0\), where \(g_0\)
is the zeroth Legendre coefficient. The bordered density-source multiplier
measures a discrete solvability defect. It is not a physical particle source.

At zero electric field a resolved closed-field-line harmonic can leave
additional density null modes. For example, \(q=4/5\) and three field
periods give \(D_\parallel\cos(12\theta-15\zeta)=0\). A single density
constraint does not remove that nullspace. The bordered matrix is rejected
when its estimated reciprocal one-norm condition number is at most
\((n+1)\epsilon_{\rm mach}\), where \(n\) is the number of surface points
and \(\epsilon_{\rm mach}\) machine precision. A successful factorization
alone does not establish a unique solution.

Use the inner product
\((u,w)=\langle\int_{-1}^1u^*w\,d\xi/2\rangle\), where the star denotes
complex conjugation. Let \(s_1=-\kappa(D_\perp b)(1+\xi^2)/(2b)\),
\(s_3=\xi b\), and \(d_{ij}=(s_i,A^{-1}s_j)\) on the constrained space.
With pitch reversal \(Pg(\xi)=g(-\xi)\), streaming and mirror terms are odd,
whereas collisions and electric advection are even. Therefore

\[
A(\omega_E)^\dagger=P A(-\omega_E)P. \tag{12}
\]
\[
d_{11}(\omega_E)=d_{11}(-\omega_E),\quad
d_{33}(\omega_E)=d_{33}(-\omega_E),\quad
d_{13}(\omega_E)=-d_{31}(-\omega_E). \tag{13}
\]

These relations use an even radial drive and an odd parallel drive with the
same flux pairing as `mono_coefficients_3d`. They are continuum identities.
Fourier product aliasing can leave a discrete defect in streaming and mirror
terms. Electric skew symmetry itself is exact to rounding on the chosen grid.

## Discretisation and derivatives

The electric term is diagonal in Legendre degree. Each angular diagonal block
acquires \(\omega_E E\), where \(E\) is the unit-electric operator. The same
term enters the bordered degree-zero block and the residual action. The
existing elimination and factor reuse remain valid.

The pointwise coefficient and its forward and reverse derivatives come from
`gen/gen_dk_mono.f90`. At fixed geometry the solution derivative \(g_E\) is

\[
A g_E+\mu_E e_0=-E g,\qquad\langle(g_E)_0\rangle=0. \tag{14}
\]

Here \(e_0\) is a spatially constant degree-zero source, and \(\mu_E\)
the derivative of its bordered source multiplier. The bordered factorisation
is reused. This derivative holds the drive, speed,
collision frequency and magnetic geometry fixed. It is not a total derivative
of a thermal flux with respect to a physical electric field.

## Verification gates

| Observable | Independent check | Acceptance |
|---|---|---|
| Operator action | Analytic Fourier-Legendre manufactured distribution | Relative maximum error below \(10^{-12}\) |
| Bordered inverse | Same analytic continuum distribution | Maximum error below \(2\times10^{-11}\) |
| Electric weighted skew symmetry and constant nullspace | Matrix identities | Relative maximum defect below \(10^{-13}\) |
| Electric derivative | Central difference, step \(10^{-6}\,\mathrm{m}^{-1}\) | Componentwise relative error below \(2\times10^{-6}\), derivative floor \(10^{-5}\) |
| Onsager-Casimir | Reversed electric field on a non-symmetric surface | Relative defect below \(2\times10^{-4}\) |
| Zero electric field | Existing default and explicit zero | Relative defect below \(10^{-13}\) |
| Closed field lines | Analytic additional density null mode at \(q=4/5\) | Null-mode action below \(10^{-12}\) and singular status |

## Independent coefficient reference

Let \(H=G+\iota I\), and let \(d_{ij}\) be the pitch-averaged coefficients
returned by MILCA. SFINCS's reference deflection frequency gives

\[
c_D=\frac{3\sqrt\pi}{4}
\left[\operatorname{erf}(1)-\frac{\operatorname{erf}(1)-2e^{-1}/\sqrt\pi}{2}\right],
\quad\nu=c_D\mathtt{nuPrime}\,B_0/H,
\quad\omega_E=\kappa\mathtt{EStar}. \tag{15}
\]

Its density-gradient and parallel-electric-field right-hand sides and moment
integrals give its native Fortran matrix \(M\) as

\[
M_{11}=-\frac{16H}{\sqrt\pi\iota^2B_0}d_{11},\quad
M_{21}=-\frac8{\sqrt\pi\iota}d_{31},\quad
M_{12}=\frac8{\sqrt\pi\iota}d_{13},\quad
M_{22}=\frac{4B_0}{\sqrt\pi H}d_{33}. \tag{16}
\]

The pitch integral is unaveraged in SFINCS and averaged in MILCA. The
unit-speed Maxwellian, speed quadrature and gyrophase integral contribute
\(2/\sqrt\pi\) per unaveraged pitch moment. The two right-hand-side scale
factors cancel their diagnostic normalization factors, leaving (16).
SFINCS's HDF5 matrix is transposed relative to its printed Fortran matrix.
The source definitions are `doc/manual/version3/runs.tex`, `sfincs_main.F90`,
`populateMatrix.F90`, `solver.F90` and `diagnostics.F90` of the public
revision recorded in the specification.

For both nonzero electric fields in that synthetic case, all four entries
agree within \(2\times10^{-9}\) relative. MILCA's independent angular and
pitch refinements from \(15\times9\times32\) to \(41\times23\times80\),
and SFINCS's combined refinement to \(31\times17\times64\), change the
coefficients by less than \(10^{-9}\) relative. The
[comparison record](../../test/data/sfincs_mono_exb_refinement.csv) retains
both codes' matrices and the independent refinements. SFINCS reports true relative
residuals below \(5\times10^{-8}\) in these runs. The regression tolerance is one
percent, fixed before inspecting the comparison.

At zero field the resolved rational harmonic makes both fine discretizations
singular. Their finite LU outputs are rejected as coefficient evidence.
This does not establish existence or uniqueness at rational surfaces.

`test_dk_mono_exb manufactured trace.csv` also writes the analytic action,
the computed action and the action without electric advection on one angular
cut, for a reproducible before/after figure. The units are inverse metres.

## Open physics gates

Low-collisionality resolution ladders remain required. Full electric velocity-space transport,
tangential magnetic drift, radial orbit coupling, self-consistent potential and physical wall losses are
separate extensions. Numerical certification of this new solve is unavailable.

[The multispecies surface solver](multispecies-electric-drift.md) adds
energy-coupled collisions in the same incompressible electric approximation.
