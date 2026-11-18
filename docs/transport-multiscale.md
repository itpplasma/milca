# Transport and multiscale coupling

Status: canonical transport and multiscale specification. MILCA implements
the flux-shell toroidal momentum balance below.

## Principle

Transport is a **macro evolution problem** whose closure can range from an
empirical coefficient to the direct flux of a global nonlinear kinetic solve.
MILCA should support this entire ladder without forcing all microphysics into a
scalar diffusion coefficient.

For transported macro variables \(u_a\),

\[
\partial_t u_a + \nabla\cdot\Gamma_a = S_a.
\]

The transport state may contain, as demanded by the chosen model,

\[
n_s,\ T_s,\ \text{rotation/momentum},\ j_\parallel,\ \psi,\
n_Z,\ \text{fast-particle moments},\ldots
\]

and their 1D, 2D or 3D spatial dependence.

## Flux-shell toroidal momentum balance

The balance advances the mechanical toroidal angular momentum inside fixed
shells. Let \(\mathbf x\) be Cartesian position, \(\hat{\mathbf z}\) the
chosen rotation axis and \(\mathbf k=\hat{\mathbf z}\times\mathbf x\) the
generator of positive toroidal rotation. Let \(\mathbf p\) be mechanical
momentum per unit volume, \(\boldsymbol{\Pi}\) its symmetric momentum-flux
tensor and \(\mathbf f\) the force per unit volume on the modeled species.
The momentum moment gives

\[
\begin{aligned}
\partial_t\mathbf p+\nabla\cdot\boldsymbol{\Pi}&=\mathbf f, \tag{1}\\
\boldsymbol{\Pi}:\nabla\mathbf k&=0, \tag{2}\\
\partial_t(\mathbf k\cdot\mathbf p)
  +\nabla\cdot(\boldsymbol{\Pi}\mathbf k)&=\mathbf k\cdot\mathbf f. \tag{3}
\end{aligned}
\]

Equation (2) follows from symmetry of the tensor and antisymmetry of
\(\nabla\mathbf k\). The species momentum equation and stress convention are
given in [R. Fitzpatrick, *Plasma Fluid Theory*, equations (4.94) and (4.97)](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node53.html).

For shell \(i\), let \(L_i\) be its volume integral of
\(\mathbf k\cdot\mathbf p\), in \(\mathrm{kg\,m^2\,s^{-1}}\), and \(T_i\)
its integral of \(\mathbf k\cdot\mathbf f\), in \(\mathrm{N\,m}\). At face
\(j\), let \(F_j\) be the integral of
\(\mathbf n\cdot\boldsymbol{\Pi}\mathbf k\), in \(\mathrm{N\,m}\), with the
normal \(\mathbf n\) pointing toward the exterior on every face. Faces
\(i\) and \(i+1\) are the inner and outer faces of shell \(i\).
For a time step \(\Delta t\) the overbar denotes a time average over the step.

\[
\begin{aligned}
L_i^{n+1}&=L_i^n+\Delta t\,(\overline F_i-\overline F_{i+1}
  +\overline T_i), \tag{4}\\
\sum_{i=1}^{N}(L_i^{n+1}-L_i^n)
  &=\Delta t\,\left(\overline F_1-\overline F_{N+1}
  +\sum_{i=1}^{N}\overline T_i\right). \tag{5}
\end{aligned}
\]

Here \(n\) is the step index and \(N\) the number of shells. One shared flux
per interior face gives cancellation in (5). A shell that starts on a regular
magnetic axis has zero inner face flux. A positive outer edge flux removes
positive toroidal angular momentum. Reversing the toroidal direction reverses
\(L_i\), \(F_j\) and \(T_i\). It leaves the face normals directed outward.

In a flux coordinate \(s\) increasing outward, integrate a flux-surface-averaged
torque density over each shell with the volume measure \((dV/ds)\,ds\), where
\(V(s)\) is the enclosed volume. For neoclassical toroidal viscosity with
finite orbit width, the input is the real-space deposited torque. When the
momentum tensor includes a stress that carries an exchange through the shell
faces, the associated torque enters the ledger through that flux.
For angular velocity constant on a shell, \(\Omega_i=L_i/I_i\), with
\(I_i=\int_i \rho R^2\,dV\), mass density \(\rho\) and cylindrical distance
\(R\) from the chosen axis. Changes in inertia alter \(\Omega_i\) even in a
shell with conserved \(L_i\).

The axisymmetric kinetic closure gives the local toroidal angular velocity
\(\Omega_\phi(R)=\omega+K/R^2\), where \(\omega\) is its rigid rotation term
and \(K\) its parallel-flow correction. If both coefficients are held constant
across shell \(i\), then \(L_i=I_i\omega+M_iK\) with shell mass
\(M_i=\int_i\rho\,dV\). The rotation \(L_i/I_i\) then equals
\(\Omega_\phi\) at \(R=\sqrt{I_i/M_i}\). [Axisymmetric profile](../app/ntv_axisymmetric_profile.f90)

For a species of charge \(q_a\), number density \(n_a\), pressure \(p_a\)
and mean velocity \(\mathbf v_a\), the radial projection of its momentum
equation gives

\[
E_r=\frac{\hat{\mathbf e}_r\cdot\nabla p_a}{q_a n_a}
    -\hat{\mathbf e}_r\cdot(\mathbf v_a\times\mathbf B), \tag{6}
\]

where \(\hat{\mathbf e}_r\) is the outward unit normal, \(\mathbf B\) the
magnetic field and \(E_r\) the electric field along the normal. This model
neglects radial inertia, viscous forces and collisional friction.
Use the same toroidal direction for the momentum ledger and the rotation
given to `ntv_omega_te_from_rotation`. It returns
\(\Omega_{tE}=d\Phi/d\psi_{\rm pol}\), with electrostatic potential \(\Phi\)
and poloidal flux per radian \(\psi_{\rm pol}\) in the native convention.
The outward field is \(E_r=-\Omega_{tE}\hat{\mathbf e}_r\cdot\nabla\psi_{\rm pol}\).
[Axisymmetric force balance](../src/apps/ntv_device.f90)

`toroidal_momentum_step` in `milca_toroidal_momentum` takes \(N\) shell
inventories, \(N+1\) face fluxes and \(N\) torques in inside-to-outside order.
The face flux and torque averages
are supplied by the caller. Equation (4) integrates exact averages without
time discretization error. Start-of-step inputs give a forward Euler step and require the
stability limit of the chosen flux closure. The FortSym definition supplies
the value, Jacobian-vector and transpose-Jacobian-vector products. The
interval overload encloses the discrete update for enclosed inputs. Errors
from calculating the averages, from the closure and from the geometry enter
through the input enclosures.

## Transport dimensionality

Transport dimensionality is independent of the microphysics model:

| Level | Meaning |
|---|---|
| 0D | global inventory/burn/circuit-like balance |
| 1D | flux-surface or other single-coordinate transport |
| 1.5D | 1D transport with repeatedly solved 2D equilibrium |
| 2D | transport retaining a second physical/angular coordinate |
| 3D | full configuration-space macro transport |
| 5D/6D | direct GK/FK evolution with no separate macro closure |

The 2D coordinates are not hard-coded. They may be \((R,Z)\),
\((\rho,\theta)\), or another mapped pair. Three-dimensional transport must not
assume flux surfaces and is important for stellarators, RMPs, islands,
stochastic fields and edge/divertor geometry.

## Closure-provider contract

The core abstraction is not a scalar diffusion-coefficient getter. A
microphysics provider maps the macro state and geometry to whatever closure
data it can justify:

\[
\mathcal M[U]
\mapsto
\{
\Gamma,\ Q,\ \Pi,\ J,\
L_{ij},\
K(x,x'),\
C_{\alpha\beta}(x,k,\omega),\
D_U\Gamma,\
\epsilon
\}.
\]

A provider need only supply the supported subset. Examples:

- prescribed \(D,\chi,V\);
- neoclassical Onsager matrix;
- quasilinear mode spectrum and fluxes;
- nonlinear local flux-tube fluxes/Jacobians;
- nonlocal orbit/global transport kernel;
- direct global gyrofluid/GK/FK flux fields;
- electromagnetic cross-spectral statistics.

The transport solver consumes a mathematically declared closure, not a
code-specific file format.

## Closure-fidelity ladder

The intended progression is

\[
\text{empirical/prescribed}
\rightarrow
\text{neoclassical/local}
\rightarrow
\text{linear/QL}
\rightarrow
\text{nonlinear local}
\rightarrow
\text{nonlocal orbit/global kernel}
\rightarrow
\text{global nonlinear microphysics}
\rightarrow
\text{fully coupled evolution}.
\]

A local Fick-like law

\[
\Gamma=-D\nabla n + Vn
\]

is the cheapest member, not the universal interface.

A nonlocal model may instead use

\[
\Gamma(x)=\int K(x,x';U)X(x')\,dx',
\]

where \(X\) denotes thermodynamic forces or profile perturbations.

## Traditional integrated modelling as an early mode

MILCA should support the historically successful equilibrium/transport
separation directly. Here "equilibrium" means a selected stationary forward
state; it may be a native Grad--Shafranov/general-MHD solve or, early on, a
validated imported candidate. The transport architecture must not assume that
the state owns nested flux surfaces unless the selected 1D/1.5D representation
explicitly requires them:

\[
\text{equilibrium}
\rightarrow
\text{transport step}
\rightarrow
\text{updated profiles/current}
\rightarrow
\text{new equilibrium}.
\]

Fidelity can then increase continuously:

\[
\text{frozen equilibrium}
\rightarrow
\text{periodically re-solved equilibrium}
\rightarrow
\text{implicitly coupled quasi-static equilibrium+transport}
\rightarrow
\text{dynamic MHD/kinetic evolution}.
\]

This cheap mode is important for scenarios, scans, control and teaching and
should not be delayed until global GK/FK exists.

ASTRA-like 1D/1.5D transport is therefore a valid MILCA configuration, not an
external workflow concept.

The first planned 1.5D instance is the resistive background of a tokamak
core: surface-averaged evolution of q, pressure and density under the
visco-resistive Ohm's law with a Grad-Shafranov constraint, rung L0 of
[the resistive model ladder](nonlinear-resistive-mhd.md#l0-surface-averaged-background),
which states its equations, sources and gates. Its mean-field extension
(rung LM there) is a closure provider under the contract above and is a
hypothesis until checked.

## KAMEL/QL-Balance heritage

The current KAMEL repository already provides a concrete internal example:

\[
\text{profiles}
\rightarrow
\text{KiLCA/KIM linear kinetic response}
\rightarrow
\text{quasilinear transport tensor}
\rightarrow
\text{QL-Balance profile evolution}.
\]

MILCA should first reproduce this pattern in memory through common state and
closure interfaces, preserving QL-Balance/KAMEL as independent regression
oracles.

Then generalize beyond one-dimensional profiles and beyond one response model.

## Stability-to-transport coupling

A linear microstability solve can provide

\[
\{\gamma_k,\omega_k,\xi_k\}
\]

plus a saturation rule and cross-phases to produce transport.

The closure should retain richer information than \(D,\chi\) whenever useful.
At minimum, support complex mode structures and cross-spectral quantities
needed for particle, heat, momentum and electromagnetic transport.

This connects [stability](stability.md) directly to transport and lets the same
linearized operator serve both science and integrated modelling.

## Electromagnetic spectrum/cross-spectrum closure

A turbulence provider may return spectral statistics such as

\[
C_{\alpha\beta}(x,k,\omega)
=
\left\langle
\delta A_\alpha(x,k,\omega)
\delta A_\beta^*(x,k,\omega)
\right\rangle,
\]

with \(A_\alpha\) drawn from the selected fluctuating field/distribution
variables.

Cross-spectra and phases matter. A scalar power spectrum is insufficient for
fluxes such as \(\langle\tilde n\tilde v_E\rangle\).

Possible providers include:

- linear modes + saturation model;
- nonlinear local flux tube;
- sparse samples from a global simulation;
- a global nonlinear gyrofluid/GK/FK state.

The consumer sees one closure contract with declared statistical/temporal
averaging and uncertainty.

## Global orbit transport

Guiding-center/full-orbit transport provides an important intermediate regime
between local 1D coefficients and global 5D/6D turbulence.

A GORILLA-like orbit engine can accumulate transition probabilities, response
kernels or direct particle/energy/momentum transfer between macro cells:

\[
K(x,x';E,\mu,\ldots)
\]

or an equivalent conservative transport operator.

This is especially important for 2D/3D stellarator, RMP, island, stochastic and
finite-orbit-width transport where a radial local coefficient is inadequate.

The public GORILLA code is currently an orbit/characteristic engine; MILCA
should treat richer global-transport/spectrum coupling as a planned extension,
not claim it as existing public GORILLA functionality.

## Turbulence spectrum along global orbits

A particularly valuable multiscale mode is

\[
\text{local/global turbulence statistics}
\rightarrow
\text{orbit-resolved resonant interaction}
\rightarrow
\text{2D/3D macro transport}.
\]

Schematically, for action variables \(J\),

\[
D^{ij}(J)
=
\sum_k\int d\omega\,
\mathcal R^{ij}(J;k,\omega)\,
C(k,\omega),
\]

where \(\mathcal R\) is computed using the actual global orbit/geometry.

This preserves geometry and resonance information that is lost when
microturbulence is collapsed immediately into a radial scalar diffusivity.

The first implementation should start with a simple prescribed spectrum and an
analytic orbit/resonance limit before consuming turbulence generated elsewhere.

## Coupling-strength ladder

MILCA should support progressively tighter coupling:

### C0 — offline/table closure
Precomputed coefficients/tables feed transport.

### C1 — in-memory loose coupling
At macro step \(n\),

\[
U_n\rightarrow \mathcal M[U_n]\rightarrow U_{n+1}.
\]

### C2 — iterated self-consistent closure
Macro state and micro closure are iterated to convergence at each macro step or
steady transport solve.

### C3 — separated-timescale multiscale
Fast local/global microproblems are embedded in slow macro transport, using
spatial/time scale separation.

### C4 — sparse global-micro sampling
Global nonlinear microphysics is sampled only where/when the macro solver or
adjoint error indicator requires it; interpolation/surrogates carry explicit
validity/error status.

### C5 — fully coupled evolution
When scale separation/locality fails, evolve the global nonlinear model
directly.

TRINITY- and Tango-like strategies are conceptual precedents for C3/C4.
GENE-X-like full-f global evolution is a reference for the C5 end of the
spectrum. These are comparison architectures, not code modes.

## Locality and timescale are separate questions

True multiscale modelling must distinguish:

1. whether the fast and slow timescales separate;
2. whether the microphysics is spatially local.

A transport calculation may exploit

\[
\tau_{\rm micro}\ll\tau_{\rm transport}
\]

while still using **global** gyrokinetic/orbit microphysics rather than local
flux tubes.

Thus the hierarchy is two-dimensional:

\[
\text{local}\leftrightarrow\text{nonlocal/global},
\qquad
\text{loose}\leftrightarrow\text{fully coupled}.
\]

Do not hard-wire multiscale to mean local flux tubes.

## Local and temporal adaptive microphysics

Microphysics fidelity need not be fixed for the whole macro domain or macro
step. The transport controller may request richer closure evaluations only on
selected cells/patches, phase-space regions and time intervals.

A microprovider can be activated because

- an adjoint-weighted closure defect is large;
- a surrogate/interpolant leaves its validity domain;
- a local asymptotic ordering fails;
- an instability/resonance is entered;
- stochastic/UQ variance is dominated by a particular region.

The microproblem may be neoclassical DK, GK/FK turbulence, global orbits,
Maxwell--Lorentz particles, Wigner/HF/QMC, detailed radiation or
state-resolved chemistry.

Local timestepping/subcycling/multirate integration is allowed when exchange
conservation and coupling error are controlled. When a transient event ends,
the rich model may be demoted only through a defined projection/reconstruction
with explicit information-loss error.

Cache or surrogate reuse requires a validity key in state/parameter/model
space and an uncertainty envelope. If the error indicator says locality or
timescale separation has failed, enlarge the microdomain or move toward C5
fully coupled evolution.

This is the generalization of C3/C4 from "run a microcode sometimes" to a
space-time adaptive heterogeneous multiscale method. See
[adaptive-multiscale-uq.md](adaptive-multiscale-uq.md).

## Differentiable coupled transport

A closure provider should expose, when available,

\[
D_U\Gamma\,v
\quad\text{and/or}\quad
D_U\Gamma^\dagger w.
\]

Then the macro transport balance can be solved implicitly:

\[
R_{\rm tr}(U,\mathcal M[U])=0
\]

rather than by fragile explicit profile iteration.

For cheap symbolic closures these derivatives are generated directly. For
expensive nonlinear microphysics they may come from tangent/adjoint solves,
validated secants or a surrogate with an explicit error/validity envelope.

Adjoints also determine which closure evaluations matter most to a requested
observable and therefore where expensive microphysics should be recomputed.

## Verification and conservation

Every transport model declares conservation of:

- particles/nuclei;
- charge;
- momentum;
- energy;
- current/flux where relevant.

Internal micro-to-macro transfers cancel in the total ledger. Sources, sinks,
radiation and walls are explicit ports.

Verification starts with exact diffusion/advection/source problems in 1D,
then separable 2D/3D manufactured cases, Onsager/neoclassical limits, and
cross-code integrated cases.

The first 2D/3D transport implementation should use the same FortFEM mapped
spaces as the rest of MILCA rather than introduce an independent transport mesh
stack.

## Certification

Transport certificates separate:

- macro discretization/time error;
- closure-provider numerical error;
- closure/model-reduction error;
- interpolation/surrogate error;
- statistical error for sampled microphysics;
- equilibrium/geometry uncertainty.

A quasilinear coefficient with a tiny numerical error is not a rigorous
turbulence-model error bound.

Goal-oriented adjoints should eventually drive where and **when** to refine the
macro grid/time step, recompute microphysics, increase spectrum resolution,
promote correlation/quantum/radiation fidelity or replace a cheap closure by a
richer one. Multifidelity/MLMC/control-variate UQ should spend many samples on
cheap closures and sparse samples on expensive corrections while retaining
separate deterministic model-error status.

Chris&AI
