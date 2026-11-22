# Exact-force-balance families and Boozer/Hamada coordinate consistency

**Status (2026-11-22): exploratory research plan; no implementation or new numerical results.** This is not a replacement for the [equilibrium roadmap](equilibrium-roadmap.md) or its P0--P3 Cartesian-flux-map versus Lagrangian-map comparison. [ROADMAP.md](../ROADMAP.md) is the only active priority list. Do not change the ideal-equilibrium core or build a generic neural/GP framework before a small falsifiable gate has passed.

## Question and exactness contract

Can one generate an expressive family of static isotropic-MHD equilibria whose *interior differential equations* hold identically, and spend approximation error on matching **boundary geometry, pressure/current/flux profiles, coil currents and diagnostic signals**?

The continuum model is

$$
\nabla\cdot B=0,\qquad
(\nabla\times B)\times B/\mu_0=\nabla p.
$$

Different meanings of "exact" must be kept separate:

1. **Analytically exact family:** identities follow for all admissible parameters, independent of the fit. A symbolic identity, regularity/domain hypotheses, and an independent evaluation are required. Floating-point evaluation still has rounding error.
2. **Algebraically constrained family:** basis coefficients satisfy exact polynomial/Fourier or other operator identities only if the full coefficient constraints are solved and, for a mathematical claim, certified. A small algebraic residual is not an identity.
3. **Nonlinear manifold corrector:** a discrete equilibrium residual is driven to solver precision. This is not automatically pointwise continuum force balance between grid points; report strong/weak residual and discretisation/certification status independently.

Inputs must be classified **hard** or **soft**. Pressure, current, rotational transform, edge flux, magnetic topology, conductor model, coil currents and moving boundaries must never change silently. Failure to fit incompatible data is an informative physical result, not a reason to weaken the declared equations.

The basis is not restricted to GPs: analytical functions, polynomial/Trefftz, Fourier, B-spline, rational, RBF/kernel, GP, or neural representations are admissible *if* the chosen field construction preserves the actual differential identities. Monge/physics-consistent GPs parametrise kernels of **linear** differential operators; arbitrary superpositions of nonlinear MHD equilibria are not equilibria.

## Track E0: exact axisymmetric Trefftz family (first low-risk gate)

For $R>0$, use the restricted Grad--Shafranov profiles

$$
p'(\psi)=a,\quad F F'(\psi)=b,\quad
\Delta^*\psi=-\mu_0 a R^2-b,\quad
\Delta^*=R\partial_R(R^{-1}\partial_R)+\partial_Z^2.
$$

An exact family is

$$
\psi=-\mu_0 a R^4/8-b Z^2/2+\sum_k c_k H_k(R,Z),
\qquad \Delta^*H_k=0,
$$

with $p=p_0+a\psi$ and $F^2=F_0^2+2b\psi$, subject to $F^2\geq0$, domain regularity and chosen topology. Optimize the coefficients against a *target* boundary ($\psi=\psi_b$), profiles and synthetic diagnostics, keeping the exact interior GS equation.

**Acceptance E0:** fit an ellipse, a rounded triangular target, and an intentionally incompatible boundary/profile target; verify the complete continuous identities symbolically and independently at off-grid points to roundoff; report separate boundary, profile, coil and plasma--vacuum stress mismatches, condition estimates and time-to-fit. Enrich the Trefftz family and verify that interior exactness is unchanged even when data fit is poor. An exact *interior* solution alone does not make a coil-driven free-boundary equilibrium.

## Track E1: exact three-dimensional finite-pressure families

Use independently specified published analytic nonsymmetric nested equilibria (e.g. Landreman's 2026 families cited in [equilibria-inference.md](equilibria-inference.md)). Identify the family's actual free parameters, flux/pressure restrictions and regularity; do not claim arbitrary boundary or profile completeness. Fit perturbed target surfaces, profiles and synthetic signals *within* that exact solution family. Compare with both numerical MILCA map directions as independent references.

**Acceptance E1:** one analytic 3D finite-$p'$ state, independently checked in physical space for $B$, $J$, $p$, $\nabla p$, $\nabla\cdot B$, $J\times B-\nabla p$, flux and $\iota$; then a deliberate target mismatch with a measured approximation floor. An analytic case with $p'=0$ alone is insufficient for the Hamada-current tests below.

## Track E2: stress-potential construction (high risk, not yet a method)

The symmetric MHD stress

$$
T=\left(p+\frac{B^2}{2\mu_0}\right)\mathbf 1-\frac{B\otimes B}{\mu_0}
$$

satisfies $\nabla\cdot T=0$ exactly when $\nabla\cdot B=0$ and MHD force balance holds. A local Beltrami double-curl stress potential on a suitable contractible domain can impose the **linear** differential identity:

$$
T_{ij}=\epsilon_{ikl}\epsilon_{jmn}\partial_k\partial_m H_{ln},
\qquad \partial_j T_{ij}=0,
\qquad B=\nabla\times A.
$$

Crucial missing constraint: almost every divergence-free symmetric tensor is **not** an MHD stress. Enforce the *entire* nonlinear algebraic identity linking $T$, $B$ and $p$, positivity of the magnetic eigenvalue gap, smooth global orientation, divergence freedom and flux/topology constraints. A double stress eigenvalue alone is not sufficient; global potential completeness is topology dependent. Coefficient matching may make the differential part structural, but the nonlinear algebra may be as hard as the original PDE.

**Acceptance E2:** represent and independently verify a known finite-pressure solution by exact polynomial/coefficient identities, reject non-MHD stresses, and only then attempt a genuinely non-axisymmetric finite-pressure family. Report solution-space dimension, conditioning, algebraic solve cost and any unsatisfied identities. Do not claim general pointwise exactness based on a truncated collocation solve.

## Track E3: nonlinear equilibrium-manifold corrector (practical alternative)

Let a generator $U_0(a)$ use any useful basis, including PCGP/kernel or learned latent coordinates, and allow explicitly selected boundary/profile/coil parameters $\vartheta$ to move:

$$
R_h(U,\vartheta)=0,\qquad
\min_{U,a,\vartheta}\;
\|U-U_0(a)\|^2_W
+\|\mathcal O(U,\vartheta)-d\|^2_{C_d^{-1}}
+\|\vartheta-\vartheta_{\rm target}\|^2_{C_\vartheta^{-1}} .
$$

The equality means **discrete force balance** unless an exact analytical function family is supplied. Tangents satisfy

$$
R_U\delta U+R_\vartheta\delta\vartheta=0.
$$

Differentiate the converged equations implicitly (JVP/VJP), not Newton iterations. A linearized physics-consistent kernel can generate tangent directions, followed by nonlinear correction; its linear nullspace alone is not an exact nonlinear manifold.

**Acceptance E3:** first exact Solov'ev, then analytic 3D, then coil-driven free-boundary fixtures. Compare the baseline method at matched *independently evaluated physical residual*, input mismatch, convergence basin, nonlinear/linear cost and sensitivities against re-solved finite differences. Preserve branches and distinguish solver residual from continuum accuracy. Use the selected P3 map or current proven Lagrangian map; do not introduce a third production core.

## Track C: "perfect" Boozer/Hamada via independent compatible routes

For regular nested isotropic-MHD equilibria, independent constructions of a **given** straight-field-line coordinate system should agree up to normalized angular/flux gauges in the resolved continuum limit:

- **Boozer covariant:** tangential covariant $B_\theta=I(\psi)$ and $B_\zeta=G(\psi)$.
- **Boozer Jacobian:** $\sqrt{g_B}B^2=G+\iota I$ in a consistent flux convention.
- **Hamada I:** $\sqrt{g_H}$ is a flux function.
- **Hamada II:** electric-current lines are straight.

**Boozer and Hamada coordinates are not generally the same coordinates.** Equality of the **two Boozer routes** probes current tangency, and equality of the **two Hamada routes** probes corresponding current/force-balance compatibility, but neither proves *all* components of $J\times B=\nabla p$. Always retain an independently computed vector force-balance residual and averaged radial balance.

Use [libneo's existing differentiable routes](https://github.com/itpplasma/libneo/blob/47b1988cddf5e6544d92eb7418c87f83be3d7ea4/src/boozer_hybrid.f90) and [quality diagnostics](https://github.com/itpplasma/libneo/blob/47b1988cddf5e6544d92eb7418c87f83be3d7ea4/src/boozer_hybrid_quality.f90); do not duplicate them in MILCA. Match $\psi$, $\iota$, field period, handedness/sign, Fourier range, angle gauge, interpolation and radial position. Compare coordinate maps/round trips, transformed fields, spectra and Jacobians, not only $B_{00}$.

**Non-negotiable exceptions:**
- Magnetic differential equations have $m\iota-nN_{\rm fp}$ denominators. At rational surfaces, solvability conditions are necessary; near rationals, small divisors amplify errors. Report the resonant forcing and excluded modes. A numerically zeroed mode must not be claimed as a perfect coordinate solution.
- libneo Hamada II has a denominator proportional to $p'$ under radial equilibrium. At $p'=0$ it is underdetermined/singular even for a perfect vacuum or force-free state: use Hamada I and mark Hamada-II comparison unavailable.
- Exact continuum force balance **does not guarantee machine-precision numerical coordinate transforms** at finite Fourier/truncation resolution. Convergence, regularity and conditioning must be reported.
- Current tangency $J\cdot\nabla\psi=0$ or route agreement is **necessary but not sufficient** for general isotropic force balance.

**Gate C0:** exact axisymmetric Solov'ev with $p'\ne0$: covariant versus Jacobian Boozer and Hamada I versus current-based Hamada II, gauge-normalized angular differences and spectra, and direct force balance. Add a $p'=0$ case to verify Hamada-II failure is classified as a degeneracy, not an equilibrium defect. Record libneo's radial balance, $J^s$, lambda consistency, Boozer family mismatch and resonant forcing diagnostics. Check interpolation/refinement towards roundoff where well-conditioned.

**Gate C1:** repeat on a known analytic nonsymmetric 3D finite-pressure state and on a divergence-free, surface-tangent perturbation which violates force balance. Treat route agreement as an additional **independent consistency check**, not a substitute for force balance. Test rational and near-rational surfaces separately (including an honest unavailable status).

**Gate C2 (optional):** if C0/C1 show useful discrimination, add Boozer/Hamada *route mismatch* as a supplementary fit observable for E0/E1/E3. Never silently adjust profiles to obtain agreement or claim that zero route difference establishes existence of an ideal equilibrium.

## Common result ledger and dependencies

| Ledger | Report for every proposed gate |
|---|---|
| Physical identity | Exact symbolic/coefficient identity status separately from independent off-grid $\|\nabla\cdot B\|$ and $\|J\times B-\nabla p\|$; units, normalization, conditioning, rounding/certification status |
| Fit to inputs | Boundary displacement and normal-field defect, $p/I/F/\iota$ profile mismatch, coil-current deviation, plasma--vacuum traction/trace error and diagnostic misfit; hard versus soft |
| Topology | Valid domain, orientation/nonfolding, flux/winding, axis, surfaces/islands and excluded stochastic/separatrix regions, interface current sheets if admitted |
| Coordinates | Covariant/Jacobian Boozer discrepancy, Hamada I/II where valid, $J^s$, radial balance, resonance/zero-mode diagnostics, spectral/coordinate round trip |
| Algorithm | Basis size, algebraic rank/condition, steps, wall time, memory, parameter sensitivities and comparison at matched *physical* error |
| Evidence | candidate, symbolic exact family, certified algebraic/discrete/continuous, or unavailable/failed. Distinguish numerical, fit/input, model, and statistical uncertainty |

**Order:** preserve P0--P3. Then small offline E0 and C0 (reusing analytic fixtures), E1/C1 if supported; E2 only after exact coefficient construction is justified; E3 when a verified equilibrium residual is available. These tracks are research alternatives, not active MILCA implementation instructions until promoted in ROADMAP.md. A result is useful if it reduces *time to verified force-balanced fit* at a specified boundary/observable tolerance, not merely because it moves residuals into another variable.

**Ownership:** MILCA generates/solves equilibria and publishes physical observables and error status; FortSym/FortFEM/FortNum/FortSparse own reusable symbolic/numerical primitives; libneo owns SFL coordinate construction and JVP/VJP; TIAGO owns diagnostic likelihoods, priors and reconstruction; stellarator optimizers consume the differentiable provider. No new general-purpose GP/ML framework is a prerequisite.

## Source ledger and open literature checks

This is a plan derived from equations and existing group repositories; **no tests were run here**. Before implementation, verify the precise hypotheses of each published construction, source licences and independent derivations.

- Solov'ev (1968) exact restricted Grad--Shafranov family and MILCA's [existing GS reference](grad-shafranov.md).
- Landreman, analytic nonsymmetric equilibria, arXiv:2609.26742 (2026), recorded in [equilibria-inference.md](equilibria-inference.md).
- Classical Beltrami elasticity stress potentials: review the local/global tensor-potential theorem and MHD algebraic compatibility before deriving software from it.
- Moser, Albert and Ranftl, *Parameter learning with physics-consistent Gaussian Processes*, and group-owned [PCGP](https://github.com/krystophny/PCGP) as conceptual inspiration only; no code copied.
- Group-owned libneo [Boozer/Hamada constructions](https://github.com/itpplasma/libneo/blob/47b1988cddf5e6544d92eb7418c87f83be3d7ea4/src/boozer_hybrid_routes.f90) and [independent tests](https://github.com/itpplasma/libneo/blob/47b1988cddf5e6544d92eb7418c87f83be3d7ea4/test/source/test_sfl_routes.f90), snapshot 47b1988cddf5e6544d92eb7418c87f83be3d7ea4; no code copied.
