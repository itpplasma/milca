# Kinetic layer response to resonant perturbations (KL)

Status: plan and first stage, 2026-11-15. The linear response of a magnetised plasma to a
static helical perturbation near its resonant surface, with the ion Larmor radius treated by
integral kernels instead of an expansion, first in a straight periodic cylinder and then in
the torus. The cylinder model is that of Markl, Kasilov and Albert (J. Plasma Phys., accepted
2026; the code KIM in KAMEL), which also outlines the toroidal extension.

## Stages

1. **KL0, susceptibility functions** of the Ornstein-Uhlenbeck collision model (this page,
   `milca_ou_susceptibility`).
2. **KL1, cylinder response**: charge and current kernels in radial wavenumber with the
   gyroaverage, Poisson solve, shielding current. Compared with KIM (integral and second-order
   Larmor-radius runs) on a device case and an ITER q = 3 surface.
3. **KL2, toroidal response**: ions in canonical action-angle variables with the effective
   toroidal radius psi* (resonances of bounce, transit and precession frequencies), electrons
   from the drift-kinetic operator with the misalignment-field source, the finite-orbit kernel
   of [drift-kinetics-extensions.md](drift-kinetics-extensions.md) as independent ion check.
   Limits: KL1 in the cylinder, NEO-RT for thin orbits, POTATO at finite orbit width,
   NEO-2-QL for the electrons.

## KL0: susceptibility functions

In units of the thermal speed v_T and the collision frequency nu, with x1 = k_par v_T/nu and
x2 = (omega - omega_E)/nu, the velocity-space Green's function of
(-i x2 + i x1 u - L) G = delta(u - u'), L g = (g' + u g)', gives

    I^{kl}(x1, x2) = int du du' u^k G(u, u') u'^l exp(-u'^2/2)/sqrt(2 pi).

L generates the stationary Ornstein-Uhlenbeck process (unit variance, correlation
exp(-|t|)), so by Feynman-Kac I^{kl} = int_0^inf dt exp(i x2 t) E[u(t)^k u(0)^l
exp(-i x1 X)] with X = int_0^t u ds. The process is Gaussian: Var X = 2 (t - 1 + exp(-t)),
and the tilt exp(-i x1 X) shifts u(0) and u(t) by mu = -i x1 (1 - exp(-t)) at unchanged
correlation. Hence

    I^{kl} = int_0^inf dt exp(i x2 t - x1^2 (t - 1 + exp(-t))) E[(b + mu)^k (a + mu)^l]

with (a, b) a standard normal pair of correlation exp(-t). Substituting s = exp(-t), each power
s^j of the moment polynomial gives a Kummer integral M(1, p + j + 1, x1^2)/(p + j),
p = x1^2 - i x2. This sum is used for |x1| < 3; above, where its terms cancel like x1^(k + l),
the time integral is evaluated directly on adaptive Gauss-Legendre panels with the moments
computed from mu (no cancellation). The singular point is x1 = x2 = 0 (no damping).

The integral part of the collision operator, nu phi(u) (u^2 - 1) <u^2 - 1, f> for energy and
nu phi(u) u <u, f> for momentum, is a rank-one source; the conserving functions follow from the
Woodbury identity on I (`ou_conserving`). Particles are conserved by L itself.

Checks (`ou_susceptibility`):
- an independent solution of the same equation in Hermite functions, where L is diagonal and u
  tridiagonal, for all four conservation models at x1 from 0.05 to 12 and x2 from -5 to 400,
  within 2e-12 times the condition number of the conserving correction;
- the collisionless limit -i Z(x2/(sqrt(2) x1))/(sqrt(2) x1) of I^{00} at x1 = 3000;
- the collisional limit x1 = 0 in closed form, and symmetry I^{kl} = I^{lk};
- the closed forms of KiLCA (calc_Imn_array). On a grid of 24 points they agree to 1e-12 or
  better except at small x1 with large x2, where the KiLCA forms lose up to 6e-9 and the
  Hermite solution agrees with this module to 2e-12.

## KL1: cylinder response (`milca_kinetic_layer_cylinder`)

Gaussian units. The background is given on one period of a radially periodised
window (equidistant, endpoint excluded): k_par, k_s, omega_E, B0, E_r and per
species n, T, n', T' and the collision frequency. For one harmonic the kernels
in the Fourier basis exp(i k_l r), |l| <= M, are

    K(l, l') = (2 pi/N) sum_j exp(-i (k_l - k_l') r_j)/(8 pi^2) sum_species core,

with the gyroaverage harmonic only. With b+ = rho_T^2 (k_perp^2 + k_perp'^2)/2,
bx = rho_T^2 k_perp k_perp', S_n = exp(-b+) I_n(bx) and
F(I^a, I^b) = I^a (S0 (A1 + A2 (1 - b+)) + A2 bx S1) + A2 S0 I^b/2:

| kernel | core |
|---|---|
| rho-Phi | (1/lambda_D^2) (-G_B + i v_T^2 k_s/(omega_c nu) F(I00, I02)) |
| rho-B_r | -(1/lambda_D^2) v_T^3/(omega_c nu c) F(I01, I03) |
| j_par-Phi | (1/lambda_D^2) i v_T^3 k_s/(omega_c nu) F(I10, I12) |
| j_par-B_r | -(1/lambda_D^2) v_T^4/(omega_c nu c) F(I11, I13) |

The susceptibility indices follow from moment counting: the output moment
supplies k, the source (u_par for B_r and the u_par^2 of the energy factor of the
A2 drive) supplies l. KIM uses I^{21} in place of I^{03} for rho-B_r; the two
are equal for the energy-conserving operator (checked numerically), not for
the bare Ornstein-Uhlenbeck operator. G_B = exp(-rho_T^2 (k - k')^2/2) is the
gyroaverage of the Boltzmann term at the particle position (Eq. 2.41 of the
paper); KIM sets it to 1 (`boltzmann_gyroaverage = .false.`). Poisson's
equation (4 pi K^rhoPhi - diag(k_l^2)) Phi = -4 pi B_r K^rhoB(:, 0) is solved for
a radially constant B_r, and j_par = K^jPhi Phi + B_r K^jB(:, 0).

Checks (`kinetic_layer_cylinder`): the overflow-free Bessel products against the
direct ones; the local rho-Phi core against the second-order Larmor-radius form
(Eq. 2.66) with an O(b^2) remainder; diagonal kernels on a uniform background;
and the ideal limit Phi = -i E_r B_r/(k_par B0) far from the resonance, whose
remaining deviation is the ion Larmor-radius correction (rho_i k_s)^2 (3.4e-3 at
k_s = 0.12/cm, 3.4e-5 at 0.012/cm).

`kl_cylinder` solves one window from a text background. Given the window background
that KIM writes for an a tokamak pedestal case (harmonic (6, 2), 512 points,
M = 225, inputs not public), Phi, the electron and ion parallel currents and the
integrated shielding current agree with KIM's local periodic branch to 4e-5
within 1.5 cm of the resonant surface; the Boltzmann gyroaverage that KIM omits
changes the integrated current by 1.6e-4 there.

`milca_kinetic_layer_window` builds the window from radial profiles (n_e, T_e,
T_i, q, E_r against the effective radius): cubic splines of the profiles, the
cylinder force balance du/dr = -2 r u/(q^2 R0^2 + r^2) - 8 pi dp/dr for
u = B0^2 from the vacuum value at the first radius, the geometry k_par, k_s,
omega_E and B0 of that true equilibrium, and the periodic continuation of
profiles and geometry as functions (`fortnum_periodize`), so that the core of
the window is the true background; integrating the force balance through the
blended profiles would add the artificial pressure gradient of the
transition bands to B0. Collision frequencies follow KIM after KAMEL#399.
`kl_profiles` runs the whole chain from profile files. Checks
(`kinetic_layer_window`): closed-form B0 and wavenumbers for constant q
without pressure, the given profiles on the core, continuity at the window
edge and the diamagnetic drop 8 pi Delta p of B0^2. On the a tokamak case
the window matches KIM's to 1e-8 in the geometry and 1e-11 in the profiles
and collision frequencies, and Phi, j_par and the integrated shielding
current from the profiles agree with KIM to 2e-5. Next: the kinetic layer
closure of the per-n torus response and an ITER case.

## KL closure: the kinetic layer in the per-n torus response (plan)

The per-n FEEC response (`mhd_ampere_n`, `mhd_response_n`) represents each
rational surface by two edge sheets whose amplitudes follow from linear
constraints (`layer_closure_t`); the resistive closure L5 sets the net layer
current to Delta_layer times the reconnected flux A (`resistive_closure`).
The kinetic closure uses the same constraints with Delta from the cylinder
response of the rational surface, built from the flux-surface profiles mapped
to the effective radius (KL1). In MILCA's conventions a sheet of amplitude
lambda carries the surface current K = -i m lambda/r_s along the toroidal
direction, and the resonant flux of B^r_m is Phi_m = 2 R0 r_s B^r_m (peak
amplitude). With div B = 0 and the jump of B_theta across the sheet,

    Delta_r = [dPhi_m/dr]/Phi_m = i m mu0 K/(r_s B^r_m),
    Delta_psin = Delta_r dr/dpsin,

where K = I_par/(2 pi r_s) is the shielding current per unit length of the
cylinder response for the resonant B^r_m (constant psi), converted to SI.
The orientation of the toroidal field and of the poloidal harmonic enter the
sign of Delta. The ideal limit (|Delta| large) does not fix it; the energy
balance does: the electromagnetic torque on the layer, proportional to
Im(Delta) |A|^2, must brake the plasma rotation relative to the static field,
with the sign that L5's resistive-inertial layer has at the same rotation.
Checks before landing, in this order:

1. Cylinder, end to end (done): KIM's electromagnetic run type solves
   Poisson's and Ampere's laws together with the kinetic current for a
   prescribed B^r at the outer boundary. On a device case, (m, n) =
   (6, 2), the jump of d(r B^r)/dr across the layer in that run equals
   -2 i k_s I_par/(c r_s B^r) from its own current to 0.3 %, which fixes sign
   and normalisation of Delta without orientation conventions (it required
   the fix of KIM's Ampere operator, itpplasma/KAMEL#414). The constant-psi
   Delta from the electrostatic current (B^r uniform across the layer) agrees
   with the electromagnetic one to 2 % in the real and 18 % in the imaginary
   part. This gap is not the variation of psi across the strongly shielded
   layer (B^r(r_s) is 4 % of its edge value): MILCA's kernels driven by the
   electromagnetic run's own B^r(r) give a net current within 2 % of the
   constant-psi one, but 14 % below the current of that run in magnitude.
   The electromagnetic run uses KIM's non-periodic kernels, the periodic
   electrostatic branch (which MILCA reproduces to 2e-5) its periodic ones;
   the difference between the two KIM branches on this case is open.
   `milca_kinetic_layer_delta` (`kl_delta`) returns Delta and I_par per gauss
   for the profiles of one rational surface.
2. Cylinder to torus (implemented as `closure = 'kinetic'` of
   `gs_response_3d`): the closure itself is the L5 matching with a
   prescribed complex Delta, which the L5 suite checks against the cylinder
   boundary-value problem at large aspect ratio. New are the orientation and
   the radial scale of Delta:
   - MILCA expands a per-n field as e^{i (n zeta - m theta_s)} with theta_s
     counter-clockwise in (R, Z) and (R, phi, Z) right handed
     (`milca_vacuum_surfaces`). The cylinder (r, theta, z) of the kinetic
     layer is right handed, so with theta = theta_s its axis points along
     -phi: the harmonic is (-m, -n) in the cylinder's e^{i (m theta +
     n z/R0)}, q = -q_s and B_z = -F/R. Both descriptions then have the
     resonance at m = n q_s and the same phase, so Delta carries over
     without conjugation.
   - The rational sits where |q| of the kinetic profiles equals |q_s|, and
     Delta per psin is Delta_r |dq/dpsin|/|dq/dr|, so the local shear of
     the two descriptions matches and no definition of the effective radius
     enters.
   - Sign of the time convention: in a flat collisional plasma, Im Delta of
     the kinetic layer at the Doppler frequency omega_E has the sign of the
     resistive-inertial layer of L5 with gamma = -i (0 - n Omega), n Omega =
     omega_E (`test_kinetic_layer_delta`), so KIM's and MILCA's conventions
     agree.
   On the A = 10 test torus the solution meets the matching condition to
   1e-14 (`gs_response_3d_kinetic_closure`).
3. Ideal limit (shielding grows with |Delta|) and the comparison with
   MEPHIT's kinetic hybrid on a device case; the ITER demonstration case.

## Layers matched to a toroidal outer solution (KIM-1L)

`milca_kinetic_layer_match` (`klm_solve`) closes N layers with an outer
region given by its tearing matrix E (dimensionless) and coil drive C
[G cm], psi = r B^r at the rational surfaces,

    sum_k (r_j Delta_j delta_jk - E_jk) psi_k = C_j,
    T_j = -pi R0 r_j n Re(sigma_j) |psi_j/r_j|^2/(m_j c),

with Delta_j and sigma_j = I_par/B^r of `kl_delta` at r_j. E and C come from
a toroidal ideal outer solution (GPEC/DCON or the per-n response of
`gs_response_3d`). Checks (`test_kinetic_layer_match`): a vacuum cylinder
(E = -2|m|) returns the single-surface torque of `klt_torque`, with a
difference falling as R0^-2 because `klt_torque` approximates the binormal
wavenumber by m s/r_s; |Delta| -> infinity removes flux and torque;
Delta = 0 passes psi = -E^{-1} C, and a purely reactive layer carries no
torque.
