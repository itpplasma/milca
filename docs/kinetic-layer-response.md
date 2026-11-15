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
   Larmor-radius runs) on ASDEX Upgrade 33353 and an ITER q = 3 surface.
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
k_s = 0.12/cm, 3.4e-5 at 0.012/cm). Next: the periodised background from
profiles and the comparison with KIM-1L on ASDEX Upgrade 33353.

