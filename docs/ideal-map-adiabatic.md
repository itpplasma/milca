# Adiabatic pressure in the map lane

With an adiabatic index gamma > 0 the mass per flux surface is fixed instead of
the pressure. The equilibrium is a stationary point of the Kruskal-Kulsrud
energy over maps X with frozen flux,

    W[X] = int |B|^2/(2 mu0) dV + int_0^1 K(s) f(V'(s)) ds,             (1)
    f(v) = (v^(1-gamma) - 1)/(gamma - 1),   f'(v) = -v^(-gamma),        (2)
    p(s) = K(s)/V'(s)^gamma,                                            (3)

with s the normalised toroidal flux of the surface, V' = dV/ds over the full
torus and K the mass function. At gamma = 1, f(v) = -ln v (isothermal).

## Discrete V'

Every reference point y carries its surface label s(y), and the coarea formula
gives, for any test function phi(s),

    int phi(s) V'(s) ds = int phi(s(y)) J dy dphi-hat,                   (4)

J the map Jacobian. V' is the L2 projection on shifted Legendre polynomials
P_k(2s - 1), k <= N (`adiabatic_degree`), with moments per field period

    b_k = int P_k(s(y)) J dy dphi-hat,                                   (5)
    V'_h(s) = nfp sum_k (2k + 1) b_k P_k(2s - 1),                        (6)

so V'_h is a smooth function of the map coefficients and (5) uses the same
quadrature as the magnetic energy. The s integral in (1) uses 2(N + 1) + 16
Gauss-Legendre points.

## Gradient and Hessian

With v = V'_h/nfp, p_h = K/(nfp v)^gamma and the projected pressure
P_h = sum pi_k P_k, pi_k = (2k + 1) int p_h P_k ds,

    dW/dX = dW_frozen/dX at p = P_h,                                     (7)
    d2W/dX2 = H_frozen(P_h) + G Q G^T,                                   (8)
    G = (db_0/dX, ..., db_N/dX),
    Q_kl = (2k + 1)(2l + 1) int gamma p_h/v P_k P_l ds.                  (9)

(7) holds because the pressure term of the frozen energy is -int P_h J =
-pi . b, which cancels the dependence of the thermal energy on b. Q is
positive definite while V'_h > 0, so compression stiffens every mode. The
volume gradient dJ/dX is the generated energy leaf at b = 0, p = -1. The
generator proves that the density is affine in p with slope -J and quadratic
in b (`gen/gen_ideal_map.f90`).

The Newton step and the current-matching sensitivity factor the sparse
frozen part A = H_frozen + shift. The rank N + 1 term enters through Woodbury,

    (A + U Q U^T)^-1 r = y - Z (I + Q U^T Z)^-1 Q U^T y,
    y = A^-1 r,  Z = A^-1 U,  U = R^T G,                                  (10)

with N + 1 extra solves per factorisation. The nested omega elimination, the
Krylov step and the lagged quasi-Newton step reject gamma > 0.

## Verification

`test/test_ideal_map_adiabatic.f90` on a sheared 3D seed with gamma = 5/3,
N = 6:

| check | result |
|---|---|
| gradient against energy differences, seed and after two steps | 1.1e-9, 1.2e-8 |
| Hessian (8) against gradient differences | 2.0e-8, 2.4e-8 |
| P_h against p* on a fixed-pressure equilibrium with K = p* V'^gamma | 6.5e-16 |
| adiabatic gradient there against the fixed-pressure one | equal (2.5e-12) |
| Newton (10) from the seed back to that state | 9 iterations, energy 2.9e-16 |

## VMEC normalisation

VMEC holds mass(s) = mu0 PRES_SCALE AM(s) (|Phi'| R00)^gamma, Phi' = phiedge/(2 pi), and
p = mass/vp^gamma with vp = V'/(4 pi^2). In the full-torus V' of (3) this is
K(s) = mu0 P(s) (2 pi |phiedge| R00)^gamma, which `nested_equilibrium` uses for
its `pressure` profile when `gamma > 0`.
