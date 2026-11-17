# SFINCS reference for SI surface moments

`surface_fluxes_sfincs.txt` contains two synthetic transport matrices used as
inputs to `surface_fluxes` and independent SFINCS moments as expected outputs.
It tests the SI conversion, including its factor of two, geometric factors,
thermal force and parallel-field normalization.

The SFINCS source is the public BSD 2-Clause
[Landreman repository](https://github.com/landreman/sfincs/tree/8df5453472e982df0f6ae005243ce38d57a83711).
The reference uses PETSc 3.18 and MUMPS. The full matrix is factored with
`useIterativeLinearSolver=.false.`. PETSc GMRES with that direct
preconditioner reports the true residual. The companion
[`full`](surface_fluxes_sfincs_full.namelist) and
[`Lorentz`](surface_fluxes_sfincs_lorentz.namelist) namelists retain all
physical and discretization settings.

## Field and plasma

Use Boozer angles \(\theta,\zeta\), field scale \(B_0\), and
\(b=B/B_0\). Both codes use

\[
b=1-0.12\cos\theta-0.04\cos(\theta-3\zeta)
+0.02\sin(2\theta-3\zeta).
\]

The number of field periods is 3, rotational transform \(\iota=1.25\),
safety factor \(q=0.8\), covariant field components \(G=6\) and \(I=0.3\)
in tesla metres, and parallel coefficient
\(\kappa=1.25/6.375\) inverse metres. The common incompressible electric
rate is 4000 inverse seconds. Neither code includes electric speed/pitch
terms, tangential magnetic drift or a varying surface potential here.

The species is deuterium, with mass
\(3.3435837768\times10^{-27}\) kilograms, charge number 1, density
\(3\times10^{19}\) inverse cubic metres, and temperature 2 kiloelectronvolts.
MILCA uses Coulomb logarithm 136.04778278483653 and
`collop_make_grid(6,12,6,...)`, with 11 poloidal and 7 toroidal points and
24 Legendre coefficients. The matrix inputs were calculated with
`surface_transport` and these parameters. The field values are sampled on
its uniform angular grid.

SFINCS uses reference length 1 metre, magnetic field 1 tesla, the above
density, mass and temperature, and reference potential 2000 volts. Its
collision parameter is `nu_n=0.005`, with `Delta=0.0091365219824869892`
and `alpha=1`. Constraint scheme 2 uses a density source at each speed.
`Nxi_for_x_option=0` keeps 24 Legendre coefficients at every speed,
matching MILCA's angular-pitch discretization.

## Fixture format and units

The first line is the number of cases. Each case contains:

1. MILCA collision model, then \(\langle b^2\rangle\).
2. Nine MILCA matrix entries in Fortran order, first matrix index fastest.
3. Nine expected SI moments in the same order, with rows particle flux,
   kinetic-energy flux, and \(\langle b u_\parallel\rangle\).

The average uses weight \(b^{-2}\). If \(a=10^{-6}\) is the fixture's
dimensionless drive multiplier, the three tests use

| Drive | \(n'/n\) [1/Wb] | \(T'/T\) [1/Wb] | \(E_0\) [V/m] |
|---|---|---|---|
| Density | \(a\) | 0 | 0 |
| Temperature | \(3a/2\) | \(a\) | 0 |
| Parallel electric | 0 | 0 | \(2000a/\langle b^2\rangle\) |

Primes mean derivatives in toroidal flux divided by \(2\pi\).
The force input \(\Phi'\) is zero in these basis tests. SFINCS
`RHSMode=2` similarly separates the three force columns from electric
advection in its operator. These are matrix-basis responses, rather than a
self-consistent plasma state.

Let \(n_0,T_0,v_{t0}\) be the reference density, temperature in joules and
thermal speed. The reference length is \(L_0=1\,\mathrm m\) and the
reference flux is \(\psi_0=B_0L_0^2=1\,\mathrm{Wb}\).
From SFINCS HDF5 diagnostics the expected rows are

\[
\Gamma=a n_0v_{t0}\frac{\psi_0}{L_0}\,\mathtt{particleFlux\_vm\_psiHat},\quad
Q=2a n_0T_0v_{t0}\frac{\psi_0}{L_0}\,\mathtt{heatFlux\_vm\_psiHat},\quad
U=a v_{t0}\,\mathtt{FSABFlow}.
\]

The factor of two in the energy conversion follows SFINCS's reference
energy \(m_0v_{t0}^2=2T_0\). These factors were derived from its manual and
`diagnostics.F90` before comparing coefficients.

## Refinement and acceptance

| Check at fixed angular-pitch grid | Lorentz | Full self-collisions |
|---|---|---|
| SFINCS speed refinement | 24 to 32 points, maximum change 0.0629% | 12 to 16 points, maximum change 0.0325% |
| MILCA speed refinement | 8 to 12 elements, maximum change 0.0144% | 8 to 12 elements, maximum change 0.00448% |
| Maximum final true SFINCS relative residual | \(3.35\times10^{-8}\) | \(2.07\times10^{-11}\) |
| Maximum SI moment discrepancy | 0.0160% | 0.0370% |

The benchmark protocol requires speed-refinement changes below 0.3%, true
relative residuals below \(5\times10^{-8}\), and a cross-code discrepancy
below 1%. The native normalization check uses that 1% criterion.
Angular and pitch resolution are held equal in this fixture. This check
does not claim a continuum-converged machine prediction.

The raw HDF5 SHA-256 values are
`1753630c6c8bc23e5bc3997bb5c6cb7af9dd6c7912475ae4adae6e887fca23d7`
for Lorentz and
`0aec64c4eaf97564b75bad194e9a01590bb959022e1d9c6699af2d7084902bd6`
for full self-collisions.
