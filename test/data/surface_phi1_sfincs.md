# SFINCS reference for the in-surface potential drive

`surface_phi1_sfincs.txt` holds the change of the parallel flow caused by a
prescribed potential \(\Phi_1(\theta,\zeta)\) at a fixed density or temperature
gradient, computed with SFINCS in its linear external-\(\Phi_1\) mode, as the
reference for `test_surface_phi1 sfincs`.

## Code and solver

[SFINCS](https://github.com/landreman/sfincs/tree/8df5453472e982df0f6ae005243ce38d57a83711)
at 8df5453 (BSD 2-Clause), PETSc 3.18 with MUMPS, direct solve. SFINCS sets the
MUMPS pivoting threshold to \(10^{-6}\), with which these solves leave relative
residuals up to 1. The runs therefore pass `-mat_mumps_cntl_1 0.01` (the MUMPS
default), which gives residuals of \(10^{-11}\) to \(10^{-15}\). The same change is
the default of the itpplasma fork
([fix/mumps-pivot-threshold](https://github.com/itpplasma/sfincs/tree/fix/mumps-pivot-threshold)).

## Case

Field, species, resolution and collision settings are those of
[surface_fluxes_sfincs.md](surface_fluxes_sfincs.md) (three-helicity Boozer
field, deuterium at 2 keV and \(3\times10^{19}\,\mathrm m^{-3}\), 11 x 7
angles, 24 Legendre modes, 16 speed nodes), with `RHSMode=1`,
`dPhiHatdpsiHat=0`, `constraintScheme=2`, `includePhi1=readExternalPhi1=.true.`,
`includePhi1InKineticEquation=.true.`, `includePhi1InCollisionOperator=.false.`.
The potential is

\[
\hat\Phi_1 = 10^{-4}\left[\cos\theta + 0.6\sin\theta + 0.5\cos(2\theta-3\zeta+0.3)\right],
\]

in units of 2000 V. Drive 1 has \(d\hat n/d\hat\psi = 1\); drive 2 has
\(d\hat T/d\hat\psi = 1\) and \(d\hat n/d\hat\psi = 1.5\), so that the forces are
\(f_1 = 1, f_2 = 0\) and \(f_1 = 0, f_2 = 1\) per weber.

Each entry is the part linear in \(\Phi_1\) of runs at \(0\), \(\Phi_1\) and
\(2\Phi_1\), \(L = [4(u_1-u_0) - (u_2-u_0)]/2\). The quadratic part at
\(\Phi_1\) is below \(2\times10^{-3}\) of \(L\) except for the full operator with
drive 2, where \(L\) nearly cancels (0.0114 T m/s against a quadratic part of
0.0100 T m/s).

## Format and units

The first line is the number of cases. Each case has three lines: the MILCA
collision model (1 pitch-angle scattering, SFINCS `collisionOperator=1`; 3 full
operator, `collisionOperator=0`) and the force index; \(\langle B u_\parallel\rangle\)
in T m/s; \(u_\parallel\) in m/s at the 77 grid points, theta fastest. SFINCS
`flow` and `FSABFlow` are converted with \(\bar v = (2\cdot 2000\,\mathrm{eV}/m_D)^{1/2}\).
