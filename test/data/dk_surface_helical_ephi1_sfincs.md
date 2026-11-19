# SFINCS reference for the potential response of the inductive drive

`dk_surface_helical_ephi1_sfincs.txt` holds the part of the parallel flow of
the inductive drive that is linear in a prescribed potential \(\Phi_1\), on
the helical-core file of
[dk_surface_helical_sfincs.md](dk_surface_helical_sfincs.md) with
\(\delta = 0.03\), for `test_dk_surface_device sfincs_phi1`.

## Code and solver

[SFINCS](https://github.com/landreman/sfincs/tree/8df5453472e982df0f6ae005243ce38d57a83711)
at 8df5453 (BSD 2-Clause), PETSc 3.18 with MUMPS, direct solve with
`-mat_mumps_cntl_1 0.01 -mat_mumps_icntl_10 10` (iterative refinement).
Relative residuals: \(5\times10^{-10}\) or less. Without the refinement
they were \(2\times10^{-7}\) to \(7\times10^{-7}\) for the full operator.

## Case

The geometry, species, resolution and collisionality are those of the
\(\nu_n = 0.005\) cases of
[dk_surface_helical_sfincs.md](dk_surface_helical_sfincs.md). The density and
temperature gradients are zero (`dNHatdpsiNs=0`, `dTHatdpsiNs=0`).
`includePhi1=readExternalPhi1=.true.`, `includePhi1InKineticEquation=.true.`,
`includePhi1InCollisionOperator=.false.`. In MILCA's reading of the file,
\(Z=1\) and \(T\) equal to the reference temperature,

\[
\varepsilon = \frac{Ze\Phi_1}{T} = 10^{-4}\left[0.4 + \cos\theta + 0.6\sin\theta
+ 0.5\cos(2\theta-\zeta+0.3)\right].
\]

SFINCS reads the file as its mirror image, so `Phi1Hat` is given at
\((\theta, -\zeta)\). Its `BHat` agrees with MILCA's \(|B|\) at \((\theta,-\zeta)\)
to \(4\times10^{-14}\).

The inductive drive is proportional to \(F = f_M(\psi, H)\), so it carries the
Boltzmann factor \(e^{-\varepsilon}\). SFINCS 8df5453 builds this drive from
`EParallelHat` without that factor. The runs therefore set `EParallelHat=0`,
`EParallelHatSpec=1` and pass the profile \(1 - \varepsilon\) through
`EParallelHatSpec_bcdatFile`. That file has the layout of the `.bc` file, with
data lines `m n cos sin` in the file's convention, which SFINCS mirrors like
\(|B|\). The profile path with the profile 1 reproduces `EParallelHat=1` exactly.

Each entry is the part linear in \(\Phi_1\) of runs at \(0\), \(\Phi_1\) and
\(2\Phi_1\), \(L = [4(u_1-u_0) - (u_2-u_0)]/2\). It is divided by
`FSABFlow` of the inductive drive without \(\Phi_1\), so the normalisation of
the force cancels. The quadratic part is below \(2.4\times10^{-4}\) of \(L\).

| collision operator | \(L/\langle Bu_\parallel\rangle\) | parallel-field term alone |
|---|---|---|
| pitch-angle scattering | \(9.9760\times10^{-5}\) | \(1.3976\times10^{-4}\) |
| full operator | \(8.3845\times10^{-4}\) | \(8.7845\times10^{-4}\) |

The parallel-field term alone (runs with `EParallelHat=1` and no profile) is
larger by the mean of \(\varepsilon\), \(4\times10^{-5}\). The varying part of
\(\varepsilon B\) has \(\langle\varepsilon B^2\rangle = 0\) in the Boozer measure
\(\propto 1/B^2\), so it is a parallel gradient. Its response is a Boltzmann
density without flow: without the constant, both runs agreed to \(10^{-8}\).

## Format

The first line is the number of cases. Each case has three lines: the MILCA
collision model (1 pitch-angle scattering, SFINCS `collisionOperator=1`; 3 full
operator, `collisionOperator=0`); \(L/\langle Bu_\parallel\rangle\) for
`FSABFlow`; and the same ratio for `flow` at the 77 points of MILCA's grid
(theta fastest), in units of 1/T.
