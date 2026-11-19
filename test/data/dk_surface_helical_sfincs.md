# SFINCS reference on a helical-core Boozer file

`dk_surface_helical_sfincs.txt` holds SI fluxes for a density gradient on the
synthetic one-period file of `test_dk_surface_device` with helical amplitude
\(\delta = 0.03\), computed with
[SFINCS](https://github.com/landreman/sfincs/tree/8df5453472e982df0f6ae005243ce38d57a83711)
at 8df5453 (BSD 2-Clause), `geometryScheme=12`, PETSc 3.18 with MUMPS and
`-mat_mumps_cntl_1 0.01` (relative residuals \(10^{-10}\) to \(10^{-14}\)).

## File and case

\(|B| = B_0[1 - \varepsilon r\cos\theta + \delta r(\cos(\theta-\varphi_B) +
0.5\sin(\theta-\varphi_B))]\), \(r=\sqrt s\), \(B_0 = 2\) T, \(\varepsilon = 0.15\),
\(\iota = 1.2\), \(G = 3.4\) T m, \(I = 0\), toroidal flux 0.785 Wb, one field
period, surface \(s = 0.5\) (`psiN_wish=0.5`). SFINCS reads lines into a
200-character buffer, so the file for SFINCS uses 19-character fields; the
values are those the test writes. SFINCS converts the file to a right-handed
system and reads it as its mirror image (\(\psi_a = -0.1249\), \(\iota = -1.2\));
the fluxes per \(\psi\) and \(\langle Bu_\parallel\rangle\) below agree in sign
with MILCA's reading.

Deuterium at 2 keV and \(3\times10^{19}\,\mathrm m^{-3}\), `dNHatdpsiNs=-1`
(\(dn/ds = -n\)), no temperature gradient or radial electric field,
`constraintScheme=2`, 11 x 7 angles, 24 Legendre modes, 16 speed nodes,
`nu_n` = 0.005 and 0.05, which MILCA represents by the Coulomb-logarithm input
136.0478 and 1360.478 of `surface_fluxes_sfincs.md`.

## Format and units

The first line is the number of cases. Each case: MILCA collision model
(1 pitch-angle scattering = SFINCS `collisionOperator=1`, 3 full operator =
`collisionOperator=0`), MILCA collisionality input, particle flux
\(\langle\Gamma\cdot\nabla\psi\rangle\) [Wb m\(^{-3}\) s\(^{-1}\)],
kinetic-energy flux \(\langle Q\cdot\nabla\psi\rangle\) [W Wb m\(^{-3}\)], with
\(\psi\) the toroidal flux over \(2\pi\), and \(\langle Bu_\parallel\rangle\) [T m/s]. From SFINCS:
\(\Gamma = n_0 v_0|\psi_a|\,\mathtt{particleFlux\_vm\_psiN}\),
\(Q = 2 n_0 T_0 v_0|\psi_a|\,\mathtt{heatFlux\_vm\_psiN}\),
\(\langle Bu_\parallel\rangle = v_0\,\mathtt{FSABFlow}\), \(v_0 = (2T_0/m_D)^{1/2}\).
