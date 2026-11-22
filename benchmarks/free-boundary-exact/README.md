# Exact free-boundary benchmark (axisymmetric Solov'ev)

`solovev_exact.py OUT` writes an exact free-boundary problem for every
VMEC-class code: the Solov'ev equilibrium psi = alpha (R^2 - r0^2)^2 +
beta R^2 Z^2 with constant F inside D = {psi < psi_b}, held by ring coils
outside the mgrid box whose currents are fitted to the exact external field
B_ext = B_Solovev - B_plasma (B_plasma by the ring Green's function over D,
polar quadrature about each point). The fit is exact to about 1e-7 of
|B_pol| on independent points inside D and on its boundary (64 rings), and
the field is smooth on the grid, so a coarse mgrid (41 x 41) suffices.

Outputs: `mgrid_solovev_exact.nc` (one coil group, unit current),
`input.solovev_exact` (VMEC INDATA: the exact iota and pressure as dense
cubic splines in s, PHIEDGE, an elliptic starting boundary),
`coils.solovev_exact` (R, Z, I of the rings) and `reference.json` (volume,
wb, betatotal, the axis R = r0, iota at axis and edge, psi_b).

`compare.py reference.json RESULT.nc [label]` prints the errors of a VMEC
wout or a MILCA `nested_equilibrium` result: relative volume, wb, beta,
edge iota, the axis distance and, for a wout, the largest distance of the
s = 1 surface from the exact boundary.

Results (beta 1.7 %, iota 0.39 to 0.28, 41 x 41 grid, 64 rings; the grid
spacing does not change them):

| code | resolution | time | volume | axis | boundary |
|---|---|---|---|---|---|
| VMEC2000 | mpol 8, ns 64 | 0.48 s | 7.4e-4 | 3.4 mm | 4.0 mm |
| VMEC2000 | mpol 12, ns 128 | 1.3 s | 3.6e-4 | 1.4 mm | 2.2 mm |
| VMEC++ | mpol 8, ns 64 | 2.1 s | 7.4e-4 | 3.4 mm | 4.0 mm |
| MILCA, vacuum bem | automatic | 11 s | 2.1e-2 | 41 mm | |
| MILCA, vacuum shell | automatic | fails (status 6) | | | |

VMEC converges at first order in 1/ns. MILCA's axisymmetric BEM vacuum
converges to a wrong state (its own mesh-error estimate does not fall);
the shell vacuum does not converge at this beta. Both are open.
