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

`ladder.sh CASE XVMEC PYTHON MILCA_BUILD THREADS` runs the speed and accuracy
ladders of VMEC2000, VMEC++ (Python package) and MILCA on a case directory.

Results (beta 1.7 %, iota 0.39 to 0.28, 41 x 41 grid, 64 rings; the grid
spacing does not change them). Errors against the exact solution; wall time
on an Apple M-class machine, VMEC2000 one process, VMEC++ and MILCA 8
threads; VMEC mpol 8 up to ns 64, 12 up to ns 256, 16 at ns 512; MILCA
degree 3 with the given boundary nodes:

| code | resolution | time | volume | axis | wb |
|---|---|---|---|---|---|
| VMEC2000 | ns 32 | 0.32 s | 1.0e-3 | 7.5 mm | 2.1e-5 |
| VMEC++ | ns 32 | 1.9 s | 1.0e-3 | 7.5 mm | 2.1e-5 |
| VMEC2000 | ns 64 | 0.54 s | 7.4e-4 | 3.4 mm | 9.5e-7 |
| VMEC++ | ns 64 | 2.9 s | 7.4e-4 | 3.4 mm | 9.5e-7 |
| VMEC2000 | ns 128 | 1.45 s | 3.6e-4 | 1.4 mm | 2.5e-6 |
| VMEC++ | ns 128 | 4.4 s | 3.6e-4 | 1.4 mm | 2.5e-6 |
| VMEC2000 | ns 256 | 3.9 s | 3.4e-4 | 0.75 mm | 3.9e-6 |
| VMEC++ | ns 256 | 6.5 s | 3.5e-4 | 0.76 mm | 3.9e-6 |
| VMEC2000 | ns 512 | 17.8 s | 2.1e-4 | 0.36 mm | 2.5e-6 |
| VMEC++ | ns 512 | 10.9 s | 2.1e-4 | 0.36 mm | 2.5e-6 |
| MILCA | 12 nodes | 1.67 s | 2.8e-4 | 0.53 mm | 2.2e-4 |
| MILCA | 16 nodes | 1.64 s | 3.5e-4 | 0.83 mm | 6.2e-5 |
| MILCA | 24 nodes | 3.5 s | 4.9e-4 | 1.1 mm | 3.1e-6 |
| MILCA | 32 nodes | 6.2 s | 2.7e-4 | 0.51 mm | 1.6e-6 |
| MILCA | 48 nodes | 20 s | 5.0e-5 | 0.11 mm | 2.6e-7 |
| MILCA fixed boundary on the exact D | automatic | 0.9 s | 6e-6 | 0.14 mm | 2e-6 |

VMEC converges at first order in 1/ns. At equal time (about 1.5 s) MILCA
matches VMEC2000 in volume, is better in the axis and worse in wb (open); at
equal accuracy MILCA at 32 nodes matches VMEC at ns 512 in a third to half
of the time, and at 48 nodes it is beyond what VMEC reaches with ns <= 512.
MILCA's error is not yet monotone in the mesh (24 nodes), which makes a cheap
level unreliable (open).

How MILCA got there (the first run from the elliptic start was 2.1 % off in
volume and 4 cm in the axis at every mesh):

- The frozen field of the seed carried its field-line labelling through the
  free-boundary solve. Axisymmetric decks now seed again from the solved
  boundary until wb and the axis change by less than 1e-5, on 16 boundary
  nodes first, then once or twice on the requested mesh.
- The plasma-boundary single layer was over-resolved: near/far Gauss orders
  4/3 and 8 graded radial points at up to 24 nodes, 6/4/12 up to 48, 8/6/20
  above (MILCA_BEM_ORDERS overrides). FortFEM assembles each separated panel
  pair once.
- The Newton step cap is 5 % of the plasma extent (was 2 %; 22 instead of
  40 iterations from the elliptic start; MILCA_STEP_CAP overrides).
- MILCA_FB_PROFILE=1 prints the wall time of the free-boundary pieces.

Open: wb on coarse meshes, the non-monotone mesh error, the shell vacuum on
this case, and the 3D free boundary of VMEC++'s CTH-like deck (imported, but
the solve meets a folded map; VMEC++ solves it in 0.8 s).
