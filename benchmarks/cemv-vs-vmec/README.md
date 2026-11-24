# Fixed-boundary equilibrium benchmark

MILCA's VMEC mode maps flux coordinates to physical space. Its CEMV mode
maps physical space back to flux coordinates. Compare both with stock
[VMEC2000](https://github.com/ORNL-Fusion/VMEC2000) on the same boundary,
toroidal flux, rotational transform and pressure. Publication-quality
comparisons also require a converged reference for the magnetic axis, flux
surfaces and magnetic field. See the [equilibrium roadmap](../../docs/equilibrium-roadmap.md).

## Cases

| Case | Input | Independent check |
|---|---|---|
| H1 | Integer-transform and sheared Landreman equilibria from `test_equilibrium_3d deck` | Analytic magnetic field for both cases. Analytic pressure and flux label for the integer case. |
| H2 | Quasi-axisymmetric vacuum case in [label-qa](../label-qa/) | Radial force-balance bands and resolution of the wall |
| H3 | Compact Toroidal Hybrid case in [flux-map-standard](../flux-map-standard/cth/input.cth-iota) | Radial force-balance bands with finite pressure and shaped geometry |

The H1 decks use Fourier approximations to the analytic boundaries. Generate
the integer case with `test_equilibrium_3d deck INPUT` and the sheared case
with `test_equilibrium_3d deck INPUT landreman-sheared`. Set
`MILCA_FM3_EXACT` to `landreman-integer` or `landreman-sheared` to report
volume magnetic-field error. The `equilibrium_d8_*` tests use equal-weight
samples for the integer case.

The sheared CEMV benchmark starts from polar labels. Subsequent stages
place their outboard section anchors at the coarse magnetic axis height
and solve for stationarity.

With magnetic field \(B\), pressure \(p\) and vacuum permeability \(\mu_0\),
let \(s\) be normalized toroidal flux. The common force diagnostic is

\[
\eta_{s>0.1} =
\frac{\| (\nabla\times B)\times B - \nabla(\mu_0 p)\|_{L^2(s>0.1)}}
     {\|\nabla(B^2/2)\|_{L^2(s>0.1)}} .
\]

The drivers also report the whole-volume norm and ten flux bands. The
finite-element integral samples cell interiors. For CEMV, also inspect the
separate block-interface jump and wall distance. Stock VMEC's whole-volume
norm is sensitive to radial differentiation near the axis. The
`test_flux_map_qa ... vmec` check reconstructs the field from Wout
geometry and straight-field-line angle. An independent check using the
saved magnetic-field coefficients gives different force estimates on the
same physical points. The Lagrangian diagnostic matches poloidal orientation
to the frozen-in field and rejects the result if quadrature reveals a fold
in the geometry or flux labels.

## Runs

From the checkout root, set `build` to the absolute CMake build path, `deck`
to the prescribed-transform VMEC input, `nested_deck` to the same case in
MILCA's `nested_input` syntax, and `xvmec` to the stock VMEC2000 executable.
Use the same CPU allocation and binding for comparisons, retain the exact
decks, and record the compiler, library, thread and rank settings.

```sh
cmake --build "$build" --target test_equilibrium_3d test_flux_map_qa nested_equilibrium milca_import
benchmarks/cemv-vs-vmec/run_cemv.sh "$build" "$deck" cemv 3:4:21 8
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MILCA_NESTED_EXPORT="$PWD/lagrangian" \
  "$build/nested_equilibrium" --config "$nested_deck"
benchmarks/cemv-vs-vmec/run_stock.sh "$xvmec" "$deck" stock mpirun -np 8 --bind-to core
MILCA_FM3_VMEC="$PWD/stock/wout_benchmark.nc" \
  "$build/test_flux_map_qa" "$deck" stock vmec
```

Run the sheared CEMV case at seven planes per field period.

```sh
"$build/test_equilibrium_3d" deck input.landreman-sheared landreman-sheared
MILCA_FM3_EXACT=landreman-sheared \
  benchmarks/cemv-vs-vmec/run_cemv.sh "$build" input.landreman-sheared cemv-sheared 2:4:7 8 2:2:7
```

For finer CEMV output, refine the in-plane mesh before the toroidal grid.

```sh
MILCA_FM3_EXACT=landreman-sheared \
  benchmarks/cemv-vs-vmec/run_cemv.sh "$build" input.landreman-sheared cemv-sheared-fine 2:8:13 8 2:2:7,2:4:7,2:8:7
```

The sheared forward map uses the stationary solver to reach an energy
saddle at the prescribed rotational transform. Import the same deck.

```sh
"$build/milca_import" --format vmec input.landreman-sheared landreman-sheared.nml
```

Replace the importer's `angle_mpol = 4, max_iterations = 400,` line
inside the `nested_input` group with these controls.

```fortran
stationary = .true.,
degree = 4, boundary_nodes = 16, map_ntor = 8, angle_mpol = 6,
ladder_nodes = 12, ladder_degree = 3,
tolerance = 1e-5, max_iterations = 48,
```

```sh
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MILCA_NESTED_EXPORT="$PWD/lagrangian-sheared" \
  "$build/nested_equilibrium" --config landreman-sheared.nml
```

Pass the local MPI launcher and its arguments after the stock output directory.
With no launcher the script runs one rank. The runner checks the normal exit
message and Wout exit flag and rejects `Try increasing NITER`. Compare the
final `FSQR`, `FSQZ` and `FSQL` in `run.log` with the final `FTOL_ARRAY`
entry in the deck. Record each MILCA Newton status. The geometry comparison
rejects forward output flagged as unconverged. Check accuracy before comparing
times. Measure setup through converged output and record postprocessing
separately.

For geometry and field comparisons, export a converged symmetric Wout with
`export_wout.py` (NumPy and netCDF4). Compare two Wout refinements before
using the finer one as a reference.

```sh
python benchmarks/cemv-vs-vmec/export_wout.py stock/wout_benchmark.nc stock
python benchmarks/label-qa/compare.py cemv stock
python benchmarks/label-qa/compare.py lagrangian stock
```

For the sheared comparison, use `cemv-sheared-fine` and
`lagrangian-sheared` with a converged reference from the same sheared deck.

Set `MILCA_FM3_RADIAL` to increase the number of integration points in the
Wout force check and quantify radial sampling error.

`MILCA_FM3_GRADING=1` selects curvature-graded sections for CEMV.
The meshes must pass the geometry check. The CEMV stage syntax is poloidal
degree, cells per block direction, and number of toroidal planes.
CEMV uses the input boundary's winding to convert the sign of the
rotational transform to its counterclockwise poloidal convention.
`MILCA_NESTED_BENCHMARK_LAGGED=1` experiments with reusing the Lagrangian
Newton Hessian. The default uses the original solver.
