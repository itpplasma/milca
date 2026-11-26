# Flux map on standard fixed-boundary decks

Opt-in solves for the 3D flux map (`docs/flux-map.md`, "Section meshes by
harmonic maps"). Each directory holds the deck
with the rotational transform prescribed (`input.<name>-iota`, written by
`make_reference.py` from the VMEC++ state of the deck as given) and the
VMEC++ reference in `reference/` (axis, surfaces s = 0.25, 0.5, 1 with |B|,
`summary.json`).

| case | deck | content |
|---|---|---|
| `qa` | Landreman–Paul QA | vacuum, NFP 2, weak shear |
| `qh` | `input.nfp4_QH_warm_start` of vmecpp | vacuum, NFP 4 |
| `cth` | `input.cth_like_fixed_bdy` of vmecpp | pressure, current, strong shear, NFP 5 |
| `w7x` | `input.w7x` of vmecpp | W7-X standard, pressure, NFP 5; bean and triangle sections |
| `li383` | `input.li383_low_res` of VMEC2000 | NCSX-like, pressure, current, NFP 3 |

Each stage starts from the previous solution. Specify a stage as
`degree:cells:planes[:radial]`. The central block has `cells` cells in each
direction. Each outer block has `cells` along the wall and `radial` across
it. Omit `radial` to use `cells` in both directions.

```sh
cmake --build build --target test_flux_map_qa
OMP_NUM_THREADS=8 build/test_flux_map_qa benchmarks/flux-map-standard/w7x/input.w7x-iota OUT 2:2:9,2:4:13
python benchmarks/label-qa/compare.py OUT benchmarks/flux-map-standard/w7x/reference
```

The driver prints per stage the mesh report (wall samples, residual of the
harmonic maps, wall distance, area error, certified Jacobian ratio), the
Newton steps, the energy, the number of negative Hessian eigenvalues, the
force balance J×B − ∇p of the state with J = ∇×B_h (relative to
|∇(B²/2)| and |∇p|, over the volume, over s > 0.1 and per band of s), the
range of ∇φ·(∇x×∇y) over the quadrature points, the winding of the wall
angle and the jump of B²/2 at the block interfaces,
and writes `axis.txt`, `surface_s*.txt`, `current.txt` (circulation of B
along flux surfaces: the enclosed toroidal current) and `summary.json`.
Diagnostic switches: `MILCA_FM3_PSCALE` (pressure factor),
`MILCA_FM3_FREE_AREA` (no area correction of the mesh wall),
`MILCA_FM3_CORE` (core fraction of the disc O-grid), `MILCA_FM3_VMEC`,
`MILCA_FM3_HESSIAN`, `MILCA_FM3_GRADING` (curvature grading of the section
meshes, 0 = arc length), `MILCA_FM3_BANDS` (bands of s of the VMEC report).

The harmonic-disc boundary angles follow the deck's poloidal sample spacing
by default. This fixes the poloidal angle origin across the toroidal planes,
so resampling the wall adds no toroidal harmonics. `MILCA_FM3_SAMPLE_ANGLE=0`
places the wall by arc length, and `MILCA_FM3_GRADING=c` with c > 0 by
curvature grading. The standard stages are `2:2:9,3:4:21,4:6:21`, and the
finer 17-plane configuration uses `2:2:9,3:4:17,5:8:17:4`.

With `MILCA_FM3_VMEC=<wout>` the force balance of the VMEC state is printed
first, on the same measure and without a mesh; the stages argument `vmec`
stops after it:

```sh
MILCA_FM3_VMEC=benchmarks/flux-map-standard/w7x/work/wout_w7x_iota.nc \
  build/test_flux_map_qa benchmarks/flux-map-standard/w7x/input.w7x-iota OUT vmec
```

`checks/`: `desc_axis.py` (third-code axis with DESC),
`vmec_iota_scale.py` (VMEC++ with a scaled transform),
`hessian_toroidal_content.py` (lowest Hessian modes and their toroidal
harmonics), `harmonic_map_check.py` (the harmonic map in numpy, independent
of the Fortran implementation). Measured numbers: `docs/flux-map.md`.
