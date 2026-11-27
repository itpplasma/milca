# Label lane on a real stellarator: Landreman–Paul QA reference

Reference data for the ideal label lane (docs/equilibrium-roadmap.md, step 4)
on the fixed-boundary, vacuum (p = 0) Landreman–Paul QA deck
`benchmarks/nested/input.landreman-paul-qa` (NFP 2, provenance in
`benchmarks/nested/README.md`).

The deck prescribes zero net toroidal current (NCURR = 1, CURTOR = 0). The
label lane prescribes the flux profiles, so the benchmark problem is the
same boundary with the rotational transform of the current-free VMEC++ state,
fitted by a power series of degree 6 in s (normalized toroidal flux) and
prescribed: `input.landreman-paul-qa-iota` (NCURR = 0, AI = fit). In label
form: chi = s, F(chi) = -PHIEDGE chi/(2 pi), G'(chi) = -iota(chi)
PHIEDGE/(2 pi), p = 0 (signs as in MILCA's E4 convention; check against
VMEC's signgs = -1 of this deck).

## Commands

```sh
python make_reference.py --xvmec ~/code/VMEC2000/build-cm/build/bin/xvmec   # 2.3 min
python compare.py MILCA_OUTPUT_DIR
```

`make_reference.py` runs VMEC++ (Python package) and VMEC2000 on both decks
(ns 16, 51, 101; FTOL 1e-14; MPOL 16, NTOR 12 as in the deck; VMEC++ 8
threads), writes the prescribed-iota deck and, from the prescribed-iota
VMEC++ state, the reference files to `benchmarks/label-qa/reference/` of an
external data directory (`MILCA_DATA_DIR`, or `--out`); wout files go to the
sibling `work/`. `compare.py` reads the reference from there.

## Reference files

| file | content |
|---|---|
| `iota_coefficients.txt` | c_0..c_6 of iota(s) = sum c_k s^k; PHIEDGE, NFP in the header |
| `iota.txt` | s, iota of VMEC++ (full grid, ns 101) |
| `axis.txt` | phi, R, Z of the magnetic axis on 64 cylindrical planes phi = 2 pi k/(64 NFP), one period |
| `surface_s0.25.txt`, `_s0.50`, `_s1.00` | phi, VMEC theta, R, Z, \|B\| on the same planes, 128 equally spaced theta |
| `summary.json` | versions, times, residuals, energy, the differences below |

Coefficients at s between full-grid surfaces are linear in s; \|B\| (half
grid) is linear between half-grid surfaces and extrapolated linearly to
s = 1. Energy W = int \|B\|^2/2 dV over the torus with B in tesla and without
mu0: 0.281794 T^2 m^3 (W/mu0 = 2.2424e5 J); it equals VMEC's wb (2 pi)^2 to
2e-11. Volume 0.564712 m^3.

`compare.py DIR` reads a MILCA result in the same layout (`axis.txt`,
`surface_s*.txt` with any theta on the 64 planes, optional `iota.txt` and
`summary.json` with `energy_T2m3`) and prints axis max \|dR\|, \|dZ\|, the
distance of each point to the reference curve of its plane (trigonometric
interpolant, Newton on the closest point; max and rms), the relative \|B\|
difference at the closest point, iota and energy differences.

## Measured (Apple M-class, 2026-11-22)

| quantity | value |
|---|---|
| VMEC++ 0.9.3, ns 16/51/101 | 9.1 s (current-free), 9.0 s (prescribed iota) |
| VMEC2000 (`build-cm`, one process) | 53.8 s, 52.8 s |
| final residuals (fsqr, fsqz, fsql), prescribed iota | 1.0e-14, 4.8e-15, 1.9e-15 |
| iota, axis to edge | 0.42304 to 0.41577 |
| degree-6 fit, max error | 4.4e-5 |
| prescribed vs current-free (same VMEC angles) | axis 1.6e-7 m, s 0.25: 8.5e-5 m, s 0.5: 7.4e-5 m, iota 4.5e-5 |
| VMEC2000 vs VMEC++, prescribed iota (same angles) | axis 1.1e-6 m, s 0.5: 2.6e-5 m |
| VMEC2000 vs VMEC++, s 0.5, normal distance (`compare.py`) | max 8.4e-6 m, \|B\| 6.5e-6 |
| VMEC++ radial resolution, max \|rmnc(ns 151) - rmnc(ns 101)\| | axis 4.4e-5 m, s 0.25: 9.3e-6 m, s 0.5: 7.4e-6 m |

The iota fit changes the surfaces by less than 1e-4 m (minor radius about
0.17 m), so the prescribed-iota problem is the QA equilibrium to that level.
VMEC converges at first order in 1/ns; from ns 101 and 151 the reference
error is about 1.3e-4 m at the axis and 2e-5 m at s = 0.25 and 0.5, the
level below which a MILCA difference is not resolved. GVEC and DESC are not
installed in this environment (`benchmarks/nested/run_reference.py` has
adapters for both).

## Provenance

VMEC++ 0.9.3 (Python package, github.com/proximafusion/vmecpp, MIT);
VMEC2000 (github.com/itpplasma/VMEC2000, fork of hiddenSymmetries/VMEC2000,
MIT, commit fff954b025e4d0934a011c3a4c376d4d619807fa, local build
`build-cm`).
Only their outputs are used; no code of either enters MILCA.
