# Foreign fixed-boundary equilibrium inputs

`milca_import [--format auto|gvec|desc|spectre|vmec|wout] INPUT NAMELIST`
writes the `nested_input` namelist of
[nested_equilibrium](nested-equilibrium-input.md) for an input of another
equilibrium code. Exit status 0 means the namelist is written, 3 that the
input is read but describes physics `nested_equilibrium` does not model (the
reason goes to stderr), 1 that the input cannot be read. Notes on the
conversion and its errors follow the namelist group in the written file.

## Conventions

Every reader produces one intermediate form, the convention of VMEC:

    R = sum_{m,n} rbc(n,m) cos(m t - n N phi) + rbs(n,m) sin(m t - n N phi),      (1)

likewise Z with zbc, zbs, where phi is the cylindrical angle counterclockwise
seen from above, N the number of field periods and t a poloidal angle of
either orientation. iota is measured in (t, phi), the enclosed toroidal
current I_V has the sign of VMEC's CURTOR, and the toroidal flux Phi keeps its
sign. `nested_equilibrium` wants t counterclockwise in (R, Z), Phi > 0 and
its own current I, so the last step is

    A = (1/2) oint (R dZ - Z dR) at phi = 0,  and t -> -t if A < 0,              (2)
    t -> -t:  (n, m) -> (-n, m) for m > 0,  rbs, zbs -> -rbs, -zbs,  iota -> -iota,   (3)
    Phi -> |Phi|,  I = -sign(Phi) I_V.                                             (4)

The current keeps its sign under (3) because it is physical. (4) matches the
benchmark converter and the PARVMEC comparison of the
[input page](nested-equilibrium-input.md). Modes beyond m <= 64, |n| <= 32 or
below 1e-12 of the minor radius are dropped, and the sum of their |c| bounds
the displacement of the boundary. Above 1e-3 of the minor radius the input is
refused.

A profile is a power series in s, the normalised toroidal flux, or a
piecewise cubic Hermite profile. Several Hermite profiles share the union of
their knots, which leaves each unchanged.

## Formats

| Source | Read | Boundary | Profiles |
|---|---|---|---|
| GVEC parameter file (INI, TOML, YAML), hmap 1 | X1_b_*, X2_b_* within X?_mn_max and X?_sin_cos, plus X?pert_b_* when boundary_perturb | (5) | polynomial exact, interpolation as the not-a-knot spline exact, B-spline of degree <= 3 exact, higher on 257 added knots, error noted |
| GVEC, boundary from file | boundary/X, Y on (t_i, zeta_j), projected on m <= m_max, abs(n) <= n_max | (5) | as above |
| GVEC whichInitEquilibrium = 1 | the VMEC wout file it names | as wout | as wout |
| GVEC hmap 20 (Frenet), 21 (G-frame) | axis and frame, the G-frame from axis/xyz, Nxyz, Bxyz | (6) | as above, with the factor s_zeta |
| VMEC wout | rmnc, zmns (rmns, zmnc) of the last surface | direct | iotaf, presf as the not-a-knot spline on the uniform s grid, as GVEC builds them |
| VMEC &INDATA | RBC ... ZBS, NCURR, AM, AI, AC, CURTOR | direct | power_series only, exact |
| VMEC &VSEQ | its last EXTENSION | as INDATA | as INDATA |
| DESC text input | R1, Z1 on G_m(t) G_n(N phi), p, i, c in powers of rho | grid fit (6) | even powers of rho exact, odd ones on 257 knots, error noted |
| SPECTRE TOML, igeometry 3 | physics.rbc ... keyed (m, n), tflux, iota, ivolume, isurf | direct, then (2) | lconstraint 1: iota, 3: enclosed current, not-a-knot spline through the interface values |

GVEC maps (X1, X2, zeta) to x = (X1 cos zeta, -X1 sin zeta, X2) for hmap 1,
so phi = -zeta and

    rbc(-n, m) = X1_b_cos(m; n),  zbs(-n, m) = X2_b_sin(m; n),
    iota_V = -iota_G,  Phi_V = -Phi_G,  I_V = -I_G,                                (5)

the flips of GVEC's own `flip_parameters_zeta`. For hmap 20 and 21 the
surface x(t, zeta) = X0(zeta) + X1 N(zeta) + X2 B(zeta) is cut by the planes
phi = phi_j: on each line t = t_i, zeta(phi) is found by regula falsi on a
table of phi(zeta), and

    R(t_i, phi_j), Z(t_i, phi_j)  ->  (1) by a discrete Fourier fit,               (6)

with m up to 2 m_max + 8 and |n| up to 64. iota, Phi and I take the factor
s_zeta = sign(dphi/dzeta). The surface is admitted only if phi is monotone in
zeta on every line t = t_i and turns once, which excludes knots. The largest
distance between the fit and the cut at the midpoints of the fit grid is the
stated geometry error. DESC's product basis is evaluated on a grid and fitted
by (6) without a cut. Its current is c = signgs I_V with signgs = -1 for a
counterclockwise boundary (DESC's VMEC reader).

## Refused inputs

GVEC hmap 3 (a straight periodic cylinder) and hmap 10 (a knot), G-frame axes
winding more than once, boundaries not resolved by 64 toroidal modes over phi
(fit error above 1e-3 of the minor radius), SPECTRE slab or cylinder geometry,
free boundary, a stepped pressure (not C1) or lconstraint other than 1 and 3,
VMEC free-boundary decks (LFREEB true with an MGRID_FILE), VMEC profile types
other than power_series, and VMEX &MIRROR (an open mirror).

## Evidence

CTest `equilibrium_import` (label fast) writes its fixtures and checks:

- GVEC hmap 1 in INI, TOML and YAML spelling against (5) by hand, and that
  the written namelist reads with `nested_equilibrium`'s declarations,
- a clamped cubic B-spline against Cox-de Boor to 1e-14, a not-a-knot spline
  through a line,
- a Frenet circle R0 = 3 with a helical X1 mode: (6) and (2) give
  rbc(0,1) = -0.5, rbc(-1,1) = -0.05, zbs(0,1) = -0.6 to 1e-12,
- a G-frame with N tilted out of the phi planes: every raw grid point lies on
  the imported surface to 1e-8 m,
- a clockwise wout with cubic iotaf, presf: flipped boundary and exact
  profiles,
- DESC's DSHAPE input against the VMEC DSHAPE deck it came from, a 3-D DESC
  boundary and the sign of its current,
- SPECTRE interface currents and the flip, VSEQ, INDATA repeat counts, and
  the refusals.

On the benchmark inputs the importer agrees with the benchmark's VMEC
converter on the VMEC twin of each case:

| Case | Boundary distance (m) | Profiles |
|---|---|---|
| GVEC ellipstell_lowres, w7x, current_profile, bspline_profiles, advanced_control | <= 1.8e-15 | equal |
| GVEC tilted_elliptok_hmap_RZ (boundary file) | 1.8e-11 | not compared |
| GVEC wout of test-signgs-{neg,pos}-asym{F,T}, d23p4_tm | <= 1.6e-15 | wout profiles |
| DESC DSHAPE, QSC_r2_5.5 | 1.8e-15, 5.5e-9 | equal, current sign included |

The GVEC G-frame N2-63 boundary lies within 3.2e-7 m of the imported surface
at every point of the file grid. The QH nfp 7 G-frame and Frenet boundaries
are refused: an independent cut shows their R(t, phi) spectrum decaying only
like n^-2 (3.6e-5 m at n = 32) because the strongly helical axis makes the
surface nearly tangent to the phi planes.

## Open

- VMEC profile types other than power_series are converted by the benchmark
  converter, not here.
- SPECTRE's relaxed volumes are replaced by one nested ideal state whose
  profile between interfaces is a model choice, not an error estimate.

## Sources

- GVEC documentation ("Coordinate conventions", parameter list) and source
  (`src/hmap`, `src/vmec/vmec.F90`, `python/gvec/util/geometry.py`),
  [itpplasma/gvec](https://github.com/itpplasma/gvec) 1ca0c4b30404b2f1cc2b76c918d8fb0bb95b7ef1,
  MIT licence.
- DESC `desc/input_reader.py`,
  [itpplasma/DESC](https://github.com/itpplasma/DESC) de3529ef75693c49a04ba6aecf77b96528889bc6,
  MIT licence.
- SPECTRE `spectre/converters/vmec_helpers.py` and `tests/compare` inputs,
  [krystophny/spectre](https://gitlab.com/krystophny/spectre) 0d8e45c4ed5ab807ec78b18222ce736b5084e2bd,
  MIT licence.
- S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983), the DSHAPE
  case.
