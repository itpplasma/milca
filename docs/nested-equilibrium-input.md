# Nested equilibrium from a boundary and profiles

`nested_equilibrium --config FILE` solves the fixed-boundary nested 3D map
equilibrium (the D8 boundary/profile lane, [ideal-seed-3d](ideal-seed-3d.md))
for a boundary, a toroidal flux and prescribed profiles, or with
`free_boundary` the free-boundary equilibrium in the field of coil groups
([free boundary](#free-boundary)), and writes the solved
state's global quantities to `output/nested_equilibrium.nc`. It is the native
entry for cross-code comparison with VMEC-class codes. `milca_import` writes
this namelist from GVEC, DESC, SPECTRE and VMEC inputs
([equilibrium-import](equilibrium-import.md)).

## Input: namelist `nested_input`

| Name | Meaning | Default |
|---|---|---|
| `nfp` | field periods | 1 |
| `mpol`, `ntor` | largest poloidal and toroidal boundary mode, at most 64 and 256 | 1, 0 |
| `map_ntor` | toroidal modes of the interior map, `ntor` to 256; -1 chooses them by `fourier_tolerance` | -1 |
| `fourier_tolerance` | largest accepted `fourier_error` with `map_ntor = -1` | 1e-4 |
| `rbc(n,m)`, `rbs(n,m)`, `zbc(n,m)`, `zbs(n,m)` | boundary R, Z = sum C cos(m theta - n nfp phi) + S sin(m theta - n nfp phi), indexed (n, m) as in VMEC | 0 |
| `phiedge` | total toroidal flux (Wb), positive | 1 |
| `iota(k)`, `pressure(k)` | coefficients of s^k, s the normalized toroidal flux; pressure in Pa | 0 |
| `iota_type`, `pressure_type`, `current_type` | `power_series` (the coefficients above and `current(k)`) or `hermite` | `power_series` |
| `n_knots`, `profile_knots(:)` | knots 0 = s_1 < ... < s_K = 1 of the `hermite` profiles, K at most 1024 | 0 |
| `iota_values(:)`, `iota_slopes(:)`, `pressure_values(:)`, `pressure_slopes(:)`, `current_values(:)`, `current_slopes(:)` | values and s-derivatives at the knots of a piecewise cubic Hermite profile; pressure in Pa, current in A | 0 |
| `boundary_nodes`, `degree` | disc mesh of the final level; `boundary_nodes = -1` chooses it by `mesh_tolerance` | -1, 3 |
| `mesh_tolerance` | largest accepted `mesh_error` with `boundary_nodes = -1` | 1e-4 |
| `ladder_nodes(:)`, `ladder_degree(:)` | coarser levels solved first; each corrects the next seed (`ideal_map_prolong_correction`) | none |
| `angle_mpol` | poloidal modes of the boundary field-line angle; -1 pins the boundary angle | -1 |
| `tolerance`, `max_iterations` | free reduced-gradient reduction and Newton limit per level | 1e-8, 200 |
| `stationary` | solve for a stationary point at fixed boundary and prescribed rotational transform | false |
| `prescribe_current` | match the enclosed toroidal current `current(k)` (A, coefficients of s^k) instead of keeping `iota` | false |
| `max_current_iterations`, `current_tolerance` | Newton steps on iota and the largest iota increment of the last step; with automatic resolution at most a tenth of the active `fourier_tolerance`, `mesh_tolerance` and `iota_tolerance` | 30, 1e-5 |
| `iota_nodes`, `iota_tolerance` | Chebyshev-Gauss surfaces of the matched iota polynomial, 2 to 32, or -1 to raise them from 6 under error control | -1, 1e-4 |
| `free_boundary` | the boundary is the initial guess of a free-boundary solve | false |
| `mgrid_file`, `extcur(k)` | MAKEGRID mgrid file of the coil groups and the current of group k (VMEC's EXTCUR) | none, 0 |
| `mgrid_degree` | interpolation degree of the mgrid in (R, Z), 1, 3 or 5 | 5 |
| `coupling_margin` | clearance of the coupling surface around the guess, in units of the minor radius | 0.5 |
| `vacuum` | `auto`: `bem` for an axisymmetric deck, otherwise `shell`; `bem`: sheet current on the plasma boundary, no vacuum mesh (axisymmetric decks only); `shell`: vacuum shell between the plasma and the coupling surface | `auto` |

A `hermite` profile is C1 and cubic on each knot interval; the iota integral
from the axis is exact on it. Any other profile form (VMEC's splines, sums of
arctangents, pedestals) is passed as values and slopes at enough knots; the
interpolation error of that representation belongs to the converter.

With a prescribed current, `iota` is the starting profile. Iota is the
polynomial through its values c_j at n = `iota_nodes` Chebyshev-Gauss
surfaces s_j = (1 - cos((2j - 1) pi/(2n)))/2, which exclude the axis, stored
by values and slopes on 129 uniform Hermite knots. The current is matched at
the s_j:

    F_j(c) = I(s_j; x(c), c) - I_T(s_j) = 0,                                  (1)

with x(c) the map in equilibrium and I(s) = (1/mu0) oint B . dl along the
image of the reference circle s at phi = 0. The circle is cut at its
crossings with the straight cell edges, which a quadratic gives exactly, and
every arc lies in one cell, where Gauss-Legendre with 8 points converges
spectrally (`milca_map_3d_loop`). The arcs depend only on the reference mesh
and are built once per level.

On a fixed map, B on the surface s depends on iota(s) alone and linearly, so
I(s_j) = A_j c_j + C_j. Newton on (1) uses

    dF/dc = diag(A) + (dI/dx) dx/dc,    dx/dc = -H^-1 dg/dc,                  (2)

with g the energy gradient and H its Hessian (`ideal_map_sensitivity`, one
factorisation for all columns). A comes exactly from one shift of every c_j
by 1. The energy is quadratic in c on a fixed map (third and fourth
differences are at roundoff), so dg/dc_k from the unit steps +-e_k is exact.
(dI/dx) dx_k is a central difference of the loop integrals along dx_k, whose
O(h^2) error changes only the rate, not the matched current. The iota
increment of the last Newton step is `current_error`. Intermediate maps are
solved to max(tolerance, min(1e-4, 1e-3 * increment)) against the first
seed's gradient, and the last step is repeated on the fully solved map. A
step may reverse the sign of iota on every surface (a current-free
stellarator started from the other handedness). An iota that changes sign
inside [0, 1] is not admitted by the frozen map
([ideal-seed-3d](ideal-seed-3d.md)), so a step that would produce one is
halved until iota keeps one sign. An iterate whose |iota| comes within
`current_tolerance` of zero has no resolved sign and ends the solve with
`current_status = "iota vanishes"` and a `failure` reason. On an axisymmetric boundary Ampere's
law on the surface s gives mu0 I(s) = psi'(s) oint |grad s|/R dl with a
positive integral of order s near the axis, so iota has the zeros of
I_T(s)/s. A prescribed current that vanishes on (0, 1], or has dI_T/ds = 0
on the axis, is therefore reported as a `failure` before the seed.

The current is matched pointwise because the energy weighs iota with s near
the axis: for a large-aspect-ratio circle the poloidal field energy is
proportional to int s iota^2 ds. dW/dc_k = mu0 phiedge int I l_k ds (l_k the iota
of c = e_k) holds on the discrete map (to 4e-5 against the loop integral on
the circular case), but a Galerkin match of these weak currents fixes iota
only in that weighted norm. On the circular case its iota on axis moved
between 0.59 and 0.72 from mesh to mesh, and raising the polynomial degree on
a fixed mesh changed the observables by 2e-2 to 3e-1 instead of less.

With `iota_nodes = -1` the nodes start at 6 and are raised by 2 until the
change of the observables, `iota_error`, is at most `iota_tolerance`.

The resolution is chosen by the change of the observables between levels:
wb (relative), raxis/aminor_p and iota on axis and edge. With
`boundary_nodes = -1` the mesh runs through 16, 24, 32, 48, 64, 96 and 128
boundary nodes at fixed `degree` until that change, `mesh_error`, is at most
`mesh_tolerance`. The interior map needs more toroidal harmonics than the
boundary: the equilibrium inside a boundary with `ntor` modes is not limited
to them. With `map_ntor = -1` the map's modes are raised by 2 until the
change, `fourier_error`, is at most `fourier_tolerance`. Each level is seeded
with the correction of the solved previous one
(`ideal_map_prolong_correction`), and the current is matched to a tenth of the
tolerances so that the changes measure truncation. With both automatic, the
mesh ladder runs first at the boundary's modes, then the toroidal ladder and
the iota-node ladder, and the active ladders alternate until each is met by
its first refinement. Axisymmetric data
(`ntor = 0`) keep `map_ntor = 0`.

With `stationary = .true.`, the Newton step solves \(H d = -g\), where
\(H\) is the undamped Hessian and \(g\) is the free energy gradient. The
line search reduces \(\|H^{-1} g_{\mathrm{trial}}\|_2\) with \(H\) held
at the current iterate. This can reach a saddle of the energy. Each
refinement transfers the solved geometry and field-line angle and retains
the target level's boundary. The stopping tolerance is a fraction of the
free-gradient norm of the uncorrected seed at that level.

The cost grows steeply with the interior modes. For the rotating ellipse
nfp = 5, R0 = 5, 16 boundary nodes, at 8 threads:

| map_ntor | 2 | 4 | 8 | 12 | 16 |
|---|---|---|---|---|---|
| wall time (s) | 2.3 | 4.9 | 36 | 109 | 114 |
| peak memory (GB) | 0.11 | 0.31 | 1.0 | 2.1 | 3.6 |

Memory grows like map_ntor^2, about 14 GB at 32 and 58 GB at 64 modes on this
mesh. The energy sees the map only at its
quadrature points, so the continued start of a level can reach a stationary
point that folds between them; the volume-integral points detect the fold,
and that level is solved again from its Winslow seed. A level whose map still
folds does not resolve the map: a fold on a plane of the energy's toroidal
quadrature lies between in-plane points and raises the mesh, a fold only
between those planes raises the toroidal modes (the interior of a boundary
with a large axis excursion needs more harmonics than the boundary). With
that resolution fixed or at its limit the run ends with a `failure`.

The seed of every toroidal plane is the Winslow map of its cross-section
(`milca_winslow_seed`, [ideal-seed-3d](ideal-seed-3d.md)): the interior that
minimises the Dirichlet energy of the inverse map, which is injective for any
simply connected section and needs no initial axis. A boundary without an
injective interior map, such as a curve traversed twice, ends the run with
`failure = "the seed is not admitted ..."`.

The boundary must run counterclockwise in (R, Z) for increasing theta. A
boundary given in an arbitrary poloidal angle (every VMEC input) needs
`angle_mpol >= 0`: with the boundary angle pinned, the field lines on the wall
are tied to that angle and the solve converges to a different, constrained
state (on the circular test case the axis is then 8-20 mm off).

## Free boundary

With `free_boundary` the plasma map of the
[3D free boundary](free-boundary-nested-map-3d.md#10-free-boundary-decks-of-the-nested-equilibrium-app)
is solved on the cylindrical chart with the per-n exterior. The vacuum field
is sum_k `extcur(k)` b_k with b_k group k of the mgrid file at unit current
(`mgrid_mode` 'S'; a raw grid 'R' is divided by `raw_coil_cur`). The grid is
interpolated in (R, Z) by the tensor Lagrange polynomial of degree
`mgrid_degree` (error O(h^(degree + 1))), and across its planes by the exact
discrete Fourier series. The toroidal flux is positive: a deck with a negative
toroidal field reverses the signs of all `extcur` and of the current.

An axisymmetric deck takes `vacuum = 'bem'` by default: the vacuum is the
coil field plus a sheet current on the plasma boundary, whose energy is a
boundary integral of the Grad-Shafranov single layer
(theory/derivations/free-boundary-plasma-boundary-bem.md), solved by the
inertia-controlled Newton. The coil flux is fitted by ring sources to the
grid nodes inside Gamma_c, which only selects those nodes; the misfit is
printed. An up-down symmetric deck holds the vertical position. The result
records the vacuum used in the attribute `vacuum`.

The coupling surface Gamma_c is an axisymmetric curve that keeps
`coupling_margin` times the minor radius from every cross-section of the guess
and stays inside the grid and off the coils that lie inside it (grid nodes
where the curl of the field is not small). Where they are closer, the
clearance shrinks to 70 % of the room. The reference disc of the frozen field
is the ellipse inscribed in the mean bounding box of the sections.

A prescribed current fixes the circulation c = mu0 I(1) of the vacuum field.
Before the first free-boundary solve the interior of the guess is brought to
fixed-boundary equilibrium with the current matched there, so that the
plasma's edge current equals c from the start; the free-boundary current is
then matched as in (1) and (2), with dx/dc from the slide-gauged KKT matrix of
the free boundary (`ideal_fb_sensitivity`). With prescribed iota, c is set to
mu0 times the solved edge current by secant steps. Further mesh and toroidal
levels start from the solved map interpolated to the new mesh. `angle_mpol`
is not used: the boundary nodes move freely.

The free-boundary output adds `coupling_margin` (the smallest clearance used),
`coil_closure` (a bound of the n = 0 closure of the coil traces, zero when no
coil lies inside Gamma_c), `circulation` and `vacuum_f` (R B_phi of the coil
field); `force_residual` is the free-boundary residual sqrt(merit/(sigma V)).

## Output

`volume_p`, `aminor_p` (sqrt of the mean cross-section area / pi), `rmajor_p`
(volume / (2 pi mean area)), `aspect`, `raxis_cc` (toroidal mean of the
O-point's major radius), `iotaf_axis`, `iotaf_edge`, `pressure_axis` (Pa),
`wb` (int B^2/2 dV / (2 pi)^2 with mu0 = 1, VMEC's normalisation), `betatotal`,
`force_residual` (free reduced-gradient norm relative to the final level's
uncorrected seed), `tolerance`, `iterations`, `ctor` (enclosed toroidal
current at the edge, (1/mu0) of B.dl along the discrete wall at phi = 0),
with a prescribed current also
`current_error`, `current_iterations`, `iota_nodes` and the attribute
`current_status`, with automatic iota nodes `iota_error`, the
map's toroidal modes `map_ntor` and `boundary_nodes`, with automatic
resolution `fourier_error` and `mesh_error`, the final solve's wall seconds
(`energy_seconds`, `hessian_seconds`, `factor_seconds`) and counts
(`factorizations`, `lu_fallbacks`, `analyses`), and the global attribute
`converged`. A seed that is not admitted, a folded map met by a readback or a
failed correction ends the run with exit status 0, `converged = "false"` and
the global attribute `failure` giving the reason; invalid input stops with an
error. Volume integrals use collapsed Gauss-Legendre on each triangle and
the trapezoidal rule in phi; the mean area is int dV/R / (2 pi).

These are candidate values: `mesh_error` and `fourier_error` estimate the
truncation from two levels each; no bound is reported.

## Checks

`equilibrium_nested_circular-tokamak` solves
`examples/nested-equilibrium/circular-tokamak.nml` (R0 = 3 m, a = 1 m,
iota = 0.4 + 0.2 s, p = 0) and requires the exact volume 2 pi^2 R0 a^2 to 1e-4,
the axis within 1.5 mm of 3.0366 m and `wb` within 2e-5 of 0.0750949; PARVMEC
(ns = 101, mpol = 8) gives 3.03634 m and 0.0750949, GVEC gives 3.03683 m.
`equilibrium_nested_circular-tokamak-current` prescribes I(s) = 340876 A * s on
the same boundary, starting from iota = 0.3; PARVMEC with CURTOR = -340876 A and
AC = 1 (its current sign convention is opposite) gives raxis 3.04266 m and wb
0.0755050, MILCA 3.04281 m and 0.0755055 at 6 iota nodes (0.4 s); with
`boundary_nodes = -1` the ladders stop at 24 boundary nodes and 8 iota nodes
with 3.04281 m and 0.0755051 (1.1 s).
`equilibrium_nested_qi-nfp2-current-free` (label slow) solves the boundary of
`input.nfp2_QI` of the VMEX examples at zero net toroidal current. The matched
iota is negative from an iota = 0.3 start; with 16 boundary nodes and 8
iota nodes the toroidal levels 4 to 12 give fourier_error 0.19, 1.9e-2,
6.4e-3, 4.4e-4 and 5.1e-5 (188 s at 8 threads). VMEC2000 at NTOR 6 gives
raxis 0.98461 m, wb 2.276210e-3 and edge iota -0.44174; MILCA 0.98450 m,
2.276223e-3 and -0.44229, a difference within VMEC2000's own truncation at six
toroidal modes as estimated by the levels above. The iota nodes are fixed
because on this mesh raising them from 10 to 30 changes the observables by
2e-4 to 6e-4 without decreasing, the mesh's limit on pointwise iota; a ladder
whose change does not shrink stops there.
