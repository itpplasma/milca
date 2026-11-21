# Lessons: free-boundary equilibrium against VMEC

Physics area: axisymmetric free-boundary equilibrium with external coils given
as an mgrid (MAKEGRID's tabulated coil-group fields), solved by MILCA's
Eulerian Grad-Shafranov solver (`gs_free`). The solver can take VMEC's
constraint set: pressure p(s), enclosed toroidal current I(s) in the
normalised toroidal flux s = Phi/Phi_edge, and the edge toroidal flux
Phi_edge.

Reference: VMEC2000 on the DIII-D deck of its example set (`input.DIII-D`,
`lasym = T`, `mgrid_d3d_ef.nc` with 126 x 201 nodes). Campaign of November
2026.

Outcome:
- The deck as distributed is not a well-posed reference (cause 3.1).
- At 0.95 Phi_edge MILCA agrees with VMEC2000 to:
  - 1e-5 in I_p;
  - 1.1e-4 in volume and 0.11 mm in the axis against VMEC2000 extrapolated
    in its mgrid and ns resolution (cause 4.1);
  - 4e-3 in q(0.99).

Each cause below says whether the reference can be fixed (implementation) or
is limited in principle (representation). MILCA reproduces neither kind.

## 3. Different model or convention in the reference

### 3.1 VMEC's boundary at or beyond the separatrix (limit in principle)

The VMEC2000 DIII-D deck puts its boundary on the near-double-null separatrix
of its own coil field. With VMEC's coil currents, p'(psi_N) and FF'(psi_N):
- MILCA's separatrix encloses 3.654 Wb against PHIEDGE = 3.661 Wb;
- the X-points lie at (1.29, +-1.12) m, just outside VMEC's top and bottom
  corners.

A Fourier boundary cannot form an X-point, so VMEC converges to a boundary
whose corners are rounded and which encloses more flux than closed surfaces
allow.

Evidence, from the mpol scan at ns 128 (`scluster:~/bench/vmpol`):

| PHIEDGE | boundary change, mpol 12 to 16 / 16 to 20 | q(0.99) at mpol 12 / 16 / 20 | corner curvature (1/m) at mpol 12 / 16 / 20 |
|---|---|---|---|
| full | 19 mm / 4.4 mm | 5.00 / 5.30 / 5.40 | 9.1 / 13.9 / 15.2 |
| 0.95 | 3.8 mm / 0.6 mm | 3.952 / 3.973 / 3.975 | 6.6 / 7.7 / 7.8 |

At the full flux the corners keep sharpening with mpol. mpol 24 did not
converge at ns 128 in either case. The exact closed-surface limit of the flux
is open: within about 1 % of PHIEDGE.

Rule: a free-boundary reference must have its boundary inside the separatrix
of its coil field. Check this with MILCA's X-point search before comparing.
The benchmark case is the same deck at 0.95 PHIEDGE.

### 3.2 Edge pressure against the vacuum (VMEC: balanced through F)

VMEC adds the edge pressure to the vacuum pressure balance
(`funct3d.f`: `pgcon = bsqvac + presf_ns`). Lowering p(1) from 1127 Pa to 0
(AM0 shifted) changes:
- the boundary by 1e-8 m and R_axis by 6e-9 relative;
- beta by 5 %.

The balance is taken up by B_phi, through the plasma's F, not by the shape.
MILCA's Grad-Shafranov equation sees only p', so a pressure jump at Gamma_p
carries no current sheet. Physically, [[p + B^2/2 mu0]] = 0 requires a sheet,
poloidal (F jump), toroidal (B_pol jump) or both. Open item for MILCA: the
sheet as an interface term in the solver, not a choice made by convention.

## 4. Resolution of MILCA or of the reference

### 4.1 VMEC2000 interpolates the mgrid bilinearly (fixable)

`NESTOR_vacuum/becoil.f` interpolates B_R, B_phi and B_Z of the mgrid
bilinearly. The error is O(h^2). On the plasma boundary of the DIII-D grid it
is 2.2e-4 of |B_pol| (rms 6e-5), and the equilibrium responds strongly.

Effect, measured by refining the same mgrid 2x with local 6-point Lagrange
interpolation (exact at the original nodes,
`scluster:~/bench/mgrid_d3d_ef_x2.nc`), at ns 128, mpol 16, 0.95 PHIEDGE:

| | V (m^3) | R_axis (m) | q(0.5) |
|---|---|---|---|
| VMEC2000, mgrid as given | 18.569152 | 1.714226 | 1.775731 |
| VMEC2000, mgrid refined 2x | 18.560785 | 1.713856 | 1.776178 |
| MILCA | 18.5544 | 1.71350 | 1.77709 |

The original grid costs VMEC 4.5e-4 in volume and 0.37 mm in the axis. With
grids refined 2x and 4x and ns 128 and 256, both errors extrapolated (h^2 in
the mgrid, ns^-2 in ns, mpol 16), VMEC2000 converges to V = 18.5565 m^3,
R_axis = 1.71361 m, q(0.5) = 1.77635 and q(0.9) = 3.35261. MILCA (96 boundary
nodes, P2, VMEC's input profiles) differs from that by -1.1e-4 in V, -0.11 mm
in R_axis and +4e-4 and +6e-4 in q, inside its own mesh spread of 1e-4 in V.
The 0.65 mm axis difference to VMEC2000 on the mgrid as given is almost
entirely VMEC's interpolation error. A
global spline (scipy `RegularGridInterpolator`, quintic) is not a valid
refinement: the near-singular nodes next to coil filaments make it ring by
1e-4 across the grid.

MILCA uses the coil fields differently. It integrates the poloidal flux
psi_g per group from B by least squares over coil-free grid edges, with
fourth-order edge integrals; nodes next to a coil are excluded. It then
interpolates psi_g with degree-5 tensor Lagrange polynomials (`gs_coils`,
kind `grid`). The flux error is below 1e-5 in the plasma region, checked by
halving the grid. Rule: never interpolate a tabulated vacuum field at less
than fourth order.

### 4.2 VMEC's radial discretisation of the current (resolution)

VMEC's `ctor` differs from the input CURTOR:
- 1.9e-4 at ns 128;
- 4.6e-5 at ns 256.

R_axis moves by 0.08 mm from ns 128 to ns 256. MILCA meets the input current
to 3e-6. Rule: compare against the input constraint, not against VMEC's
output of it.

### 4.3 MILCA's q readback on uniform rho levels (MILCA, fixed)

The readback read q on about 45 levels uniform in rho = sqrt(psi_N). Near an
edge where q grows quickly, the cubic integral Phi = int 2 pi q dpsi is then
too coarse:
- s_tor came out 1 % low near the edge;
- q(s_tor = 0.95) was 12 % high on the DIII-D fixed boundary, which looked
  like a 4-7 % solver error.

With levels packed towards the edge (rho_j = sin(pi x_j/2)), `gs_phase1` and
`gs_free` agree with VMEC ns 512 at q(0.5, 0.9, 0.95, 0.99) to 5e-5, 8e-5,
0.14 % and 0.36 %. That is inside VMEC's own mpol and ns spread.

Check: `grad_shafranov_phase1_edge_readback_app` requires the contour and
region toroidal fluxes to agree to 8e-6 on a steep-edge D shape. Packed
levels give 2.9e-6, uniform ones 2.0e-5.

## 6. Inconsistent input data

### 6.1 Profiles rebuilt from the wout file (comparison procedure)

I(s) rebuilt from `buco` (half grid, extrapolated to the edge) and p(s) from
`presf` agree with VMEC's input polynomials inside, but not at the edge:
- dI/ds(1) = 5.46e4 against 2.91e4 A;
- dp/ds(1) = 1.39e4 against 2.01e4 Pa.

The rebuilt table also kept two nodes 1e-16 apart. Cubic Lagrange
interpolation through them blew up FF' at one knot (the step from -0.456 to
43.4), and the outer iteration showed a false floor of 1e-3.

Rules:
- Take a reference's profiles from its input (AM, AC with CURTOR for VMEC),
  not from its output.
- `gs_read_current_target` rejects nodes closer than 1e-9.

## 7. Implementation and toolchain defects (MILCA)

- **Picard on the profiles is unstable.** It amplifies a deviation 2-7 fold
  per pass, even when started at the fixed point, because the flux depth and
  I_p respond strongly to a law in physical units. Anderson acceleration
  (memory 6) converges in 6-7 solves.
- **Double null.** The X-point mode followed one hint. With two nearly equal
  X-points the bounding one changes between iterations. `two_xpoints` takes
  the saddle with the larger axis_sign*psi.
- **Thread oversubscription.** `gs_free` with 32 OpenMP threads under a
  workstation load of 26 stalls completely. Set OMP_NUM_THREADS to the free
  cores.

## Rules from this campaign

1. **Well-posed references only.** Before a free-boundary comparison, check
   that the reference boundary lies inside its separatrix.
2. **Reference convergence includes its auxiliary data.** Converge the
   reference in ns and mpol, and also in its coil-field interpolation.
3. **Compare against input constraints, not output.** Use the reference's
   input profiles and currents, not quantities rebuilt from its output.
4. **Do not tune MILCA towards an unconverged reference.** A remaining
   difference is attributed only after the reference's own errors are
   measured.
