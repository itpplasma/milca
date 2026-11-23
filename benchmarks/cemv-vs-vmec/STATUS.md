# P1 status (stopped 2026-11-23, work in progress)

P1 is the head-to-head comparison of CEMV, the Lagrangian map and VMEC++ (see `docs/equilibrium-roadmap.md`). It was stopped at the user's request before anything landed. Nothing here is on main yet.

## Branches

| Repo | Branch | State |
|---|---|---|
| FortFEM | `p1/wall-corner-tangent` | Commit c53c754e. The harmonic-disc mesh wall is now one periodic spline, so it is smooth at the four block corners. Its own tests pass. Not landed |
| MILCA | `eq/p1-lag` | Lagrangian map: force balance J×B − ∇p, measured the same way as for CEMV (`milca_map_3d_force`); benchmark export of `nested_equilibrium`; perturbed starts for the basin study (`MILCA_NESTED_START_SHIFT`, `_ANGLE`). The energy in summary.json now includes the pressure term. Not landed |
| MILCA | `eq/p1` | Optional outputs of `flux_map3d_force_point` (label, ∇(B²/2), volume element); band and corner diagnostics in `test_flux_map_qa`; run scripts in `benchmarks/cemv-vs-vmec/`. Not landed |

## P1.0: the QA wall layer, s > 0.9

- **Finding (from the run logs, not yet written into flux-map.md):** 72–85 % of the squared force-balance residual in the wall layer sits in the cells next to the four block corners of the section mesh.
- **Cause:** the mesh wall was interpolated block by block, so it had a kink at each corner, up to 11.6° on a coarse, strongly shaped section. The FortFEM commit above removes the kink.
- **Not yet measured:** whether the fix makes the wall-layer residual on QA converge.

## Lagrangian-map force balance

The check on the circular tokamak gives a residual over s > 0.1 of 5.7e-3, 2.3e-3 and 1.3e-3 at 16, 24 and 32 boundary nodes (degree 3, angle_mpol 6). With angle_mpol 2 it stalls at 2e-2.

## Next steps

1. Land the FortFEM fix, then pin it in MILCA.
2. Rerun QA at 3:4:21 and 3:8:21 and classify the wall-corner hypothesis.
3. Run H2 (QA) completely, then H1 and H3, as specified in the roadmap's P1.

Scratch logs, which may not survive: `/private/tmp/claude-502/-Users-markl-m/e0271dd0-2407-4315-80ea-4ff4d4d3436e/scratchpad/p1/`.

## Lagrangian-map column (a sub-agent finished after the stop)

Branch `eq/p1-lag`, logs in `scratchpad/p1/lag/`. The force balance is computed with the same generated leaf as for CEMV. B agrees with `map_3d_sample` to 1e-14.

**H2, QA with prescribed ι** (map_ntor 5, angle_mpol 4; every run converged with one sign of the Jacobian):

| deg | nodes | unknowns | Newton | 8 thr (s) | 1 thr (s) | peak mem | energy vs VMEC++ | η_B, s > 0.1 | η_B, s > 0.9 |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 16 | 5,532 | 17 | 4.1 | 13.6 | 0.43 GB | +7.2e-6 | 6.2e-2 | 0.15 |
| 3 | 24 | 13,584 | 50 | 21.6 | 92.9 | 0.97 GB | +1.4e-6 | 2.0e-2 | 5.3e-2 |
| 3 | 32 | 24,012 | 62 | 49.8 | 203.6 | 1.6 GB | +4.4e-7 | 1.1e-2 | 3.0e-2 |
| 4 | 12 | 4,168 | 39 | 10.0 | 42.5 | 0.57 GB | +8.0e-8 | 1.4e-2 | 3.1e-2 |
| 4 | 16 | 9,976 | 126 | 70.4 | 330.5 | 1.0 GB | +4.1e-9 | 5.5e-3 | 1.3e-2 |
| 4 | 24 | 24,408 | 234 | 283 | – | 2.1 GB | −1.1e-8 | 3.7e-3 | 8.6e-3 |

- **Geometry:** axis 0.14–0.16 mm, surfaces 0.03–0.06 mm rms in every row. This is at the uncertainty of the ns-101 reference (about 0.13 mm), so a finer VMEC++ reference is needed to rank the runs.
- **Interior:** η_B for s < 0.5 is 1.2e-3 to 2e-3, about the level of VMEC++ (1.8e-3).
- **Cost:** Newton steps grow steeply with resolution.

**H3, CTH-like deck** (fixed boundary, cylindrical chart; all converged, no fold):
- Axis 0.14–0.16 mm.
- Energy −3.5e-7 to −6.8e-7 against ∫(B²/2 − μ0 p).
- η_B does not converge (9.0e-2 to 6.6e-2, order about 0.4). Resonant surfaces are suspected but not isolated.
- The curve chart is not wired.

**Basins on H2** (degree 3, 16 nodes):
- Axis shifts of ±1 to 3 cm reach the same state; ±6 cm and beyond fold.
- Boundary-angle perturbations of up to 1 rad reach the same state.

**Not run:** H1; implicit sensitivities (no machinery exists); a registered test for `milca_map_3d_force`.
