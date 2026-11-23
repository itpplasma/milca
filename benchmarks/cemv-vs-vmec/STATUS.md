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
