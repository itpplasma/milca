# CEMV or VMEC direction: which one suits which downstream code

Recorded 2026-11-23. This is guidance for planning. It has not been measured.

## Two directions

- **CEMV, the flux map.** The unknowns are the labels as functions of position, u(x) = (s, θ*). This is `milca_flux_map_3d`.
- **VMEC direction, the Lagrangian map.** The unknowns are positions as functions of labels, X(s, θ, φ). This is the Lagrangian map, and also VMEC, DESC and GVEC.

The exact solutions are inverse to each other. Each direction gives one question directly and needs an inversion for the other, usually a Newton solve per point:

- "Where is the point with these labels?" is direct in the VMEC direction.
- "What are the labels at this point in space?" is direct in CEMV.

## Which direction suits which code

| Area | Codes | What they ask of the equilibrium | Preferable | Why |
|---|---|---|---|---|
| Guiding-centre orbits, confined region | SIMPLE, VENUS-LEVIS | B and its derivatives in flux or canonical coordinates; exact surfaces | VMEC direction | Integration in (s, θ, φ) keeps surfaces exact and supports symplectic schemes |
| Guiding-centre orbits with losses, walls, fast-ion loads | ASCOT5, BEAMS3D, KORC | B on a cylindrical grid, the wall, the region outside the last surface | CEMV | The orbit lives in x, so s(x) is direct; the wall and the vacuum are native. These codes resample VMEC onto (R, φ, Z) today |
| Ideal linear stability, ballooning, Alfvén gaps | CAS3D, COBRAVMEC, TERPSICHORE, STELLGAP/AE3D, FAR3d | Flux or Boozer coordinates; metric on surfaces; resonant surfaces | VMEC direction | The modes are written in (s, θ, φ) |
| Nonlinear or extended MHD, hybrid | M3D-C1 (stellarator), JOREK (stellarator), NIMROD, MEGA | Initial state on a fixed real-space mesh, often beyond the last surface | CEMV | Same mesh and representation, no mapping; non-nested regions possible later |
| Gyrokinetic turbulence, core | GENE (GIST), stella, GX, EUTERPE, GENE-3D | Field-aligned geometry | VMEC direction | Field-aligned coordinates are built from flux coordinates |
| Gyrokinetic turbulence, edge or X-point | XGC-S and similar | Real-space mesh across the separatrix | CEMV | Nested surfaces end there |
| Neoclassical theory | DKES, NEO-2, SFINCS, MONKES, KNOSOS, NEO (ε_eff) | Boozer spectrum on each surface | VMEC direction | Boozer transformation from X (BOOZ_XFORM) |
| 1D transport, profile evolution | Trinity3D, TANGO, NTSS, ASTRA (stellarator) | V′(s), surface averages, metric terms | VMEC direction, slightly | Surface averages are direct in X; CEMV gets them by integrating over a level set (co-area formula) |
| Edge and scrape-off layer, divertor | EMC3-EIRENE | Field lines in real space, islands, stochastic layer, targets | CEMV, or real-space codes such as HINT | Real space and open field lines, which the VMEC direction cannot represent |
| Heating, ray tracing | TRAVIS, ECRH and ICRH ray tracers | n and T at points x along rays, and ∇s | CEMV | n(s(x)) is direct; with the VMEC direction the code inverts per point |
| Profile reconstruction (Thomson, ECE, interferometry) | Minerva, V3FIT, TIAGO | s along lines of sight; profiles as f(s(x)); derivatives with respect to the parameters | CEMV | Lines of sight live in x; s(x) and its derivatives are direct, with no Newton inversion inside the fit |
| Magnetic reconstruction (coils, loops) | V3FIT, Minerva, STELLOPT reconstruction | Plasma current in space; virtual casing | CEMV or free-boundary VMEC direction | Plasma current in space plus Biot–Savart is natural in x; free boundary with virtual casing is the established route |
| Optimisation, physics targets | SIMSOPT, STELLOPT, DESC | Quasi-symmetry residual, ι, magnetic well, ε_eff | VMEC direction | The targets are defined in flux or Boozer coordinates |
| Optimisation, coils and engineering | FOCUS, REGCOIL, SIMSOPT stage 2 | Plasma boundary, B·n on it, coil–plasma distance | CEMV | Everything is in real space, and the free boundary is natural |

The equilibrium solvers themselves split as follows:

- **Inverse, VMEC direction:** VMEC and VMEC++, DESC, GVEC, NSTAB.
- **Real space or non-nested:** HINT, PIES, SIESTA.
- **Multi-region, inverse within each volume:** SPEC.

## Summary

- **VMEC direction:** suits work inside the nested core in flux coordinates, namely stability, turbulence, neoclassical theory, transport and the physics targets of optimisation.
- **CEMV:** suits codes that work at points in space, namely walls, the edge, coils, lines of sight, heating rays, losses and nonlinear MHD.

Reconstruction for experiments and edge work lean towards CEMV, while optimisation and core theory lean towards the VMEC direction. Keeping both, tied together by the inverse check (P2 in [equilibrium-roadmap.md](equilibrium-roadmap.md)), gives every consumer its natural direction without an inversion of its own.

The entries rest on the published descriptions of these codes, not on a check of their current input options.
