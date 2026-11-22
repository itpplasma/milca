# Flux-map prototypes (throwaway)

Python prototypes behind [docs/flux-map.md](../../docs/flux-map.md): `fluxmap.py` (axisymmetric Solov'ev), `fluxmap3d.py` and `fluxmap3d_diag.py` (Landreman–Paul QA, C0 Lagrange hexahedra).
They need a virtual environment with numpy, scipy and jax and are not production code; the solver is `milca_flux_map` and `milca_flux_map_3d`.
