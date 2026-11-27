# MILCA roadmap — active work first

Status: 2026-11-20. **This is the only active priority list.** Issues not listed here are
backlog, not instructions to implement. Physics specifications and long-range
scope remain in [DESIGN.md](DESIGN.md), [the capability catalogue](docs/capability-catalogue.md),
[the model graph](docs/model-graph.md), and the
[previous detailed milestone plan](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md).

## Mission and division of work

MILCA is the group's fast-moving, independently verifiable **physics and numerics
laboratory**, not a requirement that everyone use one executable. Established
specialized codes (GORILLA, SIMPLE, NEO-2, MARS, BEAMS3D, VMEC++, etc.) are
independent products, competitors, consumers and scientific comparators.
Improve them through focused, independently evidenced upstream PRs; their
maintainers decide adoption and review. Do not transfer MILCA implementations
into comparator codes just to obtain agreement.

The scientific target is **time and cost to a named observable at matched
physical model, accuracy and uncertainty**, not amount of code, number of
features or raw runtime alone. Keep the group's human effort on physics.

## Current work (ordered by priority)

| Priority | Work package and next gate | Tracked work |
|---|---|---|
| **P0** | **Restore trustworthy gates:** finish the adaptive NTV mesh error/size-limit defects (refinement and band joins are capped since #372 and #382; the initial mesh of an energy can still exceed the limit); investigate open nightly and equilibrium timeout failures. Do not accept an accuracy or resource claim when its gate was bypassed. | [#363](https://github.com/itpplasma/milca/issues/363), [#362](https://github.com/itpplasma/milca/issues/362), [#259](https://github.com/itpplasma/milca/issues/259), [#145](https://github.com/itpplasma/milca/issues/145) |
| **P1a** | **NTV/neoclassical physics:** converge device-case observables (data outside the repository) with explicit model/convention match, resolution, error and cost; focus on unresolved resonant-ion discrepancies. Finish the currently exercised collision/source work: the energy-coupled kernel with finite orbit width on the hp orbit-class mesh (stacked branches `ntv/a01` to `ntv/a05`) and the torque density in real space. | [#233](https://github.com/itpplasma/milca/issues/233), [#197](https://github.com/itpplasma/milca/issues/197), [#198](https://github.com/itpplasma/milca/issues/198), [#369](https://github.com/itpplasma/milca/issues/369) |
| **P1d** | **Drift-kinetic transport at NEO-2 feature parity, helical core first (Chris, 2026-11-19):** local 3D transport, bootstrap current and conductivity on general 3D surfaces without stellarator symmetry, the prescribed in-surface potential drive with local parallel current, helical-core geometry, quasineutral Phi_1, then the rest of the tracker in dependency order. Each slice is checked against an analytic limit, a published result or SFINCS. | [#198](https://github.com/itpplasma/milca/issues/198), [#194](https://github.com/itpplasma/milca/issues/194), [#220](https://github.com/itpplasma/milca/issues/220), [#221](https://github.com/itpplasma/milca/issues/221), [#222](https://github.com/itpplasma/milca/issues/222), [#224](https://github.com/itpplasma/milca/issues/224) |
| **P1b** | **3D equilibrium, nested and with islands, fixed and free boundary:** one residual and one state for every topology, in at most the time of the fastest comparator at the same accuracy: the VMEC implementations, GVEC and DESC on nested cases, SPECTRE and SIESTA on island cases. Gates E1 to E10 of the [P1b procedure](#p1b-procedure). Next, by the [equilibrium roadmap](docs/equilibrium-roadmap.md) (P0–P6): finish the verification of the 3D flux map, compare it head to head with the Lagrangian map (its inverse) against VMEC++, choose the ideal core from that, then couple it to models without imposed surfaces, first in axisymmetry against Grad–Shafranov. Every equilibrium provides its Jacobian and kinetic-energy mass to P1e. | [plan](docs/equilibrium-roadmap.md), [benchmark](https://github.com/itpplasma/benchmark_vmec), [free-boundary plan](docs/free-boundary-ale.md), [islands](docs/island-equilibrium.md) |
| **P1c** | **Drift-wave turbulence for the resonant-layer response:** spectra, correlation scales and transport coefficients at one rational surface ([drift-turbulence.md](docs/drift-turbulence.md)). The linear fluid roots (stage L), Hasegawa-Wakatani turbulence with its saturated-flux benchmark (stage HW) and the homogeneous two-dimensional step of the drift-Alfven model (stage DA) are in; next is stage DA beyond the homogeneous periodic box. | [#371](https://github.com/itpplasma/milca/pull/371), [#379](https://github.com/itpplasma/milca/pull/379), [#383](https://github.com/itpplasma/milca/pull/383) |
| **P1e** | **Global 3D MHD stability with and without islands:** a fixed- and free-boundary stability capability for 3D equilibria with nested surfaces, islands and stochastic regions that matches or exceeds the published capabilities of CAS3D, CASTOR3D, TERPSICHORE and DCON (unstable global modes, marginal points, stable spectra, wall and vacuum response), with certified errors, derivatives with respect to equilibrium and design parameters, and cost per observable at accepted accuracy as the figure of merit. States with islands use the Jacobian of the same equilibrium residual under a declared constraint class ([islands](docs/stability.md#stability-of-equilibria-with-islands)). Built under the procedure below; each rung starts when the P1b gate of its equilibrium has passed, as listed in [stages](docs/stability.md#stages). | [specification](docs/stability.md#global-3d-ideal-mhd-stability-p1e) |
| **P1f** | **Visco-resistive MHD on the resistive time scale:** how the q profile of a driven tokamak core evolves and which stationary states it reaches, on the ladder of symmetry reductions of the one visco-resistive residual. Each rung is run as time evolution, Newton steady state and spectrum, and checked against its parent and a published reference. Order: the source-free axisymmetric torus (L2, gate RM-B11: uniform loop voltage, steady state with flows, the low-Hartmann-number flow of Bates and Montgomery 1998); the loop-voltage ledger of the helical cylinder (L1, RM-B13); the surface-averaged background (L0, RM-B12); the torus with n = 0, 1 (L4). The few-harmonic rung L3 and the mean-field closure LM wait for L4. Resistive spectra are the Jacobian eigenproblem of P1e. | [model ladder](docs/nonlinear-resistive-mhd.md#model-ladder-on-the-resistive-time-scale), [torus case](docs/vrmhd-torus.md#steady-ohmic-state) |
| **P1g** | **Radiation collisions and QED (started 2026-11-15):** relaxation and power balance of a hot magnetized electron-ion plasma with radiation at p-¹¹B and D-³He conditions and at runaway energies. Order: tree-level QED kernels (electron-ion and electron-electron bremsstrahlung, Møller, Compton, pair creation and annihilation) generated by FortSym and checked against closed forms; relativistic bremsstrahlung loss rates against published fits; radiation friction, then the full radiation collision operator with a declared photon distribution and the measured defect of the friction-only reduction; thermal pair balance against Svensson (1982); the p-¹¹B and D-³He power balance with synchrotron wall reflection through a declared reflection coefficient. | [QED and radiation](docs/qed-radiation.md#first-slice-p1g) |

P1g started on 2026-11-15; its gate is in [its first slice](docs/qed-radiation.md#first-slice-p1g).
P1d is active since 2026-11-19. P1a to P1c and P1f are **existing active campaigns**, not a proposal to restart them. Each
campaign closes a specific observable-level verification gate before adding
adjacent scope. Document unresolved discrepancies explicitly.

### P1b procedure

The comparators are the VMEC implementations VMEC++, VMEX, PARVMEC and
VMEC2000, the nested codes GVEC and DESC, the island codes SPEC and SPECTRE
(multi-region relaxed MHD) and SIESTA, all run on the same input, and the
island and stochastic-region codes HINT and PIES through their published
results. The scoreboard is the
[equilibrium benchmark](https://github.com/itpplasma/benchmark_vmec), which
runs every code on the test inputs of all their repositories, plus the island
cases below. Accuracy is the deviation of physical outputs (volume, magnetic
energy, axis, iota profile, B on the boundary, island width and phase) from a
resolution-converged MILCA reference of the same case, with each comparator's
deviation from it recorded; residuals of different codes are not compared.
Time is wall time including setup, on the same node and thread count. A case
MILCA cannot take is classified as unsupported input, conversion failure,
solver failure or timeout. G-gates are those of
[islands](docs/island-equilibrium.md#verification-ladder-before-machines).

E1 to E5 and E11 (nested) and E6 to E10 and E12 (islands) are two chains,
each in order. E4 needs only E2, E11 needs E1, and E12 needs E7 and E11. A failed gate is recorded with its cause, not relaxed, and the
existing convergence gates of the non-nested (NN1, NN2), island and Beltrami
lanes stay.

| Gate | Case | Comparators | Passed when |
|---|---|---|---|
| E1 | every fixed-boundary deck of the benchmark: prescribed iota or current, all VMEC profile types, adiabatic pressure | VMEC++, VMEX, PARVMEC, DESC, GVEC | each deck converges or carries a classified failure; outputs within the reference's estimated error |
| E2 | W7-X, Landreman-Paul QA and QH, a shaped tokamak | VMEC++, VMEX, PARVMEC, GVEC, DESC | G6; time at most the fastest comparator's at the same accuracy (G-speed) |
| E3 | prescribed toroidal current I(s), bootstrap-like on QA, zero on W7-X | VMEC++, DESC | iota and current within the reference's estimated error |
| E4 | nested free boundary: finite pressure from a displaced seed, W7-X standard case with coil currents, then the coil currents reconstructed | VMEC++ with refined mgrid, GVEC with coils; Biot-Savart and virtual casing | Newton converges without stalls; the vacuum-field error converges with the plasma mesh; G8 |
| E5 | every free-boundary deck of the benchmark | VMEC++, VMEX, PARVMEC | as E1 |
| E6 | islands, fixed boundary: classical stellarator with 1/2 and 3/5 islands; the fixed-boundary QA vacuum case of SPECTRE's paper | SIESTA, SPECTRE | G7; O- and X-points within the estimated error |
| E7 | ideal limit and current constraint with islands: the ideal-limit sequence of [islands](docs/island-equilibrium.md#execution-order-and-status-9-november-2026), prescribed I(s) | SIESTA; SPEC with prescribed current (Baillod et al., J. Plasma Phys. 87, 905870403, 2021) | residual floor reported along the sequence; island observables within the estimated error |
| E8 | islands, free boundary: W7-X vacuum, classical stellarator at finite beta | SPECTRE | G8 and G7 |
| E9 | island divertor and stochastic edge: W7-X 5/5, 5/4, 5/6 at finite beta with open field lines and sheath boundary data; LHD n = 1 resonant perturbation; the stochastic edge of the CFQS stellarator | SPECTRE on closed field lines; HINT and PIES published results | G10 for W7-X and LHD; CFQS against published HINT field-line plots; open-field-line cases by internal convergence only |
| E10 | machine resolution: matrix-free actions, mode-block preconditioning, hp refinement at islands | SPECTRE, SIESTA | time and memory at matched island observables at most those of the fastest comparator |
| E11 | profiles as the frozen flux form (0, chi'(s), Phi'(s)) with chi'^2 + Phi'^2 > 0 instead of iota and Phi: axisymmetric reversed-field pinch through its reversal surface Phi' = 0, then a single-helical-axis pinch state | VMEC2000 with LRFP (poloidal-flux label); SPEC for the helical state | as E1; reversal surface and helical axis within the estimated error |
| E12 | n = 0 resonance: tokamak with iota through zero, where every n = 0 displacement is resonant; inertia of the constrained Hessian at the iota = 0 surface under refinement, then the state with n = 0 islands (ideal current sheet and relaxed) | PARVMEC (nested state); SPEC | inertia reported under refinement; island observables as E6 |

### P1e procedure

1. **Clean room.** All P1e work follows "Source provenance (strict)" in
   [AGENTS.md](AGENTS.md): fresh sessions, publications and licensed
   open-source code only, a declaration in every pull request.
2. **Specification first.** Extend [docs/stability.md](docs/stability.md)
   with a capability matrix of the four reference codes taken from their
   publications, the model hierarchy to be covered, the target observables
   and their acceptance tolerances. Freeze the tolerances before measuring an
   implementation.
3. **Shared base.** Pointwise mathematics is generated by FortSym; spaces and
   assembly use FortFEM, factorization and inertia FortSparse, generic
   numerics FortNum and derivatives FortAD. Generic gaps are repaired in the
   owning library. The MILCA-owned ideal-MHD energy derivation that replaces
   the pinned external generator (see AGENTS.md) comes first.
4. **Verification ladder.** Analytic and asymptotic limits first, then
   published tokamak and stellarator results. Open-source TERPSICHORE and
   DCON/GPEC serve as running comparators; closed codes enter through
   published results only. Each rung is a frozen gate, and discrepancies are
   recorded rather than fitted.
5. **Performance.** Measured on the clusters at accepted accuracy, with the
   CPU path as reference and CUDA optional; generic speedups go to the owning
   library.
6. **Sources.** Each capability lists the publications and open-source
   repositories (commit and licence) it was built and validated from.

## Next work (start small, after a P0 gate is stable)

1. **P1d — First independent Monte Carlo/continuum check:**
   [#366](https://github.com/itpplasma/milca/issues/366). Use a common
   Lorentz or Ornstein–Uhlenbeck collision model, an analytic limit, a MILCA
   continuum implementation and an *independent* GORILLA_APPLETS stochastic
   implementation. Converge mesh, timestep and sample count; compare the same
   weak observable with deterministic and statistical errors kept separate.
   Do not assume test-particle trajectories implement a full linearized
   multispecies Landau operator.
2. **P2 — Evidence-to-fix automation:** make a failed analytic, cross-method or
   cross-code case produce a minimal reproducer, discrepancy classification,
   evidence record, targeted test and optionally a proposed upstream PR.
   Reuse the [benchmark atlas](docs/benchmarks.md), fo and Fortyard's
   [testing-at-scale design](https://github.com/lazy-fortran/fortyard/blob/main/doc/testing-at-scale.md);
   avoid a new universal verification framework.
3. **P2 — Expand comparisons only where physics needs them:** NTV regime limits
   [#232](https://github.com/itpplasma/milca/issues/232), established transport
   codes [#227](https://github.com/itpplasma/milca/issues/227)–[#235](https://github.com/itpplasma/milca/issues/235),
   collision hierarchy [#188](https://github.com/itpplasma/milca/issues/188),
   and native-field orbit coupling [#331](https://github.com/itpplasma/milca/issues/331).
   Promote each independently validated slice on scientific demand.
4. **After P1g, in this order:** nuclear reaction rates N3-N4 for D-T, D-³He
   and p-¹¹B ([nuclear reactions](docs/nuclear-reactions.md)); electron-positron
   plasmas, with positronium formation and the positron species in the kinetic
   models ([pair plasmas](docs/pair-plasmas.md)); the real-time
   QED reference in Fock space (Q0, [QED](docs/qed-radiation.md)); the
   self-gravitating many-body slice G2-G3 on the shared fast-multipole code,
   then the dark-matter hypotheses (cold, warm, self-interacting, fuzzy,
   QUMOND), post-Newtonian dynamics, cosmology and general relativity
   ([gravity](docs/gravity.md)). MILCA is the differentiable forward model of
   a later, separate inverse code for Bayesian comparison of dark-matter and
   gravity hypotheses with observations.

## Verification and integration rules

These habits make a landing informative; none of them is a gate. MILCA lands
fast and finds the remaining errors when it is applied in depth.

1. **Different evidence, not duplicate implementations:** challenge a result by
   analytic/asymptotic limits, conservation/adjoint/reciprocity identities,
   independent symbolic derivation, discretization, Monte Carlo versus
   continuum, established codes and ultimately observations as applicable.
   Shared equations, generated kernels, inputs and historical derivations are
   *correlated evidence*. Record exactly which physical models intersect.
2. **A test must be able to fail:** use failing-before/passing-after
   reproducers, convergence and resolution sweeps, mutation/counterexample
   searches, boundary and sign checks. A golden record or two solvers using
   the same symbolic expression is not independent truth.
3. **Benchmark observable-first:** pin model, conventions, solver revision,
   input, geometry, tolerance, random seed, hardware and actual result status.
   Report continuous error bounds, empirical discretization estimates,
   statistical uncertainty and model mismatch separately. A disagreement is
   allowed to remain *unresolved*.
4. **AI investigates continuously:** generate tests, minimize failures,
   check paper formulas and errata, compare independent CAS systems (including
   FortSym/SymPy/Mathics), diagnose regressions and propose upstream repairs.
   No agent may alter production code *and* silently redefine its authoritative
   benchmark or mathematical acceptance condition.
5. **Review scientific risk, not lines or AI authorship:** automatically
   integrate low-risk, well-tested mechanical changes; audit samples. Demand
   independently checkable evidence for numerical changes. Escalate new physics,
   unresolved high-impact result changes, reference/oracle changes and
   publication claims. External maintainers always choose their own process.
6. **Continuous verification, bounded cost:** affected tests immediately,
   cross-method/downstream tests after integration, expensive reference and
   device tests in a resource-limited background queue, and complete
   evidence for published claims. A passing cached result is valid only for its
   correctly identified full input/toolchain closure. Track stale and failed
   tests; repair the test harness rather than silently relaxing the gate.

## Admission and delivery

Implement the **smallest complete scientific slice** with a named question and
an explicit mathematical model, and land it. Independent oracles, convergence
and error status, certificates and time/memory at matched accuracy are added as
the slice is applied; record what is missing rather than waiting for it. Keep
derivations, reference inputs, reproducible commands and provenance.

Generic numerics belong in FortNum/FortSparse/FortFEM/FortSym/FortIO and
interchange in libneo as defined in [AGENTS.md](AGENTS.md). Maintain independent
specialized reference implementations; do not require established codes to
depend on MILCA.

## Future, not active

All other roadmap families (comprehensive 3D kinetic
response, gyrokinetics, edge/neutrals, reactor integration, generalized orbit
framework, turbulence beyond the drift-turbulence branch above, radiation transport beyond P1g, microscopic/quantum extensions and broad
multifidelity automation) remain **future capability choices**, even where
pieces already run. Their requirements and milestone detail remain in the
[previous staged plan](https://github.com/itpplasma/milca/blob/94caec95e474bd23640278c5e03b50204943f369/ROADMAP.md) and domain
specifications. Start one only when a live physics question and available
verification evidence justify moving it into the current table.

Chris&AI
