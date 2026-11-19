# Model ladder plan

Status: plan, 2026-11-19. Planned items are not support claims.

The [graph of plasma models](model-graph.md) runs from speculative unified
theories to trivial models. This plan states how every rung is implemented and
checked. The middle of the ladder, from many-particle mechanics to MHD and
neoclassical transport, has its own pages. It appears here only where the new
rungs attach to it.

## Contract of a rung

Every rung, from string theory toy models to the ideal gas law, meets the same
contract, derived from the [design principles](../DESIGN.md):

1. **Model page.** The equations with every symbol defined, the state, the
   conserved quantities and the range of validity, on one documentation page.
2. **Reduction arrow.** For each parent, the ordering parameter epsilon, the
   derivation as numbered steps with every dropped term, and its entry in the
   [theory map](theory-map.md) with a checked status.
3. **One residual.** The model is one state and one residual (or one action,
   energy or Hamiltonian) in the canonical equation interface, so that primal,
   JVP, VJP and error bounds come from the same definition.
4. **Generated kernels.** Pointwise mathematics is derived and generated with
   FortSym, with the identities of the model (conservation, symmetry, limits)
   proved symbolically in the generator.
5. **Structure-preserving discretisation.** The best available method for the
   model: symplectic or variational integrators for Hamiltonian systems,
   compatible finite elements for fields, conservative and positivity-preserving
   schemes for kinetics, exact gauge invariance on lattices.
6. **Error control.** Every result carries its discretisation error estimate,
   and a certificate where the mathematics permits. Model, numerical and
   sampling errors stay separate.
7. **Ladder check.** The rung is compared with its parent in their common range
   of validity, and the model defect is measured as a function of epsilon with
   the order the derivation predicts. This is the defining test of the ladder.
8. **Benchmark.** Where other codes or tabulated results exist, accuracy and
   time at equal accuracy are recorded.
9. **Teaching artefact.** One figure or interactive example per rung that shows
   the model, its limit and the comparison with its parent.

A rung is "running" when 1 to 7 hold for at least one case, "designed" when its
page and arrow exist, "speculative" when the theory itself is unconfirmed.

## Order of work

The ladder is built from both ends towards the running middle, cheapest and
most useful rungs first:

| Stage | Rungs | Reason |
|---|---|---|
| L1 | triviality, dimensional analysis, formulary, thermodynamic equilibrium | oracles for every other rung, hours each |
| L2 | analytic equilibria and limits, 0-D balance, circuit models | checks of the running equilibrium and stability rungs, teaching |
| L3 | reduced profile transport, surrogates | consumers of neoclassical and later gyrokinetic results |
| L4 | thermal field theory, QED kinetics, electroweak rates | complete the running QED, collision and nuclear rungs |
| L5 | QCD toy models, quark-gluon plasma kinetics and hydrodynamics | a second plasma with the same ladder structure |
| L6 | semiclassical gravity, general relativity rungs | after GRMHD exists in the middle |
| L7 | speculative toy models | teaching and structural checks only |

<a id="trivial-bottom"></a>
## Trivial bottom

### Triviality

States that every model must reproduce exactly, used as oracles everywhere:

- vacuum: no particles, vacuum fields satisfy the field equations alone,
- a uniform plasma at rest: zero current, uniform density, no force,
- quasineutrality, sum_s q_s n_s = 0, as the limit lambda_D / L -> 0,
- the ideal gas p = n T for every species.

Implementation: a test module that runs every registered model on its trivial
state and requires the residual to vanish to roundoff and the state to stay
unchanged under time stepping. Ladder check: none (no parent below). Every
other rung checks against this one.

### Dimensional analysis

The dimensionless numbers of a case (beta, collisionality nu*, rho* = rho/L,
Lundquist number S, Mach number, lambda_D/L, Coulomb logarithm, plasma
parameter N_D) and the similarity groups of a family of cases (Buckingham pi).

Implementation: the unit-tagged parameter record of every input
([interface](interface-reproducibility.md)) evaluates these numbers, writes them
to every result and states which rungs are valid for them (for example
gyrokinetics for rho* << 1, MHD for nu* >> 1 along the field line). Ladder
check: two cases with equal dimensionless numbers and different units give equal
dimensionless results to roundoff in every model that claims the similarity.

### Formulary

Debye length, plasma and cyclotron frequencies, gyroradii, thermal speeds,
collision times, Coulomb logarithm, mean free path, resistivity (Spitzer), with
their quantum and relativistic corrections as separate terms.

Implementation: one module of FortSym-generated pure functions with intervals,
so that every value carries its bound. Ladder check: each formula against its
parent rung, for example the collision time against the relaxation rate of the
running Landau operator, the Debye length against the screened potential of the
Debye-sphere many-body run, the Spitzer resistivity against the drift-kinetic
conductivity at zero trapping.

### Thermodynamic equilibrium

Maxwell-Boltzmann, Fermi-Dirac and Maxwell-Juttner distributions, Saha and
Boltzmann populations, Planck radiation, the equation of state of an ideal and a
weakly coupled (Debye-Huckel) plasma.

Implementation: exact moments and entropies by FortSym, the Saha solve as a
convex minimisation of the free energy (unique, so Newton with a certificate).
Ladder check: the stationary state of the 1X1V Vlasov-Fokker-Planck laboratory,
the jellium and Fermi-gas runs at finite temperature, and the Debye-Huckel
correction against the one-component-plasma many-body energy at small coupling.

<a id="reduced-bottom"></a>
## Reduced bottom

### Analytic equilibria and limits

Solov'ev and Cerfon-Freidberg equilibria, the large-aspect-ratio circular
tokamak with its Shafranov shift, the cylindrical screw pinch, the
reversed-field-pinch Taylor state. Stability limits as algebraic criteria
(Kruskal-Shafranov, Suydam, Mercier, Troyon beta limit, Greenwald density
limit, Hugill diagram).

Implementation: closed forms generated by FortSym with their expansion
parameter (inverse aspect ratio, beta) and the first neglected term as an error
estimate. Ladder check: each closed form against the running Grad-Shafranov,
3D-equilibrium and ideal-stability rungs. The measured difference must scale
with the predicted power of epsilon (for example the Shafranov shift with
epsilon squared beta_p). The Troyon and Greenwald limits are checked as
empirical fits, so the check reports where the fit and the ideal limit part.

### Circuit and filament models

The plasma as an inductive circuit (current, internal inductance, loop
voltage, coupled to coils and passive structures) and as rigid current
filaments for the vertical instability and its growth rate.

Implementation: a small ODE system with exact inductances from the running
free-boundary Green's functions, integrated by an implicit symplectic or
energy-conserving scheme. Ladder check: growth rate against the n = 0
eigenvalue of the free-boundary Hessian, flux balance against the time-dependent
free-boundary equilibrium.

### 0-D balance

Global power and particle balance of a plasma (heating, fusion power,
bremsstrahlung and line radiation, transport losses with a confinement time),
confinement-time scaling laws as fits, the Lawson triple product, POPCON
operating maps and a systems model of a power plant (net electric power from
plasma, magnets, blanket and recirculating power).

Implementation: algebraic and ODE models with every closure named (scaling
law, radiation model, profile peaking factors) and differentiable end to end, so
that sensitivities and operating-point optimisation use the same residual.
Ladder check: the 0-D energy against the volume integral of the 1.5-D
transport solution with the same sources, the radiation terms against the
radiation rung, the fusion power against the reactivity rung. Scaling laws are
reported as fits with their database range, never extrapolated silently.

### Reduced profile transport

1.5-D flux-surface-averaged transport of density, temperatures, rotation and
current on running equilibria, with neoclassical fluxes from the drift-kinetic
rung and anomalous fluxes from critical-gradient (stiff) models or surrogates.

Implementation: implicit Newton on the coupled profile equations with the
equilibrium in the same residual (no operator splitting), conservative finite
volumes or finite elements in the flux label, adaptive time steps under error
control. Ladder check: fluxes against the drift-kinetic and later gyrokinetic
rungs on the same profiles. The 0-D balance as its volume average.

### Surrogates

Neural networks, Gaussian processes or tables trained on the rungs above
(transport fluxes from gyrokinetics, neoclassical coefficients from drift
kinetics, equilibrium responses), used where their speed is needed.

Implementation: training data generated by MILCA with error bars, the training
range stored with the surrogate, inference in Fortran (dense layers as
generated kernels), and an out-of-range flag that falls back to the parent
rung. Ladder check: error against held-out parent runs over the whole training
range and at its edge. A surrogate is never reported outside its range.

<a id="fundamental-top"></a>
## Fundamental top

### Thermal field theory

Finite-temperature QED: the photon self-energy (Debye screening and plasma
frequency from first principles), damping rates, the equation of state of a
relativistic electron-positron-photon plasma.

Implementation: imaginary-time (Matsubara) and real-time (Keldysh) one-loop
expressions generated by FortSym, momentum integrals with certified quadrature,
hard-thermal-loop resummation where soft scales enter. Ladder check: the
nonrelativistic, weak-coupling limit reproduces the formulary (Debye length,
plasma frequency) and the classical dielectric function of the kinetic rung. The equation of state matches the thermodynamic-equilibrium rung.

### QED kinetics completion

The designed QED rungs ([QED and radiation](qed-radiation.md)) with Boltzmann
collision integrals from tree-level cross sections (Compton, Moller, Bhabha,
bremsstrahlung, pair creation and annihilation), conservative discretisations
of the collision integrals, and the radiation collision operator.

Ladder check: the Landau operator as the small-angle limit of the Moller and
Coulomb integrals, Kompaneets as the small-energy-transfer limit of Compton,
and detailed balance at the Planck and Maxwell-Juttner equilibria to roundoff.

### Electroweak theory

The weak interaction at low energy (Fermi theory and its electroweak parent):
beta decay and electron capture rates, neutrino-electron and neutrino-nucleon
scattering, and neutrino transport by the same kinetic machinery as photons.

Implementation: tree-level amplitudes by FortSym, phase-space integrals with
certified quadrature, Fermi functions with Coulomb corrections. Ladder check:
the pp-chain weak rate of the running nuclear rung, neutron lifetime and
measured decay rates. Neutrino transport against the diffusion limit.

### QCD and chiral effective field theory

Lattice gauge theory as the computable definition of QCD: pure gauge SU(2) and
SU(3) with the Wilson action and hybrid Monte Carlo, then light quarks. Chiral
effective field theory for the nucleon-nucleon force as the bridge to the
nuclear rung.

Implementation: gauge links on a periodic lattice with exact gauge invariance,
symplectic molecular dynamics inside hybrid Monte Carlo, autocorrelation-aware
error bars, continuum extrapolation as error control. The production scale of
lattice QCD is out of reach. MILCA implements the small lattices that show the
method and its limits (confinement, string tension, deconfinement transition).
Ladder check: the pure-gauge deconfinement temperature against published values
at the same lattice spacings. The chiral nucleon-nucleon potential against the
phase shifts that the nuclear rung uses.

### Quark-gluon plasma

The deconfined phase as a plasma: hard-thermal-loop screening and damping,
QCD kinetic theory with elastic and radiative processes, viscous relativistic
hydrodynamics (second-order, causal) for its expansion.

Implementation: the kinetic and hydrodynamic machinery of the electromagnetic
rungs with colour factors, so that both plasmas share code and tests. Ladder
check: the kinetic shear viscosity against the hydrodynamic transport
coefficient. The abelian limit (colour factors to their QED values) reproduces
the electromagnetic rungs exactly.

### Semiclassical gravity

Quantum fields on curved backgrounds: Schwinger pair creation in strong
electric fields, Unruh and Hawking radiation, particle creation in an expanding
universe.

Implementation: mode functions of a scalar and a Dirac field on the given
background by certified ODE integration, Bogoliubov coefficients as the
observable. Ladder check: the flat-space Schwinger rate of the QED rung. The
Hawking spectrum against the Planck rung at the Hawking temperature. The
classical limit of the backreaction against general relativity.

<a id="speculative-top"></a>
## Speculative top

These theories are unconfirmed or incomplete. MILCA implements only tractable
toy limits that reduce to the confirmed rungs, as structural checks and for
teaching. Each page states clearly that the theory is speculative.

| Theory | Tractable limit in MILCA | Check against the confirmed rung |
|---|---|---|
| grand unification | running gauge couplings at one and two loops in the Standard Model and SU(5) | the measured couplings at the Z mass |
| supersymmetry | supersymmetric quantum mechanics (partner potentials, Witten index) | exact spectra of the many-electron rung's one-body problems |
| axions and dark sectors | axion electrodynamics and axion-photon conversion in a plasma | Maxwell's equations at zero axion coupling, plasma dispersion of the kinetic rung |
| string theory | classical relativistic strings (Nambu-Goto), the Veneziano amplitude, Kaluza-Klein reduction of five-dimensional gravity | Maxwell plus general relativity as the zero-mode sector |
| loop quantum gravity | the area and volume spectra of spin networks, spin-foam vertex amplitudes for small graphs | the semiclassical (large-spin) limit against classical geometry |
| asymptotic safety | the Einstein-Hilbert truncation of the functional renormalisation group flow | Newton's constant and the cosmological constant at low scale |
| causal sets | sprinkling of Minkowski and de Sitter spacetime, the d'Alembertian on causal sets | the continuum wave operator as the density grows |

The ladder check of a speculative rung tests only the consistency of the
implementation with its confirmed limits. It is never presented as evidence for
the theory.

## Middle rungs still designed

The designed rungs between the many-particle parents and MHD have their own
pages and stay on those plans: gyrokinetics and gyrofluids
([model family](model-family.md), [moment models](moment-models.md)), extended
MHD, edge physics ([edge](edge-physics.md)), reacting and radiating plasmas
([reacting](reacting-partially-ionized.md), [radiation](radiation.md)), macro
transport ([transport](transport-multiscale.md)) and gravity
([gravity](gravity.md)). Each gets the same contract when it starts.

## Ladder chains

The ladder becomes evidence when a chain of rungs is run end to end and every
link is measured. Planned chains:

1. formulary collision time <- Landau operator <- Lenard-Balescu <- BBGKY <-
   classical many-body, with the Coulomb logarithm from the quantum rung,
2. Spitzer resistivity <- drift kinetics <- Landau operator, then the loop
   voltage of the circuit model on a running equilibrium,
3. Lawson triple product <- 0-D balance <- 1.5-D transport <- neoclassical
   fluxes <- drift kinetics, with the fusion power from the reactivity rung,
4. Debye length <- thermodynamic equilibrium <- thermal QED,
5. Shafranov shift <- Grad-Shafranov <- 3D frozen-in equilibrium in the
   axisymmetric limit,
6. Maxwell plus general relativity <- Kaluza-Klein zero modes <- five-dimensional
   gravity (speculative top, structural only).

Each chain is one figure in the teaching material: the observable at every rung
against epsilon, with the measured model defect and its predicted order.
