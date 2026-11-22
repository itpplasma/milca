# Coils, walls and passive structures

Status: plan of 2026-11-22. One model of the conductors outside the plasma for
every MILCA application that needs them: free-boundary equilibrium
([free-boundary-ale.md](free-boundary-ale.md), [eulerian-label-plan.md](eulerian-label-plan.md)
E6), equilibrium reconstruction (TIAGO), linear stability with walls
([stability.md](stability.md), P1e) and time-dependent free-boundary
evolution. Nothing here is implemented beyond what the gates of those plans
record.

## 1. Principles

- **One model, four uses.** The same conductor description and the same
  operators serve equilibrium, reconstruction, stability and evolution; the
  uses differ only in which quantities are given and which are unknown.
- **Currents as parameters, not as solver state.** A static equilibrium is a
  function u(I) of the external currents I = (coil circuits, passive
  elements). The forward model takes I as input (prescribed, or from a circuit
  evolution); the inverse problem (TIAGO) treats I as unknowns with priors.
  MILCA returns u and the derivatives du/dI = -J^-1 dF/dI on the Newton
  factorisation (extra right-hand sides), and the JVP/VJP of observables, as
  the map lane's circuit interface already does
  ([free-boundary-nested-map-3d.md](free-boundary-nested-map-3d.md)).
- **Symbolic and differentiable.** The pointwise terms (the field and energy
  of a circuit element on a quadrature point, the thin-wall magnetic and
  dissipation densities, the Ohm's-law and circuit couplings) are FortSym
  leaves with their first and second derivatives and proved identities, like
  the label-lane energy. Conductivities, geometry parameters and currents are
  differentiable inputs.
- **Verified against limits.** Every representation has an analytic or
  independent oracle (section 6) before it enters a physics gate.
- **Geometry from CAD offline.** GECKO (itpplasma, MIT, Python, on libneo)
  turns STEP winding packs into coil filaments and STL vessel meshes into
  wall surfaces and R-Z slices. MILCA reads its neutral output files in
  Fortran; no Python enters the executable.

## 2. Representations

| Id | Conductor | Unknowns | Field in the plasma and the shell | Use |
|---|---|---|---|---|
| W0 | coils: 3D polylines (`milca_coils_3d`), axisymmetric rectangles or filaments (`milca_gs_coils`), grouped into circuits by an incidence matrix | circuit currents I_c | Biot-Savart / Green function per unit current, a basis of source fields | all |
| W1 | axisymmetric passive structures (vessel, stabilising plates): toroidal filaments or segments with resistance | element currents I_p, optionally reduced to the eigenmodes of L^-1 R | as W0 | equilibrium (n = 0 eddy currents at a time instant), reconstruction, evolution, vertical stability |
| W2 | thin resistive 3D wall Gamma_w: sheet current K = n x grad_s Psi_w with a stream function Psi_w on the surface, sheet resistivity eta/d, holes and ports through harmonic currents (homology) | Psi_w (and the net currents round the holes) | boundary integral of K (single layer), or the jump [n x B] = K inside the shell | RWM stability, 3D eddy currents, evolution with n != 0 |
| W3 | ideal wall: B.n = 0 frozen at Gamma_w | none | natural condition of the shell potential | the ideal-wall limit of stability, flux-conserving equilibria (fast events, superconducting shells) |
| W4 | thick or volumetric conductors | volume currents | volume integral | only if W2 is shown insufficient (thick blankets, halo currents) |
| L | first wall (limiter surface) | none | none | limiter contact of the plasma boundary |

Structure in the free-boundary shell: Gamma_p (plasma) - shell - Gamma_w
(wall, W2 or W3, an interior surface of the shell) - shell - Gamma_c (fixed
coupling surface) - exterior through the Steklov-Poincare map per toroidal
mode (axisymmetric Gamma_c, FortFEM `tm_bem`) or a 3D BEM for a shaped
Gamma_c. Gamma_c stays a computational surface; the wall is never forced
onto it. The vacuum field is split, B = sum_k I_k B_k + c B_fil + grad Phi:
the W0/W1 basis fields, the circulation c (plasma current) of a filament
inside the plasma, and the single-valued Phi (E6c).

## 3. Uses

### 3.1 Static free-boundary equilibrium

- Default wall model **exterior**: coils and passive currents as given
  parameters, the vacuum continued to infinity. A resistive wall carries no
  current in a steady state; at a time instant of a discharge the passive
  currents (W1, W2 projected) are additional given sources.
- **ideal**: a conducting Gamma_w (W3), for flux-conserving states and the
  ideal-wall limit.
- Limiter contact (L): the plasma boundary touches the first wall; the
  boundary flux is then fixed by the contact point instead of a prescribed
  enclosed flux. A complementarity condition (contact or not) between solves,
  not inside Newton.
- Literature: coil-only free boundary is the stellarator standard (VMEC with
  NESTOR [N], DESC [A], SPEC through virtual casing [N]); axisymmetric codes
  carry passive circuits with prescribed currents (FreeGSNKE [A], FEEQS.M and
  CEDRES++ [A], TokaMaker [A]).

### 3.2 Reconstruction (TIAGO)

- MILCA: the forward map u(theta, I) and its derivatives with respect to the
  profile parameters theta and all external currents I.
- TIAGO: I_c from measured coil currents with their uncertainty; I_p of the
  passive elements as unknowns with a physics prior, the circuit model driven
  by measured loop voltages, I_p ~ N(I_circuit, Sigma_circuit), or the
  eigenmode amplitudes of W1 with a prior on their decay. Flux loops,
  pick-up coils and Rogowskis enter through their geometry as observables
  with exact derivatives.
- Literature: EFIT fits vessel segment currents or takes them from loop
  voltages [N]; LIUQE estimates segment currents from per-segment flux-loop
  voltages and calibrated resistances [A]; CLISTE and IDE fit passive
  conductor currents [N]. Bayesian current tomography (Minerva) models
  plasma and coil currents with Gaussian priors [A]; passive currents as
  inferred variables with circuit priors were not found in the survey and are
  a candidate contribution of the MILCA-TIAGO pair.

### 3.3 Linear stability with walls (P1e)

- No-wall and ideal-wall limits: the vacuum energy of the plasma
  displacement with the exterior map (no wall) and with W3 at Gamma_w.
- Resistive wall modes: the extended eigenproblem of plasma and wall (W2),
  growth rates on the wall time; feedback and error-field coils as W0
  circuits with given gains. Thin wall first; 3D wall structure (ports)
  through the hole currents of W2.
- Literature: DCON with VACUUM gives the ideal-wall and no-wall limits [N];
  MARS-F/K: a thin resistive wall, feedback coils and kinetic damping [N];
  CarMa couples the plasma response on a surface to 3D volumetric
  conductors [N]; VALEN: 3D thin-shell circuits [N]; STARWALL: a triangulated
  thin wall that couples toroidal modes, and CASTOR3D solves the extended
  eigenproblem of plasma inertia and wall for 2D and 3D equilibria [A];
  M3D-C1 meshes wall and vacuum as regions, so wall thickness is
  arbitrary [A].

### 3.4 Time-dependent free boundary

- Quasi-static evolution: at each step the equilibrium with the current
  external currents, and the circuit equations M dI/dt + R I = V - d Psi_plasma/dt
  for coils, W1 and W2, implicit in time (BDF), one Newton(-Krylov) system
  for plasma, shell and circuits; vertical stability control as a W0
  circuit with a controller.
- Later: nonlinear MHD with a wall through the NN lane, where halo currents
  require the wall inside the domain or a dedicated coupling.
- Literature: CREATE-L/NL (vessel eigenmodes, controller design), DINA, TSC,
  CORSICA (vessel filaments, feedback) [N]; FreeGSNKE's evolutive solver
  (circuits coupled through Poynting's theorem, implicit Euler with
  Newton-Krylov) [A]; JOREK-STARWALL evolves wall currents implicitly from
  STARWALL's response matrices [A], JOREK-CARIDDI couples volumetric
  conductors through virtual casing [A]; M3D-C1, NIMROD and JOREK agree on
  wall forces in an NSTX-like benchmark [A].

## 4. Numerics

- Dense blocks (coil and wall mutual inductances, BEM single layers, the
  exterior map) belong to FortFEM; above a few thousand unknowns they are
  compressed (hierarchical matrices / HODLR, as ThinCurr [A] and the
  JOREK-STARWALL work [A]) or summed by FMM (F6 of free-boundary-ale.md).
- The W2 stream function lives in the spline trace space of Gamma_w
  (FortFEM surface splines) with the hole currents as extra unknowns; the
  magnetic energy (1/2) int int K(x).K(y)/(4 pi |x - y|) and the dissipation
  int eta |K|^2/d are the FortSym leaves, with the singular quadrature of the
  BEM layer.
- Eigenmode reduction of W1/W2 (the generalised eigenproblem of L and R) is
  optional and exact up to truncation; the reduced basis is a change of
  variables, so derivatives pass through.

## 5. Order

1. E6c (exterior map, plasma-current circulation) and the solver module with
   a circuit source basis (W0); axisymmetric coil cases against the 2D spline
   free boundary (`gs_spline_newton_free`) on the gs_free fixture.
2. W1 as additional circuits; limiter contact (L); GECKO readers for coil
   filaments and first-wall slices.
3. Derivatives with respect to the circuit currents (JVP/VJP) for TIAGO.
4. W3 ideal wall at an interior Gamma_w (no-wall and ideal-wall limits for
   P1e).
5. W2 thin resistive wall with holes; the RWM extended eigenproblem.
6. Time-dependent circuits with quasi-static equilibrium; vertical
   displacement growth rate.

## 6. Gates

| Gate | Case | Oracle | Passed when |
|---|---|---|---|
| WG0 | coil basis | Biot-Savart of a circular loop and of a straight segment (closed forms), `coils_3d` against `gs_coils` for an axisymmetric loop | fields agree to quadrature accuracy; derivatives with respect to currents exact (linear) |
| WG1 | axisymmetric vessel eigenmodes | L/R decay times of a thin circular torus or cylinder shell (closed form for the cylinder) | decay times at the discretisation order |
| WG2 | thin-wall RWM in a cylinder | analytic growth rate of the ideal-wall/no-wall window for a cylindrical plasma with a thin resistive shell | growth rate at the discretisation order; gamma tau_w finite between the limits |
| WG3 | ideal and no-wall limits | published DCON limits for a reference equilibrium | beta limits within the stated error |
| WG4 | vertical displacement | analytic n = 0 rigid-displacement growth rate with a resistive wall (elongated cylinder model) | growth rate on the wall time at the discretisation order |
| WG5 | 3D wall with holes | a toroidal shell with a port: hole currents and decay times against an independent open-source thin-shell code (ThinCurr, OFT, with commit and licence recorded) | decay times within the discretisation error |
| WG6 | wall forces | the published NSTX-like benchmark (Krebs et al.) | forces within the spread of the published codes |

## 7. Sources

Reading state of this survey: [A] abstract or documentation only, [N] not
checked in this survey (confirm before a claim rests on it). No non-public
code was used.

- Lao et al., Nucl. Fusion 25 (1985) 1611 (EFIT) [N].
- Moret et al., Fusion Eng. Des. 91 (2015) 1 (LIUQE) [A].
- McCarthy, Phys. Plasmas 6 (1999) 3554 (CLISTE) [N]; Fischer et al., Fusion
  Sci. Technol. 69 (2016) 526 (IDE) [N].
- Svensson and Werner, Plasma Phys. Control. Fusion 50 (2008) 085002 [A];
  Kwak et al., Nucl. Fusion 62 (2022) 126069 [A]; von Nessi and Hole, Plasma
  Phys. Control. Fusion 56 (2014) 114011 [A].
- Hirshman, van Rij and Merkel, Comput. Phys. Commun. 43 (1986) 143 (VMEC);
  Merkel, J. Comput. Phys. 66 (1986) 83 (NESTOR) [N].
- Conlin et al., arXiv:2412.05680 (DESC free boundary) [A].
- Pentland et al., Phys. Scr. 100 (2025) 025608, arXiv:2407.12432
  (FreeGSNKE) [A]; Amorisco et al., Phys. Plasmas 31 (2024) 042517 [N].
- Heumann et al., J. Plasma Phys. 81 (2015) 905810301 (FEEQS.M, CEDRES++) [A].
- Hansen et al., arXiv:2311.07719 (TokaMaker) [A]; arXiv:2412.14962
  (ThinCurr) [A].
- Albanese and Villone, Nucl. Fusion 38 (1998) 723 (CREATE-L) [N];
  Khayrutdinov and Lukash, J. Comput. Phys. 109 (1993) 193 (DINA) [N];
  Jardin et al., J. Comput. Phys. 66 (1986) 481 (TSC) [N].
- Chance, Phys. Plasmas 4 (1997) 2161 (VACUUM) [N]; Glasser, Phys. Plasmas 23
  (2016) 072505 (DCON) [N]; Liu et al., Phys. Plasmas 7 (2000) 3681 (MARS-F)
  [N]; Portone et al., Plasma Phys. Control. Fusion 50 (2008) 085004 (CarMa)
  [N]; Bialek et al., Phys. Plasmas 8 (2001) 2170 (VALEN) [N].
- Merkel and Strumberger, arXiv:1508.04911 (STARWALL) [A]; Strumberger and
  Günter, Nucl. Fusion 57 (2017) 016032 (CASTOR3D) [A].
- Ferraro et al., Phys. Plasmas 23 (2016) 056114 (M3D-C1 resistive wall) [A].
- Hölzl et al., J. Phys. Conf. Ser. 401 (2012) 012010 (JOREK-STARWALL) [A];
  arXiv:2501.05956 (full-MHD coupling) [A]; Cipolletta et al.,
  arXiv:2404.16546 (matrix compression) [A]; Isernia et al., EPS 2022
  P1a.107 (JOREK-CARIDDI) [A].
- Krebs et al., arXiv:1908.02387 and arXiv:2011.04523 (wall-force
  benchmark) [A].
- Malhotra et al., Plasma Phys. Control. Fusion 62 (2020) 024004 (BIEST) [N].
