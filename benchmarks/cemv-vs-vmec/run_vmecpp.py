"""One VMEC++ run of the head-to-head (docs/equilibrium-roadmap.md, P1).

Runs VMEC++ on a prescribed-iota deck at one resolution (multigrid in ns up
to --ns, MPOL, NTOR), in this process, and writes to OUT:

- wout.nc (for MILCA's force-balance check of the state, test_flux_map_qa
  with MILCA_FM3_VMEC and the stages argument "vmec");
- axis.txt, surface_s0.25.txt, surface_s0.50.txt, summary.json in the layout
  benchmarks/label-qa/compare.py reads.

summary.json holds the wall time of deck reading plus the run, the peak
resident memory of the process (python and numpy included), the iteration
count, the unknowns (the Fourier coefficients of R, Z and lambda on all
surfaces) and the energy int (|B|^2/2 - mu0 p) dV.

--axis-offset D starts from the axis guess of the deck moved by D metres in
R (convergence basin).

Usage: python run_vmecpp.py DECK OUT --ns 101 --mpol 16 --ntor 12 [--threads 8]
       (python with vmecpp and netCDF4, e.g. ~/code/bench-venv)
"""
import argparse
import json
import pathlib
import resource
import sys
import time

t_start = time.perf_counter()
import numpy as np
import vmecpp

from export_wout import export_wout

LADDER = (16, 51, 101, 151, 201)

ap = argparse.ArgumentParser()
ap.add_argument("deck", type=pathlib.Path)
ap.add_argument("out", type=pathlib.Path)
ap.add_argument("--ns", type=int, required=True)
ap.add_argument("--mpol", type=int, required=True)
ap.add_argument("--ntor", type=int, required=True)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--ftol", type=float, default=1e-14)
ap.add_argument("--niter", type=int, default=30000)
ap.add_argument("--axis-offset", type=float, default=0.0)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

t0 = time.perf_counter()
vin = vmecpp.VmecInput.from_file(a.deck)
# the boundary of the deck, truncated or padded to MPOL, NTOR
vin.rbc = vmecpp.VmecInput.resize_2d_coeff(vin.rbc, a.mpol, a.ntor)
vin.zbs = vmecpp.VmecInput.resize_2d_coeff(vin.zbs, a.mpol, a.ntor)
vin.mpol, vin.ntor = a.mpol, a.ntor
for name in ("raxis_c", "zaxis_s"):
    v = np.zeros(a.ntor + 1)
    old = np.array(getattr(vin, name))
    v[: min(len(old), a.ntor + 1)] = old[: a.ntor + 1]
    setattr(vin, name, v)
if a.axis_offset:
    v = np.array(vin.raxis_c)
    v[0] += a.axis_offset
    vin.raxis_c = v
ns = [n for n in LADDER if n < a.ns] + [a.ns]
vin.ns_array = np.array(ns)
vin.ftol_array = np.full(len(ns), a.ftol)
vin.niter_array = np.full(len(ns), a.niter)
vin.return_outputs_even_if_not_converged = False
status, reason = "converged", ""
try:
    res = vmecpp.run(vin, max_threads=a.threads, verbose=False)
except Exception as e:  # the failure mode is the result
    status, reason = "failed", str(e).strip().splitlines()[-1][:300]
dt = time.perf_counter() - t0
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
peak_mb = peak / (1024**2 if sys.platform == "darwin" else 1024)
summary = {"code": "VMEC++", "deck": a.deck.name, "ns": a.ns, "mpol": a.mpol, "ntor": a.ntor,
           "threads": a.threads, "status": status, "reason": reason, "time_s": dt,
           "time_with_imports_s": time.perf_counter() - t_start, "peak_memory_MB": peak_mb,
           "axis_offset_m": a.axis_offset}
if status != "converged":
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    sys.exit(0)
res.wout.save(str(a.out / "wout.nc"))

summary.update(export_wout(a.out / "wout.nc", a.out))
(a.out / "summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary))
