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
import netCDF4
import numpy as np
import vmecpp

NPHI, NTHETA = 64, 128
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

with netCDF4.Dataset(a.out / "wout.nc") as f:
    g = lambda k: np.array(f.variables[k][:])
    nsf, nfp = int(g("ns")), int(g("nfp"))
    xm, xn, xmn, xnn = g("xm"), g("xn"), g("xm_nyq"), g("xn_nyq")
    rmnc, zmns, bmnc, gmnc = g("rmnc"), g("zmns"), g("bmnc"), g("gmnc")
    presf = g("presf")
    fsq = [float(g(k)) for k in ("fsqr", "fsqz", "fsql")]
    itfsq = int(g("itfsq")) if "itfsq" in f.variables else -1
    signgs = int(g("signgs"))
sf = np.linspace(0, 1, nsf)
sh = (np.arange(1, nsf) - 0.5) / (nsf - 1)


def full(arr, s):
    return np.array([np.interp(s, sf, arr[:, i]) for i in range(arr.shape[1])])


def half(arr, s):
    if s < sh[0]:
        lo, hi = 0, 1
    elif s > sh[-1]:
        lo, hi = -2, -1
    else:
        return np.array([np.interp(s, sh, arr[1:, i]) for i in range(arr.shape[1])])
    y0, y1 = arr[1:][lo], arr[1:][hi]
    return y0 + (y1 - y0) * (s - sh[lo]) / (sh[hi] - sh[lo])


def surface(s, theta, phi):
    ang = theta[:, None] * xm[None, :] - phi * xn[None, :]
    angb = theta[:, None] * xmn[None, :] - phi * xnn[None, :]
    return np.cos(ang) @ full(rmnc, s), np.sin(ang) @ full(zmns, s), np.cos(angb) @ half(bmnc, s)


th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
rows = []
for k in range(NPHI):
    phi = 2 * np.pi * k / (NPHI * nfp)
    r, z, _ = surface(0.0, np.zeros(1), phi)
    rows.append((phi, r[0], z[0]))
np.savetxt(a.out / "axis.txt", rows, header="phi R Z")
for s in (0.25, 0.5):
    rows = []
    for k in range(NPHI):
        phi = 2 * np.pi * k / (NPHI * nfp)
        r, z, b = surface(s, th, phi)
        rows += [(phi, th[i], r[i], z[i], b[i]) for i in range(NTHETA)]
    np.savetxt(a.out / f"surface_s{s:.2f}.txt", rows, header="phi angle R Z |B|")

# energy int (|B|^2/2 - mu0 p) dV, midpoint rule in s on the half grid
mu0 = 4e-7 * np.pi
tq = np.linspace(0, 2 * np.pi, 65)[:-1]
zq = np.linspace(0, 2 * np.pi / nfp, 65)[:-1]
T, Z = np.meshgrid(tq, zq, indexing="ij")
c = np.cos(T[..., None] * xmn - Z[..., None] * xnn).reshape(-1, len(xmn))
ph = np.interp(sh, sf, presf)
w, jmin, jmax = 0.0, np.inf, -np.inf
for j in range(1, nsf):
    b, gj = c @ bmnc[j], c @ gmnc[j]
    w += np.mean((0.5 * b**2 - mu0 * ph[j - 1]) * np.abs(gj))
    jmin, jmax = min(jmin, gj.min()), max(jmax, gj.max())
summary.update(energy_T2m3=w * (2 * np.pi) ** 2 / (nsf - 1), fsq=fsq, iterations=itfsq,
               unknowns=3 * nsf * len(xm), mnmax=len(xm), signgs=signgs,
               jacobian_range=[float(jmin), float(jmax)], time_s=dt)
(a.out / "summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary))
