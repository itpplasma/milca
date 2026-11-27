"""VMEC++ reference for the flux map on a standard fixed-boundary deck.

The flux map prescribes the rotational transform and the pressure as
functions of s (normalized toroidal flux). For a deck with prescribed current
(NCURR = 1) this script

1. runs VMEC++ on the deck as given;
2. fits its iota(s) by a power series of degree --degree and writes the deck
   with that iota prescribed (NCURR = 0), keeping boundary and pressure;
3. runs VMEC++ on the prescribed-iota deck (the problem the flux map solves);
4. writes, for that state, plain-text references in the layout of
   benchmarks/label-qa (axis, surfaces s = 0.25, 0.5, 1 with |B|, iota,
   summary.json) in OUT/reference and keeps both wout files in OUT/work for
   MILCA's residual check of the VMEC state. OUT defaults to
   $MILCA_DATA_DIR/benchmarks/flux-map-standard/NAME.

Usage: python make_reference.py DECK NAME [--out DIR] [--ns ...] [--degree 6]
       (python with vmecpp and netCDF4, e.g. ~/code/bench-venv)
"""
import argparse
import importlib.metadata as md
import json
import os
import pathlib
import re
import sys
import time

import netCDF4
import numpy as np
import vmecpp

NPHI, NTHETA, SURFACES = 64, 128, (0.25, 0.5, 1.0)

ap = argparse.ArgumentParser()
ap.add_argument("deck", type=pathlib.Path)
ap.add_argument("name")
ap.add_argument("--out", type=pathlib.Path, default=None)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--ns", type=int, nargs="+", default=[16, 51, 101])
ap.add_argument("--ftol", type=float, default=1e-14)
ap.add_argument("--niter", type=int, default=30000)
ap.add_argument("--degree", type=int, default=6)
ap.add_argument("--reuse", action="store_true", help="reuse wout files present in OUT/work")
ap.add_argument("--mpol", type=int, default=None, help="override MPOL of the deck (reference resolution)")
ap.add_argument("--ntor", type=int, default=None, help="override NTOR of the deck")
a = ap.parse_args()


def milca_data(relative):
    root = os.environ.get("MILCA_DATA_DIR")
    if not root:
        sys.exit(f"set MILCA_DATA_DIR to the external MILCA data directory ({relative} is read or written there)")
    return pathlib.Path(root) / relative


OUT = a.out or milca_data(f"benchmarks/flux-map-standard/{a.name}")
WORK = OUT / "work"
REF = OUT / "reference"
for d in (OUT, WORK, REF):
    d.mkdir(parents=True, exist_ok=True)


def run_vmecpp(vin, tag):
    path = WORK / f"wout_{a.name}_{tag}.nc"
    if a.reuse and path.exists():
        return path, None
    vin.ns_array = np.array(a.ns)
    vin.ftol_array = np.full(len(a.ns), a.ftol)
    vin.niter_array = np.full(len(a.ns), a.niter)
    vin.return_outputs_even_if_not_converged = False
    t0 = time.perf_counter()
    res = vmecpp.run(vin, max_threads=a.threads, verbose=False)
    dt = time.perf_counter() - t0
    path = WORK / f"wout_{a.name}_{tag}.nc"
    res.wout.save(str(path))
    return path, dt


class Wout:
    """The wout fields used here (VMEC++ netCDF)."""

    def __init__(self, path):
        with netCDF4.Dataset(path) as f:
            g = lambda k: np.array(f.variables[k][:])
            self.ns, self.nfp = int(g("ns")), int(g("nfp"))
            self.xm, self.xn = g("xm"), g("xn")
            self.xmn, self.xnn = g("xm_nyq"), g("xn_nyq")
            self.rmnc, self.zmns = g("rmnc"), g("zmns")
            self.bmnc, self.gmnc = g("bmnc"), g("gmnc")
            self.iotaf, self.presf = g("iotaf"), g("presf")
            self.phiedge = float(g("phi")[-1])
            self.wb, self.wp = float(g("wb")), float(g("wp"))
            self.volume = float(g("volume_p"))
            self.mpol, self.ntor = int(g("mpol")), int(g("ntor"))
            self.fsq = [float(g(k)) for k in ("fsqr", "fsqz", "fsql")]
            self.lasym = bool(int(g("lasym__logical__"))) if "lasym__logical__" in f.variables else False
        self.s = np.linspace(0, 1, self.ns)

    def coef(self, arr, s, half=False):
        if half:
            sh = (np.arange(1, self.ns) - 0.5) / (self.ns - 1)
            lo, hi = (0, 1) if s < sh[0] else ((-2, -1) if s > sh[-1] else (None, None))
            if lo is None:
                return np.array([np.interp(s, sh, arr[1:, i]) for i in range(arr.shape[1])])
            x0, x1 = sh[lo], sh[hi]
            y0, y1 = arr[1:][lo], arr[1:][hi]
            return y0 + (y1 - y0) * (s - x0) / (x1 - x0)
        return np.array([np.interp(s, self.s, arr[:, i]) for i in range(arr.shape[1])])

    def surface(self, s, theta, phi):
        rc, zs = self.coef(self.rmnc, s), self.coef(self.zmns, s)
        bc = self.coef(self.bmnc, s, half=True)
        ang = theta[:, None] * self.xm[None, :] - phi * self.xn[None, :]
        angb = theta[:, None] * self.xmn[None, :] - phi * self.xnn[None, :]
        return (np.cos(ang) @ rc, np.sin(ang) @ zs, np.cos(angb) @ bc)

    def energy(self):
        """int (|B|^2/2 - mu0 p) dV over the torus, B in T, p in Pa: the
        gamma = 0 energy in the flux map's units (T^2 m^3), midpoint in s."""
        mu0 = 4e-7 * np.pi
        th = np.linspace(0, 2 * np.pi, 65)[:-1]
        ze = np.linspace(0, 2 * np.pi / self.nfp, 65)[:-1]
        T, Z = np.meshgrid(th, ze, indexing="ij")
        c = np.cos(T[..., None] * self.xmn - Z[..., None] * self.xnn).reshape(-1, len(self.xmn))
        sh = (np.arange(1, self.ns) - 0.5) / (self.ns - 1)
        ph = np.interp(sh, self.s, self.presf)
        w = 0.0
        for j in range(1, self.ns):
            b, gj = c @ self.bmnc[j], c @ self.gmnc[j]
            w += np.mean((0.5 * b**2 - mu0 * ph[j - 1]) * np.abs(gj))
        return w * (2 * np.pi) ** 2 / (self.ns - 1)


def compare(w1, w2):
    th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
    out = {}
    for s in (0.0,) + SURFACES:
        d = 0.0
        for k in range(NPHI):
            phi = 2 * np.pi * k / (NPHI * w1.nfp)
            r1, z1, _ = w1.surface(s, th, phi)
            r2, z2, _ = w2.surface(s, th, phi)
            d = max(d, np.max(np.hypot(r1 - r2, z1 - z2)))
        out[f"s={s}"] = d
    out["iota"] = float(np.max(np.abs(np.interp(w2.s, w1.s, w1.iotaf) - w2.iotaf)))
    return out


record = {"deck": str(a.deck), "vmecpp": md.version("vmecpp"), "ns": a.ns, "ftol": a.ftol,
          "threads": a.threads}
text = a.deck.read_text()
if a.mpol:
    text = re.sub(r"(?im)^(\s*)MPOL\s*=[^,\n]*,?", rf"\g<1>MPOL = {a.mpol}", text)
if a.ntor:
    text = re.sub(r"(?im)^(\s*)NTOR\s*=[^,\n]*,?", rf"\g<1>NTOR = {a.ntor}", text)
deck_given = OUT / f"input.{a.name}"
deck_given.write_text(text)
record["mpol_ntor_override"] = [a.mpol, a.ntor]
vin = vmecpp.VmecInput.from_file(deck_given)
if vin.ncurr == 1:
    p_free, t_free = run_vmecpp(vin, "as_given")
    wf = Wout(p_free)
    coef = np.polynomial.polynomial.polyfit(wf.s, wf.iotaf, a.degree)
    fit = float(np.max(np.abs(np.polynomial.polynomial.polyval(wf.s, coef) - wf.iotaf)))
    # the prescribed-iota deck as INDATA text (what MILCA's import_vmec reads)
    t = re.sub(r"(?im)^\s*(AI|PIOTA_TYPE)\s*=.*\n", "", text)
    t = re.sub(r"(?im)^(\s*)NCURR\s*=\s*\d+\s*,?",
               lambda m: f"{m.group(1)}NCURR = 0\n{m.group(1)}PIOTA_TYPE = 'power_series'\n"
                         f"{m.group(1)}AI = " + " ".join(f"{c:.16e}" for c in coef) + "\n", t, count=1)
    deck_iota = OUT / f"input.{a.name}-iota"
    deck_iota.write_text(t)
    vin = vmecpp.VmecInput.from_file(deck_iota)
    record.update(time_as_given_s=t_free, iota_fit_degree=a.degree, iota_fit_max_error=fit,
                  fsq_as_given=wf.fsq)
else:
    coef, wf = None, None
p_iota, t_iota = run_vmecpp(vin, "iota")
wi = Wout(p_iota)
record.update(time_iota_s=t_iota, mpol=wi.mpol, ntor=wi.ntor, ns_final=wi.ns, fsq_iota=wi.fsq,
              iota_axis=float(wi.iotaf[0]), iota_edge=float(wi.iotaf[-1]),
              pressure_axis_Pa=float(wi.presf[0]), lasym=wi.lasym)
if wf is not None:
    record["iota_vs_as_given"] = compare(wi, wf)

# radial resolution of the reference: ns 151 against a final ns below it
if max(a.ns) < 151:
    p151 = WORK / f"wout_{a.name}_iota_ns151.nc"
    if not (a.reuse and p151.exists()):
        vin.ns_array = np.array(a.ns + [151])
        vin.ftol_array = np.full(len(a.ns) + 1, a.ftol)
        vin.niter_array = np.full(len(a.ns) + 1, a.niter)
        res = vmecpp.run(vin, max_threads=a.threads, verbose=False)
        res.wout.save(str(p151))
    record["iota_ns151_vs_final"] = compare(Wout(p151), wi)

W = wi.energy()
record.update(energy_T2m3=W, wb_times_4pi2=wi.wb * (2 * np.pi) ** 2, volume=wi.volume,
              phiedge=wi.phiedge, nfp=wi.nfp, wout_iota=str(p_iota))
if coef is not None:
    np.savetxt(REF / "iota_coefficients.txt", coef,
               header=f"iota(s) = sum_k c_k s^k; s normalized toroidal flux; "
                      f"PHIEDGE {wi.phiedge:.12e} Wb, NFP {wi.nfp}")
th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
rows = []
for k in range(NPHI):
    phi = 2 * np.pi * k / (NPHI * wi.nfp)
    r, z, _ = wi.surface(0.0, np.zeros(1), phi)
    rows.append((phi, r[0], z[0]))
np.savetxt(REF / "axis.txt", rows, header="phi R Z (m), cylindrical phi over one period")
for s in SURFACES:
    rows = []
    for k in range(NPHI):
        phi = 2 * np.pi * k / (NPHI * wi.nfp)
        r, z, b = wi.surface(s, th, phi)
        rows += [(phi, t, ri, zi, bi) for t, ri, zi, bi in zip(th, r, z, b)]
    np.savetxt(REF / f"surface_s{s:.2f}.txt", rows, fmt="%.10e",
               header=f"s = {s}: phi theta_VMEC R Z (m) |B| (T)")
np.savetxt(REF / "iota.txt", np.c_[wi.s, wi.iotaf], header="s iota (VMEC++ full grid)")
np.savetxt(REF / "pressure.txt", np.c_[wi.s, wi.presf], header="s p (Pa, VMEC++ full grid)")
(REF / "summary.json").write_text(json.dumps(record, indent=1, default=float))
print(json.dumps(record, indent=1, default=float))
