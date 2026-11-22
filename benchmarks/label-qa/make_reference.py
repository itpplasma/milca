"""VMEC reference for the label-lane benchmark on the Landreman-Paul QA deck.

1. VMEC++ on benchmarks/nested/input.landreman-paul-qa as given (fixed
   boundary, p = 0, zero net toroidal current: NCURR = 1, CURTOR = 0).
2. Power-series fit of its iota(s), s the normalized toroidal flux; the deck
   with that iota prescribed (NCURR = 0) is written to
   input.landreman-paul-qa-iota and run again. That problem is the one the
   label lane solves: chi = s, F = -PHIEDGE chi / (2 pi),
   G' = -iota PHIEDGE / (2 pi), p = 0.
3. VMEC2000 (xvmec, optional: --xvmec PATH) on both decks.
4. Plain-text reference of the prescribed-iota VMEC++ state in reference/.

Usage: python make_reference.py [--xvmec PATH] [--threads N] [--ns ...]
"""
import argparse
import importlib.metadata as md
import json
import pathlib
import re
import subprocess
import time

import netCDF4
import numpy as np
import vmecpp

HERE = pathlib.Path(__file__).resolve().parent
DECK = HERE.parent / "nested" / "input.landreman-paul-qa"
OUT = HERE / "reference"
WORK = HERE / "work"
NPHI, NTHETA, SURFACES, DEGREE = 64, 128, (0.25, 0.5, 1.0), 6

ap = argparse.ArgumentParser()
ap.add_argument("--xvmec", default=None)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--ns", type=int, nargs="+", default=[16, 51, 101])
ap.add_argument("--ftol", type=float, default=1e-14)
a = ap.parse_args()
OUT.mkdir(exist_ok=True)
WORK.mkdir(exist_ok=True)


def run_vmecpp(vin, name):
    vin.ns_array = np.array(a.ns)
    vin.ftol_array = np.full(len(a.ns), a.ftol)
    vin.niter_array = np.full(len(a.ns), 30000)
    vin.return_outputs_even_if_not_converged = False
    t0 = time.perf_counter()
    res = vmecpp.run(vin, max_threads=a.threads, verbose=False)
    dt = time.perf_counter() - t0
    path = WORK / f"wout_{name}.nc"
    res.wout.save(str(path))
    return path, dt


def run_xvmec(deck_text, name):
    """VMEC2000 on deck_text with the ladder of --ns and --ftol."""
    d = WORK / f"x_{name}"
    d.mkdir(exist_ok=True)
    k = len(a.ns)
    t = re.sub(r"NS_ARRAY\s*=.*", "NS_ARRAY = " + " ".join(map(str, a.ns)), deck_text)
    t = re.sub(r"NITER_ARRAY\s*=.*", "NITER_ARRAY = " + " ".join(["30000"] * k), t)
    t = re.sub(r"FTOL_ARRAY\s*=.*", "FTOL_ARRAY = " + " ".join([f"{a.ftol:.1e}"] * k), t)
    (d / f"input.{name}").write_text(t)
    t0 = time.perf_counter()
    subprocess.run([a.xvmec, f"input.{name}"], cwd=d, check=True,
                   stdout=open(d / "x.log", "w"), stderr=subprocess.STDOUT)
    return d / f"wout_{name}.nc", time.perf_counter() - t0


class Wout:
    """The few wout fields used here, from netCDF (VMEC++ and VMEC2000)."""

    def __init__(self, path):
        with netCDF4.Dataset(path) as f:
            g = lambda k: np.array(f.variables[k][:])
            self.ns, self.nfp = int(g("ns")), int(g("nfp"))
            self.xm, self.xn = g("xm"), g("xn")
            self.xmn, self.xnn = g("xm_nyq"), g("xn_nyq")
            self.rmnc, self.zmns = g("rmnc"), g("zmns")        # (ns, mnmax), full grid
            self.bmnc, self.gmnc = g("bmnc"), g("gmnc")        # half grid, row 0 unused
            self.iotaf, self.phiedge = g("iotaf"), float(g("phi")[-1])
            self.wb, self.volume = float(g("wb")), float(g("volume_p"))
            self.mpol, self.ntor = int(g("mpol")), int(g("ntor"))
            self.fsq = [float(g(k)) for k in ("fsqr", "fsqz", "fsql")]
        self.s = np.linspace(0, 1, self.ns)

    def coef(self, arr, s, half=False):
        """Coefficients at s, linear in s between grid surfaces."""
        if half:
            sh = (np.arange(1, self.ns) - 0.5) / (self.ns - 1)
            return np.array([np.interp(s, sh, arr[1:, i], left=None, right=None)
                             if sh[0] <= s <= sh[-1] else
                             np.polyval(np.polyfit(sh[-2:] if s > sh[-1] else sh[:2],
                                                   arr[-2:, i] if s > sh[-1] else arr[1:3, i],
                                                   1), s)
                             for i in range(arr.shape[1])])
        return np.array([np.interp(s, self.s, arr[:, i]) for i in range(arr.shape[1])])

    def surface(self, s, theta, phi):
        rc, zs = self.coef(self.rmnc, s), self.coef(self.zmns, s)
        bc = self.coef(self.bmnc, s, half=True)
        ang = theta[:, None] * self.xm[None, :] - phi * self.xn[None, :]
        angb = theta[:, None] * self.xmn[None, :] - phi * self.xnn[None, :]
        return (np.cos(ang) @ rc, np.sin(ang) @ zs, np.cos(angb) @ bc)

    def energy(self):
        """W = int |B|^2 / 2 dV over the torus, B in T (no mu0), m^3 T^2,
        midpoint rule in s on the half grid, 64 x 64 points per period."""
        th = np.linspace(0, 2 * np.pi, 65)[:-1]
        ze = np.linspace(0, 2 * np.pi / self.nfp, 65)[:-1]
        T, Z = np.meshgrid(th, ze, indexing="ij")
        ang = T[..., None] * self.xmn - Z[..., None] * self.xnn
        c = np.cos(ang).reshape(-1, len(self.xmn))
        w = 0.0
        for j in range(1, self.ns):
            b, gj = c @ self.bmnc[j], c @ self.gmnc[j]
            w += np.mean(0.5 * b**2 * np.abs(gj))
        # mean over (theta, zeta) -> times (2 pi)^2, ds = 1/(ns - 1)
        return w * (2 * np.pi) ** 2 / (self.ns - 1)


def compare(w1, w2):
    """Max distance of the axis and of the s = 0.25, 0.5, 1 surfaces (per
    point, same VMEC angles), and max |iota| difference."""
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
    out["iota"] = float(np.max(np.abs(w1.iotaf - w2.iotaf)))
    return out


record = {"vmecpp": md.version("vmecpp")}
record["ns"], record["ftol"], record["threads"] = a.ns, a.ftol, a.threads

# 1. current-free, as given
vin = vmecpp.VmecInput.from_file(DECK)
p_free, t_free = run_vmecpp(vin, "current_free")
wf = Wout(p_free)
# 2. iota fit and the prescribed-iota deck
coef = np.polynomial.polynomial.polyfit(wf.s, wf.iotaf, DEGREE)
fit_err = float(np.max(np.abs(np.polynomial.polynomial.polyval(wf.s, coef) - wf.iotaf)))
text = DECK.read_text()
text = re.sub(r"NCURR\s*=\s*\d+", "NCURR = 0", text)
text = text.replace("PIOTA_TYPE = 'power_series'",
                    "PIOTA_TYPE = 'power_series'\n  AI = " +
                    " ".join(f"{c:.16e}" for c in coef))
(HERE / "input.landreman-paul-qa-iota").write_text(text)
vin = vmecpp.VmecInput.from_file(HERE / "input.landreman-paul-qa-iota")
p_iota, t_iota = run_vmecpp(vin, "iota")
wi = Wout(p_iota)
record.update(time_current_free_s=t_free, time_iota_s=t_iota, iota_fit_degree=DEGREE,
              iota_fit_max_error=fit_err, iota_axis=float(wf.iotaf[0]),
              iota_edge=float(wf.iotaf[-1]), mpol=wi.mpol, ntor=wi.ntor, ns_final=wi.ns,
              fsq_iota=wi.fsq, fsq_current_free=wf.fsq,
              iota_vs_current_free=compare(wi, wf))
# 3. VMEC2000
if a.xvmec:
    px_free, tx_free = run_xvmec(DECK.read_text(), "qa")
    px_iota, tx_iota = run_xvmec(text, "qa_iota")
    record.update(xvmec=str(a.xvmec), xvmec_time_current_free_s=tx_free,
                  xvmec_time_iota_s=tx_iota,
                  xvmec_vs_vmecpp_current_free=compare(Wout(px_free), wf),
                  xvmec_vs_vmecpp_iota=compare(Wout(px_iota), wi))

# 4. reference files of the prescribed-iota VMEC++ state
W = wi.energy()
mu0 = 4e-7 * np.pi
record.update(energy_T2m3=W, energy_J=W / mu0, wb=wi.wb,
              wb_times_4pi2=wi.wb * (2 * np.pi) ** 2, volume=wi.volume,
              phiedge=wi.phiedge, nfp=wi.nfp)
np.savetxt(OUT / "iota_coefficients.txt", coef,
           header=f"iota(s) = sum_k c_k s^k, k = 0..{DEGREE}; s normalized toroidal flux; "
                  f"PHIEDGE {wi.phiedge:.12e} Wb, NFP {wi.nfp}")
th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
rows = []
for k in range(NPHI):
    phi = 2 * np.pi * k / (NPHI * wi.nfp)
    r, z, _ = wi.surface(0.0, np.zeros(1), phi)
    rows.append((phi, r[0], z[0]))
np.savetxt(OUT / "axis.txt", rows, header="phi R Z (m), cylindrical phi over one period")
for s in SURFACES:
    rows = []
    for k in range(NPHI):
        phi = 2 * np.pi * k / (NPHI * wi.nfp)
        r, z, b = wi.surface(s, th, phi)
        rows += [(phi, t, ri, zi, bi) for t, ri, zi, bi in zip(th, r, z, b)]
    np.savetxt(OUT / f"surface_s{s:.2f}.txt", rows, fmt="%.10e",
               header=f"s = {s}: phi theta_VMEC R Z (m) |B| (T)")
np.savetxt(OUT / "iota.txt", np.c_[wi.s, wi.iotaf], header="s iota (VMEC++ full grid)")
(OUT / "summary.json").write_text(json.dumps(record, indent=1, default=float))
print(json.dumps(record, indent=1, default=float))
