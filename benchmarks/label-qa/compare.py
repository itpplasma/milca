"""Compare a MILCA label-lane result with the VMEC++ reference in reference/
(or in the directory given as second argument).

MILCA writes, in a directory DIR, the layout of reference/:
  axis.txt             phi R Z           on the 64 reference planes
  surface_s0.25.txt    phi theta R Z |B| (any theta, the same planes)
  surface_s0.50.txt, surface_s1.00.txt
  iota.txt             s iota            (optional)
  summary.json         {"energy_T2m3": W} (optional, W = int |B|^2/2 dV, B in T)

Printed: axis max |dR|, |dZ|; per surface the distance of each MILCA point to
the reference curve of its plane (max and rms, m), the relative |B| difference
at the closest reference point (max and rms); max |iota difference|; relative
energy difference.

Usage: python compare.py DIR
"""
import json
import pathlib
import sys

import numpy as np

REF = pathlib.Path(__file__).resolve().parent / "reference"
SURFACES = ("0.25", "0.50", "1.00")


def planes(rows):
    return {p: rows[np.isclose(rows[:, 0], p, rtol=0, atol=1e-7)]
            for p in np.unique(np.round(rows[:, 0], 7))}


def plane_of(table, phi):
    """The entry of table whose plane angle is phi (to 1e-6 rad)."""
    key = min(table, key=lambda q: abs(q - phi))
    if abs(key - phi) > 1e-6:
        raise KeyError(f"no reference plane at phi = {phi}")
    return table[key]


def curve_distance(ref, pts, nsub=16):
    """Distance of pts (R, Z, |B|) to the closed reference curve ref (R, Z,
    |B| at equally spaced VMEC theta, trigonometric interpolant), and |B| of
    ref at the closest point: dense search, then Newton on the squared
    distance in the curve parameter."""
    n = len(ref)
    c = np.fft.fft(ref, axis=0) / n
    k = np.fft.fftfreq(n, 1.0 / n)
    if n % 2 == 0:
        c[n // 2] *= 0.5                     # split the Nyquist mode
        c, k = np.vstack([c, c[n // 2:n // 2 + 1]]), np.append(k, n // 2)

    def ev(t, d=0):
        e = np.exp(1j * np.outer(t, k)) * (1j * k) ** d
        return np.real(e @ c)

    td = np.linspace(0, 2 * np.pi, n * nsub, endpoint=False)
    dense = ev(td)
    d = np.hypot(pts[:, None, 0] - dense[None, :, 0], pts[:, None, 1] - dense[None, :, 1])
    t = td[np.argmin(d, axis=1)]
    for _ in range(8):
        x, x1, x2 = ev(t), ev(t, 1), ev(t, 2)
        r = x[:, :2] - pts[:, :2]
        g = np.sum(r * x1[:, :2], axis=1)
        h = np.sum(x1[:, :2] ** 2 + r * x2[:, :2], axis=1)
        t = t - g / np.where(h > 0, h, np.inf)
    x = ev(t)
    return np.hypot(*(x[:, :2] - pts[:, :2]).T), x[:, 2]


def read_summary(directory):
    path = directory / "summary.json"
    if not path.exists():
        return {}
    result = json.loads(path.read_text())
    if result.get("converged") is False:
        raise SystemExit(f"{directory}: equilibrium solve did not converge")
    return result


def main(out):
    out = pathlib.Path(out)
    summary, reference = read_summary(out), read_summary(REF)
    ref_axis, axis = np.loadtxt(REF / "axis.txt"), np.loadtxt(out / "axis.txt")
    ra = {p: (r, z) for p, r, z in ref_axis}
    dr = max(abs(r - plane_of(ra, p)[0]) for p, r, z in axis)
    dz = max(abs(z - plane_of(ra, p)[1]) for p, r, z in axis)
    print(f"axis: max |dR| {dr:.3e} m, max |dZ| {dz:.3e} m")
    for s in SURFACES:
        path = out / f"surface_s{s}.txt"
        if not path.exists():
            continue
        ref, mil = planes(np.loadtxt(REF / f"surface_s{s}.txt")), planes(np.loadtxt(path))
        dist, db = [], []
        for p, rows in mil.items():
            d, bref = curve_distance(plane_of(ref, p)[:, 2:5], rows[:, 2:5])
            dist.append(d)
            db.append(np.abs(rows[:, 4] - bref) / bref)
        dist, db = np.concatenate(dist), np.concatenate(db)
        print(f"s = {s}: distance max {dist.max():.3e} m, rms {np.sqrt(np.mean(dist**2)):.3e} m;"
              f" |B| rel. max {db.max():.3e}, rms {np.sqrt(np.mean(db**2)):.3e}")
    if (out / "iota.txt").exists():
        mi, ri = np.loadtxt(out / "iota.txt"), np.loadtxt(REF / "iota.txt")
        print(f"iota: max |diff| {np.max(np.abs(mi[:, 1] - np.interp(mi[:, 0], *ri.T))):.3e}")
    w, wr = summary.get("energy_T2m3"), reference.get("energy_T2m3")
    if w is not None and wr is not None:
        print(f"energy: {w:.10e} T^2 m^3, relative difference {(w - wr) / wr:.3e}")


if __name__ == "__main__":
    # optional second argument: another reference directory in the same
    # layout (e.g. benchmarks/flux-map-standard/<case>/reference)
    if len(sys.argv) > 2:
        REF = pathlib.Path(sys.argv[2]).resolve()
    main(sys.argv[1])
