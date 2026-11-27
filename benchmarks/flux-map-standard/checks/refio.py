"""Shared part of the third-code references (DESC, GVEC): the deck's
profiles as power series in s, the boundary of the deck, and the reference
layout of make_reference.py (axis, surfaces s = 0.25, 0.5, 1 on 64 planes
per period, iota, summary.json) for any code that evaluates R, Z, |B| at
(s, own poloidal angle, own toroidal angle).

The toroidal angle of a code may run against VMEC's cylindrical phi
(GVEC with hmap 1, possibly DESC). The orientation is taken from the
boundary: at s = 1 every code reproduces the deck's boundary, and only one
sign of the toroidal angle does so."""
import json
import math
import pathlib
import re

import numpy as np

NPHI, NTHETA, SURFACES = 64, 128, (0.25, 0.5, 1.0)


def _numbers(text, key):
    m = re.search(rf"(?im)^\s*{key}\s*=\s*([^\n!]*)", text)
    if not m:
        return None
    return [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", m.group(1).replace("d", "e").replace("D", "e"))]


def _word(text, key, default):
    m = re.search(rf"(?im)^\s*{key}\s*=\s*'([^']*)'", text)
    return m.group(1).lower() if m else default


class Deck:
    """The parts of a VMEC INDATA deck that fix the prescribed-iota problem:
    NFP, PHIEDGE, iota(s) and p(s) as power series in s, boundary RBC/ZBS."""

    def __init__(self, path):
        t = pathlib.Path(path).read_text()
        self.nfp = int(_numbers(t, "NFP")[0])
        self.phiedge = _numbers(t, "PHIEDGE")[0]
        if _word(t, "PIOTA_TYPE", "power_series") != "power_series":
            raise SystemExit("only PIOTA_TYPE = 'power_series'")
        nc = _numbers(t, "NCURR")
        if nc and int(nc[0]) != 0:
            raise SystemExit("deck must prescribe iota (NCURR = 0)")
        self.iota = np.array(_numbers(t, "AI"))
        am = np.array(_numbers(t, "AM") or [0.0])
        scale = (_numbers(t, "PRES_SCALE") or [1.0])[0]
        kind = _word(t, "PMASS_TYPE", "power_series")
        if kind == "power_series":
            self.pres = scale * am
        elif kind == "two_power":
            # p = AM0 (1 - s^AM1)^AM2, integer exponents expanded binomially
            a0, a1, a2 = am[:3]
            if a1 != int(a1) or a2 != int(a2):
                raise SystemExit("two_power with non-integer exponents")
            c = np.zeros(int(a1) * int(a2) + 1)
            for k in range(int(a2) + 1):
                c[k * int(a1)] = math.comb(int(a2), k) * (-1) ** k
            self.pres = scale * a0 * c
        else:
            raise SystemExit(f"PMASS_TYPE {kind} not supported")
        sp = _numbers(t, "SPRES_PED")
        if sp and sp[0] < 1.0:
            raise SystemExit("SPRES_PED < 1 not supported")
        self.rbc, self.zbs = {}, {}
        for key, tab in (("RBC", self.rbc), ("ZBS", self.zbs)):
            for n, m, v in re.findall(rf"(?i){key}\(\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*\)\s*=\s*([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)", t):
                tab[(int(m), int(n))] = float(v.replace("d", "e").replace("D", "e"))

    def iota_of(self, s):
        return np.polynomial.polynomial.polyval(s, self.iota)

    def pres_of(self, s):
        return np.polynomial.polynomial.polyval(s, self.pres)

    def boundary(self, theta, phi):
        """R, Z of the deck boundary at VMEC theta and cylindrical phi."""
        r = sum(v * np.cos(m * theta - n * self.nfp * phi) for (m, n), v in self.rbc.items())
        z = sum(v * np.sin(m * theta - n * self.nfp * phi) for (m, n), v in self.zbs.items())
        return r, z


def curve_distance(ref, pts, nsub=16):
    """Distance of pts (R, Z) to the closed curve ref (R, Z at equally spaced
    angle, trigonometric interpolant); dense search then Newton."""
    n = len(ref)
    c = np.fft.fft(ref, axis=0) / n
    k = np.fft.fftfreq(n, 1.0 / n)
    if n % 2 == 0:
        c[n // 2] *= 0.5
        c, k = np.vstack([c, c[n // 2:n // 2 + 1]]), np.append(k, n // 2)

    def ev(t, d=0):
        return np.real((np.exp(1j * np.outer(t, k)) * (1j * k) ** d) @ c)

    td = np.linspace(0, 2 * np.pi, n * nsub, endpoint=False)
    dense = ev(td)
    t = td[np.argmin(np.hypot(pts[:, None, 0] - dense[None, :, 0], pts[:, None, 1] - dense[None, :, 1]), axis=1)]
    for _ in range(8):
        x, x1, x2 = ev(t), ev(t, 1), ev(t, 2)
        r = x - pts
        h = np.sum(x1**2 + r * x2, axis=1)
        t = t - np.sum(r * x1, axis=1) / np.where(h > 0, h, np.inf)
    return np.hypot(*(ev(t) - pts).T)


def orientation(deck, surface):
    """+1 or -1: the sign s_z with code toroidal angle = s_z * phi, from the
    boundary. surface(s, theta, zeta) -> R, Z, B arrays over theta."""
    th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
    best = []
    for sz in (1, -1):
        d = 0.0
        for k in range(8):
            phi = 2 * np.pi * (k + 0.37) / (8 * deck.nfp)
            r, z, _ = surface(1.0, th, sz * phi)
            rb, zb = deck.boundary(th, phi)
            d = max(d, curve_distance(np.c_[rb, zb], np.c_[r, z]).max())
        best.append((d, sz))
    best.sort()
    print(f"boundary check: distance {best[0][0]:.2e} m (sign {best[0][1]:+d}), other sign {best[1][0]:.2e} m", flush=True)
    if best[0][0] > 1e-6:
        raise SystemExit("no orientation reproduces the deck boundary")
    return best[0][1], best[0][0]


def write_reference(out, deck, surface, axis, record, s_iota=None, iota=None):
    """Write the make_reference.py layout into out/reference."""
    ref = pathlib.Path(out) / "reference"
    ref.mkdir(parents=True, exist_ok=True)
    sz, dbnd = orientation(deck, surface)
    th = np.linspace(0, 2 * np.pi, NTHETA, endpoint=False)
    phis = 2 * np.pi * np.arange(NPHI) / (NPHI * deck.nfp)
    rows = []
    for phi in phis:
        r, z = axis(sz * phi)
        rows.append((phi, r, z))
    np.savetxt(ref / "axis.txt", rows, header="phi R Z (m), cylindrical phi over one period")
    for s in SURFACES:
        rows = []
        for phi in phis:
            r, z, b = surface(s, th, sz * phi)
            rows += [(phi, t, ri, zi, bi) for t, ri, zi, bi in zip(th, r, z, b)]
        np.savetxt(ref / f"surface_s{s:.2f}.txt", rows, fmt="%.10e",
                   header=f"s = {s}: phi theta R Z (m) |B| (T), theta the code's poloidal angle")
    if s_iota is not None:
        np.savetxt(ref / "iota.txt", np.c_[s_iota, iota], header="s iota")
    record.update(toroidal_sign=sz, boundary_distance_m=dbnd)
    (ref / "summary.json").write_text(json.dumps(record, indent=1, default=float))
    print(json.dumps(record, indent=1, default=float), flush=True)
