"""Export a symmetric VMEC Wout in the label-qa comparison layout."""
import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np

NPHI, NTHETA = 64, 128


def export_wout(path, out):
    """Write geometry and field samples and return Wout energy and residuals."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path) as f:
        if "ier_flag" in f.variables and int(f["ier_flag"][:]) != 0:
            raise ValueError("Wout is from an unsuccessful solve")
        if "rmns" in f.variables or "zmnc" in f.variables:
            raise ValueError("only stellarator-symmetric Wout is supported")
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
    np.savetxt(out / "axis.txt", rows, header="phi R Z")
    for s in (0.25, 0.5):
        rows = []
        for k in range(NPHI):
            phi = 2 * np.pi * k / (NPHI * nfp)
            r, z, b = surface(s, th, phi)
            rows += [(phi, th[i], r[i], z[i], b[i]) for i in range(NTHETA)]
        np.savetxt(out / f"surface_s{s:.2f}.txt", rows, header="phi angle R Z |B|")

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
    return dict(energy_T2m3=w * (2 * np.pi) ** 2 / (nsf - 1), fsq=fsq, iterations=itfsq,
                unknowns=3 * nsf * len(xm), mnmax=len(xm), signgs=signgs,
                jacobian_range=[float(jmin), float(jmax)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wout", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    summary = export_wout(args.wout, args.out)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1))
