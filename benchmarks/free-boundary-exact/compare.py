"""Errors of a free-boundary result against the exact Solov'ev reference.

Usage: compare.py REFERENCE.json RESULT.nc [label]

RESULT is a VMEC wout file or MILCA's nested_equilibrium output. Reported:
relative errors of volume, wb, betatotal, iota at the edge, the axis
position, and for a wout the largest distance of the s = 1 surface from the
exact boundary psi = psi_b (|psi - psi_b| / |grad psi|).
"""
import json
import sys

import netCDF4
import numpy as np


def main():
    ref = json.load(open(sys.argv[1]))
    d = netCDF4.Dataset(sys.argv[2])
    label = sys.argv[3] if len(sys.argv) > 3 else sys.argv[2]
    v = d.variables

    def get(name):
        if name in v:
            return float(np.ravel(v[name][:])[0] if v[name].ndim else v[name][:])
        return None

    out = {}
    out["volume"] = get("volume_p")
    out["wb"] = get("wb")
    out["beta"] = get("betatotal")
    if "iotaf" in v:
        io = np.ravel(v["iotaf"][:])
        out["iota_edge"] = abs(float(io[-1]))
    else:
        ie = get("iotaf_edge")
        out["iota_edge"] = abs(ie) if ie is not None else None
    out["raxis"] = float(np.ravel(v["raxis_cc"][:])[0]) if "raxis_cc" in v else None
    errs = {
        "volume": rel(out["volume"], ref["volume_p"]),
        "wb": rel(out["wb"], ref["wb"]),
        "beta": rel(out["beta"], ref["betatotal"]),
        "iota_edge": rel(out["iota_edge"], ref["iota_edge"]),
        "axis_m": abs(out["raxis"] - ref["raxis_cc"]) if out["raxis"] is not None else None,
    }
    if "rmnc" in v:
        errs["boundary_m"] = boundary_distance(ref, v)
    line = f"{label:28s}" + "".join(
        f" {k} {val:9.2e}" if val is not None else f" {k} {'-':>9s}" for k, val in errs.items())
    print(line)
    if "ier_flag" in v and int(v["ier_flag"][:]) != 0:
        print(f"  ier_flag {int(v['ier_flag'][:])}")


def rel(a, b):
    return None if a is None else abs(a - b) / abs(b)


def boundary_distance(ref, v):
    rmnc = np.array(v["rmnc"][:])[-1]
    zmns = np.array(v["zmns"][:])[-1]
    xm = np.array(v["xm"][:])
    th = np.linspace(0, 2 * np.pi, 721)
    r = (rmnc[None, :] * np.cos(np.outer(th, xm))).sum(1)
    z = (zmns[None, :] * np.sin(np.outer(th, xm))).sum(1)
    a, b, r0 = ref["alpha"], ref["beta"], ref["raxis_cc"]
    psi = a * (r**2 - r0**2) ** 2 + b * r**2 * z**2
    dr = 4 * a * r * (r**2 - r0**2) + 2 * b * r * z**2
    dz = 2 * b * r**2 * z
    return float(np.max(np.abs(psi - ref["psi_b"]) / np.hypot(dr, dz)))


if __name__ == "__main__":
    main()
