"""GVEC reference for a prescribed-iota deck at one resolution.

Conventions (boundary, axis guess, PHIEDGE sign, sign of iota) come from
GVEC's own VMEC converter applied to a VMEC++ wout of the deck; the profiles
are then replaced by the deck's exact power series in s and the start is
GVEC's initial guess from boundary and axis, not the VMEC++ interior.
Writes OUT/reference in the layout of make_reference.py.

Usage: python gvec_reference.py DECK WOUT OUT M N NELEMS [DEG] [TOL]
       (python with gvec 1.5.0)"""
import pathlib
import sys
import time

import numpy as np
import gvec
from gvec.scripts.convert_wout import convert_vmec_wout

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import refio  # noqa: E402

deck_path, wout, out = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
M, N, nel = (int(v) for v in sys.argv[4:7])
deg = int(sys.argv[7]) if len(sys.argv) > 7 else 5
tol = float(sys.argv[8]) if len(sys.argv) > 8 else 1e-10
deck = refio.Deck(deck_path)
out.mkdir(parents=True, exist_ok=True)

p = convert_vmec_wout(wout.absolute())
conv = p["iota"]
sign = float(np.sign(np.interp(0.5, conv["rho2"], conv["vals"]) / deck.iota_of(0.5)))
print(f"gvec {gvec.__version__}; iota sign relative to the deck {sign:+.0f}; phiedge {p['phiedge']}", flush=True)
p["iota"] = {"type": "polynomial", "coefs": list(sign * deck.iota)}
p["pres"] = {"type": "polynomial", "coefs": list(deck.pres)}
p.update(ProjectName="ref", whichInitEquilibrium=0, X1_mn_max=[M, N], X2_mn_max=[M, N], LA_mn_max=[M, N],
         X1X2_deg=deg, LA_deg=deg, sgrid={"nElems": nel}, minimize_tol=tol, maxIter=400000, totalIter=400000)
for k in ("VMECwoutfile", "VMECwoutfile_format", "init_LA"):
    p.pop(k, None)

t0 = time.perf_counter()
run = gvec.run(p, runpath=out / "run", keep_intermediates=None)
dt = time.perf_counter() - t0
st = run.state


def surface(s, th, zeta):
    d = st.evaluate("X1", "X2", "mod_B", rho=np.array([np.sqrt(max(s, 0.0))]), theta=th, zeta=np.array([zeta]))
    f = lambda q: d[q].transpose("rad", "pol", "tor").values[0, :, 0]
    return f("X1"), f("X2"), f("mod_B")


def axis(zeta):
    r, z, _ = surface(0.0, np.zeros(1), zeta)
    return r[0], z[0]


s = np.linspace(0, 1, 101)
ev = st.evaluate("iota", rho=np.sqrt(s), theta=None, zeta=None)
io = sign * ev.iota.values
record = {"code": "GVEC", "version": gvec.__version__, "deck": str(deck_path), "M": M, "N": N, "nElems": nel,
          "deg": deg, "minimize_tol": tol, "time_s": dt, "max_force": float(run.max_force),
          "iterations": int(run.GVEC_iter_used), "iota_sign": sign,
          "iota_vs_deck": float(np.max(np.abs(io - deck.iota_of(s))))}
refio.write_reference(out, deck, surface, axis, record, s, io)
