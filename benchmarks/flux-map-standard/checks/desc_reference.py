"""DESC reference for a prescribed-iota deck, as a resolution ladder.

Conventions (boundary, Psi, sign of iota) come from DESC's own VMEC reader
applied to a VMEC++ wout of the deck; the profiles are then replaced by the
deck's exact power series (s = rho^2), the interior is reset to DESC's
initial guess from the boundary, and DESC solves by continuation. Each level
raises (L, M, N) and re-solves from the previous level; each writes
OUT/L<l>/reference in the layout of make_reference.py.

Usage: python desc_reference.py DECK WOUT OUT L:M:N [L:M:N ...]
       (python with desc-opt 0.17.4)"""
import json
import pathlib
import sys
import time

import numpy as np
import desc
from desc.continuation import solve_continuation_automatic
from desc.grid import LinearGrid
from desc.profiles import PowerSeriesProfile
from desc.vmec import VMECIO

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import refio  # noqa: E402

deck_path, wout, out = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3])
levels = [tuple(int(v) for v in x.split(":")) for x in sys.argv[4:]]
deck = refio.Deck(deck_path)

L0, M0, N0 = levels[0]
eq = VMECIO.load(wout, L=L0, M=M0, N=N0)
sign = float(np.sign(eq.iota(np.array([np.sqrt(0.5)]))[0] / deck.iota_of(0.5)))
print(f"desc {desc.__version__}; iota sign relative to the deck {sign:+.0f}; Psi {float(eq.Psi):.6e}", flush=True)
iota = PowerSeriesProfile(params=sign * deck.iota, modes=2 * np.arange(len(deck.iota)), sym=True)
pres = PowerSeriesProfile(params=deck.pres, modes=2 * np.arange(len(deck.pres)), sym=True)
eq.iota, eq.pressure = iota, pres
eq.set_initial_guess()  # interior from the boundary only, not from VMEC++

TOL = dict(ftol=1e-12, xtol=1e-12, gtol=1e-12)
for i, (L, M, N) in enumerate(levels):
    t0 = time.perf_counter()
    if i == 0:
        eq = solve_continuation_automatic(eq, verbose=1, maxiter=500, **TOL)[-1]
    else:
        eq.change_resolution(L=L, M=M, N=N, L_grid=2 * L, M_grid=2 * M, N_grid=2 * N)
        eq.solve(verbose=1, maxiter=500, **TOL)
    dt = time.perf_counter() - t0
    f = eq.compute(["<|F|>_vol", "<|grad(|B|^2)|/2mu0>_vol", "W", "V"])
    fnorm = float(f["<|F|>_vol"]) / float(f["<|grad(|B|^2)|/2mu0>_vol"])

    def surface(s, th, zeta, eq=eq):
        g = LinearGrid(rho=np.array([np.sqrt(s)]), theta=th, zeta=np.array([zeta]), NFP=1, sym=False)
        d = eq.compute(["R", "Z", "|B|"], grid=g)
        return np.asarray(d["R"]), np.asarray(d["Z"]), np.asarray(d["|B|"])

    def axis(zeta, eq=eq):
        r, z, _ = surface(0.0, np.zeros(1), zeta)
        return r[0], z[0]

    s = np.linspace(0, 1, 101)
    gi = LinearGrid(rho=np.sqrt(s), M=0, N=0, NFP=eq.NFP)
    io = np.asarray(eq.compute("iota", grid=gi)["iota"])
    record = {"code": "DESC", "version": desc.__version__, "deck": str(deck_path), "L": L, "M": M, "N": N,
              "time_s": dt, "force_normalised_vol": fnorm, "iota_sign": sign,
              "iota_vs_deck": float(np.max(np.abs(sign * io - deck.iota_of(s)))),
              "energy_T2m3": float(f["W"]) * 4e-7 * np.pi, "volume": float(f["V"])}
    print(f"L M N {L} {M} {N}: {dt:.0f} s, normalised force {fnorm:.2e}", flush=True)
    refio.write_reference(out / f"L{i}", deck, surface, axis, record, s, sign * io)
    eq.save(str(out / f"L{i}" / "eq.h5"))
