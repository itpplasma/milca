"""Lowest eigenvalues of the flux-map Hessian at a solved state.

Reads the stream files written by test_flux_map_3d with MILCA_FM3_HESSIAN
set, removes the fixed slots (y of the wall functions) and reports the lowest
eigenvalues of the Jacobi-scaled Hessian D^-1/2 H D^-1/2. The scaling is a
congruence, so the number of negative eigenvalues (the inertia) is that of H.
For each mode it reports the overlap with the rigid rotation of u (the exact
symmetry x + i y -> e^{i c}(x + i y), Theta -> Theta + c) and with uniform
interior shifts of x and y (a proxy for a displacement of the axis, the zero
of u), and how the mode is spread over the planes.

Usage: python hessian_spectrum.py FILE [FILE ...] [--k K] [--dense-max N]
"""
import argparse
import sys

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def read(path):
    with open(path, "rb") as f:
        nrow, nnz, n0, nt, nfix, ngauge = np.fromfile(f, np.int32, 6)
        row_ptr = np.fromfile(f, np.int32, nrow + 1) - 1
        col = np.fromfile(f, np.int32, nnz) - 1
        val = np.fromfile(f, np.float64, nnz)
        z = np.fromfile(f, np.float64, nrow)
        g = np.fromfile(f, np.float64, nrow)
        fixed = np.fromfile(f, np.int32, nfix) - 1
        wall = np.fromfile(f, np.int32, n0).astype(bool)
        w = np.fromfile(f, np.float64, 1)[0]
        mass = np.fromfile(f, np.float64, n0 * n0 * nt).reshape(nt, n0, n0)
        base = np.fromfile(f, np.float64, 2 * n0 * nt).reshape(nt, n0, 2)
    h = sp.csr_matrix((val, col, row_ptr), shape=(nrow, nrow))
    blocks = [sp.csr_matrix(mass[k]) for k in range(nt)]
    m = sp.block_diag(blocks + blocks, format="csr")
    return dict(h=h, m=m, z=z, g=g, fixed=fixed, wall=wall, n0=n0, nt=nt, w=w,
                ngauge=ngauge, volume=mass.sum(), base=base)


def generators(d):
    n0, nt, z, wall = d["n0"], d["nt"], d["z"], d["wall"]
    x = z[: n0 * nt].reshape(nt, n0)
    y = z[n0 * nt:].reshape(nt, n0)
    wall2 = np.broadcast_to(wall, (nt, n0))
    rot = np.concatenate([np.where(wall2, 1.0, -y).ravel(), np.where(wall2, 0.0, x).ravel()])
    one = np.where(wall2, 0.0, 1.0).ravel()
    tx = np.concatenate([one, 0 * one])
    ty = np.concatenate([0 * one, one])
    return dict(rotation=rot, shift_x=tx, shift_y=ty)


def analyse_l2(d, path, k):
    """Generalized problem H v = mu M v with the L2 metric M of delta u:
    eigenvalues mu that compare between meshes."""
    h, m = d["h"], d["m"]
    free = np.ones(h.shape[0], bool)
    free[d["fixed"]] = False
    hf, mf = h[free][:, free], m[free][:, free]
    hf = 0.5 * (hf + hf.T)
    lam, vec = spla.eigsh(hf.tocsc(), k=k, M=mf.tocsc(), sigma=-1.0, which="LM")
    order = np.argsort(lam)
    print(f"  L2 metric (volume {d['volume']:.6f}): lowest mu " + " ".join(f"{x:.4e}" for x in lam[order]))
    return lam[order]


def analyse(path, k, dense_max):
    d = read(path)
    h = d["h"]
    n = h.shape[0]
    asym = abs(h - h.T).max() / abs(h).max()
    free = np.ones(n, bool)
    free[d["fixed"]] = False
    hf = h[free][:, free]
    dg = np.abs(hf.diagonal())
    s = sp.diags(1 / np.sqrt(dg)) @ hf @ sp.diags(1 / np.sqrt(dg))
    s = 0.5 * (s + s.T)
    gscaled = np.linalg.norm(d["g"][free] / np.sqrt(dg))
    m = s.shape[0]
    if m <= dense_max:
        lam, vec = sla.eigh(s.toarray(), subset_by_index=[0, k - 1])
        method = "dense"
    else:
        lam, vec = spla.eigsh(s.tocsc(), k=k, sigma=-1e-2, which="LM")
        order = np.argsort(lam)
        lam, vec = lam[order], vec[:, order]
        method = "shift-invert about -1e-2"
    gens = {name: np.sqrt(dg) * v[free] for name, v in generators(d).items()}
    gens = {name: v / np.linalg.norm(v) for name, v in gens.items()}
    nt, n0 = d["nt"], d["n0"]
    print(f"{path}: {m} free unknowns, W {d['w']:.9e}, scaled |g| {gscaled:.1e}, "
          f"gauge rows {d['ngauge']}, Hessian asymmetry {asym:.1e}, {method}")
    print("   eigenvalue   rotation  shift_x  shift_y   plane weights (rms per plane / max)")
    full = np.zeros(n)
    for j in range(k):
        full[:] = 0
        full[free] = vec[:, j] / np.sqrt(dg)
        per_plane = np.sqrt((full[: n0 * nt].reshape(nt, n0) ** 2).sum(1)
                            + (full[n0 * nt:].reshape(nt, n0) ** 2).sum(1))
        per_plane /= per_plane.max()
        ov = [abs(vec[:, j] @ gens[name]) for name in ("rotation", "shift_x", "shift_y")]
        print(f"  {lam[j]: .4e}  {ov[0]:8.3f} {ov[1]:8.3f} {ov[2]:8.3f}   "
              + " ".join(f"{p:.2f}" for p in per_plane))
    print(f"  negative among these: {(lam < 0).sum()}")
    analyse_l2(d, path, k)
    return lam


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--dense-max", type=int, default=6000)
    a = ap.parse_args()
    for path in a.files:
        analyse(path, a.k, a.dense_max)
    sys.exit(0)
