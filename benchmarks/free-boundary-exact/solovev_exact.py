"""Exact axisymmetric free-boundary benchmark from the Solov'ev equilibrium.

Inside the plasma D = {psi < psi_b} the field is the Solov'ev equilibrium

    psi = alpha (R^2 - r0^2)^2 + beta R^2 Z^2,   B = grad psi x grad phi + F grad phi,

with F constant, so the plasma current is toroidal, mu0 J_phi = -(8 alpha +
2 beta) R, and the pressure p = (8 alpha + 2 beta)(psi_b - psi)/mu0 vanishes on
the boundary. The external field is B_ext = B_Solovev - B_plasma, B_plasma the
field of J_phi in D by the ring-current Green's function. B_ext is exact on and
inside D, so every free-boundary code given B_ext (one mgrid coil group at unit
current), the toroidal flux and the profiles p(s), iota(s) must return D, the
axis R = r0, Z = 0 and the integrals below. Outside D, B_ext is the field of a
distributed coil (the Solov'ev current continued), which a free-boundary code
never sees except through B_ext.

B_plasma at each grid node is integrated in polar coordinates about the node,
which removes the 1/r singularity; along each ray psi is a quartic, so the
chord of D is found from its exact roots. Gauss-Legendre in angle and radius.

Writes: mgrid_solovev_exact.nc, input.solovev_exact (VMEC INDATA with exact
spline profiles), reference.json (exact values) and the boundary curve.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import ellipe, ellipk

MU0 = 4e-7 * np.pi


class Solovev:
    def __init__(self, r0=4.0, kappa=1.5, a=1.0, b0=1.5, beta_target=0.02):
        self.r0, self.f = r0, r0 * b0
        # psi ~ 4 alpha r0^2 dR^2 + beta r0^2 Z^2 near the axis: elongation
        # kappa = 2 sqrt(alpha/beta); psi_b from the outboard midplane at a
        self.alpha = 1.0
        self.beta = 4 * self.alpha / kappa**2
        self.psi_b = self.psi(r0 + a, 0.0)
        # scale psi so that the volume-average beta is near beta_target
        self._scale_to_beta(b0, beta_target)

    def psi(self, r, z):
        return self.alpha * (r**2 - self.r0**2) ** 2 + self.beta * r**2 * z**2

    def grad(self, r, z):
        dr = 4 * self.alpha * r * (r**2 - self.r0**2) + 2 * self.beta * r * z**2
        dz = 2 * self.beta * r**2 * z
        return dr, dz

    def pprime(self):
        return -(8 * self.alpha + 2 * self.beta) / MU0  # dp/dpsi

    def pressure(self, psi):
        return self.pprime() * (psi - self.psi_b)

    def jphi(self, r):
        return -(8 * self.alpha + 2 * self.beta) * r / MU0

    def field(self, r, z):
        dr, dz = self.grad(r, z)
        return -dz / r, self.f / r, dr / r  # B_R, B_phi, B_Z

    def _scale_to_beta(self, b0, beta_target):
        # p0 = -p' psi_b scales like c^2 with psi -> c psi; B_phi fixed
        p0 = -self.pprime() * self.psi_b
        c2 = beta_target * b0**2 / (2 * MU0) / (0.5 * p0)
        c = np.sqrt(c2)
        self.alpha *= c
        self.beta *= c
        self.psi_b *= c

    def ray_chord(self, x0, z0, ang):
        """Parameters t1 < t2 >= 0 where the ray x0 + t e enters and leaves D."""
        e_r, e_z = np.cos(ang), np.sin(ang)
        # psi(x0 + t e) - psi_b as a quartic in t
        r = np.polynomial.Polynomial([x0, e_r])
        z = np.polynomial.Polynomial([z0, e_z])
        q = self.alpha * (r**2 - self.r0**2) ** 2 + self.beta * r**2 * z**2 - self.psi_b
        roots = q.roots()
        roots = np.sort(roots[np.abs(roots.imag) < 1e-12].real)
        inside0 = self.psi(x0, z0) < self.psi_b
        if inside0:
            pos = roots[roots > 0]
            return (0.0, pos[0]) if len(pos) else None
        pos = roots[roots > 0]
        # D is convex for these parameters: the chord is the first pair
        for k in range(len(pos) - 1):
            tm = 0.5 * (pos[k] + pos[k + 1])
            if self.psi(x0 + tm * e_r, z0 + tm * e_z) < self.psi_b:
                return pos[k], pos[k + 1]
        return None


def ring_field(r, z, rs, zs, current):
    """(B_R, B_Z) at (r, z) of coaxial rings at (rs, zs) carrying current."""
    dz = z - zs
    s2 = (r + rs) ** 2 + dz**2
    d2 = (r - rs) ** 2 + dz**2
    m = 4 * r * rs / s2
    kk, ee = ellipk(m), ellipe(m)
    pref = MU0 * current / (2 * np.pi * np.sqrt(s2))
    br = pref * dz / r * (-kk + (rs**2 + r**2 + dz**2) / d2 * ee)
    bz = pref * (kk + (rs**2 - r**2 - dz**2) / d2 * ee)
    return br, bz


def quartic_coefficients(sol, x0, z0, er, ez):
    """Coefficients c0..c4 (arrays over rays) of psi(x0 + t e) - psi_b."""
    # r = x0 + t er, z = z0 + t ez; u = r^2 - r0^2 = (x0^2 - r0^2) + 2 x0 er t + er^2 t^2
    u0, u1, u2 = x0**2 - sol.r0**2, 2 * x0 * er, er**2
    r2 = (x0**2, 2 * x0 * er, er**2)
    z2 = (z0**2, 2 * z0 * ez, ez**2)
    c = np.zeros((5,) + np.shape(er))
    # alpha u^2
    for i, ui in enumerate((u0, u1, u2)):
        for j, uj in enumerate((u0, u1, u2)):
            c[i + j] += sol.alpha * ui * uj
    # beta r^2 z^2
    for i, ri in enumerate(r2):
        for j, zj in enumerate(z2):
            c[i + j] += sol.beta * ri * zj
    c[0] -= sol.psi_b
    return c


def batched_roots(c):
    """Real roots (nan elsewhere) of c0 + c1 t + ... + c4 t^4 per column."""
    n = c.shape[1]
    comp = np.zeros((n, 4, 4))
    comp[:, 1:, :3] = np.eye(3)
    comp[:, :, 3] = -(c[:4] / c[4]).T
    ev = np.linalg.eigvals(comp)
    real = np.where(np.abs(ev.imag) < 1e-9 * (1 + np.abs(ev.real)), ev.real, np.nan)
    return np.sort(real, axis=1)


def plasma_field(sol, x0, z0, nang, nrad, panels=16):
    """B_plasma (R, Z) at (x0, z0): polar Gauss rule about the point (the
    area element t dt removes the 1/t singularity of the ring field)."""
    ga, wa = np.polynomial.legendre.leggauss(nang)
    gr, wr = np.polynomial.legendre.leggauss(nrad)
    a0 = 2 * np.pi * np.arange(panels)[:, None] / panels
    angs = (a0 + np.pi / panels * (ga[None, :] + 1)).ravel()
    wang = np.tile(np.pi / panels * wa, panels)
    er, ez = np.cos(angs), np.sin(angs)
    roots = batched_roots(quartic_coefficients(sol, x0, z0, er, ez))
    inside0 = sol.psi(x0, z0) < sol.psi_b
    t1 = np.full(angs.shape, np.nan)
    t2 = np.full(angs.shape, np.nan)
    for k in range(len(angs)):
        rt = roots[k]
        pos = rt[np.isfinite(rt) & (rt > 0)]
        if inside0:
            if len(pos):
                t1[k], t2[k] = 0.0, pos[0]
        else:
            for m in range(len(pos) - 1):
                tm = 0.5 * (pos[m] + pos[m + 1])
                if sol.psi(x0 + tm * er[k], z0 + tm * ez[k]) < sol.psi_b:
                    t1[k], t2[k] = pos[m], pos[m + 1]
                    break
    ok = np.isfinite(t1)
    t1, t2, er, ez, wang = t1[ok], t2[ok], er[ok], ez[ok], wang[ok]
    t = 0.5 * (t2 - t1)[:, None] * gr[None, :] + 0.5 * (t2 + t1)[:, None]
    wt = 0.5 * (t2 - t1)[:, None] * wr[None, :]
    rs = x0 + t * er[:, None]
    zs = z0 + t * ez[:, None]
    cur = sol.jphi(rs) * t * wt * wang[:, None]
    fr, fz = ring_field(x0, z0, rs, zs, cur)
    return fr.sum(), fz.sum()


def level_points(sol, level, npts):
    """Points on psi = level by rays from the axis."""
    th = 2 * np.pi * (np.arange(npts) + 0.5) / npts
    er, ez = np.cos(th), np.sin(th)
    c = quartic_coefficients(sol, sol.r0, 0.0, er, ez)
    c[0] += sol.psi_b - level
    roots = batched_roots(c)
    rho = np.array([rt[np.isfinite(rt) & (rt > 0)].min() for rt in roots])
    return sol.r0 + rho * er, rho * ez


def fit_coils(sol, args):
    """Ring currents on an ellipse outside the grid box fitted to B_ext."""
    ncoil = args.coils
    th = 2 * np.pi * (np.arange(ncoil) + 0.5) / ncoil
    coils_r = sol.r0 + args.coil_a * np.cos(th)
    coils_z = args.coil_b * np.sin(th)

    def exact(levels, npts):
        pts_r, pts_z = [], []
        for lv in levels:
            r, z = level_points(sol, lv * sol.psi_b, npts)
            pts_r.append(r)
            pts_z.append(z)
        pts_r, pts_z = np.concatenate(pts_r), np.concatenate(pts_z)
        b = np.array([plasma_field(sol, r, z, args.nang, args.nrad) for r, z in zip(pts_r, pts_z)])
        sr, _, sz = sol.field(pts_r, pts_z)
        return pts_r, pts_z, np.concatenate([sr - b[:, 0], sz - b[:, 1]])

    def matrix(pr, pz):
        m = np.zeros((2 * len(pr), ncoil))
        for k in range(ncoil):
            fr, fz = ring_field(pr, pz, coils_r[k], coils_z[k], 1.0)
            m[:len(pr), k] = fr
            m[len(pr):, k] = fz
        return m

    pr, pz, rhs = exact([0.05, 0.3, 0.6, 0.9, 0.999], 2 * ncoil)
    a = matrix(pr, pz)
    cur, *_ = np.linalg.lstsq(a, rhs, rcond=None)
    scale = np.max(np.abs(rhs))
    fit = dict(fit=np.max(np.abs(a @ cur - rhs)) / scale)
    pr, pz, rhs = exact([0.15, 0.45, 0.75], 2 * ncoil + 7)
    fit["check"] = np.max(np.abs(matrix(pr, pz) @ cur - rhs)) / scale
    pr, pz, rhs = exact([0.9995], 4 * ncoil + 3)
    fit["boundary"] = np.max(np.abs(matrix(pr, pz) @ cur - rhs)) / scale
    return coils_r, coils_z, cur, fit


def write_coils(path, rs, zs, cur):
    """The rings as (R, Z, I), one per line."""
    with open(path, "w") as f:
        f.write("# R Z I of the coaxial ring coils (A, +phi)\n")
        for r, z, i in zip(rs, zs, cur):
            f.write(f"{r:.16e} {z:.16e} {i:.16e}\n")


def flux_profiles(sol, nsurf, ntheta=512):
    """Toroidal flux Phi(psi), iota(psi) on psi_k = psi_b s-ish levels, and
    the volume, wb, beta of D."""
    # flux surfaces by rays from the axis: radius rho(theta, psi)
    th = 2 * np.pi * (np.arange(ntheta) + 0.5) / ntheta
    psis = sol.psi_b * (np.arange(1, nsurf + 1) / nsurf) ** 2

    def radius(level, ang):
        r = np.polynomial.Polynomial([sol.r0, np.cos(ang)])
        z = np.polynomial.Polynomial([0.0, np.sin(ang)])
        q = sol.alpha * (r**2 - sol.r0**2) ** 2 + sol.beta * r**2 * z**2 - level
        rt = q.roots()
        rt = rt[(np.abs(rt.imag) < 1e-10) & (rt.real > 0)].real
        return rt.min()

    phi = np.zeros(nsurf)
    q_safety = np.zeros(nsurf)
    for k, lev in enumerate(psis):
        rho = np.array([radius(lev, a) for a in th])
        rr = sol.r0 + rho * np.cos(th)
        zz = rho * np.sin(th)
        # Phi = F int dA/R = F int_0^2pi int_0^rho (1/R) rho' drho' dth: inner by Gauss
        g, w = np.polynomial.legendre.leggauss(32)
        acc = 0.0
        for j, a in enumerate(th):
            rp = 0.5 * rho[j] * (g + 1)
            acc += np.sum(w * 0.5 * rho[j] * rp / (sol.r0 + rp * np.cos(a)))
        phi[k] = sol.f * acc * 2 * np.pi / ntheta
        # q = F/(2 pi) oint dl/(R |grad psi|) on the surface
        drho = np.gradient(rho, th, edge_order=2)  # periodic enough at 512 points
        dr = drho * np.cos(th) - rho * np.sin(th)
        dz = drho * np.sin(th) + rho * np.cos(th)
        dl = np.hypot(dr, dz)
        gr_, gz_ = sol.grad(rr, zz)
        q_safety[k] = sol.f / (2 * np.pi) * np.sum(dl / (rr**2 * np.hypot(gr_, gz_) / rr)) * \
            2 * np.pi / ntheta
    return psis, phi, q_safety


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--grid", type=int, default=81, help="mgrid nodes per direction")
    ap.add_argument("--nang", type=int, default=24)
    ap.add_argument("--nrad", type=int, default=24)
    ap.add_argument("--coils", type=int, default=48)
    ap.add_argument("--beta", type=float, default=0.02, help="target volume-average beta")
    ap.add_argument("--b0", type=float, default=1.5, help="vacuum B_phi at r0 (T)")
    ap.add_argument("--coil-a", type=float, default=2.2)
    ap.add_argument("--coil-b", type=float, default=2.8)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    sol = Solovev(beta_target=args.beta, b0=args.b0)

    # ring coils on an ellipse outside the grid box, currents by least squares
    # against the exact B_ext = B_Solovev - B_plasma on levels inside D
    coils_r, coils_z, coils_i, fit = fit_coils(sol, args)
    print(f"coil fit: {len(coils_r)} rings, max |B_fit - B_ext| / max |B_pol| on independent "
          f"points: fit {fit['fit']:.2e}, check {fit['check']:.2e}, boundary {fit['boundary']:.2e}")
    rmin, rmax = sol.r0 - 1.5, sol.r0 + 1.5
    zmin, zmax = -2.0, 2.0
    n = args.grid
    rg = np.linspace(rmin, rmax, n)
    zg = np.linspace(zmin, zmax, n)
    rr, zz = np.meshgrid(rg, zg)
    br = np.zeros((n, n))
    bz = np.zeros((n, n))
    for rc, zc, ic in zip(coils_r, coils_z, coils_i):
        fr, fz = ring_field(rr, zz, rc, zc, ic)
        br += fr
        bz += fz
    bp = sol.f / rr
    write_mgrid(out / "mgrid_solovev_exact.nc", rg, zg, br, bp, bz)
    write_coils(out / "coils.solovev_exact", coils_r, coils_z, coils_i)

    psis, phi, qs = flux_profiles(sol, 64)
    phi_b = phi[-1]
    s = phi / phi_b
    iota = 1 / qs
    pres = sol.pressure(psis)
    # axis values by extrapolation in s (smooth)
    s_all = np.concatenate([[0.0], s])
    iota0 = np.polyval(np.polyfit(s[:6], iota[:6], 3), 0.0)
    iota_all = np.concatenate([[iota0], iota])
    p_all = np.concatenate([[sol.pressure(0.0)], pres])
    write_indata(out / "input.solovev_exact", sol, phi_b, s_all, p_all, iota_all)

    # exact integrals over D: volume, wb (VMEC: int B^2/(2 mu0) dV / (2pi)^2 * mu0?), beta
    ref = exact_integrals(sol)
    ref.update(raxis_cc=sol.r0, zaxis=0.0, phiedge=phi_b, iota_axis=iota0, iota_edge=iota[-1],
               psi_b=sol.psi_b, alpha=sol.alpha, beta=sol.beta, F=sol.f)
    (out / "reference.json").write_text(json.dumps(ref, indent=1))
    print(json.dumps(ref, indent=1))


def exact_integrals(sol, ntheta=1024, nrad=48):
    th = 2 * np.pi * (np.arange(ntheta) + 0.5) / ntheta
    g, w = np.polynomial.legendre.leggauss(nrad)
    vol = wmag = wp = area = 0.0
    for a in th:
        ch = sol.ray_chord(sol.r0, 0.0, a)
        rho = ch[1]
        rp = 0.5 * rho * (g + 1)
        wt = 0.5 * rho * w * rp * (2 * np.pi / ntheta)  # dA
        rr = sol.r0 + rp * np.cos(a)
        zz = rp * np.sin(a)
        b_r, b_p, b_z = sol.field(rr, zz)
        dv = 2 * np.pi * rr * wt
        vol += dv.sum()
        area += wt.sum()
        wmag += np.sum((b_r**2 + b_p**2 + b_z**2) / (2 * MU0) * dv)
        wp += np.sum(sol.pressure(sol.psi(rr, zz)) * dv)
    return dict(volume_p=vol, area=area, wb=wmag / (2 * np.pi) ** 2 * MU0,
                wp=wp / (2 * np.pi) ** 2 * MU0, betatotal=wp / wmag)


def write_mgrid(path, rg, zg, br, bp, bz):
    import netCDF4
    ds = netCDF4.Dataset(path, "w")
    ds.createDimension("stringsize", 30)
    ds.createDimension("external_coil_groups", 1)
    ds.createDimension("dim_00001", 1)
    ds.createDimension("external_coils", 1)
    ds.createDimension("rad", len(rg))
    ds.createDimension("zee", len(zg))
    ds.createDimension("phi", 1)
    for name, val in [("ir", len(rg)), ("jz", len(zg)), ("kp", 1), ("nfp", 1), ("nextcur", 1)]:
        v = ds.createVariable(name, "i4")
        v.assignValue(val)
    for name, val in [("rmin", rg[0]), ("rmax", rg[-1]), ("zmin", zg[0]), ("zmax", zg[-1])]:
        v = ds.createVariable(name, "f8")
        v.assignValue(val)
    cg = ds.createVariable("coil_group", "S1", ("external_coil_groups", "stringsize"))
    cg[0, :] = netCDF4.stringtochar(np.array(["solovev_exact".ljust(30)], "S30"))[0]
    mm = ds.createVariable("mgrid_mode", "S1", ("dim_00001",))
    mm[0] = b"S"
    rc = ds.createVariable("raw_coil_cur", "f8", ("external_coils",))
    rc[0] = 1.0
    for name, arr in [("br_001", br), ("bp_001", bp), ("bz_001", bz)]:
        v = ds.createVariable(name, "f8", ("phi", "zee", "rad"))
        v[0, :, :] = arr
    ds.close()


def write_indata(path, sol, phi_b, s, p, iota):
    def arr(name, vals):
        return f"  {name} = " + " ".join(f"{v:.16e}" for v in vals) + "\n"

    # initial boundary guess: an ellipse (the free boundary must find D)
    a_guess = 0.9
    txt = "&INDATA\n"
    txt += "  LFREEB = T, MGRID_FILE = 'mgrid_solovev_exact.nc', EXTCUR(1) = 1.0,\n"
    txt += "  NFP = 1, MPOL = 8, NTOR = 0, NZETA = 1, LASYM = F,\n"
    txt += f"  PHIEDGE = {phi_b:.16e},\n  NCURR = 0, GAMMA = 0.0, CURTOR = 0.0,\n"
    txt += "  NS_ARRAY = 16 32 64, FTOL_ARRAY = 1e-10 1e-12 1e-14, NITER_ARRAY = 5000 10000 20000,\n"
    txt += "  NSTEP = 200, DELT = 0.9, NVACSKIP = 6,\n"
    txt += "  PMASS_TYPE = 'cubic_spline', PRES_SCALE = 1.0,\n"
    txt += arr("AM_AUX_S", s) + arr("AM_AUX_F", p)
    txt += "  PIOTA_TYPE = 'cubic_spline',\n"
    txt += arr("AI_AUX_S", s) + arr("AI_AUX_F", iota)
    txt += f"  RAXIS_CC = {sol.r0:.6f}, ZAXIS_CS = 0.0,\n"
    txt += f"  RBC(0,0) = {sol.r0:.6f}, RBC(0,1) = {a_guess:.6f}, ZBS(0,1) = {1.4 * a_guess:.6f},\n"
    txt += "/\n"
    path.write_text(txt)


if __name__ == "__main__":
    main()
