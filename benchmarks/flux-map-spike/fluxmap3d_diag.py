"""Throwaway spike: Cartesian flux map u=(x,y) in 3D, fixed boundary, p = 0, prescribed iota(s).
B = c0 [2 grad x X grad y + iota(s) grad zeta X grad s], s = x^2+y^2 (normalised toroidal flux),
c0 = PHIEDGE/(2 pi), zeta = cylindrical angle (spike only).  One field period, hex Q_p, Cartesian vectors.
usage: fluxmap3d.py DECK REFDIR p proj n nz [kappa list, e.g. 8,4,2,1] [-v]
"""
import sys, re, time, json
import numpy as np
import jax, jax.numpy as jnp
import scipy.sparse as sp, scipy.sparse.linalg as spla
jax.config.update("jax_enable_x64", True)

def read_deck(path):
    txt = open(path).read()
    nfp = int(re.search(r'NFP\s*=\s*(\d+)', txt).group(1))
    phiedge = float(re.search(r'PHIEDGE\s*=\s*([-\d.eE+]+)', txt).group(1))
    rbc = {(int(a), int(b)): float(c) for a, b, c in re.findall(r'RBC\(\s*(-?\d+),\s*(-?\d+)\)\s*=\s*([-\d.eE+]+)', txt)}
    zbs = {(int(a), int(b)): float(c) for a, b, c in re.findall(r'ZBS\(\s*(-?\d+),\s*(-?\d+)\)\s*=\s*([-\d.eE+]+)', txt)}
    return nfp, phiedge, rbc, zbs

def boundary(rbc, zbs, nfp, th, phi):
    R = sum(c*np.cos(m*th - n*nfp*phi) for (n, m), c in rbc.items())
    Z = sum(c*np.sin(m*th - n*nfp*phi) for (n, m), c in zbs.items())
    return R, Z

def disc_mesh(n, p):
    N = n*p; a = 0.45/np.sqrt(2); lat = np.linspace(-1, 1, N + 1)
    blocks = []
    X, Y = np.meshgrid(a*lat, a*lat, indexing='ij'); blocks.append(np.stack([X, Y], -1))
    T, S = np.meshgrid(np.linspace(0, 1, N + 1), lat, indexing='ij')
    for k in range(4):
        sq = np.stack([a*np.ones_like(S), a*S], -1)
        arc = np.stack([np.cos(S*np.pi/4), np.sin(S*np.pi/4)], -1)
        P = (1 - T[..., None])*sq + T[..., None]*arc
        c, s_ = np.cos(k*np.pi/2), np.sin(k*np.pi/2)
        blocks.append(np.stack([c*P[..., 0] - s_*P[..., 1], s_*P[..., 0] + c*P[..., 1]], -1))
    ids = {}; disc = []; conn = []; centre = []
    def nid(pt):
        key = (round(pt[0]*1e8), round(pt[1]*1e8))
        if key not in ids: ids[key] = len(disc); disc.append(pt)
        return ids[key]
    for ib, B in enumerate(blocks):
        gid = np.array([[nid(B[i, j]) for j in range(N + 1)] for i in range(N + 1)])
        for i in range(n):
            for j in range(n):
                conn.append(gid[i*p:(i + 1)*p + 1, j*p:(j + 1)*p + 1].ravel()); centre.append(ib == 0)
    return np.array(disc), np.array(conn), np.array(centre)

def section(disc, rbc, zbs, nfp, phi):
    """(R, Z) of the disc nodes on the plane phi: centre = m=0 part, inner ellipse = m=1 part of the boundary."""
    r = np.hypot(disc[:, 0], disc[:, 1]); th = np.arctan2(disc[:, 1], disc[:, 0])
    tg = 2*np.pi*np.arange(2048)/2048
    Rg, Zg = boundary(rbc, zbs, nfp, tg, phi)
    # area centroid; rays from it (the section must be star-shaped about it)
    dR = np.gradient(Rg, tg, edge_order=2); dZ = np.gradient(Zg, tg, edge_order=2)
    area = 0.5*np.sum(Rg*dZ - Zg*dR); c = np.array([np.sum(Rg**2*dZ)/(2*area), -np.sum(Zg**2*dR)/(2*area)])
    # affine normalisation by the second moments of the section, then rays from the centroid
    x_, z_ = Rg - c[0], Zg - c[1]
    Ixx = np.sum(x_**3*dZ)/3/area; Izz = -np.sum(z_**3*dR)/3/area; Ixz = np.sum(x_**2*z_*dZ)/2/area
    ev, U = np.linalg.eigh(np.array([[Ixx, Ixz], [Ixz, Izz]])); S = U @ np.diag(np.sqrt(ev)) @ U.T; Si = np.linalg.inv(S)
    xn = Si[0, 0]*x_ + Si[0, 1]*z_; zn = Si[1, 0]*x_ + Si[1, 1]*z_
    om = np.unwrap(np.arctan2(zn, xn)); rho = np.hypot(xn, zn)
    sgn = np.sign(om[-1] - om[0]); assert (sgn*np.diff(om) > 0).all(), 'section not star-shaped about its centroid'
    if sgn < 0: om = om[::-1]; rho = rho[::-1]
    rb = np.interp(th, om, rho, period=2*np.pi); rbar = rho.min()
    scale = rbar + r**2*(rb - rbar)
    return c[None, :] + (scale[:, None]*disc) @ S.T

def lagrange(p, xq):
    xn = np.linspace(-1, 1, p + 1); V = np.zeros((len(xq), p + 1)); D = np.zeros_like(V)
    for a_ in range(p + 1):
        c = np.zeros(p + 1); c[a_] = 1; poly = np.polyfit(xn, c, p)
        V[:, a_] = np.polyval(poly, xq); D[:, a_] = np.polyval(np.polyder(poly), xq)
    return V, D

def main():
    deck, refdir = sys.argv[1], sys.argv[2]
    p, proj, n, nz = [int(a) for a in sys.argv[3:7]]
    kappas = [float(a) for a in sys.argv[7].split(',')] if len(sys.argv) > 7 and not sys.argv[7].startswith('-') else [1.0]
    verbose = '-v' in sys.argv
    nfp, phiedge, rbc, zbs = read_deck(deck)
    c0 = phiedge/(2*np.pi)
    iota_c = np.loadtxt(refdir + '/iota_coefficients.txt')[::-1].copy()
    t00 = time.time()
    disc, conn2, centre2 = disc_mesh(n, p)
    nd = len(disc); K = nz*p; per = 2*np.pi/nfp
    rd = np.hypot(disc[:, 0], disc[:, 1]); isb2 = rd > 1 - 1e-9
    phis = per*np.arange(K + 1)/K
    xyz = np.zeros((K + 1, nd, 3))
    for k, ph in enumerate(phis):
        RZ = section(disc, rbc, zbs, nfp, ph)
        xyz[k] = np.stack([RZ[:, 0]*np.cos(ph), RZ[:, 0]*np.sin(ph), RZ[:, 1]], -1)
    # elements
    ne2, nen2 = conn2.shape; nen = nen2*(p + 1); nw = (p + 1)**2
    gconn = []; dnode = []
    for kz in range(nz):
        pl = kz*p + np.arange(p + 1)
        for e2 in range(ne2):
            gconn.append((pl[None, :]*nd + conn2[e2][:, None]).ravel())
            dnode.append(((pl[None, :] % K)*nd + conn2[e2][:, None]).ravel())
    gconn = np.array(gconn); dnode = np.array(dnode); ne = len(gconn)
    # tables
    nq = p + 2
    xq, wq = np.polynomial.legendre.leggauss(nq); V, D = lagrange(p, xq)
    e3 = lambda A, B, C: np.einsum('ia,jb,kc->ijkabc', A, B, C).reshape(nq**3, -1)
    Nq = e3(V, V, V); dN = np.stack([e3(D, V, V), e3(V, D, V), e3(V, V, D)], 0)      # (r, q, a)
    wref = np.einsum('i,j,k->ijk', wq, wq, wq).ravel()
    q1, q2, q3 = [a.ravel() for a in np.meshgrid(np.arange(nq), np.arange(nq), np.arange(nq), indexing='ij')]
    Nb, dNb = V[q1, p], D[q1, p]
    Lw = np.einsum('qb,qc->qbc', V[q2], V[q3]).reshape(nq**3, -1)
    d2L = np.einsum('qb,qc->qbc', D[q2], V[q3]).reshape(nq**3, -1)
    d3L = np.einsum('qb,qc->qbc', V[q2], D[q3]).reshape(nq**3, -1)
    ref = np.stack([xq[q1], xq[q2], xq[q3]], -1)
    def projector(d):
        if d is None: return np.eye(nq**3)
        Pm = np.stack([ref[:, 0]**i*ref[:, 1]**j*ref[:, 2]**k for i in range(d[0] + 1) for j in range(d[1] + 1) for k in range(d[2] + 1)], -1)
        return Pm @ np.linalg.solve(Pm.T @ (wref[:, None]*Pm), Pm.T*wref[None, :])
    q = p - 1
    degs = {1: [(q, q, q)]*3, 2: [None, None, (q, q, p)], 3: [(p, q, q), (q, p, q), (q, q, p)], 4: [None, None, (q, q, q)],
            5: [(q, q, q), (q, q, q), (q, q, p)], 6: [(p, q, p), (q, p, p), (q, q, p)]}.get(proj, [None]*3)
    Pi = np.stack([projector(d) for d in degs], 0)
    xe = xyz.reshape(-1, 3)[gconn]                                   # (ne, nen, 3)
    J = np.einsum('rqa,ead->eqdr', dN, xe)                           # dX_d/dxi_r
    detJ = np.linalg.det(J); print('# detJ min max', detJ.min(), detJ.max()); assert detJ.min()*detJ.max() > 0, 'folded mesh'
    A = np.linalg.inv(J)                                             # (r, d)
    G = np.einsum('rqa,eqrd->eqad', dN, A)
    pos = np.einsum('qa,ead->eqd', Nq, xe)
    R2 = pos[..., 0]**2 + pos[..., 1]**2
    gz = np.stack([-pos[..., 1]/R2, pos[..., 0]/R2, np.zeros_like(R2)], -1)
    wdet = wref[None, :]*np.abs(detJ)
    print(f'# mesh: {ne} elements, {nd} nodes/plane, {K} planes, volume x nfp {nfp*wdet.sum():.6f}, tables {time.time() - t00:.1f}s', flush=True)

    Nqj, Nbj, dNbj, Lwj, d2Lj, d3Lj, Pij = [jnp.array(a) for a in (Nq, Nb, dNb, Lw, d2L, d3L, Pi)]
    ic = jnp.array(iota_c)
    def bfield(z, G, A, J, detJ, gz, flag, kappa):
        ue = z[:2*nen].reshape(-1, 2); thn = z[2*nen:]
        th = Lwj @ thn; c, s_ = jnp.cos(th), jnp.sin(th)
        rg = jnp.stack([dNbj, Nbj*(d2Lj @ thn), Nbj*(d3Lj @ thn)], -1)      # ref gradient pieces
        gxw = jnp.stack([rg[:, 0]*c, -rg[:, 1]*s_, -rg[:, 2]*s_], -1)
        gyw = jnp.stack([rg[:, 0]*s_, rg[:, 1]*c, rg[:, 2]*c], -1)
        x = Nqj @ ue[:, 0] + flag*Nbj*c; y = Nqj @ ue[:, 1] + flag*Nbj*s_
        gx = jnp.einsum('qad,a->qd', G, ue[:, 0]) + flag*jnp.einsum('qr,qrd->qd', gxw, A)
        gy = jnp.einsum('qad,a->qd', G, ue[:, 1]) + flag*jnp.einsum('qr,qrd->qd', gyw, A)
        T = 2*jnp.cross(gx, gy)
        if proj:
            That = detJ[:, None]*jnp.einsum('qrd,qd->qr', A, T)
            T = jnp.einsum('qdr,qr->qd', J, jnp.einsum('rqk,kr->qr', Pij, That))/detJ[:, None]
        s = x*x + y*y; gs = 2*(x[:, None]*gx + y[:, None]*gy)
        return c0*(T + kappa*jnp.polyval(ic, s)[:, None]*jnp.cross(gz, gs)), s
    def energy(z, G, A, J, detJ, gz, w, flag, kappa):
        B, _ = bfield(z, G, A, J, detJ, gz, flag, kappa)
        return 0.5*jnp.sum(w*jnp.sum(B*B, -1))
    args = [jnp.array(a) for a in (G, A, J, detJ, gz, wdet)]
    E = jax.jit(jax.vmap(energy, in_axes=(0, 0, 0, 0, 0, 0, 0, 0, None)))
    dE = jax.jit(jax.vmap(jax.grad(energy), in_axes=(0, 0, 0, 0, 0, 0, 0, 0, None)))
    d2E = jax.jit(jax.vmap(jax.hessian(energy), in_axes=(0, 0, 0, 0, 0, 0, 0, 0, None)))

    # unknowns
    nn = K*nd; isb = np.tile(isb2, K)
    col = np.zeros((nn, 2), int); k_ = 0
    for i in range(nn):
        if isb[i]: col[i] = k_; k_ += 1
        else: col[i] = (k_, k_ + 1); k_ += 2
    nu = k_
    dd = np.tile(disc, (K, 1))
    qv = np.zeros(nu + 1)
    qv[col[isb, 0]] = np.arctan2(dd[isb, 1], dd[isb, 0])
    qv[col[~isb, 0]] = dd[~isb, 0]; qv[col[~isb, 1]] = dd[~isb, 1]
    nzl = 2*nen + nw
    dofs = np.full((ne, nzl), nu); off = np.zeros((ne, nzl)); flag = np.zeros(ne)
    bl = np.arange(p*nw, (p + 1)*nw)
    for e in range(ne):
        nodes = dnode[e]; ib = isb[nodes]
        dofs[e, 0:2*nen:2] = np.where(ib, nu, col[nodes, 0]); dofs[e, 1:2*nen:2] = np.where(ib, nu, col[nodes, 1])
        if ib.any():
            assert ib[bl].all() and ib.sum() == nw
            flag[e] = 1; th = qv[col[nodes[bl], 0]]
            dofs[e, 2*nen:] = col[nodes[bl], 0]; off[e, 2*nen:] = 2*np.pi*np.round((th[nw//2] - th)/(2*np.pi))
    flagj = jnp.array(flag)
    # gauge: stellarator symmetry gives y = 0 on Z = 0 of the plane phi = 0
    fixed = np.zeros(nu, bool)
    for j in range(nd):
        if abs(disc[j, 1]) < 1e-12:
            if isb2[j]: fixed[col[j, 0]] = True
            else: fixed[col[j, 1]] = True; qv[col[j, 1]] = 0
    free = np.where(~fixed)[0]
    rows = np.repeat(dofs, nzl, axis=1).ravel(); cols = np.tile(dofs, (1, nzl)).ravel()
    zloc = lambda qv: jnp.array(qv[dofs] + off)
    print(f'# unknowns {nu}', flush=True)

    def axis_error(qv):
        refax = np.loadtxt(refdir + '/axis.txt'); err = []; sgn = []
        xq1 = np.linspace(-1, 1, p + 1)
        for k in range(0, K, max(1, K//8)):
            best = None
            for e2 in np.where(centre2)[0]:
                nodes = k*nd + conn2[e2]
                ue = np.stack([qv[col[nodes, 0]], qv[col[nodes, 1]]], -1)      # centre block: interior nodes
                xi = np.zeros(2)
                for _ in range(20):
                    V1, D1 = lagrange(p, xi[:1]); V2, D2 = lagrange(p, xi[1:])
                    N_ = np.outer(V1[0], V2[0]).ravel()
                    dN_ = np.stack([np.outer(D1[0], V2[0]).ravel(), np.outer(V1[0], D2[0]).ravel()], -1)
                    f = N_ @ ue; Jm = ue.T @ dN_
                    try: xi = xi - np.linalg.solve(Jm, f)
                    except np.linalg.LinAlgError: break
                    if np.abs(xi).max() > 3: break
                if np.abs(xi).max() <= 1 + 1e-9 and np.linalg.norm(f) < 1e-10:
                    V1, _ = lagrange(p, xi[:1]); V2, _ = lagrange(p, xi[1:])
                    X = np.outer(V1[0], V2[0]).ravel() @ xyz[k][conn2[e2]]
                    best = (np.hypot(X[0], X[1]), X[2]); break
            if best is None: err.append(np.nan); continue
            def trig(f):
                c = np.fft.rfft(f)/len(f); m = np.arange(len(c))
                return (c[0] + 2*np.sum(c[1:]*np.exp(1j*m[1:]*nfp*phis[k]))).real
            err.append(np.hypot(best[0] - trig(refax[:, 1]), best[1] - trig(refax[:, 2])))
            sgn.append((phis[k], best[0] - trig(refax[:, 1]), best[1] - trig(refax[:, 2])))
        print('   axis dR dZ per plane [mm]: ' + ' '.join(f'({a:.2f}: {1e3*b:+.1f} {1e3*c:+.1f})' for a, b, c in sgn))
        return np.array(err)

    def label_on_reference(qv):
        # interior elements only (wall elements skipped): s_h at VMEC surface points of the plane phi = 0
        out = []
        for sv in ('0.25', '0.50'):
            d = np.loadtxt(refdir + f'/surface_s{sv}.txt'); d = d[np.abs(d[:, 0]) < 1e-12][::4]
            vals = []; bvals = []
            for Rp, Zp, Bp in d[:, 2:5]:
                for e2 in range(ne2):
                    nodes = conn2[e2]
                    if isb2[nodes].any(): continue
                    xe2 = np.stack([np.hypot(xyz[0][nodes, 0], xyz[0][nodes, 1]), xyz[0][nodes, 2]], -1)
                    if not (xe2[:, 0].min() - 1e-9 <= Rp <= xe2[:, 0].max() + 1e-9 and xe2[:, 1].min() - 1e-9 <= Zp <= xe2[:, 1].max() + 1e-9): continue
                    xi = np.zeros(2); ok = False
                    for _ in range(15):
                        V1, D1 = lagrange(p, xi[:1]); V2, D2 = lagrange(p, xi[1:])
                        N_ = np.outer(V1[0], V2[0]).ravel()
                        dN_ = np.stack([np.outer(D1[0], V2[0]).ravel(), np.outer(V1[0], D2[0]).ravel()], -1)
                        f = N_ @ xe2 - np.array([Rp, Zp])
                        if np.linalg.norm(f) < 1e-12: ok = True; break
                        xi = xi - np.linalg.solve(xe2.T @ dN_, f)
                        if np.abs(xi).max() > 2: break
                    if ok and np.abs(xi).max() <= 1 + 1e-9:
                        vals.append((N_ @ qv[col[nodes, 0]])**2 + (N_ @ qv[col[nodes, 1]])**2)
                        # |B| at the point: 3D element of the first layer, xi3 = -1
                        e = e2; gn = gconn[e]; dn = dnode[e]
                        V3, D3 = lagrange(p, np.array([-1.0]))
                        N3 = np.einsum('a,b,c->abc', V1[0], V2[0], V3[0]).ravel()
                        dN3 = np.stack([np.einsum('a,b,c->abc', D1[0], V2[0], V3[0]).ravel(), np.einsum('a,b,c->abc', V1[0], D2[0], V3[0]).ravel(),
                                        np.einsum('a,b,c->abc', V1[0], V2[0], D3[0]).ravel()], 0)
                        X3 = xyz.reshape(-1, 3)[gn]; J3 = (dN3 @ X3).T; G3 = np.linalg.solve(J3.T, dN3)       # (d, a)
                        xn_, yn_ = qv[col[dn, 0]], qv[col[dn, 1]]
                        x_, y_ = N3 @ xn_, N3 @ yn_; gx_, gy_ = G3 @ xn_, G3 @ yn_
                        s_ = x_*x_ + y_*y_; gs_ = 2*(x_*gx_ + y_*gy_)
                        Bv = c0*(2*np.cross(gx_, gy_) + np.polyval(iota_c, s_)*np.cross(np.array([0.0, 1.0/Rp, 0.0]), gs_))
                        bvals.append(np.linalg.norm(Bv)/Bp - 1); break
            vals = np.array(vals)
            out.append(f's={sv}: {len(vals)} pts, mean s_h {vals.mean():.4f}, rms dev {np.sqrt(np.mean((vals - float(sv))**2)):.2e}, |B| rel dev rms {np.sqrt(np.mean(np.array(bvals)**2)):.2e} max {np.abs(bvals).max():.2e}' if len(vals) else f's={sv}: no pts')
        print('   label on VMEC surfaces (phi = 0): ' + '; '.join(out))

    refW = json.load(open(refdir + '/summary.json')).get('energy_T2m3', 0.281794)
    for kappa in kappas:
        total = lambda qv: float(jnp.sum(E(zloc(qv), *args, flagj, kappa)))
        def derivs(qv):
            z = zloc(qv)
            ge = np.array(dE(z, *args, flagj, kappa)); He = np.array(d2E(z, *args, flagj, kappa))
            g = np.zeros(nu + 1); np.add.at(g, dofs.ravel(), ge.ravel())
            return g, sp.csr_matrix((He.ravel(), (rows, cols)), shape=(nu + 1, nu + 1))
        lam = 1e-9; W = total(qv); steps = 0; t0 = time.time(); W0 = W
        g, H = derivs(qv)
        for it in range(60):
            gf = g[free]; Hf = H[free][:, free]; dg = np.abs(Hf.diagonal()) + 1e-30
            gn = np.linalg.norm(gf/np.sqrt(dg))
            if verbose: print(f'   it {it:3d} W {nfp*W:.10e} |g| {gn:.2e} lam {lam:.1e}', flush=True)
            if gn < 1e-9*np.sqrt(abs(W)): break
            ok = False
            while lam < 1e8:
                d = spla.spsolve((Hf + lam*dg.min()*sp.identity(len(free)) ).tocsc(), -gf)
                pred = -(gf @ d + 0.5*d @ (Hf @ d))
                qt = qv.copy(); qt[free] += d; Wt = total(qt)
                if pred > 0 and np.isfinite(Wt) and W - Wt > 1e-4*pred - 1e-13*abs(W):
                    qv = qt; W = Wt; steps += 1; lam = max(lam/100, 1e-12); ok = True; break
                lam *= 10
            if not ok: break
            g, H = derivs(qv)
        ax = axis_error(qv); label_on_reference(qv)
        print(f'kappa {kappa:5.2f} p {p} proj {proj} n {n} nz {nz} dofs {nu} steps {steps:3d} |g| {gn:.1e} '
              f'W {nfp*W:.6f} (start {nfp*W0:.4f}; VMEC at kappa 1: {refW:.6f}) axis err max {np.nanmax(ax):.2e} '
              f'mean {np.nanmean(ax):.2e} nan {int(np.isnan(ax).sum())}  {time.time() - t0:.1f}s', flush=True)

main()
