"""Throwaway spike: Cartesian flux map u=(x,y): Omega -> unit disc, axisymmetric Solov'ev.
B_phi = tau(s) det(grad u), B_pol = psib grad(s) x grad(phi), s = x^2+y^2 = psi_N.
W = int [B_phi^2/2 + B_pol^2/2 - p(s)] R dR dZ.  Fixed O-grid mesh NOT centred on the axis.
usage: fluxmap.py f0 p proj(0/1) n1 n2 ...
"""
import sys, time
import numpy as np
import jax, jax.numpy as jnp
import scipy.sparse as sp, scipy.sparse.linalg as spla
jax.config.update("jax_enable_x64", True)

R0, AL, BE, PSIB = 1.0, 0.25, 0.5, 0.09
EPS = 1.0
PP = 8*AL + 2*BE                      # p = PP*PSIB*(1-s)

def psi(R, Z): return AL*(R**2 - R0**2)**2 + BE*R**2*Z**2
def dpsi(R, Z): return 4*AL*R*(R**2 - R0**2) + 2*BE*R*Z**2, 2*BE*R**2*Z

def ray_radius(c, th, level, rmax=0.9):
    lo = np.zeros_like(th); hi = np.full_like(th, rmax)
    for _ in range(60):
        mid = 0.5*(lo + hi)
        v = psi(c[0] + mid*np.cos(th), c[1] + mid*np.sin(th)) - level
        lo = np.where(v < 0, mid, lo); hi = np.where(v < 0, hi, mid)
    return 0.5*(lo + hi)

def qhat_fit(deg=12):
    """q/f0 as a polynomial in t = 2 s/1.1 - 1."""
    t = np.cos(np.pi*(np.arange(48) + 0.5)/48); s = 0.55*(t + 1)
    th = 2*np.pi*np.arange(400)/400
    q = []
    for sv in s:
        rho = ray_radius((R0, 0.0), th, sv*PSIB)
        R = R0 + rho*np.cos(th); Z = rho*np.sin(th)
        gr, gz = dpsi(R, Z)
        q.append(np.mean(rho/(R*(gr*np.cos(th) + gz*np.sin(th)))))
    c = np.polynomial.chebyshev.chebfit(t, np.array(q), deg)
    res = np.max(np.abs(np.polynomial.chebyshev.chebval(t, c) - q))
    return np.polynomial.chebyshev.cheb2poly(c)[::-1], res   # polyval order in t

def mesh(n, p, centre):
    """5-block O-grid of the unit disc mapped to the Solov'ev section; Q_p Lagrange nodes."""
    N = n*p; a = 0.45/np.sqrt(2)
    lat = np.linspace(-1, 1, N + 1)
    blocks = []
    X, Y = np.meshgrid(a*lat, a*lat, indexing='ij'); blocks.append(np.stack([X, Y], -1))
    T, S = np.meshgrid(np.linspace(0, 1, N + 1), lat, indexing='ij')
    for k in range(4):
        sq = np.stack([a*np.ones_like(S), a*S], -1)
        arc = np.stack([np.cos(S*np.pi/4), np.sin(S*np.pi/4)], -1)
        P = (1 - T[..., None])*sq + T[..., None]*arc
        c, s_ = np.cos(k*np.pi/2), np.sin(k*np.pi/2)
        blocks.append(np.stack([c*P[..., 0] - s_*P[..., 1], s_*P[..., 0] + c*P[..., 1]], -1))
    ids = {}; disc = []; conn = []
    def nid(pt):
        key = (round(pt[0]*1e8), round(pt[1]*1e8))
        if key not in ids: ids[key] = len(disc); disc.append(pt)
        return ids[key]
    for B in blocks:
        gid = np.array([[nid(B[i, j]) for j in range(N + 1)] for i in range(N + 1)])
        for i in range(n):
            for j in range(n):
                conn.append(gid[i*p:(i + 1)*p + 1, j*p:(j + 1)*p + 1].ravel())
    disc = np.array(disc); conn = np.array(conn)
    r = np.hypot(disc[:, 0], disc[:, 1]); th = np.arctan2(disc[:, 1], disc[:, 0])
    rb = ray_radius(centre, th, PSIB); rbar = ray_radius(centre, 2*np.pi*np.arange(64)/64, PSIB).min()
    scale = rbar + r**2*(rb - rbar)
    xy = np.stack([centre[0] + scale*disc[:, 0], centre[1] + scale*disc[:, 1]], -1)
    bnd = np.where(r > 1 - 1e-9)[0]
    return disc, xy, conn, bnd

def lagrange(p, xq):
    xn = np.linspace(-1, 1, p + 1)
    V = np.zeros((len(xq), p + 1)); D = np.zeros_like(V)
    for a_ in range(p + 1):
        c = np.zeros(p + 1); c[a_] = 1
        poly = np.polyfit(xn, c, p)
        V[:, a_] = np.polyval(poly, xq); D[:, a_] = np.polyval(np.polyder(poly), xq)
    return V, D

def tabulate(xy, conn, p, nq):
    xq, wq = np.polynomial.legendre.leggauss(nq)
    V, D = lagrange(p, xq)
    Nq = np.einsum('ia,jb->ijab', V, V).reshape(nq*nq, -1)
    d1 = np.einsum('ia,jb->ijab', D, V).reshape(nq*nq, -1)
    d2 = np.einsum('ia,jb->ijab', V, D).reshape(nq*nq, -1)
    w = np.outer(wq, wq).ravel()
    ref = np.stack(np.meshgrid(xq, xq, indexing='ij'), -1).reshape(-1, 2)
    xe = xy[conn]                                              # (ne, nen, 2)
    J = np.stack([np.einsum('qa,ead->eqd', d1, xe), np.einsum('qa,ead->eqd', d2, xe)], -2)  # (ne,nq,ref,phys)
    det = J[..., 0, 0]*J[..., 1, 1] - J[..., 0, 1]*J[..., 1, 0]
    assert det.min() > 0, 'folded mesh %g' % det.min()
    Ji = np.linalg.inv(J)                                      # (phys, ref)
    dref = np.stack([d1, d2], 0)                               # (ref, q, a)
    G = np.einsum('eqdr,rqa->eqad', Ji, dref)                  # dN_a/dx_d
    pos = np.einsum('qa,ead->eqd', Nq, xe)
    i1 = np.repeat(np.arange(nq), nq); i2 = np.tile(np.arange(nq), nq)
    tab = dict(Nb=V[i1, p], dNb=D[i1, p], L=V[i2], dL=D[i2])
    return Nq, G, pos, w[None, :]*det, ref, Ji, tab

def make_energy(Nq, ref, qpoly, f0, p, proj, tab):
    Nq = jnp.array(Nq); qpoly = jnp.array(qpoly)
    if proj >= 2: terms = [ref[:, 0]**i*ref[:, 1]**j for i in range(p) for j in range(p)]   # Q_{p-1}, dA measure
    else: terms = [ref[:, 0]**i*ref[:, 1]**j for i in range(p) for j in range(p - i)]   # P_{p-1}
    Pm = jnp.array(np.stack(terms, -1))
    Nb, dNb, L, dL = [jnp.array(tab[k]) for k in ('Nb', 'dNb', 'L', 'dL')]
    nen = Nq.shape[1]
    def fields(z, G, Ji, flag):
        ue = z[:2*nen].reshape(-1, 2); th = L @ z[2*nen:]; dth = dL @ z[2*nen:]
        c, s_ = jnp.cos(th), jnp.sin(th)
        x = Nq @ ue[:, 0] + flag*Nb*c; y = Nq @ ue[:, 1] + flag*Nb*s_
        gx = jnp.einsum('qad,a->qd', G, ue[:, 0]) + flag*(Ji[:, :, 0]*(dNb*c)[:, None] + Ji[:, :, 1]*(-Nb*s_*dth)[:, None])
        gy = jnp.einsum('qad,a->qd', G, ue[:, 1]) + flag*(Ji[:, :, 0]*(dNb*s_)[:, None] + Ji[:, :, 1]*(Nb*c*dth)[:, None])
        det = gx[:, 0]*gy[:, 1] - gx[:, 1]*gy[:, 0]
        s = x*x + y*y
        gs = 2*(x[:, None]*gx + y[:, None]*gy)
        tau = 2*f0*PSIB*jnp.polyval(qpoly, 2*s/1.1 - 1)
        return s, gs, tau*det
    def bbar(bphi, wR):
        M = Pm.T @ (wR[:, None]*Pm)
        return Pm @ jnp.linalg.solve(M, Pm.T @ (wR*bphi))
    def energy(z, G, R, w, Ji, flag):
        s, gs, bphi = fields(z, G, Ji, flag)
        wR = w*R
        hg = 0.0
        if proj:
            bb = bbar(bphi, w if proj >= 2 else wR)
            if proj == 3: hg = 0.5*(EPS/f0**2)*(bphi - bb)**2     # stabilisation at the poloidal scale
            bphi = bb
        return jnp.sum(wR*(hg + 0.5*bphi**2 + 0.5*PSIB**2*jnp.sum(gs**2, -1)/R**2 - PP*PSIB*(1 - s)))
    def diag(z, G, R, w, Ji, flag):
        s, gs, bphi = fields(z, G, Ji, flag)
        return s, gs, bphi, bbar(bphi, w*R)
    return (jax.jit(jax.vmap(energy)), jax.jit(jax.vmap(jax.grad(energy))),
            jax.jit(jax.vmap(jax.hessian(energy))), jax.jit(jax.vmap(diag)))

def run(f0s, p, proj, n, qpoly, verbose=False, pin=False, lm=True, guard=False):
    centre = (0.5*(np.sqrt(1 - 0.6) + np.sqrt(1 + 0.6)), 0.0)
    disc, xy, conn, bnd = mesh(n, p, centre)
    Nq, G, pos, w, ref, Ji, tab = tabulate(xy, conn, p, p + 3)
    G = jnp.array(G); R = jnp.array(pos[..., 0]); wj = jnp.array(w); Ji = jnp.array(Ji)
    nn = len(xy); ne, nen = conn.shape
    isb = np.zeros(nn, bool); isb[bnd] = True
    col = np.zeros((nn, 2), int); k = 0
    for i in range(nn):
        if isb[i]: col[i] = k; k += 1
        else: col[i] = (k, k + 1); k += 2
    nu = k
    qv = np.zeros(nu + 1)                       # last entry: dummy
    for i in range(nn):
        if isb[i]: qv[col[i, 0]] = np.arctan2(disc[i, 1], disc[i, 0])
        else: qv[col[i]] = disc[i]
    # element dof map: 2 nen nodal entries, then p+1 wall angles (with 2 pi offsets)
    nz = 2*nen + p + 1
    dofs = np.full((ne, nz), nu); off = np.zeros((ne, nz)); flag = np.zeros(ne)
    bl = np.arange(p*(p + 1), (p + 1)**2)       # local wall nodes of a wall element
    for e in range(ne):
        wall = isb[conn[e, bl]].all()
        for a_ in range(nen):
            i = conn[e, a_]
            if not isb[i]: dofs[e, 2*a_:2*a_ + 2] = col[i]
            else: assert wall and a_ in bl
        if wall:
            flag[e] = 1; th = qv[col[conn[e, bl], 0]]
            dofs[e, 2*nen:] = col[conn[e, bl], 0]
            off[e, 2*nen:] = 2*np.pi*np.round((th[p//2] - th)/(2*np.pi))
    flagj = jnp.array(flag)
    fixed = np.zeros(nu, bool)
    if pin:   # gauge by up-down symmetry: y = 0 on Z = 0 (interior nodes)
        for i in range(nn):
            if abs(xy[i, 1]) < 1e-12 and not isb[i]: fixed[col[i, 1]] = True; qv[col[i, 1]] = 0
            if abs(xy[i, 1]) < 1e-12 and isb[i]: fixed[col[i, 0]] = True
    free = np.where(~fixed)[0]
    rows = np.repeat(dofs, nz, axis=1).ravel(); cols = np.tile(dofs, (1, nz)).ravel()
    for f0 in f0s:
        E, dE, d2E, diag = make_energy(Nq, ref, qpoly, f0, p, proj, tab)
        zloc = lambda qv: jnp.array(qv[dofs] + off)
        def total(qv): return float(jnp.sum(E(zloc(qv), G, R, wj, Ji, flagj)))
        def derivs(qv):
            z = zloc(qv)
            ge = np.array(dE(z, G, R, wj, Ji, flagj)); He = np.array(d2E(z, G, R, wj, Ji, flagj))
            g = np.zeros(nu + 1); np.add.at(g, dofs.ravel(), ge.ravel())
            H = sp.csr_matrix((He.ravel(), (rows, cols)), shape=(nu + 1, nu + 1))
            return g, H
        lam = 1e-2 if lm else 1e-9; W = total(qv); steps = 0; t0 = time.time()
        g, H = derivs(qv)
        for it in range(400):
            gf = g[free]; Hf = H[free][:, free]; D = sp.diags(np.abs(Hf.diagonal()) + 1e-30) if lm else sp.identity(len(free))*np.abs(Hf.diagonal()).min()
            gn = np.linalg.norm(gf/np.sqrt(np.abs(Hf.diagonal()) + 1e-30))
            if verbose: print(f'   it {it:3d} W {W:.12e} |g| {gn:.2e} lam {lam:.1e}')
            if gn < 1e-9*max(1.0, f0): break
            while True:
                d = spla.spsolve((Hf + lam*D).tocsc(), -gf)
                pred = -(gf @ d + 0.5*d @ (Hf @ d))
                qt = qv.copy(); qt[free] += d
                Wt = total(qt)
                if guard and float(jnp.min(diag(zloc(qt), G, R, wj, Ji, flagj)[2])) <= 0: Wt = np.inf
                if pred > 0 and np.isfinite(Wt) and W - Wt > 1e-4*pred - 1e-14*abs(W):
                    rho = (W - Wt)/pred
                    qv = qt; W = Wt; steps += 1
                    lam = max(lam/10 if rho > 0.75 else (lam if rho > 0.25 else lam*4), 1e-12) if lm else max(lam/100, 1e-12)
                    break
                lam *= 10
                if lam > 1e12: break
            if lam > 1e12: break
            g, H = derivs(qv)
        # errors at the quadrature points
        s, gs, bphi, bb = [np.array(a) for a in diag(zloc(qv), G, R, wj, Ji, flagj)]
        Rq, Zq = pos[..., 0], pos[..., 1]
        sx = psi(Rq, Zq)/PSIB; gr, gz = dpsi(Rq, Zq)
        wR = w*Rq
        l2 = lambda e, r: np.sqrt(np.sum(wR*e)/np.sum(wR*r))
        e_s = np.sqrt(np.sum(wR*(s - sx)**2)/np.sum(wR))
        e_bp = l2((gs[..., 0] - gr/PSIB)**2 + (gs[..., 1] - gz/PSIB)**2, (gr**2 + gz**2)/PSIB**2)
        e_bt = l2((bphi - f0/Rq)**2, (f0/Rq)**2); e_bb = l2((bb - f0/Rq)**2, (f0/Rq)**2)
        # axis: Newton for u_h = 0 is overkill; take argmin of s_h on the quadrature points
        k_ = np.unravel_index(np.argmin(s), s.shape)
        print(f'f0 {f0:6.1f} p {p} proj {int(proj)} n {n:2d} dofs {nu:5d} steps {steps:3d} |g| {gn:.1e} '
              f's_L2 {e_s:.3e} Bpol {e_bp:.3e} Bphi {e_bt:.3e} Bphi_bar {e_bb:.3e} '
              f'min s {s.min():.1e} at R {Rq[k_]:.4f} Z {Zq[k_]:+.4f}  {time.time() - t0:.1f}s', flush=True)

if __name__ == '__main__':
    for a in sys.argv:
        if a.startswith('eps='): EPS = float(a[4:])
    f0s = [float(a) for a in sys.argv[1].split(',')]; p = int(sys.argv[2]); proj = int(sys.argv[3])
    ns = [int(a) for a in sys.argv[4:] if a.isdigit()]
    qpoly, res = qhat_fit()
    print(f'# q/f0 fit residual {res:.1e}; q/f0 axis {np.polyval(qpoly, -1):.4f} edge {np.polyval(qpoly, 2/1.1 - 1):.4f}')
    for n in ns: run(f0s, p, proj, n, qpoly, verbose='-v' in sys.argv, pin='-pin' in sys.argv, lm='-newton' not in sys.argv, guard='-guard' in sys.argv)
