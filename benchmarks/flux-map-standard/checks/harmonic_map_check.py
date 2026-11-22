"""Inverse-harmonic (Winslow) map of a section onto the unit disc by a boundary integral
(double layer as the real part of a Cauchy integral, trapezoidal rule, barycentric
interior evaluation), and the pullback of a disc O-grid. Own derivation; numpy only."""
import re, sys, numpy as np
def deck(path):
    t = open(path).read(); nfp = int(re.search(r'(?i)NFP\s*=\s*(\d+)', t).group(1)); c = {}
    for k, n, m, v in re.findall(r'(?i)(RBC|ZBS)\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*=\s*([-+0-9.eEdD]+)', t):
        c[(k.upper(), int(n), int(m))] = float(v.replace('d','e').replace('D','e'))
    return nfp, c
def wall(nfp, c, phi, th):
    R = np.zeros_like(th); Z = np.zeros_like(th); dR = R.copy(); dZ = R.copy(); ddR=R.copy(); ddZ=R.copy()
    for (k, n, m), v in c.items():
        a = m*th - n*nfp*phi
        if k == 'RBC': R += v*np.cos(a); dR -= m*v*np.sin(a); ddR -= m*m*v*np.cos(a)
        else: Z += v*np.sin(a); dZ += m*v*np.cos(a); ddZ -= m*m*v*np.sin(a)
    return R+1j*Z, dR+1j*dZ, ddR+1j*ddZ
class HMap:
    def __init__(s, nfp, c, phi, N=512):
        th = 2*np.pi*np.arange(N)/N; y, dy, ddy = wall(nfp, c, phi, th)
        area = np.sum(np.imag(np.conj(y)*dy))*np.pi/N
        s.sgn = 1.0
        if area < 0:   # make counter-clockwise
            th = -th; y, dy, ddy = wall(nfp, c, phi, th); dy = -dy; s.sgn = -1.0
        h = 2*np.pi/N; s.y, s.w, s.N = y, dy*h, N
        cen = np.sum(y*np.imag(np.conj(y)*dy))*h/3/abs(area)   # centroid: (2/3) int y dA-ish (polygon formula in the limit)
        s.cen = cen; s.area = abs(area)
        # arc length from the outboard point at centroid height (smooth in phi)
        sp = np.abs(dy); k = np.fft.fftfreq(N, 1/N); f = np.fft.fft(sp)
        L = f[0].real*h; s.L = L
        def arc(t):  # int_0^t |y'|
            t = np.atleast_1d(t); kk = k[1:]
            return (f[0].real*t + np.real((f[1:]/(1j*kk))[None,:]*(np.exp(1j*np.outer(t,kk))-1)).sum(1))/N
        # origin: root of Im(y - cen) = 0 on the outboard side
        cand = [i for i in range(N) if (y[i].imag-cen.imag)*(y[(i+1)%N].imag-cen.imag) <= 0 and y[i].real > cen.real]
        i0 = max(cand, key=lambda i: y[i].real); lo, hi = i0*h, (i0+1)*h
        for _ in range(60):
            mid = (lo+hi)/2; ym = wall(nfp, c, phi, s.sgn*np.array([mid]))[0][0]
            ylo = wall(nfp, c, phi, s.sgn*np.array([lo]))[0][0]
            if (ym.imag-cen.imag)*(ylo.imag-cen.imag) <= 0: hi = mid
            else: lo = mid
        s.t0 = (lo+hi)/2; tt = 2*np.pi*np.arange(N)/N
        s.ang = 2*np.pi*(arc(tt)-arc(s.t0))/L       # disc angle of the boundary nodes
        g = np.exp(1j*s.ang)
        D = y[None,:]-y[:,None]; np.fill_diagonal(D, 1)
        C = s.w[None,:]/D/(2j*np.pi); np.fill_diagonal(C, ddy/(2*dy)*h/(2j*np.pi))
        A = 0.5*np.eye(N) + C.real
        s.v = []
        for rhs in (g.real, g.imag):
            mu = np.linalg.solve(A, rhs)
            dmu = np.real(np.fft.ifft(1j*k*np.fft.fft(mu)))
            Cd = s.w[None,:]/D/(2j*np.pi); np.fill_diagonal(Cd, 0)
            v = mu + Cd@mu - mu*Cd.sum(1) + dmu*h/(2j*np.pi)
            s.v.append(v)
        s.bc_err = max(abs(s.v[0].real-g.real).max(), abs(s.v[1].real-g.imag).max())
        s.nfp, s.c, s.phi = nfp, c, phi
    def G(s, z):
        """xi+i eta and the 2x2 Jacobian rows (d xi/dR, d xi/dZ),(d eta/dR, d eta/dZ) at interior points z."""
        z = np.atleast_1d(z); q = s.w[None,:]/(s.y[None,:]-z[:,None]); den = q.sum(1); out = []; der = []
        for v in s.v:
            f = (q@v)/den; d = ((q/(s.y[None,:]-z[:,None]))*(v[None,:]-f[:,None])).sum(1)/den
            out.append(f.real); der.append(d)
        return out[0]+1j*out[1], der
    def bpoint(s, a):
        """wall point at disc angle a (inverse of the arc-length correspondence)."""
        a = np.atleast_1d(a) % (2*np.pi); t = s.t0 + a  # Newton on arc(t)
        tt = 2*np.pi*np.arange(s.N)/s.N; ext = np.concatenate([s.ang, s.ang[:1]+2*np.pi])
        base = (s.ang - s.ang[0]) ; base = np.concatenate([base, [2*np.pi]])
        t = np.interp((a - s.ang[0]) % (2*np.pi), base, np.concatenate([tt, [2*np.pi]]))
        return wall(s.nfp, s.c, s.phi, s.sgn*t)[0]
    def inv(s, zeta, z0):
        """z with G(z) = zeta, Newton from z0 (arrays)."""
        z = np.array(z0, complex)
        for it in range(50):
            g, d = s.G(z); r = g - zeta
            a, b, cc, dd = d[0].real, -d[0].imag, d[1].real, -d[1].imag
            det = a*dd - b*cc
            step = -((dd*r.real - b*r.imag) + 1j*(-cc*r.real + a*r.imag))/det
            z = z + step
            if abs(step).max() < 1e-13: break
        return z, det
    def table(s, nr=48, na=256):
        """z on the polar grid rho_j = 1 - j/nr (j = 0..nr-1), marching inward."""
        a = 2*np.pi*np.arange(na)/na; z = s.bpoint(a); rows = [z]; dets = []
        zc = s.cen
        z = z + 0.5/nr*(zc - z)   # first guess a little inside
        for j in range(1, nr):
            z, det = s.inv((1-j/nr)*np.exp(1j*a), z); rows.append(z); dets.append(det)
        return a, np.array(rows), np.array(dets)
if __name__ == '__main__':
    nfp, c = deck(sys.argv[1]); N = int(sys.argv[3]) if len(sys.argv) > 3 else 512
    for frac in [float(x) for x in sys.argv[2].split(',')]:
        phi = frac*2*np.pi/nfp
        h = HMap(nfp, c, phi, N); a, Z, det = h.table()
        # det of G along rings: > 0 everywhere <=> diffeomorphism; ratio = quality of the pullback
        zc, _ = h.inv(np.array([0j]), np.array([Z[-1].mean()]))
        print('phi/period %.3f  N %d  bc err %.1e  det G min %.3e max %.3e (ratio %.4f)  image of disc centre %.4f %+.4f  centroid %.4f %+.4f' % (frac, N, h.bc_err, det.min(), det.max(), det.min()/det.max(), zc[0].real, zc[0].imag, h.cen.real, h.cen.imag))
