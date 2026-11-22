"""Lowest eigenvalues of H v = mu M v of flux-map Hessian dumps (MILCA_FM3_HESSIAN)
with the two dominant toroidal harmonics of each mode. Usage: python hessian_toroidal_content.py FILE..."""
import pathlib, sys, numpy as np, scipy.linalg as sla
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "flux-map-spike"))
import hessian_spectrum as hs
for path in sys.argv[1:]:
    d = hs.read(path); n = d['h'].shape[0]; free = np.setdiff1d(np.arange(n), d['fixed'])
    H = d['h'].toarray(); H = 0.5*(H+H.T)[np.ix_(free, free)]; M = d['m'].toarray()[np.ix_(free, free)]
    mu, V = sla.eigh(H, M, subset_by_index=[0, 11])
    n0, nt = d['n0'], d['nt']
    out = []
    for j in range(len(mu)):
        v = np.zeros(n); v[free] = V[:, j]
        x = v[:n0*nt].reshape(nt, n0); y = v[n0*nt:].reshape(nt, n0); inter = ~d['wall']
        c = np.fft.fft((x + 1j*y)[:, inter], axis=0); pw = (abs(c)**2).sum(1); k = np.fft.fftfreq(nt, 1/nt).astype(int)
        top = np.argsort(pw)[::-1][:2]
        out.append('%+.2e[n=%+d %.0f%%, %+d %.0f%%]' % (mu[j], k[top[0]], 100*pw[top[0]]/pw.sum(), k[top[1]], 100*pw[top[1]]/pw.sum()))
    print(path.split('/')[-1], 'negative:', int((mu < 0).sum()), '\n   ', '\n    '.join(out[:8]))
