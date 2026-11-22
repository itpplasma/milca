"""Third-code check of the vacuum QA axis: DESC (zero current, zero pressure) from the VMEC
deck, resolution ladder, axis against a reference axis.txt.
Usage: python desc_axis.py DECK REFERENCE/axis.txt OUTPREFIX L:M:N,L:M:N,...  (python with desc-opt)"""
import sys, time, numpy as np
import desc
from desc.geometry import FourierRZToroidalSurface
from desc.equilibrium import Equilibrium
from desc.profiles import PowerSeriesProfile
from desc.continuation import solve_continuation_automatic
from desc.grid import LinearGrid
deck, refp, out = sys.argv[1], sys.argv[2], sys.argv[3]
res = [tuple(int(v) for v in x.split(':')) for x in sys.argv[4].split(',')]
ref = np.loadtxt(refp); phi = ref[:,0]
surf = FourierRZToroidalSurface.from_input_file(deck)
def report(eq, tag):
    ga = LinearGrid(rho=np.array([0.0]), theta=np.array([0.0]), zeta=phi, NFP=1)
    da = eq.compute(['R','Z','iota'], grid=ga); a = np.c_[phi, np.array(da['R']), np.array(da['Z'])]
    np.savetxt(f'{out}_{tag}.txt', a)
    f = eq.compute(['<|F|>_vol','<|grad(|B|^2)|/2mu0>_vol','W_B'])
    best = None
    for nm, b in (('as is', a), ('phi flipped', np.roll(a[::-1], 1, axis=0))):
        dR = b[:,1]-ref[:,1]; dZ = b[:,2]-ref[:,2]; m = np.hypot(dR,dZ).max()
        if best is None or m < best[0]: best = (m, dR.mean(), nm)
    print('%s: DESC axis - VMEC ref: mean dR %+.4f mm  max|d| %.4f mm (%s); iota0 %.6f  normalised force %.2e' % (tag, 1e3*best[1], 1e3*best[0], best[2], float(np.array(da['iota'])[0]), float(f['<|F|>_vol'])/float(f['<|grad(|B|^2)|/2mu0>_vol'])), flush=True)
eq = None
for (L, M, N) in res:
    t0 = time.time()
    if eq is None:
        eq = Equilibrium(surface=surf, L=L, M=M, N=N, NFP=surf.NFP, Psi=0.08385727554, current=PowerSeriesProfile([0.0]), pressure=PowerSeriesProfile([0.0]), sym=True)
        eq = solve_continuation_automatic(eq, verbose=0, ftol=1e-12, xtol=1e-12, gtol=1e-12, maxiter=300)[-1]
    else:
        eq.change_resolution(L=L, M=M, N=N, L_grid=2*L, M_grid=2*M, N_grid=2*N)
        eq.solve(verbose=0, ftol=1e-12, xtol=1e-12, gtol=1e-12, maxiter=300)
    print('L M N', L, M, N, 'time %.0f s' % (time.time()-t0), flush=True)
    report(eq, f'{L}_{M}_{N}')
