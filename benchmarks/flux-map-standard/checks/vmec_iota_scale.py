"""VMEC++ axes of a prescribed-iota deck with iota scaled by the given factors (ns 16/51/101).
Usage: python vmec_iota_scale.py DECK OUTDIR KAPPA..."""
import sys, re, pathlib, numpy as np, vmecpp, netCDF4
deck, out = sys.argv[1], pathlib.Path(sys.argv[2])
for kap in [float(x) for x in sys.argv[3:]]:
    vin = vmecpp.VmecInput.from_file(deck)
    vin.ai = np.array(vin.ai)*kap
    vin.ns_array = np.array([16,51,101]); vin.ftol_array = np.full(3,1e-14); vin.niter_array = np.full(3,40000)
    r = vmecpp.run(vin, max_threads=4, verbose=False); p = out/f'wout_k{kap:.2f}.nc'; r.wout.save(str(p))
    f = netCDF4.Dataset(p); xn=f['xn'][:]; nfp=int(f['nfp'][:]); phi = 2*np.pi*np.arange(64)/(64*nfp); a=-np.outer(phi,xn)
    d = out/f'ref_k{kap:.2f}'; d.mkdir(exist_ok=True)
    np.savetxt(d/'axis.txt', np.c_[phi, np.cos(a)@np.array(f['rmnc'][0]), np.sin(a)@np.array(f['zmns'][0])])
    print(kap, 'iota0', float(f['iotaf'][0]), 'wb', float(f['wb'][:])*4*np.pi**2, flush=True)
