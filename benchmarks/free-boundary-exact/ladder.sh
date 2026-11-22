#!/usr/bin/env bash
# Speed and accuracy ladders of VMEC2000, VMEC++ and MILCA on one exact case.
# Usage: ladder.sh CASE_DIR (from solovev_exact.py) XVMEC PYTHON MILCA_BUILD THREADS
set -u
dir=$(cd "$1" && pwd); xvmec=$2; py=$3; build=$(cd "$4" && pwd); threads=${5:-8}
here=$(cd "$(dirname "$0")" && pwd)
ref=$dir/reference.json
row() { printf '%-10s %-22s %8s  ' "$1" "$2" "$3"; "$py" "$here/compare.py" "$ref" "$4" "" | cut -c30-; }
if [ -z "${MILCA_ONLY:-}" ]; then
for spec in "8:16 32" "8:16 32 64" "12:16 32 64 128" "12:16 32 64 128 256" "16:16 32 64 128 256 512"; do
    mp=${spec%%:*}; ns=${spec#*:}; last=${ns##* }; nn=$(echo $ns | wc -w)
    w=$dir/vmec-m$mp-ns$last; mkdir -p "$w"; cp "$dir/mgrid_solovev_exact.nc" "$w/"
    ft=$(for i in $(seq $nn); do printf '1e-14 '; done); it=$(for i in $(seq $nn); do printf '50000 '; done)
    sed -e "s/MPOL = 8/MPOL = $mp/" -e "s/NS_ARRAY = .*/NS_ARRAY = $ns, FTOL_ARRAY = $ft, NITER_ARRAY = $it,/" \
        "$dir/input.solovev_exact" > "$w/input.solovev_exact"
    t0=$(python3 -c "import time;print(time.time())")
    (cd "$w" && "$xvmec" input.solovev_exact > x.log 2>&1)
    t1=$(python3 -c 'import time;print(time.time())')
    row VMEC2000 "mpol $mp ns $last" "$(python3 -c "print(f'{$t1-$t0:.2f}s')")" "$w/wout_solovev_exact.nc"
    (cd "$w" && env -u CPATH "$py" - > vpp.log 2>&1 <<PY
import time, vmecpp
inp = vmecpp.VmecInput.from_file("input.solovev_exact")
t = time.time(); out = vmecpp.run(inp, max_threads=$threads, verbose=False)
open("vpp.time", "w").write(f"{time.time() - t:.2f}s"); out.wout.save("wout_vmecpp.nc")
PY
    )
    row VMEC++ "mpol $mp ns $last" "$(cat "$w/vpp.time" 2>/dev/null)" "$w/wout_vmecpp.nc"
done
fi
"$build/milca_import" --format vmec "$dir/input.solovev_exact" "$dir/milca.nml" > /dev/null 2>&1
for nb in 12 16 24 32 48; do
    w=$dir/milca-nb$nb; mkdir -p "$w/output"; cp "$dir/mgrid_solovev_exact.nc" "$w/"
    sed "s/  angle_mpol = 4/  vacuum = 'bem', boundary_nodes = $nb, angle_mpol = 8/" "$dir/milca.nml" > "$w/in.nml"
    t0=$(python3 -c 'import time;print(time.time())')
    (cd "$w" && OMP_NUM_THREADS=$threads "$build/nested_equilibrium" --config in.nml > log 2>&1)
    t1=$(python3 -c 'import time;print(time.time())')
    row MILCA "$nb boundary nodes" "$(python3 -c "print(f'{$t1-$t0:.2f}s')")" "$w/output/nested_equilibrium.nc"
done
