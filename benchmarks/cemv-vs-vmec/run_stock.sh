#!/bin/sh
# Run the prescribed-transform deck verbatim and reject unsuccessful VMEC output.
set -eu
if [ "$#" -lt 3 ]; then
    echo "usage: $0 XVMECPATH DECK OUT [LAUNCHER ARGS ...]" >&2
    exit 2
fi
vmec=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
deck=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"
cd "$3"
shift 3
case $(uname -s) in
    Darwin) time_option=-l ;;
    *) time_option=-v ;;
esac
# Run in a fresh directory so a stale Wout cannot qualify a failed invocation.
if [ -e wout_benchmark.nc ]; then
    echo "Wout already exists in $PWD" >&2
    exit 2
fi
cp "$deck" input.benchmark
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/time "$time_option" \
    "$@" "$vmec" benchmark > run.log 2> time.txt
# Stock VMEC2000 can write ier_flag=0 after exhausting the iteration budget.
if grep -qi 'Try increasing NITER' run.log ||
    ! grep -q 'EXECUTION TERMINATED NORMALLY' run.log ||
    ! ncdump -v ier_flag wout_benchmark.nc | grep -q 'ier_flag = 0 *;'; then
    echo "VMEC did not converge. See $PWD/run.log" >&2
    exit 1
fi
