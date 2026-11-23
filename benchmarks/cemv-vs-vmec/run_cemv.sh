#!/bin/sh
# One CEMV (flux map) run of the head-to-head: the driver test_flux_map_qa on a
# prescribed-iota deck, from the coarse stage 2:2:9 (start u = G, no axis
# information) to the stage DEGREE:CELLS:PLANES. Writes OUT/log.txt (stage
# lines with unknowns, Newton steps, stage times, the force balance per band),
# OUT/time.txt (/usr/bin/time -l of the whole process: peak memory) and the
# files benchmarks/label-qa/compare.py reads.
# Usage: run_cemv.sh BUILD DECK OUT DEGREE:CELLS:PLANES THREADS [start stage]
set -e
build=$1; deck=$2; out=$3; stage=$4; threads=$5; first=${6:-2:2:9}
mkdir -p "$out"
OMP_NUM_THREADS=$threads /usr/bin/time -l "$build/test_flux_map_qa" "$deck" "$out" \
    "$first,$stage" > "$out/log.txt" 2> "$out/time.txt"
