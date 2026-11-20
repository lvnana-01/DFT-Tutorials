#!/bin/bash
# Full DFPT workflow for Si: scf -> ph.x (4x4x4 q grid) -> q2r.x -> matdyn.x
# Usage: bash run_si.sh            (prints the commands only, dry run)
#        RUN=1 bash run_si.sh      (actually executes them)
# Override launchers, e.g. PW="srun pw.x -nk 4" PH="srun ph.x -nk 4".
set -euo pipefail
PW=${PW:-"mpirun -np 8 pw.x -nk 4"}
PH=${PH:-"mpirun -np 8 ph.x -nk 4"}
SERIAL=${SERIAL:-""}            # e.g. "mpirun -np 1" if needed on your system

run () {
    echo "+ $*"
    if [ "${RUN:-0}" = "1" ]; then eval "$@"; fi
}

run "$PW -in si.scf.in > si.scf.out"
run "$PH -in si.ph.in > si.ph.out"
run "grep -e 'Convergence has been achieved' -e 'No convergence' si.ph.out | sort | uniq -c"
run "$SERIAL q2r.x -in si.q2r.in > si.q2r.out"
run "$SERIAL matdyn.x -in si.matdyn.in > si.matdyn.out"
run "$SERIAL matdyn.x -in si.matdyn_dos.in > si.matdyn_dos.out"
run "python ../../scripts/qe_phonon_plot.py -h"
