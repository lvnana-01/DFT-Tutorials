#!/bin/bash
# Charge-density difference and Bader workflow for CO/Cu(111).
# Usage: bash run_diff.sh            (prints the commands only, dry run)
#        RUN=1 bash run_diff.sh      (actually executes them)
# Override the launcher, e.g. PW="srun pw.x -nk 4" PP="srun pp.x".
set -euo pipefail
PW=${PW:-"mpirun -np 8 pw.x -nk 4"}
PP=${PP:-"mpirun -np 8 pp.x"}
PJ=${PJ:-"mpirun -np 8 projwfc.x"}

run () {
    echo "+ $*"
    if [ "${RUN:-0}" = "1" ]; then eval "$@"; fi
}

for s in ab slab co; do
    run "$PW -in $s.scf.in > $s.scf.out"
done

# All three runs must use the same FFT grid
run "grep 'Dense  grid' ab.scf.out slab.scf.out co.scf.out"

for s in ab slab co; do
    run "$PP -in rho_$s.pp.in > rho_$s.pp.out"
done
run "$PP -in diff_plot.pp.in > diff_plot.pp.out"

# Bader inputs (PAW): valence density + all-electron reference
run "$PP -in ab_rho.pp.in > ab_rho.pp.out"
run "$PP -in ab_ae.pp.in > ab_ae.pp.out"
run "bader ab_rho.cube -ref ab_ae.cube"

# Optional: STM image and Lowdin charges
run "$PP -in stm.pp.in > stm.pp.out"
run "$PJ -in ab.projwfc.in > ab.projwfc.out"
