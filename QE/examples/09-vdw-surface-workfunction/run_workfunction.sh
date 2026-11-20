#!/usr/bin/env bash
#
# run_workfunction.sh - Al(111) slab: relax -> scf -> pp.x (plot_num=11)
#                       -> average.x (planar + macroscopic average)
#                       -> work function = V_vacuum - E_Fermi
#
# Usage:
#   bash run_workfunction.sh            # full chain
#   SKIP_RELAX=1 bash run_workfunction.sh
#   PW_CMD="srun pw.x -nk 4" PP_CMD="srun -n 4" bash run_workfunction.sh
#
# average.x must run on ONE process (it stops otherwise).
# avg.dat columns: z (bohr), planar average (Ry), macroscopic average (Ry).
# The tutorial repository owner does not run calculations.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
PP_CMD=${PP_CMD-mpirun -np 4}

if [ -z "${SKIP_RELAX:-}" ]; then
  $PW_CMD -in al111_slab.relax.in > al111_slab.relax.out
  grep -e 'bfgs converged' -e 'End final coordinates' al111_slab.relax.out || true
  echo "Paste the relaxed ATOMIC_POSITIONS into al111_slab.scf.in if you want"
  echo "the work function of the relaxed slab (this script does not edit it)."
fi

$PW_CMD -in al111_slab.scf.in > al111_slab.scf.out
grep -e '^!' -e 'the Fermi energy' al111_slab.scf.out | tail -2
grep -A6 'Computed dipole along edir' al111_slab.scf.out | tail -7 || true

$PP_CMD pp.x -in pp_v.in > pp_v.out
average.x < average.in > average.out            # serial on purpose
ls -l avg.dat

EF=$(grep 'the Fermi energy' al111_slab.scf.out | tail -1 | awk '{print $5}')
echo "E_Fermi (eV) = $EF"
# Rough estimate of the two vacuum levels: planar average (Ry) at the
# middle of each vacuum region. For this cell (c = 30 A) the slab spans
# fractional z = 0.07-0.46 and the dipole sawtooth drops at 0.75-0.80, so
# the vacuum middles are at fractional z ~ 0.60 (top) and ~ 0.93 (bottom).
# For production use: python ../../scripts/qe_workfunction.py -h
awk -v ef="$EF" '
  { z[NR]=$1; v[NR]=$2 }
  END {
    L = z[NR] + (z[2]-z[1]);              # cell length in bohr
    for (k=1; k<=2; k++) {
      f = (k==1) ? 0.60 : 0.93;
      best=1; dmin=1e9;
      for (i=1; i<=NR; i++) { d=z[i]-f*L; if (d<0) d=-d; if (d<dmin){dmin=d; best=i} }
      vac = v[best]*13.605693;            # Ry -> eV
      printf "side %d (frac z=%.2f): V_vac = %9.4f eV   Phi = %8.4f eV\n", k, f, vac, vac-ef
    }
  }' avg.dat
