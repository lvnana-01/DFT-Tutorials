#!/usr/bin/env bash
#
# run_uv.sh - DFT+U+V linear response for LiCoO2 (one cycle).
#   1) pw.x licoo2.uv.scf.in  (U, V = 0.0001 eV)
#   2) hp.x licoo2.hp.in      -> licoo2.Hubbard_parameters.dat, HUBBARD.dat
# HUBBARD.dat contains a ready-to-paste HUBBARD card (U on Co-3d as the
# Co-Co "V" entry, inter-site V Co-3d/O-2p with atom-pair indices).
#
# Usage:  bash run_uv.sh
#         PW_CMD="srun pw.x -nk 4" HP_CMD="srun hp.x -nk 4" bash run_uv.sh
#
# NOTE: the tutorial repository owner does not run calculations.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
HP_CMD=${HP_CMD:-mpirun -np 4 hp.x}

rm -rf tmp_uv
$PW_CMD -in licoo2.uv.scf.in > licoo2.uv.scf.out
grep -q 'convergence has been achieved' licoo2.uv.scf.out
$HP_CMD -in licoo2.hp.in > licoo2.hp.out

echo "== on-site U"
awk '/Hubbard U parameters/{f=1} f' licoo2.Hubbard_parameters.dat | head -8
echo "== HUBBARD card for the next pw.x run"
cat HUBBARD.dat
