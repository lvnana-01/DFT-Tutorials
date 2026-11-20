#!/usr/bin/env bash
#
# run_all.sh - Chapter 6 workflow for Si:
#   scf -> bands -> bands.x -> nscf -> dos.x -> projwfc.x -> sumpdos.x
#
# Usage:
#   bash run_all.sh                 # run everything with $PW_CMD / $PP_CMD
#   PW_CMD="mpirun -np 8 pw.x -nk 4" PP_CMD="mpirun -np 8" bash run_all.sh
#
# Environment variables:
#   PW_CMD   default "mpirun -np 4 pw.x"
#   PP_CMD   launcher for post-processing codes, default "mpirun -np 4"
#            (use PP_CMD="" to run bands.x/dos.x/projwfc.x serially)
#
# Pseudopotential Si.pbe-n-rrkjus_psl.1.0.0.UPF must be in ./pseudo.
# NOTE: the tutorial repository owner does not run calculations.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
PP_CMD=${PP_CMD-mpirun -np 4}

$PW_CMD -in si.scf.in   > si.scf.out
grep -e '^!' -e 'highest occupied' si.scf.out

$PW_CMD -in si.bands.in > si.bands.out
$PP_CMD bands.x -in si.bands_pp.in > si.bands_pp.out
grep 'high-symmetry point' si.bands_pp.out

$PW_CMD -in si.nscf.in  > si.nscf.out
grep -e 'the Fermi energy is' -e 'highest occupied' si.nscf.out || true

$PP_CMD dos.x -in si.dos.in > si.dos.out
head -1 si.dos

$PP_CMD projwfc.x -in si.projwfc.in > si.projwfc.out
grep -A4 'Lowdin Charges' si.projwfc.out | head -5 || true
grep 'Spilling Parameter' si.projwfc.out

# sum all s and all p projections of every Si atom (sumpdos.x writes to stdout)
sumpdos.x si.pdos_atm#*\(Si\)_wfc#*\(s\) > si.pdos_Si_s.dat
sumpdos.x si.pdos_atm#*\(Si\)_wfc#*\(p\) > si.pdos_Si_p.dat
ls si.pdos_*
