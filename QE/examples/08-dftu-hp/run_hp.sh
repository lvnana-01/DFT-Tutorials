#!/usr/bin/env bash
#
# run_hp.sh - linear-response (DFPT) Hubbard U for AFM NiO with hp.x,
# optionally iterated to self-consistency (U_in -> U_out -> U_in ...).
#
# One cycle = three steps (following QE HP/examples/example02, option 1):
#   1) pw.x  nio.lr_step1.in : "fake metal" scf (tiny smearing) at U_in
#   2) pw.x  nio.lr_step2.in : insulating scf, occupations='fixed',
#                              integer tot_magnetization, same nbnd
#   3) hp.x  nio.hp.in       : U_out written to nio_lr.Hubbard_parameters.dat
#
# Usage:
#   bash run_hp.sh               # one cycle starting from U_in = 0.0001 eV
#   NCYCLE=5 bash run_hp.sh      # up to 5 cycles, stop when |U_out-U_in| < TOL
#   PW_CMD="srun pw.x -nk 4" HP_CMD="srun hp.x -nk 4" bash run_hp.sh
#
# Environment variables:
#   PW_CMD  default "mpirun -np 4 pw.x"
#   HP_CMD  default "mpirun -np 4 hp.x"
#   U0      starting U in eV, default 0.0001 (i.e. start from plain PBE)
#   NCYCLE  maximum number of cycles, default 1
#   TOL     convergence threshold on U in eV, default 0.01
#
# NOTE: the tutorial repository owner does not run calculations; this
# script is provided for readers. Check every step's output by hand the
# first time you use it.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
HP_CMD=${HP_CMD:-mpirun -np 4 hp.x}
U_IN=${U0:-0.0001}
NCYCLE=${NCYCLE:-1}
TOL=${TOL:-0.01}

echo "# cycle   U_in(eV)   U_out(eV)" > u_history.dat

for cycle in $(seq 1 "$NCYCLE"); do
  rm -rf tmp_lr
  tag="c${cycle}"

  # --- step 1: fake-metal scf --------------------------------------------
  sed "s/@U_IN@/${U_IN}/g" nio.lr_step1.in > "nio.${tag}.step1.in"
  $PW_CMD -in "nio.${tag}.step1.in" > "nio.${tag}.step1.out"
  grep -q 'convergence has been achieved' "nio.${tag}.step1.out" || {
    echo "step 1 did not converge (cycle $cycle)" >&2; exit 1; }

  nbnd=$(grep 'number of Kohn-Sham states' "nio.${tag}.step1.out" | awk '{print $5}')
  totmag=$(grep 'total magnetization' "nio.${tag}.step1.out" | tail -1 | awk '{print $4}')
  # tot_magnetization must be an integer for fixed occupations
  totmag_int=$(awk -v m="$totmag" 'BEGIN{printf "%d", (m<0)? m-0.5 : m+0.5}')
  echo "cycle $cycle: nbnd=$nbnd total magnetization=$totmag -> $totmag_int"

  # --- step 2: insulating scf -------------------------------------------
  sed -e "s/@U_IN@/${U_IN}/g" -e "s/@NBND@/${nbnd}/" -e "s/@TOTMAG@/${totmag_int}/" \
      nio.lr_step2.in > "nio.${tag}.step2.in"
  $PW_CMD -in "nio.${tag}.step2.in" > "nio.${tag}.step2.out"
  grep -q 'convergence has been achieved' "nio.${tag}.step2.out" || {
    echo "step 2 did not converge (cycle $cycle)" >&2; exit 1; }

  # --- step 3: hp.x ------------------------------------------------------
  $HP_CMD -in nio.hp.in > "nio.${tag}.hp.out"
  cp nio_lr.Hubbard_parameters.dat "nio.${tag}.Hubbard_parameters.dat"

  # U of site 1 (label Ni1): last column of its line in the U table
  U_OUT=$(awk '/Hubbard U parameters/{f=1} f && $3=="Ni1"{print $NF; exit}' \
          nio_lr.Hubbard_parameters.dat)
  echo "$cycle  $U_IN  $U_OUT" | tee -a u_history.dat

  done_flag=$(awk -v a="$U_IN" -v b="$U_OUT" -v t="$TOL" \
              'BEGIN{d=a-b; if(d<0)d=-d; print (d<t)?1:0}')
  U_IN=$U_OUT
  if [ "$done_flag" = "1" ]; then
    echo "U converged: $U_OUT eV"
    break
  fi
done

echo "History in u_history.dat; last table:"
awk '/Hubbard U parameters/{f=1} f' nio_lr.Hubbard_parameters.dat | head -8
