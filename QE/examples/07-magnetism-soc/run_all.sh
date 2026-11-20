#!/usr/bin/env bash
#
# run_all.sh - Chapter 7 workflow (magnetism, noncollinear, SOC).
#
# Usage:
#   bash run_all.sh [fe|ni|nio|soc|mae|nc|all]     (default: all)
#   PW_CMD="mpirun -np 8 pw.x -nk 4" bash run_all.sh mae
#
# Environment variables:
#   PW_CMD   default "mpirun -np 4 pw.x"
#   PP_CMD   launcher for bands.x, default "mpirun -np 4" (use "" for serial)
#
# Pseudopotentials must be placed in ./pseudo (see README.md).
# NOTE: the tutorial repository owner does not run calculations; this
# script is provided for readers.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
PP_CMD=${PP_CMD-mpirun -np 4}
TASK=${1:-all}

mag_summary () {
  # print final total/absolute magnetization and per-site moments
  local out=$1
  grep '^!' "$out" | tail -1
  grep -e 'total magnetization' -e 'absolute magnetization' "$out" | tail -2
  grep -A3 'Magnetic moment per site' "$out" | tail -4 || true
}

run_fe () {
  $PW_CMD -in fe.scf.in    > fe.scf.out
  $PW_CMD -in fe.nm.scf.in > fe.nm.scf.out
  echo "== bcc Fe FM";  mag_summary fe.scf.out
  echo "== bcc Fe NM";  grep '^!' fe.nm.scf.out | tail -1
}

run_ni () {
  $PW_CMD -in ni.scf.in > ni.scf.out
  echo "== fcc Ni FM"; mag_summary ni.scf.out
}

run_nio () {
  $PW_CMD -in nio_afm.scf.in > nio_afm.scf.out
  echo "== NiO AFM-II (plain PBE)"; mag_summary nio_afm.scf.out
  grep -e 'the Fermi energy' -e 'highest occupied' nio_afm.scf.out | tail -1 || true
}

run_soc () {
  $PW_CMD -in gaas_soc.scf.in   > gaas_soc.scf.out
  $PW_CMD -in gaas_soc.bands.in > gaas_soc.bands.out
  $PP_CMD bands.x -in gaas_soc.bands_pp.in > gaas_soc.bands_pp.out
  echo "== GaAs SOC"
  grep -e 'Noncollinear calculation' -e 'Non magnetic calculation' gaas_soc.scf.out || true
  grep 'highest occupied' gaas_soc.scf.out | tail -1
}

run_mae () {
  $PW_CMD -in fe_mae.lsda.in > fe_mae.lsda.out
  for d in 001 100 110; do
    # each direction works on its own copy of the LSDA save directory,
    # so that one nscf run cannot modify the files read by the next one
    rm -rf "tmp_${d}"; mkdir -p "tmp_${d}"
    cp -r tmp/fe_mae.save "tmp_${d}/"
    $PW_CMD -in fe_mae.nscf_${d}.in > fe_mae.nscf_${d}.out
  done
  echo "== force-theorem band energies (eV)"
  for d in 001 100 110; do
    printf '%s  ' "$d"; grep 'eband, Ef (eV)' fe_mae.nscf_${d}.out | tail -1
  done
}

run_nc () {
  $PW_CMD -in fe_noncolin.scf.in > fe_noncolin.scf.out
  echo "== bcc Fe noncollinear (no SOC)"
  grep '^!' fe_noncolin.scf.out | tail -1
  grep -A6 'Magnetic moment per site' fe_noncolin.scf.out | tail -7 || \
    grep -B1 -A4 'magnetization :' fe_noncolin.scf.out | tail -6
}

case "$TASK" in
  fe)  run_fe ;;
  ni)  run_ni ;;
  nio) run_nio ;;
  soc) run_soc ;;
  mae) run_mae ;;
  nc)  run_nc ;;
  all) run_fe; run_ni; run_nio; run_soc; run_mae; run_nc ;;
  *)   echo "unknown task: $TASK" >&2; exit 1 ;;
esac
