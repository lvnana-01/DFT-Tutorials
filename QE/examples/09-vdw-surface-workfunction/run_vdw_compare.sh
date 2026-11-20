#!/usr/bin/env bash
#
# run_vdw_compare.sh - relax bulk graphite with several dispersion schemes
# and tabulate the interlayer distance c/2 and the final energy.
#
# Each method replaces the line marked "@VDW@" in graphite.vdw.in:
#   none      : plain PBE (no dispersion; graphite is barely bound)
#   d2        : vdw_corr = 'grimme-d2'
#   d3bj      : vdw_corr = 'grimme-d3' with dftd3_version = 4 (BJ damping)
#   d3zero    : vdw_corr = 'grimme-d3' with dftd3_version = 3 (zero damping)
#   ts        : vdw_corr = 'ts'  -- NOT in the default list: pw.x stops with
#               'Tkatchenko-Scheffler not implemented with USPP/PAW', and the
#               SSSP carbon pseudopotential used here is PAW. Use a
#               norm-conserving C pseudopotential before running 'ts'.
#   rvv10     : input_dft = 'rvv10'          (nonlocal functional)
#   df2b86r   : input_dft = 'vdw-df2-b86r'   (rev-vdW-DF2, nonlocal)
#
# Usage:
#   bash run_vdw_compare.sh                   # all methods
#   bash run_vdw_compare.sh d3bj rvv10        # selected methods
#   PW_CMD="srun pw.x -nk 4" bash run_vdw_compare.sh
#
# NOTE: energies of different functionals are NOT comparable with each
# other; compare geometries (c, interlayer distance) instead.
# The tutorial repository owner does not run calculations.

set -euo pipefail
PW_CMD=${PW_CMD:-mpirun -np 4 pw.x}
METHODS=${*:-none d2 d3bj d3zero rvv10 df2b86r}

line_for () {
  case "$1" in
    none)    echo "  vdw_corr    = 'none'" ;;
    d2)      echo "  vdw_corr    = 'grimme-d2'" ;;
    d3bj)    echo "  vdw_corr    = 'grimme-d3'" ;;
    d3zero)  echo "  vdw_corr    = 'grimme-d3'" ;;
    ts)      echo "  vdw_corr    = 'ts'" ;;
    rvv10)   echo "  input_dft   = 'rvv10'" ;;
    df2b86r) echo "  input_dft   = 'vdw-df2-b86r'" ;;
    *) echo "unknown method $1" >&2; exit 1 ;;
  esac
}

printf '%-9s %12s %18s\n' method c_final_A E_final_Ry > vdw_summary.dat
for m in $METHODS; do
  d="vdw_$m"; mkdir -p "$d"
  newline=$(line_for "$m")
  ver=4; [ "$m" = "d3zero" ] && ver=3
  sed -e "s|^  vdw_corr    = 'grimme-d3'      ! @VDW@|${newline}|" \
      -e "s|^  dftd3_version = 4 .*|  dftd3_version = ${ver}|" \
      -e "s|outdir        = './tmp'|outdir        = './tmp_${m}'|" \
      -e "s|pseudo_dir    = './pseudo'|pseudo_dir    = '../pseudo'|" \
      graphite.vdw.in > "$d/graphite.in"
  ( cd "$d" && $PW_CMD -in graphite.in > graphite.out )
  # final cell: third CELL_PARAMETERS row after "Begin final coordinates".
  # The header is "CELL_PARAMETERS (alat= <bohr>)", "(bohr)" or "(angstrom)";
  # convert the z component of the third vector to Angstrom.
  c=$(awk '/Begin final coordinates/{f=1}
           f && /CELL_PARAMETERS/ {
             s = 1.0
             if ($0 ~ /alat=/) { h = $0; sub(/.*alat=/, "", h); sub(/\).*/, "", h); s = h * 0.529177 }
             else if ($0 ~ /bohr/) { s = 0.529177 }
             getline; getline; getline; printf "%.4f", $3 * s; exit }' "$d/graphite.out")
  e=$(grep '^!' "$d/graphite.out" | tail -1 | awk '{print $5}')
  printf '%-9s %12s %18s\n' "$m" "${c:-NA}" "${e:-NA}" | tee -a vdw_summary.dat
done
echo "c_final_A = z component of the third final cell vector, converted to Angstrom."
