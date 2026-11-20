#!/usr/bin/env bash
# run_all.sh - run the three chapter-2 examples and print the key output lines.
# Nothing is executed unless RUN=1 is set; by default it only checks inputs.
#
# Usage:
#   bash run_all.sh                       # check that inputs and pseudopotentials exist
#   RUN=1 NP=4 bash run_all.sh            # run pw.x for si, si (ibrav=0), al, h2o
#   RUN=1 PW="srun -n 4 pw.x" bash run_all.sh   # inside a Slurm allocation

set -u
NP=${NP:-1}
RUN=${RUN:-0}
PW=${PW:-"mpirun -np $NP pw.x"}
INPUTS="si.scf.in si.scf.ibrav0.in al.scf.in h2o.scf.in"
PSEUDOS="Si.pbe-n-rrkjus_psl.1.0.0.UPF Al.pbe-n-kjpaw_psl.1.0.0.UPF O.pbe-n-kjpaw_psl.0.1.UPF H.pbe-rrkjus_psl.1.0.0.UPF"

missing=0
for p in $PSEUDOS; do
    if [ -f "pseudo/$p" ]; then echo "[ OK ] pseudo/$p"; else echo "[MISS] pseudo/$p"; missing=1; fi
done
for f in $INPUTS; do
    if grep -q $'\r' "$f"; then echo "[FAIL] $f has CRLF line endings"; missing=1; fi
done

if [ "$RUN" != "1" ]; then
    echo "Static checks done. Set RUN=1 to run pw.x."
    exit $missing
fi
[ $missing -ne 0 ] && { echo "Fix missing files first."; exit 1; }

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
for f in $INPUTS; do
    out=${f%.in}.out
    echo "=== $f -> $out"
    $PW -i "$f" > "$out" 2>&1
    grep -E "^!" "$out" || echo "  no final energy line"
    grep -E "highest occupied|Fermi energy" "$out"
    grep -E "convergence has been achieved|convergence NOT achieved" "$out" | tail -n 1
    grep -A1 "total   stress" "$out" | head -n 1
    grep -q "JOB DONE" "$out" && echo "  JOB DONE" || echo "  JOB DONE missing"
done
