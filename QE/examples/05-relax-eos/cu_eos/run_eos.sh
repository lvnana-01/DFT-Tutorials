#!/usr/bin/env bash
#
# run_eos.sh - fcc Cu E(a) scan for Chapter 5, then fit with ev.x.
#
# Usage:
#   bash run_eos.sh            # write inputs only (a_<value>/pw.in)
#   RUN=1 bash run_eos.sh      # also run pw.x and ev.x
#   bash run_eos.sh collect    # only collect energies and run ev.x
#
# Environment variables:
#   RUN=1      run pw.x (default 0)
#   PW_CMD     default "pw.x", e.g. "mpirun -np 8 pw.x -nk 4"
#   ALATS      list of lattice parameters in Angstrom
#              default "3.50 3.54 3.58 3.62 3.66 3.70 3.74"
#   EOS        ev.x equation of state: 1=birch1 2=birch2 3=keane 4=murnaghan
#              (default 2)
#
# Output: cu_eos.dat (two columns: a in Angstrom, E in Ry) and cu_eos.ev.out
# (ev.x fit). ev.x is interactive; its five answers are given on stdin below.
# NOTE: the tutorial repository owner does not run calculations.

set -euo pipefail
RUN=${RUN:-0}
PW_CMD=${PW_CMD:-pw.x}
ALATS=${ALATS:-3.50 3.54 3.58 3.62 3.66 3.70 3.74}
EOS=${EOS:-2}

if [ "${1:-}" != "collect" ]; then
    for a in $ALATS; do
        d=a_$a; mkdir -p "$d"
        [ -e "$d/pseudo" ] || ln -s ../pseudo "$d/pseudo"
        sed "s/@ALAT@/$a/" cu_eos_template.in > "$d/pw.in"
        if [ "$RUN" = "1" ]; then
            (cd "$d" && $PW_CMD -in pw.in > pw.out)
        fi
    done
fi

: > cu_eos.dat
for a in $ALATS; do
    f=a_$a/pw.out
    [ -f "$f" ] || { echo "missing $f"; continue; }
    # '|| true': an unfinished run must not abort the script (set -e/pipefail)
    e=$(grep '^!' "$f" | tail -1 | awk '{print $5}' || true)
    p=$(grep 'P=' "$f" | tail -1 | awk '{print $NF}' || true)
    [ -n "$e" ] || { echo "no final energy in $f (not converged?)"; continue; }
    echo "$a $e" >> cu_eos.dat
    echo "a = $a A  E = $e Ry  P = $p kbar"
done

if [ -s cu_eos.dat ] && command -v ev.x > /dev/null; then
    # answers: units, lattice type, EOS type, input file, output file
    ev.x << EOI
Ang
fcc
$EOS
cu_eos.dat
cu_eos.ev.out
EOI
    head -3 cu_eos.ev.out
fi
