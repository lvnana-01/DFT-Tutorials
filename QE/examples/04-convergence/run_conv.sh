#!/usr/bin/env bash
#
# run_conv.sh - generate (and optionally run) convergence-test inputs for
#               Chapter 4 of the QE tutorial.
#
# Usage:
#   bash run_conv.sh ecut            # Si: scan ecutwfc at fixed dual and k grid
#   bash run_conv.sh kpt             # Si: scan k grid at fixed cutoff
#   bash run_conv.sh alsmear         # Al: scan k grid x degauss (2D table)
#   bash run_conv.sh collect DIR     # only summarise an existing DIR
#
# Environment variables (all optional):
#   RUN=1          actually call pw.x (default 0: only write inputs)
#   PW_CMD         pw.x command, default "pw.x"  (e.g. "mpirun -np 4 pw.x -nk 2")
#   ECUTS          list for ecut mode, default "25 30 35 40 45 50 60 70 80"
#   DUAL           ecutrho/ecutwfc ratio, default 8 (ultrasoft / PAW)
#   KGRID          fixed k grid for ecut mode, default 8
#   ECUT_FIX       fixed ecutwfc for kpt mode, default 40
#   KS             list of N for N x N x N grids, default "2 4 6 8 10 12 14 16"
#                  (kpt mode) or "6 8 10 12 14 16 20 24" (alsmear mode)
#   DEGAUSS_LIST   list for alsmear mode, default "0.005 0.01 0.02 0.04"
#   SMEAR          smearing type for alsmear mode, default "mv"
#
# Each calculation lives in its own sub-directory and writes pw.out there.
# The summary file summary.dat contains:
#   label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf_iter  status
# (dE is relative to the LAST line of the scan = the most converged one.)
#
# NOTE: the tutorial repository owner does not run calculations; this script
# is provided for readers. Pseudopotentials must be placed in ./pseudo first.

set -euo pipefail

MODE=${1:-}
RUN=${RUN:-0}
PW_CMD=${PW_CMD:-pw.x}
RY2MEV=13605.693

usage() { sed -n '3,27p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

# collect DIR : one line per sub-directory (in creation order = list order)
collect() {
    local dir=$1 nat=$2
    local out=$dir/summary.dat
    : > "$out.tmp"
    while read -r label; do
        local f=$dir/$label/pw.out e p n st
        if [ -f "$f" ]; then
            # '|| true': an unfinished run must not abort the summary (set -e/pipefail)
            e=$(grep '^!' "$f" | tail -1 | awk '{print $5}' || true)
            p=$(grep 'P=' "$f" | tail -1 | awk '{print $NF}' || true)
            n=$(grep 'convergence has been achieved in' "$f" | tail -1 | awk '{print $6}' || true)
            if grep -q 'JOB DONE' "$f"; then st=done; else st=unfinished; fi
        else
            e=""; p=""; n=""; st=not_run
        fi
        echo "$label ${e:-NA} ${p:-NA} ${n:-NA} $st" >> "$out.tmp"
    done < "$dir/labels.txt"
    local eref
    eref=$(awk '$2!="NA"{e=$2} END{print e}' "$out.tmp")
    {
        echo "# label  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status"
        awk -v r="${eref:-NA}" -v c=$RY2MEV -v nat="$nat" '{
            if ($2=="NA" || r=="NA") d="NA"; else d=sprintf("%.3f", ($2-r)*c/nat);
            printf "%-14s %18s %12s %10s %6s %s\n", $1, $2, d, $3, $4, $5 }' "$out.tmp"
    } > "$out"
    rm -f "$out.tmp"
    cat "$out"
}

run_one() {
    local d=$1
    if [ "$RUN" = "1" ]; then
        (cd "$d" && $PW_CMD -in pw.in > pw.out)
    fi
}

link_pseudo() {
    local d=$1
    mkdir -p "$d"
    [ -e "$d/pseudo" ] || ln -s ../../pseudo "$d/pseudo"
}

case "$MODE" in
    -h|--help|"") usage 0 ;;
    ecut)
        DIR=conv_ecut; mkdir -p $DIR; : > $DIR/labels.txt
        DUAL=${DUAL:-8}; KGRID=${KGRID:-8}
        for e in ${ECUTS:-25 30 35 40 45 50 60 70 80}; do
            d=$DIR/ecut_$e; link_pseudo "$d"
            rho=$(awk -v a="$e" -v b="$DUAL" 'BEGIN{print a*b}')
            sed -e "s/@ECUT@/$e/" -e "s/@ECUTRHO@/$rho/" \
                -e "s/@K@/$KGRID $KGRID $KGRID/" si_ecut_template.in > "$d/pw.in"
            echo "ecut_$e" >> $DIR/labels.txt
            run_one "$d"
        done
        collect $DIR 2 ;;
    kpt)
        DIR=conv_kpt; mkdir -p $DIR; : > $DIR/labels.txt
        DUAL=${DUAL:-8}; E=${ECUT_FIX:-40}
        rho=$(awk -v a="$E" -v b="$DUAL" 'BEGIN{print a*b}')
        for k in ${KS:-2 4 6 8 10 12 14 16}; do
            d=$DIR/k_$k; link_pseudo "$d"
            sed -e "s/@ECUT@/$E/" -e "s/@ECUTRHO@/$rho/" \
                -e "s/@K@/$k $k $k/" si_ecut_template.in > "$d/pw.in"
            echo "k_$k" >> $DIR/labels.txt
            run_one "$d"
        done
        collect $DIR 2 ;;
    alsmear)
        SMEAR=${SMEAR:-mv}
        for s in ${DEGAUSS_LIST:-0.005 0.01 0.02 0.04}; do
            DIR=conv_al_${SMEAR}_$s; mkdir -p $DIR; : > $DIR/labels.txt
            for k in ${KS:-6 8 10 12 14 16 20 24}; do
                d=$DIR/k_$k; link_pseudo "$d"
                sed -e "s/@K@/$k $k $k/" -e "s/@DEGAUSS@/$s/" \
                    -e "s/@SMEAR@/$SMEAR/" al_kpt_template.in > "$d/pw.in"
                echo "k_$k" >> $DIR/labels.txt
                run_one "$d"
            done
            echo "== smearing=$SMEAR degauss=$s Ry =="
            collect $DIR 1
        done ;;
    collect)
        [ $# -ge 2 ] || usage 1
        nat=2; case "$2" in conv_al_*) nat=1 ;; esac
        collect "$2" "$nat" ;;
    *) usage 1 ;;
esac
