#!/usr/bin/env bash
#
# qe_conv_test.sh - template-based convergence scan (ecutwfc / k-grid / degauss) for pw.x
#
# Usage:
#   bash qe_conv_test.sh -t TEMPLATE.in -m ecut  -v "30 40 50 60" [-d 8] [options]
#   bash qe_conv_test.sh -t TEMPLATE.in -m kpt   -v "4 6 8 10"    [options]
#   bash qe_conv_test.sh -t TEMPLATE.in -m smear -v "0.005 0.01 0.02" -s "K=12 12 12" -s SMEAR=mv
#   bash qe_conv_test.sh -t TEMPLATE.in -m ecut  -v "30 40 50" --run   # also run pw.x
#   bash qe_conv_test.sh --collect conv_ecut                          # only summarise
#   bash qe_conv_test.sh -h
#
# Placeholders in TEMPLATE.in (replaced by sed):
#   @ECUT@   ecutwfc in Ry               (mode ecut)
#   @ECUTRHO@ ecutrho = dual * ecutwfc   (mode ecut, dual from -d, default 4)
#   @K@      k grid "N N N"  or @K1@ @K2@ @K3@ for "N N N2" with -z N3 fixed (mode kpt)
#   @DEGAUSS@ smearing width in Ry        (mode smear)
#   Any other placeholder @KEY@ is filled with -s KEY=VALUE (e.g. -s "K=8 8 8" -s SMEAR=mv);
#   an unreplaced @...@ left in a generated input is an error.
#   Example lines:   ecutwfc = @ECUT@ , ecutrho = @ECUTRHO@
#                    K_POINTS automatic
#                      @K@ 0 0 0
#
# Options:
#   -t FILE    template pw.x input (required unless --collect)
#   -m MODE    ecut | kpt | smear
#   -v LIST    space-separated values (quoted)
#   -d DUAL    ecutrho/ecutwfc ratio for @ECUTRHO@ (default 4; 8-12 for US/PAW)
#   -z N3      fixed third k value for slabs / 2D (uses @K1@ @K2@ @K3@ -> N N N3)
#   -o DIR     work directory (default conv_<mode>)
#   -s KEY=VAL fixed value for placeholder @KEY@ (repeatable; applied after the scanned one)
#   --run      run pw.x in each sub-directory (default: only generate inputs).
#              Command from $PW_CMD, default "pw.x"; e.g. PW_CMD="mpirun -np 4 pw.x -nk 2"
#   --collect DIR  only collect results of an existing DIR
#
# Output: DIR/<value>/pw.in (+ pw.out after running) and DIR/summary.dat with
#   value  E_total(Ry)  dE_vs_last(meV/atom)  P(kbar)  n_scf  status
#   (n_scf = total number of SCF iterations in the run)
# The "!" line ("!    total energy = ... Ry") of pw.x is used; dE is relative to the
# LAST (most converged) value, divided by nat (read from "number of atoms/cell").
# Typical criterion: |dE| < 1 meV/atom (energies), stress < 0.5 kbar change.
#
# Dependencies: bash, sed, awk, grep; pw.x only with --run.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# qe_conv_test.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

set -euo pipefail

TEMPLATE=""; MODE=""; VALUES=""; DUAL=4; KZ=""; WORKDIR=""; RUN=0; COLLECT=""
SETS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) usage 0 ;;
        -t) TEMPLATE=$2; shift 2 ;;
        -m) MODE=$2; shift 2 ;;
        -v) VALUES=$2; shift 2 ;;
        -d) DUAL=$2; shift 2 ;;
        -z) KZ=$2; shift 2 ;;
        -o) WORKDIR=$2; shift 2 ;;
        -s) SETS+=("$2"); shift 2 ;;
        --run) RUN=1; shift ;;
        --collect) COLLECT=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; usage 1 ;;
    esac
done

PW_CMD=${PW_CMD:-pw.x}

collect() {
    local dir=$1 sum="$1/summary.dat"
    [ -d "$dir" ] || { echo "ERROR: $dir not found" >&2; exit 1; }
    local tmp; tmp=$(mktemp)
    for sub in "$dir"/*/; do
        sub=${sub%/}; local v; v=$(basename "$sub")
        local out="$sub/pw.out"
        if [ ! -f "$out" ]; then
            echo "$v NaN NaN 0 not_run" >> "$tmp"; continue
        fi
        local e p n st nat
        e=$(grep '^!' "$out" | tail -n 1 | awk '{print $(NF-1)}' || true)
        p=$(grep 'P=' "$out" | tail -n 1 | awk '{print $NF}' || true)
        n=$(grep -c 'iteration #' "$out" || true)
        nat=$(grep 'number of atoms/cell' "$out" | head -n 1 | awk '{print $NF}' || true)
        if grep -q 'JOB DONE' "$out"; then st=done; else st=unfinished; fi
        grep -q 'convergence NOT achieved' "$out" && st=scf_not_converged
        echo "$v ${e:-NaN} ${p:-NaN} ${n:-0} $st ${nat:-1}" >> "$tmp"
    done
    sort -g "$tmp" | awk '
        { v[NR]=$1; e[NR]=$2; p[NR]=$3; n[NR]=$4; s[NR]=$5; nat[NR]=($6==""?1:$6) }
        END {
            last=e[NR]
            printf "# %-10s %18s %22s %10s %6s %s\n", "value", "E_total(Ry)", "dE_vs_last(meV/atom)", "P(kbar)", "n_scf", "status"
            for (i=1;i<=NR;i++) {
                if (e[i]=="NaN" || last=="NaN") d="NaN"
                else d=sprintf("%.3f", (e[i]-last)*13605.693122994/nat[i])
                printf "  %-10s %18s %22s %10s %6s %s\n", v[i], e[i], d, p[i], n[i], s[i]
            }
        }' > "$sum"
    rm -f "$tmp"
    cat "$sum"
    echo "summary written to $sum"
}

if [ -n "$COLLECT" ]; then collect "$COLLECT"; exit 0; fi

[ -n "$TEMPLATE" ] && [ -f "$TEMPLATE" ] || { echo "ERROR: template (-t) missing" >&2; usage 1; }
case "$MODE" in ecut|kpt|smear) ;; *) echo "ERROR: -m must be ecut|kpt|smear" >&2; usage 1 ;; esac
[ -n "$VALUES" ] || { echo "ERROR: -v values missing" >&2; usage 1; }
if grep -q $'\r' "$TEMPLATE"; then
    echo "WARNING: $TEMPLATE has CRLF line endings; converting copies with tr -d '\\r'" >&2
fi
WORKDIR=${WORKDIR:-conv_$MODE}
mkdir -p "$WORKDIR"

for v in $VALUES; do
    d="$WORKDIR/$v"; mkdir -p "$d"
    case "$MODE" in
        ecut)
            rho=$(awk -v a="$v" -v b="$DUAL" 'BEGIN{printf "%g", a*b}')
            sed -e "s/@ECUT@/$v/g" -e "s/@ECUTRHO@/$rho/g" "$TEMPLATE" ;;
        kpt)
            k3=${KZ:-$v}
            sed -e "s/@K@/$v $v $v/g" -e "s/@K1@/$v/g" -e "s/@K2@/$v/g" -e "s/@K3@/$k3/g" "$TEMPLATE" ;;
        smear)
            sed -e "s/@DEGAUSS@/$v/g" "$TEMPLATE" ;;
    esac | tr -d '\r' > "$d/pw.in"
    for kv in "${SETS[@]+"${SETS[@]}"}"; do
        key=${kv%%=*}; val=${kv#*=}
        sed -i.bak "s/@${key}@/${val}/g" "$d/pw.in" && rm -f "$d/pw.in.bak"
    done
    if grep -q '@[A-Z0-9]*@' "$d/pw.in"; then
        echo "ERROR: unreplaced placeholder in $d/pw.in:" >&2
        grep -n '@[A-Z0-9]*@' "$d/pw.in" >&2; exit 1
    fi
    echo "generated $d/pw.in"
done

if [ "$RUN" -eq 1 ]; then
    for v in $VALUES; do
        d="$WORKDIR/$v"
        echo "running: (cd $d && $PW_CMD -in pw.in > pw.out)"
        # shellcheck disable=SC2086
        (cd "$d" && $PW_CMD -in pw.in > pw.out 2> pw.err) || echo "WARNING: pw.x failed in $d" >&2
    done
    collect "$WORKDIR"
else
    echo "inputs generated only. Run them yourself (or add --run), then:"
    echo "  bash $0 --collect $WORKDIR"
fi
