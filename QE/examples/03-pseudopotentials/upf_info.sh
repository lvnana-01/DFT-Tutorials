#!/usr/bin/env bash
# upf_info.sh - print the key header fields of UPF files and flag common problems.
# Read-only: it never modifies or runs anything.
#
# Usage:
#   bash upf_info.sh pseudo/*.UPF pseudo/*.upf
#   bash upf_info.sh ~/pseudo/SSSP_1.3.0_PBE_efficiency/*
#
# Handles both header styles:
#   UPF v2 : <PP_HEADER pseudo_type="PAW" z_valence="6.0" functional="PBE" .../>
#   UPF v1 : plain-text lines inside <PP_HEADER> ... </PP_HEADER>, e.g.
#            "    4.00000000000      Z valence"   (used by GBRV files in SSSP)

set -u
[ $# -eq 0 ] && { echo "usage: $0 file.UPF [file2.UPF ...]"; exit 1; }

get() {  # get <attribute> <file> -> value of attribute="..." (UPF v2)
    grep -o -m1 -E "$1=\"[^\"]*\"" "$2" | head -n 1 \
        | sed -E "s/^$1=\"[[:space:]]*//; s/[[:space:]]*\"$//"
}

get_v1() {  # get_v1 <label regex> <file> -> text before the label (UPF v1)
    sed -n '/<PP_HEADER>/,/<\/PP_HEADER>/p' "$2" | grep -m1 -E "$1" \
        | sed -E "s/($1).*//; s/^[[:space:]]+//; s/[[:space:]]+$//"
}

printf "%-42s %-6s %-8s %-6s %-6s %-20s %s\n" FILE TYPE REL NLCC ZVAL FUNCTIONAL SUGGESTED_WFC_RHO
funcs=""
for f in "$@"; do
    [ -f "$f" ] || continue
    if ! head -c 2000 "$f" | grep -q -E "<UPF|<PP_INFO|<PP_HEADER"; then
        echo "$(basename "$f"): does not look like a UPF file (HTML page or other format?)"
        continue
    fi
    if grep -q -m1 'pseudo_type="' "$f"; then
        typ=$(get pseudo_type "$f")
        rel=$(get relativistic "$f")
        nlcc=$(get core_correction "$f")
        zval=$(get z_valence "$f")
        fun=$(get functional "$f" | tr -s ' ')
    else
        typ=$(get_v1 "(Norm - Conserving|Ultrasoft|US|PAW).*(pseudopotential|PP)" "$f" | awk '{print $1}')
        rel="(v1)"
        nlcc=$(get_v1 "Nonlinear Core Correction" "$f")
        zval=$(get_v1 "Z valence" "$f")
        fun=$(get_v1 "Exchange-Correlation functional" "$f" | tr -s ' ')
    fi
    sug=$(grep -i -m2 "Suggested" "$f" | grep -o -E "[0-9]+\.?[0-9]* *Ry" | tr -d ' ' | paste -sd/ -)
    zfmt=$(printf '%.2f' "${zval:-0}" 2>/dev/null) || zfmt="${zval:-?}"
    printf "%-42s %-6s %-8s %-6s %-6s %-20s %s\n" "$(basename "$f")" "${typ:-?}" "${rel:-?}" \
        "${nlcc:-?}" "$zfmt" "${fun:-?}" "${sug:--}"
    funcs="$funcs|$fun"
    [ "$rel" = "full" ] && echo "   note: fully relativistic file (use with noncolin/lspinorb)"
done

# PBE may be written as "PBE", "SLA PW PBX PBC" (v2) or "SLA PW PBX PBC PBE" (v1);
# normalise these spellings before comparing. Other functionals are compared as-is.
nf=$(echo "$funcs" | tr '|' '\n' | sed '/^$/d' | tr '[:lower:]' '[:upper:]' \
     | sed -E 's/^SLA PW PBX PBC( PBE)?$/PBE/' | sort -u | wc -l)
if [ "$nf" -gt 1 ]; then
    echo "WARNING: more than one functional spelling found. If they really differ,"
    echo "         pw.x stops with 'inconsistent DFT read from PP files'"
    echo "         (input_dft can override this, at your own risk)."
fi
