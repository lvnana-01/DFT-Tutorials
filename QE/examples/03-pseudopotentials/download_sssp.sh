#!/usr/bin/env bash
# download_sssp.sh - download an SSSP 1.3.0 configuration from the Materials
# Cloud Archive, extract it and verify md5 checksums against the JSON file.
#
# This script is provided for documentation; run it yourself if you want it.
# UPF files must NOT be committed to this repository.
#
# Usage:
#   bash download_sssp.sh                       # PBE efficiency into ~/pseudo
#   FUNC=PBE PROT=precision bash download_sssp.sh
#   FUNC=PBEsol PROT=efficiency DEST=/scratch/pseudo bash download_sssp.sh
#   LINK="Mg O" bash download_sssp.sh           # also symlink Mg and O into ./pseudo
#
# Archive record (latest version as of 2026-11): rcyfm-68h65
# File names follow SSSP_<version>_<functional>_<protocol>.{tar.gz,json}

set -euo pipefail

VERSION=${VERSION:-1.3.0}
FUNC=${FUNC:-PBE}            # PBE or PBEsol
PROT=${PROT:-efficiency}     # efficiency or precision
DEST=${DEST:-$HOME/pseudo}
RECORD=${RECORD:-rcyfm-68h65}
LINK=${LINK:-}

NAME="SSSP_${VERSION}_${FUNC}_${PROT}"
URL="https://archive.materialscloud.org/api/records/${RECORD}/files"
OUTDIR="$DEST/$NAME"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

fetch() {  # fetch <remote-name> <local-name>
    if command -v wget >/dev/null 2>&1; then
        wget -q -O "$2" "$URL/$1/content"
    else
        curl -sSL -o "$2" "$URL/$1/content"
    fi
}

[ -f "$NAME.json" ]   || fetch "$NAME.json"   "$NAME.json"
[ -f "$NAME.tar.gz" ] || fetch "$NAME.tar.gz" "$NAME.tar.gz"
tar -xzf "$NAME.tar.gz"

echo "== verifying md5 checksums listed in $NAME.json =="
python3 - "$NAME.json" <<'PY'
import hashlib, json, os, sys
meta = json.load(open(sys.argv[1]))
bad = 0
for el, v in sorted(meta.items()):
    fn = v["filename"]
    if not os.path.exists(fn):
        print(f"MISSING {el:3s} {fn}")
        bad += 1
        continue
    md5 = hashlib.md5(open(fn, "rb").read()).hexdigest()
    if md5 != v["md5"]:
        print(f"BADMD5  {el:3s} {fn}")
        bad += 1
print(f"{len(meta)} elements checked, {bad} problems")
sys.exit(1 if bad else 0)
PY

echo "== recommended cutoffs (Ry) =="
python3 - "$NAME.json" <<'PY'
import json, sys
meta = json.load(open(sys.argv[1]))
for el in ["H", "C", "O", "Mg", "Al", "Si"]:
    if el in meta:
        v = meta[el]
        print(f"{el:3s} {v['filename']:42s} wfc={v['cutoff_wfc']:6.1f} rho={v['cutoff_rho']:7.1f}")
PY

if [ -n "$LINK" ]; then
    cd - >/dev/null
    mkdir -p pseudo
    for el in $LINK; do
        fn=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]]['filename'])" \
             "$OUTDIR/$NAME.json" "$el")
        ln -sf "$OUTDIR/$fn" "pseudo/$fn"
        echo "linked pseudo/$fn"
    done
fi

echo "Done: $OUTDIR"
