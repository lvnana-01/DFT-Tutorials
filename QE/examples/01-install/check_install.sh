#!/usr/bin/env bash
# check_install.sh - quick sanity check of a Quantum ESPRESSO installation
#
# By default this script does NOT run any calculation: it only checks that
# executables are on PATH, that the input file and pseudopotential exist, and
# that the input has LF line endings. Set RUN=1 to also run the Si smoke test.
#
# Usage:
#   bash check_install.sh                 # static checks only
#   RUN=1 NP=2 bash check_install.sh      # also run pw.x on NP MPI ranks
#   RUN=1 NP=4 NK=2 MPI_LAUNCHER=srun bash check_install.sh

set -u

INPUT=${INPUT:-si_test.in}
OUTPUT=${INPUT%.in}.out
PSEUDO_DIR=${PSEUDO_DIR:-./pseudo}
UPF=${UPF:-Si.pbe-n-rrkjus_psl.1.0.0.UPF}
NP=${NP:-1}
NK=${NK:-1}
MPI_LAUNCHER=${MPI_LAUNCHER:-mpirun}
RUN=${RUN:-0}

status=0
ok()   { echo "[ OK ] $*"; }
warn() { echo "[WARN] $*"; }
fail() { echo "[FAIL] $*"; status=1; }

echo "== 1. executables =="
for exe in pw.x ph.x pp.x dos.x projwfc.x bands.x; do
    if command -v "$exe" >/dev/null 2>&1; then
        ok "$exe -> $(command -v "$exe")"
    else
        if [ "$exe" = "pw.x" ]; then fail "pw.x not found in PATH"; else warn "$exe not found"; fi
    fi
done
if command -v "$MPI_LAUNCHER" >/dev/null 2>&1; then
    ok "launcher $MPI_LAUNCHER -> $(command -v "$MPI_LAUNCHER")"
else
    warn "launcher $MPI_LAUNCHER not found (serial run only)"
fi

echo "== 2. linked libraries (informational) =="
if command -v pw.x >/dev/null 2>&1 && command -v ldd >/dev/null 2>&1; then
    ldd "$(command -v pw.x)" 2>/dev/null \
        | grep -Ei "mkl|openblas|blas|lapack|fftw|scalapack|elpa|mpi|cuda" \
        | awk '{print "       " $1}' | sort -u
fi

echo "== 3. input and pseudopotential =="
if [ -f "$INPUT" ]; then
    ok "input $INPUT found"
    if grep -q $'\r' "$INPUT"; then
        fail "$INPUT has CRLF line endings; fix with: sed -i 's/\\r\$//' $INPUT"
    else
        ok "$INPUT has LF line endings"
    fi
else
    fail "input $INPUT not found"
fi
if [ -f "$PSEUDO_DIR/$UPF" ]; then
    ok "pseudopotential $PSEUDO_DIR/$UPF found"
    grep -o 'pseudo_type="[^"]*"\|functional="[^"]*"\|z_valence="[^"]*"' \
        "$PSEUDO_DIR/$UPF" | head -n 3 | sed 's/^/       /'
else
    fail "missing $PSEUDO_DIR/$UPF (see QE tutorial chapter 3 for download)"
fi

if [ "$RUN" != "1" ]; then
    echo "== static checks finished (set RUN=1 to run pw.x) =="
    exit $status
fi

[ $status -ne 0 ] && { echo "Fix the failures above before running."; exit 1; }

echo "== 4. running pw.x (NP=$NP, NK=$NK) =="
rm -rf ./tmp "$OUTPUT"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
if [ "$NP" -gt 1 ]; then
    # srun uses -n; mpirun/mpiexec accept -np
    case "$(basename "$MPI_LAUNCHER")" in
        srun) NPFLAG="-n" ;;
        *)    NPFLAG="-np" ;;
    esac
    "$MPI_LAUNCHER" "$NPFLAG" "$NP" pw.x -nk "$NK" -i "$INPUT" > "$OUTPUT" 2>&1
else
    pw.x -i "$INPUT" > "$OUTPUT" 2>&1
fi

echo "== 5. checking $OUTPUT =="
grep -E "Program PWSCF" "$OUTPUT" | head -n 1
grep -E "Serial version|Parallel version|Number of MPI processes|Threads/MPI process|K-points division" "$OUTPUT"
grep -q "convergence has been achieved" "$OUTPUT" && ok "SCF converged" || fail "SCF did not converge"
grep -E "^!" "$OUTPUT" || fail "no final total energy line (starting with '!')"
grep -E "highest occupied|Fermi energy" "$OUTPUT"
grep -q "JOB DONE" "$OUTPUT" && ok "JOB DONE found" || fail "JOB DONE missing"
[ -f ./tmp/si_test.save/data-file-schema.xml ] && ok "XML data file written" \
    || warn "tmp/si_test.save/data-file-schema.xml not found"
# A mismatched launcher starts NP independent serial runs: the header then
# appears NP times. ("Number of MPI processes" is printed only by MPI+OpenMP builds.)
nhead=$(grep -c "Program PWSCF" "$OUTPUT")
[ "$nhead" -gt 1 ] && fail "header printed $nhead times: MPI launcher does not match pw.x MPI"

if [ $status -eq 0 ]; then echo "RESULT: PASS"; else echo "RESULT: FAIL"; fi
exit $status
