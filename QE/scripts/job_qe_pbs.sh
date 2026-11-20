#!/usr/bin/env bash
#PBS -N qe
#PBS -l select=1:ncpus=32:mpiprocs=32:ompthreads=1
#PBS -l walltime=24:00:00
#PBS -j oe
##PBS -q normal                        # <- remove one '#' and edit for your cluster
##PBS -A myproject
## Torque / OpenPBS without 'select': use instead
##PBS -l nodes=1:ppn=32
#
# job_qe_pbs.sh - PBS Pro / OpenPBS / Torque job template for Quantum ESPRESSO (pw.x, ph.x, ...)
#
# Usage:
#   qsub -v INP=si.scf.in job_qe_pbs.sh                 # PBS does not pass script arguments
#   qsub -v INP=si.scf.in,NK=4 job_qe_pbs.sh
#   qsub -v EXE=ph.x,INP=si.ph.in,NK=4 job_qe_pbs.sh
#   bash job_qe_pbs.sh -h                               # this help
#
# MPI ranks are taken from $PBS_NODEFILE (one line per rank with PBS Pro mpiprocs; one line
# per core with Torque ppn) divided by NTHREADS (default: $OMP_NUM_THREADS set by
# ompthreads, else 1).  Override with NRANKS=... if your site differs.
#
# Variables (environment or edit below):
#   INP   input file (default: first argument, else pw.in);  OUT  output (default <INP>.out)
#   EXE   QE executable (pw.x)          NK  pools  -nk  (k-point parallelisation)
#   ND    -nd (ScaLAPACK diag. group, must be a square number)   NT  -nt task groups
#   NI    -ni images (neb.x, ph.x images)                         EXTRA  more flags
#   ESPRESSO_TMPDIR  scratch for outdir if outdir is not set in the input
#   ESPRESSO_PSEUDO  pseudopotential directory if pseudo_dir is not set in the input
# Rules of thumb: NK should divide the number of MPI ranks and not exceed the number of
# irreducible k points; ranks per pool = ranks / NK.  With no -nk/-nd QE >= 7.1 chooses
# defaults itself.  OpenMP: OMP_NUM_THREADS = ompthreads (PBS Pro) or NTHREADS.
#
# What the script does: loads the QE environment (edit block below), sets OpenMP variables,
# keeps an existing output as <out>.<n>.bak, runs "mpirun -np N EXE <flags> -in INP > OUT", then checks
# for "JOB DONE." and prints convergence problems ("convergence NOT achieved", "%%%%" errors).
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# job_qe_pbs.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
INP=${INP:-${1:-pw.in}}
OUT=${OUT:-${INP%.*}.out}
EXE=${EXE:-pw.x}
NK=${NK:-}            # e.g. 4
ND=${ND:-}            # e.g. 16 (square number)
NT=${NT:-}
NI=${NI:-}
EXTRA=${EXTRA:-}
LAUNCHER=${LAUNCHER:-mpirun}          # mpirun / mpiexec

# ---- environment: keep exactly ONE block that matches your installation ----
# module load quantum-espresso/7.5                 # cluster module (check: module avail espresso)
# source "$HOME/miniconda3/bin/activate" qe        # conda-forge package (conda install -c conda-forge qe)
# export PATH=$HOME/q-e-qe-7.5/bin:$PATH           # own source build
# export ESPRESSO_PSEUDO=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency
export ESPRESSO_TMPDIR=${ESPRESSO_TMPDIR:-${TMPDIR:-${PBS_O_WORKDIR:-$PWD}/tmp}}
# ==========================================================================

cd "${PBS_O_WORKDIR:-$PWD}"
[ -f "$INP" ] || { echo "ERROR: input file '$INP' not found in $(pwd)" >&2; exit 1; }
command -v "$EXE" >/dev/null 2>&1 || {
    echo "ERROR: $EXE not in PATH - edit the environment block of $0" >&2; exit 1; }
if grep -q $'\r' "$INP"; then
    echo "ERROR: $INP has Windows (CRLF) line endings - run: dos2unix $INP" >&2; exit 1
fi
mkdir -p "$ESPRESSO_TMPDIR"

NTHREADS=${NTHREADS:-${OMP_NUM_THREADS:-1}}
export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close
if [ -z "${NRANKS:-}" ]; then
    if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
        NRANKS=$(( $(wc -l < "$PBS_NODEFILE") / NTHREADS ))
    else
        NRANKS=1
    fi
fi
[ "$NRANKS" -ge 1 ] || NRANKS=1

FLAGS=""
[ -n "$NI" ] && FLAGS="$FLAGS -ni $NI"
[ -n "$NK" ] && FLAGS="$FLAGS -nk $NK"
[ -n "$NT" ] && FLAGS="$FLAGS -nt $NT"
[ -n "$ND" ] && FLAGS="$FLAGS -nd $ND"
FLAGS="$FLAGS $EXTRA"

if [ -f "$OUT" ]; then
    n=1; while [ -e "$OUT.$n.bak" ]; do n=$((n + 1)); done
    mv "$OUT" "$OUT.$n.bak"
    echo "previous output moved to $OUT.$n.bak"
fi

echo "========================================================"
echo "Job ID        : ${PBS_JOBID:-none}   Host: $(hostname)"
echo "Executable    : $(command -v "$EXE")"
echo "Input/Output  : $INP -> $OUT"
echo "MPI x OpenMP  : $NRANKS x $OMP_NUM_THREADS   flags:${FLAGS}"
echo "TMPDIR        : $ESPRESSO_TMPDIR"
echo "Start         : $(date)"
echo "========================================================"

set +e
if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
    # shellcheck disable=SC2086
    $LAUNCHER -np "$NRANKS" -machinefile "$PBS_NODEFILE" "$EXE" $FLAGS -in "$INP" > "$OUT"
else
    # shellcheck disable=SC2086
    $LAUNCHER -np "$NRANKS" "$EXE" $FLAGS -in "$INP" > "$OUT"
fi
rc=$?
set -e
echo "End           : $(date)   exit code $rc"

status=0
if grep -q "JOB DONE" "$OUT"; then
    echo "OK: $EXE finished ($OUT)"
    grep '^!' "$OUT" | tail -n 1 || true
    grep -E "the Fermi energy is|highest occupied" "$OUT" | tail -n 1 || true
    grep -E "bfgs converged|End of BFGS|The maximum number of steps" "$OUT" | tail -n 2 || true
    if grep -q "convergence NOT achieved" "$OUT"; then
        echo "WARNING: some SCF cycle did NOT converge (see $OUT)"; status=1
    fi
else
    status=1
    echo "PROBLEM: no 'JOB DONE.' in $OUT (crash, error or wall-time limit)." >&2
    grep -A4 "^ %%%%%%" "$OUT" | head -n 8 >&2 || true
    grep -q "convergence NOT achieved" "$OUT" && \
        echo "-> SCF not converged: lower mixing_beta, more electron_maxstep, smearing (chapter 4/18)." >&2
    echo "-> wall time: set max_seconds in &CONTROL so pw.x stops cleanly, then rerun with" >&2
    echo "   restart_mode='restart' (same prefix/outdir); summary: python qe_parse_out.py $OUT" >&2
fi
[ "$rc" -ne 0 ] && status=1
exit "$status"
