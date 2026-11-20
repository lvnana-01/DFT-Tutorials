#!/usr/bin/env bash
#SBATCH --job-name=qe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32          # MPI ranks per node
#SBATCH --cpus-per-task=1             # OpenMP threads per rank (QE built with OpenMP only)
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal           # <- remove one '#' and edit for your cluster
##SBATCH --account=myproject
##SBATCH --gpus-per-node=4            # GPU build: one MPI rank per GPU (ntasks-per-node=4)
#
# job_qe_slurm.sh - Slurm job template for Quantum ESPRESSO (pw.x, ph.x, pp.x, ...)
#
# Usage:
#   sbatch job_qe_slurm.sh si.scf.in                         # input as argument (pw.x)
#   sbatch --export=ALL,INP=si.scf.in,NK=4 job_qe_slurm.sh   # or as variables
#   sbatch --export=ALL,EXE=ph.x,INP=si.ph.in,NK=4 job_qe_slurm.sh
#   sbatch --export=ALL,EXE=neb.x,INP=neb.in,NI=4 job_qe_slurm.sh   # images (-ni)
#   bash job_qe_slurm.sh -h                                  # this help
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
# defaults itself.  OpenMP: export OMP_NUM_THREADS = cpus-per-task.
#
# What the script does: loads the QE environment (edit block below), sets OpenMP variables,
# keeps an existing output as <out>.<n>.bak, runs "srun EXE <flags> -in INP > OUT", then checks
# for "JOB DONE." and prints convergence problems ("convergence NOT achieved", "%%%%" errors).
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# job_qe_slurm.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
INP=${1:-${INP:-pw.in}}
OUT=${OUT:-${INP%.*}.out}
EXE=${EXE:-pw.x}
NK=${NK:-}            # e.g. 4
ND=${ND:-}            # e.g. 16 (square number)
NT=${NT:-}
NI=${NI:-}
EXTRA=${EXTRA:-}
LAUNCHER=${LAUNCHER:-srun}            # srun (recommended) or mpirun

# ---- environment: keep exactly ONE block that matches your installation ----
# module load quantum-espresso/7.5                 # cluster module (check: module avail espresso)
# source "$HOME/miniconda3/bin/activate" qe        # conda-forge package (conda install -c conda-forge qe)
# export PATH=$HOME/q-e-qe-7.5/bin:$PATH           # own source build
# export ESPRESSO_PSEUDO=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency
export ESPRESSO_TMPDIR=${ESPRESSO_TMPDIR:-${SLURM_TMPDIR:-$PWD/tmp}}
# ==========================================================================

cd "${SLURM_SUBMIT_DIR:-$PWD}"
[ -f "$INP" ] || { echo "ERROR: input file '$INP' not found in $(pwd)" >&2; exit 1; }
command -v "$EXE" >/dev/null 2>&1 || {
    echo "ERROR: $EXE not in PATH - edit the environment block of $0" >&2; exit 1; }
if grep -q $'\r' "$INP"; then
    echo "ERROR: $INP has Windows (CRLF) line endings - run: dos2unix $INP" >&2; exit 1
fi
mkdir -p "$ESPRESSO_TMPDIR"

NTHREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export SRUN_CPUS_PER_TASK=$NTHREADS     # Slurm >= 22.05: srun does not inherit --cpus-per-task

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

NRANKS=${SLURM_NTASKS:-1}
echo "========================================================"
echo "Job ID        : ${SLURM_JOB_ID:-none}   Nodes: ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "Executable    : $(command -v "$EXE")"
echo "Input/Output  : $INP -> $OUT"
echo "MPI x OpenMP  : $NRANKS x $OMP_NUM_THREADS   flags:${FLAGS}"
echo "TMPDIR        : $ESPRESSO_TMPDIR"
echo "Start         : $(date)"
echo "========================================================"

set +e
if [ "$LAUNCHER" = "srun" ]; then
    # shellcheck disable=SC2086
    srun --cpu-bind=cores "$EXE" $FLAGS -in "$INP" > "$OUT"
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
