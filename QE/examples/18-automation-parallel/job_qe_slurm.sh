#!/bin/bash
# Annotated Slurm template for pw.x (CPU, MPI + optional OpenMP) used in chapter 18.
# The repository-wide template is QE/scripts/job_qe_slurm.sh; this copy is kept
# next to the chapter so every line can be discussed in the text.
#
# Usage: sbatch job_qe_slurm.sh si.scf.in
#SBATCH --job-name=qe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16       # MPI ranks per node
#SBATCH --cpus-per-task=2          # OpenMP threads per rank (1 = pure MPI)
#SBATCH --time=01:00:00
#SBATCH --signal=B:USR1@300        # send USR1 to this script 5 min before the limit
#SBATCH --output=%x-%j.log

set -euo pipefail
INPUT=${1:?usage: sbatch job_qe_slurm.sh input.in}
OUTPUT=${INPUT%.in}.out
PREFIX=$(sed -n "s/^ *prefix *= *'\([^']*\)'.*/\1/p" "${INPUT}" | head -1)
PREFIX=${PREFIX:-pwscf}

# module load quantum-espresso/7.6           # site specific
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
# export ESPRESSO_PSEUDO=$HOME/pseudo        # default pseudo_dir if not set in the input
# export ESPRESSO_TMPDIR=/scratch/$USER/$SLURM_JOB_ID   # default outdir if not set in the input

NK=${NK:-4}                                  # k-point pools: divisor of ranks, <= number of k points
NP=${SLURM_NTASKS:-16}

cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p tmp
rm -f "tmp/${PREFIX}.EXIT" "${PREFIX}.EXIT"   # a stale EXIT file stops the run immediately

# Clean stop: pw.x checks for prefix.EXIT in outdir and in the working directory.
trap 'echo "USR1: requesting clean stop"; touch "tmp/${PREFIX}.EXIT"' USR1

srun -n "${NP}" --cpu-bind=cores pw.x -nk "${NK}" -in "${INPUT}" > "${OUTPUT}" &
PID=$!
until wait "${PID}"; do
  kill -0 "${PID}" 2>/dev/null || break
done

grep -q "JOB DONE" "${OUTPUT}" && echo "finished: ${OUTPUT}" || echo "check ${OUTPUT}: no JOB DONE"
