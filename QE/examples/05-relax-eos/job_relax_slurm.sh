#!/bin/bash
#SBATCH --job-name=qe-relax
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log
#
# Slurm job for the Chapter 5 examples. Edit the module line and the
# partition/account options for your cluster. See also
# QE/scripts/job_qe_slurm.sh for a more general template.

set -euo pipefail
# module load quantum-espresso/7.6      # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

INPUT=${INPUT:-si.vc-relax.in}
OUTPUT=${INPUT%.in}.out

# k-point pools: Si 8x8x8 has enough irreducible k points for -nk 4.
# Use -nk 1 for the Gamma-only CO molecule.
NK=${NK:-4}

srun pw.x -nk "$NK" -in "$INPUT" > "$OUTPUT"

grep -e 'bfgs converged' -e 'The maximum number of steps' -e 'JOB DONE' "$OUTPUT" || true
