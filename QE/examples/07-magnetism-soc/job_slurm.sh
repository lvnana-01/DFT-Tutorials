#!/bin/bash
#SBATCH --job-name=qe-mag
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log
#
# Slurm wrapper for run_all.sh (Chapter 7). Edit the module line and the
# partition/account options for your cluster. A more general template is
# QE/scripts/job_qe_slurm.sh.
#
# Usage:  sbatch job_slurm.sh            # all tasks
#         sbatch --export=ALL,TASK=mae job_slurm.sh

set -euo pipefail
# module load quantum-espresso/7.6     # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

NK=${NK:-4}
export PW_CMD="srun pw.x -nk ${NK}"
export PP_CMD="srun"
bash run_all.sh "${TASK:-all}"
