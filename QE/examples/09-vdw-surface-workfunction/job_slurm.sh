#!/bin/bash
#SBATCH --job-name=qe-surf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.log
#
# Slurm wrapper for the Chapter 9 examples (edit module/partition lines).
# General template: QE/scripts/job_qe_slurm.sh
#
# Usage:
#   sbatch job_slurm.sh                              # work function chain
#   sbatch --export=ALL,WHAT=vdw job_slurm.sh        # graphite vdW comparison
#   sbatch --export=ALL,WHAT=mos2 job_slurm.sh       # MoS2 monolayer 2Dxy relax

set -euo pipefail
# module load quantum-espresso/7.6     # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

NK=${NK:-4}
export PW_CMD="srun pw.x -nk ${NK}"
export PP_CMD="srun"

case "${WHAT:-wf}" in
  wf)   bash run_workfunction.sh ;;
  vdw)  bash run_vdw_compare.sh ;;
  mos2) $PW_CMD -in mos2_mono.vc-relax.in > mos2_mono.vc-relax.out
        grep -e 'bfgs converged' -e 'P=' mos2_mono.vc-relax.out | tail -2 ;;
  *)    echo "unknown WHAT=${WHAT}" >&2; exit 1 ;;
esac
