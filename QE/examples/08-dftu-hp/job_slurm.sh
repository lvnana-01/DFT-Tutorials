#!/bin/bash
#SBATCH --job-name=qe-hp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.log
#
# Slurm wrapper for the Chapter 8 examples. Edit the module line and the
# partition/account options for your cluster (general template:
# QE/scripts/job_qe_slurm.sh).
#
# Usage:
#   sbatch job_slurm.sh                                   # NiO, 1 cycle
#   sbatch --export=ALL,NCYCLE=5 job_slurm.sh             # self-consistent U
#   sbatch --export=ALL,WHAT=uv job_slurm.sh              # LiCoO2 U+V
#   sbatch --export=ALL,WHAT=dftu job_slurm.sh            # NiO DFT+U scf only

set -euo pipefail
# module load quantum-espresso/7.6     # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

NK=${NK:-4}
export PW_CMD="srun pw.x -nk ${NK}"
export HP_CMD="srun hp.x -nk ${NK}"

case "${WHAT:-hp}" in
  hp)   bash run_hp.sh ;;
  uv)   bash run_uv.sh ;;
  dftu) $PW_CMD -in nio.dftu.scf.in > nio.dftu.scf.out
        grep -e '^!' -e 'Hubbard energy' -e 'highest occupied' nio.dftu.scf.out || true ;;
  *)    echo "unknown WHAT=${WHAT}" >&2; exit 1 ;;
esac
