#!/bin/bash
#SBATCH --job-name=qe-neb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log
#
# Slurm job for the Chapter 11 NEB example.
#   STEP=relax sbatch job_neb_slurm.sh   -> relax both end points
#   STEP=neb   sbatch job_neb_slurm.sh   -> CI-NEB with image parallelization
# NI = number of image groups (-ni); the ranks are split evenly among them.
# Edit the module line and partition/account options for your cluster.
# A more general template is QE/scripts/job_qe_slurm.sh.

set -euo pipefail
# module load quantum-espresso/7.6      # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

STEP=${STEP:-neb}
NI=${NI:-5}

case "$STEP" in
  relax)
    for s in initial final; do
      srun pw.x -nk 4 -in h_al100.$s.relax.in > h_al100.$s.relax.out
      grep -e 'bfgs converged' -e 'JOB DONE' h_al100.$s.relax.out || true
    done
    ;;
  neb)
    # neb.x does NOT read standard input: the input must be given with -inp
    srun neb.x -ni "$NI" -inp h_al100.neb.in > h_al100.neb.out
    grep -e 'activation energy' -e 'neb: convergence achieved' \
         -e 'reached the maximum number' h_al100.neb.out || true
    ;;
  *)
    echo "unknown STEP=$STEP (use relax or neb)"; exit 1 ;;
esac
