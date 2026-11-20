#!/bin/bash
#SBATCH --job-name=qe-ph
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log
#
# Slurm job for the Chapter 12 examples.
#   SYS=si sbatch job_ph_slurm.sh            -> scf + ph.x + q2r.x + matdyn.x
#   SYS=al sbatch job_ph_slurm.sh
#   NI=4 SYS=si sbatch job_ph_slurm.sh       -> ph.x with 4 images (-ni), then
#                                               a short recover run collects the
#                                               dynamical matrices (see Chapter 12)
# Edit the module line and partition/account options for your cluster.
# A more general template is QE/scripts/job_qe_slurm.sh.

set -euo pipefail
# module load quantum-espresso/7.6      # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SYS=${SYS:-si}
NK=${NK:-4}
NI=${NI:-1}

srun pw.x -nk "$NK" -in "$SYS.scf.in" > "$SYS.scf.out"

if [ "$NI" -gt 1 ]; then
    srun ph.x -ni "$NI" -nk "$NK" -in "$SYS.ph.in" > "$SYS.ph.out"
    # Collect the results of all images: same input with recover = .true.
    sed 's/^&INPUTPH/\&INPUTPH\n  recover = .true./' "$SYS.ph.in" > "$SYS.ph_collect.in"
    srun ph.x -nk "$NK" -in "$SYS.ph_collect.in" > "$SYS.ph_collect.out"
    PHOUT="$SYS.ph_collect.out"
else
    srun ph.x -nk "$NK" -in "$SYS.ph.in" > "$SYS.ph.out"
    PHOUT="$SYS.ph.out"
fi

q2r.x -in "$SYS.q2r.in" > "$SYS.q2r.out"
matdyn.x -in "$SYS.matdyn.in" > "$SYS.matdyn.out"
# One 'Diagonalizing' block per irreducible q point is expected
grep -c 'Diagonalizing the dynamical matrix' "$PHOUT" || true
