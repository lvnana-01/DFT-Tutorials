#!/bin/bash
#SBATCH --job-name=hca2-orca
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000M
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
# Run one ORCA input.  Usage:  sbatch sub_orca.slurm orca_opt_reactant.inp
# --ntasks must match "%pal nprocs"; --mem-per-cpu must exceed %maxcore (MB)
# by ~25% (the templates use %maxcore 3000 for def2-SVP and 4000 for def2-TZVP).
# ORCA must be called with its FULL PATH for parallel runs.
set -euo pipefail
INP=${1:?usage: sbatch sub_orca.slurm input.inp}
NAME=${INP%.inp}

# module load orca openmpi     # adapt to your cluster
ORCA_BIN=$(command -v orca)

WORK=${SLURM_TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}
mkdir -p "$WORK"
cp "$INP" ./*.xyz ./*.hess "$WORK"/ 2>/dev/null || true
cd "$WORK"
"$ORCA_BIN" "$INP" > "$SLURM_SUBMIT_DIR/$NAME.out"
cp ./*.xyz ./*.hess ./*.allxyz ./*.gbw "$SLURM_SUBMIT_DIR"/ 2>/dev/null || true

cd "$SLURM_SUBMIT_DIR"
grep "ORCA TERMINATED NORMALLY" "$NAME.out" || echo "[!!] abnormal termination"
grep -A12 "VIBRATIONAL FREQUENCIES" "$NAME.out" | grep -i "imaginary" || true
