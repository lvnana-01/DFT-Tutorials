#!/bin/bash
#SBATCH --job-name=hca2-g16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
# Run one Gaussian 16 input.  Usage:  sbatch sub_g16.slurm opt_reactant.gjf
# --cpus-per-task must match %nprocshared; --mem must exceed %mem by a few GB
# (the templates use %mem=24GB for opt/TS/IRC/freq and 32GB for def2-TZVP/ONIOM).
set -euo pipefail
INP=${1:?usage: sbatch sub_g16.slurm input.gjf}
NAME=${INP%.gjf}

# module load gaussian/16      # adapt to your cluster
export GAUSS_SCRDIR=${SLURM_TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}
mkdir -p "$GAUSS_SCRDIR"

g16 < "$INP" > "$NAME.log"

grep -c "Normal termination" "$NAME.log" || true
grep "imaginary frequencies" "$NAME.log" || echo "no imaginary frequencies reported"
