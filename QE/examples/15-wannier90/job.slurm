#!/bin/bash
#SBATCH --job-name=qe15-w90
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log
#
# Chapter 15: pw.x scf -> pw.x nscf (full grid) -> wannier90.x -pp -> pw2wannier90.x -> wannier90.x
# SEEDS selects which Wannier runs to do (both reuse the same nscf):
#   sbatch --export=ALL,SEEDS="si_val" job.slurm

set -euo pipefail
# module load quantum-espresso/7.6 wannier90/3.1     # site specific

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
RUN="srun -n ${NP}"            # or: mpirun -np ${NP}
SEEDS=${SEEDS:-"si_val si"}

ls pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF > /dev/null
command -v wannier90.x > /dev/null
command -v pw2wannier90.x > /dev/null

${RUN} pw.x -nk 4 -in si.scf.in  > si.scf.out
${RUN} pw.x -nk 4 -in si.nscf.in > si.nscf.out

for seed in ${SEEDS}; do
  wannier90.x -pp ${seed}                                   # -> ${seed}.nnkp
  ${RUN} pw2wannier90.x -in ${seed}.pw2wan.in > ${seed}.pw2wan.out
  wannier90.x ${seed}                                       # -> ${seed}.wout, ${seed}_band.dat
  grep -A6 "Final State" ${seed}.wout | tail -4 || true
done
