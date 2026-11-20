#!/bin/bash
#SBATCH --job-name=qe14-hybrid
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.log
#
# Chapter 14 workflow. Edit module names for your cluster.
#   sbatch --export=ALL,STEPS="pbe hse" job.slurm
# Steps:
#   pbe     : Si PBE reference
#   hse     : Si HSE on the 4x4x4 mesh (gap on the mesh only)
#   pbe0    : MgO PBE0
#   r2scan  : Si r2SCAN (needs a Libxc-enabled build)
#   wannier : Si HSE on explicit full grid -> wannier90 -pp -> pw2wannier90 -> wannier90
#             (needs wannier90.x and pw2wannier90.x in PATH)

set -euo pipefail
# module load quantum-espresso/7.6 wannier90/3.1     # site specific

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
RUN="srun -n ${NP}"            # or: mpirun -np ${NP}
STEPS=${STEPS:-"pbe hse pbe0"}

for step in ${STEPS}; do
  case ${step} in
    pbe)    ${RUN} pw.x -nk 4 -in si.pbe.scf.in    > si.pbe.scf.out ;;
    hse)    ${RUN} pw.x -nk 4 -in si.hse.scf.in    > si.hse.scf.out ;;
    pbe0)   ${RUN} pw.x -nk 4 -in mgo.pbe0.scf.in  > mgo.pbe0.scf.out ;;
    r2scan) ${RUN} pw.x -nk 4 -in si.r2scan.scf.in > si.r2scan.scf.out ;;
    wannier)
      ${RUN} pw.x -nk 4 -in si.hse.fullgrid.in > si.hse.fullgrid.out
      wannier90.x -pp si_hse
      ${RUN} pw2wannier90.x -in si_hse.pw2wan.in > si_hse.pw2wan.out
      wannier90.x si_hse
      ;;
    *) echo "unknown step ${step}"; exit 1 ;;
  esac
done
grep -H -e "highest occupied" -e "EXX self-consistency reached" *.out || true
