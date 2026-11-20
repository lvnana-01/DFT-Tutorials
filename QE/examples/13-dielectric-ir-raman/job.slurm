#!/bin/bash
#SBATCH --job-name=qe13-dielec
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log
#
# Chapter 13 workflow (GaAs + H2O). Edit module names for your cluster.
# STEPS selects what to run, e.g.  sbatch --export=ALL,STEPS="gaas_gamma" job.slurm
#   gaas_gamma : scf -> ph.x at Gamma (epsil + Z* + Raman) -> dynmat.x (LO/IR/Raman, then eps0)
#   gaas_disp  : ph.x on a q grid -> q2r.x -> matdyn.x (LO-TO in the dispersion)
#   gaas_eps   : nscf with empty bands -> epsilon.x (overwrites gaas.save wavefunctions)
#   h2o        : relax -> ph.x at Gamma -> dynmat.x (IR + Raman of a molecule)

set -euo pipefail
# module load quantum-espresso/7.6     # site specific

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
RUN="srun -n ${NP}"            # or: mpirun -np ${NP}
STEPS=${STEPS:-"gaas_gamma gaas_disp gaas_eps h2o"}

ls pseudo/Ga.upf pseudo/As.upf pseudo/O.upf pseudo/H.upf > /dev/null

for step in ${STEPS}; do
  case ${step} in
    gaas_gamma)
      ${RUN} pw.x -nk 4 -in gaas.scf.in > gaas.scf.out
      ${RUN} ph.x -nk 4 -in gaas.ph_epsil.in > gaas.ph_epsil.out
      dynmat.x -in gaas.dynmat.in      > gaas.dynmat.out        # LO along q=(1,0,0), IR, Raman
      dynmat.x -in gaas.dynmat_eps0.in > gaas.dynmat_eps0.out   # eps0 from TO modes (lperm)
      ;;
    gaas_disp)
      [ -d tmp/gaas.save ] || ${RUN} pw.x -nk 4 -in gaas.scf.in > gaas.scf.out
      ${RUN} ph.x -nk 4 -in gaas.ph_disp.in > gaas.ph_disp.out
      q2r.x -in gaas.q2r.in > gaas.q2r.out
      matdyn.x -in gaas.matdyn.in > gaas.matdyn.out
      ;;
    gaas_eps)
      [ -d tmp/gaas.save ] || ${RUN} pw.x -nk 4 -in gaas.scf.in > gaas.scf.out
      ${RUN} pw.x -nk 4 -in gaas.nscf_eps.in > gaas.nscf_eps.out
      ${RUN} epsilon.x -in gaas.eps.in > gaas.eps.out
      ;;
    h2o)
      ${RUN} pw.x -in h2o.relax.in > h2o.relax.out
      ${RUN} ph.x -in h2o.ph_gamma.in > h2o.ph_gamma.out
      dynmat.x -in h2o.dynmat.in > h2o.dynmat.out
      ;;
    *) echo "unknown step ${step}"; exit 1 ;;
  esac
done
grep -c "JOB DONE" *.out || true
