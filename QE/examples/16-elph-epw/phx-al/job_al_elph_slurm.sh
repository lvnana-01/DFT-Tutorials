#!/bin/bash
# Slurm job: ph.x electron-phonon route for fcc Al (pw.x -> ph.x -> q2r.x -> matdyn.x -> lambda.x).
# Adjust partition/account/modules to your cluster. Submit with: sbatch job_al_elph_slurm.sh
#SBATCH --job-name=al-elph
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
# module load quantum-espresso/7.6      # site specific
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-16}
RUN="srun -n ${NP}"                      # or: mpirun -np ${NP}

cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p tmp

# 1) dense-k SCF that writes a2Fsave, 2) phonon-grid SCF (same prefix/outdir)
${RUN} pw.x -nk 4 -in al.scf.fit.in > al.scf.fit.out
${RUN} pw.x -nk 4 -in al.scf.in     > al.scf.out
# 3) phonons + el-ph on the 4x4x4 q grid (el-ph does not support -ni images)
${RUN} ph.x -nk 4 -in al.elph.in    > al.elph.out
# 4) real-space IFCs and interpolated linewidths (serial is enough)
q2r.x    -in al.q2r.in          > al.q2r.out
matdyn.x -in al.matdyn.disp.in  > al.matdyn.disp.out
matdyn.x -in al.matdyn.dos.in   > al.matdyn.dos.out
# 5) lambda, omega_log and Allen-Dynes Tc from the per-q files in elph_dir/
lambda.x < al.lambda.in > al.lambda.out

grep -c "JOB DONE" al.scf.fit.out al.scf.out al.elph.out || true
