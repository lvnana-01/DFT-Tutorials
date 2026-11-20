#!/bin/bash
# Slurm job: EPW route for fcc Pb (scf -> ph.x -> pp.py -> q2r -> nscf -> epw1 -> epw2).
# Requires epw.x (built with "make epw"), EPW/bin/pp.py and wannier90/utility/kmesh.pl.
# Adjust paths, partition and modules to your cluster.
#SBATCH --job-name=pb-epw
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
# module load quantum-espresso/7.6      # site specific, must provide epw.x
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-32}
RUN="srun -n ${NP}"
QE_SRC=${QE_SRC:-$HOME/src/q-e-qe-7.6}  # QE source tree that contains EPW/bin/pp.py
KMESH=${KMESH:-$HOME/src/wannier90/utility/kmesh.pl}  # kmesh.pl from the Wannier90 distribution

cd "${SLURM_SUBMIT_DIR:-.}"
# All inputs use outdir = './': pp.py expects _ph0/ and pb.save/ in this directory.

${RUN} pw.x -nk 8 -in pb.scf.in > pb.scf.out
${RUN} ph.x -nk 8 -in pb.ph.in  > pb.ph.out

# Collect dyn, dvscf and patterns into ./save (pp.py asks for the prefix on stdin)
echo "pb" | python3 "${QE_SRC}/EPW/bin/pp.py"

q2r.x -in pb.q2r.in > pb.q2r.out

# Build the nscf input with the full 6x6x6 k list (216 points, no symmetry reduction)
cat pb.nscf.in.head > pb.nscf.in
perl "${KMESH}" 6 6 6 >> pb.nscf.in
${RUN} pw.x -nk 8 -in pb.nscf.in > pb.nscf.out

# EPW: -nk pools over k points; nprocs must equal npool x nimage
${RUN} epw.x -nk ${NP} -in pb.epw1.in > pb.epw1.out
${RUN} epw.x -nk ${NP} -in pb.epw2.in > pb.epw2.out
