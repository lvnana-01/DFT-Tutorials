#!/bin/bash
# Slurm job for chapter 17: pw.x BOMD of Si64 (NVT) and cp.x water (damp -> verlet).
# Usage: sbatch job_md_slurm.sh pw      # run si64.md.in
#        sbatch job_md_slurm.sh pwnve   # run si64.nve.in
#        sbatch job_md_slurm.sh cp      # run h2o.cp_damp.in then h2o.cp_verlet.in
#SBATCH --job-name=qe-md
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --output=%x-%j.log

set -euo pipefail
# module load quantum-espresso/7.6      # site specific
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-32}
RUN="srun -n ${NP}"
MODE=${1:-pw}

cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p tmp
# A leftover prefix.EXIT file makes pw.x/cp.x stop right after initialisation: remove it.
rm -f tmp/si64.EXIT tmp/si64nve.EXIT tmp/h2o.EXIT si64.EXIT si64nve.EXIT h2o.EXIT

# Ten minutes before the time limit Slurm sends USR1 to this script; we create
# the EXIT file so that pw.x/cp.x stop cleanly and can be restarted.
trap 'echo "USR1 received, requesting clean stop"; touch tmp/si64.EXIT tmp/si64nve.EXIT tmp/h2o.EXIT' USR1

# Run a QE program in the background so that the USR1 trap can fire while it runs.
# "wait" returns early when a signal arrives, so keep waiting until the process exits.
run_qe() {
  local out=$1; shift
  "$@" > "${out}" &
  local pid=$!
  until wait "${pid}"; do
    kill -0 "${pid}" 2>/dev/null || break
  done
}

case "${MODE}" in
  pw)
    # Gamma-point run: no k-point pools; use plane-wave parallelization only.
    run_qe si64.md.out ${RUN} pw.x -in si64.md.in
    ;;
  pwnve)
    run_qe si64.nve.out ${RUN} pw.x -in si64.nve.in
    ;;
  cp)
    run_qe h2o.cp_damp.out   ${RUN} cp.x -in h2o.cp_damp.in
    run_qe h2o.cp_verlet.out ${RUN} cp.x -in h2o.cp_verlet.in
    ;;
  *)
    echo "unknown mode ${MODE} (use pw, pwnve or cp)"; exit 1
    ;;
esac
