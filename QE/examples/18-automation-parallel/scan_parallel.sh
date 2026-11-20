#!/bin/bash
# Generate (and optionally submit) a small grid of pw.x runs with different
# -nk / OMP_NUM_THREADS settings for the same input (add -nd by hand for large
# systems), then summarise the
# wall time and memory estimate. It never runs pw.x by itself unless SUBMIT=1.
#
# Usage:
#   bash scan_parallel.sh si.scf.in              # create scan_*/ directories and job files
#   SUBMIT=1 bash scan_parallel.sh si.scf.in     # also sbatch every job
#   bash scan_parallel.sh collect                # summarise finished runs
set -euo pipefail

NTOT=${NTOT:-16}                # total cores per job (ranks x threads)
NK_LIST=${NK_LIST:-"1 2 4 8"}
NT_LIST=${NT_LIST:-"1 2"}      # OpenMP threads per rank

if [[ "${1:-}" == "collect" ]]; then
  printf "%-22s %12s %16s\n" "run" "PWSCF_WALL" "max_RAM_estimate"
  for d in scan_*/; do
    out="${d}pw.out"
    [[ -f "${out}" ]] || continue
    wall=$(grep "PWSCF *:.*WALL" "${out}" | tail -1 | sed "s/.*CPU *//; s/ *WALL.*//")
    ram=$(grep "Estimated max dynamical RAM per process" "${out}" | awk '{print $(NF-1), $NF}')
    printf "%-22s %12s %16s\n" "${d%/}" "${wall:-n/a}" "${ram:-n/a}"
  done
  exit 0
fi

INPUT=${1:?usage: bash scan_parallel.sh input.in | collect}
for nt in ${NT_LIST}; do
  nranks=$(( NTOT / nt ))
  for nk in ${NK_LIST}; do
    (( nranks % nk == 0 )) || continue
    dir="scan_nk${nk}_omp${nt}"
    mkdir -p "${dir}"
    cp "${INPUT}" "${dir}/pw.in"
    [[ -d pseudo ]] && ln -sfn ../pseudo "${dir}/pseudo"
    cat > "${dir}/job.sh" <<JOB
#!/bin/bash
#SBATCH --job-name=${dir}
#SBATCH --nodes=1
#SBATCH --ntasks=${nranks}
#SBATCH --cpus-per-task=${nt}
#SBATCH --time=00:30:00
export OMP_NUM_THREADS=${nt}
srun --cpu-bind=cores pw.x -nk ${nk} -in pw.in > pw.out
JOB
    if [[ "${SUBMIT:-0}" == "1" ]]; then
      (cd "${dir}" && sbatch job.sh)
    fi
    echo "prepared ${dir} (ranks=${nranks}, -nk ${nk}, OMP=${nt})"
  done
done
