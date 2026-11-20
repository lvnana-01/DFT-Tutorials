#!/bin/bash
#SBATCH --job-name=hca2-xtb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# GFN2-xTB pre-optimisation of the hCA II cluster models.
# Usage:  sbatch run_xtb.sh            (or: bash run_xtb.sh on a workstation)
# Expects ../00-structure/build/cluster.xyz        (charge +1, no CO2)
#     and ../00-structure/build_co2/cluster.xyz    (charge +1, with CO2)
# Edit CHRG to 0 if the models were built with --with-glu106.
set -euo pipefail

# module load xtb        # adapt to your cluster
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_STACKSIZE=4G
ulimit -s unlimited || true

CHRG=1
UHF=0           # number of unpaired electrons; Zn(II) d10 -> 0
SOLV=ether      # ALPB solvent with eps ~4, a crude "protein-like" medium
HERE=$(pwd)

run_opt () {   # $1 = input xyz, $2 = work dir, $3 = xcontrol file
    mkdir -p "$2"
    cp "$1" "$2/cluster.xyz"
    cp "$3" "$2/xcontrol.inp"
    cd "$2"
    xtb cluster.xyz --gfn 2 --chrg "$CHRG" --uhf "$UHF" --alpb "$SOLV" \
        --opt tight --input xcontrol.inp > xtb.out 2>&1 || true
    # with "set -e" a failed xtb run would otherwise stop the whole script
    if grep -q "GEOMETRY OPTIMIZATION CONVERGED" xtb.out && [ ! -e NOT_CONVERGED ]; then
        echo "[OK] $2 converged"
    else
        echo "[!!] $2 not converged - check xtb.out"
    fi
    grep "TOTAL ENERGY" xtb.out | tail -1 || true
    cd "$HERE"
}

# 1) cluster without substrate (resting state Zn-OH-)
run_opt ../00-structure/build/cluster.xyz     xtb_resting  xcontrol_fix.inp
# 2) cluster with CO2 (reactant / Michaelis complex)
run_opt ../00-structure/build_co2/cluster.xyz xtb_reactant xcontrol_fix.inp

# 3) optional cheap relaxed scan O...C, started from the pre-optimised reactant
#    (atom numbers in xcontrol_scan.inp refer to the default CO2 model)
if [ ! -s xtb_reactant/xtbopt.xyz ]; then
    echo "[!!] no xtb_reactant/xtbopt.xyz - skipping the scan"
    exit 0
fi
mkdir -p xtb_scan
cp xtb_reactant/xtbopt.xyz xtb_scan/start.xyz
cp xcontrol_scan.inp xtb_scan/xcontrol.inp
cd xtb_scan
xtb start.xyz --gfn 2 --chrg "$CHRG" --uhf "$UHF" --alpb "$SOLV" \
    --opt --input xcontrol.inp > xtb_scan.out 2>&1 || true
ls -l xtbscan.log 2>/dev/null || echo "[!!] no xtbscan.log - check xtb_scan.out"
cd "$HERE"
