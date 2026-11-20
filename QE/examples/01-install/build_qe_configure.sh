#!/usr/bin/env bash
# build_qe_configure.sh - example configure + make recipe for Quantum ESPRESSO 7.6
#
# READ AND EDIT BEFORE USE. This script was written as documentation and has
# never been executed by the repository authors.
#
# Usage:
#   bash build_qe_configure.sh                    # CPU, MPI + OpenMP
#   MODE=gpu CUDA_CC=90 bash build_qe_configure.sh  # NVIDIA GPU (NVHPC SDK)

set -euo pipefail

QE_VERSION=${QE_VERSION:-7.6}
SRC_DIR=${SRC_DIR:-$HOME/src/qe-${QE_VERSION}}
PREFIX=${PREFIX:-$HOME/opt/qe-${QE_VERSION}}
MODE=${MODE:-cpu}
CUDA_CC=${CUDA_CC:-80}
NJOBS=${NJOBS:-8}
TARGETS=${TARGETS:-"pw ph pp"}

export LC_ALL=C     # QE scripts assume the C locale

if [ ! -d "$SRC_DIR" ]; then
    git clone --branch "qe-${QE_VERSION}" --depth 1 \
        https://gitlab.com/QEF/q-e.git "$SRC_DIR"
fi

# Out-of-source build (supported by autoconf since QE 7.6)
BUILD_DIR="$SRC_DIR/build-make-${MODE}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

if [ "$MODE" = "gpu" ]; then
    # Requires the NVIDIA HPC SDK environment (e.g. 'module load nvhpc')
    ../configure --prefix="$PREFIX" --enable-openmp \
        --with-cuda="${NVHPC_CUDA_HOME:?set NVHPC_CUDA_HOME}" \
        --with-cuda-cc="$CUDA_CC" \
        --with-cuda-mpi=yes
else
    ../configure --prefix="$PREFIX" --enable-openmp --with-scalapack=yes
fi

echo "=== configure summary ==="
cat configure.msg || true
echo "=== key make.inc lines ==="
grep -E "^(DFLAGS|MPIF90|FFLAGS|BLAS_LIBS|LAPACK_LIBS|SCALAPACK_LIBS|FFT_LIBS)" make.inc || true

# Stop here if you want to inspect make.inc first:
if [ "${CONFIGURE_ONLY:-0}" = "1" ]; then
    exit 0
fi

# shellcheck disable=SC2086
make -j "$NJOBS" $TARGETS
make install

ls "$BUILD_DIR/bin"
echo "Done. Add $PREFIX/bin to PATH."
