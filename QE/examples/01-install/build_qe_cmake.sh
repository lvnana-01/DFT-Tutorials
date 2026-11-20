#!/usr/bin/env bash
# build_qe_cmake.sh - example CMake build recipe for Quantum ESPRESSO 7.6
#
# READ AND EDIT BEFORE USE. This script was written as documentation and has
# never been executed by the repository authors. Module names, compilers and
# library paths are site-specific.
#
# Usage:
#   MODE=cpu bash build_qe_cmake.sh          # CPU build (MPI + OpenMP + ScaLAPACK)
#   MODE=gpu GPU_ARCH=sm_90 bash build_qe_cmake.sh   # NVIDIA GPU build (NVHPC SDK)
#
# Notes for QE 7.6:
#   * CMake >= 3.20 is required.
#   * QE_ENABLE_CUDA / QE_ENABLE_OPENACC / QE_ENABLE_ROCM / QE_ENABLE_ONEAPI /
#     QE_ENABLE_OFFLOAD were removed. Use QE_GPU and QE_GPU_ARCHS instead.
#     Any of the old options, even set to OFF, stops CMake with FATAL_ERROR,
#     so always configure in a fresh build directory (or delete CMakeCache.txt).

set -euo pipefail

QE_VERSION=${QE_VERSION:-7.6}
SRC_DIR=${SRC_DIR:-$HOME/src/qe-${QE_VERSION}}
PREFIX=${PREFIX:-$HOME/opt/qe-${QE_VERSION}-${MODE:-cpu}}
MODE=${MODE:-cpu}
GPU_ARCH=${GPU_ARCH:-sm_80}       # A100 = sm_80, H100 = sm_90
NJOBS=${NJOBS:-8}
TARGETS=${TARGETS:-"pw ph pp"}    # e.g. "pwall", "pw ph pp hp neb", "epw"

export LC_ALL=C

# ---------------------------------------------------------------------------
# 1) Get the source (skipped if the directory already exists)
# ---------------------------------------------------------------------------
if [ ! -d "$SRC_DIR" ]; then
    git clone --branch "qe-${QE_VERSION}" --depth 1 \
        https://gitlab.com/QEF/q-e.git "$SRC_DIR"
fi
cd "$SRC_DIR"

# ---------------------------------------------------------------------------
# 2) Load the toolchain (EDIT: module names differ between clusters)
# ---------------------------------------------------------------------------
if command -v module >/dev/null 2>&1; then
    if [ "$MODE" = "gpu" ]; then
        : # module load nvhpc            # NVIDIA HPC SDK >= 21.7
    else
        : # module load gcc openmpi openblas fftw scalapack
        : # or: module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
    fi
fi

cmake --version | head -n 1

# ---------------------------------------------------------------------------
# 3) Configure in a fresh build directory
# ---------------------------------------------------------------------------
BUILD_DIR="build-${MODE}"
rm -rf "$BUILD_DIR"

if [ "$MODE" = "gpu" ]; then
    cmake -B "$BUILD_DIR" \
        -DCMAKE_Fortran_COMPILER=mpif90 \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DQE_GPU="openacc;cuda" \
        -DQE_GPU_ARCHS="$GPU_ARCH" \
        -DQE_ENABLE_MPI_GPU_AWARE=ON \
        -DQE_ENABLE_OPENMP=ON
else
    cmake -B "$BUILD_DIR" \
        -DCMAKE_Fortran_COMPILER=mpif90 \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DQE_ENABLE_OPENMP=ON \
        -DQE_ENABLE_SCALAPACK=ON \
        -DQE_FFTW_VENDOR=AUTO
fi

# ---------------------------------------------------------------------------
# 4) Build and install
# ---------------------------------------------------------------------------
# shellcheck disable=SC2086
cmake --build "$BUILD_DIR" -j "$NJOBS" --target $TARGETS
cmake --install "$BUILD_DIR"

echo "Executables:"
ls "$BUILD_DIR/bin"

# ---------------------------------------------------------------------------
# 5) Optional: run the pw.x system tests (labels look like system--pw)
# ---------------------------------------------------------------------------
if [ "${RUN_TESTS:-0}" = "1" ]; then
    (cd "$BUILD_DIR" && ctest -L system--pw -j 4 --output-on-failure)
fi

echo "Done. Add $PREFIX/bin to PATH."
