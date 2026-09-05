#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SRC="${ROOT}/src/opposedflow_rmt3h.c"
BIN="${ROOT}/build/opposedflow_rmt3h"

mkdir -p "${ROOT}/build"

echo "Compiler:"
gcc --version | head -1

echo
echo "Building:"
echo "${SRC}"

gcc \
    -O2 \
    -std=c11 \
    -fopenmp \
    -Wall \
    -Wextra \
    "${SRC}" \
    -lm \
    -o "${BIN}"

echo
echo "Built:"
ls -lh "${BIN}"

echo
echo "OpenMP linkage:"
ldd "${BIN}" | grep -E 'gomp|omp' || true

echo
echo "Build complete."
