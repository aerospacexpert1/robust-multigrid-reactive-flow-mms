#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

echo "======================================================================"
echo "RMT FINAL build"
echo "======================================================================"
echo "Compiler: $(gcc --version | head -n1)"

rm -rf build
mkdir -p build

echo "[1/5] Generate pressure-solver replacement from frozen V95 source"
python3 tools/generate_rmt_final.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_rmt_final.c

echo "[2/5] Audit frozen V95 physics and physical timestep"
python3 tools/audit_fairness.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_rmt_final.c

echo "[3/5] Compile RMT mesh-ladder self-test"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
  tests/rmt_pressure_selftest.c -lm -o build/rmt_pressure_selftest

echo "[4/5] Run boundary-CV + 108/216/432/864 mesh-ladder self-test"
RMT_TEST_SWEEPS=${RMT_TEST_SWEEPS:-4} \
  ./build/rmt_pressure_selftest | tee build/rmt_pressure_selftest.log
grep -q "BOUNDARY_CV_RESTRICTION PASS" build/rmt_pressure_selftest.log
grep -q "SELFTEST PASS" build/rmt_pressure_selftest.log

echo "[5/5] Compile frozen-V95 + RMT FINAL reacting-flow solver"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp -Isrc \
  build/opposedflow_rmt_final.c -lm -o build/opposedflow_rmt_final

sha256sum \
  src/opposedflow_v95_frozen_sg_reference.c \
  src/rmt_final2d_impl.h \
  tools/generate_rmt_final.py \
  tools/audit_fairness.py \
  tests/rmt_pressure_selftest.c \
  build/opposedflow_rmt_final.c \
  build/opposedflow_rmt_final \
  > build/SHA256SUMS.txt

echo "RMT FINAL BUILD PASS"
