#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

echo "======================================================================"
echo "FINALVOL2 build (optimized implementation, calibration not applied)"
echo "======================================================================"
echo "Compiler: $(gcc --version | head -n1)"

rm -rf build
mkdir -p build

echo "[1/7] Generate solver from frozen V95 reference"
python3 tools/generate_finalvol2.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_finalvol2.c

echo "[2/7] Audit frozen V95 physics/timestep"
python3 tools/audit_fairness.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_finalvol2.c

echo "[3/7] Audit calibration/production separation"
python3 tools/audit_calibration_separation.py

echo "[4/9] Compile optimized RMT self-test"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
  tests/rmt_pressure_selftest.c -lm -o build/rmt_pressure_selftest

echo "[5/9] Validate robust baseline setting (16 sweeps)"
RMT_TEST_SWEEPS=16 ./build/rmt_pressure_selftest | tee build/rmt_pressure_selftest.log
grep -q "BOUNDARY_CV_RESTRICTION PASS" build/rmt_pressure_selftest.log
grep -q "SELFTEST PASS" build/rmt_pressure_selftest.log

echo "[6/9] Compile optimization-equivalence test against RMT FINAL"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
  tests/optimization_equivalence.c -lm -o build/optimization_equivalence

echo "[7/9] Prove implementation optimizations preserve the RMT correction"
OMP_NUM_THREADS=4 OMP_DYNAMIC=false ./build/optimization_equivalence | tee build/optimization_equivalence.log
grep -q "OPTIMIZATION_EQUIVALENCE PASS" build/optimization_equivalence.log

echo "[8/9] Compile independent calibration executable"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
  calibration/rmt_smoothing_calibration.c -lm -o build/rmt_smoothing_calibration

echo "[9/9] Compile frozen-V95 + FINALVOL2 reacting-flow solver"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp -Isrc \
  build/opposedflow_finalvol2.c -lm -o build/opposedflow_finalvol2

sha256sum \
  src/opposedflow_v95_frozen_sg_reference.c \
  src/rmt_finalvol2_impl.h \
  tools/generate_finalvol2.py \
  tools/audit_fairness.py \
  tools/audit_calibration_separation.py \
  tests/rmt_pressure_selftest.c \
  tests/reference_rmt_final2d_impl.h \
  tests/optimization_equivalence.c \
  calibration/rmt_smoothing_calibration.c \
  build/opposedflow_finalvol2.c build/opposedflow_finalvol2 \
  > build/SHA256SUMS.txt

echo "FINALVOL2 BUILD PASS"
