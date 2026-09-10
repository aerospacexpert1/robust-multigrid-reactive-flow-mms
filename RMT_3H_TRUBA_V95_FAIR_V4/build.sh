#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

echo "======================================================================"
echo "RMT V95 FAIR V4 build"
echo "======================================================================"
echo "Compiler: $(gcc --version | head -n1)"

rm -rf build
mkdir -p build

echo "[1/6] Generate RMT source from frozen V95 SG reference"
python3 tools/generate_rmt_v95.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_rmt3h_v95_fair.c

echo "[2/6] Finalize fail-fast convergence and restore byte-identical timestep labels"
python3 tools/finalize_generated_v95.py build/opposedflow_rmt3h_v95_fair.c

echo "[3/6] Audit that non-pressure V95 kernels/timestep are unchanged"
python3 tools/audit_fairness.py \
  src/opposedflow_v95_frozen_sg_reference.c \
  build/opposedflow_rmt3h_v95_fair.c

echo "[4/6] Compile pressure RMT self-test"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
  tests/rmt_pressure_selftest.c -lm -o build/rmt_pressure_selftest

echo "[5/6] Run pressure RMT self-test"
./build/rmt_pressure_selftest | tee build/rmt_pressure_selftest.log
grep -q "SELFTEST PASS" build/rmt_pressure_selftest.log

echo "[6/6] Compile frozen-V95 + RMT reacting-flow solver"
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp -Isrc \
  build/opposedflow_rmt3h_v95_fair.c -lm -o build/opposedflow_rmt3h_v95_fair

sha256sum src/opposedflow_v95_frozen_sg_reference.c src/rmt_book2d_impl.h \
  tools/generate_rmt_v95.py tools/finalize_generated_v95.py \
  build/opposedflow_rmt3h_v95_fair.c build/opposedflow_rmt3h_v95_fair \
  > build/SHA256SUMS.txt

echo "RMT V95 FAIR V4 BUILD PASS"
