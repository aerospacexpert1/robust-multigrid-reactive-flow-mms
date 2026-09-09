#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LEGACY="${ROOT}/src/opposedflow_rmt3h_LEGACY_V2.c"
GENERATED="${ROOT}/build/opposedflow_rmt3h_BOOKFIX_V3.c"
BIN="${ROOT}/build/opposedflow_rmt3h_bookfix_v3"
SELFTEST="${ROOT}/build/rmt_pressure_selftest"
mkdir -p "${ROOT}/build"

echo "======================================================================"
echo "RMT BOOKFIX V3 build"
echo "======================================================================"
echo "Compiler: $(gcc --version | head -1)"

echo "[1/4] Generate corrected source from frozen legacy source"
python3 "${ROOT}/tools/patch_bookfix_v3.py" "${LEGACY}" "${GENERATED}"

echo "[2/4] Compile Martynenko-style pressure self-test"
gcc -O2 -std=c11 -fopenmp -Wall -Wextra -I"${ROOT}/src" \
    "${ROOT}/tests/rmt_pressure_selftest.c" -lm -o "${SELFTEST}"

echo "[3/4] Run pressure self-test"
OMP_NUM_THREADS=1 "${SELFTEST}" | tee "${ROOT}/build/rmt_pressure_selftest.log"
grep -q "SELFTEST PASS" "${ROOT}/build/rmt_pressure_selftest.log"

echo "[4/4] Compile corrected reacting-flow solver"
gcc -O2 -std=c11 -fopenmp -Wall -Wextra -I"${ROOT}/src" \
    "${GENERATED}" -lm -o "${BIN}"

echo
sha256sum "${LEGACY}" "${GENERATED}" "${BIN}" | tee "${ROOT}/build/SHA256SUMS.txt"
echo
ls -lh "${BIN}"
ldd "${BIN}" | grep -E 'gomp|omp' || true

echo "RMT BOOKFIX V3 BUILD PASS"
