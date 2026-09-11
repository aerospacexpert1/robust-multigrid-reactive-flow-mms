#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
[[ -x build/opposedflow_rmt_final ]] || ./build.sh

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
RMT_SWEEPS=${RMT_SWEEPS:-16}

rm -rf first_pressure_logs
mkdir -p first_pressure_logs

run_one() {
  local tag="$1" nx="$2" ny="$3" vel="$4" A="$5" Ta="$6"
  local log="first_pressure_logs/${tag}.log"
  echo "===== FIRST PRESSURE: $tag Nx=$nx Ny=$ny v=$vel A=$A Ta=$Ta sweeps=$RMT_SWEEPS ====="
  ./build/opposedflow_rmt_final \
    -Nx "$nx" -Ny "$ny" \
    -end 1e-7 -dt 1e-7 -dtMax 1e-4 \
    -vFuel "$vel" -vOx "-$vel" \
    -threads "$OMP_NUM_THREADS" -caseId "$tag" \
    -rmtLevels 0 -rmtSweeps "$RMT_SWEEPS" -rmtOmega 1.0 \
    -pCycles 500 -pRelTol 1e-4 -pAbsTol 1e-6 \
    -sCycles 4 -sSweeps 5 \
    -writeOutput 0 -progressEvery 1 \
    -A "$A" -beta 5 -Ta "$Ta" \
    -perfectGas 0 -variableCp 0 -sutherland 0 -rhoRelax 1 \
    -rho 1 -mu 2e-5 -cp 1000 -Pr 0.7 -Sc 1 -HfF 3.571428571e6 \
    >"$log" 2>&1

  grep -q '^Done\. Final time=1e-07' "$log"
  grep -q 'RMT_FINAL_DIAGNOSTICS' "$log"
  ! grep -q 'RMT_PRESSURE_NOT_CONVERGED' "$log"
  grep -E 'step=|Final pressure:|Pressure totals:|RMT_FINAL_DIAGNOSTICS' "$log" | tail -8
}

run_one M1_LOW       108  36  0.05 1e7 1000
run_one M4_HIGH_FLOW 864 288  0.20 1e7 1000
run_one M4_STIFF     864 288  0.10 1e9 100

echo "FIRST_PRESSURE_MATRIX PASS"
