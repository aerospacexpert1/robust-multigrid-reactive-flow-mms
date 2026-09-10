#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
[[ -x build/opposedflow_rmt3h_v95_fair ]] || ./build.sh
rm -rf smoke_output
mkdir -p smoke_output
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores

./build/opposedflow_rmt3h_v95_fair \
  -Nx 108 -Ny 36 \
  -end 0.02 -dt 1e-7 -dtMax 1e-4 \
  -vFuel 0.05 -vOx -0.05 \
  -threads "$OMP_NUM_THREADS" -caseId RMT_V95_FAIR_SMOKE \
  -rmtLevels 0 -rmtPost 8 -rmtCoarseSweeps 16 -rmtOmega 1.0 \
  -pCycles 500 -pRelTol 1e-4 -pAbsTol 1e-6 \
  -sCycles 4 -sSweeps 5 \
  -writeOutput 1 -write 0.02 -outputDir smoke_output \
  -progressEvery 100 \
  -A 1e7 -beta 5 -Ta 1000 \
  -perfectGas 0 -variableCp 0 -sutherland 0 -rhoRelax 1 \
  -rho 1.0 -mu 2e-5 -cp 1000 -Pr 0.7 -Sc 1.0 -HfF 3.571428571e6 \
  2>&1 | tee smoke_output/smoke.log

! grep -q "RMT_PRESSURE_REJECT" smoke_output/smoke.log
grep -q '^Done\. Final time=0.02' smoke_output/smoke.log
grep -q 'RMT_BOOK_DIAGNOSTICS' smoke_output/smoke.log
[[ -s smoke_output/summary.csv ]]

python3 - <<'PY'
import csv, math
p='smoke_output/summary.csv'
with open(p,newline='') as f:
    row=next(csv.DictReader(f))

def num(k): return float(row[k])
assert abs(num('final_time_s')-0.02) < 1e-12, row['final_time_s']
assert int(float(row['time_steps'])) < 10000, row['time_steps']
assert math.isfinite(num('mass_residual_max'))
assert math.isfinite(num('Tmax_K')) and 200.0 <= num('Tmax_K') <= 5000.0
for k in ('YFmax','YOmax','YPmax','YN2max'):
    assert math.isfinite(num(k)) and -1e-12 <= num(k) <= 1.000001, (k,row[k])
for k in ('umin_m_per_s','umax_m_per_s','vmin_m_per_s','vmax_m_per_s'):
    assert math.isfinite(num(k)) and abs(num(k)) < 100.0, (k,row[k])
assert num('final_pressure_rel_residual') <= 1.01e-4 or num('final_pressure_abs_residual') <= 1.01e-6
print('RMT_V95_FAIR_SMOKE PASS')
print('steps=',row['time_steps'],'pressure_time_s=',row['pressure_time_s'],'Tmax=',row['Tmax_K'],'massMax=',row['mass_residual_max'])
PY
