#!/bin/bash
set -euo pipefail
ROOT=${ROOT:-$PWD}
INDEX=${1:?usage: run_case.sh ARRAY_INDEX}
MANIFEST="$ROOT/campaign_manifest.csv"
EXE="$ROOT/build/opposedflow_rmt_final"
[[ -x "$EXE" ]] || { echo "ERROR: executable not found: $EXE" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 2; }

line=$(sed -n "$((INDEX + 2))p" "$MANIFEST")
[[ -n "$line" ]] || { echo "ERROR: manifest has no row for index $INDEX" >&2; exit 3; }
IFS=',' read -r array_index case_id study_family group_dir physical_case mesh_name Nx Ny velocity stiffness_class A beta Ta threads core_label <<< "$line"
[[ "$array_index" == "$INDEX" ]] || { echo "ERROR: manifest/index mismatch" >&2; exit 4; }

run_dir="$ROOT/results/$group_dir/$physical_case/$mesh_name/$core_label"
out_dir="$run_dir/output"
mkdir -p "$run_dir"
if [[ -f "$run_dir/RUN_COMPLETE" ]]; then
  echo "SKIP completed: $case_id"
  exit 0
fi
rm -rf "$out_dir"
rm -f "$run_dir/RUN_FAILED" "$run_dir/RUN_COMPLETE"
mkdir -p "$out_dir"

# EXACT matched V95 campaign physics/numerics used by SG/MG2V/MG2W/MG3V.
END_TIME=${END_TIME:-2.0}
WRITE_INTERVAL=${WRITE_INTERVAL:-2.0}
DT0=${DT0:-1e-7}
DTMAX=${DTMAX:-1e-4}
PROGRESS_EVERY=${PROGRESS_EVERY:-1000}
RHO0=${RHO0:-1.0}
MU0=${MU0:-2.0e-5}
CP0=${CP0:-1000.0}
PR0=${PR0:-0.7}
SC0=${SC0:-1.0}
HFF=${HFF:-3.571428571e6}
P_CYCLES=${P_CYCLES:-500}
P_REL_TOL=${P_REL_TOL:-1e-4}
P_ABS_TOL=${P_ABS_TOL:-1e-6}

# RMT FINAL has one RMT tuning knob: postsmoothing sweep count.
RMT_LEVELS=${RMT_LEVELS:-0}
RMT_SMOOTH_SWEEPS=${RMT_SMOOTH_SWEEPS:-4}
RMT_OMEGA=${RMT_OMEGA:-1.0}

cat > "$run_dir/parameters.txt" <<PARAMS
array_index=$array_index
case_id=$case_id
study_family=$study_family
group_dir=$group_dir
physical_case=$physical_case
mesh_name=$mesh_name
Nx=$Nx
Ny=$Ny
cells=$((Nx*Ny))
velocity_m_per_s=$velocity
stiffness_class=$stiffness_class
A=$A
beta=$beta
Ta_K=$Ta
threads=$threads
core_label=$core_label
solver=RMT_FINAL
rmt_levels=$RMT_LEVELS
rmt_smoothing_sweeps=$RMT_SMOOTH_SWEEPS
rmt_coarsest_solver=dense_direct
rmt_correction_factor=$RMT_OMEGA
rmt_cycle_policy=full_correction_no_monotonic_line_search_no_hidden_fallback
pressure_max_cycles=$P_CYCLES
pressure_rel_tol=$P_REL_TOL
pressure_abs_tol=$P_ABS_TOL
end_time=$END_TIME
write_interval=$WRITE_INTERVAL
dt0=$DT0
dtMax=$DTMAX
rho=$RHO0
mu=$MU0
cp=$CP0
Pr=$PR0
Sc=$SC0
HfF=$HFF
variableCp=0
variableDensity=0
sutherland=0
job_id=${SLURM_JOB_ID:-manual}
array_job_id=${SLURM_ARRAY_JOB_ID:-manual}
array_task_id=${SLURM_ARRAY_TASK_ID:-manual}
node=$(hostname)
start=$(date --iso-8601=seconds)
PARAMS

export OMP_NUM_THREADS="$threads"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores

cmd=("$EXE"
  -Nx "$Nx" -Ny "$Ny"
  -end "$END_TIME" -dt "$DT0" -dtMax "$DTMAX"
  -vFuel "$velocity" -vOx "-$velocity"
  -threads "$threads" -caseId "$case_id"
  -rmtLevels "$RMT_LEVELS" -rmtSweeps "$RMT_SMOOTH_SWEEPS" -rmtOmega "$RMT_OMEGA"
  -pCycles "$P_CYCLES" -pRelTol "$P_REL_TOL" -pAbsTol "$P_ABS_TOL"
  -sCycles 4 -sSweeps 5
  -writeOutput 1 -write "$WRITE_INTERVAL" -outputDir "$out_dir"
  -progressEvery "$PROGRESS_EVERY"
  -A "$A" -beta "$beta" -Ta "$Ta"
  -perfectGas 0 -variableCp 0 -sutherland 0 -rhoRelax 1
  -rho "$RHO0" -mu "$MU0" -cp "$CP0" -Pr "$PR0" -Sc "$SC0" -HfF "$HFF")

start_ns=$(date +%s%N)
set +e
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
  srun --nodes=1 --ntasks=1 --cpus-per-task="$threads" --exclusive --cpu-bind=cores \
    "${cmd[@]}" > "$run_dir/run.log" 2>&1
  rc=$?
else
  "${cmd[@]}" > "$run_dir/run.log" 2>&1
  rc=$?
fi
set -e
end_ns=$(date +%s%N)
launcher_wall=$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN {printf "%.9f", (e-s)/1.0e9}')
cat > "$run_dir/launcher_timing.txt" <<TIMING
launcher_wall_seconds_including_launcher_and_output=$launcher_wall
exit_code=$rc
finish=$(date --iso-8601=seconds)
TIMING

valid=1
[[ $rc -eq 0 ]] || valid=0
[[ -s "$out_dir/summary.csv" ]] || valid=0
[[ -s "$out_dir/fields.pvd" ]] || valid=0
[[ -s "$out_dir/performance/pressure_step_history.csv" ]] || valid=0
[[ -s "$out_dir/performance/final_pressure_cycle_history.csv" ]] || valid=0
grep -q '^Done\.' "$run_dir/run.log" || valid=0
grep -q 'RMT_FINAL_DIAGNOSTICS' "$run_dir/run.log" || valid=0
if grep -q 'RMT_PRESSURE_NOT_CONVERGED' "$run_dir/run.log"; then valid=0; fi

if [[ -s "$out_dir/summary.csv" ]]; then
  if ! python3 - "$out_dir/summary.csv" "$END_TIME" "$P_REL_TOL" "$P_ABS_TOL" <<'PY'
import csv, math, sys
path,end,rtol,atol=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])
with open(path,newline='') as f:
    row=next(csv.DictReader(f))
def x(k): return float(row[k])
assert row['solver']=='RMT_FINAL', row['solver']
assert abs(x('final_time_s')-end) <= 1e-11*max(1.0,abs(end))
assert math.isfinite(x('mass_residual_max'))
assert math.isfinite(x('Tmax_K'))
assert x('final_pressure_rel_residual') <= 1.01*rtol or x('final_pressure_abs_residual') <= 1.01*atol
for k in ('YFmax','YOmax','YPmax','YN2max'):
    assert math.isfinite(x(k)) and -1e-10 <= x(k) <= 1.00001, (k,row[k])
PY
  then
    valid=0
  fi
fi

if [[ $valid -eq 1 ]]; then
  touch "$run_dir/RUN_COMPLETE"
  echo "COMPLETED: $case_id launcher_wall=${launcher_wall}s"
else
  cat > "$run_dir/RUN_FAILED" <<FAIL
case_id=$case_id
exit_code=$rc
finish=$(date --iso-8601=seconds)
FAIL
  echo "FAILED: $case_id exit=$rc; inspect $run_dir/run.log" >&2
  if [[ $rc -ne 0 ]]; then exit "$rc"; else exit 20; fi
fi
