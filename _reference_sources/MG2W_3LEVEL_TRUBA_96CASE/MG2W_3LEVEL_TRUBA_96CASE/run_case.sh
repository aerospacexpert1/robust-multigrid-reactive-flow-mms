#!/bin/bash
set -eo pipefail

ROOT=${ROOT:-$PWD}
INDEX=${1:?usage: run_case.sh ARRAY_INDEX}
MANIFEST="$ROOT/campaign_manifest.csv"
EXE="$ROOT/build/opposedflow_mg2w"

[[ -x "$EXE" ]] || { echo "ERROR: executable not found: $EXE" >&2; exit 2; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 2; }

line=$(sed -n "$((INDEX + 2))p" "$MANIFEST")
[[ -n "$line" ]] || { echo "ERROR: manifest has no row for index $INDEX" >&2; exit 3; }

IFS=',' read -r array_index case_id study_family group_dir physical_case mesh_name Nx Ny velocity stiffness_class A beta Ta threads core_label <<< "$line"
for _v in array_index case_id study_family group_dir physical_case mesh_name Nx Ny velocity stiffness_class A beta Ta threads core_label; do
    printf -v "$_v" '%s' "${!_v//$'\r'/}"
done
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

# Matched physical/numerical settings used by the frozen V95/OpenFOAM campaign.
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

# Classical three-level W-cycle controls (gamma=2).
MG_PRE=${MG_PRE:-3}
MG_POST=${MG_POST:-3}
MG_COARSE_SWEEPS=${MG_COARSE_SWEEPS:-40}
MG_OMEGA=${MG_OMEGA:-1.0}
P_CYCLES=${P_CYCLES:-500}
P_REL_TOL=${P_REL_TOL:-1e-4}
P_ABS_TOL=${P_ABS_TOL:-1e-6}

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
solver=MG2W_3LEVEL
level0=${Nx}x${Ny}
level1=$((Nx/2))x$((Ny/2))
level2=$((Nx/4))x$((Ny/4))
mg_pre_sweeps=$MG_PRE
mg_post_sweeps=$MG_POST
mg_coarse_sweeps=$MG_COARSE_SWEEPS
mg_omega=$MG_OMEGA
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

start_ns=$(date +%s%N)
set +e
srun --nodes=1 --ntasks=1 --cpus-per-task="$threads" --exclusive --cpu-bind=cores \
    "$EXE" \
    -Nx "$Nx" -Ny "$Ny" \
    -end "$END_TIME" -dt "$DT0" -dtMax "$DTMAX" \
    -vFuel "$velocity" -vOx "-$velocity" \
    -threads "$threads" -caseId "$case_id" \
    -mgPre "$MG_PRE" -mgPost "$MG_POST" \
    -mgCoarseSweeps "$MG_COARSE_SWEEPS" -mgOmega "$MG_OMEGA" \
    -pCycles "$P_CYCLES" -pRelTol "$P_REL_TOL" -pAbsTol "$P_ABS_TOL" \
    -writeOutput 1 -write "$WRITE_INTERVAL" -outputDir "$out_dir" \
    -progressEvery "$PROGRESS_EVERY" \
    -A "$A" -beta "$beta" -Ta "$Ta" \
    -perfectGas 0 -variableCp 0 -sutherland 0 -rhoRelax 1 \
    -rho "$RHO0" -mu "$MU0" -cp "$CP0" -Pr "$PR0" -Sc "$SC0" -HfF "$HFF" \
    > "$run_dir/run.log" 2>&1
rc=$?
set -e
end_ns=$(date +%s%N)
launcher_wall=$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN {printf "%.9f", (e-s)/1.0e9}')

cat > "$run_dir/launcher_timing.txt" <<TIMING
launcher_wall_seconds_including_srun_and_output=$launcher_wall
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

if [[ $valid -eq 1 ]]; then
    touch "$run_dir/RUN_COMPLETE"
    echo "COMPLETED: $case_id launcher_wall=${launcher_wall}s"
else
    cat > "$run_dir/RUN_FAILED" <<FAIL
case_id=$case_id
exit_code=$rc
finish=$(date --iso-8601=seconds)
FAIL
    echo "FAILED: $case_id exit=$rc" >&2
    exit ${rc:-1}
fi
