#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASK_ID="${1:?usage: ./run_case.sh TASK_ID}"
MANIFEST="${ROOT}/campaign_manifest.csv"
BIN="${ROOT}/build/opposedflow_rmt3h_bookfix_v3"

[[ -f "${MANIFEST}" ]] || { echo "ERROR: missing ${MANIFEST}; run ./rebuild_manifest.sh"; exit 1; }
[[ -x "${BIN}" ]] || { echo "ERROR: missing binary ${BIN}; run ./build.sh"; exit 1; }
[[ "${TASK_ID}" =~ ^[0-9]+$ ]] || { echo "TASK_ID must be integer"; exit 1; }
(( TASK_ID >= 0 && TASK_ID <= 95 )) || { echo "TASK_ID must be 0..95"; exit 1; }
ROW="$(sed -n "$((TASK_ID+2))p" "${MANIFEST}")"
[[ -n "${ROW}" ]] || { echo "No manifest row for task ${TASK_ID}"; exit 1; }
IFS=',' read -r task_id case_id family level mesh Nx Ny threads vIn arrA arrBeta arrTa endTime <<< "${ROW}"
[[ "${task_id}" == "${TASK_ID}" ]] || { echo "Manifest mismatch"; exit 1; }

OUTDIR="${ROOT}/results/${family}/${level}/${mesh}/T${threads}"
mkdir -p "${OUTDIR}"
rm -f "${OUTDIR}/RUN_COMPLETE" "${OUTDIR}/RUN_FAILED"
cat > "${OUTDIR}/metadata.txt" <<META
rmt_version=BOOKFIX_V3
legacy_reference=LEGACY_V2
comparison_policy=same_physics_mesh_endTime_dt0_maxCo_pressure_tolerance
algorithm_change=Martynenko_style_index_mapping_shifted_boundary_independent_xy_levels_full_strength_correction
fallback_policy=omega_1_then_halving_only_if_residual_does_not_decrease;fallback_count_reported
task_id=${task_id}
case_id=${case_id}
family=${family}
level=${level}
mesh=${mesh}
Nx=${Nx}
Ny=${Ny}
threads=${threads}
vIn=${vIn}
arrA=${arrA}
arrBeta=${arrBeta}
arrTa=${arrTa}
endTime=${endTime}
host=$(hostname)
date_start=$(date --iso-8601=seconds)
slurm_job_id=${SLURM_JOB_ID:-none}
slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-none}
slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-none}
META

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    allocated="${SLURM_CPUS_PER_TASK:-0}"
    if ! [[ "${allocated}" =~ ^[0-9]+$ ]] || (( allocated < threads )); then
        echo "ERROR: Slurm allocation (${allocated}) < benchmark threads (${threads})" >&2
        exit 20
    fi
fi
export OMP_NUM_THREADS="${threads}"
export OMP_THREAD_LIMIT="${threads}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "======================================================================"
echo "RMT BOOKFIX V3 case"
echo "======================================================================"
echo "Task      : ${task_id}"
echo "Case      : ${case_id}"
echo "Family    : ${family}"
echo "Level     : ${level}"
echo "Mesh      : ${mesh}"
echo "Nx Ny     : ${Nx} ${Ny}"
echo "Threads   : ${threads}"
echo "Allocated : ${SLURM_CPUS_PER_TASK:-manual} CPUs/task"
echo "vIn       : ${vIn}"
echo "Arr A/b/Ta: ${arrA} ${arrBeta} ${arrTa}"
echo "End time  : ${endTime}"
echo "RMT       : auto independent x/y levels; 8 smooth; 16 coarsest; omega start 1.0"
echo "Fairness  : dt0=1e-7 maxCo=0.05 pCycles(max)=4 sCycles(max)=4, unchanged physics"
echo "======================================================================"

cd "${OUTDIR}"
set +e
if command -v stdbuf >/dev/null 2>&1; then
    stdbuf -oL -eL "${BIN}" \
        -Nx "${Nx}" -Ny "${Ny}" -end "${endTime}" \
        -dt 1.0e-7 -maxCo 0.05 \
        -vIn "${vIn}" -arrA "${arrA}" -arrBeta "${arrBeta}" -arrTa "${arrTa}" \
        -pCycles 4 -sCycles 4 \
        -rmtLevels 0 -rmtPost 8 -rmtCoarseSweeps 16 -rmtCoarseRepeats 1 -rmtOmega 1.0 \
        -progressEvery 100 -fieldOutput 0 2>&1 | tee solver.log
else
    "${BIN}" \
        -Nx "${Nx}" -Ny "${Ny}" -end "${endTime}" \
        -dt 1.0e-7 -maxCo 0.05 \
        -vIn "${vIn}" -arrA "${arrA}" -arrBeta "${arrBeta}" -arrTa "${arrTa}" \
        -pCycles 4 -sCycles 4 \
        -rmtLevels 0 -rmtPost 8 -rmtCoarseSweeps 16 -rmtCoarseRepeats 1 -rmtOmega 1.0 \
        -progressEvery 100 -fieldOutput 0 2>&1 | tee solver.log
fi
RC=${PIPESTATUS[0]}
set -e

echo "date_end=$(date --iso-8601=seconds)" >> metadata.txt
echo "exit_code=${RC}" >> metadata.txt
if (( RC != 0 )); then touch RUN_FAILED; exit "${RC}"; fi
if ! grep -q "RMT_RUN_SUMMARY" solver.log; then
    echo "ERROR: RMT_RUN_SUMMARY missing"; touch RUN_FAILED; exit 50
fi
if ! grep -q "RMT_BOOK_DIAGNOSTICS" solver.log; then
    echo "ERROR: RMT_BOOK_DIAGNOSTICS missing"; touch RUN_FAILED; exit 51
fi
python3 "${ROOT}/tools/extract_one_summary.py" "${OUTDIR}" "${ROW}"
touch RUN_COMPLETE
echo "RUN COMPLETE: ${case_id}"
