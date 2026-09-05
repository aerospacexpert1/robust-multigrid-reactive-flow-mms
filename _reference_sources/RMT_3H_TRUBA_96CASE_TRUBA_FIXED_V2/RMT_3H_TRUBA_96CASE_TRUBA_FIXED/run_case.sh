#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TASK_ID="${1:?usage: ./run_case.sh TASK_ID}"

MANIFEST="${ROOT}/campaign_manifest.csv"
BIN="${ROOT}/build/opposedflow_rmt3h"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: ${MANIFEST} does not exist."
    echo "Run:"
    echo "  ./rebuild_manifest.sh"
    exit 1
fi

if [[ ! -x "${BIN}" ]]; then
    echo "ERROR: binary not found:"
    echo "  ${BIN}"
    echo
    echo "Run ./build.sh first."
    exit 1
fi

if ! [[ "${TASK_ID}" =~ ^[0-9]+$ ]]; then
    echo "TASK_ID must be an integer."
    exit 1
fi

if (( TASK_ID < 0 || TASK_ID > 95 )); then
    echo "TASK_ID must be 0..95"
    exit 1
fi

ROW="$(sed -n "$((TASK_ID + 2))p" "${MANIFEST}")"

if [[ -z "${ROW}" ]]; then
    echo "No manifest row for task ${TASK_ID}"
    exit 1
fi

IFS=',' read -r \
    task_id \
    case_id \
    family \
    level \
    mesh \
    Nx \
    Ny \
    threads \
    vIn \
    arrA \
    arrBeta \
    arrTa \
    endTime \
    <<< "${ROW}"

if [[ "${task_id}" != "${TASK_ID}" ]]; then
    echo "Manifest task mismatch."
    echo "Requested : ${TASK_ID}"
    echo "Manifest  : ${task_id}"
    exit 1
fi

OUTDIR="${ROOT}/results/${family}/${level}/${mesh}/T${threads}"

mkdir -p "${OUTDIR}"

rm -f \
    "${OUTDIR}/RUN_COMPLETE" \
    "${OUTDIR}/RUN_FAILED"

cat > "${OUTDIR}/metadata.txt" <<META
task_id=${task_id}
case_id=${case_id}
solver=RMT_3H_9SHIFT_PRESSURE
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

# The Orfoz Slurm allocation is deliberately larger (56 CPUs) than the
# benchmark thread count.  The manifest value is the authoritative OpenMP size.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    allocated="${SLURM_CPUS_PER_TASK:-0}"
    if ! [[ "${allocated}" =~ ^[0-9]+$ ]] || (( allocated < threads )); then
        echo "ERROR: Slurm allocation (${allocated}) is smaller than benchmark threads (${threads})." >&2
        exit 20
    fi
fi

export OMP_NUM_THREADS="${threads}"
export OMP_THREAD_LIMIT="${threads}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "======================================================================"
echo "RMT case"
echo "======================================================================"
echo "Task      : ${task_id}"
echo "Case      : ${case_id}"
echo "Family    : ${family}"
echo "Level     : ${level}"
echo "Mesh      : ${mesh}"
echo "Nx Ny     : ${Nx} ${Ny}"
echo "Threads   : ${threads} (actual OpenMP benchmark threads)"
echo "Allocated : ${SLURM_CPUS_PER_TASK:-manual} CPUs/task (Slurm allocation)"
echo "vIn       : ${vIn}"
echo "Arr A     : ${arrA}"
echo "End time  : ${endTime}"
echo "OUTDIR    : ${OUTDIR}"
echo "======================================================================"

cd "${OUTDIR}"

set +e

"${BIN}" \
    -Nx "${Nx}" \
    -Ny "${Ny}" \
    -end "${endTime}" \
    -dt 1.0e-7 \
    -maxCo 0.05 \
    -vIn "${vIn}" \
    -arrA "${arrA}" \
    -arrBeta "${arrBeta}" \
    -arrTa "${arrTa}" \
    -pCycles 4 \
    -sCycles 4 \
    -rmtLevels 4 \
    -rmtPost 3 \
    -rmtCoarseSweeps 48 \
    -rmtCoarseRepeats 1 \
    -rmtOmega 0.005 \
    -fieldOutput 0 \
    2>&1 | tee solver.log

RC=${PIPESTATUS[0]}

set -e

echo "date_end=$(date --iso-8601=seconds)" >> metadata.txt
echo "exit_code=${RC}" >> metadata.txt

if (( RC != 0 )); then
    echo "Run failed with exit code ${RC}"
    touch RUN_FAILED
    exit "${RC}"
fi

if ! grep -q "RMT_RUN_SUMMARY" solver.log; then
    echo "ERROR: solver exited but RMT_RUN_SUMMARY was not found."
    touch RUN_FAILED
    exit 50
fi

python3 "${ROOT}/tools/extract_one_summary.py" \
    "${OUTDIR}" \
    "${task_id}" \
    "${case_id}" \
    "${family}" \
    "${level}" \
    "${mesh}" \
    "${Nx}" \
    "${Ny}" \
    "${threads}" \
    "${vIn}" \
    "${arrA}" \
    "${arrBeta}" \
    "${arrTa}"

touch RUN_COMPLETE

echo
echo "======================================================================"
echo "RUN COMPLETE"
echo "${case_id}"
echo "======================================================================"
