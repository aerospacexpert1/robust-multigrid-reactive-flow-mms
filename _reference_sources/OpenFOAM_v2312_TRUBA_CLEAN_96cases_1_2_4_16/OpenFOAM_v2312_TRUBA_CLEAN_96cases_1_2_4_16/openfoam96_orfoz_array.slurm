#!/bin/bash
#SBATCH --job-name=OF2312_96
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=3-00:00:00
#SBATCH --array=0-23%3
#SBATCH --output=slurm-OF2312_96_%A_%a.out
#SBATCH --error=slurm-OF2312_96_%A_%a.err

# Do not use `set -u` around the TRUBA OpenFOAM module environment.
set -o pipefail

# TRUBA-native OpenFOAM environment verified on an Orfoz compute node.
set +e
module purge
module load apps/openfoam/v2312
MODULE_RC=$?
set -e

echo "module load return code: ${MODULE_RC}"

# Validate the environment by variables/executables.  The foamVersion helper is
# not used here because it is not reliably exported into non-interactive jobs.
if [[ "${WM_PROJECT_VERSION:-}" != "v2312" ]]; then
    echo "ERROR: expected WM_PROJECT_VERSION=v2312, got '${WM_PROJECT_VERSION:-unset}'" >&2
    exit 10
fi

required=(blockMesh checkMesh decomposePar reconstructPar rhoReactingFoam mpirun)
for exe in "${required[@]}"; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "ERROR: required executable not found after module load: $exe" >&2
        exit 11
    }
done

# Intel MPI settings for the TRUBA v2312 module.
export I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export OMP_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR}"
RUNNER="${ROOT}/run_case_orfoz.sh"

if [[ ! -d "$ROOT/cases" ]]; then
    echo "ERROR: cases/ not found under SLURM_SUBMIT_DIR=$ROOT" >&2
    echo "Submit this script from the campaign root directory." >&2
    exit 12
fi
if [[ ! -f "$RUNNER" ]]; then
    echo "ERROR: run_case_orfoz.sh not found under $ROOT" >&2
    exit 13
fi
if [[ ! -f "$ROOT/physical_configurations.txt" ]]; then
    echo "ERROR: physical_configurations.txt not found under $ROOT" >&2
    exit 14
fi

cd "$ROOT"

LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$ROOT/physical_configurations.txt")"
if [[ -z "$LINE" ]]; then
    echo "ERROR: no physical configuration for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 15
fi

read -r GROUP PHYSICAL_CASE MESH <<< "$LINE"
CORES=(C1 C2 C4 C16)

echo "============================================================"
echo "OPENFOAM V2312 ORFOZ CAMPAIGN"
echo "============================================================"
echo "Date:               $(date --iso-8601=seconds)"
echo "Node:               $(hostname)"
echo "Submit directory:   $ROOT"
echo "Job ID:             ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}"
echo "Array task:         ${SLURM_ARRAY_TASK_ID}"
echo "Configuration:      ${GROUP}/${PHYSICAL_CASE}/${MESH}"
echo "Reserved tasks:     ${SLURM_NTASKS}"
echo "OpenFOAM version:   ${WM_PROJECT_VERSION}"
echo "WM_PROJECT_DIR:     ${WM_PROJECT_DIR:-unknown}"
echo "WM_MPLIB:           ${WM_MPLIB:-unknown}"
echo "MPI:                $(command -v mpirun)"
echo "MPI sizes:          1, 2, 4, 16 (sequential)"
echo "Excluded node:      orfoz359"
echo "============================================================"

for C in "${CORES[@]}"; do
    CASE_DIR="$ROOT/cases/$GROUP/$PHYSICAL_CASE/$MESH/$C"

    if [[ ! -d "$CASE_DIR" ]]; then
        echo "ERROR: case directory is missing: $CASE_DIR" >&2
        exit 16
    fi

    if [[ -f "$CASE_DIR/RUN_COMPLETE" ]]; then
        echo "SKIP completed: $GROUP/$PHYSICAL_CASE/$MESH/$C"
        continue
    fi

    echo
    echo "============================================================"
    echo "RUN $GROUP/$PHYSICAL_CASE/$MESH/$C"
    echo "START $(date --iso-8601=seconds)"
    echo "============================================================"

    bash "$RUNNER" "$CASE_DIR"
done

echo
echo "ARRAY TASK COMPLETE: $GROUP/$PHYSICAL_CASE/$MESH"
echo "END $(date --iso-8601=seconds)"
