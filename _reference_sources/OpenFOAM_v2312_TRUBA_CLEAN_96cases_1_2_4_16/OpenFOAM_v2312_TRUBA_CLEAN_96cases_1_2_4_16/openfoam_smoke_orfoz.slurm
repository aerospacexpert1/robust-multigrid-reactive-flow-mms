#!/bin/bash
#SBATCH --job-name=OF2312_smoke
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=slurm-OF2312_smoke_%j.out
#SBATCH --error=slurm-OF2312_smoke_%j.err

set -o pipefail
set +e
module purge
module load apps/openfoam/v2312
MODULE_RC=$?
set -e

echo "module load return code: ${MODULE_RC}"
[[ "${WM_PROJECT_VERSION:-}" == "v2312" ]] || { echo "ERROR: OpenFOAM v2312 environment missing" >&2; exit 10; }
for exe in blockMesh checkMesh decomposePar reconstructPar rhoReactingFoam mpirun; do
    command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: missing $exe" >&2; exit 11; }
done

export I_MPI_HYDRA_BOOTSTRAP=slurm
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export OMP_NUM_THREADS=1

ROOT="${SLURM_SUBMIT_DIR}"
TESTROOT="${ROOT}/_v2312_campaign_smoke_${SLURM_JOB_ID}"
cd "$ROOT"

[[ -d "$ROOT/cases" ]] || { echo "ERROR: submit from campaign root" >&2; exit 12; }
rm -rf "$TESTROOT"
mkdir -p "$TESTROOT"

cp -a "$ROOT/cases/FLOW/BC2_NOMINAL_v010/M1_108x36/C1" \
      "$TESTROOT/serial_nominal_M1_C1"
cp -a "$ROOT/cases/STIFFNESS/S3_HIGH_A1e9_Ta100/M1_108x36/C16" \
      "$TESTROOT/parallel16_high_stiffness_M1_C16"

for case_dir in \
    "$TESTROOT/serial_nominal_M1_C1" \
    "$TESTROOT/parallel16_high_stiffness_M1_C16"
do
    sed -Ei 's/^(startFrom[[:space:]]+)[^;]+;/\1startTime;/' "$case_dir/system/controlDict"
    sed -Ei 's/^(startTime[[:space:]]+)[^;]+;/\10;/' "$case_dir/system/controlDict"
    sed -Ei 's/^(endTime[[:space:]]+)[^;]+;/\11e-5;/' "$case_dir/system/controlDict"
    sed -Ei 's/^(writeInterval[[:space:]]+)[^;]+;/\11e-5;/' "$case_dir/system/controlDict"
    bash "$ROOT/run_case_orfoz.sh" "$case_dir"
done

test -f "$TESTROOT/serial_nominal_M1_C1/RUN_COMPLETE"
test -f "$TESTROOT/parallel16_high_stiffness_M1_C16/RUN_COMPLETE"
test -f "$TESTROOT/parallel16_high_stiffness_M1_C16/log.decomposePar"
test -f "$TESTROOT/parallel16_high_stiffness_M1_C16/log.reconstructPar"

echo "============================================================"
echo "V2312 ORFOZ PREFLIGHT PASSED"
echo "============================================================"
echo "Temporary results: $TESTROOT"
