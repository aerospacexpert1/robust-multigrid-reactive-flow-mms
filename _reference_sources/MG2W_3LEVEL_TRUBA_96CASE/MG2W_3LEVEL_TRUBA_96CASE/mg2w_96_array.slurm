#!/bin/bash
#SBATCH --job-name=MG2W_96
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=3-00:00:00
#SBATCH --array=0-95%3
#SBATCH --output=slurm-MG2W_96_%A_%a.out
#SBATCH --error=slurm-MG2W_96_%A_%a.err

set -eo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"

[[ -x build/opposedflow_mg2w ]] || { echo "ERROR: build/opposedflow_mg2w is missing" >&2; exit 10; }
[[ -f campaign_manifest.csv ]] || { echo "ERROR: campaign_manifest.csv is missing" >&2; exit 11; }
[[ -x run_case.sh ]] || { echo "ERROR: run_case.sh is missing or not executable" >&2; exit 12; }

export ROOT="$SLURM_SUBMIT_DIR"

echo "============================================================"
echo "MG2W THREE-LEVEL W-CYCLE 96-CASE CAMPAIGN"
echo "============================================================"
echo "job=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID} task=$SLURM_ARRAY_TASK_ID"
echo "node=$(hostname) date=$(date --iso-8601=seconds)"
echo "Each array task is one physical/mesh/thread case with its own 3-day allocation."
echo "============================================================"

./run_case.sh "$SLURM_ARRAY_TASK_ID"
