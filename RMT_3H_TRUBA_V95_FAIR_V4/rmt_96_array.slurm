#!/bin/bash
#SBATCH --job-name=RMTV95_96
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=3-00:00:00
#SBATCH --array=0-95%3
#SBATCH --output=slurm-RMTV95_96_%A_%a.out
#SBATCH --error=slurm-RMTV95_96_%A_%a.err
set -euo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"
[[ -x build/opposedflow_rmt3h_v95_fair ]] || { echo "ERROR: build/opposedflow_rmt3h_v95_fair is missing" >&2; exit 10; }
[[ -f campaign_manifest.csv ]] || { echo "ERROR: campaign_manifest.csv is missing; run ./rebuild_manifest.sh" >&2; exit 11; }
[[ -x run_case.sh ]] || { echo "ERROR: run_case.sh not executable" >&2; exit 12; }
export ROOT="$SLURM_SUBMIT_DIR"
echo "RMT V95 FAIR V4 job=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID} task=$SLURM_ARRAY_TASK_ID node=$(hostname)"
./run_case.sh "$SLURM_ARRAY_TASK_ID"
