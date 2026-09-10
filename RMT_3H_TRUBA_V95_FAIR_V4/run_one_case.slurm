#!/bin/bash
#SBATCH --job-name=RMTV95_ONE
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --output=slurm-RMTV95_ONE_%j.out
#SBATCH --error=slurm-RMTV95_ONE_%j.err
set -euo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"
[[ -x build/opposedflow_rmt3h_v95_fair ]] || { echo "ERROR: build first" >&2; exit 10; }
[[ -f campaign_manifest.csv ]] || ./rebuild_manifest.sh
TASK_ID=${TASK_ID:-3}
echo "Running RMT V95 FAIR single case TASK_ID=$TASK_ID"
./run_case.sh "$TASK_ID"
