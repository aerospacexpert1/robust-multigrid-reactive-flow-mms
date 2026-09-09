#!/bin/bash
#SBATCH --job-name=RMTV3_ONE
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm-RMTV3-one-%j.out
#SBATCH --error=logs/slurm-RMTV3-one-%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT}"
mkdir -p logs
TASK_ID="${TASK_ID:-0}"
echo "Running BOOKFIX_V3 task ${TASK_ID}"
./run_case.sh "${TASK_ID}"
