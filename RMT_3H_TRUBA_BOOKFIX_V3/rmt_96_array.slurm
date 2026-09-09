#!/bin/bash
#SBATCH --job-name=RMTV3_96
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-95%3
#SBATCH --output=logs/slurm-RMTV3_96-%A_%a.out
#SBATCH --error=logs/slurm-RMTV3_96-%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${ROOT}"
mkdir -p logs
./run_case.sh "${SLURM_ARRAY_TASK_ID}"
