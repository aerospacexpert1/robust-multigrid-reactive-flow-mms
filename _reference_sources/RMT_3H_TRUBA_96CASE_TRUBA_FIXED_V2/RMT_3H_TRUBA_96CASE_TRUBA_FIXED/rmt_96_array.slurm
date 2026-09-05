#!/bin/bash
#SBATCH --job-name=RMT3H_96
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-95%3
#SBATCH --output=logs/slurm-RMT3H_96-%A_%a.out
#SBATCH --error=logs/slurm-RMT3H_96-%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "${ROOT}"

# Orfoz requires a 56/112 CPU-per-node allocation.  The allocation is
# intentionally independent from the benchmark OpenMP thread count.
# run_case.sh reads the manifest and sets OMP_NUM_THREADS to 1, 2, 4, or 16.

echo "======================================================================"
echo "SLURM job"
echo "======================================================================"
echo "JOB ID        : ${SLURM_JOB_ID}"
echo "ARRAY JOB ID  : ${SLURM_ARRAY_JOB_ID:-none}"
echo "ARRAY TASK    : ${SLURM_ARRAY_TASK_ID}"
echo "NODE           : $(hostname)"
echo "CPUS allocated : ${SLURM_CPUS_PER_TASK:-unknown}"
echo "NOTE           : benchmark OpenMP threads come from campaign_manifest.csv"
echo "======================================================================"

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
echo "Working dir    : $(pwd)"
./run_case.sh "${SLURM_ARRAY_TASK_ID}"
