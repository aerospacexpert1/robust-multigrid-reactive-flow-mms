#!/bin/bash
#SBATCH --job-name=RMT3H_BUILD
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/slurm-RMT3H_BUILD-%j.out
#SBATCH --error=logs/slurm-RMT3H_BUILD-%j.err

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "${ROOT}"

# Build jobs must obey the same Orfoz 56/112 CPU-per-node allocation rule.
# Compilation itself remains a normal single compiler process.
./build.sh
