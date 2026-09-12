#!/bin/bash
#SBATCH --job-name=RMTV2_BUILD
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=slurm-RMTV2_BUILD_%j.out
#SBATCH --error=slurm-RMTV2_BUILD_%j.err
set -euo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"
./build.sh
