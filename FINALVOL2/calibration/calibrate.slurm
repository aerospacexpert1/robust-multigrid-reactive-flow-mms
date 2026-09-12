#!/bin/bash
#SBATCH --job-name=RMTV2_CAL
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=slurm-RMTV2_CAL_%j.out
#SBATCH --error=slurm-RMTV2_CAL_%j.err
set -euo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"
./build.sh
python3 calibration/run_smoothing_calibration.py --threads 16 --candidates 8,10,12,14,16
