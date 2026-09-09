#!/bin/bash
#SBATCH --job-name=build_RMTV3
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=build-RMTV3-%j.out
#SBATCH --error=build-RMTV3-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
./rebuild_manifest.sh
./verify_campaign.sh
./build.sh
