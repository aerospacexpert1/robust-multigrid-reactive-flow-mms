#!/bin/bash
#SBATCH --job-name=MG2W_BUILD
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --output=slurm-MG2W_BUILD_%j.out
#SBATCH --error=slurm-MG2W_BUILD_%j.err

set -eo pipefail
module purge
module load comp/gcc/12.3.0
cd "$SLURM_SUBMIT_DIR"

mkdir -p build
start_ns=$(date +%s%N)
gcc -O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp \
    src/opposedflow_mg2w.c -lm -o build/opposedflow_mg2w
end_ns=$(date +%s%N)
compile_s=$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN {printf "%.9f", (e-s)/1.0e9}')

cat > build/build_info.txt <<INFO
build_job_id=${SLURM_JOB_ID}
node=$(hostname)
partition=${SLURM_JOB_PARTITION}
compiler=$(gcc --version | head -n 1)
flags=-O3 -march=native -std=c11 -Wall -Wextra -Wpedantic -fopenmp
compile_seconds=${compile_s}
source_sha256=$(sha256sum src/opposedflow_mg2w.c | awk '{print $1}')
build_time=$(date --iso-8601=seconds)
INFO

echo "Build complete: build/opposedflow_mg2w"
echo "Compile time: ${compile_s} s"
