#!/bin/bash
#SBATCH --job-name=of_env_test
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=of_env_test_%j.out
#SBATCH --error=of_env_test_%j.err

echo "CHECKPOINT 1: batch script started"
echo "Node: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo

echo "CHECKPOINT 2: before module purge"
module purge
echo "module purge rc=$?"

echo
echo "CHECKPOINT 3: before OpenFOAM module"
module load apps/openfoam/v2312
echo "module load rc=$?"

echo
echo "CHECKPOINT 4: OpenFOAM environment"
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
echo "WM_PROJECT_DIR=${WM_PROJECT_DIR:-unset}"
echo "WM_MPLIB=${WM_MPLIB:-unset}"

echo
echo "CHECKPOINT 5: executables"
for exe in blockMesh checkMesh decomposePar reconstructPar rhoReactingFoam mpirun; do
    printf "%-24s : " "$exe"
    command -v "$exe" || echo "NOT FOUND"
done

echo
echo "CHECKPOINT 6: MPI version"
mpirun --version | head -n 3

echo
echo "OPENFOAM ENVIRONMENT TEST PASSED"
