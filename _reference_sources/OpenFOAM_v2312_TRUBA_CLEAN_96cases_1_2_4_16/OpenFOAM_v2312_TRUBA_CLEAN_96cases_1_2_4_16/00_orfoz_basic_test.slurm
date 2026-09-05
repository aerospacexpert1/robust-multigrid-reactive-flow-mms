#!/bin/bash
#SBATCH --job-name=orfoz_basic
#SBATCH --partition=orfoz
#SBATCH --exclude=orfoz359
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=orfoz_basic_%j.out
#SBATCH --error=orfoz_basic_%j.err

echo "BASIC ORFOZ TEST STARTED"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Job: $SLURM_JOB_ID"
echo "Tasks: $SLURM_NTASKS"
hostname
sleep 5
echo "BASIC ORFOZ TEST PASSED"
