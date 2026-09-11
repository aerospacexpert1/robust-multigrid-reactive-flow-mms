#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
[[ -f campaign_manifest.csv ]] || ./rebuild_manifest.sh

BUILD_JOB=$(sbatch --parsable build_rmt.slurm)
J47=$(sbatch --parsable --dependency=afterok:${BUILD_JOB} --export=ALL,TASK_ID=47 run_one_case.slurm)
J95=$(sbatch --parsable --dependency=afterok:${BUILD_JOB} --export=ALL,TASK_ID=95 run_one_case.slurm)
J96=$(sbatch --parsable --dependency=afterok:${J47}:${J95} rmt_96_array.slurm)

{
  echo "BUILD_JOB=$BUILD_JOB"
  echo "TASK47_JOB=$J47"
  echo "TASK95_JOB=$J95"
  echo "FULL96_JOB=$J96"
  date
} | tee campaign_jobs.txt

echo "Submitted gated campaign: build -> (TASK47,TASK95) -> 96-case array"
