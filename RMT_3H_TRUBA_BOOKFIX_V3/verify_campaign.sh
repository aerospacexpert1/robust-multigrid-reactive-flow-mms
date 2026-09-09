#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT}"
[[ -f src/opposedflow_rmt3h_LEGACY_V2.c ]]
[[ -f src/rmt_book2d_impl.h ]]
[[ -f tools/patch_bookfix_v3.py ]]
[[ -f campaign_manifest.csv ]] || ./rebuild_manifest.sh
rows=$(( $(wc -l < campaign_manifest.csv) - 1 ))
[[ "${rows}" -eq 96 ]] || { echo "FAIL manifest rows=${rows}"; exit 2; }
# Known task mapping from the legacy campaign; these must not move.
grep -q '^28,FLOW_NOMINAL_M4_T1,FLOW,NOMINAL,M4,864,288,1,' campaign_manifest.csv
grep -q '^40,FLOW_HIGH_M3_T1,FLOW,HIGH,M3,432,144,1,' campaign_manifest.csv
grep -q '^44,FLOW_HIGH_M4_T1,FLOW,HIGH,M4,864,288,1,' campaign_manifest.csv
grep -q '^9,FLOW_LOW_M3_T2,FLOW,LOW,M3,432,144,2,' campaign_manifest.csv
grep -q -- '-dt 1.0e-7 -maxCo 0.05' run_case.sh
grep -q -- '-pCycles 4 -sCycles 4' run_case.sh
grep -q -- '-rmtLevels 0 -rmtPost 8 -rmtCoarseSweeps 16 -rmtCoarseRepeats 1 -rmtOmega 1.0' run_case.sh
grep -q -- '-fieldOutput 0' run_case.sh
python3 -m py_compile tools/patch_bookfix_v3.py tools/extract_one_summary.py tools/collect_summary.py
bash -n build.sh run_case.sh rebuild_manifest.sh run_one_case.slurm
printf 'VERIFY PASS: 96-case mapping preserved; physics/timestep controls unchanged; BOOKFIX controls explicit.\n'
