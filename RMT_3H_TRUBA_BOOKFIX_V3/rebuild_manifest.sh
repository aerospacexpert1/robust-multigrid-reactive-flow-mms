#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${ROOT}/campaign_manifest.csv"
printf '%s\n' 'task_id,case_id,family,level,mesh,Nx,Ny,threads,vIn,arrA,arrBeta,arrTa,endTime' > "${OUT}"
task=0
emit_block() {
    local family="$1" level="$2" vin="$3" A="$4" beta="$5" Ta="$6"
    local mesh nx ny th
    for spec in 'M1 108 36' 'M2 216 72' 'M3 432 144' 'M4 864 288'; do
        read -r mesh nx ny <<< "${spec}"
        for th in 1 2 4 16; do
            printf '%d,%s_%s_%s_T%d,%s,%s,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,2\n' \
                "${task}" "${family}" "${level}" "${mesh}" "${th}" \
                "${family}" "${level}" "${mesh}" "${nx}" "${ny}" "${th}" \
                "${vin}" "${A}" "${beta}" "${Ta}" >> "${OUT}"
            task=$((task+1))
        done
    done
}
# FLOW sweep: chemistry fixed at S1; inlet speed varied.
emit_block FLOW LOW     0.05 1e7 5 1000
emit_block FLOW NOMINAL 0.10 1e7 5 1000
emit_block FLOW HIGH    0.20 1e7 5 1000
# CHEM sweep: flow fixed at nominal inlet speed; Arrhenius stiffness varied.
emit_block CHEM S1 0.10 1e7 5 1000
emit_block CHEM S2 0.10 1e8 5 500
emit_block CHEM S3 0.10 1e9 5 100
[[ "${task}" -eq 96 ]] || { echo "manifest generation error: ${task}"; exit 2; }
echo "Wrote ${OUT} with ${task} cases"
