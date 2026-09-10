#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
OUT="$ROOT/campaign_manifest.csv"
echo 'array_index,case_id,study_family,group_dir,physical_case,mesh_name,Nx,Ny,velocity_m_per_s,stiffness_class,A,beta,Ta,threads,core_label' > "$OUT"
idx=0
emit() {
  local study="$1" group="$2" phys="$3" vel="$4" stiff="$5" A="$6" beta="$7" Ta="$8"
  local mesh nx ny th
  for spec in 'M1_108x36 108 36' 'M2_216x72 216 72' 'M3_432x144 432 144' 'M4_864x288 864 288'; do
    read -r mesh nx ny <<< "$spec"
    for th in 1 2 4 16; do
      local core="C${th}"
      local cid="${study}__${phys}__${mesh}__${core}"
      printf '%d,%s,%s,%s,%s,%s,%d,%d,%.17g,%s,%.17g,%.17g,%.17g,%d,%s\n' \
        "$idx" "$cid" "$study" "$group" "$phys" "$mesh" "$nx" "$ny" "$vel" "$stiff" "$A" "$beta" "$Ta" "$th" "$core" >> "$OUT"
      idx=$((idx+1))
    done
  done
}
emit FLOW_INTENSITY FLOW BC1_LOW_v005 0.05 nominal_chemistry 1e7 5 1000
emit FLOW_INTENSITY FLOW BC2_NOMINAL_v010 0.10 nominal_chemistry 1e7 5 1000
emit FLOW_INTENSITY FLOW BC3_HIGH_v020 0.20 nominal_chemistry 1e7 5 1000
emit REACTION_STIFFNESS STIFFNESS S1_LOW_A1e7_Ta1000 0.10 low 1e7 5 1000
emit REACTION_STIFFNESS STIFFNESS S2_INTERMEDIATE_A1e8_Ta500 0.10 intermediate 1e8 5 500
emit REACTION_STIFFNESS STIFFNESS S3_HIGH_A1e9_Ta100 0.10 high 1e9 5 100
[[ "$idx" -eq 96 ]] || { echo "manifest generation failed: $idx" >&2; exit 2; }
echo "Wrote $OUT with $idx matched cases"
