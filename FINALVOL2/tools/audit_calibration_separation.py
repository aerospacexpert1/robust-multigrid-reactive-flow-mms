#!/usr/bin/env python3
"""Guard against accidental benchmark tuning or physics drift in FINALVOL2."""
from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
run=(root/"run_case.sh").read_text()
gen=(root/"tools"/"generate_finalvol2.py").read_text()
cal=(root/"calibration"/"run_smoothing_calibration.py").read_text()

must_run=[
    "P_REL_TOL=${P_REL_TOL:-1e-4}",
    "P_ABS_TOL=${P_ABS_TOL:-1e-6}",
    "RMT_SMOOTH_SWEEPS=${RMT_SMOOTH_SWEEPS:-16}",
    "-perfectGas 0 -variableCp 0 -sutherland 0 -rhoRelax 1",
    "-sCycles 4 -sSweeps 5",
]
for x in must_run:
    if x not in run: raise SystemExit(f"CALIBRATION_SEPARATION FAIL missing production invariant: {x}")

# Inspect assignment/control lines only; do not use DOTALL across the whole runner.
for line in run.splitlines():
    if "RMT_SMOOTH_SWEEPS" not in line:
        continue
    low=line.lower()
    if any(tok in low for tok in ("physical_case","mesh_name","case_id","bc1","bc2","bc3","s1_","s2_","s3_")):
        raise SystemExit(f"CALIBRATION_SEPARATION FAIL case-dependent tuning line: {line}")

for forbidden in ["write_text(run_case", "open('../run_case.sh','w'", "sed -i", "git commit", "git push"]:
    if forbidden in cal:
        raise SystemExit(f"CALIBRATION_SEPARATION FAIL calibration mutation token: {forbidden}")

if gen.count("cfg.smoothSweeps =") != 1:
    raise SystemExit("CALIBRATION_SEPARATION FAIL expected one global smoothSweeps assignment")
print("CALIBRATION_SEPARATION PASS: calibration is advisory; production parameter remains global and frozen at 16.")
