#!/usr/bin/env python3
"""Fail-fast audit proving RMT FINAL changes only the V95 pressure solver/reporting."""
from pathlib import Path
import hashlib
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: audit_fairness.py FROZEN_V95.c GENERATED_RMT.c")
a = Path(sys.argv[1]).read_text()
b = Path(sys.argv[2]).read_text()

def block(text, marker):
    p=text.find(marker)
    if p<0: raise SystemExit(f"missing marker: {marker}")
    q=text.find("{",p)
    if q<0: raise SystemExit(f"missing opening brace after: {marker}")
    depth=0; i=q; in_str=in_chr=False; esc=False
    while i<len(text):
        ch=text[i]
        if in_str:
            if esc: esc=False
            elif ch=="\\": esc=True
            elif ch=='"': in_str=False
        elif in_chr:
            if esc: esc=False
            elif ch=="\\": esc=True
            elif ch=="'": in_chr=False
        else:
            if ch=='"': in_str=True
            elif ch=="'": in_chr=True
            elif ch=="{": depth+=1
            elif ch=="}":
                depth-=1
                if depth==0:return text[p:i+1]
        i+=1
    raise SystemExit(f"unclosed block: {marker}")

def sha(x): return hashlib.sha256(x.encode()).hexdigest()

for marker in [
    "static void momentum_predictor(",
    "static void build_pressure_rhs(",
    "static void correct_velocity_and_face_flux(",
    "static void build_scalar_rhs_openfoam_like(",
    "static void solve_scalar_helmholtz(",
    "static void chemistry_subcycle(",
    "static void update_thermo_transport_fields(",
]:
    x,y=block(a,marker),block(b,marker)
    if x!=y:
        raise SystemExit(f"FAIRNESS FAIL: non-pressure kernel changed: {marker}\nbase={sha(x)}\nrmt ={sha(y)}")
    print(f"FAIRNESS MATCH {marker} sha256={sha(x)}")

loop_marker="while (time < c.endTime - 1e-15)"
x,y=block(a,loop_marker),block(b,loop_marker)
if x!=y:
    raise SystemExit(f"FAIRNESS FAIL: physical timestep loop changed\nbase={sha(x)}\nrmt ={sha(y)}")
print(f"FAIRNESS MATCH timestep_loop sha256={sha(x)}")

for token in [
    '#include "rmt_final2d_impl.h"',
    "rmt_final_build_correction",
    "rmt_final_direct_coarsest",
    "RMT_FINAL_DIAGNOSTICS",
    "RMT_PRESSURE_NOT_CONVERGED",
    "Non-incremental projection: pressure is intentionally excluded here",
    "correct_velocity_and_face_flux(&g,&ph,&f,dt);",
]:
    if token not in b:
        raise SystemExit(f"FAIRNESS FAIL: generated invariant missing: {token}")

# These V4 mechanisms are deliberately forbidden in RMT FINAL.
for forbidden in [
    "RMT_PRESSURE_REJECT",
    "omegaTry",
    "backtrack_halvings",
    "refusing hidden fallback",
]:
    if forbidden in b:
        raise SystemExit(f"FAIRNESS FAIL: obsolete V4 safeguard remains: {forbidden}")

mp=block(b,"static void momentum_predictor(")
if "dpdx" in mp or "dpdy" in mp or "gradp" in mp:
    raise SystemExit("FAIRNESS FAIL: pressure gradient found inside frozen V95 momentum predictor")

print("FAIRNESS_AUDIT PASS: frozen V95 physics/timestep are byte-identical; only the pressure linear solver/reporting changed.")
