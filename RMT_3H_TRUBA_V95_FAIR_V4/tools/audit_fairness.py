#!/usr/bin/env python3
"""Fail-fast audit proving V95 physical/timestep kernels are unchanged by RMT generation."""
from pathlib import Path
import hashlib, sys

if len(sys.argv) != 3:
    raise SystemExit("usage: audit_fairness.py FROZEN_V95.c GENERATED_RMT.c")
a = Path(sys.argv[1]).read_text()
b = Path(sys.argv[2]).read_text()

def block(text, marker):
    p = text.find(marker)
    if p < 0:
        raise SystemExit(f"missing marker: {marker}")
    q = text.find('{', p)
    if q < 0:
        raise SystemExit(f"missing opening brace after: {marker}")
    depth = 0
    i = q
    in_str = in_chr = False
    esc = False
    while i < len(text):
        ch = text[i]
        if in_str:
            if esc: esc = False
            elif ch == '\\': esc = True
            elif ch == '"': in_str = False
        elif in_chr:
            if esc: esc = False
            elif ch == '\\': esc = True
            elif ch == "'": in_chr = False
        else:
            if ch == '"': in_str = True
            elif ch == "'": in_chr = True
            elif ch == '{': depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0:
                    return text[p:i+1]
        i += 1
    raise SystemExit(f"unclosed block: {marker}")

def sha(x): return hashlib.sha256(x.encode()).hexdigest()

functions = [
    "static void momentum_predictor(",
    "static void build_pressure_rhs(",
    "static void correct_velocity_and_face_flux(",
    "static void build_scalar_rhs_openfoam_like(",
    "static void solve_scalar_helmholtz(",
    "static void chemistry_subcycle(",
    "static void update_thermo_transport_fields(",
]
for marker in functions:
    x, y = block(a, marker), block(b, marker)
    if x != y:
        raise SystemExit(f"FAIRNESS FAIL: non-pressure kernel changed: {marker}\nbase={sha(x)}\nrmt ={sha(y)}")
    print(f"FAIRNESS MATCH {marker} sha256={sha(x)}")

# Compare the complete timed physical timestep loop exactly.
loop_marker = "while (time < c.endTime - 1e-15)"
x, y = block(a, loop_marker), block(b, loop_marker)
if x != y:
    raise SystemExit(f"FAIRNESS FAIL: physical timestep loop changed\nbase={sha(x)}\nrmt ={sha(y)}")
print(f"FAIRNESS MATCH timestep_loop sha256={sha(x)}")

# Positive and negative invariants for the generated solver.
for token in [
    '#include "rmt_book2d_impl.h"',
    "rmt_book_build_correction",
    "RMT_PRESSURE_REJECT",
    "RMT_BOOK_DIAGNOSTICS",
    "Non-incremental projection: pressure is intentionally excluded here",
    "correct_velocity_and_face_flux(&g,&ph,&f,dt);",
]:
    if token not in b:
        raise SystemExit(f"FAIRNESS FAIL: generated invariant missing: {token}")

# The old failed-RMT coupling must not reappear: predictor stays pressure-free.
mp = block(b, "static void momentum_predictor(")
if "dpdx" in mp or "dpdy" in mp or "gradp" in mp:
    raise SystemExit("FAIRNESS FAIL: pressure gradient found inside V95 momentum predictor")

print("FAIRNESS_AUDIT PASS: V95 physical kernels and timestep loop are byte-identical; pressure solver is the isolated change.")
