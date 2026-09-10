#!/usr/bin/env python3
"""Final fail-fast edits after deterministic V95->RMT generation.

This restores the frozen timestep-loop progress label so the complete physical
while-loop remains byte-identical to the V95 baseline, installs the production RMT
smoothing defaults established by the exact V95 first-pressure diagnostic, and makes
failure to meet the frozen V95 pressure tolerances a hard error.
"""
from pathlib import Path
import sys

if len(sys.argv) != 2:
    raise SystemExit("usage: finalize_generated_v95.py GENERATED.c")
p = Path(sys.argv[1])
s = p.read_text()

# Restore the exact frozen progress-loop label; explicit RMT diagnostics are printed
# after the timed physical solve.
s = s.replace("rhoRMT", "rhoRBGS")

# Production default selected from the exact V95 pressure diagnostic:
# post=8 rejected; post=12 only barely met relTol=1e-4; post=16 reached 3.22e-5
# in 7 cycles with omega=1 and no backtracking.  Coarsest sweeps are 2*post.
default_old = "c->mgPostSmooth = 8;\n    c->mgCoarseSweeps = 16;"
default_new = "c->mgPostSmooth = 16;\n    c->mgCoarseSweeps = 32;"
if default_old not in s:
    raise SystemExit("finalization refused: generated RMT default block missing")
s = s.replace(default_old, default_new, 1)

anchor = '''    if (stats.cycles > 0 && stats.initialResidual > 0.0 && stats.absResidual > 0.0)
        stats.convergenceFactor = pow(stats.absResidual/stats.initialResidual, 1.0/(double)stats.cycles);
'''
if anchor not in s:
    raise SystemExit("finalization refused: pressure convergence anchor missing")
check = '''    if (stats.absResidual > c->pressureAbsTol && stats.relResidual > c->pressureRelTol) {
        fprintf(stderr, "RMT_PRESSURE_NOT_CONVERGED cycles=%d absResidual=%.17g relResidual=%.17g absTol=%.17g relTol=%.17g\\n",
                stats.cycles, stats.absResidual, stats.relResidual,
                c->pressureAbsTol, c->pressureRelTol);
        die("RMT pressure solve hit max cycles without meeting the frozen V95 tolerance");
    }

'''
s = s.replace(anchor, check + anchor, 1)
p.write_text(s)
print("Finalized generated V95 RMT source: post=16/coarse=32; frozen timestep loop restored; pressure nonconvergence is fatal")
