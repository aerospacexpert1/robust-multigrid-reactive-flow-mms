#!/usr/bin/env python3
"""Final fail-fast edits after deterministic V95->RMT generation.

This intentionally restores the frozen timestep-loop progress label so the complete
physical while-loop remains byte-identical to the V95 baseline. It also makes failure
to meet the frozen V95 pressure tolerances a hard error rather than a silent success.
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
print("Finalized generated V95 RMT source: frozen timestep-loop labels restored; pressure nonconvergence is fatal")
