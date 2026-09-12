#!/usr/bin/env python3
from pathlib import Path

ROOT=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')
for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()
    # The solver's `it` counter is a work-unit estimate (MG3V adds 46 per cycle),
    # not a literal outer-cycle counter.  24000 therefore caps MG3V at only about
    # 522 cycles on refined grids.  Preserve the 1e-7 algebraic tolerance and
    # increase only the safety budget so the MMS error is not contaminated by an
    # unconverged algebraic solve.
    s=s.replace('solve(&s,solver,1e-7,24000)', 'solve(&s,solver,1e-7,240000)')
    c.write_text(s)
print('V2_BUDGETFIX_COMPLETE')
