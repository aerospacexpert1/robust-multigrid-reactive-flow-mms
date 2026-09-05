#!/usr/bin/env python3
from pathlib import Path
c=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2/Benchmark_B_Momentum/src/mms_solver.c')
s=c.read_text()
old="s->b[k]=srcvar(var,xx,yy,st)*vol;if(BENCH_ID=='B')s->b[k]-=(var=='u'?mms_pdx(xx,yy,st):mms_pdy(xx,yy,st))*vol; (void)rho;(void)ua;(void)va;(void)dif;"
new="s->b[k]=srcvar(var,xx,yy,st)*vol; (void)rho;(void)ua;(void)va;(void)dif;"
if old not in s:
    raise SystemExit('duplicate B pressure subtraction not found')
s=s.replace(old,new)
c.write_text(s)

# The postprocess stage already implements the production predictor split:
# srcvar(B) = continuous total MMS source - grad(p).  The FV assembly must not
# subtract grad(p) a second time.
doc=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2/V2_METHOD_CHANGES.md')
d=doc.read_text()
d += '''\n## Benchmark B pressure-gradient split\nThe continuous momentum MMS source contains the analytic pressure gradient. `postprocess_package.py` moves that known pressure-gradient term out of the predictor transport operator exactly once. V2 removes the former second subtraction in FV assembly, which otherwise produced a non-vanishing refinement error.\n'''
doc.write_text(d)
print('V2_BFIX_COMPLETE')
