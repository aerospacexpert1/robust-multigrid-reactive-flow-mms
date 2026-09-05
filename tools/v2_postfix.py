#!/usr/bin/env python3
from pathlib import Path
ROOT=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')
if not ROOT.exists(): raise SystemExit('V2 root missing')
for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()
    # re.sub replacement processing in v2_patch turns C \\n escapes into literal newlines.
    s=s.replace('x,y,variable,exact,numerical,error\n");for', 'x,y,variable,exact,numerical,error\\n");for')
    s=s.replace('%.17g,%.17g,%s,%.17g,%.17g,%.17g\n",x,y,vv', '%.17g,%.17g,%s,%.17g,%.17g,%.17g\\n",x,y,vv')
    s=s.replace('T_rel_L2,rho_rel_L2\n");fprintf', 'T_rel_L2,rho_rel_L2\\n");fprintf')
    s=s.replace('%.12e,%.12e\n",BENCH_ID', '%.12e,%.12e\\n",BENCH_ID')
    s=s.replace('cont=%.3e time=%.6f\n",BENCH_ID', 'cont=%.3e time=%.6f\\n",BENCH_ID')
    c.write_text(s)
# Common engine parses the B branch at compile time for every BENCH_ID, so provide
# harmless zero pressure-gradient stubs outside B.
for h in ROOT.glob('Benchmark_*/src/mms_generated.h'):
    if 'Benchmark_B_Momentum' not in str(h):
        s=h.read_text()
        if 'mms_pdx(double' not in s:
            s += '\nstatic inline double mms_pdx(double x,double y,int s){(void)x;(void)y;(void)s;return 0.0;}\n'
            s += 'static inline double mms_pdy(double x,double y,int s){(void)x;(void)y;(void)s;return 0.0;}\n'
        h.write_text(s)
print('V2_POSTFIX_COMPLETE')
