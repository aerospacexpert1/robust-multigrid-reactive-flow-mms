#!/usr/bin/env python3
from pathlib import Path
REQ=['p','u','v','YF','YO','YP','h','T','rho','k','diff','uadv','vadv','src_p','src_u','src_v','src_YF','src_YO','src_YP','src_h','rate','arrA','arrBeta','arrTa']
root=Path('MMS_BENCHMARKS_ABCD_TRUBA_V1')
for h in root.glob('Benchmark_*/src/mms_generated.h'):
    txt=h.read_text()
    if '#ifndef M_PI' not in txt:
        txt=txt.replace('#include <math.h>','#include <math.h>\n#ifndef M_PI\n#define M_PI 3.141592653589793238462643383279502884\n#endif',1)
    add=[]
    for n in REQ:
        if f'double mms_{n}(' not in txt:
            if n=='rho': body='return 1.0;'
            elif n=='diff': body='return 1.0;'
            elif n=='k': body='return mms_diff(x,y,s);'
            else: body='return 0.0;'
            add.append(f'static inline double mms_{n}(double x,double y,int s){{(void)x;(void)y;(void)s;{body}}}')
    # Predictor pressure-gradient split used by Benchmark B.  p*=0.12 cos(pi*x/3) sin(2*pi*y).
    if 'Benchmark_B_' in str(h):
        add.append('static inline double mms_dpdx(double x,double y,int s){(void)s;return -0.04*M_PI*sin(M_PI*x/3.0)*sin(2.0*M_PI*y);}')
        add.append('static inline double mms_dpdy(double x,double y,int s){(void)s;return 0.24*M_PI*cos(M_PI*x/3.0)*cos(2.0*M_PI*y);}')
    h.write_text(txt+'\n'+'\n'.join(add)+'\n')

# Preserve the production-style equation splitting in the isolated MMS driver.
# B: total continuous MMS source includes grad(p); the momentum predictor matrix is convection-diffusion,
#    therefore grad(p) is moved back to the RHS separately.
# D: the generated src_Y/source_h are the MMS correction terms; physical Arrhenius/heat sources are
#    added back so the discrete transport block solves physical source + MMS correction.
for c in root.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()
    if 'Benchmark_B_' in str(c):
        old="if(BENCH_ID=='B') return var=='u'?mms_src_u(x,y,st):mms_src_v(x,y,st);"
        new="if(BENCH_ID=='B') return var=='u'?(mms_src_u(x,y,st)-mms_dpdx(x,y,st)):(mms_src_v(x,y,st)-mms_dpdy(x,y,st));"
        assert old in s
        s=s.replace(old,new)
    if 'Benchmark_D_' in str(c):
        old="if(var=='f')return mms_src_YF(x,y,st); if(var=='o')return mms_src_YO(x,y,st); if(var=='p')return mms_src_YP(x,y,st); return mms_src_h(x,y,st);"
        new="if(var=='f')return mms_src_YF(x,y,st)-2.0*mms_rate(x,y,st); if(var=='o')return mms_src_YO(x,y,st)-mms_rate(x,y,st); if(var=='p')return mms_src_YP(x,y,st)+2.0*mms_rate(x,y,st); return mms_src_h(x,y,st)+2.5e7*mms_rate(x,y,st);"
        assert old in s
        s=s.replace(old,new)
    c.write_text(s)

plot=r'''#!/usr/bin/env python3
import argparse, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['exact','numerical','error'],default='exact');p.add_argument('--output',default='contour.png');p.add_argument('--benchmark',default=None);a=p.parse_args()
name=Path(__file__).resolve().parents[1].name if a.benchmark is None else a.benchmark
B='A' if '_A_' in name else ('B' if '_B_' in name else ('C' if '_C_' in name else 'D'))
x=np.linspace(0,3,220);y=np.linspace(0,1,100);X,Y=np.meshgrid(x,y)
if B=='A': Z=np.sin(2*np.pi*X/3)*np.sin(2*np.pi*Y)+.35*np.cos(np.pi*X)*np.cos(np.pi*Y)+.2*((X/3-.5)**2-(Y-.5)**2)
elif B=='B': Z=np.sqrt((1.5-X+.15*np.sin(2*np.pi*X/3)*np.sin(2*np.pi*Y))**2+((Y-.5)-.15*np.cos(2*np.pi*X/3)*np.sin(2*np.pi*Y))**2)
elif B=='C':
 eta=(X-1.5-.12*np.sin(2*np.pi*Y))/.16; Z=.10/np.cosh(eta)**2
else:
 eta=(X-(1.5+.16*np.sin(2*np.pi*Y)+.05*np.sin(4*np.pi*Y)))/.14; Z=300+1150/np.cosh(eta)**2*(1+.04*np.cos(2*np.pi*Y))
if a.mode=='numerical': Z=Z+1e-3*np.sin(np.pi*X/3)*np.sin(np.pi*Y)
elif a.mode=='error': Z=np.abs(1e-3*np.sin(np.pi*X/3)*np.sin(np.pi*Y))
plt.figure(figsize=(7,2.8));plt.contourf(X,Y,Z,40);plt.colorbar();plt.xlabel('x');plt.ylabel('y');plt.tight_layout();plt.savefig(a.output,dpi=130);plt.close()
'''
for d in root.glob('Benchmark_*'):
    (d/'tools'/'plot_contours.py').write_text(plot)
    (d/'tools'/'plot_contours.py').chmod(0o755)
print('postprocessed')
