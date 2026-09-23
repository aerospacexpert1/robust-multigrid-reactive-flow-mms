#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')

RMT = r'''
/* V2 RMT production-flow correction:
   1) form the fine-grid residual of the current solution;
   2) solve a zero-initialized correction equation recursively with the factor-three
      3x3 shifted coarse-grid family;
   3) inject coarse corrections without per-level damping;
   4) apply omega=0.005 only at the top level with an 8-step descent line search;
   5) perform the production six post-correction RBGS sweeps.
   The benchmark boundary traces are Dirichlet, so the correction BC is homogeneous
   Dirichlet rather than the mixed physical pressure BC of the opposed-flow case. */
static void rmt_error_recursive(Sys*s,int level){
    if(level>=5 || s->nx<12 || s->ny<8){rbgs(s,48);return;}
    double*r=residual_array(s);
    for(int rep=0;rep<2;rep++){
        for(int sx=0;sx<3;sx++)for(int sy=0;sy<3;sy++){
            Sys c; make_coarse(s,&c,3);
            for(int ic=1;ic<c.nx-1;ic++)for(int jc=1;jc<c.ny-1;jc++){
                int i0=sx+1+3*(ic-1),j0=sy+1+3*(jc-1);
                double sum=0.0; int cnt=0;
                int ilo=i0-1,ihi=i0+1,jlo=j0-1,jhi=j0+1;
                if(ilo<1)ilo=1;if(jlo<1)jlo=1;
                if(ihi>s->nx-2)ihi=s->nx-2;if(jhi>s->ny-2)jhi=s->ny-2;
                for(int fi=ilo;fi<=ihi;fi++)for(int fj=jlo;fj<=jhi;fj++){
                    sum+=r[IDX(fi,fj,s->nx)];cnt++;
                }
                c.b[IDX(ic,jc,c.nx)]=cnt?sum/(double)cnt:0.0;
            }
            correction_bc_zero(&c);
            rmt_error_recursive(&c,level+1);
            for(int ic=1;ic<c.nx-1;ic++)for(int jc=1;jc<c.ny-1;jc++){
                int fi=sx+1+3*(ic-1),fj=sy+1+3*(jc-1);
                if(fi>0&&fj>0&&fi<s->nx-1&&fj<s->ny-1)
                    s->q[IDX(fi,fj,s->nx)] += c.q[IDX(ic,jc,c.nx)];
            }
            freesys(&c);
        }
    }
    rbgs(s,3);
    free(r);
}

static void rmt_cycle_ref(Sys*s){
    double before=residual(s);
    double*r=residual_array(s);
    Sys corr; make_coarse(s,&corr,1);
    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
        corr.b[IDX(i,j,corr.nx)] = r[IDX(i,j,s->nx)];
    correction_bc_zero(&corr);
    rmt_error_recursive(&corr,0);

    double omegaTry=0.005;
    int accepted=0;
    for(int ls=0;ls<8;ls++){
        for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
            s->q[IDX(i,j,s->nx)] += omegaTry*corr.q[IDX(i,j,corr.nx)];
        double after=residual(s);
        if(isfinite(after) && (before<=0.0 || after<before)){accepted=1;break;}
        for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
            s->q[IDX(i,j,s->nx)] -= omegaTry*corr.q[IDX(i,j,corr.nx)];
        omegaTry*=0.5;
    }
    (void)accepted;
    rbgs(s,6);
    freesys(&corr);
    free(r);
}
'''

for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()
    pat=r'static void rmt_error_recursive\(Sys\*s,int level\)\{.*?\nstatic int solve\(Sys\*s,const char\*solver,double tol,int maxit\)'
    m=re.search(pat,s,flags=re.S)
    if not m:
        raise SystemExit(f'RMT block not found in {c}')
    s=re.sub(pat,lambda _: RMT+'\nstatic int solve(Sys*s,const char*solver,double tol,int maxit)',s,flags=re.S)
    s=s.replace('else {rmt_cycle_ref(s);it+=99;}','else {rmt_cycle_ref(s);it+=1;}')
    c.write_text(s)

print('V2_RMTFIX_COMPLETE')
