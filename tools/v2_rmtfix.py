#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')

RMT = r'''
/* Production-aligned RMT benchmark adaptation.
   Invariants intentionally matched to the reacting-flow pressure RMT:
     - factor-three, 3x3 shifted coarse families;
     - no presmoothing;
     - control-volume averaging of the correction defect;
     - direct solve on every coarsest shifted grid;
     - coarsest-to-finest post-smoothing;
     - full correction added without interpolation;
     - one globally fixed 16-sweep postsmoothing count;
     - no damping, monotonic line search, rollback, or hidden fallback.
   The MMS benchmark uses homogeneous Dirichlet correction traces, whereas the
   opposed-flow pressure solver uses its production mixed correction BCs. */

static void rmt_direct_dense(Sys*s){
    const int ni=s->nx-2,nj=s->ny-2,n=ni*nj;
    if(n<=0)return;
    double *A=calloc((size_t)n*(size_t)n,sizeof(double));
    double *b=calloc((size_t)n,sizeof(double));
    if(!A||!b)die("RMT direct allocation");

    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++){
        const int k=IDX(i,j,s->nx);
        const int row=(j-1)*ni+(i-1);
        A[(size_t)row*n+row]=s->ap[k];
        b[row]=s->b[k];
        if(i>1)       A[(size_t)row*n+(row-1)] = -s->aw[k];
        if(i<s->nx-2) A[(size_t)row*n+(row+1)] = -s->ae[k];
        if(j>1)       A[(size_t)row*n+(row-ni)] = -s->as[k];
        if(j<s->ny-2) A[(size_t)row*n+(row+ni)] = -s->an[k];
    }

    for(int k=0;k<n;k++){
        int piv=k;
        double vmax=fabs(A[(size_t)k*n+k]);
        for(int r=k+1;r<n;r++){
            double v=fabs(A[(size_t)r*n+k]);
            if(v>vmax){vmax=v;piv=r;}
        }
        if(vmax<1e-30)die("RMT direct singular coarsest system");
        if(piv!=k){
            for(int c=k;c<n;c++){
                double t=A[(size_t)k*n+c];
                A[(size_t)k*n+c]=A[(size_t)piv*n+c];
                A[(size_t)piv*n+c]=t;
            }
            double tb=b[k];b[k]=b[piv];b[piv]=tb;
        }
        const double akk=A[(size_t)k*n+k];
        for(int r=k+1;r<n;r++){
            const double m=A[(size_t)r*n+k]/akk;
            if(m==0.0)continue;
            A[(size_t)r*n+k]=0.0;
            for(int c=k+1;c<n;c++)A[(size_t)r*n+c]-=m*A[(size_t)k*n+c];
            b[r]-=m*b[k];
        }
    }

    for(int r=n-1;r>=0;r--){
        double z=b[r];
        for(int c=r+1;c<n;c++)z-=A[(size_t)r*n+c]*b[c];
        b[r]=z/A[(size_t)r*n+r];
    }

    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++){
        const int row=(j-1)*ni+(i-1);
        s->q[IDX(i,j,s->nx)]=b[row];
    }
    correction_bc_zero(s);
    free(A);free(b);
}

static void rmt_error_recursive(Sys*s,int level){
    if(level>=5 || s->nx<12 || s->ny<8){
        rmt_direct_dense(s);
        return;
    }

    double*r=residual_array(s);
    for(int sy=0;sy<3;sy++)for(int sx=0;sx<3;sx++){
        Sys c; make_coarse(s,&c,3);
        for(int jc=1;jc<c.ny-1;jc++)for(int ic=1;ic<c.nx-1;ic++){
            const int i0=sx+1+3*(ic-1),j0=sy+1+3*(jc-1);
            int ilo=i0-1,ihi=i0+1,jlo=j0-1,jhi=j0+1;
            if(ilo<1)ilo=1;if(jlo<1)jlo=1;
            if(ihi>s->nx-2)ihi=s->nx-2;
            if(jhi>s->ny-2)jhi=s->ny-2;
            double sum=0.0;int cnt=0;
            for(int fj=jlo;fj<=jhi;fj++)for(int fi=ilo;fi<=ihi;fi++){
                sum+=r[IDX(fi,fj,s->nx)];
                cnt++;
            }
            c.b[IDX(ic,jc,c.nx)]=cnt?sum/(double)cnt:0.0;
        }
        correction_bc_zero(&c);
        rmt_error_recursive(&c,level+1);

        /* Index-space transfer: each shifted family writes directly to the
           corresponding fine-grid locations; no interpolation is used. */
        for(int jc=1;jc<c.ny-1;jc++)for(int ic=1;ic<c.nx-1;ic++){
            const int fi=sx+1+3*(ic-1),fj=sy+1+3*(jc-1);
            if(fi>0&&fj>0&&fi<s->nx-1&&fj<s->ny-1)
                s->q[IDX(fi,fj,s->nx)] += c.q[IDX(ic,jc,c.nx)];
        }
        freesys(&c);
    }
    free(r);

    /* Single global smoothing choice, matched to the calibrated production RMT. */
    rbgs(s,16);
}

static void rmt_cycle_ref(Sys*s){
    size_t n=(size_t)s->nx*s->ny;
    double*r=residual_array(s);
    Sys corr; make_coarse(s,&corr,1);
    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
        corr.b[IDX(i,j,corr.nx)] = r[IDX(i,j,s->nx)];
    correction_bc_zero(&corr);

    rmt_error_recursive(&corr,0);

    /* Full correction: no damping and no residual-dependent acceptance test. */
    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
        s->q[IDX(i,j,s->nx)] += corr.q[IDX(i,j,corr.nx)];

    freesys(&corr);
    free(r);
    (void)n;
}
'''

SOLVE = r'''
static int solve(Sys*s,const char*solver,double tol,int maxit){
    int it=0;
    if(strcmp(solver,"SG_RBGS")==0){
        while(it<maxit&&residual(s)>tol){
            rbgs(s,2);
            it+=2; /* SG iteration metric is RBGS sweeps. */
        }
        return it;
    }

    while(it<maxit&&residual(s)>tol){
        if(strcmp(solver,"MG2V")==0)
            mg_cycle_ref(s,2,2,1);
        else if(strcmp(solver,"MG2W")==0)
            mg_cycle_ref(s,2,2,2);
        else if(strcmp(solver,"MG3V")==0)
            mg_cycle_ref(s,3,2,1);
        else
            rmt_cycle_ref(s);

        /* For every multilevel solver, one iteration is one complete cycle.
           No solver is rescued by an SG fallback or rollback. */
        it+=1;
        if(!isfinite(residual(s)))break;
    }
    return it;
}
'''

for c in ROOT.glob('Benchmark_*/src/mms_solver.c'):
    s=c.read_text()

    # Replace only the RMT implementation while leaving the following solve()
    # declaration as a stable anchor.
    pat_rmt=r'static void rmt_error_recursive\(Sys\*s,int level\)\{.*?\nstatic int solve\(Sys\*s,const char\*solver,double tol,int maxit\)'
    if not re.search(pat_rmt,s,flags=re.S):
        raise SystemExit(f'RMT block not found in {c}')
    s=re.sub(
        pat_rmt,
        lambda _: RMT+'\nstatic int solve(Sys*s,const char*solver,double tol,int maxit)',
        s,flags=re.S
    )

    # Replace the complete solver dispatcher up to the next stable function.
    pat_solve=r'static int solve\(Sys\*s,const char\*solver,double tol,int maxit\).*?\nstatic void norms'
    if not re.search(pat_solve,s,flags=re.S):
        raise SystemExit(f'solve block not found in {c}')
    s=re.sub(pat_solve,lambda _: SOLVE.rstrip()+'\nstatic void norms',s,flags=re.S)

    c.write_text(s)


# Update generated correspondence notes so the package does not describe the
# superseded damped/line-searched benchmark RMT.
for d in ROOT.glob('Benchmark_*'):
    doc=d/'docs'/'PRODUCTION_SOLVER_CORRESPONDENCE.md'
    if doc.exists():
        q=doc.read_text()
        q=re.sub(
            r'RMT3H uses factor-three nine shifted coarse families.*?(?=\n\n|\Z)',
            'RMT3H uses factor-three nine shifted coarse families, no presmoothing, '
            'control-volume defect restriction, direct solution on every coarsest '
            'shifted grid, a single globally fixed 16-sweep postsmoothing count, '
            'and full index-space correction without interpolation, damping, '
            'line search, rollback, or hidden fallback.',
            q,flags=re.S
        )
        doc.write_text(q)

note=ROOT/'V2_METHOD_CHANGES.md'
if note.exists():
    q=note.read_text()
    q=re.sub(
        r'- RMT3H is recursively factor-three.*',
        '- RMT3H is a production-aligned factor-three multiple-shifted-grid '
        'benchmark adaptation: no presmoothing, control-volume defect restriction, '
        'direct coarsest solves, 16 post-smoothing sweeps, and full correction '
        'without interpolation, damping, line search, rollback, or hidden fallback.',
        q
    )
    note.write_text(q)
    for d in ROOT.glob('Benchmark_*'):
        p=d/'docs'/'V2_FIDELITY_NOTE.md'
        if p.exists(): p.write_text(q)

print('V3_RMT_PRODUCTION_ALIGNMENT_COMPLETE')
