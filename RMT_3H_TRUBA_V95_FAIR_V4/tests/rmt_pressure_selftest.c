#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define IDX(i,j,ny) ((i)*((ny)+2)+(j))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
typedef struct { int Nx,Ny; double Lx,Ly,dx,dy,x1,x2; } Grid;
static void die(const char*m){fprintf(stderr,"FATAL: %s\n",m);exit(2);}
#include "../src/rmt_book2d_impl.h"

static double norm_inf(const Grid*g,const double*r){
    double m=0.0;
    for(int i=1;i<=g->Nx;++i)
        for(int j=1;j<=g->Ny;++j)
            m=MAX(m,fabs(r[IDX(i,j,g->Ny)]));
    return m;
}

static void residual(const Grid*g,const double*u,const double*b,double*Au,double*r){
    rmt_book_apply_fine_A(g,u,Au);
    for(int i=1;i<=g->Nx;++i)
        for(int j=1;j<=g->Ny;++j)
            r[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-Au[IDX(i,j,g->Ny)];
}

static void add_scaled(const Grid*g,double*u,const double*c,double omega){
    for(int i=1;i<=g->Nx;++i)
        for(int j=1;j<=g->Ny;++j)
            u[IDX(i,j,g->Ny)] += omega*c[IDX(i,j,g->Ny)];
}

int main(void){
    Grid g={81,27,1,1,1.0/81.0,1.0/27.0,0,0};
    size_t n=(size_t)(g.Nx+2)*(size_t)(g.Ny+2);
    double*b=calloc(n,sizeof(double));
    double*u=calloc(n,sizeof(double));
    double*Au=calloc(n,sizeof(double));
    double*r=calloc(n,sizeof(double));
    double*c=calloc(n,sizeof(double));
    if(!b||!u||!Au||!r||!c)die("alloc");

    double pi=acos(-1.0);
    for(int i=1;i<=g.Nx;++i){
        double x=(i-.5)*g.dx;
        for(int j=1;j<=g.Ny;++j){
            double y=(j-.5)*g.dy;
            b[IDX(i,j,g.Ny)]=2*pi*pi*sin(pi*x)*cos(pi*y);
        }
    }

    residual(&g,u,b,Au,r);
    double r0=norm_inf(&g,r);
    printf("DIAG initial %.12e\n",r0);

    /* Level isolation: these are raw, full-strength corrections. */
    int caps[]={1,2,3,0};
    double level_ratio[4]={0};
    for(int q=0;q<4;++q){
        memset(u,0,n*sizeof(double));
        residual(&g,u,b,Au,r);
        RMTBookConfig cfg={caps[q],8,16,5};
        RMTBookDiag d={0};
        rmt_book_build_correction(&g,r,c,&cfg,&d);
        add_scaled(&g,u,c,1.0);
        residual(&g,u,b,Au,r);
        double rr=norm_inf(&g,r);
        level_ratio[q]=rr/r0;
        printf("LEVEL_DIAG cap=%d levels=%d/%d ratio=%.12e updates=%lld\n",
               caps[q],d.levelXMax,d.levelYMax,level_ratio[q],d.pointUpdates);
    }

    /* RAW_RMT is deliberately reported without line search.  This prevents the
       safeguard from hiding the intrinsic full-strength RMT behaviour. */
    memset(u,0,n*sizeof(double));
    double prev=r0;
    int raw_decreases=0;
    RMTBookConfig cfg={0,8,16,5};
    for(int q=0;q<8;++q){
        residual(&g,u,b,Au,r);
        RMTBookDiag d={0};
        rmt_book_build_correction(&g,r,c,&cfg,&d);
        add_scaled(&g,u,c,1.0);
        residual(&g,u,b,Au,r);
        double rr=norm_inf(&g,r);
        printf("RAW_RMT q=%d ratio_step=%.12e relative=%.12e\n",q,rr/prev,rr/r0);
        if(rr<prev)++raw_decreases;
        prev=rr;
    }
    double raw_final=prev/r0;
    printf("RAW_RMT_RESULT relative=%.12e decreases=%d/8\n",raw_final,raw_decreases);

    /* Production acceptance policy: omega=1 is always tried first.  Halving is
       allowed only when the full correction does not reduce the exact fine-grid
       residual.  Every backtrack is printed and counted, so this cannot silently
       improve the thesis timing/convergence comparison. */
    memset(u,0,n*sizeof(double));
    prev=r0;
    int full_accepts=0, backtracked_cycles=0, total_halvings=0, rejects=0;
    for(int q=0;q<8;++q){
        residual(&g,u,b,Au,r);
        double oldr=norm_inf(&g,r);
        RMTBookDiag d={0};
        rmt_book_build_correction(&g,r,c,&cfg,&d);

        double omega=1.0;
        int accepted=0, halvings=0;
        double rr=oldr;
        for(int ls=0;ls<8;++ls){
            add_scaled(&g,u,c,omega);
            residual(&g,u,b,Au,r);
            rr=norm_inf(&g,r);
            if(isfinite(rr) && rr<oldr){
                accepted=1;
                halvings=ls;
                break;
            }
            add_scaled(&g,u,c,-omega);
            omega*=0.5;
        }

        if(!accepted){
            ++rejects;
            omega=0.0;
            rr=oldr;
        }else if(halvings==0){
            ++full_accepts;
        }else{
            ++backtracked_cycles;
            total_halvings+=halvings;
        }

        printf("SAFEGUARDED_RMT q=%d omega=%.12e halvings=%d ratio_step=%.12e relative=%.12e\n",
               q,omega,halvings,rr/prev,rr/r0);
        prev=rr;
    }

    double safe_final=prev/r0;
    printf("SAFEGUARDED_RMT_RESULT relative=%.12e full_accepts=%d/8 backtracked_cycles=%d total_halvings=%d rejects=%d\n",
           safe_final,full_accepts,backtracked_cycles,total_halvings,rejects);

    /* Acceptance is mathematical, not cosmetic:
       1) finest-only GS must reduce the exact production-consistent residual;
       2) one full raw RMT correction must reduce it;
       3) eight raw full-strength RMT iterations must give net convergence;
       4) the explicit production safeguard must never reject a correction and
          must reduce the residual by at least one order of magnitude. */
    int pass=1;
    if(!(level_ratio[0] < 1.0)) pass=0;
    if(!(level_ratio[3] < 1.0)) pass=0;
    if(!(raw_final < 0.20)) pass=0;
    if(rejects != 0) pass=0;
    if(!(safe_final < 0.10)) pass=0;

    free(b);free(u);free(Au);free(r);free(c);
    rmt_book_workspace_release();
    if(!pass){fprintf(stderr,"SELFTEST FAIL\n");return 3;}
    puts("SELFTEST PASS");
    return 0;
}
