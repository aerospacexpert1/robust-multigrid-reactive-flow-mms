#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define IDX(i,j,ny) ((i)*((ny)+2)+(j))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
typedef struct { int Nx,Ny; double Lx,Ly,dx,dy,x1,x2; } Grid;
static void die(const char*m){fprintf(stderr,"FATAL: %s\n",m);exit(2);}

#include "reference_rmt_final2d_impl.h"
#include "../src/rmt_finalvol2_impl.h"

static int one(int nx,int ny){
    Grid g={nx,ny,0.02,0.02,0,0,0,0};
    g.dx=g.Lx/g.Nx;g.dy=g.Ly/g.Ny;
    const size_t n=(size_t)(nx+2)*(size_t)(ny+2);
    double*d=calloc(n,sizeof(double)),*a=calloc(n,sizeof(double)),*b=calloc(n,sizeof(double));
    if(!d||!a||!b)die("alloc equivalence");
    for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j)
        d[IDX(i,j,ny)]=sin(0.071*i)+cos(0.113*j)+0.05*sin(0.013*i*j);

    RMTFinalConfig c0={0,16,0,5};
    RMTFinalDiag q0={0};
    rmt_final_build_correction(&g,d,a,&c0,&q0);

    RMTV2Config c1={0,16,5};
    RMTV2Diag q1={0};
    rmt_v2_build_correction(&g,d,b,&c1,&q1);

    double md=0.0,mr=0.0;
    for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j){
        md=MAX(md,fabs(a[IDX(i,j,ny)]-b[IDX(i,j,ny)]));
        mr=MAX(mr,fabs(a[IDX(i,j,ny)]));
    }
    const double rel=md/MAX(mr,1e-300);
    printf("OPT_EQ Nx=%d Ny=%d old_updates=%lld new_updates=%lld rel_max_diff=%.12e\n",
           nx,ny,q0.pointUpdates,q1.pointUpdates,rel);
    free(d);free(a);free(b);
    return isfinite(rel)&&rel<5e-11&&q0.pointUpdates==q1.pointUpdates;
}

int main(void){
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    int ok=1;
    ok&=one(108,36);
    ok&=one(864,288);
    rmt_final_workspace_release();
    rmt_v2_workspace_release();
    if(!ok){fprintf(stderr,"OPTIMIZATION_EQUIVALENCE FAIL\n");return 5;}
    puts("OPTIMIZATION_EQUIVALENCE PASS");
    return 0;
}
