#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define IDX(i,j,ny) ((i)*((ny)+2)+(j))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
typedef struct { int Nx,Ny; double Lx,Ly,dx,dy,x1,x2; } Grid;
static void die(const char*m){fprintf(stderr,"FATAL: %s\n",m);exit(2);}
#include "../src/rmt_finalvol2_impl.h"

static double wall_now(void){
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;
#endif
}
static double norm_inf(const Grid*g,const double*r){
    double m=0.0;
    for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)m=MAX(m,fabs(r[IDX(i,j,g->Ny)]));
    return m;
}
static void residual(const Grid*g,const double*u,const double*b,double*Au,double*r){
    rmt_v2_apply_fine_A(g,u,Au);
    for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)
        r[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-Au[IDX(i,j,g->Ny)];
}
static double exact_mode(int mode,double x,double y,double Lx,double Ly){
    const double pi=acos(-1.0);
    if(mode==1) return sin(pi*x/Lx)*cos(pi*y/Ly);
    if(mode==2) return sin(7*pi*x/Lx)*cos(5*pi*y/Ly);
    return 0.70*sin(pi*x/Lx)*cos(pi*y/Ly)
         + 0.20*sin(7*pi*x/Lx)*cos(5*pi*y/Ly)
         + 0.10*sin(13*pi*x/Lx)*cos(4*pi*y/Ly);
}
static const char*mode_name(int mode){return mode==1?"LOW":mode==2?"MIXED":"COMBINED";}

static int run_case(int nx,int ny,int mode,int sweeps){
    Grid g={nx,ny,0.02,0.02,0,0,0,0};
    g.dx=g.Lx/g.Nx;g.dy=g.Ly/g.Ny;
    const size_t n=(size_t)(nx+2)*(size_t)(ny+2);
    double*exact=calloc(n,sizeof(double)),*b=calloc(n,sizeof(double)),*u=calloc(n,sizeof(double));
    double*Au=calloc(n,sizeof(double)),*r=calloc(n,sizeof(double)),*c=calloc(n,sizeof(double));
    if(!exact||!b||!u||!Au||!r||!c)die("alloc calibration");
    for(int i=1;i<=nx;++i){
        const double x=(i-.5)*g.dx;
        for(int j=1;j<=ny;++j){
            const double y=(j-.5)*g.dy;
            exact[IDX(i,j,ny)]=exact_mode(mode,x,y,g.Lx,g.Ly);
        }
    }
    rmt_v2_apply_fine_A(&g,exact,b);
    residual(&g,u,b,Au,r);
    const double r0=norm_inf(&g,r);
    double rr=r0;
    RMTV2Config cfg={0,sweeps,5};
    const long long upd0=g_rmtV2PointUpdatesTotal;
    const long long lu0=g_rmtV2LUFactorizationsTotal;
    int cyc=0;
    const double t0=wall_now();
    for(;cyc<30 && rr/MAX(r0,1e-300)>1e-4;++cyc){
        residual(&g,u,b,Au,r);
        RMTV2Diag d={0};
        rmt_v2_build_correction(&g,r,c,&cfg,&d);
        for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j)u[IDX(i,j,ny)]+=c[IDX(i,j,ny)];
        residual(&g,u,b,Au,r);
        rr=norm_inf(&g,r);
        if(!isfinite(rr))break;
    }
    const double elapsed=wall_now()-t0;
    const double rel=rr/MAX(r0,1e-300);
    const double rho=(cyc>0&&rel>0)?pow(rel,1.0/(double)cyc):rel;
    const long long updates=g_rmtV2PointUpdatesTotal-upd0;
    const long long lufact=g_rmtV2LUFactorizationsTotal-lu0;
    const int pass=isfinite(rel)&&rel<=1e-4&&cyc<=30;
    printf("CAL_RESULT sweeps=%d mode=%s Nx=%d Ny=%d cycles=%d relative=%.12e rho=%.9f "
           "point_updates=%lld lu_factorizations=%lld elapsed_s=%.9f pass=%d\n",
           sweeps,mode_name(mode),nx,ny,cyc,rel,rho,updates,lufact,elapsed,pass);
    free(exact);free(b);free(u);free(Au);free(r);free(c);
    return pass;
}
int main(int argc,char**argv){
    if(argc!=2){fprintf(stderr,"usage: %s SWEEPS\n",argv[0]);return 2;}
    const int sweeps=atoi(argv[1]);
    if(sweeps<1||sweeps>64){fprintf(stderr,"invalid sweeps\n");return 2;}
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    const int meshes[4][2]={{108,36},{216,72},{432,144},{864,288}};
    int ok=1;
    for(int mode=1;mode<=3;++mode)
        for(int k=0;k<4;++k)
            if(!run_case(meshes[k][0],meshes[k][1],mode,sweeps))ok=0;
    rmt_v2_workspace_release();
    printf("CAL_CANDIDATE sweeps=%d pass=%d\n",sweeps,ok);
    return ok?0:4;
}
