#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define IDX(i,j,ny) ((i)*((ny)+2)+(j))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
typedef struct { int Nx,Ny; double Lx,Ly,dx,dy,x1,x2; } Grid;
static void die(const char*m){fprintf(stderr,"FATAL: %s\n",m);exit(2);}
#include "../src/rmt_final2d_impl.h"

static double norm_inf(const Grid*g,const double*r){
    double m=0.0;
    for(int i=1;i<=g->Nx;++i)
        for(int j=1;j<=g->Ny;++j)
            m=MAX(m,fabs(r[IDX(i,j,g->Ny)]));
    return m;
}

static void residual(const Grid*g,const double*u,const double*b,double*Au,double*r){
    rmt_final_apply_fine_A(g,u,Au);
    for(int i=1;i<=g->Nx;++i)
        for(int j=1;j<=g->Ny;++j)
            r[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-Au[IDX(i,j,g->Ny)];
}

static int restriction_conservation_test(void){
    Grid g={23,17,1,1,1.0/23.0,1.0/17.0,0,0};
    size_t n=(size_t)(g.Nx+2)*(size_t)(g.Ny+2);
    double*r=calloc(n,sizeof(double));
    if(!r)die("alloc restriction test");
    RMTFinalWorkspace*w=rmt_final_workspace(&g);
    double fine=0.0;
    for(int i=1;i<=g.Nx;++i)for(int j=1;j<=g.Ny;++j){
        double v=0.31*i-0.27*j+0.013*i*j+sin(0.2*i+0.11*j);
        r[IDX(i,j,g.Ny)]=v;
        fine+=v;
    }

    int strides[]={3,9};
    for(int ss=0;ss<2;++ss){
        int sx=strides[ss], sy=strides[ss];
        rmt_final_restrict_level(&g,r,w->rhsLevel,w->prefix,sx,sy);
        for(int ox=1;ox<=MIN(sx,g.Nx);++ox)for(int oy=1;oy<=MIN(sy,g.Ny);++oy){
            double coarseIntegral=0.0;
            for(int i=ox;i<=g.Nx;i+=sx)for(int j=oy;j<=g.Ny;j+=sy){
                int i1,i2,j1,j2;
                rmt_final_cv_bounds_1d(g.Nx,i,sx,&i1,&i2);
                rmt_final_cv_bounds_1d(g.Ny,j,sy,&j1,&j2);
                double cells=(double)(i2-i1+1)*(double)(j2-j1+1);
                coarseIntegral+=w->rhsLevel[IDX(i,j,g.Ny)]*cells;
            }
            double scale=MAX(1.0,fabs(fine));
            if(fabs(coarseIntegral-fine)>5e-12*scale){
                fprintf(stderr,"RESTRICTION FAIL stride=%d offset=%d,%d fine=%.17g coarse=%.17g\n",
                        sx,ox,oy,fine,coarseIntegral);
                free(r); return 0;
            }
        }
    }
    free(r);
    puts("BOUNDARY_CV_RESTRICTION PASS");
    return 1;
}

static int solve_manufactured(int nx,int ny,int sweeps,int*cyclesOut,double*rhoOut){
    Grid g={nx,ny,0.02,0.02,0,0,0,0};
    g.dx=g.Lx/g.Nx; g.dy=g.Ly/g.Ny;
    size_t n=(size_t)(nx+2)*(size_t)(ny+2);
    double*exact=calloc(n,sizeof(double));
    double*b=calloc(n,sizeof(double));
    double*u=calloc(n,sizeof(double));
    double*Au=calloc(n,sizeof(double));
    double*r=calloc(n,sizeof(double));
    double*c=calloc(n,sizeof(double));
    if(!exact||!b||!u||!Au||!r||!c)die("alloc ladder");

    const double pi=acos(-1.0);
    for(int i=1;i<=nx;++i){
        double x=(i-.5)*g.dx;
        for(int j=1;j<=ny;++j){
            double y=(j-.5)*g.dy;
            exact[IDX(i,j,ny)]=sin(pi*x/g.Lx)*cos(pi*y/g.Ly);
        }
    }
    rmt_final_apply_fine_A(&g,exact,b);
    residual(&g,u,b,Au,r);
    double r0=norm_inf(&g,r), rr=r0;
    RMTFinalConfig cfg={0,sweeps,5};
    int cyc=0;
    for(;cyc<30 && rr/r0>1e-4;++cyc){
        residual(&g,u,b,Au,r);
        RMTFinalDiag d={0};
        rmt_final_build_correction(&g,r,c,&cfg,&d);
        for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j)
            u[IDX(i,j,ny)]+=c[IDX(i,j,ny)];
        residual(&g,u,b,Au,r);
        rr=norm_inf(&g,r);
        if(!isfinite(rr))break;
    }
    double rel=rr/r0;
    double rho=(cyc>0 && rel>0)?pow(rel,1.0/cyc):rel;
    printf("MESH_LADDER Nx=%d Ny=%d levels=%d/%d sweeps=%d cycles=%d relative=%.12e rho=%.6f\n",
           nx,ny,g_rmtFinalLastLevelX,g_rmtFinalLastLevelY,sweeps,cyc,rel,rho);
    if(cyclesOut)*cyclesOut=cyc;
    if(rhoOut)*rhoOut=rho;
    free(exact);free(b);free(u);free(Au);free(r);free(c);
    return isfinite(rel) && rel<=1.0e-4;
}

int main(void){
    int sweeps=4;
    const char*e=getenv("RMT_TEST_SWEEPS");
    if(e&&*e)sweeps=atoi(e);
    if(sweeps<1)sweeps=1;
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif

    if(!restriction_conservation_test())return 3;

    const int meshes[4][2]={{108,36},{216,72},{432,144},{864,288}};
    int cycles[4]={0};
    double rho[4]={0};
    int pass=1;
    for(int k=0;k<4;++k)
        if(!solve_manufactured(meshes[k][0],meshes[k][1],sweeps,&cycles[k],&rho[k]))
            pass=0;

    int cmin=cycles[0],cmax=cycles[0];
    for(int k=1;k<4;++k){cmin=MIN(cmin,cycles[k]);cmax=MAX(cmax,cycles[k]);}
    printf("MESH_INDEPENDENCE cycles_min=%d cycles_max=%d spread=%d\n",cmin,cmax,cmax-cmin);
    if(cmax>20 || cmax-cmin>6)pass=0;

    rmt_final_workspace_release();
    if(!pass){fprintf(stderr,"SELFTEST FAIL\n");return 4;}
    puts("SELFTEST PASS");
    return 0;
}
