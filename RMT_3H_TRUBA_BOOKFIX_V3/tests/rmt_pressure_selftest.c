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
static double norm(const Grid*g,const double*r){double m=0;for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)m=MAX(m,fabs(r[IDX(i,j,g->Ny)]));return m;}
static void residual(const Grid*g,const double*u,const double*b,double*Au,double*r){rmt_book_apply_fine_A(g,u,Au);for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)r[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-Au[IDX(i,j,g->Ny)];}
int main(void){
 Grid g={81,27,1,1,1.0/81.0,1.0/27.0,0,0};size_t n=(size_t)(g.Nx+2)*(g.Ny+2);
 double*b=calloc(n,sizeof(double)),*u=calloc(n,sizeof(double)),*Au=calloc(n,sizeof(double)),*r=calloc(n,sizeof(double)),*c=calloc(n,sizeof(double));if(!b||!u||!Au||!r||!c)die("alloc");
 double pi=acos(-1.0);for(int i=1;i<=g.Nx;++i){double x=(i-.5)*g.dx;for(int j=1;j<=g.Ny;++j){double y=(j-.5)*g.dy;b[IDX(i,j,g.Ny)]=2*pi*pi*sin(pi*x)*cos(pi*y);}}
 residual(&g,u,b,Au,r);double r0=norm(&g,r);printf("DIAG initial %.12e\n",r0);
 int caps[]={1,2,3,0};
 for(int q=0;q<4;++q){memset(u,0,n*sizeof(double));residual(&g,u,b,Au,r);RMTBookConfig cfg={caps[q],8,16,5};RMTBookDiag d={0};rmt_book_build_correction(&g,r,c,&cfg,&d);for(int i=1;i<=g.Nx;++i)for(int j=1;j<=g.Ny;++j)u[IDX(i,j,g.Ny)]+=c[IDX(i,j,g.Ny)];residual(&g,u,b,Au,r);double rr=norm(&g,r);printf("LEVEL_DIAG cap=%d levels=%d/%d ratio=%.12e updates=%lld\n",caps[q],d.levelXMax,d.levelYMax,rr/r0,d.pointUpdates);}
 /* Acceptance gate: finest-only GS must work. The full RMT is deliberately required to
    decrease residual too; otherwise packaging stops. */
 memset(u,0,n*sizeof(double));double prev=r0;int dec=0;RMTBookConfig cfg={0,8,16,5};for(int q=0;q<8;++q){residual(&g,u,b,Au,r);RMTBookDiag d={0};rmt_book_build_correction(&g,r,c,&cfg,&d);for(int i=1;i<=g.Nx;++i)for(int j=1;j<=g.Ny;++j)u[IDX(i,j,g.Ny)]+=c[IDX(i,j,g.Ny)];residual(&g,u,b,Au,r);double rr=norm(&g,r);printf("SELFTEST q=%d ratio_step=%.12e relative=%.12e\n",q,rr/prev,rr/r0);if(rr<prev)++dec;prev=rr;}
 printf("SELFTEST_RESULT relative=%.12e decreases=%d/8\n",prev/r0,dec);free(b);free(u);free(Au);free(r);free(c);rmt_book_workspace_release();if(dec<7||prev/r0>=0.2){fprintf(stderr,"SELFTEST FAIL\n");return 3;}puts("SELFTEST PASS");return 0;
}
