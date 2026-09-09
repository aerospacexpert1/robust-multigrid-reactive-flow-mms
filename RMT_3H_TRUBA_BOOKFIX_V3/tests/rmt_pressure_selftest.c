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

typedef struct {
    int Nx, Ny;
    double Lx, Ly;
    double dx, dy;
    double x1, x2;
} Grid;

static void die(const char *m) { fprintf(stderr,"FATAL: %s\n",m); exit(2); }

#include "../src/rmt_book2d_impl.h"

static double max_norm(const Grid *g, const double *r) {
    double m=0.0;
    for (int i=1;i<=g->Nx;++i) for (int j=1;j<=g->Ny;++j) {
        double a=fabs(r[IDX(i,j,g->Ny)]);
        if (a>m) m=a;
    }
    return m;
}

int main(void) {
    Grid g={0};
    g.Nx=81; g.Ny=27; g.Lx=1.0; g.Ly=1.0;
    g.dx=g.Lx/g.Nx; g.dy=g.Ly/g.Ny;
    size_t n=(size_t)(g.Nx+2)*(size_t)(g.Ny+2);
    double *u=(double*)calloc(n,sizeof(double));
    double *b=(double*)calloc(n,sizeof(double));
    double *Au=(double*)calloc(n,sizeof(double));
    double *r=(double*)calloc(n,sizeof(double));
    double *c=(double*)calloc(n,sizeof(double));
    if(!u||!b||!Au||!r||!c) die("selftest allocation");

    const double pi=acos(-1.0);
    for(int i=1;i<=g.Nx;++i) {
        double x=(i-0.5)*g.dx;
        for(int j=1;j<=g.Ny;++j) {
            double y=(j-0.5)*g.dy;
            /* -Laplace[sin(pi x) cos(pi y)] = 2*pi^2*u;
               x is homogeneous Dirichlet, y is homogeneous Neumann. */
            b[IDX(i,j,g.Ny)]=2.0*pi*pi*sin(pi*x)*cos(pi*y);
        }
    }

    RMTBookConfig cfg={0,8,16,5};
    double r0=0.0, prev=0.0;
    int full_decrease_count=0;
    for(int q=0;q<8;++q) {
        rmt_book_apply_fine_A(&g,u,Au);
        for(int i=1;i<=g.Nx;++i) for(int j=1;j<=g.Ny;++j)
            r[IDX(i,j,g.Ny)]=b[IDX(i,j,g.Ny)]-Au[IDX(i,j,g.Ny)];
        double rn=max_norm(&g,r);
        if(q==0){r0=rn; prev=rn;}
        printf("SELFTEST before q=%d residual=%.12e relative=%.12e\n",q,rn,rn/r0);
        RMTBookDiag d={0,0,0};
        rmt_book_build_correction(&g,r,c,&cfg,&d);
        for(int i=1;i<=g.Nx;++i) for(int j=1;j<=g.Ny;++j)
            u[IDX(i,j,g.Ny)] += c[IDX(i,j,g.Ny)];
        rmt_book_apply_fine_A(&g,u,Au);
        for(int i=1;i<=g.Nx;++i) for(int j=1;j<=g.Ny;++j)
            r[IDX(i,j,g.Ny)]=b[IDX(i,j,g.Ny)]-Au[IDX(i,j,g.Ny)];
        double after=max_norm(&g,r);
        printf("SELFTEST after  q=%d residual=%.12e ratio=%.6f levels=%d/%d updates=%lld\n",
               q,after,after/MAX(rn,1e-300),d.levelXMax,d.levelYMax,d.pointUpdates);
        if(after < rn) ++full_decrease_count;
        prev=after;
    }

    rmt_book_apply_fine_A(&g,u,Au);
    for(int i=1;i<=g.Nx;++i) for(int j=1;j<=g.Ny;++j)
        r[IDX(i,j,g.Ny)]=b[IDX(i,j,g.Ny)]-Au[IDX(i,j,g.Ny)];
    double rf=max_norm(&g,r);
    double rel=rf/MAX(r0,1e-300);
    printf("SELFTEST_RESULT initial=%.12e final=%.12e relative=%.12e full_decreases=%d/8\n",
           r0,rf,rel,full_decrease_count);

    free(u);free(b);free(Au);free(r);free(c);
    rmt_book_workspace_release();

    if(full_decrease_count < 7) {
        fprintf(stderr,"SELFTEST FAIL: full-strength RMT correction was not consistently decreasing residual\n");
        return 3;
    }
    if(!(rel < 2.0e-1)) {
        fprintf(stderr,"SELFTEST FAIL: insufficient residual reduction\n");
        return 4;
    }
    puts("SELFTEST PASS");
    return 0;
}
