#ifndef RMT_BOOK2D_IMPL_H
#define RMT_BOOK2D_IMPL_H

/*
 * 2-D Martynenko-style multiple-coarse-grid correction for the pressure block.
 *
 * Design goals:
 *   - factor-three coarsening in each coordinate direction;
 *   - independent maximum levels in x and y;
 *   - all 3^L shifted grid families represented by index mapping on finest-grid storage;
 *   - no presmoothing: sawtooth traversal from coarsest level to finest level;
 *   - coarse RHS obtained by control-volume-like box averaging of the current fine residual;
 *   - shifted-grid boundary stencils use the real distance from the mapped point to the
 *     physical boundary, rather than pretending every shifted-grid endpoint is a boundary;
 *   - homogeneous correction BCs matching the production pressure equation:
 *       x: Dirichlet correction = 0 at the two pressure/outlet faces,
 *       y: zero normal-gradient correction at top/bottom;
 *   - no recursive per-child calloc/free. A reusable workspace is allocated once per mesh.
 *
 * The pressure equation itself, physics, time stepping and stopping tolerance are not changed.
 */

typedef struct {
    int requestedLevels;   /* <=0 -> automatic, otherwise cap number of level indices */
    int smoothSweeps;      /* point-Seidel sweeps on non-coarsest levels */
    int coarseSweeps;      /* point-Seidel sweeps on the coarsest level */
    int nCoarsest;         /* target maximum points per 1-D shifted grid at deepest level */
} RMTBookConfig;

typedef struct {
    int levelXMax;
    int levelYMax;
    long long pointUpdates;
} RMTBookDiag;

typedef struct {
    int nx, ny;
    size_t n;
    double *remaining;
    double *rhsLevel;
    double *prefix;
    double *corr;
    double *res;
} RMTBookWorkspace;

static RMTBookWorkspace g_rmtBookWS = {0,0,0,NULL,NULL,NULL,NULL,NULL};
static long long g_rmtBookPointUpdatesTotal = 0;
static long long g_rmtBookFullAccepts = 0;
static long long g_rmtBookBacktracks = 0;
static long long g_rmtBookRejects = 0;
static double g_rmtBookLastAcceptedOmega = 0.0;
static int g_rmtBookLastLevelX = 0;
static int g_rmtBookLastLevelY = 0;

static void rmt_book_workspace_release(void) {
    free(g_rmtBookWS.remaining); g_rmtBookWS.remaining = NULL;
    free(g_rmtBookWS.rhsLevel); g_rmtBookWS.rhsLevel = NULL;
    free(g_rmtBookWS.prefix); g_rmtBookWS.prefix = NULL;
    free(g_rmtBookWS.corr); g_rmtBookWS.corr = NULL;
    free(g_rmtBookWS.res); g_rmtBookWS.res = NULL;
    g_rmtBookWS.nx = g_rmtBookWS.ny = 0;
    g_rmtBookWS.n = 0;
}

static RMTBookWorkspace *rmt_book_workspace(const Grid *g) {
    const size_t n = (size_t)(g->Nx+2) * (size_t)(g->Ny+2);
    const size_t np = (size_t)(g->Nx+1) * (size_t)(g->Ny+1);
    if (g_rmtBookWS.nx == g->Nx && g_rmtBookWS.ny == g->Ny && g_rmtBookWS.corr) {
        return &g_rmtBookWS;
    }
    rmt_book_workspace_release();
    g_rmtBookWS.nx = g->Nx;
    g_rmtBookWS.ny = g->Ny;
    g_rmtBookWS.n = n;
    g_rmtBookWS.remaining = (double*)calloc(n, sizeof(double));
    g_rmtBookWS.rhsLevel = (double*)calloc(n, sizeof(double));
    g_rmtBookWS.prefix = (double*)calloc(np, sizeof(double));
    g_rmtBookWS.corr = (double*)calloc(n, sizeof(double));
    g_rmtBookWS.res = (double*)calloc(n, sizeof(double));
    if (!g_rmtBookWS.remaining || !g_rmtBookWS.rhsLevel || !g_rmtBookWS.prefix ||
        !g_rmtBookWS.corr || !g_rmtBookWS.res) {
        die("alloc RMT book workspace");
    }
    return &g_rmtBookWS;
}

static int rmt_book_auto_level_1d(int n, int nCoarsest) {
    int level = 0;
    long long stride = 1;
    if (nCoarsest < 2) nCoarsest = 2;
    while ((n + (int)stride - 1)/(int)stride > nCoarsest) {
        if (stride > 100000000/3) break;
        stride *= 3;
        ++level;
        if (level >= 12) break;
    }
    return level;
}

static int rmt_book_ipow3(int level) {
    int s = 1;
    for (int k=0; k<level; ++k) s *= 3;
    return s;
}

/* A = -Laplacian on the finest cell-centred grid with homogeneous correction BCs. */
static void rmt_book_apply_fine_A(const Grid *g, const double *x, double *Ax) {
    const int nx=g->Nx, ny=g->Ny;
    const double ax=1.0/(g->dx*g->dx), ay=1.0/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double ap=0.0, sum=0.0;
        if (i > 1) { ap += ax; sum += ax*x[IDX(i-1,j,ny)]; }
        else       { ap += 2.0*ax; } /* cell centre is dx/2 from Dirichlet face */
        if (i < nx){ ap += ax; sum += ax*x[IDX(i+1,j,ny)]; }
        else       { ap += 2.0*ax; }
        if (j > 1) { ap += ay; sum += ay*x[IDX(i,j-1,ny)]; }
        /* missing south neighbour is homogeneous Neumann: no flux contribution */
        if (j < ny){ ap += ay; sum += ay*x[IDX(i,j+1,ny)]; }
        Ax[IDX(i,j,ny)] = ap*x[IDX(i,j,ny)] - sum;
    }
}

static void rmt_book_remaining_residual(const Grid *g, const double *b,
                                        const double *corr, double *remaining) {
    rmt_book_apply_fine_A(g, corr, remaining);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j)
        remaining[IDX(i,j,g->Ny)] = b[IDX(i,j,g->Ny)] - remaining[IDX(i,j,g->Ny)];
}

static void rmt_book_prefix_build(const Grid *g, const double *a, double *pfx) {
    const int nx=g->Nx, ny=g->Ny;
    const int pitch=ny+1;
    memset(pfx, 0, (size_t)(nx+1)*(size_t)(ny+1)*sizeof(double));
    for (int i=1; i<=nx; ++i) {
        double row=0.0;
        for (int j=1; j<=ny; ++j) {
            row += a[IDX(i,j,ny)];
            pfx[(size_t)i*pitch+j] = pfx[(size_t)(i-1)*pitch+j] + row;
        }
    }
}

static double rmt_book_box_average(const Grid *g, const double *pfx,
                                   int i, int j, int strideX, int strideY) {
    const int pitch=g->Ny+1;
    const int hx=(strideX-1)/2, hy=(strideY-1)/2;
    const int i1=MAX(1,i-hx), i2=MIN(g->Nx,i+hx);
    const int j1=MAX(1,j-hy), j2=MIN(g->Ny,j+hy);
    const double sum = pfx[(size_t)i2*pitch+j2]
                     - pfx[(size_t)(i1-1)*pitch+j2]
                     - pfx[(size_t)i2*pitch+(j1-1)]
                     + pfx[(size_t)(i1-1)*pitch+(j1-1)];
    const double cnt=(double)(i2-i1+1)*(double)(j2-j1+1);
    return cnt > 0.0 ? sum/cnt : 0.0;
}

static void rmt_book_restrict_level(const Grid *g, const double *remaining,
                                    double *rhsLevel, double *pfx,
                                    int strideX, int strideY) {
    rmt_book_prefix_build(g, remaining, pfx);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j)
        rhsLevel[IDX(i,j,g->Ny)] = rmt_book_box_average(g,pfx,i,j,strideX,strideY);
}

/*
 * Point-Seidel update on one mapped shifted grid.
 * The coarse-grid endpoint is NOT treated as the physical boundary unless it actually lies
 * there. The coefficient to a missing x neighbour is formed from the real distance from the
 * mapped cell centre to the physical Dirichlet face. For y, the physical BC is zero flux.
 */
static inline double rmt_book_update_point(const Grid *g, double *corr, const double *rhs,
                                           int i, int j, int strideX, int strideY) {
    const int nx=g->Nx, ny=g->Ny;
    const double Hx=(double)strideX*g->dx;
    const double Hy=(double)strideY*g->dy;
    const int hasW=(i-strideX >= 1), hasE=(i+strideX <= nx);
    const int hasS=(j-strideY >= 1), hasN=(j+strideY <= ny);

    const double dW=(i-0.5)*g->dx;
    const double dE=(nx-i+0.5)*g->dx;
    const double dS=(j-0.5)*g->dy;
    const double dN=(ny-j+0.5)*g->dy;

    const double xW=hasW ? 0.5*Hx : dW;
    const double xE=hasE ? 0.5*Hx : dE;
    const double yS=hasS ? 0.5*Hy : dS;
    const double yN=hasN ? 0.5*Hy : dN;
    const double Vx=MAX(xW+xE,1.0e-30);
    const double Vy=MAX(yS+yN,1.0e-30);

    double ap=0.0, sum=0.0;
    if (hasW) {
        const double a=1.0/(Hx*Vx); ap += a; sum += a*corr[IDX(i-strideX,j,ny)];
    } else {
        ap += 1.0/(MAX(dW,1.0e-30)*Vx); /* homogeneous Dirichlet correction at x=0 */
    }
    if (hasE) {
        const double a=1.0/(Hx*Vx); ap += a; sum += a*corr[IDX(i+strideX,j,ny)];
    } else {
        ap += 1.0/(MAX(dE,1.0e-30)*Vx); /* homogeneous Dirichlet correction at x=Lx */
    }
    if (hasS) {
        const double a=1.0/(Hy*Vy); ap += a; sum += a*corr[IDX(i,j-strideY,ny)];
    }
    /* no south neighbour -> homogeneous Neumann, zero flux */
    if (hasN) {
        const double a=1.0/(Hy*Vy); ap += a; sum += a*corr[IDX(i,j+strideY,ny)];
    }
    /* no north neighbour -> homogeneous Neumann, zero flux */

    return (rhs[IDX(i,j,ny)] + sum)/MAX(ap,1.0e-30);
}

static long long rmt_book_smooth_level(const Grid *g, double *corr, const double *rhs,
                                       int strideX, int strideY, int sweeps) {
    const int nGridX=MIN(strideX,g->Nx);
    const int nGridY=MIN(strideY,g->Ny);
    const int nGrids=nGridX*nGridY;
    long long updates=0;
    for (int sw=0; sw<sweeps; ++sw) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:updates) if(nGrids > 1)
#endif
        for (int gid=0; gid<nGrids; ++gid) {
            const int ox=1 + gid/nGridY;
            const int oy=1 + gid%nGridY;
            for (int i=ox; i<=g->Nx; i+=strideX) {
                for (int j=oy; j<=g->Ny; j+=strideY) {
                    corr[IDX(i,j,g->Ny)] = rmt_book_update_point(g,corr,rhs,i,j,strideX,strideY);
                    ++updates;
                }
            }
        }
    }
    return updates;
}

static void rmt_book_build_correction(const Grid *g, const double *b, double *corr,
                                      const RMTBookConfig *cfg, RMTBookDiag *diag) {
    RMTBookWorkspace *ws=rmt_book_workspace(g);
    int lx=rmt_book_auto_level_1d(g->Nx,cfg->nCoarsest);
    int ly=rmt_book_auto_level_1d(g->Ny,cfg->nCoarsest);
    if (cfg->requestedLevels > 0) {
        const int cap=MAX(0,cfg->requestedLevels-1);
        lx=MIN(lx,cap);
        ly=MIN(ly,cap);
    }
    const int lm=MAX(lx,ly);
    long long updates=0;
    memset(corr,0,ws->n*sizeof(double));

    /* Martynenko sawtooth direction: coarsest -> finest, no presmoothing. */
    for (int levelC=lm; levelC>=0; --levelC) {
        const int levelX=MIN(lx,levelC);
        const int levelY=MIN(ly,levelC);
        const int strideX=rmt_book_ipow3(levelX);
        const int strideY=rmt_book_ipow3(levelY);
        rmt_book_remaining_residual(g,b,corr,ws->remaining);
        rmt_book_restrict_level(g,ws->remaining,ws->rhsLevel,ws->prefix,strideX,strideY);
        const int sweeps=(levelC==lm) ? MAX(1,cfg->coarseSweeps) : MAX(1,cfg->smoothSweeps);
        updates += rmt_book_smooth_level(g,corr,ws->rhsLevel,strideX,strideY,sweeps);
    }

    if (diag) {
        diag->levelXMax=lx;
        diag->levelYMax=ly;
        diag->pointUpdates=updates;
    }
    g_rmtBookPointUpdatesTotal += updates;
    g_rmtBookLastLevelX=lx;
    g_rmtBookLastLevelY=ly;
}

#endif /* RMT_BOOK2D_IMPL_H */
