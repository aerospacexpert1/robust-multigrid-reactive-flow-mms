#ifndef RMT_FINAL2D_IMPL_H
#define RMT_FINAL2D_IMPL_H

/*
 * RMT FINAL: Martynenko-style multiple-coarse-grid correction for the
 * cell-centred V95 pressure equation.
 *
 * Key invariants:
 *  - factor-three multiple coarse grids; the union of all shifted grids on a
 *    level is exactly the finest-grid unknown set;
 *  - sawtooth schedule, coarsest -> finest, no presmoothing;
 *  - correction remains stored on the finest-grid array, so transfer between
 *    levels is an index-remapping operation and introduces no interpolation;
 *  - coarse RHS is the control-volume average of the exact fine-grid defect;
 *  - the control volume used for restriction is exactly the same control
 *    volume used by the shifted coarse-grid operator, including boundary grids;
 *  - coarsest shifted-grid systems are solved directly;
 *  - finer levels use parallel red-black Gauss-Seidel smoothing;
 *  - the resulting RMT correction is added in full.  No residual-monotonic
 *    line search or hidden fallback is part of the algorithm.
 */

typedef struct {
    int requestedLevels;   /* 0 = automatic independent x/y hierarchy */
    int smoothSweeps;      /* postsmoothing sweeps on every noncoarsest level */
    int nCoarsest;         /* target max points/direction on coarsest grids */
} RMTFinalConfig;

typedef struct {
    int levelXMax, levelYMax;
    long long pointUpdates;
    long long directSolves;
    long long directUnknowns;
} RMTFinalDiag;

typedef struct {
    int nx, ny;
    size_t n;
    double *rhsLevel;
    double *prefix;
    double *corr;
    double *res;
} RMTFinalWorkspace;

static RMTFinalWorkspace g_rmtFinalWS = {0,0,0,NULL,NULL,NULL,NULL};
static long long g_rmtFinalPointUpdatesTotal = 0;
static long long g_rmtFinalDirectSolvesTotal = 0;
static long long g_rmtFinalDirectUnknownsTotal = 0;
static long long g_rmtFinalNonmonotoneCycles = 0;
static double g_rmtFinalLastCycleRatio = 0.0;
static double g_rmtFinalMaxCycleRatio = 0.0;
static int g_rmtFinalLastLevelX = 0;
static int g_rmtFinalLastLevelY = 0;

static void rmt_final_workspace_release(void) {
    free(g_rmtFinalWS.rhsLevel);
    free(g_rmtFinalWS.prefix);
    free(g_rmtFinalWS.corr);
    free(g_rmtFinalWS.res);
    memset(&g_rmtFinalWS, 0, sizeof(g_rmtFinalWS));
}

static RMTFinalWorkspace *rmt_final_workspace(const Grid *g) {
    const size_t n = (size_t)(g->Nx+2)*(size_t)(g->Ny+2);
    const size_t np = (size_t)(g->Nx+1)*(size_t)(g->Ny+1);
    if (g_rmtFinalWS.nx == g->Nx && g_rmtFinalWS.ny == g->Ny && g_rmtFinalWS.corr)
        return &g_rmtFinalWS;

    rmt_final_workspace_release();
    g_rmtFinalWS.nx = g->Nx;
    g_rmtFinalWS.ny = g->Ny;
    g_rmtFinalWS.n = n;
    g_rmtFinalWS.rhsLevel = calloc(n, sizeof(double));
    g_rmtFinalWS.prefix = calloc(np, sizeof(double));
    g_rmtFinalWS.corr = calloc(n, sizeof(double));
    g_rmtFinalWS.res = calloc(n, sizeof(double));
    if (!g_rmtFinalWS.rhsLevel || !g_rmtFinalWS.prefix ||
        !g_rmtFinalWS.corr || !g_rmtFinalWS.res)
        die("alloc RMT FINAL workspace");
    return &g_rmtFinalWS;
}

static int rmt_final_auto_level_1d(int n, int nc) {
    int l = 0, s = 1;
    if (nc < 2) nc = 2;
    while ((n+s-1)/s > nc) {
        if (s > INT_MAX/3) break;
        s *= 3;
        ++l;
        if (l >= 12) break;
    }
    return l;
}

static int rmt_final_ipow3(int l) {
    int s = 1;
    while (l-- > 0) s *= 3;
    return s;
}

/* Fine-grid algebraic operator A=-Laplace with the exact V95 correction BCs:
 * homogeneous Dirichlet in x and homogeneous Neumann in y. */
static void rmt_final_apply_fine_A(const Grid *g, const double *x, double *Ax) {
    const int nx = g->Nx, ny = g->Ny;
    const double ax = 1.0/(g->dx*g->dx);
    const double ay = 1.0/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double ap = 0.0, sum = 0.0;
        if (i > 1) {
            ap += ax;
            sum += ax*x[IDX(i-1,j,ny)];
        } else {
            ap += 2.0*ax;
        }
        if (i < nx) {
            ap += ax;
            sum += ax*x[IDX(i+1,j,ny)];
        } else {
            ap += 2.0*ax;
        }
        if (j > 1) {
            ap += ay;
            sum += ay*x[IDX(i,j-1,ny)];
        }
        if (j < ny) {
            ap += ay;
            sum += ay*x[IDX(i,j+1,ny)];
        }
        Ax[IDX(i,j,ny)] = ap*x[IDX(i,j,ny)] - sum;
    }
}

static void rmt_final_prefix_build(const Grid *g, const double *a, double *p) {
    const int ny = g->Ny, pitch = ny+1;
    memset(p, 0, (size_t)(g->Nx+1)*(size_t)(ny+1)*sizeof(double));
    for (int i=1; i<=g->Nx; ++i) {
        double row = 0.0;
        for (int j=1; j<=ny; ++j) {
            row += a[IDX(i,j,ny)];
            p[(size_t)i*pitch+j] = p[(size_t)(i-1)*pitch+j] + row;
        }
    }
}

/*
 * Exact control-volume bounds for one shifted factor-three grid.
 *
 * This is the boundary bug fixed relative to V4.  A boundary coarse cell is
 * NOT generally a symmetric stride-wide box around its representative point.
 * Example: stride=3, first representative i=3 has centre x=2.5h and its first
 * east neighbour at 5.5h, hence the coarse CV is [0,4h] and contains fine
 * cells 1..4, not 2..4.  The same geometry is used below by the coarse operator.
 */
static inline void rmt_final_cv_bounds_1d(int n, int idx, int stride,
                                          int *lo, int *hi) {
    const int half = (stride-1)/2;
    const int hasMinus = (idx-stride >= 1);
    const int hasPlus  = (idx+stride <= n);
    *lo = hasMinus ? idx-half : 1;
    *hi = hasPlus  ? idx+half : n;
}

static double rmt_final_cv_average(const Grid *g, const double *p,
                                   int i, int j, int sx, int sy) {
    const int pitch = g->Ny+1;
    int i1, i2, j1, j2;
    rmt_final_cv_bounds_1d(g->Nx, i, sx, &i1, &i2);
    rmt_final_cv_bounds_1d(g->Ny, j, sy, &j1, &j2);
    const double sum = p[(size_t)i2*pitch+j2]
                     - p[(size_t)(i1-1)*pitch+j2]
                     - p[(size_t)i2*pitch+j1-1]
                     + p[(size_t)(i1-1)*pitch+j1-1];
    const double cnt = (double)(i2-i1+1)*(double)(j2-j1+1);
    return cnt > 0.0 ? sum/cnt : 0.0;
}

static void rmt_final_restrict_level(const Grid *g, const double *r,
                                     double *rhs, double *prefix,
                                     int sx, int sy) {
    rmt_final_prefix_build(g, r, prefix);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=g->Nx; ++i)
        for (int j=1; j<=g->Ny; ++j)
            rhs[IDX(i,j,g->Ny)] = rmt_final_cv_average(g,prefix,i,j,sx,sy);
}

/*
 * One-dimensional finite-volume contribution on a shifted cell-centred grid.
 * The physical boundary is a control-volume face whenever a same-residue
 * neighbour is missing.  V is therefore the true boundary-to-midpoint width,
 * exactly matching rmt_final_cv_bounds_1d().
 */
static inline void rmt_final_dim_coeff(int n, int idx, int stride, double h,
                                       int neumann,
                                       double *ap, double *am, double *az) {
    const int hm = (idx-stride >= 1);
    const int hp = (idx+stride <= n);
    const double H = stride*h;
    *am = 0.0;
    *az = 0.0;

    if (hm && hp) {
        const double a = 1.0/(H*H);
        *ap += 2.0*a;
        *am = a;
        *az = a;
        return;
    }

    if (!hm && hp) {
        const double d = (idx-0.5)*h;
        const double V = d + 0.5*H;
        const double aNbr = 1.0/(H*V);
        const double aB = neumann ? 0.0 : 1.0/(d*V);
        *ap += aNbr + aB;
        *az = aNbr;
        return;
    }

    if (hm && !hp) {
        const double d = (n-idx+0.5)*h;
        const double V = d + 0.5*H;
        const double aNbr = 1.0/(H*V);
        const double aB = neumann ? 0.0 : 1.0/(d*V);
        *ap += aNbr + aB;
        *am = aNbr;
        return;
    }

    /* One point in this shifted direction: the CV spans the whole domain. */
    {
        const double V = n*h;
        if (!neumann) {
            const double d0 = (idx-0.5)*h;
            const double d1 = (n-idx+0.5)*h;
            *ap += 1.0/(d0*V) + 1.0/(d1*V);
        }
    }
}

static inline double rmt_final_update_point(const Grid *g, double *c,
                                            const double *rhs,
                                            int i, int j, int sx, int sy) {
    const int ny = g->Ny;
    const int hw = (i-sx >= 1), he = (i+sx <= g->Nx);
    const int hs = (j-sy >= 1), hn = (j+sy <= g->Ny);
    double ap=0.0, aw=0.0, ae=0.0, as=0.0, an=0.0;

    rmt_final_dim_coeff(g->Nx,i,sx,g->dx,0,&ap,&aw,&ae);
    rmt_final_dim_coeff(g->Ny,j,sy,g->dy,1,&ap,&as,&an);

    double sum = 0.0;
    if (hw) sum += aw*c[IDX(i-sx,j,ny)];
    if (he) sum += ae*c[IDX(i+sx,j,ny)];
    if (hs) sum += as*c[IDX(i,j-sy,ny)];
    if (hn) sum += an*c[IDX(i,j+sy,ny)];
    return (rhs[IDX(i,j,ny)] + sum)/MAX(ap,1.0e-300);
}

/*
 * Parallel RBGS on all independent shifted grids of a level simultaneously.
 * Coarse-grid neighbours differ by +/-stride and therefore always have the
 * opposite checkerboard colour.  For stride=1 this is ordinary fine-grid RBGS.
 */
static long long rmt_final_smooth_level(const Grid *g, double *c,
                                        const double *rhs,
                                        int sx, int sy, int sweeps) {
    long long updates = 0;
    const long long nxy = (long long)g->Nx*(long long)g->Ny;
    for (int sw=0; sw<sweeps; ++sw) {
        for (int color=0; color<2; ++color) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) reduction(+:updates) if(nxy>2048)
#endif
            for (int i=1; i<=g->Nx; ++i) for (int j=1; j<=g->Ny; ++j) {
                const int ox = (i-1)%sx + 1;
                const int oy = (j-1)%sy + 1;
                const int qi = (i-ox)/sx;
                const int qj = (j-oy)/sy;
                if (((qi+qj)&1) != color) continue;
                c[IDX(i,j,g->Ny)] = rmt_final_update_point(g,c,rhs,i,j,sx,sy);
                ++updates;
            }
        }
    }
    return updates;
}

#define RMT_FINAL_DIRECT_MAX 64

static int rmt_final_dense_solve(double *A, double *b, int n) {
    for (int k=0; k<n; ++k) {
        int piv = k;
        double pv = fabs(A[(size_t)k*n+k]);
        for (int i=k+1; i<n; ++i) {
            const double v = fabs(A[(size_t)i*n+k]);
            if (v > pv) { pv=v; piv=i; }
        }
        if (!(pv > 1.0e-300) || !isfinite(pv)) return 0;
        if (piv != k) {
            for (int j=k; j<n; ++j) {
                const double t=A[(size_t)k*n+j];
                A[(size_t)k*n+j]=A[(size_t)piv*n+j];
                A[(size_t)piv*n+j]=t;
            }
            const double t=b[k]; b[k]=b[piv]; b[piv]=t;
        }
        const double akk=A[(size_t)k*n+k];
        for (int i=k+1; i<n; ++i) {
            const double f=A[(size_t)i*n+k]/akk;
            A[(size_t)i*n+k]=0.0;
            for (int j=k+1; j<n; ++j)
                A[(size_t)i*n+j]-=f*A[(size_t)k*n+j];
            b[i]-=f*b[k];
        }
    }
    for (int i=n-1; i>=0; --i) {
        double s=b[i];
        for (int j=i+1; j<n; ++j) s-=A[(size_t)i*n+j]*b[j];
        const double aii=A[(size_t)i*n+i];
        if (!(fabs(aii)>1.0e-300) || !isfinite(aii)) return 0;
        b[i]=s/aii;
    }
    return 1;
}

/* Exact solve of every independent shifted grid on the coarsest level. */
static void rmt_final_direct_coarsest(const Grid *g, double *c,
                                      const double *rhs,
                                      int sx, int sy,
                                      long long *solveCount,
                                      long long *unknownCount) {
    const int ngx = MIN(sx,g->Nx);
    const int ngy = MIN(sy,g->Ny);
    const int ng = ngx*ngy;
    long long solves=0, unknowns=0;
    int failed=0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,32) reduction(+:solves,unknowns) reduction(|:failed) if(ng>8)
#endif
    for (int gid=0; gid<ng; ++gid) {
        const int ox = 1 + gid/ngy;
        const int oy = 1 + gid%ngy;
        const int nxg = (ox<=g->Nx) ? 1+(g->Nx-ox)/sx : 0;
        const int nyg = (oy<=g->Ny) ? 1+(g->Ny-oy)/sy : 0;
        const int m = nxg*nyg;
        if (m <= 0) continue;
        if (m > RMT_FINAL_DIRECT_MAX) {
            failed = 1;
            continue;
        }

        double A[RMT_FINAL_DIRECT_MAX*RMT_FINAL_DIRECT_MAX];
        double b[RMT_FINAL_DIRECT_MAX];
        memset(A,0,(size_t)m*(size_t)m*sizeof(double));

        for (int a=0; a<nxg; ++a) for (int q=0; q<nyg; ++q) {
            const int i=ox+a*sx, j=oy+q*sy;
            const int row=a*nyg+q;
            double ap=0.0, aw=0.0, ae=0.0, as=0.0, an=0.0;
            rmt_final_dim_coeff(g->Nx,i,sx,g->dx,0,&ap,&aw,&ae);
            rmt_final_dim_coeff(g->Ny,j,sy,g->dy,1,&ap,&as,&an);
            A[(size_t)row*m+row]=ap;
            if (a>0)       A[(size_t)row*m+(row-nyg)] = -aw;
            if (a+1<nxg)   A[(size_t)row*m+(row+nyg)] = -ae;
            if (q>0)       A[(size_t)row*m+(row-1)]   = -as;
            if (q+1<nyg)   A[(size_t)row*m+(row+1)]   = -an;
            b[row]=rhs[IDX(i,j,g->Ny)];
        }

        if (!rmt_final_dense_solve(A,b,m)) {
            failed = 1;
            continue;
        }

        for (int a=0; a<nxg; ++a) for (int q=0; q<nyg; ++q) {
            const int i=ox+a*sx, j=oy+q*sy;
            c[IDX(i,j,g->Ny)]=b[a*nyg+q];
        }
        ++solves;
        unknowns += m;
    }

    if (failed) die("RMT FINAL coarsest direct solve failed or exceeded direct block limit");
    if (solveCount) *solveCount += solves;
    if (unknownCount) *unknownCount += unknowns;
}

static void rmt_final_build_correction(const Grid *g, const double *defect,
                                       double *c,
                                       const RMTFinalConfig *cfg,
                                       RMTFinalDiag *diag) {
    RMTFinalWorkspace *w = rmt_final_workspace(g);
    int lx = rmt_final_auto_level_1d(g->Nx,cfg->nCoarsest);
    int ly = rmt_final_auto_level_1d(g->Ny,cfg->nCoarsest);
    if (cfg->requestedLevels > 0) {
        const int cap = MAX(0,cfg->requestedLevels-1);
        lx = MIN(lx,cap);
        ly = MIN(ly,cap);
    }
    const int lm = MAX(lx,ly);
    long long updates=0, directSolves=0, directUnknowns=0;
    memset(c,0,w->n*sizeof(double));

    for (int lc=lm; lc>=0; --lc) {
        const int lX=MIN(lx,lc), lY=MIN(ly,lc);
        const int sx=rmt_final_ipow3(lX), sy=rmt_final_ipow3(lY);
        rmt_final_restrict_level(g,defect,w->rhsLevel,w->prefix,sx,sy);

        if (lc == lm) {
            rmt_final_direct_coarsest(g,c,w->rhsLevel,sx,sy,
                                      &directSolves,&directUnknowns);
        } else {
            updates += rmt_final_smooth_level(g,c,w->rhsLevel,sx,sy,
                                              MAX(1,cfg->smoothSweeps));
        }
    }

    g_rmtFinalPointUpdatesTotal += updates;
    g_rmtFinalDirectSolvesTotal += directSolves;
    g_rmtFinalDirectUnknownsTotal += directUnknowns;
    g_rmtFinalLastLevelX = lx;
    g_rmtFinalLastLevelY = ly;

    if (diag) {
        diag->levelXMax=lx;
        diag->levelYMax=ly;
        diag->pointUpdates=updates;
        diag->directSolves=directSolves;
        diag->directUnknowns=directUnknowns;
    }
}

#endif
