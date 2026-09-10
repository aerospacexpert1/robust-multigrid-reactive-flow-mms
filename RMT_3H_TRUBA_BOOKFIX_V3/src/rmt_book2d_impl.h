#ifndef RMT_BOOK2D_IMPL_H
#define RMT_BOOK2D_IMPL_H

typedef struct {
    int requestedLevels, smoothSweeps, coarseSweeps, nCoarsest;
} RMTBookConfig;

typedef struct {
    int levelXMax, levelYMax;
    long long pointUpdates;
} RMTBookDiag;

typedef struct {
    int nx, ny;
    size_t n;
    double *remaining, *rhsLevel, *prefix, *corr, *res;
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
    free(g_rmtBookWS.remaining);
    free(g_rmtBookWS.rhsLevel);
    free(g_rmtBookWS.prefix);
    free(g_rmtBookWS.corr);
    free(g_rmtBookWS.res);
    memset(&g_rmtBookWS, 0, sizeof(g_rmtBookWS));
}

static RMTBookWorkspace *rmt_book_workspace(const Grid *g) {
    size_t n = (size_t)(g->Nx+2)*(size_t)(g->Ny+2);
    size_t np = (size_t)(g->Nx+1)*(size_t)(g->Ny+1);
    if (g_rmtBookWS.nx == g->Nx && g_rmtBookWS.ny == g->Ny && g_rmtBookWS.corr)
        return &g_rmtBookWS;

    rmt_book_workspace_release();
    g_rmtBookWS.nx = g->Nx;
    g_rmtBookWS.ny = g->Ny;
    g_rmtBookWS.n = n;
    g_rmtBookWS.remaining = calloc(n, sizeof(double));
    g_rmtBookWS.rhsLevel = calloc(n, sizeof(double));
    g_rmtBookWS.prefix = calloc(np, sizeof(double));
    g_rmtBookWS.corr = calloc(n, sizeof(double));
    g_rmtBookWS.res = calloc(n, sizeof(double));
    if (!g_rmtBookWS.remaining || !g_rmtBookWS.rhsLevel ||
        !g_rmtBookWS.prefix || !g_rmtBookWS.corr || !g_rmtBookWS.res)
        die("alloc RMT book workspace");
    return &g_rmtBookWS;
}

static int rmt_book_auto_level_1d(int n, int nc) {
    int l = 0, s = 1;
    if (nc < 2) nc = 2;
    while ((n+s-1)/s > nc) {
        s *= 3;
        ++l;
        if (l >= 12) break;
    }
    return l;
}

static int rmt_book_ipow3(int l) {
    int s = 1;
    while (l-- > 0) s *= 3;
    return s;
}

static void rmt_book_apply_fine_A(const Grid *g, const double *x, double *Ax) {
    int nx = g->Nx, ny = g->Ny;
    double ax = 1.0/(g->dx*g->dx);
    double ay = 1.0/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i=1; i<=nx; ++i) for (int j=1; j<=ny; ++j) {
        double ap = 0.0, sum = 0.0;
        if (i > 1) { ap += ax; sum += ax*x[IDX(i-1,j,ny)]; }
        else       { ap += 2.0*ax; }
        if (i < nx){ ap += ax; sum += ax*x[IDX(i+1,j,ny)]; }
        else       { ap += 2.0*ax; }
        if (j > 1) { ap += ay; sum += ay*x[IDX(i,j-1,ny)]; }
        if (j < ny){ ap += ay; sum += ay*x[IDX(i,j+1,ny)]; }
        Ax[IDX(i,j,ny)] = ap*x[IDX(i,j,ny)] - sum;
    }
}

static void rmt_book_prefix_build(const Grid *g, const double *a, double *p) {
    int ny = g->Ny, pitch = ny+1;
    memset(p, 0, (size_t)(g->Nx+1)*(size_t)(ny+1)*sizeof(double));
    for (int i=1; i<=g->Nx; ++i) {
        double row = 0.0;
        for (int j=1; j<=ny; ++j) {
            row += a[IDX(i,j,ny)];
            p[(size_t)i*pitch+j] = p[(size_t)(i-1)*pitch+j] + row;
        }
    }
}

static double rmt_book_box_average(const Grid *g, const double *p,
                                   int i, int j, int sx, int sy) {
    int pitch = g->Ny+1;
    int hx = (sx-1)/2, hy = (sy-1)/2;
    int i1 = MAX(1, i-hx), i2 = MIN(g->Nx, i+hx);
    int j1 = MAX(1, j-hy), j2 = MIN(g->Ny, j+hy);
    double sum = p[(size_t)i2*pitch+j2]
               - p[(size_t)(i1-1)*pitch+j2]
               - p[(size_t)i2*pitch+j1-1]
               + p[(size_t)(i1-1)*pitch+j1-1];
    double cnt = (double)(i2-i1+1)*(double)(j2-j1+1);
    return cnt > 0.0 ? sum/cnt : 0.0;
}

static void rmt_book_restrict_level(const Grid *g, const double *r,
                                    double *rhs, double *p, int sx, int sy) {
    rmt_book_prefix_build(g, r, p);
    for (int i=1; i<=g->Nx; ++i)
        for (int j=1; j<=g->Ny; ++j)
            rhs[IDX(i,j,g->Ny)] = rmt_book_box_average(g,p,i,j,sx,sy);
}

static inline void rmt_book_dim_coeff_cell_fv(double H, double dist,
                                               int hm, int hp, int neumann,
                                               double *ap, double *am, double *az) {
    *am = 0.0;
    *az = 0.0;
    if (hm && hp) {
        double q = 1.0/(H*H);
        *ap += 2.0*q;
        *am = q;
        *az = q;
        return;
    }
    dist = MAX(dist, 1.0e-14*H);
    double V = dist + 0.5*H;
    double aNbr = 1.0/(H*V);
    double aBoundary = neumann ? 0.0 : 1.0/(dist*V);
    *ap += aNbr + aBoundary;
    if (!hm && hp) *az = aNbr;
    else if (hm && !hp) *am = aNbr;
    else *ap += aBoundary;
}

static inline double rmt_book_update_point(const Grid *g, double *c,
                                           const double *rhs,
                                           int i, int j, int sx, int sy) {
    int nx = g->Nx, ny = g->Ny;
    int hw = (i-sx >= 1), he = (i+sx <= nx);
    int hs = (j-sy >= 1), hn = (j+sy <= ny);
    double Hx = sx*g->dx, Hy = sy*g->dy;
    double ap = 0.0, aw = 0.0, ae = 0.0, as = 0.0, an = 0.0;
    double dxBoundary = !hw ? (i-0.5)*g->dx : (nx-i+0.5)*g->dx;
    double dyBoundary = !hs ? (j-0.5)*g->dy : (ny-j+0.5)*g->dy;
    rmt_book_dim_coeff_cell_fv(Hx,dxBoundary,hw,he,0,&ap,&aw,&ae);
    rmt_book_dim_coeff_cell_fv(Hy,dyBoundary,hs,hn,1,&ap,&as,&an);
    double sum = 0.0;
    if (hw) sum += aw*c[IDX(i-sx,j,ny)];
    if (he) sum += ae*c[IDX(i+sx,j,ny)];
    if (hs) sum += as*c[IDX(i,j-sy,ny)];
    if (hn) sum += an*c[IDX(i,j+sy,ny)];
    return (rhs[IDX(i,j,ny)] + sum)/MAX(ap,1.0e-30);
}

static long long rmt_book_smooth_level(const Grid *g, double *c,
                                       const double *rhs,
                                       int sx, int sy, int sweeps) {
    int ngx = MIN(sx,g->Nx), ngy = MIN(sy,g->Ny), ng = ngx*ngy;
    long long updates = 0;
    for (int sw=0; sw<sweeps; ++sw) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:updates) if(ng>1)
#endif
        for (int gid=0; gid<ng; ++gid) {
            int ox = 1 + gid/ngy;
            int oy = 1 + gid%ngy;
            for (int i=ox; i<=g->Nx; i+=sx)
                for (int j=oy; j<=g->Ny; j+=sy) {
                    c[IDX(i,j,g->Ny)] = rmt_book_update_point(g,c,rhs,i,j,sx,sy);
                    ++updates;
                }
        }
    }
    return updates;
}

static void rmt_book_build_correction(const Grid *g, const double *b, double *c,
                                      const RMTBookConfig *cfg, RMTBookDiag *d) {
    RMTBookWorkspace *w = rmt_book_workspace(g);
    int lx = rmt_book_auto_level_1d(g->Nx,cfg->nCoarsest);
    int ly = rmt_book_auto_level_1d(g->Ny,cfg->nCoarsest);
    if (cfg->requestedLevels > 0) {
        int cap = MAX(0,cfg->requestedLevels-1);
        lx = MIN(lx,cap);
        ly = MIN(ly,cap);
    }
    int lm = MAX(lx,ly);
    long long updates = 0;
    memset(c,0,w->n*sizeof(double));
    for (int lc=lm; lc>=0; --lc) {
        int lX = MIN(lx,lc), lY = MIN(ly,lc);
        int sx = rmt_book_ipow3(lX), sy = rmt_book_ipow3(lY);
        rmt_book_restrict_level(g,b,w->rhsLevel,w->prefix,sx,sy);
        int sw = (lc == lm) ? MAX(1,cfg->coarseSweeps) : MAX(1,cfg->smoothSweeps);
        updates += rmt_book_smooth_level(g,c,w->rhsLevel,sx,sy,sw);
    }
    if (d) {
        d->levelXMax = lx;
        d->levelYMax = ly;
        d->pointUpdates = updates;
    }
    g_rmtBookPointUpdatesTotal += updates;
    g_rmtBookLastLevelX = lx;
    g_rmtBookLastLevelY = ly;
}

#endif
