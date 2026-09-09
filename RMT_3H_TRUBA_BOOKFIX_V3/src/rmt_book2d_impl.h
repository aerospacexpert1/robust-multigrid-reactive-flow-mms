#ifndef RMT_BOOK2D_IMPL_H
#define RMT_BOOK2D_IMPL_H

/* Martynenko-style 2-D pressure RMT on finest-grid storage. */
typedef struct {
    int requestedLevels;
    int smoothSweeps;
    int coarseSweeps;
    int nCoarsest;
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
    const size_t n=(size_t)(g->Nx+2)*(size_t)(g->Ny+2);
    const size_t np=(size_t)(g->Nx+1)*(size_t)(g->Ny+1);
    if (g_rmtBookWS.nx==g->Nx && g_rmtBookWS.ny==g->Ny && g_rmtBookWS.corr) return &g_rmtBookWS;
    rmt_book_workspace_release();
    g_rmtBookWS.nx=g->Nx; g_rmtBookWS.ny=g->Ny; g_rmtBookWS.n=n;
    g_rmtBookWS.remaining=(double*)calloc(n,sizeof(double));
    g_rmtBookWS.rhsLevel=(double*)calloc(n,sizeof(double));
    g_rmtBookWS.prefix=(double*)calloc(np,sizeof(double));
    g_rmtBookWS.corr=(double*)calloc(n,sizeof(double));
    g_rmtBookWS.res=(double*)calloc(n,sizeof(double));
    if(!g_rmtBookWS.remaining||!g_rmtBookWS.rhsLevel||!g_rmtBookWS.prefix||!g_rmtBookWS.corr||!g_rmtBookWS.res)
        die("alloc RMT book workspace");
    return &g_rmtBookWS;
}

static int rmt_book_auto_level_1d(int n,int nCoarsest) {
    int level=0; long long stride=1;
    if(nCoarsest<2) nCoarsest=2;
    while((n+(int)stride-1)/(int)stride>nCoarsest) {
        stride*=3; ++level; if(level>=12) break;
    }
    return level;
}

static int rmt_book_ipow3(int level) { int s=1; for(int k=0;k<level;++k)s*=3; return s; }

/* Finest-grid pressure correction operator: x homogeneous Dirichlet, y homogeneous Neumann. */
static void rmt_book_apply_fine_A(const Grid *g,const double *x,double *Ax) {
    const int nx=g->Nx,ny=g->Ny;
    const double ax=1.0/(g->dx*g->dx),ay=1.0/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for(int i=1;i<=nx;++i) for(int j=1;j<=ny;++j) {
        double ap=0.0,sum=0.0;
        if(i>1){ap+=ax;sum+=ax*x[IDX(i-1,j,ny)];} else ap+=2.0*ax;
        if(i<nx){ap+=ax;sum+=ax*x[IDX(i+1,j,ny)];} else ap+=2.0*ax;
        if(j>1){ap+=ay;sum+=ay*x[IDX(i,j-1,ny)];}
        if(j<ny){ap+=ay;sum+=ay*x[IDX(i,j+1,ny)];}
        Ax[IDX(i,j,ny)]=ap*x[IDX(i,j,ny)]-sum;
    }
}

static void rmt_book_remaining_residual(const Grid *g,const double *b,const double *corr,double *remaining) {
    rmt_book_apply_fine_A(g,corr,remaining);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for(int i=1;i<=g->Nx;++i) for(int j=1;j<=g->Ny;++j)
        remaining[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-remaining[IDX(i,j,g->Ny)];
}

static void rmt_book_prefix_build(const Grid *g,const double *a,double *pfx) {
    const int nx=g->Nx,ny=g->Ny,pitch=ny+1;
    memset(pfx,0,(size_t)(nx+1)*(size_t)(ny+1)*sizeof(double));
    for(int i=1;i<=nx;++i){double row=0.0;for(int j=1;j<=ny;++j){
        row+=a[IDX(i,j,ny)]; pfx[(size_t)i*pitch+j]=pfx[(size_t)(i-1)*pitch+j]+row;
    }}
}

/* Book Eq. (2.16): arithmetic averaging of the finest-grid residuals inside the real CV part. */
static double rmt_book_box_average(const Grid *g,const double *pfx,int i,int j,int sx,int sy) {
    const int pitch=g->Ny+1,hx=(sx-1)/2,hy=(sy-1)/2;
    const int i1=MAX(1,i-hx),i2=MIN(g->Nx,i+hx),j1=MAX(1,j-hy),j2=MIN(g->Ny,j+hy);
    const double sum=pfx[(size_t)i2*pitch+j2]-pfx[(size_t)(i1-1)*pitch+j2]
                    -pfx[(size_t)i2*pitch+(j1-1)]+pfx[(size_t)(i1-1)*pitch+(j1-1)];
    const double cnt=(double)(i2-i1+1)*(double)(j2-j1+1);
    return cnt>0.0?sum/cnt:0.0;
}

static void rmt_book_restrict_level(const Grid *g,const double *remaining,double *rhsLevel,double *pfx,int sx,int sy) {
    rmt_book_prefix_build(g,remaining,pfx);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for(int i=1;i<=g->Nx;++i) for(int j=1;j<=g->Ny;++j)
        rhsLevel[IDX(i,j,g->Ny)]=rmt_book_box_average(g,pfx,i,j,sx,sy);
}

/*
 * Add one coordinate contribution to A*c = rhs, eliminating a virtual endpoint with
 * Martynenko Eqs. (2.17)-(2.20).  For a missing neighbour:
 *   ghost = alpha*c0 + beta*c1 + boundary-source.
 * The correction boundary data are homogeneous, therefore the source term is zero.
 */
static inline void rmt_book_dim_coeff(double H,double dist,int hasMinus,int hasPlus,int isNeumann,
                                      double *ap,double *aMinus,double *aPlus) {
    const double invH2=1.0/(H*H);
    *aMinus=0.0; *aPlus=0.0;
    if(hasMinus && hasPlus) {
        *ap += 2.0*invH2; *aMinus=invH2; *aPlus=invH2; return;
    }
    /* The shifted endpoint is at xi coarse spacings from the physical boundary. */
    double xi=dist/H;
    if(xi<1.0e-12) xi=1.0e-12;
    double alpha,beta;
    if(isNeumann) {
        /* Eq. (2.18)/(2.20), homogeneous derivative. */
        alpha=4.0*xi/(2.0*xi+1.0);
        beta=-(2.0*xi-1.0)/(2.0*xi+1.0);
    } else {
        /* Eq. (2.17)/(2.19), homogeneous correction value. */
        alpha=(2.0*xi-1.0)/xi;
        beta=-(xi-1.0)/(xi+1.0);
    }
    *ap += (2.0-alpha)*invH2;
    const double an=(1.0+beta)*invH2;
    if(!hasMinus && hasPlus) *aPlus=an;
    else if(hasMinus && !hasPlus) *aMinus=an;
}

static inline double rmt_book_update_point(const Grid *g,double *corr,const double *rhs,int i,int j,int sx,int sy) {
    const int nx=g->Nx,ny=g->Ny;
    const double Hx=(double)sx*g->dx,Hy=(double)sy*g->dy;
    const int hasW=(i-sx>=1),hasE=(i+sx<=nx),hasS=(j-sy>=1),hasN=(j+sy<=ny);
    double ap=0.0,aW=0.0,aE=0.0,aS=0.0,aN=0.0;
    const double distX=!hasW?(i-0.5)*g->dx:(nx-i+0.5)*g->dx;
    const double distY=!hasS?(j-0.5)*g->dy:(ny-j+0.5)*g->dy;
    rmt_book_dim_coeff(Hx,distX,hasW,hasE,0,&ap,&aW,&aE);
    rmt_book_dim_coeff(Hy,distY,hasS,hasN,1,&ap,&aS,&aN);
    double sum=0.0;
    if(hasW)sum+=aW*corr[IDX(i-sx,j,ny)];
    if(hasE)sum+=aE*corr[IDX(i+sx,j,ny)];
    if(hasS)sum+=aS*corr[IDX(i,j-sy,ny)];
    if(hasN)sum+=aN*corr[IDX(i,j+sy,ny)];
    return (rhs[IDX(i,j,ny)]+sum)/MAX(ap,1.0e-30);
}

static long long rmt_book_smooth_level(const Grid *g,double *corr,const double *rhs,int sx,int sy,int sweeps) {
    const int ngx=MIN(sx,g->Nx),ngy=MIN(sy,g->Ny),ngrids=ngx*ngy;
    long long updates=0;
    for(int sw=0;sw<sweeps;++sw) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:updates) if(ngrids>1)
#endif
        for(int gid=0;gid<ngrids;++gid) {
            const int ox=1+gid/ngy,oy=1+gid%ngy;
            for(int i=ox;i<=g->Nx;i+=sx) for(int j=oy;j<=g->Ny;j+=sy) {
                corr[IDX(i,j,g->Ny)]=rmt_book_update_point(g,corr,rhs,i,j,sx,sy); ++updates;
            }
        }
    }
    return updates;
}

static void rmt_book_build_correction(const Grid *g,const double *b,double *corr,const RMTBookConfig *cfg,RMTBookDiag *diag) {
    RMTBookWorkspace *ws=rmt_book_workspace(g);
    int lx=rmt_book_auto_level_1d(g->Nx,cfg->nCoarsest),ly=rmt_book_auto_level_1d(g->Ny,cfg->nCoarsest);
    if(cfg->requestedLevels>0){const int cap=MAX(0,cfg->requestedLevels-1);lx=MIN(lx,cap);ly=MIN(ly,cap);}
    const int lm=MAX(lx,ly); long long updates=0;
    memset(corr,0,ws->n*sizeof(double));
    for(int levelC=lm;levelC>=0;--levelC) {
        const int levelX=MIN(lx,levelC),levelY=MIN(ly,levelC);
        const int sx=rmt_book_ipow3(levelX),sy=rmt_book_ipow3(levelY);
        rmt_book_remaining_residual(g,b,corr,ws->remaining);
        rmt_book_restrict_level(g,ws->remaining,ws->rhsLevel,ws->prefix,sx,sy);
        const int sweeps=(levelC==lm)?MAX(1,cfg->coarseSweeps):MAX(1,cfg->smoothSweeps);
        updates+=rmt_book_smooth_level(g,corr,ws->rhsLevel,sx,sy,sweeps);
    }
    if(diag){diag->levelXMax=lx;diag->levelYMax=ly;diag->pointUpdates=updates;}
    g_rmtBookPointUpdatesTotal+=updates;g_rmtBookLastLevelX=lx;g_rmtBookLastLevelY=ly;
}

#endif
