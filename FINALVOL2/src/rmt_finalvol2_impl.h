#ifndef RMT_FINALVOL2_IMPL_H
#define RMT_FINALVOL2_IMPL_H

/*
 * FINALVOL2
 *
 * Martynenko-style factor-three multiple-coarse-grid RMT for the frozen V95
 * cell-centred pressure equation.
 *
 * Numerical algorithm is intentionally unchanged from RMT FINAL:
 *   - independent factor-three x/y hierarchy;
 *   - coarsest -> finest sawtooth, no presmoothing;
 *   - exact fine-grid defect restricted by control-volume averaging;
 *   - shifted-grid finite-volume operator with the same mixed correction BCs;
 *   - direct solve on every coarsest shifted grid;
 *   - RBGS postsmoothing on all finer levels;
 *   - correction stored on the fine index space and added in full;
 *   - no interpolation, no monotonic line-search and no hidden fallback.
 *
 * FINALVOL2 only removes implementation overhead:
 *   1) one fine-defect prefix integral per RMT cycle (not per level);
 *   2) precomputed shifted-grid coefficients for every level;
 *   3) precomputed red/black point maps (no modulo/division in hot smoother);
 *   4) cached LU factors for the invariant coarsest shifted-grid matrices.
 *
 * The smoothing count remains a single global Martynenko-style NSIL analogue.
 * Production default is deliberately NOT retuned here; calibration is separate.
 */

#define RMTV2_MAX_LEVELS 16
#define RMTV2_DIRECT_MAX 64

typedef struct {
    int requestedLevels;
    int smoothSweeps;
    int nCoarsest;
} RMTV2Config;

typedef struct {
    int levelXMax, levelYMax;
    long long pointUpdates;
    long long directSolves;
    long long directUnknowns;
    long long luFactorizations;
} RMTV2Diag;

typedef struct { int i,j; } RMTV2Point;

typedef struct {
    int sx, sy;
    double *xAp,*xM,*xP;
    double *yAp,*yM,*yP;
    RMTV2Point *red,*black;
    int nRed,nBlack;
} RMTV2LevelPlan;

typedef struct {
    int sx,sy,ng;
    int totalUnknowns;
    size_t totalLU;
    int *ox,*oy,*nxg,*nyg,*m;
    size_t *luOffset;
    int *pivOffset;
    double *luPool;
    int *pivPool;
} RMTV2CoarsePlan;

typedef struct {
    int nx,ny;
    size_t n;
    double *rhsLevel;
    double *prefix;
    double *corr;
    double *res;

    int planReady;
    int planRequestedLevels;
    int planNCoarsest;
    int lx,ly,lm;
    RMTV2LevelPlan level[RMTV2_MAX_LEVELS];
    RMTV2CoarsePlan coarse;
} RMTV2Workspace;

static RMTV2Workspace g_rmtV2WS = {0};
static long long g_rmtV2PointUpdatesTotal = 0;
static long long g_rmtV2DirectSolvesTotal = 0;
static long long g_rmtV2DirectUnknownsTotal = 0;
static long long g_rmtV2LUFactorizationsTotal = 0;
static long long g_rmtV2NonmonotoneCycles = 0;
static double g_rmtV2LastCycleRatio = 0.0;
static double g_rmtV2MaxCycleRatio = 0.0;
static int g_rmtV2LastLevelX = 0;
static int g_rmtV2LastLevelY = 0;

static void rmt_v2_level_release(RMTV2LevelPlan *p){
    free(p->xAp); free(p->xM); free(p->xP);
    free(p->yAp); free(p->yM); free(p->yP);
    free(p->red); free(p->black);
    memset(p,0,sizeof(*p));
}

static void rmt_v2_coarse_release(RMTV2CoarsePlan *p){
    free(p->ox); free(p->oy); free(p->nxg); free(p->nyg); free(p->m);
    free(p->luOffset); free(p->pivOffset);
    free(p->luPool); free(p->pivPool);
    memset(p,0,sizeof(*p));
}

static void rmt_v2_workspace_release(void){
    for(int l=0;l<RMTV2_MAX_LEVELS;++l) rmt_v2_level_release(&g_rmtV2WS.level[l]);
    rmt_v2_coarse_release(&g_rmtV2WS.coarse);
    free(g_rmtV2WS.rhsLevel);
    free(g_rmtV2WS.prefix);
    free(g_rmtV2WS.corr);
    free(g_rmtV2WS.res);
    memset(&g_rmtV2WS,0,sizeof(g_rmtV2WS));
}

static RMTV2Workspace *rmt_v2_workspace(const Grid *g){
    const size_t n=(size_t)(g->Nx+2)*(size_t)(g->Ny+2);
    const size_t np=(size_t)(g->Nx+1)*(size_t)(g->Ny+1);
    if(g_rmtV2WS.nx==g->Nx && g_rmtV2WS.ny==g->Ny && g_rmtV2WS.corr)
        return &g_rmtV2WS;
    rmt_v2_workspace_release();
    g_rmtV2WS.nx=g->Nx; g_rmtV2WS.ny=g->Ny; g_rmtV2WS.n=n;
    g_rmtV2WS.rhsLevel=calloc(n,sizeof(double));
    g_rmtV2WS.prefix=calloc(np,sizeof(double));
    g_rmtV2WS.corr=calloc(n,sizeof(double));
    g_rmtV2WS.res=calloc(n,sizeof(double));
    if(!g_rmtV2WS.rhsLevel||!g_rmtV2WS.prefix||!g_rmtV2WS.corr||!g_rmtV2WS.res)
        die("alloc FINALVOL2 workspace");
    return &g_rmtV2WS;
}

static int rmt_v2_auto_level_1d(int n,int nc){
    int l=0,s=1;
    if(nc<2)nc=2;
    while((n+s-1)/s>nc){
        if(s>INT_MAX/3)break;
        s*=3; ++l;
        if(l>=RMTV2_MAX_LEVELS-1)break;
    }
    return l;
}
static int rmt_v2_ipow3(int l){int s=1;while(l-->0)s*=3;return s;}

static void rmt_v2_apply_fine_A(const Grid*g,const double*x,double*Ax){
    const int nx=g->Nx,ny=g->Ny;
    const double ax=1.0/(g->dx*g->dx),ay=1.0/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j){
        double ap=0.0,sum=0.0;
        if(i>1){ap+=ax;sum+=ax*x[IDX(i-1,j,ny)];}else ap+=2.0*ax;
        if(i<nx){ap+=ax;sum+=ax*x[IDX(i+1,j,ny)];}else ap+=2.0*ax;
        if(j>1){ap+=ay;sum+=ay*x[IDX(i,j-1,ny)];}
        if(j<ny){ap+=ay;sum+=ay*x[IDX(i,j+1,ny)];}
        Ax[IDX(i,j,ny)]=ap*x[IDX(i,j,ny)]-sum;
    }
}

static void rmt_v2_prefix_build(const Grid*g,const double*a,double*p){
    const int ny=g->Ny,pitch=ny+1;
    memset(p,0,(size_t)(g->Nx+1)*(size_t)(ny+1)*sizeof(double));
    for(int i=1;i<=g->Nx;++i){
        double row=0.0;
        for(int j=1;j<=ny;++j){
            row+=a[IDX(i,j,ny)];
            p[(size_t)i*pitch+j]=p[(size_t)(i-1)*pitch+j]+row;
        }
    }
}

static inline void rmt_v2_cv_bounds_1d(int n,int idx,int stride,int*lo,int*hi){
    const int half=(stride-1)/2;
    *lo=(idx-stride>=1)?idx-half:1;
    *hi=(idx+stride<=n)?idx+half:n;
}

static inline double rmt_v2_cv_average(const Grid*g,const double*p,int i,int j,int sx,int sy){
    const int pitch=g->Ny+1;
    int i1,i2,j1,j2;
    rmt_v2_cv_bounds_1d(g->Nx,i,sx,&i1,&i2);
    rmt_v2_cv_bounds_1d(g->Ny,j,sy,&j1,&j2);
    const double sum=p[(size_t)i2*pitch+j2]-p[(size_t)(i1-1)*pitch+j2]
                    -p[(size_t)i2*pitch+j1-1]+p[(size_t)(i1-1)*pitch+j1-1];
    const double cnt=(double)(i2-i1+1)*(double)(j2-j1+1);
    return sum/cnt;
}

/* Prefix is built once per RMT cycle; each level only samples it. */
static void rmt_v2_restrict_from_prefix(const Grid*g,const double*prefix,double*rhs,int sx,int sy){
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)
        rhs[IDX(i,j,g->Ny)]=rmt_v2_cv_average(g,prefix,i,j,sx,sy);
}

static inline void rmt_v2_dim_coeff(int n,int idx,int stride,double h,int neumann,
                                    double*ap,double*am,double*az){
    const int hm=(idx-stride>=1),hp=(idx+stride<=n);
    const double H=stride*h;
    *am=0.0;*az=0.0;
    if(hm&&hp){
        const double a=1.0/(H*H);
        *ap+=2.0*a;*am=a;*az=a;return;
    }
    if(!hm&&hp){
        const double d=(idx-0.5)*h,V=d+0.5*H;
        const double aN=1.0/(H*V),aB=neumann?0.0:1.0/(d*V);
        *ap+=aN+aB;*az=aN;return;
    }
    if(hm&&!hp){
        const double d=(n-idx+0.5)*h,V=d+0.5*H;
        const double aN=1.0/(H*V),aB=neumann?0.0:1.0/(d*V);
        *ap+=aN+aB;*am=aN;return;
    }
    {
        const double V=n*h;
        if(!neumann){
            const double d0=(idx-0.5)*h,d1=(n-idx+0.5)*h;
            *ap+=1.0/(d0*V)+1.0/(d1*V);
        }
    }
}

static void rmt_v2_build_level_plan(const Grid*g,RMTV2LevelPlan*p,int sx,int sy,int needMap){
    p->sx=sx;p->sy=sy;
    const size_t nx1=(size_t)g->Nx+1,ny1=(size_t)g->Ny+1;
    p->xAp=calloc(nx1,sizeof(double));p->xM=calloc(nx1,sizeof(double));p->xP=calloc(nx1,sizeof(double));
    p->yAp=calloc(ny1,sizeof(double));p->yM=calloc(ny1,sizeof(double));p->yP=calloc(ny1,sizeof(double));
    if(!p->xAp||!p->xM||!p->xP||!p->yAp||!p->yM||!p->yP)die("alloc FINALVOL2 level coefficients");

    for(int i=1;i<=g->Nx;++i){
        double ap=0,am=0,az=0;
        rmt_v2_dim_coeff(g->Nx,i,sx,g->dx,0,&ap,&am,&az);
        p->xAp[i]=ap;p->xM[i]=am;p->xP[i]=az;
    }
    for(int j=1;j<=g->Ny;++j){
        double ap=0,am=0,az=0;
        rmt_v2_dim_coeff(g->Ny,j,sy,g->dy,1,&ap,&am,&az);
        p->yAp[j]=ap;p->yM[j]=am;p->yP[j]=az;
    }

    if(!needMap)return;
    const int nxy=g->Nx*g->Ny;
    p->red=malloc((size_t)nxy*sizeof(*p->red));
    p->black=malloc((size_t)nxy*sizeof(*p->black));
    if(!p->red||!p->black)die("alloc FINALVOL2 point maps");
    int nr=0,nb=0;
    for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j){
        const int ox=(i-1)%sx+1,oy=(j-1)%sy+1;
        const int qi=(i-ox)/sx,qj=(j-oy)/sy;
        RMTV2Point q={i,j};
        if(((qi+qj)&1)==0)p->red[nr++]=q;else p->black[nb++]=q;
    }
    p->nRed=nr;p->nBlack=nb;
}

static inline double rmt_v2_update_planned(const Grid*g,double*c,const double*rhs,
                                            const RMTV2LevelPlan*p,int i,int j){
    const int ny=g->Ny,sx=p->sx,sy=p->sy;
    const double aw=p->xM[i],ae=p->xP[i],as=p->yM[j],an=p->yP[j];
    const double ap=p->xAp[i]+p->yAp[j];
    double sum=0.0;
    if(aw!=0.0)sum+=aw*c[IDX(i-sx,j,ny)];
    if(ae!=0.0)sum+=ae*c[IDX(i+sx,j,ny)];
    if(as!=0.0)sum+=as*c[IDX(i,j-sy,ny)];
    if(an!=0.0)sum+=an*c[IDX(i,j+sy,ny)];
    return (rhs[IDX(i,j,ny)]+sum)/MAX(ap,1.0e-300);
}

static long long rmt_v2_smooth_level(const Grid*g,double*c,const double*rhs,
                                     const RMTV2LevelPlan*p,int sweeps){
    long long updates=0;
    for(int sw=0;sw<sweeps;++sw){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:updates) if(p->nRed>2048)
#endif
        for(int k=0;k<p->nRed;++k){
            const int i=p->red[k].i,j=p->red[k].j;
            c[IDX(i,j,g->Ny)]=rmt_v2_update_planned(g,c,rhs,p,i,j);
            ++updates;
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:updates) if(p->nBlack>2048)
#endif
        for(int k=0;k<p->nBlack;++k){
            const int i=p->black[k].i,j=p->black[k].j;
            c[IDX(i,j,g->Ny)]=rmt_v2_update_planned(g,c,rhs,p,i,j);
            ++updates;
        }
    }
    return updates;
}

/* LU factorization with stored pivot sequence. */
static int rmt_v2_lu_factor(double*A,int*piv,int n){
    for(int k=0;k<n;++k){
        int pk=k;double pv=fabs(A[(size_t)k*n+k]);
        for(int i=k+1;i<n;++i){const double v=fabs(A[(size_t)i*n+k]);if(v>pv){pv=v;pk=i;}}
        if(!(pv>1e-300)||!isfinite(pv))return 0;
        piv[k]=pk;
        if(pk!=k)for(int j=0;j<n;++j){double t=A[(size_t)k*n+j];A[(size_t)k*n+j]=A[(size_t)pk*n+j];A[(size_t)pk*n+j]=t;}
        const double akk=A[(size_t)k*n+k];
        for(int i=k+1;i<n;++i){
            A[(size_t)i*n+k]/=akk;
            const double f=A[(size_t)i*n+k];
            for(int j=k+1;j<n;++j)A[(size_t)i*n+j]-=f*A[(size_t)k*n+j];
        }
    }
    return 1;
}
static int rmt_v2_lu_solve(const double*LU,const int*piv,double*b,int n){
    for(int k=0;k<n;++k)if(piv[k]!=k){double t=b[k];b[k]=b[piv[k]];b[piv[k]]=t;}
    for(int i=0;i<n;++i)for(int j=0;j<i;++j)b[i]-=LU[(size_t)i*n+j]*b[j];
    for(int i=n-1;i>=0;--i){
        for(int j=i+1;j<n;++j)b[i]-=LU[(size_t)i*n+j]*b[j];
        const double d=LU[(size_t)i*n+i];
        if(!(fabs(d)>1e-300)||!isfinite(d))return 0;
        b[i]/=d;
    }
    return 1;
}

static void rmt_v2_build_coarse_plan(const Grid*g,RMTV2CoarsePlan*p,const RMTV2LevelPlan*lp){
    const int sx=lp->sx,sy=lp->sy;
    const int ngx=MIN(sx,g->Nx),ngy=MIN(sy,g->Ny),ng=ngx*ngy;
    p->sx=sx;p->sy=sy;p->ng=ng;
    p->ox=calloc((size_t)ng,sizeof(int));p->oy=calloc((size_t)ng,sizeof(int));
    p->nxg=calloc((size_t)ng,sizeof(int));p->nyg=calloc((size_t)ng,sizeof(int));p->m=calloc((size_t)ng,sizeof(int));
    p->luOffset=calloc((size_t)ng,sizeof(size_t));p->pivOffset=calloc((size_t)ng,sizeof(int));
    if(!p->ox||!p->oy||!p->nxg||!p->nyg||!p->m||!p->luOffset||!p->pivOffset)
        die("alloc FINALVOL2 coarse metadata");

    size_t totalLU=0;int totalU=0;
    for(int gid=0;gid<ng;++gid){
        const int ox=1+gid/ngy,oy=1+gid%ngy;
        const int nxg=(ox<=g->Nx)?1+(g->Nx-ox)/sx:0;
        const int nyg=(oy<=g->Ny)?1+(g->Ny-oy)/sy:0;
        const int m=nxg*nyg;
        if(m<=0||m>RMTV2_DIRECT_MAX)die("FINALVOL2 coarsest block exceeds direct limit");
        p->ox[gid]=ox;p->oy[gid]=oy;p->nxg[gid]=nxg;p->nyg[gid]=nyg;p->m[gid]=m;
        p->luOffset[gid]=totalLU;p->pivOffset[gid]=totalU;
        totalLU+=(size_t)m*(size_t)m;totalU+=m;
    }
    p->totalLU=totalLU;p->totalUnknowns=totalU;
    p->luPool=calloc(totalLU,sizeof(double));p->pivPool=calloc((size_t)totalU,sizeof(int));
    if(!p->luPool||!p->pivPool)die("alloc FINALVOL2 cached LU");

    int failed=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,32) reduction(|:failed) if(ng>8)
#endif
    for(int gid=0;gid<ng;++gid){
        const int ox=p->ox[gid],oy=p->oy[gid],nxg=p->nxg[gid],nyg=p->nyg[gid],m=p->m[gid];
        double*A=p->luPool+p->luOffset[gid];
        int*pv=p->pivPool+p->pivOffset[gid];
        for(int a=0;a<nxg;++a)for(int q=0;q<nyg;++q){
            const int i=ox+a*sx,j=oy+q*sy,row=a*nyg+q;
            const double ap=lp->xAp[i]+lp->yAp[j];
            A[(size_t)row*m+row]=ap;
            if(a>0)A[(size_t)row*m+(row-nyg)]=-lp->xM[i];
            if(a+1<nxg)A[(size_t)row*m+(row+nyg)]=-lp->xP[i];
            if(q>0)A[(size_t)row*m+(row-1)]=-lp->yM[j];
            if(q+1<nyg)A[(size_t)row*m+(row+1)]=-lp->yP[j];
        }
        if(!rmt_v2_lu_factor(A,pv,m))failed=1;
    }
    if(failed)die("FINALVOL2 cached coarsest LU factorization failed");
}

static void rmt_v2_prepare_hierarchy(const Grid*g,const RMTV2Config*cfg){
    RMTV2Workspace*w=rmt_v2_workspace(g);
    if(w->planReady&&w->planRequestedLevels==cfg->requestedLevels&&w->planNCoarsest==cfg->nCoarsest)return;

    for(int l=0;l<RMTV2_MAX_LEVELS;++l)rmt_v2_level_release(&w->level[l]);
    rmt_v2_coarse_release(&w->coarse);
    int lx=rmt_v2_auto_level_1d(g->Nx,cfg->nCoarsest);
    int ly=rmt_v2_auto_level_1d(g->Ny,cfg->nCoarsest);
    if(cfg->requestedLevels>0){
        const int cap=MAX(0,cfg->requestedLevels-1);
        lx=MIN(lx,cap);ly=MIN(ly,cap);
    }
    const int lm=MAX(lx,ly);
    if(lm>=RMTV2_MAX_LEVELS)die("FINALVOL2 hierarchy exceeds level limit");
    for(int lc=0;lc<=lm;++lc){
        const int lX=MIN(lx,lc),lY=MIN(ly,lc);
        const int sx=rmt_v2_ipow3(lX),sy=rmt_v2_ipow3(lY);
        rmt_v2_build_level_plan(g,&w->level[lc],sx,sy,lc<lm);
    }
    rmt_v2_build_coarse_plan(g,&w->coarse,&w->level[lm]);
    w->lx=lx;w->ly=ly;w->lm=lm;w->planReady=1;
    w->planRequestedLevels=cfg->requestedLevels;w->planNCoarsest=cfg->nCoarsest;
    g_rmtV2LUFactorizationsTotal+=w->coarse.ng;
}

static void rmt_v2_direct_coarsest(const Grid*g,double*c,const double*rhs,
                                    const RMTV2CoarsePlan*p,
                                    long long*solveCount,long long*unknownCount){
    int failed=0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,32) reduction(|:failed) if(p->ng>8)
#endif
    for(int gid=0;gid<p->ng;++gid){
        const int ox=p->ox[gid],oy=p->oy[gid],nxg=p->nxg[gid],nyg=p->nyg[gid],m=p->m[gid];
        double b[RMTV2_DIRECT_MAX];
        for(int a=0;a<nxg;++a)for(int q=0;q<nyg;++q){
            const int i=ox+a*p->sx,j=oy+q*p->sy;
            b[a*nyg+q]=rhs[IDX(i,j,g->Ny)];
        }
        const double*LU=p->luPool+p->luOffset[gid];
        const int*pv=p->pivPool+p->pivOffset[gid];
        if(!rmt_v2_lu_solve(LU,pv,b,m)){failed=1;continue;}
        for(int a=0;a<nxg;++a)for(int q=0;q<nyg;++q){
            const int i=ox+a*p->sx,j=oy+q*p->sy;
            c[IDX(i,j,g->Ny)]=b[a*nyg+q];
        }
    }
    if(failed)die("FINALVOL2 cached coarsest LU solve failed");
    if(solveCount)*solveCount+=p->ng;
    if(unknownCount)*unknownCount+=p->totalUnknowns;
}

static void rmt_v2_build_correction(const Grid*g,const double*defect,double*c,
                                    const RMTV2Config*cfg,RMTV2Diag*diag){
    RMTV2Workspace*w=rmt_v2_workspace(g);
    rmt_v2_prepare_hierarchy(g,cfg);
    long long updates=0,directSolves=0,directUnknowns=0;
    memset(c,0,w->n*sizeof(double));

    /* Same defect for every Martynenko mini-cycle level: integrate once. */
    rmt_v2_prefix_build(g,defect,w->prefix);

    for(int lc=w->lm;lc>=0;--lc){
        RMTV2LevelPlan*lp=&w->level[lc];
        rmt_v2_restrict_from_prefix(g,w->prefix,w->rhsLevel,lp->sx,lp->sy);
        if(lc==w->lm)
            rmt_v2_direct_coarsest(g,c,w->rhsLevel,&w->coarse,&directSolves,&directUnknowns);
        else
            updates+=rmt_v2_smooth_level(g,c,w->rhsLevel,lp,MAX(1,cfg->smoothSweeps));
    }

    g_rmtV2PointUpdatesTotal+=updates;
    g_rmtV2DirectSolvesTotal+=directSolves;
    g_rmtV2DirectUnknownsTotal+=directUnknowns;
    g_rmtV2LastLevelX=w->lx;g_rmtV2LastLevelY=w->ly;
    if(diag){
        diag->levelXMax=w->lx;diag->levelYMax=w->ly;diag->pointUpdates=updates;
        diag->directSolves=directSolves;diag->directUnknowns=directUnknowns;
        diag->luFactorizations=w->coarse.ng;
    }
}

#endif
