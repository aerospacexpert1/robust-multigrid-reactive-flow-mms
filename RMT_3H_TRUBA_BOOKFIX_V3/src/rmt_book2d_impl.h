#ifndef RMT_BOOK2D_IMPL_H
#define RMT_BOOK2D_IMPL_H

typedef struct { int requestedLevels,smoothSweeps,coarseSweeps,nCoarsest; } RMTBookConfig;
typedef struct { int levelXMax,levelYMax; long long pointUpdates; } RMTBookDiag;
typedef struct { int nx,ny; size_t n; double *remaining,*rhsLevel,*prefix,*corr,*res; } RMTBookWorkspace;

static RMTBookWorkspace g_rmtBookWS={0,0,0,NULL,NULL,NULL,NULL,NULL};
static long long g_rmtBookPointUpdatesTotal=0,g_rmtBookFullAccepts=0,g_rmtBookBacktracks=0,g_rmtBookRejects=0;
static double g_rmtBookLastAcceptedOmega=0.0;
static int g_rmtBookLastLevelX=0,g_rmtBookLastLevelY=0;

static void rmt_book_workspace_release(void){
 free(g_rmtBookWS.remaining);free(g_rmtBookWS.rhsLevel);free(g_rmtBookWS.prefix);free(g_rmtBookWS.corr);free(g_rmtBookWS.res);
 memset(&g_rmtBookWS,0,sizeof(g_rmtBookWS));
}
static RMTBookWorkspace *rmt_book_workspace(const Grid *g){
 size_t n=(size_t)(g->Nx+2)*(g->Ny+2),np=(size_t)(g->Nx+1)*(g->Ny+1);
 if(g_rmtBookWS.nx==g->Nx&&g_rmtBookWS.ny==g->Ny&&g_rmtBookWS.corr)return &g_rmtBookWS;
 rmt_book_workspace_release();g_rmtBookWS.nx=g->Nx;g_rmtBookWS.ny=g->Ny;g_rmtBookWS.n=n;
 g_rmtBookWS.remaining=calloc(n,sizeof(double));g_rmtBookWS.rhsLevel=calloc(n,sizeof(double));g_rmtBookWS.prefix=calloc(np,sizeof(double));g_rmtBookWS.corr=calloc(n,sizeof(double));g_rmtBookWS.res=calloc(n,sizeof(double));
 if(!g_rmtBookWS.remaining||!g_rmtBookWS.rhsLevel||!g_rmtBookWS.prefix||!g_rmtBookWS.corr||!g_rmtBookWS.res)die("alloc RMT book workspace");return &g_rmtBookWS;
}
static int rmt_book_auto_level_1d(int n,int nc){int l=0,s=1;if(nc<2)nc=2;while((n+s-1)/s>nc){s*=3;++l;if(l>=12)break;}return l;}
static int rmt_book_ipow3(int l){int s=1;while(l-->0)s*=3;return s;}

static void rmt_book_apply_fine_A(const Grid *g,const double *x,double *Ax){
 int nx=g->Nx,ny=g->Ny;double ax=1/(g->dx*g->dx),ay=1/(g->dy*g->dy);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
 for(int i=1;i<=nx;++i)for(int j=1;j<=ny;++j){double ap=0,sum=0;if(i>1){ap+=ax;sum+=ax*x[IDX(i-1,j,ny)];}else ap+=2*ax;if(i<nx){ap+=ax;sum+=ax*x[IDX(i+1,j,ny)];}else ap+=2*ax;if(j>1){ap+=ay;sum+=ay*x[IDX(i,j-1,ny)];}if(j<ny){ap+=ay;sum+=ay*x[IDX(i,j+1,ny)];}Ax[IDX(i,j,ny)]=ap*x[IDX(i,j,ny)]-sum;}
}
static void rmt_book_remaining_residual(const Grid *g,const double*b,const double*c,double*r){rmt_book_apply_fine_A(g,c,r);for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)r[IDX(i,j,g->Ny)]=b[IDX(i,j,g->Ny)]-r[IDX(i,j,g->Ny)];}
static void rmt_book_prefix_build(const Grid*g,const double*a,double*p){int ny=g->Ny,pitch=ny+1;memset(p,0,(size_t)(g->Nx+1)*(ny+1)*sizeof(double));for(int i=1;i<=g->Nx;++i){double row=0;for(int j=1;j<=ny;++j){row+=a[IDX(i,j,ny)];p[(size_t)i*pitch+j]=p[(size_t)(i-1)*pitch+j]+row;}}}
static double rmt_book_box_average(const Grid*g,const double*p,int i,int j,int sx,int sy){int pitch=g->Ny+1,hx=(sx-1)/2,hy=(sy-1)/2,i1=MAX(1,i-hx),i2=MIN(g->Nx,i+hx),j1=MAX(1,j-hy),j2=MIN(g->Ny,j+hy);double sum=p[(size_t)i2*pitch+j2]-p[(size_t)(i1-1)*pitch+j2]-p[(size_t)i2*pitch+j1-1]+p[(size_t)(i1-1)*pitch+j1-1],cnt=(double)(i2-i1+1)*(j2-j1+1);return cnt?sum/cnt:0;}
static void rmt_book_restrict_level(const Grid*g,const double*r,double*rhs,double*p,int sx,int sy){rmt_book_prefix_build(g,r,p);for(int i=1;i<=g->Nx;++i)for(int j=1;j<=g->Ny;++j)rhs[IDX(i,j,g->Ny)]=rmt_book_box_average(g,p,i,j,sx,sy);}

/* Exact elimination of Martynenko virtual endpoints, Eqs. (2.17)-(2.20). */
static inline void rmt_book_dim_coeff(double H,double dist,int hm,int hp,int neu,double*ap,double*am,double*az){
 double q=1/(H*H);*am=*az=0;if(hm&&hp){*ap+=2*q;*am=*az=q;return;}double xi=MAX(dist/H,1e-12),alpha,beta;if(neu){alpha=4*xi/(2*xi+1);beta=-(2*xi-1)/(2*xi+1);}else{alpha=(2*xi-1)/xi;beta=-(xi-1)/(xi+1);}*ap+=(2-alpha)*q;double an=(1+beta)*q;if(!hm&&hp)*az=an;else if(hm&&!hp)*am=an;
}
static inline double rmt_book_update_point(const Grid*g,double*c,const double*rhs,int i,int j,int sx,int sy){
 int nx=g->Nx,ny=g->Ny,hw=i-sx>=1,he=i+sx<=nx,hs=j-sy>=1,hn=j+sy<=ny;double Hx=sx*g->dx,Hy=sy*g->dy,ap=0,aw=0,ae=0,as=0,an=0;double dx=!hw?(i-.5)*g->dx:(nx-i+.5)*g->dx,dy=!hs?(j-.5)*g->dy:(ny-j+.5)*g->dy;rmt_book_dim_coeff(Hx,dx,hw,he,0,&ap,&aw,&ae);rmt_book_dim_coeff(Hy,dy,hs,hn,1,&ap,&as,&an);double sum=0;if(hw)sum+=aw*c[IDX(i-sx,j,ny)];if(he)sum+=ae*c[IDX(i+sx,j,ny)];if(hs)sum+=as*c[IDX(i,j-sy,ny)];if(hn)sum+=an*c[IDX(i,j+sy,ny)];return(rhs[IDX(i,j,ny)]+sum)/MAX(ap,1e-30);
}
static long long rmt_book_smooth_level(const Grid*g,double*c,const double*rhs,int sx,int sy,int sweeps){int ngx=MIN(sx,g->Nx),ngy=MIN(sy,g->Ny),ng=ngx*ngy;long long u=0;for(int sw=0;sw<sweeps;++sw){
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:u) if(ng>1)
#endif
 for(int gid=0;gid<ng;++gid){int ox=1+gid/ngy,oy=1+gid%ngy;for(int i=ox;i<=g->Nx;i+=sx)for(int j=oy;j<=g->Ny;j+=sy){c[IDX(i,j,g->Ny)]=rmt_book_update_point(g,c,rhs,i,j,sx,sy);++u;}}}return u;}

static void rmt_book_build_correction(const Grid*g,const double*b,double*c,const RMTBookConfig*cfg,RMTBookDiag*d){
 RMTBookWorkspace*w=rmt_book_workspace(g);int lx=rmt_book_auto_level_1d(g->Nx,cfg->nCoarsest),ly=rmt_book_auto_level_1d(g->Ny,cfg->nCoarsest);if(cfg->requestedLevels>0){int cap=MAX(0,cfg->requestedLevels-1);lx=MIN(lx,cap);ly=MIN(ly,cap);}int lm=MAX(lx,ly);long long updates=0;memset(c,0,w->n*sizeof(double));
 /* Linear sawtooth: the same defect equation is represented on every level.  The coarse
    correction is retained as the initial guess and refined on successively finer levels. */
 for(int lc=lm;lc>=0;--lc){int lX=MIN(lx,lc),lY=MIN(ly,lc),sx=rmt_book_ipow3(lX),sy=rmt_book_ipow3(lY);rmt_book_restrict_level(g,b,w->rhsLevel,w->prefix,sx,sy);int sw=lc==lm?MAX(1,cfg->coarseSweeps):MAX(1,cfg->smoothSweeps);updates+=rmt_book_smooth_level(g,c,w->rhsLevel,sx,sy,sw);}
 if(d){d->levelXMax=lx;d->levelYMax=ly;d->pointUpdates=updates;}g_rmtBookPointUpdatesTotal+=updates;g_rmtBookLastLevelX=lx;g_rmtBookLastLevelY=ly;
}
#endif
