/* MMS V3 RMT kernel: production-structure adaptation of RMT FINAL / FINALVOL2.
 * This file is injected into the generated MMS solver by tools/v3_patch.py.
 */

typedef struct {
    int nx, ny;
    size_t n;
    double *corr;
    double *rhs;
    double *prefix;
} RMTV3Workspace;

static RMTV3Workspace g_rmtv3 = {0};

static void rmt_v3_workspace_release(void){
    free(g_rmtv3.corr);
    free(g_rmtv3.rhs);
    free(g_rmtv3.prefix);
    memset(&g_rmtv3,0,sizeof(g_rmtv3));
}

static RMTV3Workspace *rmt_v3_workspace(Sys *s){
    size_t n=(size_t)s->nx*(size_t)s->ny;
    if(g_rmtv3.nx==s->nx && g_rmtv3.ny==s->ny &&
       g_rmtv3.corr && g_rmtv3.rhs && g_rmtv3.prefix)
        return &g_rmtv3;
    rmt_v3_workspace_release();
    g_rmtv3.nx=s->nx; g_rmtv3.ny=s->ny; g_rmtv3.n=n;
    g_rmtv3.corr=calloc(n,sizeof(double));
    g_rmtv3.rhs=calloc(n,sizeof(double));
    g_rmtv3.prefix=calloc((size_t)(s->nx+1)*(size_t)(s->ny+1),sizeof(double));
    if(!g_rmtv3.corr||!g_rmtv3.rhs||!g_rmtv3.prefix) die("RMT V3 workspace allocation");
    return &g_rmtv3;
}

static int rmt_v3_ipow3(int l){
    int s=1;
    while(l-- > 0) s*=3;
    return s;
}

static int rmt_v3_auto_level_1d(int nInterior,int nCoarsest){
    int l=0,s=1;
    if(nCoarsest<2)nCoarsest=2;
    while((nInterior+s-1)/s > nCoarsest){
        if(s > 100000000/3) break;
        s*=3; l++;
        if(l>=10) break;
    }
    return l;
}

static void rmt_v3_prefix_build_density(Sys*s,const double*r,double*p){
    int nx=s->nx,ny=s->ny,pitch=ny+1;
    const double vol=s->dx*s->dy;
    memset(p,0,(size_t)(nx+1)*(size_t)(ny+1)*sizeof(double));
    for(int i=1;i<nx-1;i++){
        double row=0.0;
        for(int j=1;j<ny-1;j++){
            row += r[IDX(i,j,nx)]/vol;
            p[(size_t)i*pitch+j]=p[(size_t)(i-1)*pitch+j]+row;
        }
    }
}

static inline void rmt_v3_cv_bounds_1d(int n,int idx,int stride,int*lo,int*hi){
    const int last=n-2;
    const int half=(stride-1)/2;
    const int hasMinus=(idx-stride>=1);
    const int hasPlus=(idx+stride<=last);
    *lo=hasMinus?idx-half:1;
    *hi=hasPlus?idx+half:last;
    if(*lo<1)*lo=1;
    if(*hi>last)*hi=last;
}

static double rmt_v3_cv_average(Sys*s,const double*p,int i,int j,int sx,int sy){
    int i1,i2,j1,j2,pitch=s->ny+1;
    rmt_v3_cv_bounds_1d(s->nx,i,sx,&i1,&i2);
    rmt_v3_cv_bounds_1d(s->ny,j,sy,&j1,&j2);
    double sum=p[(size_t)i2*pitch+j2]
              -p[(size_t)(i1-1)*pitch+j2]
              -p[(size_t)i2*pitch+j1-1]
              +p[(size_t)(i1-1)*pitch+j1-1];
    double cnt=(double)(i2-i1+1)*(double)(j2-j1+1);
    return cnt>0.0?sum/cnt:0.0;
}

static void rmt_v3_fill_rhs(Sys*s,const double*prefix,double*rhs,int sx,int sy){
    int nx=s->nx,ny=s->ny;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++)
        rhs[IDX(i,j,nx)]=rmt_v3_cv_average(s,prefix,i,j,sx,sy);
}

static inline void rmt_v3_laplace_dim(int n,int idx,int stride,double h,
                                      double*ap,double*am,double*az){
    int last=n-2;
    int hm=(idx-stride>=1), hp=(idx+stride<=last);
    double H=stride*h;
    *am=0.0; *az=0.0;
    if(hm&&hp){
        double a=1.0/(H*H);
        *ap+=2.0*a; *am=a; *az=a; return;
    }
    if(!hm&&hp){
        double d=idx*h;
        double V=d+0.5*H;
        double an=1.0/(H*V), ab=1.0/(d*V);
        *ap+=an+ab; *az=an; return;
    }
    if(hm&&!hp){
        double d=(n-1-idx)*h;
        double V=d+0.5*H;
        double an=1.0/(H*V), ab=1.0/(d*V);
        *ap+=an+ab; *am=an; return;
    }
    {
        double L=(n-1)*h;
        double d0=idx*h, d1=(n-1-idx)*h;
        if(d0>0.0)*ap+=1.0/(d0*L);
        if(d1>0.0)*ap+=1.0/(d1*L);
    }
}

static inline double rmt_v3_clip(double x,double a,double b){
    return x<a?a:(x>b?b:x);
}

static void rmt_v3_coeff(Sys*s,int i,int j,int sx,int sy,
                         double*ap,double*aw,double*ae,double*as,double*an){
    *ap=*aw=*ae=*as=*an=0.0;
    if(BENCH_ID=='A'){
        rmt_v3_laplace_dim(s->nx,i,sx,s->dx,ap,aw,ae);
        rmt_v3_laplace_dim(s->ny,j,sy,s->dy,ap,as,an);
        return;
    }

    double Hx=sx*s->dx, Hy=sy*s->dy;
    double x=i*s->dx, y=j*s->dy;
    double xe=rmt_v3_clip(x+0.5*Hx,0.0,3.0);
    double xw=rmt_v3_clip(x-0.5*Hx,0.0,3.0);
    double yn=rmt_v3_clip(y+0.5*Hy,0.0,1.0);
    double ys=rmt_v3_clip(y-0.5*Hy,0.0,1.0);

    double re=mms_rho(xe,y,s->st), rw=mms_rho(xw,y,s->st);
    double rn=mms_rho(x,yn,s->st), rs=mms_rho(x,ys,s->st);
    double ue=mms_uadv(xe,y,s->st), uw=mms_uadv(xw,y,s->st);
    double vn=mms_vadv(x,yn,s->st), vs=mms_vadv(x,ys,s->st);
    double de=mms_diff(xe,y,s->st), dw=mms_diff(xw,y,s->st);
    double dn=mms_diff(x,yn,s->st), ds=mms_diff(x,ys,s->st);

    double Fe=re*ue/Hx, Fw=rw*uw/Hx, Fn=rn*vn/Hy, Fs=rs*vs/Hy;
    double De=de/(Hx*Hx), Dw=dw/(Hx*Hx), Dn=dn/(Hy*Hy), Ds=ds/(Hy*Hy);

    *ae=De+fmax(-Fe,0.0);
    *aw=Dw+fmax( Fw,0.0);
    *an=Dn+fmax(-Fn,0.0);
    *as=Ds+fmax( Fs,0.0);
    *ap=*ae+*aw+*an+*as+(Fe-Fw+Fn-Fs);
    if(*ap<1e-30)*ap=1e-30;
}

static inline double rmt_v3_update_point(Sys*s,double*c,const double*rhs,
                                         int i,int j,int sx,int sy){
    int nx=s->nx,ny=s->ny;
    double ap,aw,ae,as,an;
    rmt_v3_coeff(s,i,j,sx,sy,&ap,&aw,&ae,&as,&an);
    double sum=0.0;
    if(i-sx>=1)     sum+=aw*c[IDX(i-sx,j,nx)];
    if(i+sx<=nx-2)  sum+=ae*c[IDX(i+sx,j,nx)];
    if(j-sy>=1)     sum+=as*c[IDX(i,j-sy,nx)];
    if(j+sy<=ny-2)  sum+=an*c[IDX(i,j+sy,nx)];
    return (rhs[IDX(i,j,nx)]+sum)/ap;
}

static long long rmt_v3_smooth_level(Sys*s,double*c,const double*rhs,
                                     int sx,int sy,int sweeps){
    int nx=s->nx,ny=s->ny;
    long long updates=0;
    long long nxy=(long long)(nx-2)*(long long)(ny-2);
    for(int sw=0;sw<sweeps;sw++){
        for(int color=0;color<2;color++){
            #pragma omp parallel for collapse(2) schedule(static) reduction(+:updates) if(nxy>2048)
            for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++){
                int ox=(i-1)%sx, oy=(j-1)%sy;
                int qi=(i-1-ox)/sx, qj=(j-1-oy)/sy;
                if(((qi+qj)&1)!=color)continue;
                c[IDX(i,j,nx)]=rmt_v3_update_point(s,c,rhs,i,j,sx,sy);
                updates++;
            }
        }
    }
    return updates;
}

#define RMT_V3_DIRECT_MAX 64

static int rmt_v3_dense_solve(double*A,double*b,int n){
    for(int k=0;k<n;k++){
        int piv=k; double pv=fabs(A[(size_t)k*n+k]);
        for(int i=k+1;i<n;i++){
            double v=fabs(A[(size_t)i*n+k]);
            if(v>pv){pv=v;piv=i;}
        }
        if(!(pv>1e-300)||!isfinite(pv))return 0;
        if(piv!=k){
            for(int j=k;j<n;j++){
                double t=A[(size_t)k*n+j];
                A[(size_t)k*n+j]=A[(size_t)piv*n+j];
                A[(size_t)piv*n+j]=t;
            }
            double t=b[k]; b[k]=b[piv]; b[piv]=t;
        }
        double akk=A[(size_t)k*n+k];
        for(int i=k+1;i<n;i++){
            double f=A[(size_t)i*n+k]/akk;
            A[(size_t)i*n+k]=0.0;
            for(int j=k+1;j<n;j++)A[(size_t)i*n+j]-=f*A[(size_t)k*n+j];
            b[i]-=f*b[k];
        }
    }
    for(int i=n-1;i>=0;i--){
        double z=b[i];
        for(int j=i+1;j<n;j++)z-=A[(size_t)i*n+j]*b[j];
        double aii=A[(size_t)i*n+i];
        if(!(fabs(aii)>1e-300)||!isfinite(aii))return 0;
        b[i]=z/aii;
    }
    return 1;
}

static void rmt_v3_direct_coarsest(Sys*s,double*c,const double*rhs,
                                   int sx,int sy){
    int nx=s->nx,ny=s->ny,nxi=nx-2,nyi=ny-2;
    int ngx=sx<nxi?sx:nxi, ngy=sy<nyi?sy:nyi;
    int ng=ngx*ngy,failed=0;

    #pragma omp parallel for schedule(dynamic,32) reduction(|:failed) if(ng>8)
    for(int gid=0;gid<ng;gid++){
        int ox=1+gid/ngy, oy=1+gid%ngy;
        int nxg=(ox<=nx-2)?1+(nx-2-ox)/sx:0;
        int nyg=(oy<=ny-2)?1+(ny-2-oy)/sy:0;
        int m=nxg*nyg;
        if(m<=0)continue;
        if(m>RMT_V3_DIRECT_MAX){failed=1;continue;}

        double A[RMT_V3_DIRECT_MAX*RMT_V3_DIRECT_MAX];
        double b[RMT_V3_DIRECT_MAX];
        memset(A,0,(size_t)m*(size_t)m*sizeof(double));

        for(int a=0;a<nxg;a++)for(int q=0;q<nyg;q++){
            int i=ox+a*sx,j=oy+q*sy,row=a*nyg+q;
            double ap,aw,ae,as,an;
            rmt_v3_coeff(s,i,j,sx,sy,&ap,&aw,&ae,&as,&an);
            A[(size_t)row*m+row]=ap;
            if(a>0)       A[(size_t)row*m+(row-nyg)]=-aw;
            if(a+1<nxg)   A[(size_t)row*m+(row+nyg)]=-ae;
            if(q>0)       A[(size_t)row*m+(row-1)]=-as;
            if(q+1<nyg)   A[(size_t)row*m+(row+1)]=-an;
            b[row]=rhs[IDX(i,j,nx)];
        }

        if(!rmt_v3_dense_solve(A,b,m)){failed=1;continue;}
        for(int a=0;a<nxg;a++)for(int q=0;q<nyg;q++){
            int i=ox+a*sx,j=oy+q*sy;
            c[IDX(i,j,nx)]=b[a*nyg+q];
        }
    }
    if(failed)die("RMT V3 direct coarsest solve failed");
}

static void rmt_v3_build_correction(Sys*s,const double*defect,double*c){
    RMTV3Workspace*w=rmt_v3_workspace(s);
    int nxi=s->nx-2,nyi=s->ny-2;
    int lx=rmt_v3_auto_level_1d(nxi,5);
    int ly=rmt_v3_auto_level_1d(nyi,5);
    int lm=lx>ly?lx:ly;

    memset(c,0,w->n*sizeof(double));
    rmt_v3_prefix_build_density(s,defect,w->prefix);

    for(int lc=lm;lc>=0;--lc){
        int lX=lx<lc?lx:lc, lY=ly<lc?ly:lc;
        int sx=rmt_v3_ipow3(lX), sy=rmt_v3_ipow3(lY);
        rmt_v3_fill_rhs(s,w->prefix,w->rhs,sx,sy);
        if(lc==lm)
            rmt_v3_direct_coarsest(s,c,w->rhs,sx,sy);
        else
            rmt_v3_smooth_level(s,c,w->rhs,sx,sy,16);
    }
}

static void rmt_cycle_ref(Sys*s){
    RMTV3Workspace*w=rmt_v3_workspace(s);
    double*r=residual_array(s);
    rmt_v3_build_correction(s,r,w->corr);
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)
        s->q[IDX(i,j,s->nx)] += w->corr[IDX(i,j,s->nx)];
    free(r);
}
