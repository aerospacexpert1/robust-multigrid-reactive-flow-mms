/* MMS V4 conventional-MG kernel.
 *
 * Fair-comparison rules:
 *   - the same fine-grid discrete system, initial guess and stopping tolerance are
 *     used by SG_RBGS, MG2V, MG2W, MG3V and RMT3H;
 *   - conventional MG residual transfer is dimensionally consistent with the
 *     integrated finite-volume fine-grid equations;
 *   - coarse grids are exactly nested in physical coordinates;
 *   - transfer operators are tensor-product full weighting / linear interpolation;
 *   - no solver-specific line search, rollback or hidden RBGS fallback is used.
 *
 * The logical campaign Nx,Ny values are interval/cell counts.  V4 allocates
 * Nx+1 by Ny+1 nodal points, so the campaign sizes 108x36, 216x72, ... are
 * exactly divisible by both factor two and factor three.
 */

static double *mg4_residual_density(Sys *s){
    size_t n=(size_t)s->nx*(size_t)s->ny;
    double *r=calloc(n,sizeof(double));
    if(!r) die("MG V4 residual allocation");
    const double vol=s->dx*s->dy;
    int nx=s->nx,ny=s->ny;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++){
        int k=IDX(i,j,nx);
        double Aq=s->ap[k]*s->q[k]
                 -s->ae[k]*s->q[k+1]-s->aw[k]*s->q[k-1]
                 -s->an[k]*s->q[k+nx]-s->as[k]*s->q[k-nx];
        r[k]=(s->b[k]-Aq)/vol;
    }
    return r;
}

static void mg4_zero_correction(Sys *s){
    size_t n=(size_t)s->nx*(size_t)s->ny;
    memset(s->q,0,n*sizeof(double));
    memset(s->b,0,n*sizeof(double));
}

static void mg4_make_coarse(Sys *f,Sys *c,int ratio){
    if(ratio!=2 && ratio!=3) die("MG V4 unsupported coarsening ratio");
    if(((f->nx-1)%ratio)!=0 || ((f->ny-1)%ratio)!=0)
        die("MG V4 requires exactly nested interval counts");
    int nx=(f->nx-1)/ratio+1;
    int ny=(f->ny-1)/ratio+1;
    if(nx<3 || ny<3) die("MG V4 coarse grid too small");
    allocsys(c,nx,ny);
    assemble(c,f->var,f->st);
    mg4_zero_correction(c);
}

static inline double mg4_rweight_1d(int ratio,int off){
    int a=off<0?-off:off;
    if(ratio==2){
        if(a==0)return 0.5;
        if(a==1)return 0.25;
        return 0.0;
    }
    if(a==0)return 1.0/3.0;
    if(a==1)return 2.0/9.0;
    if(a==2)return 1.0/9.0;
    return 0.0;
}

/* rFine is residual density.  Full weighting produces a coarse residual
 * density; multiplication by the coarse control-volume area converts it to
 * the integrated RHS expected by the rediscretized coarse Sys operator.
 */
static void mg4_restrict_density(Sys *f,const double *rFine,Sys *c,int ratio){
    const double cvol=c->dx*c->dy;
    const int rad=ratio-1;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int jc=1;jc<c->ny-1;jc++)for(int ic=1;ic<c->nx-1;ic++){
        int fi=ratio*ic,fj=ratio*jc;
        double z=0.0,ws=0.0;
        for(int dj=-rad;dj<=rad;dj++){
            double wy=mg4_rweight_1d(ratio,dj);
            if(wy==0.0)continue;
            for(int di=-rad;di<=rad;di++){
                double wx=mg4_rweight_1d(ratio,di);
                if(wx==0.0)continue;
                int ii=fi+di,jj=fj+dj;
                if(ii<1||ii>=f->nx-1||jj<1||jj>=f->ny-1)continue;
                double w=wx*wy;
                z+=w*rFine[IDX(ii,jj,f->nx)];
                ws+=w;
            }
        }
        if(ws>0.0) z/=ws;
        c->b[IDX(ic,jc,c->nx)]=z*cvol;
    }
}

/* Exact nested-grid bilinear prolongation. */
static void mg4_prolong_add(Sys *c,Sys *f,int ratio){
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j=1;j<f->ny-1;j++)for(int i=1;i<f->nx-1;i++){
        int ic=i/ratio,jc=j/ratio;
        int ri=i%ratio,rj=j%ratio;
        double tx=(double)ri/(double)ratio;
        double ty=(double)rj/(double)ratio;
        int ie=ic+1,jn=jc+1;
        if(ie>=c->nx)ie=c->nx-1;
        if(jn>=c->ny)jn=c->ny-1;
        double c00=c->q[IDX(ic,jc,c->nx)];
        double c10=c->q[IDX(ie,jc,c->nx)];
        double c01=c->q[IDX(ic,jn,c->nx)];
        double c11=c->q[IDX(ie,jn,c->nx)];
        double e=(1.0-tx)*(1.0-ty)*c00
                +tx*(1.0-ty)*c10
                +(1.0-tx)*ty*c01
                +tx*ty*c11;
        f->q[IDX(i,j,f->nx)]+=e;
    }
}

static void mg4_cycle(Sys *s,int ratio,int levels,int visits){
    rbgs(s,3);

    double *r=mg4_residual_density(s);
    Sys c;
    mg4_make_coarse(s,&c,ratio);
    mg4_restrict_density(s,r,&c,ratio);
    free(r);

    if(levels>1 && c.nx>=7 && c.ny>=7){
        for(int q=0;q<visits;q++)
            mg4_cycle(&c,ratio,levels-1,visits);
    }else{
        rbgs(&c,40);
    }

    mg4_prolong_add(&c,s,ratio);
    freesys(&c);
    rbgs(s,3);
}
