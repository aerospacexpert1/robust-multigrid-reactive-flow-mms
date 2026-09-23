#!/usr/bin/env python3
from pathlib import Path
import re, shutil
import sympy as sp

OLD=Path('MMS_BENCHMARKS_ABCD_TRUBA_V1')
ROOT=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2')
if ROOT.exists(): shutil.rmtree(ROOT)
if not OLD.exists(): raise SystemExit('V1 generated root missing')
OLD.rename(ROOT)

# Rename package references everywhere.
for p in ROOT.rglob('*'):
    if p.is_file():
        try: s=p.read_text()
        except UnicodeDecodeError: continue
        s=s.replace('MMS_BENCHMARKS_ABCD_TRUBA_V1','MMS_BENCHMARKS_ABCD_TRUBA_V2')
        p.write_text(s)

# ---- Continuous MMS corrections / extra analytic fields ----
x,y=sp.symbols('x y', real=True); pi=sp.pi; X=x/3; Y=y

def cexpr(e): return sp.ccode(sp.simplify(e), standard='C99')
def swfun(name, exprs):
    z=[f'static inline double {name}(double x,double y,int s){{','  switch(s){']
    for i,e in enumerate(exprs,1): z.append(f'    case {i}: return {cexpr(e)};')
    z += ['    default: return 0.0;','  }','}']
    return '\n'.join(z)+'\n'

# A is the production pressure-correction Laplacian. Variable rho enters predictor/corrector,
# but rho cancels from rho*(dt/rho)*grad(p), leaving a constant pressure Laplacian.
ha=ROOT/'Benchmark_A_Pressure_Projection/src/mms_generated.h'
sa=ha.read_text()
ps=[]; rhos=[]; us=[]; vs=[]; sm=[]; dpx=[]; dpy=[]; src=[]
for lev,(amp,eps) in enumerate([(0.8,.08),(1.2,.20),(1.8,.35)],1):
    p=amp*(sp.sin(2*pi*X)*sp.sin(2*pi*Y)+sp.Rational(7,20)*sp.cos(3*pi*X)*sp.cos(pi*Y)+sp.Rational(1,5)*((X-sp.Rational(1,2))**2-(Y-sp.Rational(1,2))**2))
    rho=1+eps*sp.sin(2*pi*X)*sp.cos(2*pi*Y)
    u=.18*(sp.Rational(3,2)-x)+.06*sp.sin(2*pi*X)*sp.sin(2*pi*Y)
    v=.18*(y-sp.Rational(1,2))-.06*sp.cos(2*pi*X)*sp.sin(2*pi*Y)
    ps.append(p);rhos.append(rho);us.append(u);vs.append(v)
    dpx.append(sp.diff(p,x));dpy.append(sp.diff(p,y))
    src.append(-(sp.diff(p,x,2)+sp.diff(p,y,2)))
    sm.append(sp.diff(rho*u,x)+sp.diff(rho*v,y))
# Replace source function; append projection fields.
sa=re.sub(r'static inline double mms_src_p\(double x,double y,int s\)\{.*?\n\}',swfun('mms_src_p',src).rstrip(),sa,flags=re.S)
sa += '\n/* V2 pressure-projection manufactured fields. */\n'
sa += swfun('mms_u_corr',us)+swfun('mms_v_corr',vs)+swfun('mms_mass_src',sm)+swfun('mms_dpx',dpx)+swfun('mms_dpy',dpy)
ha.write_text(sa)

# B: source includes the known pressure gradient; expose its analytic derivatives so the
# numerical momentum operator does not silently omit/cancel that term.
hb=ROOT/'Benchmark_B_Momentum/src/mms_generated.h'; sb=hb.read_text()
p=.12*sp.cos(pi*X)*sp.sin(2*pi*Y)
sb += '\n/* V2 known pressure-gradient forcing for the momentum predictor. */\n'
sb += swfun('mms_pdx',[sp.diff(p,x)]*3)+swfun('mms_pdy',[sp.diff(p,y)]*3)
hb.write_text(sb)

# D production stiffness: A=1e7/1e8/1e9, beta=5, Ta=1000/500/100.
hd=ROOT/'Benchmark_D_Thermochemistry/src/mms_generated.h'; sd=hd.read_text()
# Generator V1 used 350 for S3. Replace only the generated Arrhenius Ta branch.
sd=sd.replace('case 3: return 350.0;', 'case 3: return 100.0;')
hd.write_text(sd)

RMT_RECURSIVE=r'''
/* V2: production-structure RMT. The error equation recursively receives the same
   factor-three, 3x3 shifted-family treatment; it does not terminate each shifted
   grid immediately in a direct coarse RBGS solve as V1 did. */
static void rmt_error_recursive(Sys*s,int level){
    if(level>=5 || s->nx<12 || s->ny<8){rbgs(s,48);return;}
    for(int rep=0;rep<2;rep++){
        double*r=residual_array(s);
        for(int sy=0;sy<3;sy++)for(int sx=0;sx<3;sx++){
            Sys c; make_coarse(s,&c,3);
            for(int jc=1;jc<c.ny-1;jc++)for(int ic=1;ic<c.nx-1;ic++){
                int i0=sx+1+3*(ic-1),j0=sy+1+3*(jc-1);double sum=0;int cnt=0;
                for(int dj=0;dj<3;dj++)for(int di=0;di<3;di++){
                    int fi=i0+di,fj=j0+dj;
                    if(fi>0&&fj>0&&fi<s->nx-1&&fj<s->ny-1){sum+=r[IDX(fi,fj,s->nx)];cnt++;}
                }
                c.b[IDX(ic,jc,c.nx)]=cnt?sum/cnt:0.0;
            }
            correction_bc_zero(&c);
            rmt_error_recursive(&c,level+1);
            for(int jc=1;jc<c.ny-1;jc++)for(int ic=1;ic<c.nx-1;ic++){
                int fi=sx+1+3*(ic-1),fj=sy+1+3*(jc-1);
                if(fi>0&&fj>0&&fi<s->nx-1&&fj<s->ny-1)
                    s->q[IDX(fi,fj,s->nx)] += 0.005*c.q[IDX(ic,jc,c.nx)];
            }
            freesys(&c);
        }
        free(r);
    }
    rbgs(s,3);
}
static void rmt_cycle_ref(Sys*s){
    double before=residual(s); size_t n=(size_t)s->nx*s->ny;
    double*save=malloc(n*sizeof(double)); if(!save)die("rmt save"); memcpy(save,s->q,n*sizeof(double));
    rmt_error_recursive(s,0);
    double after=residual(s);
    /* Production pressure solver also safeguards the multiple correction.  If a
       manufactured operator makes the recursive correction non-descent, retain
       the post-smoothing path rather than accepting a divergent update. */
    if(!isfinite(after)||after>before){memcpy(s->q,save,n*sizeof(double));rbgs(s,6);}
    free(save);
}
'''

TAIL=r'''
typedef struct {double l1,l2,li,t,r;int it;} VR;
static double thermo_h(double T,double YF,double YO,double YP){double YN=1.0-YF-YO-YP;double C=YF*(1800.0+.25*T)+YO*(950.0+.12*T)+YP*(1200.0+.20*T)+YN*(1000.0+.10*T);return C*(T-298.15);}
static double thermo_T(double h,double YF,double YO,double YP,double guess){double T=(guess>200&&guess<2500)?guess:600.0;for(int it=0;it<20;it++){double YN=1.0-YF-YO-YP;double A=YF*1800.0+YO*950.0+YP*1200.0+YN*1000.0;double B=YF*.25+YO*.12+YP*.20+YN*.10;double f=(A+B*T)*(T-298.15)-h;double df=A+B*(2*T-298.15);double d=f/(fabs(df)>1e-12?df:1.0);T-=d;if(T<200)T=200;if(T>2500)T=2500;if(fabs(d)<1e-10*fmax(1.0,T))break;}return T;}
static double thermo_rho(double T,double YF,double YO,double YP){double YN=1.0-YF-YO-YP;double invW=YF/28.01055+YO/31.99880+YP/44.00995+YN/28.01340;double W=1.0/fmax(invW,1e-15);double R=8314.46261815324/W;return 101325.0/fmax(R*T,1.0);}
static double qchem_num(double h,double x,double y,int st){double YF=mms_YF(x,y,st),YO=mms_YO(x,y,st),YP=mms_YP(x,y,st);double T=thermo_T(h,YF,YO,YP,mms_T(x,y,st));double A=mms_arrA(x,y,st),be=mms_arrBeta(x,y,st),Ta=mms_arrTa(x,y,st);double rate=A*pow(fmax(T,200.0),be)*exp(-Ta/fmax(T,200.0))*YF*YF*YO*1e-16;return 2.5e7*rate;}
static void d_rhs_nonlinear(Sys*s,int st){for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++){int k=IDX(i,j,s->nx);double x=i*s->dx,y=j*s->dy;s->b[k]=(mms_src_h(x,y,st)+qchem_num(s->q[k],x,y,st))*s->dx*s->dy;}}
static void thermo_norms(Sys*s,int st,double*T_L2,double*rho_L2,double*T_rel,double*rho_rel){double et=0,er=0,nt=0,nr=0;long n=0;for(int j=0;j<s->ny;j++)for(int i=0;i<s->nx;i++){int k=IDX(i,j,s->nx);double x=i*s->dx,y=j*s->dy,YF=mms_YF(x,y,st),YO=mms_YO(x,y,st),YP=mms_YP(x,y,st);double Te=mms_T(x,y,st),Tn=thermo_T(s->q[k],YF,YO,YP,Te),re=mms_rho(x,y,st),rn=thermo_rho(Tn,YF,YO,YP);et+=(Tn-Te)*(Tn-Te);er+=(rn-re)*(rn-re);nt+=Te*Te;nr+=re*re;n++;}*T_L2=sqrt(et/n);*rho_L2=sqrt(er/n);*T_rel=sqrt(et/(nt+1e-300));*rho_rel=sqrt(er/(nr+1e-300));}
static void onevar_v2(int nx,int ny,int st,char var,const char*solver,VR*o,double*TL2,double*rL2,double*Trel,double*rrel){Sys s;allocsys(&s,nx,ny);assemble(&s,var,st);init_zero(&s);double t=wall();int it=0;if(BENCH_ID=='D'&&var=='h'){for(int outer=0;outer<8;outer++){d_rhs_nonlinear(&s,st);it+=solve(&s,solver,1e-7,24000);if(residual(&s)<1e-7&&outer>=2)break;}}else it=solve(&s,solver,1e-7,24000);o->t=wall()-t;o->it=it;o->r=residual(&s);norms(&s,&o->l1,&o->l2,&o->li);if(BENCH_ID=='D'&&var=='h')thermo_norms(&s,st,TL2,rL2,Trel,rrel);freesys(&s);}
static double projection_continuity(int nx,int ny,int st,const char*solver){Sys s;allocsys(&s,nx,ny);assemble(&s,'a',st);init_zero(&s);solve(&s,solver,1e-7,24000);double dt=.02,sum=0,den=0;long n=0;for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++){double x=i*s.dx,y=j*s.dy;double xe=x+.5*s.dx,xw=x-.5*s.dx,yn=y+.5*s.dy,ys=y-.5*s.dy;double mpe=mms_rho(xe,y,st)*mms_u_corr(xe,y,st)+dt*mms_dpx(xe,y,st);double mpw=mms_rho(xw,y,st)*mms_u_corr(xw,y,st)+dt*mms_dpx(xw,y,st);double mpn=mms_rho(x,yn,st)*mms_v_corr(x,yn,st)+dt*mms_dpy(x,yn,st);double mps=mms_rho(x,ys,st)*mms_v_corr(x,ys,st)+dt*mms_dpy(x,ys,st);int k=IDX(i,j,nx);double pxe=(s.q[k+1]-s.q[k])/s.dx,pxw=(s.q[k]-s.q[k-1])/s.dx,pyn=(s.q[k+nx]-s.q[k])/s.dy,pys=(s.q[k]-s.q[k-nx])/s.dy;double div=((mpe-dt*pxe)-(mpw-dt*pxw))/s.dx+((mpn-dt*pyn)-(mps-dt*pys))/s.dy;double sm=mms_mass_src(x,y,st);sum+=(div-sm)*(div-sm);den+=sm*sm;n++;}freesys(&s);return sqrt(sum/(den+1e-300));}
static void write_real_field_v2(int nx,int ny,int st,char var,const char*solver,const char*path){Sys s;allocsys(&s,nx,ny);assemble(&s,var,st);init_zero(&s);if(BENCH_ID=='D'&&var=='h'){for(int q=0;q<8;q++){d_rhs_nonlinear(&s,st);solve(&s,solver,1e-7,24000);}}else solve(&s,solver,1e-7,24000);FILE*f=fopen(path,"w");if(!f)die("field output");fprintf(f,"x,y,variable,exact,numerical,error\n");for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){int k=IDX(i,j,nx);double x=i*s.dx,y=j*s.dy,ex=s.exact[k],nu=s.q[k];const char*vv="scalar";if(BENCH_ID=='D'&&var=='h'){double YF=mms_YF(x,y,st),YO=mms_YO(x,y,st),YP=mms_YP(x,y,st);ex=mms_T(x,y,st);nu=thermo_T(s.q[k],YF,YO,YP,ex);vv="T";}fprintf(f,"%.17g,%.17g,%s,%.17g,%.17g,%.17g\n",x,y,vv,ex,nu,nu-ex);}fclose(f);freesys(&s);}
int main(int argc,char**argv){int nx=108,ny=36,st=1;const char*solver="SG_RBGS";const char*out="case_summary.csv";const char*fieldout=NULL;for(int a=1;a<argc;a++){if(!strcmp(argv[a],"--Nx")&&a+1<argc)nx=atoi(argv[++a]);else if(!strcmp(argv[a],"--Ny")&&a+1<argc)ny=atoi(argv[++a]);else if(!strcmp(argv[a],"--stiffness")&&a+1<argc){const char*z=argv[++a];st=z[1]-'0';}else if(!strcmp(argv[a],"--solver")&&a+1<argc)solver=argv[++a];else if(!strcmp(argv[a],"--out")&&a+1<argc)out=argv[++a];else if(!strcmp(argv[a],"--field-out")&&a+1<argc)fieldout=argv[++a];}if(nx<12||ny<8||st<1||st>3)die("bad args");if(!finite_sources(st))die("nonfinite MMS source");char vars[3];int nv=1;if(BENCH_ID=='A')vars[0]='a';else if(BENCH_ID=='B'){vars[0]='u';vars[1]='v';nv=2;}else if(BENCH_ID=='C'){vars[0]='f';vars[1]='o';vars[2]='p';nv=3;}else vars[0]='h';VR vr[3]={{0}};double maxL1=0,maxL2=0,maxLi=0,tot=0,rr=0,TL2=0,rL2=0,Trel=0,rrel=0;int its=0;for(int z=0;z<nv;z++){onevar_v2(nx,ny,st,vars[z],solver,&vr[z],&TL2,&rL2,&Trel,&rrel);maxL1=fmax(maxL1,vr[z].l1);maxL2=fmax(maxL2,vr[z].l2);maxLi=fmax(maxLi,vr[z].li);rr=fmax(rr,vr[z].r);tot+=vr[z].t;its+=vr[z].it;}double cont=(BENCH_ID=='A')?projection_continuity(nx,ny,st,solver):0.0;double p2=(BENCH_ID=='A')?vr[0].l2:0,u2=(BENCH_ID=='B')?vr[0].l2:0,v2=(BENCH_ID=='B')?vr[1].l2:0,YF2=(BENCH_ID=='C')?vr[0].l2:0,YO2=(BENCH_ID=='C')?vr[1].l2:0,YP2=(BENCH_ID=='C')?vr[2].l2:0,h2=(BENCH_ID=='D')?vr[0].l2:0;FILE*f=fopen(out,"w");if(!f)die("output");fprintf(f,"benchmark,stiffness,Nx,Ny,threads,solver,converged,iterations,relative_residual,L1,L2,Linf,solve_wall_time,total_relevant_wall_time,continuity_metric,p_L2,u_L2,v_L2,YF_L2,YO_L2,YP_L2,h_L2,T_L2,rho_L2,T_rel_L2,rho_rel_L2\n");fprintf(f,"%c,S%d,%d,%d,%d,%s,%d,%d,%.12e,%.12e,%.12e,%.12e,%.9f,%.9f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",BENCH_ID,st,nx,ny,omp_get_max_threads(),solver,rr<1e-6,its,rr,maxL1,maxL2,maxLi,tot,tot,cont,p2,u2,v2,YF2,YO2,YP2,h2,TL2,rL2,Trel,rrel);fclose(f);if(fieldout)write_real_field_v2(nx,ny,st,vars[0],solver,fieldout);printf("BENCH=%c solver=%s residual=%.3e L2=%.3e cont=%.3e time=%.6f\n",BENCH_ID,solver,rr,maxL2,cont,tot);return rr<1e-6?0:3;}
'''

for d in ROOT.glob('Benchmark_*'):
    c=d/'src/mms_solver.c'; s=c.read_text()
    # Pressure correction operator is constant Laplacian in production low-Mach split.
    s=re.sub(r"if\(BENCH_ID=='A'\)\{double ke=.*?continue;\}","if(BENCH_ID=='A'){double De=dy/dx,Dw=dy/dx,Dn=dx/dy,Ds=dx/dy;s->ae[k]=De;s->aw[k]=Dw;s->an[k]=Dn;s->as[k]=Ds;s->ap[k]=De+Dw+Dn+Ds;s->b[k]=srcvar(var,xx,yy,st)*vol;continue;}",s,flags=re.S)
    # Momentum source contains +grad(p); the known predictor pressure gradient belongs on RHS subtraction.
    s=s.replace("s->b[k]=srcvar(var,xx,yy,st)*vol; (void)rho;(void)ua;(void)va;(void)dif;", "s->b[k]=srcvar(var,xx,yy,st)*vol;if(BENCH_ID=='B')s->b[k]-=(var=='u'?mms_pdx(xx,yy,st):mms_pdy(xx,yy,st))*vol; (void)rho;(void)ua;(void)va;(void)dif;")
    s=re.sub(r'static void rmt_cycle_ref\(Sys\*s\)\{.*?\nstatic int solve',RMT_RECURSIVE+'\nstatic int solve',s,flags=re.S)
    s=re.sub(r'typedef struct \{double l1,l2,li,t,r;int it;\} VR;.*\Z',TAIL,s,flags=re.S)
    c.write_text(s)
    # Real campaign mesh is mandatory smoke, not the old 36x12 convenience mesh.
    sm=d/'smoke_test.sh'
    if sm.exists():
        z=sm.read_text().replace('--Nx 36 --Ny 12','--Nx 108 --Ny 36').replace('Nx=36','Nx=108').replace('Ny=12','Ny=36')
        sm.write_text(z)

# Documentation: state exactly what is production and what is solver-extension.
(ROOT/'V2_METHOD_CHANGES.md').write_text('''# V2 production-fidelity corrections\n\n- **A** uses the production low-Mach pressure-correction Laplacian and reconstructs predictor mass flux, numerical pressure correction, corrected flux and a discrete continuity defect. Density is manufactured in the predictor/corrector chain; it is not used to invent a variable-k Poisson equation.\n- **B** keeps the known pressure gradient in the momentum predictor and verifies the conservative advection-diffusion operator for both u and v.\n- **C** verifies chemistry-off YF/YO/YP transport.\n- **D** is an energy/thermochemistry isolation MMS: YF/YO/YP are prescribed manufactured composition fields (their transport is verified in C), while sensible enthalpy is the numerical unknown. Each Picard update recomputes T(h,Y), low-Mach perfect-gas density and Arrhenius heat release. Raw h error and derived T/rho absolute and relative errors are reported.\n- RMT3H is recursively factor-three on all eligible coarse levels with the nine shifted (sx,sy) families, two repeated multiple corrections, omega=0.005 and post-smoothing, matching the structural invariants of the uploaded production RMT implementation.\n- A is the direct production pressure-solver comparison. Applying the five solver kernels to B/C/D is explicitly an **operator-solver extension** for numerical-method study, not a claim that the production opposed-flow executable switches those transport blocks among five solver families.\n\nContinuous MMS sources are derived before discretization; no `b_h=A_h phi_exact` construction is used.\n''')
for d in ROOT.glob('Benchmark_*'):
    docs=d/'docs'; docs.mkdir(exist_ok=True)
    (docs/'V2_FIDELITY_NOTE.md').write_text((ROOT/'V2_METHOD_CHANGES.md').read_text())
print('V2_PATCH_COMPLETE')
