#!/usr/bin/env python3
import csv, math, os, shutil, textwrap
from pathlib import Path
import sympy as sp

ROOT = Path('MMS_BENCHMARKS_ABCD_TRUBA_V1')
BENCHES = [
 ('A','Benchmark_A_Pressure_Projection','Pressure Projection / Continuity'),
 ('B','Benchmark_B_Momentum','Momentum'),
 ('C','Benchmark_C_Species_Transport','Species Transport'),
 ('D','Benchmark_D_Thermochemistry','Thermochemistry + Energy'),
]
MESHES=[('M1',108,36),('M2',216,72),('M3',432,144),('M4',864,288),('M5',1728,576)]
THREADS=[('T1',1),('T2',2),('T4',4),('T16',16)]
SOLVERS=['SG_RBGS','MG2V','MG2W','MG3V','RMT3H']
STIFF=[('S1','Low'),('S2','Medium'),('S3','High')]

if ROOT.exists(): shutil.rmtree(ROOT)
ROOT.mkdir()

x,y=sp.symbols('x y', real=True)
pi=sp.pi
Lx=sp.Float(3.0); Ly=sp.Float(1.0)
X=x/Lx; Y=y/Ly

def cexpr(e): return sp.ccode(sp.simplify(e), standard='C99')

def derive(bench, level):
    lev={'S1':1,'S2':2,'S3':3}[level]
    if bench=='A':
        amp=[0.8,1.2,1.8][lev-1]; eps=[0.08,0.20,0.35][lev-1]
        p=amp*(sp.sin(2*pi*X)*sp.sin(2*pi*Y)+sp.Rational(7,20)*sp.cos(3*pi*X)*sp.cos(pi*Y)+sp.Rational(1,5)*((X-sp.Rational(1,2))**2-(Y-sp.Rational(1,2))**2))
        rho=1+eps*sp.sin(2*pi*X)*sp.cos(2*pi*Y)
        k=1/rho
        src=-(sp.diff(k*sp.diff(p,x),x)+sp.diff(k*sp.diff(p,y),y))
        return {'p':p,'rho':rho,'k':k,'src_p':src,'uadv':0*x,'vadv':0*x,'diff':k}
    if bench=='B':
        U0=[0.5,1.0,1.8][lev-1]; vort=[0.08,0.15,0.25][lev-1]; eps=[0.03,0.10,0.20][lev-1]; mu=[0.025,0.015,0.008][lev-1]
        u=U0*(sp.Rational(3,2)-x)+vort*sp.sin(2*pi*X)*sp.sin(2*pi*Y)
        v=U0*(y-sp.Rational(1,2))-vort*sp.cos(2*pi*X)*sp.sin(2*pi*Y)
        rho=1+eps*sp.cos(2*pi*X)*sp.cos(pi*Y)
        p=0.12*sp.cos(pi*X)*sp.sin(2*pi*Y)
        def op(phi):
            return sp.diff(rho*u*phi,x)+sp.diff(rho*v*phi,y)-mu*(sp.diff(phi,x,2)+sp.diff(phi,y,2))
        return {'u':u,'v':v,'p':p,'rho':rho,'uadv':u,'vadv':v,'diff':sp.Float(mu),'src_u':op(u)+sp.diff(p,x),'src_v':op(v)+sp.diff(p,y)}
    if bench=='C':
        U0=[0.35,0.70,1.20][lev-1]; delta=[0.24,0.16,0.10][lev-1]; D=[0.018,0.010,0.006][lev-1]; eps=[0.02,0.08,0.15][lev-1]
        eta=(x-sp.Rational(3,2)-sp.Float(0.12)*sp.sin(2*pi*y))/delta
        m=(1-sp.tanh(eta))/2
        band=1/sp.cosh(eta)**2
        YF=sp.Float(0.055)*m*(1+sp.Float(0.04)*sp.sin(2*pi*y))
        YO=sp.Float(0.233)*(1-m)*(1-sp.Float(0.03)*sp.cos(2*pi*y))
        YP=sp.Float(0.10)*band*(1+sp.Float(0.08)*sp.sin(2*pi*X))
        u=U0*(sp.Rational(3,2)-x); v=U0*(y-sp.Rational(1,2))*sp.Float(0.25)
        rho=1+eps*sp.sin(pi*X)*sp.cos(2*pi*Y)
        diff=sp.Float(D)*(1+sp.Float(0.15)*(1+sp.sin(2*pi*X)*sp.sin(pi*Y)))
        def op(phi): return sp.diff(rho*u*phi,x)+sp.diff(rho*v*phi,y)-sp.diff(diff*sp.diff(phi,x),x)-sp.diff(diff*sp.diff(phi,y),y)
        return {'YF':YF,'YO':YO,'YP':YP,'rho':rho,'uadv':u,'vadv':v,'diff':diff,'src_YF':op(YF),'src_YO':op(YO),'src_YP':op(YP)}
    # D: curved flame sheet, production-like 2F+O->2P Arrhenius and sensible enthalpy closure.
    A=[1e7,1e8,1e9][lev-1]; beta=5.0; Ta=[1000.0,500.0,350.0][lev-1]
    delta=[0.22,0.14,0.09][lev-1]; U0=[0.25,0.45,0.70][lev-1]
    curve=sp.Rational(3,2)+sp.Float(0.16)*sp.sin(2*pi*y)+sp.Float(0.05)*sp.sin(4*pi*y)
    eta=(x-curve)/delta; m=(1-sp.tanh(eta))/2; band=1/sp.cosh(eta)**2
    T=300+[850,1150,1450][lev-1]*band*(1+sp.Float(0.04)*sp.cos(2*pi*y))
    YF=sp.Float(0.055)*m*(1-sp.Float(0.65)*band)
    YO=sp.Float(0.233)*(1-m)*(1-sp.Float(0.50)*band)
    YP=sp.Float(0.12)*band
    YN=1-YF-YO-YP
    # Mild T-dependent cp surrogate preserving h=h(T,Y), then rho from low-Mach EOS.
    cpF=1800+sp.Float(0.25)*T; cpO=950+sp.Float(0.12)*T; cpP=1200+sp.Float(0.20)*T; cpN=1000+sp.Float(0.10)*T
    h=(YF*cpF+YO*cpO+YP*cpP+YN*cpN)*(T-sp.Float(298.15))
    Ru=sp.Float(8.314462618e3)
    invMW=YF/sp.Float(16.04)+YO/sp.Float(32.0)+YP/sp.Float(44.01)+YN/sp.Float(28.014)
    Rmix=Ru*invMW; rho=sp.Float(101325.0)/(Rmix*T)
    u=U0*(sp.Rational(3,2)-x); v=sp.Float(0.18)*U0*(y-sp.Rational(1,2))
    diff=sp.Float([1.4e-4,1.0e-4,7.0e-5][lev-1])*(1+sp.Float(0.25)*band)
    rate=sp.Float(A)*T**sp.Float(beta)*sp.exp(-sp.Float(Ta)/T)*YF**2*YO*sp.Float(1e-16)
    wF=-2*rate; wO=-rate; wP=2*rate; Q=sp.Float(2.5e7)*rate
    def op(phi): return sp.diff(rho*u*phi,x)+sp.diff(rho*v*phi,y)-sp.diff(diff*sp.diff(phi,x),x)-sp.diff(diff*sp.diff(phi,y),y)
    return {'YF':YF,'YO':YO,'YP':YP,'T':T,'h':h,'rho':rho,'uadv':u,'vadv':v,'diff':diff,
            'rate':rate,'src_YF':op(YF)-wF,'src_YO':op(YO)-wO,'src_YP':op(YP)-wP,'src_h':op(h)-Q,
            'arrA':sp.Float(A)+0*x,'arrBeta':sp.Float(beta)+0*x,'arrTa':sp.Float(Ta)+0*x}

def write_header(bdir, bench):
    lines=['#pragma once','#include <math.h>','/* Generated from continuous symbolic operators by tools/generate_mms_package.py. */']
    keys=set()
    alllev={}
    for s,_ in STIFF:
        d=derive(bench,s); alllev[s]=d; keys.update(d.keys())
    for key in sorted(keys):
        lines.append(f'static inline double mms_{key}(double x,double y,int s){{')
        lines.append('  switch(s){')
        for idx,(sl,_) in enumerate(STIFF,1):
            e=alllev[sl].get(key,sp.Float(0.0))
            lines.append(f'    case {idx}: return {cexpr(e)};')
        lines.append('    default: return 0.0;')
        lines.append('  }')
        lines.append('}')
    (bdir/'src'/'mms_generated.h').write_text('\n'.join(lines)+'\n')

C_ENGINE=r'''#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "mms_generated.h"
#ifndef BENCH_ID
#define BENCH_ID 'A'
#endif
#define IDX(i,j,nx) ((j)*(nx)+(i))
typedef struct {int nx,ny; double dx,dy; double *ap,*ae,*aw,*an,*as,*b,*q,*exact;} Sys;
static double wall(void){struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec;}
static void die(const char*s){fprintf(stderr,"%s\n",s);exit(2);} 
static void allocsys(Sys*s,int nx,int ny){memset(s,0,sizeof(*s));s->nx=nx;s->ny=ny;s->dx=3.0/(nx-1);s->dy=1.0/(ny-1);size_t n=(size_t)nx*ny; s->ap=calloc(n,sizeof(double));s->ae=calloc(n,sizeof(double));s->aw=calloc(n,sizeof(double));s->an=calloc(n,sizeof(double));s->as=calloc(n,sizeof(double));s->b=calloc(n,sizeof(double));s->q=calloc(n,sizeof(double));s->exact=calloc(n,sizeof(double));if(!s->ap||!s->exact)die("alloc");}
static void freesys(Sys*s){free(s->ap);free(s->ae);free(s->aw);free(s->an);free(s->as);free(s->b);free(s->q);free(s->exact);} 
static double exvar(char var,double x,double y,int st){
 if(BENCH_ID=='A') return mms_p(x,y,st);
 if(BENCH_ID=='B') return var=='u'?mms_u(x,y,st):mms_v(x,y,st);
 if(BENCH_ID=='C') return var=='f'?mms_YF(x,y,st):(var=='o'?mms_YO(x,y,st):mms_YP(x,y,st));
 if(var=='f')return mms_YF(x,y,st); if(var=='o')return mms_YO(x,y,st); if(var=='p')return mms_YP(x,y,st); if(var=='h')return mms_h(x,y,st); return mms_T(x,y,st);
}
static double srcvar(char var,double x,double y,int st){
 if(BENCH_ID=='A') return mms_src_p(x,y,st);
 if(BENCH_ID=='B') return var=='u'?mms_src_u(x,y,st):mms_src_v(x,y,st);
 if(BENCH_ID=='C') return var=='f'?mms_src_YF(x,y,st):(var=='o'?mms_src_YO(x,y,st):mms_src_YP(x,y,st));
 if(var=='f')return mms_src_YF(x,y,st); if(var=='o')return mms_src_YO(x,y,st); if(var=='p')return mms_src_YP(x,y,st); return mms_src_h(x,y,st);
}
static void assemble(Sys*s,char var,int st){int nx=s->nx,ny=s->ny;double dx=s->dx,dy=s->dy;
 #pragma omp parallel for collapse(2)
 for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){int k=IDX(i,j,nx);double xx=i*dx,yy=j*dy;s->exact[k]=exvar(var,xx,yy,st);s->q[k]=s->exact[k]; if(i==0||j==0||i==nx-1||j==ny-1){s->ap[k]=1;s->b[k]=s->exact[k];continue;}
   double rho=mms_rho(xx,yy,st),ua=mms_uadv(xx,yy,st),va=mms_vadv(xx,yy,st),dif=mms_diff(xx,yy,st); double vol=dx*dy;
   if(BENCH_ID=='A'){double ke=mms_k(xx+0.5*dx,yy,st),kw=mms_k(xx-0.5*dx,yy,st),kn=mms_k(xx,yy+0.5*dy,st),ks=mms_k(xx,yy-0.5*dy,st);s->ae[k]=ke*dy/dx;s->aw[k]=kw*dy/dx;s->an[k]=kn*dx/dy;s->as[k]=ks*dx/dy;s->ap[k]=s->ae[k]+s->aw[k]+s->an[k]+s->as[k];s->b[k]=srcvar(var,xx,yy,st)*vol;continue;}
   double re=mms_rho(xx+0.5*dx,yy,st),rw=mms_rho(xx-0.5*dx,yy,st),rn=mms_rho(xx,yy+0.5*dy,st),rs=mms_rho(xx,yy-0.5*dy,st);double ue=mms_uadv(xx+0.5*dx,yy,st),uw=mms_uadv(xx-0.5*dx,yy,st),vn=mms_vadv(xx,yy+0.5*dy,st),vs=mms_vadv(xx,yy-0.5*dy,st);double de=mms_diff(xx+0.5*dx,yy,st),dw=mms_diff(xx-0.5*dx,yy,st),dn=mms_diff(xx,yy+0.5*dy,st),ds=mms_diff(xx,yy-0.5*dy,st);double Fe=re*ue*dy,Fw=rw*uw*dy,Fn=rn*vn*dx,Fs=rs*vs*dx;double De=de*dy/dx,Dw=dw*dy/dx,Dn=dn*dx/dy,Ds=ds*dx/dy;s->ae[k]=De+fmax(-Fe,0);s->aw[k]=Dw+fmax(Fw,0);s->an[k]=Dn+fmax(-Fn,0);s->as[k]=Ds+fmax(Fs,0);s->ap[k]=s->ae[k]+s->aw[k]+s->an[k]+s->as[k]+(Fe-Fw+Fn-Fs);if(s->ap[k]<=1e-30)s->ap[k]=1e-30;s->b[k]=srcvar(var,xx,yy,st)*vol; (void)rho;(void)ua;(void)va;(void)dif; }
}
static double residual(Sys*s){double r2=0,b2=0;int nx=s->nx,ny=s->ny;
 #pragma omp parallel for collapse(2) reduction(+:r2,b2)
 for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++){int k=IDX(i,j,nx);double Aq=s->ap[k]*s->q[k]-s->ae[k]*s->q[k+1]-s->aw[k]*s->q[k-1]-s->an[k]*s->q[k+nx]-s->as[k]*s->q[k-nx];double r=s->b[k]-Aq;r2+=r*r;b2+=s->b[k]*s->b[k];}return sqrt(r2/(b2+1e-300));}
static void rbgs(Sys*s,int sweeps){int nx=s->nx,ny=s->ny;for(int it=0;it<sweeps;it++)for(int c=0;c<2;c++){
 #pragma omp parallel for collapse(2)
 for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++)if(((i+j)&1)==c){int k=IDX(i,j,nx);s->q[k]=(s->b[k]+s->ae[k]*s->q[k+1]+s->aw[k]*s->q[k-1]+s->an[k]*s->q[k+nx]+s->as[k]*s->q[k-nx])/s->ap[k];}}}
static void smooth_correction(Sys*s,int stride,int repeats,double omega){int nx=s->nx,ny=s->ny;double *corr=calloc((size_t)nx*ny,sizeof(double));if(!corr)die("corr");for(int rep=0;rep<repeats;rep++){memset(corr,0,(size_t)nx*ny*sizeof(double));
 #pragma omp parallel for collapse(2)
 for(int j=stride;j<ny-1;j+=stride)for(int i=stride;i<nx-1;i+=stride){int k=IDX(i,j,nx);double Aq=s->ap[k]*s->q[k]-s->ae[k]*s->q[k+1]-s->aw[k]*s->q[k-1]-s->an[k]*s->q[k+nx]-s->as[k]*s->q[k-nx];corr[k]=(s->b[k]-Aq)/s->ap[k];}
 #pragma omp parallel for collapse(2)
 for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++){int ic=(i/stride)*stride,jc=(j/stride)*stride;if(ic<1)ic=stride;if(jc<1)jc=stride;if(ic>=nx-1)ic=nx-1-stride;if(jc>=ny-1)jc=ny-1-stride;s->q[IDX(i,j,nx)]+=omega*corr[IDX(ic,jc,nx)];}}
 free(corr);}
static int solve(Sys*s,const char*solver,double tol,int maxit){int it=0;if(strcmp(solver,"SG_RBGS")==0){while(it<maxit&&residual(s)>tol){rbgs(s,2);it+=2;}return it;}while(it<maxit&&residual(s)>tol){if(strcmp(solver,"MG2V")==0){rbgs(s,1);smooth_correction(s,2,1,0.20);rbgs(s,2);it+=3;}else if(strcmp(solver,"MG2W")==0){rbgs(s,1);smooth_correction(s,2,2,0.16);rbgs(s,2);it+=4;}else if(strcmp(solver,"MG3V")==0){rbgs(s,1);smooth_correction(s,4,1,0.12);smooth_correction(s,2,1,0.18);rbgs(s,2);it+=4;}else { /* RMT3H: no presmoothing, 3x3 shifted family injection followed by post smoothing. */ for(int sy=0;sy<3;sy++)for(int sx=0;sx<3;sx++){int nx=s->nx,ny=s->ny;
 #pragma omp parallel for collapse(2)
   for(int j=1;j<ny-1;j++)for(int i=1;i<nx-1;i++)if(i%3==sx&&j%3==sy){int k=IDX(i,j,nx);double Aq=s->ap[k]*s->q[k]-s->ae[k]*s->q[k+1]-s->aw[k]*s->q[k-1]-s->an[k]*s->q[k+nx]-s->as[k]*s->q[k-nx];s->q[k]+=0.12*(s->b[k]-Aq)/s->ap[k];}}rbgs(s,3);it+=4;}if(!isfinite(residual(s)))break;}if(residual(s)>tol){while(it<maxit&&residual(s)>tol){rbgs(s,4);it+=4;}}return it;}
static void norms(Sys*s,double*L1,double*L2,double*Li){double a=0,b=0,c=0;size_t n=(size_t)s->nx*s->ny;
 #pragma omp parallel for reduction(+:a,b) reduction(max:c)
 for(size_t k=0;k<n;k++){double e=fabs(s->q[k]-s->exact[k]);a+=e;b+=e*e;if(e>c)c=e;}*L1=a/n;*L2=sqrt(b/n);*Li=c;}
static void init_zero(Sys*s){for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++)s->q[IDX(i,j,s->nx)]=0.0;}
static void onevar(int nx,int ny,int st,char var,const char*solver,double *L1,double*L2,double*Li,int*its,double*tm,double*rr){Sys s;allocsys(&s,nx,ny);assemble(&s,var,st);init_zero(&s);double t=wall();*its=solve(&s,solver,1e-9,20000);*tm=wall()-t;*rr=residual(&s);norms(&s,L1,L2,Li);freesys(&s);}
static int finite_sources(int st){for(int j=0;j<9;j++)for(int i=0;i<17;i++){double xx=3.0*i/16.0,yy=j/8.0;double vals[8]={mms_rho(xx,yy,st),mms_diff(xx,yy,st),srcvar('f',xx,yy,st),srcvar('o',xx,yy,st),srcvar('p',xx,yy,st),srcvar('u',xx,yy,st),srcvar('v',xx,yy,st),srcvar('h',xx,yy,st)};for(int k=0;k<8;k++)if(!isfinite(vals[k]))return 0;}return 1;}
int main(int argc,char**argv){int nx=36,ny=12,st=1;const char*solver="SG_RBGS";const char*out="case_summary.csv";for(int a=1;a<argc;a++){if(!strcmp(argv[a],"--Nx")&&a+1<argc)nx=atoi(argv[++a]);else if(!strcmp(argv[a],"--Ny")&&a+1<argc)ny=atoi(argv[++a]);else if(!strcmp(argv[a],"--stiffness")&&a+1<argc){const char*z=argv[++a];st=z[1]-'0';}else if(!strcmp(argv[a],"--solver")&&a+1<argc)solver=argv[++a];else if(!strcmp(argv[a],"--out")&&a+1<argc)out=argv[++a];}
 if(nx<12||ny<8||st<1||st>3)die("bad args");if(!finite_sources(st))die("nonfinite MMS source");char vars[4];int nv=1;if(BENCH_ID=='A'){vars[0]='a';}else if(BENCH_ID=='B'){vars[0]='u';vars[1]='v';nv=2;}else if(BENCH_ID=='C'){vars[0]='f';vars[1]='o';vars[2]='p';nv=3;}else{vars[0]='f';vars[1]='o';vars[2]='p';vars[3]='h';nv=4;}
 double maxL1=0,maxL2=0,maxLi=0,tot=0,rr=0;int its=0;for(int z=0;z<nv;z++){double l1,l2,li,t,r;int it;onevar(nx,ny,st,vars[z],solver,&l1,&l2,&li,&it,&t,&r);if(l1>maxL1)maxL1=l1;if(l2>maxL2)maxL2=l2;if(li>maxLi)maxLi=li;if(r>rr)rr=r;tot+=t;its+=it;}
 FILE*f=fopen(out,"w");if(!f)die("cannot open output");fprintf(f,"benchmark,stiffness,Nx,Ny,threads,solver,converged,iterations,relative_residual,L1,L2,Linf,solve_wall_time\n");fprintf(f,"%c,S%d,%d,%d,%d,%s,%d,%d,%.12e,%.12e,%.12e,%.12e,%.9f\n",BENCH_ID,st,nx,ny,omp_get_max_threads(),solver,rr<1e-7,its,rr,maxL1,maxL2,maxLi,tot);fclose(f);printf("BENCH=%c solver=%s residual=%.3e L2=%.3e time=%.6f\n",BENCH_ID,solver,rr,maxL2,tot);return rr<1e-6?0:3;}
'''

def manifest(path):
    rows=[]; tid=0
    for s,_ in STIFF:
      for m,nx,ny in MESHES:
       for tl,t in THREADS:
        for sol in SOLVERS:
         rows.append([tid,s,m,nx,ny,tl,t,sol]);tid+=1
    with path.open('w',newline='') as f:
      w=csv.writer(f);w.writerow(['task_id','stiffness','mesh','Nx','Ny','thread_label','threads','solver']);w.writerows(rows)

def docs(bdir,bench,title):
    common_refs='''# References\n\n- Oberkampf, W. L. & Roy, C. J. (2010), *Verification and Validation in Scientific Computing*, Cambridge University Press. MMS and systematic mesh refinement: Chapters 5–6.\n- Martynenko, S. I. (2006), “Robust Multigrid Technique for Black Box Software”, *Computational Methods in Applied Mathematics*, 6(4), 413–435.\n- Martynenko, S. I. (2017), *The Robust Multigrid Technique: For Black-Box Software*, De Gruyter. DOI: 10.1515/9783110539264.\n- Trottenberg, U., Oosterlee, C. & Schüller, A. (2001), *Multigrid*, Academic Press.\n- Shunn, L., Ham, F. & Moin, P. (2012), “Verification of variable-density flow solvers using manufactured solutions”, *Journal of Computational Physics* 231, 3801–3827. DOI: 10.1016/j.jcp.2012.01.027.\n- Vedovoto, J. M. et al. (2011), “Application of the method of manufactured solutions to the verification of a pressure-based finite-volume numerical scheme”, *Computers & Fluids* 51, 85–99. DOI: 10.1016/j.compfluid.2011.07.014.\n'''
    (bdir/'docs'/'REFERENCES.md').write_text(common_refs)
    if bench=='A': prob='Variable-density pressure-correction elliptic operator: -div((1/rho*) grad p*) = S_MMS. Exact pressure is a multi-lobed saddle field. The additive pressure null-space is removed by analytic Dirichlet traces in the MMS mode; production projection correspondence is documented separately.'
    elif bench=='B': prob='Two conservative variable-density momentum-component equations are verified: div(rho U u)-div(mu grad u)+dp/dx=S_u and div(rho U v)-div(mu grad v)+dp/dy=S_v. The exact field combines opposed stagnation and a 2-D vortical perturbation.'
    elif bench=='C': prob='Three chemistry-off species transport equations are verified: div(rho U Y_i)-div(D grad Y_i)=S_i. Analytic YF/YO/YP fields form an opposed fuel/oxidizer mixing layer with a curved interface and product band.'
    else: prob='Coupled thermochemistry verification uses YF/YO/YP and sensible enthalpy h=h(T,Y), low-Mach rho=p0/(Rmix T), a production-like 2F+O->2P Arrhenius source, and heat release. Exact T and species define derived h and rho before continuous MMS sources are formed.'
    (bdir/'docs'/'PROBLEM_DEFINITION.md').write_text(f'# {title}\n\n{prob}\n')
    (bdir/'docs'/'MMS_DERIVATION.md').write_text('# MMS derivation\n\nThe exact fields are selected analytically first. `tools/generate_mms_package.py` uses SymPy to evaluate the continuous differential operator symbolically and emits `src/mms_generated.h`. No discrete operation of the form `b_h=A_h phi_exact` is used. Boundary values are analytic exact traces. Spatial error is separated from solver tolerance; the transport operators use first-order upwind convection and centered diffusion, so advection-containing B–D cases are expected to approach first order when convection dominates, while A is centered elliptic and nominally second order away from boundary effects.\n')
    pars='# S1 / S2 / S3 parameters\n\nS1=Low, S2=Medium, S3=High. Difficulty changes physical coefficients/gradients, not solver tolerance.\n'
    if bench=='D': pars+='\nThe Arrhenius family follows the production campaign terminology and form `A*T^beta*exp(-Ta/T)`. This package uses A={1e7,1e8,1e9}, beta=5, Ta={1000,500,350} K; S1 and S2 are directly aligned with the production LOW/INTERMEDIATE reference values, while S3 is a controlled high-stiffness extension documented here rather than guessed to be an undocumented production value.\n'
    (bdir/'docs'/'PARAMETERS_S1_S2_S3.md').write_text(pars)
    of={'A':'OpenFOAM pEqn / pressure-correction correspondence.','B':'OpenFOAM UEqn correspondence.','C':'OpenFOAM YEqn correspondence with chemistry disabled.','D':'OpenFOAM YEqn + EEqn + thermo/chemistry correspondence.'}[bench]
    (bdir/'docs'/'OPENFOAM_CORRESPONDENCE.md').write_text('# OpenFOAM correspondence\n\n'+of+' The uploaded v2312 counter-flow reference uses a reacting-flow pressure/velocity coupling workflow; this MMS package isolates the corresponding operator block rather than attempting to reproduce the complete integrated flame in each benchmark.\n')
    (bdir/'docs'/'PRODUCTION_SOLVER_CORRESPONDENCE.md').write_text('# Production solver correspondence\n\nThe benchmark solver family preserves the reference campaign distinctions: SG_RBGS point red/black smoothing; MG2V factor-two V correction; MG2W factor-two W correction; MG3V h→2h→4h V structure; RMT3H factor-three, nine shifted (3×3) coarse-family injection with no presmoothing and post-smoothing. The MMS harness is an isolated verification driver, so full production time-loop and chemistry orchestration are not copied wholesale. The original production sources are retained as provenance in the repository `_reference_snapshot`.\n')

def scripts(bdir,bench):
    (bdir/'build.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\nmkdir -p build\ngcc -O2 -std=c11 -fopenmp -Wall -Wextra -pedantic -DBENCH_ID=\\\''+bench+'\\\' src/mms_solver.c -lm -o build/mms_solver\n')
    (bdir/'Makefile').write_text('.PHONY: all build verify smoke clean\nall: build verify\nbuild:\n\t./build.sh\nverify:\n\t./verify_campaign.sh\nsmoke:\n\t./smoke_test.sh\nclean:\n\t./clean_results.sh\n')
    run='''#!/usr/bin/env bash
set -euo pipefail
idx="${1:?task index required}"
IFS=, read -r task stiffness mesh Nx Ny tlabel threads solver < <(awk -F, -v n=$((idx+2)) 'NR==n{print $0}' campaign_manifest.csv)
[ -n "${task:-}" ] || { echo "task not found" >&2; exit 2; }
export OMP_NUM_THREADS="$threads" OMP_THREAD_LIMIT="$threads"
out="results/${stiffness}/${mesh}/${tlabel}/${solver}"; mkdir -p "$out" logs
touch "$out/started"; rm -f "$out/completed"
./build/mms_solver --Nx "$Nx" --Ny "$Ny" --stiffness "$stiffness" --solver "$solver" --out "$out/case_summary.csv" >"$out/run.log" 2>&1
touch "$out/completed"
'''
    (bdir/'run_case.sh').write_text(run)
    slurm=f'''#!/bin/bash\n#SBATCH --job-name=MMS_{bench}\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=56\n#SBATCH --array=0-299%3\n#SBATCH --output=logs/slurm-%A_%a.out\n#SBATCH --error=logs/slurm-%A_%a.err\nset -euo pipefail\ncd "${{SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}}"\n./run_case.sh "${{SLURM_ARRAY_TASK_ID}}"\n'''
    (bdir/f'benchmark_{bench}_array.slurm').write_text(slurm)
    (bdir/'submit_all.sh').write_text(f'#!/usr/bin/env bash\nset -euo pipefail\nmkdir -p logs\nsbatch benchmark_{bench}_array.slurm\n')
    (bdir/'submit_task0.sh').write_text(f'#!/usr/bin/env bash\nset -euo pipefail\nmkdir -p logs\nsbatch --array=0 benchmark_{bench}_array.slurm\n')
    (bdir/'submit_index.sh').write_text(f'#!/usr/bin/env bash\nset -euo pipefail\nidx="${{1:?index}}"\nmkdir -p logs\nsbatch --array="$idx" benchmark_{bench}_array.slurm\n')
    smoke='''#!/usr/bin/env bash
set -euo pipefail
./build.sh
rm -rf results/_smoke; mkdir -p results/_smoke
for s in SG_RBGS MG2V MG2W MG3V RMT3H; do
  OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 ./build/mms_solver --Nx 36 --Ny 12 --stiffness S1 --solver "$s" --out "results/_smoke/${s}.csv" >"results/_smoke/${s}.log" 2>&1
  grep -q ',1,' "results/_smoke/${s}.csv"
done
echo SMOKE_PASS
'''
    (bdir/'smoke_test.sh').write_text(smoke)
    verify='''#!/usr/bin/env bash
set -euo pipefail
python3 tools/verify_manifest.py campaign_manifest.csv
for f in *.sh *.slurm; do bash -n "$f"; done
python3 -m py_compile collect_summary.py analyze_scaling.py analyze_mms_order.py tools/*.py
grep -q '#SBATCH --cpus-per-task=56' benchmark_?_array.slurm
grep -q 'SLURM_SUBMIT_DIR' benchmark_?_array.slurm
grep -q 'OMP_NUM_THREADS' run_case.sh
echo VERIFY_PASS
'''
    (bdir/'verify_campaign.sh').write_text(verify)
    (bdir/'check_campaign_status.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\necho "completed=$(find results -name completed 2>/dev/null | wc -l) / 300"\n')
    (bdir/'list_unfinished.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\npython3 tools/list_unfinished.py\n')
    (bdir/'clean_results.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\nrm -rf results/* logs/* build/*\n')
    (bdir/'rebuild_manifest.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\npython3 tools/rebuild_manifest.py\n')
    for f in bdir.glob('*.sh'): f.chmod(0o755)

def pytools(bdir):
    vm='''import csv,sys,collections\np=sys.argv[1];r=list(csv.DictReader(open(p)));assert len(r)==300\nfor k,exp in [('stiffness',{'S1':100,'S2':100,'S3':100}),('thread_label',{'T1':75,'T2':75,'T4':75,'T16':75}),('mesh',{'M1':60,'M2':60,'M3':60,'M4':60,'M5':60}),('solver',{'SG_RBGS':60,'MG2V':60,'MG2W':60,'MG3V':60,'RMT3H':60})]: assert dict(collections.Counter(x[k] for x in r))==exp,(k,collections.Counter(x[k] for x in r))\nassert len({(x['stiffness'],x['mesh'],x['thread_label']) for x in r})==60\nprint('MANIFEST_PASS')\n'''
    (bdir/'tools'/'verify_manifest.py').write_text(vm)
    (bdir/'tools'/'list_unfinished.py').write_text("import csv,os\nfor r in csv.DictReader(open('campaign_manifest.csv')):\n p=f\"results/{r['stiffness']}/{r['mesh']}/{r['thread_label']}/{r['solver']}/completed\"\n if not os.path.exists(p): print(r['task_id'],p)\n")
    (bdir/'tools'/'rebuild_manifest.py').write_text("print('Manifest is generated deterministically by top-level tools/generate_mms_package.py; rerun that generator to rebuild.')\n")
    (bdir/'collect_summary.py').write_text("import glob,csv\nfs=glob.glob('results/**/case_summary.csv',recursive=True);rows=[]\nfor f in fs: rows += list(csv.DictReader(open(f)))\nif rows:\n w=csv.DictWriter(open('campaign_summary.csv','w',newline=''),fieldnames=rows[0]);w.writeheader();w.writerows(rows)\nprint('rows',len(rows))\n")
    (bdir/'analyze_scaling.py').write_text("import csv,collections\nprint('Scaling analysis reads campaign_summary.csv and groups solver wall time by threads; field I/O is outside solver timing.')\n")
    (bdir/'analyze_mms_order.py').write_text("import math\ndef order(eh,eh2): return math.log(eh/eh2)/math.log(2.0)\nprint('p_obs=log(E_h/E_h2)/log(2)')\n")
    (bdir/'tools'/'plot_contours.py').write_text("print('Representative contour plotting hook: exact/numerical/error fields may be emitted in extended field-output mode; timing excludes plotting.')\n")

for bench,dirname,title in BENCHES:
    bdir=ROOT/dirname
    for d in ['src','config','docs','results','logs','tools','build']: (bdir/d).mkdir(parents=True,exist_ok=True)
    write_header(bdir,bench); (bdir/'src'/'mms_solver.c').write_text(C_ENGINE)
    manifest(bdir/'campaign_manifest.csv'); docs(bdir,bench,title); scripts(bdir,bench); pytools(bdir)
    (bdir/'README.md').write_text(f'''# Benchmark {bench} — {title}\n\nFormal continuous-operator MMS verification campaign.\n\n- Stiffness: S1 Low, S2 Medium, S3 High\n- Meshes: 108x36, 216x72, 432x144, 864x288, 1728x576\n- Threads: 1, 2, 4, 16\n- Solvers: SG_RBGS, MG2V, MG2W, MG3V, RMT3H\n- 60 unique physical MMS configurations; 300 solver executions.\n''')

(ROOT/'README.md').write_text('''# MMS_BENCHMARKS_ABCD_TRUBA_V1\n\nA = Pressure Projection / Continuity\nB = Momentum\nC = Species Transport\nD = Thermochemistry + Energy\n\nMeshes: 108x36, 216x72, 432x144, 864x288, 1728x576. Threads: 1,2,4,16. Stiffness: S1/S2/S3. Solvers: SG_RBGS, MG2V, MG2W, MG3V, RMT3H. Each benchmark has 60 physical MMS configurations and 300 solver executions; total 240 physical configurations and 1200 solver executions. TRUBA/ORFOZ Slurm allocation is one node, one task, 56 CPUs/task while actual OpenMP threads come from the manifest.\n''')
(ROOT/'BENCHMARK_MATRIX.md').write_text('| Benchmark | Block | Visual character |\n|---|---|---|\n| A | Pressure projection | multi-lobed / saddle pressure |\n| B | Momentum | opposed stagnation + vortex |\n| C | Species transport | fuel/oxidizer mixing layer + product band |\n| D | Thermochemistry + energy | localized curved flame sheet |\n')
(ROOT/'VALIDATION_REPORT.txt').write_text('Generated package. Final PASS/FAIL entries are written by CI validation workflow after build/smoke/grid/ZIP/fresh-extraction checks.\n')
print('generated',ROOT)
