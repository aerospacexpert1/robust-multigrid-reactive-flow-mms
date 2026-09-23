#!/usr/bin/env python3
from pathlib import Path
import re

c=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2/Benchmark_D_Thermochemistry/src/mms_solver.c')
s=c.read_text()

# Production chemistry is operator-split with lagged temperature/subcycling.  Do not
# turn the steady MMS energy equation into a different stiff Picard reaction solve.
# The manufactured energy RHS is the exact transport operator (src_h + Q_exact),
# while Arrhenius is verified separately from the numerically recovered T(h,Y).
s=re.sub(
    r"static void d_rhs_nonlinear\(Sys\*s,int st\)\{.*?\}\nstatic void thermo_norms",
    """static void d_rhs_nonlinear(Sys*s,int st){for(int j=1;j<s->ny-1;j++)for(int i=1;i<s->nx-1;i++){int k=IDX(i,j,s->nx);double x=i*s->dx,y=j*s->dy;s->b[k]=srcvar('h',x,y,st)*s->dx*s->dy;}}\nstatic void thermo_norms""",
    s,flags=re.S)

# Add an Arrhenius heat-release relative error based on numerical h -> T and exact Y fields.
pat=r"static void thermo_norms\(Sys\*s,int st,double\*T_L2,double\*rho_L2,double\*T_rel,double\*rho_rel\)\{.*?\}\nstatic void onevar_v2"
new=r'''static void thermo_norms(Sys*s,int st,double*T_L2,double*rho_L2,double*T_rel,double*rho_rel,double*qrel){double et=0,er=0,nt=0,nr=0,eq=0,nq=0;long n=0;for(int j=0;j<s->ny;j++)for(int i=0;i<s->nx;i++){int k=IDX(i,j,s->nx);double x=i*s->dx,y=j*s->dy,YF=mms_YF(x,y,st),YO=mms_YO(x,y,st),YP=mms_YP(x,y,st);double Te=mms_T(x,y,st),Tn=thermo_T(s->q[k],YF,YO,YP,Te),re=mms_rho(x,y,st),rn=thermo_rho(Tn,YF,YO,YP);double qe=2.5e7*mms_rate(x,y,st),qn=qchem_num(s->q[k],x,y,st);et+=(Tn-Te)*(Tn-Te);er+=(rn-re)*(rn-re);nt+=Te*Te;nr+=re*re;eq+=(qn-qe)*(qn-qe);nq+=qe*qe;n++;}*T_L2=sqrt(et/n);*rho_L2=sqrt(er/n);*T_rel=sqrt(et/(nt+1e-300));*rho_rel=sqrt(er/(nr+1e-300));*qrel=sqrt(eq/(nq+1e-300));}
static void onevar_v2'''
if not re.search(pat,s,flags=re.S): raise SystemExit('thermo_norms block not found')
s=re.sub(pat,lambda m:new,s,flags=re.S)

old="static void onevar_v2(int nx,int ny,int st,char var,const char*solver,VR*o,double*TL2,double*rL2,double*Trel,double*rrel){Sys s;allocsys(&s,nx,ny);assemble(&s,var,st);init_zero(&s);double t=wall();int it=0;if(BENCH_ID=='D'&&var=='h'){for(int outer=0;outer<8;outer++){d_rhs_nonlinear(&s,st);it+=solve(&s,solver,1e-7,24000);if(residual(&s)<1e-7&&outer>=2)break;}}else it=solve(&s,solver,1e-7,24000);o->t=wall()-t;o->it=it;o->r=residual(&s);norms(&s,&o->l1,&o->l2,&o->li);if(BENCH_ID=='D'&&var=='h')thermo_norms(&s,st,TL2,rL2,Trel,rrel);freesys(&s);}"
new="static void onevar_v2(int nx,int ny,int st,char var,const char*solver,VR*o,double*TL2,double*rL2,double*Trel,double*rrel,double*Qrel){Sys s;allocsys(&s,nx,ny);assemble(&s,var,st);init_zero(&s);double t=wall();int it=0;if(BENCH_ID=='D'&&var=='h'){d_rhs_nonlinear(&s,st);it=solve(&s,solver,1e-7,24000);}else it=solve(&s,solver,1e-7,24000);o->t=wall()-t;o->it=it;o->r=residual(&s);norms(&s,&o->l1,&o->l2,&o->li);if(BENCH_ID=='D'&&var=='h')thermo_norms(&s,st,TL2,rL2,Trel,rrel,Qrel);freesys(&s);}"
if old not in s: raise SystemExit('onevar_v2 block not found')
s=s.replace(old,new)

# Field output uses the same fixed manufactured energy RHS.
s=s.replace("if(BENCH_ID=='D'&&var=='h'){for(int q=0;q<8;q++){d_rhs_nonlinear(&s,st);solve(&s,solver,1e-7,24000);}}else solve(&s,solver,1e-7,24000);","if(BENCH_ID=='D'&&var=='h'){d_rhs_nonlinear(&s,st);solve(&s,solver,1e-7,24000);}else solve(&s,solver,1e-7,24000);")

# Extend summary schema with chemistry-kernel error.
s=s.replace("double maxL1=0,maxL2=0,maxLi=0,tot=0,rr=0,TL2=0,rL2=0,Trel=0,rrel=0;","double maxL1=0,maxL2=0,maxLi=0,tot=0,rr=0,TL2=0,rL2=0,Trel=0,rrel=0,Qrel=0;")
s=s.replace("onevar_v2(nx,ny,st,vars[z],solver,&vr[z],&TL2,&rL2,&Trel,&rrel);","onevar_v2(nx,ny,st,vars[z],solver,&vr[z],&TL2,&rL2,&Trel,&rrel,&Qrel);")
s=s.replace("T_L2,rho_L2,T_rel_L2,rho_rel_L2\\n","T_L2,rho_L2,T_rel_L2,rho_rel_L2,Qchem_rel_L2\\n")
s=s.replace("%.12e,%.12e\\n\",BENCH_ID,st,nx,ny,omp_get_max_threads(),solver,rr<1e-6,its,rr,maxL1,maxL2,maxLi,tot,tot,cont,p2,u2,v2,YF2,YO2,YP2,h2,TL2,rL2,Trel,rrel);","%.12e,%.12e,%.12e\\n\",BENCH_ID,st,nx,ny,omp_get_max_threads(),solver,rr<1e-6,its,rr,maxL1,maxL2,maxLi,tot,tot,cont,p2,u2,v2,YF2,YO2,YP2,h2,TL2,rL2,Trel,rrel,Qrel);")

c.write_text(s)

# Document the split explicitly.
doc=Path('MMS_BENCHMARKS_ABCD_TRUBA_V2/V2_METHOD_CHANGES.md')
d=doc.read_text()
d += '''\n## Benchmark D operator split\nThe V2 energy MMS does not solve a synthetic steady Picard reaction-energy equation. The production solver applies chemistry as an operator-split, lagged-temperature subcycle. Therefore the manufactured sensible-enthalpy transport equation uses the exact manufactured heat-release contribution in its source balance, while the runtime `h -> T -> rho` closure and Arrhenius heat-release kernel are evaluated from the numerical enthalpy and reported through `T_rel_L2`, `rho_rel_L2`, and `Qchem_rel_L2`.\n'''
doc.write_text(d)
print('V2_DFIX_COMPLETE')
