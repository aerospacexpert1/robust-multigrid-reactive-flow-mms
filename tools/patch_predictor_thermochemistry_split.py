#!/usr/bin/env python3
from pathlib import Path
p=Path('tools/generate_mms_package.py')
s=p.read_text()
old="return {'u':u,'v':v,'p':p,'rho':rho,'uadv':u,'vadv':v,'diff':sp.Float(mu),'src_u':op(u)+sp.diff(p,x),'src_v':op(v)+sp.diff(p,y)}"
new="return {'u':u,'v':v,'p':p,'rho':rho,'uadv':u,'vadv':v,'diff':sp.Float(mu),'dpdx':sp.diff(p,x),'dpdy':sp.diff(p,y),'src_u':op(u)+sp.diff(p,x),'src_v':op(v)+sp.diff(p,y)}"
assert old in s
s=s.replace(old,new)
old2="if(BENCH_ID=='B') return var=='u'?mms_src_u(x,y,st):mms_src_v(x,y,st);"
new2="if(BENCH_ID=='B') return var=='u'?(mms_src_u(x,y,st)-mms_dpdx(x,y,st)):(mms_src_v(x,y,st)-mms_dpdy(x,y,st));"
assert old2 in s
s=s.replace(old2,new2)
old3="if(var=='f')return mms_src_YF(x,y,st); if(var=='o')return mms_src_YO(x,y,st); if(var=='p')return mms_src_YP(x,y,st); return mms_src_h(x,y,st);"
new3="if(var=='f')return mms_src_YF(x,y,st)-2.0*mms_rate(x,y,st); if(var=='o')return mms_src_YO(x,y,st)-mms_rate(x,y,st); if(var=='p')return mms_src_YP(x,y,st)+2.0*mms_rate(x,y,st); return mms_src_h(x,y,st)+2.5e7*mms_rate(x,y,st);"
assert old3 in s
s=s.replace(old3,new3)
p.write_text(s)

q=Path('tools/postprocess_package.py')
t=q.read_text()
t=t.replace("'arrA','arrBeta','arrTa']","'arrA','arrBeta','arrTa','dpdx','dpdy']")
q.write_text(t)
print('predictor and thermochemistry split patched')
