#!/usr/bin/env python3
from pathlib import Path
import csv, re, sys

if len(sys.argv) != 3:
    raise SystemExit('usage: extract_one_summary.py OUTDIR MANIFEST_ROW')
out=Path(sys.argv[1])
row=sys.argv[2].split(',')
keys=['task_id','case_id','family','level','mesh','Nx','Ny','threads','vIn','arrA','arrBeta','arrTa','endTime']
meta=dict(zip(keys,row))
text=(out/'solver.log').read_text(errors='replace')
sm=re.findall(r'RMT_RUN_SUMMARY\s+([^\n]+)',text)
dg=re.findall(r'RMT_BOOK_DIAGNOSTICS\s+([^\n]+)',text)
if not sm or not dg:
    raise SystemExit('missing RMT summary/diagnostics')

def kv(line):
    d={}
    for k,v in re.findall(r'(\w+)=([^\s]+)',line): d[k]=v
    return d
meta.update(kv(sm[-1])); meta.update(kv(dg[-1]))
fields=keys+['final_time','total_wall_s','pressure_wall_s','pressure_calls','Tmax','YPmax','YOmax','YFmax',
             'levelsX','levelsY','full_accepts','backtracks','rejects','last_omega','point_updates']
with (out/'summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerow({k:meta.get(k,'') for k in fields})
print(out/'summary.csv')
