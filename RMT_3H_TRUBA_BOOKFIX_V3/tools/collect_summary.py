#!/usr/bin/env python3
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for p in ROOT.glob('results/*/*/*/T*/summary.csv'):
    with p.open(newline='') as f:
        r=list(csv.DictReader(f))
        if r: rows.append(r[0])
rows.sort(key=lambda r:int(r['task_id']))
out=ROOT/'RMT_BOOKFIX_V3_campaign_summary.csv'
if rows:
    fields=list(rows[0].keys())
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
else:
    out.write_text('')
print(f'Wrote {out} with {len(rows)} rows')
