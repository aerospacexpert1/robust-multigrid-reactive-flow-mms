#!/usr/bin/env python3
"""Compute matched 1/2/4/16-thread scaling for RMT FINAL."""

from __future__ import annotations
import argparse,csv,math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
DEFAULT_INPUT=ROOT/"RMT_FINAL_campaign_summary.csv"
DEFAULT_OUTPUT=ROOT/"RMT_FINAL_scaling_summary.csv"

def fval(row,key):
    try:return float(row.get(key,""))
    except (TypeError,ValueError):return math.nan

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,default=DEFAULT_INPUT)
    ap.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    args=ap.parse_args()
    if not args.input.is_file():
        print(f"ERROR: missing {args.input}",file=sys.stderr)
        print("Run python3 collect_summary.py first.",file=sys.stderr)
        return 2

    with args.input.open(newline="") as f:
        rows=list(csv.DictReader(f))
    complete=[r for r in rows if r.get("status")=="complete"]
    groups={}
    for r in complete:
        groups.setdefault((r["study_family"],r["physical_case"],r["mesh_name"]),[]).append(r)

    out=[]
    for _,grows in sorted(groups.items()):
        by={int(r["threads"]):r for r in grows}
        if 1 not in by: continue
        t1=fval(by[1],"compute_time_s")
        p1=fval(by[1],"pressure_time_s")
        for th in sorted(by):
            r=by[th]; tc=fval(r,"compute_time_s"); tp=fval(r,"pressure_time_s")
            sp=t1/tc if t1>0 and tc>0 else math.nan
            eff=sp/th if math.isfinite(sp) else math.nan
            psp=p1/tp if p1>0 and tp>0 else math.nan
            out.append({
                "study_family":r["study_family"],
                "physical_case":r["physical_case"],
                "mesh_name":r["mesh_name"],
                "Nx":r["Nx"],"Ny":r["Ny"],"threads":str(th),
                "compute_time_s":r.get("compute_time_s",""),
                "pressure_time_s":r.get("pressure_time_s",""),
                "pressure_fraction_pct":r.get("pressure_fraction_pct",""),
                "speedup_vs_T1":f"{sp:.12g}" if math.isfinite(sp) else "",
                "parallel_efficiency":f"{eff:.12g}" if math.isfinite(eff) else "",
                "pressure_speedup_vs_T1":f"{psp:.12g}" if math.isfinite(psp) else "",
                "avg_rmt_cycles_per_pressure_solve":r.get("avg_rmt_cycles_per_pressure_solve",""),
                "max_rmt_cycles_per_pressure_solve":r.get("max_rmt_cycles_per_pressure_solve",""),
                "mean_pressure_convergence_factor":r.get("mean_pressure_convergence_factor",""),
                "final_pressure_rel_residual":r.get("final_pressure_rel_residual",""),
                "mass_residual_max":r.get("mass_residual_max",""),
                "Tmax_K":r.get("Tmax_K",""),
                "diag_nonmonotone_cycles":r.get("diag_nonmonotone_cycles",""),
                "diag_max_cycle_ratio":r.get("diag_max_cycle_ratio",""),
                "diag_direct_solves":r.get("diag_direct_solves",""),
            })
    fields=list(out[0].keys()) if out else [
        "study_family","physical_case","mesh_name","Nx","Ny","threads",
        "compute_time_s","pressure_time_s","pressure_fraction_pct","speedup_vs_T1",
        "parallel_efficiency","pressure_speedup_vs_T1"
    ]
    with args.output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(out)
    print(f"Wrote {args.output} with {len(out)} rows")
    print(f"Complete campaign rows used: {len(complete)} / {len(rows)}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
