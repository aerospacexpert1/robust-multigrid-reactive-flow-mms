#!/usr/bin/env python3
"""Collect all 96 RMT FINAL case summaries and RMT diagnostics."""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
DEFAULT_MANIFEST=ROOT/"campaign_manifest.csv"
DEFAULT_OUTPUT=ROOT/"RMT_FINAL_campaign_summary.csv"

def read_key_value_file(path:Path)->dict[str,str]:
    data={}
    if not path.is_file():
        return data
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            k,v=line.split("=",1); data[k.strip()]=v.strip()
    return data

def read_rmt_diag(path:Path)->dict[str,str]:
    if not path.is_file():
        return {}
    last=None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("RMT_FINAL_DIAGNOSTICS "):
            last=line
    if not last:
        return {}
    out={}
    for tok in last.split()[1:]:
        if "=" in tok:
            k,v=tok.split("=",1)
            out[f"diag_{k}"]=v
    return out

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,default=DEFAULT_MANIFEST)
    ap.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument("--require-complete",action="store_true")
    args=ap.parse_args()

    if not args.manifest.is_file():
        print(f"ERROR: missing manifest: {args.manifest}",file=sys.stderr)
        print("Run ./rebuild_manifest.sh first.",file=sys.stderr)
        return 2

    with args.manifest.open(newline="") as f:
        manifest=list(csv.DictReader(f))
    if len(manifest)!=96:
        print(f"ERROR: expected 96 manifest rows, found {len(manifest)}",file=sys.stderr)
        return 3

    rows=[]
    extra_fields=[]
    counts={"complete":0,"failed":0,"partial":0,"missing":0}
    for m in manifest:
        run_dir=ROOT/"results"/m["group_dir"]/m["physical_case"]/m["mesh_name"]/m["core_label"]
        summary=run_dir/"output"/"summary.csv"
        complete=(run_dir/"RUN_COMPLETE").is_file()
        failed=(run_dir/"RUN_FAILED").is_file()
        has_summary=summary.is_file()
        status="complete" if complete and has_summary else "failed" if failed else "partial" if has_summary else "missing"
        counts[status]+=1

        row=dict(m)
        row["status"]=status
        row["run_dir"]=str(run_dir.relative_to(ROOT))
        timing=read_key_value_file(run_dir/"launcher_timing.txt")
        row["launcher_wall_seconds"]=timing.get("launcher_wall_seconds_including_launcher_and_output","")
        row["exit_code"]=timing.get("exit_code","")

        if has_summary:
            with summary.open(newline="") as f:
                sr=next(csv.DictReader(f),None)
            if sr:
                for k,v in sr.items():
                    outk=k if k not in row else f"summary_{k}"
                    row[outk]=v
                    if outk not in extra_fields: extra_fields.append(outk)

        for k,v in read_rmt_diag(run_dir/"run.log").items():
            row[k]=v
            if k not in extra_fields: extra_fields.append(k)
        rows.append(row)

    base=list(manifest[0].keys())+["status","run_dir","launcher_wall_seconds","exit_code"]
    fields=base+[k for k in extra_fields if k not in base]
    with args.output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    incomplete=counts["failed"]+counts["partial"]+counts["missing"]
    print(f"Wrote {args.output} with {len(rows)} rows")
    print("Status: "+" ".join(f"{k}={counts[k]}" for k in ("complete","failed","partial","missing")))
    print(f"Missing summaries: {sum(r['status']=='missing' for r in rows)}")
    return 4 if args.require_complete and incomplete else 0

if __name__=="__main__":
    raise SystemExit(main())
