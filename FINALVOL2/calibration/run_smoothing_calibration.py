#!/usr/bin/env python3
"""Independent RMT smoothing calibration for FINALVOL2.

This script never edits production source or run_case.sh.  It evaluates one
single uniform postsmoothing count across deterministic manufactured pressure
problems on the four campaign mesh sizes.  Production reacting-flow cases are
not used for parameter selection.

Recommendation rule:
  1. candidate must pass every manufactured mode/mesh;
  2. max cycles <= --max-cycles and max observed rho <= --max-rho;
  3. among robust candidates choose the smallest total point-update work;
  4. wall time is reported, not used to select the numerical parameter.
"""
from __future__ import annotations
import argparse,csv,math,os,re,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CAL=Path(__file__).resolve().parent
BIN=ROOT/"build"/"rmt_smoothing_calibration"
PAT=re.compile(
    r"CAL_RESULT sweeps=(?P<sweeps>\d+) mode=(?P<mode>\S+) "
    r"Nx=(?P<Nx>\d+) Ny=(?P<Ny>\d+) cycles=(?P<cycles>\d+) "
    r"relative=(?P<relative>\S+) rho=(?P<rho>\S+) "
    r"point_updates=(?P<updates>\d+) lu_factorizations=(?P<lu>\d+) "
    r"elapsed_s=(?P<elapsed>\S+) pass=(?P<pass>[01])"
)

def compile_binary():
    (ROOT/"build").mkdir(exist_ok=True)
    cmd=["gcc","-O3","-march=native","-std=c11","-Wall","-Wextra","-Wpedantic",
         "-fopenmp",str(CAL/"rmt_smoothing_calibration.c"),"-lm","-o",str(BIN)]
    subprocess.run(cmd,check=True,cwd=ROOT)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",default="8,10,12,14,16")
    ap.add_argument("--threads",type=int,default=16)
    ap.add_argument("--max-cycles",type=int,default=20)
    ap.add_argument("--max-rho",type=float,default=0.75)
    args=ap.parse_args()
    cand=sorted({int(x) for x in args.candidates.split(",") if x.strip()})
    if not cand or min(cand)<1: raise SystemExit("invalid candidates")
    compile_binary()

    env=os.environ.copy()
    env.update({
        "OMP_NUM_THREADS":str(args.threads),
        "OMP_DYNAMIC":"false",
        "OMP_PROC_BIND":"close",
        "OMP_PLACES":"cores",
    })
    rows=[]
    for sw in cand:
        print(f"===== CALIBRATION sweeps={sw} threads={args.threads} =====",flush=True)
        p=subprocess.run([str(BIN),str(sw)],cwd=ROOT,env=env,text=True,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        print(p.stdout,end="")
        parsed=0
        for line in p.stdout.splitlines():
            m=PAT.search(line)
            if not m: continue
            d=m.groupdict()
            rows.append({
                "sweeps":int(d["sweeps"]),"mode":d["mode"],
                "Nx":int(d["Nx"]),"Ny":int(d["Ny"]),
                "cycles":int(d["cycles"]),"relative":float(d["relative"]),
                "rho":float(d["rho"]),"point_updates":int(d["updates"]),
                "lu_factorizations":int(d["lu"]),"elapsed_s":float(d["elapsed"]),
                "pass":int(d["pass"]),"process_rc":p.returncode,
                "threads":args.threads,
            })
            parsed+=1
        if parsed!=12:
            raise SystemExit(f"expected 12 CAL_RESULT rows for sweeps={sw}, got {parsed}")

    results=CAL/"calibration_results.csv"
    with results.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

    summary=[]
    for sw in cand:
        rr=[r for r in rows if r["sweeps"]==sw]
        allpass=all(r["pass"]==1 for r in rr)
        maxcy=max(r["cycles"] for r in rr)
        maxrho=max(r["rho"] for r in rr)
        work=sum(r["point_updates"] for r in rr)
        elapsed=sum(r["elapsed_s"] for r in rr)
        robust=allpass and maxcy<=args.max_cycles and maxrho<=args.max_rho
        summary.append({
            "sweeps":sw,"all_12_pass":int(allpass),"max_cycles":maxcy,
            "max_rho":maxrho,"total_point_updates":work,
            "total_elapsed_s":elapsed,"robust_candidate":int(robust),
        })

    sumpath=CAL/"calibration_summary.csv"
    with sumpath.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]))
        w.writeheader();w.writerows(summary)

    valid=[r for r in summary if r["robust_candidate"]]
    rec=min(valid,key=lambda r:(r["total_point_updates"],r["max_cycles"],r["sweeps"])) if valid else None
    txt=[
        "FINALVOL2 RMT SMOOTHING CALIBRATION",
        f"threads_for_timing_only={args.threads}",
        f"acceptance_max_cycles={args.max_cycles}",
        f"acceptance_max_rho={args.max_rho}",
        "selection_metric=minimum_total_point_updates_among_robust_candidates",
        "wall_time_used_for_selection=no",
        "production_reacting_cases_used_for_selection=no",
        "",
    ]
    for r in summary:
        txt.append(
            f"sweeps={r['sweeps']} robust={r['robust_candidate']} "
            f"max_cycles={r['max_cycles']} max_rho={r['max_rho']:.6f} "
            f"point_updates={r['total_point_updates']} elapsed_s={r['total_elapsed_s']:.6f}"
        )
    txt+=[""]
    if rec:
        txt.append(f"RECOMMENDED_SMOOTH_SWEEPS={rec['sweeps']}")
        txt.append("NOTE=Recommendation is advisory only; production source was NOT modified.")
    else:
        txt.append("RECOMMENDED_SMOOTH_SWEEPS=NONE")
        txt.append("NOTE=No candidate satisfied the predeclared robustness criteria.")
    (CAL/"CALIBRATION_RECOMMENDATION.txt").write_text("\n".join(txt)+"\n")
    print("\n".join(txt))
    print(f"Wrote {results}")
    print(f"Wrote {sumpath}")

if __name__=="__main__":
    main()
