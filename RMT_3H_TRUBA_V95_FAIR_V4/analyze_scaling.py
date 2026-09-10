#!/usr/bin/env python3
"""Compute per-mesh thread scaling from RMT_V95_campaign_summary.csv."""

from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "RMT_V95_campaign_summary.csv"
DEFAULT_OUTPUT = ROOT / "RMT_V95_scaling_summary.csv"


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return math.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()

    if not args.input.is_file():
        print(f"ERROR: missing {args.input}", file=sys.stderr)
        print("Run python3 collect_summary.py first.", file=sys.stderr)
        return 2

    with args.input.open(newline="") as f:
        rows = list(csv.DictReader(f))

    complete = [r for r in rows if r.get("status") == "complete"]
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for r in complete:
        key = (r["study_family"], r["physical_case"], r["mesh_name"])
        groups.setdefault(key, []).append(r)

    out: list[dict[str, str]] = []
    for key, grows in sorted(groups.items()):
        by_threads = {int(r["threads"]): r for r in grows}
        if 1 not in by_threads:
            continue
        base = by_threads[1]
        t1 = as_float(base, "compute_time_s")
        p1 = as_float(base, "pressure_time_s")
        for th in sorted(by_threads):
            r = by_threads[th]
            tc = as_float(r, "compute_time_s")
            tp = as_float(r, "pressure_time_s")
            speedup = t1 / tc if t1 > 0.0 and tc > 0.0 else math.nan
            efficiency = speedup / th if math.isfinite(speedup) else math.nan
            pspeedup = p1 / tp if p1 > 0.0 and tp > 0.0 else math.nan
            out.append(
                {
                    "study_family": r["study_family"],
                    "physical_case": r["physical_case"],
                    "mesh_name": r["mesh_name"],
                    "Nx": r["Nx"],
                    "Ny": r["Ny"],
                    "threads": str(th),
                    "compute_time_s": r.get("compute_time_s", ""),
                    "pressure_time_s": r.get("pressure_time_s", ""),
                    "pressure_fraction_pct": r.get("pressure_fraction_pct", ""),
                    "speedup_vs_T1": f"{speedup:.12g}" if math.isfinite(speedup) else "",
                    "parallel_efficiency": f"{efficiency:.12g}" if math.isfinite(efficiency) else "",
                    "pressure_speedup_vs_T1": f"{pspeedup:.12g}" if math.isfinite(pspeedup) else "",
                    "final_pressure_rel_residual": r.get("final_pressure_rel_residual", ""),
                    "mass_residual_max": r.get("mass_residual_max", ""),
                    "Tmax_K": r.get("Tmax_K", ""),
                }
            )

    fields = [
        "study_family",
        "physical_case",
        "mesh_name",
        "Nx",
        "Ny",
        "threads",
        "compute_time_s",
        "pressure_time_s",
        "pressure_fraction_pct",
        "speedup_vs_T1",
        "parallel_efficiency",
        "pressure_speedup_vs_T1",
        "final_pressure_rel_residual",
        "mass_residual_max",
        "Tmax_K",
    ]

    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(out)

    print(f"Wrote {args.output} with {len(out)} rows")
    print(f"Complete campaign rows used: {len(complete)} / {len(rows)}")
    print(f"Scaling groups with a T1 reference: {len(set((r['study_family'], r['physical_case'], r['mesh_name']) for r in out))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
