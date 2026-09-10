#!/usr/bin/env python3
"""Collect the 96 V95 FAIR RMT case summaries into one campaign CSV."""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "campaign_manifest.csv"
DEFAULT_OUTPUT = ROOT / "RMT_V95_campaign_summary.csv"


def read_key_value_file(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    if not path.is_file():
        return data
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            data[key.strip()] = value.strip()
    return data


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument(
        "--require-complete",
        action="store_true",
        help="return non-zero unless every manifest row has RUN_COMPLETE and summary.csv",
    )
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"ERROR: missing manifest: {args.manifest}", file=sys.stderr)
        print("Run ./rebuild_manifest.sh first.", file=sys.stderr)
        return 2

    with args.manifest.open(newline="") as f:
        manifest_rows = list(csv.DictReader(f))

    if len(manifest_rows) != 96:
        print(f"ERROR: expected 96 manifest rows, found {len(manifest_rows)}", file=sys.stderr)
        return 3

    merged: list[dict[str, str]] = []
    summary_fields: list[str] = []
    counts = {"complete": 0, "failed": 0, "partial": 0, "missing": 0}

    for m in manifest_rows:
        run_dir = (
            ROOT
            / "results"
            / m["group_dir"]
            / m["physical_case"]
            / m["mesh_name"]
            / m["core_label"]
        )
        summary_path = run_dir / "output" / "summary.csv"
        timing = read_key_value_file(run_dir / "launcher_timing.txt")

        has_complete = (run_dir / "RUN_COMPLETE").is_file()
        has_failed = (run_dir / "RUN_FAILED").is_file()
        has_summary = summary_path.is_file()

        if has_complete and has_summary:
            status = "complete"
        elif has_failed:
            status = "failed"
        elif has_summary:
            status = "partial"
        else:
            status = "missing"
        counts[status] += 1

        row: dict[str, str] = dict(m)
        row["status"] = status
        row["run_dir"] = str(run_dir.relative_to(ROOT))
        row["launcher_wall_seconds"] = timing.get(
            "launcher_wall_seconds_including_launcher_and_output", ""
        )
        row["exit_code"] = timing.get("exit_code", "")

        if has_summary:
            with summary_path.open(newline="") as sf:
                reader = csv.DictReader(sf)
                srow = next(reader, None)
            if srow is None:
                status = "partial"
                row["status"] = status
            else:
                for key, value in srow.items():
                    if key not in summary_fields:
                        summary_fields.append(key)
                    # Keep manifest metadata authoritative when names overlap.
                    out_key = key if key not in row else f"summary_{key}"
                    row[out_key] = value

        merged.append(row)

    base_fields = list(manifest_rows[0].keys())
    meta_fields = ["status", "run_dir", "launcher_wall_seconds", "exit_code"]
    output_fields = base_fields + meta_fields
    for key in summary_fields:
        out_key = key if key not in output_fields else f"summary_{key}"
        if out_key not in output_fields:
            output_fields.append(out_key)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=output_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    incomplete = counts["failed"] + counts["partial"] + counts["missing"]
    print(f"Wrote {args.output} with {len(merged)} rows")
    print(
        "Status: "
        f"complete={counts['complete']} "
        f"failed={counts['failed']} "
        f"partial={counts['partial']} "
        f"missing={counts['missing']}"
    )
    print(f"Missing summaries: {sum(1 for r in merged if r['status'] == 'missing')}")

    if args.require_complete and incomplete:
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
