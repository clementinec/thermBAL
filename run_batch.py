#!/usr/bin/env python3
"""Run headless thermBAL batch studies."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List

from floor_plan.io import load_template_json
from floor_plan.study import ScenarioResult, load_study_definition, run_study


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run thermBAL batch studies from exported geometry or simulator state JSON.",
    )
    parser.add_argument(
        "--plan",
        help="Path to geometry.json or simulator plan_state.json. Overrides study.base_plan.",
    )
    parser.add_argument(
        "--study",
        required=True,
        help="Path to the study definition JSON.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output directory for results.",
    )
    return parser.parse_args()


def resolve_plan_path(cli_plan: str | None, study_base_plan: str | None, study_path: Path) -> Path:
    if cli_plan:
        return Path(cli_plan).expanduser().resolve()
    if study_base_plan:
        candidate = Path(study_base_plan).expanduser()
        if not candidate.is_absolute():
            candidate = (study_path.parent / candidate).resolve()
        return candidate
    raise SystemExit("A plan path is required. Pass --plan or set base_plan in the study JSON.")


def ensure_output_dir(path: str | Path) -> Path:
    out_dir = Path(path).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def flatten_record(record) -> Dict[str, object]:
    snap = record.snapshot
    position = snap.get("position", {})
    environment = snap.get("environment", {})
    computed = snap.get("computed", {})
    dynamic = snap.get("dynamic_state", {})
    coords = position.get("coordinates", [None, None])

    return {
        "scenario_id": record.scenario_id,
        "time_min": record.time_min,
        "agent_id": record.agent_id,
        "room_id": position.get("room_id"),
        "room_name": position.get("room_name"),
        "x": coords[0] if isinstance(coords, list) and len(coords) > 0 else None,
        "y": coords[1] if isinstance(coords, list) and len(coords) > 1 else None,
        "activity": dynamic.get("activity"),
        "clo": dynamic.get("clo"),
        "duration_mins": dynamic.get("duration_mins"),
        "pmv": computed.get("pmv"),
        "ppd": computed.get("ppd"),
        "met": computed.get("met"),
        "air_temp": environment.get("air_temp"),
        "mean_radiant_temp": environment.get("mean_radiant_temp"),
        "humidity": environment.get("humidity"),
        "air_velocity": environment.get("air_velocity"),
        "lux": environment.get("lux"),
        "noise_db": environment.get("noise_db"),
        "position_error": position.get("error"),
        "environment_error": environment.get("error"),
        "computed_error": computed.get("error"),
    }


def write_summary_csv(path: Path, results: Iterable[ScenarioResult]) -> None:
    rows = [result.summary for result in results]
    fieldnames = [
        "scenario_id",
        "description",
        "records",
        "agents",
        "times",
        "mean_pmv",
        "min_pmv",
        "max_pmv",
        "worst_ppd",
        "error_records",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            normalized = dict(row)
            normalized["times"] = ",".join(str(item) for item in row.get("times", []))
            writer.writerow(normalized)


def write_agent_timeseries_csv(path: Path, results: Iterable[ScenarioResult]) -> None:
    rows = [flatten_record(record) for result in results for record in result.records]
    fieldnames = [
        "scenario_id",
        "time_min",
        "agent_id",
        "room_id",
        "room_name",
        "x",
        "y",
        "activity",
        "clo",
        "duration_mins",
        "pmv",
        "ppd",
        "met",
        "air_temp",
        "mean_radiant_temp",
        "humidity",
        "air_velocity",
        "lux",
        "noise_db",
        "position_error",
        "environment_error",
        "computed_error",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_snapshots_jsonl(path: Path, results: Iterable[ScenarioResult]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for result in results:
            for record in result.records:
                handle.write(json.dumps(record.to_jsonable()))
                handle.write("\n")


def write_run_manifest(
    path: Path,
    plan_path: Path,
    study_path: Path,
    results: List[ScenarioResult],
) -> None:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "plan_path": str(plan_path),
        "study_path": str(study_path),
        "scenario_count": len(results),
        "scenarios": [result.summary for result in results],
        "outputs": {
            "summary_csv": "summary.csv",
            "agent_timeseries_csv": "agent_timeseries.csv",
            "snapshots_jsonl": "snapshots.jsonl",
        },
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def main() -> None:
    args = parse_args()
    study_path = Path(args.study).expanduser().resolve()
    study = load_study_definition(study_path)
    plan_path = resolve_plan_path(args.plan, study.base_plan, study_path)
    out_dir = ensure_output_dir(args.out)

    state = load_template_json(plan_path)
    results = run_study(state, study)

    write_summary_csv(out_dir / "summary.csv", results)
    write_agent_timeseries_csv(out_dir / "agent_timeseries.csv", results)
    write_snapshots_jsonl(out_dir / "snapshots.jsonl", results)
    write_run_manifest(out_dir / "run_manifest.json", plan_path, study_path, results)

    print("thermBAL batch run complete")
    print("  plan:   {0}".format(plan_path))
    print("  study:  {0}".format(study_path))
    print("  out:    {0}".format(out_dir))
    print("  files:  summary.csv, agent_timeseries.csv, snapshots.jsonl, run_manifest.json")


if __name__ == "__main__":
    main()
