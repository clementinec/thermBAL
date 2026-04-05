# Batch Workflow

This is the non-GUI study pipeline for thermBAL.

The browser tools and Rhino exporter are for authoring. The batch runner is for reproducible scenario studies, sweeps, and scheduled replay.

## Environment

Use normal terminal Python, not RhinoPython and not the browser runtime.

Recommended:

- Python `3.10+`
- no extra pip dependencies
- standard library only

Optional setup:

```bash
python -m venv .venv
source .venv/bin/activate
```

No `pip install` step is required for the current repo layout.

## Inputs

The batch runner accepts either:

- `geometry.json`
  Geometry only. Good when agents are defined entirely in the study file.
- `plan_state.json`
  Exported from [simulator.html](../simulator.html) using `Export State`. This includes geometry, room environments, and any current agents.

Typical upstream paths:

1. Rhino -> [thermbal_rhino_exporter.py](../rhino_export/thermbal_rhino_exporter.py) -> `geometry.json`
2. Browser geometry importer -> simulator -> `Export State` -> `plan_state.json`

## Run A Study

Use [run_batch.py](../run_batch.py):

```bash
python run_batch.py \
  --plan path/to/thermBAL-plan-state.json \
  --study batch/example_study.json \
  --out out/example_run
```

Outputs:

- `summary.csv`
- `agent_timeseries.csv`
- `snapshots.jsonl`
- `run_manifest.json`

## Study File Shape

See [example_study.json](example_study.json).

Top-level keys:

- `base_plan`
  Optional path to the plan JSON. If omitted, pass `--plan`.
- `scenarios`
  List of scenarios to run.

Scenario keys:

- `id`
- `description`
- `use_base_agents`
  Default `true`. If the plan JSON already contains agents, keep them unless you turn this off.
- `agents`
  Add or override agents.
- `room_overrides`
  Override room environmental values before the scenario starts.
- `schedule`
  Timed events such as moves or activity changes.
- `snapshot_times`
  Additional times to capture even if no events happen there.

## Agent Specs

Agent specs can use built-in presets:

- `elder_walker`
- `young_visitor`
- `seated_elder`

Example:

```json
{
  "id": "persona_10",
  "profile": "elder_walker",
  "position": [4, 7],
  "activity": "walking",
  "clo": 0.8
}
```

You can also override demographics and preferences directly:

```json
{
  "id": "persona_12",
  "demographics": {
    "age": 81,
    "mobility": "uses_walker"
  },
  "preferences": {
    "noise_tolerance_db": 48
  },
  "position": [3, 6]
}
```

## Schedule Events

Example:

```json
{
  "t_min": 30,
  "agent_id": "persona_01",
  "move_to": [10, 6],
  "activity": "walking"
}
```

Supported event fields:

- `t_min`
- `agent_id`
- `move_to`
- `activity`
- `clo`
- `demographics`
- `preferences`
- `room_overrides`

Events at the same `t_min` are applied together, then snapshots are recorded.

## How To Think About The Workflow

Recommended split:

1. Author geometry in Rhino or the browser importer.
2. If needed, review room conditions and occupant setup in the browser simulator.
3. Export `plan_state.json`.
4. Define study scenarios in JSON.
5. Run `run_batch.py`.
6. Analyze CSV/JSONL outputs in Python, Excel, Grasshopper, or another downstream workflow.

## Current Scope

What the batch runner does now:

- loads current thermBAL geometry or simulator-state JSON
- applies room overrides
- places or overrides agents
- replays scheduled movement / activity / clothing changes
- records deterministic snapshots over time

What it does not do yet:

- carried-over room thermal loads from occupants
- HVAC / ventilation dynamics
- automatic placement sweeps from probabilistic distributions

Those should be added as a Python engine layer later, not as GUI logic.
