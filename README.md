# thermBAL

thermBAL is a floor-plan-based thermal comfort toolkit for authoring simplified building layouts, assigning room conditions and occupant profiles, and studying agent comfort across space.

It combines:

- a traceable PMV/PPD comfort engine in Python
- a browser-based floor plan simulator
- a PDF-backed geometry importer
- a Rhino exporter for `geometry.json`
- a headless batch runner for repeatable scenario studies

This repository is aimed at fast spatial comfort studies on simplified plans, not full building energy simulation.

## What Is In This Repo

- `comfort_engine/`
  Python PMV/PPD core with psychrometrics, metabolic helpers, and detailed trace output.
- `floor_plan/`
  Grid-based floor plan models, agent runtime, JSON I/O, presets, and study execution helpers.
- `simulator.html`
  Interactive browser simulator for editing footprint, walls, openings, rooms, and occupant states.
- `import-geometry.html`
  Assisted PDF-backed tracing workflow that outputs simulator-ready geometry.
- `rhino_export/`
  RhinoPython exporter for taking simplified Rhino geometry straight to `geometry.json`.
- `run_batch.py`
  Headless CLI for scenario replay and batch output generation.
- `batch/`
  Batch workflow docs and example study JSON.
- `index.html`
  Browser comfort explainer / formula viewer for the comfort core.

## Main Data Artifacts

The repo is easiest to use if you keep these as separate layers:

- `geometry.json`
  Geometry only: footprint cells, walls, openings, and optional room seeds.
- `plan_state.json`
  Geometry plus room environmental settings and optional placed agents.
- `study.json`
  Batch scenario definitions: room overrides, agent distributions, and timed events.

Recommended split:

1. Author geometry.
2. Set room conditions and baseline agents.
3. Run scenarios separately.

## Environment

Python workflows:

- Python `3.10+`
- no extra pip dependencies for the current Python tools

Browser workflows:

- modern browser
- serving the repo over HTTP is recommended, especially for bundled PDF loading
- `import-geometry.html` uses `pdf.js` from a CDN
- `index.html` uses KaTeX from a CDN

Simple local server:

```bash
python -m http.server 8000
```

Then open:

- `http://localhost:8000/simulator.html`
- `http://localhost:8000/import-geometry.html`
- `http://localhost:8000/index.html`

Optional upstream authoring:

- Rhino 7/8 with RhinoPython for the Rhino exporter workflow

## Quick Starts

### 1. Python Comfort Core

```python
from comfort_engine import ComfortState, Occupant, compute_pmv

state = ComfortState(Ta=24.0, Tr=24.0, RH=50.0, v=0.1)
occupant = Occupant(m_wm2=70.0, clo=0.7, wme=0.0)

result = compute_pmv(state, occupant)

print(result.PMV, result.PPD)
print(result.q_conv, result.q_rad, result.e_sk, result.q_res)
```

If you want metabolic rate from demographics rather than specifying `M` directly:

```python
from comfort_engine import estimate_metabolic_profile

profile = estimate_metabolic_profile(
    sex="male",
    weight_kg=60,
    height_cm=175,
    age=65,
    rmr_factor=1.1,
    activity_multiplier=1.2,
)

print(profile["M_active_Wm2"])
```

### 2. Browser Geometry And Simulator Workflow

1. Open `import-geometry.html`.
2. Load a PDF or the sample plan.
3. Auto-detect the candidate footprint.
4. Refine footprint, interior walls, windows, doors, and exits.
5. `Download Geometry JSON` or `Send To Simulator`.
6. In `simulator.html`, review rooms, set room conditions, add agents, and export state.

The simulator can also load geometry directly through `Load Geometry JSON`.

### 3. Rhino Workflow

1. Simplify the plan in Rhino to the level the simulator needs.
2. Put geometry on the expected layers:
   `tb_footprint`, `tb_voids`, `tb_walls`, `tb_windows`, `tb_doors`, `tb_exits`, `tb_rooms`, `tb_origin`
3. Run `rhino_export/thermbal_rhino_exporter.py` in Rhino.
4. Save the generated `geometry.json`.
5. Load that JSON in `simulator.html`.

See [`rhino_export/README.md`](rhino_export/README.md) for the layer contract and room metadata keys.

### 4. Headless Batch Workflow

The batch runner is the non-GUI study path.

1. Export `plan_state.json` from the simulator, or start directly from `geometry.json`.
2. Define scenarios in JSON.
3. Run the CLI.
4. Analyze the outputs downstream in Python, Excel, Grasshopper, or another pipeline.

Example:

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

See [`batch/README.md`](batch/README.md) for the study schema and workflow details.

## Python Package Structure

### `comfort_engine`

Primary exports:

- `ComfortState`
- `Occupant`
- `compute_pmv`
- `estimate_metabolic_profile`
- `mifflin_bmr`
- `du_bois_bsa`
- `dew_point_from_rh`
- `humidity_ratio`
- `sat_vapor_pressure`

This layer is the traceable thermal comfort core.

### `floor_plan`

Primary exports:

- `FloorPlan`, `Room`, `RoomEnvironment`
- `Agent`, `AgentDemographics`, `AgentPreferences`
- `Simulator`
- `load_geometry_json`, `load_plan_state_json`, `save_plan_state_json`
- `Scenario`, `ScheduledEvent`, `StudyDefinition`
- `run_scenario`, `run_study`
- preset helpers like `default_floor_plan()` and `demo_simulator()`

Quick demo:

```bash
python demo_floorplan.py
```

## Current Scope

What thermBAL does now:

- traceable PMV/PPD comfort calculation in Python
- psychrometric helper functions
- metabolic-rate estimation from demographics
- browser editing of footprint, walls, openings, rooms, and occupants
- assisted PDF-backed geometry tracing
- Rhino export to simulator-ready geometry JSON
- static comfort snapshots for positioned agents
- scheduled batch replay with room overrides, movement, clothing changes, and activity changes

## Important Scope Boundary On Latent Behavior

This repo already includes latent and evaporative terms inside the comfort calculation itself through the PMV heat-balance model and its psychrometric inputs.

What it does not yet include is transient room-level latent-load carryover over time, such as occupants gradually changing room moisture or thermal state and that updated room state feeding the next timestep automatically.

That future layer belongs in the Python simulation pipeline, not in the browser GUI first.

## Not Yet Implemented

- carried-over room thermal or moisture loads from occupants over time
- HVAC or ventilation dynamics
- probabilistic agent placement sweeps as a first-class study primitive
- robust automatic understanding of arbitrary messy floor plan PDFs
- clean support for highly angled or curved geometry without simplification to the working grid

## Suggested Workflow

For most projects, the practical path is:

1. Author geometry in Rhino or `import-geometry.html`.
2. Review room conditions and baseline occupants in `simulator.html`.
3. Export `plan_state.json`.
4. Define scenario sets in batch JSON.
5. Run `run_batch.py`.
6. Use the browser tools again for inspection and communication, not as the batch engine.

## Related Docs

- [`batch/README.md`](batch/README.md)
- [`rhino_export/README.md`](rhino_export/README.md)
