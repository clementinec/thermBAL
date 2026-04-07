# Apartment Daytime Cooled Case

This case turns the traced apartment template into a reproducible population sweep for the batch runner.

## Study Shape

- Plan: [`floor_plan/apartment_template.json`](../../floor_plan/apartment_template.json)
- Operation: daytime cooled apartment
- Orientation assumption: bottom of plan is south
- Cohorts:
  - `young_mixed`
  - `young_male`
  - `older_shift`
  - `higher_clo_sedentary`
  - `lighter_clothing_mobile`
  - `higher_bmi_warm_sensitive`
- Climate sweep:
  - `Ta = 20.0, 22.0, 23.5, 25.0, 28.0 C`
  - `RH = 35, 55, 75, 80 %`
- Spatial differentiation:
  - north/interior rooms use `Tr = Ta`
  - south-facing rooms use `Tr = Ta + 1.5 C`
  - east/west-exposed rooms use `Tr = Ta + 0.75 C`
- Occupancy rule:
  - one occupant in each enclosed room except `Entry Gallery`, `Service Hall`, and `East Passage`

This produces `120` scenarios total with `21` placed occupants per cohort.

## Files

- [`generate_study.py`](./generate_study.py)
  Builds the deterministic study JSON and cohort manifest.
- [`study.json`](./study.json)
  Generated batch study definition.
- [`study_manifest.json`](./study_manifest.json)
  Cohort metadata, occupied-room map, and scenario catalog.
- [`render_report.py`](./render_report.py)
  Renders the stored outputs into a dependency-free HTML report.
- [`out/summary.csv`](./out/summary.csv)
- [`out/agent_timeseries.csv`](./out/agent_timeseries.csv)
- [`out/snapshots.jsonl`](./out/snapshots.jsonl)
- [`out/run_manifest.json`](./out/run_manifest.json)
- [`out/report.html`](./out/report.html)

## Rebuild

From the repo root:

```bash
python cases/apartment_daytime_cooled/generate_study.py
python run_batch.py \
  --study cases/apartment_daytime_cooled/study.json \
  --out cases/apartment_daytime_cooled/out
python cases/apartment_daytime_cooled/render_report.py
```

## Interpretation Note

The current six-cohort rerun spans both cold and slightly warm cases:

- global mean PMV now ranges from roughly `-5.17` to `0.91`
- worst-room PPD ranges from `19.1%` to `100.0%`
- `older_shift` remains the coldest cohort overall and never reaches positive mean PMV
- `higher_clo_sedentary` is currently the warmest cohort in the sweep, peaking at `higher_clo_sedentary_t28p0_rh80`
- the apartment still shows strong cold-stress sensitivity at `20 C / 35 % RH`
- report colors now diverge by PMV sign:
  - blue for cold-side stress (`PMV < 0`)
  - beige for near-neutral
  - warm orange/red for hot-side stress (`PMV > 0`)

That gives a better heterogeneity story than the earlier narrow sweep because the case now crosses from strong cold stress into near-neutral or mildly warm conditions depending on cohort, while keeping hot and cold discomfort visually distinct in the report.
