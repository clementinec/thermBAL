# Apartment Daytime Cooled Case

This case turns the traced apartment template into a reproducible population sweep for the batch runner.

## Study Shape

- Plan: [`floor_plan/apartment_template.json`](../../floor_plan/apartment_template.json)
- Operation: daytime cooled apartment
- Orientation assumption: bottom of plan is south
- Occupancy rule: one occupant in each occupied habitable room
- Cohorts:
  - `young_mixed`
  - `young_male`
  - `older_shift`
- Climate sweep:
  - `Ta = 22.0, 23.5, 25.0 C`
  - `RH = 35, 55, 75 %`
- Spatial differentiation:
  - north/interior rooms use `Tr = Ta`
  - south-facing rooms use `Tr = Ta + 1.5 C`
  - east/west-exposed rooms use `Tr = Ta + 0.75 C`

This produces `27` scenarios total.

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

The current first-pass sweep behaves mostly like an overcooling study: all `27` scenarios remain on the cool side in aggregate, with the `older_shift` cohort showing the strongest cold-stress response. That is still useful for heterogeneity analysis, but if you want a broader neutral-to-warm crossover later, the next adjustment should be clothing/activity assumptions or a slightly warmer `Ta` sweep.
