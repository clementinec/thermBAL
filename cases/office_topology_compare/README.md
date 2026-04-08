# Office Topology Comparison

This case reuses the apartment footprint as an office floor and compares two derived templates:

- `floor_plan/office_cellular_template.json`
  The apartment rooms are re-labeled as enclosed office rooms.
- `floor_plan/office_openplan_template.json`
  The same footprint is cleared into two major open-office zones with a retained core/wet spine and a west corner office.

The comparison is intentionally controlled:

- same outer footprint
- same window/perimeter geometry
- same deterministic desk/sample points
- same six office cohorts with sex-weighted clothing and age ranges
- same `Ta x RH` sweep used in the apartment case

What changes is how the thermal field is resolved:

- `cellular_office`
  The perimeter-to-window field is averaged back to each enclosed room.
- `open_office`
  The same field is preserved positionally across the open floor.

That makes this case useful for showing how topology changes the legibility of spatial thermal heterogeneity.

## Generate

From the repo root:

```bash
python cases/office_topology_compare/generate_study.py
python cases/office_topology_compare/render_report.py
```

## Outputs

- `study_manifest.json`
  Case assumptions, template references, sample points, and cohort metadata.
- `out/summary.csv`
  Topology/cohort/scenario summaries for the sampled desk points.
- `out/comparison_summary.csv`
  Open-minus-cellular deltas for PMV range, spatial spread, and perimeter-core gap.
- `out/sample_results.csv`
  Per-sample comfort results for all scenario rows.
- `out/field_maps.json`
  Full-cell PMV maps for selected reference cases.
- `out/report.html`
  Visual comparison report with topology diagrams and field maps.

## Current Scope

This case uses a custom position-based field model for `Ta` and `MRT`. It does not model occupant-to-room thermal feedback or transient load carry-over.

The highlighted stress case in the report is `22 C / 80% RH`, where the mixed, male-formal, and female-light office cohorts separate most clearly.
