# Pathway Discomfort

This folder is a severe-discomfort extension of [pathway_watts](/Users/hongshanguo/projs/thermBAL/cases/pathway_watts).

The goal is simpler and harsher:

- start from one representative office occupant
- construct environmental states where the occupant is plainly uncomfortable
- require each baseline to satisfy `|PMV| >= 2`
- test which single environmental lever recovers the state most directly
- inspect the recovery in pathway watts, not just in PMV

This is still a mechanism study, not a spatial study.

## What This Prototype Does

1. Uses the repo's existing ISO 7730 PMV engine for one seated office occupant.
2. Builds six clearly severe discomfort archetypes:
   - severe cold bulk
   - cold radiant sink
   - cold high-airflow draft
   - severe hot bulk
   - hot radiant gain
   - hot humid still air
3. Verifies that every baseline sits at `|PMV| >= 2`.
4. Solves bounded **single-lever** recoveries for:
   - air temperature only
   - mean radiant temperature only
   - airflow only
   - humidity only
5. Records:
   - baseline PMV / PPD
   - baseline pathway decomposition
   - shortest feasible recovery for each scenario by control displacement
   - lowest estimated system-exergy recovery for each scenario
   - dominant recovery channel in watts

## Important Boundary

This branch is intentionally simpler than `pathway_watts`.

It does **not** add local radiant-asymmetry or ankle-draft dissatisfaction models. Instead, it shifts the entire experiment to stronger baseline states where the PMV signal itself is already unambiguous. The point is to understand recovery under clearly severe cold and hot conditions before layering local criteria back in.

So the logic here is:

- baseline discomfort: `|PMV| >= 2`
- recovered state: `|PMV| <= 0.5`

It now also adds a second layer of explicit system-side cost surrogates:

- all-air sensible control: `Q = ρ c_p V̇ |ΔTa|`
- radiant control: `Q = h_r A |ΔTr|`
- airflow control: `P ∝ v^3`
- humidity control: `Q = ṁ h_fg |ΔW|`

These are converted into approximate electric / exergy input using fixed COP assumptions for a canonical office control zone. This is not a full HVAC model, but it avoids treating `ΔTa` and `ΔTr` as if they had the same plant implication.

## Outputs

Running the analysis writes:

- `out/reference_summary.json`
- `out/scenario_baselines.csv`
- `out/scenario_recoveries.csv`
- `out/scenario_winners.csv`
- `out/scenario_exergy_winners.csv`
- `out/report.html`
- `out/figures/extreme_baselines.png`
- `out/figures/recovery_matrix.png`
- `out/figures/system_cost_matrix.png`

## Run

Recommended interpreter in this repo:

```bash
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pathway_discomfort/run_analysis.py
```

## Why This Exists

`pathway_watts` proved the basic mechanism logic and the value of local criteria.

This folder asks a narrower question:

- when the occupant is already far into clear discomfort,
- which environmental lever gives the shortest path back,
- which heat-transfer channel actually carries that recovery,
- and whether the same ranking still holds once a simple system exergy layer is added.

That makes it easier to read the physics without mild cases blurring the result.
