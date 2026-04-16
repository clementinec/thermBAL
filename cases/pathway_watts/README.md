# Pathway Watts

This folder is a fresh exploratory branch for the `Watts-first` paper direction discussed around [pathway_01.pdf](/Users/hongshanguo/projs/thermBAL/pathway_01.pdf).

The intent is to test a simpler paper core than `thermBAL`:

- start from one representative office occupant
- benchmark the reference heat-balance pathways
- perturb the environment through a small set of canonical discomfort states
- test which single environmental lever most directly recovers PMV and local discomfort acceptability
- inspect the recovery in pathway watts, not just in PMV

This is **not** a spatial paper. It is a mechanism paper scaffold.

## What This Prototype Does

1. Computes a reference office benchmark using the repo's existing comfort engine.
2. Computes a local finite-difference sensitivity table for:
   - `Ta`
   - `Tr`
   - `RH`
   - `v`
   - `clo`
   - `met`
3. Evaluates six canonical discomfort proxies:
   - whole room too cold
   - cold radiant asymmetry
   - ankle draft
   - mild warm bulk
   - hot radiant asymmetry
   - warm humid still air
4. Solves bounded **single-lever** recoveries for:
   - air temperature only
   - mean radiant temperature only
   - airflow only
   - humidity only
5. Records which pathway channel changes the most in each feasible recovery:
   - convection
   - radiation
   - evaporation
   - respiration

## Important Boundary

This prototype uses the current PMV engine plus explicit local discomfort checks.

That means:

- bulk scenarios are judged by `|PMV| <= 0.5`
- radiant asymmetry uses ISO-style dissatisfied-percentage curves for:
  - `cool wall`
  - `warm wall`
- draft uses an ankle-draft dissatisfied-percentage model with whole-body thermal sensation approximated by PMV

So this is now a stronger mechanism prototype, while still remaining intentionally low-order and occupant-side.

## Outputs

Running the analysis writes:

- `out/reference_summary.json`
- `out/pathway_sensitivities.csv`
- `out/scenario_baselines.csv`
- `out/scenario_recoveries.csv`
- `out/scenario_winners.csv`
- `out/report.html`
- `out/figures/reference_pathways.png`
- `out/figures/pathway_sensitivities.png`
- `out/figures/recovery_matrix.png`

## Run

Recommended interpreter in this repo:

```bash
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pathway_watts/run_analysis.py
```

## Why This Exists

The question here is not whether the solver works. It already does.

The question is whether the paper is really about:

- a scenario catalog, or
- a cleaner claim that discomfort and recovery belong to different heat-transfer pathways

This folder is meant to answer that second question quickly.
