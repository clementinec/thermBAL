# PMV Balance Relationship

This study pulls apart four things that are often bundled together too casually:

- `PMV` as a thermal-load mapping
- `PPD` as a dissatisfied-fraction mapping from PMV
- pathway decomposition as the way the environment changes body heat balance
- system-side exergy as the cost of producing those environmental changes

The aim is not to replace `pathway_watts` or `pathway_discomfort`. It is to clarify what PMV and PPD are actually doing in the model before we infer too much from them.

## Core Questions

1. In this implementation, what is the real relationship between PMV and thermal balance?
2. Does `PMV = 0` correspond to a unique environment, or to a manifold of possible environments?
3. What does PPD mean here, and what does it *not* mean?
4. If air and radiant systems both move PMV toward zero, is one more effective:
   - per body-load watt?
   - per system exergy watt?

## What This Case Computes

1. A sampled `PMV` vs `thermal load (H - L)` dataset across many environmental states and several metabolic rates.
2. The deterministic `PPD(PMV)` curve.
3. `PMV = 0` contour manifolds in `Ta-Tr` space under several humidity / air-speed conditions.
4. An `air vs radiant` recovery comparison on the same severe discomfort states used in `pathway_discomfort`, with:
   - recovered PMV
   - recovered PPD
   - dominant pathway shift
   - estimated system exergy surrogate
5. A profile-sensitive `directed 1 C` response study that asks how much a one-degree air or radiant move changes:
   - distance to neutrality
   - total thermal load
   - convective and radiative pathway watts in `W/person`
6. A one-PMV-step air-response study that asks what it takes to move:
   - `-3 -> -2`
   - `-2 -> -1`
   - `-1 -> 0`
   - and the hot-side mirrors
   across several environmental families, so the distinction between equal `ΔPMV`, equal `Δ(H-L)`, and unequal `ΔTa` becomes explicit.

## Important Boundary

This branch is not a validation study against measured TSV datasets.

The default office reference occupant in this branch is now set to `0.9 met`, and the sampled sensitivity band is `0.8 / 0.9 / 1.0 / 1.2 met`, to avoid quietly treating `1.0 met` as the uncontested sedentary office default.

It is an internal model-logic study that uses the repo's existing PMV implementation to clarify:

- what is mathematically exact inside the model
- what is empirical in the PMV / PPD mapping
- what only appears once system-side assumptions are added

## Outputs

Running the analysis writes:

- `out/pmv_load_samples.csv`
- `out/ppd_curve.csv`
- `out/zero_pmv_grid.csv`
- `out/zero_pmv_conditions.json`
- `out/air_radiant_effectiveness.csv`
- `out/profile_response_1c.csv`
- `out/pmv_step_air_response.csv`
- `out/report.html`
- `out/figures/pmv_vs_load.png`
- `out/figures/ppd_curve.png`
- `out/figures/zero_pmv_manifold.png`
- `out/figures/air_vs_radiant_effectiveness.png`
- `out/figures/profile_response_heatmaps.png`
- `out/figures/profile_pathway_heatmaps.png`
- `out/figures/pmv_step_required_air_shift.png`
- `out/figures/pmv_step_load_change.png`

## Run

Recommended interpreter in this repo:

```bash
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pmv_balance_relationship/run_analysis.py
```

## Why This Exists

There is a real conceptual distinction between:

- `PMV` moving because total thermal load changed
- `PMV` moving because a specific pathway was targeted
- a system being more or less effective in actual plant work / exergy

There is also a body-scale distinction that matters here:

- PMV in the standard formulation is built on `W/m²` heat balance, so body area does not directly change PMV once everything is normalized by area.
- But body area does change the total `W/person` attached to each pathway, which is why profile-sensitive `1 C` response maps are useful alongside the PMV plots.

And there is a label-step distinction:

- equal PMV steps imply equal thermal-load steps within a fixed profile, because PMV is linear in load for fixed `M`
- but they do not imply equal `ΔTa`, equal `ΔTr`, or equal pathway redistribution

This folder exists to make that separation explicit.
