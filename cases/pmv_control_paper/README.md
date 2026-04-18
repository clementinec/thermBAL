# PMV Control Paper

This folder contains a new journal-style manuscript pipeline for the PMV/pathway/control paper.

The paper is intentionally separate from the ACADIA material. It uses the current mechanism studies as its evidence core:

- `cases/pathway_discomfort`
- `cases/pmv_balance_relationship`

## Paper Claim

The manuscript argues that `PMV/PPD` are useful as bounded steady-state comfort indices near neutrality, but are physically underdetermined and physiologically incomplete as sole control objectives.

The paper does **not** claim:

- human-subject validation
- a dynamic physiology model
- a full HVAC simulation

It is a model-logic and control-interpretation paper.

## Build

1. Regenerate the supporting analyses if needed:

```bash
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pathway_discomfort/run_analysis.py
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pmv_balance_relationship/run_analysis.py
```

2. Build the manuscript assets:

```bash
/Users/hongshanguo/anaconda3/bin/python3.10 cases/pmv_control_paper/build_paper.py
```

3. Compile the LaTeX draft:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Run the LaTeX command from inside `cases/pmv_control_paper`.

## Outputs

The build script writes:

- `out/figures/*.png`
- `out/tables/*.tex`
- `out/generated/macros.tex`
- `out/generated/summary.json`

The manuscript source is:

- `main.tex`
- `references.bib`

## Current Structure

- `Figure 1`: PMV linearity and PPD mapping
- `Figure 2`: zero-PMV manifold
- `Figure 3`: one-PMV-step required air shifts, fixed load steps, and PPD changes
- `Figure 4`: profile-sensitive directed `1 C` PMV/PPD response and pathway watt shifts
- `Figure 5`: severe recovery control-move versus exergy comparison

Tables are generated from the current CSV outputs and are intended to be traceable to code rather than manually typed into the manuscript.
