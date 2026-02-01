# Comfort Engine (Python core)

Transparent PMV computation with traceable equations. This core expects **metabolic rate in W/m²** as input and exposes detailed intermediate values for auditability.

## What this package includes

- `ComfortState`: environmental inputs (Ta, Tr, RH, v, patm)
- `Occupant`: occupant inputs (**m_wm2**, clo, wme)
- `compute_pmv`: PMV/PPD solver with component breakdown and trace
- Metabolic helpers: `estimate_metabolic_profile`, `mifflin_bmr`, `du_bois_bsa`
- Psychrometrics: `dew_point_from_rh`, `humidity_ratio`, `sat_vapor_pressure`

## Quick start

```python
from comfort_engine import ComfortState, Occupant, compute_pmv

state = ComfortState(Ta=24.0, Tr=24.0, RH=50.0, v=0.1)
occupant = Occupant(m_wm2=70.0, clo=0.7, wme=0.0)

result = compute_pmv(state, occupant)

print(result.PMV, result.PPD)
print(result.q_conv, result.q_rad, result.e_sk, result.q_res)
```

## Using demographics to get M (W/m²)

```python
from comfort_engine import estimate_metabolic_profile, ComfortState, Occupant, compute_pmv

profile = estimate_metabolic_profile(
    sex="male",
    weight_kg=60,
    height_cm=175,
    age=65,
    rmr_factor=1.1,
    activity_multiplier=1.2,
)

# Use the derived M in W/m² directly
m_wm2 = profile["M_active_Wm2"]

state = ComfortState(Ta=24.0, Tr=24.0, RH=50.0, v=0.1)
occupant = Occupant(m_wm2=m_wm2, clo=0.7)
result = compute_pmv(state, occupant)

print(result.M, result.PMV, result.PPD)
```

## Met-equivalent (for reference only)

If you have **met units** and want W/m²:

```python
from comfort_engine import MET_TO_WM2

m_wm2 = 1.2 * MET_TO_WM2
```

The solver also provides a derived `met_equiv` value in `result.inputs` and the trace for QA checks.

## Sanity check (typical baseline)

At Ta = Tr = 25 C, RH = 50%, v = 0.1 m/s, and M = 58.15 W/m^2 (1 met), PMV depends strongly on clothing. This implementation yields approximately:

- `clo = 0.0` -> PMV ~ -2.32
- `clo = 0.5` -> PMV ~ -0.40
- `clo = 0.7` -> PMV ~ 0.01

Use this as a quick check that units and clothing inputs are sensible.

## Trace and diagnostics

```python
for step in result.trace:
    print(step.id, step.equation)

print(result.warnings)
print(result.solver)
```

## Notes

- PMV inputs are per ASHRAE/ISO conventions, but **M is taken in W/m²**.
- `wme` is external work in W/m² (defaults to 0).
- `Occupant.body_area_m2` (optional) enables component totals in W.
