#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
OUT = CASE / "out"
FIGS = OUT / "figures"
MPLCONFIGDIR = OUT / ".mplconfig"

os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from comfort_engine import ComfortState, Occupant, compute_pmv, du_bois_bsa, met_to_wm2


@dataclass(frozen=True)
class Scenario:
    id: str
    label: str
    description: str
    state: ComfortState


@dataclass(frozen=True)
class Lever:
    id: str
    label: str
    unit: str
    values: Tuple[float, ...]
    apply: Callable[[ComfortState, float], ComfortState]


REFERENCE_OCCUPANT = {
    "label": "Representative office occupant",
    "sex": "male",
    "age": 40,
    "height_cm": 170.0,
    "weight_kg": 70.0,
    "clo": 0.5,
    "met": 1.0,
}

REFERENCE_STATE = ComfortState(Ta=25.5, Tr=25.5, RH=50.0, v=0.10)

SENSITIVITY_STEPS = {
    "Ta": 0.5,
    "Tr": 0.5,
    "RH": 5.0,
    "v": 0.05,
    "clo": 0.05,
    "met": 0.10,
}

SCENARIOS: Tuple[Scenario, ...] = (
    Scenario(
        id="cool_bulk",
        label="Whole room too cold",
        description="Bulk cool proxy: air and radiant temperature both reduced.",
        state=ComfortState(Ta=21.5, Tr=21.5, RH=50.0, v=0.10),
    ),
    Scenario(
        id="cold_radiant_skew",
        label="Cold-radiant skew",
        description="Whole-body MRT proxy for a cold-wall condition; not a local asymmetry model.",
        state=ComfortState(Ta=25.5, Tr=21.5, RH=50.0, v=0.10),
    ),
    Scenario(
        id="draft",
        label="Draft / elevated air movement",
        description="Cool-side discomfort induced by elevated air speed.",
        state=ComfortState(Ta=26.0, Tr=26.0, RH=50.0, v=0.75),
    ),
    Scenario(
        id="warm_bulk",
        label="Mild warm bulk",
        description="Uniform warm-bulk state where multiple levers may recover similar PMV.",
        state=ComfortState(Ta=28.0, Tr=28.0, RH=50.0, v=0.10),
    ),
    Scenario(
        id="warm_humid_still",
        label="Warm humid still air",
        description="Warm-humid state with low air movement and evaporative constraint.",
        state=ComfortState(Ta=28.5, Tr=28.5, RH=70.0, v=0.05),
    ),
)


def _frange(start: float, stop: float, step: float) -> Tuple[float, ...]:
    n = int(round((stop - start) / step))
    return tuple(round(start + i * step, 6) for i in range(n + 1))


LEVERS: Tuple[Lever, ...] = (
    Lever(
        id="air_temp",
        label="All-air (Ta only)",
        unit="C",
        values=_frange(18.0, 32.0, 0.1),
        apply=lambda state, value: ComfortState(Ta=value, Tr=state.Tr, RH=state.RH, v=state.v),
    ),
    Lever(
        id="radiant",
        label="Radiant (Tr only)",
        unit="C",
        values=_frange(18.0, 32.0, 0.1),
        apply=lambda state, value: ComfortState(Ta=state.Ta, Tr=value, RH=state.RH, v=state.v),
    ),
    Lever(
        id="airflow",
        label="Airflow control (v only)",
        unit="m/s",
        values=_frange(0.0, 1.0, 0.01),
        apply=lambda state, value: ComfortState(Ta=state.Ta, Tr=state.Tr, RH=state.RH, v=value),
    ),
    Lever(
        id="humidity",
        label="Latent control (RH only)",
        unit="%",
        values=_frange(20.0, 80.0, 1.0),
        apply=lambda state, value: ComfortState(Ta=state.Ta, Tr=state.Tr, RH=value, v=state.v),
    ),
)


def build_occupant(*, met: float | None = None, clo: float | None = None) -> Occupant:
    bsa = du_bois_bsa(REFERENCE_OCCUPANT["weight_kg"], REFERENCE_OCCUPANT["height_cm"])
    return Occupant(
        m_wm2=met_to_wm2(met if met is not None else REFERENCE_OCCUPANT["met"]),
        clo=clo if clo is not None else REFERENCE_OCCUPANT["clo"],
        posture="seated",
        body_area_m2=bsa,
    )


def evaluate(state: ComfortState, occupant: Occupant) -> Dict[str, float]:
    result = compute_pmv(state, occupant)
    return {
        "Ta": state.Ta,
        "Tr": state.Tr,
        "RH": state.RH,
        "v": state.v,
        "PMV": result.PMV,
        "PPD": result.PPD,
        "Q_conv": result.q_conv,
        "Q_rad": result.q_rad,
        "E_sk": result.e_sk,
        "Q_res": result.q_res,
        "L_total": result.L,
        "h_c": result.h_c,
        "Tcl": result.Tcl,
        "met": result.inputs["met_equiv"],
        "clo": result.inputs["clo"],
    }


def finite_difference(var: str, step: float, state: ComfortState, occupant: Occupant) -> Dict[str, float]:
    if var in {"Ta", "Tr", "RH", "v"}:
        kwargs_plus = {"Ta": state.Ta, "Tr": state.Tr, "RH": state.RH, "v": state.v}
        kwargs_minus = dict(kwargs_plus)
        kwargs_plus[var] += step
        kwargs_minus[var] -= step
        plus = evaluate(ComfortState(**kwargs_plus), occupant)
        minus = evaluate(ComfortState(**kwargs_minus), occupant)
    elif var == "clo":
        plus = evaluate(state, build_occupant(clo=REFERENCE_OCCUPANT["clo"] + step))
        minus = evaluate(state, build_occupant(clo=REFERENCE_OCCUPANT["clo"] - step))
    elif var == "met":
        plus = evaluate(state, build_occupant(met=REFERENCE_OCCUPANT["met"] + step))
        minus = evaluate(state, build_occupant(met=REFERENCE_OCCUPANT["met"] - step))
    else:
        raise ValueError(f"Unknown sensitivity variable: {var}")

    fields = ["PMV", "PPD", "Q_conv", "Q_rad", "E_sk", "Q_res", "L_total"]
    return {field: (plus[field] - minus[field]) / (2.0 * step) for field in fields}


def dominant_channel(delta: Dict[str, float]) -> str:
    channel_map = {
        "convection": abs(delta["dQ_conv"]),
        "radiation": abs(delta["dQ_rad"]),
        "evaporation": abs(delta["dE_sk"]),
        "respiration": abs(delta["dQ_res"]),
    }
    return max(channel_map, key=channel_map.get)


def recover_with_lever(
    scenario: Scenario,
    lever: Lever,
    occupant: Occupant,
    pmv_target: float = 0.5,
) -> Dict[str, object]:
    baseline = evaluate(scenario.state, occupant)
    feasible: List[Tuple[float, Dict[str, float], float]] = []

    scenario_value = {
        "air_temp": scenario.state.Ta,
        "radiant": scenario.state.Tr,
        "airflow": scenario.state.v,
        "humidity": scenario.state.RH,
    }[lever.id]

    for value in lever.values:
        candidate_state = lever.apply(scenario.state, value)
        candidate = evaluate(candidate_state, occupant)
        if abs(candidate["PMV"]) <= pmv_target:
            feasible.append((abs(value - scenario_value), candidate, value))

    if not feasible:
        return {
            "scenario_id": scenario.id,
            "scenario_label": scenario.label,
            "lever_id": lever.id,
            "lever_label": lever.label,
            "lever_unit": lever.unit,
            "feasible": False,
            "baseline_value": scenario_value,
            "recovered_value": None,
            "delta_value": None,
            "baseline_pmv": baseline["PMV"],
            "recovered_pmv": None,
            "baseline_ppd": baseline["PPD"],
            "recovered_ppd": None,
            "dominant_channel": None,
            "dQ_conv": None,
            "dQ_rad": None,
            "dE_sk": None,
            "dQ_res": None,
        }

    delta_abs, recovered, chosen_value = min(feasible, key=lambda item: item[0])
    pathway_delta = {
        "dQ_conv": recovered["Q_conv"] - baseline["Q_conv"],
        "dQ_rad": recovered["Q_rad"] - baseline["Q_rad"],
        "dE_sk": recovered["E_sk"] - baseline["E_sk"],
        "dQ_res": recovered["Q_res"] - baseline["Q_res"],
    }
    return {
        "scenario_id": scenario.id,
        "scenario_label": scenario.label,
        "lever_id": lever.id,
        "lever_label": lever.label,
        "lever_unit": lever.unit,
        "feasible": True,
        "baseline_value": scenario_value,
        "recovered_value": chosen_value,
        "delta_value": chosen_value - scenario_value,
        "baseline_pmv": baseline["PMV"],
        "recovered_pmv": recovered["PMV"],
        "baseline_ppd": baseline["PPD"],
        "recovered_ppd": recovered["PPD"],
        "dominant_channel": dominant_channel(pathway_delta),
        **pathway_delta,
    }


def write_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def save_reference_plot(reference: Dict[str, float]) -> Path:
    labels = ["Convection", "Radiation", "Evaporation", "Respiration"]
    values = [reference["Q_conv"], reference["Q_rad"], reference["E_sk"], reference["Q_res"]]
    colors = ["#5B8FF9", "#F6A04D", "#74C0A3", "#9AA5B1"]

    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=180)
    bars = ax.bar(labels, values, color=colors, edgecolor="#24313A", linewidth=0.8)
    ax.set_ylabel("W/m²")
    ax.set_title("Reference Office Pathway Balance")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.5, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    out = FIGS / "reference_pathways.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_sensitivity_heatmap(rows: List[Dict[str, object]]) -> Path:
    variables = [row["variable"] for row in rows]
    metrics = ["PMV", "Q_conv", "Q_rad", "E_sk", "Q_res", "L_total"]
    matrix = np.array([[float(row[m]) for m in metrics] for row in rows], dtype=float)

    # Normalize per column to retain legibility across units.
    norm = np.zeros_like(matrix)
    for j in range(matrix.shape[1]):
        denom = max(1e-9, np.max(np.abs(matrix[:, j])))
        norm[:, j] = matrix[:, j] / denom

    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=180)
    im = ax.imshow(norm, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(metrics)), metrics)
    ax.set_yticks(range(len(variables)), variables)
    ax.set_title("Local Pathway Sensitivity at the Reference Office State")
    for i in range(len(variables)):
        for j in range(len(metrics)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=7, color="#1F1F1F")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="column-normalized sign / magnitude")
    fig.tight_layout()
    out = FIGS / "pathway_sensitivities.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def render_report(
    reference: Dict[str, float],
    sensitivity_rows: List[Dict[str, object]],
    baselines: List[Dict[str, object]],
    recoveries: List[Dict[str, object]],
) -> None:
    recovery_groups: Dict[str, List[Dict[str, object]]] = {}
    for row in recoveries:
        recovery_groups.setdefault(str(row["scenario_id"]), []).append(row)

    def html_escape(text: object) -> str:
        return (
            str(text)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
        )

    style = """
    body { font-family: Georgia, 'Times New Roman', serif; margin: 32px auto; max-width: 1080px; color: #222; line-height: 1.45; }
    h1, h2, h3 { font-family: Helvetica, Arial, sans-serif; color: #1f2b33; }
    .lede { font-size: 18px; max-width: 900px; }
    .note { color: #5b6670; font-size: 14px; }
    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; margin: 24px 0 32px; align-items: start; }
    .card { border: 1px solid #d8dde1; border-radius: 10px; padding: 18px; background: #fbfcfd; }
    table { border-collapse: collapse; width: 100%; margin: 10px 0 24px; font-size: 14px; }
    th, td { border-bottom: 1px solid #e1e6ea; padding: 8px 10px; text-align: left; vertical-align: top; }
    th { background: #f3f6f8; font-family: Helvetica, Arial, sans-serif; }
    img { width: 100%; border: 1px solid #d8dde1; border-radius: 10px; background: white; }
    .mono { font-family: 'SFMono-Regular', Consolas, monospace; font-size: 13px; }
    .ok { color: #1d7f44; font-weight: 600; }
    .no { color: #9a2c2c; font-weight: 600; }
    """

    baseline_rows_html = "\n".join(
        f"<tr><td>{html_escape(row['label'])}</td><td>{row['Ta']:.1f}</td><td>{row['Tr']:.1f}</td><td>{row['RH']:.0f}</td><td>{row['v']:.2f}</td><td>{row['PMV']:.2f}</td><td>{row['PPD']:.1f}</td><td>{row['dominant_loss']}</td></tr>"
        for row in baselines
    )

    sens_rows_html = "\n".join(
        f"<tr><td>{row['variable']}</td><td>{row['step']}</td><td>{row['unit']}</td><td>{row['PMV']:.3f}</td><td>{row['Q_conv']:.3f}</td><td>{row['Q_rad']:.3f}</td><td>{row['E_sk']:.3f}</td><td>{row['Q_res']:.3f}</td></tr>"
        for row in sensitivity_rows
    )

    scenario_sections = []
    for base in baselines:
        rows = recovery_groups.get(str(base["scenario_id"]), [])
        table_rows = []
        for row in rows:
            if row["feasible"]:
                delta_text = f"{row['delta_value']:+.2f} {row['lever_unit']}"
                pmv_text = f"{row['baseline_pmv']:.2f} → {row['recovered_pmv']:.2f}"
                channel_text = html_escape(row["dominant_channel"])
            else:
                delta_text = "—"
                pmv_text = f"{row['baseline_pmv']:.2f} → —"
                channel_text = "—"
            table_rows.append(
                "<tr>"
                f"<td>{html_escape(row['lever_label'])}</td>"
                f"<td class=\"{'ok' if row['feasible'] else 'no'}\">{'yes' if row['feasible'] else 'no'}</td>"
                f"<td>{delta_text}</td>"
                f"<td>{pmv_text}</td>"
                f"<td>{channel_text}</td>"
                "</tr>"
            )
        scenario_sections.append(
            "<div class='card'>"
            f"<h3>{html_escape(base['label'])}</h3>"
            f"<p class='note'>{html_escape(base['description'])}</p>"
            f"<p><span class='mono'>Ta={base['Ta']:.1f} C, Tr={base['Tr']:.1f} C, RH={base['RH']:.0f}%, v={base['v']:.2f} m/s, PMV={base['PMV']:.2f}</span></p>"
            "<table>"
            "<thead><tr><th>Lever</th><th>Feasible</th><th>Minimum control move</th><th>PMV recovery</th><th>Dominant recovery channel</th></tr></thead>"
            "<tbody>"
            + "\n".join(table_rows)
            + "</tbody></table></div>"
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Pathway Watts Prototype</title>
  <style>{style}</style>
</head>
<body>
  <h1>Pathway Watts Prototype</h1>
  <p class="lede">This is a deliberately stripped mechanism study built on the existing comfort engine. The question is whether a <em>Watts-first</em> paper core emerges once the problem is reduced to one office occupant, a few canonical discomfort states, and single-lever recoveries judged by PMV acceptability.</p>

  <div class="grid">
    <div class="card">
      <h2>Reference benchmark</h2>
      <p><span class="mono">Ta={reference['Ta']:.1f} C, Tr={reference['Tr']:.1f} C, RH={reference['RH']:.0f}%, v={reference['v']:.2f} m/s, met={reference['met']:.2f}, clo={reference['clo']:.2f}</span></p>
      <table>
        <tbody>
          <tr><th>PMV</th><td>{reference['PMV']:.2f}</td></tr>
          <tr><th>PPD</th><td>{reference['PPD']:.1f}%</td></tr>
          <tr><th>Convective loss</th><td>{reference['Q_conv']:.2f} W/m²</td></tr>
          <tr><th>Radiative loss</th><td>{reference['Q_rad']:.2f} W/m²</td></tr>
          <tr><th>Evaporative loss</th><td>{reference['E_sk']:.2f} W/m²</td></tr>
          <tr><th>Respiratory loss</th><td>{reference['Q_res']:.2f} W/m²</td></tr>
        </tbody>
      </table>
      <p class="note">At this reference office state, the interesting question is not just whether PMV is neutral, but which pathway is already carrying the larger share of the occupant heat balance.</p>
    </div>
    <div>
      <img src="figures/reference_pathways.png" alt="Reference pathway plot">
    </div>
  </div>

  <div class="grid">
    <div>
      <img src="figures/pathway_sensitivities.png" alt="Sensitivity heatmap">
    </div>
    <div class="card">
      <h2>Local pathway sensitivities</h2>
      <p class="note">Values below are centered finite-difference gradients at the reference state. They tell us which variable perturbs which heat-transfer channel most strongly, before any scenario logic is added.</p>
      <table>
        <thead><tr><th>Variable</th><th>Step</th><th>Unit</th><th>dPMV</th><th>dQ_conv</th><th>dQ_rad</th><th>dE_sk</th><th>dQ_res</th></tr></thead>
        <tbody>{sens_rows_html}</tbody>
      </table>
    </div>
  </div>

  <h2>Scenario baselines</h2>
  <table>
    <thead><tr><th>Scenario</th><th>Ta</th><th>Tr</th><th>RH</th><th>v</th><th>PMV</th><th>PPD</th><th>Largest baseline pathway</th></tr></thead>
    <tbody>{baseline_rows_html}</tbody>
  </table>

  <h2>Single-lever recovery tests</h2>
  <p class="note">These are deliberately minimal. Each row asks whether one control variable alone can return the state to <span class="mono">|PMV| ≤ 0.5</span> within bounded search ranges, and which heat-transfer pathway changes the most in the successful recovery.</p>
  {' '.join(scenario_sections)}

  <h2>Working interpretation</h2>
  <div class="card">
    <ul>
      <li>This prototype is useful if the paper is really about <strong>pathway ownership</strong> rather than about a larger scenario catalog.</li>
      <li>PMV is acting here as an acceptability boundary, not as the main object of interpretation.</li>
      <li>The next decision is whether to keep this as a single-occupant mechanism paper or couple it back to demographic heterogeneity and spatial variation later.</li>
    </ul>
  </div>
</body>
</html>
"""

    (OUT / "report.html").write_text(html, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)

    occupant = build_occupant()
    reference = evaluate(REFERENCE_STATE, occupant)
    reference["label"] = REFERENCE_OCCUPANT["label"]

    sensitivity_rows: List[Dict[str, object]] = []
    units = {"Ta": "C", "Tr": "C", "RH": "%", "v": "m/s", "clo": "clo", "met": "met"}
    for variable, step in SENSITIVITY_STEPS.items():
        row: Dict[str, object] = {"variable": variable, "step": step, "unit": units[variable]}
        row.update(finite_difference(variable, step, REFERENCE_STATE, occupant))
        sensitivity_rows.append(row)

    baselines: List[Dict[str, object]] = []
    for scenario in SCENARIOS:
        baseline = evaluate(scenario.state, occupant)
        baseline.update(
            {
                "scenario_id": scenario.id,
                "label": scenario.label,
                "description": scenario.description,
                "dominant_loss": dominant_channel(
                    {
                        "dQ_conv": baseline["Q_conv"],
                        "dQ_rad": baseline["Q_rad"],
                        "dE_sk": baseline["E_sk"],
                        "dQ_res": baseline["Q_res"],
                    }
                ),
            }
        )
        baselines.append(baseline)

    recoveries = [recover_with_lever(scenario, lever, occupant) for scenario in SCENARIOS for lever in LEVERS]

    write_csv(OUT / "pathway_sensitivities.csv", sensitivity_rows)
    write_csv(OUT / "scenario_baselines.csv", baselines)
    write_csv(OUT / "scenario_recoveries.csv", recoveries)
    (OUT / "reference_summary.json").write_text(json.dumps(reference, indent=2), encoding="utf-8")

    save_reference_plot(reference)
    save_sensitivity_heatmap(sensitivity_rows)
    render_report(reference, sensitivity_rows, baselines, recoveries)

    print(f"Wrote report to {OUT / 'report.html'}")


if __name__ == "__main__":
    main()
