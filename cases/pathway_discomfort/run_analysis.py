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

from comfort_engine import (
    ComfortState,
    Occupant,
    compute_pmv,
    du_bois_bsa,
    humidity_ratio,
    met_to_wm2,
    sat_vapor_pressure,
)


@dataclass(frozen=True)
class Scenario:
    id: str
    label: str
    description: str
    mechanism_class: str
    state: ComfortState


@dataclass(frozen=True)
class Lever:
    id: str
    label: str
    unit: str
    values: Tuple[float, ...]
    base_value: Callable[[Scenario], float]
    apply: Callable[[Scenario, float], ComfortState]


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
PMV_DISCOMFORT_LIMIT = 2.0
PMV_RECOVERY_TARGET = 0.5

SYSTEM_ASSUMPTIONS = {
    "rho_air": 1.2,
    "cp_air": 1006.0,
    "supply_air_flow_m3s": 0.05,
    "radiant_active_area_m2": 12.0,
    "h_r_panel_w_m2k": 4.7,
    "cop_air_heating": 3.0,
    "cop_air_cooling": 3.0,
    "cop_radiant_heating": 4.5,
    "cop_radiant_cooling": 5.5,
    "cop_dehumidification": 1.8,
    "cop_humidification": 1.0,
    "latent_heat_vaporization_j_kg": 2.45e6,
    "fan_power_max_w": 60.0,
    "fan_v_max_m_s": 1.2,
}

SCENARIOS: Tuple[Scenario, ...] = (
    Scenario(
        id="severe_cold_bulk",
        label="Severe cold bulk",
        description="Uniform cold room where both air and radiant temperatures are depressed.",
        mechanism_class="bulk cold",
        state=ComfortState(Ta=19.0, Tr=19.0, RH=50.0, v=0.10),
    ),
    Scenario(
        id="cold_radiant_sink",
        label="Cold radiant sink",
        description="Moderate air temperature with markedly colder surrounding surfaces.",
        mechanism_class="radiant cold",
        state=ComfortState(Ta=23.0, Tr=17.0, RH=50.0, v=0.10),
    ),
    Scenario(
        id="cold_draft",
        label="Cold high-airflow draft",
        description="Cool room where elevated air speed pushes the occupant into strong convective heat loss.",
        mechanism_class="convective cold",
        state=ComfortState(Ta=20.0, Tr=20.0, RH=50.0, v=0.60),
    ),
    Scenario(
        id="severe_hot_bulk",
        label="Severe hot bulk",
        description="Uniform hot room where both air and radiant temperatures are elevated.",
        mechanism_class="bulk hot",
        state=ComfortState(Ta=32.0, Tr=32.0, RH=50.0, v=0.10),
    ),
    Scenario(
        id="hot_radiant_gain",
        label="Hot radiant gain",
        description="Moderate air temperature with strongly overheated surrounding surfaces causing radiant heat gain.",
        mechanism_class="radiant hot",
        state=ComfortState(Ta=29.0, Tr=35.0, RH=50.0, v=0.10),
    ),
    Scenario(
        id="hot_humid_still",
        label="Hot humid still air",
        description="Hot-humid state with low air movement and suppressed evaporative relief.",
        mechanism_class="humid hot",
        state=ComfortState(Ta=31.0, Tr=31.0, RH=80.0, v=0.05),
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
        values=_frange(14.0, 36.0, 0.1),
        base_value=lambda scenario: scenario.state.Ta,
        apply=lambda scenario, value: ComfortState(Ta=value, Tr=scenario.state.Tr, RH=scenario.state.RH, v=scenario.state.v),
    ),
    Lever(
        id="radiant",
        label="Radiant (Tr only)",
        unit="C",
        values=_frange(14.0, 36.0, 0.1),
        base_value=lambda scenario: scenario.state.Tr,
        apply=lambda scenario, value: ComfortState(Ta=scenario.state.Ta, Tr=value, RH=scenario.state.RH, v=scenario.state.v),
    ),
    Lever(
        id="airflow",
        label="Airflow control (v only)",
        unit="m/s",
        values=_frange(0.0, 1.20, 0.01),
        base_value=lambda scenario: scenario.state.v,
        apply=lambda scenario, value: ComfortState(Ta=scenario.state.Ta, Tr=scenario.state.Tr, RH=scenario.state.RH, v=value),
    ),
    Lever(
        id="humidity",
        label="Latent control (RH only)",
        unit="%",
        values=_frange(20.0, 90.0, 1.0),
        base_value=lambda scenario: scenario.state.RH,
        apply=lambda scenario, value: ComfortState(Ta=scenario.state.Ta, Tr=scenario.state.Tr, RH=value, v=scenario.state.v),
    ),
)


def build_occupant() -> Occupant:
    bsa = du_bois_bsa(REFERENCE_OCCUPANT["weight_kg"], REFERENCE_OCCUPANT["height_cm"])
    return Occupant(
        m_wm2=met_to_wm2(REFERENCE_OCCUPANT["met"]),
        clo=REFERENCE_OCCUPANT["clo"],
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


def actual_vapor_pressure_pa(state: ComfortState) -> float:
    return sat_vapor_pressure(state.Ta) * (state.RH / 100.0)


def humidity_ratio_from_state(state: ComfortState) -> float:
    return humidity_ratio(actual_vapor_pressure_pa(state))


def dominant_channel(delta: Dict[str, float]) -> str:
    channel_map = {
        "convection": abs(delta["dQ_conv"]),
        "radiation": abs(delta["dQ_rad"]),
        "evaporation": abs(delta["dE_sk"]),
        "respiration": abs(delta["dQ_res"]),
    }
    return max(channel_map, key=channel_map.get)


def dominant_baseline_pathway(result: Dict[str, float]) -> str:
    channel_map = {
        "convection": abs(result["Q_conv"]),
        "radiation": abs(result["Q_rad"]),
        "evaporation": abs(result["E_sk"]),
        "respiration": abs(result["Q_res"]),
    }
    return max(channel_map, key=channel_map.get)


def signed_pathway_note(result: Dict[str, float]) -> str:
    values = {
        "Q_conv": result["Q_conv"],
        "Q_rad": result["Q_rad"],
        "E_sk": result["E_sk"],
        "Q_res": result["Q_res"],
    }
    key = max(values, key=lambda name: abs(values[name]))
    sign = "loss" if values[key] >= 0.0 else "gain"
    label = {
        "Q_conv": "convective",
        "Q_rad": "radiative",
        "E_sk": "evaporative",
        "Q_res": "respiratory",
    }[key]
    return f"{label} {sign}"


def recover_with_lever(scenario: Scenario, lever: Lever, occupant: Occupant) -> Dict[str, object]:
    baseline = evaluate(scenario.state, occupant)
    feasible: List[Tuple[float, Dict[str, float], float]] = []

    scenario_value = lever.base_value(scenario)
    lever_span = max(1e-9, max(lever.values) - min(lever.values))

    for value in lever.values:
        candidate_state = lever.apply(scenario, value)
        candidate = evaluate(candidate_state, occupant)
        if abs(candidate["PMV"]) <= PMV_RECOVERY_TARGET:
            effort = abs(value - scenario_value) / lever_span
            feasible.append((effort, candidate, value))

    if not feasible:
        return {
            "scenario_id": scenario.id,
            "scenario_label": scenario.label,
            "scenario_mechanism": scenario.mechanism_class,
            "lever_id": lever.id,
            "lever_label": lever.label,
            "lever_unit": lever.unit,
            "feasible": False,
            "baseline_value": scenario_value,
            "recovered_value": None,
            "delta_value": None,
            "normalized_effort": None,
            "control_move": None,
            "baseline_pmv": baseline["PMV"],
            "recovered_pmv": None,
            "baseline_ppd": baseline["PPD"],
            "recovered_ppd": None,
            "dominant_channel": None,
            "system_mode": None,
            "system_thermal_w": None,
            "system_exergy_w": None,
            "dQ_conv": None,
            "dQ_rad": None,
            "dE_sk": None,
            "dQ_res": None,
        }

    effort, recovered, chosen_value = min(feasible, key=lambda item: item[0])
    candidate_state = lever.apply(scenario, chosen_value)
    pathway_delta = {
        "dQ_conv": recovered["Q_conv"] - baseline["Q_conv"],
        "dQ_rad": recovered["Q_rad"] - baseline["Q_rad"],
        "dE_sk": recovered["E_sk"] - baseline["E_sk"],
        "dQ_res": recovered["Q_res"] - baseline["Q_res"],
    }
    system_cost = estimate_system_cost(lever, scenario.state, candidate_state)
    return {
        "scenario_id": scenario.id,
        "scenario_label": scenario.label,
        "scenario_mechanism": scenario.mechanism_class,
        "lever_id": lever.id,
        "lever_label": lever.label,
        "lever_unit": lever.unit,
        "feasible": True,
        "baseline_value": scenario_value,
        "recovered_value": chosen_value,
        "delta_value": chosen_value - scenario_value,
        "normalized_effort": effort,
        "control_move": abs(chosen_value - scenario_value),
        "baseline_pmv": baseline["PMV"],
        "recovered_pmv": recovered["PMV"],
        "baseline_ppd": baseline["PPD"],
        "recovered_ppd": recovered["PPD"],
        "dominant_channel": dominant_channel(pathway_delta),
        "system_mode": system_cost["mode"],
        "system_thermal_w": system_cost["thermal_w"],
        "system_exergy_w": system_cost["exergy_w"],
        **pathway_delta,
    }


def select_winning_recoveries(recoveries: List[Dict[str, object]]) -> List[Dict[str, object]]:
    winners: Dict[str, Dict[str, object]] = {}
    for row in recoveries:
        if not row["feasible"]:
            continue
        sid = str(row["scenario_id"])
        current = winners.get(sid)
        if current is None or float(row["normalized_effort"]) < float(current["normalized_effort"]):
            winners[sid] = row
    return list(winners.values())


def select_min_exergy_recoveries(recoveries: List[Dict[str, object]]) -> List[Dict[str, object]]:
    winners: Dict[str, Dict[str, object]] = {}
    for row in recoveries:
        if not row["feasible"]:
            continue
        sid = str(row["scenario_id"])
        current = winners.get(sid)
        if current is None or float(row["system_exergy_w"]) < float(current["system_exergy_w"]):
            winners[sid] = row
    return list(winners.values())


def estimate_system_cost(lever: Lever, baseline: ComfortState, recovered: ComfortState) -> Dict[str, float | str]:
    rho = SYSTEM_ASSUMPTIONS["rho_air"]
    cp = SYSTEM_ASSUMPTIONS["cp_air"]
    vdot = SYSTEM_ASSUMPTIONS["supply_air_flow_m3s"]
    mdot = rho * vdot

    if lever.id == "air_temp":
        d_t = recovered.Ta - baseline.Ta
        thermal_w = mdot * cp * abs(d_t)
        if d_t >= 0.0:
            return {
                "mode": "air heating",
                "thermal_w": thermal_w,
                "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_air_heating"],
            }
        return {
            "mode": "air cooling",
            "thermal_w": thermal_w,
            "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_air_cooling"],
        }

    if lever.id == "radiant":
        d_t = recovered.Tr - baseline.Tr
        thermal_w = SYSTEM_ASSUMPTIONS["h_r_panel_w_m2k"] * SYSTEM_ASSUMPTIONS["radiant_active_area_m2"] * abs(d_t)
        if d_t >= 0.0:
            return {
                "mode": "radiant heating",
                "thermal_w": thermal_w,
                "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_radiant_heating"],
            }
        return {
            "mode": "radiant cooling",
            "thermal_w": thermal_w,
            "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_radiant_cooling"],
        }

    if lever.id == "airflow":
        p_max = SYSTEM_ASSUMPTIONS["fan_power_max_w"]
        v_max = SYSTEM_ASSUMPTIONS["fan_v_max_m_s"]
        base_power = p_max * (max(baseline.v, 0.0) / v_max) ** 3
        new_power = p_max * (max(recovered.v, 0.0) / v_max) ** 3
        return {
            "mode": "fan power",
            "thermal_w": abs(new_power - base_power),
            "exergy_w": new_power - base_power,
        }

    if lever.id == "humidity":
        w0 = humidity_ratio_from_state(baseline)
        w1 = humidity_ratio_from_state(recovered)
        thermal_w = mdot * SYSTEM_ASSUMPTIONS["latent_heat_vaporization_j_kg"] * abs(w1 - w0)
        if w1 >= w0:
            return {
                "mode": "humidification",
                "thermal_w": thermal_w,
                "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_humidification"],
            }
        return {
            "mode": "dehumidification",
            "thermal_w": thermal_w,
            "exergy_w": thermal_w / SYSTEM_ASSUMPTIONS["cop_dehumidification"],
        }

    raise ValueError(f"Unknown lever id: {lever.id}")


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


def save_baseline_pathway_plot(baselines: List[Dict[str, object]]) -> Path:
    labels = [row["label"] for row in baselines]
    x = np.arange(len(labels))
    width = 0.18
    components = [
        ("Q_conv", "Convection", "#5B8FF9"),
        ("Q_rad", "Radiation", "#F6A04D"),
        ("E_sk", "Evaporation", "#74C0A3"),
        ("Q_res", "Respiration", "#9AA5B1"),
    ]

    fig, ax = plt.subplots(figsize=(10.6, 5.4), dpi=180)
    for idx, (field, title, color) in enumerate(components):
        values = [float(row[field]) for row in baselines]
        ax.bar(x + (idx - 1.5) * width, values, width=width, label=title, color=color, edgecolor="#27343C", linewidth=0.7)

    ax.axhline(0.0, color="#27343C", linewidth=0.9)
    ax.set_ylabel("W/m²")
    ax.set_title("Extreme Discomfort Baselines: Occupant Heat-Transfer Pathways")
    ax.set_xticks(x, labels, rotation=18, ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)

    for idx, row in enumerate(baselines):
        pmv = float(row["PMV"])
        top = max(float(row["Q_conv"]), float(row["Q_rad"]), float(row["E_sk"]), float(row["Q_res"]))
        bottom = min(float(row["Q_conv"]), float(row["Q_rad"]), float(row["E_sk"]), float(row["Q_res"]))
        y = top + 2.0 if pmv >= 0.0 else bottom - 3.5
        va = "bottom" if pmv >= 0.0 else "top"
        ax.text(idx, y, f"PMV {pmv:+.2f}", ha="center", va=va, fontsize=8, color="#1F1F1F")

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=4, frameon=False, fontsize=8)
    fig.tight_layout()
    out = FIGS / "extreme_baselines.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_recovery_matrix(recoveries: List[Dict[str, object]]) -> Path:
    scenario_order = [s.id for s in SCENARIOS]
    scenario_labels = [s.label for s in SCENARIOS]
    lever_order = [l.id for l in LEVERS]
    lever_labels = [l.label.replace(" (", "\n(") for l in LEVERS]

    color_map = {
        "convection": "#5B8FF9",
        "radiation": "#F6A04D",
        "evaporation": "#74C0A3",
        "respiration": "#9AA5B1",
        "infeasible": "#E3E6E8",
    }

    lookup = {(str(r["scenario_id"]), str(r["lever_id"])): r for r in recoveries}
    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=180)

    for i, sid in enumerate(scenario_order):
        for j, lid in enumerate(lever_order):
            row = lookup[(sid, lid)]
            channel = row["dominant_channel"] if row["feasible"] else "infeasible"
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor=color_map[channel], edgecolor="white", linewidth=1.5)
            ax.add_patch(rect)
            text = f"{float(row['normalized_effort']):.2f}" if row["feasible"] else "—"
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color="#1F1F1F")

    ax.set_xlim(-0.5, len(lever_order) - 0.5)
    ax.set_ylim(len(scenario_order) - 0.5, -0.5)
    ax.set_xticks(range(len(lever_labels)), lever_labels)
    ax.set_yticks(range(len(scenario_labels)), scenario_labels)
    ax.set_title("Severe Discomfort Recovery Matrix\ncell color = dominant recovery channel, text = normalized control move")
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[:].set_visible(False)

    legend_items = [
        ("convection", "Convection"),
        ("radiation", "Radiation"),
        ("evaporation", "Evaporation"),
        ("respiration", "Respiration"),
        ("infeasible", "Infeasible"),
    ]
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[k]) for k, _ in legend_items]
    labels = [label for _, label in legend_items]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout()
    out = FIGS / "recovery_matrix.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_system_cost_matrix(recoveries: List[Dict[str, object]]) -> Path:
    scenario_order = [s.id for s in SCENARIOS]
    scenario_labels = [s.label for s in SCENARIOS]
    lever_order = [l.id for l in LEVERS]
    lever_labels = [l.label.replace(" (", "\n(") for l in LEVERS]
    lookup = {(str(r["scenario_id"]), str(r["lever_id"])): r for r in recoveries}

    feasible_values = [
        float(row["system_exergy_w"])
        for row in recoveries
        if row["feasible"] and row["system_exergy_w"] is not None and float(row["system_exergy_w"]) >= 0.0
    ]
    vmax = max(feasible_values) if feasible_values else 1.0
    cmap = plt.get_cmap("YlGnBu")

    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=180)
    for i, sid in enumerate(scenario_order):
        for j, lid in enumerate(lever_order):
            row = lookup[(sid, lid)]
            if row["feasible"]:
                value = float(row["system_exergy_w"])
                if value < 0.0:
                    color = "#B9E3C6"
                else:
                    color = cmap(min(1.0, value / max(vmax, 1e-9)))
                text = f"{value:.0f}"
            else:
                color = "#E3E6E8"
                text = "—"
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor=color, edgecolor="white", linewidth=1.5)
            ax.add_patch(rect)
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color="#1F1F1F")

    ax.set_xlim(-0.5, len(lever_order) - 0.5)
    ax.set_ylim(len(scenario_order) - 0.5, -0.5)
    ax.set_xticks(range(len(lever_labels)), lever_labels)
    ax.set_yticks(range(len(scenario_labels)), scenario_labels)
    ax.set_title("Estimated System Exergy Matrix\ncell text = electric / exergy surrogate W")
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[:].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=vmax))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="estimated system exergy W")
    fig.tight_layout()
    out = FIGS / "system_cost_matrix.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def render_report(
    reference: Dict[str, float],
    baselines: List[Dict[str, object]],
    recoveries: List[Dict[str, object]],
    winners: List[Dict[str, object]],
    exergy_winners: List[Dict[str, object]],
) -> None:
    recovery_groups: Dict[str, List[Dict[str, object]]] = {}
    for row in recoveries:
        recovery_groups.setdefault(str(row["scenario_id"]), []).append(row)

    winner_counts: Dict[str, int] = {}
    for row in winners:
        label = str(row["lever_label"])
        winner_counts[label] = winner_counts.get(label, 0) + 1

    exergy_winner_counts: Dict[str, int] = {}
    for row in exergy_winners:
        label = str(row["lever_label"])
        exergy_winner_counts[label] = exergy_winner_counts.get(label, 0) + 1

    feasible_counts: Dict[str, int] = {}
    for row in recoveries:
        if row["feasible"]:
            label = str(row["lever_label"])
            feasible_counts[label] = feasible_counts.get(label, 0) + 1

    def html_escape(text: object) -> str:
        return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    style = """
    body { font-family: Georgia, 'Times New Roman', serif; margin: 32px auto; max-width: 1120px; color: #222; line-height: 1.45; }
    h1, h2, h3 { font-family: Helvetica, Arial, sans-serif; color: #1f2b33; }
    .lede { font-size: 18px; max-width: 930px; }
    .note { color: #5b6670; font-size: 14px; }
    .grid { display: grid; grid-template-columns: 1.15fr 1fr; gap: 24px; margin: 24px 0 32px; align-items: start; }
    .card { border: 1px solid #d8dde1; border-radius: 10px; padding: 18px; background: #fbfcfd; }
    table { border-collapse: collapse; width: 100%; margin: 10px 0 24px; font-size: 14px; }
    th, td { border-bottom: 1px solid #e1e6ea; padding: 8px 10px; text-align: left; vertical-align: top; }
    th { background: #f3f6f8; font-family: Helvetica, Arial, sans-serif; }
    img { width: 100%; border: 1px solid #d8dde1; border-radius: 10px; background: white; }
    .mono { font-family: 'SFMono-Regular', Consolas, monospace; font-size: 13px; }
    .ok { color: #1d7f44; font-weight: 600; }
    .no { color: #9a2c2c; font-weight: 600; }
    .section { margin-top: 34px; }
    ul { margin-top: 6px; }
    """

    baseline_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(row['label'])}</td>"
        f"<td>{html_escape(row['mechanism_class'])}</td>"
        f"<td>{row['Ta']:.1f}</td><td>{row['Tr']:.1f}</td><td>{row['RH']:.0f}</td><td>{row['v']:.2f}</td>"
        f"<td>{row['PMV']:.2f}</td><td>{row['PPD']:.1f}</td>"
        f"<td>{html_escape(row['dominant_loss'])}</td>"
        f"<td>{html_escape(row['dominant_note'])}</td>"
        "</tr>"
        for row in baselines
    )

    winner_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(row['scenario_label'])}</td>"
        f"<td>{html_escape(row['lever_label'])}</td>"
        f"<td>{row['normalized_effort']:.3f}</td>"
        f"<td>{html_escape(row['dominant_channel'])}</td>"
        f"<td>{row['baseline_pmv']:.2f} → {row['recovered_pmv']:.2f}</td>"
        f"<td>{row['baseline_ppd']:.1f}% → {row['recovered_ppd']:.1f}%</td>"
        "</tr>"
        for row in winners
    )

    exergy_winner_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(row['scenario_label'])}</td>"
        f"<td>{html_escape(row['lever_label'])}</td>"
        f"<td>{row['system_exergy_w']:.1f}</td>"
        f"<td>{html_escape(row['system_mode'])}</td>"
        f"<td>{row['normalized_effort']:.3f}</td>"
        f"<td>{row['baseline_pmv']:.2f} → {row['recovered_pmv']:.2f}</td>"
        "</tr>"
        for row in exergy_winners
    )

    winner_map = {str(row["scenario_id"]): row for row in winners}
    exergy_winner_map = {str(row["scenario_id"]): row for row in exergy_winners}
    scenario_sections = []
    for base in baselines:
        rows = recovery_groups.get(str(base["scenario_id"]), [])
        winner = winner_map.get(str(base["scenario_id"]))
        exergy_winner = exergy_winner_map.get(str(base["scenario_id"]))
        table_rows = []
        for row in rows:
            if row["feasible"]:
                delta_text = f"{row['delta_value']:+.2f} {row['lever_unit']}"
                pmv_text = f"{row['baseline_pmv']:.2f} → {row['recovered_pmv']:.2f}"
                ppd_text = f"{row['baseline_ppd']:.1f}% → {row['recovered_ppd']:.1f}%"
                channel_text = html_escape(row["dominant_channel"])
                effort_text = f"{row['normalized_effort']:.3f}"
                exergy_text = f"{row['system_exergy_w']:.1f} W"
                thermal_text = f"{row['system_thermal_w']:.1f} W"
            else:
                delta_text = "—"
                pmv_text = f"{row['baseline_pmv']:.2f} → —"
                ppd_text = f"{row['baseline_ppd']:.1f}% → —"
                channel_text = "—"
                effort_text = "—"
                exergy_text = "—"
                thermal_text = "—"
            table_rows.append(
                "<tr>"
                f"<td>{html_escape(row['lever_label'])}</td>"
                f"<td class=\"{'ok' if row['feasible'] else 'no'}\">{'yes' if row['feasible'] else 'no'}</td>"
                f"<td>{delta_text}</td>"
                f"<td>{effort_text}</td>"
                f"<td>{exergy_text}</td>"
                f"<td>{thermal_text}</td>"
                f"<td>{pmv_text}</td>"
                f"<td>{ppd_text}</td>"
                f"<td>{channel_text}</td>"
                "</tr>"
            )
        winner_text = ""
        if winner is not None:
            winner_text = (
                f"Shortest feasible path: <strong>{html_escape(winner['lever_label'])}</strong> "
                f"({winner['normalized_effort']:.3f} normalized control move, dominant channel {html_escape(winner['dominant_channel'])})."
            )
        exergy_text = ""
        if exergy_winner is not None:
            exergy_text = (
                f" Lowest estimated system exergy: <strong>{html_escape(exergy_winner['lever_label'])}</strong> "
                f"({exergy_winner['system_exergy_w']:.1f} W electric/exergy surrogate, mode {html_escape(exergy_winner['system_mode'])})."
            )
        scenario_sections.append(
            "<div class='card'>"
            f"<h3>{html_escape(base['label'])}</h3>"
            f"<p class='note'>{html_escape(base['description'])}</p>"
            f"<p><span class='mono'>Ta={base['Ta']:.1f} C, Tr={base['Tr']:.1f} C, RH={base['RH']:.0f}%, v={base['v']:.2f} m/s, PMV={base['PMV']:.2f}, PPD={base['PPD']:.1f}%</span></p>"
            f"<p class='note'>Mechanism class: <strong>{html_escape(base['mechanism_class'])}</strong>. Largest absolute pathway = <strong>{html_escape(base['dominant_note'])}</strong>. {winner_text}{exergy_text}</p>"
            "<table>"
            "<thead><tr><th>Lever</th><th>Feasible</th><th>Minimum control move</th><th>Normalized control move</th><th>Estimated exergy</th><th>Estimated thermal delivery</th><th>PMV recovery</th><th>PPD recovery</th><th>Dominant recovery channel</th></tr></thead>"
            "<tbody>"
            + "\n".join(table_rows)
            + "</tbody></table></div>"
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Pathway Discomfort</title>
  <style>{style}</style>
</head>
<body>
  <h1>Pathway Discomfort</h1>
  <p class="lede">This branch pushes the mechanism study into clearly severe states. Every baseline scenario is selected so the representative office occupant begins at <span class="mono">|PMV| ≥ {PMV_DISCOMFORT_LIMIT:.1f}</span>, which makes the discomfort signal unambiguous before any recovery logic is tested. The question is then practical rather than abstract: once the occupant is plainly too cold or too hot, which single environmental lever recovers comfort most directly, and through which pathway watts?</p>

  <div class="card">
    <h2>Method in One Paragraph</h2>
    <p>One representative seated office occupant is evaluated with the repo's existing ISO 7730 PMV engine. Six severe discomfort archetypes are constructed so that each baseline state sits at or beyond <span class="mono">|PMV| ≥ {PMV_DISCOMFORT_LIMIT:.1f}</span>: bulk cold, radiant cold, convective cold, bulk hot, radiant hot, and hot-humid still air. Each state is then tested against four bounded single-lever interventions: air temperature only, mean radiant temperature only, airflow only, and humidity only. Recovery is defined by <span class="mono">|PMV| ≤ {PMV_RECOVERY_TARGET:.1f}</span>. A second layer then estimates system-side actuation using a canonical office control-zone surrogate: 50 L/s supply air for the all-air lever, a 12 m² radiant panel with <span class="mono">h_r = 4.7 W/m²K</span> for the radiant lever, cubic fan power for airflow, and latent load from humidity-ratio change for the humidity lever. Reported exergy is the estimated electric/input work after applying simplified COP assumptions.</p>
  </div>

  <div class="card section">
    <h2>System-Side Surrogate Assumptions</h2>
    <table>
      <tbody>
        <tr><th>All-air sensible control</th><td><span class="mono">Q = ρ c_p V̇ |ΔTa|</span> with <span class="mono">V̇ = 0.05 m³/s</span>, heating/cooling COP = 3.0</td></tr>
        <tr><th>Radiant sensible control</th><td><span class="mono">Q = h_r A |ΔTr|</span> with <span class="mono">A = 12 m²</span>, <span class="mono">h_r = 4.7 W/m²K</span>, heating COP = 4.5, cooling COP = 5.5</td></tr>
        <tr><th>Airflow control</th><td>Fan power surrogate <span class="mono">P ∝ v³</span> with <span class="mono">60 W</span> at <span class="mono">1.2 m/s</span></td></tr>
        <tr><th>Humidity control</th><td><span class="mono">Q = ṁ h_fg |ΔW|</span> from humidity-ratio change, dehumidification COP = 1.8, humidification COP = 1.0</td></tr>
      </tbody>
    </table>
    <p class="note">This is intentionally a system surrogate, not a full HVAC plant model. It exists to keep numeric input displacement from being misread as true plant effort.</p>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/extreme_baselines.png" alt="Extreme discomfort baseline pathways">
    </div>
    <div class="card">
      <h2>Reference and Scenario Framing</h2>
      <p><span class="mono">Reference office = Ta={reference['Ta']:.1f} C, Tr={reference['Tr']:.1f} C, RH={reference['RH']:.0f}%, v={reference['v']:.2f} m/s, PMV={reference['PMV']:.2f}</span></p>
      <p class="note">The figure at left shows the actual heat-balance decomposition for each severe baseline. This matters because the states are not interchangeable. Some are dominated by elevated radiative loss, others by excessive convective loss, and hot radiant gain is the only case where the radiative term reverses sign and becomes net heat gain to the body.</p>
      <ul>
        <li><strong>Severe cold bulk</strong> is split almost evenly between convective and radiative loss.</li>
        <li><strong>Cold radiant sink</strong> is strongly radiation-led.</li>
        <li><strong>Cold high-airflow draft</strong> is strongly convection-led.</li>
        <li><strong>Hot radiant gain</strong> carries an explicit negative radiative term, meaning the environment is heating the body radiatively.</li>
      </ul>
    </div>
  </div>

  <div class="section">
    <h2>Severe baseline states</h2>
    <table>
      <thead><tr><th>Scenario</th><th>Mechanism</th><th>Ta</th><th>Tr</th><th>RH</th><th>v</th><th>PMV</th><th>PPD</th><th>Largest absolute pathway</th><th>Signed note</th></tr></thead>
      <tbody>{baseline_rows_html}</tbody>
    </table>
  </div>

  <div class="section">
    <h2>Shortest feasible recovery by scenario</h2>
    <table>
      <thead><tr><th>Scenario</th><th>Winning lever</th><th>Normalized control move</th><th>Dominant channel</th><th>PMV</th><th>PPD</th></tr></thead>
      <tbody>{winner_rows_html}</tbody>
    </table>
  </div>

  <div class="section">
    <h2>Lowest estimated system exergy by scenario</h2>
    <table>
      <thead><tr><th>Scenario</th><th>Winning lever</th><th>Estimated exergy</th><th>System mode</th><th>Normalized control move</th><th>PMV</th></tr></thead>
      <tbody>{exergy_winner_rows_html}</tbody>
    </table>
  </div>

  <div class="section">
    <h2>Recovery ownership matrix</h2>
    <div class="grid">
      <div>
        <img src="figures/recovery_matrix.png" alt="Recovery ownership matrix">
      </div>
      <div class="card">
        <p>The matrix condenses the single-lever recovery study. Every cell asks whether one control variable alone can drag the occupant from a clearly severe state back to <span class="mono">|PMV| ≤ {PMV_RECOVERY_TARGET:.1f}</span>. The cell color identifies the dominant pathway change in that recovery, while the overlaid number gives the smallest normalized control move found in the bounded search.</p>
        <ul>
          <li><strong>{winner_counts.get('All-air (Ta only)', 0)} of {len(winners)}</strong> scenarios are won by all-air recovery when the target is PMV-only comfort restoration.</li>
          <li>Radiant correction remains feasible in <strong>{feasible_counts.get('Radiant (Tr only)', 0)}</strong> of {len(SCENARIOS)} cases and is the consistent second path in every radiative archetype.</li>
          <li>Airflow and humidity alone do not recover any of these severe PMV states to <span class="mono">|PMV| ≤ {PMV_RECOVERY_TARGET:.1f}</span> within the bounded search.</li>
          <li>So the matrix is doing two jobs at once: it shows PMV's bias toward dry-bulb correction, and it still preserves the distinct pathway signatures behind each baseline state.</li>
        </ul>
      </div>
    </div>
  </div>

  <div class="section">
    <h2>Estimated system exergy matrix</h2>
    <div class="grid">
      <div>
        <img src="figures/system_cost_matrix.png" alt="Estimated system exergy matrix">
      </div>
      <div class="card">
        <p>This second matrix uses the same feasible recoveries, but ranks them by estimated system exergy rather than by normalized input displacement. That is the distinction the original report was missing: changing <span class="mono">Tr</span> by one degree is not assumed to have the same plant implication as changing <span class="mono">Ta</span> by one degree.</p>
        <ul>
          <li><strong>{exergy_winner_counts.get('Radiant (Tr only)', 0)} of {len(exergy_winners)}</strong> scenarios are exergy-optimal under the current surrogate assumptions.</li>
          <li>Radiant correction becomes the lower-work option in most feasible hot and cold cases, even where all-air still wins on normalized control displacement.</li>
          <li>The divergence between the two winner tables is the practical reason to separate occupant-side effectiveness from system-side cost.</li>
        </ul>
      </div>
    </div>
  </div>

  <div class="section">
    <h2>Detailed recovery tests</h2>
    <p class="note">Each scenario below keeps three things visible at once: the severe baseline state, which single-lever controls are actually feasible, and which pathway does the real work when recovery is possible.</p>
    {' '.join(scenario_sections)}
  </div>

  <div class="section">
    <h2>Working interpretation</h2>
    <div class="card">
      <ul>
        <li>Severe discomfort states separate much more cleanly than mild PMV deviations. Once the baseline crosses <span class="mono">|PMV| ≥ {PMV_DISCOMFORT_LIMIT:.1f}</span>, the heat-balance signature of each archetype is visually obvious in the decomposition plot.</li>
        <li>At the same time, the shortest PMV-only recovery collapses toward dry-bulb control: all-air wins every scenario in this branch, even when the baseline state is visibly radiative or convective in origin.</li>
        <li>Once the simple system surrogate is added, the story changes. Radiant correction often becomes the lower-exergy option even when all-air still wins on normalized input displacement.</li>
        <li>That tension is the main result. Occupant-side PMV recovery and system-side exergy cost are not the same ranking problem, and collapsing them into one scalar “effort” would be misleading.</li>
      </ul>
    </div>
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

    baselines: List[Dict[str, object]] = []
    for scenario in SCENARIOS:
        baseline = evaluate(scenario.state, occupant)
        baseline.update(
            {
                "scenario_id": scenario.id,
                "label": scenario.label,
                "description": scenario.description,
                "mechanism_class": scenario.mechanism_class,
                "dominant_loss": dominant_baseline_pathway(baseline),
                "dominant_note": signed_pathway_note(baseline),
            }
        )
        baselines.append(baseline)

    recoveries = [recover_with_lever(scenario, lever, occupant) for scenario in SCENARIOS for lever in LEVERS]
    winners = select_winning_recoveries(recoveries)
    exergy_winners = select_min_exergy_recoveries(recoveries)

    write_csv(OUT / "scenario_baselines.csv", baselines)
    write_csv(OUT / "scenario_recoveries.csv", recoveries)
    write_csv(OUT / "scenario_winners.csv", winners)
    write_csv(OUT / "scenario_exergy_winners.csv", exergy_winners)
    (OUT / "reference_summary.json").write_text(json.dumps(reference, indent=2), encoding="utf-8")

    save_baseline_pathway_plot(baselines)
    save_recovery_matrix(recoveries)
    save_system_cost_matrix(recoveries)
    render_report(reference, baselines, recoveries, winners, exergy_winners)

    print(f"Wrote report to {OUT / 'report.html'}")


if __name__ == "__main__":
    main()
