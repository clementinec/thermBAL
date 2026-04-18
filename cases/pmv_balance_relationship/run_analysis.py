#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

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
    state: ComfortState


@dataclass(frozen=True)
class Lever:
    id: str
    label: str
    unit: str
    values: Tuple[float, ...]


@dataclass(frozen=True)
class Profile:
    id: str
    label: str
    age: int
    sex: str
    height_cm: float
    weight_kg: float
    clo: float
    met: float


@dataclass(frozen=True)
class StepFamily:
    id: str
    label: str
    RH: float
    v: float
    mode: str
    radiant_offset_c: float = 0.0


REFERENCE_OCCUPANT = {
    "height_cm": 170.0,
    "weight_kg": 70.0,
    "clo": 0.5,
    "met": 0.9,
}

PROFILES: Tuple[Profile, ...] = (
    Profile(
        id="smaller_female",
        label="Smaller female, lower-met",
        age=30,
        sex="female",
        height_cm=160.0,
        weight_kg=55.0,
        clo=0.45,
        met=0.8,
    ),
    Profile(
        id="reference_office",
        label="Reference office",
        age=40,
        sex="male",
        height_cm=170.0,
        weight_kg=70.0,
        clo=0.5,
        met=0.9,
    ),
    Profile(
        id="larger_male",
        label="Larger male, formal",
        age=35,
        sex="male",
        height_cm=182.0,
        weight_kg=88.0,
        clo=0.7,
        met=1.0,
    ),
    Profile(
        id="older_layered",
        label="Older layered occupant",
        age=70,
        sex="female",
        height_cm=165.0,
        weight_kg=64.0,
        clo=0.7,
        met=0.8,
    ),
)

PMV_STEP_STARTS: Tuple[int, ...] = (-3, -2, -1, 1, 2, 3)

STEP_FAMILIES: Tuple[StepFamily, ...] = (
    StepFamily("neutral_still", "Neutral still bulk", 50.0, 0.10, "bulk"),
    StepFamily("humid_still", "Humid still bulk", 80.0, 0.05, "bulk"),
    StepFamily("moving_air", "Moving-air bulk", 50.0, 0.20, "bulk"),
    StepFamily("radiant_bias", "Radiant-biased bulk", 50.0, 0.10, "signed_offset", 4.0),
)

REFERENCE_STATE = ComfortState(Ta=25.5, Tr=25.5, RH=50.0, v=0.10)
PMV_TARGET = 0.5

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
}

SCENARIOS: Tuple[Scenario, ...] = (
    Scenario("severe_cold_bulk", "Severe cold bulk", ComfortState(Ta=19.0, Tr=19.0, RH=50.0, v=0.10)),
    Scenario("cold_radiant_sink", "Cold radiant sink", ComfortState(Ta=23.0, Tr=17.0, RH=50.0, v=0.10)),
    Scenario("cold_draft", "Cold high-airflow draft", ComfortState(Ta=20.0, Tr=20.0, RH=50.0, v=0.60)),
    Scenario("severe_hot_bulk", "Severe hot bulk", ComfortState(Ta=32.0, Tr=32.0, RH=50.0, v=0.10)),
    Scenario("hot_radiant_gain", "Hot radiant gain", ComfortState(Ta=29.0, Tr=35.0, RH=50.0, v=0.10)),
    Scenario("hot_humid_still", "Hot humid still air", ComfortState(Ta=31.0, Tr=31.0, RH=80.0, v=0.05)),
)

LEVERS: Tuple[Lever, ...] = (
    Lever("air_temp", "All-air (Ta only)", "C", tuple(round(14.0 + 0.1 * i, 6) for i in range(int((36.0 - 14.0) / 0.1) + 1))),
    Lever("radiant", "Radiant (Tr only)", "C", tuple(round(14.0 + 0.1 * i, 6) for i in range(int((36.0 - 14.0) / 0.1) + 1))),
)


def build_occupant(*, met: float | None = None, clo: float | None = None) -> Occupant:
    bsa = du_bois_bsa(REFERENCE_OCCUPANT["weight_kg"], REFERENCE_OCCUPANT["height_cm"])
    return Occupant(
        m_wm2=met_to_wm2(met if met is not None else REFERENCE_OCCUPANT["met"]),
        clo=clo if clo is not None else REFERENCE_OCCUPANT["clo"],
        posture="seated",
        body_area_m2=bsa,
    )


def build_profile_occupant(profile: Profile) -> Occupant:
    bsa = du_bois_bsa(profile.weight_kg, profile.height_cm)
    return Occupant(
        m_wm2=met_to_wm2(profile.met),
        clo=profile.clo,
        posture="seated",
        body_area_m2=bsa,
    )


def pmv_factor_for_occupant(occupant: Occupant) -> float:
    return 0.303 * math.exp(-0.036 * occupant.m_wm2) + 0.028


def evaluate(state: ComfortState, occupant: Occupant) -> Dict[str, float]:
    result = compute_pmv(state, occupant)
    sensible_generation = occupant.m_wm2 - occupant.wme
    thermal_load = sensible_generation - result.L
    area = occupant.body_area_m2 or 1.0
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
        "thermal_load": thermal_load,
        "Q_conv_w": result.q_conv * area,
        "Q_rad_w": result.q_rad * area,
        "E_sk_w": result.e_sk * area,
        "Q_res_w": result.q_res * area,
        "thermal_load_w": thermal_load * area,
        "pmv_factor": pmv_factor_for_occupant(occupant),
        "met": result.inputs["met_equiv"],
        "clo": result.inputs["clo"],
        "body_area_m2": area,
    }


def ppd_from_pmv(pmv: float) -> float:
    return 100.0 - 95.0 * math.exp(-0.03353 * (pmv ** 4) - 0.2179 * (pmv ** 2))


def actual_vapor_pressure_pa(state: ComfortState) -> float:
    return sat_vapor_pressure(state.Ta) * (state.RH / 100.0)


def humidity_ratio_from_state(state: ComfortState) -> float:
    return humidity_ratio(actual_vapor_pressure_pa(state))


def apply_lever(scenario: Scenario, lever: Lever, value: float) -> ComfortState:
    if lever.id == "air_temp":
        return ComfortState(Ta=value, Tr=scenario.state.Tr, RH=scenario.state.RH, v=scenario.state.v)
    if lever.id == "radiant":
        return ComfortState(Ta=scenario.state.Ta, Tr=value, RH=scenario.state.RH, v=scenario.state.v)
    raise ValueError(f"Unknown lever: {lever.id}")


def step_target_label(start_label: int) -> int:
    if start_label < 0:
        return start_label + 1
    return start_label - 1


def format_step_label(start_label: int, target_label: int) -> str:
    def fmt(value: int) -> str:
        if value > 0:
            return f"+{value}"
        return str(value)

    return f"{fmt(start_label)}→{fmt(target_label)}"


def build_step_family_state(family: StepFamily, ta: float, start_label: int) -> ComfortState:
    if family.mode == "bulk":
        tr = ta
    elif family.mode == "signed_offset":
        sign = 1.0 if start_label > 0 else -1.0
        tr = ta + sign * family.radiant_offset_c
    else:
        raise ValueError(f"Unknown step family mode: {family.mode}")
    return ComfortState(Ta=ta, Tr=tr, RH=family.RH, v=family.v)


def solve_scalar_target(
    evaluator,
    target: float,
    low: float,
    high: float,
    samples: int = 500,
) -> float | None:
    xs = np.linspace(low, high, samples)
    ys = [float(evaluator(x)) - target for x in xs]

    best_idx = int(np.argmin(np.abs(ys)))
    if abs(ys[best_idx]) <= 1e-4:
        return float(xs[best_idx])

    for i in range(len(xs) - 1):
        y1 = ys[i]
        y2 = ys[i + 1]
        if y1 == 0.0:
            return float(xs[i])
        if y1 * y2 > 0:
            continue
        a = float(xs[i])
        b = float(xs[i + 1])
        fa = y1
        fb = y2
        for _ in range(50):
            m = 0.5 * (a + b)
            fm = float(evaluator(m)) - target
            if abs(fm) <= 1e-6:
                return m
            if fa * fm <= 0:
                b = m
                fb = fm
            else:
                a = m
                fa = fm
        return 0.5 * (a + b)
    return None


def solve_family_baseline_state(profile: Profile, family: StepFamily, start_label: int) -> ComfortState | None:
    occupant = build_profile_occupant(profile)

    def evaluator(ta: float) -> float:
        return evaluate(build_step_family_state(family, ta, start_label), occupant)["PMV"]

    ta = solve_scalar_target(evaluator, float(start_label), 5.0, 40.0)
    if ta is None:
        return None
    return build_step_family_state(family, ta, start_label)


def solve_air_shift_to_target(profile: Profile, baseline_state: ComfortState, target_pmv: float) -> ComfortState | None:
    occupant = build_profile_occupant(profile)

    def evaluator(ta: float) -> float:
        trial = ComfortState(Ta=ta, Tr=baseline_state.Tr, RH=baseline_state.RH, v=baseline_state.v)
        return evaluate(trial, occupant)["PMV"]

    ta = solve_scalar_target(evaluator, target_pmv, 5.0, 40.0)
    if ta is None:
        return None
    return ComfortState(Ta=ta, Tr=baseline_state.Tr, RH=baseline_state.RH, v=baseline_state.v)


def estimate_system_exergy(lever: Lever, baseline: ComfortState, recovered: ComfortState) -> Tuple[str, float, float]:
    if lever.id == "air_temp":
        rho = SYSTEM_ASSUMPTIONS["rho_air"]
        cp = SYSTEM_ASSUMPTIONS["cp_air"]
        mdot = rho * SYSTEM_ASSUMPTIONS["supply_air_flow_m3s"]
        thermal_w = mdot * cp * abs(recovered.Ta - baseline.Ta)
        if recovered.Ta >= baseline.Ta:
            return "air heating", thermal_w, thermal_w / SYSTEM_ASSUMPTIONS["cop_air_heating"]
        return "air cooling", thermal_w, thermal_w / SYSTEM_ASSUMPTIONS["cop_air_cooling"]

    if lever.id == "radiant":
        thermal_w = (
            SYSTEM_ASSUMPTIONS["h_r_panel_w_m2k"]
            * SYSTEM_ASSUMPTIONS["radiant_active_area_m2"]
            * abs(recovered.Tr - baseline.Tr)
        )
        if recovered.Tr >= baseline.Tr:
            return "radiant heating", thermal_w, thermal_w / SYSTEM_ASSUMPTIONS["cop_radiant_heating"]
        return "radiant cooling", thermal_w, thermal_w / SYSTEM_ASSUMPTIONS["cop_radiant_cooling"]

    raise ValueError(f"Unknown lever: {lever.id}")


def dominant_pathway_delta(before: Dict[str, float], after: Dict[str, float]) -> str:
    deltas = {
        "convection": abs(after["Q_conv"] - before["Q_conv"]),
        "radiation": abs(after["Q_rad"] - before["Q_rad"]),
        "evaporation": abs(after["E_sk"] - before["E_sk"]),
        "respiration": abs(after["Q_res"] - before["Q_res"]),
    }
    return max(deltas, key=deltas.get)


def find_recovery(scenario: Scenario, lever: Lever, occupant: Occupant) -> Dict[str, object] | None:
    baseline = evaluate(scenario.state, occupant)
    low = min(lever.values)
    high = max(lever.values)
    span = max(1e-9, high - low)
    base_value = scenario.state.Ta if lever.id == "air_temp" else scenario.state.Tr

    feasible: List[Tuple[float, ComfortState, Dict[str, float]]] = []
    for value in lever.values:
        candidate_state = apply_lever(scenario, lever, value)
        candidate = evaluate(candidate_state, occupant)
        if abs(candidate["PMV"]) <= PMV_TARGET:
            effort = abs(value - base_value) / span
            feasible.append((effort, candidate_state, candidate))

    if not feasible:
        return None

    effort, candidate_state, candidate = min(feasible, key=lambda item: item[0])
    mode, thermal_w, exergy_w = estimate_system_exergy(lever, scenario.state, candidate_state)
    baseline = evaluate(scenario.state, occupant)
    delta_pmv = candidate["PMV"] - baseline["PMV"]
    delta_load = candidate["thermal_load"] - baseline["thermal_load"]
    return {
        "scenario_id": scenario.id,
        "scenario_label": scenario.label,
        "lever_id": lever.id,
        "lever_label": lever.label,
        "baseline_pmv": baseline["PMV"],
        "recovered_pmv": candidate["PMV"],
        "baseline_ppd": baseline["PPD"],
        "recovered_ppd": candidate["PPD"],
        "baseline_load": baseline["thermal_load"],
        "recovered_load": candidate["thermal_load"],
        "delta_pmv": delta_pmv,
        "delta_load": delta_load,
        "control_move": abs((candidate_state.Ta if lever.id == "air_temp" else candidate_state.Tr) - base_value),
        "normalized_control_move": effort,
        "system_mode": mode,
        "system_thermal_w": thermal_w,
        "system_exergy_w": exergy_w,
        "pmv_per_100w_exergy": abs(delta_pmv) / exergy_w * 100.0 if exergy_w > 0 else float("nan"),
        "load_per_100w_exergy": abs(delta_load) / exergy_w * 100.0 if exergy_w > 0 else float("nan"),
        "dominant_pathway": dominant_pathway_delta(baseline, candidate),
    }


def write_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_pmv_load_samples() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for met in (0.8, 0.9, 1.0, 1.2):
        occupant = build_occupant(met=met)
        for Ta in (18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0):
            for Tr in (18.0, 22.0, 26.0, 30.0, 34.0):
                for RH in (30.0, 50.0, 70.0):
                    for v in (0.05, 0.10, 0.20):
                        result = evaluate(ComfortState(Ta=Ta, Tr=Tr, RH=RH, v=v), occupant)
                        rows.append(
                            {
                                "met": met,
                                "Ta": Ta,
                                "Tr": Tr,
                                "RH": RH,
                                "v": v,
                                "thermal_load": result["thermal_load"],
                                "PMV": result["PMV"],
                                "PPD": result["PPD"],
                                "pmv_factor": result["pmv_factor"],
                            }
                        )
    return rows


def save_pmv_vs_load_plot(rows: List[Dict[str, object]]) -> Path:
    fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=180)
    colors = {0.8: "#4C78A8", 0.9: "#72B7B2", 1.0: "#F58518", 1.2: "#54A24B"}
    for met, color in colors.items():
        sub = [row for row in rows if float(row["met"]) == met]
        x = np.array([float(row["thermal_load"]) for row in sub])
        y = np.array([float(row["PMV"]) for row in sub])
        ax.scatter(x, y, s=12, alpha=0.55, color=color, label=f"{met:.1f} met samples")
        k = float(sub[0]["pmv_factor"])
        xx = np.linspace(x.min(), x.max(), 200)
        ax.plot(xx, k * xx, color=color, linewidth=1.8, linestyle="--", label=f"{met:.1f} met theory")

    ax.axhline(0.0, color="#2B2B2B", linewidth=0.9)
    ax.axvline(0.0, color="#2B2B2B", linewidth=0.9)
    ax.set_xlabel("Thermal load H - L (W/m²)")
    ax.set_ylabel("PMV")
    ax.set_title("PMV Is Linear in Thermal Load for a Fixed Metabolic Rate")
    ax.grid(color="#D8D8D8", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")
    fig.tight_layout()
    out = FIGS / "pmv_vs_load.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def build_ppd_curve() -> List[Dict[str, object]]:
    rows = []
    for i in range(601):
        pmv = -3.0 + 6.0 * i / 600.0
        rows.append({"PMV": pmv, "PPD": ppd_from_pmv(pmv)})
    return rows


def save_ppd_curve_plot(rows: List[Dict[str, object]]) -> Path:
    fig, ax = plt.subplots(figsize=(6.8, 4.4), dpi=180)
    x = np.array([float(row["PMV"]) for row in rows])
    y = np.array([float(row["PPD"]) for row in rows])
    ax.plot(x, y, color="#4C78A8", linewidth=2.2)
    ax.axvspan(-0.5, 0.5, color="#E8F2F8", alpha=0.85, zorder=0)
    key_points = [(-1.0, ppd_from_pmv(-1.0)), (-0.5, ppd_from_pmv(-0.5)), (0.0, ppd_from_pmv(0.0)), (0.5, ppd_from_pmv(0.5)), (1.0, ppd_from_pmv(1.0))]
    ax.scatter([point[0] for point in key_points], [point[1] for point in key_points], color="#D62728", s=22, zorder=3)
    ax.annotate("minimum = 5% at PMV = 0", xy=(0.0, 5.0), xytext=(0.62, 17.0), arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8)
    ax.annotate("about 10% at |PMV| = 0.5", xy=(0.5, ppd_from_pmv(0.5)), xytext=(1.15, 28.0), arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8)
    ax.annotate("about 26% at |PMV| = 1.0", xy=(1.0, ppd_from_pmv(1.0)), xytext=(1.55, 43.0), arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8)
    ax.set_xlabel("PMV")
    ax.set_ylabel("PPD (%)")
    ax.set_title("PPD Restores Nonlinear Acceptability Meaning, but Not New Physics")
    ax.grid(color="#D8D8D8", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    out = FIGS / "ppd_curve.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def build_zero_pmv_grid() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    conditions = [
        {"id": "dry_still", "label": "RH 30%, v 0.05 m/s", "RH": 30.0, "v": 0.05},
        {"id": "neutral_still", "label": "RH 50%, v 0.10 m/s", "RH": 50.0, "v": 0.10},
        {"id": "humid_still", "label": "RH 80%, v 0.05 m/s", "RH": 80.0, "v": 0.05},
        {"id": "neutral_moving", "label": "RH 50%, v 0.20 m/s", "RH": 50.0, "v": 0.20},
    ]
    rows: List[Dict[str, object]] = []
    for condition in conditions:
        for Ta in np.linspace(18.0, 32.0, 45):
            for Tr in np.linspace(18.0, 36.0, 45):
                pmv = evaluate(ComfortState(Ta=float(Ta), Tr=float(Tr), RH=condition["RH"], v=condition["v"]), build_occupant())["PMV"]
                rows.append(
                    {
                        "condition_id": condition["id"],
                        "condition_label": condition["label"],
                        "Ta": float(Ta),
                        "Tr": float(Tr),
                        "RH": condition["RH"],
                        "v": condition["v"],
                        "PMV": pmv,
                    }
                )
    return rows, conditions


def save_zero_pmv_plot(rows: List[Dict[str, object]], conditions: List[Dict[str, object]]) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.4), dpi=180, sharex=True, sharey=True)
    levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    cmap = plt.get_cmap("RdYlBu_r")
    for ax, condition in zip(axes.flat, conditions):
        sub = [row for row in rows if row["condition_id"] == condition["id"]]
        ta_values = sorted({float(row["Ta"]) for row in sub})
        tr_values = sorted({float(row["Tr"]) for row in sub})
        z = np.array([[next(float(r["PMV"]) for r in sub if float(r["Ta"]) == ta and float(r["Tr"]) == tr) for ta in ta_values] for tr in tr_values])
        ax.contourf(ta_values, tr_values, z, levels=np.linspace(-2.5, 2.5, 21), cmap=cmap, vmin=-2.5, vmax=2.5)
        cs = ax.contour(ta_values, tr_values, z, levels=levels, colors=["#666666", "#999999", "#111111", "#999999", "#666666"], linewidths=[0.8, 0.8, 1.8, 0.8, 0.8])
        ax.clabel(cs, fmt={-1.0: "-1", -0.5: "-0.5", 0.0: "0", 0.5: "0.5", 1.0: "1"}, fontsize=7)
        ax.set_title(condition["label"], fontsize=9)
        ax.set_xlabel("Ta (°C)")
        ax.set_ylabel("Tr (°C)")
    fig.suptitle("Zero-PMV Manifolds in Ta-Tr Space", y=0.98, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIGS / "zero_pmv_manifold.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def build_air_radiant_effectiveness() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    occupant = build_occupant()
    for scenario in SCENARIOS:
        for lever in LEVERS:
            row = find_recovery(scenario, lever, occupant)
            if row is not None:
                rows.append(row)
    return rows


def build_profile_response_rows() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for profile in PROFILES:
        occupant = build_profile_occupant(profile)
        for scenario in SCENARIOS:
            baseline = evaluate(scenario.state, occupant)
            if abs(baseline["PMV"]) < 0.5:
                continue
            directed_step = 1.0 if baseline["PMV"] < 0.0 else -1.0
            for lever_id in ("air_temp", "radiant"):
                if lever_id == "air_temp":
                    shifted = ComfortState(
                        Ta=scenario.state.Ta + directed_step,
                        Tr=scenario.state.Tr,
                        RH=scenario.state.RH,
                        v=scenario.state.v,
                    )
                    lever_label = "Directed 1°C air shift"
                else:
                    shifted = ComfortState(
                        Ta=scenario.state.Ta,
                        Tr=scenario.state.Tr + directed_step,
                        RH=scenario.state.RH,
                        v=scenario.state.v,
                    )
                    lever_label = "Directed 1°C radiant shift"
                after = evaluate(shifted, occupant)
                rows.append(
                    {
                        "profile_id": profile.id,
                        "profile_label": profile.label,
                        "age": profile.age,
                        "sex": profile.sex,
                        "height_cm": profile.height_cm,
                        "weight_kg": profile.weight_kg,
                        "body_area_m2": baseline["body_area_m2"],
                        "clo": profile.clo,
                        "met": profile.met,
                        "scenario_id": scenario.id,
                        "scenario_label": scenario.label,
                        "lever_id": lever_id,
                        "lever_label": lever_label,
                        "baseline_pmv": baseline["PMV"],
                        "after_pmv": after["PMV"],
                        "pmv_change": after["PMV"] - baseline["PMV"],
                        "abs_pmv_improvement": abs(baseline["PMV"]) - abs(after["PMV"]),
                        "baseline_ppd": baseline["PPD"],
                        "after_ppd": after["PPD"],
                        "ppd_improvement": baseline["PPD"] - after["PPD"],
                        "baseline_load": baseline["thermal_load"],
                        "after_load": after["thermal_load"],
                        "load_change": after["thermal_load"] - baseline["thermal_load"],
                        "baseline_load_w": baseline["thermal_load_w"],
                        "after_load_w": after["thermal_load_w"],
                        "load_change_w": after["thermal_load_w"] - baseline["thermal_load_w"],
                        "dQ_conv": after["Q_conv"] - baseline["Q_conv"],
                        "dQ_rad": after["Q_rad"] - baseline["Q_rad"],
                        "dE_sk": after["E_sk"] - baseline["E_sk"],
                        "dQ_res": after["Q_res"] - baseline["Q_res"],
                        "dQ_conv_w": after["Q_conv_w"] - baseline["Q_conv_w"],
                        "dQ_rad_w": after["Q_rad_w"] - baseline["Q_rad_w"],
                        "dE_sk_w": after["E_sk_w"] - baseline["E_sk_w"],
                        "dQ_res_w": after["Q_res_w"] - baseline["Q_res_w"],
                    }
                )
    return rows


def save_air_radiant_effectiveness_plot(rows: List[Dict[str, object]]) -> Path:
    scenario_order = [scenario.id for scenario in SCENARIOS]
    scenario_labels = [scenario.label for scenario in SCENARIOS]
    air_vals = []
    rad_vals = []
    air_load_vals = []
    rad_load_vals = []
    for sid in scenario_order:
        air = next((row for row in rows if row["scenario_id"] == sid and row["lever_id"] == "air_temp"), None)
        rad = next((row for row in rows if row["scenario_id"] == sid and row["lever_id"] == "radiant"), None)
        air_vals.append(float(air["pmv_per_100w_exergy"]) if air is not None else np.nan)
        rad_vals.append(float(rad["pmv_per_100w_exergy"]) if rad is not None else np.nan)
        air_load_vals.append(float(air["load_per_100w_exergy"]) if air is not None else np.nan)
        rad_load_vals.append(float(rad["load_per_100w_exergy"]) if rad is not None else np.nan)

    x = np.arange(len(scenario_order))
    width = 0.36
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.0), dpi=180, sharex=True)

    axes[0].bar(x - width / 2, air_vals, width=width, color="#4C78A8", label="All-air")
    axes[0].bar(x + width / 2, rad_vals, width=width, color="#F58518", label="Radiant")
    axes[0].set_ylabel("|ΔPMV| per 100 W exergy")
    axes[0].set_title("Air vs Radiant: PMV Recovery per Unit Estimated System Exergy")
    axes[0].grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    axes[0].set_axisbelow(True)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(x - width / 2, air_load_vals, width=width, color="#4C78A8", label="All-air")
    axes[1].bar(x + width / 2, rad_load_vals, width=width, color="#F58518", label="Radiant")
    axes[1].set_ylabel("|Δ(H-L)| per 100 W exergy (W/m²)")
    axes[1].set_title("Air vs Radiant: Thermal-Load Change per Unit Estimated System Exergy")
    axes[1].grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    axes[1].set_axisbelow(True)
    axes[1].set_xticks(x, scenario_labels, rotation=18, ha="right")
    axes[1].legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out = FIGS / "air_vs_radiant_effectiveness.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_profile_response_heatmaps(rows: List[Dict[str, object]]) -> Path:
    scenario_order = [scenario.id for scenario in SCENARIOS]
    scenario_labels = [scenario.label for scenario in SCENARIOS]
    profile_order = [profile.id for profile in PROFILES]
    profile_labels = [profile.label for profile in PROFILES]

    def matrix(metric: str, lever_id: str) -> np.ndarray:
        arr = np.zeros((len(profile_order), len(scenario_order)))
        for i, pid in enumerate(profile_order):
            for j, sid in enumerate(scenario_order):
                row = next(r for r in rows if r["profile_id"] == pid and r["scenario_id"] == sid and r["lever_id"] == lever_id)
                arr[i, j] = float(row[metric])
        return arr

    pmv_air = matrix("abs_pmv_improvement", "air_temp")
    pmv_rad = matrix("abs_pmv_improvement", "radiant")
    ppd_air = matrix("ppd_improvement", "air_temp")
    ppd_rad = matrix("ppd_improvement", "radiant")

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.8), dpi=180, constrained_layout=True)

    pmv_max = max(np.max(np.abs(pmv_air)), np.max(np.abs(pmv_rad)), 1e-9)
    ppd_max = max(np.max(np.abs(ppd_air)), np.max(np.abs(ppd_rad)), 1e-9)

    def draw(ax, data, title, cmap, vmax, fmt):
        im = ax.imshow(data, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(scenario_labels)), scenario_labels, rotation=18, ha="right")
        ax.set_yticks(range(len(profile_labels)), profile_labels)
        ax.set_title(title, fontsize=10)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, fmt.format(data[i, j]), ha="center", va="center", fontsize=7, color="#111111")
        return im

    im1 = draw(axes[0, 0], pmv_air, "Directed 1°C Air Shift\nImprovement in |PMV|", "RdYlBu", pmv_max, "{:.2f}")
    im2 = draw(axes[0, 1], pmv_rad, "Directed 1°C Radiant Shift\nImprovement in |PMV|", "RdYlBu", pmv_max, "{:.2f}")
    im3 = draw(axes[1, 0], ppd_air, "Directed 1°C Air Shift\nReduction in PPD (percentage points)", "RdYlBu", ppd_max, "{:.1f}")
    im4 = draw(axes[1, 1], ppd_rad, "Directed 1°C Radiant Shift\nReduction in PPD (percentage points)", "RdYlBu", ppd_max, "{:.1f}")

    cbar1 = fig.colorbar(im1, ax=axes[0, :], fraction=0.025, pad=0.02)
    cbar1.set_label("positive = closer to PMV 0")
    cbar2 = fig.colorbar(im3, ax=axes[1, :], fraction=0.025, pad=0.02)
    cbar2.set_label("positive = lower predicted dissatisfied fraction")

    out = FIGS / "profile_response_heatmaps.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_profile_load_heatmaps(rows: List[Dict[str, object]]) -> Path:
    scenario_order = [scenario.id for scenario in SCENARIOS]
    scenario_labels = [scenario.label for scenario in SCENARIOS]
    profile_order = [profile.id for profile in PROFILES]
    profile_labels = [profile.label for profile in PROFILES]

    def matrix(metric: str, lever_id: str) -> np.ndarray:
        arr = np.zeros((len(profile_order), len(scenario_order)))
        for i, pid in enumerate(profile_order):
            for j, sid in enumerate(scenario_order):
                row = next(r for r in rows if r["profile_id"] == pid and r["scenario_id"] == sid and r["lever_id"] == lever_id)
                arr[i, j] = float(row[metric])
        return arr

    load_air = matrix("load_change_w", "air_temp")
    load_rad = matrix("load_change_w", "radiant")
    vmax = max(np.max(np.abs(load_air)), np.max(np.abs(load_rad)), 1e-9)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), dpi=180, constrained_layout=True)

    for ax, data, title in zip(
        axes,
        (load_air, load_rad),
        ("Directed 1°C Air Shift\nChange in Thermal Load (W/person)", "Directed 1°C Radiant Shift\nChange in Thermal Load (W/person)"),
    ):
        im = ax.imshow(data, cmap="RdYlBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(scenario_labels)), scenario_labels, rotation=18, ha="right")
        ax.set_yticks(range(len(profile_labels)), profile_labels)
        ax.set_title(title, fontsize=10)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, f"{data[i, j]:.1f}", ha="center", va="center", fontsize=7, color="#111111")

    cbar = fig.colorbar(im, ax=axes[:], fraction=0.025, pad=0.02)
    cbar.set_label("positive = higher thermal load W/person")

    out = FIGS / "profile_load_heatmaps.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_profile_pathway_heatmaps(rows: List[Dict[str, object]]) -> Path:
    scenario_order = [scenario.id for scenario in SCENARIOS]
    scenario_labels = [scenario.label for scenario in SCENARIOS]
    profile_order = [profile.id for profile in PROFILES]
    profile_labels = [profile.label for profile in PROFILES]

    def matrix(metric: str, lever_id: str) -> np.ndarray:
        arr = np.zeros((len(profile_order), len(scenario_order)))
        for i, pid in enumerate(profile_order):
            for j, sid in enumerate(scenario_order):
                row = next(r for r in rows if r["profile_id"] == pid and r["scenario_id"] == sid and r["lever_id"] == lever_id)
                arr[i, j] = float(row[metric])
        return arr

    air_conv = matrix("dQ_conv_w", "air_temp")
    air_rad = matrix("dQ_rad_w", "air_temp")
    rad_conv = matrix("dQ_conv_w", "radiant")
    rad_rad = matrix("dQ_rad_w", "radiant")

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.8), dpi=180, constrained_layout=True)
    vmax = max(
        np.max(np.abs(air_conv)),
        np.max(np.abs(air_rad)),
        np.max(np.abs(rad_conv)),
        np.max(np.abs(rad_rad)),
        1e-9,
    )

    def draw(ax, data, title):
        im = ax.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(scenario_labels)), scenario_labels, rotation=18, ha="right")
        ax.set_yticks(range(len(profile_labels)), profile_labels)
        ax.set_title(title, fontsize=10)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, f"{data[i, j]:.1f}", ha="center", va="center", fontsize=7, color="#111111")
        return im

    im = draw(axes[0, 0], air_conv, "Directed 1°C Air Shift\nΔ Convective Heat Loss (W/person)")
    draw(axes[0, 1], air_rad, "Directed 1°C Air Shift\nΔ Radiative Heat Loss (W/person)")
    draw(axes[1, 0], rad_conv, "Directed 1°C Radiant Shift\nΔ Convective Heat Loss (W/person)")
    draw(axes[1, 1], rad_rad, "Directed 1°C Radiant Shift\nΔ Radiative Heat Loss (W/person)")

    cbar = fig.colorbar(im, ax=axes[:, :], fraction=0.025, pad=0.02)
    cbar.set_label("positive = higher body heat loss in that pathway")

    out = FIGS / "profile_pathway_heatmaps.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def build_pmv_step_rows() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for profile in PROFILES:
        occupant = build_profile_occupant(profile)
        k = pmv_factor_for_occupant(occupant)
        theoretical_step_load = 1.0 / k
        for family in STEP_FAMILIES:
            for start_label in PMV_STEP_STARTS:
                target_label = step_target_label(start_label)
                baseline_state = solve_family_baseline_state(profile, family, start_label)
                if baseline_state is None:
                    continue
                baseline = evaluate(baseline_state, occupant)
                after_state = solve_air_shift_to_target(profile, baseline_state, float(target_label))
                if after_state is None:
                    continue
                after = evaluate(after_state, occupant)
                rows.append(
                    {
                        "profile_id": profile.id,
                        "profile_label": profile.label,
                        "family_id": family.id,
                        "family_label": family.label,
                        "start_label": start_label,
                        "target_label": target_label,
                        "transition_label": format_step_label(start_label, target_label),
                        "baseline_Ta": baseline_state.Ta,
                        "baseline_Tr": baseline_state.Tr,
                        "baseline_RH": baseline_state.RH,
                        "baseline_v": baseline_state.v,
                        "baseline_pmv": baseline["PMV"],
                        "baseline_ppd": baseline["PPD"],
                        "after_Ta": after_state.Ta,
                        "after_Tr": after_state.Tr,
                        "after_pmv": after["PMV"],
                        "after_ppd": after["PPD"],
                        "delta_ta": after_state.Ta - baseline_state.Ta,
                        "abs_delta_ta": abs(after_state.Ta - baseline_state.Ta),
                        "load_change": after["thermal_load"] - baseline["thermal_load"],
                        "load_change_w": after["thermal_load_w"] - baseline["thermal_load_w"],
                        "ppd_improvement": baseline["PPD"] - after["PPD"],
                        "theoretical_step_load_wm2": math.copysign(theoretical_step_load, target_label - start_label),
                        "theoretical_step_load_w": math.copysign(theoretical_step_load * occupant.body_area_m2, target_label - start_label),
                        "dQ_conv_w": after["Q_conv_w"] - baseline["Q_conv_w"],
                        "dQ_rad_w": after["Q_rad_w"] - baseline["Q_rad_w"],
                        "dE_sk_w": after["E_sk_w"] - baseline["E_sk_w"],
                        "dQ_res_w": after["Q_res_w"] - baseline["Q_res_w"],
                    }
                )
    return rows


def save_pmv_step_heatmap(rows: List[Dict[str, object]], metric: str, title: str, cbar_label: str, fmt: str, cmap: str) -> Path:
    transition_order = [format_step_label(start, step_target_label(start)) for start in PMV_STEP_STARTS]
    family_order = [family.id for family in STEP_FAMILIES]
    family_labels = [family.label for family in STEP_FAMILIES]
    profile_order = [profile.id for profile in PROFILES]
    profile_labels = [profile.label for profile in PROFILES]

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.2), dpi=180, constrained_layout=True)
    vmax = max(abs(float(row[metric])) for row in rows) if rows else 1.0
    vmax = max(vmax, 1e-9)

    for ax, family_id, family_label in zip(axes.flat, family_order, family_labels):
        data = np.zeros((len(profile_order), len(transition_order)))
        for i, pid in enumerate(profile_order):
            for j, tid in enumerate(transition_order):
                row = next(r for r in rows if r["profile_id"] == pid and r["family_id"] == family_id and r["transition_label"] == tid)
                data[i, j] = float(row[metric])
        im = ax.imshow(data, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(transition_order)), transition_order, rotation=20, ha="right")
        ax.set_yticks(range(len(profile_labels)), profile_labels)
        ax.set_title(family_label, fontsize=10)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, fmt.format(data[i, j]), ha="center", va="center", fontsize=7, color="#111111")

    cbar = fig.colorbar(im, ax=axes[:, :], fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label)
    out_name = {
        "delta_ta": "pmv_step_required_air_shift.png",
        "load_change_w": "pmv_step_load_change.png",
        "ppd_improvement": "pmv_step_ppd_change.png",
    }[metric]
    out = FIGS / out_name
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def render_report(
    pmv_load_rows: List[Dict[str, object]],
    ppd_rows: List[Dict[str, object]],
    zero_rows: List[Dict[str, object]],
    effectiveness_rows: List[Dict[str, object]],
    profile_response_rows: List[Dict[str, object]],
    pmv_step_rows: List[Dict[str, object]],
) -> None:
    def html_escape(text: object) -> str:
        return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    style = """
    body { font-family: Georgia, 'Times New Roman', serif; margin: 32px auto; max-width: 1120px; color: #222; line-height: 1.45; }
    h1, h2, h3 { font-family: Helvetica, Arial, sans-serif; color: #1f2b33; }
    .lede { font-size: 18px; max-width: 920px; }
    .note { color: #5b6670; font-size: 14px; }
    .grid { display: grid; grid-template-columns: 1.08fr 1fr; gap: 24px; margin: 24px 0 32px; align-items: start; }
    .card { border: 1px solid #d8dde1; border-radius: 10px; padding: 18px; background: #fbfcfd; }
    table { border-collapse: collapse; width: 100%; margin: 10px 0 24px; font-size: 14px; }
    th, td { border-bottom: 1px solid #e1e6ea; padding: 8px 10px; text-align: left; vertical-align: top; }
    th { background: #f3f6f8; font-family: Helvetica, Arial, sans-serif; }
    img { width: 100%; border: 1px solid #d8dde1; border-radius: 10px; background: white; }
    .mono { font-family: 'SFMono-Regular', Consolas, monospace; font-size: 13px; }
    .section { margin-top: 34px; }
    ul { margin-top: 6px; }
    """

    pmv_factor_rows = []
    for met in (0.8, 0.9, 1.0, 1.2):
        row = next(r for r in pmv_load_rows if float(r["met"]) == met)
        pmv_factor_rows.append(f"<tr><td>{met:.1f} met</td><td>{float(row['pmv_factor']):.4f}</td></tr>")
    pmv_factor_rows_html = "\n".join(pmv_factor_rows)

    profile_table_rows = []
    for profile in PROFILES:
        area = du_bois_bsa(profile.weight_kg, profile.height_cm)
        total_w = met_to_wm2(profile.met) * area
        profile_table_rows.append(
            "<tr>"
            f"<td>{html_escape(profile.label)}</td>"
            f"<td>{profile.age}</td>"
            f"<td>{html_escape(profile.sex)}</td>"
            f"<td>{profile.height_cm:.0f}</td>"
            f"<td>{profile.weight_kg:.0f}</td>"
            f"<td>{area:.2f}</td>"
            f"<td>{profile.clo:.2f}</td>"
            f"<td>{profile.met:.1f}</td>"
            f"<td>{total_w:.1f}</td>"
            "</tr>"
        )
    profile_table_rows_html = "\n".join(profile_table_rows)

    profile_summary_rows = []
    for profile in PROFILES:
        air = [row for row in profile_response_rows if row["profile_id"] == profile.id and row["lever_id"] == "air_temp"]
        rad = [row for row in profile_response_rows if row["profile_id"] == profile.id and row["lever_id"] == "radiant"]
        profile_summary_rows.append(
            "<tr>"
            f"<td>{html_escape(profile.label)}</td>"
            f"<td>{np.mean([float(row['abs_pmv_improvement']) for row in air]):.2f}</td>"
            f"<td>{np.mean([float(row['abs_pmv_improvement']) for row in rad]):.2f}</td>"
            f"<td>{np.mean([float(row['ppd_improvement']) for row in air]):.1f}</td>"
            f"<td>{np.mean([float(row['ppd_improvement']) for row in rad]):.1f}</td>"
            f"<td>{np.mean([abs(float(row['dQ_conv_w'])) for row in air]):.1f}</td>"
            f"<td>{np.mean([abs(float(row['dQ_rad_w'])) for row in rad]):.1f}</td>"
            "</tr>"
        )
    profile_summary_rows_html = "\n".join(profile_summary_rows)

    pmv_step_summary_rows = []
    for profile in PROFILES:
        sub = [row for row in pmv_step_rows if row["profile_id"] == profile.id]
        pmv_step_summary_rows.append(
            "<tr>"
            f"<td>{html_escape(profile.label)}</td>"
            f"<td>{np.mean([abs(float(row['delta_ta'])) for row in sub]):.2f}</td>"
            f"<td>{np.mean([abs(float(row['load_change'])) for row in sub]):.2f}</td>"
            f"<td>{np.mean([abs(float(row['load_change_w'])) for row in sub]):.1f}</td>"
            f"<td>{np.mean([float(row['ppd_improvement']) for row in sub]):.1f}</td>"
            f"<td>{np.mean([abs(float(row['dQ_conv_w'])) for row in sub]):.1f}</td>"
            f"<td>{np.mean([abs(float(row['dQ_rad_w'])) for row in sub]):.1f}</td>"
            "</tr>"
        )
    pmv_step_summary_rows_html = "\n".join(pmv_step_summary_rows)

    effectiveness_table_rows = []
    for scenario in SCENARIOS:
        air = next(row for row in effectiveness_rows if row["scenario_id"] == scenario.id and row["lever_id"] == "air_temp")
        rad = next((row for row in effectiveness_rows if row["scenario_id"] == scenario.id and row["lever_id"] == "radiant"), None)
        exergy_winner = air if rad is None or float(air["system_exergy_w"]) <= float(rad["system_exergy_w"]) else rad
        rad_control = "—" if rad is None else f"{rad['control_move']:.1f} C"
        rad_exergy = "—" if rad is None else f"{rad['system_exergy_w']:.1f} W"
        rad_pmv_eff = "—" if rad is None else f"{rad['pmv_per_100w_exergy']:.3f}"
        rad_pathway = "—" if rad is None else html_escape(rad["dominant_pathway"])
        effectiveness_table_rows.append(
            "<tr>"
            f"<td>{html_escape(scenario.label)}</td>"
            f"<td>{air['control_move']:.1f} C</td>"
            f"<td>{air['system_exergy_w']:.1f} W</td>"
            f"<td>{air['pmv_per_100w_exergy']:.3f}</td>"
            f"<td>{html_escape(air['dominant_pathway'])}</td>"
            f"<td>{rad_control}</td>"
            f"<td>{rad_exergy}</td>"
            f"<td>{rad_pmv_eff}</td>"
            f"<td>{rad_pathway}</td>"
            f"<td>{html_escape(exergy_winner['lever_label'])}</td>"
            "</tr>"
        )
    effectiveness_table_rows_html = "\n".join(effectiveness_table_rows)

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>PMV Balance Relationship</title>
  <style>{style}</style>
</head>
<body>
  <h1>PMV, PPD, Thermal Load, and Pathway Effectiveness</h1>
  <p class="lede">This study separates four things that are often collapsed together: <strong>PMV</strong> as a deterministic mapping from thermal load, <strong>PPD</strong> as a deterministic mapping from PMV, <strong>pathway decomposition</strong> as the way the environment achieves that load, and <strong>system-side exergy</strong> as the cost of producing that environmental change. That distinction matters because PMV itself does not “know” whether a recovery came through air or radiation once the total thermal load has been changed.</p>

  <div class="card">
    <h2>What the Model Actually Says</h2>
    <p>In this codebase, PMV is computed exactly as <span class="mono">PMV = (0.303 * exp(-0.036 * M) + 0.028) * (H - L)</span>, where <span class="mono">H = M - W</span> is internal sensible generation and <span class="mono">L = Q_conv + Q_rad + E_sk + Q_res</span> is the total heat loss term. PPD is then computed as <span class="mono">100 - 95 * exp(-0.03353 * PMV^4 - 0.2179 * PMV^2)</span>. So for a fixed metabolic rate, PMV is linear in thermal load, and PPD is a deterministic transform of PMV. PPD is <strong>not</strong> model accuracy, confidence, or probability that a TSV prediction is correct; it is the model's predicted fraction dissatisfied under the PMV framework. The default office reference in this branch is set to <span class="mono">0.9 met</span>, and the comparison band is <span class="mono">0.8 / 0.9 / 1.0 / 1.2 met</span>, rather than quietly assuming <span class="mono">1.0 met</span> as the universal sedentary office default.</p>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/pmv_vs_load.png" alt="PMV versus thermal load">
    </div>
    <div class="card">
      <h2>PMV vs Thermal Load</h2>
      <p>The scatter at left is sampled directly from the comfort engine across many environmental states. The dashed lines are the exact theoretical relation for each metabolic rate. They lie on top of the samples because, for a fixed <span class="mono">M</span>, PMV is literally linear in <span class="mono">H - L</span>.</p>
      <table>
        <thead><tr><th>Occupant</th><th>PMV factor</th></tr></thead>
        <tbody>{pmv_factor_rows_html}</tbody>
      </table>
      <p class="note">Implication: once total thermal load is known, PMV mostly hides pathway identity. Pathways matter because they are how the environment changes <span class="mono">H - L</span>, not because PMV assigns separate meaning to convective and radiative watts.</p>
    </div>
  </div>

  <div class="grid section">
    <div class="card">
      <h2>PPD vs PMV</h2>
      <p>The PPD curve has a fixed minimum of 5% at PMV = 0. That minimum is built into the model form. This is why PPD should be read as a dissatisfied-fraction surrogate under the PMV framework, not as an empirical accuracy score for the PMV predictor.</p>
      <ul>
        <li>PMV = 0 does not imply 0% dissatisfaction.</li>
        <li>PPD grows symmetrically with |PMV|.</li>
        <li>The comfort rim at <span class="mono">|PMV| ≤ 0.5</span> still spans roughly 5% to 10% dissatisfied.</li>
        <li>PPD adds acceptability curvature, but it still does not add pathway or physiological meaning.</li>
      </ul>
    </div>
    <div>
      <img src="figures/ppd_curve.png" alt="PPD curve">
    </div>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/zero_pmv_manifold.png" alt="Zero PMV manifold">
    </div>
    <div class="card">
      <h2>Zero-PMV Manifolds</h2>
      <p>The contour panels show that PMV = 0 is not a single state. It is a manifold in environmental input space. Different combinations of air temperature, radiant temperature, humidity, and air speed can all settle at PMV = 0 for the same occupant.</p>
      <ul>
        <li>Raising air speed shifts the neutral contour.</li>
        <li>Humidity changes matter, but less than Ta and Tr in this office range.</li>
        <li>Ta and Tr trade against each other along the PMV = 0 contour, which is why air and radiant systems can both hit neutrality through different pathways.</li>
      </ul>
    </div>
  </div>

  <div class="section">
    <h2>One-PMV-Step Air Moves Across TSV Labels</h2>
    <div class="card">
      <p>This layer answers the more precise question behind TSV labels: if you ask the model to move exactly one PMV label toward neutrality, such as <span class="mono">-3→-2</span> or <span class="mono">-1→0</span>, is the required change the same? The answer splits in two. For a fixed profile, the required change in <span class="mono">H-L</span> is nearly identical by construction because PMV is linear in thermal load. But the required <span class="mono">ΔTa</span> is <strong>not</strong> identical, because the mapping from air temperature to thermal load depends on the starting environmental family and pathway balance.</p>
      <table>
        <thead><tr><th>Profile</th><th>Avg |ΔTa| for one PMV step (°C)</th><th>Avg |Δ(H-L)| (W/m²)</th><th>Avg |Δ(H-L)| (W/person)</th><th>Avg ΔPPD (points)</th><th>Avg |ΔQ_conv| (W/person)</th><th>Avg |ΔQ_rad| (W/person)</th></tr></thead>
        <tbody>{pmv_step_summary_rows_html}</tbody>
      </table>
      <p class="note">This is where PMV and PPD split cleanly. The model says one PMV step is a constant thermal-load step within a profile, but the same one-step transition does not carry constant air-temperature, pathway, or dissatisfaction meaning.</p>
    </div>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/pmv_step_required_air_shift.png" alt="Required air shift per PMV step">
    </div>
    <div class="card">
      <h2>Required ΔTa Depends on the Starting State</h2>
      <p>Each cell shows the air-temperature move needed to go one PMV label closer to zero under a given environmental family. Positive values mean warming; negative values mean cooling. This is where the nonlinearity lives.</p>
      <ul>
        <li><span class="mono">-3→-2</span> and <span class="mono">-1→0</span> are not the same required <span class="mono">ΔTa</span>.</li>
        <li>The required air move also changes under humid, moving-air, and radiant-biased conditions.</li>
        <li>So equal PMV steps do not imply equal control displacement.</li>
      </ul>
    </div>
  </div>

  <div class="grid section">
    <div class="card">
      <h2>But the Thermal-Load Step Stays Nearly Constant</h2>
      <p>These heatmaps show the signed change in total occupant thermal load in <span class="mono">W/person</span> for the same one-step PMV move. Within each profile the rows are almost uniform, because a one-unit PMV step corresponds to a fixed load step set by that profile's PMV factor. What changes with profile is the <strong>size</strong> of that fixed load step, and what changes with state is how the air system has to redistribute pathway watts to achieve it.</p>
      <ul>
        <li>Lower-met profiles need a smaller load correction per PMV unit.</li>
        <li>Larger bodies show larger <span class="mono">W/person</span> for the same per-area PMV step.</li>
        <li>The interesting variation is therefore in <span class="mono">ΔTa</span> and pathway shares, not in the PMV-to-load algebra itself.</li>
      </ul>
    </div>
    <div>
      <img src="figures/pmv_step_load_change.png" alt="Load change per PMV step">
    </div>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/pmv_step_ppd_change.png" alt="PPD change per PMV step">
    </div>
    <div class="card">
      <h2>PPD Gives Unequal Meaning to Equal PMV Steps</h2>
      <p>The same one-label PMV transition does not produce the same change in predicted dissatisfaction. That is the useful thing PPD adds near the comfort rim. A step from <span class="mono">-2→-1</span> or <span class="mono">+2→+1</span> can remove much more dissatisfaction than a step from <span class="mono">-1→0</span> or <span class="mono">+1→0</span>, even though all are one-unit PMV moves.</p>
      <ul>
        <li>PMV steps are equal in thermal-load space for a fixed profile.</li>
        <li>PPD changes are not equal across the same PMV steps.</li>
        <li>PPD therefore restores operational comfort meaning near neutrality, but still does not restore pathway or physiology.</li>
      </ul>
    </div>
  </div>

  <div class="section">
    <h2>Air vs Radiant Recovery Under the Same PMV Target</h2>
    <table>
      <thead>
        <tr>
          <th rowspan="2">Scenario</th>
          <th colspan="4">All-Air</th>
          <th colspan="4">Radiant</th>
          <th rowspan="2">Lower exergy winner</th>
        </tr>
        <tr>
          <th>Δ control</th><th>Exergy</th><th>|ΔPMV| / 100W</th><th>Dominant pathway</th>
          <th>Δ control</th><th>Exergy</th><th>|ΔPMV| / 100W</th><th>Dominant pathway</th>
        </tr>
      </thead>
      <tbody>{effectiveness_table_rows_html}</tbody>
    </table>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/air_vs_radiant_effectiveness.png" alt="Air versus radiant effectiveness">
    </div>
    <div class="card">
      <h2>What Changes When System Work Is Added</h2>
      <p>This comparison answers the practical question directly. If you compare air and radiant recovery only by PMV target, both can often hit the same recovered PMV. If you compare them by <strong>system exergy</strong>, the ranking changes.</p>
      <ul>
        <li>Per occupant thermal-load logic, both levers are just different routes to changing <span class="mono">H - L</span>.</li>
        <li>Per estimated system exergy, radiant recovery is usually more efficient in the radiative scenarios and in most of the hot/cold bulk scenarios here.</li>
        <li>So the statement “radiant is more effective” is only meaningful once you specify the denominator: per PMV target, per body-load watt, or per system exergy watt.</li>
      </ul>
    </div>
  </div>

  <div class="section">
    <h2>Profile-Sensitive 1°C Response</h2>
    <div class="card">
      <p>This layer asks the next question directly: for different occupant archetypes, how much does a <span class="mono">directed 1°C</span> change in <span class="mono">Ta</span> or <span class="mono">Tr</span> move PMV toward neutrality, and which pathway watts actually move? The four profiles below change not only <span class="mono">met</span> and <span class="mono">clo</span>, but also body size. In the standard PMV formulation, body size does <strong>not</strong> directly change PMV once everything is normalized by body area, but it <strong>does</strong> change the total wattage when those pathway terms are converted back into <span class="mono">W/person</span>. Age and sex are therefore only entering here through assumed body size, clothing, and metabolic-rate choices; they are not explicit terms in Fanger's PMV equation.</p>
      <table>
        <thead><tr><th>Profile</th><th>Age</th><th>Sex</th><th>Height (cm)</th><th>Weight (kg)</th><th>BSA (m²)</th><th>Clo</th><th>Met</th><th>Metabolic rate (W/person)</th></tr></thead>
        <tbody>{profile_table_rows_html}</tbody>
      </table>
    </div>
  </div>

  <div class="grid section">
    <div>
      <img src="figures/profile_response_heatmaps.png" alt="Profile response heatmaps">
    </div>
    <div class="card">
      <h2>Moving PMV Toward Zero</h2>
      <p>The heatmaps at left show the reduction in <span class="mono">|PMV|</span> from a <strong>directed</strong> 1°C change: colder-than-neutral baselines get a warming move; hotter-than-neutral baselines get a cooling move. The lower row translates the same perturbation into total occupant thermal-load watts.</p>
      <table>
        <thead><tr><th>Profile</th><th>Avg |PMV| improvement, air</th><th>Avg |PMV| improvement, radiant</th><th>Avg ΔPPD, air</th><th>Avg ΔPPD, radiant</th><th>Avg |ΔQ_conv|, air (W/person)</th><th>Avg |ΔQ_rad|, radiant (W/person)</th></tr></thead>
        <tbody>{profile_summary_rows_html}</tbody>
      </table>
      <p class="note">This is the point your original question was pushing toward: the same <span class="mono">1°C</span> control move does not have the same PMV or PPD leverage across profiles and discomfort origins, even before pathway watts are compared.</p>
    </div>
  </div>

  <div class="grid section">
    <div class="card">
      <h2>Total Thermal-Load Motion Still Matters</h2>
      <p>The PPD view is useful because it restores acceptability meaning, but it should not erase the body-load interpretation. These panels keep the total thermal-load movement visible in <span class="mono">W/person</span>. In the current cases, directed air shifts usually move more total load than directed radiant shifts for the same nominal <span class="mono">1°C</span> input, because the air move perturbs multiple terms at once.</p>
    </div>
    <div>
      <img src="figures/profile_load_heatmaps.png" alt="Profile load heatmaps">
    </div>
  </div>

  <div class="grid section">
    <div class="card">
      <h2>Which Watts Actually Move</h2>
      <p>The pathway heatmaps separate <span class="mono">ΔQ_conv</span> from <span class="mono">ΔQ_rad</span> in <span class="mono">W/person</span>. They make the mechanism visible instead of leaving it hidden behind a single PMV scalar.</p>
      <ul>
        <li>A directed <span class="mono">Ta</span> change primarily moves convective watts, with only secondary radiative spillover.</li>
        <li>A directed <span class="mono">Tr</span> change primarily moves radiative watts, while still nudging convection through the changed surface-air balance.</li>
        <li>Larger bodies do not automatically get different PMV for the same per-area load, but they do show larger total watt shifts in <span class="mono">W/person</span>.</li>
      </ul>
    </div>
    <div>
      <img src="figures/profile_pathway_heatmaps.png" alt="Profile pathway heatmaps">
    </div>
  </div>

  <div class="section">
    <h2>Working Interpretation</h2>
    <div class="card">
      <ul>
        <li>PMV is not a free-floating empirical vote regressor in this implementation. It is a thermal-load expression wrapped in Fanger's empirical scaling factor.</li>
        <li>PPD is not prediction confidence. It is a dissatisfaction mapping from PMV, with a built-in floor of 5%.</li>
        <li>PPD is useful for interpreting the comfort rim, but it is still not an edge-physiology model.</li>
        <li>PMV = 0 means the modeled thermal load has been balanced for the chosen occupant and conditions, not that the environment is unique or uniformly comfortable for everyone.</li>
        <li>Pathway decomposition is still useful because it explains how the environment reaches that load balance.</li>
        <li>Radiant vs air is therefore a system question layered on top of the PMV model: per body-load watt the comparison is similar, but per system exergy watt the levers can diverge sharply.</li>
        <li>Changing demographics in the current model works mostly by changing <span class="mono">met</span>, <span class="mono">clo</span>, and body area. That changes both the slope from thermal load to PMV and the total pathway watts required to move a person back toward neutrality.</li>
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
    (OUT / ".mplconfig").mkdir(parents=True, exist_ok=True)

    pmv_load_rows = build_pmv_load_samples()
    ppd_rows = build_ppd_curve()
    zero_rows, conditions = build_zero_pmv_grid()
    effectiveness_rows = build_air_radiant_effectiveness()
    profile_response_rows = build_profile_response_rows()
    pmv_step_rows = build_pmv_step_rows()

    write_csv(OUT / "pmv_load_samples.csv", pmv_load_rows)
    write_csv(OUT / "ppd_curve.csv", ppd_rows)
    write_csv(OUT / "zero_pmv_grid.csv", zero_rows)
    write_csv(OUT / "air_radiant_effectiveness.csv", effectiveness_rows)
    write_csv(OUT / "profile_response_1c.csv", profile_response_rows)
    write_csv(OUT / "pmv_step_air_response.csv", pmv_step_rows)
    (OUT / "zero_pmv_conditions.json").write_text(json.dumps(conditions, indent=2), encoding="utf-8")

    save_pmv_vs_load_plot(pmv_load_rows)
    save_ppd_curve_plot(ppd_rows)
    save_zero_pmv_plot(zero_rows, conditions)
    save_air_radiant_effectiveness_plot(effectiveness_rows)
    save_profile_response_heatmaps(profile_response_rows)
    save_profile_load_heatmaps(profile_response_rows)
    save_profile_pathway_heatmaps(profile_response_rows)
    save_pmv_step_heatmap(
        pmv_step_rows,
        "delta_ta",
        "Required Air Shift for One PMV Step Toward Neutrality",
        "positive = warmer air, negative = cooler air",
        "{:.2f}",
        "RdBu_r",
    )
    save_pmv_step_heatmap(
        pmv_step_rows,
        "load_change_w",
        "Total Occupant Load Change for One PMV Step Toward Neutrality",
        "positive = higher occupant thermal load W/person",
        "{:.1f}",
        "RdBu_r",
    )
    save_pmv_step_heatmap(
        pmv_step_rows,
        "ppd_improvement",
        "Change in Predicted Dissatisfied Fraction for One PMV Step",
        "positive = lower predicted dissatisfied fraction",
        "{:.1f}",
        "RdYlBu",
    )
    render_report(pmv_load_rows, ppd_rows, zero_rows, effectiveness_rows, profile_response_rows, pmv_step_rows)

    print(f"Wrote report to {OUT / 'report.html'}")


if __name__ == "__main__":
    main()
