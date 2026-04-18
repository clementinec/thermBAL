#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
OUT = CASE / "out"
FIGS = OUT / "figures"
TABLES = OUT / "tables"
GENERATED = OUT / "generated"
MPLCONFIGDIR = OUT / ".mplconfig"
CACHE_DIR = OUT / ".cache"

PMV_CASE = ROOT / "cases" / "pmv_balance_relationship"
DISCOMFORT_CASE = ROOT / "cases" / "pathway_discomfort"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)
import os
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from comfort_engine import du_bois_bsa, met_to_wm2


def load_module(module_path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


PMV_MOD = load_module(PMV_CASE / "run_analysis.py", "pmv_balance_relationship_run_analysis")
DISCOMFORT_MOD = load_module(DISCOMFORT_CASE / "run_analysis.py", "pathway_discomfort_run_analysis")


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def latex_escape(text: Any) -> str:
    s = str(text)
    replacements = {
        "\\": "\\textbackslash{}",
        "&": "\\&",
        "%": "\\%",
        "$": "\\$",
        "#": "\\#",
        "_": "\\_",
        "{": "\\{",
        "}": "\\}",
        "~": "\\textasciitilde{}",
        "^": "\\textasciicircum{}",
        "→": "$\\rightarrow$",
        "°": "$^\\circ$",
        "≤": "$\\leq$",
        "≥": "$\\geq$",
        "|": "$|$",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    return s


def float_range(values: Iterable[float]) -> tuple[float, float]:
    values = list(values)
    return (min(values), max(values))


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else float("nan")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_table(path: Path, caption: str, label: str, col_spec: str, header: List[str], rows: List[List[str]], notes: str | None = None) -> None:
    body = "\n".join(" & ".join(row) + r" \\" for row in rows)
    note_block = f"\n\\begin{{flushleft}}\\footnotesize {notes}\\end{{flushleft}}" if notes else ""
    text = rf"""\begin{{table}}[t]
\centering
\caption{{{caption}}}
\label{{{label}}}
\begin{{adjustbox}}{{max width=\textwidth}}
\begin{{tabular}}{{{col_spec}}}
\toprule
{' & '.join(header)} \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}{note_block}
\end{{table}}
"""
    write_text(path, text)


def save_two_panel(left_path: Path, right_path: Path, left_title: str, right_title: str, out_path: Path, figsize=(13.5, 6.0)) -> None:
    fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=220, constrained_layout=True)
    for ax, path, title, label in zip(axes, (left_path, right_path), (left_title, right_title), ("a", "b")):
        image = mpimg.imread(path)
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(title, fontsize=12)
        ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", ha="left", fontsize=14, fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2})
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def save_panel_strip(panels: List[tuple[Path, str, str]], out_path: Path, figsize=(18.0, 5.8)) -> None:
    fig, axes = plt.subplots(1, len(panels), figsize=figsize, dpi=220, constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (path, title, label) in zip(axes, panels):
        image = mpimg.imread(path)
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(title, fontsize=12)
        ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", ha="left", fontsize=14, fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2})
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def build_figures() -> Dict[str, str]:
    figures = {}
    save_two_panel(
        PMV_CASE / "out" / "figures" / "pmv_vs_load.png",
        PMV_CASE / "out" / "figures" / "ppd_curve.png",
        "PMV as a linear function of thermal load",
        "PPD as a deterministic transform of PMV",
        FIGS / "fig01_pmv_foundations.png",
        figsize=(13.5, 5.8),
    )
    figures["pmv_foundations"] = "out/figures/fig01_pmv_foundations.png"

    shutil.copy2(PMV_CASE / "out" / "figures" / "zero_pmv_manifold.png", FIGS / "fig02_zero_manifold.png")
    figures["zero_manifold"] = "out/figures/fig02_zero_manifold.png"

    save_panel_strip(
        [
            (
                PMV_CASE / "out" / "figures" / "pmv_step_required_air_shift.png",
                "Required air-temperature move for a one-label PMV transition",
                "a",
            ),
            (
                PMV_CASE / "out" / "figures" / "pmv_step_load_change.png",
                "Occupant thermal-load change for the same PMV transition",
                "b",
            ),
            (
                PMV_CASE / "out" / "figures" / "pmv_step_ppd_change.png",
                "Predicted dissatisfied-fraction change for the same PMV transition",
                "c",
            ),
        ],
        FIGS / "fig03_step_response.png",
        figsize=(18.8, 5.8),
    )
    figures["step_response"] = "out/figures/fig03_step_response.png"

    save_two_panel(
        PMV_CASE / "out" / "figures" / "profile_response_heatmaps.png",
        PMV_CASE / "out" / "figures" / "profile_pathway_heatmaps.png",
        "Directed 1$^\\circ$C air/radiant PMV and PPD response across profiles",
        "Directed 1$^\\circ$C pathway watt shifts across profiles",
        FIGS / "fig04_profile_pathways.png",
        figsize=(13.8, 8.0),
    )
    figures["profile_pathways"] = "out/figures/fig04_profile_pathways.png"

    save_two_panel(
        DISCOMFORT_CASE / "out" / "figures" / "recovery_matrix.png",
        DISCOMFORT_CASE / "out" / "figures" / "system_cost_matrix.png",
        "Shortest PMV recovery by normalized control move",
        "Lowest estimated system exergy for the same severe states",
        FIGS / "fig05_recovery_exergy.png",
        figsize=(13.8, 6.5),
    )
    figures["recovery_exergy"] = "out/figures/fig05_recovery_exergy.png"
    return figures


def build_tables_and_macros() -> Dict[str, Any]:
    pmv_load_rows = read_csv(PMV_CASE / "out" / "pmv_load_samples.csv")
    step_rows = read_csv(PMV_CASE / "out" / "pmv_step_air_response.csv")
    profile_rows = read_csv(PMV_CASE / "out" / "profile_response_1c.csv")
    discomfort_baselines = read_csv(DISCOMFORT_CASE / "out" / "scenario_baselines.csv")
    discomfort_recoveries = read_csv(DISCOMFORT_CASE / "out" / "scenario_recoveries.csv")
    discomfort_winners = read_csv(DISCOMFORT_CASE / "out" / "scenario_winners.csv")
    discomfort_exergy = read_csv(DISCOMFORT_CASE / "out" / "scenario_exergy_winners.csv")

    profiles_table_rows: List[List[str]] = []
    pmv_step_summary_rows: List[List[str]] = []
    macros: Dict[str, str] = {}

    one_step_loads: List[float] = []
    air_improvements: List[float] = []
    rad_improvements: List[float] = []
    air_ppd_improvements: List[float] = []
    rad_ppd_improvements: List[float] = []
    step_ppd_changes: List[float] = []

    for profile in PMV_MOD.PROFILES:
        bsa = du_bois_bsa(profile.weight_kg, profile.height_cm)
        metabolic_w = met_to_wm2(profile.met) * bsa
        profiles_table_rows.append([
            latex_escape(profile.label),
            str(profile.age),
            latex_escape(profile.sex),
            f"{profile.height_cm:.0f}",
            f"{profile.weight_kg:.0f}",
            f"{bsa:.2f}",
            f"{profile.clo:.2f}",
            f"{profile.met:.1f}",
            f"{metabolic_w:.1f}",
        ])

        sub_steps = [row for row in step_rows if row["profile_id"] == profile.id]
        sub_air = [row for row in profile_rows if row["profile_id"] == profile.id and row["lever_id"] == "air_temp"]
        sub_rad = [row for row in profile_rows if row["profile_id"] == profile.id and row["lever_id"] == "radiant"]

        step_load = mean(abs(float(row["load_change_w"])) for row in sub_steps)
        step_std = math.sqrt(mean((abs(float(row["load_change_w"])) - step_load) ** 2 for row in sub_steps))
        delta_ta_min, delta_ta_max = float_range(abs(float(row["delta_ta"])) for row in sub_steps)
        step_ppd = mean(float(row["ppd_improvement"]) for row in sub_steps)
        avg_air = mean(float(row["abs_pmv_improvement"]) for row in sub_air)
        avg_rad = mean(float(row["abs_pmv_improvement"]) for row in sub_rad)
        avg_air_ppd = mean(float(row["ppd_improvement"]) for row in sub_air)
        avg_rad_ppd = mean(float(row["ppd_improvement"]) for row in sub_rad)

        one_step_loads.append(step_load)
        air_improvements.append(avg_air)
        rad_improvements.append(avg_rad)
        air_ppd_improvements.append(avg_air_ppd)
        rad_ppd_improvements.append(avg_rad_ppd)
        step_ppd_changes.extend(float(row["ppd_improvement"]) for row in sub_steps)

        prefix = {
            "smaller_female": "Smaller",
            "reference_office": "Reference",
            "larger_male": "Larger",
            "older_layered": "Older",
        }[profile.id]
        macros[f"{prefix}OneStepLoad"] = f"{step_load:.1f}"
        macros[f"{prefix}DeltaTaMin"] = f"{delta_ta_min:.2f}"
        macros[f"{prefix}DeltaTaMax"] = f"{delta_ta_max:.2f}"
        macros[f"{prefix}AirImprove"] = f"{avg_air:.2f}"
        macros[f"{prefix}RadImprove"] = f"{avg_rad:.2f}"
        macros[f"{prefix}StepPPD"] = f"{step_ppd:.1f}"
        macros[f"{prefix}AirPPDImprove"] = f"{avg_air_ppd:.1f}"
        macros[f"{prefix}RadPPDImprove"] = f"{avg_rad_ppd:.1f}"

    severe_scenario_rows: List[List[str]] = []
    for scenario in DISCOMFORT_MOD.SCENARIOS:
        baseline = next(row for row in discomfort_baselines if row["scenario_id"] == scenario.id)
        severe_scenario_rows.append([
            latex_escape(scenario.label),
            latex_escape(scenario.mechanism_class),
            f"{float(baseline['Ta']):.1f}",
            f"{float(baseline['Tr']):.1f}",
            f"{float(baseline['RH']):.0f}",
            f"{float(baseline['v']):.2f}",
            f"{float(baseline['PMV']):.2f}",
            f"{float(baseline['PPD']):.1f}",
        ])

    system_rows = [
        ["All-air sensible control", r"$Q = \rho c_p \dot V |\Delta T_a|$", f"{DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['supply_air_flow_m3s']:.2f} m$^3$/s, COP = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_air_heating']:.1f}/{DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_air_cooling']:.1f}"],
        ["Radiant control", r"$Q = h_r A |\Delta T_r|$", f"$h_r$ = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['h_r_panel_w_m2k']:.1f} W/m$^2$K, $A$ = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['radiant_active_area_m2']:.0f} m$^2$, COP = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_radiant_heating']:.1f}/{DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_radiant_cooling']:.1f}"],
        ["Airflow control", r"$P \propto v^3$", f"fan power max = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['fan_power_max_w']:.0f} W at {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['fan_v_max_m_s']:.1f} m/s"],
        ["Latent control", r"$Q = \dot m h_{{fg}} |\Delta W|$", f"$h_{{fg}}$ = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['latent_heat_vaporization_j_kg']:.2e} J/kg, COP = {DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_dehumidification']:.1f}/{DISCOMFORT_MOD.SYSTEM_ASSUMPTIONS['cop_humidification']:.1f}"],
    ]

    step_example_rows: List[List[str]] = []
    for profile in PMV_MOD.PROFILES:
        if profile.id not in {"reference_office", "larger_male"}:
            continue
        sub_steps = [row for row in step_rows if row["profile_id"] == profile.id and row["family_id"] == "neutral_still"]
        sub_steps.sort(key=lambda row: float(row["start_label"]))
        for row in sub_steps:
            step_example_rows.append([
                latex_escape(profile.label),
                latex_escape(row["transition_label"]),
                f"{float(row['delta_ta']):.2f}",
                f"{float(row['load_change_w']):.1f}",
                f"{float(row['ppd_improvement']):.1f}",
                f"{float(row['dQ_conv_w']):.1f}",
                f"{float(row['dQ_rad_w']):.1f}",
            ])

    outcome_rows: List[List[str]] = []
    control_winners = 0
    radiant_exergy_winners = 0
    for scenario in DISCOMFORT_MOD.SCENARIOS:
        control = next(row for row in discomfort_winners if row["scenario_id"] == scenario.id)
        exergy = next(row for row in discomfort_exergy if row["scenario_id"] == scenario.id)
        if control["lever_label"] == "All-air (Ta only)":
            control_winners += 1
        if exergy["lever_label"] == "Radiant (Tr only)":
            radiant_exergy_winners += 1
        air = next(row for row in discomfort_recoveries if row["scenario_id"] == scenario.id and row["lever_id"] == "air_temp")
        radiant_candidates = [row for row in discomfort_recoveries if row["scenario_id"] == scenario.id and row["lever_id"] == "radiant"]
        radiant = radiant_candidates[0] if radiant_candidates else None
        outcome_rows.append([
            latex_escape(scenario.label),
            f"{float(air['control_move']):.1f}",
            f"{float(air['recovered_pmv']):.2f}",
            f"{float(air['system_exergy_w']):.1f}",
            latex_escape(air["dominant_channel"]),
            "—" if radiant is None or radiant.get("feasible") == "False" else f"{float(radiant['control_move']):.1f}",
            "—" if radiant is None or radiant.get("feasible") == "False" else f"{float(radiant['recovered_pmv']):.2f}",
            "—" if radiant is None or radiant.get("feasible") == "False" else f"{float(radiant['system_exergy_w']):.1f}",
            "—" if radiant is None or radiant.get("feasible") == "False" else latex_escape(radiant["dominant_channel"]),
        ])

    pmv_factor_ref = next(float(row["pmv_factor"]) for row in pmv_load_rows if abs(float(row["met"]) - 0.9) < 1e-9)

    macros["ScenarioCount"] = str(len(DISCOMFORT_MOD.SCENARIOS))
    macros["ProfileCount"] = str(len(PMV_MOD.PROFILES))
    macros["RadiantExergyWins"] = str(radiant_exergy_winners)
    macros["AllAirControlWins"] = str(control_winners)
    macros["ReferencePMVFactor"] = f"{pmv_factor_ref:.4f}"
    macros["OneStepLoadMin"] = f"{min(one_step_loads):.1f}"
    macros["OneStepLoadMax"] = f"{max(one_step_loads):.1f}"
    macros["OneStepPPDMin"] = f"{min(step_ppd_changes):.1f}"
    macros["OneStepPPDMax"] = f"{max(step_ppd_changes):.1f}"
    macros["AirImproveMin"] = f"{min(air_improvements):.2f}"
    macros["AirImproveMax"] = f"{max(air_improvements):.2f}"
    macros["RadImproveMin"] = f"{min(rad_improvements):.2f}"
    macros["RadImproveMax"] = f"{max(rad_improvements):.2f}"
    macros["AirPPDImproveMin"] = f"{min(air_ppd_improvements):.1f}"
    macros["AirPPDImproveMax"] = f"{max(air_ppd_improvements):.1f}"
    macros["RadPPDImproveMin"] = f"{min(rad_ppd_improvements):.1f}"
    macros["RadPPDImproveMax"] = f"{max(rad_ppd_improvements):.1f}"
    macros["PPDAtHalf"] = f"{PMV_MOD.ppd_from_pmv(0.5):.1f}"
    macros["PPDAtOne"] = f"{PMV_MOD.ppd_from_pmv(1.0):.1f}"

    write_table(
        TABLES / "table_profiles.tex",
        "Normative occupant profiles used in the profile-sensitive PMV step and 1$^\\circ$C perturbation analyses.",
        "tab:profiles",
        "p{3.1cm}cccccccc",
        ["Profile", "Age", "Sex", "Height", "Weight", "BSA", "Clo", "Met", "W/person"],
        profiles_table_rows,
        notes="Profiles are illustrative control-design archetypes. Age and sex are not explicit PMV inputs; they enter only through the assumed body size, clothing insulation, and metabolic rate.",
    )
    write_table(
        TABLES / "table_scenarios.tex",
        "Severe discomfort archetypes used for recovery and exergy comparisons.",
        "tab:scenarios",
        "p{3.3cm}p{1.8cm}cccccc",
        ["Scenario", "Class", "$T_a$", "$T_r$", "RH", "$v$", "PMV", "PPD"],
        severe_scenario_rows,
        notes="All severe scenarios satisfy $|PMV| \\geq 2$ for the representative office occupant before recovery.",
    )
    write_table(
        TABLES / "table_system.tex",
        "System-side surrogate assumptions used to convert control moves into approximate exergy/electric input.",
        "tab:system",
        "p{3.0cm}p{4.5cm}p{6.0cm}",
        ["Lever", "Surrogate relation", "Assumptions"],
        system_rows,
        notes="These are low-order plant surrogates used only for ranking, not full HVAC models.",
    )
    write_table(
        TABLES / "table_step_summary.tex",
        "Representative one-label PMV transitions under the neutral-still family for two office profiles.",
        "tab:stepsummary",
        "p{3.0cm}p{1.0cm}cccccc",
        ["Profile", "Step", "$\\Delta T_a$", "$\\Delta(H-L)$", "$\\Delta PPD$", "$\\Delta Q_{conv}$", "$\\Delta Q_{rad}$"],
        step_example_rows,
        notes="All rows are taken from the neutral-still environmental family to keep the evidence readable. The full family-to-family variation appears in Figure~\\ref{fig:stepresponse}. Equal PMV steps retain the same body-load step within a profile, but not the same air-temperature move or dissatisfaction change.",
    )
    write_table(
        TABLES / "table_outcomes.tex",
        "Air-versus-radiant comparison for the six severe recovery archetypes.",
        "tab:outcomes",
        "p{2.5cm}cccccccc",
        ["Scenario", "Air move", "Air PMV", "Air exergy", "Air dom.", "Rad move", "Rad PMV", "Rad exergy", "Rad dom."],
        outcome_rows,
        notes="Air and radiant are shown side by side rather than only as winners. The recovery-ranking reversal is meaningful because both levers reach nearly the same recovered PMV in five of the six scenarios, while their estimated exergy requirements differ substantially.",
    )

    macro_lines = []
    for key, value in sorted(macros.items()):
        macro_lines.append(rf"\newcommand{{\{key}}}{{{value}}}")
    write_text(GENERATED / "macros.tex", "\n".join(macro_lines) + "\n")

    summary = {
        "macros": macros,
        "profile_step_summary": pmv_step_summary_rows,
        "outcome_rows": outcome_rows,
    }
    write_text(GENERATED / "summary.json", json.dumps(summary, indent=2))
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    GENERATED.mkdir(parents=True, exist_ok=True)
    figures = build_figures()
    summary = build_tables_and_macros()
    print(f"Built figures: {len(figures)}")
    print(f"Built summary tables and macros at {GENERATED / 'macros.tex'}")
    print(json.dumps(summary["macros"], indent=2))


if __name__ == "__main__":
    main()
