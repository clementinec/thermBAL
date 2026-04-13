#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CASE_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "thermbal_mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

from build_templates import CELLULAR_TEMPLATE, OPEN_TEMPLATE
from generate_study import (
    WORK_ROOM_ORDER,
    WORK_ROOM_ABBREVIATIONS,
    build_agents,
    compute_comfort,
    farthest_point_sample,
    make_field,
    room_average_envs,
    sample_count_for_room,
)


OUT_DIR = CASE_DIR / "out"
FIG_DIR = OUT_DIR / "figures"
FIELD_MAPS_PATH = OUT_DIR / "field_maps.json"
SUMMARY_PATH = OUT_DIR / "summary.csv"
SAMPLES_PATH = OUT_DIR / "sample_results.csv"

FOCUS_TA = 22.0
FOCUS_RH = 80.0
DISPLAY_COHORTS = [
    ("default_male_35", "Default Male 35"),
    ("realistic_mixed", "Realistic Mixed"),
    ("female_light", "Female Light"),
]
FIELD_COHORT = ("realistic_mixed", "Realistic Mixed")
DENSITY_COHORTS = [
    ("default_male_35", "Default Male 35"),
    ("realistic_male", "Realistic Male"),
    ("realistic_mixed", "Realistic Mixed"),
    ("female_light", "Female Light"),
]

WORK_FILL = "#efe7dc"
GRID_EDGE = "#f7f2ea"
BOUNDARY = "#7e7469"
TEXT = "#2d2926"
MUTED = "#6d645b"
RULE = "#d7cab9"
CELLULAR_TOPOLOGY = "#b48b60"
OPEN_TOPOLOGY = "#567fae"
CELLULAR_COMPARE = "#6b6158"
OPEN_COMPARE = "#5f7d5f"


def read_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: object, default: float = 0.0) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def build_field_lookup(rows: Iterable[Dict[str, object]]) -> Dict[Tuple[str, str, float, float], Dict[str, object]]:
    out: Dict[Tuple[str, str, float, float], Dict[str, object]] = {}
    for row in rows:
        out[
            (
                str(row["topology"]),
                str(row["cohort"]),
                as_float(row["base_ta_c"]),
                as_float(row["base_rh_pct"]),
            )
        ] = row
    return out


def build_summary_lookup(rows: Iterable[Dict[str, str]]) -> Dict[Tuple[str, str, float, float], Dict[str, str]]:
    out: Dict[Tuple[str, str, float, float], Dict[str, str]] = {}
    for row in rows:
        out[
            (
                str(row["topology"]),
                str(row["cohort"]),
                as_float(row["base_ta_c"]),
                as_float(row["base_rh_pct"]),
            )
        ] = row
    return out


def cell_map(entry: Dict[str, object]) -> Dict[Tuple[int, int], Dict[str, object]]:
    return {(int(cell["x"]), int(cell["y"])): cell for cell in entry["cells"]}


def field_range(entry: Dict[str, object]) -> float:
    values = [as_float(cell["pmv"]) for cell in entry["cells"] if bool(cell["is_work_cell"])]
    return max(values) - min(values)


def boundary_segments(cells_by_coord: Dict[Tuple[int, int], object]) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    segments = set()
    for (col, row), cell in cells_by_coord.items():
        room_name = str(cell["room_name"]) if isinstance(cell, dict) else str(cell)
        edges = [
            ((col, row), (col + 1, row), (col, row - 1)),
            ((col + 1, row), (col + 1, row + 1), (col + 1, row)),
            ((col, row + 1), (col + 1, row + 1), (col, row + 1)),
            ((col, row), (col, row + 1), (col - 1, row)),
        ]
        for start, end, neighbor_coord in edges:
            neighbor = cells_by_coord.get(neighbor_coord)
            if neighbor is None:
                segments.add((start, end))
                continue
            neighbor_name = str(neighbor["room_name"]) if isinstance(neighbor, dict) else str(neighbor)
            if neighbor_name != room_name:
                segments.add((start, end))
    return sorted(segments)


def plan_cells_from_payload(payload: Dict[str, object]) -> Dict[Tuple[int, int], str]:
    out: Dict[Tuple[int, int], str] = {}
    for room in payload["rooms"]:
        room_name = str(room["name"])
        for col, row in room["cells"]:
            out[(int(col), int(row))] = room_name
    return out


def work_rooms_for_payload(payload: Dict[str, object]) -> set[str]:
    room_names = {str(room["name"]) for room in payload["rooms"]}
    if "North Open Office" in room_names:
        return {"North Open Office", "South Open Office", "West Corner Office"}
    return {
        "Reception",
        "Admin Bay",
        "Team Office West",
        "Team Office East",
        "Executive Office",
        "Focus Room",
        "Print Hub",
        "Collaboration West",
        "Collaboration East",
        "Meeting Room West",
        "Meeting Room East",
        "South Studio West",
        "South Studio East",
        "Project Bay West",
        "Project Bay East",
    }


def room_index(payload: Dict[str, object]) -> Dict[str, set[Tuple[int, int]]]:
    return {
        str(room["name"]): {tuple(cell) for cell in room["cells"]}
        for room in payload["rooms"]
    }


def interpolate_metric(
    records: List[Dict[str, str]],
    payload: Dict[str, object],
    metric_key: str,
) -> Dict[Tuple[int, int], float]:
    work_rooms = work_rooms_for_payload(payload)
    cell_to_values: Dict[Tuple[int, int], List[float]] = {}
    for record in records:
        cell = (int(float(record["x"])), int(float(record["y"])))
        cell_to_values.setdefault(cell, []).append(as_float(record[metric_key]))

    sample_points = [
        (int(float(record["x"])), int(float(record["y"])), as_float(record[metric_key]))
        for record in records
    ]

    out: Dict[Tuple[int, int], float] = {}
    for room in payload["rooms"]:
        room_name = str(room["name"])
        if room_name not in work_rooms:
            continue
        for cell in (tuple(item) for item in room["cells"]):
            if cell in cell_to_values:
                out[cell] = statistics.mean(cell_to_values[cell])
                continue

            weighted = 0.0
            weights = 0.0
            for sx, sy, value in sample_points:
                dist_sq = (cell[0] - sx) ** 2 + (cell[1] - sy) ** 2
                weight = 1.0 / max(0.6, dist_sq) ** 1.15
                weighted += weight * value
                weights += weight
            out[cell] = weighted / weights if weights else 0.0
    return out


def draw_field_panel(
    ax: plt.Axes,
    payload: Dict[str, object],
    entry: Dict[str, object],
    *,
    cmap: mcolors.Colormap,
    norm: mcolors.Normalize,
    row_label: str | None = None,
) -> None:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    payload_cells = plan_cells_from_payload(payload)
    cells = cell_map(entry)

    ax.set_xlim(0, cols)
    ax.set_ylim(rows, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    for (col, row), room_name in payload_cells.items():
        cell = cells.get((col, row))
        if cell is None:
            facecolor = "#fbf7f0"
        elif bool(cell["is_work_cell"]):
            facecolor = cmap(norm(as_float(cell["pmv"])))
        else:
            facecolor = WORK_FILL
        ax.add_patch(
            Rectangle(
                (col, row),
                1.0,
                1.0,
                facecolor=facecolor,
                edgecolor=GRID_EDGE,
                linewidth=0.35,
            )
        )

    lines = LineCollection(boundary_segments(cells), colors=BOUNDARY, linewidths=0.8, zorder=5)
    ax.add_collection(lines)

    if row_label:
        ax.text(
            -1.55,
            rows / 2.0,
            row_label,
            rotation=90,
            va="center",
            ha="center",
            fontsize=10.5,
            color=TEXT,
            fontweight="bold",
        )


def metric_arrays(
    payload: Dict[str, object],
    pmv_map: Dict[Tuple[int, int], float],
    ppd_map: Dict[Tuple[int, int], float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    pmv = np.full((rows, cols), np.nan)
    ppd = np.full((rows, cols), np.nan)
    support = np.zeros((rows, cols), dtype=bool)

    room_cells = room_index(payload)
    work_rooms = work_rooms_for_payload(payload)
    for room in payload["rooms"]:
        room_name = str(room["name"])
        for col, row in room_cells[room_name]:
            if room_name not in work_rooms:
                support[row, col] = True
                continue
            pmv[row, col] = pmv_map.get((col, row), np.nan)
            ppd[row, col] = ppd_map.get((col, row), np.nan)
    return pmv, ppd, support


def draw_smooth_field_panel(
    ax: plt.Axes,
    payload: Dict[str, object],
    pmv_map: Dict[Tuple[int, int], float],
    ppd_map: Dict[Tuple[int, int], float],
    summary: Dict[str, str],
    span_value: float,
    *,
    cmap: mcolors.Colormap,
    norm: mcolors.Normalize,
    row_label: str | None = None,
) -> None:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    payload_cells = plan_cells_from_payload(payload)
    pmv, ppd, support = metric_arrays(payload, pmv_map, ppd_map)

    ax.set_xlim(0, cols)
    ax.set_ylim(rows, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    for (col, row) in payload_cells:
        fill = "#ebe2d6" if support[row, col] else "#f7f0e6"
        ax.add_patch(
            Rectangle(
                (col, row),
                1.0,
                1.0,
                facecolor=fill,
                edgecolor="none",
                linewidth=0.0,
                zorder=0,
            )
        )

    heat = np.ma.masked_invalid(pmv)
    heat_cmap = cmap.copy()
    heat_cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    ax.imshow(
        heat,
        cmap=heat_cmap,
        norm=norm,
        interpolation="bicubic",
        extent=(0, cols, rows, 0),
        zorder=1,
    )

    lines = LineCollection(boundary_segments(payload_cells), colors=BOUNDARY, linewidths=0.75, zorder=4)
    ax.add_collection(lines)

    stats = (
        f"PMV {as_float(summary['mean_pmv']):.2f}   "
        f"PPD {as_float(summary['mean_ppd']):.0f}%   "
        f"span {span_value:.3f}"
    )
    ax.text(
        0.03,
        0.03,
        stats,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.24", facecolor="#fffaf3", edgecolor="#d8ccbc", linewidth=0.6),
        zorder=6,
    )

    if row_label:
        ax.text(
            -1.45,
            rows / 2.0,
            row_label,
            rotation=90,
            va="center",
            ha="center",
            fontsize=10.2,
            color=TEXT,
            fontweight="bold",
        )


def make_pmv_cmap() -> mcolors.Colormap:
    return mcolors.LinearSegmentedColormap.from_list(
        "thermbal_pmv",
        [
            (0.00, "#355c8a"),
            (0.23, "#6d95be"),
            (0.47, "#b9d0e1"),
            (0.63, "#ebe3d4"),
            (0.80, "#d9c7a8"),
            (1.00, "#c68654"),
        ],
    )


def make_delta_cmap() -> mcolors.Colormap:
    return mcolors.LinearSegmentedColormap.from_list(
        "thermbal_delta",
        [
            (0.00, "#4d78a3"),
            (0.48, "#dce3e9"),
            (0.50, "#f5f1ea"),
            (0.72, "#ddb189"),
            (1.00, "#c66b3d"),
        ],
    )


def kde_curve(values: List[float], x_grid: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.zeros_like(x_grid)
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.12
    bandwidth = max(0.10, 1.06 * std * (arr.size ** (-1.0 / 5.0)))
    kernels = np.exp(-0.5 * ((x_grid[:, None] - arr[None, :]) / bandwidth) ** 2)
    density = kernels.sum(axis=1) / (arr.size * bandwidth * math.sqrt(2.0 * math.pi))
    return density


def build_work_points(cellular_payload_state) -> List[Dict[str, object]]:
    room_cells = {room.name: set(room.cells) for room in cellular_payload_state.plan.rooms}
    work_points: List[Dict[str, object]] = []
    for room_name in WORK_ROOM_ORDER:
        cells = sorted(room_cells[room_name], key=lambda item: (item[1], item[0]))
        sample_count = sample_count_for_room(len(cells))
        for point_index, cell in enumerate(farthest_point_sample(cells, sample_count), start=1):
            work_points.append(
                {
                    "sample_id": f"{WORK_ROOM_ABBREVIATIONS[room_name]}_{point_index:02d}",
                    "source_room": room_name,
                    "cell": list(cell),
                }
            )
    return work_points


def room_anchored_delta_data():
    from floor_plan.io import load_template_json

    cellular_state = load_template_json(CELLULAR_TEMPLATE)
    active_cells = set(cellular_state.active_cells)
    cellular_rooms = {room.name: set(room.cells) for room in cellular_state.plan.rooms}
    field = make_field(active_cells, cellular_state.features, FOCUS_TA, FOCUS_RH)
    cellular_room_env = room_average_envs(cellular_rooms, field)
    work_points = build_work_points(cellular_state)
    agents_by_cohort = build_agents(work_points)

    density_values: Dict[str, List[float]] = {}
    map_values: Dict[str, Dict[Tuple[int, int], float]] = {}
    for cohort_key, agents in agents_by_cohort.items():
        deltas: List[float] = []
        per_cell_values: Dict[Tuple[int, int], List[float]] = defaultdict(list)
        for agent in agents:
            source_room = str(agent["source_room"])
            source_cells = cellular_rooms[source_room]
            cell_env = cellular_room_env[source_room]
            cellular_pmv = compute_comfort(agent, cell_env)["pmv"]
            for cell in source_cells:
                open_pmv = compute_comfort(agent, field[cell])["pmv"]
                delta = open_pmv - cellular_pmv
                deltas.append(delta)
                per_cell_values[cell].append(delta)
        density_values[cohort_key] = deltas
        map_values[cohort_key] = {
            cell: statistics.mean(values)
            for cell, values in per_cell_values.items()
        }
    return cellular_state, density_values, map_values


def draw_delta_map_panel(
    ax: plt.Axes,
    payload: Dict[str, object],
    delta_map: Dict[Tuple[int, int], float],
    *,
    cmap: mcolors.Colormap,
    norm: mcolors.Normalize,
) -> None:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    payload_cells = plan_cells_from_payload(payload)
    room_cells = room_index(payload)
    work_rooms = work_rooms_for_payload(payload)

    ax.set_xlim(0, cols)
    ax.set_ylim(rows, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    for room in payload["rooms"]:
        room_name = str(room["name"])
        for col, row in room_cells[room_name]:
            if room_name in work_rooms and (col, row) in delta_map:
                facecolor = cmap(norm(delta_map[(col, row)]))
            elif room_name in work_rooms:
                facecolor = "#f3ecdf"
            else:
                facecolor = "#ece3d7"
            ax.add_patch(
                Rectangle(
                    (col, row),
                    1.0,
                    1.0,
                    facecolor=facecolor,
                    edgecolor="none",
                    linewidth=0.0,
                    zorder=0,
                )
            )

    lines = LineCollection(boundary_segments(payload_cells), colors=BOUNDARY, linewidths=0.75, zorder=3)
    ax.add_collection(lines)


def render_distribution_story_figure(
    records_lookup: Dict[Tuple[str, str, float, float], List[Dict[str, str]]],
    summary_lookup: Dict[Tuple[str, str, float, float], Dict[str, str]],
    cellular_payload: Dict[str, object],
    open_payload: Dict[str, object],
) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    focus_key = FIELD_COHORT[0]
    cell_records = records_lookup[("cellular_office", focus_key, FOCUS_TA, FOCUS_RH)]
    open_records = records_lookup[("open_office", focus_key, FOCUS_TA, FOCUS_RH)]
    cell_summary = summary_lookup[("cellular_office", focus_key, FOCUS_TA, FOCUS_RH)]
    open_summary = summary_lookup[("open_office", focus_key, FOCUS_TA, FOCUS_RH)]
    cell_interp_pmv = interpolate_metric(cell_records, cellular_payload, "pmv")
    open_interp_pmv = interpolate_metric(open_records, open_payload, "pmv")
    cell_interp_ppd = interpolate_metric(cell_records, cellular_payload, "ppd")
    open_interp_ppd = interpolate_metric(open_records, open_payload, "ppd")

    all_field_values = list(cell_interp_pmv.values()) + list(open_interp_pmv.values())
    pmv_norm = mcolors.Normalize(
        vmin=math.floor((min(all_field_values) - 0.05) * 10.0) / 10.0,
        vmax=math.ceil((max(all_field_values) + 0.05) * 10.0) / 10.0,
    )
    pmv_cmap = make_pmv_cmap()

    all_sample_values = []
    for cohort_key, _ in DENSITY_COHORTS:
        for topology in ("cellular_office", "open_office"):
            values = [
                as_float(row["pmv"])
                for row in records_lookup[(topology, cohort_key, FOCUS_TA, FOCUS_RH)]
            ]
            all_sample_values.extend(values)
    x_min = math.floor((min(all_sample_values) - 0.15) * 10.0) / 10.0
    x_max = math.ceil((max(all_sample_values) + 0.15) * 10.0) / 10.0
    x_grid = np.linspace(x_min, x_max, 280)

    fig = plt.figure(figsize=(11.6, 8.2), facecolor="white")
    grid = fig.add_gridspec(
        nrows=2,
        ncols=4,
        left=0.06,
        right=0.97,
        top=0.84,
        bottom=0.14,
        hspace=0.34,
        wspace=0.22,
        height_ratios=[1.10, 0.82],
    )

    fig.text(
        0.06,
        0.965,
        "Office Topology Counterfactual: Spatial Field and Population Spread",
        fontsize=15.2,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )
    fig.text(
        0.06,
        0.922,
        "Top: interpolated PMV field for the realistic mixed cohort on the shared footprint at 22 C / 80% RH.\n"
        "Bottom: sampled worker PMV distributions for four cohorts. Open plan changes mean PMV very little but broadens the visible span.",
        fontsize=8.9,
        color=MUTED,
        ha="left",
        va="top",
    )

    top_left = fig.add_subplot(grid[0, 0:2])
    top_right = fig.add_subplot(grid[0, 2:4])
    draw_smooth_field_panel(
        top_left,
        cellular_payload,
        cell_interp_pmv,
        cell_interp_ppd,
        cell_summary,
        max(cell_interp_pmv.values()) - min(cell_interp_pmv.values()),
        cmap=pmv_cmap,
        norm=pmv_norm,
    )
    draw_smooth_field_panel(
        top_right,
        open_payload,
        open_interp_pmv,
        open_interp_ppd,
        open_summary,
        max(open_interp_pmv.values()) - min(open_interp_pmv.values()),
        cmap=pmv_cmap,
        norm=pmv_norm,
    )
    top_left.set_title(
        f"Cellular Office  |  {FIELD_COHORT[1]}",
        fontsize=9.9,
        color=TEXT,
        pad=6,
    )
    top_right.set_title(
        f"Open Office  |  {FIELD_COHORT[1]}",
        fontsize=9.9,
        color=TEXT,
        pad=6,
    )

    cax = fig.add_axes([0.40, 0.535, 0.22, 0.022])
    sm = plt.cm.ScalarMappable(norm=pmv_norm, cmap=pmv_cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#a29281")
    cbar.ax.tick_params(labelsize=8.2, colors=MUTED)
    cbar.set_label("Interpolated PMV", fontsize=9.5, color=TEXT)

    legend_handles = None
    for idx, (cohort_key, cohort_label) in enumerate(DENSITY_COHORTS):
        ax = fig.add_subplot(grid[1, idx])
        cell_values = [as_float(row["pmv"]) for row in records_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]]
        open_values = [as_float(row["pmv"]) for row in records_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]]
        cell_density = kde_curve(cell_values, x_grid)
        open_density = kde_curve(open_values, x_grid)

        ax.axvspan(-0.5, 0.5, color="#ece4d7", alpha=0.55, zorder=0)
        ax.axvline(0.0, color=RULE, linewidth=0.7, zorder=1)
        cell_fill = ax.fill_between(x_grid, cell_density, color=CELLULAR_TOPOLOGY, alpha=0.18, zorder=2)
        open_fill = ax.fill_between(x_grid, open_density, color=OPEN_TOPOLOGY, alpha=0.18, zorder=2)
        cell_line = ax.plot(x_grid, cell_density, color=CELLULAR_TOPOLOGY, linewidth=1.8, zorder=3)[0]
        open_line = ax.plot(x_grid, open_density, color=OPEN_TOPOLOGY, linewidth=1.8, zorder=3)[0]
        ax.axvline(statistics.mean(cell_values), color=CELLULAR_TOPOLOGY, linewidth=1.0, linestyle="--", zorder=4)
        ax.axvline(statistics.mean(open_values), color=OPEN_TOPOLOGY, linewidth=1.0, linestyle="--", zorder=4)

        if legend_handles is None:
            legend_handles = [cell_line, open_line]

        max_density = max(float(cell_density.max()), float(open_density.max()))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0.0, max_density * 1.16)
        ax.set_title(cohort_label, fontsize=9.7, color=TEXT, pad=8)
        ax.set_xlabel("PMV", fontsize=8.2, color=TEXT)
        if idx == 0:
            ax.set_ylabel("Density", fontsize=8.2, color=TEXT)
        else:
            ax.set_yticklabels([])
        ax.tick_params(axis="both", labelsize=7.8, colors=MUTED)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(RULE)
        ax.spines["bottom"].set_color(RULE)
        ax.grid(axis="y", color="#ece3d7", linewidth=0.6, alpha=0.7)

        cell_span = max(cell_values) - min(cell_values)
        open_span = max(open_values) - min(open_values)
        ax.text(
            0.02,
            0.96,
            f"C {statistics.mean(cell_values):+.2f} | span {cell_span:.3f}\n"
            f"O {statistics.mean(open_values):+.2f} | span {open_span:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.22", facecolor="#fffaf3", edgecolor="#ddd1c2", linewidth=0.6),
        )

    if legend_handles is not None:
        fig.legend(
            legend_handles,
            ["Cellular office samples", "Open office samples"],
            loc="lower center",
            bbox_to_anchor=(0.52, 0.055),
            ncol=2,
            frameon=False,
            fontsize=8.6,
        )

    fig.text(
        0.06,
        0.030,
        "Bottom-row densities use the shared desk/sample points at 22 C / 80% RH. Dashed markers indicate topology-specific cohort means; span values are max-minus-min sampled PMV.",
        fontsize=8.6,
        color=MUTED,
        ha="left",
    )

    svg_path = FIG_DIR / "fig_office_topology_distribution_story.svg"
    png_path = FIG_DIR / "fig_office_topology_distribution_story_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_delta_distribution_figure(cellular_payload: Dict[str, object]) -> None:
    cellular_state, density_values, map_values = room_anchored_delta_data()

    all_deltas = [value for cohort_key, _ in DENSITY_COHORTS for value in density_values[cohort_key]]
    delta_min = min(all_deltas)
    delta_max = max(all_deltas)
    lim = max(abs(delta_min), abs(delta_max))
    delta_norm = mcolors.TwoSlopeNorm(
        vmin=round(-lim - 0.02, 2),
        vcenter=0.0,
        vmax=round(lim + 0.02, 2),
    )
    delta_cmap = make_delta_cmap()
    x_grid = np.linspace(delta_norm.vmin, delta_norm.vmax, 280)

    fig = plt.figure(figsize=(11.6, 7.9), facecolor="white")
    grid = fig.add_gridspec(
        nrows=2,
        ncols=4,
        left=0.06,
        right=0.97,
        top=0.86,
        bottom=0.12,
        hspace=0.34,
        wspace=0.18,
        height_ratios=[1.02, 0.86],
    )

    fig.text(
        0.06,
        0.965,
        "Open minus Cellular PMV: Room-Anchored Spatial Discrepancy",
        fontsize=15.2,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )
    fig.text(
        0.06,
        0.922,
        "Each agent is evaluated across every cell in its original cellular room footprint, then compared against the same coordinates under the open-office field.\n"
        "Top: mean open-minus-cellular delta for the realistic mixed cohort. Bottom: delta-density span for four cohorts on the same axis.",
        fontsize=8.9,
        color=MUTED,
        ha="left",
        va="top",
    )

    top_ax = fig.add_subplot(grid[0, :])
    draw_delta_map_panel(
        top_ax,
        cellular_payload,
        map_values["realistic_mixed"],
        cmap=delta_cmap,
        norm=delta_norm,
    )
    top_ax.set_title("Realistic Mixed  |  mean PMV delta across original room footprints", fontsize=10.2, color=TEXT, pad=8)
    mixed_vals = density_values["realistic_mixed"]
    top_ax.text(
        0.015,
        0.03,
        f"mean Δ {statistics.mean(mixed_vals):+.3f}   min {min(mixed_vals):+.3f}   max {max(mixed_vals):+.3f}   span {max(mixed_vals)-min(mixed_vals):.3f}",
        transform=top_ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.2,
        color=TEXT,
        bbox=dict(boxstyle="round,pad=0.24", facecolor="#fffaf3", edgecolor="#d8ccbc", linewidth=0.6),
    )

    cax = fig.add_axes([0.40, 0.515, 0.22, 0.022])
    sm = plt.cm.ScalarMappable(norm=delta_norm, cmap=delta_cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#a29281")
    cbar.ax.tick_params(labelsize=8.2, colors=MUTED)
    cbar.set_label("Open minus cellular PMV", fontsize=9.5, color=TEXT)

    cohort_line_colors = {
        "default_male_35": "#8d6a44",
        "realistic_male": "#9a7a56",
        "realistic_mixed": "#567fae",
        "female_light": "#385f8d",
    }

    for idx, (cohort_key, cohort_label) in enumerate(DENSITY_COHORTS):
        ax = fig.add_subplot(grid[1, idx])
        values = density_values[cohort_key]
        density = kde_curve(values, x_grid)
        color = cohort_line_colors[cohort_key]
        ax.axvspan(-0.05, 0.05, color="#ede5d9", alpha=0.70, zorder=0)
        ax.axvline(0.0, color=RULE, linewidth=0.8, zorder=1)
        ax.fill_between(x_grid, density, color=color, alpha=0.18, zorder=2)
        ax.plot(x_grid, density, color=color, linewidth=1.9, zorder=3)
        ax.axvline(statistics.mean(values), color=color, linewidth=1.0, linestyle="--", zorder=4)
        ax.set_xlim(delta_norm.vmin, delta_norm.vmax)
        ax.set_ylim(0.0, float(density.max()) * 1.18)
        ax.set_title(cohort_label, fontsize=9.7, color=TEXT, pad=8)
        ax.set_xlabel("Open minus cellular PMV", fontsize=8.0, color=TEXT)
        if idx == 0:
            ax.set_ylabel("Density", fontsize=8.2, color=TEXT)
        else:
            ax.set_yticklabels([])
        ax.tick_params(axis="both", labelsize=7.8, colors=MUTED)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color(RULE)
        ax.spines["bottom"].set_color(RULE)
        ax.grid(axis="y", color="#ece3d7", linewidth=0.6, alpha=0.7)
        ax.text(
            0.02,
            0.96,
            f"mean {statistics.mean(values):+.3f}\n[{min(values):+.3f}, {max(values):+.3f}]",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.22", facecolor="#fffaf3", edgecolor="#ddd1c2", linewidth=0.6),
        )

    fig.text(
        0.06,
        0.038,
        "Positive values mean the open-office interpretation reads warmer or less cold than the cellular office at the same original room coordinates.",
        fontsize=8.6,
        color=MUTED,
        ha="left",
    )

    svg_path = FIG_DIR / "fig_office_topology_delta_distribution.svg"
    png_path = FIG_DIR / "fig_office_topology_delta_distribution_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_story_with_density_sidecar_figure() -> None:
    story_png = FIG_DIR / "fig_office_topology_story_hires.png"
    if not story_png.exists():
        return

    _, density_values, _ = room_anchored_delta_data()
    all_deltas = [value for cohort_key, _ in DENSITY_COHORTS for value in density_values[cohort_key]]
    delta_min = min(all_deltas)
    delta_max = max(all_deltas)
    lim = max(abs(delta_min), abs(delta_max))
    y_min = round(-lim - 0.02, 2)
    y_max = round(lim + 0.02, 2)
    y_grid = np.linspace(y_min, y_max, 320)

    fig = plt.figure(figsize=(14.8, 7.2), facecolor="white")
    grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.03,
        right=0.98,
        top=0.98,
        bottom=0.08,
        width_ratios=[3.4, 1.15],
        wspace=0.08,
    )

    ax_story = fig.add_subplot(grid[0, 0])
    image = plt.imread(story_png)
    ax_story.imshow(image)
    ax_story.axis("off")

    ax_density = fig.add_subplot(grid[0, 1])
    ax_density.set_title("Open minus cellular\nroom-anchored PMV span", fontsize=10.8, color=TEXT, pad=10)
    ax_density.axhline(0.0, color=RULE, linewidth=0.9, zorder=1)
    ax_density.grid(axis="y", color="#ece3d7", linewidth=0.6, alpha=0.8)
    ax_density.set_axisbelow(True)
    ax_density.set_xlim(0.35, len(DENSITY_COHORTS) + 0.65)
    ax_density.set_ylim(y_min, y_max)
    ax_density.set_ylabel("Open minus cellular PMV", fontsize=8.8, color=TEXT)
    ax_density.tick_params(axis="y", labelsize=8.0, colors=MUTED)
    ax_density.tick_params(axis="x", labelsize=7.8, colors=MUTED)
    ax_density.spines["top"].set_visible(False)
    ax_density.spines["right"].set_visible(False)
    ax_density.spines["left"].set_color(RULE)
    ax_density.spines["bottom"].set_color(RULE)

    cohort_line_colors = {
        "default_male_35": "#8d6a44",
        "realistic_male": "#9a7a56",
        "realistic_mixed": "#567fae",
        "female_light": "#385f8d",
    }
    max_density = 0.0
    densities = {}
    for cohort_key, _ in DENSITY_COHORTS:
        density = kde_curve(density_values[cohort_key], y_grid)
        densities[cohort_key] = density
        max_density = max(max_density, float(density.max()))

    scale = 0.34 / max_density if max_density else 0.0
    xticks = []
    xlabels = []
    for idx, (cohort_key, cohort_label) in enumerate(DENSITY_COHORTS, start=1):
        xticks.append(idx)
        if cohort_key == "default_male_35":
            xlabels.append("Default\nMale")
        elif cohort_key == "realistic_male":
            xlabels.append("Realistic\nMale")
        elif cohort_key == "realistic_mixed":
            xlabels.append("Realistic\nMixed")
        else:
            xlabels.append("Female\nLight")

        color = cohort_line_colors[cohort_key]
        values = density_values[cohort_key]
        density = densities[cohort_key]
        width = density * scale
        ax_density.fill_betweenx(
            y_grid,
            idx - width,
            idx + width,
            facecolor=color,
            alpha=0.22,
            linewidth=0.0,
            zorder=2,
        )
        ax_density.plot(idx + width, y_grid, color=color, linewidth=1.6, zorder=3)
        ax_density.plot(idx - width, y_grid, color=color, linewidth=1.6, zorder=3)
        mean_delta = statistics.mean(values)
        ax_density.plot([idx - 0.18, idx + 0.18], [mean_delta, mean_delta], color=color, linewidth=2.0, zorder=4)
        ax_density.text(
            idx,
            y_max - 0.04 * (y_max - y_min),
            f"{mean_delta:+.3f}",
            ha="center",
            va="top",
            fontsize=7.1,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="#fffaf3", edgecolor="#ddd1c2", linewidth=0.5),
        )

    ax_density.set_xticks(xticks)
    ax_density.set_xticklabels(xlabels)
    ax_density.text(
        0.5,
        0.01,
        "Top label = mean delta",
        transform=ax_density.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.8,
        color=MUTED,
    )

    svg_path = FIG_DIR / "fig_office_topology_story_with_delta_sidecar.svg"
    png_path = FIG_DIR / "fig_office_topology_story_with_delta_sidecar_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_story_with_value_sidecar_figure(
    field_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
) -> None:
    story_png = FIG_DIR / "fig_office_topology_story_hires.png"
    if not story_png.exists():
        return

    fig = plt.figure(figsize=(14.8, 7.2), facecolor="white")
    grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.03,
        right=0.98,
        top=0.98,
        bottom=0.08,
        width_ratios=[3.4, 1.20],
        wspace=0.08,
    )

    ax_story = fig.add_subplot(grid[0, 0])
    image = plt.imread(story_png)
    ax_story.imshow(image)
    ax_story.axis("off")

    ax_violin = fig.add_subplot(grid[0, 1])
    ax_violin.set_title("Actual PMV span\ncellular vs open", fontsize=10.8, color=TEXT, pad=10)

    cohort_keys = [
        ("default_male_35", "Default\nMale"),
        ("realistic_male", "Realistic\nMale"),
        ("realistic_mixed", "Realistic\nMixed"),
        ("female_light", "Female\nLight"),
    ]
    all_values: List[float] = []
    values_by_key: Dict[Tuple[str, str], List[float]] = {}
    for cohort_key, _ in cohort_keys:
        for topology in ("cellular_office", "open_office"):
            entry = field_lookup[(topology, cohort_key, FOCUS_TA, FOCUS_RH)]
            values = [as_float(cell["pmv"]) for cell in entry["cells"] if bool(cell["is_work_cell"])]
            values_by_key[(topology, cohort_key)] = values
            all_values.extend(values)

    y_min = math.floor((min(all_values) - 0.05) * 10.0) / 10.0
    y_max = math.ceil((max(all_values) + 0.05) * 10.0) / 10.0
    y_grid = np.linspace(y_min, y_max, 320)

    ax_violin.axhspan(-0.5, 0.5, color="#ece4d7", alpha=0.45, zorder=0)
    ax_violin.axhline(0.0, color=RULE, linewidth=0.9, zorder=1)
    ax_violin.grid(axis="y", color="#ece3d7", linewidth=0.6, alpha=0.8)
    ax_violin.set_axisbelow(True)
    ax_violin.set_xlim(0.35, len(cohort_keys) + 0.65)
    ax_violin.set_ylim(y_min, y_max)
    ax_violin.set_ylabel("PMV", fontsize=8.8, color=TEXT)
    ax_violin.tick_params(axis="y", labelsize=8.0, colors=MUTED)
    ax_violin.tick_params(axis="x", labelsize=7.8, colors=MUTED)
    ax_violin.spines["top"].set_visible(False)
    ax_violin.spines["right"].set_visible(False)
    ax_violin.spines["left"].set_color(RULE)
    ax_violin.spines["bottom"].set_color(RULE)

    max_density = 0.0
    density_lookup: Dict[Tuple[str, str], np.ndarray] = {}
    for cohort_key, _ in cohort_keys:
        for topology in ("cellular_office", "open_office"):
            density = kde_curve(values_by_key[(topology, cohort_key)], y_grid)
            density_lookup[(topology, cohort_key)] = density
            max_density = max(max_density, float(density.max()))

    half_width = 0.26 / max_density if max_density else 0.0
    for idx, (cohort_key, short_label) in enumerate(cohort_keys, start=1):
        cell_values = values_by_key[("cellular_office", cohort_key)]
        open_values = values_by_key[("open_office", cohort_key)]
        cell_density = density_lookup[("cellular_office", cohort_key)] * half_width
        open_density = density_lookup[("open_office", cohort_key)] * half_width

        ax_violin.fill_betweenx(
            y_grid,
            idx - cell_density,
            idx,
            facecolor=CELLULAR_TOPOLOGY,
            alpha=0.24,
            linewidth=0.0,
            zorder=2,
        )
        ax_violin.fill_betweenx(
            y_grid,
            idx,
            idx + open_density,
            facecolor=OPEN_TOPOLOGY,
            alpha=0.24,
            linewidth=0.0,
            zorder=2,
        )
        ax_violin.plot(idx - cell_density, y_grid, color=CELLULAR_TOPOLOGY, linewidth=1.5, zorder=3)
        ax_violin.plot(idx + open_density, y_grid, color=OPEN_TOPOLOGY, linewidth=1.5, zorder=3)
        ax_violin.plot([idx - 0.17, idx - 0.01], [statistics.mean(cell_values)] * 2, color=CELLULAR_TOPOLOGY, linewidth=2.0, zorder=4)
        ax_violin.plot([idx + 0.01, idx + 0.17], [statistics.mean(open_values)] * 2, color=OPEN_TOPOLOGY, linewidth=2.0, zorder=4)
        ax_violin.text(
            idx,
            y_max - 0.04 * (y_max - y_min),
            f"{statistics.mean(cell_values):+.2f} | {statistics.mean(open_values):+.2f}",
            ha="center",
            va="top",
            fontsize=7.0,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.16", facecolor="#fffaf3", edgecolor="#ddd1c2", linewidth=0.5),
        )

    ax_violin.set_xticks(range(1, len(cohort_keys) + 1))
    ax_violin.set_xticklabels([label for _, label in cohort_keys])
    ax_violin.text(
        0.02,
        0.02,
        "Left half = cellular\nRight half = open",
        transform=ax_violin.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color=MUTED,
    )

    svg_path = FIG_DIR / "fig_office_topology_story_with_value_sidecar.svg"
    png_path = FIG_DIR / "fig_office_topology_story_with_value_sidecar_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_story_with_box_sidecar_figure(
    field_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
) -> None:
    story_png = FIG_DIR / "fig_office_topology_story_hires.png"
    if not story_png.exists():
        return

    fig = plt.figure(figsize=(14.8, 7.2), facecolor="white")
    grid = fig.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.03,
        right=0.98,
        top=0.98,
        bottom=0.08,
        width_ratios=[3.4, 1.28],
        wspace=0.08,
    )

    ax_story = fig.add_subplot(grid[0, 0])
    image = plt.imread(story_png)
    ax_story.imshow(image)
    ax_story.axis("off")

    ax_box = fig.add_subplot(grid[0, 1])
    ax_box.set_title("Actual PMV spread\ncellular vs open", fontsize=10.8, color=TEXT, pad=10)

    cohort_keys = [
        ("default_male_35", "Default\nMale"),
        ("realistic_male", "Realistic\nMale"),
        ("realistic_mixed", "Realistic\nMixed"),
        ("female_light", "Female\nLight"),
    ]
    all_values: List[float] = []
    groups: List[List[float]] = []
    positions: List[float] = []
    cell_pos: List[float] = []
    open_pos: List[float] = []
    for idx, (cohort_key, _) in enumerate(cohort_keys, start=1):
        cell_values = [
            as_float(cell["pmv"])
            for cell in field_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]["cells"]
            if bool(cell["is_work_cell"])
        ]
        open_values = [
            as_float(cell["pmv"])
            for cell in field_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]["cells"]
            if bool(cell["is_work_cell"])
        ]
        groups.extend([cell_values, open_values])
        cell_pos.append(idx - 0.16)
        open_pos.append(idx + 0.16)
        positions.extend([idx - 0.16, idx + 0.16])
        all_values.extend(cell_values)
        all_values.extend(open_values)

    y_min = math.floor((min(all_values) - 0.05) * 10.0) / 10.0
    y_max = math.ceil((max(all_values) + 0.05) * 10.0) / 10.0

    ax_box.axhspan(-0.5, 0.5, color="#ece4d7", alpha=0.45, zorder=0)
    ax_box.axhline(0.0, color=RULE, linewidth=0.9, zorder=1)
    ax_box.grid(axis="y", color="#ece3d7", linewidth=0.6, alpha=0.8)
    ax_box.set_axisbelow(True)
    ax_box.set_xlim(0.4, len(cohort_keys) + 0.6)
    ax_box.set_ylim(y_min, y_max)
    ax_box.set_ylabel("PMV", fontsize=8.8, color=TEXT)
    ax_box.tick_params(axis="y", labelsize=8.0, colors=MUTED)
    ax_box.tick_params(axis="x", labelsize=7.8, colors=MUTED)
    ax_box.spines["top"].set_visible(False)
    ax_box.spines["right"].set_visible(False)
    ax_box.spines["left"].set_color(RULE)
    ax_box.spines["bottom"].set_color(RULE)

    bp = ax_box.boxplot(
        groups,
        positions=positions,
        widths=0.22,
        patch_artist=True,
        showfliers=False,
        whis=(5, 95),
        medianprops=dict(color=TEXT, linewidth=1.5),
        whiskerprops=dict(color=RULE, linewidth=1.0),
        capprops=dict(color=RULE, linewidth=1.0),
        boxprops=dict(linewidth=1.1),
    )

    for i, box in enumerate(bp["boxes"]):
        is_cell = i % 2 == 0
        color = CELLULAR_TOPOLOGY if is_cell else OPEN_TOPOLOGY
        color = CELLULAR_COMPARE if is_cell else OPEN_COMPARE
        box.set_facecolor(color)
        box.set_alpha(0.26)
        box.set_edgecolor(color)

    for i, median in enumerate(bp["medians"]):
        median.set_color(CELLULAR_COMPARE if i % 2 == 0 else OPEN_COMPARE)
        median.set_linewidth(1.9)

    # Draw full-range thin lines behind the box to make tail differences explicit.
    for idx, values in enumerate(groups):
        x = positions[idx]
        color = CELLULAR_COMPARE if idx % 2 == 0 else OPEN_COMPARE
        ax_box.plot([x, x], [min(values), max(values)], color=color, linewidth=0.8, alpha=0.55, zorder=1.5)

    ax_box.set_xticks(range(1, len(cohort_keys) + 1))
    ax_box.set_xticklabels([label for _, label in cohort_keys])

    for idx, (cohort_key, _) in enumerate(cohort_keys, start=1):
        cell_values = groups[(idx - 1) * 2]
        open_values = groups[(idx - 1) * 2 + 1]
        ax_box.text(
            idx,
            y_max - 0.045 * (y_max - y_min),
            f"{(max(cell_values)-min(cell_values)):.2f} | {(max(open_values)-min(open_values)):.2f}",
            ha="center",
            va="top",
            fontsize=7.0,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.16", facecolor="#fffaf3", edgecolor="#ddd1c2", linewidth=0.5),
        )

    ax_box.text(
        0.02,
        0.02,
        "Graphite = cellular\nMoss = open\nTop tag = full span C | O",
        transform=ax_box.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color=MUTED,
    )

    svg_path = FIG_DIR / "fig_office_topology_story_with_box_sidecar.svg"
    png_path = FIG_DIR / "fig_office_topology_story_with_box_sidecar_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_story_figure(
    field_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
    summary_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
    cellular_payload: Dict[str, object],
    open_payload: Dict[str, object],
) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    pmv_cmap = make_pmv_cmap()
    pmv_norm = mcolors.Normalize(vmin=-3.0, vmax=0.1)

    fig = plt.figure(figsize=(10.8, 7.2), facecolor="white")
    grid = fig.add_gridspec(
        nrows=2,
        ncols=3,
        left=0.08,
        right=0.97,
        top=0.72,
        bottom=0.17,
        hspace=0.18,
        wspace=0.10,
    )

    fig.text(
        0.08,
        0.965,
        "Office Counterfactual at 22 C / 80% RH",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.92,
        "Same footprint, two topologies, and three representative cohorts.\n"
        "Default-male workers stay least cold, while mixed and female-light cohorts shift colder.\n"
        "Opening the plan barely changes the cohort mean but stretches the work-cell field.",
        fontsize=10.2,
        color=MUTED,
        ha="left",
        va="top",
    )

    note_cell_summary = summary_lookup[("cellular_office", "realistic_male", FOCUS_TA, FOCUS_RH)]
    note_open_summary = summary_lookup[("open_office", "realistic_male", FOCUS_TA, FOCUS_RH)]
    note_cell_field = field_lookup[("cellular_office", "realistic_male", FOCUS_TA, FOCUS_RH)]
    note_open_field = field_lookup[("open_office", "realistic_male", FOCUS_TA, FOCUS_RH)]
    fig.text(
        0.08,
        0.84,
        (
            "Realistic male closely tracks the default baseline: mean PMV "
            f"{as_float(note_cell_summary['mean_pmv']):.2f}\u2192{as_float(note_open_summary['mean_pmv']):.2f}, "
            f"field range {field_range(note_cell_field):.3f}\u2192{field_range(note_open_field):.3f}."
        ),
        fontsize=9.4,
        color=MUTED,
        ha="left",
        va="top",
    )

    cax = fig.add_axes([0.31, 0.09, 0.38, 0.022])
    sm = plt.cm.ScalarMappable(norm=pmv_norm, cmap=pmv_cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#a29281")
    cbar.ax.tick_params(labelsize=8.5, colors=MUTED)
    cbar.set_label("PMV", fontsize=10, color=TEXT)

    for col_index, (cohort_key, cohort_label) in enumerate(DISPLAY_COHORTS):
        cell_entry = field_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        open_entry = field_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        cell_summary = summary_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        open_summary = summary_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]

        top_ax = fig.add_subplot(grid[0, col_index])
        bottom_ax = fig.add_subplot(grid[1, col_index])
        draw_field_panel(
            top_ax,
            cellular_payload,
            cell_entry,
            cmap=pmv_cmap,
            norm=pmv_norm,
        )
        draw_field_panel(
            bottom_ax,
            open_payload,
            open_entry,
            cmap=pmv_cmap,
            norm=pmv_norm,
        )

        header = (
            f"{cohort_label}\n"
            f"mean PMV {as_float(cell_summary['mean_pmv']):.2f}\u2192{as_float(open_summary['mean_pmv']):.2f}   "
            f"mean PPD {as_float(cell_summary['mean_ppd']):.0f}%\u2192{as_float(open_summary['mean_ppd']):.0f}%\n"
            f"field range {field_range(cell_entry):.3f}\u2192{field_range(open_entry):.3f}"
        )
        top_ax.set_title(header, fontsize=9.2, color=TEXT, pad=12, loc="center")

    fig.text(0.035, 0.565, "Cellular Office", rotation=90, va="center", ha="center", fontsize=10.5, color=TEXT, fontweight="bold")
    fig.text(0.035, 0.285, "Open Office", rotation=90, va="center", ha="center", fontsize=10.5, color=TEXT, fontweight="bold")

    svg_path = FIG_DIR / "fig_office_topology_story.svg"
    png_path = FIG_DIR / "fig_office_topology_story_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def render_delta_figure(
    field_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
    open_payload: Dict[str, object],
) -> None:
    delta_cmap = make_delta_cmap()
    delta_norm = mcolors.TwoSlopeNorm(vmin=-0.30, vcenter=0.0, vmax=0.70)

    fig = plt.figure(figsize=(10.8, 4.5), facecolor="white")
    grid = fig.add_gridspec(
        nrows=1,
        ncols=3,
        left=0.08,
        right=0.97,
        top=0.66,
        bottom=0.20,
        wspace=0.10,
    )

    fig.text(
        0.08,
        0.96,
        "Open Office Minus Cellular Office PMV",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.89,
        "Positive values mean the open-plan reinterpretation reads warmer or less cold\n"
        "than the cellular office at the same location. The strongest gains sit along the\n"
        "perimeter and window-facing work bands rather than in the core.",
        fontsize=10.2,
        color=MUTED,
        ha="left",
        va="top",
    )

    cax = fig.add_axes([0.33, 0.10, 0.34, 0.024])
    sm = plt.cm.ScalarMappable(norm=delta_norm, cmap=delta_cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#a29281")
    cbar.ax.tick_params(labelsize=8.5, colors=MUTED)
    cbar.set_label("Open minus cellular PMV", fontsize=10, color=TEXT)

    payload_cells = plan_cells_from_payload(open_payload)
    for col_index, (cohort_key, cohort_label) in enumerate(DISPLAY_COHORTS):
        cell_entry = field_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        open_entry = field_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        cell_cells = cell_map(cell_entry)
        open_cells = cell_map(open_entry)

        ax = fig.add_subplot(grid[0, col_index])
        cols = int(open_payload["meta"]["grid_cols"])
        rows = int(open_payload["meta"]["grid_rows"])
        ax.set_xlim(0, cols)
        ax.set_ylim(rows, 0)
        ax.set_aspect("equal")
        ax.axis("off")

        for coord in payload_cells:
            open_cell = open_cells.get(coord)
            cell_cell = cell_cells.get(coord)
            if open_cell is None:
                facecolor = "#fbf7f0"
            elif not bool(open_cell["is_work_cell"]) and not bool(cell_cell and cell_cell["is_work_cell"]):
                facecolor = WORK_FILL
            elif cell_cell is None:
                facecolor = "#fbf7f0"
            else:
                delta = as_float(open_cell["pmv"]) - as_float(cell_cell["pmv"])
                facecolor = delta_cmap(delta_norm(delta))
            ax.add_patch(
                Rectangle(
                    coord,
                    1.0,
                    1.0,
                    facecolor=facecolor,
                    edgecolor=GRID_EDGE,
                    linewidth=0.35,
                )
            )

        lines = LineCollection(boundary_segments(open_cells), colors=BOUNDARY, linewidths=0.8, zorder=5)
        ax.add_collection(lines)

        overlaps = []
        for coord, open_cell in open_cells.items():
            cell_cell = cell_cells.get(coord)
            if cell_cell is None or not bool(open_cell["is_work_cell"]) or not bool(cell_cell["is_work_cell"]):
                continue
            overlaps.append(as_float(open_cell["pmv"]) - as_float(cell_cell["pmv"]))

        ax.set_title(
            f"{cohort_label}\nmean \u0394 {sum(overlaps) / len(overlaps):+.3f}   max \u0394 {max(overlaps):+.3f}",
            fontsize=9.5,
            color=TEXT,
            pad=12,
        )

    svg_path = FIG_DIR / "fig_office_topology_delta.svg"
    png_path = FIG_DIR / "fig_office_topology_delta_hires.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


def main() -> None:
    field_maps = read_json(FIELD_MAPS_PATH)
    summary_rows = read_csv(SUMMARY_PATH)
    sample_rows = read_csv(SAMPLES_PATH)
    field_lookup = build_field_lookup(field_maps["maps"])
    summary_lookup = build_summary_lookup(summary_rows)
    cellular_payload = read_json(CELLULAR_TEMPLATE)
    open_payload = read_json(OPEN_TEMPLATE)
    records_lookup: Dict[Tuple[str, str, float, float], List[Dict[str, str]]] = {}
    for row in sample_rows:
        key = (
            row["topology"],
            row["cohort"],
            as_float(row["base_ta_c"]),
            as_float(row["base_rh_pct"]),
        )
        records_lookup.setdefault(key, []).append(row)

    render_delta_distribution_figure(cellular_payload)
    render_distribution_story_figure(records_lookup, summary_lookup, cellular_payload, open_payload)
    render_story_figure(field_lookup, summary_lookup, cellular_payload, open_payload)
    render_story_with_density_sidecar_figure()
    render_story_with_value_sidecar_figure(field_lookup)
    render_story_with_box_sidecar_figure(field_lookup)
    render_delta_figure(field_lookup, open_payload)


if __name__ == "__main__":
    main()
