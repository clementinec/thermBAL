#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import statistics
import tempfile
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

WORK_FILL = "#efe7dc"
GRID_EDGE = "#f7f2ea"
BOUNDARY = "#7e7469"
TEXT = "#2d2926"
MUTED = "#6d645b"


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


def render_story_figure(
    records_lookup: Dict[Tuple[str, str, float, float], List[Dict[str, str]]],
    summary_lookup: Dict[Tuple[str, str, float, float], Dict[str, object]],
    cellular_payload: Dict[str, object],
    open_payload: Dict[str, object],
) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    pmv_cmap = make_pmv_cmap()
    panel_data: List[Dict[str, object]] = []
    for cohort_key, _ in DISPLAY_COHORTS:
        cell_records = records_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        open_records = records_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        cell_interp_pmv = interpolate_metric(cell_records, cellular_payload, "pmv")
        open_interp_pmv = interpolate_metric(open_records, open_payload, "pmv")
        cell_interp_ppd = interpolate_metric(cell_records, cellular_payload, "ppd")
        open_interp_ppd = interpolate_metric(open_records, open_payload, "ppd")
        panel_data.append(
            {
                "cohort_key": cohort_key,
                "cell_pmv": cell_interp_pmv,
                "open_pmv": open_interp_pmv,
                "cell_ppd": cell_interp_ppd,
                "open_ppd": open_interp_ppd,
            }
        )

    display_pmvs = [
        value
        for panel in panel_data
        for mapping in (panel["cell_pmv"], panel["open_pmv"])
        for value in mapping.values()
    ]
    pmv_min = min(display_pmvs)
    pmv_max = max(display_pmvs)
    pmv_norm = mcolors.Normalize(
        vmin=math.floor((pmv_min - 0.05) * 10.0) / 10.0,
        vmax=math.ceil((pmv_max + 0.05) * 10.0) / 10.0,
    )

    fig = plt.figure(figsize=(10.9, 7.8), facecolor="white")
    grid = fig.add_gridspec(
        nrows=2,
        ncols=3,
        left=0.08,
        right=0.97,
        top=0.80,
        bottom=0.16,
        hspace=0.16,
        wspace=0.08,
    )

    fig.text(
        0.08,
        0.965,
        "Office Topology Counterfactual at 22 C / 80% RH",
        fontsize=15.5,
        fontweight="bold",
        color=TEXT,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.925,
        "Interpolated work-area PMV is shown on the shared footprint; support spaces stay muted so the spatial gradient stays legible.",
        fontsize=9.2,
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.text(
        0.08,
        0.900,
        "The open-plan reinterpretation changes cohort means very little, but it preserves a longer and more continuous cold gradient from perimeter to core.",
        fontsize=9.2,
        color=MUTED,
        ha="left",
        va="top",
    )

    cax = fig.add_axes([0.37, 0.08, 0.27, 0.022])
    sm = plt.cm.ScalarMappable(norm=pmv_norm, cmap=pmv_cmap)
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.outline.set_linewidth(0.6)
    cbar.outline.set_edgecolor("#a29281")
    cbar.ax.tick_params(labelsize=8.5, colors=MUTED)
    cbar.set_label("PMV", fontsize=10, color=TEXT)

    for col_index, (cohort_key, cohort_label) in enumerate(DISPLAY_COHORTS):
        panel = next(item for item in panel_data if item["cohort_key"] == cohort_key)
        cell_summary = summary_lookup[("cellular_office", cohort_key, FOCUS_TA, FOCUS_RH)]
        open_summary = summary_lookup[("open_office", cohort_key, FOCUS_TA, FOCUS_RH)]

        top_ax = fig.add_subplot(grid[0, col_index])
        bottom_ax = fig.add_subplot(grid[1, col_index])
        draw_smooth_field_panel(
            top_ax,
            cellular_payload,
            panel["cell_pmv"],
            panel["cell_ppd"],
            cell_summary,
            max(panel["cell_pmv"].values()) - min(panel["cell_pmv"].values()),
            cmap=pmv_cmap,
            norm=pmv_norm,
            row_label="Cellular Office" if col_index == 0 else None,
        )
        draw_smooth_field_panel(
            bottom_ax,
            open_payload,
            panel["open_pmv"],
            panel["open_ppd"],
            open_summary,
            max(panel["open_pmv"].values()) - min(panel["open_pmv"].values()),
            cmap=pmv_cmap,
            norm=pmv_norm,
            row_label="Open Office" if col_index == 0 else None,
        )

        top_ax.set_title(cohort_label, fontsize=10.2, color=TEXT, pad=10, loc="center")

    fig.text(
        0.08,
        0.045,
        "Panel chips list mean PMV, mean PPD, and interpolated work-area span for each topology.",
        fontsize=8.8,
        color=MUTED,
        ha="left",
    )

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

    render_story_figure(records_lookup, summary_lookup, cellular_payload, open_payload)
    render_delta_figure(field_lookup, open_payload)


if __name__ == "__main__":
    main()
