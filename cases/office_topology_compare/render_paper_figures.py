#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CASE_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "thermbal_mplconfig"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

from build_templates import CELLULAR_TEMPLATE, OPEN_TEMPLATE


OUT_DIR = CASE_DIR / "out"
FIG_DIR = OUT_DIR / "figures"
FIELD_MAPS_PATH = OUT_DIR / "field_maps.json"
SUMMARY_PATH = OUT_DIR / "summary.csv"

FOCUS_TA = 22.0
FOCUS_RH = 80.0
DISPLAY_COHORTS = [
    ("default_male_35", "Default Male 35"),
    ("realistic_mixed", "Realistic Mixed"),
    ("female_light", "Female Light"),
]
NOTE_COHORT = ("realistic_male", "Realistic Male")

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


def boundary_segments(cells_by_coord: Dict[Tuple[int, int], Dict[str, object]]) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    segments = set()
    for (col, row), cell in cells_by_coord.items():
        room_name = str(cell["room_name"])
        edges = [
            ((col, row), (col + 1, row), (col, row - 1)),
            ((col + 1, row), (col + 1, row + 1), (col + 1, row)),
            ((col, row + 1), (col + 1, row + 1), (col, row + 1)),
            ((col, row), (col, row + 1), (col - 1, row)),
        ]
        for start, end, neighbor_coord in edges:
            neighbor = cells_by_coord.get(neighbor_coord)
            if neighbor is None or str(neighbor["room_name"]) != room_name:
                segments.add((start, end))
    return sorted(segments)


def plan_cells_from_payload(payload: Dict[str, object]) -> Dict[Tuple[int, int], str]:
    out: Dict[Tuple[int, int], str] = {}
    for room in payload["rooms"]:
        room_name = str(room["name"])
        for col, row in room["cells"]:
            out[(int(col), int(row))] = room_name
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

    note_cell_summary = summary_lookup[("cellular_office", NOTE_COHORT[0], FOCUS_TA, FOCUS_RH)]
    note_open_summary = summary_lookup[("open_office", NOTE_COHORT[0], FOCUS_TA, FOCUS_RH)]
    note_cell_field = field_lookup[("cellular_office", NOTE_COHORT[0], FOCUS_TA, FOCUS_RH)]
    note_open_field = field_lookup[("open_office", NOTE_COHORT[0], FOCUS_TA, FOCUS_RH)]
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

    fig.text(
        0.08,
        0.035,
        "Work cells are colored by mean cohort PMV. Support rooms stay muted. The figure is intentionally not a comfort-average chart: the point is that the open office preserves a longer spatial cold gradient even when the cohort mean barely changes.",
        fontsize=9.2,
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
    field_lookup = build_field_lookup(field_maps["maps"])
    summary_lookup = build_summary_lookup(summary_rows)
    cellular_payload = read_json(CELLULAR_TEMPLATE)
    open_payload = read_json(OPEN_TEMPLATE)

    render_story_figure(field_lookup, summary_lookup, cellular_payload, open_payload)
    render_delta_figure(field_lookup, open_payload)


if __name__ == "__main__":
    main()
