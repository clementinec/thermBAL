#!/usr/bin/env python3
from __future__ import annotations

import csv
import html
import json
import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from build_templates import CELLULAR_TEMPLATE, OPEN_TEMPLATE
from generate_study import CELLULAR_WORK_ROOMS, OPEN_WORK_ROOMS


CASE_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = CASE_DIR / "study_manifest.json"
SAMPLES_PATH = CASE_DIR / "out" / "sample_results.csv"
COMPARISON_PATH = CASE_DIR / "out" / "comparison_summary.csv"
FIELD_MAPS_PATH = CASE_DIR / "out" / "field_maps.json"
REPORT_PATH = CASE_DIR / "out" / "report.html"

CELL_PX = 11


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


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def color_for_pmv(pmv: float | None) -> str:
    if pmv is None:
        return "#ece5d9"
    if abs(pmv) < 0.15:
        return "#e6ddcf"
    severity = clamp(abs(pmv) / 3.0, 0.0, 1.0)
    if pmv < 0:
        hue = 208.0
        sat = 24.0 + 56.0 * severity
        light = 95.0 - 50.0 * severity
    else:
        hue = 18.0
        sat = 28.0 + 58.0 * severity
        light = 95.0 - 50.0 * severity
    return "hsl({0:.0f} {1:.0f}% {2:.0f}%)".format(hue, sat, light)


def color_for_ppd(ppd: float | None) -> str:
    if ppd is None:
        return "#ece5d9"
    severity = clamp(ppd / 100.0, 0.0, 1.0)
    hue = 24.0
    sat = 18.0 + 58.0 * severity
    light = 96.0 - 52.0 * severity
    return "hsl({0:.0f} {1:.0f}% {2:.0f}%)".format(hue, sat, light)


def gender_stroke(gender: str) -> str:
    return "#2d2926" if gender == "male" else "#5f7f99"


def room_centroid(cells: Iterable[Tuple[int, int]]) -> Tuple[float, float]:
    cells = list(cells)
    return (
        sum(col + 0.5 for col, _ in cells) / len(cells),
        sum(row + 0.5 for _, row in cells) / len(cells),
    )


def room_fill_color(name: str) -> str:
    if "Open Office" in name or "Office" in name or "Studio" in name or "Bay" in name or "Corner" in name:
        return "#dce8f1"
    if "Collaboration" in name or "Meeting" in name:
        return "#ebdcc2"
    if "Reception" in name or "Admin" in name or "Focus" in name:
        return "#e6eee4"
    return "#eee6db"


def work_rooms_for_payload(payload: Dict[str, object]) -> set[str]:
    room_names = {str(room["name"]) for room in payload["rooms"]}
    if "North Open Office" in room_names:
        return set(OPEN_WORK_ROOMS)
    return set(CELLULAR_WORK_ROOMS)


def room_index(payload: Dict[str, object]) -> Dict[str, set[Tuple[int, int]]]:
    return {
        str(room["name"]): {tuple(cell) for cell in room["cells"]}
        for room in payload["rooms"]
    }


def draw_plan_background(
    fragments: List[str],
    payload: Dict[str, object],
    *,
    work_fill: str = "#f7f0e4",
    support_fill: str = "#ece3d7",
) -> None:
    rooms = payload["rooms"]
    room_cells = room_index(payload)
    work_rooms = work_rooms_for_payload(payload)

    for room in rooms:
        name = str(room["name"])
        fill = work_fill if name in work_rooms else support_fill
        for col, row in room_cells[name]:
            fragments.append(
                '<rect x="{0}" y="{1}" width="{2}" height="{2}" fill="{3}" stroke="#ebe1d3" stroke-width="0.85"/>'.format(
                    col * CELL_PX,
                    row * CELL_PX,
                    CELL_PX,
                    fill,
                )
            )

    for room in rooms:
        name = str(room["name"])
        cells = room_cells[name]
        for col, row in cells:
            if (col, row - 1) not in cells:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#7a7064" stroke-width="1.2"/>'.format(
                        col * CELL_PX,
                        row * CELL_PX,
                        (col + 1) * CELL_PX,
                    )
                )
            if (col + 1, row) not in cells:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#7a7064" stroke-width="1.2"/>'.format(
                        (col + 1) * CELL_PX,
                        row * CELL_PX,
                        (row + 1) * CELL_PX,
                    )
                )
            if (col, row + 1) not in cells:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#7a7064" stroke-width="1.2"/>'.format(
                        col * CELL_PX,
                        (row + 1) * CELL_PX,
                        (col + 1) * CELL_PX,
                    )
                )
            if (col - 1, row) not in cells:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#7a7064" stroke-width="1.2"/>'.format(
                        col * CELL_PX,
                        row * CELL_PX,
                        (row + 1) * CELL_PX,
                    )
                )


def render_topology_svg(payload: Dict[str, object], abbreviations: Dict[str, str]) -> str:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    width = cols * CELL_PX
    height = rows * CELL_PX

    fragments = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Topology diagram">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fbf7f0" rx="10" ry="10"/>',
    ]
    draw_plan_background(fragments, payload, work_fill="#f3eee4", support_fill="#efe7db")

    for room in payload["rooms"]:
        name = str(room["name"])
        cells = [tuple(cell) for cell in room["cells"]]
        for col, row in cells:
            fragments.append(
                '<rect x="{0}" y="{1}" width="{2}" height="{2}" fill="{3}" fill-opacity="0.55" stroke="none"/>'.format(
                    col * CELL_PX,
                    row * CELL_PX,
                    CELL_PX,
                    room_fill_color(name),
                )
            )
        cx, cy = room_centroid(cells)
        fragments.append(
            '<text x="{0}" y="{1}" text-anchor="middle" dominant-baseline="middle" font-size="9" font-family="Avenir Next, Helvetica Neue, sans-serif" fill="#2d2926">{2}</text>'.format(
                cx * CELL_PX,
                cy * CELL_PX,
                html.escape(abbreviations.get(name, name)),
            )
        )

    fragments.append("</svg>")
    return "".join(fragments)


def render_sample_scatter_svg(payload: Dict[str, object], records: List[Dict[str, str]]) -> str:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    width = cols * CELL_PX
    height = rows * CELL_PX
    fragments = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Sampled PMV scatter">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fbf7f0" rx="10" ry="10"/>',
    ]
    draw_plan_background(fragments, payload)

    for record in records:
        x = (as_float(record["x"]) + 0.5) * CELL_PX
        y = (as_float(record["y"]) + 0.5) * CELL_PX
        pmv = as_float(record["pmv"])
        ppd = as_float(record["ppd"])
        radius = 2.0 + 3.2 * clamp(ppd / 100.0, 0.0, 1.0)
        fragments.append("<g>")
        fragments.append(
            "<title>{0}</title>".format(
                html.escape(
                    "{0} · {1} · age {2} · PMV {3:.2f} · PPD {4:.1f}%".format(
                        record["sample_id"],
                        record["gender"],
                        int(float(record["age"])),
                        pmv,
                        ppd,
                    )
                )
            )
        )
        fragments.append(
            '<circle cx="{0}" cy="{1}" r="{2}" fill="{3}" stroke="{4}" stroke-width="1.4" fill-opacity="0.95"/>'.format(
                x,
                y,
                radius,
                color_for_pmv(pmv),
                gender_stroke(str(record["gender"])),
            )
        )
        fragments.append("</g>")

    fragments.append("</svg>")
    return "".join(fragments)


def interpolate_metric(records: List[Dict[str, str]], payload: Dict[str, object], metric_key: str) -> Dict[Tuple[int, int], float]:
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


def render_interpolated_svg(payload: Dict[str, object], metric_map: Dict[Tuple[int, int], float], metric: str) -> str:
    cols = int(payload["meta"]["grid_cols"])
    rows = int(payload["meta"]["grid_rows"])
    width = cols * CELL_PX
    height = rows * CELL_PX
    fragments = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Interpolated metric field">',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#fbf7f0" rx="10" ry="10"/>',
    ]
    draw_plan_background(fragments, payload, work_fill="#f2ebde", support_fill="#ebe2d6")

    room_cells = room_index(payload)
    work_rooms = work_rooms_for_payload(payload)
    for room in payload["rooms"]:
        room_name = str(room["name"])
        if room_name not in work_rooms:
            continue
        for cell in room_cells[room_name]:
            value = metric_map.get(cell)
            color = color_for_pmv(value) if metric == "pmv" else color_for_ppd(value)
            fragments.append(
                '<rect x="{0}" y="{1}" width="{2}" height="{2}" fill="{3}" fill-opacity="0.96" stroke="#e8ded0" stroke-width="0.6"/>'.format(
                    cell[0] * CELL_PX,
                    cell[1] * CELL_PX,
                    CELL_PX,
                    color,
                )
            )

    fragments.append("</svg>")
    return "".join(fragments)


def metric_summary(records: List[Dict[str, str]]) -> Dict[str, float]:
    pmvs = [as_float(record["pmv"]) for record in records]
    ppds = [as_float(record["ppd"]) for record in records]
    return {
        "mean_pmv": statistics.mean(pmvs),
        "worst_ppd": max(ppds),
        "mean_ppd": statistics.mean(ppds),
    }


def focus_section(
    label: str,
    cohort_key: str,
    records_by_key: Dict[Tuple[str, str, float, float], List[Dict[str, str]]],
    comparison_lookup: Dict[Tuple[str, float, float], Dict[str, str]],
    cellular_payload: Dict[str, object],
    open_payload: Dict[str, object],
    *,
    base_ta_c: float,
    base_rh_pct: float,
) -> str:
    cell_records = records_by_key[("cellular_office", cohort_key, base_ta_c, base_rh_pct)]
    open_records = records_by_key[("open_office", cohort_key, base_ta_c, base_rh_pct)]

    demo_counts = Counter(record["gender"] for record in open_records)
    mean_age = statistics.mean(as_float(record["age"]) for record in open_records)
    comparison = comparison_lookup[(cohort_key, base_ta_c, base_rh_pct)]

    cell_interp_pmv = interpolate_metric(cell_records, cellular_payload, "pmv")
    open_interp_pmv = interpolate_metric(open_records, open_payload, "pmv")
    cell_interp_ppd = interpolate_metric(cell_records, cellular_payload, "ppd")
    open_interp_ppd = interpolate_metric(open_records, open_payload, "ppd")

    cell_pmv_values = list(cell_interp_pmv.values())
    open_pmv_values = list(open_interp_pmv.values())

    return """
    <section class="focus-section">
      <div class="section-copy">
        <h3>{label}</h3>
        <p><strong>{female}</strong> female and <strong>{male}</strong> male sampled workers, with mean age <strong>{age:.1f}</strong>. At <strong>{ta:.0f} C / {rh:.0f}% RH</strong>, the open office increases sampled PMV range by <strong>{range_delta:+.3f}</strong> and perimeter-core gap by <strong>{gap_delta:+.3f}</strong>.</p>
        <div class="chip-row">
          <span class="chip">circle fill = PMV</span>
          <span class="chip">circle radius = PPD</span>
          <span class="chip">male stroke = charcoal</span>
          <span class="chip">female stroke = blue-grey</span>
        </div>
      </div>
      <div class="plot-pair">
        <article class="plot-card">
          <h4>Cellular Office</h4>
          <div class="metric-strip">
            <span>mean PMV {cell_mean_pmv:.2f}</span>
            <span>mean PPD {cell_mean_ppd:.1f}%</span>
            <span>worst PPD {cell_worst_ppd:.1f}%</span>
            <span>interpolated range {cell_interp_range:.3f}</span>
          </div>
          <div class="figure-block">
            <div class="figure-card">
              <h5>Sampled PMV / PPD Points</h5>
              {cell_scatter}
            </div>
            <div class="mini-grid">
              <div class="figure-card">
                <h5>Interpolated PMV</h5>
                {cell_pmv_map}
              </div>
              <div class="figure-card">
                <h5>Interpolated PPD</h5>
                {cell_ppd_map}
              </div>
            </div>
          </div>
        </article>
        <article class="plot-card">
          <h4>Open Office</h4>
          <div class="metric-strip">
            <span>mean PMV {open_mean_pmv:.2f}</span>
            <span>mean PPD {open_mean_ppd:.1f}%</span>
            <span>worst PPD {open_worst_ppd:.1f}%</span>
            <span>interpolated range {open_interp_range:.3f}</span>
          </div>
          <div class="figure-block">
            <div class="figure-card">
              <h5>Sampled PMV / PPD Points</h5>
              {open_scatter}
            </div>
            <div class="mini-grid">
              <div class="figure-card">
                <h5>Interpolated PMV</h5>
                {open_pmv_map}
              </div>
              <div class="figure-card">
                <h5>Interpolated PPD</h5>
                {open_ppd_map}
              </div>
            </div>
          </div>
        </article>
      </div>
    </section>
    """.format(
        label=html.escape(label),
        female=demo_counts.get("female", 0),
        male=demo_counts.get("male", 0),
        age=mean_age,
        ta=base_ta_c,
        rh=base_rh_pct,
        range_delta=as_float(comparison["pmv_range_delta_open_minus_cellular"]),
        gap_delta=as_float(comparison["perimeter_core_gap_delta_open_minus_cellular"]),
        cell_mean_pmv=metric_summary(cell_records)["mean_pmv"],
        cell_mean_ppd=metric_summary(cell_records)["mean_ppd"],
        cell_worst_ppd=metric_summary(cell_records)["worst_ppd"],
        cell_interp_range=max(cell_pmv_values) - min(cell_pmv_values),
        cell_scatter=render_sample_scatter_svg(cellular_payload, cell_records),
        cell_pmv_map=render_interpolated_svg(cellular_payload, cell_interp_pmv, "pmv"),
        cell_ppd_map=render_interpolated_svg(cellular_payload, cell_interp_ppd, "ppd"),
        open_mean_pmv=metric_summary(open_records)["mean_pmv"],
        open_mean_ppd=metric_summary(open_records)["mean_ppd"],
        open_worst_ppd=metric_summary(open_records)["worst_ppd"],
        open_interp_range=max(open_pmv_values) - min(open_pmv_values),
        open_scatter=render_sample_scatter_svg(open_payload, open_records),
        open_pmv_map=render_interpolated_svg(open_payload, open_interp_pmv, "pmv"),
        open_ppd_map=render_interpolated_svg(open_payload, open_interp_ppd, "ppd"),
    )


def main() -> None:
    manifest = read_json(MANIFEST_PATH)
    sample_rows = read_csv(SAMPLES_PATH)
    comparison_rows = read_csv(COMPARISON_PATH)
    field_maps = read_json(FIELD_MAPS_PATH)
    cellular_payload = read_json(CELLULAR_TEMPLATE)
    open_payload = read_json(OPEN_TEMPLATE)

    records_by_key: Dict[Tuple[str, str, float, float], List[Dict[str, str]]] = {}
    for row in sample_rows:
        key = (
            row["topology"],
            row["cohort"],
            as_float(row["base_ta_c"]),
            as_float(row["base_rh_pct"]),
        )
        records_by_key.setdefault(key, []).append(row)

    comparison_lookup = {
        (row["cohort"], as_float(row["base_ta_c"]), as_float(row["base_rh_pct"])): row
        for row in comparison_rows
    }

    focus_rows = [
        ("Realistic Mixed / 22 C / 80% RH", "realistic_mixed", 22.0, 80.0),
        ("Default Male 35 / 22 C / 80% RH", "default_male_35", 22.0, 80.0),
        ("Realistic Male / 22 C / 80% RH", "realistic_male", 22.0, 80.0),
        ("Female Light / 22 C / 80% RH", "female_light", 22.0, 80.0),
    ]

    realistic_open = records_by_key[("open_office", "realistic_mixed", 22.0, 80.0)]
    realistic_counts = Counter(record["gender"] for record in realistic_open)
    realistic_age = statistics.mean(as_float(record["age"]) for record in realistic_open)

    range_gains = [as_float(row["pmv_range_delta_open_minus_cellular"]) for row in comparison_rows]
    mean_range_gain = statistics.mean(range_gains)
    max_range_gain = max(range_gains)

    mrt_ref = next(
        entry
        for entry in field_maps["maps"]
        if entry["topology"] == "open_office"
        and entry["cohort"] == "realistic_mixed"
        and as_float(entry["base_ta_c"]) == 25.0
        and as_float(entry["base_rh_pct"]) == 55.0
    )
    peak_mrt = max(as_float(cell["mean_radiant_temp"]) for cell in mrt_ref["cells"] if cell["is_work_cell"])

    focus_sections = "".join(
        focus_section(
            label,
            cohort_key,
            records_by_key,
            comparison_lookup,
            cellular_payload,
            open_payload,
            base_ta_c=ta,
            base_rh_pct=rh,
        )
        for label, cohort_key, ta, rh in focus_rows
    )

    html_doc = """
    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>thermBAL Office Spatial Comparison</title>
      <style>
        :root {{
          --paper: rgba(255, 252, 247, 0.92);
          --paper-strong: #fffaf2;
          --ink: #2d2926;
          --muted: #7d7367;
          --rule: #d8cab8;
          --accent: #c7875f;
        }}
        * {{ box-sizing: border-box; }}
        body {{
          margin: 0;
          color: var(--ink);
          font-family: "Avenir Next", "Helvetica Neue", Helvetica, Arial, sans-serif;
          background:
            radial-gradient(circle at top left, rgba(214, 194, 166, 0.28), transparent 34%),
            linear-gradient(180deg, #f7f1e8 0%, #efe5d8 100%);
        }}
        main {{
          max-width: 1440px;
          margin: 0 auto;
          padding: 38px 26px 56px;
        }}
        h1, h2, h3, h4, h5, p {{ margin: 0; }}
        .hero {{
          display: grid;
          grid-template-columns: 1.15fr 0.85fr;
          gap: 18px;
          margin-bottom: 22px;
        }}
        .hero-copy, .stats, .panel, .section-copy, .plot-card, .footnote {{
          background: var(--paper);
          border: 1px solid var(--rule);
          border-radius: 24px;
        }}
        .hero-copy {{
          padding: 28px;
        }}
        .eyebrow {{
          font-size: 11px;
          letter-spacing: 0.18em;
          text-transform: uppercase;
          color: var(--accent);
          margin-bottom: 10px;
        }}
        .hero-copy h1 {{
          font-size: clamp(34px, 4vw, 56px);
          line-height: 0.95;
          max-width: 12ch;
        }}
        .hero-copy p {{
          margin-top: 14px;
          color: var(--muted);
          line-height: 1.55;
          max-width: 66ch;
        }}
        .stats {{
          padding: 18px;
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 12px;
          align-content: start;
        }}
        .stat-card {{
          padding: 16px;
          border-radius: 18px;
          background: var(--paper-strong);
          border: 1px solid #e2d5c6;
        }}
        .stat-card span {{
          display: block;
          font-size: 12px;
          text-transform: uppercase;
          letter-spacing: 0.08em;
          color: var(--muted);
          margin-bottom: 8px;
        }}
        .stat-card strong {{
          display: block;
          font-size: 28px;
          line-height: 1;
        }}
        .stat-card small {{
          display: block;
          margin-top: 8px;
          line-height: 1.45;
          color: var(--muted);
        }}
        .topologies {{
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 16px;
          margin-bottom: 22px;
        }}
        .panel {{
          padding: 18px;
        }}
        .panel p {{
          margin: 8px 0 14px;
          line-height: 1.5;
          color: var(--muted);
        }}
        .panel svg, .figure-card svg {{
          width: 100%;
          height: auto;
          display: block;
        }}
        .legend, .chip-row {{
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
          margin-top: 12px;
        }}
        .chip {{
          padding: 6px 10px;
          border-radius: 999px;
          border: 1px solid #e3d7c8;
          background: #faf5ed;
          color: var(--muted);
          font-size: 12px;
        }}
        .focus-section {{
          display: grid;
          grid-template-columns: 320px 1fr;
          gap: 16px;
          align-items: start;
          margin-bottom: 20px;
        }}
        .section-copy {{
          padding: 18px;
          position: sticky;
          top: 20px;
        }}
        .section-copy p {{
          margin-top: 10px;
          color: var(--muted);
          line-height: 1.55;
        }}
        .plot-pair {{
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 16px;
        }}
        .plot-card {{
          padding: 16px;
        }}
        .metric-strip {{
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
          margin: 10px 0 14px;
        }}
        .metric-strip span {{
          padding: 6px 10px;
          border-radius: 999px;
          background: #f7f0e5;
          border: 1px solid #e5d8c8;
          font-size: 12px;
          color: var(--muted);
        }}
        .figure-block {{
          display: grid;
          gap: 12px;
        }}
        .mini-grid {{
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 12px;
        }}
        .figure-card {{
          padding: 12px;
          border-radius: 18px;
          border: 1px solid #e6dacb;
          background: var(--paper-strong);
        }}
        .figure-card h5 {{
          margin-bottom: 10px;
          font-size: 13px;
        }}
        .footnote {{
          padding: 18px;
          color: var(--muted);
          line-height: 1.55;
        }}
        @media (max-width: 1180px) {{
          .hero, .topologies, .focus-section, .plot-pair, .mini-grid {{
            grid-template-columns: 1fr;
          }}
          .section-copy {{
            position: static;
          }}
        }}
      </style>
    </head>
    <body>
      <main>
        <section class="hero">
          <div class="hero-copy">
            <div class="eyebrow">Office Counterfactual</div>
            <h1>Spatial PMV and PPD Become Legible When the Plan Opens Up</h1>
            <p>This report focuses on the actual sampled office results rather than aggregate spread tables. The emphasis is the cool-humid stress case at <strong>22 C / 80% RH</strong>, using the same footprint in two topologies: a cellular office and an open office. Sampled PMV/PPD points are shown directly on plan, then interpolated into spatial PMV and PPD fields.</p>
          </div>
          <div class="stats">
            <div class="stat-card">
              <span>Scenario Combos</span>
              <strong>{scenario_count}</strong>
              <small>2 topologies x 6 cohorts x 5 air temperatures x 4 humidity states.</small>
            </div>
            <div class="stat-card">
              <span>Shared Desk Samples</span>
              <strong>{sample_count}</strong>
              <small>Fixed global sample points reused across both topologies.</small>
            </div>
            <div class="stat-card">
              <span>Realistic Mixed Cohort</span>
              <strong>{female_count}F / {male_count}M</strong>
              <small>Mean age {mean_age:.1f}. This now respects the intended 55% female / 45% male split.</small>
            </div>
            <div class="stat-card">
              <span>Peak Window MRT</span>
              <strong>{peak_mrt:.2f} C</strong>
              <small>Open-office work-cell maximum at the `25 C / 55% RH` reference state.</small>
            </div>
            <div class="stat-card">
              <span>Average Open Gain</span>
              <strong>{mean_gain:+.3f}</strong>
              <small>Mean open-minus-cellular sampled PMV range delta across every comparison row.</small>
            </div>
            <div class="stat-card">
              <span>Largest Gain</span>
              <strong>{max_gain:+.3f}</strong>
              <small>Maximum sampled PMV range delta across the office comparison matrix.</small>
            </div>
          </div>
        </section>

        <section class="topologies">
          <article class="panel">
            <h2>Cellular Office Template</h2>
            <p>The apartment footprint is re-labeled as a cellular office with enclosed team rooms, studios, meeting rooms, a retained core, and distributed service spaces.</p>
            {cellular_svg}
            <div class="legend">
              <span class="chip">Blue = enclosed work rooms</span>
              <span class="chip">Green = reception / admin / focus</span>
              <span class="chip">Tan = collaboration / meeting</span>
              <span class="chip">Taupe = support / wet rooms</span>
            </div>
          </article>
          <article class="panel">
            <h2>Open Office Template</h2>
            <p>The same footprint is simplified into two large open-office fields plus a west corner office, while the internal core and wet rooms remain enclosed.</p>
            {open_svg}
            <div class="legend">
              <span class="chip">NO / SO = open office fields</span>
              <span class="chip">WC = west corner office</span>
              <span class="chip">CL = retained core lobby</span>
              <span class="chip">NW / EW / CW = wet rooms</span>
            </div>
          </article>
        </section>

        {focus_sections}

        <section class="footnote">
          The point of the comparison is not that open plan always makes conditions better. It is that open plan makes spatial gradients more legible: the perimeter band and the core separate visibly in the interpolated PMV and PPD fields, while the cellular case flattens those same gradients back into room averages.
        </section>
      </main>
    </body>
    </html>
    """.format(
        scenario_count=int(manifest["scenario_grid"]["scenario_count"]),
        sample_count=len(manifest["sample_points"]),
        female_count=realistic_counts.get("female", 0),
        male_count=realistic_counts.get("male", 0),
        mean_age=realistic_age,
        peak_mrt=peak_mrt,
        mean_gain=mean_range_gain,
        max_gain=max_range_gain,
        cellular_svg=render_topology_svg(
            cellular_payload,
            {
                room["name"]: room["name"].replace("Office", "Ofc").replace("Collaboration", "Collab").replace("Project", "Proj").replace("Connector", "Conn")
                for room in cellular_payload["rooms"]
            },
        ),
        open_svg=render_topology_svg(
            open_payload,
            {
                "North Open Office": "NO",
                "South Open Office": "SO",
                "West Corner Office": "WC",
                "Core Lobby": "CL",
                "North WC": "NW",
                "East WC": "EW",
                "Wellness Room": "WR",
                "Core WC": "CW",
            },
        ),
        focus_sections=focus_sections,
    )

    REPORT_PATH.write_text(html_doc, encoding="utf-8")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
