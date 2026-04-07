#!/usr/bin/env python3
from __future__ import annotations

import csv
import html
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Tuple


CASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = CASE_DIR.parents[1]
MANIFEST_PATH = CASE_DIR / "study_manifest.json"
SUMMARY_PATH = CASE_DIR / "out" / "summary.csv"
TIMESERIES_PATH = CASE_DIR / "out" / "agent_timeseries.csv"
REPORT_PATH = CASE_DIR / "out" / "report.html"

CELL_PX = 12


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


def color_for_thermal_state(pmv: float | None, ppd: float | None) -> str:
    if pmv is None:
        return "#f2eee7"

    risk = clamp((ppd or 0.0) / 100.0, 0.0, 1.0)
    magnitude = clamp(abs(pmv) / 3.0, 0.0, 1.0)
    severity = clamp(0.58 * magnitude + 0.42 * risk, 0.0, 1.0)

    if abs(pmv) < 0.18 and risk < 0.12:
        return "#e8e1d4"

    if pmv < 0:
        hue = 208.0
        sat = 24.0 + 52.0 * severity
        light = 95.0 - 50.0 * severity
    else:
        hue = 18.0
        sat = 26.0 + 56.0 * severity
        light = 95.0 - 50.0 * severity
    return "hsl({0:.0f} {1:.0f}% {2:.0f}%)".format(hue, sat, light)


def exposure_color(exposure: str) -> str:
    palette = {
        "south": "#f4c7b0",
        "east_west": "#f0ddb1",
        "north_interior": "#dbe8f6",
        "unoccupied": "#f2eee7",
    }
    return palette.get(exposure, palette["unoccupied"])


def room_centroid(cells: Iterable[Tuple[int, int]]) -> Tuple[float, float]:
    cells = list(cells)
    return (
        sum(col + 0.5 for col, _ in cells) / len(cells),
        sum(row + 0.5 for _, row in cells) / len(cells),
    )


def room_bbox(cells: Iterable[Tuple[int, int]]) -> Tuple[int, int, int, int]:
    cells = list(cells)
    cols = [col for col, _ in cells]
    rows = [row for _, row in cells]
    return min(cols), min(rows), max(cols), max(rows)


def scenario_sort_key(meta: Dict[str, object]) -> Tuple[int, float, float]:
    cohort_order = {
        "young_mixed": 0,
        "young_male": 1,
        "older_shift": 2,
        "higher_clo_sedentary": 3,
        "lighter_clothing_mobile": 4,
        "higher_bmi_warm_sensitive": 5,
    }
    return (
        cohort_order.get(str(meta["cohort"]), 99),
        float(meta["base_ta_c"]),
        float(meta["base_rh_pct"]),
    )


def render_plan_svg(
    rooms: List[Dict[str, object]],
    cols: int,
    rows: int,
    metric_by_room: Dict[str, Dict[str, float]],
    occupied_rooms: set[str],
    abbreviations: Dict[str, str],
    exposure_by_room: Dict[str, str],
    *,
    label_rooms: bool = False,
    exposure_mode: bool = False,
) -> str:
    width = cols * CELL_PX
    height = rows * CELL_PX
    fragments: List[str] = [
        '<svg viewBox="0 0 {0} {1}" role="img" aria-label="Apartment plan">'.format(width, height),
        '<rect x="0" y="0" width="{0}" height="{1}" fill="#fcfaf6" rx="8" ry="8"/>'.format(width, height),
    ]

    room_index = {
        room["name"]: {tuple(cell) for cell in room["cells"]}
        for room in rooms
    }

    for room in rooms:
        name = str(room["name"])
        cells = [tuple(cell) for cell in room["cells"]]
        metric = metric_by_room.get(name, {})
        ppd = metric.get("ppd")
        pmv = metric.get("pmv")
        if exposure_mode:
            fill = exposure_color(exposure_by_room.get(name, "unoccupied"))
        else:
            fill = color_for_thermal_state(
                pmv if name in occupied_rooms else None,
                ppd if name in occupied_rooms else None,
            )

        tooltip_parts = [name]
        if name in occupied_rooms and ppd is not None:
            tooltip_parts.append("PPD {0:.1f}%".format(ppd))
        if name in occupied_rooms and pmv is not None:
            tooltip_parts.append("PMV {0:.2f}".format(pmv))
        if exposure_mode:
            tooltip_parts.append("Exposure {0}".format(exposure_by_room.get(name, "unoccupied")))

        fragments.append("<g>")
        fragments.append("<title>{0}</title>".format(html.escape(" · ".join(tooltip_parts))))
        for col, row in cells:
            fragments.append(
                '<rect x="{0}" y="{1}" width="{2}" height="{2}" fill="{3}" stroke="#e7dfd2" stroke-width="1"/>'.format(
                    col * CELL_PX,
                    row * CELL_PX,
                    CELL_PX,
                    fill,
                )
            )
        for col, row in cells:
            if (col, row - 1) not in room_index[name]:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#7f7365" stroke-width="1.5"/>'.format(
                        col * CELL_PX,
                        row * CELL_PX,
                        (col + 1) * CELL_PX,
                    )
                )
            if (col, row + 1) not in room_index[name]:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#7f7365" stroke-width="1.5"/>'.format(
                        col * CELL_PX,
                        (row + 1) * CELL_PX,
                        (col + 1) * CELL_PX,
                    )
                )
            if (col - 1, row) not in room_index[name]:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#7f7365" stroke-width="1.5"/>'.format(
                        col * CELL_PX,
                        row * CELL_PX,
                        (row + 1) * CELL_PX,
                    )
                )
            if (col + 1, row) not in room_index[name]:
                fragments.append(
                    '<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#7f7365" stroke-width="1.5"/>'.format(
                        (col + 1) * CELL_PX,
                        row * CELL_PX,
                        (row + 1) * CELL_PX,
                    )
                )
        if label_rooms and name in abbreviations:
            cx, cy = room_centroid(cells)
            fragments.append(
                '<text x="{0}" y="{1}" text-anchor="middle" dominant-baseline="middle" '
                'font-size="11" font-family="SF Mono, Fira Code, monospace" fill="#2d2d2d">{2}</text>'.format(
                    cx * CELL_PX,
                    cy * CELL_PX,
                    html.escape(abbreviations[name]),
                )
            )
        fragments.append("</g>")

    fragments.append("</svg>")
    return "".join(fragments)


def render_metric_table(
    cohort_label: str,
    cohort_key: str,
    temperatures: List[float],
    humidities: List[float],
    summary_by_scenario: Dict[str, Dict[str, str]],
    mean_pmv_by_scenario: Dict[str, float],
    mean_ppd_by_scenario: Dict[str, float],
    scenario_meta: Dict[str, Dict[str, object]],
) -> str:
    rows = [
        '<div class="table-card"><h4>{0}</h4><table class="matrix"><thead><tr><th>Ta \\ RH</th>'.format(
            html.escape(cohort_label)
        )
    ]
    for rh in humidities:
        rows.append("<th>{0:.0f}%</th>".format(rh))
    rows.append("</tr></thead><tbody>")

    for ta in temperatures:
        rows.append("<tr><th>{0:.1f} C</th>".format(ta))
        for rh in humidities:
            scenario_id = next(
                key
                for key, meta in scenario_meta.items()
                if meta["cohort"] == cohort_key and float(meta["base_ta_c"]) == ta and float(meta["base_rh_pct"]) == rh
            )
            summary = summary_by_scenario[scenario_id]
            mean_pmv = mean_pmv_by_scenario.get(scenario_id, 0.0)
            mean_ppd = mean_ppd_by_scenario.get(scenario_id, 0.0)
            worst_ppd = as_float(summary["worst_ppd"])
            rows.append(
                '<td style="background:{0}"><strong>{1:.1f}%</strong><span>Worst {2:.1f}% · PMV {3:.2f}</span></td>'.format(
                    color_for_thermal_state(mean_pmv, mean_ppd),
                    mean_ppd,
                    worst_ppd,
                    mean_pmv,
                )
            )
        rows.append("</tr>")
    rows.append("</tbody></table></div>")
    return "".join(rows)


def main() -> None:
    manifest = read_json(MANIFEST_PATH)
    summary_rows = read_csv(SUMMARY_PATH)
    timeseries_rows = read_csv(TIMESERIES_PATH)

    plan_path = ROOT_DIR / str(manifest["plan"]["path"])
    plan_payload = read_json(plan_path)
    rooms = list(plan_payload["rooms"])
    grid_cols = int(plan_payload["meta"]["grid_cols"])
    grid_rows = int(plan_payload["meta"]["grid_rows"])

    occupied_rooms = set(manifest["plan"]["occupied_rooms"])
    abbreviations = dict(manifest["plan"]["room_abbreviations"])

    exposure_by_room = {
        room_name: "south" if room_name in {
            "Service Vestibule",
            "Guest Room",
            "Family Lounge East",
            "Bath 2",
            "Bedroom 3",
            "Laundry",
            "Bedroom 1",
            "Bedroom 2",
        } else "east_west" if room_name in {"Primary Bedroom"} else "north_interior"
        for room_name in occupied_rooms
    }

    scenario_meta = {item["id"]: item for item in manifest["scenarios"]}
    summary_by_scenario = {row["scenario_id"]: row for row in summary_rows}

    records_by_scenario_room: Dict[str, Dict[str, List[Dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    ppd_values_by_scenario: Dict[str, List[float]] = defaultdict(list)
    for row in timeseries_rows:
        records_by_scenario_room[row["scenario_id"]][row["room_name"]].append(row)
        if row["ppd"] not in ("", None):
            ppd_values_by_scenario[row["scenario_id"]].append(as_float(row["ppd"]))

    mean_ppd_by_scenario = {
        scenario_id: mean(values) if values else 0.0
        for scenario_id, values in ppd_values_by_scenario.items()
    }
    mean_pmv_by_scenario = {
        row["scenario_id"]: as_float(row["mean_pmv"])
        for row in summary_rows
    }

    room_metrics_by_scenario: Dict[str, Dict[str, Dict[str, float]]] = {}
    for scenario_id, by_room in records_by_scenario_room.items():
        room_metrics_by_scenario[scenario_id] = {}
        for room_name, rows_for_room in by_room.items():
            room_metrics_by_scenario[scenario_id][room_name] = {
                "ppd": mean(as_float(item["ppd"]) for item in rows_for_room if item["ppd"] != ""),
                "pmv": mean(as_float(item["pmv"]) for item in rows_for_room if item["pmv"] != ""),
                "air_temp": mean(as_float(item["air_temp"]) for item in rows_for_room if item["air_temp"] != ""),
                "mean_radiant_temp": mean(
                    as_float(item["mean_radiant_temp"]) for item in rows_for_room if item["mean_radiant_temp"] != ""
                ),
                "humidity": mean(as_float(item["humidity"]) for item in rows_for_room if item["humidity"] != ""),
            }

    room_sensitivity: Dict[str, List[Tuple[str, float]]] = {}
    for cohort_key in manifest["cohorts"]:
        values: Dict[str, List[float]] = defaultdict(list)
        for scenario_id, meta in scenario_meta.items():
            if meta["cohort"] != cohort_key:
                continue
            for room_name, metric in room_metrics_by_scenario.get(scenario_id, {}).items():
                if room_name in occupied_rooms:
                    values[room_name].append(metric["ppd"])
        room_sensitivity[cohort_key] = sorted(
            ((room_name, mean(ppds)) for room_name, ppds in values.items()),
            key=lambda item: (-item[1], item[0]),
        )

    cohort_agents = {
        cohort_key: data["agents"]
        for cohort_key, data in manifest["cohorts"].items()
    }

    temperatures = sorted({float(item["base_ta_c"]) for item in manifest["scenarios"]})
    humidities = sorted({float(item["base_rh_pct"]) for item in manifest["scenarios"]})

    ordered_scenarios = sorted(manifest["scenarios"], key=scenario_sort_key)

    findings = []
    for cohort_key, meta in manifest["cohorts"].items():
        cohort_summaries = [
            summary_by_scenario[item["id"]]
            for item in ordered_scenarios
            if item["cohort"] == cohort_key
        ]
        hottest = max(cohort_summaries, key=lambda row: as_float(row["mean_pmv"]))
        riskiest = max(cohort_summaries, key=lambda row: as_float(row["worst_ppd"]))
        top_room = room_sensitivity[cohort_key][0]
        findings.append(
            "<li><strong>{0}</strong>: highest mean PMV at <code>{1}</code>, highest single-room PPD at "
            "<code>{2}</code>, most consistently stressed room <strong>{3}</strong> ({4:.1f}% mean PPD).</li>".format(
                html.escape(meta["label"]),
                html.escape(hottest["scenario_id"]),
                html.escape(riskiest["scenario_id"]),
                html.escape(top_room[0]),
                top_room[1],
            )
        )

    reference_plan = render_plan_svg(
        rooms=rooms,
        cols=grid_cols,
        rows=grid_rows,
        metric_by_room={},
        occupied_rooms=occupied_rooms,
        abbreviations=abbreviations,
        exposure_by_room=exposure_by_room,
        label_rooms=True,
        exposure_mode=True,
    )

    cohort_sections: List[str] = []
    for cohort_key, cohort_meta in manifest["cohorts"].items():
        table_html = render_metric_table(
            cohort_label=str(cohort_meta["label"]),
            cohort_key=cohort_key,
            temperatures=temperatures,
            humidities=humidities,
            summary_by_scenario=summary_by_scenario,
            mean_pmv_by_scenario=mean_pmv_by_scenario,
            mean_ppd_by_scenario=mean_ppd_by_scenario,
            scenario_meta=scenario_meta,
        )

        cards = []
        for scenario in ordered_scenarios:
            if scenario["cohort"] != cohort_key:
                continue
            scenario_id = scenario["id"]
            summary = summary_by_scenario[scenario_id]
            mean_ppd = mean_ppd_by_scenario.get(scenario_id, 0.0)
            svg = render_plan_svg(
                rooms=rooms,
                cols=grid_cols,
                rows=grid_rows,
                metric_by_room=room_metrics_by_scenario.get(scenario_id, {}),
                occupied_rooms=occupied_rooms,
                abbreviations=abbreviations,
                exposure_by_room=exposure_by_room,
            )
            cards.append(
                '<article class="scenario-card"><h5>Ta {0:.1f} C · RH {1:.0f}%</h5>'
                '<div class="scenario-meta">Mean PMV {2:.2f} · Mean PPD {3:.1f}% · Worst room PPD {4:.1f}%</div>'
                '{5}</article>'.format(
                    float(scenario["base_ta_c"]),
                    float(scenario["base_rh_pct"]),
                    as_float(summary["mean_pmv"]),
                    mean_ppd,
                    as_float(summary["worst_ppd"]),
                    svg,
                )
            )

        roster_rows = []
        for agent in cohort_agents[cohort_key]:
            roster_rows.append(
                "<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td><td>{4} kg / {5} cm</td><td>{6}</td></tr>".format(
                    html.escape(agent["room_abbrev"]),
                    html.escape(agent["room_name"]),
                    int(agent["demographics"]["age"]),
                    html.escape(agent["demographics"]["gender"]),
                    agent["demographics"]["weight_kg"],
                    agent["demographics"]["height_cm"],
                    html.escape(agent["activity"]),
                )
            )

        stress_rows = []
        for room_name, avg_ppd in room_sensitivity[cohort_key][:6]:
            stress_rows.append(
                "<tr><td>{0}</td><td>{1:.1f}%</td></tr>".format(
                    html.escape(room_name),
                    avg_ppd,
                )
            )

        cohort_sections.append(
            """
            <section class="cohort-section">
              <div class="section-head">
                <h3>{label}</h3>
                <p>{description}</p>
              </div>
              <div class="table-grid">
                {table_html}
                <div class="table-card">
                  <h4>Occupant Roster</h4>
                  <table class="roster">
                    <thead><tr><th>Tag</th><th>Room</th><th>Age</th><th>Gender</th><th>Body</th><th>Activity</th></tr></thead>
                    <tbody>{roster_rows}</tbody>
                  </table>
                </div>
                <div class="table-card">
                  <h4>Most Sensitive Rooms</h4>
                  <table class="roster">
                    <thead><tr><th>Room</th><th>Mean PPD</th></tr></thead>
                    <tbody>{stress_rows}</tbody>
                  </table>
                </div>
              </div>
              <div class="scenario-grid">{cards}</div>
            </section>
            """.format(
                label=html.escape(str(cohort_meta["label"])),
                description=html.escape(str(cohort_meta["description"])),
                table_html=table_html,
                roster_rows="".join(roster_rows),
                stress_rows="".join(stress_rows),
                cards="".join(cards),
            )
        )

    legend_rows = []
    for room_name in manifest["plan"]["occupied_rooms"]:
        legend_rows.append(
            "<tr><td>{0}</td><td>{1}</td><td>{2}</td></tr>".format(
                html.escape(abbreviations[room_name]),
                html.escape(room_name),
                html.escape(exposure_by_room.get(room_name, "north_interior")),
            )
        )

    html_out = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Apartment Daytime Cooled Study</title>
  <style>
    :root {{
      --bg: #f7f2ea;
      --panel: #fffdfa;
      --ink: #211d1a;
      --muted: #665f58;
      --border: #ddd2c2;
      --shadow: 0 10px 24px rgba(33, 29, 26, 0.08);
      --accent: #8b7355;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "SF Mono", "IBM Plex Mono", "Consolas", monospace;
      color: var(--ink);
      background:
        radial-gradient(1100px 700px at 0% 0%, #fffaf0 0%, transparent 55%),
        linear-gradient(180deg, var(--bg), #efe6d8);
    }}
    .page {{
      max-width: 1480px;
      margin: 0 auto;
      padding: 28px 20px 40px;
    }}
    .hero, .panel, .cohort-section {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 18px;
      box-shadow: var(--shadow);
    }}
    .hero {{
      padding: 22px 24px;
      margin-bottom: 20px;
    }}
    h1, h2, h3, h4, h5 {{
      margin: 0;
      font-weight: 600;
    }}
    .hero p, .section-head p, .muted {{
      color: var(--muted);
      line-height: 1.5;
    }}
    .hero-grid {{
      margin-top: 16px;
      display: grid;
      grid-template-columns: 380px 1fr;
      gap: 18px;
      align-items: start;
    }}
    .legend-table, .roster, .matrix {{
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
    }}
    .legend-table th, .legend-table td,
    .roster th, .roster td,
    .matrix th, .matrix td {{
      border-bottom: 1px solid #eee4d5;
      padding: 8px 7px;
      text-align: left;
      vertical-align: top;
    }}
    .matrix td {{
      min-width: 94px;
    }}
    .matrix td span {{
      display: block;
      margin-top: 4px;
      color: var(--muted);
      font-size: 11px;
    }}
    .thermal-legend {{
      display: grid;
      gap: 6px;
      margin: 10px 0 14px;
      padding: 10px 12px;
      border: 1px solid var(--border);
      border-radius: 12px;
      background: #fffcf8;
      font-size: 12px;
    }}
    .swatch-row {{
      display: flex;
      align-items: center;
      gap: 8px;
    }}
    .swatch {{
      width: 14px;
      height: 14px;
      border-radius: 999px;
      border: 1px solid rgba(33, 29, 26, 0.16);
      flex: 0 0 auto;
    }}
    .legend-note {{
      color: var(--muted);
      line-height: 1.45;
    }}
    .panel {{
      padding: 18px 20px;
      margin-bottom: 20px;
    }}
    .findings ul {{
      margin: 10px 0 0 18px;
      padding: 0;
    }}
    .findings li {{
      margin-bottom: 8px;
      line-height: 1.5;
    }}
    .cohort-section {{
      padding: 18px 18px 22px;
      margin-bottom: 20px;
    }}
    .section-head {{
      display: flex;
      justify-content: space-between;
      gap: 20px;
      align-items: baseline;
      margin-bottom: 16px;
    }}
    .table-grid {{
      display: grid;
      grid-template-columns: 1.2fr 1.1fr 0.9fr;
      gap: 14px;
      margin-bottom: 18px;
    }}
    .table-card {{
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 14px;
      background: #fffcf8;
    }}
    .table-card h4 {{
      margin-bottom: 10px;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      color: var(--accent);
    }}
    .scenario-grid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 14px;
    }}
    .scenario-card {{
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 12px;
      background: #fff;
    }}
    .scenario-card h5 {{
      font-size: 13px;
      margin-bottom: 6px;
    }}
    .scenario-meta {{
      font-size: 11px;
      color: var(--muted);
      margin-bottom: 10px;
    }}
    svg {{
      width: 100%;
      height: auto;
      display: block;
    }}
    code {{
      background: #f3ece2;
      padding: 1px 4px;
      border-radius: 4px;
    }}
    @media (max-width: 1200px) {{
      .hero-grid, .table-grid, .scenario-grid {{
        grid-template-columns: 1fr;
      }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <section class="hero">
      <h1>{title}</h1>
      <p>Deterministic apartment population sweep using the traced apartment template, multiple occupant cohorts, and a daytime cooled room-environment matrix with orientation-based MRT offsets.</p>
      <div class="hero-grid">
        <div>
          {reference_plan}
        </div>
        <div>
          <div class="panel" style="margin:0; box-shadow:none;">
            <h2 style="font-size:14px; text-transform:uppercase; letter-spacing:0.08em; color:var(--accent);">Case Setup</h2>
            <p class="muted">{orientation} {operation}</p>
            <p class="muted">{thermal_model}</p>
            <div class="thermal-legend">
              <div class="swatch-row"><span class="swatch" style="background:#6f9dcb"></span><span>Cold-side stress: PMV &lt; 0</span></div>
              <div class="swatch-row"><span class="swatch" style="background:#e8e1d4"></span><span>Near-neutral: PMV around 0</span></div>
              <div class="swatch-row"><span class="swatch" style="background:#e38c69"></span><span>Hot-side stress: PMV &gt; 0</span></div>
              <div class="legend-note">Plan fills use PMV sign for hue, with shade deepening as discomfort severity increases.</div>
            </div>
            <table class="legend-table">
              <thead><tr><th>Tag</th><th>Occupied Room</th><th>Exposure Rule</th></tr></thead>
              <tbody>{legend_rows}</tbody>
            </table>
          </div>
        </div>
      </div>
    </section>

    <section class="panel findings">
      <h2 style="font-size:14px; text-transform:uppercase; letter-spacing:0.08em; color:var(--accent);">Auto Findings</h2>
      <ul>{findings}</ul>
    </section>

    {cohort_sections}
  </div>
</body>
</html>
""".format(
        title=html.escape(str(manifest["title"])),
        reference_plan=reference_plan,
        orientation=html.escape(str(manifest["assumptions"]["orientation"])),
        operation=html.escape(str(manifest["assumptions"]["operation"])),
        thermal_model=html.escape(str(manifest["assumptions"]["thermal_model"])),
        legend_rows="".join(legend_rows),
        findings="".join(findings),
        cohort_sections="".join(cohort_sections),
    )

    REPORT_PATH.write_text(html_out, encoding="utf-8")
    print("Wrote report: {0}".format(REPORT_PATH))


if __name__ == "__main__":
    main()
