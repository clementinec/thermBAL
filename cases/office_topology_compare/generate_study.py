#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import random
import statistics
import sys
from collections import deque
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

CASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = CASE_DIR.parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from comfort_engine import ComfortState, Occupant, compute_pmv, estimate_metabolic_profile
from floor_plan.io import load_template_json

from build_templates import (
    CELLULAR_ROOM_RENAMES,
    CELLULAR_TEMPLATE,
    OPEN_TEMPLATE,
    derive_templates,
)

Coord = Tuple[int, int]

STUDY_PATH = CASE_DIR / "study_manifest.json"
OUT_DIR = CASE_DIR / "out"
SUMMARY_PATH = OUT_DIR / "summary.csv"
COMPARISON_PATH = OUT_DIR / "comparison_summary.csv"
SAMPLES_PATH = OUT_DIR / "sample_results.csv"
FIELD_MAPS_PATH = OUT_DIR / "field_maps.json"

SEED = 20260408
BASE_TEMPERATURES = [20.0, 22.0, 23.5, 25.0, 28.0]
BASE_HUMIDITIES = [35.0, 55.0, 75.0, 80.0]

ACTIVITY_MULTIPLIER = {
    "seated": 1.0,
    "standing": 1.2,
    "walking": 1.6,
}

ACTIVITY_POSTURE = {
    "seated": "seated",
    "standing": "standing",
    "walking": "standing",
}

COHORTS = {
    "realistic_mixed": {
        "label": "Realistic Mixed Office",
        "description": "Mixed office cohort with male/female dress differences, ages concentrated in the 20s to 50s.",
        "sex_mode": "mixed",
        "male_share": 0.45,
        "age_mean": 40.0,
        "age_sd": 6.0,
        "age_min": 28.0,
        "age_max": 58.0,
        "male_clo": 1.00,
        "female_clo": 0.72,
        "bmi_shift": 0.0,
    },
    "default_male_35": {
        "label": "Default Male Baseline",
        "description": "Male-only baseline cohort fixed around age 35 with typical long-sleeve office clothing.",
        "sex_mode": "male",
        "male_share": 1.0,
        "age_mean": 35.0,
        "age_sd": 0.0,
        "age_min": 35.0,
        "age_max": 35.0,
        "male_clo": 1.00,
        "female_clo": 1.00,
        "bmi_shift": 0.0,
    },
    "realistic_male": {
        "label": "Realistic Male Office",
        "description": "Male-only office cohort with realistic age and body variation but the same typical long-sleeve office clothing.",
        "sex_mode": "male",
        "male_share": 1.0,
        "age_mean": 40.0,
        "age_sd": 6.0,
        "age_min": 28.0,
        "age_max": 58.0,
        "male_clo": 1.00,
        "female_clo": 1.00,
        "bmi_shift": 0.1,
    },
    "female_light": {
        "label": "Female Light Office",
        "description": "Female office cohort with lighter summer office attire and lower clothing insulation than the male formal group.",
        "sex_mode": "female",
        "male_share": 0.0,
        "age_mean": 34.0,
        "age_sd": 6.0,
        "age_min": 24.0,
        "age_max": 52.0,
        "male_clo": 0.72,
        "female_clo": 0.68,
        "bmi_shift": -0.1,
    },
    "older_mixed": {
        "label": "Older Mixed Office",
        "description": "Mixed office cohort shifted into late-career ages while preserving the same sex-weighted clothing logic.",
        "sex_mode": "mixed",
        "male_share": 0.55,
        "age_mean": 61.0,
        "age_sd": 4.5,
        "age_min": 54.0,
        "age_max": 72.0,
        "male_clo": 1.04,
        "female_clo": 0.78,
        "bmi_shift": 0.35,
    },
    "older_male_formal": {
        "label": "Older Male Formal",
        "description": "Older male office cohort in higher-clo formal dress, representing the cold-sensitive end of office populations.",
        "sex_mode": "male",
        "male_share": 1.0,
        "age_mean": 64.0,
        "age_sd": 4.0,
        "age_min": 58.0,
        "age_max": 74.0,
        "male_clo": 1.06,
        "female_clo": 1.06,
        "bmi_shift": 0.45,
    },
    "older_female_light": {
        "label": "Older Female Light",
        "description": "Older female office cohort in lighter business-casual dress, contrasting against the older male formal profile.",
        "sex_mode": "female",
        "male_share": 0.0,
        "age_mean": 63.0,
        "age_sd": 4.0,
        "age_min": 56.0,
        "age_max": 74.0,
        "male_clo": 0.78,
        "female_clo": 0.74,
        "bmi_shift": 0.20,
    },
}

CELLULAR_WORK_ROOMS = {
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

OPEN_WORK_ROOMS = {
    "North Open Office",
    "South Open Office",
    "West Corner Office",
}

WORK_ROOM_ORDER = [
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
]

WORK_ROOM_ABBREVIATIONS = {
    "Reception": "RE",
    "Admin Bay": "AB",
    "Team Office West": "TW",
    "Team Office East": "TE",
    "Executive Office": "EO",
    "Focus Room": "FR",
    "Print Hub": "PH",
    "Collaboration West": "CW",
    "Collaboration East": "CE",
    "Meeting Room West": "MW",
    "Meeting Room East": "ME",
    "South Studio West": "SW",
    "South Studio East": "SE",
    "Project Bay West": "PW",
    "Project Bay East": "PE",
    "North Open Office": "NO",
    "South Open Office": "SO",
    "Core Lobby": "CL",
    "West Corner Office": "WC",
    "North WC": "NW",
    "East WC": "EW",
    "Wellness Room": "WR",
    "Core WC": "CW",
}

OFFICE_ROOM_ACTIVITY = {
    "Reception": "seated",
    "Admin Bay": "seated",
    "Team Office West": "seated",
    "Team Office East": "seated",
    "Executive Office": "seated",
    "Focus Room": "seated",
    "Print Hub": "seated",
    "Collaboration West": "seated",
    "Collaboration East": "seated",
    "Meeting Room West": "seated",
    "Meeting Room East": "seated",
    "South Studio West": "seated",
    "South Studio East": "seated",
    "Project Bay West": "seated",
    "Project Bay East": "seated",
}

REFERENCE_CASES = [
    {"cohort": "realistic_mixed", "base_ta_c": 22.0, "base_rh_pct": 80.0, "label": "Realistic Mixed / 22 C / 80% RH"},
    {"cohort": "default_male_35", "base_ta_c": 22.0, "base_rh_pct": 80.0, "label": "Default Male 35 / 22 C / 80% RH"},
    {"cohort": "realistic_male", "base_ta_c": 22.0, "base_rh_pct": 80.0, "label": "Realistic Male / 22 C / 80% RH"},
    {"cohort": "female_light", "base_ta_c": 22.0, "base_rh_pct": 80.0, "label": "Female Light / 22 C / 80% RH"},
    {"cohort": "realistic_mixed", "base_ta_c": 25.0, "base_rh_pct": 55.0, "label": "Realistic Mixed / 25 C / 55% RH"},
    {"cohort": "realistic_mixed", "base_ta_c": 28.0, "base_rh_pct": 75.0, "label": "Realistic Mixed / 28 C / 75% RH"},
]

ORIGINAL_TO_CELLULAR = {source: target for source, target in CELLULAR_ROOM_RENAMES.items()}


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def temp_label(temp_c: float) -> str:
    return str(temp_c).replace(".", "p")


def neighbor_cells(cell: Coord) -> List[Coord]:
    col, row = cell
    return [(col + 1, row), (col - 1, row), (col, row + 1), (col, row - 1)]


def bfs_depth(active_cells: Set[Coord], seed_cells: Iterable[Coord]) -> Dict[Coord, int]:
    depths = {cell: 999 for cell in active_cells}
    queue: deque[Coord] = deque()
    for cell in seed_cells:
        if cell in active_cells and depths[cell] > 0:
            depths[cell] = 0
            queue.append(cell)
    while queue:
        current = queue.popleft()
        current_depth = depths[current]
        for nxt in neighbor_cells(current):
            if nxt in active_cells and depths[nxt] > current_depth + 1:
                depths[nxt] = current_depth + 1
                queue.append(nxt)
    return depths


def centroid_cell(cells: List[Coord]) -> Coord:
    cx = sum(col + 0.5 for col, _ in cells) / len(cells)
    cy = sum(row + 0.5 for _, row in cells) / len(cells)
    return min(cells, key=lambda cell: ((cell[0] + 0.5 - cx) ** 2 + (cell[1] + 0.5 - cy) ** 2, cell[1], cell[0]))


def farthest_point_sample(cells: List[Coord], count: int) -> List[Coord]:
    cells = sorted(cells, key=lambda item: (item[1], item[0]))
    if count >= len(cells):
        return cells
    selected = [centroid_cell(cells)]
    while len(selected) < count:
        next_cell = max(
            cells,
            key=lambda cell: (
                min((cell[0] - other[0]) ** 2 + (cell[1] - other[1]) ** 2 for other in selected),
                -cell[1],
                -cell[0],
            ),
        )
        if next_cell in selected:
            break
        selected.append(next_cell)
    return sorted(selected, key=lambda item: (item[1], item[0]))


def sample_latents(sample_points: List[Dict[str, object]]) -> Dict[str, Dict[str, float]]:
    rng = random.Random(SEED)
    latents: Dict[str, Dict[str, float]] = {}
    for sample in sample_points:
        latents[str(sample["sample_id"])] = {
            "sex_roll": rng.random(),
            "age_z": rng.gauss(0.0, 1.0),
            "height_z": rng.gauss(0.0, 1.0),
            "bmi_z": rng.gauss(0.0, 1.0),
        }
    return latents


def sample_count_for_room(cell_count: int) -> int:
    if cell_count >= 24:
        return 4
    if cell_count >= 14:
        return 3
    if cell_count >= 8:
        return 2
    return 1


def build_agents(sample_points: List[Dict[str, object]]) -> Dict[str, List[Dict[str, object]]]:
    latents = sample_latents(sample_points)
    out: Dict[str, List[Dict[str, object]]] = {}
    for cohort_key, cohort in COHORTS.items():
        mixed_sex_by_sample: Dict[str, str] = {}
        if cohort["sex_mode"] == "mixed":
            target_males = int(round(len(sample_points) * float(cohort["male_share"])))
            ordered_ids = [
                str(sample["sample_id"])
                for sample in sorted(sample_points, key=lambda sample: latents[str(sample["sample_id"])]["sex_roll"])
            ]
            for index, sample_id in enumerate(ordered_ids):
                mixed_sex_by_sample[sample_id] = "male" if index < target_males else "female"

        agents: List[Dict[str, object]] = []
        for index, sample in enumerate(sample_points, start=1):
            latent = latents[str(sample["sample_id"])]
            if cohort["sex_mode"] == "male":
                sex = "male"
            elif cohort["sex_mode"] == "female":
                sex = "female"
            else:
                sex = mixed_sex_by_sample[str(sample["sample_id"])]

            age = clamp(
                float(cohort["age_mean"]) + float(cohort["age_sd"]) * latent["age_z"],
                float(cohort["age_min"]),
                float(cohort["age_max"]),
            )

            if sex == "male":
                height_cm = clamp(173.0 + 6.0 * latent["height_z"], 160.0, 190.0)
                bmi = clamp(23.6 + 2.4 * latent["bmi_z"] + float(cohort["bmi_shift"]), 19.0, 32.0)
                clo = float(cohort["male_clo"])
            else:
                height_cm = clamp(161.0 + 6.0 * latent["height_z"], 150.0, 178.0)
                bmi = clamp(22.1 + 2.4 * latent["bmi_z"] + float(cohort["bmi_shift"]), 18.0, 31.0)
                clo = float(cohort["female_clo"])

            activity = OFFICE_ROOM_ACTIVITY.get(str(sample["source_room"]), "seated")

            weight_kg = bmi * (height_cm / 100.0) ** 2
            clo = clamp(clo + 0.03 * latent["bmi_z"], 0.55, 1.10)

            agents.append(
                {
                    "id": f"worker_{index:02d}",
                    "sample_id": sample["sample_id"],
                    "position": list(sample["cell"]),
                    "source_room": sample["source_room"],
                    "activity": activity,
                    "clo": round(clo, 2),
                    "demographics": {
                        "age": int(round(age)),
                        "gender": sex,
                        "weight_kg": round(weight_kg, 1),
                        "height_cm": int(round(height_cm)),
                    },
                }
            )
        out[cohort_key] = agents
    return out


def classify_window(feature: Dict[str, object], active_cells: Set[Coord]) -> Tuple[Coord | None, str | None]:
    orient = str(feature["orient"])
    row = int(feature["r"])
    col = int(feature["c"])
    if orient == "h":
        above = (col, row - 1)
        below = (col, row)
        if above in active_cells and below not in active_cells:
            return above, "south"
        if below in active_cells and above not in active_cells:
            return below, "north"
        return None, None
    left = (col - 1, row)
    right = (col, row)
    if left in active_cells and right not in active_cells:
        return left, "east"
    if right in active_cells and left not in active_cells:
        return right, "west"
    return None, None


def make_field(active_cells: Set[Coord], features: List[Dict[str, object]], base_ta: float, base_rh: float) -> Dict[Coord, Dict[str, float]]:
    perimeter_seed = []
    for cell in active_cells:
        if any(neighbor not in active_cells for neighbor in neighbor_cells(cell)):
            perimeter_seed.append(cell)

    window_seed: Set[Coord] = set()
    south_window_seed: Set[Coord] = set()
    for feature in features:
        if feature.get("type") != "window":
            continue
        cell, orientation = classify_window(feature, active_cells)
        if cell is None:
            continue
        window_seed.add(cell)
        if orientation == "south":
            south_window_seed.add(cell)

    perimeter_depth = bfs_depth(active_cells, perimeter_seed)
    window_depth = bfs_depth(active_cells, window_seed)
    south_depth = bfs_depth(active_cells, south_window_seed)

    field: Dict[Coord, Dict[str, float]] = {}
    for cell in active_cells:
        p_depth = float(perimeter_depth[cell])
        w_depth = float(window_depth[cell])
        s_depth = float(south_depth[cell])

        perimeter_factor = math.exp(-p_depth / 1.45)
        window_factor = 0.0 if w_depth >= 999 else math.exp(-w_depth / 0.85)
        south_factor = 0.0 if s_depth >= 999 else math.exp(-s_depth / 0.95)

        air_temp = base_ta + 0.18 * perimeter_factor + 0.12 * window_factor + 0.08 * south_factor
        mean_radiant_temp = base_ta + 0.55 * perimeter_factor + 1.65 * window_factor + 0.95 * south_factor
        air_velocity = 0.07 + 0.02 * perimeter_factor + 0.012 * window_factor

        field[cell] = {
            "air_temp": round(air_temp, 3),
            "mean_radiant_temp": round(mean_radiant_temp, 3),
            "humidity": round(base_rh, 3),
            "air_velocity": round(air_velocity, 3),
            "perimeter_depth": int(perimeter_depth[cell]),
            "window_depth": int(window_depth[cell]) if window_depth[cell] < 999 else 999,
            "south_window_depth": int(south_depth[cell]) if south_depth[cell] < 999 else 999,
        }
    return field


def room_average_envs(rooms: Dict[str, Set[Coord]], field: Dict[Coord, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    for room_name, cells in rooms.items():
        envs = [field[cell] for cell in cells if cell in field]
        out[room_name] = {
            "air_temp": round(sum(env["air_temp"] for env in envs) / len(envs), 3),
            "mean_radiant_temp": round(sum(env["mean_radiant_temp"] for env in envs) / len(envs), 3),
            "humidity": round(sum(env["humidity"] for env in envs) / len(envs), 3),
            "air_velocity": round(sum(env["air_velocity"] for env in envs) / len(envs), 3),
        }
    return out


def compute_comfort(agent: Dict[str, object], env: Dict[str, float]) -> Dict[str, float]:
    demo = agent["demographics"]
    activity = str(agent["activity"])
    profile = estimate_metabolic_profile(
        sex=str(demo["gender"]),
        weight_kg=float(demo["weight_kg"]),
        height_cm=float(demo["height_cm"]),
        age=int(demo["age"]),
        activity_multiplier=ACTIVITY_MULTIPLIER[activity],
    )
    state = ComfortState(
        Ta=float(env["air_temp"]),
        Tr=float(env["mean_radiant_temp"]),
        RH=float(env["humidity"]),
        v=float(env["air_velocity"]),
    )
    occupant = Occupant(
        m_wm2=profile["M_active_Wm2"],
        clo=float(agent["clo"]),
        posture=ACTIVITY_POSTURE[activity],
        body_area_m2=profile["BSA_m2"],
    )
    result = compute_pmv(state, occupant)
    return {
        "pmv": round(result.PMV, 3),
        "ppd": round(result.PPD, 3),
    }


def summarize_records(records: List[Dict[str, object]]) -> Dict[str, float]:
    pmvs = [float(record["pmv"]) for record in records]
    ppds = [float(record["ppd"]) for record in records]
    perimeter = [float(record["pmv"]) for record in records if int(record["perimeter_depth"]) <= 1]
    core = [float(record["pmv"]) for record in records if int(record["perimeter_depth"]) >= 3]
    cold_risk = [record for record in records if float(record["ppd"]) >= 20.0 and float(record["pmv"]) < 0.0]
    hot_risk = [record for record in records if float(record["ppd"]) >= 20.0 and float(record["pmv"]) > 0.0]

    return {
        "point_count": len(records),
        "mean_pmv": round(statistics.mean(pmvs), 3),
        "min_pmv": round(min(pmvs), 3),
        "max_pmv": round(max(pmvs), 3),
        "pmv_range": round(max(pmvs) - min(pmvs), 3),
        "spatial_std_pmv": round(statistics.pstdev(pmvs), 3),
        "mean_ppd": round(statistics.mean(ppds), 3),
        "worst_ppd": round(max(ppds), 3),
        "cold_risk_share_pct": round(100.0 * len(cold_risk) / len(records), 2),
        "hot_risk_share_pct": round(100.0 * len(hot_risk) / len(records), 2),
        "perimeter_core_gap": round((statistics.mean(perimeter) - statistics.mean(core)) if perimeter and core else 0.0, 3),
    }


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def field_map_for_case(
    topology: str,
    rooms_by_name: Dict[str, Set[Coord]],
    active_cells: Set[Coord],
    cohort_agents: List[Dict[str, object]],
    field: Dict[Coord, Dict[str, float]],
) -> List[Dict[str, object]]:
    room_avg = room_average_envs(rooms_by_name, field)
    cell_to_room = {cell: room_name for room_name, cells in rooms_by_name.items() for cell in cells}
    work_cells = set().union(*(cells for name, cells in rooms_by_name.items() if (name in OPEN_WORK_ROOMS or name in CELLULAR_WORK_ROOMS)))
    cells: List[Dict[str, object]] = []
    for cell in sorted(active_cells, key=lambda item: (item[1], item[0])):
        room_name = cell_to_room[cell]
        env = field[cell] if topology == "open_office" else room_avg[room_name]
        pmvs = [compute_comfort(agent, env)["pmv"] for agent in cohort_agents]
        cells.append(
            {
                "x": cell[0],
                "y": cell[1],
                "room_name": room_name,
                "is_work_cell": cell in work_cells,
                "pmv": round(statistics.mean(pmvs), 3),
                "air_temp": env["air_temp"],
                "mean_radiant_temp": env["mean_radiant_temp"],
            }
        )
    return cells


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    derive_templates()

    cellular_state = load_template_json(CELLULAR_TEMPLATE)
    open_state = load_template_json(OPEN_TEMPLATE)

    room_cells = {room.name: set(room.cells) for room in cellular_state.plan.rooms}
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

    agents_by_cohort = build_agents(work_points)
    active_cells = set(cellular_state.active_cells)
    cellular_rooms = {room.name: set(room.cells) for room in cellular_state.plan.rooms}
    open_rooms = {room.name: set(room.cells) for room in open_state.plan.rooms}

    sample_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    comparison_rows: List[Dict[str, object]] = []
    summary_index: Dict[Tuple[str, str, float, float], Dict[str, object]] = {}
    field_cache: Dict[Tuple[float, float], Dict[Coord, Dict[str, float]]] = {}

    for base_ta in BASE_TEMPERATURES:
        for base_rh in BASE_HUMIDITIES:
            field = make_field(active_cells, cellular_state.features, base_ta, base_rh)
            field_cache[(base_ta, base_rh)] = field
            cellular_room_env = room_average_envs(cellular_rooms, field)

            for cohort_key, cohort_meta in COHORTS.items():
                for topology in ("cellular_office", "open_office"):
                    records: List[Dict[str, object]] = []
                    for agent in agents_by_cohort[cohort_key]:
                        cell = tuple(agent["position"])
                        cell_field = field[cell]
                        if topology == "cellular_office":
                            topology_room = cellular_state.plan.room_at(cell)
                            env = cellular_room_env[topology_room.name]
                        else:
                            topology_room = open_state.plan.room_at(cell)
                            env = {
                                "air_temp": cell_field["air_temp"],
                                "mean_radiant_temp": cell_field["mean_radiant_temp"],
                                "humidity": cell_field["humidity"],
                                "air_velocity": cell_field["air_velocity"],
                            }
                        comfort = compute_comfort(agent, env)
                        record = {
                            "scenario_id": f"{topology}__{cohort_key}__t{temp_label(base_ta)}_rh{int(base_rh)}",
                            "topology": topology,
                            "topology_label": "Cellular Office" if topology == "cellular_office" else "Open Office",
                            "cohort": cohort_key,
                            "cohort_label": cohort_meta["label"],
                            "sample_id": agent["sample_id"],
                            "source_room": agent["source_room"],
                            "topology_room": topology_room.name if topology_room else "Outside",
                            "x": cell[0],
                            "y": cell[1],
                            "base_ta_c": base_ta,
                            "base_rh_pct": base_rh,
                            "age": agent["demographics"]["age"],
                            "gender": agent["demographics"]["gender"],
                            "clo": agent["clo"],
                            "activity": agent["activity"],
                            "air_temp": env["air_temp"],
                            "mean_radiant_temp": env["mean_radiant_temp"],
                            "humidity": env["humidity"],
                            "air_velocity": env["air_velocity"],
                            "perimeter_depth": cell_field["perimeter_depth"],
                            "window_depth": cell_field["window_depth"],
                            "south_window_depth": cell_field["south_window_depth"],
                            "pmv": comfort["pmv"],
                            "ppd": comfort["ppd"],
                        }
                        records.append(record)
                        sample_rows.append(record)

                    summary = {
                        "scenario_id": f"{topology}__{cohort_key}__t{temp_label(base_ta)}_rh{int(base_rh)}",
                        "topology": topology,
                        "topology_label": "Cellular Office" if topology == "cellular_office" else "Open Office",
                        "cohort": cohort_key,
                        "cohort_label": cohort_meta["label"],
                        "base_ta_c": base_ta,
                        "base_rh_pct": base_rh,
                        **summarize_records(records),
                    }
                    summary_rows.append(summary)
                    summary_index[(topology, cohort_key, base_ta, base_rh)] = summary

                cellular_summary = summary_index[("cellular_office", cohort_key, base_ta, base_rh)]
                open_summary = summary_index[("open_office", cohort_key, base_ta, base_rh)]
                comparison_rows.append(
                    {
                        "cohort": cohort_key,
                        "cohort_label": cohort_meta["label"],
                        "base_ta_c": base_ta,
                        "base_rh_pct": base_rh,
                        "mean_pmv_delta_open_minus_cellular": round(open_summary["mean_pmv"] - cellular_summary["mean_pmv"], 3),
                        "pmv_range_delta_open_minus_cellular": round(open_summary["pmv_range"] - cellular_summary["pmv_range"], 3),
                        "spatial_std_delta_open_minus_cellular": round(open_summary["spatial_std_pmv"] - cellular_summary["spatial_std_pmv"], 3),
                        "worst_ppd_delta_open_minus_cellular": round(open_summary["worst_ppd"] - cellular_summary["worst_ppd"], 3),
                        "perimeter_core_gap_delta_open_minus_cellular": round(open_summary["perimeter_core_gap"] - cellular_summary["perimeter_core_gap"], 3),
                    }
                )

    field_maps = {
        "case_id": "office_topology_compare",
        "reference_cases": REFERENCE_CASES,
        "maps": [],
    }
    for reference in REFERENCE_CASES:
        cohort_key = reference["cohort"]
        base_ta = reference["base_ta_c"]
        base_rh = reference["base_rh_pct"]
        field = field_cache[(base_ta, base_rh)]
        for topology, rooms in (("cellular_office", cellular_rooms), ("open_office", open_rooms)):
            field_maps["maps"].append(
                {
                    "topology": topology,
                    "cohort": cohort_key,
                    "base_ta_c": base_ta,
                    "base_rh_pct": base_rh,
                    "label": reference["label"],
                    "cells": field_map_for_case(topology, rooms, active_cells, agents_by_cohort[cohort_key], field),
                }
            )

    study_manifest = {
        "case_id": "office_topology_compare",
        "title": "Office Topology Comparison on Shared Apartment Footprint",
        "seed": SEED,
        "templates": {
            "cellular_office": str(CELLULAR_TEMPLATE.relative_to(ROOT_DIR)),
            "open_office": str(OPEN_TEMPLATE.relative_to(ROOT_DIR)),
        },
        "assumptions": {
            "orientation": "Bottom of plan is south.",
            "operation": "Daytime cooled office condition.",
            "comparison_logic": (
                "The same perimeter-to-window thermal field is applied to both topologies. "
                "Cellular office rooms receive room-mean field values; the open office keeps the field positionally resolved."
            ),
            "field_model": (
                "Perimeter, window, and south-window depth decay control Ta and MRT. "
                "At 25 C air temperature, near-window MRT can rise above 27 C at peak exposed edges, while the core stays closer to air temperature."
            ),
            "occupancy_rule": (
                "Deterministic desk/sample points are placed across the same shared work cells in both topologies so the comparison isolates topology and field resolution."
            ),
        },
        "cohorts": COHORTS,
        "sample_points": work_points,
        "scenario_grid": {
            "base_temperatures_c": BASE_TEMPERATURES,
            "base_relative_humidity_pct": BASE_HUMIDITIES,
            "scenario_count": len(BASE_TEMPERATURES) * len(BASE_HUMIDITIES) * len(COHORTS) * 2,
        },
    }

    STUDY_PATH.write_text(json.dumps(study_manifest, indent=2), encoding="utf-8")
    FIELD_MAPS_PATH.write_text(json.dumps(field_maps, indent=2), encoding="utf-8")

    write_csv(
        SUMMARY_PATH,
        summary_rows,
        [
            "scenario_id",
            "topology",
            "topology_label",
            "cohort",
            "cohort_label",
            "base_ta_c",
            "base_rh_pct",
            "point_count",
            "mean_pmv",
            "min_pmv",
            "max_pmv",
            "pmv_range",
            "spatial_std_pmv",
            "mean_ppd",
            "worst_ppd",
            "cold_risk_share_pct",
            "hot_risk_share_pct",
            "perimeter_core_gap",
        ],
    )
    write_csv(
        COMPARISON_PATH,
        comparison_rows,
        [
            "cohort",
            "cohort_label",
            "base_ta_c",
            "base_rh_pct",
            "mean_pmv_delta_open_minus_cellular",
            "pmv_range_delta_open_minus_cellular",
            "spatial_std_delta_open_minus_cellular",
            "worst_ppd_delta_open_minus_cellular",
            "perimeter_core_gap_delta_open_minus_cellular",
        ],
    )
    write_csv(
        SAMPLES_PATH,
        sample_rows,
        [
            "scenario_id",
            "topology",
            "topology_label",
            "cohort",
            "cohort_label",
            "sample_id",
            "source_room",
            "topology_room",
            "x",
            "y",
            "base_ta_c",
            "base_rh_pct",
            "age",
            "gender",
            "clo",
            "activity",
            "air_temp",
            "mean_radiant_temp",
            "humidity",
            "air_velocity",
            "perimeter_depth",
            "window_depth",
            "south_window_depth",
            "pmv",
            "ppd",
        ],
    )

    print(f"Wrote {STUDY_PATH}")
    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote {COMPARISON_PATH}")
    print(f"Wrote {SAMPLES_PATH}")
    print(f"Wrote {FIELD_MAPS_PATH}")
    print(f"Scenarios: {len(summary_rows)} summaries")
    print(f"Sample points: {len(work_points)}")


if __name__ == "__main__":
    main()
