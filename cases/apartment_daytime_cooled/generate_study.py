#!/usr/bin/env python3
from __future__ import annotations

import json
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple

CASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = CASE_DIR.parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from floor_plan.io import load_template_json

PLAN_PATH = ROOT_DIR / "floor_plan" / "apartment_template.json"
STUDY_PATH = CASE_DIR / "study.json"
MANIFEST_PATH = CASE_DIR / "study_manifest.json"
SEED = 20260407

BASE_TEMPERATURES = [20.0, 22.0, 23.5, 25.0, 28.0]
BASE_HUMIDITIES = [35.0, 55.0, 75.0, 80.0]

COHORTS = {
    "young_mixed": {
        "label": "Young Scholars Mixed",
        "description": "Mixed male/female cohort, ages 25-45, one occupant per occupied room.",
    },
    "young_male": {
        "label": "Young Scholars Male Only",
        "description": "Male-only control cohort using the same room assignment and latent demographic seeds.",
    },
    "older_shift": {
        "label": "Age-Shifted +50 Years",
        "description": "Same mixed cohort as young_mixed with every agent age shifted upward by 50 years.",
    },
    "higher_clo_sedentary": {
        "label": "Higher Clo, Sedentary",
        "description": "Mixed cohort with heavier clothing and reduced activity relative to the young mixed baseline.",
    },
    "lighter_clothing_mobile": {
        "label": "Light Clothing, Mobile",
        "description": "Mixed cohort with lighter clothing and elevated activity relative to the young mixed baseline.",
    },
    "higher_bmi_warm_sensitive": {
        "label": "Higher BMI, Warm Sensitive",
        "description": "Mixed cohort with higher BMI and slightly lighter clothing to test a heavier, warmer-running occupant profile.",
    },
}

OCCUPIED_ROOM_ORDER = [
    "Living Room",
    "Dining Room",
    "Primary Bedroom",
    "Kitchen",
    "Breakfast Nook",
    "Powder Room",
    "Study",
    "Family Lounge West",
    "Family Lounge East",
    "Entry Closet",
    "Service Vestibule",
    "Guest Room",
    "Bath 2",
    "Bedroom 1",
    "Bedroom 2",
    "Bedroom 3",
    "West Bedroom Hall",
    "East Bedroom Hall",
    "Laundry",
    "Ensuite Bath",
    "Shared Bath",
]

ROOM_ACTIVITY = {
    "Living Room": "seated",
    "Dining Room": "seated",
    "Primary Bedroom": "seated",
    "Kitchen": "standing",
    "Breakfast Nook": "seated",
    "Powder Room": "standing",
    "Study": "seated",
    "Family Lounge West": "seated",
    "Family Lounge East": "seated",
    "Entry Closet": "standing",
    "Service Vestibule": "standing",
    "Guest Room": "seated",
    "Bath 2": "standing",
    "Bedroom 1": "seated",
    "Bedroom 2": "seated",
    "Bedroom 3": "seated",
    "West Bedroom Hall": "standing",
    "East Bedroom Hall": "standing",
    "Laundry": "standing",
    "Ensuite Bath": "standing",
    "Shared Bath": "standing",
}

ROOM_CLO = {
    "Kitchen": 0.45,
    "Powder Room": 0.45,
    "Entry Closet": 0.45,
    "Service Vestibule": 0.45,
    "Bath 2": 0.45,
    "Laundry": 0.45,
    "Ensuite Bath": 0.45,
    "Shared Bath": 0.45,
    "default": 0.55,
}

SOUTH_EXPOSED = {
    "Service Vestibule",
    "Guest Room",
    "Family Lounge East",
    "Bath 2",
    "Bedroom 3",
    "Laundry",
    "Bedroom 1",
    "Bedroom 2",
}

EAST_WEST_EXPOSED = {
    "Primary Bedroom",
}

BEDROOM_SPACES = {
    "Primary Bedroom",
    "Guest Room",
    "Bedroom 1",
    "Bedroom 2",
    "Bedroom 3",
}

DAY_SPACES = {
    "Living Room",
    "Dining Room",
    "Breakfast Nook",
    "Study",
    "Family Lounge West",
    "Family Lounge East",
}

BATH_AND_SERVICE = {
    "Bath 2",
    "Ensuite Bath",
    "Shared Bath",
    "Laundry",
}

ROOM_ABBREVIATIONS = {
    "Living Room": "LR",
    "Dining Room": "DR",
    "Primary Bedroom": "PB",
    "Kitchen": "K",
    "Breakfast Nook": "BN",
    "Powder Room": "PR",
    "Study": "ST",
    "Family Lounge West": "FW",
    "Family Lounge East": "FE",
    "Entry Closet": "EC",
    "Service Vestibule": "SV",
    "Guest Room": "GR",
    "Bath 2": "B2A",
    "Bedroom 1": "B1",
    "Bedroom 2": "B2",
    "Bedroom 3": "B3",
    "West Bedroom Hall": "WBH",
    "East Bedroom Hall": "EBH",
    "Laundry": "L",
    "Ensuite Bath": "EB",
    "Shared Bath": "SB",
}

ACTIVITY_LEVELS = ["seated", "standing", "walking"]


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def round_half(value: float) -> float:
    return round(float(value) * 2.0) / 2.0


def temp_label(temp_c: float) -> str:
    return str(temp_c).replace(".", "p")


def representative_cell(cells: List[Tuple[int, int]]) -> List[int]:
    cx = sum(col for col, _ in cells) / len(cells)
    cy = sum(row for _, row in cells) / len(cells)
    chosen = min(
        cells,
        key=lambda cell: ((cell[0] - cx) ** 2 + (cell[1] - cy) ** 2, cell[1], cell[0]),
    )
    return [chosen[0], chosen[1]]


def shift_activity(activity: str, steps: int) -> str:
    try:
        index = ACTIVITY_LEVELS.index(activity)
    except ValueError:
        index = 1
    next_index = max(0, min(len(ACTIVITY_LEVELS) - 1, index + steps))
    return ACTIVITY_LEVELS[next_index]


def sample_latents(room_names: List[str]) -> Dict[str, Dict[str, float]]:
    rng = random.Random(SEED)
    out: Dict[str, Dict[str, float]] = {}
    for room_name in room_names:
        out[room_name] = {
            "sex_roll": rng.random(),
            "age_z": rng.gauss(0.0, 1.0),
            "height_z": rng.gauss(0.0, 1.0),
            "bmi_z": rng.gauss(0.0, 1.0),
        }
    return out


def build_agent_payload(
    cohort_key: str,
    room_name: str,
    position: List[int],
    latent: Dict[str, float],
    agent_index: int,
) -> Dict[str, object]:
    mixed_sex = "male" if latent["sex_roll"] < 0.5 else "female"
    if cohort_key == "young_male":
        sex = "male"
    else:
        sex = mixed_sex

    young_age = clamp(35.0 + 5.0 * latent["age_z"], 25.0, 45.0)
    age = young_age + 50.0 if cohort_key == "older_shift" else young_age
    age = clamp(age, 25.0, 95.0)

    if sex == "male":
        height_cm = clamp(173.0 + 6.0 * latent["height_z"], 160.0, 188.0)
        bmi = clamp(23.2 + 2.4 * latent["bmi_z"], 18.5, 30.0)
    else:
        height_cm = clamp(161.0 + 6.0 * latent["height_z"], 150.0, 176.0)
        bmi = clamp(21.8 + 2.4 * latent["bmi_z"], 18.5, 30.0)

    weight_kg = bmi * (height_cm / 100.0) ** 2
    activity = ROOM_ACTIVITY[room_name]
    clo = ROOM_CLO.get(room_name, ROOM_CLO["default"])

    if cohort_key == "higher_clo_sedentary":
        activity = shift_activity(activity, -1)
        clo += 0.35
    elif cohort_key == "lighter_clothing_mobile":
        activity = shift_activity(activity, 1)
        clo -= 0.20
    elif cohort_key == "higher_bmi_warm_sensitive":
        bmi = clamp(bmi + 4.0, 21.0, 36.0)
        weight_kg = bmi * (height_cm / 100.0) ** 2
        clo -= 0.10

    clo = clamp(clo, 0.30, 1.20)

    return {
        "id": "resident_{0:02d}".format(agent_index),
        "position": position,
        "activity": activity,
        "clo": clo,
        "demographics": {
            "age": int(round(age)),
            "gender": sex,
            "weight_kg": round(weight_kg, 1),
            "height_cm": int(round(height_cm)),
            "mobility": "normal",
            "hearing": "normal",
            "vision": "normal",
        },
        "preferences": {
            "noise_tolerance_db": 60 if cohort_key not in {"older_shift", "higher_bmi_warm_sensitive"} else 55,
            "light_preference": "moderate",
            "social_density": "moderate",
        },
        "room_name": room_name,
        "room_abbrev": ROOM_ABBREVIATIONS[room_name],
    }


def build_room_overrides(room_names: List[str], base_ta: float, base_rh: float) -> Dict[str, Dict[str, float]]:
    overrides: Dict[str, Dict[str, float]] = {}
    for room_name in room_names:
        air_temp = float(base_ta)
        mean_radiant_temp = float(base_ta)
        humidity = float(base_rh)
        air_velocity = 0.10

        if room_name in SOUTH_EXPOSED:
            mean_radiant_temp = base_ta + 1.5
        elif room_name in EAST_WEST_EXPOSED:
            mean_radiant_temp = base_ta + 0.75

        if room_name in BEDROOM_SPACES:
            air_velocity = 0.08
        elif room_name in DAY_SPACES:
            air_velocity = 0.12
        elif room_name == "Kitchen":
            air_velocity = 0.15
        elif room_name in BATH_AND_SERVICE:
            air_velocity = 0.08

        if room_name in BATH_AND_SERVICE:
            humidity = min(90.0, base_rh + 10.0)

        overrides[room_name] = {
            "air_temp": round_half(air_temp),
            "mean_radiant_temp": round_half(mean_radiant_temp),
            "humidity": round_half(humidity),
            "air_velocity": round(air_velocity, 2),
        }
    return overrides


def main() -> None:
    state = load_template_json(PLAN_PATH)
    rooms_by_name = {room.name: room for room in state.plan.rooms}
    missing = [room_name for room_name in OCCUPIED_ROOM_ORDER if room_name not in rooms_by_name]
    if missing:
        raise SystemExit("Missing expected apartment rooms: {0}".format(", ".join(missing)))

    all_room_names = [room.name for room in state.plan.rooms]
    occupied_positions = {
        room_name: representative_cell(sorted(rooms_by_name[room_name].cells))
        for room_name in OCCUPIED_ROOM_ORDER
    }
    latents = sample_latents(OCCUPIED_ROOM_ORDER)

    cohort_agents: Dict[str, List[Dict[str, object]]] = {}
    for cohort_key in COHORTS:
        agents: List[Dict[str, object]] = []
        for index, room_name in enumerate(OCCUPIED_ROOM_ORDER, start=1):
            agents.append(
                build_agent_payload(
                    cohort_key=cohort_key,
                    room_name=room_name,
                    position=occupied_positions[room_name],
                    latent=latents[room_name],
                    agent_index=index,
                )
            )
        cohort_agents[cohort_key] = agents

    scenarios = []
    manifest_scenarios = []
    for cohort_key, cohort_meta in COHORTS.items():
        for base_ta in BASE_TEMPERATURES:
            for base_rh in BASE_HUMIDITIES:
                scenario_id = "{0}_t{1}_rh{2}".format(cohort_key, temp_label(base_ta), int(base_rh))
                room_overrides = build_room_overrides(all_room_names, base_ta, base_rh)
                scenarios.append(
                    {
                        "id": scenario_id,
                        "description": (
                            "{label} at Ta={ta:.1f} C, RH={rh:.0f}% with south-facing rooms at Tr=Ta+1.5 C, "
                            "east/west-exposed rooms at Tr=Ta+0.75 C, north/interior rooms at Tr=Ta."
                        ).format(
                            label=cohort_meta["label"],
                            ta=base_ta,
                            rh=base_rh,
                        ),
                        "use_base_agents": False,
                        "agents": [
                            {
                                "id": agent["id"],
                                "position": agent["position"],
                                "activity": agent["activity"],
                                "clo": agent["clo"],
                                "demographics": agent["demographics"],
                                "preferences": agent["preferences"],
                            }
                            for agent in cohort_agents[cohort_key]
                        ],
                        "room_overrides": room_overrides,
                    }
                )
                manifest_scenarios.append(
                    {
                        "id": scenario_id,
                        "cohort": cohort_key,
                        "cohort_label": cohort_meta["label"],
                        "base_ta_c": base_ta,
                        "base_rh_pct": base_rh,
                    }
                )

    study_payload = {
        "base_plan": "../../floor_plan/apartment_template.json",
        "scenarios": scenarios,
    }

    manifest_payload = {
        "case_id": "apartment_daytime_cooled",
        "title": "Apartment Daytime Cooled Population Sweep",
        "seed": SEED,
        "assumptions": {
            "orientation": "Bottom of plan is south.",
            "operation": "Daytime cooled apartment condition.",
            "thermal_model": (
                "North/interior rooms use Tr=Ta, south-facing rooms use Tr=Ta+1.5 C, "
                "east/west-exposed rooms use Tr=Ta+0.75 C."
            ),
            "cohort_design": (
                "Six cohorts share the same room placement: young mixed, young male-only, age-shifted +50 years, "
                "higher-clo sedentary, lighter-clothing mobile, and higher-BMI warm-sensitive."
            ),
            "occupancy_rule": (
                "One occupant is placed in each enclosed room except Entry Gallery, Service Hall, and East Passage."
            ),
        },
        "plan": {
            "path": str(PLAN_PATH.relative_to(ROOT_DIR)),
            "template_name": state.meta.get("name"),
            "occupied_rooms": OCCUPIED_ROOM_ORDER,
            "room_abbreviations": ROOM_ABBREVIATIONS,
        },
        "cohorts": {
            cohort_key: {
                **COHORTS[cohort_key],
                "agents": cohort_agents[cohort_key],
            }
            for cohort_key in COHORTS
        },
        "scenarios": manifest_scenarios,
    }

    STUDY_PATH.write_text(json.dumps(study_payload, indent=2), encoding="utf-8")
    MANIFEST_PATH.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")

    print("Wrote study: {0}".format(STUDY_PATH))
    print("Wrote manifest: {0}".format(MANIFEST_PATH))
    print("Scenarios: {0}".format(len(scenarios)))
    print("Agents per cohort: {0}".format(len(OCCUPIED_ROOM_ORDER)))


if __name__ == "__main__":
    main()
