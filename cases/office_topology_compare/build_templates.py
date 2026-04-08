#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


CASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = CASE_DIR.parents[1]
SOURCE_TEMPLATE = ROOT_DIR / "floor_plan" / "apartment_template.json"
CELLULAR_TEMPLATE = ROOT_DIR / "floor_plan" / "office_cellular_template.json"
OPEN_TEMPLATE = ROOT_DIR / "floor_plan" / "office_openplan_template.json"

Coord = Tuple[int, int]

CELLULAR_ROOM_RENAMES = {
    "Living Room": "Team Office West",
    "Dining Room": "Team Office East",
    "Primary Bedroom": "Executive Office",
    "Kitchen": "Reception",
    "Breakfast Nook": "Admin Bay",
    "Powder Room": "North WC",
    "Study": "Focus Room",
    "Entry Gallery": "Core Lobby",
    "Entry Closet": "Storage",
    "Service Vestibule": "Print Hub",
    "Service Hall": "West Connector",
    "Family Lounge West": "Collaboration West",
    "Family Lounge East": "Collaboration East",
    "East Passage": "East Connector",
    "Guest Room": "Meeting Room West",
    "Bath 2": "East WC",
    "Bedroom 3": "Meeting Room East",
    "Laundry": "Utility",
    "Ensuite Bath": "Wellness Room",
    "West Bedroom Hall": "South Studio West",
    "East Bedroom Hall": "South Studio East",
    "Bedroom 1": "Project Bay West",
    "Bedroom 2": "Project Bay East",
    "Shared Bath": "Core WC",
}

OPEN_GROUPS = {
    "North Open Office": {
        "Kitchen",
        "Breakfast Nook",
        "Living Room",
        "Dining Room",
        "Study",
        "Primary Bedroom",
    },
    "South Open Office": {
        "Family Lounge West",
        "Family Lounge East",
        "Guest Room",
        "Bedroom 3",
        "West Bedroom Hall",
        "East Bedroom Hall",
        "Bedroom 1",
        "Bedroom 2",
    },
    "Core Lobby": {
        "Entry Gallery",
        "Service Hall",
        "East Passage",
    },
    "West Corner Office": {
        "Entry Closet",
        "Service Vestibule",
        "Laundry",
    },
    "North WC": {"Powder Room"},
    "East WC": {"Bath 2"},
    "Wellness Room": {"Ensuite Bath"},
    "Core WC": {"Shared Bath"},
}


def load_payload() -> Dict[str, object]:
    return json.loads(SOURCE_TEMPLATE.read_text(encoding="utf-8"))


def room_cells_by_name(payload: Dict[str, object]) -> Dict[str, Set[Coord]]:
    out: Dict[str, Set[Coord]] = {}
    for room in payload.get("rooms", []):
        if not isinstance(room, dict):
            continue
        out[str(room["name"])] = {tuple(cell) for cell in room.get("cells", [])}
    return out


def room_payloads_by_name(payload: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for room in payload.get("rooms", []):
        if isinstance(room, dict):
            out[str(room["name"])] = room
    return out


def edge_neighbors(orient: str, row: int, col: int) -> Tuple[Coord, Coord]:
    if orient == "h":
        return (col, row - 1), (col, row)
    return (col - 1, row), (col, row)


def edge_is_exterior(orient: str, row: int, col: int, active_cells: Set[Coord]) -> bool:
    a, b = edge_neighbors(orient, row, col)
    return (a in active_cells) ^ (b in active_cells)


def build_walls(active_cells: Set[Coord], room_by_cell: Dict[Coord, str], cols: int, rows: int) -> Tuple[List[List[bool]], List[List[bool]]]:
    h_walls = [[False for _ in range(cols)] for _ in range(rows + 1)]
    v_walls = [[False for _ in range(cols + 1)] for _ in range(rows)]

    for col, row in active_cells:
        top = (col, row - 1)
        bottom = (col, row + 1)
        left = (col - 1, row)
        right = (col + 1, row)

        if top not in active_cells or room_by_cell.get(top) != room_by_cell[(col, row)]:
            h_walls[row][col] = True
        if bottom not in active_cells or room_by_cell.get(bottom) != room_by_cell[(col, row)]:
            h_walls[row + 1][col] = True
        if left not in active_cells or room_by_cell.get(left) != room_by_cell[(col, row)]:
            v_walls[row][col] = True
        if right not in active_cells or room_by_cell.get(right) != room_by_cell[(col, row)]:
            v_walls[row][col + 1] = True

    return h_walls, v_walls


def wall_exists(feature: Dict[str, object], h_walls: List[List[bool]], v_walls: List[List[bool]]) -> bool:
    orient = str(feature["orient"])
    row = int(feature["r"])
    col = int(feature["c"])
    if orient == "h":
        return 0 <= row < len(h_walls) and 0 <= col < len(h_walls[row]) and h_walls[row][col]
    return 0 <= row < len(v_walls) and 0 <= col < len(v_walls[row]) and v_walls[row][col]


def filter_open_features(features: Iterable[Dict[str, object]], active_cells: Set[Coord], h_walls: List[List[bool]], v_walls: List[List[bool]]) -> List[Dict[str, object]]:
    kept: List[Dict[str, object]] = []
    for feature in features:
        orient = str(feature.get("orient", ""))
        row = int(feature.get("r", -1))
        col = int(feature.get("c", -1))
        kind = str(feature.get("type", "window"))
        if orient not in {"h", "v"}:
            continue
        if kind not in {"window", "door", "exit"}:
            continue
        if not edge_is_exterior(orient, row, col, active_cells):
            continue
        if not wall_exists(feature, h_walls, v_walls):
            continue
        if kind == "window" or kind == "exit":
            kept.append({"type": kind, "orient": orient, "r": row, "c": col})
    return kept


def cellular_payload(source: Dict[str, object]) -> Dict[str, object]:
    payload = copy.deepcopy(source)
    payload["meta"] = {
        **payload.get("meta", {}),
        "name": "Office Cellular Template",
        "mode": "geometry_template",
        "source_template": "apartment_template.json",
        "derived_case": "office_topology_compare",
    }
    rooms: List[Dict[str, object]] = []
    for room in payload.get("rooms", []):
        if not isinstance(room, dict):
            continue
        renamed = copy.deepcopy(room)
        original_name = str(room["name"])
        renamed["name"] = CELLULAR_ROOM_RENAMES[original_name]
        rooms.append(renamed)
    payload["rooms"] = rooms
    return payload


def open_payload(source: Dict[str, object]) -> Dict[str, object]:
    room_cells = room_cells_by_name(source)
    room_payloads = room_payloads_by_name(source)
    active_cells = {tuple(cell) for cell in source.get("active_cells", [])}
    meta = dict(source.get("meta", {}))
    cols = int(meta["grid_cols"])
    rows = int(meta["grid_rows"])

    merged_rooms: List[Dict[str, object]] = []
    room_by_cell: Dict[Coord, str] = {}
    for index, (group_name, member_names) in enumerate(OPEN_GROUPS.items(), start=1):
        cells: Set[Coord] = set()
        for member_name in member_names:
            cells.update(room_cells[member_name])
        for cell in cells:
            room_by_cell[cell] = group_name
        source_room = room_payloads[next(iter(member_names))]
        merged_rooms.append(
            {
                "id": f"office_open_{index:02d}",
                "name": group_name,
                "cells": [list(cell) for cell in sorted(cells, key=lambda item: (item[1], item[0]))],
                "env": copy.deepcopy(source_room.get("env", {})),
                "ceilH": float(source_room.get("ceilH", 2.8)),
            }
        )

    if set(room_by_cell) != active_cells:
        missing = active_cells.difference(room_by_cell)
        raise RuntimeError(f"Open-office grouping does not cover all active cells: {sorted(missing)[:10]}")

    h_walls, v_walls = build_walls(active_cells, room_by_cell, cols, rows)
    features = filter_open_features(source.get("features", []), active_cells, h_walls, v_walls)

    return {
        "meta": {
            **meta,
            "name": "Office Open-Plan Template",
            "mode": "geometry_template",
            "source_template": "apartment_template.json",
            "derived_case": "office_topology_compare",
        },
        "pdf_reference": copy.deepcopy(source.get("pdf_reference")),
        "active_cells": [list(cell) for cell in sorted(active_cells, key=lambda item: (item[1], item[0]))],
        "walls": {
            "horizontal": h_walls,
            "vertical": v_walls,
        },
        "features": features,
        "rooms": merged_rooms,
    }


def derive_templates() -> Tuple[Path, Path]:
    source = load_payload()
    CELLULAR_TEMPLATE.write_text(json.dumps(cellular_payload(source), indent=2), encoding="utf-8")
    OPEN_TEMPLATE.write_text(json.dumps(open_payload(source), indent=2), encoding="utf-8")
    return CELLULAR_TEMPLATE, OPEN_TEMPLATE


def main() -> None:
    cellular_path, open_path = derive_templates()
    print(f"Wrote {cellular_path}")
    print(f"Wrote {open_path}")


if __name__ == "__main__":
    main()
