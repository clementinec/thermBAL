"""JSON IO helpers for thermBAL floor-plan studies.

These adapters bridge the browser-exported JSON payloads and the Python
floor-plan models used by the batch runner.
"""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .models import (
    Agent,
    AgentDemographics,
    AgentPreferences,
    Coord,
    Door,
    FloorPlan,
    Room,
    RoomEnvironment,
    Window,
)

DEFAULT_ENV = {
    "air_temp": 24.0,
    "mean_radiant_temp": 24.0,
    "humidity": 50.0,
    "air_velocity": 0.1,
    "lux": 300.0,
    "noise_db": 40.0,
}
DEFAULT_CEIL_H = 2.8


@dataclass
class LoadedPlanState:
    """Loaded floor-plan state plus the raw grid topology."""

    plan: FloorPlan
    agents: List[Agent] = field(default_factory=list)
    meta: Dict[str, Any] = field(default_factory=dict)
    active_cells: Set[Coord] = field(default_factory=set, repr=False)
    horizontal_walls: List[List[bool]] = field(default_factory=list, repr=False)
    vertical_walls: List[List[bool]] = field(default_factory=list, repr=False)
    features: List[Dict[str, Any]] = field(default_factory=list, repr=False)

    @classmethod
    def from_floor_plan(
        cls,
        plan: FloorPlan,
        agents: Optional[Iterable[Agent]] = None,
        meta: Optional[Dict[str, Any]] = None,
    ) -> "LoadedPlanState":
        """Build a serializable plan state from Python floor-plan objects."""
        plan_copy = copy.deepcopy(plan)
        agents_copy = [copy.deepcopy(agent) for agent in (agents or [])]
        active_cells, h_walls, v_walls = _derive_grid_topology(plan_copy)
        features = _serialize_plan_features(plan_copy)
        return cls(
            plan=plan_copy,
            agents=agents_copy,
            meta={
                "mode": "simulator_state",
                "grid_cols": plan_copy.grid_cols,
                "grid_rows": plan_copy.grid_rows,
                "cell_size_m": plan_copy.cell_size_m,
                **(meta or {}),
            },
            active_cells=active_cells,
            horizontal_walls=h_walls,
            vertical_walls=v_walls,
            features=features,
        )


def load_json_payload(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_geometry_json(path: str | Path) -> LoadedPlanState:
    return load_template_json(path)


def load_plan_state_json(path: str | Path) -> LoadedPlanState:
    return load_template_json(path)


def load_template_json(path: str | Path) -> LoadedPlanState:
    payload = load_json_payload(path)
    return load_template_payload(payload)


def load_template_payload(payload: Dict[str, Any]) -> LoadedPlanState:
    if not isinstance(payload, dict):
        raise TypeError("Plan payload must be a JSON object.")

    meta = payload.get("meta") or {}
    walls = payload.get("walls") or {}
    grid_cols = int(round(meta.get("grid_cols", payload.get("grid_cols", 20))))
    grid_rows = int(round(meta.get("grid_rows", payload.get("grid_rows", 15))))
    cell_size_m = float(meta.get("cell_size_m", payload.get("cell_size_m", 1.0)))

    active_cells = _sanitize_active_cells(payload.get("active_cells"), grid_cols, grid_rows)
    h_walls = _sanitize_wall_matrix(walls.get("horizontal"), grid_rows + 1, grid_cols)
    v_walls = _sanitize_wall_matrix(walls.get("vertical"), grid_rows, grid_cols + 1)
    features = _sanitize_features(payload.get("features"), grid_cols, grid_rows)

    rooms = _build_rooms(payload.get("rooms"), active_cells, h_walls, v_walls)
    plan = FloorPlan(
        grid_cols=grid_cols,
        grid_rows=grid_rows,
        cell_size_m=cell_size_m,
        rooms=rooms,
    )
    plan.windows, plan.doors = _build_openings(features, plan, active_cells)
    plan.build_index()

    agents = [_agent_from_payload(item) for item in payload.get("agents", []) or []]

    return LoadedPlanState(
        plan=plan,
        agents=agents,
        meta=dict(meta),
        active_cells=active_cells,
        horizontal_walls=h_walls,
        vertical_walls=v_walls,
        features=features,
    )


def build_plan_state_payload(state: LoadedPlanState, include_agents: bool = True) -> Dict[str, Any]:
    payload = {
        "meta": {
            "mode": state.meta.get("mode", "simulator_state"),
            "grid_cols": state.plan.grid_cols,
            "grid_rows": state.plan.grid_rows,
            "cell_size_m": state.plan.cell_size_m,
            **{k: v for k, v in state.meta.items() if k not in {"grid_cols", "grid_rows", "cell_size_m"}},
        },
        "active_cells": [[col, row] for col, row in sorted(state.active_cells)],
        "walls": {
            "horizontal": [row[:] for row in state.horizontal_walls],
            "vertical": [row[:] for row in state.vertical_walls],
        },
        "features": [dict(feature) for feature in state.features],
        "rooms": [
            {
                "id": room.id,
                "name": room.name,
                "cells": [list(cell) for cell in room.cells],
                "env": {
                    "air_temp": room.environment.air_temp,
                    "mean_radiant_temp": room.environment.mean_radiant_temp,
                    "humidity": room.environment.humidity,
                    "air_velocity": room.environment.air_velocity,
                    "lux": room.environment.lux,
                    "noise_db": room.environment.noise_db,
                },
                "ceilH": room.ceiling_height,
            }
            for room in state.plan.rooms
        ],
    }
    if include_agents:
        payload["agents"] = [_serialize_agent(agent) for agent in state.agents]
    return payload


def save_plan_state_json(
    path: str | Path,
    state: LoadedPlanState,
    include_agents: bool = True,
) -> None:
    payload = build_plan_state_payload(state, include_agents=include_agents)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _sanitize_active_cells(cells: Any, cols: int, rows: int) -> Set[Coord]:
    active: Set[Coord] = set()
    if not isinstance(cells, list):
        return active
    for cell in cells:
        if not isinstance(cell, (list, tuple)) or len(cell) < 2:
            continue
        col = int(round(cell[0]))
        row = int(round(cell[1]))
        if 0 <= col < cols and 0 <= row < rows:
            active.add((col, row))
    return active


def _sanitize_wall_matrix(matrix: Any, rows: int, cols: int) -> List[List[bool]]:
    out = [[False for _ in range(cols)] for _ in range(rows)]
    if not isinstance(matrix, list):
        return out
    for row_idx in range(min(rows, len(matrix))):
        row = matrix[row_idx]
        if not isinstance(row, list):
            continue
        for col_idx in range(min(cols, len(row))):
            out[row_idx][col_idx] = bool(row[col_idx])
    return out


def _sanitize_features(features: Any, cols: int, rows: int) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    if not isinstance(features, list):
        return out
    for feature in features:
        if not isinstance(feature, dict):
            continue
        orient = feature.get("orient")
        if orient not in {"h", "v"}:
            continue
        row = int(round(feature.get("r", -1)))
        col = int(round(feature.get("c", -1)))
        if orient == "h" and not (0 <= row <= rows and 0 <= col < cols):
            continue
        if orient == "v" and not (0 <= row < rows and 0 <= col <= cols):
            continue
        out.append(
            {
                "type": feature.get("type", "window") if feature.get("type") in {"window", "door", "exit"} else "window",
                "orient": orient,
                "r": row,
                "c": col,
            }
        )
    return out


def _room_environment_from_payload(payload: Optional[Dict[str, Any]]) -> RoomEnvironment:
    data = {**DEFAULT_ENV, **(payload or {})}
    return RoomEnvironment(
        air_temp=float(data["air_temp"]),
        mean_radiant_temp=float(data["mean_radiant_temp"]),
        humidity=float(data["humidity"]),
        air_velocity=float(data["air_velocity"]),
        lux=float(data["lux"]),
        noise_db=float(data["noise_db"]),
    )


def _build_rooms(
    payload_rooms: Any,
    active_cells: Set[Coord],
    h_walls: List[List[bool]],
    v_walls: List[List[bool]],
) -> List[Room]:
    rooms: List[Room] = []
    covered: Set[Coord] = set()

    if isinstance(payload_rooms, list):
        for index, room_payload in enumerate(payload_rooms):
            if not isinstance(room_payload, dict):
                continue
            cells = []
            for cell in room_payload.get("cells", []) or []:
                if not isinstance(cell, (list, tuple)) or len(cell) < 2:
                    continue
                coord = (int(round(cell[0])), int(round(cell[1])))
                if coord in active_cells:
                    cells.append(coord)
            if not cells:
                continue
            room = Room(
                id=str(room_payload.get("id") or "room_{0}".format(index)),
                name=str(room_payload.get("name") or "Room {0}".format(_room_label(index))),
                cells=sorted(set(cells)),
                environment=_room_environment_from_payload(room_payload.get("env")),
                ceiling_height=float(room_payload.get("ceilH", DEFAULT_CEIL_H)),
            )
            rooms.append(room)
            covered.update(room.cells)

    if covered != active_cells:
        detected = _detect_room_components(active_cells - covered, h_walls, v_walls)
        offset = len(rooms)
        for index, component in enumerate(detected):
            rooms.append(
                Room(
                    id="room_{0}".format(offset + index),
                    name="Room {0}".format(_room_label(offset + index)),
                    cells=[tuple(cell) for cell in sorted(component)],
                    environment=_room_environment_from_payload(None),
                    ceiling_height=DEFAULT_CEIL_H,
                )
            )

    return rooms


def _detect_room_components(
    active_cells: Set[Coord],
    h_walls: List[List[bool]],
    v_walls: List[List[bool]],
) -> List[Set[Coord]]:
    visited: Set[Coord] = set()
    components: List[Set[Coord]] = []

    for start in sorted(active_cells):
        if start in visited:
            continue
        queue = [start]
        visited.add(start)
        component: Set[Coord] = set()

        while queue:
            col, row = queue.pop(0)
            component.add((col, row))
            for next_col, next_row in (
                (col - 1, row),
                (col + 1, row),
                (col, row - 1),
                (col, row + 1),
            ):
                next_coord = (next_col, next_row)
                if next_coord not in active_cells or next_coord in visited:
                    continue
                if _wall_between((col, row), next_coord, h_walls, v_walls):
                    continue
                visited.add(next_coord)
                queue.append(next_coord)

        components.append(component)

    return components


def _wall_between(
    left: Coord,
    right: Coord,
    h_walls: List[List[bool]],
    v_walls: List[List[bool]],
) -> bool:
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        wall_row = max(y1, y2)
        return bool(h_walls[wall_row][x1])
    if y1 == y2:
        wall_col = max(x1, x2)
        return bool(v_walls[y1][wall_col])
    return True


def _build_openings(
    features: List[Dict[str, Any]],
    plan: FloorPlan,
    active_cells: Set[Coord],
) -> Tuple[List[Window], List[Door]]:
    windows: List[Window] = []
    doors: List[Door] = []

    for feature in _group_features(features):
        room_id = _feature_room_id(feature, plan, active_cells)
        if feature["type"] == "window":
            windows.append(Window(p1=feature["p1"], p2=feature["p2"], room_id=room_id))
        else:
            doors.append(
                Door(
                    p1=feature["p1"],
                    p2=feature["p2"],
                    is_exit=feature["type"] == "exit",
                    room_id=room_id,
                )
            )

    return windows, doors


def _group_features(features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: List[Dict[str, Any]] = []

    for feature_type in ("window", "door", "exit"):
        for orient in ("h", "v"):
            by_line: Dict[int, List[int]] = {}
            for feature in features:
                if feature["type"] != feature_type or feature["orient"] != orient:
                    continue
                line_key = feature["r"] if orient == "h" else feature["c"]
                step = feature["c"] if orient == "h" else feature["r"]
                by_line.setdefault(line_key, []).append(step)

            for line_key, steps in by_line.items():
                steps = sorted(set(steps))
                if not steps:
                    continue
                start = steps[0]
                end = start
                for step in steps[1:]:
                    if step == end + 1:
                        end = step
                        continue
                    grouped.append(_grouped_feature_payload(feature_type, orient, line_key, start, end))
                    start = step
                    end = step
                grouped.append(_grouped_feature_payload(feature_type, orient, line_key, start, end))

    return grouped


def _grouped_feature_payload(feature_type: str, orient: str, line_key: int, start: int, end: int) -> Dict[str, Any]:
    if orient == "h":
        return {
            "type": feature_type,
            "orient": orient,
            "p1": (start, line_key),
            "p2": (end + 1, line_key),
        }
    return {
        "type": feature_type,
        "orient": orient,
        "p1": (line_key, start),
        "p2": (line_key, end + 1),
    }


def _feature_room_id(feature: Dict[str, Any], plan: FloorPlan, active_cells: Set[Coord]) -> str:
    coords = _adjacent_cells_for_segment(feature["p1"], feature["p2"], feature["orient"])
    for coord in coords:
        if coord not in active_cells:
            continue
        room = plan.room_at(coord)
        if room is not None:
            return room.id
    return ""


def _adjacent_cells_for_segment(p1: Coord, p2: Coord, orient: str) -> List[Coord]:
    cells: List[Coord] = []
    if orient == "h":
        row = p1[1]
        for col in range(min(p1[0], p2[0]), max(p1[0], p2[0])):
            cells.append((col, row - 1))
            cells.append((col, row))
    else:
        col = p1[0]
        for row in range(min(p1[1], p2[1]), max(p1[1], p2[1])):
            cells.append((col - 1, row))
            cells.append((col, row))
    return cells


def _agent_from_payload(payload: Dict[str, Any]) -> Agent:
    demo = payload.get("demo") or payload.get("demographics") or {}
    prefs = payload.get("prefs") or payload.get("preferences") or {}
    position = payload.get("pos") or payload.get("position") or [0, 0]
    return Agent(
        id=str(payload.get("id") or "agent"),
        demographics=AgentDemographics(
            age=int(demo.get("age", 30)),
            gender=str(demo.get("gender", "male")),
            weight_kg=float(demo.get("weight_kg", 70.0)),
            height_cm=float(demo.get("height_cm", 170.0)),
            mobility=str(demo.get("mobility", "normal")),
            hearing=str(demo.get("hearing", "normal")),
            vision=str(demo.get("vision", "normal")),
        ),
        preferences=AgentPreferences(
            noise_tolerance_db=float(prefs.get("noise_tolerance_db", 55.0)),
            light_preference=str(prefs.get("light_preference", "bright")),
            social_density=str(prefs.get("social_density", "moderate")),
        ),
        clo=float(payload.get("clo", 1.0)),
        activity=str(payload.get("activity", "seated")),
        position=(int(position[0]), int(position[1])),
    )


def _serialize_agent(agent: Agent) -> Dict[str, Any]:
    return {
        "id": agent.id,
        "demo": {
            "age": agent.demographics.age,
            "gender": agent.demographics.gender,
            "weight_kg": agent.demographics.weight_kg,
            "height_cm": agent.demographics.height_cm,
            "mobility": agent.demographics.mobility,
            "hearing": agent.demographics.hearing,
            "vision": agent.demographics.vision,
        },
        "prefs": {
            "noise_tolerance_db": agent.preferences.noise_tolerance_db,
            "light_preference": agent.preferences.light_preference,
            "social_density": agent.preferences.social_density,
        },
        "clo": agent.clo,
        "activity": agent.activity,
        "pos": list(agent.position),
    }


def _derive_grid_topology(plan: FloorPlan) -> Tuple[Set[Coord], List[List[bool]], List[List[bool]]]:
    active_cells: Set[Coord] = set()
    room_by_cell: Dict[Coord, str] = {}
    for room in plan.rooms:
        for coord in room.cells:
            cell = (int(coord[0]), int(coord[1]))
            active_cells.add(cell)
            room_by_cell[cell] = room.id

    h_walls = [[False for _ in range(plan.grid_cols)] for _ in range(plan.grid_rows + 1)]
    v_walls = [[False for _ in range(plan.grid_cols + 1)] for _ in range(plan.grid_rows)]

    for col, row in sorted(active_cells):
        room_id = room_by_cell[(col, row)]

        north = (col, row - 1)
        south = (col, row + 1)
        west = (col - 1, row)
        east = (col + 1, row)

        if north not in active_cells or room_by_cell.get(north) != room_id:
            h_walls[row][col] = True
        if south not in active_cells or room_by_cell.get(south) != room_id:
            h_walls[row + 1][col] = True
        if west not in active_cells or room_by_cell.get(west) != room_id:
            v_walls[row][col] = True
        if east not in active_cells or room_by_cell.get(east) != room_id:
            v_walls[row][col + 1] = True

    return active_cells, h_walls, v_walls


def _serialize_plan_features(plan: FloorPlan) -> List[Dict[str, Any]]:
    features: List[Dict[str, Any]] = []
    for window in plan.windows:
        features.extend(_segment_to_feature_edges(window.p1, window.p2, "window"))
    for door in plan.doors:
        feature_type = "exit" if door.is_exit else "door"
        features.extend(_segment_to_feature_edges(door.p1, door.p2, feature_type))
    return features


def _segment_to_feature_edges(p1: Coord, p2: Coord, feature_type: str) -> List[Dict[str, Any]]:
    x1, y1 = p1
    x2, y2 = p2
    edges: List[Dict[str, Any]] = []
    if y1 == y2:
        row = y1
        for col in range(min(x1, x2), max(x1, x2)):
            edges.append({"type": feature_type, "orient": "h", "r": row, "c": col})
    elif x1 == x2:
        col = x1
        for row in range(min(y1, y2), max(y1, y2)):
            edges.append({"type": feature_type, "orient": "v", "r": row, "c": col})
    return edges


def _room_label(index: int) -> str:
    base = chr(65 + (index % 26))
    suffix = index // 26
    return base if suffix == 0 else "{0}{1}".format(base, suffix)
