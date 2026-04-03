"""Spatial computations on the floor-plan grid.

Given an agent's position, compute distances to walls, windows, doors/exits,
enclosure ratio, and count of visible neighbouring agents.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from .models import Agent, Coord, Door, FloorPlan, Room, Window


def _cell_centre(coord: Coord, cell_size: float) -> Tuple[float, float]:
    """Return the physical centre of a grid cell."""
    return (coord[0] + 0.5) * cell_size, (coord[1] + 0.5) * cell_size


def _seg_midpoint(p1: Coord, p2: Coord, cell_size: float) -> Tuple[float, float]:
    """Midpoint of an element defined by two grid points."""
    return (p1[0] + p2[0]) / 2.0 * cell_size, (p1[1] + p2[1]) / 2.0 * cell_size


def _dist(a: Tuple[float, float], b: Tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


# ── Distance to nearest room boundary (wall) ───────────────────────

def dist_to_wall(coord: Coord, room: Room, cell_size: float) -> float:
    """Minimum distance from cell centre to the room perimeter.

    We compute this by finding the closest cell that is *in* the room
    whose neighbour is *outside* the room (i.e. a boundary cell), then
    taking the distance.  For rectangular rooms this simplifies to
    min distance to any edge.
    """
    cx, cy = _cell_centre(coord, cell_size)
    room_set = set(room.cells)

    min_d = float("inf")
    for cell in room.cells:
        # check if this cell is on the perimeter
        x, y = cell
        neighbours = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
        if any(n not in room_set for n in neighbours):
            bx, by = _cell_centre(cell, cell_size)
            d = _dist((cx, cy), (bx, by))
            # the wall is at the far edge of the boundary cell
            # approximate: distance minus half a cell gets you to the wall
            wall_d = max(0.0, d - cell_size * 0.0)  # keep centre-to-centre for now
            if cell == coord:
                # agent IS on a boundary cell → distance is to nearest edge
                wall_d = cell_size * 0.5
            min_d = min(min_d, wall_d)

    return round(min_d, 2) if min_d < float("inf") else 0.0


def dist_to_nearest_window(coord: Coord, plan: FloorPlan) -> float:
    """Distance from cell centre to the nearest window midpoint."""
    cx, cy = _cell_centre(coord, plan.cell_size_m)
    if not plan.windows:
        return -1.0  # no windows
    return round(
        min(_dist((cx, cy), _seg_midpoint(w.p1, w.p2, plan.cell_size_m))
            for w in plan.windows),
        2,
    )


def dist_to_nearest_door(coord: Coord, plan: FloorPlan, exits_only: bool = False) -> float:
    """Distance from cell centre to the nearest door (or exit)."""
    cx, cy = _cell_centre(coord, plan.cell_size_m)
    doors = [d for d in plan.doors if (not exits_only or d.is_exit)]
    if not doors:
        return -1.0
    return round(
        min(_dist((cx, cy), _seg_midpoint(d.p1, d.p2, plan.cell_size_m))
            for d in doors),
        2,
    )


# ── Enclosure ratio ────────────────────────────────────────────────

def enclosure_ratio(coord: Coord, room: Room, plan: FloorPlan) -> float:
    """Fraction of the 4 cardinal neighbours that are walls / outside the room.

    0.0 = completely open (all 4 sides are room interior)
    1.0 = fully enclosed (all 4 sides are walls / out-of-room)
    """
    x, y = coord
    room_set = set(room.cells)
    blocked = 0
    for nx, ny in [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]:
        if (nx, ny) not in room_set or not plan.is_inside((nx, ny)):
            blocked += 1
    return round(blocked / 4.0, 2)


# ── Visible agents ──────────────────────────────────────────────────

def visible_agents(
    agent: Agent,
    all_agents: List[Agent],
    room: Room,
    max_range_cells: int = 10,
) -> int:
    """Count other agents in the same room within line-of-sight range."""
    room_set = set(room.cells)
    count = 0
    ax, ay = agent.position
    for other in all_agents:
        if other.id == agent.id:
            continue
        if other.position not in room_set:
            continue
        ox, oy = other.position
        if abs(ox - ax) + abs(oy - ay) <= max_range_cells:
            count += 1
    return count


# ── Aggregate spatial snapshot ──────────────────────────────────────

def compute_spatial(
    agent: Agent,
    plan: FloorPlan,
    all_agents: List[Agent],
) -> Dict:
    """Return the full spatial block for an agent's current position."""
    room = plan.room_at(agent.position)
    if room is None:
        return {"error": f"Agent {agent.id} at {agent.position} is not inside any room"}

    return {
        "coordinates": list(agent.position),
        "room_id": room.id,
        "room_name": room.name,
        "dist_to_wall": dist_to_wall(agent.position, room, plan.cell_size_m),
        "dist_to_window": dist_to_nearest_window(agent.position, plan),
        "dist_to_exit": dist_to_nearest_door(agent.position, plan, exits_only=True),
        "dist_to_door": dist_to_nearest_door(agent.position, plan, exits_only=False),
        "ceiling_h": room.ceiling_height,
        "enclosure_ratio": enclosure_ratio(agent.position, room, plan),
        "visible_agents": visible_agents(agent, all_agents, room),
    }
