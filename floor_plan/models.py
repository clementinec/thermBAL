"""Data models for the floor-plan simulator.

A FloorPlan is a 2-D grid of cells.  Rooms own rectangular (or arbitrary)
sets of cells and carry fixed environmental attributes.  Agents are
movable occupants placed on the grid whose comfort JSON is recomputed
whenever they move.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

Coord = Tuple[int, int]


# ── Environment (per-room, static) ──────────────────────────────────

@dataclass
class RoomEnvironment:
    """Fixed environmental conditions inside a room."""
    air_temp: float = 24.0       # °C
    mean_radiant_temp: float = 24.0  # °C  (defaults to air_temp if unset)
    humidity: float = 50.0       # %RH
    air_velocity: float = 0.1    # m/s
    lux: float = 300.0           # lx
    noise_db: float = 40.0       # dB


# ── Spatial elements ────────────────────────────────────────────────

@dataclass
class Window:
    """A window defined by two grid-edge endpoints."""
    p1: Coord
    p2: Coord
    room_id: str = ""


@dataclass
class Door:
    """A door / exit defined by two grid-edge endpoints."""
    p1: Coord
    p2: Coord
    is_exit: bool = False
    room_id: str = ""


@dataclass
class Room:
    """A named zone owning a set of grid cells."""
    id: str
    name: str
    cells: List[Coord]
    environment: RoomEnvironment = field(default_factory=RoomEnvironment)
    ceiling_height: float = 2.8  # m


# ── Floor plan ──────────────────────────────────────────────────────

@dataclass
class FloorPlan:
    """The stage: a grid with rooms, walls, windows, and doors."""
    grid_cols: int                # X dimension (width in cells)
    grid_rows: int                # Y dimension (depth in cells)
    cell_size_m: float = 1.0     # physical size of one cell edge
    rooms: List[Room] = field(default_factory=list)
    windows: List[Window] = field(default_factory=list)
    doors: List[Door] = field(default_factory=list)

    # ── lookups (built lazily) ──
    _cell_to_room: Dict[Coord, Room] = field(
        default_factory=dict, repr=False, init=False
    )

    def build_index(self) -> None:
        """Rebuild the cell → room lookup."""
        self._cell_to_room = {}
        for room in self.rooms:
            for c in room.cells:
                self._cell_to_room[c] = room

    def room_at(self, coord: Coord) -> Optional[Room]:
        if not self._cell_to_room:
            self.build_index()
        return self._cell_to_room.get(coord)

    def is_inside(self, coord: Coord) -> bool:
        return 0 <= coord[0] < self.grid_cols and 0 <= coord[1] < self.grid_rows


# ── Agent / persona ─────────────────────────────────────────────────

@dataclass
class AgentDemographics:
    age: int = 30
    gender: str = "male"
    weight_kg: float = 70.0
    height_cm: float = 170.0
    mobility: str = "normal"         # normal | uses_cane | uses_walker | wheelchair
    hearing: str = "normal"          # normal | impaired | deaf
    vision: str = "normal"           # normal | mild_impairment | severe_impairment | blind


@dataclass
class AgentPreferences:
    """Subjective preferences (feed into LLM layer later)."""
    noise_tolerance_db: float = 55.0
    light_preference: str = "bright"   # dim | moderate | bright
    social_density: str = "moderate"   # low_den | moderate | high_den


@dataclass
class Agent:
    """A movable occupant on the floor plan."""
    id: str
    demographics: AgentDemographics = field(default_factory=AgentDemographics)
    preferences: AgentPreferences = field(default_factory=AgentPreferences)
    clo: float = 1.0               # clothing insulation
    activity: str = "seated"       # seated | standing | walking
    position: Coord = (0, 0)       # current grid cell
