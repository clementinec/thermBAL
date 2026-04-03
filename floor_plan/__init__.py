"""Floor-plan experience simulator.

The Sims-like engine for agent-based environmental experience modelling.
Drop agents onto a floor plan, move them around, and get rich JSON
snapshots of their computed comfort and spatial context.
"""

from .engine import Simulator
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
from .presets import (
    default_floor_plan,
    demo_simulator,
    elder_walker,
    seated_elder,
    young_visitor,
)
from .spatial import compute_spatial

__all__ = [
    "Agent",
    "AgentDemographics",
    "AgentPreferences",
    "Coord",
    "Door",
    "FloorPlan",
    "Room",
    "RoomEnvironment",
    "Simulator",
    "Window",
    "compute_spatial",
    "default_floor_plan",
    "demo_simulator",
    "elder_walker",
    "seated_elder",
    "young_visitor",
]
