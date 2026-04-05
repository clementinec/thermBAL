"""Floor-plan experience simulator.

The Sims-like engine for agent-based environmental experience modelling.
Drop agents onto a floor plan, move them around, and get rich JSON
snapshots of their computed comfort and spatial context.
"""

from .engine import Simulator
from .io import (
    LoadedPlanState,
    build_plan_state_payload,
    load_geometry_json,
    load_plan_state_json,
    load_template_json,
    load_template_payload,
    save_plan_state_json,
)
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
from .study import (
    AgentSpec,
    Scenario,
    ScenarioResult,
    ScheduledEvent,
    SnapshotRecord,
    StudyDefinition,
    load_study_definition,
    run_scenario,
    run_study,
)

__all__ = [
    "Agent",
    "AgentDemographics",
    "AgentPreferences",
    "AgentSpec",
    "Coord",
    "Door",
    "FloorPlan",
    "LoadedPlanState",
    "Room",
    "RoomEnvironment",
    "Scenario",
    "ScenarioResult",
    "ScheduledEvent",
    "Simulator",
    "SnapshotRecord",
    "StudyDefinition",
    "Window",
    "build_plan_state_payload",
    "compute_spatial",
    "default_floor_plan",
    "demo_simulator",
    "elder_walker",
    "load_geometry_json",
    "load_plan_state_json",
    "load_study_definition",
    "load_template_json",
    "load_template_payload",
    "run_scenario",
    "run_study",
    "save_plan_state_json",
    "seated_elder",
    "young_visitor",
]
