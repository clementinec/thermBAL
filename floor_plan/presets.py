"""Ready-made floor plans and agent templates.

These let you spin up a simulation instantly:

    from floor_plan import presets
    sim = presets.demo_simulator()
    print(sim.snapshot_json("persona_01"))
"""

from __future__ import annotations

from .models import (
    Agent,
    AgentDemographics,
    AgentPreferences,
    Door,
    FloorPlan,
    Room,
    RoomEnvironment,
    Window,
)


# ── Helper: rectangular room cells ──────────────────────────────────

def _rect_cells(x0: int, y0: int, w: int, h: int):
    """Generate cell list for a rectangular room."""
    return [(x, y) for x in range(x0, x0 + w) for y in range(y0, y0 + h)]


# ── Default floor plan ──────────────────────────────────────────────
# A simple clinic / community-centre layout:
#
#   0 1 2 3 4 5 6 7 8 9 10 11
# 0 +-----------+------------+
# 1 |           |            |
# 2 |   LOBBY   | ACTIVITY   |
# 3 |           |   ROOM     |
# 4 |     D     |            |
# 5 +-----+----+----+-------+
# 6 |CORR |         |       |
# 7 |IDOR | DINING  | QUIET |
# 8 |     |         | ROOM  |
# 9 +-----+---------+-------+
#
# Grid: 12 cols × 10 rows, 1 m per cell

def default_floor_plan() -> FloorPlan:
    """A 12×10 community-centre floor plan with 4 rooms + corridor."""

    lobby = Room(
        id="lobby",
        name="Lobby",
        cells=_rect_cells(0, 0, 6, 5),
        environment=RoomEnvironment(
            air_temp=23.0, mean_radiant_temp=23.0,
            humidity=50, air_velocity=0.15,
            lux=350, noise_db=55,
        ),
        ceiling_height=3.2,
    )

    activity = Room(
        id="activity",
        name="Activity Room",
        cells=_rect_cells(6, 0, 6, 5),
        environment=RoomEnvironment(
            air_temp=25.0, mean_radiant_temp=25.5,
            humidity=55, air_velocity=0.2,
            lux=400, noise_db=62,
        ),
        ceiling_height=3.2,
    )

    corridor = Room(
        id="corridor",
        name="Corridor",
        cells=_rect_cells(0, 5, 3, 5),
        environment=RoomEnvironment(
            air_temp=22.0, mean_radiant_temp=21.5,
            humidity=45, air_velocity=0.3,
            lux=180, noise_db=48,
        ),
        ceiling_height=2.8,
    )

    dining = Room(
        id="dining",
        name="Dining Hall",
        cells=_rect_cells(3, 5, 6, 5),
        environment=RoomEnvironment(
            air_temp=24.0, mean_radiant_temp=24.0,
            humidity=55, air_velocity=0.1,
            lux=300, noise_db=58,
        ),
        ceiling_height=2.8,
    )

    quiet = Room(
        id="quiet",
        name="Quiet Room",
        cells=_rect_cells(9, 5, 3, 5),
        environment=RoomEnvironment(
            air_temp=23.5, mean_radiant_temp=23.5,
            humidity=50, air_velocity=0.05,
            lux=200, noise_db=32,
        ),
        ceiling_height=2.8,
    )

    # Windows along the outer walls
    windows = [
        Window(p1=(0, 1), p2=(0, 3), room_id="lobby"),       # lobby west wall
        Window(p1=(7, 0), p2=(10, 0), room_id="activity"),    # activity north wall
        Window(p1=(11, 6), p2=(11, 8), room_id="quiet"),      # quiet east wall
        Window(p1=(4, 9), p2=(7, 9), room_id="dining"),       # dining south wall
    ]

    # Doors
    doors = [
        Door(p1=(2, 0), p2=(3, 0), is_exit=True, room_id="lobby"),   # main entrance
        Door(p1=(5, 2), p2=(6, 2), room_id="lobby"),                  # lobby ↔ activity
        Door(p1=(1, 5), p2=(2, 5), room_id="corridor"),               # lobby ↔ corridor
        Door(p1=(3, 7), p2=(3, 8), room_id="dining"),                 # corridor ↔ dining
        Door(p1=(9, 7), p2=(9, 8), room_id="quiet"),                  # dining ↔ quiet
    ]

    return FloorPlan(
        grid_cols=12,
        grid_rows=10,
        cell_size_m=1.0,
        rooms=[lobby, activity, corridor, dining, quiet],
        windows=windows,
        doors=doors,
    )


# ── Agent templates ─────────────────────────────────────────────────

def elder_walker(agent_id: str = "persona_01") -> Agent:
    """78-year-old female using a walker, hearing impaired."""
    return Agent(
        id=agent_id,
        demographics=AgentDemographics(
            age=78, gender="female",
            weight_kg=62.0, height_cm=158.0,
            mobility="uses_walker", hearing="impaired", vision="mild_impairment",
        ),
        preferences=AgentPreferences(
            noise_tolerance_db=55, light_preference="bright", social_density="low_den",
        ),
        clo=0.8,
        activity="walking",
        position=(4, 7),
    )


def young_visitor(agent_id: str = "persona_02") -> Agent:
    """28-year-old male, fully able-bodied."""
    return Agent(
        id=agent_id,
        demographics=AgentDemographics(
            age=28, gender="male",
            weight_kg=78.0, height_cm=178.0,
        ),
        preferences=AgentPreferences(
            noise_tolerance_db=70, light_preference="moderate", social_density="moderate",
        ),
        clo=0.6,
        activity="standing",
        position=(8, 2),
    )


def seated_elder(agent_id: str = "persona_03") -> Agent:
    """85-year-old male, wheelchair user, severe vision impairment."""
    return Agent(
        id=agent_id,
        demographics=AgentDemographics(
            age=85, gender="male",
            weight_kg=72.0, height_cm=168.0,
            mobility="wheelchair", hearing="normal", vision="severe_impairment",
        ),
        preferences=AgentPreferences(
            noise_tolerance_db=50, light_preference="bright", social_density="low_den",
        ),
        clo=1.2,
        activity="seated",
        position=(10, 7),
    )


# ── Quick demo ──────────────────────────────────────────────────────

def demo_simulator():
    """Return a ready-to-use Simulator with the default plan and 3 agents."""
    from .engine import Simulator

    plan = default_floor_plan()
    sim = Simulator(plan)
    sim.add_agent(elder_walker())
    sim.add_agent(young_visitor())
    sim.add_agent(seated_elder())
    return sim
