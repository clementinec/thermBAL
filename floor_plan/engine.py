"""Floor-plan experience engine.

The Simulator is the main entry point:
  1.  Load a FloorPlan (rooms, windows, doors).
  2.  Add Agent(s).
  3.  Move agents around.
  4.  Call  snapshot(agent_id)  to get the full JSON for that agent at their
      current position — environment, spatial, computed PMV/PPD, and
      dynamic state — ready for downstream LLM consumption.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from comfort_engine import (
    ComfortState,
    Occupant,
    compute_pmv,
    du_bois_bsa,
    estimate_metabolic_profile,
)

from .models import Agent, Coord, FloorPlan
from .spatial import compute_spatial


# Activity → approximate met multiplier over RMR
_ACTIVITY_MET_MULT: Dict[str, float] = {
    "seated": 1.0,
    "standing": 1.2,
    "walking": 1.6,
}

# Activity → posture (for PMV model)
_ACTIVITY_POSTURE: Dict[str, str] = {
    "seated": "seated",
    "standing": "standing",
    "walking": "standing",
}


@dataclass
class _AgentState:
    """Mutable runtime state tracked per agent."""
    placed_at: float = 0.0           # epoch when agent entered current cell
    last_snapshot: Optional[Dict] = None


class Simulator:
    """The Sims-like experience engine.

    Usage
    -----
    >>> sim = Simulator(plan)
    >>> sim.add_agent(agent)
    >>> sim.move_agent("persona_01", (4, 7))
    >>> print(json.dumps(sim.snapshot("persona_01"), indent=2))
    """

    def __init__(self, plan: FloorPlan) -> None:
        self.plan = plan
        self.plan.build_index()
        self._agents: Dict[str, Agent] = {}
        self._state: Dict[str, _AgentState] = {}

    # ── Agent management ────────────────────────────────────────────

    def add_agent(self, agent: Agent) -> None:
        self._agents[agent.id] = agent
        self._state[agent.id] = _AgentState(placed_at=time.time())

    def remove_agent(self, agent_id: str) -> None:
        self._agents.pop(agent_id, None)
        self._state.pop(agent_id, None)

    def move_agent(self, agent_id: str, new_pos: Coord) -> None:
        agent = self._agents[agent_id]
        if agent.position != new_pos:
            agent.position = new_pos
            self._state[agent_id].placed_at = time.time()

    def get_agent(self, agent_id: str) -> Agent:
        return self._agents[agent_id]

    @property
    def agents(self) -> List[Agent]:
        return list(self._agents.values())

    # ── Comfort computation ─────────────────────────────────────────

    def _compute_comfort(self, agent: Agent) -> Dict[str, Any]:
        """Run PMV/PPD for the agent at their current position."""
        room = self.plan.room_at(agent.position)
        if room is None:
            return {"error": "agent not in any room"}

        env = room.environment
        d = agent.demographics

        # Metabolic rate from demographics + activity
        activity_mult = _ACTIVITY_MET_MULT.get(agent.activity, 1.2)
        profile = estimate_metabolic_profile(
            sex=d.gender,
            weight_kg=d.weight_kg,
            height_cm=d.height_cm,
            age=d.age,
            activity_multiplier=activity_mult,
        )

        state = ComfortState(
            Ta=env.air_temp,
            Tr=env.mean_radiant_temp,
            RH=env.humidity,
            v=env.air_velocity,
        )
        occupant = Occupant(
            m_wm2=profile["M_active_Wm2"],
            clo=agent.clo,
            posture=_ACTIVITY_POSTURE.get(agent.activity, "seated"),
            body_area_m2=profile["BSA_m2"],
        )

        result = compute_pmv(state, occupant)

        return {
            "pmv": round(result.PMV, 2),
            "ppd": round(result.PPD, 1),
            "metabolic_w_m2": round(profile["M_active_Wm2"], 1),
            "met": round(profile["met_active"], 2),
            "clothing_temp_c": round(result.Tcl, 1),
            "heat_loss": {
                "convection": round(result.q_conv, 1),
                "radiation": round(result.q_rad, 1),
                "evaporative_skin": round(result.e_sk, 1),
                "respiratory": round(result.q_res, 1),
            },
        }

    # ── Full snapshot ───────────────────────────────────────────────

    def snapshot(self, agent_id: str) -> Dict[str, Any]:
        """Produce the complete JSON record for an agent.

        This is the primary output — structured to match the schema the
        user specified, ready for LLM narrative generation.
        """
        agent = self._agents[agent_id]
        st = self._state[agent_id]
        room = self.plan.room_at(agent.position)
        all_agents = self.agents

        # ── agent block
        agent_block = {
            "id": agent.id,
            "age": agent.demographics.age,
            "gender": agent.demographics.gender,
            "mobility": agent.demographics.mobility,
            "hearing": agent.demographics.hearing,
            "vision": agent.demographics.vision,
        }

        # ── position / spatial block
        spatial = compute_spatial(agent, self.plan, all_agents)

        # ── environment block (from the room the agent is in)
        if room:
            env = room.environment
            env_block = {
                "lux": env.lux,
                "noise_db": env.noise_db,
                "air_temp": env.air_temp,
                "mean_radiant_temp": env.mean_radiant_temp,
                "humidity": env.humidity,
                "air_velocity": env.air_velocity,
            }
        else:
            env_block = {"error": "outside all rooms"}

        # ── computed block (PMV / PPD + heat breakdown)
        computed = self._compute_comfort(agent)

        # ── preferences block
        prefs = agent.preferences
        prefs_block = {
            "noise_tolerance_db": prefs.noise_tolerance_db,
            "light_preference": prefs.light_preference,
            "social_density": prefs.social_density,
        }

        # ── dynamic state
        duration_mins = round((time.time() - st.placed_at) / 60.0, 1)
        dynamic_block = {
            "duration_mins": duration_mins,
            "activity": agent.activity,
            "clo": agent.clo,
        }

        snapshot = {
            "agent": agent_block,
            "preferences": prefs_block,
            "position": spatial,
            "environment": env_block,
            "computed": computed,
            "dynamic_state": dynamic_block,
        }

        st.last_snapshot = snapshot
        return snapshot

    def snapshot_all(self) -> List[Dict[str, Any]]:
        """Produce snapshots for every agent on the plan."""
        return [self.snapshot(aid) for aid in self._agents]

    def snapshot_json(self, agent_id: str, indent: int = 2) -> str:
        """Convenience: snapshot as a formatted JSON string."""
        return json.dumps(self.snapshot(agent_id), indent=indent)
