"""Headless batch-study helpers for thermBAL."""

from __future__ import annotations

import copy
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean
from typing import Any, Dict, Iterable, List, Optional

from .engine import Simulator
from .io import LoadedPlanState
from .models import Agent, AgentDemographics, AgentPreferences
from .presets import elder_walker, seated_elder, young_visitor


PROFILE_BUILDERS = {
    "elder_walker": elder_walker,
    "young_visitor": young_visitor,
    "seated_elder": seated_elder,
}


@dataclass
class AgentSpec:
    id: Optional[str] = None
    profile: Optional[str] = None
    demographics: Dict[str, Any] = field(default_factory=dict)
    preferences: Dict[str, Any] = field(default_factory=dict)
    clo: Optional[float] = None
    activity: Optional[str] = None
    position: Optional[List[int]] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AgentSpec":
        return cls(
            id=data.get("id"),
            profile=data.get("profile"),
            demographics=dict(data.get("demographics") or {}),
            preferences=dict(data.get("preferences") or {}),
            clo=data.get("clo"),
            activity=data.get("activity"),
            position=list(data.get("position")) if data.get("position") is not None else None,
        )


@dataclass(order=True)
class ScheduledEvent:
    t_min: float
    agent_id: Optional[str] = None
    move_to: Optional[List[int]] = None
    activity: Optional[str] = None
    clo: Optional[float] = None
    room_overrides: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    demographics: Dict[str, Any] = field(default_factory=dict)
    preferences: Dict[str, Any] = field(default_factory=dict)
    note: str = ""

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ScheduledEvent":
        return cls(
            t_min=float(data.get("t_min", 0.0)),
            agent_id=data.get("agent_id"),
            move_to=list(data.get("move_to") or data.get("position")) if (data.get("move_to") or data.get("position")) is not None else None,
            activity=data.get("activity"),
            clo=data.get("clo"),
            room_overrides=dict(data.get("room_overrides") or {}),
            demographics=dict(data.get("demographics") or {}),
            preferences=dict(data.get("preferences") or {}),
            note=str(data.get("note") or ""),
        )


@dataclass
class Scenario:
    id: str
    description: str = ""
    use_base_agents: bool = True
    agents: List[AgentSpec] = field(default_factory=list)
    room_overrides: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    schedule: List[ScheduledEvent] = field(default_factory=list)
    snapshot_times: List[float] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Scenario":
        return cls(
            id=str(data.get("id") or "scenario"),
            description=str(data.get("description") or ""),
            use_base_agents=bool(data.get("use_base_agents", True)),
            agents=[AgentSpec.from_dict(item) for item in data.get("agents", []) or []],
            room_overrides=dict(data.get("room_overrides") or {}),
            schedule=sorted(
                [ScheduledEvent.from_dict(item) for item in data.get("schedule", []) or []],
                key=lambda event: event.t_min,
            ),
            snapshot_times=[float(item) for item in data.get("snapshot_times", []) or []],
        )


@dataclass
class StudyDefinition:
    base_plan: Optional[str] = None
    scenarios: List[Scenario] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StudyDefinition":
        return cls(
            base_plan=data.get("base_plan"),
            scenarios=[Scenario.from_dict(item) for item in data.get("scenarios", []) or []],
        )


@dataclass
class SnapshotRecord:
    scenario_id: str
    time_min: float
    agent_id: str
    snapshot: Dict[str, Any]

    def to_jsonable(self) -> Dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "time_min": self.time_min,
            "agent_id": self.agent_id,
            "snapshot": self.snapshot,
        }


@dataclass
class ScenarioResult:
    scenario: Scenario
    records: List[SnapshotRecord] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)


class ManualClock:
    """Simple controllable clock for deterministic batch runs."""

    def __init__(self) -> None:
        self._seconds = 0.0

    def now(self) -> float:
        return self._seconds

    def set_minutes(self, minutes: float) -> None:
        self._seconds = float(minutes) * 60.0


def load_study_definition(path: str | Path) -> StudyDefinition:
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return StudyDefinition.from_dict(payload)


def run_study(base_state: LoadedPlanState, study: StudyDefinition) -> List[ScenarioResult]:
    return [run_scenario(base_state, scenario) for scenario in study.scenarios]


def run_scenario(base_state: LoadedPlanState, scenario: Scenario) -> ScenarioResult:
    plan = copy.deepcopy(base_state.plan)
    agents = _materialize_agents(base_state.agents, scenario)
    _apply_room_overrides(plan, scenario.room_overrides)

    clock = ManualClock()
    sim = Simulator(plan, clock=clock.now)
    for agent in agents:
        sim.add_agent(copy.deepcopy(agent))

    records: List[SnapshotRecord] = []
    events_by_time: Dict[float, List[ScheduledEvent]] = defaultdict(list)
    for event in scenario.schedule:
        events_by_time[event.t_min].append(event)

    timeline = sorted(set([0.0] + scenario.snapshot_times + list(events_by_time.keys())))
    for time_min in timeline:
        clock.set_minutes(time_min)
        for event in events_by_time.get(time_min, []):
            _apply_event(sim, event)
        records.extend(_capture_records(sim, scenario.id, time_min))

    summary = _summarize_records(scenario, records)
    return ScenarioResult(scenario=scenario, records=records, summary=summary)


def _materialize_agents(base_agents: Iterable[Agent], scenario: Scenario) -> List[Agent]:
    agents_by_id: Dict[str, Agent] = {}
    if scenario.use_base_agents:
        for agent in base_agents:
            agents_by_id[agent.id] = copy.deepcopy(agent)

    next_index = len(agents_by_id) + 1
    for spec in scenario.agents:
        agent, next_index = _agent_from_spec(spec, agents_by_id, next_index)
        agents_by_id[agent.id] = agent

    return [agents_by_id[agent_id] for agent_id in sorted(agents_by_id)]


def _agent_from_spec(spec: AgentSpec, existing: Dict[str, Agent], next_index: int) -> tuple[Agent, int]:
    if spec.id and spec.id in existing:
        agent = copy.deepcopy(existing[spec.id])
    elif spec.profile:
        if spec.profile not in PROFILE_BUILDERS:
            raise KeyError("Unknown agent profile: {0}".format(spec.profile))
        agent_id = spec.id or _generated_agent_id(spec.profile, existing, next_index)
        agent = PROFILE_BUILDERS[spec.profile](agent_id)
        next_index += 1
    else:
        agent_id = spec.id or _generated_agent_id("agent", existing, next_index)
        agent = Agent(id=agent_id)
        next_index += 1

    _merge_demographics(agent.demographics, spec.demographics)
    _merge_preferences(agent.preferences, spec.preferences)
    if spec.clo is not None:
        agent.clo = float(spec.clo)
    if spec.activity is not None:
        agent.activity = str(spec.activity)
    if spec.position is not None:
        agent.position = (int(spec.position[0]), int(spec.position[1]))
    return agent, next_index


def _generated_agent_id(prefix: str, existing: Dict[str, Agent], next_index: int) -> str:
    while True:
        candidate = "{0}_{1:02d}".format(prefix, next_index)
        if candidate not in existing:
            return candidate
        next_index += 1


def _merge_demographics(target: AgentDemographics, updates: Dict[str, Any]) -> None:
    for key, value in updates.items():
        if not hasattr(target, key):
            continue
        if key in {"age"}:
            setattr(target, key, int(value))
        elif key in {"weight_kg", "height_cm"}:
            setattr(target, key, float(value))
        else:
            setattr(target, key, str(value))


def _merge_preferences(target: AgentPreferences, updates: Dict[str, Any]) -> None:
    for key, value in updates.items():
        if not hasattr(target, key):
            continue
        if key == "noise_tolerance_db":
            setattr(target, key, float(value))
        else:
            setattr(target, key, str(value))


def _apply_room_overrides(plan, overrides: Dict[str, Dict[str, Any]]) -> None:
    if not overrides:
        return
    room_map = {}
    for room in plan.rooms:
        room_map[room.id] = room
        room_map[room.id.lower()] = room
        room_map[room.name] = room
        room_map[room.name.lower()] = room

    for key, values in overrides.items():
        room = room_map.get(key) or room_map.get(str(key).lower())
        if room is None:
            raise KeyError("Unknown room override target: {0}".format(key))
        for attr, value in values.items():
            if hasattr(room.environment, attr):
                setattr(room.environment, attr, float(value))
            elif attr == "ceilH":
                room.ceiling_height = float(value)
            else:
                raise KeyError("Unknown room attribute override: {0}".format(attr))


def _apply_event(sim: Simulator, event: ScheduledEvent) -> None:
    if event.room_overrides:
        _apply_room_overrides(sim.plan, event.room_overrides)

    if event.agent_id is None:
        return

    agent = sim.get_agent(event.agent_id)
    if event.move_to is not None:
        sim.move_agent(event.agent_id, (int(event.move_to[0]), int(event.move_to[1])))
        agent = sim.get_agent(event.agent_id)
    if event.activity is not None:
        agent.activity = str(event.activity)
    if event.clo is not None:
        agent.clo = float(event.clo)
    if event.demographics:
        _merge_demographics(agent.demographics, event.demographics)
    if event.preferences:
        _merge_preferences(agent.preferences, event.preferences)


def _capture_records(sim: Simulator, scenario_id: str, time_min: float) -> List[SnapshotRecord]:
    records = []
    for agent in sorted(sim.agents, key=lambda item: item.id):
        records.append(
            SnapshotRecord(
                scenario_id=scenario_id,
                time_min=time_min,
                agent_id=agent.id,
                snapshot=sim.snapshot(agent.id),
            )
        )
    return records


def _summarize_records(scenario: Scenario, records: List[SnapshotRecord]) -> Dict[str, Any]:
    pmvs = []
    ppds = []
    errors = 0
    for record in records:
        computed = record.snapshot.get("computed", {})
        pmv = computed.get("pmv")
        ppd = computed.get("ppd")
        if isinstance(pmv, (int, float)):
            pmvs.append(float(pmv))
        if isinstance(ppd, (int, float)):
            ppds.append(float(ppd))
        if "error" in computed or "error" in record.snapshot.get("position", {}) or "error" in record.snapshot.get("environment", {}):
            errors += 1

    return {
        "scenario_id": scenario.id,
        "description": scenario.description,
        "records": len(records),
        "agents": len({record.agent_id for record in records}),
        "times": sorted({record.time_min for record in records}),
        "mean_pmv": round(mean(pmvs), 3) if pmvs else None,
        "min_pmv": round(min(pmvs), 3) if pmvs else None,
        "max_pmv": round(max(pmvs), 3) if pmvs else None,
        "worst_ppd": round(max(ppds), 2) if ppds else None,
        "error_records": errors,
    }
