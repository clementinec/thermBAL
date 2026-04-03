#!/usr/bin/env python3
"""Demo: The Sims-like floor-plan comfort simulator.

Run:  python demo_floorplan.py

Shows how to:
  1.  Load a default floor plan (rooms with fixed environments).
  2.  Add agent personas.
  3.  Move agents around.
  4.  Generate the full JSON snapshot at each position.
"""

import json

from floor_plan import demo_simulator, Agent, AgentDemographics, AgentPreferences


def divider(title: str) -> None:
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}\n")


def main():
    # ── 1. Spin up the simulator with the default plan + 3 agents ───
    sim = demo_simulator()

    divider("ALL AGENTS — initial positions")
    for snap in sim.snapshot_all():
        aid = snap["agent"]["id"]
        room = snap["position"].get("room_name", "?")
        pmv = snap["computed"].get("pmv", "?")
        ppd = snap["computed"].get("ppd", "?")
        print(f"  {aid:15s}  room={room:15s}  PMV={pmv:>6}  PPD={ppd:>5}%")

    # ── 2. Full JSON for persona_01 (elder walker in dining hall) ────
    divider("FULL SNAPSHOT — persona_01 at (4, 7)")
    print(sim.snapshot_json("persona_01"))

    # ── 3. Move persona_01 to the quiet room ────────────────────────
    sim.move_agent("persona_01", (10, 6))
    divider("AFTER MOVE — persona_01 → Quiet Room (10, 6)")
    print(sim.snapshot_json("persona_01"))

    # ── 4. Move persona_01 to the lobby ─────────────────────────────
    sim.move_agent("persona_01", (2, 2))
    divider("AFTER MOVE — persona_01 → Lobby (2, 2)")
    print(sim.snapshot_json("persona_01"))

    # ── 5. Add a brand-new agent on the fly ─────────────────────────
    divider("ADD NEW AGENT — child visitor in Activity Room")
    child = Agent(
        id="persona_04",
        demographics=AgentDemographics(
            age=10, gender="female",
            weight_kg=32.0, height_cm=138.0,
            mobility="normal", hearing="normal", vision="normal",
        ),
        preferences=AgentPreferences(
            noise_tolerance_db=75, light_preference="bright", social_density="high_den",
        ),
        clo=0.5,
        activity="walking",
        position=(8, 3),
    )
    sim.add_agent(child)
    print(sim.snapshot_json("persona_04"))

    # ── 6. Change agent attributes and recompute ────────────────────
    divider("CHANGE ATTRIBUTES — persona_02 sits down, puts on jacket")
    agent2 = sim.get_agent("persona_02")
    agent2.activity = "seated"
    agent2.clo = 1.0
    print(sim.snapshot_json("persona_02"))

    # ── 7. Summary table after all moves ────────────────────────────
    divider("FINAL SUMMARY")
    print(f"  {'Agent':15s}  {'Room':15s}  {'PMV':>6}  {'PPD':>5}  {'Pos'}")
    print(f"  {'-'*15}  {'-'*15}  {'-'*6}  {'-'*5}  {'-'*8}")
    for snap in sim.snapshot_all():
        a = snap["agent"]["id"]
        r = snap["position"].get("room_name", "?")
        pmv = snap["computed"].get("pmv", "?")
        ppd = snap["computed"].get("ppd", "?")
        pos = snap["position"].get("coordinates", "?")
        print(f"  {a:15s}  {r:15s}  {pmv:>6}  {ppd:>5}  {pos}")


if __name__ == "__main__":
    main()
