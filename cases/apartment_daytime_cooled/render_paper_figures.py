#!/usr/bin/env python3
"""Generate ACADIA paper figures for the apartment daytime cooled case.

Produces PDF vector figures in ./out/figures/.

Run:
    python3.10 cases/apartment_daytime_cooled/render_paper_figures.py
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

# ═══════════════════════════════════════════════════════════════════
# Paths
# ═══════════════════════════════════════════════════════════════════

CASE_DIR  = Path(__file__).resolve().parent
ROOT_DIR  = CASE_DIR.parents[1]
TEMPL     = ROOT_DIR / "floor_plan" / "apartment_template.json"
MANIFEST  = CASE_DIR / "study_manifest.json"
SNAPS     = CASE_DIR / "out" / "snapshots.jsonl"
FIG_DIR   = CASE_DIR / "out" / "figures"

if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))


# ═══════════════════════════════════════════════════════════════════
# Visual language
# ═══════════════════════════════════════════════════════════════════

INK   = "#2d2926"
MUTED = "#91887e"
RULE  = "#cdc4b7"
BG    = "none"
PANEL = "#ffffff"

FONT_BODY  = "Avenir Next"
FONT_MONO  = "Andale Mono"
FONT_TITLE = "Avenir Next"

EXP_FILLS = {
    "south":          "#eacdba",
    "east_west":      "#e6d9ba",
    "north_interior": "#c8d5e5",
    "unoccupied":     "#e9e4db",
}

# Three key cohorts (the focused comparison set)
KEY3 = ["higher_clo_sedentary", "young_mixed", "older_shift"]
KEY3_LABELS = ["Higher Clo, Sedentary", "Young Mixed", "Age-Shifted +50 y"]

COHORT_SHORT = {
    "young_mixed":               "Young Mixed",
    "young_male":                "Young Male",
    "older_shift":               "Age +50 y",
    "higher_clo_sedentary":      "Higher Clo",
    "lighter_clothing_mobile":   "Light Clo",
    "higher_bmi_warm_sensitive": "Higher BMI",
}
COHORT_CLR = {
    "young_mixed":               "#4a7fa5",
    "young_male":                "#6ba0be",
    "older_shift":               "#8064a0",
    "higher_clo_sedentary":      "#c47d3a",
    "lighter_clothing_mobile":   "#5a8a5e",
    "higher_bmi_warm_sensitive": "#b55a5a",
}

SOUTH = {"Service Vestibule", "Guest Room", "Family Lounge East",
         "Bath 2", "Bedroom 3", "Laundry", "Bedroom 1", "Bedroom 2"}
EW    = {"Primary Bedroom"}
OCC   = ["Living Room", "Dining Room", "Primary Bedroom", "Kitchen",
         "Breakfast Nook", "Powder Room", "Study", "Family Lounge West",
         "Family Lounge East", "Entry Closet", "Service Vestibule",
         "Guest Room", "Bath 2", "Bedroom 1", "Bedroom 2", "Bedroom 3",
         "West Bedroom Hall", "East Bedroom Hall", "Laundry",
         "Ensuite Bath", "Shared Bath"]

# PMV clamped to ISO 7730 validity range on the cold side
PMV_FLOOR = -3.0

# ── Colormap ───────────────────────────────────────────────────────

PMV_CMAP = LinearSegmentedColormap.from_list("pmv", [
    (0.000, "#2e5a7e"), (0.143, "#4a7fa5"), (0.286, "#86b1cc"),
    (0.400, "#bdd4e5"), (0.500, "#ede5d6"), (0.600, "#e4c9a5"),
    (0.714, "#d09460"), (0.857, "#b86b3a"), (1.000, "#8f4422"),
], N=512)


# ═══════════════════════════════════════════════════════════════════
# Style
# ═══════════════════════════════════════════════════════════════════

def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FONT_BODY, "Avenir", "Helvetica Neue", "DejaVu Sans"],
        "font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 10,
        "axes.titleweight": "medium", "axes.linewidth": 0.5,
        "axes.edgecolor": RULE, "axes.facecolor": PANEL,
        "axes.labelcolor": INK,
        "mathtext.default": "regular",
        "xtick.major.width": 0.4, "ytick.major.width": 0.4,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "figure.facecolor": "white", "figure.dpi": 200,
        "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.15,
        "legend.fontsize": 7.5, "legend.frameon": False,
        "grid.linewidth": 0.3, "grid.color": RULE, "grid.alpha": 0.6,
    })


# ═══════════════════════════════════════════════════════════════════
# Data
# ═══════════════════════════════════════════════════════════════════

def fl(v, d=0.0):
    try:    return d if v in (None, "") else float(v)
    except: return d


def exp_of(rn):
    if rn in SOUTH: return "south"
    if rn in EW:    return "east_west"
    return "north_interior"


def clamp_pmv(v):
    return max(PMV_FLOOR, v)


def load_rooms():
    from floor_plan.io import load_template_json
    st = load_template_json(TEMPL)
    return (
        [{"id": r.id, "name": r.name,
          "cells": [tuple(c) for c in r.cells]} for r in st.plan.rooms],
        st.plan.grid_cols, st.plan.grid_rows,
    )


def load_snapshots():
    out = defaultdict(dict)
    with open(SNAPS, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            snap = rec["snapshot"]
            pos  = snap.get("position", {})
            comp = snap.get("computed", {})
            env  = snap.get("environment", {})
            dyn  = snap.get("dynamic_state", {})
            rn   = pos.get("room_name", "")
            if not rn or "error" in comp:
                continue
            hl = comp.get("heat_loss", {})
            out[rec["scenario_id"]][rn] = {
                "pmv": comp.get("pmv", 0), "ppd": comp.get("ppd", 0),
                "met": comp.get("met", 0), "m_wm2": comp.get("metabolic_w_m2", 0),
                "tcl": comp.get("clothing_temp_c", 0),
                "q_conv": hl.get("convection", 0), "q_rad": hl.get("radiation", 0),
                "e_sk": hl.get("evaporative_skin", 0), "q_res": hl.get("respiratory", 0),
                "ta": env.get("air_temp", 0), "tr": env.get("mean_radiant_temp", 0),
                "rh": env.get("humidity", 0), "v": env.get("air_velocity", 0),
                "clo": dyn.get("clo", 0), "activity": dyn.get("activity", ""),
            }
    return dict(out)


def sid_for(cohort, ta, rh=55):
    return f"{cohort}_t{str(ta).replace('.', 'p')}_rh{rh}"


# ═══════════════════════════════════════════════════════════════════
# Plan renderer
# ═══════════════════════════════════════════════════════════════════

def draw_plan(ax, rooms, gcols, grows, fills, *,
              ext_lw=1.3, int_lw=0.5):
    """Draw apartment plan on *ax*. No room labels."""
    cell_room = {}
    for rm in rooms:
        for c in rm["cells"]:
            cell_room[tuple(c)] = rm["name"]
    active = set(cell_room)

    for rm in rooms:
        fill = fills.get(rm["name"], EXP_FILLS["unoccupied"])
        for col, row in rm["cells"]:
            ax.add_patch(Rectangle((col, row), 1, 1, fc=fill, ec="none", lw=0))

    ext_h, ext_v, int_h, int_v = set(), set(), set(), set()
    for rm in rooms:
        nm = rm["name"]
        for col, row in rm["cells"]:
            for dc, dr, is_h in [(0, -1, True), (0, 1, True),
                                  (-1, 0, False), (1, 0, False)]:
                nb = (col + dc, row + dr)
                if is_h:
                    edge = (col, row if dr == -1 else row + 1)
                else:
                    edge = (col if dc == -1 else col + 1, row)
                if nb not in active:
                    (ext_h if is_h else ext_v).add(edge)
                elif cell_room[nb] != nm:
                    (int_h if is_h else int_v).add(edge)

    def _draw(edges, orient, lw, color):
        by_line = defaultdict(list)
        for a, b in edges:
            (by_line[b] if orient == "h" else by_line[a]).append(
                a if orient == "h" else b)
        for key, vals in by_line.items():
            vals = sorted(set(vals))
            s = vals[0]; e = s + 1
            for v in vals[1:]:
                if v == e: e = v + 1
                else:
                    xy = ([s, e], [key, key]) if orient == "h" else ([key, key], [s, e])
                    ax.plot(*xy, color=color, lw=lw, solid_capstyle="round", zorder=3)
                    s = v; e = v + 1
            xy = ([s, e], [key, key]) if orient == "h" else ([key, key], [s, e])
            ax.plot(*xy, color=color, lw=lw, solid_capstyle="round", zorder=3)

    _draw(int_h, "h", int_lw, RULE)
    _draw(int_v, "v", int_lw, RULE)
    _draw(ext_h, "h", ext_lw, INK)
    _draw(ext_v, "v", ext_lw, INK)

    ax.set_xlim(0, gcols); ax.set_ylim(grows, 0)
    ax.set_aspect("equal"); ax.axis("off")


# ═══════════════════════════════════════════════════════════════════
# Figure 1 — Apartment plan with exposure classification
# ═══════════════════════════════════════════════════════════════════

def figure_1(rooms, gcols, grows):
    fig = plt.figure(figsize=(7.5, 4.6))
    ax = fig.add_axes([0.02, 0.06, 0.64, 0.84])

    occ_set = set(OCC)
    fills = {rm["name"]: EXP_FILLS[exp_of(rm["name"])]
             if rm["name"] in occ_set else EXP_FILLS["unoccupied"]
             for rm in rooms}
    draw_plan(ax, rooms, gcols, grows, fills)
    ax.annotate("S", xy=(gcols / 2, grows + 0.3), ha="center", va="top",
                fontsize=8, fontfamily=FONT_BODY, fontweight="bold", color=MUTED)

    # Legend panel (exposure only, no room directory)
    lg = fig.add_axes([0.70, 0.30, 0.28, 0.45]); lg.axis("off")
    y = 0.95
    lg.text(0, y, "EXPOSURE  CLASSIFICATION", fontsize=7.5, fontweight="bold",
            fontfamily=FONT_BODY, color=MUTED, transform=lg.transAxes, va="top")
    y -= 0.10
    for lab, key in [
        ("South-facing\n$T_r = T_a + 1.5$ K",   "south"),
        ("East / West\n$T_r = T_a + 0.75$ K",    "east_west"),
        ("North / Interior\n$T_r = T_a$",         "north_interior"),
        ("Circulation (unoccupied)",               "unoccupied"),
    ]:
        lg.add_patch(FancyBboxPatch((0, y - 0.04), 0.07, 0.06,
            boxstyle="round,pad=0.003", fc=EXP_FILLS[key], ec=RULE, lw=0.5,
            transform=lg.transAxes, clip_on=False))
        lg.text(0.12, y, lab, fontsize=7.0, fontfamily=FONT_BODY, color=INK,
                transform=lg.transAxes, va="top", linespacing=1.4)
        y -= 0.20

    y -= 0.05
    lg.text(0, y, "21 occupied rooms\n8 south-facing, 1 east/west,\n12 north/interior",
            fontsize=6.5, fontfamily=FONT_BODY, color=MUTED,
            transform=lg.transAxes, va="top", linespacing=1.5)

    fig.text(0.02, 0.96, "Apartment Plan", fontsize=12, fontweight="medium",
             fontfamily=FONT_TITLE, color=INK, va="top")
    fig.text(0.02, 0.935,
             r"24 rooms  $\cdot$  32 $\times$ 24 m grid  $\cdot$  "
             r"1 m cell  $\cdot$  south at bottom",
             fontsize=7, fontfamily=FONT_BODY, color=MUTED, va="top")
    _save(fig, "fig1_apartment_plan.pdf")


# ═══════════════════════════════════════════════════════════════════
# Figure 2 — Same plan, different bodies × different setpoints
# ═══════════════════════════════════════════════════════════════════

def figure_2(rooms, gcols, grows, data):
    temps = [20.0, 25.0, 28.0]
    rh = 55
    occ_set = set(OCC)

    norm = TwoSlopeNorm(vcenter=0.0, vmin=PMV_FLOOR, vmax=1.5)

    fig, axes = plt.subplots(3, 3, figsize=(7.5, 8.0))
    fig.subplots_adjust(left=0.02, right=0.88, bottom=0.03, top=0.86,
                        wspace=0.06, hspace=0.18)

    for ri, ta in enumerate(temps):
        for ci, (co, co_label) in enumerate(zip(KEY3, KEY3_LABELS)):
            ax = axes[ri, ci]
            sid = sid_for(co, ta, rh)
            rd = data.get(sid, {})

            fills = {}
            for rm in rooms:
                n = rm["name"]
                if n in rd and n in occ_set:
                    pmv = clamp_pmv(rd[n]["pmv"])
                    fills[n] = matplotlib.colors.to_hex(PMV_CMAP(norm(pmv)))
                else:
                    fills[n] = EXP_FILLS["unoccupied"]
            draw_plan(ax, rooms, gcols, grows, fills,
                      ext_lw=0.8, int_lw=0.3)

            # Column titles (top row only)
            if ri == 0:
                ax.set_title(co_label, fontsize=9, fontfamily=FONT_BODY,
                             fontweight="medium", color=INK, pad=8)

            # Row labels (left column only)
            if ci == 0:
                ax.text(-0.04, 0.5, f"$T_a$ = {ta:.0f} °C",
                        transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.5,
                        fontfamily=FONT_BODY, fontweight="medium", color=INK)

            # Mean / range annotation
            pmvs = [clamp_pmv(rd[n]["pmv"]) for n in occ_set if n in rd]
            if pmvs:
                ax.text(0.5, -0.03,
                        f"mean {mean(pmvs):+.1f}   [{min(pmvs):+.1f}, {max(pmvs):+.1f}]",
                        transform=ax.transAxes, ha="center", va="top",
                        fontsize=5.5, fontfamily=FONT_MONO, color=MUTED)

    # Shared colorbar
    cax = fig.add_axes([0.91, 0.08, 0.015, 0.74])
    sm = plt.cm.ScalarMappable(cmap=PMV_CMAP, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("PMV", fontsize=8, fontfamily=FONT_BODY, color=INK)
    cb.ax.tick_params(labelsize=6.5, colors=MUTED, width=0.4)
    cb.outline.set_linewidth(0.4); cb.outline.set_edgecolor(RULE)
    for v in (-0.5, 0.5):
        cb.ax.axhline(v, color=INK, lw=0.4, ls="--", alpha=0.5)

    fig.text(0.02, 0.94,
             "Same Spatial Design, Different Conditions and Bodies",
             fontsize=12, fontweight="medium", fontfamily=FONT_TITLE, color=INK)
    fig.text(0.02, 0.895,
             r"RH 55%  $\cdot$  rows = air temperature  $\cdot$  "
             r"columns = occupant cohort  $\cdot$  "
             r"dashed lines on colorbar = $\pm$0.5 comfort zone",
             fontsize=6.5, fontfamily=FONT_BODY, color=MUTED)
    _save(fig, "fig2_design_matrix.pdf")


# ═══════════════════════════════════════════════════════════════════
# Figure 3 — Comfort span across design conditions (strip chart)
# ═══════════════════════════════════════════════════════════════════

def figure_3(data):
    temps = [20.0, 22.0, 23.5, 25.0, 28.0]
    rh = 55

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    fig.subplots_adjust(left=0.10, right=0.74, bottom=0.14, top=0.86)

    # Comfort zone band
    ax.axhspan(-0.5, 0.5, color="#e8e2d4", alpha=0.65, zorder=0)
    ax.axhline(0, color=RULE, lw=0.4, zorder=1)

    width = 0.55
    offsets = np.linspace(-width / 2, width / 2, len(KEY3))
    jscale = 0.08
    rng = np.random.RandomState(42)

    exp_mk = {"south": "^", "east_west": "D", "north_interior": "o"}
    exp_sz = {"south": 16, "east_west": 14, "north_interior": 14}

    for ci, co in enumerate(KEY3):
        color = COHORT_CLR[co]
        for ti, ta in enumerate(temps):
            rd = data.get(sid_for(co, ta, rh), {})
            pmvs = [clamp_pmv(rd[n]["pmv"]) for n in OCC if n in rd]
            exps = [exp_of(n) for n in OCC if n in rd]
            if not pmvs: continue

            x_base = ti + offsets[ci]
            jit = rng.uniform(-jscale, jscale, len(pmvs))

            for pmv, exp, j in zip(pmvs, exps, jit):
                ax.scatter(x_base + j, pmv, s=exp_sz[exp], marker=exp_mk[exp],
                           color=color, alpha=0.80, edgecolors="white",
                           linewidths=0.3, zorder=4)

            m = mean(pmvs)
            ax.plot([x_base - 0.09, x_base + 0.09], [m, m],
                    color=color, lw=2.0, solid_capstyle="round", zorder=5)

    ax.set_xticks(range(len(temps)))
    ax.set_xticklabels([f"{t:.1f}" for t in temps])
    ax.set_xlabel(r"Air temperature  $T_a$  (°C)", fontsize=9,
                  fontfamily=FONT_BODY, color=INK, labelpad=6)
    ax.set_ylabel("Predicted Mean Vote  (PMV)", fontsize=9,
                  fontfamily=FONT_BODY, color=INK, labelpad=6)
    ax.set_xlim(-0.6, len(temps) - 0.4)
    ax.set_ylim(-3.4, 2.2)
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.grid(axis="y", zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    handles = [Line2D([], [], color=COHORT_CLR[c], lw=2.0,
                      label=COHORT_SHORT[c]) for c in KEY3]
    handles += [
        Line2D([], [], marker="o", color=MUTED, lw=0, markersize=4,
               label="North / Interior"),
        Line2D([], [], marker="^", color=MUTED, lw=0, markersize=4,
               label="South-facing"),
        Line2D([], [], marker="D", color=MUTED, lw=0, markersize=3.5,
               label="East / West"),
    ]
    leg = ax.legend(handles=handles, loc="upper left",
                    bbox_to_anchor=(1.03, 1.0), fontsize=7,
                    handlelength=1.5, labelspacing=0.6,
                    title="Cohort / Exposure", title_fontsize=7.5)
    leg.get_title().set_fontweight("medium")
    leg.get_title().set_color(INK)

    fig.text(0.10, 0.94,
             "Comfort Span Across Design Conditions",
             fontsize=12, fontweight="medium", fontfamily=FONT_TITLE, color=INK)
    fig.text(0.10, 0.895,
             r"RH 55%  $\cdot$  each dot = one occupant  $\cdot$  "
             r"horizontal bars = cohort mean  $\cdot$  "
             r"shaded band = ASHRAE comfort zone",
             fontsize=6.5, fontfamily=FONT_BODY, color=MUTED)
    _save(fig, "fig3_comfort_span.pdf")


# ═══════════════════════════════════════════════════════════════════
# Figure 4 — Heat loss decomposition (no room names, exposure chips)
# ═══════════════════════════════════════════════════════════════════

def figure_4(data):
    sid = sid_for("young_mixed", 25.0)
    rd = data.get(sid, {})

    entries = []
    for rn in OCC:
        if rn not in rd: continue
        d = rd[rn]
        entries.append({
            "exp": exp_of(rn),
            "q_conv": d["q_conv"], "q_rad": d["q_rad"],
            "e_sk": d["e_sk"], "q_res": d["q_res"],
            "total": d["q_conv"] + d["q_rad"] + d["e_sk"] + d["q_res"],
            "v": d["v"], "tr": d["tr"],
        })
    entries.sort(key=lambda e: -e["total"])

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    fig.subplots_adjust(left=0.10, right=0.82, bottom=0.12, top=0.86)

    y = np.arange(len(entries))
    bar_h = 0.64

    comps = [
        ("q_conv", "Convective",  "#5b93b8"),
        ("q_rad",  "Radiative",   "#d4956c"),
        ("e_sk",   "Evaporative", "#7bb08a"),
        ("q_res",  "Respiratory", "#b0a89c"),
    ]

    left = np.zeros(len(entries))
    for key, label, color in comps:
        vals = np.array([e[key] for e in entries])
        ax.barh(y, vals, height=bar_h, left=left, color=color,
                edgecolor="white", linewidth=0.3, label=label, zorder=3)
        left += vals

    # Total annotation
    for i, e in enumerate(entries):
        ax.text(e["total"] + 0.4, i, f'{e["total"]:.0f}',
                va="center", ha="left", fontsize=6, fontfamily=FONT_MONO,
                color=MUTED)

    # Exposure chips as Y labels
    ax.set_yticks(y)
    ax.set_yticklabels([""] * len(entries))
    for i, e in enumerate(entries):
        chip_color = EXP_FILLS[e["exp"]]
        ax.add_patch(Rectangle((-3.2, i - bar_h / 2), 2.0, bar_h,
                     fc=chip_color, ec=RULE, lw=0.3, zorder=2, clip_on=False))

    ax.set_xlabel(r"Heat loss from body  (W/m$^2$)", fontsize=9,
                  fontfamily=FONT_BODY, color=INK, labelpad=6)
    ax.set_xlim(0, max(e["total"] for e in entries) + 6)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", zorder=0)

    # Right-side: Tr and v
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(y)
    tags = [f'$T_r$={e["tr"]:.0f}   v={e["v"]:.2f}' for e in entries]
    ax2.set_yticklabels(tags, fontsize=6, fontfamily=FONT_MONO, color=MUTED)
    ax2.tick_params(axis="y", length=0, pad=4)
    for sp in ax2.spines.values():
        sp.set_visible(False)

    leg = ax.legend(loc="lower right", fontsize=7.5, ncol=2,
                    handlelength=1.2, columnspacing=1.0)
    for t in leg.get_texts():
        t.set_color(INK)

    fig.text(0.10, 0.94,
             "Heat Loss Decomposition by Room",
             fontsize=12, fontweight="medium", fontfamily=FONT_TITLE, color=INK)
    fig.text(0.10, 0.895,
             r"Young Mixed  $\cdot$  $T_a$ = 25 °C  $\cdot$  RH 55%  $\cdot$  "
             r"color chips = exposure zone  $\cdot$  "
             r"right labels = $T_r$ and air velocity",
             fontsize=6.5, fontfamily=FONT_BODY, color=MUTED)
    _save(fig, "fig4_heat_loss_decomposition.pdf")


# ═══════════════════════════════════════════════════════════════════
# Figure 5 — Convective–radiative scatter (no room annotations)
# ═══════════════════════════════════════════════════════════════════

def figure_5(data):
    rh = 55
    ta = 25.0

    fig, ax = plt.subplots(figsize=(5.8, 5.3))
    fig.subplots_adjust(left=0.14, right=0.72, bottom=0.12, top=0.88)

    exp_mk = {"north_interior": "o", "south": "^", "east_west": "D"}
    exp_sz = {"north_interior": 28, "south": 32, "east_west": 28}

    # Diagonal
    diag = np.linspace(10, 35, 50)
    ax.plot(diag, diag, color=RULE, lw=0.6, ls="--", zorder=1)
    ax.text(31, 30, r"$Q_c = Q_r$", fontsize=6.5, fontfamily=FONT_MONO,
            color="#b5ada2", rotation=40, va="bottom", ha="right")

    for co in KEY3:
        sid = sid_for(co, ta, rh)
        rd = data.get(sid, {})
        color = COHORT_CLR[co]
        for rn in OCC:
            if rn not in rd: continue
            d = rd[rn]
            exp = exp_of(rn)
            ax.scatter(d["q_conv"], d["q_rad"],
                       s=exp_sz[exp], marker=exp_mk[exp],
                       color=color, alpha=0.82,
                       edgecolors="white", linewidths=0.4, zorder=4)

    ax.set_xlabel(r"Convective heat loss  $Q_c$  (W/m$^2$)", fontsize=9,
                  fontfamily=FONT_BODY, color=INK, labelpad=6)
    ax.set_ylabel(r"Radiative heat loss  $Q_r$  (W/m$^2$)", fontsize=9,
                  fontfamily=FONT_BODY, color=INK, labelpad=6)
    ax.set_xlim(15, 33); ax.set_ylim(13, 30)
    ax.set_aspect("equal")
    ax.grid(zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Region labels
    ax.text(0.04, 0.96, "convection-\ndominant",
            transform=ax.transAxes, fontsize=6, fontfamily=FONT_BODY,
            color="#b5ada2", va="top", linespacing=1.3)
    ax.text(0.96, 0.04, "radiation-\ndominant",
            transform=ax.transAxes, fontsize=6, fontfamily=FONT_BODY,
            color="#b5ada2", ha="right", va="bottom", linespacing=1.3)

    # Legend
    coh_h = [Line2D([], [], marker="o", color=COHORT_CLR[c], lw=0,
                    markersize=5, label=COHORT_SHORT[c]) for c in KEY3]
    exp_h = [
        Line2D([], [], marker="o", color=MUTED, lw=0, markersize=4,
               label="North / Interior"),
        Line2D([], [], marker="^", color=MUTED, lw=0, markersize=4.5,
               label="South-facing"),
        Line2D([], [], marker="D", color=MUTED, lw=0, markersize=3.5,
               label="East / West"),
    ]
    leg = ax.legend(handles=coh_h + exp_h, loc="upper left",
                    bbox_to_anchor=(1.03, 1.0), fontsize=7,
                    handlelength=1.2, labelspacing=0.6,
                    title="Cohort / Exposure", title_fontsize=7.5)
    leg.get_title().set_fontweight("medium")
    leg.get_title().set_color(INK)

    fig.text(0.14, 0.95,
             "Convective\u2013Radiative Heat Loss",
             fontsize=12, fontweight="medium", fontfamily=FONT_TITLE, color=INK)
    fig.text(0.14, 0.905,
             r"$T_a$ = 25 °C  $\cdot$  RH 55%  $\cdot$  "
             r"each point = one room $\times$ one cohort  $\cdot$  "
             r"marker shape = solar exposure",
             fontsize=6.5, fontfamily=FONT_BODY, color=MUTED)
    _save(fig, "fig5_conv_rad_scatter.pdf")


# ═══════════════════════════════════════════════════════════════════

def _save(fig, name):
    out = FIG_DIR / name
    fig.savefig(out, facecolor="none", transparent=True)
    plt.close(fig)
    print(f"  \u2713 {name}")


def main():
    apply_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data\u2026")
    rooms, gcols, grows = load_rooms()
    data = load_snapshots()
    print(f"  {len(rooms)} rooms  \u00b7  {gcols}\u00d7{grows} grid"
          f"  \u00b7  {len(data)} scenarios")

    print("\nRendering figures\u2026")
    figure_1(rooms, gcols, grows)
    figure_2(rooms, gcols, grows, data)
    figure_3(data)
    figure_4(data)
    figure_5(data)
    print(f"\nAll figures \u2192 {FIG_DIR}/")


if __name__ == "__main__":
    main()
