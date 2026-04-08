#!/usr/bin/env python3
"""Render a publication-oriented thermBAL workflow diagram.

The figure focuses on the actual paper pipeline:
plan authoring -> study definition -> deterministic batch outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


CASE_DIR = Path(__file__).resolve().parent
FIG_DIR = CASE_DIR / "out" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = FIG_DIR / "fig0_framework_pipeline.pdf"
OUT_PNG = FIG_DIR / "fig0_framework_pipeline_hires.png"


INK = "#2d2926"
MUTED = "#8f867c"
RULE = "#d4cabd"
CARD = "#fffdf9"
CARD_ALT = "#f7f1e8"
ACCENT_WARM = "#c8875f"
ACCENT_BLUE = "#6f94b3"
ACCENT_GREEN = "#6f9e85"
ACCENT_GOLD = "#d0a464"


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Avenir Next", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.facecolor": "none",
    }
)


def rounded_box(ax, x, y, w, h, *, fc=CARD, ec=RULE, lw=0.8, rs=0.12, z=1):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.02,rounding_size={rs}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=z,
    )
    ax.add_patch(patch)
    return patch


def label(
    ax,
    x,
    y,
    text,
    *,
    size=8,
    weight="regular",
    color=INK,
    ha="left",
    va="center",
    z=5,
    clip_patch=None,
    linespacing=1.05,
):
    txt = ax.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        zorder=z,
        linespacing=linespacing,
        clip_on=clip_patch is not None,
    )
    if clip_patch is not None:
        txt.set_clip_path(clip_patch)
    return txt


def draw_chip(ax, x, y, w, text, *, fc=CARD_ALT, ec=RULE, color=MUTED):
    patch = rounded_box(ax, x, y, w, 0.34, fc=fc, ec=ec, lw=0.6, rs=0.09, z=3)
    label(ax, x + w / 2, y + 0.17, text, size=4.6, color=color, ha="center", clip_patch=patch)


def draw_pdf_card(ax, x, y, w, h, title, subtitle):
    patch = rounded_box(ax, x, y, w, h, fc=CARD, ec=RULE, lw=0.7, rs=0.12, z=2)
    ax.add_patch(Rectangle((x + 0.12, y + 0.18), 0.34, 0.48, facecolor="#fff8f2", edgecolor=RULE, linewidth=0.6, zorder=3))
    ax.add_patch(Rectangle((x + 0.18, y + 0.58), 0.12, 0.06, facecolor=ACCENT_WARM, edgecolor="none", zorder=4))
    for i in range(3):
        ax.plot([x + 0.18, x + 0.38], [y + 0.48 - 0.09 * i, y + 0.48 - 0.09 * i], color=MUTED, lw=0.55, zorder=4)
    label(ax, x + 0.56, y + 0.56, title, size=5.5, weight="medium", clip_patch=patch)
    label(ax, x + 0.56, y + 0.28, subtitle, size=4.0, color=MUTED, clip_patch=patch)


def draw_grid_card(ax, x, y, w, h):
    patch = rounded_box(ax, x, y, w, h, fc=CARD, ec=RULE, lw=0.7, rs=0.12, z=2)
    label(ax, x + 0.16, y + h - 0.22, "Grid Template", size=5.8, weight="medium", va="top", clip_patch=patch)
    label(ax, x + 0.16, y + h - 0.50, "cells, walls, openings,\nroom IDs", size=3.9, color=MUTED, va="top", clip_patch=patch)
    gx = x + 0.18
    gy = y + 0.18
    cell = 0.16
    cols = 6
    rows = 4
    active = {(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (1, 2), (2, 2), (3, 2), (4, 2), (2, 3), (3, 3)}
    warm = {(0, 1), (1, 1)}
    for row in range(rows):
        for col in range(cols):
            fc = "#f5f1ea"
            if (col, row) in active:
                fc = "#dfe8f0"
            if (col, row) in warm:
                fc = "#efcfb6"
            ax.add_patch(
                Rectangle(
                    (gx + col * cell, gy + row * cell),
                    cell,
                    cell,
                    facecolor=fc,
                    edgecolor=RULE,
                    linewidth=0.45,
                    zorder=3,
                )
            )
    ax.plot([gx + cell, gx + cell * 5], [gy + cell * 2, gy + cell * 2], color=INK, lw=1.2, zorder=4)
    ax.plot([gx + cell * 4, gx + cell * 4], [gy + cell, gy + cell * 3], color=INK, lw=1.0, zorder=4)
    # geometry.json chip removed — implementation detail


def draw_study_card(ax, x, y, w, h):
    patch = rounded_box(ax, x, y, w, h, fc=CARD, ec=RULE, lw=0.7, rs=0.12, z=2)
    label(ax, x + 0.16, y + h - 0.22, "Study Definition", size=5.8, weight="medium", va="top", clip_patch=patch)
    draw_chip(ax, x + 0.18, y + 0.82, 1.10, "room overrides", fc="#eef3f0", color=ACCENT_GREEN)
    draw_chip(ax, x + 1.38, y + 0.82, 0.92, "21 rooms", fc="#f7f1e8", color=ACCENT_GOLD)
    draw_chip(ax, x + 0.18, y + 0.40, 0.92, "6 cohorts", fc="#eef3f0", color=ACCENT_GREEN)
    draw_chip(ax, x + 1.20, y + 0.40, 1.18, "5 Ta x 4 RH", fc="#f7f1e8", color=ACCENT_GOLD)
    # study.json chip removed — implementation detail


def draw_engine_card(ax, x, y, w, h):
    patch = rounded_box(ax, x, y, w, h, fc="#f8f5ef", ec=RULE, lw=0.8, rs=0.12, z=2)
    label(ax, x + 0.16, y + h - 0.22, "Heat-Balance Engine", size=5.6, weight="medium", va="top", clip_patch=patch)
    label(ax, x + 0.16, y + h - 0.50, "deterministic thermal metrics", size=3.8, color=MUTED, va="top", clip_patch=patch)
    draw_chip(ax, x + 0.18, y + 0.94, 0.68, "PMV", fc="#edf2f7", color=ACCENT_BLUE)
    draw_chip(ax, x + 0.92, y + 0.94, 0.68, "PPD", fc="#edf2f7", color=ACCENT_BLUE)
    draw_chip(ax, x + 1.66, y + 0.94, 0.78, "met", fc="#edf2f7", color=ACCENT_BLUE)
    draw_chip(ax, x + 0.18, y + 0.50, 0.90, "convective", fc="#fff5ee", color=ACCENT_WARM)
    draw_chip(ax, x + 1.14, y + 0.50, 0.86, "radiative", fc="#fff5ee", color=ACCENT_WARM)
    draw_chip(ax, x + 0.18, y + 0.12, 1.00, "evaporative", fc="#eef3f0", color=ACCENT_GREEN)
    draw_chip(ax, x + 1.24, y + 0.12, 0.92, "respiratory", fc="#f0ece6", color=MUTED)


def draw_output_card(ax, x, y, w, h, title, accent, kind):
    patch = rounded_box(ax, x, y, w, h, fc=CARD, ec=RULE, lw=0.7, rs=0.12, z=2)
    label(ax, x + 0.12, y + h - 0.18, title, size=5.4, weight="medium", clip_patch=patch)
    if kind == "map":
        cell = 0.12
        colors = ["#dbe7f1", "#bdd1e1", "#f2dfc8", "#d1905d"]
        for row in range(3):
            for col in range(4):
                idx = (row * 2 + col) % len(colors)
                ax.add_patch(Rectangle((x + 0.16 + col * cell, y + 0.18 + row * cell), cell, cell, facecolor=colors[idx], edgecolor=RULE, linewidth=0.35, zorder=3))
    elif kind == "span":
        ax.plot([x + 0.14, x + 0.94], [y + 0.26, y + 0.26], color=RULE, lw=0.5, zorder=3)
        ax.plot([x + 0.14, x + 0.94], [y + 0.48, y + 0.48], color=RULE, lw=0.5, zorder=3)
        pts = [(0.22, 0.30), (0.38, 0.46), (0.56, 0.38), (0.72, 0.54), (0.86, 0.62)]
        for px, py in pts:
            ax.scatter([x + px], [y + py], s=10, color=accent, edgecolors="white", linewidths=0.3, zorder=4)
    elif kind == "bars":
        starts = [0.0, 0.24, 0.52]
        fills = [ACCENT_BLUE, ACCENT_WARM, ACCENT_GREEN]
        for i in range(4):
            left = x + 0.16
            yy = y + 0.18 + i * 0.12
            widths = [0.22, 0.28, 0.20]
            for j, width in enumerate(widths):
                ax.add_patch(Rectangle((left, yy), width, 0.08, facecolor=fills[j], edgecolor="white", linewidth=0.3, zorder=4))
                left += width
    elif kind == "json":
        # Simple record/document icon — stacked lines suggesting structured data
        for i in range(5):
            lw = 0.7 if i == 0 else 0.5
            c = ACCENT_BLUE if i < 2 else (ACCENT_WARM if i < 4 else ACCENT_GREEN)
            w_line = 0.52 - i * 0.04
            ax.plot([x + 0.16, x + 0.16 + w_line],
                    [y + h - 0.38 - i * 0.10, y + h - 0.38 - i * 0.10],
                    color=c, lw=lw, solid_capstyle="round", zorder=3)
        label(ax, x + 0.16, y + 0.14, "per agent", size=3.8, color=MUTED, clip_patch=patch)


def connector(ax, x0, y0, x1, y1, *, rad=0.0, dashed=False):
    patch = FancyArrowPatch(
        (x0, y0),
        (x1, y1),
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=0.8,
        color=MUTED,
        linestyle="--" if dashed else "-",
        connectionstyle=f"arc3,rad={rad}",
        zorder=1,
    )
    ax.add_patch(patch)


def main():
    fig, ax = plt.subplots(figsize=(8.2, 3.7))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # Panel headings
    label(ax, 0.48, 6.18, "1  Spatial Authoring", size=8.4, weight="medium")
    label(ax, 5.00, 6.18, "2  Study Assembly", size=8.4, weight="medium")
    label(ax, 9.38, 6.18, "3  Batch Outputs", size=8.4, weight="medium")

    # Inputs
    draw_pdf_card(ax, 0.48, 4.70, 1.52, 0.94, "PDF Plan", "trace / simplify")
    draw_pdf_card(ax, 2.12, 4.70, 1.52, 0.94, "Rhino", "direct export")
    draw_grid_card(ax, 0.48, 2.38, 3.28, 1.68)

    # Study and engine
    draw_study_card(ax, 5.00, 3.70, 2.94, 1.66)
    draw_engine_card(ax, 5.00, 1.28, 3.46, 1.98)

    # Output cards
    draw_output_card(ax, 9.38, 4.18, 1.72, 1.28, "PMV Maps", ACCENT_WARM, "map")
    draw_output_card(ax, 11.32, 4.18, 1.72, 1.28, "Comfort Span", ACCENT_BLUE, "span")
    draw_output_card(ax, 9.38, 2.44, 1.72, 1.28, "Heat Loss", ACCENT_GREEN, "bars")
    draw_output_card(ax, 11.32, 2.44, 1.72, 1.28, "Snapshots", ACCENT_BLUE, "json")

    draw_chip(ax, 9.42, 1.16, 1.06, "summary.csv", fc="#f6efe7", color=ACCENT_WARM)
    draw_chip(ax, 10.60, 1.16, 1.34, "snapshots.jsonl", fc="#eef3f0", color=ACCENT_GREEN)
    draw_chip(ax, 12.08, 1.16, 0.84, "report", fc="#edf2f7", color=ACCENT_BLUE)

    # Connectors — PDF/Rhino down to Grid Template
    connector(ax, 1.24, 4.70, 1.24, 4.10)
    connector(ax, 2.88, 4.70, 2.88, 4.10)

    # "or" label above the gap between PDF and Rhino
    label(ax, 2.06, 5.72, "or", size=5.5, color=MUTED, ha="center",
          va="center")

    # Grid Template → Study Definition (right edge to left edge)
    connector(ax, 3.76, 3.40, 5.00, 4.53)

    # Study Definition → Engine (bottom to top)
    connector(ax, 6.47, 3.70, 6.73, 3.26)

    # Engine → output cards: top-row arrows curve UP, bottom-row curve DOWN
    # This prevents crossing. Order: near-top, far-top, near-bottom, far-bottom.
    connector(ax, 8.46, 2.65, 9.38, 4.82, rad=0.10)    # → PMV Maps
    connector(ax, 8.46, 2.50, 11.32, 4.82, rad=0.20)   # → Comfort Span
    connector(ax, 8.46, 2.30, 9.38, 3.08, rad=0.02)    # → Heat Loss
    connector(ax, 8.46, 1.95, 11.32, 3.08, rad=-0.06)   # → Agent JSON

    # Bottom note
    rounded_box(ax, 4.84, 0.54, 4.64, 0.38, fc="#fffdf9", ec=RULE, lw=0.55, rs=0.10, z=2)
    label(
        ax,
        7.16,
        0.73,
        "deterministic static snapshots  ·  no occupant-to-room thermal feedback",
        size=4.2,
        color=MUTED,
        ha="center",
    )

    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.05, transparent=True)
    fig.savefig(OUT_PNG, bbox_inches="tight", pad_inches=0.05, dpi=300, transparent=True)
    plt.close(fig)
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
