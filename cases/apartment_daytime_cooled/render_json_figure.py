#!/usr/bin/env python3
"""Render a publication-quality JSON agent snapshot figure for ACADIA 2026.

Produces a monospace code-block-style image with syntax-highlighting colours.

Run:
    python3.10 cases/apartment_daytime_cooled/render_json_figure.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ═══════════════════════════════════════════════════════════════════
# Paths
# ═══════════════════════════════════════════════════════════════════

FIG_DIR = Path(__file__).resolve().parent / "out" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = FIG_DIR / "fig_json_snippet.pdf"
OUT_PNG = FIG_DIR / "fig_json_snippet_hires.png"

# ═══════════════════════════════════════════════════════════════════
# Visual language (matches render_paper_figures.py)
# ═══════════════════════════════════════════════════════════════════

BG_FIGURE = "#faf8f4"
BG_CODE   = "#f9f6f1"
BORDER    = "#cdc4b7"

CLR_KEY   = "#5b7a8a"   # muted blue – keys
CLR_STR   = "#8b6b4a"   # warm brown – string values
CLR_NUM   = "#2d2926"   # dark – numbers
CLR_PUNCT = "#91887e"   # muted gray – braces, brackets, commas, colons

FONT_MONO  = "Andale Mono"
FONT_TITLE = "Avenir Next"

# ═══════════════════════════════════════════════════════════════════
# JSON data
# ═══════════════════════════════════════════════════════════════════

SNAPSHOT = {
    "agent": {
        "age": 29, "gender": "female",
        "mobility": "normal"
    },
    "position": {
        "coordinates": [13, 18],
        "dist_to_wall": 2.0,
        "dist_to_window": 3.5,
        "enclosure_ratio": 0.25,
        "visible_agents": 0
    },
    "environment": {
        "air_temp": 25.0,
        "mean_radiant_temp": 26.5,
        "humidity": 55.0,
        "air_velocity": 0.08,
        "lux": 300, "noise_db": 40
    },
    "computed": {
        "pmv": -0.96, "ppd": 24.8,
        "metabolic_w_m2": 47.1,
        "clothing_temp_c": 30.7,
        "heat_loss": {
            "convection": 23.1,
            "radiation": 20.2,
            "evaporative_skin": 11.2,
            "respiratory": 3.9
        }
    },
    "dynamic_state": {
        "activity": "seated",
        "clo": 0.55,
        "duration_mins": 0.0
    }
}

# ═══════════════════════════════════════════════════════════════════
# Tokeniser – split JSON text into (token, colour) pairs per line
# ═══════════════════════════════════════════════════════════════════

def tokenise_json(obj: dict, indent: int = 2) -> list[list[tuple[str, str]]]:
    """Return a list of lines, each line a list of (text, colour) spans."""
    raw = json.dumps(obj, indent=indent)
    lines_out: list[list[tuple[str, str]]] = []

    for line in raw.split("\n"):
        spans: list[tuple[str, str]] = []
        i = 0
        n = len(line)
        while i < n:
            ch = line[i]
            # leading whitespace
            if ch == " " and (not spans or spans[-1][1] == CLR_PUNCT):
                ws = ""
                while i < n and line[i] == " ":
                    ws += " "
                    i += 1
                spans.append((ws, CLR_PUNCT))
                continue
            # punctuation
            if ch in "{}[],:":
                spans.append((ch, CLR_PUNCT))
                i += 1
                # space after colon / comma
                if i < n and line[i] == " ":
                    spans.append((" ", CLR_PUNCT))
                    i += 1
                continue
            # quoted string
            if ch == '"':
                j = i + 1
                while j < n and line[j] != '"':
                    if line[j] == '\\':
                        j += 1
                    j += 1
                token = line[i:j + 1]
                # decide colour: key if followed by ':', else string value
                rest = line[j + 1:].lstrip()
                if rest.startswith(":"):
                    spans.append((token, CLR_KEY))
                else:
                    spans.append((token, CLR_STR))
                i = j + 1
                continue
            # number / literal (true, false, null)
            m = re.match(r'-?[\d.eE+\-]+|true|false|null', line[i:])
            if m:
                spans.append((m.group(), CLR_NUM))
                i += m.end()
                continue
            # fallback
            spans.append((ch, CLR_PUNCT))
            i += 1

        lines_out.append(spans)
    return lines_out


# ═══════════════════════════════════════════════════════════════════
# Render
# ═══════════════════════════════════════════════════════════════════

def render():
    lines = tokenise_json(SNAPSHOT, indent=2)
    n_lines = len(lines)

    fig_w, fig_h = 4.0, 5.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_alpha(0.0)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Title
    ax.text(
        0.5, 0.975,
        "Agent Snapshot (one occupant record)",
        ha="center", va="top",
        fontfamily=FONT_TITLE, fontsize=9.5, fontweight="medium",
        color="#4a4440", style="italic",
        transform=ax.transAxes,
    )

    # Code block bounding box (in axes coords)
    pad_x = 0.05
    pad_y_top = 0.06
    pad_y_bot = 0.02
    box_x = pad_x
    box_y = pad_y_bot
    box_w = 1.0 - 2 * pad_x
    box_h = 1.0 - pad_y_top - pad_y_bot - 0.03  # leave room for title

    # Draw background rectangle
    rect = FancyBboxPatch(
        (box_x, box_y), box_w, box_h,
        boxstyle="round,pad=0.012",
        facecolor=BG_CODE,
        edgecolor=BORDER,
        linewidth=0.8,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(rect)

    # Compute text positions
    font_size = 8.0
    # inner margins within the code box
    mx = 0.040
    my = 0.025
    text_x = box_x + mx
    text_top = box_y + box_h - my
    text_bot = box_y + my
    line_height = (text_top - text_bot) / max(n_lines, 1)

    # Render each line with coloured spans.
    # For monospace fonts every character has the same advance width.
    # We place each span at the correct *character* offset from the left.
    #
    # Measure one character's width using the renderer so we don't guess.
    fig.canvas.draw()  # initialise renderer
    renderer = fig.canvas.get_renderer()
    _probe = ax.text(
        0, 0, "M",
        fontfamily=FONT_MONO, fontsize=font_size,
        transform=ax.transAxes,
    )
    _bb = _probe.get_window_extent(renderer=renderer)
    # Convert pixel width to axes-fraction width
    ax_bbox = ax.get_window_extent(renderer=renderer)
    char_w = _bb.width / ax_bbox.width
    _probe.remove()

    for li, spans in enumerate(lines):
        y = text_top - li * line_height
        char_offset = 0  # running character count from start of line
        for text_fragment, colour in spans:
            x_pos = text_x + char_offset * char_w
            ax.text(
                x_pos, y,
                text_fragment,
                ha="left", va="top",
                fontfamily=FONT_MONO,
                fontsize=font_size,
                color=colour,
                transform=ax.transAxes,
                clip_on=True,
            )
            char_offset += len(text_fragment)

    # ─── Save ──────────────────────────────────────────────────────
    fig.savefig(
        str(OUT_PDF),
        format="pdf",
        bbox_inches="tight",
        pad_inches=0.05,
        transparent=True,
    )
    fig.savefig(
        str(OUT_PNG),
        format="png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.05,
        transparent=True,
    )
    plt.close(fig)
    print(f"Saved  {OUT_PDF}")
    print(f"Saved  {OUT_PNG}")


if __name__ == "__main__":
    render()
