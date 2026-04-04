"""thermBAL Rhino exporter.

Export an orthogonal Rhino floor plan into thermBAL geometry JSON.

Designed to be easy to copy into the Rhino Python editor:
  1. Put geometry on the layers listed in LAYER_NAMES below.
  2. Run this script in Rhino.
  3. Enter the grid cell size in model units.
  4. Save the JSON file.

The output can be loaded directly in simulator.html via
"Load Geometry JSON".

Notes
-----
- Rhino 7/8 compatible RhinoPython style: avoids modern Python-only syntax.
- Best for rectilinear floor plans aligned to World XY.
- Exterior walls are auto-generated from the footprint.
- Interior wall/opening curves should already be snapped to the intended grid.
"""

from __future__ import print_function

import json
import math
import os
from datetime import datetime

import Rhino
import rhinoscriptsyntax as rs
import scriptcontext as sc

try:
    basestring
except NameError:
    basestring = str


LAYER_NAMES = {
    "footprint": "tb_footprint",
    "voids": "tb_voids",
    "walls": "tb_walls",
    "windows": "tb_windows",
    "doors": "tb_doors",
    "exits": "tb_exits",
    "rooms": "tb_rooms",
    "origin": "tb_origin",
}

DEFAULT_ENV = {
    "air_temp": 24.0,
    "mean_radiant_temp": 24.0,
    "humidity": 50.0,
    "air_velocity": 0.1,
    "lux": 300.0,
    "noise_db": 40.0,
}

DEFAULT_CEIL_H = 2.8
STICKY_CELL_SIZE = "thermBAL.rhino_export.cell_size"
STICKY_TEMPLATE_NAME = "thermBAL.rhino_export.template_name"

SNAP_TOL_FRAC = 0.20
DIR_TOL_FRAC = 0.02


def log(message):
    print(message)
    try:
        Rhino.RhinoApp.WriteLine(str(message))
    except Exception:
        pass


def warn_once(warnings, message):
    if message not in warnings:
        warnings.append(message)


def object_ids_on_layer(layer_name):
    if not rs.IsLayer(layer_name):
        return []
    return list(rs.ObjectsByLayer(layer_name, select=False) or [])


def curve_infos_on_layer(layer_name, require_closed, warnings):
    infos = []
    for obj_id in object_ids_on_layer(layer_name):
        if not rs.IsCurve(obj_id):
            warn_once(
                warnings,
                "Skipped non-curve object on layer '{0}'.".format(layer_name),
            )
            continue

        curve = rs.coercecurve(obj_id)
        if curve is None:
            warn_once(
                warnings,
                "Could not read a curve on layer '{0}'.".format(layer_name),
            )
            continue

        if require_closed and not curve.IsClosed:
            warn_once(
                warnings,
                "Skipped open curve on layer '{0}'. Closed planar curves are required.".format(layer_name),
            )
            continue

        ok, plane = curve.TryGetPlane(sc.doc.ModelAbsoluteTolerance)
        if not ok:
            warn_once(
                warnings,
                "Skipped non-planar curve on layer '{0}'.".format(layer_name),
            )
            continue

        infos.append(
            {
                "id": obj_id,
                "curve": curve,
                "plane": plane,
                "name": rs.ObjectName(obj_id) or "",
                "user_text": read_user_text(obj_id),
            }
        )
    return infos


def read_user_text(obj_id):
    data = {}
    keys = rs.GetUserText(obj_id)
    if not keys:
        return data
    if isinstance(keys, basestring):
        keys = [keys]
    for key in keys:
        data[str(key)] = rs.GetUserText(obj_id, key)
    return data


def point_in_curve_info(x, y, curve_info):
    plane = curve_info["plane"]
    point = Rhino.Geometry.Point3d(x, y, plane.Origin.Z)
    containment = curve_info["curve"].Contains(
        point,
        plane,
        sc.doc.ModelAbsoluteTolerance,
    )
    return containment != Rhino.Geometry.PointContainment.Outside


def point_in_any_curve(x, y, curve_infos):
    for info in curve_infos:
        if point_in_curve_info(x, y, info):
            return True
    return False


def bbox_min_max(obj_ids):
    bbox = rs.BoundingBox(obj_ids)
    if not bbox:
        return None
    min_x = min(pt.X for pt in bbox)
    min_y = min(pt.Y for pt in bbox)
    min_z = min(pt.Z for pt in bbox)
    max_x = max(pt.X for pt in bbox)
    max_y = max(pt.Y for pt in bbox)
    max_z = max(pt.Z for pt in bbox)
    return (
        Rhino.Geometry.Point3d(min_x, min_y, min_z),
        Rhino.Geometry.Point3d(max_x, max_y, max_z),
    )


def get_origin_point(all_obj_ids):
    origin_layer = LAYER_NAMES["origin"]
    for obj_id in object_ids_on_layer(origin_layer):
        if rs.IsPoint(obj_id):
            return rs.PointCoordinates(obj_id)
        if rs.IsTextDot(obj_id):
            return rs.TextDotPoint(obj_id)

    bounds = bbox_min_max(all_obj_ids)
    if not bounds:
        return None
    return bounds[0]


def model_unit_scale_to_meters():
    try:
        return Rhino.RhinoMath.UnitScale(
            sc.doc.ModelUnitSystem,
            Rhino.UnitSystem.Meters,
        )
    except Exception:
        return 1.0


def default_template_name():
    doc_name = sc.doc.Name or "rhino-model"
    return os.path.splitext(os.path.basename(doc_name))[0]


def prompt_cell_size():
    default_size = sc.sticky.get(STICKY_CELL_SIZE, 1.0)
    value = rs.GetReal(
        "Grid cell size in model units",
        default_size,
        0.000001,
    )
    if value is not None:
        sc.sticky[STICKY_CELL_SIZE] = value
    return value


def prompt_template_name():
    default_name = sc.sticky.get(STICKY_TEMPLATE_NAME, default_template_name())
    value = rs.GetString("Template name", default_name)
    if value:
        sc.sticky[STICKY_TEMPLATE_NAME] = value
    return value


def prompt_output_path(default_name):
    filename = "{0}-geometry.json".format(default_name)
    return rs.SaveFileName(
        "Save thermBAL geometry JSON",
        "JSON Files (*.json)|*.json||",
        None,
        filename,
    )


def compute_grid(origin, max_point, cell_size):
    width = max_point.X - origin.X
    height = max_point.Y - origin.Y
    if width < 0 or height < 0:
        return None
    cols = int(math.ceil(width / float(cell_size) - 1e-9))
    rows = int(math.ceil(height / float(cell_size) - 1e-9))
    return max(1, cols), max(1, rows)


def sample_active_cells(cols, rows, origin, cell_size, footprint_infos, void_infos):
    active = set()
    for row in range(rows):
        y = origin.Y + (row + 0.5) * cell_size
        for col in range(cols):
            x = origin.X + (col + 0.5) * cell_size
            inside = point_in_any_curve(x, y, footprint_infos)
            if not inside:
                continue
            if void_infos and point_in_any_curve(x, y, void_infos):
                continue
            active.add((col, row))
    return active


def empty_wall_matrices(cols, rows):
    return (
        [[False for _ in range(cols)] for _ in range(rows + 1)],
        [[False for _ in range(cols + 1)] for _ in range(rows)],
    )


def build_exterior_walls(active_cells, h_walls, v_walls):
    active = set(active_cells)
    for col, row in active:
        if (col, row - 1) not in active:
            h_walls[row][col] = True
        if (col, row + 1) not in active:
            h_walls[row + 1][col] = True
        if (col - 1, row) not in active:
            v_walls[row][col] = True
        if (col + 1, row) not in active:
            v_walls[row][col + 1] = True


def edge_has_wall(edge, h_walls, v_walls):
    orient, row, col = edge
    if orient == "h":
        return bool(h_walls[row][col])
    return bool(v_walls[row][col])


def edge_touches_active(edge, active_cells):
    orient, row, col = edge
    if orient == "h":
        return ((col, row - 1) in active_cells) or ((col, row) in active_cells)
    return ((col - 1, row) in active_cells) or ((col, row) in active_cells)


def edge_is_exterior(edge, active_cells):
    orient, row, col = edge
    if orient == "h":
        count = int((col, row - 1) in active_cells) + int((col, row) in active_cells)
    else:
        count = int((col - 1, row) in active_cells) + int((col, row) in active_cells)
    return count == 1


def segment_curves(obj_id):
    curve = rs.coercecurve(obj_id)
    if curve is None:
        return []
    try:
        segments = curve.DuplicateSegments()
    except Exception:
        segments = None
    if segments:
        return list(segments)
    return [curve]


def snap_index(value, origin_value, cell_size):
    raw = (value - origin_value) / float(cell_size)
    snapped = int(round(raw))
    if abs(raw - snapped) > SNAP_TOL_FRAC:
        return None
    return snapped


def segment_to_edges(pt_a, pt_b, origin, cell_size, cols, rows, label, warnings):
    dx = pt_b.X - pt_a.X
    dy = pt_b.Y - pt_a.Y
    abs_dx = abs(dx)
    abs_dy = abs(dy)
    dir_tol = max(sc.doc.ModelAbsoluteTolerance * 2.0, cell_size * DIR_TOL_FRAC)

    if abs_dx <= dir_tol and abs_dy <= dir_tol:
        return []

    if abs_dx > dir_tol and abs_dy > dir_tol:
        warn_once(
            warnings,
            "Skipped non-orthogonal {0} segment. Only horizontal/vertical segments export cleanly.".format(label),
        )
        return []

    edges = []

    if abs_dy <= dir_tol:
        row = snap_index((pt_a.Y + pt_b.Y) * 0.5, origin.Y, cell_size)
        col_0 = snap_index(min(pt_a.X, pt_b.X), origin.X, cell_size)
        col_1 = snap_index(max(pt_a.X, pt_b.X), origin.X, cell_size)
        if row is None or col_0 is None or col_1 is None or col_1 <= col_0:
            warn_once(
                warnings,
                "Skipped {0} segment that does not land cleanly on the grid.".format(label),
            )
            return []
        for col in range(col_0, col_1):
            if 0 <= row <= rows and 0 <= col < cols:
                edges.append(("h", row, col))
        return edges

    col = snap_index((pt_a.X + pt_b.X) * 0.5, origin.X, cell_size)
    row_0 = snap_index(min(pt_a.Y, pt_b.Y), origin.Y, cell_size)
    row_1 = snap_index(max(pt_a.Y, pt_b.Y), origin.Y, cell_size)
    if col is None or row_0 is None or row_1 is None or row_1 <= row_0:
        warn_once(
            warnings,
            "Skipped {0} segment that does not land cleanly on the grid.".format(label),
        )
        return []
    for row in range(row_0, row_1):
        if 0 <= row < rows and 0 <= col <= cols:
            edges.append(("v", row, col))
    return edges


def add_wall_segments(wall_ids, h_walls, v_walls, origin, cell_size, cols, rows, warnings):
    for obj_id in wall_ids:
        for segment in segment_curves(obj_id):
            edges = segment_to_edges(
                segment.PointAtStart,
                segment.PointAtEnd,
                origin,
                cell_size,
                cols,
                rows,
                "wall",
                warnings,
            )
            for orient, row, col in edges:
                if orient == "h":
                    h_walls[row][col] = True
                else:
                    v_walls[row][col] = True


def feature_edges_from_layer(layer_name, feature_type, origin, cell_size, cols, rows, warnings):
    edge_set = set()
    for obj_id in object_ids_on_layer(layer_name):
        if not rs.IsCurve(obj_id):
            warn_once(
                warnings,
                "Skipped non-curve object on layer '{0}'.".format(layer_name),
            )
            continue
        for segment in segment_curves(obj_id):
            edges = segment_to_edges(
                segment.PointAtStart,
                segment.PointAtEnd,
                origin,
                cell_size,
                cols,
                rows,
                feature_type,
                warnings,
            )
            for orient, row, col in edges:
                edge_set.add((feature_type, orient, row, col))
    return edge_set


def room_name(room_info, index):
    if room_info["name"]:
        return room_info["name"]
    value = room_info["user_text"].get("name")
    if value:
        return value
    return "Room {0}".format(chr(65 + (index % 26)))


def read_float(text, default_value):
    try:
        return float(text)
    except Exception:
        return default_value


def room_payload(room_info, index, cells):
    user_text = room_info["user_text"]
    env = dict(DEFAULT_ENV)
    for key in DEFAULT_ENV.keys():
        if key in user_text:
            env[key] = read_float(user_text[key], DEFAULT_ENV[key])
    ceil_h = read_float(user_text.get("ceilH"), DEFAULT_CEIL_H)
    return {
        "id": "room_{0}".format(index),
        "name": room_name(room_info, index),
        "cells": [[col, row] for col, row in sorted(cells)],
        "env": env,
        "ceilH": ceil_h,
    }


def build_room_payloads(room_infos, active_cells, origin, cell_size):
    if not room_infos:
        return []

    room_payloads = []
    sorted_active = sorted(active_cells)
    for index, room_info in enumerate(room_infos):
        cells = []
        for col, row in sorted_active:
            x = origin.X + (col + 0.5) * cell_size
            y = origin.Y + (row + 0.5) * cell_size
            if point_in_curve_info(x, y, room_info):
                cells.append((col, row))
        if cells:
            room_payloads.append(room_payload(room_info, index, cells))
    return room_payloads


def filter_features(raw_features, h_walls, v_walls, active_cells):
    payload = []
    seen = set()
    for feature_type, orient, row, col in sorted(raw_features):
        edge = (orient, row, col)
        if not edge_has_wall(edge, h_walls, v_walls):
            continue
        if not edge_touches_active(edge, active_cells):
            continue
        if feature_type == "exit" and not edge_is_exterior(edge, active_cells):
            continue
        key = (feature_type, orient, row, col)
        if key in seen:
            continue
        seen.add(key)
        payload.append(
            {
                "type": feature_type,
                "orient": orient,
                "r": row,
                "c": col,
            }
        )
    return payload


def export_geometry_json():
    warnings = []

    footprint_infos = curve_infos_on_layer(LAYER_NAMES["footprint"], True, warnings)
    if not footprint_infos:
        rs.MessageBox(
            "No closed footprint curves were found on layer '{0}'.".format(LAYER_NAMES["footprint"]),
            0,
            "thermBAL Rhino Export",
        )
        return

    void_infos = curve_infos_on_layer(LAYER_NAMES["voids"], True, warnings)
    room_infos = curve_infos_on_layer(LAYER_NAMES["rooms"], True, warnings)

    all_ids = []
    for layer_name in LAYER_NAMES.values():
        all_ids.extend(object_ids_on_layer(layer_name))

    origin = get_origin_point(all_ids)
    if origin is None:
        rs.MessageBox(
            "Could not determine an origin. Add geometry or place a point on layer '{0}'.".format(LAYER_NAMES["origin"]),
            0,
            "thermBAL Rhino Export",
        )
        return

    cell_size = prompt_cell_size()
    if not cell_size:
        return

    template_name = prompt_template_name()
    if not template_name:
        return

    bounds = bbox_min_max(all_ids)
    if not bounds:
        rs.MessageBox("No exportable geometry found.", 0, "thermBAL Rhino Export")
        return

    min_point, max_point = bounds
    if min_point.X < origin.X - sc.doc.ModelAbsoluteTolerance or min_point.Y < origin.Y - sc.doc.ModelAbsoluteTolerance:
        rs.MessageBox(
            "The origin must be at or below the lower-left corner of the exported geometry.",
            0,
            "thermBAL Rhino Export",
        )
        return

    grid = compute_grid(origin, max_point, cell_size)
    if grid is None:
        rs.MessageBox("Could not compute grid size.", 0, "thermBAL Rhino Export")
        return

    cols, rows = grid
    active_cells = sample_active_cells(
        cols,
        rows,
        origin,
        cell_size,
        footprint_infos,
        void_infos,
    )
    if not active_cells:
        rs.MessageBox(
            "Footprint sampling found no active cells. Check the cell size and footprint layer.",
            0,
            "thermBAL Rhino Export",
        )
        return

    h_walls, v_walls = empty_wall_matrices(cols, rows)
    build_exterior_walls(active_cells, h_walls, v_walls)
    add_wall_segments(
        object_ids_on_layer(LAYER_NAMES["walls"]),
        h_walls,
        v_walls,
        origin,
        cell_size,
        cols,
        rows,
        warnings,
    )

    feature_edges = set()
    feature_edges.update(
        feature_edges_from_layer(
            LAYER_NAMES["windows"],
            "window",
            origin,
            cell_size,
            cols,
            rows,
            warnings,
        )
    )
    feature_edges.update(
        feature_edges_from_layer(
            LAYER_NAMES["doors"],
            "door",
            origin,
            cell_size,
            cols,
            rows,
            warnings,
        )
    )
    feature_edges.update(
        feature_edges_from_layer(
            LAYER_NAMES["exits"],
            "exit",
            origin,
            cell_size,
            cols,
            rows,
            warnings,
        )
    )

    rooms = build_room_payloads(room_infos, active_cells, origin, cell_size)
    features = filter_features(feature_edges, h_walls, v_walls, active_cells)

    scale_to_m = model_unit_scale_to_meters()
    payload = {
        "meta": {
            "exported_at": datetime.utcnow().isoformat() + "Z",
            "mode": "geometry_template",
            "name": template_name,
            "grid_cols": cols,
            "grid_rows": rows,
            "cell_size_m": cell_size * scale_to_m,
            "source": "rhino_exporter",
            "source_model": sc.doc.Name or "",
            "rhino_version": str(Rhino.RhinoApp.ExeVersion),
            "model_units": str(sc.doc.ModelUnitSystem),
            "origin_model_units": [origin.X, origin.Y, origin.Z],
        },
        "active_cells": [[col, row] for col, row in sorted(active_cells)],
        "walls": {
            "horizontal": h_walls,
            "vertical": v_walls,
        },
        "features": features,
        "rooms": rooms,
    }

    path = prompt_output_path(template_name)
    if not path:
        return

    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)

    summary_lines = [
        "Exported thermBAL geometry JSON:",
        "  file: {0}".format(path),
        "  grid: {0} x {1}".format(cols, rows),
        "  cell size (model units): {0}".format(cell_size),
        "  cell size (m): {0}".format(round(cell_size * scale_to_m, 4)),
        "  active cells: {0}".format(len(active_cells)),
        "  features: {0}".format(len(features)),
        "  rooms: {0}".format(len(rooms)),
        "",
        "Load this file in simulator.html using 'Load Geometry JSON'.",
    ]

    if warnings:
        summary_lines.append("")
        summary_lines.append("Warnings:")
        for item in warnings:
            summary_lines.append("  - {0}".format(item))

    summary = "\n".join(summary_lines)
    log(summary)
    rs.MessageBox(summary, 0, "thermBAL Rhino Export")


if __name__ == "__main__":
    export_geometry_json()
