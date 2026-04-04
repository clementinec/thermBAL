# Rhino Export

Use [thermbal_rhino_exporter.py](thermbal_rhino_exporter.py) to export a Rhino floor plan directly into the `geometry.json` format that `simulator.html` can load.

## What It Exports

- `active_cells`: occupiable footprint cells
- `walls`: wall-edge matrices
- `features`: windows, doors, exits
- `rooms`: optional named room regions with default environment values

The exporter auto-generates the exterior wall ring from the footprint. The `tb_walls` layer is mainly for interior partitions.

## Recommended Rhino Layers

- `tb_footprint`
  Closed planar curves that outline the occupiable building footprint.
- `tb_voids`
  Optional closed planar curves to subtract courtyards, shafts, or holes from the footprint.
- `tb_walls`
  Optional interior wall centerlines as horizontal/vertical lines or polylines.
- `tb_windows`
  Optional horizontal/vertical lines or polylines for windows.
- `tb_doors`
  Optional horizontal/vertical lines or polylines for doors.
- `tb_exits`
  Optional horizontal/vertical lines or polylines for exits.
- `tb_rooms`
  Optional closed planar curves for room naming and metadata.
- `tb_origin`
  Optional point or text dot marking the lower-left grid origin. If omitted, the exporter uses the bounding-box minimum corner.

## Geometry Assumptions

- Best for plans aligned to World XY.
- Best for orthogonal geometry.
- Curved or angled footprint boundaries can still rasterize into the footprint.
- Curved or angled interior walls/openings are not exported cleanly and will be skipped.
- Wall/opening curves should already be snapped to the intended grid.

## How To Use In Rhino

1. Put your geometry on the layers above.
2. Open Rhino's Python editor.
3. Paste in [thermbal_rhino_exporter.py](thermbal_rhino_exporter.py), or load the file directly.
4. Run the script.
5. Enter the grid cell size in model units.
6. Enter a template name.
7. Save the output JSON.

Then open [simulator.html](../simulator.html) and click `Load Geometry JSON`.

## Room Names And Metadata

If you put closed curves on `tb_rooms`, the exporter will sample those curves and attach room payloads.

Room naming priority:

- Rhino object name
- user text key `name`
- fallback `Room A`, `Room B`, ...

Optional user text keys on room curves:

- `name`
- `ceilH`
- `air_temp`
- `mean_radiant_temp`
- `humidity`
- `air_velocity`
- `lux`
- `noise_db`

If these keys are absent, the exporter uses the same defaults as the geometry importer.

## Unit Handling

The script asks for the cell size in Rhino model units, then converts that to meters for `cell_size_m` using the Rhino document unit system.

Examples:

- model in meters, cell size `0.5` -> exports `cell_size_m = 0.5`
- model in millimeters, cell size `500` -> exports `cell_size_m = 0.5`

## Practical Workflow

1. Author or simplify geometry in Rhino.
2. Export geometry JSON with this script.
3. Load that JSON into [simulator.html](../simulator.html).
4. Set room conditions and occupant profiles in the simulator.
5. Export the simulator state for downstream batch studies.

## Caveats

- The exporter does not attempt to infer a grid rotation.
- It does not convert arbitrary Rhino Breps/meshes into walls.
- It assumes the Rhino drawing has already been simplified to the level the simulator needs.
