# Tangent Continuity in Trochoidal Toolpaths

A tangent discontinuity during material engagement means the CNC controller
must instantaneously change feed direction — causing jerk, tool deflection,
and potential machine damage.  This document describes which transitions must
be checked, which are structurally exempt, and how the Python representation
affects the measurement.

## Operation classes

The toolpath operation stream contains six types:

| Operation | Tool engaged? | Tangent-continuous with neighbours? |
|-----------|:------------:|:-----------------------------------:|
| `cut`     | yes          | depends (see below)                 |
| `lead_in` | yes          | must be, at cut boundary            |
| `lead_out`| yes          | must be, at cut boundary            |
| `link`    | no           | not required (at clearance Z)       |
| `retract` | no           | not required                        |
| `plunge`  | no           | not required                        |

## Trochoid chain structure

Within a single MAT edge, the cut sequence alternates:

```
 ╭── Arc/Circle (cutting circle, ≈360°)
 │      ╰── Line (bridge step along MAT edge)
 │              ╰── Arc/Circle (next cutting circle)
 │                      ╰── Line (bridge step) ...
```

The bridge lines march along the MAT edge direction.  The circle tangent
at the bridge junction is perpendicular to the MAT edge.  This ~90° angle
is inherent to trochoidal milling geometry — not a defect:

```
        MAT edge direction →

    bridge line        bridge line
    ──────────•──────────•──────────
              │          │
          ○   │      ○   │
         ╱ ╲  │     ╱ ╲  │
        │   │ ↑    │   │ ↑   ← circle tangent
         ╲ ╱       ╲ ╱       (perpendicular)
          ○         ○
```

CNC controllers handle these short direction changes through look-ahead
deceleration.  The bridge distance is comparable to stepover.

## Consistent winding

All cut arcs and circles use the same winding direction (CW for climb
milling).  At varying-radius stations — where adjacent trochoid circles
have different radii — the tangent arrival and departure points on a
circle do not coincide.  The toolpath emits a full circle followed by a
short repositioning arc, both at nominal winding, rather than a single
arc with flipped winding.  This trades a small amount of re-cutting
(over the gap region) for consistent climb/conventional milling
throughout the toolpath.

### Why this matters

Without consistent winding, the toolpath alternates between climb and
conventional milling at every varying-radius station.  This causes:

- **Alternating cutting forces** — climb milling pulls the tool into the
  material, conventional pushes it away.  Rapid alternation induces
  vibration and chatter.
- **Degraded surface finish** — each direction change leaves a witness
  mark.  Consistent climb milling produces a uniform finish.
- **Reduced tool life** — conventional milling rubs before cutting,
  generating heat.  Consistent climb avoids this.

The original implementation used a single arc with winding flipped to
guarantee coverage (>π sweep).  The fix splits this into full circle +
repositioning arc, both at nominal winding.  Empirically, ~60% of arcs
were flipped — nearly all were near-full circles (gap < 1°) where the
flip was caused by construction noise in the tangent departure point.

### Construction noise gate

The decision to emit a repositioning arc uses
`CGAL::compare_squared_distance` (an exact predicate) with a threshold
of `1e-18` (i.e. `1e-9²`).  Points closer than this are treated as
coincident, collapsing to a single full circle.  This avoids mixing
double-precision proximity checks with exact predicates — the gate stays
in exact arithmetic throughout.

## ToolpathOperation fields

Each `ToolpathOperation` carries the full C++ metadata:

| Field | Type | Description |
|-------|------|-------------|
| `geometry` | `Line \| Arc \| Circle` | compas geometry primitive |
| `operation` | `str` | `"cut"`, `"lead_in"`, `"lead_out"`, `"link"`, `"retract"`, `"plunge"` |
| `path_index` | `int` | MAT edge group index |
| `clockwise` | `bool` | winding direction from C++ (ground truth) |
| `start_tangent` | `ndarray` | unit tangent vector at start (winding-correct, z=0) |
| `end_tangent` | `ndarray` | unit tangent vector at end (winding-correct, z=0) |

The `start_tangent` and `end_tangent` are computed in C++ via
`Vector_2::perpendicular(orientation)` — exact CGAL kernel operations —
then normalized to unit length.  For lines, both equal the unit
direction vector.  For arcs/circles, they are perpendicular to the
radius at each endpoint, respecting the winding direction.

To compute a tangent at an arbitrary parameter `t` on an arc/circle in
Python, use `clockwise` to select the correct perpendicular:

```python
rx, ry = point.x - center.x, point.y - center.y
tx, ty = (ry, -rx) if op.clockwise else (-ry, rx)
```

## Winding representation gap

The C++ backend emits circles with explicit CW/CCW orientation.
`ToolpathOperation.clockwise` preserves this flag for downstream code.
The compas geometry types have a partial representation gap:

| C++ type | Python type | Winding preserved? |
|----------|-------------|--------------------|
| CW arc   | `Arc` (start_angle > end_angle) | yes |
| CCW arc  | `Arc` (start_angle < end_angle) | yes |
| CW circle| `Circle`    | **no** — always CCW |
| CCW circle| `Circle`   | yes |

Because compas `Circle` has no winding attribute, CW circles appear as
CCW in the Python representation.  Use `ToolpathOperation.clockwise`
(from the C++ metadata) for the ground truth.  The pre-computed
`start_tangent` and `end_tangent` unit vectors are also winding-correct.

## Testing strategy

The tangent continuity test uses three principles:

### 1. Analytical tangent computation

Tangent at any point on an arc/circle is the perpendicular to the radius
vector — exact, no discretization:

```python
rx = point.x - center.x
ry = point.y - center.y
tangent = (-ry, rx)  # or (ry, -rx) for CW; sign irrelevant with abs(dot)
```

### 2. Parallelism, not directionality

Because the winding representation may flip the tangent sign, the test
checks **parallelism** (angle between tangent lines, 0°–90°) rather than
co-directionality (angle between tangent vectors, 0°–180°):

```python
cos_a = abs(dot(exit_tangent, entry_tangent))  # abs handles winding flip
angle = acos(clamp(cos_a, 0, 1))
assert angle < 1.0  # degrees
```

This catches physical jerk (45°, 90° direction changes) while tolerating
the harmless 180° winding artifact.

### 3. Bridge exclusion

Line ↔ Arc/Circle transitions within the same path_index and both labelled
`cut` are trochoid bridge steps — structurally ~90° and exempt from the
tangent check.

### What the test catches

| Scenario | Measured angle | Test result |
|----------|:--------------:|:-----------:|
| Smooth arc → arc transition | < 1° | pass |
| Winding flip (CW → CCW repr) | 0° (via abs) | pass |
| Bridge line ↔ circle | ~90° | skipped |
| Missing bridge at chain junction | 30°–90° | **fail** |
| Lead-in misaligned with first cut | > 1° | **fail** |
| Broken arc encoding | arbitrary | **fail** |

## C++ tessellated polyline

The `trochoidal_mat_toolpath_circular` function now returns a dense 3D polyline
tessellated directly from the internal `ToolpathPrimitive` vector.  Because this
tessellation uses the `signed_sweep` from each `TrochoidArc` — not the compas
geometry representation — it preserves the correct traversal order for both CW
and CCW arcs.

The polyline bypasses the winding representation gap entirely: no Python-side
tessellation, no heuristic segment reversal, no `to_polyline()` ambiguity.
Angular resolution is controlled by `samples_per_radian` (default 10.0, giving
~63 points per full circle).

Use `result.polyline` for visualization and `result.operations` for G-code
generation.
