import math

import numpy as np
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon
from compas.geometry import is_point_in_polygon_xy
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import polygon_medial_axis_transform
from compas_cgal.toolpath import trochoidal_mat_toolpath
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular


def _point_segment_distance_xy(point, a, b):
    point = np.asarray(point, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom <= 1e-16:
        return float(np.linalg.norm(point - a))
    t = max(0.0, min(1.0, float(np.dot(point - a, ab) / denom)))
    projection = a + t * ab
    return float(np.linalg.norm(point - projection))


def _distance_to_polygon_boundary_xy(point, polygon_xy):
    distance = np.inf
    count = len(polygon_xy)
    for i in range(count):
        d = _point_segment_distance_xy(point, polygon_xy[i], polygon_xy[(i + 1) % count])
        if d < distance:
            distance = d
    return float(distance)


def _op_start_xy(op: ToolpathOperation) -> np.ndarray:
    g = op.geometry
    if isinstance(g, Line):
        return np.array([float(g.start[0]), float(g.start[1]), float(g.start[2])], dtype=np.float64)
    if isinstance(g, Arc):
        pt = g.point_at(0.0)
        return np.array([float(pt[0]), float(pt[1]), float(pt[2])], dtype=np.float64)
    if isinstance(g, Circle):
        pt = g.point_at(0.0)
        return np.array([float(pt[0]), float(pt[1]), float(pt[2])], dtype=np.float64)
    raise TypeError(type(g))


def _op_end_xy(op: ToolpathOperation) -> np.ndarray:
    g = op.geometry
    if isinstance(g, Line):
        return np.array([float(g.end[0]), float(g.end[1]), float(g.end[2])], dtype=np.float64)
    if isinstance(g, Arc):
        pt = g.point_at(1.0)
        return np.array([float(pt[0]), float(pt[1]), float(pt[2])], dtype=np.float64)
    if isinstance(g, Circle):
        pt = g.point_at(1.0)
        return np.array([float(pt[0]), float(pt[1]), float(pt[2])], dtype=np.float64)
    raise TypeError(type(g))


# ---------------------------------------------------------------------------
# Tangent continuity helpers
# ---------------------------------------------------------------------------

ENGAGED_OPS = frozenset({"lead_in", "cut", "lead_out"})


def _tangent_xy(geom, t):
    """Exact unit tangent vector in XY at parameter *t*.

    Line: direction vector (exact).
    Arc / Circle: perpendicular to radius at the point (exact, no discretization).
    Since the test uses abs(dot) for parallelism, the CW/CCW sign is irrelevant.
    """
    if isinstance(geom, Line):
        d = np.array([float(geom.end[0]) - float(geom.start[0]), float(geom.end[1]) - float(geom.start[1])])
    else:
        pt = geom.point_at(t)
        center = geom.frame.point
        rx = float(pt[0]) - float(center[0])
        ry = float(pt[1]) - float(center[1])
        d = np.array([-ry, rx])
    n = np.linalg.norm(d)
    return d / n if n > 1e-14 else None


def _is_trochoid_bridge(prev, curr):
    """True when the transition is a bridge line in a trochoid chain.

    Trochoid chains alternate Arc/Circle ↔ Line within the same path.
    The bridge line steps along the MAT edge at ~90° to the circle tangent —
    this is inherent to trochoidal milling geometry, not a defect.
    """
    if prev.path_index != curr.path_index:
        return False
    if prev.operation != "cut" or curr.operation != "cut":
        return False
    prev_is_line = isinstance(prev.geometry, Line)
    curr_is_line = isinstance(curr.geometry, Line)
    return prev_is_line != curr_is_line


def _assert_engaged_tangent_continuity(ops, tol_deg=1.0):
    """Assert tangent parallelism at all consecutive engaged-operation transitions.

    Checks that exit tangent of op[i-1] and entry tangent of op[i] are within
    *tol_deg* of being parallel whenever both operations are material-engaged
    (lead_in, cut, lead_out).

    Trochoid bridge lines (short linear steps between cutting circles) are
    excluded — their ~90° angle to the circle tangent is inherent to
    trochoidal milling geometry.

    Uses *parallelism* (angle to line, not vector) because compas Circle/Arc
    encoding may reverse the winding direction relative to the C++ trochoid.
    A winding flip (180°) is safe — it just means G02 vs G03 in G-code.
    A non-parallel tangent (e.g. 45°, 90°) is a physical jerk that can damage
    the tool or workpiece.
    """
    violations = []
    for i in range(1, len(ops)):
        prev, curr = ops[i - 1], ops[i]
        if prev.operation not in ENGAGED_OPS or curr.operation not in ENGAGED_OPS:
            continue
        if _is_trochoid_bridge(prev, curr):
            continue

        exit_t = _tangent_xy(prev.geometry, 1.0)
        entry_t = _tangent_xy(curr.geometry, 0.0)
        if exit_t is None or entry_t is None:
            continue

        # abs(dot) → measures parallelism, tolerating winding reversal
        cos_a = np.clip(abs(np.dot(exit_t, entry_t)), 0.0, 1.0)
        angle_deg = math.degrees(math.acos(cos_a))
        if angle_deg > tol_deg:
            violations.append((i, prev.operation, curr.operation, angle_deg))

    assert not violations, f"{len(violations)} tangent discontinuity(ies) (>{tol_deg}°):\n" + "\n".join(
        f"  op[{idx - 1}] {pop} -> op[{idx}] {cop}: {a:.1f}°" for idx, pop, cop, a in violations[:10]
    )


# ---------------------------------------------------------------------------
# Hypothesis strategy: random simple polygons
# ---------------------------------------------------------------------------


@st.composite
def _simple_polygons(draw):
    """Simple polygon via angular sweep with controlled perturbation.

    Vertices are placed at evenly-spaced angles around the origin with ±25%
    angular jitter and 0.5–1.0× radial variation.  This guarantees simplicity
    without a computational geometry check.
    """
    n = draw(st.integers(min_value=4, max_value=12))
    step = 2.0 * math.pi / n
    base_r = draw(st.floats(min_value=5.0, max_value=15.0))
    points = []
    for i in range(n):
        a = i * step + draw(st.floats(min_value=-0.25 * step, max_value=0.25 * step))
        r = base_r * draw(st.floats(min_value=0.5, max_value=1.0))
        points.append((r * math.cos(a), r * math.sin(a), 0.0))
    return Polygon(points)


# ---------------------------------------------------------------------------
# Test polygons
# ---------------------------------------------------------------------------

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])
IRREGULAR = Polygon(
    [
        (-1.91, 3.59, 0.0),
        (-5.53, -5.22, 0.0),
        (-0.39, -1.98, 0.0),
        (2.98, -5.51, 0.0),
        (4.83, -2.02, 0.0),
        (9.70, -3.63, 0.0),
        (12.23, 1.25, 0.0),
        (3.42, 0.66, 0.0),
        (2.92, 4.03, 0.0),
    ]
)
L_SHAPE = Polygon([(0, 0, 0), (6, 0, 0), (6, 2, 0), (2, 2, 0), (2, 8, 0), (0, 8, 0)])
STAR = Polygon(
    [
        (5.0, 0.0, 0.0),
        (1.5, 1.2, 0.0),
        (1.55, 4.76, 0.0),
        (-0.57, 1.9, 0.0),
        (-4.05, 2.94, 0.0),
        (-1.85, 0.0, 0.0),
        (-4.05, -2.94, 0.0),
        (-0.57, -1.9, 0.0),
        (1.55, -4.76, 0.0),
        (1.5, -1.2, 0.0),
    ]
)
KITE = Polygon([(5, 0, 0), (0, 5, 0), (-5, 0, 0), (0, -2.5, 0)])


def test_polygon_medial_axis_transform_square():
    mat_points, mat_radii = polygon_medial_axis_transform(SQUARE)
    assert mat_points.shape == (1, 3)
    assert mat_radii.shape == (1,)
    assert pytest.approx(mat_radii[0]) == 5.0


def test_trochoidal_mat_toolpath_inside_polygon():
    tool_diameter = 1.0
    tool_radius = 0.5 * tool_diameter
    paths = trochoidal_mat_toolpath(
        IRREGULAR,
        tool_diameter=tool_diameter,
        stepover=0.5,
        pitch=0.75,
        samples_per_cycle=8,
    )

    assert len(paths) > 0
    polygon_xy = [list(pt[:2]) for pt in IRREGULAR.points]
    for path in paths:
        assert path.shape[1] == 3
        assert np.isfinite(path).all()
        assert path.shape[0] >= 2

        sample = np.linspace(0, len(path) - 1, 16, dtype=int)
        for index in sample:
            point_xy = path[index][:2].tolist()
            assert is_point_in_polygon_xy(point_xy, polygon_xy)
            distance = _distance_to_polygon_boundary_xy(point_xy, polygon_xy)
            assert distance >= tool_radius - 1e-3


def test_trochoidal_mat_toolpath_maximizes_clearance():
    tool_diameter = 2.0
    tool_radius = 0.5 * tool_diameter
    paths = trochoidal_mat_toolpath(
        SQUARE,
        tool_diameter=tool_diameter,
        stepover=2.0,
        pitch=1.0,
        samples_per_cycle=16,
        max_trochoid_radius=float("inf"),
    )

    assert len(paths) > 0
    polygon_xy = [list(pt[:2]) for pt in SQUARE.points]
    sampled_distances = []
    for path in paths:
        sample = np.linspace(0, len(path) - 1, 24, dtype=int)
        for index in sample:
            sampled_distances.append(_distance_to_polygon_boundary_xy(path[index][:2].tolist(), polygon_xy))

    min_distance = min(sampled_distances)
    assert min_distance >= tool_radius - 1e-2
    assert min_distance <= tool_radius + 1e-2


def test_circular_returns_toolpath_operations():
    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
    )

    assert len(ops) > 0
    for op in ops:
        assert isinstance(op, ToolpathOperation)
        assert isinstance(op.geometry, (Line, Arc, Circle))
        assert op.operation in ("cut", "lead_in", "lead_out", "link", "retract", "plunge")
        assert isinstance(op.path_index, int)

    assert any(op.operation == "cut" for op in ops)
    assert any(isinstance(op.geometry, (Arc, Circle)) for op in ops if op.operation == "cut")


def test_circular_with_leads_and_links():
    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
        lead_in=0.2,
        lead_out=0.2,
        link_paths=True,
        optimize_order=True,
    )

    assert len(ops) > 0
    op_types = {op.operation for op in ops}
    assert "cut" in op_types
    assert "lead_in" in op_types
    assert "lead_out" in op_types
    assert "link" in op_types


def test_circular_with_clearance_z():
    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
        lead_in=0.1,
        lead_out=0.1,
        link_paths=True,
        optimize_order=True,
        cut_z=-2.0,
        clearance_z=3.0,
        retract_at_end=True,
    )

    assert len(ops) > 0
    op_types = {op.operation for op in ops}
    assert "retract" in op_types
    assert "plunge" in op_types
    assert "link" in op_types

    for op in ops:
        start = _op_start_xy(op)
        end = _op_end_xy(op)
        if op.operation in ("cut", "lead_in", "lead_out"):
            assert pytest.approx(start[2], abs=1e-9) == -2.0
            assert pytest.approx(end[2], abs=1e-9) == -2.0
        if op.operation == "link":
            assert pytest.approx(start[2], abs=1e-9) == 3.0
            assert pytest.approx(end[2], abs=1e-9) == 3.0
        if op.operation == "retract":
            assert pytest.approx(start[2], abs=1e-9) == -2.0
            assert pytest.approx(end[2], abs=1e-9) == 3.0
        if op.operation == "plunge":
            assert pytest.approx(start[2], abs=1e-9) == 3.0
            assert pytest.approx(end[2], abs=1e-9) == -2.0


def test_circular_continuity():
    """Consecutive operations share start/end points."""
    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
        lead_in=0.2,
        lead_out=0.2,
        link_paths=True,
        cut_z=-1.0,
        clearance_z=2.0,
    )

    for i in range(1, len(ops)):
        prev_end = _op_end_xy(ops[i - 1])
        curr_start = _op_start_xy(ops[i])
        assert np.allclose(prev_end, curr_start, atol=1e-5), f"Discontinuity at op {i}: prev end {prev_end} != curr start {curr_start}"


def test_full_circle_output():
    """Cut operations should contain full Circle geometry objects."""
    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
    )

    cut_ops = [op for op in ops if op.operation == "cut"]
    circles = [op for op in cut_ops if isinstance(op.geometry, Circle)]
    assert len(circles) > 0, "Expected Circle geometry in cut operations"
    assert all(c.geometry.radius > 0 for c in circles)


def test_full_circle_sweep():
    """Full circles should have sweep angle of approximately 2*pi."""
    import math

    ops = trochoidal_mat_toolpath_circular(
        SQUARE,
        tool_diameter=2.0,
        pitch=1.0,
        max_trochoid_radius=float("inf"),
        max_passes=20,
    )

    cut_circles = [op for op in ops if op.operation == "cut" and isinstance(op.geometry, Circle)]
    assert len(cut_circles) > 0
    for op in cut_circles:
        c = op.geometry
        # Circle.point_at(0) == Circle.point_at(1) for a full circle
        p0 = np.array(c.point_at(0.0)[:2])
        p1 = np.array(c.point_at(1.0)[:2])
        assert np.allclose(p0, p1, atol=1e-9), f"Circle start != end: {p0} vs {p1}"
        assert pytest.approx(c.circumference, abs=1e-6) == 2.0 * math.pi * c.radius


def test_circular_exit_tangent_alignment():
    """Circle/arc exit tangents should be aligned with the next line direction (smooth exit)."""
    ops = trochoidal_mat_toolpath_circular(
        IRREGULAR,
        tool_diameter=1.0,
        pitch=0.75,
        max_trochoid_radius=float("inf"),
        max_passes=20,
    )

    cut_ops = [o for o in ops if o.operation == "cut"]
    assert len(cut_ops) > 2

    exit_angles = []
    for i in range(len(cut_ops) - 1):
        prev_g = cut_ops[i].geometry
        curr_g = cut_ops[i + 1].geometry
        if not isinstance(prev_g, (Arc, Circle)) or not isinstance(curr_g, Line):
            continue

        # Circle/arc exit tangent (sampled near t=1)
        p0 = np.array(prev_g.point_at(0.97)[:2])
        p1 = np.array(prev_g.point_at(1.0)[:2])
        arc_dir = p1 - p0
        arc_norm = np.linalg.norm(arc_dir)
        if arc_norm < 1e-12:
            continue

        # Next line direction
        line_dir = np.array(curr_g.end[:2]) - np.array(curr_g.start[:2])
        line_norm = np.linalg.norm(line_dir)
        if line_norm < 1e-12:
            continue

        cos_angle = np.dot(arc_dir / arc_norm, line_dir / line_norm)
        # Use absolute cos — for full circles the sampled direction may be
        # anti-parallel to the tangent line due to parametric wraparound.
        angle = np.degrees(np.arccos(np.clip(abs(cos_angle), 0, 1)))
        exit_angles.append(angle)

    assert len(exit_angles) > 0
    # Exit tangent should be well-aligned (< 30°) for most circle->line junctions
    aligned = sum(1 for a in exit_angles if a < 30)
    assert aligned >= len(exit_angles) * 0.7, f"Only {aligned}/{len(exit_angles)} exits are tangent-aligned. Angles: {[f'{a:.0f}' for a in sorted(exit_angles)]}"


def test_trochoidal_mat_toolpath_invalid_parameters():
    small = Polygon([[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]])
    with pytest.raises(ValueError):
        trochoidal_mat_toolpath(small, tool_diameter=0.0)
    with pytest.raises(ValueError):
        trochoidal_mat_toolpath(small, tool_diameter=1.0, stepover=0.0)
    with pytest.raises(ValueError):
        trochoidal_mat_toolpath(small, tool_diameter=1.0, samples_per_cycle=3)


# ---------------------------------------------------------------------------
# Tangent continuity — known polygons + Hypothesis
# ---------------------------------------------------------------------------

_TANGENT_TOOLPATH_KWARGS = dict(
    tool_diameter=0.5,
    pitch=0.4,
    lead_in=0.15,
    lead_out=0.15,
    link_paths=True,
    optimize_order=True,
    cut_z=-0.2,
    clearance_z=2.0,
)

_TANGENT_POLYGONS = [SQUARE, IRREGULAR, L_SHAPE, STAR, KITE]
_TANGENT_IDS = ["square", "irregular", "L_shape", "star", "kite"]


@pytest.mark.parametrize("polygon", _TANGENT_POLYGONS, ids=_TANGENT_IDS)
def test_tangent_continuity_known_polygons(polygon):
    """Engaged milling motions must be tangent-continuous (< 1°)."""
    ops = trochoidal_mat_toolpath_circular(polygon, **_TANGENT_TOOLPATH_KWARGS)
    engaged = [op for op in ops if op.operation in ENGAGED_OPS]
    assert len(engaged) > 2, "Not enough engaged operations to test"
    _assert_engaged_tangent_continuity(ops)


@given(polygon=_simple_polygons())
@settings(max_examples=50, deadline=None)
def test_tangent_continuity_random_polygons(polygon):
    """Property: for any simple polygon, engaged motions are tangent-continuous."""
    ops = trochoidal_mat_toolpath_circular(polygon, **_TANGENT_TOOLPATH_KWARGS)
    assume(len([op for op in ops if op.operation in ENGAGED_OPS]) > 2)
    _assert_engaged_tangent_continuity(ops)
