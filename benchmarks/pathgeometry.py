"""A generated tool path reduced to the polylines a drawing needs of it.

Curved primitives are resampled from their OWN parametrisation at a fixed
angular step, never read off the tessellated `polyline` a generator emits
alongside them: that tessellation is a property of the generator's sampling
settings rather than of the path, so a figure drawn from it shows one
generator's chording beside another's and calls the difference geometry.

Nothing here knows about colour or about matplotlib. A result is read only for
`.operations`, and an operation only for `.geometry`, `.operation` and
`.path_index`, so a path that was replayed, deserialised or built by hand
reduces exactly like one the kernel just emitted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any
from typing import List
from typing import Literal
from typing import Optional
from typing import Sequence
from typing import Tuple

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.tolerance import TOL

from benchmarks.errors import EmptyToolpathError
from benchmarks.errors import UnknownOperationClassError
from benchmarks.errors import UnplottableBoundaryError
from benchmarks.errors import UnplottableGeometryError
from compas_cgal.helix import HELIX_CLOSING_TURNS
from compas_cgal.helix import Helix

# Angular step a curved primitive is resampled at, chosen so that the chord is
# invisible for ANY radius rather than for the radii this corpus happens to draw.
# A chord of theta subtends a sagitta of r*(1 - cos(theta/2)); at two degrees that
# is 1.52e-4 r. The worst case a figure can hold is a circle filling the page, so
# r is at most half the 6.5 inch figure width, 82 mm, giving 0.013 mm of paper --
# a third of a pixel at 600 dpi. A full circle costs 180 segments, which is what
# a vector figure of a few hundred arcs is allowed to cost.
ARC_DEGREES_PER_SAMPLE = 2.0

# Even a very short arc keeps enough samples to read as curved rather than as a
# two-point chord.
MIN_ARC_SAMPLES = 8

# Padding around the drawn geometry, as a fraction of its longer side.
MARGIN_FRACTION = 0.04

# The operation classes a drawing names, one per value of the kernel's
# `OperationType` (pinned by `tests/benchmarks/test_helix_pathgeometry.py`).
OperationName = Literal["cut", "lead_in", "lead_out", "link", "retract", "plunge", "helix"]

# Operation classes whose motion removes material. Everything else is travel. A
# helix bores its footprint disk on the way down, so it cuts.
CUTTING_OPERATIONS: frozenset[OperationName] = frozenset({"cut", "lead_in", "lead_out", "helix"})

# Every operation class a drawing handles, which is exactly the set the kernel's
# `OperationType` defines. An unknown class is an error rather than a default
# stroke: a figure that silently draws an unrecognised motion as a cut is worse
# than one that refuses to be drawn.
KNOWN_OPERATIONS: frozenset[OperationName] = frozenset({"cut", "lead_in", "lead_out", "link", "retract", "plunge", "helix"})


@dataclass(frozen=True)
class Motion:
    """One operation reduced to what a drawing needs of it.

    Attributes:
        index: Position in the operation stream, which indexes `engagement_deg`.
        name: Operation class, lower case.
        path_index: The chain the operation belongs to.
        points: The motion resampled in the XY plane, in travel order.
        length: Summed length of the resampled polyline.
        is_point: True when the motion has no extent in XY -- a plunge, a
            retract, or a link between coincident ends. A helix has the extent
            of its footprint circle.
        is_cutting: True when the motion removes material.
    """

    index: int
    name: OperationName
    path_index: int
    points: Tuple[Tuple[float, float], ...]
    length: float
    is_point: bool
    is_cutting: bool


def drawable_motions(result: Any) -> Tuple[Motion, ...]:
    """Reduce a result's operations to drawable motions.

    Args:
        result: Anything carrying `.operations`.

    Returns:
        One motion per operation, in stream order.

    Raises:
        EmptyToolpathError: The result carries no operations.
        UnknownOperationClassError: An operation names an unknown class.
        UnplottableGeometryError: An operation carries a primitive with no path.
    """
    operations = list(result.operations)
    if not operations:
        raise EmptyToolpathError("This result carries no operations, so there is nothing to draw.")

    motions: List[Motion] = []
    for index, operation in enumerate(operations):
        name = _operation_name(index, operation)
        points = _sample(index, operation.geometry)
        spread = max(math.hypot(x - points[0][0], y - points[0][1]) for x, y in points)
        motions.append(
            Motion(
                index=index,
                name=name,
                path_index=int(operation.path_index),
                points=points,
                length=_polyline_length(points),
                is_point=TOL.is_zero(spread),
                is_cutting=name in CUTTING_OPERATIONS,
            )
        )
    return tuple(motions)


def _operation_name(index: int, operation: Any) -> OperationName:
    """The operation's class as a lower-case string.

    Args:
        index: Position in the stream, for the error message.
        operation: The operation, whose `.operation` is a string or a str-enum.

    Returns:
        The class name.

    Raises:
        UnknownOperationClassError: The class is not one this module draws.
    """
    raw = operation.operation
    name = str(getattr(raw, "value", raw)).lower()
    for known in KNOWN_OPERATIONS:
        if known == name:
            return known
    raise UnknownOperationClassError(f"Operation {index} is of class {name!r}; this module draws {sorted(KNOWN_OPERATIONS)}.")


def _sample(index: int, geometry: Any) -> Tuple[Tuple[float, float], ...]:
    """Resample one primitive into an XY polyline, in travel order.

    Curved primitives are sampled from their own parametrisation, so a circle is
    drawn as a circle whatever tessellation the generator emitted alongside it.
    A helix is drawn as its plan view: every turn of its footprint circle, the
    flat closing turn included, in its own turn sense (`_sample_helix`).

    Args:
        index: Position in the stream, for the error message.
        geometry: A compas line, arc or circle, or a `compas_cgal.helix.Helix`.

    Returns:
        The sampled points.

    Raises:
        UnplottableGeometryError: The primitive carries no path to sample.
    """
    if isinstance(geometry, Helix):
        return _sample_helix(geometry)
    if isinstance(geometry, Circle):
        count = _arc_samples(2.0 * math.pi)
        return tuple(_xy(geometry.point_at(step / count)) for step in range(count + 1))
    if isinstance(geometry, Arc):
        count = _arc_samples(geometry.angle)
        return tuple(_xy(geometry.point_at(step / count)) for step in range(count + 1))
    if isinstance(geometry, Line):
        return (_xy(geometry.start), _xy(geometry.end))
    raise UnplottableGeometryError(f"Operation {index} carries geometry of type {type(geometry).__name__!r}, which has no drawable path.")


def _sample_helix(helix: Helix) -> Tuple[Tuple[float, float], ...]:
    """The plan view of a helix op, from its own parametrisation, in travel order.

    Station ``t`` in ``[0, 1]`` is the exact entry ``P`` turned about the
    footprint centre by ``t * (turns + HELIX_CLOSING_TURNS)`` full turns in the
    helix's sense: the ramp turns, then the flat closing turn at ``bottom_z``
    (user ruling Q6). The descent is along Z, which a plan view drops
    (`helix_elevation` draws it). The two ends are the exact entry and exit
    (the same point: every turn is whole), not a rounding of the last
    rotation. REPORTING only: a drawing decides nothing.

    Args:
        helix: The helix.

    Returns:
        The sampled points, ``turns + HELIX_CLOSING_TURNS`` times round the
        footprint circle.
    """
    turn = (-math.tau if helix.clockwise else math.tau) * (helix.turns + HELIX_CLOSING_TURNS)
    count = _arc_samples(turn)
    cx, cy = float(helix.footprint.centre.x), float(helix.footprint.centre.y)
    entry = (float(helix.entry_xy.x), float(helix.entry_xy.y))
    exit_ = (float(helix.exit_xy.x), float(helix.exit_xy.y))
    dx, dy = entry[0] - cx, entry[1] - cy
    interior = []
    for step in range(1, count):
        angle = turn * step / count
        cosine, sine = math.cos(angle), math.sin(angle)
        interior.append((cx + cosine * dx - sine * dy, cy + sine * dx + cosine * dy))
    return (entry, *interior, exit_)


def helix_elevation(helix: Helix) -> Tuple[Tuple[float, float], ...]:
    """The side elevation of a helix op: Z against plan arc length, in travel order.

    The ramp descends linearly from ``top_z`` to ``bottom_z`` over its plan run,
    ``turns`` footprint circumferences (the plan projection of
    `Helix.ramp_length`, which is the 3D length); the closing turn then runs
    flat at ``bottom_z`` for the helix's own `Helix.closing_turn_length` (user
    ruling Q6). Both pieces are straight in this chart, so their three corners
    are the whole profile; the heights are the op's exact top and bottom.
    REPORTING only.

    Args:
        helix: The helix.

    Returns:
        ``(arc length mm, z mm)`` at the ramp start, the ramp end and the
        closing turn's end.
    """
    ramp_plan = helix.turns * math.tau * float(helix.footprint.radius.value)
    top, bottom = float(helix.top_z.value), float(helix.bottom_z.value)
    return (
        (0.0, top),
        (ramp_plan, bottom),
        (ramp_plan + float(helix.closing_turn_length()), bottom),
    )


def _arc_samples(angle: float) -> int:
    """How many segments an arc of *angle* radians is drawn with.

    Args:
        angle: Swept angle in radians; its sign does not change the count.

    Returns:
        The segment count, at least `MIN_ARC_SAMPLES`.
    """
    return max(MIN_ARC_SAMPLES, int(math.ceil(abs(math.degrees(angle)) / ARC_DEGREES_PER_SAMPLE)))


def _xy(point: Any) -> Tuple[float, float]:
    """The XY part of a point, as plain floats.

    Args:
        point: Anything indexable by 0 and 1.

    Returns:
        The pair.
    """
    return (float(point[0]), float(point[1]))


def _polyline_length(points: Sequence[Tuple[float, float]]) -> float:
    """Summed length of a polyline.

    Args:
        points: The polyline's vertices.

    Returns:
        The length; zero for a single point.
    """
    return float(sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(points, points[1:])))


def closed_rings(boundary: Any, holes: Sequence[Any]) -> Tuple[Tuple[Tuple[float, float], ...], ...]:
    """Close the pocket boundary and its islands into drawable rings.

    Args:
        boundary: The outer boundary.
        holes: The islands.

    Returns:
        The outer ring first, then one ring per island.

    Raises:
        UnplottableBoundaryError: A ring has fewer than three points.
    """
    return tuple(_ring(polygon) for polygon in [boundary, *holes])


def _ring(polygon: Any) -> Tuple[Tuple[float, float], ...]:
    """One closed ring from a polygon or a sequence of points.

    Args:
        polygon: A compas polygon, or any sequence of points.

    Returns:
        The ring, with its first point repeated at the end.

    Raises:
        UnplottableBoundaryError: The ring has fewer than three points.
    """
    points = tuple(_xy(point) for point in getattr(polygon, "points", polygon))
    if len(points) < 3:
        raise UnplottableBoundaryError(f"A pocket boundary needs at least three points; got {len(points)}.")
    closed = TOL.is_zero(math.hypot(points[-1][0] - points[0][0], points[-1][1] - points[0][1]))
    return points if closed else points + (points[0],)


def view_bounds(
    rings: Sequence[Sequence[Tuple[float, float]]],
    motions: Sequence[Motion],
    envelope_diameter: Optional[float],
) -> Tuple[float, float, float, float]:
    """The view, padded so that no stroke touches the panel edge.

    Args:
        rings: The pocket rings.
        motions: Every motion that will be drawn.
        envelope_diameter: Tool diameter when the swept envelope is drawn, so the
            padding covers the half width the envelope adds outside the centre
            path; None otherwise.

    Returns:
        ``(xmin, xmax, ymin, ymax)``.
    """
    xs = [x for ring in rings for x, _ in ring] + [x for motion in motions for x, _ in motion.points]
    ys = [y for ring in rings for _, y in ring] + [y for motion in motions for _, y in motion.points]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    margin = MARGIN_FRACTION * max(xmax - xmin, ymax - ymin)
    if envelope_diameter is not None:
        margin += 0.5 * envelope_diameter
    return (xmin - margin, xmax + margin, ymin - margin, ymax + margin)
