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

# Operation classes whose motion removes material. Everything else is travel.
CUTTING_OPERATIONS = frozenset({"cut", "lead_in", "lead_out"})

# Every operation class a drawing handles, which is exactly the set the kernel's
# `OperationType` defines. An unknown class is an error rather than a default
# stroke: a figure that silently draws an unrecognised motion as a cut is worse
# than one that refuses to be drawn.
KNOWN_OPERATIONS = frozenset({"cut", "lead_in", "lead_out", "link", "retract", "plunge"})


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
            retract, or a link between coincident ends.
        is_cutting: True when the motion removes material.
    """

    index: int
    name: str
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


def _operation_name(index: int, operation: Any) -> str:
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
    if name not in KNOWN_OPERATIONS:
        raise UnknownOperationClassError(f"Operation {index} is of class {name!r}; this module draws {sorted(KNOWN_OPERATIONS)}.")
    return name


def _sample(index: int, geometry: Any) -> Tuple[Tuple[float, float], ...]:
    """Resample one primitive into an XY polyline, in travel order.

    Curved primitives are sampled from their own parametrisation, so a circle is
    drawn as a circle whatever tessellation the generator emitted alongside it.

    Args:
        index: Position in the stream, for the error message.
        geometry: A compas line, arc or circle.

    Returns:
        The sampled points.

    Raises:
        UnplottableGeometryError: The primitive carries no path to sample.
    """
    if isinstance(geometry, Circle):
        count = _arc_samples(2.0 * math.pi)
        return tuple(_xy(geometry.point_at(step / count)) for step in range(count + 1))
    if isinstance(geometry, Arc):
        count = _arc_samples(geometry.angle)
        return tuple(_xy(geometry.point_at(step / count)) for step in range(count + 1))
    if isinstance(geometry, Line):
        return (_xy(geometry.start), _xy(geometry.end))
    raise UnplottableGeometryError(f"Operation {index} carries geometry of type {type(geometry).__name__!r}, which has no drawable path.")


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
