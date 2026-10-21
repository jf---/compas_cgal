"""Exact gouge audit of a Held toolpath against its analytic pocket boundary.

A path element (a machining circle, or one straight piece of a connector) is
the locus of the tool centre. The swept tool disk of radius ``r`` leaves the
pocket exactly when ``dist(element, ∂P) < r``; ``dist == r`` is a tie (the
tool touches the wall, allowed) and ``dist > r`` is clear. Every such decision
here is exact: inputs are injected as rationals, distances are the extremes of
closed-form candidate sets (:mod:`benchmarks.held_gouge_boundary`), and each
comparison with ``r`` is the sign of one quadratic surd over the integers of
one dyadic grid that holds every input double exactly. There is no tolerance
and no escalation floor. Depths are REPORTED as doubles and never decide.

Precondition (stated, not checked): every element lies in the closed pocket.
The path is connected, so an element that does not reach within ``r`` of the
boundary is on the same side as every other element.

Candidate pruning is exact too: a primitive is skipped only when the
double-evaluated gap between the element's and the primitive's outward-rounded
boxes is ``> r``. IEEE subtraction is correctly rounded and therefore monotone,
so ``fl(x - y) > r`` (``r`` a double) implies ``x - y > r``: the skipped
primitive is provably farther than ``r``.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass
from decimal import Decimal
from decimal import localcontext
from enum import IntEnum
from fractions import Fraction
from typing import Literal

import numpy as np
import numpy.typing as npt

from benchmarks.errors import BenchmarkError
from benchmarks.held_gouge_boundary import Box
from benchmarks.held_gouge_boundary import DyadicGrid
from benchmarks.held_gouge_boundary import ExactPrimitive
from benchmarks.held_gouge_boundary import PrimitiveKind
from benchmarks.held_gouge_boundary import ReferencePrimitive
from benchmarks.held_gouge_boundary import Vec
from benchmarks.held_gouge_boundary import exact_box
from benchmarks.held_gouge_boundary import exact_primitive
from benchmarks.held_gouge_boundary import primitive_exponent
from benchmarks.held_gouge_boundary import round_down
from benchmarks.held_gouge_boundary import round_up
from benchmarks.held_gouge_surd import REPORT_DIGITS
from benchmarks.held_gouge_surd import Surd
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Reporting bucket only, mm: splits gouge counts the way the 2026-10-17 Decimal
# audit did (geometric sagitta gouges vs rounding-scale ones). Never a decision.
REPORTING_SPLIT_MM = 1e-9
ElementKind = Literal["circle", "connector"]


class HeldGougeAuditError(BenchmarkError):
    """Base error of the Held gouge audit."""


class EmptyBoundaryError(HeldGougeAuditError):
    """The pocket boundary has no primitives."""


class NonPositiveToolRadiusError(HeldGougeAuditError):
    """The tool radius is not positive."""


class HeldPathGougesBoundaryError(HeldGougeAuditError):
    """At least one path element sweeps the tool disk outside the pocket."""


class Clearance(IntEnum):
    """Exact sign of ``dist(element, ∂P) - r``."""

    GOUGE = -1
    TOUCH = 0
    CLEAR = 1


@dataclass(frozen=True)
class AuditCircle:
    """A machining circle: the tool centre runs on it."""

    centre: Point2[WorldXY]
    radius: Millimetre


@dataclass(frozen=True)
class AuditConnector:
    """A connector polyline: the tool centre runs along its straight pieces."""

    samples: tuple[Point2[WorldXY], ...]


@dataclass(frozen=True)
class ElementClearance:
    """The exact verdict on one non-clear path element.

    Attributes:
        element: ``"circle"`` or ``"connector"``.
        index: Index of the element in its sequence.
        clearance: Exact GOUGE or TOUCH.
        depth_mm: Reported ``r - dist`` (double); ``0.0`` exactly for TOUCH.
        primitive: Boundary primitive of the deepest (or touching) pair.
        primitive_kind: ``"line"``, ``"convex-arc"`` or ``"concave-arc"``.
        segment: For a connector, the straight piece of the deepest pair.
    """

    element: ElementKind
    index: int
    clearance: Clearance
    depth_mm: float
    primitive: int
    primitive_kind: PrimitiveKind
    segment: int | None


@dataclass(frozen=True)
class HeldGougeAuditReport:
    """Exact gouge verdicts of one path; only non-clear elements are listed."""

    tool_radius_mm: float
    primitives: int
    circles: int
    connectors: int
    elements: tuple[ElementClearance, ...]
    seconds: float

    def _of(self, element: ElementKind, clearance: Clearance) -> tuple[ElementClearance, ...]:
        return tuple(item for item in self.elements if item.element == element and item.clearance == clearance)

    @property
    def gouging_circles(self) -> tuple[ElementClearance, ...]:
        return self._of("circle", Clearance.GOUGE)

    @property
    def gouging_connectors(self) -> tuple[ElementClearance, ...]:
        return self._of("connector", Clearance.GOUGE)

    @property
    def touching_circles(self) -> tuple[ElementClearance, ...]:
        return self._of("circle", Clearance.TOUCH)

    @property
    def touching_connectors(self) -> tuple[ElementClearance, ...]:
        return self._of("connector", Clearance.TOUCH)

    @property
    def gouges(self) -> tuple[ElementClearance, ...]:
        return tuple(item for item in self.elements if item.clearance == Clearance.GOUGE)

    def max_depth_mm(self, element: ElementKind) -> float:
        """Reported deepest gouge of one element kind; ``0.0`` when none gouges."""
        return max((item.depth_mm for item in self._of(element, Clearance.GOUGE)), default=0.0)

    def above_split(self, element: ElementKind) -> int:
        """Gouges of one kind deeper than the reporting split (comparison with the Decimal audit)."""
        return sum(1 for item in self._of(element, Clearance.GOUGE) if item.depth_mm > REPORTING_SPLIT_MM)


def _decimal_root(square: Surd, scale: Decimal) -> Decimal:
    """Reporting-only √ of an exact non-negative surd in grid units², returned in mm."""
    with localcontext() as context:
        context.prec = REPORT_DIGITS
        value = square.to_decimal()
        return value.sqrt() / scale if value > 0 else Decimal(0)


@dataclass(frozen=True)
class _Pair:
    """A non-clear (element piece, primitive) pair with its reported distance in mm."""

    clearance: int
    distance: Decimal
    primitive: int


class _GridAudit:
    """One audit run: boundary, radius and path injected onto one dyadic grid."""

    def __init__(self, primitives: tuple[ExactPrimitive, ...], boxes: npt.NDArray[np.float64], tool_radius: float, grid: DyadicGrid) -> None:
        self._primitives = primitives
        self._boxes = boxes
        self._radius_double = tool_radius
        self._grid = grid
        self._radius = grid.integer(tool_radius)
        self._radius_mm = Decimal(tool_radius)  # exact: Decimal(float) converts without rounding
        with localcontext() as context:
            context.prec = REPORT_DIGITS
            self._scale = Decimal(2) ** grid.exponent

    def _near(self, box: Box) -> tuple[ExactPrimitive, ...]:
        """Primitives not provably farther than ``r`` from the box (see the module docstring)."""
        b = self._boxes
        gap = np.maximum.reduce([b[:, 0] - box.xmax, box.xmin - b[:, 2], b[:, 1] - box.ymax, box.ymin - b[:, 3]])
        return tuple(self._primitives[int(i)] for i in np.nonzero(~(gap > self._radius_double))[0])

    def _circle_pair(self, primitive: ExactPrimitive, c: Vec, radius: int) -> _Pair | None:
        # dist = max(lo - R, R - hi, 0); compare(lo - R, r) = cmp_lo, compare(R - hi, r) = -cmp_hi.
        lows = primitive.nearest(c)
        cmp_lo = min(square.compare_square(radius + self._radius) for square in lows)
        if cmp_lo > 0:
            return None
        inner = radius - self._radius
        highs = primitive.farthest(c)
        cmp_hi = 1 if inner < 0 else max(square.compare_square(inner) for square in highs)
        verdict = max(cmp_lo, -cmp_hi, -1)
        if verdict > 0:
            return None
        with localcontext() as context:
            context.prec = REPORT_DIGITS
            lo = min(_decimal_root(square, self._scale) for square in lows)
            hi = max(_decimal_root(square, self._scale) for square in highs)
            big_r = Decimal(radius) / self._scale
            return _Pair(verdict, max(lo - big_r, big_r - hi, Decimal(0)), primitive.index)

    def _segment_pair(self, primitive: ExactPrimitive, p: Vec, q: Vec) -> _Pair | None:
        squares = primitive.segment_candidates(p, q)
        if squares is None:
            return _Pair(Clearance.GOUGE, Decimal(0), primitive.index)
        verdict = min(square.compare_square(self._radius) for square in squares)
        if verdict > 0:
            return None
        return _Pair(verdict, min(_decimal_root(square, self._scale) for square in squares), primitive.index)

    def _verdict(self, element: ElementKind, index: int, pairs: Sequence[tuple[_Pair | None, int | None]]) -> ElementClearance | None:
        decided = [(pair, segment) for pair, segment in pairs if pair is not None]
        if not decided:
            return None
        worst = min(pair.clearance for pair, _ in decided)
        pair, segment = min(((p, s) for p, s in decided if p.clearance == worst), key=lambda item: item[0].distance)
        primitive = self._primitives[pair.primitive]
        depth = 0.0 if worst == Clearance.TOUCH else float(self._radius_mm - pair.distance)
        return ElementClearance(element, index, Clearance(worst), depth, pair.primitive, primitive.kind, segment)

    def circle(self, index: int, circle: AuditCircle) -> ElementClearance | None:
        c = self._grid.point(circle.centre.x, circle.centre.y)
        radius = self._grid.integer(circle.radius)
        cx, cy, big_r = Fraction(circle.centre.x), Fraction(circle.centre.y), Fraction(circle.radius)
        box = Box(round_down(cx - big_r), round_down(cy - big_r), round_up(cx + big_r), round_up(cy + big_r))
        return self._verdict("circle", index, [(self._circle_pair(p, c, radius), None) for p in self._near(box)])

    def connector(self, index: int, connector: AuditConnector) -> ElementClearance | None:
        pairs: list[tuple[_Pair | None, int | None]] = []
        samples = connector.samples
        for segment, (start, end) in enumerate(zip(samples, samples[1:])):
            p, q = self._grid.point(start.x, start.y), self._grid.point(end.x, end.y)
            box = Box(min(start.x, end.x), min(start.y, end.y), max(start.x, end.x), max(start.y, end.y))
            pairs.extend((self._segment_pair(prim, p, q), segment) for prim in self._near(box))
        if len(samples) == 1:
            p = self._grid.point(samples[0].x, samples[0].y)
            box = Box(samples[0].x, samples[0].y, samples[0].x, samples[0].y)
            pairs.extend((self._segment_pair(prim, p, p), 0) for prim in self._near(box))
        return self._verdict("connector", index, pairs)


class HeldGougeAudit:
    """Exact clearance oracle over one analytic pocket boundary and tool radius."""

    def __init__(self, primitives: tuple[ReferencePrimitive, ...], tool_radius: float) -> None:
        self._primitives = primitives
        self._tool_radius = tool_radius
        self._exponent = max([DyadicGrid.exponent_of(tool_radius), *map(primitive_exponent, primitives)])
        grid = DyadicGrid(self._exponent)
        boxes = [exact_box(exact_primitive(i, p, grid), grid) for i, p in enumerate(primitives)]
        self._boxes: npt.NDArray[np.float64] = np.array([(b.xmin, b.ymin, b.xmax, b.ymax) for b in boxes], dtype=np.float64)

    @classmethod
    def build(cls, primitives: Sequence[ReferencePrimitive], tool_radius: Millimetre) -> HeldGougeAudit:
        """Validate and index the analytic boundary and the tool radius.

        Raises:
            EmptyBoundaryError: no primitives.
            NonPositiveToolRadiusError: ``tool_radius <= 0``.
        """
        if not primitives:
            raise EmptyBoundaryError("a pocket boundary needs at least one primitive")
        if not float(tool_radius) > 0.0:
            raise NonPositiveToolRadiusError(f"tool radius {tool_radius} is not positive")
        return cls(tuple(primitives), float(tool_radius))

    @property
    def primitive_count(self) -> int:
        return len(self._primitives)

    def audit(self, circles: Sequence[AuditCircle], connectors: Sequence[AuditConnector]) -> HeldGougeAuditReport:
        """Decide every element exactly; list the gouging and touching ones."""
        started = time.perf_counter()
        values = [v for c in circles for v in (c.centre.x, c.centre.y, c.radius)] + [v for t in connectors for s in t.samples for v in (s.x, s.y)]
        grid = DyadicGrid(max([self._exponent, *map(DyadicGrid.exponent_of, values)]))
        run = _GridAudit(tuple(exact_primitive(i, p, grid) for i, p in enumerate(self._primitives)), self._boxes, self._tool_radius, grid)
        found = [run.circle(i, c) for i, c in enumerate(circles)]
        found += [run.connector(i, c) for i, c in enumerate(connectors)]
        elements = tuple(item for item in found if item is not None)
        return HeldGougeAuditReport(self._tool_radius, len(self._primitives), len(circles), len(connectors), elements, time.perf_counter() - started)


def require_gouge_free(report: HeldGougeAuditReport) -> None:
    """Enforce: no element may sweep the tool outside the pocket (touching is allowed).

    Raises:
        HeldPathGougesBoundaryError: listing every gouging element with its reported depth.
    """
    gouges = report.gouges
    if not gouges:
        return
    lines = [
        f"{item.element} {item.index}"
        + (f" piece {item.segment}" if item.segment is not None else "")
        + f": depth {item.depth_mm:.6g} mm against {item.primitive_kind} {item.primitive}"
        for item in sorted(gouges, key=lambda item: -item.depth_mm)
    ]
    raise HeldPathGougesBoundaryError(
        f"{len(report.gouging_circles)} circles and {len(report.gouging_connectors)} connectors gouge the pocket boundary (r = {report.tool_radius_mm} mm):\n" + "\n".join(lines)
    )
