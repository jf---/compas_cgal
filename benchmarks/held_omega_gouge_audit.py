"""Exact gouge audit of Held paths against the declared boundary: Ω* paths and polygon-Ω paths alike.

One oracle (:class:`benchmarks.held_declared_clearance.DeclaredClearance`)
decides every element of both paths, so the old-versus-new table compares like
with like:

- an Ω* path: binary64 circles and exact ∂Ω* connectors (segment and arc
  sub-pieces with one-root endpoints);
- a polygon-Ω path: binary64 circles and binary64 connector polylines.

A connector's verdict is its worst sub-piece; ``segment`` names that sub-piece.
Fan connectors (zero progress at one exact site) have no motion to audit. The
report type is the stored-boundary audit's
(:class:`benchmarks.held_gouge_audit.HeldGougeAuditReport`).

``backend`` selects where the decisions run: :attr:`ClearanceBackend.PYTHON`
(this module, the reference oracle) or :attr:`ClearanceBackend.NATIVE`
(:mod:`benchmarks.held_omega_gouge_audit_native`, the C++ oracle built from the
same declared boundary; its build time is part of the reported seconds). Both
produce the same report decisions.
"""

from __future__ import annotations

import dataclasses
import enum
import time
from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance import Element
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import ElementClearance
from benchmarks.held_gouge_audit import ElementKind
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connectors import OmegaConnector
from benchmarks.held_omega_gouge_audit_native import audit_omega_path_native
from benchmarks.held_omega_gouge_audit_native import audit_polygon_path_native
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


class ClearanceBackend(enum.Enum):
    """Where the exact clearance decisions run."""

    PYTHON = "python"
    NATIVE = "native"


def _native(clearance: DeclaredClearance) -> NativeDeclaredClearance:
    return NativeDeclaredClearance.build(clearance.boundary, clearance.tool_radius)


def _listed(clearance: DeclaredClearance, kind: ElementKind, index: int, verdict: ClearanceVerdict, segment: int | None) -> ElementClearance | None:
    if verdict.clearance is Clearance.CLEAR or verdict.primitive is None:
        return None
    return ElementClearance(kind, index, verdict.clearance, verdict.depth_mm, verdict.primitive, clearance.primitive_kind(verdict.primitive), segment)


def _worst(clearance: DeclaredClearance, index: int, elements: Sequence[Element]) -> ElementClearance | None:
    found = [(clearance.element(element), segment) for segment, element in enumerate(elements)]
    decided = [(verdict, segment) for verdict, segment in found if verdict.clearance is not Clearance.CLEAR]
    if not decided:
        return None
    verdict, segment = min(decided, key=lambda item: (item[0].clearance, -item[0].depth_mm))
    return _listed(clearance, "connector", index, verdict, segment)


def _piece_element(piece: CentreDomainPiece) -> Element:
    support = piece.support
    if isinstance(support, OffsetLineSupport):
        return SegmentElement(piece.source, piece.target, support.direction, support.normal, support.level)
    return ArcElement(support.centre, support.radius, support.sense, piece.source, piece.target)


def _connector_verdict(clearance: DeclaredClearance, index: int, connector: OmegaConnector, pieces: dict[int, ClearanceVerdict]) -> ElementClearance | None:
    """Worst verdict over the whole ∂Ω* pieces the connector runs on, sub-pieces of a gouging piece decided on their own.

    A sub-piece is a subset of its piece, so its clearance is at least the
    piece's: a piece that does not gouge clears every sub-piece of it. A
    gouging piece's sub-pieces are decided individually (the connector may
    run on the clear part).
    """
    worst: tuple[ClearanceVerdict, int] | None = None
    for segment, sub in enumerate(connector.pieces):
        if sub.index not in pieces:
            pieces[sub.index] = clearance.element(_piece_element(sub.piece))
        verdict = pieces[sub.index]
        if verdict.clearance is Clearance.GOUGE:
            verdict = clearance.element(sub.element())
        if verdict.clearance is not Clearance.CLEAR and (worst is None or (verdict.clearance, -verdict.depth_mm) < (worst[0].clearance, -worst[0].depth_mm)):
            worst = (verdict, segment)
    return None if worst is None else _listed(clearance, "connector", index, worst[0], worst[1])


def audit_omega_path(
    clearance: DeclaredClearance, circles: Sequence[OmegaCircle], connectors: Sequence[OmegaConnector], backend: ClearanceBackend = ClearanceBackend.PYTHON
) -> HeldGougeAuditReport:
    """Decide every circle and every connector of an Ω* path exactly.

    A connector's TOUCH means it runs on ∂Ω* pieces whose exact clearance is r
    (each touches its owner along its whole length by construction); its
    sub-pieces have clearance >= r.
    """
    started = time.perf_counter()
    if backend is ClearanceBackend.NATIVE:
        report = audit_omega_path_native(_native(clearance), circles, connectors)
        return dataclasses.replace(report, seconds=time.perf_counter() - started)
    found = [_listed(clearance, "circle", index, clearance.circle(circle.element()), None) for index, circle in enumerate(circles)]
    pieces: dict[int, ClearanceVerdict] = {}
    found += [_connector_verdict(clearance, index, connector, pieces) for index, connector in enumerate(connectors)]
    elements = tuple(item for item in found if item is not None)
    return HeldGougeAuditReport(float(clearance.tool_radius), len(clearance.primitives), len(circles), len(connectors), elements, time.perf_counter() - started)


def _exact(point: Point2[WorldXY]) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].from_binary64(float(point.x), float(point.y))


def _polyline_elements(samples: Sequence[Point2[WorldXY]]) -> list[Element]:
    points = [_exact(sample) for sample in samples]
    return [SegmentElement.between(start, end) for start, end in zip(points, points[1:]) if start != end]


def audit_polygon_path(
    clearance: DeclaredClearance,
    circles: Sequence[tuple[Point2[WorldXY], float]],
    polylines: Sequence[Sequence[Point2[WorldXY]]],
    backend: ClearanceBackend = ClearanceBackend.PYTHON,
) -> HeldGougeAuditReport:
    """Decide every binary64 circle ``(centre, radius)`` and connector polyline exactly against the declared boundary."""
    started = time.perf_counter()
    if backend is ClearanceBackend.NATIVE:
        report = audit_polygon_path_native(_native(clearance), circles, polylines)
        return dataclasses.replace(report, seconds=time.perf_counter() - started)
    found = [_listed(clearance, "circle", index, clearance.circle(CircleElement(_exact(centre), Fraction(radius))), None) for index, (centre, radius) in enumerate(circles)]
    found += [_worst(clearance, index, _polyline_elements(samples)) for index, samples in enumerate(polylines)]
    elements = tuple(item for item in found if item is not None)
    return HeldGougeAuditReport(float(clearance.tool_radius), len(clearance.primitives), len(circles), len(polylines), elements, time.perf_counter() - started)
