"""Exact gouge audit of Held paths on the native declared-boundary oracle: one batch per element kind.

The report is :func:`benchmarks.held_omega_gouge_audit.audit_omega_path`'s and
:func:`~benchmarks.held_omega_gouge_audit.audit_polygon_path`'s, element for
element: the same verdicts, the same reported primitive, the same worst
sub-piece selection. Only where the decisions run differs: every circle,
every whole ∂Ω* piece and every polyline piece is decided in one native call
per kind (:class:`benchmarks.held_declared_clearance_native.NativeDeclaredClearance`).

∂Ω* pieces carry their declared owners (LINE / ARC features) into the native
call; a sub-piece of a gouging piece is decided on its own, in a second batch,
as the Python audit does.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import Element
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_declared_clearance_wire import offset_owners
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import ElementClearance
from benchmarks.held_gouge_audit import ElementKind
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connectors import OmegaConnector
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def _listed(clearance: NativeDeclaredClearance, kind: ElementKind, index: int, verdict: ClearanceVerdict, segment: int | None) -> ElementClearance | None:
    if verdict.clearance is Clearance.CLEAR or verdict.primitive is None:
        return None
    return ElementClearance(kind, index, verdict.clearance, verdict.depth_mm, verdict.primitive, clearance.primitive_kind(verdict.primitive), segment)


def _rank(verdict: ClearanceVerdict) -> tuple[Clearance, float]:
    return verdict.clearance, -verdict.depth_mm


def _piece_element(piece: CentreDomainPiece) -> Element:
    support = piece.support
    if isinstance(support, OffsetLineSupport):
        return SegmentElement(piece.source, piece.target, support.direction, support.normal, support.level)
    return ArcElement(support.centre, support.radius, support.sense, piece.source, piece.target)


def _piece_verdicts(clearance: NativeDeclaredClearance, connectors: Sequence[OmegaConnector]) -> dict[int, ClearanceVerdict]:
    """Each whole ∂Ω* piece any connector runs on, decided once (first-seen order), with its owners."""
    pieces: dict[int, CentreDomainPiece] = {}
    for connector in connectors:
        for sub in connector.pieces:
            pieces.setdefault(sub.index, sub.piece)
    order = list(pieces)
    verdicts = clearance.elements([_piece_element(pieces[index]) for index in order], [offset_owners(pieces[index].owners) for index in order])
    return dict(zip(order, verdicts, strict=True))


def _sub_piece_verdicts(clearance: NativeDeclaredClearance, connectors: Sequence[OmegaConnector], pieces: dict[int, ClearanceVerdict]) -> dict[tuple[int, int], ClearanceVerdict]:
    """Each sub-piece of a gouging piece, decided on its own (keyed by connector and sub-piece position)."""
    keys = [(index, segment) for index, connector in enumerate(connectors) for segment, sub in enumerate(connector.pieces) if pieces[sub.index].clearance is Clearance.GOUGE]
    subs = [connectors[index].pieces[segment] for index, segment in keys]
    verdicts = clearance.elements([sub.element() for sub in subs], [offset_owners(sub.piece.owners) for sub in subs])
    return dict(zip(keys, verdicts, strict=True))


def _connector_verdict(
    clearance: NativeDeclaredClearance, index: int, connector: OmegaConnector, pieces: dict[int, ClearanceVerdict], subs: dict[tuple[int, int], ClearanceVerdict]
) -> ElementClearance | None:
    """Worst verdict over the connector's sub-pieces: the whole piece's, or the sub-piece's own where the piece gouges."""
    worst: tuple[ClearanceVerdict, int] | None = None
    for segment, sub in enumerate(connector.pieces):
        verdict = pieces[sub.index]
        if verdict.clearance is Clearance.GOUGE:
            verdict = subs[(index, segment)]
        if verdict.clearance is not Clearance.CLEAR and (worst is None or _rank(verdict) < _rank(worst[0])):
            worst = (verdict, segment)
    return None if worst is None else _listed(clearance, "connector", index, worst[0], worst[1])


def audit_omega_path_native(clearance: NativeDeclaredClearance, circles: Sequence[OmegaCircle], connectors: Sequence[OmegaConnector]) -> HeldGougeAuditReport:
    """Decide every circle and every connector of an Ω* path exactly, natively."""
    started = time.perf_counter()
    circle_verdicts = clearance.circles([circle.element() for circle in circles])
    found = [_listed(clearance, "circle", index, verdict, None) for index, verdict in enumerate(circle_verdicts)]
    pieces = _piece_verdicts(clearance, connectors)
    subs = _sub_piece_verdicts(clearance, connectors, pieces)
    found += [_connector_verdict(clearance, index, connector, pieces, subs) for index, connector in enumerate(connectors)]
    elements = tuple(item for item in found if item is not None)
    return HeldGougeAuditReport(float(clearance.tool_radius), clearance.primitive_count, len(circles), len(connectors), elements, time.perf_counter() - started)


def _exact(point: Point2[WorldXY]) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].from_binary64(float(point.x), float(point.y))


def _polyline_pieces(samples: Sequence[Point2[WorldXY]]) -> list[tuple[ExactPoint2[WorldXY], ExactPoint2[WorldXY]]]:
    points = [_exact(sample) for sample in samples]
    return [(start, end) for start, end in zip(points, points[1:]) if start != end]


def audit_polygon_path_native(
    clearance: NativeDeclaredClearance, circles: Sequence[tuple[Point2[WorldXY], float]], polylines: Sequence[Sequence[Point2[WorldXY]]]
) -> HeldGougeAuditReport:
    """Decide every binary64 circle ``(centre, radius)`` and connector polyline exactly against the declared boundary, natively."""
    started = time.perf_counter()
    circle_verdicts = clearance.circles([CircleElement(_exact(centre), Fraction(radius)) for centre, radius in circles])
    found = [_listed(clearance, "circle", index, verdict, None) for index, verdict in enumerate(circle_verdicts)]
    pieces = [_polyline_pieces(samples) for samples in polylines]
    verdicts = clearance.segments_binary64([piece for polyline in pieces for piece in polyline])
    offset = 0
    for index, polyline in enumerate(pieces):
        own = verdicts[offset : offset + len(polyline)]
        offset += len(polyline)
        decided = [(verdict, segment) for segment, verdict in enumerate(own) if verdict.clearance is not Clearance.CLEAR]
        if decided:
            verdict, segment = min(decided, key=lambda item: _rank(item[0]))
            found.append(_listed(clearance, "connector", index, verdict, segment))
    elements = tuple(item for item in found if item is not None)
    return HeldGougeAuditReport(float(clearance.tool_radius), clearance.primitive_count, len(circles), len(polylines), elements, time.perf_counter() - started)
