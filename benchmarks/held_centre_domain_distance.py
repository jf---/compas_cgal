"""Exact nearest point and squared distance from a rational point to one ∂Ω* piece (Python oracle).

The production query is native (``src/held_centre_domain_query_2.cpp``); this
pure-``Fraction`` formulation is its independent differential oracle in
``tests/benchmarks/test_held_centre_domain_queries.py``.

Every quantity is a one-root number: the distance to an offset line is
rational, the distance to a circle piece in its wedge is ``(√D − R)²`` with
``D = |p − c|²`` (one-root in ``√D``), and the nearest point on it is
``c + R(p − c)/√D`` (same extension). Comparisons across pieces are exact
cross-root comparisons (:meth:`OneRootNumber.compare`); no float decides
anything. Design: docs/held_exact_centre_domain.md, "Projection of a hypothesis
contact".
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from benchmarks.held_centre_domain_supports import CentreDomainSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from compas_cgal.units import WorldXY


@dataclass(frozen=True, eq=False)
class NearestOnPiece:
    """The nearest point of a piece to a query point and its exact squared distance.

    ``ambiguous`` is set when two distinct points of the piece are nearest
    (the query is the centre of a circle piece, or equidistant from both ends
    of a piece it lies outside the wedge of).
    """

    squared_distance: OneRootNumber
    point: OneRootPoint2
    ambiguous: bool


def _squared_distance(point: OneRootPoint2, query: ExactPoint2[WorldXY]) -> OneRootNumber:
    dx, dy = point.minus(query)
    return dx * dx + dy * dy


def _nearest_end(source: OneRootPoint2, target: OneRootPoint2, query: ExactPoint2[WorldXY]) -> NearestOnPiece:
    to_source, to_target = _squared_distance(source, query), _squared_distance(target, query)
    order = to_source.compare(to_target)
    if order < 0:
        return NearestOnPiece(to_source, source, False)
    if order > 0:
        return NearestOnPiece(to_target, target, False)
    return NearestOnPiece(to_source, source, source != target)


def nearest_on_piece(support: CentreDomainSupport, source: OneRootPoint2, target: OneRootPoint2, query: ExactPoint2[WorldXY]) -> NearestOnPiece:
    """Exact nearest point of the piece ``source → target`` on ``support`` to ``query``.

    A piece with ``source == target`` is a full circle.
    """
    if isinstance(support, OffsetLineSupport):
        along = support.direction[0] * (query.x - support.anchor.x) + support.direction[1] * (query.y - support.anchor.y)
        if support.parameter(source) <= along <= support.parameter(target):
            foot = OneRootPoint2.rational(support.anchor.x + along * support.direction[0], support.anchor.y + along * support.direction[1])
            height = support.normal[0] * query.x + support.normal[1] * query.y - support.level
            return NearestOnPiece(OneRootNumber.of(height * height), foot, False)
        return _nearest_end(source, target, query)
    dx, dy = query.minus(support.centre)
    radial_squared = Fraction(dx * dx + dy * dy)
    if radial_squared == 0:
        return NearestOnPiece(OneRootNumber.of(support.radius_squared), source, True)
    probe = OneRootPoint2.from_exact(query)
    full_circle = source == target
    if full_circle or compare_on_arc(support.centre, support.sense, source, probe, target) <= 0:
        radius = support.radius
        distance = OneRootNumber.build(radial_squared + support.radius_squared, -2 * radius, radial_squared)
        scale = radius / radial_squared
        point = OneRootPoint2.build(
            OneRootNumber.build(support.centre.x, scale * dx, radial_squared),
            OneRootNumber.build(support.centre.y, scale * dy, radial_squared),
        )
        return NearestOnPiece(distance, point, False)
    return _nearest_end(source, target, query)
