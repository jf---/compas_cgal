"""Offset contact projection: one exact AABB clearance batch per boundary."""

from __future__ import annotations

from collections.abc import Iterator
from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from compas_cgal import _circle_geometry_2
from compas_cgal.held_boundary_clearance import HeldBoundaryClearanceIndex
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Gap zero: the admissible set of a query is exactly its set of tied nearest
# sites (sqrt(d_i) - sqrt(min d) <= 0 iff d_i == min d), decided exactly.
_NEAREST_ONLY = Millimetre(0.0)
# The polygon contract of ``_circle_geometry_2.project_boundary_contact``.
_MINIMUM_BOUNDARY_VERTICES = 3


def _nearest_sites(
    contacts: Sequence[Point2[WorldXY]],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> Iterator[ProjectionBoundarySite | None]:
    """Resolve the contacts in order, lazily, so callers stop at the first failure as the native loop did."""
    bound = float(evidence_bound)
    if len(boundary) < _MINIMUM_BOUNDARY_VERTICES or not 0.0 <= bound < float("inf"):
        raise _circle_geometry_2.InvalidBoundaryProjectionError("Projection needs a polygon and nonnegative finite distance bound.")
    squared_bound = Fraction(bound) ** 2
    rows = HeldBoundaryClearanceIndex.build(boundary).distance_admissible_sites(contacts, _NEAREST_ONLY)
    for row in rows:
        nearest = row[0]
        if nearest.squared_distance > squared_bound:
            yield None
            continue
        # Distinct site identities are distinct points on a simple ring; the
        # footpoint test keeps the native rule (distinct points) regardless.
        if any(site.footpoint != nearest.footpoint for site in row[1:]):
            raise _circle_geometry_2.AmbiguousBoundaryProjectionError("Distinct offset contacts have equal minimum distance.")
        yield ProjectionBoundarySite.build(ProjectionBoundarySegmentId(nearest.segment), Fraction(float(nearest.parameter)), len(boundary))


def project_offset_contacts(
    contacts: Sequence[Point2[WorldXY]],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> tuple[ProjectionBoundarySite | None, ...]:
    """Nearest offset-boundary site of every contact, ``None`` beyond the bound.

    Same decisions as ``_circle_geometry_2.project_boundary_contact`` per
    contact, all exact: the nearest site, the tie (distinct nearest points at
    equal distance) and the bound (minimum squared distance greater than
    ``evidence_bound²`` is rejected, checked before the tie). The queries run
    as one ``HeldBoundaryClearanceIndex`` batch instead of an O(n) loop each.

    The site parameter is the correctly rounded double of the exact rational
    projection parameter, injected back exactly. A parameter that rounds to
    one canonicalises to the next vertex, as ``ProjectionBoundarySite.build``
    does for any parameter one.

    Raises:
        _circle_geometry_2.InvalidBoundaryProjectionError: fewer than three
            boundary vertices, or a negative or non-finite bound.
        _circle_geometry_2.AmbiguousBoundaryProjectionError: a contact within
            the bound has distinct nearest points at equal distance.
    """
    return tuple(_nearest_sites(contacts, boundary, evidence_bound))


def project_bounded_offset_contacts(
    contacts: Sequence[Point2[WorldXY]],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> tuple[ProjectionBoundarySite, ...]:
    """``project_offset_contacts`` where a contact beyond the bound fails.

    The first failing contact in order decides the error, as in a loop over
    ``_circle_geometry_2.project_boundary_contact``.

    Raises:
        _circle_geometry_2.BoundaryProjectionDistanceError: a contact's
            nearest site lies beyond the bound.
        _circle_geometry_2.AmbiguousBoundaryProjectionError: see
            ``project_offset_contacts``.
        _circle_geometry_2.InvalidBoundaryProjectionError: see
            ``project_offset_contacts``.
    """
    sites = []
    for site in _nearest_sites(contacts, boundary, evidence_bound):
        if site is None:
            raise _circle_geometry_2.BoundaryProjectionDistanceError("Closest offset contact exceeds the boundary evidence bound.")
        sites.append(site)
    return tuple(sites)
