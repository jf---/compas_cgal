"""Swept regions: full-circle cuts and boundary connectors."""

from __future__ import annotations

from dataclasses import dataclass

from compas_cgal.coverage_certificate.errors import DegenerateConnectorError
from compas_cgal.coverage_certificate.errors import InvalidCentreBoundarySiteError
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


@dataclass(frozen=True)
class FullCircleSweep:
    """A full-circle cut: sweeps the annulus max(R - r, 0) <= |x - c| <= R + r (a disk when R <= r)."""

    centre: Point2[WorldXY]
    path_radius: GuideRadius

    @classmethod
    def build(cls, centre: Point2[WorldXY], path_radius: GuideRadius) -> FullCircleSweep:
        """Build a cut; both values are already validated finite and positive by their unit types."""
        return cls(centre, path_radius)


@dataclass(frozen=True)
class BoundaryConnectorSweep:
    """A CCW interval of one centre-domain loop between two exact sites; sweeps polyline (+) D_r.

    The certificate enters its capsules as end disks (CAPSULE_CAP) and rectangles
    (CAPSULE_FLANK), after the exact union of collinear pieces on each edge.
    """

    start: CentreBoundarySite
    end: CentreBoundarySite

    @classmethod
    def build(cls, domain: CentreDomain, start: CentreBoundarySite, end: CentreBoundarySite) -> BoundaryConnectorSweep:
        """Build a connector CCW along one loop from ``start`` to ``end``.

        Raises:
            InvalidCentreBoundarySiteError: a site names a loop or segment the
                domain lacks, or the sites lie on different loops.
            DegenerateConnectorError: zero CCW progress (start == end).
        """
        start.check_in(domain)
        end.check_in(domain)
        if start.loop != end.loop:
            raise InvalidCentreBoundarySiteError(f"a connector follows one loop, got loops {start.loop} and {end.loop}")
        if (start.segment, start.parameter) == (end.segment, end.parameter):
            raise DegenerateConnectorError(f"connector from {start} to itself makes no CCW progress")
        return cls(start, end)
