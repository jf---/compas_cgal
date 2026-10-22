"""Exact containment of emitted circles in the polygon cutter-centre domain, with a reported protrusion."""

from __future__ import annotations

import numpy as np

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from compas_cgal import _coverage_2
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


def containment_report(circles: tuple[Figure5CounterclockwiseCircle, ...], boundary: tuple[Point2[WorldXY], ...]) -> tuple[int, Millimetre]:
    """Count circles whose centre the native cutter-centre domain rejects, and the largest reported protrusion.

    The rejection is the exact native predicate's; the protrusion (radius minus
    the double clearance to the boundary polygon) is reporting only.
    """
    polygon = np.array([(float(point.x), float(point.y), 0.0) for point in boundary])
    starts = polygon[:, :2]
    edges = np.roll(starts, -1, axis=0) - starts
    rejected = 0
    maximum_excess = 0.0
    for circle in circles:
        radius = float(circle.radius.value)
        center = np.array((float(circle.center.x), float(circle.center.y)))
        if _coverage_2.CutterCentreDomain2.build(polygon, [], radius).contains(float(center[0]), float(center[1])):
            continue
        rejected += 1
        # Reporting only: the exact native predicate above owns the rejection.
        parameters = np.clip(np.sum((center - starts) * edges, axis=1) / np.sum(edges * edges, axis=1), 0.0, 1.0)
        clearance = float(np.linalg.norm(center - starts - parameters[:, None] * edges, axis=1).min())
        maximum_excess = max(maximum_excess, radius - clearance)
    return rejected, Millimetre(maximum_excess)
