"""Native stock state from actual full-circle sweeps in a Held draft prefix.

Circle-only replay: entry clearing and connector sweeps are deliberately absent
from this state. It is not the paper's filled-disk contour approximation or an
engagement certificate. Python never computes the annulus radii.
"""

from __future__ import annotations

import numpy as np

from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from compas_cgal import _stock_2
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def replay_circle_stock(
    boundary: tuple[Point2[WorldXY], ...],
    circles: tuple[Figure5CounterclockwiseCircle, ...],
    tool_radius: ToolRadius,
) -> _stock_2.Stock2:
    """Subtract the supplied full circles from a fresh polygon pocket.

    Returns remaining material after exactly the supplied prefix. Each stored
    centre, guide radius and tool radius enters CGAL separately as binary64;
    the exact kernel constructs and unions their swept annuli. This preserves
    uncut central islands and all historical cuts, without inventing a plunge.
    """
    polygon = np.array([(float(point.x), float(point.y), 0.0) for point in boundary], dtype=np.float64)
    stock = _stock_2.Stock2(polygon, [])
    for circle in circles:
        stock.subtract_circle_sweep(float(circle.center.x), float(circle.center.y), float(circle.radius.value), float(tool_radius.value))
    return stock
