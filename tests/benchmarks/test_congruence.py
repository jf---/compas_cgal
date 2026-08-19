from __future__ import annotations

import math

import pytest

from benchmarks.congruence import PYTHAGOREAN_ROTATIONS
from benchmarks.congruence import rotate_point
from benchmarks.congruence import rotate_spec
from benchmarks.congruence import translate_spec
from benchmarks.families.analytic import stadium
from benchmarks.spec import PocketSpec
from compas_cgal import _stock_2
from compas_cgal.stock import Stock

# 4*sin^2(pi/2) == 4: the largest legal squared-chord surrogate, i.e. a 180 deg
# cap. The whole engaged rim therefore counts and total_tea is effectively
# uncapped, which makes the congruence comparison a statement about the measured
# engagement rather than about where a cap happened to clip it.
CAP_RATIO_FULL = 4.0

# 4*sin^2(0) == 0: no void gap is pre-absorbed, so the pessimistic runs coincide
# with the true runs. Any congruence violation observed is then attributable to
# the kernel's own seam handling rather than to gap-closure bookkeeping.
GAP_CLOSE_NONE = 0.0

TOOL_RADIUS = 0.5

# Congruent inputs must agree to the last bit the exact kernel reports; this is a
# strictly-positive floor standing in for exact equality, far below the ~1e-2 rad
# discrepancies a genuine seam bug produces.
ENGAGEMENT_ABS_TOL = 1e-9

# A rational rotation is exact, so the rotated unit vector's norm is limited only
# by the two double divisions a/c and b/c -- a few ulps.
UNIT_NORM_ABS_TOL = 1e-15

# The cut disk is offset from the probe centre by 1.2 tool radii along the local
# channel axis: far enough that the tool overlaps the void on one side only, so a
# rotation-sensitive CCW seam has a real run boundary to mishandle.
CUT_OFFSET = -0.6


def _engagement(spec: PocketSpec, probe_x: float, probe_y: float, cut_x: float, cut_y: float) -> float:
    """Total engagement angle of a probe circle against a once-cut stock.

    Args:
        spec: The pocket instance supplying the stock boundary.
        probe_x: Probe (cutter) centre x.
        probe_y: Probe (cutter) centre y.
        cut_x: Centre x of the disk already removed from the stock.
        cut_y: Centre y of the disk already removed from the stock.

    Returns:
        The total engaged angle in radians.
    """
    stock = Stock(spec.polygon, list(spec.holes))
    stock.subtract_disk(cut_x, cut_y, TOOL_RADIUS)
    total, _max_run, _exceeded = _stock_2.engagement_at(stock.raw, probe_x, probe_y, TOOL_RADIUS, CAP_RATIO_FULL, GAP_CLOSE_NONE)
    return float(total)


def test_translation_along_the_channel_preserves_engagement() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    base = _engagement(spec, 0.0, 0.0, CUT_OFFSET, 0.0)
    for dx in (1.0, 2.5, -3.0):
        moved = translate_spec(spec, dx, 0.0)
        assert _engagement(moved, dx, 0.0, dx + CUT_OFFSET, 0.0) == pytest.approx(base, abs=ENGAGEMENT_ABS_TOL)


def test_rational_rotation_preserves_engagement() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    base = _engagement(spec, 0.0, 0.0, CUT_OFFSET, 0.0)
    for triple in PYTHAGOREAN_ROTATIONS:
        turned = rotate_spec(spec, triple)
        px, py = rotate_point(0.0, 0.0, triple)
        cx, cy = rotate_point(CUT_OFFSET, 0.0, triple)
        assert _engagement(turned, px, py, cx, cy) == pytest.approx(base, abs=ENGAGEMENT_ABS_TOL)


def test_rotation_coefficients_are_exactly_rational() -> None:
    for a, b, c in PYTHAGOREAN_ROTATIONS:
        assert a * a + b * b == c * c
        x, y = rotate_point(1.0, 0.0, (a, b, c))
        assert math.hypot(x, y) == pytest.approx(1.0, abs=UNIT_NORM_ABS_TOL)
