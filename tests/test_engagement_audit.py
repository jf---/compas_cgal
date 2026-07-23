"""Stock wrapper and toolpath engagement-audit tests (Task 6).

Task 6a covers the typed `Stock` wrapper roundtrip; Task 6b appends the
engagement-audit replay tests that consume `engagement.py`.
"""

import math

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal.engagement import (
    EngagementReport,
    InvalidEngagementCapError,
    InvalidToolDiameterError,
    audit_toolpath_engagement,
)
from compas_cgal.stock import Stock
from compas_cgal.toolpath import ToolpathResult, trochoidal_mat_toolpath_circular

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


def test_stock_wrapper_roundtrip():
    stock = Stock(SQUARE)
    assert stock.contains(5.0, 5.0)
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)


# A full fill of SQUARE audits in >120s: every replayed cut removes material from
# the exact stock, and the depleting boolean grows more expensive per op as cuts
# accumulate. max_passes=1 caps the toolpath to a single skeleton chain (which
# truncates coverage, hence the warning filter) while still cutting virgin stock
# at full immersion -- exactly the regime these audits probe. Determinism is
# input-agnostic, so its test uses a smaller pocket to keep two audits fast.
DETERMINISM_POCKET = Polygon([[0, 0, 0], [6, 0, 0], [6, 6, 0], [0, 6, 0]])


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_audit_square_trochoid_toolpath():
    result = trochoidal_mat_toolpath_circular(SQUARE, tool_diameter=2.0, pitch=1.0, max_passes=1)
    report = audit_toolpath_engagement(SQUARE, result, tool_diameter=2.0, tea_cap=math.radians(120))
    assert isinstance(report, EngagementReport)
    assert report.engaged_ops > 0
    # Every operation is accounted for, engaged or not (skipped moves still map 1:1).
    assert len(report.operations) == len(result.operations)
    assert report.max_tea > 0.0
    # The existing generator has NO engagement regulation: its first motion cuts
    # virgin stock at full immersion, so the audit reports it truthfully rather
    # than flattering the toolpath. Full slotting far exceeds the 120 deg cap.
    assert report.max_tea > math.radians(150)
    assert report.cap_violations > 0


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_audit_is_deterministic():
    result = trochoidal_mat_toolpath_circular(DETERMINISM_POCKET, tool_diameter=2.0, pitch=1.0, max_passes=1)
    r1 = audit_toolpath_engagement(DETERMINISM_POCKET, result, tool_diameter=2.0, tea_cap=math.pi)
    r2 = audit_toolpath_engagement(DETERMINISM_POCKET, result, tool_diameter=2.0, tea_cap=math.pi)
    # The replay is a pure function of (stock, ops): exact-kernel measurements and
    # deterministic float arithmetic make two audits of one input agree exactly.
    assert r2.max_tea == pytest.approx(r1.max_tea)
    assert r2.cap_violations == r1.cap_violations


def _empty_result() -> ToolpathResult:
    return ToolpathResult(operations=[], polyline=np.empty((0, 3), dtype=np.float64))


@pytest.mark.parametrize("bad_cap", [0.0, -0.1, math.pi + 1e-6, 2.0 * math.pi])
def test_audit_rejects_out_of_range_cap(bad_cap):
    # The cap is validated at the boundary before any replay work (a run subtends
    # at most a half turn before the > pi case is an exact orientation verdict).
    with pytest.raises(InvalidEngagementCapError):
        audit_toolpath_engagement(SQUARE, _empty_result(), tool_diameter=2.0, tea_cap=bad_cap)


@pytest.mark.parametrize("bad_diameter", [0.0, -1.0])
def test_audit_rejects_nonpositive_tool_diameter(bad_diameter):
    with pytest.raises(InvalidToolDiameterError):
        audit_toolpath_engagement(SQUARE, _empty_result(), tool_diameter=bad_diameter, tea_cap=math.pi)
