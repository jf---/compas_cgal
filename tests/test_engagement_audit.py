"""Stock wrapper and toolpath engagement-audit tests (Task 6).

Task 6a covers the typed `Stock` wrapper roundtrip; Task 6b appends the
engagement-audit replay tests that consume `engagement.py`.
"""

import math

import numpy as np
import pytest
from compas.geometry import Arc, Frame, Point, Polygon

from compas_cgal.engagement import (
    EngagementReport,
    InvalidEngagementCapError,
    InvalidToolDiameterError,
    _certify_arc_engagement,
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


# --------------------------------------------------------------------------- #
# Circular-motion run-merge hole and its gap-closure repair (arc certifier)    #
# --------------------------------------------------------------------------- #

# Merge witness, mirroring the LINEAR witness in tests/test_stock.py: the lower
# half-plane (y <= 5) with a thin central void slot at x = 5 biting the cutter
# rim's bottom. A radius-0.5 cutter centred at (5, 5) engages the lower
# semicircle; the slot splits it into two symmetric ~84.26 deg runs separated by
# a ~11.5 deg void gap. As the centre slides off x = 5 the two runs MERGE into
# one ~180 deg run -- an O(1) jump the O(sqrt d) growth lemma cannot bridge, so a
# station that measures only the true (unmerged) 84 deg run FALSE-PASSES a cap
# set above 84 deg. Gap-closure pessimism pre-absorbs the ~11.5 deg gap at the
# station, exposing the ~180 deg merged run before the merge physically completes.
MERGE_LOWER_HALF = Polygon([[0, 0, 0], [10, 0, 0], [10, 5, 0], [0, 5, 0]])
MERGE_TOOL_RADIUS = 0.5


def _slotted_lower_half_stock() -> Stock:
    """Lower half-plane with the thin central rim-biting void slot (merge witness)."""
    stock = Stock(MERGE_LOWER_HALF)
    stock.subtract_capsule(5.0, 4.30, 5.0, 4.53, 0.05)
    return stock


def _bottom_anchored_arc(center_y: float, radius: float, sweep_deg: float = 36.0) -> Arc:
    """Arc whose tool-centre path STARTS at its lowest point ``(5, center_y - radius)``.

    Frame origin (centre of curvature) at ``(5, center_y)`` with ``start_angle`` at
    the -y direction (270 deg): ``point_at(0.0)`` is exactly ``(5, center_y - radius)``
    and, as ``params[0]``, is unconditionally a sampled station regardless of how
    the fixed-density ``ceil`` rounds. Sweeping ``sweep_deg`` to one side lifts the
    centre away from the material, so engagement is deepest at the anchored start.
    """
    frame = Frame(Point(5.0, center_y, 0.0))
    return Arc(radius=radius, start_angle=math.radians(270.0), end_angle=math.radians(270.0 + sweep_deg), frame=frame)


def test_arc_certifier_refuses_circular_merge_over_cap():
    """RED before the fix: a circular cut motion whose stations each split into two
    sub-cap engaged runs across a thin CLOSING void must NOT certify.

    The tool-centre arc is anchored at its lowest point exactly on the material
    boundary at ``(5, 5)`` (``params[0]``) and rises 36 deg to one side, so every
    sampled station sits in the two-run regime of the slotted lower half-plane.
    At the fixed station density the guard is ``gamma_guard ~44 deg`` (>> the
    ~11.5 deg gap) and the guarded cap is ``pi - gamma_guard ~136 deg``, which
    lies strictly between the single-run span (~84 deg) and the merged span
    (~180 deg).

    Pre-fix the arc certifier passes NO ``gap_close_ratio`` (defaults 0): each
    station sees only its larger ~84 deg run, nothing exceeds ~136 deg, and the
    motion FALSE-PASSES as ``cap_certified=True`` -- the O(1) run-merge slipping
    between stations. Post-fix the certifier threads ``gamma_guard`` so the thin
    gap is pre-absorbed, the ~180 deg merged pessimistic run exceeds ~136 deg, and
    the motion is correctly refused. The reported ``max_tea`` stays the TRUE
    (unmerged) sub-cap measure throughout -- the deciding/reporting split: the peak
    reported run is well under a half-turn, yet the motion is refused because the
    DECISION accounts for the merge the reported peak cannot show.
    """
    stock = _slotted_lower_half_stock()
    arc = _bottom_anchored_arc(center_y=5.0 + 0.05, radius=0.05)

    max_tea, cap_certified, stations = _certify_arc_engagement(stock, arc, MERGE_TOOL_RADIUS, math.pi)

    assert cap_certified is False  # pre-fix FALSE PASS (returns True); post-fix refuses
    assert stations >= 2  # a genuine multi-station circular walk, not a single probe
    # Deciding/reporting split: reported peak run stays the TRUE ~84 deg measure,
    # far below the pi cap -- the refusal comes from the merge DECISION, not max_tea.
    assert max_tea < math.radians(90.0)


def test_arc_certifier_preserves_benign_single_run_immersion():
    """Preservation: gap-closure pessimism does NOT spuriously refuse an ordinary arc.

    On the UNSLOTTED lower half-plane a radius-0.5 cutter whose lowest station sits
    at ``(5, 5.3)`` engages a single contiguous ~106 deg run (its deepest bite; the
    arc lifts away from there). There is no thin void to absorb, so ``gap_close_ratio``
    is a no-op: the guarded cap ``pi - gamma_guard ~136 deg`` clears the ~106 deg run
    at every station and the motion certifies ``True``. This is the exact fixed-density
    path that carries the merge repair, exercised on benign geometry to prove the
    repair narrows to genuine merges and leaves ordinary immersion arcs certified.
    """
    stock = Stock(MERGE_LOWER_HALF)  # lower half-plane, NO slot
    arc = _bottom_anchored_arc(center_y=5.0 + 0.30 + 0.05, radius=0.05)  # deepest station at (5, 5.3)

    max_tea, cap_certified, stations = _certify_arc_engagement(stock, arc, MERGE_TOOL_RADIUS, math.pi)

    assert cap_certified is True  # a genuine single sub-cap run is not spuriously refused
    assert max_tea < math.pi  # single-run immersion, well under the half-turn cap
    assert max_tea == pytest.approx(math.radians(106.26), abs=math.radians(1.0))  # the ~106 deg bite


def test_arc_certifier_certifies_non_engaged_arc():
    """Preservation: an arc that never touches material (rapid-style clearance) certifies.

    With the cutter centre lifted so the whole tool clears the ``y <= 5`` boundary,
    every station measures zero engagement; the guarded cap is positive at this
    density, so gap-closure runs but finds nothing to merge and the motion certifies
    ``True`` with ``max_tea == 0``.
    """
    stock = Stock(MERGE_LOWER_HALF)  # lower half-plane, NO slot
    arc = _bottom_anchored_arc(center_y=5.60 + 0.05, radius=0.05)  # tool bottom 5.1 > 5: no contact

    max_tea, cap_certified, stations = _certify_arc_engagement(stock, arc, MERGE_TOOL_RADIUS, math.pi)

    assert cap_certified is True
    assert max_tea == pytest.approx(0.0, abs=1e-9)
