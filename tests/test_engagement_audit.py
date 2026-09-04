"""Stock wrapper and toolpath engagement-audit tests (Task 6).

Task 6a covers the typed `Stock` wrapper roundtrip; Task 6b appends the
engagement-audit replay tests that consume `engagement.py`.
"""

import math

import numpy as np
import pytest
from compas.geometry import Arc, Circle, Frame, Line, Point, Polygon
from compas.tolerance import TOL

from compas_cgal import _stock_2
from compas_cgal.engagement import (
    AUDIT_ENGAGED,
    EngagementReport,
    InvalidEngagementCapError,
    InvalidToolDiameterError,
    UnexpectedToolpathGeometryError,
    _cap_chord_ratio,
    _certify_arc_engagement,
    _infer_cut_height,
    _subtract_operation,
    audit_toolpath_engagement,
)
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType, ToolpathOperation, ToolpathResult, trochoidal_mat_toolpath_circular

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


def test_audit_rejects_curve_below_line_authorized_cut_plane() -> None:
    operations = [
        ToolpathOperation(geometry=Line([1.0, 1.0, 0.0], [5.0, 1.0, 0.0]), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(
            geometry=Circle(1.0, frame=Frame([3.0, 3.0, -1.0])),
            operation=OperationType.CUT,
            path_index=0,
        ),
    ]
    result = ToolpathResult(operations=operations, polyline=np.empty((0, 3), dtype=np.float64))

    with pytest.raises(UnexpectedToolpathGeometryError):
        audit_toolpath_engagement(SQUARE, result, tool_diameter=2.0, tea_cap=math.pi)


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
    """The arc certifier's gap-closure wiring: a circular cut whose stations each split
    into two sub-cap runs across a thin void must be REFUSED (RED before the fix).

    The tool-centre arc is anchored at its lowest point exactly on the material
    boundary at ``(5, 5)`` (``params[0]``, unconditionally sampled regardless of the
    fixed-density ``ceil``) and rises 36 deg to one side, so every sampled station
    sits in the two-run regime of the slotted lower half-plane (two ~84 deg runs
    across a ~11.5 deg gap). At the fixed density the guard is ``gamma_guard ~44 deg``
    and the guarded cap ``pi - gamma_guard ~136 deg`` lies strictly between the
    single-run span (~84 deg) and the merged span (~180 deg).

    Pre-fix the certifier passes NO ``gap_close_ratio`` (defaults 0): each station
    sees only its larger ~84 deg run, nothing exceeds ~136 deg, and the motion
    certifies ``True``. Post-fix the certifier threads ``gamma_guard`` so the thin
    gap is pre-absorbed and the ~180 deg pessimistic merged run exceeds ~136 deg ->
    the motion is refused.

    HONEST SCOPE (verified by a dense exact scan): on THIS motion the two runs never
    actually fuse -- the true peak run is ~86 deg and the centre never submerges -- so
    ``False`` is the accepted OVER-CONSERVATIVE gap-closure refusal, NOT a demonstrated
    genuine merge hole (a robust certifier-level RED for a real hole needs a brittle
    sub-floor construction). Its value is as a WIRING-DISCRIMINATION GUARD: it fails
    iff the certifier stops threading ``gap_close_ratio``. It also exercises the
    deciding/reporting split -- reported ``max_tea`` stays the TRUE sub-cap measure
    while the DECISION refuses.
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


# --------------------------------------------------------------------------- #
# Rim-span normalisation: no run may report more engagement than a full turn   #
# --------------------------------------------------------------------------- #

# The 12x8 reference pocket, driven by the UNREGULATED generator at tool 2.0 --
# the instance on which the rim-span defect was measured. Its 99 cutting motions
# graze the depleting stock often enough that a rim sub-arc of ~1e-15 rad occurs.
SPAN_POCKET = Polygon([[0, 0, 0], [12, 0, 0], [12, 8, 0], [0, 8, 0]])
SPAN_TOOL_DIAMETER = 2.0
SPAN_CLEARANCE_Z = 2.0

FULL_TURN = 2.0 * math.pi

# Mirror of FULL_TURN_REPORTING_SLACK in src/engagement_2.cpp: the engaged rim
# sub-arcs partition the cutter circle, so their reported spans sum to at most a
# full turn EXACTLY, and the only error is the rounding of each span evaluation
# and of their summation -- a few ulps of 2*pi (~8.9e-16). 1e-12 rad is ~1000x
# that and still 12 orders of magnitude below the failure it guards, which
# over-reports by a WHOLE TURN. Reporting only; no verdict consults it.
FULL_TURN_REPORTING_SLACK = 1e-12

# A station reporting within this much of a full turn is claiming the ENTIRE rim
# sits in material. 1e-3 rad (0.057 deg) is far below any real feature at this
# scale and far above the ~1e-15 rad summation noise, so the classification is
# unambiguous in both directions.
BURIED_RIM_TOLERANCE_RAD = 1e-3

# Rim sample points used to refute a "whole rim is buried" claim. Eight is ample:
# the claim asserts EVERY rim point is in material, so one counter-example
# suffices, and a spurious full turn comes from a near-zero contact arc whose
# cutter sits almost entirely in cleared void.
BURIED_RIM_SAMPLES = 8

# Sample radius as a fraction of the tool radius: just inside the rim, so a
# genuinely buried rim's samples are strictly interior points and `contains` is
# never asked about a boundary point (where membership is a measure-zero coin
# flip rather than a fact about the engagement).
BURIED_RIM_PROBE_FRACTION = 0.999

# Stations per motion for the replay scan. The endpoints are always included, and
# the defect's witness station is a motion START, so a coarse walk finds it.
SPAN_SCAN_STATIONS = 12


def _span_pocket_toolpath() -> ToolpathResult:
    return trochoidal_mat_toolpath_circular(SPAN_POCKET, tool_diameter=SPAN_TOOL_DIAMETER, clearance_z=SPAN_CLEARANCE_Z)


def _motion_stations(op) -> list:
    """Cutter-centre stations along one motion, endpoints included."""
    g = op.geometry
    if isinstance(g, (Arc, Circle)):
        return [g.point_at(i / SPAN_SCAN_STATIONS) for i in range(SPAN_SCAN_STATIONS)]
    return [
        (
            float(g.start[0]) + (float(g.end[0]) - float(g.start[0])) * i / SPAN_SCAN_STATIONS,
            float(g.start[1]) + (float(g.end[1]) - float(g.start[1])) * i / SPAN_SCAN_STATIONS,
        )
        for i in range(SPAN_SCAN_STATIONS + 1)
    ]


def _rim_outside_material(stock: Stock, px: float, py: float, tool_radius: float) -> list:
    """Rim sample points at ``(px, py)`` that are NOT in *stock*'s material.

    Sampled just inside the rim so a genuinely buried cutter yields strictly
    interior points and membership is never asked about a boundary point. Must be
    read while *stock* still holds the state the station was measured against.
    """
    probe_radius = BURIED_RIM_PROBE_FRACTION * tool_radius
    outside = []
    for k in range(BURIED_RIM_SAMPLES):
        angle = 2.0 * math.pi * k / BURIED_RIM_SAMPLES
        sx = px + probe_radius * math.cos(angle)
        sy = py + probe_radius * math.sin(angle)
        if not stock.contains(sx, sy):
            outside.append((sx, sy))
    return outside


def _scan_engagement(result: ToolpathResult) -> list:
    """Replay *result*, reading the engagement at every station of every cut motion.

    Mirrors the audit's replay (measure against the CURRENT stock, then subtract)
    but keeps ``total_tea``, which the audit discards -- and which is where a
    mis-normalised rim sub-arc shows up first. A station claiming a buried rim has
    its rim sampled HERE, against the stock the measurement actually saw; the
    stock is destroyed by the very next subtraction, so the check cannot be
    deferred to the caller.

    Returns:
        One ``(op_index, x, y, total_tea, buried_claim, rim_outside)`` row per
        station, where ``rim_outside`` lists the rim samples found in void and is
        empty for every station that made no buried-rim claim.
    """
    tool_radius = 0.5 * SPAN_TOOL_DIAMETER
    ratio = _cap_chord_ratio(math.radians(120.0))
    stock = Stock(SPAN_POCKET)
    cut_z = _infer_cut_height(result.operations)
    rows = []
    for index, op in enumerate(result.operations):
        if op.operation == OperationType.RETRACT:
            continue
        g = op.geometry
        if isinstance(g, Line):
            z0, z1 = float(g.start[2]), float(g.end[2])
            if abs(z0 - z1) > TOL.absolute:
                if z1 < z0:
                    stock.subtract_disk(float(g.end[0]), float(g.end[1]), tool_radius)
                continue
            if z0 > cut_z + TOL.absolute:
                continue
        if op.operation not in AUDIT_ENGAGED:
            continue
        for p in _motion_stations(op):
            px, py = float(p[0]), float(p[1])
            total, _max_run, _exceeded = _stock_2.engagement_at(stock.raw, px, py, tool_radius, ratio, 0.0)
            buried_claim = total >= FULL_TURN - BURIED_RIM_TOLERANCE_RAD
            outside = _rim_outside_material(stock, px, py, tool_radius) if buried_claim else []
            rows.append((index, px, py, total, buried_claim, outside))
        _subtract_operation(stock, op, tool_radius)
    return rows


def test_no_operation_reports_more_than_a_full_turn():
    """No audited operation may report engagement beyond one full turn of the rim.

    The engaged sub-arcs PARTITION the cutter circle, so an assembled run covers
    it at most once. Before the rim-span repair this pocket reported
    ``report.max_tea == 6.28318530717959`` -- above ``2*pi``, with operation 32
    claiming a whole buried rim where the true contact was ``5.7e-15`` rad. It now
    reports exactly ``2*pi``, from operation 1's genuinely buried virgin-stock cut.

    Asserted against the named reporting slack rather than a bare ``<= 2*pi``
    because the reported number is a SUM of rounded spans: a run assembled from
    several sub-arcs may legitimately overshoot by an ulp. The slack is 12 orders
    of magnitude below the whole-turn error it exists to catch;
    ``test_a_full_turn_report_means_a_buried_rim`` is the discriminating half of
    this pair.
    """
    result = _span_pocket_toolpath()
    report = audit_toolpath_engagement(SPAN_POCKET, result, tool_diameter=SPAN_TOOL_DIAMETER, tea_cap=math.radians(120.0))

    assert report.engaged_ops > 0  # a vacuous pass would satisfy the bound trivially
    over = [(e.op_index, e.max_tea) for e in report.operations if e.max_tea > FULL_TURN + FULL_TURN_REPORTING_SLACK]
    assert over == [], f"operations reporting beyond a full turn: {over}"
    assert report.max_tea <= FULL_TURN + FULL_TURN_REPORTING_SLACK


def test_a_full_turn_report_means_a_buried_rim():
    """A station reporting a full turn must actually have its whole rim in material.

    The physical cross-check behind the reported angle, and the discriminating
    regression for the rim-span defect: a degenerate sub-arc promoted to a full
    turn claims a buried rim at a station whose cutter sits in cleared void, which
    a handful of rim samples refutes outright. Both branches are exercised on this
    pocket -- operation 1 cuts virgin stock with a genuinely buried rim (the claim
    holds), operation 32 grazes at ``5.7e-15`` rad (the claim must not be made).
    """
    rows = _scan_engagement(_span_pocket_toolpath())

    for index, px, py, total, _buried_claim, _outside in rows:
        assert total <= FULL_TURN + FULL_TURN_REPORTING_SLACK, f"op {index} at ({px}, {py}) reports total_tea {total!r} beyond a full turn"

    claims = [row for row in rows if row[4]]
    assert claims, "no station claimed a buried rim, so the check never ran"
    broken = [(index, px, py, total, outside) for index, px, py, total, _claim, outside in claims if outside]
    assert broken == [], f"stations claiming a buried rim whose rim is in void: {broken}"


# --------------------------------------------------------------------------- #
# Baseline audit script (Task 8) -- SP2 reference generator                    #
# --------------------------------------------------------------------------- #


def test_baseline_script_generates_report(tmp_path):
    """`scripts/engagement_baseline.py` runs end-to-end and emits a BLUF-first report + JSON.

    The full 7-pocket baseline is O(n^2) exact stock depletion (~1 hour) and its generation is
    DEFERRED to SP2 (docs/superpowers/state/sp1-gate-c-analysis.md) -- the committed artifact is
    the harness, not baseline data. This test drives the SAME generate -> audit -> report code
    path over a tiny synthetic smoke pocket via ``--quick`` (~8 s). It guards the report
    CONTRACT -- BLUF heading, the "worst TEA" phrasing, the SP2-reference sentence, and a JSON
    sidecar whose per-op numbers round-trip through ``json.loads``.
    """
    import json
    import os
    import pathlib
    import subprocess
    import sys

    repo_root = pathlib.Path(__file__).resolve().parent.parent
    env = dict(os.environ, PYTHONPATH="src")
    out = subprocess.run(
        [sys.executable, "scripts/engagement_baseline.py", "--out", str(tmp_path), "--quick"],
        capture_output=True,
        text=True,
        env=env,
        cwd=repo_root,
        timeout=120,
    )
    assert out.returncode == 0, out.stderr

    md = (tmp_path / "engagement_baseline.md").read_text()
    assert md.splitlines()[0].startswith("# Engagement baseline")
    assert "worst TEA" in md
    assert "This baseline is the reference SP2 must beat." in md

    data = json.loads((tmp_path / "engagement_baseline.json").read_text())
    assert data["pockets"], "expected at least one audited pocket"
    first = data["pockets"][0]
    for key in ("name", "tool_diameter", "max_tea_deg", "cap_violations", "stations", "wall_clock_s"):
        assert key in first, f"missing per-op field {key!r}"
    assert "total_wall_clock_s" in data
