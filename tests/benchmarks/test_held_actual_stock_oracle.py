"""Actual-stock differential oracle (unit U0): the circle view C_t, the typed part polygon P, the D2 / D4 harnesses."""

from __future__ import annotations

import math

import numpy as np
import pytest

from benchmarks import held_actual_stock_oracle as oracle_module
from benchmarks.held_actual_stock_last16 import CIRCLE
from benchmarks.held_actual_stock_last16 import LAST_M
from benchmarks.held_actual_stock_last16 import SEGMENT
from benchmarks.held_actual_stock_last16 import BoundBreakError
from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_last16 import Verdict
from benchmarks.held_actual_stock_oracle import SKIS
from benchmarks.held_actual_stock_oracle import CertificateViolationError
from benchmarks.held_actual_stock_oracle import CircleViewRefusalError
from benchmarks.held_actual_stock_oracle import EntryCircleNotRefusedError
from benchmarks.held_actual_stock_oracle import OracleWorkerNotInitializedError
from benchmarks.held_actual_stock_oracle import circle_view_vs_actual
from benchmarks.held_actual_stock_oracle import production_implies_certificate
from benchmarks.held_actual_stock_oracle import view_prefixes
from benchmarks.held_actual_stock_oracle import workload_on_part
from benchmarks.held_material_domain import CentreDomainPolygon
from benchmarks.held_material_domain import DegeneratePolygonError
from benchmarks.held_material_domain import MaterialDomainMismatchError
from benchmarks.held_material_domain import PartPolygon
from benchmarks.held_material_domain import require_part_material
from benchmarks.held_toolpath_benchmark import generate
from compas_cgal import _stock_2
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY
from tools.actual_stock_profile import CircleSweep
from tools.actual_stock_profile import ReplayWorkload
from tools.actual_stock_profile import SegmentSweep

TOOL = 1.0
CAP_DEG = 80.0
WIDE = ((-10.0, -10.0), (50.0, -10.0), (50.0, 10.0), (-10.0, 10.0))
FAR_X = 30.0


def _circle(cx: float, cy: float, rho: float) -> CircleSweep:
    return CircleSweep(cx, cy, rho, cx + rho, cy)


def _rows(corners: tuple[tuple[float, float], ...]) -> np.ndarray:
    return np.array([(x, y, 0.0) for x, y in corners], dtype=np.float64)


def _workload(circles: list[CircleSweep], segments: tuple[SegmentSweep, ...] = (), corners: tuple[tuple[float, float], ...] = WIDE) -> ReplayWorkload:
    return ReplayWorkload(
        name="synthetic",
        polygon=_rows(corners),
        tool_radius=TOOL,
        cap_chord_ratio=_stock_2.cap_chord_ratio(math.radians(CAP_DEG)),
        circles=tuple(circles),
        segments=segments,
        generation_s=0.0,
    )


def _connector_cleared() -> ReplayWorkload:
    """Circle 1 sits inside a band cleared ONLY by circle 0's two connectors: within on S_1, a plunge on C_1."""
    band = (SegmentSweep(0, 4.0, 0.5, 12.0, 0.5), SegmentSweep(0, 4.0, -0.5, 12.0, -0.5))
    return _workload([_circle(0.0, 0.0, 3.0), _circle(8.0, 0.0, 0.25)], band)


def _escalation() -> ReplayWorkload:
    """Circle 0 big; LAST_M plunges far away; circle 0 again (refused by last-16, certified on the full window)."""
    fillers = [_circle(FAR_X + 2.0 * k / LAST_M, 0.0, 0.5) for k in range(LAST_M)]
    return _workload([_circle(0.0, 0.0, 3.0), *fillers, _circle(0.0, 0.0, 3.0)])


def _overlapping() -> ReplayWorkload:
    """A plunge, then two circles that each read only the previous circle's cut: every t >= 1 is certified."""
    return _workload([_circle(0.0, 0.0, 3.0), _circle(0.0, 0.0, 3.0), _circle(0.0, 0.0, 3.0)], (SegmentSweep(0, 3.0, 0.0, 3.5, 0.7),))


def _entry_clear() -> ReplayWorkload:
    """A tiny part far inside circle 0's guide circle: no station of circle 0 reads material, so the entry circle is certified."""
    return _workload([_circle(0.0, 0.0, 5.0)], corners=((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)))


# ---- the circle view in Last16Oracle -----------------------------------------------


def test_circle_view_selects_prior_circles_and_no_connector() -> None:
    w = _connector_cleared()
    actual = Last16Oracle(w).candidates(1, None)
    circles = Last16Oracle(w, sweeps=SweepView.CIRCLES).candidates(1, None)
    assert [e[1] for e in actual].count(SEGMENT) == 2
    assert circles == [e for e in actual if e[1] == CIRCLE] == [(0, CIRCLE, 0)]


def test_circle_view_stock_equals_its_reference_and_contains_the_actual_stock() -> None:
    w = _connector_cleared()
    view = Last16Oracle(w, sweeps=SweepView.CIRCLES)
    actual = Last16Oracle(w)
    for last in (LAST_M, None):
        _v, circle_stock = view.view(1, last)
        view.audit_stock(1, circle_stock, last)
        _a, actual_stock = actual.view(1, last)
        assert actual_stock.is_subset_of(circle_stock)
        assert not circle_stock.is_subset_of(actual_stock), "the connectors cut material the circle view keeps"


def test_connector_cleared_circle_is_certified_on_the_actual_stock_and_refused_on_the_circle_view() -> None:
    w = _connector_cleared()
    assert Last16Oracle(w).certify(1).verdict is Verdict.CERTIFIED_LAST16
    assert Last16Oracle(w, sweeps=SweepView.CIRCLES).certify(1).verdict is Verdict.REFUSED


def test_view_prefixes_are_the_reference_engine_on_the_window_stock() -> None:
    w = _connector_cleared()
    circles = view_prefixes(Last16Oracle(w, sweeps=SweepView.CIRCLES), 1, None)
    actual = view_prefixes(Last16Oracle(w), 1, None)
    assert len(circles) == len(actual) == 45
    assert any(e for *_x, e in circles) and not any(e for *_x, e in actual)
    assert all(c[1] >= a[1] for c, a in zip(circles, actual)), "max run is monotone under stock inclusion"


# ---- the typed part polygon ------------------------------------------------------------


def test_domain_polygons_refuse_fewer_than_three_vertices() -> None:
    two = [Point2[WorldXY].build(0.0, 0.0), Point2[WorldXY].build(1.0, 0.0)]
    with pytest.raises(DegeneratePolygonError):
        PartPolygon.build(two)
    with pytest.raises(DegeneratePolygonError):
        CentreDomainPolygon.build(two)


def test_require_part_material_is_exact() -> None:
    part = PartPolygon.from_stock_rows(_rows(WIDE))
    require_part_material(_rows(WIDE), part, "fixture")
    nudged = _rows(WIDE)
    nudged[0, 0] = math.nextafter(nudged[0, 0], 0.0)
    with pytest.raises(MaterialDomainMismatchError):
        require_part_material(nudged, part, "fixture")
    with pytest.raises(MaterialDomainMismatchError):
        require_part_material(_rows(WIDE[:3]), part, "fixture")


def test_wall_band_is_material_on_the_part_and_absent_on_the_offset() -> None:
    """Station rims on the wall: the r-wide band beyond the offset is the only material, so the offset accepts what P refuses."""
    plunge = [CircleSweep(0.0, 0.0, 0.01, 0.01, 0.0)]
    part = _workload(plunge, corners=((0.0, -10.0), (20.0, -10.0), (20.0, 10.0), (0.0, 10.0)))
    offset = _workload(plunge, corners=((TOOL, -10.0 + TOOL), (20.0 - TOOL, -10.0 + TOOL), (20.0 - TOOL, 10.0 - TOOL), (TOOL, 10.0 - TOOL)))
    assert Last16Oracle(part).certify(0).verdict is Verdict.REFUSED
    assert Last16Oracle(offset).certify(0).verdict is Verdict.CERTIFIED_LAST16


def test_workload_on_a_generated_path_decides_on_the_part_polygon() -> None:
    path = generate(SKIS)
    workload = workload_on_part(path, 0.0)
    part = PartPolygon.from_case(path.case)
    require_part_material(workload.polygon, part, SKIS)
    centre = CentreDomainPolygon.build(path.boundary)
    assert len(centre.points) != len(part.points) or centre.points != part.points
    assert Last16Oracle(workload).material == part


# ---- D4: circle view vs actual ------------------------------------------------------------


@pytest.mark.parametrize("processes", [1, 2])
def test_d4_reports_the_connector_cleared_circle(processes: int) -> None:
    w = _connector_cleared()
    report = circle_view_vs_actual(w, range(len(w.circles)), processes)
    assert report.compared == 2
    assert report.violations == (1,)
    assert [(r.circle, r.circle_view, r.actual) for r in report.rows] == [(0, Verdict.REFUSED, Verdict.REFUSED), (1, Verdict.REFUSED, Verdict.CERTIFIED_LAST16)]
    with pytest.raises(CircleViewRefusalError):
        report.require()


def test_d4_is_zero_without_connectors_and_skips_the_actual_on_circle_view_certificates() -> None:
    w = _escalation()
    report = circle_view_vs_actual(w, range(len(w.circles)), 1)
    report.require()
    assert report.violations == ()
    last = len(w.circles) - 1
    assert report.rows[last].circle_view is Verdict.CERTIFIED_FULL_WINDOW
    assert report.rows[last].actual is None, "a circle-view certificate needs no actual decision"
    assert all(r.actual is Verdict.REFUSED for r in report.rows[:-1])


def test_d4_audit_decides_the_actual_on_sampled_circles() -> None:
    w = _overlapping()
    report = circle_view_vs_actual(w, range(len(w.circles)), 1, audit_every=2)
    assert [r.actual is not None for r in report.rows] == [True, False, True]
    assert report.audited == 2
    assert report.violations == (), "an audited circle certified on both views is no violation"


def test_d4_bound_break_is_raised_when_the_circle_view_certifies_an_actual_refusal(monkeypatch: pytest.MonkeyPatch) -> None:
    w = _connector_cleared()
    real = oracle_module._oracle
    swapped = {SweepView.CIRCLES: SweepView.ACTUAL, SweepView.ACTUAL: SweepView.CIRCLES}
    monkeypatch.setattr(oracle_module, "_oracle", lambda view: real(swapped[view]))
    with pytest.raises(BoundBreakError):
        circle_view_vs_actual(w, range(len(w.circles)), 1, audit_every=1)


def test_unstarted_oracle_worker_raises_a_named_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(oracle_module, "_WORKLOAD", None)
    with pytest.raises(OracleWorkerNotInitializedError):
        oracle_module._oracle(SweepView.ACTUAL)


# ---- D2: production implies certificate ------------------------------------------------------


def test_d2_passes_when_every_later_circle_is_certified_and_the_entry_refused() -> None:
    w = _overlapping()
    report = production_implies_certificate(w, range(len(w.circles)), 1)
    report.require()
    assert report.entry is Verdict.REFUSED and report.violations == () and report.compared == 3


def test_d2_names_refused_later_circles() -> None:
    w = _escalation()
    report = production_implies_certificate(w, range(len(w.circles)), 2)
    assert report.violations == tuple(range(1, LAST_M + 1))
    with pytest.raises(CertificateViolationError):
        report.require()


def test_d2_on_the_circle_view_refuses_the_connector_cleared_circle() -> None:
    w = _connector_cleared()
    production_implies_certificate(w, range(len(w.circles)), 1).require()
    report = production_implies_certificate(w, range(len(w.circles)), 1, sweeps=SweepView.CIRCLES)
    assert report.violations == (1,)


def test_d2_refuses_a_certified_entry_circle() -> None:
    w = _entry_clear()
    report = production_implies_certificate(w, range(len(w.circles)), 1)
    assert report.entry is Verdict.CERTIFIED_LAST16 and report.violations == ()
    with pytest.raises(EntryCircleNotRefusedError):
        report.require()
