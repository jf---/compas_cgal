import hashlib
import math
from dataclasses import replace

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _coverage_2
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.errors import CoverageTransitionError
from compas_cgal.adaptive.errors import IncompletePocketCoverageError
from compas_cgal.adaptive.errors import InvalidCoverageCertificateError
from compas_cgal.adaptive.errors import InvalidCoverageGeometryError
from compas_cgal.adaptive.errors import InvalidCoverageSweepError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock


def _array(points: tuple[tuple[float, float], ...]) -> np.ndarray:
    return np.asarray(tuple((x, y, 0.0) for x, y in points), dtype=np.float64)


def _ring(points: tuple[tuple[float, float], ...]) -> CanonicalRingV1:
    return CanonicalRingV1.build_outer(tuple(Point2[WorldXY].build(x, y) for x, y in points))


def _domain(
    points: tuple[tuple[float, float], ...],
    radius: float,
) -> ReachableDomain:
    return ReachableDomain.build(
        design_boundary=_ring(points),
        holes=(),
        tool_radius=ToolRadius.build(radius),
    )


def _native_coverage(
    points: tuple[tuple[float, float], ...],
    domain_radius: float,
    precleared_center: tuple[float, float] = (0.0, 0.0),
    precleared_radius: float = 1.0,
) -> _coverage_2.Coverage2:
    reachable = _coverage_2.ReachableDomain2(
        _array(points),
        [],
        domain_radius,
    ).reachable_material()
    return _coverage_2.Coverage2(
        reachable,
        precleared_center[0],
        precleared_center[1],
        precleared_radius,
    )


LARGE_SQUARE = ((-40.0, -40.0), (40.0, -40.0), (40.0, 40.0), (-40.0, 40.0))
COMPLETE_RECTANGLE = ((0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0))


def test_native_binding_publishes_task5_surface_and_named_errors() -> None:
    assert hasattr(_coverage_2, "CoverageSweepRecord2")
    assert hasattr(_coverage_2, "Coverage2")
    assert InvalidCoverageGeometryError is _coverage_2.InvalidCoverageGeometryError
    assert CoverageTransitionError is _coverage_2.CoverageTransitionError
    assert not hasattr(_coverage_2, "CoverageConstructionError")


def test_native_geometry_and_transition_errors_are_distinct_and_atomic() -> None:
    domain = _coverage_2.ReachableDomain2(
        _array(COMPLETE_RECTANGLE),
        [],
        1.0,
    )
    coverage = _coverage_2.Coverage2(
        domain.reachable_material(),
        2.0,
        2.0,
        1.0,
    )
    residual_before = coverage.residual()

    with pytest.raises(InvalidCoverageGeometryError, match="distinct endpoints"):
        coverage.add_segment_sweep(1.0, 1.0, 1.0, 1.0, 1.0)

    assert coverage.residual().exactly_equals(residual_before)
    assert tuple(coverage.sweep_records) == ()
    with pytest.raises(CoverageTransitionError, match="reachable-material"):
        _coverage_2.Coverage2(
            domain.design_region(),
            2.0,
            2.0,
            1.0,
        )


def test_exact_non_axis_capsule_membership_at_endpoints_sides_and_seams() -> None:
    coverage = _native_coverage(LARGE_SQUARE, 1.0, (30.0, 30.0), 1.0)

    coverage.add_segment_sweep(0.0, 0.0, 3.0, 4.0, 5.0)
    sweep = coverage.accumulated_sweeps()

    assert sweep.contains(0.0, 0.0)
    assert sweep.contains(3.0, 4.0)
    assert sweep.contains(-2.5, 5.0)
    assert sweep.contains(5.5, -1.0)
    assert not sweep.contains(np.nextafter(-2.5, -3.0), 5.0)
    assert not sweep.contains(np.nextafter(5.5, 6.0), -1.0)


def test_native_segment_record_rejects_any_different_binary64_input() -> None:
    coverage = _native_coverage(LARGE_SQUARE, 1.0, (30.0, 30.0), 1.0)

    record = coverage.add_segment_sweep(0.0, 0.0, 3.0, 4.0, 5.0)

    assert record.matches_exact_segment(0.0, 0.0, 3.0, 4.0, 5.0)
    assert not record.matches_exact_segment(
        0.0,
        0.0,
        np.nextafter(3.0, 4.0),
        4.0,
        5.0,
    )
    assert record.strategy_version == coverage.strategy_version
    assert tuple(coverage.sweep_records) == (record.structural_record,)


@pytest.mark.parametrize(
    ("phase", "tool_radius", "inside", "boundary", "outside"),
    (
        ((3.0, 4.0), 10.0, (0.0, 0.0), (9.0, 12.0), (9.0, np.nextafter(12.0, 13.0))),
        ((3.0, 4.0), 5.0, (0.0, 0.0), (6.0, 8.0), (6.0, np.nextafter(8.0, 9.0))),
        ((6.0, 8.0), 5.0, (6.0, 8.0), (9.0, 12.0), (9.0, np.nextafter(12.0, 13.0))),
    ),
)
def test_full_circle_sweep_handles_guide_radius_less_equal_greater_than_tool(
    phase: tuple[float, float],
    tool_radius: float,
    inside: tuple[float, float],
    boundary: tuple[float, float],
    outside: tuple[float, float],
) -> None:
    coverage = _native_coverage(LARGE_SQUARE, 1.0, (30.0, 30.0), 1.0)

    coverage.add_full_circle_sweep(0.0, 0.0, phase[0], phase[1], tool_radius)
    sweep = coverage.accumulated_sweeps()

    assert sweep.contains(*inside)
    assert sweep.contains(*boundary)
    assert not sweep.contains(*outside)


def test_annulus_membership_binds_inner_outer_guide_and_seam_boundaries() -> None:
    coverage = _native_coverage(LARGE_SQUARE, 1.0, (30.0, 30.0), 1.0)

    coverage.add_full_circle_sweep(0.0, 0.0, 6.0, 8.0, 5.0)
    sweep = coverage.accumulated_sweeps()

    assert sweep.contains(3.0, 4.0)
    assert sweep.contains(6.0, 8.0)
    assert sweep.contains(9.0, 12.0)
    assert sweep.contains(-6.0, -8.0)
    assert not sweep.contains(np.nextafter(3.0, 0.0), 4.0)
    assert not sweep.contains(9.0, np.nextafter(12.0, 13.0))


def test_native_full_circle_record_rejects_any_different_binary64_input() -> None:
    coverage = _native_coverage(LARGE_SQUARE, 1.0, (30.0, 30.0), 1.0)

    record = coverage.add_full_circle_sweep(0.0, 0.0, 3.0, 4.0, 5.0)

    assert record.matches_exact_full_circle(0.0, 0.0, 3.0, 4.0, 5.0)
    assert not record.matches_exact_full_circle(
        0.0,
        0.0,
        3.0,
        np.nextafter(4.0, 5.0),
        5.0,
    )
    assert record.strategy_version == coverage.strategy_version
    assert tuple(coverage.sweep_records) == (record.structural_record,)


def _complete_ledger() -> tuple[CoverageLedger, ExactCircleMotion]:
    domain = _domain(COMPLETE_RECTANGLE, 1.0)
    ledger = CoverageLedger.build(
        reachable_domain=domain,
        precleared_center=Point2[WorldXY].build(2.0, 2.0),
        precleared_radius=ToolRadius.build(1.0),
    )
    circle = ExactCircleMotion.build(
        center=Point2[WorldXY].build(2.0, 2.0),
        phase_vector=Vector2[WorldXY].build(1.0, 1.0),
        clockwise=False,
    )
    return ledger, circle


def test_exact_residual_is_nonempty_before_and_empty_after_complete_sweep() -> None:
    ledger, circle = _complete_ledger()

    assert not ledger.residual_is_empty()
    assert ledger.residual_component_count() == 1
    residual_records = ledger.residual_component_records
    assert len(residual_records) == 1

    ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert ledger.residual_is_empty()
    assert ledger.residual_component_count() == 0
    assert ledger.residual_component_records == ()
    assert ledger.certificate.exact_residual_empty


def test_removing_final_sweep_leaves_named_residual_component() -> None:
    incomplete, circle = _complete_ledger()
    complete, _ = _complete_ledger()
    complete.add_sweep(circle, ToolRadius.build(1.0))

    assert not incomplete.residual_is_empty()
    assert incomplete.residual_component_records
    assert all(record.startswith(b"residual-component-v1") for record in incomplete.residual_component_records)
    with pytest.raises(IncompletePocketCoverageError, match="residual component"):
        incomplete.require_complete()
    complete.require_complete()


def _order_fixture() -> tuple[ReachableDomain, ExactSegmentMotion, ExactCircleMotion]:
    domain = _domain(LARGE_SQUARE, 1.0)
    segment = ExactSegmentMotion.build(
        Point2[WorldXY].build(-3.0, 0.0),
        Point2[WorldXY].build(3.0, 0.0),
    )
    circle = ExactCircleMotion.build(
        center=Point2[WorldXY].build(0.0, 0.0),
        phase_vector=Vector2[WorldXY].build(3.0, 4.0),
        clockwise=True,
    )
    return domain, segment, circle


def _ledger(domain: ReachableDomain) -> CoverageLedger:
    return CoverageLedger.build(
        reachable_domain=domain,
        precleared_center=Point2[WorldXY].build(30.0, 30.0),
        precleared_radius=ToolRadius.build(1.0),
    )


def test_operation_order_and_duplicates_preserve_sets_but_bind_canonical_lineage() -> None:
    domain, segment, circle = _order_fixture()
    forward = _ledger(domain)
    repeated = _ledger(domain)
    reverse = _ledger(domain)
    duplicate = _ledger(domain)

    for ledger in (forward, repeated):
        ledger.add_sweep(segment, ToolRadius.build(1.0))
        ledger.add_sweep(circle, ToolRadius.build(1.0))
    reverse.add_sweep(circle, ToolRadius.build(1.0))
    reverse.add_sweep(segment, ToolRadius.build(1.0))
    duplicate.add_sweep(segment, ToolRadius.build(1.0))
    duplicate.add_sweep(circle, ToolRadius.build(1.0))
    duplicate.add_sweep(segment, ToolRadius.build(1.0))

    assert forward.residual.exactly_equals(repeated.residual)
    assert forward.residual.exactly_equals(reverse.residual)
    assert forward.residual.exactly_equals(duplicate.residual)
    assert forward.certificate.digest == repeated.certificate.digest
    assert forward.certificate.digest != reverse.certificate.digest
    assert forward.certificate.digest != duplicate.certificate.digest
    assert forward.lineage == repeated.lineage
    assert forward.lineage != reverse.lineage
    assert forward.lineage != duplicate.lineage


def test_clone_is_alias_safe_for_residual_and_lineage() -> None:
    ledger, circle = _complete_ledger()
    clone = ledger.clone()
    snapshot = ledger.residual
    certificate_digest = ledger.certificate.digest

    clone.add_sweep(circle, ToolRadius.build(1.0))

    assert ledger.residual.exactly_equals(snapshot)
    assert not ledger.residual_is_empty()
    assert ledger.certificate.digest == certificate_digest
    assert clone.residual_is_empty()
    assert ledger.lineage == ()
    assert len(clone.lineage) == 1


class _CircleRecordBoundary:
    def __init__(self, record: object, matches: bool) -> None:
        self.strategy_version = record.strategy_version
        self.structural_record = record.structural_record
        self._matches = matches

    def matches_exact_full_circle(
        self,
        cx: float,
        cy: float,
        phase_x: float,
        phase_y: float,
        tool_radius: float,
    ) -> bool:
        return self._matches


class _CircleTrialBoundary:
    def __init__(
        self,
        native: _coverage_2.Coverage2,
        *,
        matches: bool,
        corrupt_history: bool,
        corrupt_certificate: bool,
    ) -> None:
        self._native = native
        self._matches = matches
        self._corrupt_history = corrupt_history
        self._corrupt_certificate = corrupt_certificate

    def add_full_circle_sweep(
        self,
        cx: float,
        cy: float,
        phase_x: float,
        phase_y: float,
        tool_radius: float,
    ) -> _CircleRecordBoundary:
        return _CircleRecordBoundary(
            self._native.add_full_circle_sweep(
                cx,
                cy,
                phase_x,
                phase_y,
                tool_radius,
            ),
            self._matches,
        )

    @property
    def strategy_version(self) -> bytes:
        return self._native.strategy_version

    @property
    def sweep_records(self) -> tuple[bytes, ...]:
        records = tuple(self._native.sweep_records)
        return records + ((b"divergent-native-history",) if self._corrupt_history else ())

    @property
    def residual_component_records(self) -> object:
        if self._corrupt_certificate:
            return (b"",)
        return self._native.residual_component_records

    def residual_component_count(self) -> int:
        return self._native.residual_component_count()

    def exact_residual_relation(self) -> bool:
        return self._native.exact_residual_relation()

    def residual_is_empty(self) -> bool:
        return self._native.residual_is_empty()


class _NativeBoundary:
    def __init__(
        self,
        native: _coverage_2.Coverage2,
        *,
        matches: bool = True,
        corrupt_history: bool = False,
        corrupt_certificate: bool = False,
    ) -> None:
        self._native = native
        self._matches = matches
        self._corrupt_history = corrupt_history
        self._corrupt_certificate = corrupt_certificate

    def clone(self) -> _CircleTrialBoundary:
        return _CircleTrialBoundary(
            self._native.clone(),
            matches=self._matches,
            corrupt_history=self._corrupt_history,
            corrupt_certificate=self._corrupt_certificate,
        )

    @property
    def residual_component_records(self) -> object:
        return self._native.residual_component_records

    def residual_component_count(self) -> int:
        return self._native.residual_component_count()

    def exact_residual_relation(self) -> bool:
        return self._native.exact_residual_relation()

    def residual_is_empty(self) -> bool:
        return self._native.residual_is_empty()

    @property
    def strategy_version(self) -> bytes:
        return self._native.strategy_version


def test_rejected_native_matcher_keeps_python_ledger_atomic() -> None:
    ledger, circle = _complete_ledger()
    certificate_digest = ledger.certificate.digest
    ledger._native = _NativeBoundary(ledger._native, matches=False)

    with pytest.raises(InvalidCoverageSweepError, match="bind the exact motion inputs"):
        ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert ledger.lineage == ()
    assert ledger.certificate.digest == certificate_digest
    assert not ledger.residual_is_empty()


def test_divergent_native_history_keeps_python_ledger_atomic() -> None:
    ledger, circle = _complete_ledger()
    certificate_digest = ledger.certificate.digest
    ledger._native = _NativeBoundary(
        ledger._native,
        corrupt_history=True,
    )

    with pytest.raises(InvalidCoverageSweepError, match="history diverged"):
        ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert ledger.lineage == ()
    assert ledger.certificate.digest == certificate_digest
    assert not ledger.residual_is_empty()


def test_invalid_trial_certificate_keeps_python_ledger_atomic() -> None:
    ledger, circle = _complete_ledger()
    certificate_digest = ledger.certificate.digest
    ledger._native = _NativeBoundary(
        ledger._native,
        corrupt_certificate=True,
    )

    with pytest.raises(
        InvalidCoverageCertificateError,
        match="nonempty exact bytes",
    ):
        ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert ledger.lineage == ()
    assert ledger.certificate.digest == certificate_digest
    assert not ledger.residual_is_empty()


def test_lineage_binds_each_parent_digest_in_order() -> None:
    domain, segment, circle = _order_fixture()
    ledger = _ledger(domain)

    first = ledger.add_sweep(segment, ToolRadius.build(1.0))
    second = ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert first.parent_lineage == ()
    assert second.parent_lineage == (first.digest,)
    assert ledger.certificate.ordered_sweep_records == (
        first.canonical_bytes,
        second.canonical_bytes,
    )


def test_under_cover_stock_remains_while_full_sweep_ledger_proves_coverage() -> None:
    ledger, circle = _complete_ledger()
    stock = Stock(Polygon(tuple((x, y, 0.0) for x, y in COMPLETE_RECTANGLE)))
    stock.subtract_disk(2.0, 2.0, 1.0)
    stock.raw.subtract_exact_full_circle(
        2.0,
        2.0,
        1.0,
        1.0,
        False,
        1.0,
        0.75,
        4096,
    )

    ledger.add_sweep(circle, ToolRadius.build(1.0))

    assert not stock.is_empty()
    assert ledger.residual_is_empty()


def test_rounded_guide_radius_fails_exact_outer_boundary_fixture() -> None:
    exact = _native_coverage(COMPLETE_RECTANGLE, 1.0, (2.0, 2.0), 1.0)
    rounded = exact.clone()

    exact.add_full_circle_sweep(2.0, 2.0, 1.0, 1.0, 1.0)
    rounded.add_full_circle_sweep(
        2.0,
        2.0,
        math.nextafter(1.0, 0.0),
        1.0,
        1.0,
    )

    assert exact.residual_is_empty()
    assert not rounded.residual_is_empty()
    assert rounded.residual_component_count() > 0


def test_coverage_certificate_binds_motions_m_r_and_separate_unreachable_digest() -> None:
    ledger, circle = _complete_ledger()
    initial = ledger.certificate
    ledger.add_sweep(circle, ToolRadius.build(1.0))
    complete = ledger.certificate

    assert initial.reachable_material_digest == complete.reachable_material_digest
    assert initial.unreachable_residual_digest == complete.unreachable_residual_digest
    assert initial.ordered_sweep_records == ()
    assert len(complete.ordered_sweep_records) == 1
    assert complete.ordered_sweep_records[0] == ledger.lineage[0].canonical_bytes
    assert complete.digest == hashlib.sha256(complete.canonical_bytes).digest()
    assert complete.exact_residual_empty


def test_invalid_coverage_certificate_cannot_mint_identity() -> None:
    ledger, _ = _complete_ledger()
    invalid = replace(
        ledger.certificate,
        exact_residual_relation=False,
    )

    with pytest.raises(InvalidCoverageCertificateError, match="relation"):
        _ = invalid.digest
