import inspect
from dataclasses import FrozenInstanceError
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import ExactCenterParameterV1
from compas_cgal.adaptive.errors import ExactDepletionCenterLimitError
from compas_cgal.adaptive.errors import InvalidDepletionTraceError
from compas_cgal.adaptive.errors import InvalidDepletionWitnessError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import DepletionPolicy
import compas_cgal.adaptive.stock_area as stock_area_module
from compas_cgal.adaptive.stock_area import DepletionWitness
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock

SQUARE = np.array(
    [[-8.0, -8.0, 0.0], [8.0, -8.0, 0.0], [8.0, 8.0, 0.0], [-8.0, 8.0, 0.0]],
    dtype=np.float64,
)


def _raw_stock() -> _stock_2.Stock2:
    return _stock_2.Stock2(SQUARE, [])


def _stock() -> Stock:
    return Stock(Polygon([[-8.0, -8.0, 0.0], [8.0, -8.0, 0.0], [8.0, 8.0, 0.0], [-8.0, 8.0, 0.0]]))


def _policy(chord_bound: float = 0.5, center_count_limit: int = 4096) -> DepletionPolicy:
    return DepletionPolicy.build(
        chord_bound=ChordBound.build(chord_bound),
        center_count_limit=center_count_limit,
    )


def _segment() -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(-3.0, -1.0),
        Point2[WorldXY].build(3.0, 7.0),
    )


def _circle(*, clockwise: bool) -> ExactCircleMotion:
    return ExactCircleMotion.build(
        Point2[WorldXY].build(0.0, 0.0),
        Vector2[WorldXY].build(3.0, 4.0),
        clockwise,
    )


def test_exact_segment_trace_proves_every_construction_invariant() -> None:
    stock = _raw_stock()

    trace = stock.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 0.5, 4096)

    assert trace.center_count == len(trace.center_parameters)
    assert trace.center_parameters[0] == (-1, 0, trace.center_parameters[0][2])
    assert trace.center_parameters[-1] == (-1, trace.center_parameters[-1][2], trace.center_parameters[-1][2])
    assert trace.exact_incidence
    assert trace.exact_parameters_in_range
    assert trace.exact_anchors_present
    assert trace.exact_removal_radius_valid
    assert trace.exact_chord_bound_holds
    assert trace.matches_exact_inputs(0.75, 0.5, 4096)
    assert not trace.matches_exact_inputs(0.5, 0.5, 4096)
    assert not trace.matches_exact_inputs(0.75, 0.25, 4096)
    assert not trace.matches_exact_inputs(0.75, 0.5, trace.center_count - 1)
    assert not hasattr(trace, "max_chord")
    assert not hasattr(trace, "removal_radius")
    assert not trace.cyclic


def test_exact_circle_trace_proves_cardinal_anchors_and_cyclic_seam() -> None:
    stock = _raw_stock()

    trace = stock.subtract_exact_full_circle(0.0, 0.0, 3.0, 4.0, False, 0.75, 0.5, 4096)

    assert trace.center_count == len(trace.center_parameters)
    assert {(chart, numerator) for chart, numerator, _ in trace.center_parameters} >= {
        (0, 0),
        (1, 0),
        (2, 0),
        (3, 0),
    }
    assert trace.exact_incidence
    assert trace.exact_parameters_in_range
    assert trace.exact_anchors_present
    assert trace.exact_removal_radius_valid
    assert trace.exact_chord_bound_holds
    assert trace.exact_seam_chord_bound_holds
    assert trace.cyclic


@pytest.mark.parametrize(
    ("guide_radius", "tool_radius"),
    [(1.0, 2.0), (2.0, 2.0), (5.0, 0.75)],
)
def test_exact_circle_union_is_exact_subset_of_representable_sweep(
    guide_radius: float,
    tool_radius: float,
) -> None:
    assert _stock_2.exact_full_circle_undercover_holds(
        0.0,
        0.0,
        guide_radius,
        0.0,
        guide_radius,
        tool_radius,
        0.5,
        4096,
    )


def test_exact_segment_union_is_exact_subset_of_representable_capsule() -> None:
    assert _stock_2.exact_segment_undercover_holds(
        -3.0,
        -1.0,
        3.0,
        7.0,
        10.0,
        0.75,
        0.5,
        4096,
    )


def test_direct_exact_induction_holds_for_segment_and_circle() -> None:
    initial = _raw_stock()

    assert _stock_2.exact_segment_induction_holds(
        initial,
        -3.0,
        -1.0,
        3.0,
        7.0,
        10.0,
        0.75,
        0.5,
        4096,
    )
    assert _stock_2.exact_full_circle_induction_holds(
        initial,
        0.0,
        0.0,
        5.0,
        0.0,
        5.0,
        0.75,
        0.5,
        4096,
    )


def test_minimum_density_and_repeated_overlap_preserve_exact_induction() -> None:
    minimum_segment = _raw_stock()
    segment_trace = minimum_segment.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 10.0, 2)
    minimum_circle = _raw_stock()
    circle_trace = minimum_circle.subtract_exact_full_circle(0.0, 0.0, 1.0, 0.0, False, 0.25, 2.0, 4)

    assert segment_trace.center_count == 2
    assert circle_trace.center_count == 4

    once = _raw_stock()
    twice = once.clone()
    once.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 0.5, 4096)
    twice.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 0.5, 4096)
    twice.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 0.5, 4096)
    assert twice.exactly_equals(once)

    assert _stock_2.exact_segment_induction_holds(
        once,
        -2.5,
        -1.0,
        3.5,
        7.0,
        10.0,
        0.75,
        0.5,
        4096,
    )


def test_rounded_constructions_fail_native_exact_incidence_predicates() -> None:
    interpolation = 1.0 / 3.0
    x = 0.1 + interpolation * (0.3 - 0.1)
    y = 0.2 + interpolation * (0.7 - 0.2)
    assert not _stock_2.exact_segment_point_is_incident(0.1, 0.2, 0.3, 0.7, x, y)

    chart_parameter = 0.3
    denominator = 1.0 + chart_parameter * chart_parameter
    cosine = (1.0 - chart_parameter * chart_parameter) / denominator
    sine = 2.0 * chart_parameter / denominator
    assert not _stock_2.exact_circle_point_is_incident(0.0, 0.0, 3.0, 4.0, 3.0 * cosine - 4.0 * sine, 3.0 * sine + 4.0 * cosine)


def test_clone_is_an_exact_independent_deep_copy() -> None:
    original = _raw_stock()
    clone = original.clone()

    assert original.exactly_equals(clone)
    clone.subtract_disk(0.0, 0.0, 1.0)

    assert original.contains(0.0, 0.0)
    assert not clone.contains(0.0, 0.0)
    assert clone.is_subset_of(original)
    assert not original.exactly_equals(clone)


def test_center_limit_rejection_precedes_native_mutation() -> None:
    stock = _raw_stock()
    snapshot = stock.clone()

    with pytest.raises(ExactDepletionCenterLimitError):
        stock.subtract_exact_full_circle(0.0, 0.0, 3.0, 4.0, False, 0.75, 0.03125, 4)

    assert stock.exactly_equals(snapshot)


def test_tighter_chord_bound_only_depletes_more_modeled_stock() -> None:
    coarse = _raw_stock()
    fine = coarse.clone()

    coarse_trace = coarse.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 1.0, 4096)
    fine_trace = fine.subtract_exact_segment(-3.0, -1.0, 3.0, 7.0, 0.75, 0.25, 4096)

    assert fine_trace.center_count > coarse_trace.center_count
    assert fine.is_subset_of(coarse)


def test_cw_and_ccw_share_removal_set_but_not_ordered_witness() -> None:
    ccw = Stock2Area.build(_stock())
    cw = ccw.fork()

    ccw_witness = ccw.deplete(_circle(clockwise=False), ToolRadius.build(0.75), _policy())
    cw_witness = cw.deplete(_circle(clockwise=True), ToolRadius.build(0.75), _policy())

    assert ccw.raw.exactly_equals(cw.raw)
    assert set(ccw_witness.center_parameters) == set(cw_witness.center_parameters)
    assert ccw_witness.center_parameters != cw_witness.center_parameters
    assert ccw_witness.digest != cw_witness.digest


def test_public_area_construction_deep_copies_mutable_native_input() -> None:
    source = _raw_stock()
    area = Stock2Area(source, ())
    snapshot = area.raw.clone()

    source.subtract_disk(0.0, 0.0, 1.0)

    assert area.raw.exactly_equals(snapshot)
    assert area.lineage == ()


def test_public_raw_access_returns_an_independent_snapshot() -> None:
    area = Stock2Area.build(_stock())
    snapshot = area.raw.clone()
    outward = area.raw

    outward.subtract_disk(0.0, 0.0, 1.0)

    assert area.raw.exactly_equals(snapshot)
    assert area.lineage == ()


def test_exact_trace_mismatch_rejects_trial_without_commit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    area = Stock2Area.build(_stock())
    snapshot = area.raw

    def reject_mismatched_trace(
        trace: _stock_2.DepletionTrace,
        motion: ExactSegmentMotion | ExactCircleMotion,
        tool_radius: ToolRadius,
        policy: DepletionPolicy,
    ) -> tuple[ExactCenterParameterV1, ...]:
        del motion
        if not trace.matches_exact_inputs(
            tool_radius.value,
            policy.chord_bound.value / 2.0,
            policy.center_count_limit,
        ):
            raise InvalidDepletionTraceError("injected exact trace mismatch.")
        raise AssertionError("fault injection did not create an exact mismatch.")

    monkeypatch.setattr(stock_area_module, "_validated_center_parameters", reject_mismatched_trace)

    with pytest.raises(InvalidDepletionTraceError, match="exact trace mismatch"):
        area.deplete(_segment(), ToolRadius.build(0.75), _policy())

    assert area.raw.exactly_equals(snapshot)
    assert area.lineage == ()


def test_witness_owns_full_inputs_and_is_immutable() -> None:
    area = Stock2Area.build(_stock())
    motion = _segment()
    policy = _policy()
    radius = ToolRadius.build(0.75)

    witness = area.deplete(motion, radius, policy)

    assert witness.motion is motion
    assert witness.policy is policy
    assert witness.tool_radius is radius
    assert all(type(parameter) is ExactCenterParameterV1 for parameter in witness.center_parameters)
    assert type(witness.native_strategy_version) is bytes
    assert witness.native_strategy_version
    assert witness.parent_lineage == ()
    assert area.lineage == (witness,)
    with pytest.raises(FrozenInstanceError):
        witness.parent_lineage = ()  # type: ignore[misc]


def _rebuild_witness(
    witness: DepletionWitness,
    *,
    center_parameters: tuple[ExactCenterParameterV1, ...] | None = None,
    policy: DepletionPolicy | None = None,
    native_strategy_version: bytes | None = None,
) -> DepletionWitness:
    return DepletionWitness(
        witness.motion,
        policy if policy is not None else witness.policy,
        witness.tool_radius,
        center_parameters if center_parameters is not None else witness.center_parameters,
        native_strategy_version if native_strategy_version is not None else witness.native_strategy_version,
        witness.parent_lineage,
    )


def test_direct_witness_rejects_omitted_or_reordered_segment_centers() -> None:
    witness = Stock2Area.build(_stock()).deplete(_segment(), ToolRadius.build(0.75), _policy())
    reordered = list(witness.center_parameters)
    reordered[0], reordered[1] = reordered[1], reordered[0]

    with pytest.raises(InvalidDepletionWitnessError, match="segment"):
        _rebuild_witness(witness, center_parameters=witness.center_parameters[:-1])
    with pytest.raises(InvalidDepletionWitnessError, match="segment"):
        _rebuild_witness(witness, center_parameters=tuple(reordered))


def test_direct_witness_rejects_policy_limit_and_strategy_mutations() -> None:
    witness = Stock2Area.build(_stock()).deplete(_segment(), ToolRadius.build(0.75), _policy())
    too_small = _policy(center_count_limit=len(witness.center_parameters) - 1)

    with pytest.raises(InvalidDepletionWitnessError, match="center-count"):
        _rebuild_witness(witness, policy=too_small)
    with pytest.raises(InvalidDepletionWitnessError, match="strategy"):
        _rebuild_witness(witness, native_strategy_version=b"unsupported-exact-strategy-v0")


def test_direct_witness_rejects_circle_chart_denominator_and_order_mutations() -> None:
    witness = Stock2Area.build(_stock()).deplete(_circle(clockwise=False), ToolRadius.build(0.75), _policy())
    wrong_chart = (replace(witness.center_parameters[0], chart=1),) + witness.center_parameters[1:]
    wrong_denominator = (replace(witness.center_parameters[0], denominator=witness.center_parameters[0].denominator * 2),) + witness.center_parameters[1:]
    reordered = (witness.center_parameters[1], witness.center_parameters[0]) + witness.center_parameters[2:]

    for mutation in (wrong_chart, wrong_denominator, reordered):
        with pytest.raises(InvalidDepletionWitnessError, match="circle"):
            _rebuild_witness(witness, center_parameters=mutation)


def test_direct_valid_witness_reconstruction_preserves_identity() -> None:
    witness = Stock2Area.build(_stock()).deplete(_circle(clockwise=True), ToolRadius.build(0.75), _policy())

    rebuilt = _rebuild_witness(witness)

    assert rebuilt == witness
    assert rebuilt.digest == witness.digest


def test_fork_retains_lineage_but_owns_independent_raw_stock() -> None:
    area = Stock2Area.build(_stock())
    first = area.deplete(_segment(), ToolRadius.build(0.5), _policy())
    fork = area.fork()
    snapshot = area.raw.clone()

    fork.deplete(_circle(clockwise=False), ToolRadius.build(0.5), _policy())

    assert area.raw.exactly_equals(snapshot)
    assert area.lineage == (first,)
    assert fork.lineage[:1] == area.lineage
    assert len(fork.lineage) == 2


def test_failed_area_trial_preserves_stock_and_lineage_atomically() -> None:
    area = Stock2Area.build(_stock())
    snapshot = area.raw.clone()

    with pytest.raises(ExactDepletionCenterLimitError):
        area.deplete(_circle(clockwise=False), ToolRadius.build(0.75), _policy(0.03125, 4))

    assert area.raw.exactly_equals(snapshot)
    assert area.lineage == ()


def test_structural_parameter_canonicalization_retains_exact_odd_fraction() -> None:
    parameter = ExactCenterParameterV1.build(chart=2, numerator=1, denominator=3)

    assert parameter.fraction == Fraction(1, 3)
    assert parameter.canonical_bytes


def test_public_stock_area_api_has_google_style_contract_docs() -> None:
    assert inspect.getdoc(DepletionWitness)
    assert inspect.getdoc(Stock2Area)
    assert "Returns:" in inspect.getdoc(Stock2Area.build)
    assert "Returns:" in inspect.getdoc(Stock2Area.fork)
    deplete_docs = inspect.getdoc(Stock2Area.deplete)
    assert "Returns:" in deplete_docs
    assert "Raises:" in deplete_docs
