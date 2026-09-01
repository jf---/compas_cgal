from __future__ import annotations

import math
import sys
from fractions import Fraction

import pytest

import benchmarks.held_reference_certification as certification
import benchmarks.held_reference_geometry as geometry
from benchmarks.errors import DisconnectedPublishedBoundaryError
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import InvalidReferenceProjectionError
from benchmarks.errors import InvalidReferenceReconstructionError
from benchmarks.held_reference_certification import BINARY64_UNIT_ROUNDOFF
from benchmarks.held_reference_certification import BIARC_POLYNOMIAL_OPERATION_COUNT
from benchmarks.held_reference_certification import certify_biarc_root
from benchmarks.held_reference_geometry import MillimetresPerPdfPoint
from benchmarks.held_reference_geometry import NORMALIZED_ROUNDOFF
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_geometry import SourceCubic
from benchmarks.held_reference_geometry import SourceLine
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import _certified_line
from benchmarks.held_reference_geometry import project_boundary
from benchmarks.held_reference_geometry import reconstruct_cubic_certified
from benchmarks.held_reference_geometry import reconstruct_cubic
from benchmarks.held_reference_geometry import reconstruct_source_path
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _world(x: float, y: float) -> Point2[WorldXY]:
    return Point2[WorldXY].build(x, y)


def _identity_transform() -> SourceToWorld:
    return SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 0.0),
        world_origin=_world(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0),
    )


def _quarter_circle_source() -> SourceCubic:
    control = 4.0 * (math.sqrt(2.0) - 1.0) / 3.0
    return SourceCubic.build(
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, control),
        PdfPoint2.build(control, 1.0),
        PdfPoint2.build(0.0, 1.0),
    )


def _circle_quarter(start_index: int, end_index: int, radius: float = 1.0) -> SourceCubic:
    points = (
        (radius, 0.0),
        (0.0, radius),
        (-radius, 0.0),
        (0.0, -radius),
    )
    start = points[start_index]
    end = points[end_index]
    sweep_sign = 1.0 if (end_index - start_index) % 4 == 1 else -1.0
    tangent_start = sweep_sign * -start[1] / radius, sweep_sign * start[0] / radius
    tangent_end = sweep_sign * -end[1] / radius, sweep_sign * end[0] / radius
    handle = radius * 4.0 * (math.sqrt(2.0) - 1.0) / 3.0
    return SourceCubic.build(
        PdfPoint2.build(*start),
        PdfPoint2.build(start[0] + handle * tangent_start[0], start[1] + handle * tangent_start[1]),
        PdfPoint2.build(end[0] - handle * tangent_end[0], end[1] - handle * tangent_end[1]),
        PdfPoint2.build(*end),
    )


def _near_circle_second_quarter() -> SourceCubic:
    handle = 4.0 * (math.sqrt(2.0) - 1.0) / 3.0
    return SourceCubic.build(
        PdfPoint2.build(0.0, 1.0),
        PdfPoint2.build(-handle, 1.0),
        PdfPoint2.build(-1.01, handle),
        PdfPoint2.build(-1.01, 0.0),
    )


def _non_circular_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, 2.0),
        PdfPoint2.build(2.0, 2.0),
    )


def _asymmetric_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(2.0, 1.0),
        PdfPoint2.build(2.0, 3.0),
    )


def _semicircle_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(-1.0, 0.0),
        PdfPoint2.build(-1.0, 4.0 / 3.0),
        PdfPoint2.build(1.0, 4.0 / 3.0),
        PdfPoint2.build(1.0, 0.0),
    )


def _angular_reversal_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(1.0, 5.0),
        PdfPoint2.build(5.0, 1.0),
        PdfPoint2.build(0.0, 1.0),
    )


def _monstera_source_125() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(311.593887480182, 423.7029828915639),
        PdfPoint2.build(311.593887480182, 423.58579534974297),
        PdfPoint2.build(311.601700060447, 423.46860780792196),
        PdfPoint2.build(311.61341922138297, 423.35142026610094),
    )


def _monstera_transform() -> SourceToWorld:
    return SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 0.0),
        world_origin=_world(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0 / 2.826174326591),
    )


def _figure5_proof_hard_source() -> SourceCubic:
    return SourceCubic.build(
        PdfPoint2.build(170.574243138248, 93.339818109552),
        PdfPoint2.build(152.20705251096499, 116.695291267899),
        PdfPoint2.build(118.38673357500198, 120.738260606424),
        PdfPoint2.build(95.031260416655, 102.37106950249799),
    )


def _figure5_transform() -> SourceToWorld:
    return SourceToWorld.build(
        source_origin=PdfPoint2.build(261.69535427191397, 146.761702919259),
        world_origin=_world(0.0, 0.0),
        scale=MillimetresPerPdfPoint(1.0 / 3.226562815407),
        reflect_source_y=True,
    )


def _unit_tangent(arc: ReferenceArc, *, at_end: bool) -> tuple[float, float]:
    point = arc.end if at_end else arc.start
    radius_x = float(point.x) - float(arc.centre.x)
    radius_y = float(point.y) - float(arc.centre.y)
    direction = 1.0 if float(arc.sweep) > 0.0 else -1.0
    length = math.hypot(radius_x, radius_y)
    return direction * -radius_y / length, direction * radius_x / length


def _vector_residual(left: tuple[float, float], right: tuple[float, float]) -> float:
    return math.hypot(left[0] - right[0], left[1] - right[1])


def _rounded_rectangle_boundary() -> ReferenceBoundary:
    primitives = (
        ReferenceLine.build(_world(1.0, 0.0), _world(3.0, 0.0)),
        ReferenceArc.build(_world(3.0, 0.0), _world(4.0, 1.0), _world(3.0, 1.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(4.0, 1.0), _world(4.0, 3.0)),
        ReferenceArc.build(_world(4.0, 3.0), _world(3.0, 4.0), _world(3.0, 3.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(3.0, 4.0), _world(1.0, 4.0)),
        ReferenceArc.build(_world(1.0, 4.0), _world(0.0, 3.0), _world(1.0, 3.0), Radian(math.pi / 2.0)),
        ReferenceLine.build(_world(0.0, 3.0), _world(0.0, 1.0)),
        ReferenceArc.build(_world(0.0, 1.0), _world(1.0, 0.0), _world(1.0, 1.0), Radian(math.pi / 2.0)),
    )
    return ReferenceBoundary.build(primitives, ToolRadius.build(1.0), Millimetre(0.1))


def _unit_circle_boundary(offset: float = 0.0) -> ReferenceBoundary:
    centre = _world(offset, offset)
    cardinal_points = (
        _world(offset + 1.0, offset),
        _world(offset, offset + 1.0),
        _world(offset - 1.0, offset),
        _world(offset, offset - 1.0),
    )
    arcs = tuple(ReferenceArc.build(start, end, centre, Radian(math.pi / 2.0)) for start, end in zip(cardinal_points, (*cardinal_points[1:], cardinal_points[0])))
    return ReferenceBoundary.build(arcs, ToolRadius.build(1.0), Millimetre(0.1))


def _stored_semantics_counterexample(order: tuple[int, ...], sweep: float) -> ReferenceBoundary:
    centre = (4091241540306.3516, 3228197555522.9883)
    construction_radius = 0.5
    phase = 5.972552524945925
    points = tuple(
        _world(
            centre[0] + construction_radius * math.cos(phase + index * math.pi / 2.0),
            centre[1] + construction_radius * math.sin(phase + index * math.pi / 2.0),
        )
        for index in range(4)
    )
    world_centre = _world(*centre)
    arcs = tuple(ReferenceArc.build(points[start], points[end], world_centre, Radian(sweep)) for start, end in zip(order, (*order[1:], order[0])))
    return ReferenceBoundary.build(arcs, ToolRadius.build(1.0), Millimetre(0.1))


def _point_segment_distance(
    point: tuple[float, float],
    start: tuple[float, float],
    end: tuple[float, float],
) -> float:
    segment = end[0] - start[0], end[1] - start[1]
    parameter = ((point[0] - start[0]) * segment[0] + (point[1] - start[1]) * segment[1]) / (segment[0] ** 2 + segment[1] ** 2)
    parameter = min(1.0, max(0.0, parameter))
    closest = start[0] + parameter * segment[0], start[1] + parameter * segment[1]
    return math.hypot(point[0] - closest[0], point[1] - closest[1])


def test_pdf_point_rejects_non_finite_coordinate() -> None:
    with pytest.raises(InvalidPublishedPrimitiveError):
        PdfPoint2.build(float("nan"), 0.0)


def test_source_line_rejects_identical_endpoints() -> None:
    point = PdfPoint2.build(2.0, 3.0)
    with pytest.raises(InvalidPublishedPrimitiveError):
        SourceLine.build(point, point)


def test_reference_boundary_rejects_disconnected_primitives() -> None:
    first = ReferenceLine.build(_world(0.0, 0.0), _world(1.0, 0.0))
    second = ReferenceLine.build(_world(2.0, 0.0), _world(0.0, 0.0))

    with pytest.raises(DisconnectedPublishedBoundaryError):
        ReferenceBoundary.build(
            (first, second),
            ToolRadius.build(1.0),
            Millimetre(0.1),
        )


@pytest.mark.parametrize("offset", (0.0, 2.0**30))
def test_reference_arc_rejects_local_radius_defect_independent_of_translation(offset: float) -> None:
    radius = 2.0**-30
    relative_defect = math.sqrt(sys.float_info.epsilon)

    with pytest.raises(InvalidPublishedPrimitiveError):
        ReferenceArc.build(
            _world(offset, radius),
            _world(offset, -radius * (1.0 + relative_defect)),
            _world(offset, 0.0),
            Radian(math.pi),
        )


def test_quarter_circle_cubic_recovers_one_arc() -> None:
    arcs = reconstruct_cubic(
        _quarter_circle_source(),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(arcs) == 1
    assert float(arcs[0].centre.x) == pytest.approx(0.0)
    assert float(arcs[0].centre.y) == pytest.approx(0.0)
    assert float(arcs[0].sweep) == pytest.approx(math.pi / 2.0)


def test_non_circular_cubic_splits_into_g1_arcs() -> None:
    arcs = reconstruct_cubic(
        _non_circular_source(),
        _identity_transform(),
        Millimetre(0.01),
    )

    assert len(arcs) >= 2
    for left, right in zip(arcs, arcs[1:]):
        assert left.end == right.start
        assert _unit_tangent(left, at_end=True) == pytest.approx(_unit_tangent(right, at_end=False))


def test_asymmetric_cubic_reconstructs_without_unequal_radius_failure() -> None:
    arcs = reconstruct_cubic(
        _asymmetric_source(),
        _identity_transform(),
        Millimetre(0.01),
    )

    assert arcs[0].start == _world(0.0, 0.0)
    assert arcs[-1].end == _world(2.0, 3.0)
    for left, right in zip(arcs, arcs[1:]):
        assert _unit_tangent(left, at_end=True) == pytest.approx(_unit_tangent(right, at_end=False))


def test_antiparallel_endpoint_tangents_recover_one_semicircle() -> None:
    arcs = reconstruct_cubic(
        _semicircle_source(),
        _identity_transform(),
        Millimetre(0.02),
    )

    assert len(arcs) == 1
    assert float(arcs[0].centre.x) == pytest.approx(0.0)
    assert float(arcs[0].centre.y) == pytest.approx(0.0)
    assert float(arcs[0].sweep) == pytest.approx(-math.pi)


def test_near_parallel_biarc_preserves_both_authored_endpoint_tangents() -> None:
    tangent_x = math.sqrt(32.0 * sys.float_info.epsilon)
    start_tangent = (tangent_x, math.sqrt(1.0 - tangent_x**2))
    end_tangent = (0.0, 1.0)
    source = SourceCubic.build(
        PdfPoint2.build(-1.0, 0.0),
        PdfPoint2.build(-1.0 + start_tangent[0], start_tangent[1]),
        PdfPoint2.build(1.0 - end_tangent[0], -end_tangent[1]),
        PdfPoint2.build(1.0, 0.0),
    )

    arcs = reconstruct_cubic(source, _identity_transform(), Millimetre(2.0))

    assert _vector_residual(_unit_tangent(arcs[0], at_end=False), start_tangent) <= NORMALIZED_ROUNDOFF
    assert _vector_residual(_unit_tangent(arcs[-1], at_end=True), end_tangent) <= NORMALIZED_ROUNDOFF
    for left, right in zip(arcs, arcs[1:]):
        assert _vector_residual(_unit_tangent(left, at_end=True), _unit_tangent(right, at_end=False)) <= NORMALIZED_ROUNDOFF


def test_monstera_source_125_recovers_root_biarc_without_subdivision() -> None:
    primitives = reconstruct_cubic(
        _monstera_source_125(),
        _monstera_transform(),
        Millimetre(0.10280275256426046),
    )

    assert len(primitives) == 2
    assert all(isinstance(primitive, ReferenceArc) for primitive in primitives)


def test_figure5_proof_hard_source_closes_with_bounded_certificate_nodes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    node_limit = 2_048
    calls = 0
    original = geometry.biarc_correspondence_node_bound

    def counted(*args: object, **kwargs: object) -> float | None:
        nonlocal calls
        calls += 1
        assert calls <= node_limit
        return original(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(geometry, "biarc_correspondence_node_bound", counted)

    reconstruction = reconstruct_cubic_certified(
        _figure5_proof_hard_source(),
        _figure5_transform(),
        Millimetre(0.07386234629061081),
    )

    assert reconstruction.primitives
    assert float(reconstruction.deviation_upper_bound) <= 0.07386234629061081
    assert calls <= node_limit


def test_figure5_biarc_children_do_not_invent_rational_breakpoint_witnesses(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _figure5_proof_hard_source()
    transform = _figure5_transform()
    limit = Millimetre(0.07386234629061081)
    witnessed_intervals: list[tuple[Fraction, Fraction]] = []
    merge_witnesses: list[tuple[object | None, object | None]] = []
    original = geometry._merge_arc_entries

    def capture(first: object, second: object, deviation_limit: float) -> object:
        merge_witnesses.append(
            (
                first.source_witnesses,  # type: ignore[attr-defined]
                second.source_witnesses,  # type: ignore[attr-defined]
            )
        )
        for entry in (first, second):
            witnesses = entry.source_witnesses  # type: ignore[attr-defined]
            if witnesses is not None:
                witnessed_intervals.extend((start, end) for _, start, end in witnesses)
        return original(first, second, deviation_limit)  # type: ignore[arg-type]

    monkeypatch.setattr(geometry, "_merge_arc_entries", capture)

    merged = reconstruct_source_path((source,), transform, limit)

    assert (None, None) in merge_witnesses
    assert all(start == 0 and end == 1 for start, end in witnessed_intervals)
    assert float(merged.deviation_upper_bound) <= float(limit)


def test_partial_source_witness_cannot_authorize_adjacent_arc_merge(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    midpoint = math.sqrt(0.5)
    first_arc = ReferenceArc.build(
        _world(1.0, 0.0),
        _world(midpoint, midpoint),
        _world(0.0, 0.0),
        Radian(math.pi / 4.0),
    )
    second_arc = ReferenceArc.build(
        _world(midpoint, midpoint),
        _world(0.0, 1.0),
        _world(0.0, 0.0),
        Radian(math.pi / 4.0),
    )
    source = _quarter_circle_source()
    controls = tuple(geometry._point_xy(_identity_transform().point(point)) for point in (source.start, source.control1, source.control2, source.end))
    first = geometry._MergeEntry(first_arc, ((controls, Fraction(0), Fraction(1, 2)),), 0.0)
    second = geometry._MergeEntry(second_arc, ((controls, Fraction(1, 2), Fraction(1)),), 0.0)
    monkeypatch.setattr(geometry, "certify_cubic_interval_circle", lambda *_args: 0.0)

    assert geometry._merge_arc_entries(first, second, 0.01) is None


def test_represented_biarc_root_perturbation_exceeds_exact_backward_certificate() -> None:
    start_tangent = (1.0, 0.0)
    end_tangent = (0.0, 1.0)
    chord = (math.sqrt(0.5), -math.sqrt(0.5))
    represented_root = math.sqrt(0.5) + math.sqrt(sys.float_info.epsilon)
    exact_root = Fraction.from_float(represented_root)
    exact_chord = tuple(Fraction.from_float(value) for value in chord)
    exact_chord_squared = exact_chord[0] ** 2 + exact_chord[1] ** 2
    exact_residual = abs(2 * exact_root**2 - exact_chord_squared)
    absolute_sum = 2 * exact_root**2 + exact_chord_squared
    unit_roundoff = Fraction.from_float(BINARY64_UNIT_ROUNDOFF)
    accumulated = BIARC_POLYNOMIAL_OPERATION_COUNT * unit_roundoff
    exact_backward_bound = accumulated / (1 - accumulated) * absolute_sum

    assert exact_residual > exact_backward_bound
    assert certify_biarc_root(start_tangent, end_tangent, chord, represented_root) is None


def test_biarc_root_certificate_uses_exact_stored_vectors() -> None:
    start_tangent = (-0.4849139190664588, -0.8745618852291746)
    end_tangent = (-0.4618103023613498, -0.8869787171251172)
    chord = (0.8852837220510297, -0.4650513213307486)

    assert certify_biarc_root(start_tangent, end_tangent, chord, 74.49504729012146) is None
    assert certify_biarc_root(start_tangent, end_tangent, chord, 74.49504729012952) is not None


def test_biarc_root_certificate_refuses_subnormal_conditioning() -> None:
    assert certify_biarc_root((sys.float_info.min / 2.0, 0.0), (0.0, 1.0), (1.0, 0.0), 1.0) is None


def test_mapped_g1_bound_never_uses_infinity_as_acceptance_slack() -> None:
    local_biarc = (
        ReferenceArc.build(_world(0.0, 0.0), _world(1.0, 1.0), _world(0.0, 1.0), Radian(math.pi / 2.0)),
        ReferenceArc.build(_world(1.0, 1.0), _world(2.0, 2.0), _world(2.0, 1.0), Radian(-math.pi / 2.0)),
    )
    origin = (2.0**50, 2.0**50)
    mapped = geometry._map_local_biarc(
        *local_biarc,
        origin,
        (origin[0] + 0.5, origin[1] + 0.5),
        0.2,
    )

    bounds = geometry._mapped_biarc_tangent_bounds(local_biarc, mapped, 0.2)

    assert bounds is None or all(math.isfinite(bound) for bound in bounds)


def test_equal_distance_biarc_refuses_unclosable_mapping(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(geometry, "_mapped_biarc_tangent_bounds", lambda *_args: None)
    controls = tuple(
        (float(point.x), float(point.y))
        for point in (
            _monstera_transform().point(_monstera_source_125().start),
            _monstera_transform().point(_monstera_source_125().control1),
            _monstera_transform().point(_monstera_source_125().control2),
            _monstera_transform().point(_monstera_source_125().end),
        )
    )

    assert geometry._equal_distance_biarc(controls) is None


@pytest.mark.parametrize(
    ("scale_power", "world_origin"),
    (
        (-4, (0.0, 0.0)),
        (0, (64.0, -32.0)),
        (5, (-128.0, 256.0)),
    ),
)
def test_monstera_root_biarc_is_invariant_under_dyadic_scale_and_translation(
    scale_power: int,
    world_origin: tuple[float, float],
) -> None:
    scale = (1.0 / 2.826174326591) * 2.0**scale_power
    transform = SourceToWorld.build(
        source_origin=PdfPoint2.build(0.0, 0.0),
        world_origin=_world(*world_origin),
        scale=MillimetresPerPdfPoint(scale),
    )

    primitives = reconstruct_cubic(
        _monstera_source_125(),
        transform,
        Millimetre(0.10280275256426046 * 2.0**scale_power),
    )

    assert len(primitives) == 2


def test_source_transform_reflects_y_explicitly() -> None:
    transform = SourceToWorld.build(
        source_origin=PdfPoint2.build(10.0, 20.0),
        world_origin=_world(2.0, 3.0),
        scale=MillimetresPerPdfPoint(0.5),
        reflect_source_y=True,
    )

    assert transform.point(PdfPoint2.build(14.0, 26.0)) == _world(4.0, 0.0)


@pytest.mark.parametrize("invalid", (1, "false"))
def test_source_transform_requires_runtime_bool(invalid: object) -> None:
    with pytest.raises(InvalidPublishedPrimitiveError):
        SourceToWorld.build(
            source_origin=PdfPoint2.build(0.0, 0.0),
            world_origin=_world(0.0, 0.0),
            scale=MillimetresPerPdfPoint(1.0),
            reflect_source_y=invalid,  # type: ignore[arg-type]
        )


def test_exact_control_hull_preserves_cancellation_sensitive_turn() -> None:
    magnitude = 2.0**52
    points = ((0.0, 0.0), (magnitude, magnitude - 1.0), (magnitude + 1.0, magnitude), (0.0, 1.0))

    hull = certification.exact_convex_hull(points)

    assert (magnitude, magnitude - 1.0) in hull


@pytest.mark.parametrize(
    ("controls", "centre", "expected"),
    (
        (((1.0, 0.0), (0.0, 2.0), (-1.0, 0.0), (0.0, -2.0)), (0.0, 0.0), (Fraction(0), Fraction(4))),
        (((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)), (3.0, 0.0), (Fraction(4), Fraction(10))),
    ),
)
def test_exact_control_hull_radius_squared_oracle(
    controls: tuple[tuple[float, float], ...],
    centre: tuple[float, float],
    expected: tuple[Fraction, Fraction],
) -> None:
    assert certification.exact_control_hull_radius_squared_bounds(controls, centre) == expected


def test_circle_certificate_refuses_limit_one_ulp_below_exact_radial_bound() -> None:
    controls = ((1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0))
    bound = certification.circle_deviation_upper_bound(controls, (0.0, 0.0), controls[0])

    assert bound is not None
    assert (
        geometry._cubic_within_circle(
            controls,
            ((0.0, 0.0), math.pi / 2.0),
            math.nextafter(bound, -math.inf),
            depth=geometry.MAX_RECONSTRUCTION_SUBDIVISIONS,
        )
        is None
    )


def test_exact_cubic_midpoint_survives_large_translation() -> None:
    magnitude = 2.0**52
    controls = ((magnitude, 0.0), (magnitude + 1.0, 1.0), (magnitude + 2.0, 1.0), (magnitude + 3.0, 0.0))

    assert certification.exact_cubic_point(controls, Fraction(1, 2)) == (
        Fraction.from_float(magnitude) + Fraction(3, 2),
        Fraction(3, 4),
    )


def test_exact_biarc_breakpoint_branch_uses_squared_lengths() -> None:
    assert (
        certification.biarc_branch_at_parameter(
            Fraction(1, 4),
            Fraction(4),
            Fraction(1),
            Fraction(1),
            Fraction(2, 3),
        )
        == 2
    )


def test_biarc_certificate_refuses_limit_one_ulp_below_node_bound() -> None:
    source = _monstera_source_125()
    transform = _monstera_transform()
    controls = tuple(geometry._point_xy(transform.point(point)) for point in (source.start, source.control1, source.control2, source.end))
    biarc = geometry._equal_distance_biarc(controls)
    assert biarc is not None
    proof_biarc = tuple(
        (
            geometry._point_xy(arc.start),
            geometry._point_xy(arc.end),
            geometry._point_xy(arc.centre),
            float(arc.sweep),
        )
        for arc in biarc
    )
    bound = certification.biarc_correspondence_node_bound(
        controls,
        (proof_biarc[0], proof_biarc[1]),
        Fraction(0),
        Fraction(1),
    )

    assert bound is not None
    assert (
        geometry._cubic_within_biarc(
            controls,
            biarc,
            math.nextafter(bound, -math.inf),
            start_parameter=Fraction(0),
            end_parameter=Fraction(1),
            depth=geometry.MAX_RECONSTRUCTION_SUBDIVISIONS,
        )
        is None
    )


def test_finite_biased_trig_enlarges_audited_error(monkeypatch: pytest.MonkeyPatch) -> None:
    angle = Fraction(1, 3)
    baseline = certification.audited_sin_cos(angle)
    original_sin = math.sin
    original_cos = math.cos
    monkeypatch.setattr(certification.math, "sin", lambda value: original_sin(value) + 0.000001)
    monkeypatch.setattr(certification.math, "cos", lambda value: original_cos(value) - 0.000001)

    biased = certification.audited_sin_cos(angle)

    assert biased[1] > baseline[1]


def test_auxiliary_arc_paths_share_stored_join_exactly() -> None:
    join = (1.0, 1.0)
    first = ((0.0, 0.0), join, (0.0, 1.0), math.pi / 2.0)
    second = (join, (2.0, 2.0), (2.0, 1.0), -math.pi / 2.0)

    assert certification.closed_arc_point_box(first, Fraction(1)) == (
        (Fraction(1), Fraction(1)),
        (Fraction(1), Fraction(1)),
    )
    assert certification.closed_arc_point_box(second, Fraction(0)) == (
        (Fraction(1), Fraction(1)),
        (Fraction(1), Fraction(1)),
    )


def test_exact_polar_coefficients_detect_negative_cancellation_turn() -> None:
    magnitude = 2.0**52
    controls = (
        (magnitude, magnitude),
        (2.0 * magnitude, 2.0 * magnitude - 1.0),
        (2.0 * magnitude + 2.0, 2.0 * magnitude),
        (2.0 * magnitude + 4.0, 2.0 * magnitude + 1.0),
    )

    coefficients = certification.exact_polar_bernstein_coefficients(controls, (0.0, 0.0))

    assert coefficients[1] < 0


def test_certified_reconstruction_carries_proved_upper_bound() -> None:
    reconstruction = reconstruct_cubic_certified(
        _quarter_circle_source(),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert isinstance(reconstruction, ReferenceReconstruction)
    assert 0.0 <= float(reconstruction.deviation_upper_bound) < 0.001
    assert reconstruction.primitives == reconstruct_cubic(
        _quarter_circle_source(),
        _identity_transform(),
        Millimetre(0.001),
    )


def test_reconstruction_factory_rejects_empty_chain() -> None:
    with pytest.raises(InvalidReferenceReconstructionError):
        ReferenceReconstruction.build((), Millimetre(0.0))


def test_collinear_cubic_reconstructs_as_certified_line() -> None:
    source = SourceCubic.build(
        PdfPoint2.build(0.0, 0.0),
        PdfPoint2.build(1.0, 0.0),
        PdfPoint2.build(2.0, 0.0),
        PdfPoint2.build(3.0, 0.0),
    )

    reconstruction = reconstruct_cubic_certified(source, _identity_transform(), Millimetre(0.001))

    assert reconstruction.primitives == (ReferenceLine.build(_world(0.0, 0.0), _world(3.0, 0.0)),)
    assert float(reconstruction.deviation_upper_bound) == 0.0


def test_near_collinear_line_bound_covers_independent_curve_oracle() -> None:
    controls = ((0.0, 0.0), (1.0, 0.0001), (2.0, -0.0001), (3.0, 0.0))

    candidate = _certified_line(controls, 0.001)

    assert candidate is not None
    _, upper_bound = candidate
    measured = max(
        abs(3.0 * (1.0 - parameter) ** 2 * parameter * controls[1][1] + 3.0 * (1.0 - parameter) * parameter**2 * controls[2][1])
        for parameter in (index / 1000.0 for index in range(1001))
    )
    assert measured <= upper_bound <= 0.001


def test_line_certificate_uses_finite_segment_endpoint_regions() -> None:
    controls = ((0.0, 0.0), (-1.0, 0.0001), (2.0, 0.0), (3.0, 0.0))

    assert _certified_line(controls, 0.001) is None


def test_line_certificate_refuses_reversed_end_handle() -> None:
    controls = ((0.0, 0.0), (1.0, 0.0), (4.0, 0.0), (3.0, 0.0))

    assert _certified_line(controls, 0.001) is None


def test_source_path_merges_recertified_co_circular_arcs() -> None:
    reconstruction = reconstruct_source_path(
        (_circle_quarter(0, 1), _circle_quarter(1, 2)),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(reconstruction.primitives) == 1
    assert isinstance(reconstruction.primitives[0], ReferenceArc)
    assert float(reconstruction.primitives[0].sweep) == pytest.approx(math.pi)


def test_source_path_refuses_near_co_circular_merge() -> None:
    reconstruction = reconstruct_source_path(
        (_circle_quarter(0, 1), _near_circle_second_quarter()),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(reconstruction.primitives) >= 2


def test_source_path_refuses_opposite_sweep_merge() -> None:
    reconstruction = reconstruct_source_path(
        (_circle_quarter(0, 1), _circle_quarter(1, 0)),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(reconstruction.primitives) == 2


def test_source_path_does_not_merge_beyond_one_turn() -> None:
    reconstruction = reconstruct_source_path(
        tuple(_circle_quarter(index % 4, (index + 1) % 4) for index in range(5)),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert len(reconstruction.primitives) == 2
    assert float(reconstruction.primitives[0].sweep) == pytest.approx(math.tau)


def test_source_path_bound_is_maximum_retained_span_bound() -> None:
    first = reconstruct_cubic_certified(_circle_quarter(0, 1), _identity_transform(), Millimetre(0.001))
    second = reconstruct_cubic_certified(_circle_quarter(1, 2), _identity_transform(), Millimetre(0.001))

    reconstruction = reconstruct_source_path(
        (_circle_quarter(0, 1), _circle_quarter(1, 2)),
        _identity_transform(),
        Millimetre(0.001),
    )

    assert float(reconstruction.deviation_upper_bound) <= max(
        float(first.deviation_upper_bound),
        float(second.deviation_upper_bound),
    )


def test_circle_candidate_rejects_angular_reversal() -> None:
    arcs = reconstruct_cubic(
        _angular_reversal_source(),
        _identity_transform(),
        Millimetre(10.0),
    )

    assert len(arcs) >= 2


def test_projection_closes_and_meets_chord_bound() -> None:
    projection = project_boundary(_rounded_rectangle_boundary(), Millimetre(0.002))

    assert isinstance(projection, PolygonProjection)
    assert projection.points[0] != projection.points[-1]
    assert projection.observed_deviation <= projection.deviation_limit
    assert len(projection.points) >= 8


def test_projection_bound_covers_independently_measured_emitted_chord_midpoints() -> None:
    projection = project_boundary(_unit_circle_boundary(), Millimetre(0.002))
    segments_per_quarter = len(projection.points) // 4
    angle_step = math.pi / (2.0 * segments_per_quarter)
    oracle = 0.0
    for index, start in enumerate(projection.points):
        end = projection.points[(index + 1) % len(projection.points)]
        midpoint_angle = (index + 0.5) * angle_step
        oracle = max(
            oracle,
            _point_segment_distance(
                (math.cos(midpoint_angle), math.sin(midpoint_angle)),
                (float(start.x), float(start.y)),
                (float(end.x), float(end.y)),
            ),
        )

    assert float(projection.observed_deviation) >= oracle


def test_projection_rejects_translated_circle_when_emitted_coordinates_exceed_bound() -> None:
    offset = 2.0**40

    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_unit_circle_boundary(offset), Millimetre(math.ulp(offset)))


def test_projection_rejects_collapsed_emitted_chord() -> None:
    with pytest.raises(InvalidReferenceProjectionError, match="collapsed"):
        project_boundary(_unit_circle_boundary(2.0**50), Millimetre(0.001))


@pytest.mark.parametrize(
    ("order", "sweep"),
    (
        ((0, 1, 2, 3), math.pi / 2.0),
        ((0, 3, 2, 1), -math.pi / 2.0),
    ),
)
def test_projection_rejects_bidirectional_stored_arc_segment_counterexample(
    order: tuple[int, ...],
    sweep: float,
) -> None:
    boundary = _stored_semantics_counterexample(order, sweep)

    with pytest.raises(InvalidReferenceProjectionError, match="exceed"):
        project_boundary(boundary, Millimetre(0.0006852378679477024))


def test_projection_rejects_non_positive_bound() -> None:
    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_rounded_rectangle_boundary(), Millimetre(0.0))
