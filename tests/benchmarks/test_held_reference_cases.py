from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import benchmarks.held_reference_cases as cases_module
from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.errors import UnknownHeldReferenceCaseError
from benchmarks.errors import UnsupportedHeldReferenceVersionError
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import Figure7Observation
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceBoundary
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_geometry import PolygonProjection
from benchmarks.held_reference_geometry import project_boundary
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY


def _payload(name: str = "figure5") -> dict[str, object]:
    path = cases_module.DATA_DIRECTORY / f"{name}.json"
    value = json.loads(path.read_text())
    assert isinstance(value, dict)
    return value


def _install_payload(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    payload: dict[str, object],
    *,
    name: str = "figure5",
) -> None:
    directory = tmp_path / "held_pfeiffer_2025"
    directory.mkdir()
    (directory / f"{name}.json").write_text(json.dumps(payload))
    monkeypatch.setattr(cases_module, "DATA_DIRECTORY", directory)


def _factory_kwargs(case: HeldReferenceCase) -> dict[str, object]:
    return {
        "name": case.name,
        "authors": case.authors,
        "title": case.title,
        "doi": case.doi,
        "publication_page": case.publication_page,
        "pdf_page": case.pdf_page,
        "figure": case.figure,
        "subfigure": case.subfigure,
        "source_crop": case.source_crop,
        "source_tool_centre": case.source_tool_centre,
        "source_tool_radius": case.source_tool_radius,
        "boundary": case.boundary,
        "reconstruction": case.reconstruction,
        "projection": case.projection,
        "tool_radius": case.tool_radius,
        "tea_cap": case.tea_cap,
        "start_marker": case.start_marker,
        "start_marker_radius": case.start_marker_radius,
        "figure7_observation": case.figure7_observation,
    }


def _next_float(value: float, direction: float, steps: int) -> float:
    for _ in range(steps):
        value = math.nextafter(value, direction)
    return value


def _reconstruction_with_shifted_centre(case: HeldReferenceCase, steps: int) -> ReferenceReconstruction:
    primitives = list(case.reconstruction.primitives)
    first = primitives[0]
    assert isinstance(first, ReferenceArc)
    primitives[0] = ReferenceArc.build(
        first.start,
        first.end,
        Point2[WorldXY].build(_next_float(float(first.centre.x), math.inf, steps), float(first.centre.y)),
        first.sweep,
    )
    return ReferenceReconstruction.build(primitives, case.reconstruction.deviation_upper_bound)


def _projection_with_shifted_point(case: HeldReferenceCase, steps: int) -> PolygonProjection:
    points = list(case.projection.points)
    first = points[0]
    points[0] = Point2[WorldXY].build(_next_float(float(first.x), math.inf, steps), float(first.y))
    return PolygonProjection.build(points, case.projection.deviation_limit, case.projection.observed_deviation)


def test_reconstruction_componentwise_comparison_accepts_one_representation_step_only() -> None:
    case = load_held_reference_case("figure5")

    assert cases_module._reconstruction_matches_within_one_binary64_step(
        _reconstruction_with_shifted_centre(case, 1),
        case.reconstruction,
    )
    assert not cases_module._reconstruction_matches_within_one_binary64_step(
        _reconstruction_with_shifted_centre(case, 2),
        case.reconstruction,
    )


def test_projection_componentwise_comparison_accepts_one_representation_step_only() -> None:
    case = load_held_reference_case("figure5")

    assert cases_module._projection_matches_within_one_binary64_step(
        _projection_with_shifted_point(case, 1),
        case.projection,
    )
    assert not cases_module._projection_matches_within_one_binary64_step(
        _projection_with_shifted_point(case, 2),
        case.projection,
    )


@pytest.mark.parametrize("evidence", ["reconstruction", "projection"])
def test_evidence_upper_bounds_accept_one_upward_representation_step(evidence: str) -> None:
    case = load_held_reference_case("figure5")
    if evidence == "reconstruction":
        computed = case.reconstruction
        recorded = ReferenceReconstruction.build(
            computed.primitives,
            Millimetre(math.nextafter(float(computed.deviation_upper_bound), math.inf)),
        )
        assert cases_module._reconstruction_matches_within_one_binary64_step(recorded, computed)
    else:
        computed_projection = case.projection
        recorded_projection = PolygonProjection.build(
            computed_projection.points,
            computed_projection.deviation_limit,
            Millimetre(math.nextafter(float(computed_projection.observed_deviation), math.inf)),
        )
        assert cases_module._projection_matches_within_one_binary64_step(recorded_projection, computed_projection)


@pytest.mark.parametrize(("direction", "steps"), [(math.inf, 2), (-math.inf, 1)])
@pytest.mark.parametrize("evidence", ["reconstruction", "projection"])
def test_evidence_upper_bounds_reject_two_steps_or_downward_records(
    evidence: str,
    direction: float,
    steps: int,
) -> None:
    case = load_held_reference_case("figure5")
    if evidence == "reconstruction":
        computed = case.reconstruction
        recorded = ReferenceReconstruction.build(
            computed.primitives,
            Millimetre(_next_float(float(computed.deviation_upper_bound), direction, steps)),
        )
        assert not cases_module._reconstruction_matches_within_one_binary64_step(recorded, computed)
    else:
        computed_projection = case.projection
        recorded_projection = PolygonProjection.build(
            computed_projection.points,
            computed_projection.deviation_limit,
            Millimetre(_next_float(float(computed_projection.observed_deviation), direction, steps)),
        )
        assert not cases_module._projection_matches_within_one_binary64_step(recorded_projection, computed_projection)


@pytest.mark.parametrize(
    ("section", "field"),
    [
        ("analytic_boundary", "certified_deviation_upper_bound_mm"),
        ("polygon_projection", "observed_deviation_mm"),
    ],
)
@pytest.mark.parametrize(("direction", "steps", "accepted"), [(math.inf, 1, True), (math.inf, 2, False), (-math.inf, 1, False)])
def test_loader_applies_directional_one_step_bound_contract(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    section: str,
    field: str,
    direction: float,
    steps: int,
    accepted: bool,
) -> None:
    payload = _payload()
    evidence = payload[section]
    assert isinstance(evidence, dict)
    evidence[field] = _next_float(float(evidence[field]), direction, steps)
    _install_payload(monkeypatch, tmp_path, payload)

    if accepted:
        load_held_reference_case("figure5")
    else:
        with pytest.raises(MalformedHeldReferenceCaseError):
            load_held_reference_case("figure5")


def test_public_factory_rejects_mixed_boundary_and_reconstruction() -> None:
    figure5 = load_held_reference_case("figure5")
    upper = load_held_reference_case("figure8_upper")
    kwargs = _factory_kwargs(figure5)
    kwargs["reconstruction"] = upper.reconstruction

    with pytest.raises(MalformedHeldReferenceCaseError):
        HeldReferenceCase.build(**kwargs)  # type: ignore[arg-type]


def test_public_factory_rejects_mixed_boundary_projection() -> None:
    figure5 = load_held_reference_case("figure5")
    upper = load_held_reference_case("figure8_upper")
    kwargs = _factory_kwargs(figure5)
    kwargs["projection"] = upper.projection

    with pytest.raises(MalformedHeldReferenceCaseError):
        HeldReferenceCase.build(**kwargs)  # type: ignore[arg-type]


def test_public_factory_rejects_boundary_tool_mismatch() -> None:
    figure5 = load_held_reference_case("figure5")
    kwargs = _factory_kwargs(figure5)
    kwargs["boundary"] = ReferenceBoundary.build(
        figure5.boundary.primitives,
        ToolRadius.build(2.0),
        figure5.boundary.boundary_stroke_width,
    )

    with pytest.raises(MalformedHeldReferenceCaseError):
        HeldReferenceCase.build(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("authors", ("Wrong",)),
        ("publication_page", 741),
        ("source_crop", (cases_module.PdfPoint2.build(0.0, 0.0), cases_module.PdfPoint2.build(1.0, 1.0))),
        ("tool_radius", ToolRadius.build(2.0)),
        ("tea_cap", cases_module.Degree(79.0)),
        ("start_marker", Point2[WorldXY].build(0.0, 0.0)),
        (
            "figure7_observation",
            Figure7Observation.build(
                publication_page=743,
                pdf_page=14,
                panels=("a", "b", "c"),
                role="shape_only_tool_centre_observation",
                boundary_authority=False,
                numeric_fidelity_gate=False,
            ),
        ),
    ],
)
def test_public_factory_rejects_noncanonical_scalar_metadata(field: str, value: object) -> None:
    kwargs = _factory_kwargs(load_held_reference_case("figure5"))
    kwargs[field] = value

    with pytest.raises(MalformedHeldReferenceCaseError):
        HeldReferenceCase.build(**kwargs)  # type: ignore[arg-type]


def test_public_factory_rejects_noncanonical_publisher_observations() -> None:
    kwargs = _factory_kwargs(load_held_reference_case("figure5"))
    kwargs.update(
        source_tool_centre=cases_module.PdfPoint2.build(0.0, 0.0),
        source_tool_radius=cases_module.PdfPointUnit(3.226562815407),
        start_marker_radius=Millimetre(1.0012107515803956),
    )

    with pytest.raises(MalformedHeldReferenceCaseError):
        HeldReferenceCase.build(**kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("source_tool", "centre", [1e300, 81.00485482270125]),
        ("start_marker", "radius_mm", 999.0),
    ],
)
def test_loader_rejects_noncanonical_publisher_observations(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    section: str,
    field: str,
    value: object,
) -> None:
    payload = _payload()
    parent = payload["normalization"] if section == "source_tool" else payload["machining"]
    assert isinstance(parent, dict)
    observation = parent[section]
    assert isinstance(observation, dict)
    observation[field] = value
    _install_payload(monkeypatch, tmp_path, payload)

    with pytest.raises(MalformedHeldReferenceCaseError):
        load_held_reference_case("figure5")


def test_all_cases_load_in_canonical_order() -> None:
    loaded = load_all_held_reference_cases()

    assert tuple(case.name for case in loaded) == CANONICAL_CASE_NAMES
    assert CANONICAL_CASE_NAMES == (
        "figure5",
        "figure8_upper",
        "figure8_crossed_skis",
        "figure8_monstera",
    )


@pytest.mark.parametrize("name", ["missing", "figure7"])
def test_unknown_case_has_named_error(name: str) -> None:
    with pytest.raises(UnknownHeldReferenceCaseError):
        load_held_reference_case(name)


def test_figure7_observation_is_shape_only_and_unique() -> None:
    loaded = load_all_held_reference_cases()

    assert loaded[0].figure7_observation is not None
    observation = loaded[0].figure7_observation
    assert observation.role == "shape_only_tool_centre_observation"
    assert observation.boundary_authority is False
    assert observation.numeric_fidelity_gate is False
    assert all(case.figure7_observation is None for case in loaded[1:])


def test_cases_close_typed_geometry_and_pocket_contracts() -> None:
    expected_census = {
        "figure5": (5, 26, 65),
        "figure8_upper": (4, 66, 144),
        "figure8_crossed_skis": (2, 56, 106),
        "figure8_monstera": (103, 214, 387),
    }
    for case in load_all_held_reference_cases():
        spec = case.pocket_spec()
        assert spec.name == case.name
        assert spec.family == "held_pfeiffer_2025"
        assert spec.tool_diameter == 2.0
        assert spec.tea_cap_deg == 80.0
        assert spec.holes == ()
        assert case.analytic_signed_area > 0.0
        assert case.polygon_signed_area > 0.0
        assert case.boundary.primitives[0].start == case.boundary.primitives[-1].end
        assert len(case.projection.points) == case.projection_vertex_count
        assert (
            sum(isinstance(primitive, ReferenceLine) for primitive in case.boundary.primitives),
            sum(isinstance(primitive, ReferenceArc) for primitive in case.boundary.primitives),
            len(case.projection.points),
        ) == expected_census[case.name]


def test_projection_is_exact_deterministic_recomputation() -> None:
    for case in load_all_held_reference_cases():
        projected = project_boundary(case.boundary, case.projection.deviation_limit)
        assert projected == case.projection
        tighter = project_boundary(case.boundary, Millimetre(float(case.projection.deviation_limit) / 2.0))
        assert len(tighter.points) >= len(case.projection.points)


def test_round_trip_preserves_semantic_payload() -> None:
    payload = _payload()

    decoded = cases_module.decode_case_document(json.dumps(payload).encode(), expected_name="figure5")
    assert decoded == load_held_reference_case("figure5")


def test_unsupported_version_has_named_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    payload = _payload()
    payload["version"] = 2
    _install_payload(monkeypatch, tmp_path, payload)

    with pytest.raises(UnsupportedHeldReferenceVersionError):
        load_held_reference_case("figure5")


@pytest.mark.parametrize(
    "mutation",
    [
        lambda value: value.update(extra=True),
        lambda value: value["analytic_boundary"].update(frame="source"),  # type: ignore[union-attr]
        lambda value: value["source"].update(unit="mm"),  # type: ignore[union-attr]
        lambda value: value["analytic_boundary"]["primitives"][0].update(kind="quadratic"),  # type: ignore[index,union-attr]
        lambda value: value["polygon_projection"].pop("observed_deviation_mm"),  # type: ignore[union-attr]
        lambda value: value["analytic_boundary"].pop("certified_deviation_upper_bound_mm"),  # type: ignore[union-attr]
        lambda value: value["citation"].update(title="Wrong title"),  # type: ignore[union-attr]
        lambda value: value["publication"].update(publication_page=741),  # type: ignore[union-attr]
        lambda value: value["source"]["crop"].update(minimum=[41.0, 70.0]),  # type: ignore[index,union-attr]
        lambda value: value["normalization"].update(source_origin=[0.0, 0.0]),  # type: ignore[union-attr]
        lambda value: value["normalization"].update(world_origin=[1.0, 0.0]),  # type: ignore[union-attr]
        lambda value: value["analytic_boundary"].update(normalized_stroke_width_mm=1.0),  # type: ignore[union-attr]
        lambda value: value["polygon_projection"].update(deviation_limit_mm=0.01),  # type: ignore[union-attr]
    ],
)
def test_closed_schema_rejects_malformed_documents(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mutation: object,
) -> None:
    payload = _payload()
    mutation(payload)  # type: ignore[operator]
    _install_payload(monkeypatch, tmp_path, payload)

    with pytest.raises(MalformedHeldReferenceCaseError):
        load_held_reference_case("figure5")


def test_filename_name_mismatch_fails(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    payload = _payload()
    payload["name"] = "figure8_upper"
    _install_payload(monkeypatch, tmp_path, payload)

    with pytest.raises(MalformedHeldReferenceCaseError):
        load_held_reference_case("figure5")


@pytest.mark.parametrize(
    "raw",
    [
        b'{"schema":"held-pfeiffer-reference-case","schema":"duplicate","version":1}',
        b'{"schema":"held-pfeiffer-reference-case","version":NaN}',
        b"{",
    ],
)
def test_strict_json_rejects_duplicate_nonfinite_and_invalid_syntax(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    raw: bytes,
) -> None:
    directory = tmp_path / "held_pfeiffer_2025"
    directory.mkdir()
    (directory / "figure5.json").write_bytes(raw)
    monkeypatch.setattr(cases_module, "DATA_DIRECTORY", directory)

    with pytest.raises(MalformedHeldReferenceCaseError):
        load_held_reference_case("figure5")


def test_reversed_boundary_claiming_ccw_fails(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    payload = _payload()
    analytic = payload["analytic_boundary"]
    assert isinstance(analytic, dict)
    primitives = analytic["primitives"]
    assert isinstance(primitives, list)
    reversed_primitives: list[dict[str, object]] = []
    for primitive in reversed(primitives):
        assert isinstance(primitive, dict)
        reversed_primitive = dict(primitive)
        reversed_primitive["start"], reversed_primitive["end"] = primitive["end"], primitive["start"]
        if primitive["kind"] == "arc":
            reversed_primitive["sweep_rad"] = -float(primitive["sweep_rad"])
        reversed_primitives.append(reversed_primitive)
    analytic["primitives"] = reversed_primitives
    _install_payload(monkeypatch, tmp_path, payload)

    with pytest.raises(MalformedHeldReferenceCaseError):
        load_held_reference_case("figure5")
