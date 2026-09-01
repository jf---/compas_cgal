from __future__ import annotations

import json
from pathlib import Path

import pytest

import benchmarks.held_reference_cases as cases_module
from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.errors import UnknownHeldReferenceCaseError
from benchmarks.errors import UnsupportedHeldReferenceVersionError
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import project_boundary
from compas_cgal.adaptive.units import Millimetre


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
