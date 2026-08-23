from __future__ import annotations

import math

import pytest

from compas_cgal import _stock_2

WORLD_X = (1.0, 0.0, 0.0)
WORLD_Y = (0.0, 1.0, 0.0)


def test_line_plane_classification_is_exact_after_binary64_ingress() -> None:
    result = _stock_2.classify_audit_line(
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        0.0,
        5.0,
        "cut",
    )

    assert isinstance(result, _stock_2.AuditSegmentMotion2)

    with pytest.raises(_stock_2.AuditOffPlaneError):
        off_plane = math.nextafter(0.0, math.inf)
        _stock_2.classify_audit_line(
            (0.0, 0.0, off_plane),
            (1.0, 0.0, off_plane),
            0.0,
            5.0,
            "cut",
        )


@pytest.mark.parametrize(
    ("start", "end", "role", "expected_type"),
    [
        ((0.0, 0.0, 5.0), (0.0, 0.0, 0.0), "plunge", _stock_2.AuditVerticalPlunge2),
        ((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), "retract", _stock_2.AuditVerticalRetract2),
        ((0.0, 0.0, 5.0), (1.0, 0.0, 5.0), "link", _stock_2.AuditClearanceTransport2),
    ],
)
def test_line_role_is_derived_by_native_geometry(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    role: str,
    expected_type: type[object],
) -> None:
    result = _stock_2.classify_audit_line(start, end, 0.0, 5.0, role)

    assert isinstance(result, expected_type)


@pytest.mark.parametrize(
    ("cut_z", "clearance_z"),
    [(0.0, 0.0), (1.0, 0.0)],
)
def test_native_classifier_rejects_invalid_plane_contract(
    cut_z: float,
    clearance_z: float,
) -> None:
    with pytest.raises(_stock_2.AuditInvalidPlaneError):
        _stock_2.classify_audit_line(
            (0.0, 0.0, cut_z),
            (1.0, 0.0, cut_z),
            cut_z,
            clearance_z,
            "cut",
        )


def test_native_line_classifier_rejects_ramps_and_role_contradictions() -> None:
    with pytest.raises(_stock_2.AuditUnsupportedGeometryError):
        _stock_2.classify_audit_line(
            (0.0, 0.0, 5.0),
            (1.0, 0.0, 0.0),
            0.0,
            5.0,
            "link",
        )

    with pytest.raises(_stock_2.AuditContradictoryRoleError):
        _stock_2.classify_audit_line(
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            0.0,
            5.0,
            "retract",
        )


def test_native_circle_classifier_returns_opaque_exact_motion() -> None:
    result = _stock_2.classify_audit_circle(
        (1.0, 2.0, 0.0),
        WORLD_X,
        WORLD_Y,
        2.0,
        True,
        0.0,
        5.0,
        "cut",
    )

    assert isinstance(result, _stock_2.AuditCircleMotion2)


@pytest.mark.parametrize(
    ("xaxis", "yaxis"),
    [
        ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        ((2.0, 0.0, 0.0), (0.0, 2.0, 0.0)),
        ((1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
        ((1.0, 0.0, 0.0), (0.0, -1.0, 0.0)),
    ],
)
def test_native_circle_classifier_rejects_noncanonical_frame(
    xaxis: tuple[float, float, float],
    yaxis: tuple[float, float, float],
) -> None:
    with pytest.raises(_stock_2.AuditUnsupportedGeometryError):
        _stock_2.classify_audit_circle(
            (1.0, 2.0, 0.0),
            xaxis,
            yaxis,
            2.0,
            False,
            0.0,
            5.0,
            "cut",
        )


@pytest.mark.parametrize("radius", [0.0, -1.0])
def test_native_circle_classifier_rejects_nonpositive_radius(radius: float) -> None:
    with pytest.raises(_stock_2.AuditUnsupportedGeometryError):
        _stock_2.classify_audit_circle(
            (1.0, 2.0, 0.0),
            WORLD_X,
            WORLD_Y,
            radius,
            False,
            0.0,
            5.0,
            "cut",
        )


def test_native_arc_classifier_owns_phase_and_sweep_direction() -> None:
    result = _stock_2.classify_audit_arc(
        (1.0, 2.0, 0.0),
        WORLD_X,
        WORLD_Y,
        2.0,
        0.0,
        math.pi / 2.0,
        False,
        0.0,
        5.0,
        "cut",
    )

    assert isinstance(result, _stock_2.AuditArcMotion2)

    with pytest.raises(_stock_2.AuditContradictoryOrientationError):
        _stock_2.classify_audit_arc(
            (1.0, 2.0, 0.0),
            WORLD_X,
            WORLD_Y,
            2.0,
            0.0,
            math.pi / 2.0,
            True,
            0.0,
            5.0,
            "cut",
        )


@pytest.mark.parametrize("end_angle", [0.0, math.nextafter(math.tau, math.inf)])
def test_native_arc_classifier_rejects_degenerate_sweep(end_angle: float) -> None:
    with pytest.raises(_stock_2.AuditUnsupportedGeometryError):
        _stock_2.classify_audit_arc(
            (1.0, 2.0, 0.0),
            WORLD_X,
            WORLD_Y,
            2.0,
            0.0,
            end_angle,
            False,
            0.0,
            5.0,
            "cut",
        )


def test_native_classifier_rejects_nonfinite_input() -> None:
    with pytest.raises(_stock_2.AuditNonFiniteInputError):
        _stock_2.classify_audit_line(
            (math.nan, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            0.0,
            5.0,
            "cut",
        )


@pytest.mark.parametrize(
    "native_type",
    [
        _stock_2.AuditSegmentMotion2,
        _stock_2.AuditCircleMotion2,
        _stock_2.AuditArcMotion2,
        _stock_2.AuditVerticalPlunge2,
        _stock_2.AuditVerticalRetract2,
        _stock_2.AuditClearanceTransport2,
    ],
)
def test_native_classification_values_cannot_be_forged(
    native_type: type[object],
) -> None:
    with pytest.raises(TypeError):
        native_type()
