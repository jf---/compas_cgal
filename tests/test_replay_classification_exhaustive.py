"""Exhaustive, fail-loud replay classification contract (Task 1.2 amendment A0).

Every ``OperationType`` member must carry an explicit replay rule; a member
without one raises instead of silently replaying as a rapid. The pinned tables
below were captured from the HEAD (7b001d47) classifier for every member and
every z-configuration, so the rule table is bit-identical to the old logic.
"""

from enum import Enum
from typing import Callable
from typing import Optional

import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.tolerance import TOL

import compas_cgal.replay_classification as rc
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import Millimetre

CUT_Z = Millimetre(0.0)
AT_TOLERANCE = Millimetre(TOL.absolute)
ABOVE_TOLERANCE = Millimetre(2.0 * TOL.absolute)
BELOW_PLANE = Millimetre(-1.0)
XY_TRAVEL = Millimetre(1.0)

# (start_z, end_z, xy_travel) per named line configuration, cut plane at CUT_Z.
LINE_CONFIGS: dict[str, tuple[Millimetre, Millimetre, Millimetre]] = {
    "flat_at_cut": (CUT_Z, CUT_Z, XY_TRAVEL),
    "flat_at_tol": (AT_TOLERANCE, AT_TOLERANCE, XY_TRAVEL),
    "flat_above": (ABOVE_TOLERANCE, ABOVE_TOLERANCE, XY_TRAVEL),
    "flat_below": (BELOW_PLANE, BELOW_PLANE, XY_TRAVEL),
    "vertical_down": (ABOVE_TOLERANCE, CUT_Z, AT_TOLERANCE),
    "vertical_up": (CUT_Z, ABOVE_TOLERANCE, AT_TOLERANCE),
    "ramp": (ABOVE_TOLERANCE, CUT_Z, ABOVE_TOLERANCE),
    "flat_zero_length_at_cut": (CUT_Z, CUT_Z, Millimetre(0.0)),
    "flat_zero_length_above": (ABOVE_TOLERANCE, ABOVE_TOLERANCE, Millimetre(0.0)),
}

PLANAR_CONFIGS: dict[str, Millimetre] = {
    "at_cut": CUT_Z,
    "at_tol": AT_TOLERANCE,
    "above": ABOVE_TOLERANCE,
    "below": Millimetre(-ABOVE_TOLERANCE),
}

_ENGAGED_LINE = {
    "flat_at_cut": "motion",
    "flat_at_tol": "motion",
    "flat_above": "rapid",
    "flat_below": "motion",
    "vertical_down": "plunge",
    "vertical_up": "rapid",
    "ramp": "CutPlaneRampError",
    "flat_zero_length_at_cut": "motion",
    "flat_zero_length_above": "rapid",
}
# HEAD results of classify_line_replay, keyed by member then configuration.
HEAD_LINE: dict[OperationType, dict[str, str]] = {
    OperationType.CUT: _ENGAGED_LINE,
    OperationType.LEAD_IN: _ENGAGED_LINE,
    OperationType.LEAD_OUT: _ENGAGED_LINE,
    OperationType.LINK: _ENGAGED_LINE,
    OperationType.RETRACT: dict.fromkeys(LINE_CONFIGS, "retract"),
    OperationType.PLUNGE: {
        "flat_at_cut": "rapid",
        "flat_at_tol": "rapid",
        "flat_above": "rapid",
        "flat_below": "rapid",
        "vertical_down": "plunge",
        "vertical_up": "rapid",
        "ramp": "CutPlaneRampError",
        "flat_zero_length_at_cut": "rapid",
        "flat_zero_length_above": "rapid",
    },
    # HELIX rows (Task 1.2): a HELIX operation must carry Helix geometry.
    OperationType.HELIX: dict.fromkeys(LINE_CONFIGS, "HelixGeometryError"),
}

_CUTTING_PLANAR = {
    "at_cut": "motion",
    "at_tol": "motion",
    "above": "OffPlaneReplayCurveError",
    "below": "OffPlaneReplayCurveError",
}
# HEAD results of classify_planar_replay.
HEAD_PLANAR: dict[OperationType, dict[str, str]] = {
    OperationType.CUT: _CUTTING_PLANAR,
    OperationType.LEAD_IN: _CUTTING_PLANAR,
    OperationType.LEAD_OUT: _CUTTING_PLANAR,
    OperationType.LINK: {
        "at_cut": "motion",
        "at_tol": "motion",
        "above": "rapid",
        "below": "OffPlaneReplayCurveError",
    },
    OperationType.RETRACT: dict.fromkeys(PLANAR_CONFIGS, "retract"),
    OperationType.PLUNGE: dict.fromkeys(PLANAR_CONFIGS, "rapid"),
    OperationType.HELIX: dict.fromkeys(PLANAR_CONFIGS, "HelixGeometryError"),
}

# HEAD results of line_cut_height_anchor: "start" = start_z, "min" = min(start_z, end_z), None = no evidence.
_CUTTING_ANCHOR: dict[str, Optional[str]] = {
    "flat_at_cut": "start",
    "flat_at_tol": "start",
    "flat_above": "start",
    "flat_below": "start",
    "vertical_down": "min",
    "vertical_up": "min",
    "ramp": None,
    "flat_zero_length_at_cut": "start",
    "flat_zero_length_above": "start",
}
_NON_CUTTING_ANCHOR: dict[str, Optional[str]] = {
    "flat_at_cut": None,
    "flat_at_tol": None,
    "flat_above": None,
    "flat_below": None,
    "vertical_down": "min",
    "vertical_up": "min",
    "ramp": None,
    "flat_zero_length_at_cut": None,
    "flat_zero_length_above": None,
}
HEAD_ANCHOR: dict[OperationType, dict[str, Optional[str]]] = {
    OperationType.CUT: _CUTTING_ANCHOR,
    OperationType.LEAD_IN: _CUTTING_ANCHOR,
    OperationType.LEAD_OUT: _CUTTING_ANCHOR,
    OperationType.LINK: _NON_CUTTING_ANCHOR,
    OperationType.RETRACT: _NON_CUTTING_ANCHOR,
    OperationType.PLUNGE: _NON_CUTTING_ANCHOR,
    OperationType.HELIX: dict.fromkeys(LINE_CONFIGS, "HelixGeometryError"),
}

CURVE_Z = 3.0
LINE_Z = 5.0
# HEAD results of infer_operation_cut_height for a circle at CURVE_Z with the given member:
# (curve alone, curve + flat LINK line at LINE_Z, curve + flat CUT line at LINE_Z).
_CUTTING_INFER = (CURVE_Z, CURVE_Z, LINE_Z)
_NON_CUTTING_INFER = (0.0, 0.0, LINE_Z)
HEAD_INFER: dict[OperationType, tuple[object, object, object]] = {
    OperationType.CUT: _CUTTING_INFER,
    OperationType.LEAD_IN: _CUTTING_INFER,
    OperationType.LEAD_OUT: _CUTTING_INFER,
    OperationType.LINK: _NON_CUTTING_INFER,
    OperationType.RETRACT: _NON_CUTTING_INFER,
    OperationType.PLUNGE: _NON_CUTTING_INFER,
    OperationType.HELIX: ("HelixGeometryError",) * 3,
}

HEAD_AUDIT_ENGAGED = frozenset({OperationType.CUT, OperationType.LEAD_IN, OperationType.LEAD_OUT, OperationType.LINK})


def _outcome(call: Callable[[], object]) -> object:
    try:
        return call()
    except (rc.CutPlaneRampError, rc.OffPlaneReplayCurveError, rc.HelixGeometryError) as error:
        return type(error).__name__


def test_pinned_tables_cover_every_operation_type() -> None:
    members = set(OperationType)
    assert set(HEAD_LINE) == members
    assert set(HEAD_PLANAR) == members
    assert set(HEAD_ANCHOR) == members
    assert set(HEAD_INFER) == members


def test_every_operation_type_has_an_explicit_replay_rule() -> None:
    table = getattr(rc, "REPLAY_RULE_BY_OPERATION", None)
    assert table is not None, "replay_classification has no exhaustive REPLAY_RULE_BY_OPERATION table"
    assert set(table) == set(OperationType)


def test_audit_engaged_keeps_its_head_value() -> None:
    assert rc.AUDIT_ENGAGED == HEAD_AUDIT_ENGAGED


def test_engagement_reexport_is_the_same_object() -> None:
    from compas_cgal.engagement import AUDIT_ENGAGED

    assert AUDIT_ENGAGED is rc.AUDIT_ENGAGED


@pytest.mark.parametrize("operation", list(OperationType))
@pytest.mark.parametrize("config", list(LINE_CONFIGS))
def test_line_replay_is_bit_identical_to_head(operation: OperationType, config: str) -> None:
    start_z, end_z, xy_travel = LINE_CONFIGS[config]
    outcome = _outcome(lambda: rc.classify_line_replay(operation, start_z=start_z, end_z=end_z, xy_travel=xy_travel, cut_z=CUT_Z))
    assert outcome == HEAD_LINE[operation][config]


@pytest.mark.parametrize("operation", list(OperationType))
@pytest.mark.parametrize("config", list(PLANAR_CONFIGS))
def test_planar_replay_is_bit_identical_to_head(operation: OperationType, config: str) -> None:
    outcome = _outcome(lambda: rc.classify_planar_replay(operation, motion_z=PLANAR_CONFIGS[config], cut_z=CUT_Z))
    assert outcome == HEAD_PLANAR[operation][config]


@pytest.mark.parametrize("operation", list(OperationType))
@pytest.mark.parametrize("config", list(LINE_CONFIGS))
def test_line_anchor_is_bit_identical_to_head(operation: OperationType, config: str) -> None:
    start_z, end_z, xy_travel = LINE_CONFIGS[config]
    expected_kind = HEAD_ANCHOR[operation][config]
    expected = {None: None, "start": start_z, "min": min(start_z, end_z)}.get(expected_kind, expected_kind)
    anchor = _outcome(lambda: rc.line_cut_height_anchor(operation, start_z=start_z, end_z=end_z, xy_travel=xy_travel))
    assert anchor == expected


@pytest.fixture
def table_without_cut(monkeypatch: pytest.MonkeyPatch) -> None:
    table = getattr(rc, "REPLAY_RULE_BY_OPERATION", None)
    assert table is not None, "replay_classification has no exhaustive REPLAY_RULE_BY_OPERATION table"
    monkeypatch.setattr(
        rc,
        "REPLAY_RULE_BY_OPERATION",
        {operation: rule for operation, rule in table.items() if operation is not OperationType.CUT},
    )


@pytest.mark.usefixtures("table_without_cut")
def test_line_replay_of_an_unruled_operation_raises() -> None:
    with pytest.raises(rc.UnclassifiedOperationError):
        rc.classify_line_replay(OperationType.CUT, start_z=CUT_Z, end_z=CUT_Z, xy_travel=XY_TRAVEL, cut_z=CUT_Z)


@pytest.mark.usefixtures("table_without_cut")
def test_planar_replay_of_an_unruled_operation_raises() -> None:
    with pytest.raises(rc.UnclassifiedOperationError):
        rc.classify_planar_replay(OperationType.CUT, motion_z=CUT_Z, cut_z=CUT_Z)


@pytest.mark.usefixtures("table_without_cut")
def test_line_anchor_of_an_unruled_operation_raises() -> None:
    with pytest.raises(rc.UnclassifiedOperationError):
        rc.line_cut_height_anchor(OperationType.CUT, start_z=ABOVE_TOLERANCE, end_z=CUT_Z, xy_travel=AT_TOLERANCE)


def _curve(operation: OperationType) -> ToolpathOperation:
    return ToolpathOperation(
        geometry=Circle(1.0, frame=Frame([0.0, 0.0, CURVE_Z], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])),
        operation=operation,
        path_index=0,
    )


def _flat_line(operation: OperationType) -> ToolpathOperation:
    return ToolpathOperation(geometry=Line([0.0, 0.0, LINE_Z], [1.0, 0.0, LINE_Z]), operation=operation, path_index=0)


@pytest.mark.parametrize("operation", list(OperationType))
def test_inferred_cut_height_from_curves_is_bit_identical_to_head(operation: OperationType) -> None:
    curve = _curve(operation)
    observed = (
        _outcome(lambda: rc.infer_operation_cut_height([curve])),
        _outcome(lambda: rc.infer_operation_cut_height([curve, _flat_line(OperationType.LINK)])),
        _outcome(lambda: rc.infer_operation_cut_height([curve, _flat_line(OperationType.CUT)])),
    )
    assert observed == HEAD_INFER[operation]


@pytest.mark.usefixtures("table_without_cut")
def test_inferred_cut_height_of_an_unruled_curve_raises() -> None:
    with pytest.raises(rc.UnclassifiedOperationError):
        rc.infer_operation_cut_height([_curve(OperationType.CUT)])


class _FifthRule(Enum):
    """Stand-in for a ReplayRule value no dispatch handles yet (e.g. a future HELIX rule)."""

    HELIX = "helix"


@pytest.fixture
def cut_maps_to_a_fifth_rule(monkeypatch: pytest.MonkeyPatch) -> None:
    table = getattr(rc, "REPLAY_RULE_BY_OPERATION", None)
    assert table is not None, "replay_classification has no exhaustive REPLAY_RULE_BY_OPERATION table"
    monkeypatch.setattr(rc, "REPLAY_RULE_BY_OPERATION", {**table, OperationType.CUT: _FifthRule.HELIX})


def _raised_name(call: Callable[[], object]) -> object:
    """Return the raised exception's class name, or the returned value when nothing raised."""
    try:
        return call()
    except Exception as error:  # noqa: BLE001 - the name is asserted by the caller
        return type(error).__name__


UNHANDLED = "UnhandledReplayRuleError"

FIFTH_RULE_PROBES: dict[str, Callable[[], object]] = {
    **{
        f"line_replay[{config}]": (
            lambda start_z=start_z, end_z=end_z, xy_travel=xy_travel: rc.classify_line_replay(OperationType.CUT, start_z=start_z, end_z=end_z, xy_travel=xy_travel, cut_z=CUT_Z)
        )
        for config, (start_z, end_z, xy_travel) in LINE_CONFIGS.items()
    },
    **{
        f"line_anchor[{config}]": (
            lambda start_z=start_z, end_z=end_z, xy_travel=xy_travel: rc.line_cut_height_anchor(OperationType.CUT, start_z=start_z, end_z=end_z, xy_travel=xy_travel)
        )
        for config, (start_z, end_z, xy_travel) in LINE_CONFIGS.items()
    },
    **{
        f"planar_replay[{config}]": (lambda motion_z=motion_z: rc.classify_planar_replay(OperationType.CUT, motion_z=motion_z, cut_z=CUT_Z))
        for config, motion_z in PLANAR_CONFIGS.items()
    },
    "infer_cut_height[curve]": lambda: rc.infer_operation_cut_height([_curve(OperationType.CUT)]),
    "audit_engaged_rule": lambda: rc.rule_is_engaged(_FifthRule.HELIX),  # type: ignore[arg-type]
}


@pytest.mark.usefixtures("cut_maps_to_a_fifth_rule")
@pytest.mark.parametrize("probe", list(FIFTH_RULE_PROBES))
def test_a_rule_without_a_dispatch_branch_raises(probe: str) -> None:
    assert _raised_name(FIFTH_RULE_PROBES[probe]) == UNHANDLED
