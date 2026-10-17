"""Shared cut-plane replay classification for production and evidence consumers."""

import math
from collections.abc import Mapping
from enum import Enum
from types import MappingProxyType
from typing import Iterable
from typing import Literal
from typing import Never
from typing import NoReturn
from typing import Optional
from typing import TypeAlias

from compas.geometry import Line
from compas.tolerance import TOL

from compas_cgal.helix import Helix
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import Millimetre

ReplayCategory: TypeAlias = Literal["motion", "rapid", "plunge", "retract", "helix"]


class ReplayRule(Enum):
    """How one operation type replays under the single-plane contract.

    Attributes:
        CUTTING: Engaged; a flat line anchors the cut plane, and an arc or
            circle must lie on it.
        LINK: Engaged at the cut plane; carries no anchor, and an arc or
            circle above the plane replays as a rapid.
        AIR: Never engaged; flat lines and curves replay as rapids.
        RETRACT: Always a retract.
        HELIX: A helical entry (Task 1.2): replays as its exact bore disk,
            category ``"helix"``. Not cap-sampled (like a plunge; its next
            motion is), and its bottom Z is independent cut-plane evidence.
            It carries `compas_cgal.helix.Helix` geometry only; a
            line or planar curve under it, or a helix under any other rule,
            raises `HelixGeometryError`. Never AIR: a rapid replay would skip
            the bore.
    """

    CUTTING = "cutting"
    LINK = "link"
    AIR = "air"
    RETRACT = "retract"
    HELIX = "helix"


REPLAY_RULE_BY_OPERATION: Mapping[OperationType, ReplayRule] = MappingProxyType(
    {
        OperationType.CUT: ReplayRule.CUTTING,
        OperationType.LEAD_IN: ReplayRule.CUTTING,
        OperationType.LEAD_OUT: ReplayRule.CUTTING,
        OperationType.LINK: ReplayRule.LINK,
        OperationType.RETRACT: ReplayRule.RETRACT,
        OperationType.PLUNGE: ReplayRule.AIR,
        OperationType.HELIX: ReplayRule.HELIX,
    }
)
"""Exhaustive replay rule per operation; keys must equal ``set(OperationType)``."""


class CutPlaneRampError(ValueError):
    """A differing-Z line also travels in XY outside the cut-plane model."""


class OffPlaneReplayCurveError(ValueError):
    """An engaged arc or circle lies outside the inferred cutting plane."""


class HelixGeometryError(ValueError):
    """A HELIX operation without `Helix` geometry, or `Helix` geometry under another operation."""


class OffPlaneHelixError(ValueError):
    """A helix's bottom Z lies outside the inferred cutting plane."""


class UnclassifiedOperationError(RuntimeError):
    """An operation type has no entry in ``REPLAY_RULE_BY_OPERATION``."""


class UnhandledReplayRuleError(RuntimeError):
    """A replay rule reached a dispatch that has no branch for it."""


def _unhandled_rule(rule: Never) -> NoReturn:
    """Fail on a rule no ``match`` branch handles.

    Typed ``Never`` like ``typing.assert_never``: every dispatch is an
    exhaustive ``is`` chain ending here, so ``mypy --strict`` rejects a chain
    missing a ``ReplayRule`` member; at runtime it raises a named error
    instead of ``AssertionError``.

    Raises:
        UnhandledReplayRuleError: Always.
    """
    raise UnhandledReplayRuleError(f"no dispatch branch for replay rule {rule!r}")


def rule_is_engaged(rule: ReplayRule) -> bool:
    """Return whether operations under ``rule`` engage material at the cut plane.

    Raises:
        UnhandledReplayRuleError: ``rule`` has no branch here.
    """
    if rule is ReplayRule.CUTTING or rule is ReplayRule.LINK:
        return True
    if rule is ReplayRule.AIR or rule is ReplayRule.RETRACT or rule is ReplayRule.HELIX:
        return False
    _unhandled_rule(rule)


def rule_anchors_cut_plane(rule: ReplayRule) -> bool:
    """Return whether a flat line or a curve under ``rule`` is cut-plane evidence.

    Raises:
        HelixGeometryError: ``rule`` is HELIX, which carries no line or planar
            curve (its anchor is the helix bottom, read by
            `infer_operation_cut_height`).
        UnhandledReplayRuleError: ``rule`` has no branch here.
    """
    if rule is ReplayRule.CUTTING:
        return True
    if rule is ReplayRule.LINK or rule is ReplayRule.AIR or rule is ReplayRule.RETRACT:
        return False
    if rule is ReplayRule.HELIX:
        raise HelixGeometryError("a HELIX operation carries Helix geometry, never a line or planar curve")
    _unhandled_rule(rule)


AUDIT_ENGAGED = frozenset(operation for operation, rule in REPLAY_RULE_BY_OPERATION.items() if rule_is_engaged(rule))


def replay_rule(operation: OperationType) -> ReplayRule:
    """Return the explicit replay rule for one operation type.

    Args:
        operation: The operation type to look up.

    Returns:
        The rule registered for ``operation``.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
    """
    rule = REPLAY_RULE_BY_OPERATION.get(operation)
    if rule is None:
        raise UnclassifiedOperationError(f"no replay rule for operation {operation!r}")
    return rule


def minimum_cut_height(heights: Iterable[Millimetre]) -> Millimetre:
    """Return the lowest supplied motion height, or zero for no motions."""
    return min(heights, default=Millimetre(0.0))


def line_cut_height_anchor(
    operation: OperationType,
    *,
    start_z: Millimetre,
    end_z: Millimetre,
    xy_travel: Millimetre,
) -> Optional[Millimetre]:
    """Return independent cutting-plane evidence carried by one line.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
        HelixGeometryError: ``operation`` is HELIX.
        UnhandledReplayRuleError: The operation's rule has no dispatch branch.
    """
    anchors = rule_anchors_cut_plane(replay_rule(operation))
    if xy_travel <= TOL.absolute and abs(start_z - end_z) > TOL.absolute:
        return min(start_z, end_z)
    if abs(start_z - end_z) <= TOL.absolute and anchors:
        return start_z
    return None


def infer_operation_cut_height(operations: Iterable[ToolpathOperation]) -> Millimetre:
    """Infer a cut plane without letting a curve override independent line evidence.

    A helix bottom is independent evidence, like a plunge foot.

    Raises:
        UnclassifiedOperationError: An operation type has no registered rule.
        HelixGeometryError: A HELIX operation and its geometry disagree.
        UnhandledReplayRuleError: An operation's rule has no dispatch branch.
    """
    line_anchors: list[Millimetre] = []
    curve_heights: list[Millimetre] = []
    for operation in operations:
        geometry = operation.geometry
        if isinstance(geometry, Helix):
            _require_helix_rule(operation.operation)
            line_anchors.append(Millimetre(float(geometry.bottom_z.value)))
        elif isinstance(geometry, Line):
            anchor = line_cut_height_anchor(
                operation.operation,
                start_z=Millimetre(float(geometry.start[2])),
                end_z=Millimetre(float(geometry.end[2])),
                xy_travel=Millimetre(
                    math.hypot(
                        float(geometry.end[0]) - float(geometry.start[0]),
                        float(geometry.end[1]) - float(geometry.start[1]),
                    )
                ),
            )
            if anchor is not None:
                line_anchors.append(anchor)
        elif rule_anchors_cut_plane(replay_rule(operation.operation)):
            curve_heights.append(Millimetre(float(geometry.frame.point[2])))
    return minimum_cut_height(line_anchors if line_anchors else curve_heights)


def classify_line_replay(
    operation: OperationType,
    *,
    start_z: Millimetre,
    end_z: Millimetre,
    xy_travel: Millimetre,
    cut_z: Millimetre,
) -> ReplayCategory:
    """Classify a line under the single-plane replay contract.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
        HelixGeometryError: ``operation`` is HELIX.
        UnhandledReplayRuleError: The operation's rule has no dispatch branch.
        CutPlaneRampError: A differing-Z line also travels in XY.
    """
    rule = replay_rule(operation)
    flat_on_plane: ReplayCategory
    if rule is ReplayRule.RETRACT:
        return "retract"
    elif rule is ReplayRule.CUTTING or rule is ReplayRule.LINK:
        flat_on_plane = "motion"
    elif rule is ReplayRule.AIR:
        flat_on_plane = "rapid"
    elif rule is ReplayRule.HELIX:
        raise HelixGeometryError("a HELIX operation carries Helix geometry, never a line")
    else:
        _unhandled_rule(rule)
    if abs(start_z - end_z) > TOL.absolute:
        if xy_travel > TOL.absolute:
            raise CutPlaneRampError(f"differing-z line has nonzero XY travel (len={xy_travel:.3e})")
        return "plunge" if end_z < start_z else "rapid"
    if start_z > cut_z + TOL.absolute:
        return "rapid"
    return flat_on_plane


def classify_planar_replay(operation: OperationType, *, motion_z: Millimetre, cut_z: Millimetre) -> ReplayCategory:
    """Classify a planar arc or circle under the replay contract.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
        HelixGeometryError: ``operation`` is HELIX.
        UnhandledReplayRuleError: The operation's rule has no dispatch branch.
        OffPlaneReplayCurveError: An engaged curve lies off the cut plane.
    """
    rule = replay_rule(operation)
    if rule is ReplayRule.RETRACT:
        return "retract"
    elif rule is ReplayRule.AIR:
        return "rapid"
    elif rule is ReplayRule.HELIX:
        raise HelixGeometryError("a HELIX operation carries Helix geometry, never a planar curve")
    elif rule is ReplayRule.LINK:
        if motion_z > cut_z + TOL.absolute:
            return "rapid"
    elif rule is not ReplayRule.CUTTING:
        _unhandled_rule(rule)
    if not TOL.is_between(motion_z, cut_z, cut_z, atol=TOL.absolute):
        raise OffPlaneReplayCurveError("engaged curve lies outside the inferred cut plane")
    return "motion"


def _require_helix_rule(operation: OperationType) -> None:
    """Refuse `Helix` geometry under any rule but HELIX.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
        HelixGeometryError: ``operation``'s rule is not HELIX.
    """
    if replay_rule(operation) is not ReplayRule.HELIX:
        raise HelixGeometryError(f"Helix geometry under operation {operation!r}; only HELIX carries it")


def classify_helix_replay(operation: OperationType, *, bottom_z: Millimetre, cut_z: Millimetre) -> ReplayCategory:
    """Classify a helix: it replays as its bore when it ends on the cut plane.

    Args:
        operation: The operation type; must be HELIX.
        bottom_z: The helix's end height.
        cut_z: The inferred cut plane.

    Returns:
        ``"helix"``.

    Raises:
        UnclassifiedOperationError: ``operation`` has no registered rule.
        HelixGeometryError: ``operation``'s rule is not HELIX.
        OffPlaneHelixError: ``bottom_z`` is not on the cut plane.
    """
    _require_helix_rule(operation)
    if not TOL.is_between(bottom_z, cut_z, cut_z, atol=TOL.absolute):
        raise OffPlaneHelixError(f"helix bottom z={bottom_z!r} lies outside the inferred cut plane z={cut_z!r}")
    return "helix"


def classify_operation_replay(operation: ToolpathOperation, cut_z: Millimetre) -> ReplayCategory:
    """Classify one typed toolpath operation without stock side effects.

    Raises:
        UnclassifiedOperationError: The operation type has no registered rule.
        HelixGeometryError: A HELIX operation and its geometry disagree.
        UnhandledReplayRuleError: The operation's rule has no dispatch branch.
        CutPlaneRampError: A differing-Z line also travels in XY.
        OffPlaneReplayCurveError: An engaged curve lies off the cut plane.
        OffPlaneHelixError: A helix ends off the cut plane.
    """
    geometry = operation.geometry
    if isinstance(geometry, Helix):
        return classify_helix_replay(operation.operation, bottom_z=Millimetre(float(geometry.bottom_z.value)), cut_z=cut_z)
    if isinstance(geometry, Line):
        return classify_line_replay(
            operation.operation,
            start_z=Millimetre(float(geometry.start[2])),
            end_z=Millimetre(float(geometry.end[2])),
            xy_travel=Millimetre(
                math.hypot(
                    float(geometry.end[0]) - float(geometry.start[0]),
                    float(geometry.end[1]) - float(geometry.start[1]),
                )
            ),
            cut_z=cut_z,
        )
    return classify_planar_replay(
        operation.operation,
        motion_z=Millimetre(float(geometry.frame.point[2])),
        cut_z=cut_z,
    )
