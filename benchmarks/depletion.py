"""Replay a generated toolpath onto a fresh stock, outside any timed region.

The kernel diagnostics in `benchmarks.instrument` are only meaningful on a
DEPLETED stock: a virgin arrangement has a handful of vertices and two-digit
rationals no matter how hard the instance was, so probing before depletion pins
both readings to their virgin values and neither can ever move. The engagement
audit already depletes a stock internally, but discards it, so the corpus
re-creates that final state here -- deliberately after the timed region, because
`probe_digits` collapses the lazy-exact filter and would otherwise inflate the
very certification time it exists to explain.

The classification below is the shared semantic decision used by depletion and
quality evidence: retracts and clearance-height links are rapid travel and remove
nothing, plunges bore a disk, and cut-plane motions remove their swept area. The
public replay projection deliberately folds retracts into rapid travel; evidence
consumers retain the four-way source category.

`replay_events` is the primitive: it walks the toolpath and hands each cut-plane
motion to the caller BEFORE that motion removes its own material, which is the
only moment at which the material a motion actually meets still exists, and
reports each plunge with the exact verdict of whether its disk bored material.
`replay_cuts` is its projection onto the cut motions, and `replay_depletion` is
the degenerate consumer that wants only the end state.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import Iterator
from typing import Union

from compas.geometry import Line

from benchmarks.errors import UnreplayableOperationError
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import _infer_cut_height
from compas_cgal.engagement import _subtract_operation
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import ReplayCategory
from compas_cgal.replay_classification import classify_line_replay
from compas_cgal.replay_classification import classify_planar_replay
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import Millimetre


class ReplayKind(enum.Enum):
    """How the cut-plane model treats one toolpath operation.

    Attributes:
        RAPID: Removes nothing -- a retract, an up-move, or a link above the
            cutting plane.
        PLUNGE: Bores the tool disk at its foot.
        CUT: Removes its swept area at the cutting plane, and is the only kind the
            engagement audit measures.
    """

    RAPID = "rapid"
    PLUNGE = "plunge"
    CUT = "cut"


@dataclass(frozen=True)
class CutMotion:
    """One cut-plane motion, paired with the material it is about to meet.

    Attributes:
        index: Position of the operation in its toolpath.
        operation: The motion itself.
        stock: The depleting stock in the state this motion cuts into. Valid only
            until the replay resumes, which immediately subtracts this motion's
            swept area -- so a consumer must read what it needs before yielding
            control back, and must not mutate it.
    """

    index: int
    operation: ToolpathOperation
    stock: Stock


# A plunge bores the tool disk at its foot. As a sweep that disk is the full-turn
# sweep of a guide circle of radius ZERO: `Stock.circle_sweep_removal` with
# rho = 0 builds the disk of radius r about the foot exactly and decides its
# regularized interior intersection with the stock without mutating it. Not a
# tolerance: the guide radius of a point.
PLUNGE_GUIDE_RADIUS = 0.0


@dataclass(frozen=True)
class PlungeEvent:
    """One plunge, judged against the stock it bores into.

    Attributes:
        index: Position of the operation in its toolpath.
        x: X of the plunge's foot, as authored (binary64, injected exactly).
        y: Y of the plunge's foot, as authored.
        into_material: Whether the plunge disk meets the stock's interior, an
            EXACT verdict taken before the disk is subtracted. Touching the
            cleared boundary is not boring material.
    """

    index: int
    x: float
    y: float
    into_material: bool


ReplayEvent = Union[CutMotion, PlungeEvent]


def replay_events(spec: PocketSpec, result: ToolpathResult, stock: Stock) -> Iterator[ReplayEvent]:
    """Deplete *stock* by replaying *result*, yielding every cut motion and plunge.

    A cut motion is yielded BEFORE it removes its swept area; a plunge is
    yielded after its exact into-material verdict and before its disk is
    removed, but carries no stock, so the caller cannot observe the order.
    Rapid moves are skipped. On exhaustion *stock* holds the end state.

    Args:
        spec: The instance whose tool the toolpath was generated for.
        result: The generated toolpath to replay.
        stock: A fresh stock for *spec*, depleted in place as the replay proceeds.

    Yields:
        One `CutMotion` per cut-plane motion and one `PlungeEvent` per plunge,
        in toolpath order.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
    """
    cut_z = _infer_cut_height(result.operations)
    for index, operation in enumerate(result.operations):
        kind = _replay_kind(index, operation, cut_z)
        if kind is ReplayKind.RAPID:
            continue
        if kind is ReplayKind.PLUNGE:
            geometry = operation.geometry
            assert isinstance(geometry, Line)  # guaranteed by _replay_kind's PLUNGE branch
            x, y = float(geometry.end[0]), float(geometry.end[1])
            verdict = stock.circle_sweep_removal(x, y, PLUNGE_GUIDE_RADIUS, spec.tool_radius)
            yield PlungeEvent(index=index, x=x, y=y, into_material=verdict is SweepRemoval.REMOVES)
            stock.subtract_disk(x, y, spec.tool_radius)
            continue
        yield CutMotion(index=index, operation=operation, stock=stock)
        _subtract_operation(stock, operation, spec.tool_radius)


def replay_cuts(spec: PocketSpec, result: ToolpathResult, stock: Stock) -> Iterator[CutMotion]:
    """Deplete *stock* by replaying *result*, yielding each cut motion before it cuts.

    The projection of `replay_events` onto its cut motions: rapid moves are
    skipped and plunges bore their disk without being yielded. On exhaustion
    *stock* holds the state the audit ended with.

    Args:
        spec: The instance whose tool the toolpath was generated for.
        result: The generated toolpath to replay.
        stock: A fresh stock for *spec*, depleted in place as the replay proceeds.

    Yields:
        One `CutMotion` per material-removing cut-plane motion, in toolpath order.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
    """
    return (event for event in replay_events(spec, result, stock) if isinstance(event, CutMotion))


def replay_depletion(spec: PocketSpec, result: ToolpathResult) -> Stock:
    """Rebuild the stock the audit ended with, by re-cutting *result* from virgin.

    Args:
        spec: The instance whose boundary and tool the toolpath was generated for.
        result: The generated toolpath to replay.

    Returns:
        A stock depleted by every material-removing motion in *result*.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        InvalidPolygonError: The instance's boundary or a hole is degenerate.
    """
    stock = Stock(spec.polygon, list(spec.holes))
    for _motion in replay_cuts(spec, result, stock):
        pass
    return stock


def _replay_kind(index: int, op: ToolpathOperation, cut_z: float) -> ReplayKind:
    """Classify one operation under the cut-plane model, without touching any stock.

    Args:
        index: Position of the operation in its toolpath, for error messages.
        op: The operation to classify.
        cut_z: The single cutting-plane height inferred for the toolpath.

    Returns:
        The operation's `ReplayKind`.

    Raises:
        UnreplayableOperationError: The operation is a ramped 3D move, which the
            cut-plane model cannot represent.
    """
    geometry = op.geometry
    if isinstance(geometry, Line):
        z_start = float(geometry.start[2])
        z_end = float(geometry.end[2])
        xy_travel = math.hypot(float(geometry.end[0]) - float(geometry.start[0]), float(geometry.end[1]) - float(geometry.start[1]))
        category = _classify_replay_category(
            index,
            op.operation,
            cut_z,
            is_line=True,
            z_start=z_start,
            z_end=z_end,
            xy_travel=xy_travel,
        )
    else:
        motion_z = float(geometry.frame.point[2])
        category = _classify_replay_category(
            index,
            op.operation,
            cut_z,
            is_line=False,
            z_start=motion_z,
            z_end=0.0,
            xy_travel=0.0,
        )
    if category == "motion":
        return ReplayKind.CUT
    if category == "plunge":
        return ReplayKind.PLUNGE
    return ReplayKind.RAPID


def _classify_replay_category(
    index: int,
    operation: OperationType,
    cut_z: float,
    *,
    is_line: bool,
    z_start: float,
    z_end: float,
    xy_travel: float,
) -> ReplayCategory:
    """Adapt geometry facts to the shared cut-plane replay decision."""
    try:
        if is_line:
            return classify_line_replay(
                operation,
                start_z=Millimetre(z_start),
                end_z=Millimetre(z_end),
                xy_travel=Millimetre(xy_travel),
                cut_z=Millimetre(cut_z),
            )
        return classify_planar_replay(operation, motion_z=Millimetre(z_start), cut_z=Millimetre(cut_z))
    except (CutPlaneRampError, OffPlaneReplayCurveError) as error:
        raise UnreplayableOperationError(f"Operation {index} ({operation.value}) {error}.") from error
