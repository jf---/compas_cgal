"""Replay a generated toolpath onto a fresh stock, outside any timed region.

The kernel diagnostics in `benchmarks.instrument` are only meaningful on a
DEPLETED stock: a virgin arrangement has a handful of vertices and two-digit
rationals no matter how hard the instance was, so probing before depletion pins
both readings to their virgin values and neither can ever move. The engagement
audit already depletes a stock internally, but discards it, so the corpus
re-creates that final state here -- deliberately after the timed region, because
`probe_digits` collapses the lazy-exact filter and would otherwise inflate the
very certification time it exists to explain.

The classification below is a MIRROR of ``_replay_operation`` in
``src/compas_cgal/engagement.py`` with the measurement removed: retracts and
clearance-height links are rapid travel and remove nothing, plunges bore a disk,
and cut-plane motions remove their swept area. The geometry-to-boolean mapping
and the cut-plane inference are imported from that module rather than re-derived,
so only the branch structure is duplicated and a change to the audit's cut-plane
rule cannot silently desynchronise the two.

`replay_cuts` is the primitive: it walks the toolpath and hands each cut-plane
motion to the caller BEFORE that motion removes its own material, which is the
only moment at which the material a motion actually meets still exists.
`replay_depletion` is the degenerate consumer that wants only the end state.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import Iterator

from compas.geometry import Line
from compas.tolerance import TOL

from benchmarks.errors import UnreplayableOperationError
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import AUDIT_ENGAGED
from compas_cgal.engagement import _infer_cut_height
from compas_cgal.engagement import _subtract_operation
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult


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


def replay_cuts(spec: PocketSpec, result: ToolpathResult, stock: Stock) -> Iterator[CutMotion]:
    """Deplete *stock* by replaying *result*, yielding each cut motion before it cuts.

    Rapid moves are skipped and plunges bore their disk without being yielded;
    only the motions the engagement audit measures reach the caller. On exhaustion
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
    cut_z = _infer_cut_height(result.operations)
    for index, operation in enumerate(result.operations):
        kind = _replay_kind(index, operation, cut_z)
        if kind is ReplayKind.RAPID:
            continue
        if kind is ReplayKind.PLUNGE:
            geometry = operation.geometry
            assert isinstance(geometry, Line)  # guaranteed by _replay_kind's PLUNGE branch
            stock.subtract_disk(float(geometry.end[0]), float(geometry.end[1]), spec.tool_radius)
            continue
        yield CutMotion(index=index, operation=operation, stock=stock)
        _subtract_operation(stock, operation, spec.tool_radius)


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
    if op.operation == OperationType.RETRACT:
        # Rapid clearance-plane up-move: no material interaction, whatever its geometry.
        return ReplayKind.RAPID

    geometry = op.geometry
    if isinstance(geometry, Line):
        z_start = float(geometry.start[2])
        z_end = float(geometry.end[2])
        if abs(z_start - z_end) > TOL.absolute:
            xy_travel = math.hypot(float(geometry.end[0]) - float(geometry.start[0]), float(geometry.end[1]) - float(geometry.start[1]))
            if xy_travel > TOL.absolute:
                raise UnreplayableOperationError(
                    f"Operation {index} ({op.operation.value}) is a differing-z line with nonzero XY travel (len={xy_travel:.3e}); "
                    f"ramped 3D cutting is outside the cut-plane depletion model."
                )
            # Downward: a full-immersion bore. Upward: retract-shaped, removes nothing.
            return ReplayKind.PLUNGE if z_end < z_start else ReplayKind.RAPID
        if z_start > cut_z + TOL.absolute:
            # Horizontal move above the cutting plane: rapid travel across an
            # already-cleared corridor, not a cut.
            return ReplayKind.RAPID

    if op.operation not in AUDIT_ENGAGED:
        return ReplayKind.RAPID
    return ReplayKind.CUT
