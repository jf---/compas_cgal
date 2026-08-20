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
"""

from __future__ import annotations

import math

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
    cut_z = _infer_cut_height(result.operations)
    for index, operation in enumerate(result.operations):
        _deplete_operation(stock, index, operation, spec.tool_radius, cut_z)
    return stock


def _deplete_operation(stock: Stock, index: int, op: ToolpathOperation, tool_radius: float, cut_z: float) -> None:
    """Remove one operation's material from *stock*, or nothing if it cuts none.

    Args:
        stock: The depleting stock.
        index: Position of the operation in its toolpath, for error messages.
        op: The operation to replay.
        tool_radius: Tool radius.
        cut_z: The single cutting-plane height inferred for the toolpath.

    Raises:
        UnreplayableOperationError: The operation is a ramped 3D move, which the
            cut-plane model cannot represent.
    """
    if op.operation == OperationType.RETRACT:
        # Rapid clearance-plane up-move: no material interaction, whatever its geometry.
        return

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
            if z_end < z_start:
                # Plunge: a full-immersion bore removes the tool disk at its foot.
                stock.subtract_disk(float(geometry.end[0]), float(geometry.end[1]), tool_radius)
            return
        if z_start > cut_z + TOL.absolute:
            # Horizontal move above the cutting plane: rapid travel across an
            # already-cleared corridor, not a cut.
            return

    if op.operation not in AUDIT_ENGAGED:
        return
    _subtract_operation(stock, op, tool_radius)
