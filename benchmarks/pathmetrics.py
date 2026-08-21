"""Measure a generated toolpath, separating the entry cut from steady-state cutting.

Every generator that enters solid stock cuts a FULL SLOT on its first circle after
each plunge: the tool is surrounded by material and no stepover, advance, or
engagement control can change that. The entry strategy that does -- a pre-drilled
hole, a helical ramp -- is a different feature from the stepover control, and it
is the stepover control that Held's Figure 6 varies.

So the raw maximum engagement over a whole path is pinned near a full turn for
every generator at every setting, and carries no information about the axis being
swept. A comparison filtered on it is empty by construction: no spacing and no cap
ever complies, every row reads "no compliant path", and the figure cannot show the
effect it exists to show. Measured on the 20x12 pocket at tool 2.0, the raw
maximum is 360.00 degrees at every one of the twelve default trial spacings, while
the after-entry maximum over the same twelve ranges from 98.73 to 205.73 degrees.

This module therefore reports BOTH, and every downstream comparison names which
one it used. The after-entry maximum is an exclusion, not a certificate: the entry
cuts are still there, still over any cap worth setting, and `entry_cuts` says how
many of them a path has.

ONE REFERENCE CAP. Engagement is measured at `REFERENCE_CAP_DEG` for every path,
never at the cap the path is later compared against. The audit refines line
motions against the cap it is handed, so a per-cap measurement would make the
engagement axis of a plot move with the cap axis it is plotted against.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import FrozenSet
from typing import Sequence
from typing import Set

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from benchmarks.depletion import ReplayKind
from benchmarks.depletion import _replay_kind
from benchmarks.errors import UnmeasurableOperationLengthError
from benchmarks.exceedance import exceedance_positions
from benchmarks.spec import MAX_CAP_DEG
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import OperationEngagement
from compas_cgal.engagement import _infer_cut_height
from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# The single cap every engagement measurement in this module is taken at. The
# kernel's contract tops out at a half turn (`benchmarks.spec.MAX_CAP_DEG`), which
# is the loosest legal request and therefore the one that constrains the audit
# least -- a measurement, not a test the path could fail.
REFERENCE_CAP_DEG = MAX_CAP_DEG

# Retract height as a multiple of the tool diameter, matching
# `benchmarks.runner.CLEARANCE_Z_TOOL_DIAMETERS`: links at this height are rapid
# travel the audit records as unmeasured, and expressing it in tool diameters
# keeps every comparison scale-free.
CLEARANCE_Z_TOOL_DIAMETERS = 2.0


@dataclass(frozen=True)
class PathMetrics:
    """What one generated toolpath costs and how hard it cuts.

    Attributes:
        length: Total path length, summed analytically from the operation
            primitives rather than from the tessellated polyline.
        cut_motions: Cut-plane motions the audit measured.
        entry_cuts: Cut motions that are the first material contact after a
            plunge. Each is a full slot by construction.
        max_tea_deg: Worst engagement angle over ALL cut motions, in degrees.
            Pinned near a full turn whenever `entry_cuts` is nonzero, so it is
            reported for honesty and is never a comparison axis.
        max_tea_after_entry_deg: Worst engagement angle over the cut motions that
            are NOT entries, in degrees. The quantity the stepover control moves.
    """

    length: float
    cut_motions: int
    entry_cuts: int
    max_tea_deg: float
    max_tea_after_entry_deg: float


def entry_cut_indices(result: ToolpathResult) -> FrozenSet[int]:
    """Operation indices of the first material contact after each plunge.

    The classification is `benchmarks.depletion`'s, so a path's entry cuts cannot
    drift away from the motions the depletion replay treats as cutting. Rapid
    moves between a plunge and the first cut remove nothing and therefore do not
    disarm the plunge: the cut that follows them still meets virgin stock.

    Args:
        result: The generated toolpath to classify.

    Returns:
        The indices, possibly empty when the path never plunges.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
    """
    cut_z = _infer_cut_height(result.operations)
    entries: Set[int] = set()
    armed = False
    for index, operation in enumerate(result.operations):
        kind = _replay_kind(index, operation, cut_z)
        if kind is ReplayKind.PLUNGE:
            armed = True
        elif kind is ReplayKind.CUT and armed:
            entries.add(index)
            armed = False
    return frozenset(entries)


def path_length(result: ToolpathResult) -> float:
    """Total path length, summed exactly over the operation primitives.

    Analytic, never the tessellated polyline: a chorded circle is shorter than the
    circle, by an amount that depends on the sampling density each generator
    happens to use, which would put a bias into a length comparison between
    generators that has nothing to do with either path.

    Args:
        result: The generated toolpath.

    Returns:
        The summed length of every operation, cutting and rapid alike.

    Raises:
        UnmeasurableOperationLengthError: An operation carries a primitive whose
            length is undefined.
    """
    return sum(_primitive_length(index, operation) for index, operation in enumerate(result.operations))


def _primitive_length(index: int, operation: ToolpathOperation) -> float:
    """Exact length of one operation's primitive.

    compas exposes `length` as a property on `Line` and `Arc` but as an
    unimplemented `Curve` method on `Circle`, whose length is `circumference`;
    reading `.length` uniformly raises `NotImplementedError` on every circle.

    Args:
        index: Position of the operation, for the error message.
        operation: The operation to measure.

    Returns:
        The primitive's length.

    Raises:
        UnmeasurableOperationLengthError: The primitive is not a line, arc, or
            circle.
    """
    geometry = operation.geometry
    if isinstance(geometry, Circle):
        return float(geometry.circumference)
    if isinstance(geometry, (Arc, Line)):
        return float(geometry.length)
    raise UnmeasurableOperationLengthError(f"Operation {index} ({operation.operation.value}) carries geometry {type(geometry).__name__!r}, whose length is not defined.")


def max_tea_after_entry(operations: Sequence[OperationEngagement], entries: FrozenSet[int]) -> float:
    """Worst engaged-run angle over the audited operations that are not entries.

    Args:
        operations: The audit's per-operation records.
        entries: Indices from `entry_cut_indices`.

    Returns:
        The maximum in radians; ``0.0`` when every measured operation is an entry.
    """
    return max((op.max_tea for op in operations if op.op_index not in entries), default=0.0)


def measure_path(spec: PocketSpec, result: ToolpathResult) -> PathMetrics:
    """Audit *result* at the reference cap and report its length and engagement.

    The cap on *spec* is deliberately ignored: see the module docstring on why the
    engagement axis must not move with the cap axis.

    Args:
        spec: The instance whose boundary, holes, and tool the path was generated
            for.
        result: The generated toolpath.

    Returns:
        The path's metrics.

    Raises:
        UnmeasurableOperationLengthError: An operation carries a primitive whose
            length is undefined.
        UnreplayableOperationError: An operation lies outside the cut-plane model.
    """
    report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, math.radians(REFERENCE_CAP_DEG), holes=list(spec.holes))
    entries = entry_cut_indices(result)
    return PathMetrics(
        length=path_length(result),
        cut_motions=report.engaged_ops,
        entry_cuts=len(entries),
        max_tea_deg=math.degrees(report.max_tea),
        max_tea_after_entry_deg=math.degrees(max_tea_after_entry(report.operations, entries)),
    )


def demonstrated_exceedances_after_entry(spec: PocketSpec, result: ToolpathResult) -> int:
    """Non-entry cut motions where the exact engagement predicate actually fired.

    Sampled, so it is a LOWER BOUND on true exceedance and never a certificate --
    the peak of a motion can sit between two sampled positions. It answers "did
    the predicate ever fire away from the entry?", which is a different question
    from "could the cap be proved?" and leans the opposite way. Measured against
    the cap on *spec*, because unlike the engagement measurement this IS a
    question about a particular cap.

    Args:
        spec: The instance whose tool and cap the path is judged against.
        result: The generated toolpath.

    Returns:
        The count of non-entry cut motions with at least one position over the cap.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnsampleableMotionError: A cut motion carries no cutter-centre path.
    """
    entries = entry_cut_indices(result)
    return len({row[0] for row in exceedance_positions(spec, result)} - entries)
