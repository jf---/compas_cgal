"""Count the cut motions DEMONSTRATED to exceed the engagement cap.

`EngagementReport.cap_violations` counts operations the certifier could not
CERTIFY, which is a strictly larger and different set. For a circular motion the
audit gives up *before measuring anything* whenever the analytic growth guard
cannot close at that station density, sets `cap_certified = False`, and then
deliberately never consults the per-station `exceeded` flag. Operations with
0 of 36 sampled positions over the cap are counted there all the same. Reading
that number as "violations" is how a generator that had improved 2.4x came to
look like a regression.

The two questions are complementary and neither answers the other:

* **uncertified** -- "could the cap be *proved*?" Sound, conservative, and the
  number that decides whether the certifier is usable on geometry it did not
  author. It over-counts by construction.
* **truly_exceeding** -- "did the exact predicate ever *fire*?" Measured here,
  by sampling. It under-counts by construction: a motion can rise above the cap
  strictly between two sampled positions and go uncounted, so this is a LOWER
  BOUND on true exceedance and NEVER a certificate.

A report that carries only one of them is misleading in whichever direction that
one leans.
"""

from __future__ import annotations

import math
from typing import List
from typing import Sequence
from typing import Tuple

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line

from benchmarks.depletion import CutMotion
from benchmarks.depletion import replay_cuts
from benchmarks.errors import UnsampleableMotionError
from benchmarks.spec import PocketSpec
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.stock import Stock
from compas_cgal.toolpath import ToolpathResult

# Cutter positions sampled per cut motion: 12 around a machining circle is one
# probe per 30 degrees of tool-centre travel, which resolves the engagement swing
# of a loop without paying for the audit's 20-station walk -- this measurement
# only asks WHETHER the exact predicate ever fires on a motion, never where. The
# same count along a bridge keeps segments and loops comparable, and it is the
# density docs/local_depletion.md quotes its "truly exceeding" figures at, so the
# numbers stay comparable across documents.
#
# CONVERGENCE (measured 2026-08-21 on the 12x8 rectangle and the L-shaped pocket,
# both generators, tool 2.0, cap 120 deg): the count is unchanged at every density
# from 12 through 60 and the last count to move does so at 10. Twelve sits just
# above that knee, so the number reported is the converged one and not an artifact
# of how finely the loop was walked.
EXCEEDANCE_SAMPLES_PER_MOTION = 12

# The engagement query is asked WITHOUT gap-closure pessimism: this is a question
# about the material actually engaged at a position, not about what a certifier
# must assume could happen before the next station. Pessimism belongs to the
# certificate, and the certificate is the other column.
NO_GAP_CLOSURE = 0.0


def count_truly_exceeding(spec: PocketSpec, result: ToolpathResult) -> int:
    """Count *result*'s cut motions shown to exceed *spec*'s engagement cap.

    Replays the toolpath against a depleting stock and asks the exact predicate
    `_stock_2.engagement_at` at `EXCEEDANCE_SAMPLES_PER_MOTION` cutter positions
    along each cut motion, counting the motion once if ANY position reports
    `cap_exceeded`. Each motion is measured BEFORE it removes its own material,
    so every answer is about the material that motion truly meets.

    The per-position verdict is exact; the count over a motion is not. Sampling
    can only ever demonstrate exceedance, never absence of it, so the result is a
    LOWER BOUND on the number of genuinely violating motions and must not be
    presented as a certificate. The sound direction is the audit's
    `cap_violations` (operations that could not be certified), which over-counts.

    Args:
        spec: The instance whose tool and cap the toolpath was generated for.
        result: The generated toolpath to replay.

    Returns:
        The number of cut motions with at least one sampled position over the cap.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnsampleableMotionError: A cut motion's geometry is neither a line, an
            arc, nor a circle, so no cutter-centre positions can be read from it.
    """
    ratio = _cap_chord_ratio(spec.tea_cap_rad)
    stock = Stock(spec.polygon, list(spec.holes))
    return sum(1 for motion in replay_cuts(spec, result, stock) if _motion_exceeds(motion, spec.tool_radius, ratio))


def _motion_exceeds(motion: CutMotion, tool_radius: float, cap_chord_ratio: float) -> bool:
    """Report whether any sampled position on *motion* is over the cap.

    Args:
        motion: The cut motion and the stock it is about to cut.
        tool_radius: Tool radius.
        cap_chord_ratio: The cap's exact squared-chord surrogate.

    Returns:
        True as soon as one sampled position reports `cap_exceeded`.

    Raises:
        UnsampleableMotionError: The motion's geometry carries no cutter path.
    """
    raw = motion.stock.raw
    for x, y in _sample_positions(motion):
        _total_tea, _max_run_tea, exceeded = _stock_2.engagement_at(raw, x, y, tool_radius, cap_chord_ratio, NO_GAP_CLOSURE)
        if exceeded:
            return True
    return False


def _sample_positions(motion: CutMotion) -> Sequence[Tuple[float, float]]:
    """Cutter-centre positions to probe along one cut motion.

    A closed circle wraps, so its seam is sampled once; an open arc and a segment
    keep both endpoints, where a motion's engagement is typically extreme. This is
    the station convention `_certify_arc_engagement` uses, at a coarser density.

    Args:
        motion: The cut motion to sample.

    Returns:
        The positions, in travel order.

    Raises:
        UnsampleableMotionError: The geometry is neither a line, an arc, nor a
            circle.
    """
    n = EXCEEDANCE_SAMPLES_PER_MOTION
    geometry = motion.operation.geometry
    if isinstance(geometry, Circle):
        return [_xy(geometry.point_at(i / n)) for i in range(n)]
    if isinstance(geometry, Arc):
        return [_xy(geometry.point_at(i / n)) for i in range(n + 1)]
    if isinstance(geometry, Line):
        x0, y0 = float(geometry.start[0]), float(geometry.start[1])
        x1, y1 = float(geometry.end[0]), float(geometry.end[1])
        return [(x0 + (x1 - x0) * i / n, y0 + (y1 - y0) * i / n) for i in range(n + 1)]
    raise UnsampleableMotionError(f"Operation {motion.index} ({motion.operation.operation.value}) has geometry {type(geometry).__name__!r}, which carries no cutter-centre path.")


def _xy(point: Sequence[float]) -> Tuple[float, float]:
    """Drop a compas point's z, which the cut-plane model fixes."""
    return (float(point[0]), float(point[1]))


def exceedance_positions(spec: PocketSpec, result: ToolpathResult) -> List[Tuple[int, float, float, float]]:
    """Every sampled position over the cap, for attributing a count to geometry.

    The count `count_truly_exceeding` returns says how many motions are over; this
    says which, where, and by how much, so a suspicious number can be traced to a
    place on the pocket instead of being taken on faith.

    Args:
        spec: The instance whose tool and cap the toolpath was generated for.
        result: The generated toolpath to replay.

    Returns:
        One ``(op_index, x, y, max_run_tea_deg)`` row per sampled position that
        reports `cap_exceeded`, in toolpath order.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnsampleableMotionError: A cut motion's geometry carries no cutter path.
    """
    ratio = _cap_chord_ratio(spec.tea_cap_rad)
    stock = Stock(spec.polygon, list(spec.holes))
    rows: List[Tuple[int, float, float, float]] = []
    for motion in replay_cuts(spec, result, stock):
        raw = motion.stock.raw
        for x, y in _sample_positions(motion):
            _total_tea, max_run_tea, exceeded = _stock_2.engagement_at(raw, x, y, spec.tool_radius, ratio, NO_GAP_CLOSURE)
            if exceeded:
                rows.append((motion.index, x, y, math.degrees(max_run_tea)))
    return rows
