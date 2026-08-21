"""Neck-severity sweep: a 20x10 pocket pinched to a controllable width.

As the pinch narrows toward the tool diameter, void gaps on the cutter rim shrink
continuously toward zero. That is precisely the regime where a conservative
SAMPLED method must assume every sub-resolution gap is closed and its bound
degenerates to the vacuous full circle, while the exact certifier decides the gap.
This sweep measures where each behaviour begins.

The two threshold functions here read DIFFERENT columns and answer different
questions; the module deliberately offers both rather than one function named
"the" threshold, because a single number would be read as whichever the reader
already had in mind:

* `find_uncertified_pinch_threshold` reads `uncertified` -- the widest neck at
  which the certificate still fails to CLOSE. Sound and conservative, and the
  number a "can this certifier be pointed at geometry it did not author?"
  decision rests on. It over-counts genuine violations by construction.
* `find_exceeding_pinch_threshold` reads `truly_exceeding` -- the widest neck at
  which a sampled cutter position was DEMONSTRATED over the cap. A lower bound,
  never a certificate.

See `benchmarks.exceedance` for why neither substitutes for the other.

MEASURED (2026-08-21, 20x10 dumbbell, tool 2.0, cap 120 deg, via
`benchmarks.runner.run_spec`) -- the two columns behave completely differently
along this axis, and the sweep is only honest if that is said up front:

    pinch/diameter        1.05    1.40    2.00    3.00    4.00   | control
    uncertified / ops    0.481   0.487   0.524   0.522   0.509   |   0.471
    truly_exceeding         57      29      15      17      16   |      18
    truly_exceeding/cut  0.168   0.083   0.044   0.053   0.058   |   0.126

The control column is `analytic.rectangle(20, 10)` -- the SAME pocket with the
notches removed, so it has no neck at all. Its uncertified rate is 0.471, inside
the spread of the whole sweep. **The pinch axis does not move `uncertified`.**
Roughly half of all operations fail to certify at every neck severity, including
none, so that rate is a property of this branch's certifier and not of the
geometry, and `find_uncertified_pinch_threshold` consequently SATURATES: it
returns the widest pinch in whatever sweep it is given. Read a saturated answer
as "no neck in this range certifies", never as "necks this wide are the problem".
That saturation is the number the corpus exists to expose, and it is expected to
move when the certifier does -- which is what makes the threshold worth keeping.

`truly_exceeding` does move: 3.6x in count (57 -> 16) and 2.9x per cut motion
(0.168 -> 0.058) from the tightest neck to the most open, bottoming out at
pinch/diameter 2.0. The tightest neck runs 1.3x the neckless control and the open
ones sit below it, so the neck drives genuine exceedance even though it does not
drive certifiability.
"""

from __future__ import annotations

import math
from typing import Callable
from typing import Optional
from typing import Sequence

from compas.geometry import Polygon

from benchmarks.errors import ImpassableNeckError
from benchmarks.errors import MissingSweepParameterError
from benchmarks.errors import UnpinchedChannelError
from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec

# Pinch widths as multiples of the tool diameter. Below 1.0 the tool cannot pass
# at all; the sweep starts just above and widens to a comfortably open channel.
PINCH_SWEEP_DEFAULT: tuple[float, ...] = (1.05, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0)

# The pocket the neck is cut into. Fixed so the sweep varies neck severity only:
# every instance removes very nearly the same material through a different waist.
DUMBBELL_WIDTH = 20.0
DUMBBELL_HEIGHT = 10.0

# Half-span of the notch along x: each notch runs from the wall at x = 9 up to the
# tip at x = 10 and back down at x = 11, so the neck is a wedge one tool diameter
# long at its narrowest rather than a parallel-sided slot. A wedge is the harder
# instance -- the gap on the cutter rim closes CONTINUOUSLY as the cutter advances
# instead of jumping -- which is the behaviour the sweep exists to resolve.
NOTCH_HALF_SPAN = 1.0

# The sweep coordinate every record in this family carries.
PINCH_PARAM = "pinch"


def traversable_pinch(tool_diameter: float) -> float:
    """The narrowest pinch a tool of this diameter can pass through.

    Args:
        tool_diameter: Cutter diameter.

    Returns:
        The limiting pinch width, equal to the tool diameter.
    """
    return tool_diameter


def dumbbell(pinch: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A 20x10 pocket pinched to *pinch* at x = 10 by two facing reflex notches.

    Args:
        pinch: Channel width at the neck. Must be at least the tool diameter and
            strictly less than the pocket height.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        ImpassableNeckError: The neck is narrower than the tool.
        UnpinchedChannelError: The neck is NaN, or as wide as the pocket, so the
            notch tips have reached the walls and there is no neck left.
    """
    if not math.isfinite(pinch) or pinch >= DUMBBELL_HEIGHT:
        raise UnpinchedChannelError(f"dumbbell: pinch must be finite and below the pocket height {DUMBBELL_HEIGHT:g}, got {pinch!r}.")
    if pinch < traversable_pinch(tool_diameter):
        raise ImpassableNeckError(f"dumbbell: pinch {pinch:g} is narrower than the tool diameter {tool_diameter:g}, so no toolpath reaches through it.")
    h = 0.5 * pinch
    mid_x = 0.5 * DUMBBELL_WIDTH
    left_x, right_x = mid_x - NOTCH_HALF_SPAN, mid_x + NOTCH_HALF_SPAN
    mid_y = 0.5 * DUMBBELL_HEIGHT
    points = [
        [0.0, 0.0, 0.0],
        [left_x, 0.0, 0.0],
        [mid_x, mid_y - h, 0.0],
        [right_x, 0.0, 0.0],
        [DUMBBELL_WIDTH, 0.0, 0.0],
        [DUMBBELL_WIDTH, DUMBBELL_HEIGHT, 0.0],
        [right_x, DUMBBELL_HEIGHT, 0.0],
        [mid_x, mid_y + h, 0.0],
        [left_x, DUMBBELL_HEIGHT, 0.0],
        [0.0, DUMBBELL_HEIGHT, 0.0],
    ]
    return PocketSpec.build(
        name=f"dumbbell_p{pinch:g}",
        family="necks",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={PINCH_PARAM: pinch, "pinch_over_diameter": pinch / tool_diameter},
    )


def pinch_sweep(pinches: tuple[float, ...], tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep neck severity from just-traversable to comfortably open.

    Args:
        pinches: Pinch widths as multiples of the tool diameter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per pinch, ascending.

    Raises:
        ImpassableNeckError: A requested multiple is below 1.0.
        UnpinchedChannelError: A requested pinch reaches the pocket height.
    """
    return [dumbbell(pinch=m * tool_diameter, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for m in sorted(pinches)]


def find_uncertified_pinch_threshold(records: Sequence[MeasurementRecord]) -> Optional[float]:
    """The widest pinch at which the cap could still not be PROVED.

    Reads the `uncertified` column: operations whose certificate did not close,
    including those the growth guard abandoned before measuring anything. This is
    the sound, conservative side of the cap question and the one a certifiability
    threshold is about. It is NOT a count of violations; for the measured lower
    bound on genuine exceedance use `find_exceeding_pinch_threshold`.

    Warning:
        As measured on 2026-08-21 this SATURATES -- every instance in the sweep
        carries uncertified operations, and so does a neckless control, so the
        answer is simply the widest pinch supplied. See the module docstring: a
        saturated answer means "nothing in this range certifies", never "necks
        this wide are the cause".

    Args:
        records: Measurements from a pinch sweep, any order. Instances that
            failed to run are ignored: their zeroed columns are absence of
            measurement, not a clean certificate.

    Returns:
        The widest pinch with at least one uncertified operation, or None when
        every measured instance certified completely.

    Raises:
        MissingSweepParameterError: A record does not carry a `pinch` parameter.
    """
    return _widest_pinch_where(records, lambda record: record.uncertified > 0, "uncertified")


def find_exceeding_pinch_threshold(records: Sequence[MeasurementRecord]) -> Optional[float]:
    """The widest pinch at which a sampled cutter position was over the cap.

    Reads the `truly_exceeding` column, which is measured by sampling and is
    therefore a LOWER BOUND on genuine exceedance, never a certificate: a motion
    can rise above the cap strictly between two sampled positions and go
    uncounted. For the sound direction use `find_uncertified_pinch_threshold`.

    Args:
        records: Measurements from a pinch sweep, any order. Instances that
            failed to run are ignored.

    Returns:
        The widest pinch with at least one demonstrated exceedance, or None when
        no measured instance produced one.

    Raises:
        MissingSweepParameterError: A record does not carry a `pinch` parameter.
    """
    return _widest_pinch_where(records, lambda record: record.truly_exceeding > 0, "truly_exceeding")


def _widest_pinch_where(records: Sequence[MeasurementRecord], predicate: Callable[[MeasurementRecord], bool], column: str) -> Optional[float]:
    """The widest pinch among measured records satisfying *predicate*.

    The single implementation both thresholds share, so the ONLY difference
    between them is the column named at the call site.

    Args:
        records: Measurements from a pinch sweep, any order.
        predicate: Callable deciding whether one record counts.
        column: Column name the predicate reads, for the error message.

    Returns:
        The widest qualifying pinch, or None when none qualifies.

    Raises:
        MissingSweepParameterError: A record does not carry a `pinch` parameter.
    """
    hits: list[float] = []
    for record in records:
        if record.error is not None:
            continue
        if PINCH_PARAM not in record.params:
            raise MissingSweepParameterError(
                f"{record.name!r} carries params {sorted(record.params)} and no {PINCH_PARAM!r}, so its {column} count cannot be placed on the sweep."
            )
        if predicate(record):
            hits.append(float(record.params[PINCH_PARAM]))
    return max(hits) if hits else None
