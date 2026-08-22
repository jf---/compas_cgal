"""Which mark every motion of a tool path gets, under one colour mode.

Colour carries exactly one quantity per figure, and everything the active mode
does not encode goes muted rather than taking a hue it has no claim to.
`ColourBy.OPERATION` is the only categorical mode -- what kind of motion this is.
The other three encode a magnitude on one blue ramp: where the move sits along
the path, which chain it belongs to, or how hard it is measured to cut.
`benchmarks.palette` carries the validator evidence for both sets, including why
traversal identity is a ramp with direct labels rather than eight hues.

MOTIONS WITHOUT EXTENT BECOME POINTS. A plunge and a retract are pure Z moves;
projected into the cutting plane they are one location, not a stroke. They
therefore take a shaped marker -- plunge down, retract up -- in secondary ink, in
every colour mode. A link whose two ends coincide gets the same treatment,
because that is what it is.

Nothing here imports matplotlib: a mark is a description, and
`benchmarks.plotting` is what draws one.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict
from typing import Literal
from typing import Mapping
from typing import Optional
from typing import Sequence

from compas.tolerance import TOL

from benchmarks.errors import EngagementLengthMismatchError
from benchmarks.errors import MissingEngagementDataError
from benchmarks.errors import UnknownColourModeError
from benchmarks.palette import Palette
from benchmarks.palette import band_edges
from benchmarks.palette import band_index
from benchmarks.pathgeometry import Motion

# Stroke weights in points. Cutting is the figure; linking and leads are
# scaffolding, drawn lighter and with a dash pattern, so operation class survives
# greyscale printing and colour-vision deficiency without leaning on hue.
CUT_WIDTH_PT = 1.6
LEAD_WIDTH_PT = 1.3
LINK_WIDTH_PT = 0.9

# Size of a point event's marker, in points.
POINT_MARKER_PT = 4.5

# The line styles a mark may carry. Spelled out rather than left as `str` so that
# a typo is a type error instead of a matplotlib exception at draw time.
LineStyle = Literal["-", "--", "-.", ":", "none"]

# Legend ordering. Explicit integers rather than insertion order, so two panels of
# a comparison that start with different motions still produce one legend in one
# order.
_ORDER_CUT = 10
_ORDER_LEAD_IN = 20
_ORDER_LEAD_OUT = 21
_ORDER_LINK = 30
_ORDER_RAMP_BASE = 100
_ORDER_FOLDED = 200
_ORDER_UNENCODED = 210
_ORDER_PLUNGE = 300
_ORDER_RETRACT = 301
_ORDER_POINT_OTHER = 302

_POINT_MARKERS = {"plunge": "v", "retract": "^"}
_POINT_ORDERS = {"plunge": _ORDER_PLUNGE, "retract": _ORDER_RETRACT}


class ColourBy(str, Enum):
    """Which quantity the stroke colour encodes; str-valued so it survives a CLI flag.

    Attributes:
        OPERATION: Categorical -- what kind of motion this is.
        TRAVERSAL: Ordinal -- which chain the motion belongs to, in machining
            order. Exact identity comes from `annotate_traversals`, because the
            ramp carries only as many steps as it validated as distinguishable.
        SEQUENCE: Sequential -- how far along the whole path the motion sits.
        ENGAGEMENT: Sequential -- the measured engagement angle of the motion,
            supplied by the caller so this module never needs the exact kernel.
    """

    OPERATION = "operation"
    TRAVERSAL = "traversal"
    SEQUENCE = "sequence"
    ENGAGEMENT = "engagement"


@dataclass(frozen=True)
class Mark:
    """The complete visual specification of one motion.

    Attributes:
        label: Legend text; motions sharing a label share one legend entry.
        order: Sort key for the legend.
        colour: Stroke or marker colour.
        width_pt: Stroke weight in points.
        linestyle: The line style, or ``"none"`` for a point mark.
        marker: A matplotlib marker, or ``""`` for a stroke.
    """

    label: str
    order: int
    colour: str
    width_pt: float
    linestyle: LineStyle
    marker: str


def encode(
    motions: Sequence[Motion],
    colour_by: ColourBy,
    palette: Palette,
    engagement_deg: Optional[Sequence[Optional[float]]],
) -> Dict[int, Mark]:
    """Give every motion its mark, under one colour mode.

    Args:
        motions: The motions to encode.
        colour_by: What the colour encodes.
        palette: The theme's tokens.
        engagement_deg: One measurement per operation, for the engagement mode.

    Returns:
        Operation index to mark.

    Raises:
        MissingEngagementDataError: The engagement mode was asked for with no
            measurements.
        EngagementLengthMismatchError: Not exactly one measurement per operation.
        UnknownColourModeError: The colour mode is not one this module encodes.
    """
    if colour_by is ColourBy.OPERATION:
        return {motion.index: _operation_mark(motion, palette) for motion in motions}
    if colour_by is ColourBy.TRAVERSAL:
        return _traversal_marks(motions, palette)
    if colour_by is ColourBy.SEQUENCE:
        return _sequence_marks(motions, palette)
    if colour_by is ColourBy.ENGAGEMENT:
        return _engagement_marks(motions, palette, engagement_deg)
    raise UnknownColourModeError(f"{colour_by!r} is not a colour mode this module encodes.")


def _operation_mark(motion: Motion, palette: Palette) -> Mark:
    """The categorical mark for one motion.

    Args:
        motion: The motion.
        palette: The theme's tokens.

    Returns:
        Its mark. A plunge or retract that does have extent in the cutting plane
        is a ramp rather than a pure Z move, and is named as one instead of being
        quietly folded into the link colour.
    """
    if motion.is_point:
        return _point_mark(motion, palette)
    if motion.name == "cut":
        return Mark(label="cut", order=_ORDER_CUT, colour=palette.cut, width_pt=CUT_WIDTH_PT, linestyle="-", marker="")
    if motion.name == "lead_in":
        return Mark(label="lead in", order=_ORDER_LEAD_IN, colour=palette.lead, width_pt=LEAD_WIDTH_PT, linestyle=":", marker="")
    if motion.name == "lead_out":
        return Mark(label="lead out", order=_ORDER_LEAD_OUT, colour=palette.lead, width_pt=LEAD_WIDTH_PT, linestyle="-.", marker="")
    if motion.name == "link":
        return Mark(label="link", order=_ORDER_LINK, colour=palette.link, width_pt=LINK_WIDTH_PT, linestyle="--", marker="")
    return Mark(
        label=f"{motion.name} (ramped)",
        order=_POINT_ORDERS.get(motion.name, _ORDER_POINT_OTHER),
        colour=palette.secondary,
        width_pt=LINK_WIDTH_PT,
        linestyle="--",
        marker="",
    )


def _point_mark(motion: Motion, palette: Palette) -> Mark:
    """The mark for a motion with no extent in the cutting plane.

    Args:
        motion: The motion.
        palette: The theme's tokens.

    Returns:
        A shaped marker in secondary ink: plunge down, retract up, and a dot for
        anything else that turns out to have coincident ends.
    """
    spoken = motion.name.replace("_", " ")
    if motion.name in _POINT_MARKERS:
        return Mark(label=spoken, order=_POINT_ORDERS[motion.name], colour=palette.secondary, width_pt=0.0, linestyle="none", marker=_POINT_MARKERS[motion.name])
    return Mark(label=f"{spoken} (no extent)", order=_ORDER_POINT_OTHER, colour=palette.secondary, width_pt=0.0, linestyle="none", marker="o")


def _unencoded_mark(label: str, palette: Palette) -> Mark:
    """The muted mark a ramp mode gives to motions it does not encode.

    Args:
        label: Legend text.
        palette: The theme's tokens.

    Returns:
        The mark.
    """
    return Mark(label=label, order=_ORDER_UNENCODED, colour=palette.muted, width_pt=LINK_WIDTH_PT, linestyle="--", marker="")


def _ramp_mark(label: str, slot: int, palette: Palette) -> Mark:
    """A cutting stroke on the ordinal ramp.

    Args:
        label: Legend text for this step.
        slot: Index into the ramp.
        palette: The theme's tokens.

    Returns:
        The mark.
    """
    return Mark(label=label, order=_ORDER_RAMP_BASE + slot, colour=palette.ramp[slot], width_pt=CUT_WIDTH_PT, linestyle="-", marker="")


def _traversal_marks(motions: Sequence[Motion], palette: Palette) -> Dict[int, Mark]:
    """Colour cutting motions by which chain they belong to.

    The ramp carries `len(palette.ramp)` distinguishable steps and not one more,
    so chains past that fold into a single muted bucket the legend counts. Hues
    are never cycled: two chains sharing a colour would be a claim the reader has
    no way to check, while a fold is visible and named, and `annotate_traversals`
    still gives every chain its exact index.

    Args:
        motions: The motions to encode.
        palette: The theme's tokens.

    Returns:
        Operation index to mark.
    """
    chains = sorted({motion.path_index for motion in motions})
    slots = {chain: slot for slot, chain in enumerate(chains[: len(palette.ramp)])}
    folded = len(chains) - len(slots)
    marks: Dict[int, Mark] = {}
    for motion in motions:
        if motion.is_point:
            marks[motion.index] = _point_mark(motion, palette)
        elif not motion.is_cutting:
            marks[motion.index] = _unencoded_mark("link / rapid", palette)
        elif motion.path_index in slots:
            marks[motion.index] = _ramp_mark(f"chain {motion.path_index}", slots[motion.path_index], palette)
        else:
            marks[motion.index] = Mark(
                label=f"other ({folded} chains)",
                order=_ORDER_FOLDED,
                colour=palette.muted,
                width_pt=CUT_WIDTH_PT,
                linestyle="-",
                marker="",
            )
    return marks


def _sequence_marks(motions: Sequence[Motion], palette: Palette) -> Dict[int, Mark]:
    """Colour cutting motions by how far along the path they sit.

    Position is measured by arc length over the whole path, rapids included,
    because that is what "along the path" means to the machine.

    Args:
        motions: The motions to encode.
        palette: The theme's tokens.

    Returns:
        Operation index to mark.
    """
    total = float(sum(motion.length for motion in motions))
    bands = len(palette.ramp)
    edges = band_edges(0.0, 1.0, bands)
    marks: Dict[int, Mark] = {}
    travelled = 0.0
    for motion in motions:
        midpoint = travelled + 0.5 * motion.length
        travelled += motion.length
        if motion.is_point:
            marks[motion.index] = _point_mark(motion, palette)
            continue
        if not motion.is_cutting:
            marks[motion.index] = _unencoded_mark("link / rapid", palette)
            continue
        fraction = midpoint / total if TOL.is_positive(total) else 0.0
        slot = band_index(fraction, 0.0, 1.0, bands)
        low, high = edges[min(slot, len(edges) - 1)]
        marks[motion.index] = _ramp_mark(f"{low:.0%}–{high:.0%} along path", slot, palette)
    return marks


def _engagement_marks(
    motions: Sequence[Motion],
    palette: Palette,
    engagement_deg: Optional[Sequence[Optional[float]]],
) -> Dict[int, Mark]:
    """Colour cutting motions by their measured engagement angle.

    The measurements are the caller's: this module never computes engagement, so
    a figure can be drawn without the exact kernel, and the figure cannot
    disagree with the audit that produced the numbers.

    Args:
        motions: The motions to encode.
        palette: The theme's tokens.
        engagement_deg: One measurement per operation, None where a motion
            carries none.

    Returns:
        Operation index to mark.

    Raises:
        MissingEngagementDataError: No measurements were supplied.
        EngagementLengthMismatchError: Not exactly one measurement per operation.
    """
    if engagement_deg is None:
        raise MissingEngagementDataError("ColourBy.ENGAGEMENT draws measurements it is given; pass engagement_deg=[...], one entry per operation.")
    values = list(engagement_deg)
    if len(values) != len(motions):
        raise EngagementLengthMismatchError(f"engagement_deg has {len(values)} entries for {len(motions)} operations; it needs exactly one per operation.")

    measured = [value for motion, value in zip(motions, values) if value is not None and motion.is_cutting and not motion.is_point]
    bands = len(palette.ramp)
    low, high = (min(measured), max(measured)) if measured else (0.0, 0.0)
    edges = band_edges(low, high, bands)
    marks: Dict[int, Mark] = {}
    for motion, value in zip(motions, values):
        if motion.is_point:
            marks[motion.index] = _point_mark(motion, palette)
        elif value is None or not motion.is_cutting:
            marks[motion.index] = _unencoded_mark("not measured", palette)
        else:
            slot = band_index(value, low, high, bands)
            band_low, band_high = edges[min(slot, len(edges) - 1)]
            marks[motion.index] = _ramp_mark(f"{band_low:.0f}–{band_high:.0f}°", slot, palette)
    return marks


def legend_entries(marks: Mapping[int, Mark]) -> Dict[str, Mark]:
    """The legend entries a set of marks will produce.

    Computed before the figure exists, because the legend's row count decides how
    much of the page it needs and therefore how tall the figure is.

    Args:
        marks: Operation index to mark.

    Returns:
        Label to the first mark carrying it, in first-drawn order.
    """
    entries: Dict[str, Mark] = {}
    for mark in marks.values():
        entries.setdefault(mark.label, mark)
    return entries
