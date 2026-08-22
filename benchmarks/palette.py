"""Design tokens for every figure this corpus draws, with the evidence they passed.

Colour here is not a preference and not a per-figure decision: every value below is
a hex from the project's data-viz instance file, and every SET this module ships
was put through that skill's validator, in both modes, on the pairlist the figure
form actually needs. A tool-path plan view is a MAP -- any two marks can end up
side by side -- so the categorical sets are gated on `--pairs all`, which is the
strictly harder pairlist and the one that decided the two designs below.

WHAT THE VALIDATOR DECIDED, and what it rejected (`validate_palette.js`):

* Three categorical slots, `--pairs all`: PASS in both modes. Light worst CVD
  ΔE 9.2 (deutan), worst normal-vision ΔE 24.0; dark 9.4 and 20.9. Light mode
  carries one contrast WARN -- aqua `#1baf7a` at 2.74:1 -- which obliges the
  relief channel the drawing already ships: a legend on every figure, plus a
  stroke style per operation class so colour is never the only carrier.
* FOUR categorical slots, `--pairs all`: FAIL. Slot 4 (yellow `#eda100`) against
  slot 2 (orange `#eb6834`) measures normal-vision ΔE 13.7, under the hard floor
  of 15. So an operation palette of four hues is not available, and the fourth
  and fifth operation classes -- plunge and retract -- take a MARK SHAPE instead
  of a hue. That is the honest encoding anyway: in plan view a plunge and a
  retract are pure Z motions with no extent, so they are events at a point, and
  spending an identity hue on a zero-length mark buys nothing.
* EIGHT categorical slots, `--pairs all`: FAIL hard -- worst CVD ΔE 3.2, worst
  normal-vision ΔE 7.1. There is no ordering of eight hues that passes, which is
  why traversal identity is NOT eight hues. Traversal index is an ORDERED
  quantity (it is the machining order; swapping two chains changes what the
  figure says), so by the skill's own categorical-versus-ordinal test it takes a
  one-hue ordinal ramp, and exact identity comes from direct labels.
* The five-step blue ordinal ramp, `--ordinal`: PASS in both modes -- monotone
  lightness, every adjacent ΔL ≥ 0.06, and the step nearest the surface clears
  the 2:1 floor (light end 2.06:1 on the light surface, dark end 2.15:1 on the
  dark surface). Five is the CAPACITY, not a preference: a ten-step ramp over the
  same window fails the ΔL gate at every adjacent pair (0.047 against 0.06), so
  six or more individually distinguishable steps do not exist on this ramp.

The dark column is not an inversion of the light one. It is the same hue family
stepped for the dark surface, taken from the instance file's dark column, and the
ramp's anchor flips: on light, larger magnitude is DARKER; on dark, larger
magnitude is LIGHTER, so distance from the surface always means more.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Tuple

from compas.tolerance import TOL

from benchmarks.errors import InvalidBandCountError


class Theme(str, Enum):
    """Which surface a figure is drawn for; str-valued so it survives a CLI flag."""

    LIGHT = "light"
    DARK = "dark"


@dataclass(frozen=True)
class Palette:
    """Every token one figure may use, resolved for one theme.

    Attributes:
        theme: The theme these tokens were selected for.
        surface: Chart surface; the figure and every panel are painted with it.
        ink: Primary ink, for text and for the pocket boundary.
        secondary: Secondary ink, for subtitles and for point-event markers.
        muted: Muted ink, for motions the active colour mode does not encode.
        grid: Hairline chrome, used for the swept tool envelope.
        cut: Categorical slot 1 -- cutting motions.
        link: Categorical slot 2 -- linking motions.
        lead: Categorical slot 3 -- lead-in and lead-out motions.
        ramp: The ordinal ramp, ordered LOW magnitude to HIGH magnitude. Its
            length is the number of distinguishable steps the ramp has, and it
            is therefore the fold threshold for any discrete encoding.
    """

    theme: Theme
    surface: str
    ink: str
    secondary: str
    muted: str
    grid: str
    cut: str
    link: str
    lead: str
    ramp: Tuple[str, ...]


LIGHT = Palette(
    theme=Theme.LIGHT,
    surface="#fcfcfb",
    ink="#0b0b0b",
    secondary="#52514e",
    muted="#898781",
    grid="#e1e0d9",
    cut="#2a78d6",
    link="#eb6834",
    lead="#1baf7a",
    ramp=("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"),
)

DARK = Palette(
    theme=Theme.DARK,
    surface="#1a1a19",
    ink="#ffffff",
    secondary="#c3c2b7",
    muted="#898781",
    grid="#2c2c2a",
    cut="#3987e5",
    link="#d95926",
    lead="#199e70",
    ramp=("#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4"),
)


def palette_for(theme: Theme) -> Palette:
    """The tokens for one theme.

    Args:
        theme: The surface the figure will be read on.

    Returns:
        The resolved palette.
    """
    return DARK if theme is Theme.DARK else LIGHT


def band_edges(low: float, high: float, bands: int) -> Tuple[Tuple[float, float], ...]:
    """Split ``[low, high]`` into equal closed-open bands, one per ramp step.

    Args:
        low: Smallest value in the data.
        high: Largest value in the data.
        bands: How many bands to cut, normally the ramp length.

    Returns:
        One ``(lower, upper)`` pair per band, in ascending order. A degenerate
        range (every value equal) yields a single band holding that value.

    Raises:
        InvalidBandCountError: Fewer than one band was asked for.
    """
    if bands < 1:
        raise InvalidBandCountError(f"A ramp needs at least one band; got {bands}.")
    if not TOL.is_positive(high - low):
        return ((low, high),)
    width = (high - low) / bands
    return tuple((low + index * width, low + (index + 1) * width) for index in range(bands))


def band_index(value: float, low: float, high: float, bands: int) -> int:
    """Which band *value* falls in, clamped to the band range.

    Args:
        value: The magnitude to place.
        low: Smallest value in the data.
        high: Largest value in the data.
        bands: How many bands the range was cut into.

    Returns:
        A band index in ``[0, bands - 1]``. A degenerate range places every
        value in band zero, so a constant field reads as one flat colour rather
        than an arbitrary point on the ramp.

    Raises:
        InvalidBandCountError: Fewer than one band was asked for.
    """
    if bands < 1:
        raise InvalidBandCountError(f"A ramp needs at least one band; got {bands}.")
    if not TOL.is_positive(high - low):
        return 0
    scaled = int((value - low) / (high - low) * bands)
    return max(0, min(bands - 1, scaled))
