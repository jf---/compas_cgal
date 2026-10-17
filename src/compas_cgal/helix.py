"""The exact helix of the helical entry: its footprint loop, descent, turn count and closing turn.

A `Helix` descends the whole `CutPlane` interval on a footprint circle in an
integer number of turns, so it ends where it starts: at the footprint's exact
entry ``P``. It then closes with `HELIX_CLOSING_TURNS` flat full turn on the
footprint at the cut plane, in the same operation (ruling "Q6 shape"), again
ending at ``P``: the bore is the full-turn sweep at depth, which is what the
replay removes at the end of the operation. The footprint is authored by the pipeline
(`compas_cgal.pipeline.pencil.helix_footprint`); this core module only models
the descent, so `compas_cgal.toolpath` can carry a `Helix` as operation
geometry without importing `pipeline/`.

Deciding vs reporting: the turn count, the ramp verdict and the footprint's
incidence are exact (`fractions.Fraction`); `length` is REPORTING and feeds no
decision.

The helix is emitted as a ``ToolpathOperation`` with ``operation=HELIX``,
``geometry=<Helix>`` and ``entry_xy = P``. Everything an emitter or survey
needs (``clockwise``, the exact entry and exit, the footprint circle) is on the
helix itself, so emission can append the operation as
is, without the CUT/LINK relabel `PathState.emit_motion` applies to motions
(amendment M5), and a survey can measure the seam exit on the footprint circle
(amendment M6).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from compas.geometry import Circle
from compas.geometry import Frame

from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import LoopEntryIncidenceError
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutPlane
from compas_cgal.units import CutZ
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Flat full turns on the footprint at the cut plane after the ramp, in the same
# HELIX operation (USER 2026-10-16 Q6, ruling "Q6 shape"): one, so the bottom of
# the bore is swept at depth before the pencil leaves P.
HELIX_CLOSING_TURNS = 1


class HelixTurnsError(ValueError):
    """A helix turn count is not a positive integer."""


class HelixDescentError(ValueError):
    """A helix top Z is not strictly above its bottom Z."""


def ramp_turns(descent: ExactMillimetre, footprint_radius: GuideRadius, pitch_per_radius: HelixPitchPerRadius) -> int:
    """The helix turn-count rule: the least ``n >= 1`` with ``descent / n <= exact(K) * rho_h``, exact.

    ``ceil(descent / (exact(K) * rho_h))`` over the rationals, at least one.

    Args:
        descent: Exact Z travelled by the ramp, mm (positive).
        footprint_radius: The footprint radius ``rho_h``.
        pitch_per_radius: The exact ramp-cap surrogate ``K``.

    Returns:
        The ramp's turn count.

    Raises:
        HelixDescentError: *descent* is not positive.
    """
    if descent <= 0:
        raise HelixDescentError(f"a helix ramp descends a positive distance; got {descent!r} mm.")
    return max(1, math.ceil(descent / (pitch_per_radius.exact * Fraction(footprint_radius.value))))


@dataclass(frozen=True)
class Helix:
    """A helical descent on an exactly incident footprint loop, closed by one flat turn at the bottom.

    The operation is ``turns`` ramp turns from ``top_z`` to ``bottom_z``, then
    `HELIX_CLOSING_TURNS` flat turn at ``bottom_z``; both end at the exact
    entry ``P``. The closing turn is always present: the ramp verdict and turn
    count concern the ramp only.

    Attributes:
        footprint: The footprint loop; its exact ``entry`` is where the helix
            starts (at ``top_z``) and ends (at ``bottom_z``).
        top_z: Start height, the cut plane's clearance.
        bottom_z: End height, the cut plane.
        turns: Number of full ramp turns, a positive integer (the closing
            turn is not counted here).
        clockwise: Turn sense in world XY, the root loop's.
    """

    footprint: IncidentLoop
    top_z: ClearanceZ
    bottom_z: CutZ
    turns: int
    clockwise: bool

    def __post_init__(self) -> None:
        if isinstance(self.turns, bool) or not isinstance(self.turns, int) or self.turns < 1:
            raise HelixTurnsError(f"a helix needs a positive integer turn count; got {self.turns!r}.")
        if Fraction(self.top_z.value) <= Fraction(self.bottom_z.value):
            raise HelixDescentError(f"helix top {self.top_z.value!r} is not above its bottom {self.bottom_z.value!r}.")
        if not self.footprint.certifies_incidence():
            raise LoopEntryIncidenceError(f"helix footprint entry {self.footprint.entry!r} is not exactly on its circle.")

    @classmethod
    def build(cls, footprint: IncidentLoop, cut_plane: CutPlane, pitch_per_radius: HelixPitchPerRadius, clockwise: bool) -> Helix:
        """The helix over the whole cut-plane interval with the fewest turns the ramp cap admits.

        ``turns`` is `ramp_turns`: the least positive integer with
        ``descent / turns <= exact(K) * rho_h``, computed over the rationals.

        Args:
            footprint: The footprint loop (see `compas_cgal.pipeline.pencil.helix_footprint`).
            cut_plane: Clearance (top) and cut (bottom) heights.
            pitch_per_radius: The exact ramp-cap surrogate ``K``.
            clockwise: Turn sense, the root loop's.

        Returns:
            The helix.

        Raises:
            LoopEntryIncidenceError: *footprint* fails its incidence certificate.
        """
        descent = ExactMillimetre(Fraction(cut_plane.clearance_z.value) - Fraction(cut_plane.cut_z.value))
        return cls(
            footprint=footprint,
            top_z=cut_plane.clearance_z,
            bottom_z=cut_plane.cut_z,
            turns=ramp_turns(descent, footprint.radius, pitch_per_radius),
            clockwise=clockwise,
        )

    @property
    def descent(self) -> ExactMillimetre:
        """Exact Z travelled."""
        return ExactMillimetre(Fraction(self.top_z.value) - Fraction(self.bottom_z.value))

    @property
    def pitch(self) -> ExactMillimetre:
        """Exact descent per turn."""
        return ExactMillimetre(self.descent / self.turns)

    def ramp_within(self, limit: HelixPitchPerRadius) -> bool:
        """The exact ramp verdict ``pitch <= exact(K) * rho_h``.

        Args:
            limit: The ramp-cap surrogate ``K``.

        Returns:
            True iff the helix ramps within the surrogate cap.
        """
        return self.pitch <= limit.exact * Fraction(self.footprint.radius.value)

    @property
    def entry_xy(self) -> Point2[WorldXY]:
        """Exact start point in XY, the footprint entry ``P``."""
        return self.footprint.entry

    @property
    def exit_xy(self) -> Point2[WorldXY]:
        """Exact end point in XY: ``P`` again, because the ramp and the closing turn are whole turns."""
        return self.footprint.entry

    def ramp_length(self) -> Millimetre:
        """REPORTING: the ramp's path length, ``turns * hypot(2*pi*rho_h, pitch)``, mm. Decides nothing."""
        circumference = 2.0 * math.pi * self.footprint.radius.value
        return Millimetre(self.turns * math.hypot(circumference, float(self.pitch)))

    def closing_turn_length(self) -> Millimetre:
        """REPORTING: the flat closing turn's path length, ``HELIX_CLOSING_TURNS * 2*pi*rho_h``, mm. Decides nothing."""
        return Millimetre(HELIX_CLOSING_TURNS * 2.0 * math.pi * self.footprint.radius.value)

    def length(self) -> Millimetre:
        """REPORTING: the operation's path length, ramp plus closing turn, mm. Decides nothing."""
        return Millimetre(self.ramp_length() + self.closing_turn_length())

    def footprint_circle(self) -> Circle:
        """The footprint as a compas `Circle` on the cut plane (the closing turn).

        Returns:
            The circle; its ``point_at(0)`` is a double rounding of `exit_xy`,
            the exact exit is `exit_xy`.
        """
        circle = self.footprint.circle()
        frame = Frame([self.footprint.centre.x, self.footprint.centre.y, self.bottom_z.value], circle.frame.xaxis, circle.frame.yaxis)
        return Circle(circle.radius, frame=frame)
