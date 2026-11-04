"""Machining loops whose entry point lies EXACTLY on the guide circle.

A plunge at a loop's entry is an exact no-op on the loop's swept annulus iff
the entry is at distance exactly ``rho`` from the centre: ``disk(e, r)`` lies in
the closed annulus ``[rho - r, rho + r]`` iff ``|e - c| = rho``. A double entry
``c + rho * (cos t, sin t)`` is ulps off that circle for every non-axis phase, so
its plunge disk pokes a sliver out of the annulus. No replay can fix that: it
is an AUTHORING property, and this module authors it.

`IncidentLoop.build` chooses the geometry a generator emits -- it is not a snap
at the exact seam, which stays snap-free:

* the entry direction is the rational unit vector of the quarter-chart atlas
  (``src/exact_circle_chart_2.cpp``, ``src/exact_circle_chart_atlas_2.h``):
  quarter ``k`` rotates ``((Q^2 - j^2), 2 j Q) / (Q^2 + j^2)``, ``t = j / Q`` the
  chart parameter nearest the target phase;
* centre, entry and radius are multiples of ``2^-K`` mm, the radius a multiple
  of the reduced hypotenuse ``m`` times ``2^-K``, so ``e - c = (rho / m)(a, b)``
  is a grid vector and every coordinate is an exact double;
* the CERTIFICATE ``(ex - cx)^2 + (ey - cy)^2 == rho^2`` is checked over
  `fractions.Fraction` on the emitted doubles; a miss raises
  `LoopEntryIncidenceError`. There is no fallback.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Protocol

from compas.geometry import Circle
from compas.geometry import Frame

from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

# Chart denominator Q. The chart angle 2*atan(t) has slope <= 2 in t, and the
# chosen t = j/Q is within 1/(2Q) of the target's, so |dphase| <= 1/Q rad = 2^-10
# rad (0.056 deg).
INCIDENCE_CHART_DENOMINATOR = 2**10

# Dyadic grid exponent K: centre, entry and radius are integer multiples of
# 2^-K mm. Coordinates of magnitude <= 2^11 mm (centre bound plus one radius)
# are then integers below 2^51 times 2^-K, inside binary64's 53-bit significand.
INCIDENCE_GRID_EXPONENT = 40

# Largest supported |centre coordinate| and radius, mm: 2^10 mm = 1.024 m, the
# magnitude the grid exponent above is derived for.
INCIDENCE_COORDINATE_BOUND_MM = 2.0**10

# Bound on |authored radius - target|, mm. The radius quantum is m * 2^-K with
# m <= 2 Q^2 = 2^21, so rounding moves rho by at most 2^20 * 2^-40 = 2^-20 mm
# (9.5e-7 mm); the promise is the ruling's 1e-5 mm.
INCIDENCE_RADIUS_BOUND_MM = 1e-5

# Bound on |authored phase - target|, rad: 1/Q, derived at the chart denominator.
INCIDENCE_PHASE_BOUND_RAD = 1.0 / INCIDENCE_CHART_DENOMINATOR

# Bound on |authored centre - target| per coordinate, mm: half a grid step,
# 2^-(K+1) = 4.5e-13 mm.
INCIDENCE_CENTRE_BOUND_MM = 2.0 ** -(INCIDENCE_GRID_EXPONENT + 1)

# Number of quarter charts of the atlas: `IncidentLoop.chart` is in [0, 3].
INCIDENCE_CHART_COUNT = 4

_QUARTER_TURN = math.pi / 2.0


class NonFiniteIncidentLoopPhaseError(ValueError):
    """The target phase is not a finite angle (centre and radius are validated by their unit types)."""


class IncidentLoopRangeError(ValueError):
    """A centre coordinate or the radius exceeds `INCIDENCE_COORDINATE_BOUND_MM`."""


class IncidentLoopRadiusBelowQuantumError(ValueError):
    """The target radius rounds to zero on its chart's radius quantum."""


class GuideLoop(Protocol):
    """A full machining loop as the in-process engagement queries read it: guide circle and start point, doubles.

    `IncidentLoop` is one (its ``entry`` exactly on the circle); a published
    Held circle wrapped by ``benchmarks.held_helix_entry`` is another (its
    ``entry`` the certified nearest double of the exact entry).
    """

    @property
    def centre(self) -> Point2[WorldXY]: ...

    @property
    def radius(self) -> GuideRadius: ...

    @property
    def entry(self) -> Point2[WorldXY]: ...


class LoopEntryIncidenceError(RuntimeError):
    """The emitted doubles fail the exact incidence certificate (broken invariant)."""


class ChartParameterRangeError(ValueError):
    """A chart index is outside ``[0, INCIDENCE_CHART_COUNT)`` or a numerator outside ``[0, Q]``."""


@dataclass(frozen=True)
class IncidentLoop:
    """One machining loop authored with an exactly incident entry.

    Attributes:
        centre: Guide centre, an exact grid point.
        radius: Authored guide radius.
        entry: Entry point on the guide, exactly at distance *radius* from
            *centre*.
        phase: Direction of *entry* from *centre* -- the realised phase,
            which REPORTS the exact entry; it decides nothing.
        chart: Quarter-chart index in ``[0, 3]`` of the entry direction.
        chart_numerator: Chart parameter numerator ``j``, ``t = j / Q``.
    """

    centre: Point2[WorldXY]
    radius: GuideRadius
    entry: Point2[WorldXY]
    phase: Radian
    chart: int
    chart_numerator: int

    @classmethod
    def build(cls, centre: Point2[WorldXY], radius: GuideRadius, phase: Radian) -> IncidentLoop:
        """Author the loop nearest the target with an exactly incident entry.

        Args:
            centre: Target guide centre (finite by its type).
            radius: Target guide radius (positive by its type).
            phase: Target entry direction (any finite angle).

        Returns:
            The authored loop, within `INCIDENCE_CENTRE_BOUND_MM`,
            `INCIDENCE_RADIUS_BOUND_MM` and `INCIDENCE_PHASE_BOUND_RAD` of the
            target.

        Raises:
            NonFiniteIncidentLoopPhaseError: *phase* is not finite.
            IncidentLoopRangeError: An input exceeds the supported magnitude.
            IncidentLoopRadiusBelowQuantumError: *radius* rounds to zero.
            LoopEntryIncidenceError: The certificate fails on the emitted doubles.
        """
        cx, cy = float(centre.x), float(centre.y)
        target_radius, target_phase = float(radius.value), float(phase)
        if not math.isfinite(target_phase):
            raise NonFiniteIncidentLoopPhaseError(f"incident loop needs a finite phase; got {target_phase!r}.")
        if max(abs(cx), abs(cy), target_radius) > INCIDENCE_COORDINATE_BOUND_MM:
            raise IncidentLoopRangeError(f"incident loop supports |coordinate|, radius <= {INCIDENCE_COORDINATE_BOUND_MM} mm; got ({cx!r}, {cy!r}), {target_radius!r}.")

        chart, numerator = _chart_parameter(target_phase)
        a, b, hypotenuse = chart_triple(chart, numerator)
        grid = 2**INCIDENCE_GRID_EXPONENT
        steps = round(Fraction(target_radius) * grid / hypotenuse)
        if steps < 1:
            raise IncidentLoopRadiusBelowQuantumError(f"radius {target_radius!r} mm is below the radius quantum {hypotenuse / grid!r} mm of chart ({chart}, {numerator}).")

        exact_cx = Fraction(round(Fraction(cx) * grid), grid)
        exact_cy = Fraction(round(Fraction(cy) * grid), grid)
        exact_radius = Fraction(steps * hypotenuse, grid)
        exact_ex = exact_cx + Fraction(steps * a, grid)
        exact_ey = exact_cy + Fraction(steps * b, grid)
        emitted = tuple(exact_double(value) for value in (exact_cx, exact_cy, exact_radius, exact_ex, exact_ey))
        loop = cls(
            centre=Point2[WorldXY](Millimetre(emitted[0]), Millimetre(emitted[1])),
            radius=GuideRadius(Millimetre(emitted[2])),
            entry=Point2[WorldXY](Millimetre(emitted[3]), Millimetre(emitted[4])),
            phase=chart_phase(chart, numerator),
            chart=chart,
            chart_numerator=numerator,
        )
        if not loop.certifies_incidence():
            raise LoopEntryIncidenceError(f"entry {loop.entry!r} is not exactly at radius {loop.radius!r} from {loop.centre!r}.")
        return loop

    def certifies_incidence(self) -> bool:
        """The exact certificate ``|entry - centre|^2 == radius^2`` over the rationals.

        Returns:
            True iff the emitted doubles satisfy it exactly.
        """
        dx = Fraction(self.entry.x) - Fraction(self.centre.x)
        dy = Fraction(self.entry.y) - Fraction(self.centre.y)
        return dx * dx + dy * dy == Fraction(self.radius.value) ** 2

    def certifies_chart(self) -> bool:
        """The exact check that *chart*, *chart_numerator* and *phase* describe the entry direction.

        ``(entry - centre) * m == (a, b) * radius`` over the rationals for the
        chart triple ``(a, b, m)``, and *phase* is the double `chart_phase`
        derives from ``(chart, chart_numerator)``. The phase comparison is a
        consistency check on a DERIVED value, not a transcendental decision:
        `build` stores exactly `chart_phase`'s result, so the two doubles are
        bit-identical by construction and any difference is a contradicting
        record.

        Returns:
            True iff all three agree exactly.

        Raises:
            ChartParameterRangeError: The chart parameter is out of range.
        """
        a, b, m = chart_triple(self.chart, self.chart_numerator)
        rho = Fraction(self.radius.value)
        dx = Fraction(self.entry.x) - Fraction(self.centre.x)
        dy = Fraction(self.entry.y) - Fraction(self.centre.y)
        return dx * m == a * rho and dy * m == b * rho and self.phase == chart_phase(self.chart, self.chart_numerator)

    def circle(self) -> Circle:
        """The loop as a compas `Circle` carrying the authored centre and radius.

        Returns:
            A world-XY circle whose frame x-axis points along the realised phase.
            Its ``point_at(0)`` is a double rounding of `entry`; the exact entry
            is `entry`.
        """
        ux, uy = math.cos(self.phase), math.sin(self.phase)
        frame = Frame([self.centre.x, self.centre.y, 0.0], [ux, uy, 0.0], [-uy, ux, 0.0])
        return Circle(self.radius.value, frame=frame)


def _chart_parameter(phase: float) -> tuple[int, int]:
    """The quarter chart and numerator ``j`` whose chart angle is nearest *phase*.

    Placement intent only: the angle picks a rational direction, which is then
    exact; nothing downstream decides on the double.

    Args:
        phase: Target angle, rad.

    Returns:
        ``(chart, j)`` with chart in ``[0, 3]`` and ``j`` in ``[0, Q]``.
    """
    turned = math.fmod(phase, math.tau)
    if turned < 0.0:
        turned += math.tau
    chart = min(int(turned // _QUARTER_TURN), INCIDENCE_CHART_COUNT - 1)
    local = turned - chart * _QUARTER_TURN
    numerator = round(math.tan(local / 2.0) * INCIDENCE_CHART_DENOMINATOR)
    return chart, min(max(numerator, 0), INCIDENCE_CHART_DENOMINATOR)


def chart_triple(chart: int, numerator: int) -> tuple[int, int, int]:
    """The reduced Pythagorean triple ``(a, b, m)`` of a chart point, ``a^2 + b^2 = m^2``.

    Mirrors the atlas of ``src/exact_circle_chart_atlas_2.h``: quarter 0 is
    ``((Q^2 - j^2), 2 j Q) / (Q^2 + j^2)``, each further quarter a rotation by a
    quarter turn.

    Args:
        chart: Quarter index in ``[0, 3]``.
        numerator: Chart numerator ``j`` in ``[0, Q]``.

    Returns:
        ``(a, b, m)`` with ``gcd(a, b) == 1``.

    Raises:
        ChartParameterRangeError: *chart* or *numerator* is out of range.
    """
    if not 0 <= chart < INCIDENCE_CHART_COUNT or not 0 <= numerator <= INCIDENCE_CHART_DENOMINATOR:
        raise ChartParameterRangeError(f"chart parameter ({chart!r}, {numerator!r}) is outside [0, {INCIDENCE_CHART_COUNT}) x [0, {INCIDENCE_CHART_DENOMINATOR}].")
    q = INCIDENCE_CHART_DENOMINATOR
    a, b, m = q * q - numerator * numerator, 2 * numerator * q, q * q + numerator * numerator
    divisor = math.gcd(a, b)
    a, b, m = a // divisor, b // divisor, m // divisor
    for _ in range(chart):
        a, b = -b, a
    return a, b, m


def chart_phase(chart: int, numerator: int) -> Radian:
    """The phase of a chart point: ``atan2(b, a)`` of its `chart_triple`.

    The ONE computation of a loop's phase: `IncidentLoop.build` stores this
    result, so a phase re-derived here from ``(chart, numerator)`` is
    bit-identical to the stored one. The phase REPORTS the exact direction;
    it decides nothing.

    Args:
        chart: Quarter index in ``[0, 3]``.
        numerator: Chart numerator ``j`` in ``[0, Q]``.

    Returns:
        The direction angle, rad.

    Raises:
        ChartParameterRangeError: *chart* or *numerator* is out of range.
    """
    a, b, _ = chart_triple(chart, numerator)
    return Radian(math.atan2(b, a))


def exact_double(value: Fraction) -> float:
    """The double equal to *value*, which the grid construction guarantees exists.

    Args:
        value: A grid rational.

    Returns:
        ``float(value)``.

    Raises:
        LoopEntryIncidenceError: *value* is not exactly a double.
    """
    emitted = float(value)
    if Fraction(emitted) != value:
        raise LoopEntryIncidenceError(f"grid value {value!r} is not an exact double.")
    return emitted
