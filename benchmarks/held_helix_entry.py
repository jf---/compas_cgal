"""The Held lane's helical entry: a HELIX into the stock, the pencil ramp-out, the hand-off into the published circle 0.

Rulings (USER 2026-11-04; ledger .superpowers/sdd/2026-10-13-product-pipeline/progress.md):

* The entry is its own operation stream ahead of the published circle 0:
  `compas_cgal.helix.Helix` (ramp turns + one flat closing turn, Q6), then
  the pipeline's tangent-pencil ramp-out loops (`compas_cgal.pipeline.pencil`,
  the `enter_first_cut` machinery), then the hand-off into circle 0. The
  published circles and connectors are untouched.
* Every junction of the entry is gated by Q4's ``theta_e = cap * pi / 45`` at
  gate station 1 (`Admissibility.loop_entry_within_rise`), the last one being
  the hand-off into circle 0.
* Cut depth is the case's declared `HeldReferenceCase.cut_depth` below the
  stock top `HELD_STOCK_TOP`; the ramp is the D1 default
  `DEFAULT_HELIX_PITCH_PER_RADIUS` (2 deg); the cutter must be centre-cutting
  (`NonCentreCuttingHelixError`, fail closed).

Two footprints, decided exactly on ``rho0 <= r`` (circle 0's radius vs the tool's):

* ``rho0 <= r``: the helix descends on circle 0 itself, injected exactly
  (`HeldCircleFootprint`: circle 0's own doubles, not re-authored). The bore
  ``disk(c0, rho0 + r)`` is then circle 0's swept disk exactly, no ramp-out
  loop is needed, and the hand-off is circle 0 at its own contact point.
* ``rho0 > r``: a helix on circle 0 would leave a core pin of radius
  ``rho0 - r``. The pencil root is `inscribed_root` (an `IncidentLoop` on the
  incidence grid whose swept disk lies exactly in circle 0's), the helix
  descends on the pencil footprint (``rho_h <= r``) and `RampOut` walks the
  pencil furthest-first, every member gated by the coverage bound, theta_e and
  the cap. The root itself is not cut: once the scan admits it, circle 0 is
  decided in its place (`EntryHandOffRefusedError`), after a straight link
  from the pencil's entry to circle 0's contact (sub-micron; gated by the cap).

CONTAINMENT. Every swept disk of the entry lies exactly in circle 0's swept
disk ``disk(c0, rho0 + r)``: on circle 0's footprint the bore IS that disk; on
the pencil, `inscribed_root` certifies the root's swept disk inside it and
every member (the helix footprint included) is internally tangent to the root
with a smaller radius, so its swept disk lies in the root's. The entry
therefore inherits circle 0's containment, the published path's own, and the
pencil is taken `TangentPencil.on_grid` instead of re-certifying the root's
disk against the part polygon (`TangentPencil.build`'s
`adaptive.containment.GougeContainment`, whose `ReachableDomain` takes minutes
on a 65-vertex Held polygon, measured 2026-11-04).

EXACTNESS. Circle 0's contact ``q`` is not exactly on its circle (``|q - c0|^2 -
rho0^2`` is ``-3.3e-18`` mm^2 on Figure 5). The exact entry ``P`` is the point
of the exact circle ``(c0, rho0)`` in the exact rational direction ``q - c0``, a
one-root point; `HeldCircleFootprint.certifies_incidence` certifies over the
rationals that ``q`` is ``P``'s nearest double, so the emitted start is ``q``
and nothing is snapped.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction

from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Point
from compas.geometry import Polygon

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal.engagement import _subtract_operation
from compas_cgal.helix import Helix
from compas_cgal.helix_units import DEFAULT_HELIX_PITCH_PER_RADIUS
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.incident_loop import INCIDENCE_GRID_EXPONENT
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import LoopEntryIncidenceError
from compas_cgal.incident_loop import chart_triple
from compas_cgal.incident_loop import exact_double
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.pipeline.emit import helix_operation
from compas_cgal.pipeline.emit import line_operation
from compas_cgal.pipeline.emit import loop_operation
from compas_cgal.pipeline.pencil import PencilConditions
from compas_cgal.pipeline.pencil import RampOut
from compas_cgal.pipeline.pencil import TangentPencil
from compas_cgal.post.errors import NonCentreCuttingHelixError
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import ClearanceZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

# The stock top of every Held case, world Z, mm: the datum of the declared cut
# depth and where the helix starts ramping (the tri-dexel replay's stock top is
# this datum shifted by its own STOCK_TOP_MM).
HELD_STOCK_TOP = ClearanceZ(Millimetre(0.0))

# The Held path is one closed path: every entry operation carries path index 0.
HELD_PATH_INDEX = 0

# Held circles turn counter-clockwise (`Figure5CounterclockwiseCircle.clockwise`).
HELD_CLOCKWISE = False

_GRID = 2**INCIDENCE_GRID_EXPONENT


class HeldHelixEntryError(BenchmarkError):
    """Base error of the Held helical entry."""


class NonCanonicalEntryReprError(HeldHelixEntryError):
    """The entry record's ``repr`` carries a memory address, so its digest would differ run to run."""


class DegenerateFootprintEntryError(HeldHelixEntryError):
    """Circle 0's contact coincides with its centre, so its entry direction is undefined."""


class InscribedRootError(HeldHelixEntryError):
    """No grid loop on circle 0's entry chart is inscribed in circle 0's swept disk."""


class EntryHandOffRefusedError(HeldHelixEntryError):
    """The hand-off into the published circle 0 fails Q4's theta_e at station 1, or the cap."""


class HandOffCoverageError(HeldHelixEntryError):
    """Circle 0's uncut core is not inside the disk the entry cleared (the pencil coverage bound at the hand-off)."""


class HandOffLinkRefusedError(HeldHelixEntryError):
    """The straight link from the pencil's entry to circle 0's contact exceeds the cap."""


def _sqrt_quotient_at_most(numerator: Fraction, radicand: Fraction, bound: Fraction) -> bool:
    """Exact ``numerator / sqrt(radicand) <= bound`` for rationals, ``radicand > 0``."""
    if numerator <= 0:
        return bound >= 0 or numerator * numerator >= bound * bound * radicand
    return bound > 0 and numerator * numerator <= bound * bound * radicand


def _nearest_double_of(coordinate: float, centre: Fraction, scaled: Fraction, radicand: Fraction) -> bool:
    """Exact: *coordinate* is a nearest double of ``centre + scaled / sqrt(radicand)`` (ties included).

    The rounding interval of a double ``x`` is ``[(x- + x)/2, (x + x+)/2]``
    with ``x-``/``x+`` its neighbour doubles (asymmetric at a binade edge).
    """
    value = Fraction(coordinate)
    low = (Fraction(math.nextafter(coordinate, -math.inf)) + value) / 2 - centre
    high = (Fraction(math.nextafter(coordinate, math.inf)) + value) / 2 - centre
    return _sqrt_quotient_at_most(scaled, radicand, high) and _sqrt_quotient_at_most(-scaled, radicand, -low)


@dataclass(frozen=True)
class HeldCircleFootprint:
    """A published Held circle injected exactly as a helix footprint and guide loop (`compas_cgal.helix.HelixFootprint`).

    ``centre`` and ``radius`` are the circle's own doubles; ``entry`` is its
    contact ``q``. The exact entry is ``P = c + rho (q - c) / |q - c|``, a
    one-root point; `certifies_incidence` certifies that ``q`` is its nearest
    double.

    Attributes:
        centre: The circle's centre.
        radius: The circle's guide radius.
        entry: The circle's contact, the nearest double of the exact entry.
    """

    centre: Point2[WorldXY]
    radius: GuideRadius
    entry: Point2[WorldXY]

    @classmethod
    def build(cls, circle: Figure5CounterclockwiseCircle) -> HeldCircleFootprint:
        """Inject *circle* exactly and certify its entry.

        Raises:
            DegenerateFootprintEntryError: The contact coincides with the centre.
            LoopEntryIncidenceError: The contact is not the nearest double of the exact entry.
        """
        footprint = cls(centre=circle.center, radius=circle.radius, entry=circle.contact_point)
        if footprint.centre == footprint.entry:
            raise DegenerateFootprintEntryError(f"circle contact {circle.contact_point!r} coincides with its centre.")
        if not footprint.certifies_incidence():
            raise LoopEntryIncidenceError(f"contact {circle.contact_point!r} is not the nearest double of the exact entry of circle ({circle.center!r}, {circle.radius!r}).")
        return footprint

    @property
    def entry_direction(self) -> tuple[Fraction, Fraction]:
        """The exact rational direction ``q - c`` of the exact entry."""
        return Fraction(self.entry.x) - Fraction(self.centre.x), Fraction(self.entry.y) - Fraction(self.centre.y)

    def certifies_incidence(self) -> bool:
        """Exact: ``entry`` is a nearest double, per coordinate, of ``P = c + rho u / |u|`` with ``u = q - c != 0``."""
        ux, uy = self.entry_direction
        radicand = ux * ux + uy * uy
        if radicand == 0:
            return False
        rho = Fraction(self.radius.value)
        return _nearest_double_of(self.entry.x, Fraction(self.centre.x), rho * ux, radicand) and _nearest_double_of(self.entry.y, Fraction(self.centre.y), rho * uy, radicand)

    def circle(self) -> Circle:
        """REPORTING: the circle as a compas `Circle` whose frame x-axis points at the entry."""
        ux, uy = (float(component) for component in self.entry_direction)
        length = math.hypot(ux, uy)
        ux, uy = ux / length, uy / length
        return Circle(self.radius.value, frame=Frame([self.centre.x, self.centre.y, 0.0], [ux, uy, 0.0], [-uy, ux, 0.0]))


def inscribed_root(circle: Figure5CounterclockwiseCircle) -> IncidentLoop:
    """The pencil root standing in for circle 0: the largest grid loop on its entry chart inscribed in its swept disk.

    `IncidentLoop.build` places the chart point and the grid centre nearest
    circle 0's; the radius is then the largest grid multiple ``s m / 2^K``
    with ``|c_R - c0| <= rho0 - rho_R`` exactly, so the root's swept disk (and
    every pencil member's) lies in circle 0's. Authoring, not a snap: circle 0
    is not changed and the root is never cut.

    Raises:
        InscribedRootError: No positive grid radius is inscribed.
    """
    nearest = IncidentLoop.build(circle.center, circle.radius, Radian(math.atan2(circle.contact_point.y - circle.center.y, circle.contact_point.x - circle.center.x)))
    a, b, m = chart_triple(nearest.chart, nearest.chart_numerator)
    steps = Fraction(nearest.radius.value) * _GRID / m
    if steps.denominator != 1:
        raise InscribedRootError(f"authored radius {nearest.radius.value!r} mm is off its chart grid.")
    rho0 = Fraction(circle.radius.value)
    dx = Fraction(nearest.centre.x) - Fraction(circle.center.x)
    dy = Fraction(nearest.centre.y) - Fraction(circle.center.y)
    count = int(steps)
    while count >= 1:
        gap = rho0 - Fraction(count * m, _GRID)
        if gap >= 0 and gap * gap >= dx * dx + dy * dy:
            break
        count -= 1
    if count < 1:
        raise InscribedRootError(f"no grid loop on chart ({nearest.chart}, {nearest.chart_numerator}) is inscribed in circle ({circle.center!r}, {circle.radius!r}).")
    cx, cy = Fraction(nearest.centre.x), Fraction(nearest.centre.y)
    root = IncidentLoop(
        centre=nearest.centre,
        radius=GuideRadius(Millimetre(exact_double(Fraction(count * m, _GRID)))),
        entry=Point2[WorldXY](Millimetre(exact_double(cx + Fraction(count * a, _GRID))), Millimetre(exact_double(cy + Fraction(count * b, _GRID)))),
        phase=nearest.phase,
        chart=nearest.chart,
        chart_numerator=nearest.chart_numerator,
    )
    if not (root.certifies_incidence() and root.certifies_chart()):
        raise LoopEntryIncidenceError(f"inscribed root {root!r} fails its incidence or chart certificate.")
    return root


def part_polygon(case: HeldReferenceCase) -> Polygon:
    """The case's material polygon P (``case.projection.points``) in world XY."""
    return Polygon([Point(float(point.x), float(point.y), 0.0) for point in case.projection.points])


def entry_admissibility(case: HeldReferenceCase) -> Admissibility:
    """The in-process verdicts of the entry: the case's cap, its centre domain, and Q4's declared entry-rise limit."""
    cap = Radian(math.radians(float(case.tea_cap)))
    return Admissibility.build(
        case.tool_radius,
        cap_ratio(cap),
        centre_domain(part_polygon(case), [], case.tool_radius),
        entry_rise_ratio=entry_rise_ratio(cap, EntryRampLength.default()),
    )


@dataclass(frozen=True)
class HeldHelixEntry:
    """The decided entry of one Held path, ahead of its published circle 0.

    Attributes:
        helix: The HELIX (ramp turns + one flat closing turn).
        members: The ramp-out loops in emission order (none on circle 0's own footprint).
        root: The inscribed pencil root (`inscribed_root`), or None on circle 0's own footprint.
        hand_off: The straight link ``(start, end)`` from the last entry exit to
            circle 0's contact, or None when they coincide.
        hand_off_operation: The link's label: CUT if it removes material, else LINK.
        circle_zero: Circle 0 as decided at the hand-off.
    """

    helix: Helix
    members: tuple[IncidentLoop, ...]
    root: IncidentLoop | None
    hand_off: tuple[Point2[WorldXY], Point2[WorldXY]] | None
    hand_off_operation: OperationType | None
    circle_zero: HeldCircleFootprint

    @classmethod
    def build(cls, case: HeldReferenceCase, circles: tuple[Figure5CounterclockwiseCircle, ...]) -> HeldHelixEntry:
        """Decide the entry into *circles[0]* on the case's stock, every junction gated by theta_e.

        Args:
            case: The reference case (tool, cap, declared cut depth and cutter).
            circles: The published circles in emission order.

        Returns:
            The entry.

        Raises:
            NonCentreCuttingHelixError: The case's cutter is not centre-cutting.
            DegenerateFootprintEntryError: Circle 0's contact is its centre.
            LoopEntryIncidenceError: Circle 0's contact is not its exact entry's nearest double.
            InscribedRootError: ``rho0 > r`` and no inscribed pencil root exists.
            RampOutStuckError: No pencil member admits the junction and the cap.
            HandOffLinkRefusedError: The link to circle 0's contact exceeds the cap.
            HandOffCoverageError: Circle 0's core is not inside the cleared disk.
            EntryHandOffRefusedError: Circle 0 fails theta_e at station 1 or the cap.
        """
        if not case.centre_cutting:
            raise NonCentreCuttingHelixError("a HELIX entry ramps on the end face: the Held case's cutter must be centre-cutting.")
        circle_zero = HeldCircleFootprint.build(circles[0])
        tool = case.tool_radius
        r = float(tool.value)
        plane = case.cut_depth.cut_plane(HELD_STOCK_TOP)
        admissibility = entry_admissibility(case)
        stock = Stock(part_polygon(case))
        members: list[IncidentLoop] = []
        root: IncidentLoop | None = None
        if Fraction(circle_zero.radius.value) <= Fraction(tool.value):
            helix = Helix.build(circle_zero, plane, DEFAULT_HELIX_PITCH_PER_RADIUS, HELD_CLOCKWISE)
            _subtract_operation(stock, helix_operation(helix, HELD_PATH_INDEX), r)
            last_exit, last_centre, last_radius = circle_zero.entry, circle_zero.centre, circle_zero.radius
        else:
            root = inscribed_root(circles[0])
            pencil = TangentPencil.on_grid(root)
            ramp = RampOut.build(pencil, tool, PencilConditions(admissibility, HELD_CLOCKWISE).ordered)
            helix = Helix.build(ramp.footprint, plane, DEFAULT_HELIX_PITCH_PER_RADIUS, HELD_CLOCKWISE)
            _subtract_operation(stock, helix_operation(helix, HELD_PATH_INDEX), r)
            last: IncidentLoop = ramp.footprint
            member = ramp.next_member(stock, last)
            while not ramp.is_root(member):
                _subtract_operation(stock, loop_operation(member, plane.cut_z, HELD_PATH_INDEX, HELD_CLOCKWISE), r)
                members.append(member)
                last = member
                member = ramp.next_member(stock, last)
            last_exit, last_centre, last_radius = last.entry, last.centre, last.radius
        hand_off: tuple[Point2[WorldXY], Point2[WorldXY]] | None = None
        hand_off_operation: OperationType | None = None
        if last_exit != circle_zero.entry:
            hand_off = (last_exit, circle_zero.entry)
            if not admissibility.straight_within_cap(stock, last_exit, circle_zero.entry):
                raise HandOffLinkRefusedError(f"the link {last_exit!r} -> {circle_zero.entry!r} into circle 0 exceeds the cap.")
            link = line_operation(last_exit, circle_zero.entry, plane.cut_z.value, plane.cut_z.value, OperationType.CUT, HELD_PATH_INDEX)
            hand_off_operation = OperationType.LINK if admissibility.removal(stock, link) is SweepRemoval.NOTHING else OperationType.CUT
            _subtract_operation(stock, link, r)
        _require_core_covered(circle_zero, last_centre, last_radius, tool.value)
        if not admissibility.loop_entry_within_rise(stock, circle_zero, HELD_CLOCKWISE):
            raise EntryHandOffRefusedError("the hand-off into circle 0 exceeds theta_e at gate station 1.")
        if not admissibility.loop_within_cap(stock, circle_zero, HELD_CLOCKWISE):
            raise EntryHandOffRefusedError("circle 0 exceeds the cap after the entry.")
        return cls(helix=helix, members=tuple(members), root=root, hand_off=hand_off, hand_off_operation=hand_off_operation, circle_zero=circle_zero)

    def operations(self, path_index: int = HELD_PATH_INDEX) -> tuple[ToolpathOperation, ...]:
        """The entry as toolpath operations: HELIX, the ramp-out CUT loops, the hand-off link.

        Rebuilt from the exact record on every call (operations carry arrays and
        are mutable; the record is the value).
        """
        cut = self.helix.bottom_z
        ops = [helix_operation(self.helix, path_index)]
        ops += [loop_operation(member, cut, path_index, self.helix.clockwise) for member in self.members]
        if self.hand_off is not None and self.hand_off_operation is not None:
            ops.append(line_operation(self.hand_off[0], self.hand_off[1], cut.value, cut.value, self.hand_off_operation, path_index))
        return tuple(ops)


def _require_core_covered(circle_zero: HeldCircleFootprint, centre: Point2[WorldXY], radius: GuideRadius, tool: float) -> None:
    """Exact: circle 0's uncut core ``disk(c0, rho0 - r)`` lies in the cleared ``disk(c, rho + r)`` (vacuous for ``rho0 <= r``).

    Raises:
        HandOffCoverageError: It does not.
    """
    r = Fraction(tool)
    core = Fraction(circle_zero.radius.value) - r
    if core <= 0:
        return
    slack = Fraction(radius.value) + r - core
    dx = Fraction(circle_zero.centre.x) - Fraction(centre.x)
    dy = Fraction(circle_zero.centre.y) - Fraction(centre.y)
    if slack < 0 or slack * slack < dx * dx + dy * dy:
        raise HandOffCoverageError(f"circle 0's core (radius {float(core)!r} mm) is not inside the disk cleared by the entry's last loop.")


def entry_digest(entry: HeldHelixEntry) -> str:
    """SHA-256 of the entry record's ``repr``: frozen dataclasses of doubles (shortest round-trip) and enums, so canonical.

    Raises:
        NonCanonicalEntryReprError: The ``repr`` carries a memory address.
    """
    text = repr(entry)
    if " at 0x" in text:
        raise NonCanonicalEntryReprError("the entry repr carries a memory address.")
    return hashlib.sha256(text.encode()).hexdigest()
