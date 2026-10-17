"""Exact tangent-pencil authoring at a root loop's entry, and its furthest-first ramp-out.

Members of the pencil at ``P`` share ``P`` and the root's entry direction.
With the reduced Pythagorean direction ``(a, b, m)`` of ``P - c_root``
(``a^2 + b^2 = m^2``) and the root at ``s_root`` grid steps
(``P - c_root = s_root (a, b) / 2^K``, ``rho_root = s_root m / 2^K``), member
``s`` is ``c = P - s (a, b) / 2^K``, ``rho = s m / 2^K``: internally tangent to
the root at ``P`` for ``0 < s <= s_root``. Every coordinate is a ``2^-K`` grid
value, hence an exact double, and incidence is certified over the rationals.

The direction is read from the root's EXACT geometry, not from its chart
fields: ``(a, b) / m`` is ``(P - c_root) / rho_root`` in lowest terms, which
equals the atlas triple of `compas_cgal.incident_loop` for every loop that
module authors (its triple is reduced, so ``gcd(a, b, m) == 1``). This keeps
the module on public names only; `incident_loop` exposes no public chart-triple
helper.

THE RAMP-OUT (Task 1.2, N1). From the helix footprint to the root, one member
per `ladder.furthest_first` scan over the pencil's rungs (the root first, then
down by `RampOut.rung_steps`), never a bisection: engagement is not monotone in
the member radius. Conditions per candidate, cheapest first:

1. `PencilStepCondition` (P-C, prepended by the instance, Ruling
   1.1-B-boundinvariant): ``rho_next - rho_prev <= r``, exact in `Fraction`.
   Coverage: member ``k+1``'s uncut core ``disk(c_{k+1}, rho_{k+1} - r)`` lies in
   the cleared ``disk(c_k, rho_k + r)`` iff
   ``|c_{k+1} - c_k| + rho_{k+1} - r <= rho_k + r``, and on a pencil
   ``|c_{k+1} - c_k| = rho_{k+1} - rho_k``; the bore is member 0.
2. The junction (`PencilConditions.junction_admissible`): the candidate's gate
   station 1 is under the declared entry-rise limit (USER 2026-10-16 Q4).
3. The full cap (`PencilConditions.loop_admissible`): every in-process station
   ``k/45``, ``k >= 0`` (`Admissibility.loop_within_cap`); station 0 is asked
   here because the ramp-out has no at-rest ring to supply it, and the ring is
   never used (the AR exception stays the narrow walk's alone).

L1 (closed-form scan start): P-C is linear in ``s``, so the rungs it refuses
are exactly a prefix of the scan; the scan starts at the first rung it admits.
The bound stays prepended, so the start is a cost device, never a decision.

DOMAIN (P-D). The conditions decide nothing about the cutter-centre domain.
Members are internally tangent at ``P`` with ``rho_k <= rho_root``, so
``|c_k - c_root| = rho_root - rho_k`` and each member's swept disk
``disk(c_k, rho_k + r)`` lies in the root's swept disk
``disk(c_root, rho_root + r)``. `TangentPencil.build` certifies that ONE disk
exactly in the design pocket (`adaptive.containment.GougeContainment`, the
entry-disk certificate probe P6 used) and raises `HoledDomainPencilError`
otherwise: an island inside the pencil is refused by name until Task 1.7
supplies a per-member containment condition.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import overload

from compas.geometry import Polygon

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.incident_loop import INCIDENCE_GRID_EXPONENT
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.incident_loop import IncidentLoopRadiusBelowQuantumError
from compas_cgal.incident_loop import LoopEntryIncidenceError
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.errors import HoledDomainPencilError
from compas_cgal.pipeline.errors import RampOutStuckError
from compas_cgal.pipeline.ladder import CandidateCondition
from compas_cgal.pipeline.ladder import NoAdmissibleCandidate
from compas_cgal.pipeline.ladder import furthest_first
from compas_cgal.stock import Stock
from compas_cgal.units import EntryRadius
from compas_cgal.units import GuideRadius
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

_GRID = 2**INCIDENCE_GRID_EXPONENT

# Ramp-out rung spacing as a fraction of the coverage window: P-C admits a rise of
# at most r, and the furthest-first scan resolves that window into 100 rungs, so an
# admitted member is at most 1 % of r short of the furthest admissible one on the
# rung grid. Resolution intent, not a tolerance (it decides nothing; every rung is
# decided exactly). At r = 1 mm it is the 0.01 mm rung of probes P2 / te. One grid
# step (2^-40 m) is not usable: a refusal run of ~0.05 mm would be ~5e10 rungs.
PENCIL_RUNG_TOOL_RADII = Fraction(1, 100)


class RootOffChartGridError(ValueError):
    """A root loop is not an integer number of grid steps along its entry direction, so no pencil exists at its entry."""


class HelixFootprintExceedsRootError(ValueError):
    """The floored footprint radius exceeds the root radius: the footprint is not inside the root's pencil."""


class NotAPencilMemberError(ValueError):
    """A loop's radius is not an integer number of the pencil's grid steps, so it is no member of the pencil."""


@dataclass(frozen=True)
class TangentPencil:
    """The pencil of loops internally tangent to a root loop at its exact entry ``P``.

    Attributes:
        root: The root loop.
        a: Entry-direction numerator along x.
        b: Entry-direction numerator along y.
        hypotenuse: ``m``, with ``a^2 + b^2 = m^2``; one grid step of radius is ``m / 2^K``.
        root_steps: ``s_root``, the root's radius in grid steps.
    """

    root: IncidentLoop
    a: int
    b: int
    hypotenuse: int
    root_steps: int

    @classmethod
    def build(cls, root: IncidentLoop, tool: ToolRadius, pocket: GougeContainment) -> TangentPencil:
        """The pencil at *root*'s entry, its swept disk certified exactly in the design pocket (P-D).

        Args:
            root: An incidence-certified loop.
            tool: Tool radius ``r``.
            pocket: The design pocket's exact containment authority (`entry_containment`).

        Returns:
            The pencil.

        Raises:
            LoopEntryIncidenceError: *root* fails its incidence certificate.
            RootOffChartGridError: ``rho_root 2^K / m`` is not an integer.
            HoledDomainPencilError: ``disk(c_root, rho_root + r)`` is not in the pocket.
        """
        pencil = cls.on_grid(root)
        disk = Fraction(root.radius.value) + Fraction(tool.value)
        try:
            pocket.certify_entry_disk(center=root.centre, radius=EntryRadius.build(_exact_float(disk)), tool_radius=tool)
        except GougeContainmentError as error:
            raise HoledDomainPencilError(f"the pencil's swept disk (centre {root.centre!r}, radius {float(disk)!r} mm) is not in the design pocket: {error}") from error
        return pencil

    @classmethod
    def on_grid(cls, root: IncidentLoop) -> TangentPencil:
        """The pencil's grid at *root*'s entry, WITHOUT the pocket certificate (footprint arithmetic only; never cut from).

        Args:
            root: An incidence-certified loop.

        Returns:
            The pencil grid.

        Raises:
            LoopEntryIncidenceError: *root* fails its incidence certificate.
            RootOffChartGridError: ``rho_root 2^K / m`` is not an integer.
        """
        if not root.certifies_incidence():
            raise LoopEntryIncidenceError(f"root {root!r} fails its incidence certificate.")
        a, b, m = _entry_direction(root)
        root_steps = Fraction(root.radius.value) * _GRID / m
        if root_steps.denominator != 1:
            raise RootOffChartGridError(f"root radius {root.radius.value!r} mm is not an integer grid step {m}/2^{INCIDENCE_GRID_EXPONENT} mm.")
        return cls(root=root, a=a, b=b, hypotenuse=m, root_steps=int(root_steps))

    def member(self, steps: int) -> IncidentLoop:
        """Member *steps*: ``c = P - steps (a, b) / 2^K``, ``rho = steps m / 2^K``, entry ``P``.

        Args:
            steps: Grid steps, ``0 < steps <= root_steps``.

        Returns:
            The member, exactly incident at ``P``.

        Raises:
            LoopEntryIncidenceError: A member value is not an exact double, or
                the member fails its certificate (broken invariant).
        """
        px, py = Fraction(self.root.entry.x), Fraction(self.root.entry.y)
        loop = IncidentLoop(
            centre=Point2[WorldXY](Millimetre(_exact_float(px - Fraction(steps * self.a, _GRID))), Millimetre(_exact_float(py - Fraction(steps * self.b, _GRID)))),
            radius=GuideRadius(Millimetre(_exact_float(Fraction(steps * self.hypotenuse, _GRID)))),
            entry=self.root.entry,
            phase=self.root.phase,
            chart=self.root.chart,
            chart_numerator=self.root.chart_numerator,
        )
        if not loop.certifies_incidence():
            raise LoopEntryIncidenceError(f"pencil member {steps} entry {loop.entry!r} is not exactly on its circle.")
        return loop

    def steps_of(self, member: IncidentLoop) -> int:
        """The grid steps of a member of this pencil (its radius over one step, an integer).

        Args:
            member: A loop authored by `member`.

        Returns:
            Its ``steps``.

        Raises:
            NotAPencilMemberError: The radius is off the pencil's grid.
        """
        steps = Fraction(member.radius.value) * _GRID / self.hypotenuse
        if steps.denominator != 1:
            raise NotAPencilMemberError(f"radius {member.radius.value!r} mm is not an integer grid step {self.hypotenuse}/2^{INCIDENCE_GRID_EXPONENT} mm.")
        return int(steps)

    def footprint(self, tool: ToolRadius) -> IncidentLoop:
        """The helix footprint: the member with the largest radius ``<= r``, ``s = floor(r 2^K / m)``.

        ``round`` would put ``rho_h`` above ``r`` on 470 of 1025 chart-0
        numerators (probe P4), leaving an uncut core pin of radius
        ``rho_h - r`` inside the bore; with the floor the full bottom turn
        sweeps exactly ``disk(c_h, rho_h + r)``.

        Args:
            tool: Tool radius ``r``.

        Returns:
            The footprint loop, exactly incident at ``P``.

        Raises:
            IncidentLoopRadiusBelowQuantumError: ``r`` is below one grid step.
            HelixFootprintExceedsRootError: ``s > s_root``.
        """
        steps = Fraction(tool.value) * _GRID // self.hypotenuse
        if steps < 1:
            raise IncidentLoopRadiusBelowQuantumError(f"tool radius {tool.value!r} mm is below the grid step {self.hypotenuse}/2^{INCIDENCE_GRID_EXPONENT} mm.")
        if steps > self.root_steps:
            raise HelixFootprintExceedsRootError(f"footprint radius {float(Fraction(steps * self.hypotenuse, _GRID))!r} mm exceeds root radius {self.root.radius.value!r} mm.")
        return self.member(int(steps))


def helix_footprint(root: IncidentLoop, tool: ToolRadius) -> IncidentLoop:
    """The helix footprint of *root*'s pencil (`TangentPencil.footprint`).

    Args:
        root: The root loop the entry ramps out to.
        tool: Tool radius ``r``.

    Returns:
        The footprint loop, exactly incident at ``root.entry``.

    Raises:
        LoopEntryIncidenceError: *root* fails its incidence certificate, or a
            member value is not an exact double (broken invariant).
        RootOffChartGridError: ``rho_root 2^K / m`` is not an integer.
        IncidentLoopRadiusBelowQuantumError: ``r`` is below one grid step.
        HelixFootprintExceedsRootError: ``s > s_root``.
    """
    return TangentPencil.on_grid(root).footprint(tool)


def entry_containment(polygon: Polygon, holes: list[Polygon], tool: ToolRadius) -> GougeContainment:
    """The design pocket's exact containment authority for a pencil (P-D).

    Args:
        polygon: Outer pocket boundary in world XY.
        holes: Islands.
        tool: Tool radius.

    Returns:
        The containment authority over the pocket's reachable domain.
    """
    domain = ReachableDomain.build(
        design_boundary=CanonicalRingV1.build_outer([Point2[WorldXY].build(float(x), float(y)) for x, y, _z in polygon.points]),
        holes=[CanonicalRingV1.build_hole([Point2[WorldXY].build(float(x), float(y)) for x, y, _z in hole.points]) for hole in holes],
        tool_radius=tool,
    )
    return GougeContainment.build(domain)


@dataclass(frozen=True)
class PencilStepCondition:
    """P-C, the pencil coverage bound: a member rises at most ``r`` above the one before, exact.

    Attributes:
        tool: Tool radius ``r``.
    """

    tool: ToolRadius

    def __call__(self, stock: Stock, origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        """``rho_candidate - rho_origin <= r`` over the rationals (*stock* unused)."""
        del stock
        return Fraction(candidate.radius.value) - Fraction(origin.radius.value) <= Fraction(self.tool.value)


@dataclass(frozen=True)
class PencilConditions:
    """The ramp-out's injected conditions after the coverage bound, both in process.

    Attributes:
        admissibility: The verdicts (with the declared entry-rise limit).
        clockwise: The root's turn direction, every member's.
    """

    admissibility: Admissibility
    clockwise: bool

    def junction_admissible(self, stock: Stock, origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        """The candidate's gate station 1 is under the declared entry-rise limit (*origin* unused).

        Raises:
            MissingEntryRiseLimitError: The admissibility declares no limit.
        """
        del origin
        return self.admissibility.loop_entry_within_rise(stock, candidate, self.clockwise)

    def loop_admissible(self, stock: Stock, origin: IncidentLoop, candidate: IncidentLoop, /) -> bool:
        """Every in-process station of the candidate's full loop is under the cap (*origin* unused)."""
        del origin
        return self.admissibility.loop_within_cap(stock, candidate, self.clockwise)

    @property
    def ordered(self) -> tuple[CandidateCondition[IncidentLoop], ...]:
        """Junction (one station) before the full cap (45 stations): cheapest refusal first; the conjunction decides."""
        return (self.junction_admissible, self.loop_admissible)


class _Rungs(Sequence[IncidentLoop]):
    """The scan's candidates, authored lazily: member ``steps[i]``."""

    def __init__(self, pencil: TangentPencil, steps: range) -> None:
        self._pencil = pencil
        self._steps = steps

    def __len__(self) -> int:
        return len(self._steps)

    @overload
    def __getitem__(self, index: int) -> IncidentLoop: ...

    @overload
    def __getitem__(self, index: slice) -> Sequence[IncidentLoop]: ...

    def __getitem__(self, index: int | slice) -> IncidentLoop | Sequence[IncidentLoop]:
        if isinstance(index, slice):
            return _Rungs(self._pencil, self._steps[index])
        return self._pencil.member(self._steps[index])


@dataclass(frozen=True)
class RampOut:
    """The furthest-first ramp-out from the helix footprint to the root, one member per scan.

    Attributes:
        pencil: The root's pencil.
        tool: Tool radius ``r``.
        footprint: The helix footprint, member 0.
        rung_steps: Grid steps between consecutive rungs, ``ceil(PENCIL_RUNG_TOOL_RADII r 2^K / m)``.
        conditions: The coverage bound, then the injected conditions, in evaluation order.
    """

    pencil: TangentPencil
    tool: ToolRadius
    footprint: IncidentLoop
    rung_steps: int
    conditions: tuple[CandidateCondition[IncidentLoop], ...]

    @classmethod
    def build(cls, pencil: TangentPencil, tool: ToolRadius, conditions: tuple[CandidateCondition[IncidentLoop], ...]) -> RampOut:
        """The ramp-out over *pencil*, prepending `PencilStepCondition` to *conditions*.

        Args:
            pencil: The root's pencil.
            tool: Tool radius.
            conditions: Injected conditions in evaluation order (`PencilConditions.ordered`).

        Returns:
            The ramp-out.

        Raises:
            IncidentLoopRadiusBelowQuantumError: ``r`` is below one grid step.
            HelixFootprintExceedsRootError: The footprint exceeds the root.
        """
        return cls(
            pencil=pencil,
            tool=tool,
            footprint=pencil.footprint(tool),
            rung_steps=math.ceil(PENCIL_RUNG_TOOL_RADII * Fraction(tool.value) * _GRID / pencil.hypotenuse),
            conditions=(PencilStepCondition(tool), *conditions),
        )

    def is_root(self, member: IncidentLoop) -> bool:
        """Whether *member* is the root, the ramp-out's last member."""
        return self.pencil.steps_of(member) == self.pencil.root_steps

    def next_member(self, stock: Stock, previous: IncidentLoop) -> IncidentLoop:
        """The furthest admissible member after *previous* on *stock* (one furthest-first scan).

        Rungs run from the root down by `rung_steps`, strictly above
        *previous*; the scan starts at the first rung the coverage bound
        admits (L1), the bound itself still asked first on every candidate.

        Args:
            stock: The stock after the helix and every member before.
            previous: The member the tool is on (the footprint first).

        Returns:
            The member; `is_root` tells whether the ramp-out is complete.

        Raises:
            RampOutStuckError: No rung passes every condition.
        """
        previous_steps = self.pencil.steps_of(previous)
        reach = Fraction(self.tool.value) * _GRID // self.pencil.hypotenuse
        rungs = range(self.pencil.root_steps, previous_steps, -self.rung_steps)
        overshoot = self.pencil.root_steps - (previous_steps + int(reach))
        start = 0 if overshoot <= 0 else -(-overshoot // self.rung_steps)
        got = furthest_first(stock, previous, _Rungs(self.pencil, rungs[start:]), self.conditions)
        if isinstance(got, NoAdmissibleCandidate):
            raise RampOutStuckError(
                f"no pencil member above radius {previous.radius.value!r} mm admits the junction and the cap "
                f"({got.scanned} rung(s) scanned toward root radius {self.pencil.root.radius.value!r} mm)."
            )
        return got


def _entry_direction(root: IncidentLoop) -> tuple[int, int, int]:
    """The reduced integer direction ``(a, b, m)`` of ``P - c_root``, ``a^2 + b^2 = m^2``.

    Args:
        root: An incidence-certified loop.

    Returns:
        ``(a, b, m)`` with ``(a, b) / m == (P - c_root) / rho_root`` and ``m``
        the least common denominator.
    """
    radius = Fraction(root.radius.value)
    cos = (Fraction(root.entry.x) - Fraction(root.centre.x)) / radius
    sin = (Fraction(root.entry.y) - Fraction(root.centre.y)) / radius
    m = math.lcm(cos.denominator, sin.denominator)
    return int(cos * m), int(sin * m), m


def _exact_float(value: Fraction) -> float:
    """The double equal to a grid rational.

    Args:
        value: A ``2^-K`` grid rational within the supported magnitude.

    Returns:
        ``float(value)``.

    Raises:
        LoopEntryIncidenceError: *value* is not exactly a double.
    """
    emitted = float(value)
    if Fraction(emitted) != value:
        raise LoopEntryIncidenceError(f"grid value {value!r} is not an exact double.")
    return emitted
