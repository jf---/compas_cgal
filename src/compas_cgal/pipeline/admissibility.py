"""In-process admissibility of loops and straight moves (Rulings 1.1-C8, 1.1-I17).

Every decision is an exact native predicate, asked through the decision-only
native queries (``engagement_*_exceeded``: the same stations and exact verdicts
as ``engagement_*_prefixes``, one bool, early exit); the reporting prefix
samples are never computed for a decision:

* the cap verdict is `_stock_2`'s ``cap_exceeded`` boolean on the exact
  arrangement, evaluated IN PROCESS (against the stock minus the motion's own
  swept prefix) at sampled stations -- and, for a loop, also AT REST (against
  the stock as it stands) on rho's probe ring, the frozen rho generator's own
  rule (`loop_within_cap_at_rest`, Ruling 1.1-C7; the one named exception
  to plan R1(a)'s in-process rule, User ruling (2026-10-16), docs/pipeline.md);
* containment of a straight is exact segment containment in the cutter-centre
  domain (the start's point location plus exact edge clearance);
* removal is the gate's own three-valued bracket (`compas_cgal.engagement`).

SAMPLED, NOT BOUNDED: each verdict is exact at the stations it is asked at;
nothing here bounds engagement between stations. Nothing here mutates the stock.

The cap crosses into exact-land once, in `cap_ratio`, as the rational chord
surrogate ``4 sin^2(cap / 2)`` computed by `_stock_2.cap_chord_ratio`; the
sub-ulp gap between it and the transcendental cap is API semantics. The
declared entry-rise limit of the helical entry's pencil junctions crosses the
same way, once, in `entry_rise_ratio`, and is injected beside the cap in
`Admissibility.build` (USER 2026-10-16 Q4).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field
from typing import NewType
from typing import Self

from compas.geometry import Polygon

from compas_cgal import _coverage_2
from compas_cgal import _stock_2
from compas_cgal.engagement import _operation_removal
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.errors import InvalidCapChordRatioError
from compas_cgal.pipeline.errors import InvalidEngagementCapError
from compas_cgal.pipeline.errors import InvalidEntryRiseChordRatioError
from compas_cgal.pipeline.errors import MissingEntryRiseLimitError
from compas_cgal.pipeline.ladder import WarmStart
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import _polygon_to_ccw_vertices
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

CapChordRatio = NewType("CapChordRatio", float)

# The chord surrogate ``4 sin^2(theta_e / 2)`` of the declared entry-rise limit
# (`entry_rise_ratio`); validated in `Admissibility.build` like the cap's.
EntryRiseChordRatio = NewType("EntryRiseChordRatio", float)

# The largest chord surrogate a cap in (0, pi] has: 4 sin^2(pi / 2) = 4, the
# half-turn cap.
MAX_CAP_CHORD_RATIO = 4.0

# The gate's samples per motion (benchmarks.survey.QUALITY_SAMPLES_PER_MOTION,
# held equal by tests/pipeline/test_admissibility.py; restated because
# pipeline/ imports no benchmarks module). Every in-process decision is taken
# at exactly the gate's own prefix parameters (Ruling 1.1-I4), so a motion the
# generator admits has no cap_exceeded sample in the gate: the two share one
# stock model and one station set. A sampling density, not a certificate.
GATE_SAMPLES_PER_MOTION = 45

# Stations on a full machining circle: the gate's, k / N for k < N, entry
# first; the seam (parameter 1) is the entry again and is not resampled.
LOOP_PROBE_COUNT_INPROCESS = GATE_SAMPLES_PER_MOTION


# The gate's junction station (benchmarks.quality_observations.JUNCTION_STATION,
# restated because pipeline/ imports no benchmarks module): the step a
# junction reads is the next motion's sample 1 of `GATE_SAMPLES_PER_MOTION`.
JUNCTION_STATION = 1


def loop_parameters() -> list[float]:
    """The prefix parameters a loop is decided at: the gate's ``k / N``, ``k = 0 .. N - 1``."""
    return [index / LOOP_PROBE_COUNT_INPROCESS for index in range(LOOP_PROBE_COUNT_INPROCESS)]


def motion_parameters() -> list[float]:
    """The prefix parameters an open motion (straight or arc) is decided at: the gate's ``k / N``, ``k = 0 .. N``."""
    return [index / GATE_SAMPLES_PER_MOTION for index in range(GATE_SAMPLES_PER_MOTION + 1)]


# rho's at-rest probe ring: the loop entry, then this many positions uniform in
# angle from the advance direction (engagement_toolpath.LOOP_PROBE_COUNT, held
# equal by test_against_rho). A sampling density, not a certificate.
LOOP_PROBE_COUNT_AT_REST = 32

# A full turn in degrees: the ring is built in degrees then converted, as rho's
# LOOP_PROBE_ANGLES_DEG is, so the positions are rho's to the bit.
FULL_TURN_DEG = 360.0

# Straight-move stations are spaced this many tool radii apart, ends included,
# so consecutive tool disks overlap by a full diameter (port of the
# chain-ordered generator's link probes).
LINK_PROBE_SPACING_TOOL_RADII = 1.0

# Fewest stations on any straight move: both ends.
MIN_LINK_PROBES = 2

# No gap closure: closure is the pessimism a continuous certifier needs between
# stations, and this makes no between-station claim.
GAP_CLOSE_RATIO = 0.0


def cap_ratio(cap: Radian) -> CapChordRatio:
    """The exact rational surrogate of an engagement cap.

    Args:
        cap: Engagement cap in ``(0, pi]``.

    Returns:
        ``4 sin^2(cap / 2)`` as computed by the native kernel.

    Raises:
        InvalidEngagementCapError: *cap* is NaN or outside ``(0, pi]``.
    """
    if not (cap > 0.0 and cap <= math.pi):
        raise InvalidEngagementCapError(f"engagement cap must be in (0, pi] rad; got {cap!r}.")
    return CapChordRatio(_stock_2.cap_chord_ratio(float(cap)))


def entry_rise_ratio(cap: Radian, ramp: EntryRampLength) -> EntryRiseChordRatio:
    """The exact rational surrogate of the declared entry-rise limit at the gate's junction station.

    The limit is ``theta_e = cap * s_e / L``: engagement may rise from the seam
    to the cap no faster than linearly over the declared ramp ``L``. It is
    evaluated at the gate's station 1 of the smallest legal loop (radius
    ``r``), ``s_e = 2 pi r / GATE_SAMPLES_PER_MOTION``, so
    ``theta_e = cap * 2 pi / (GATE_SAMPLES_PER_MOTION * L / r)`` -- ``cap * pi / 45``
    at the default ``L = 2r``. ``theta_e`` is computed once as a double and its
    chord surrogate injected exactly; the sub-ulp gap between the surrogate
    and the transcendental intent is this parameter's semantics. SAMPLED
    (Ruling 0.3-C6): it bounds the rise at that one station, nothing between.

    Args:
        cap: Engagement cap in ``(0, pi]``.
        ramp: The declared entry ramp.

    Returns:
        ``4 sin^2(theta_e / 2)`` as computed by the native kernel.

    Raises:
        InvalidEngagementCapError: *cap* is NaN or outside ``(0, pi]``.
        InvalidEntryRiseChordRatioError: ``theta_e`` exceeds ``pi`` (a ramp
            shorter than the station-1 arc).
    """
    if not (cap > 0.0 and cap <= math.pi):
        raise InvalidEngagementCapError(f"engagement cap must be in (0, pi] rad; got {cap!r}.")
    theta = cap * math.tau / (GATE_SAMPLES_PER_MOTION * ramp.tool_radii)
    if not theta <= math.pi:
        raise InvalidEntryRiseChordRatioError(f"entry-rise limit {theta!r} rad exceeds pi: ramp {ramp.tool_radii!r} tool radii is shorter than the station-1 arc.")
    return EntryRiseChordRatio(_stock_2.cap_chord_ratio(theta))


def centre_domain(polygon: Polygon, holes: list[Polygon], tool: ToolRadius) -> _coverage_2.CutterCentreDomain2:
    """The exact region a cutter of radius *tool* may be centred in.

    Args:
        polygon: Outer pocket boundary.
        holes: Islands.
        tool: Tool radius.

    Returns:
        The native cutter-centre domain.
    """
    return _coverage_2.CutterCentreDomain2.build(_polygon_to_ccw_vertices(polygon), [_polygon_to_ccw_vertices(hole) for hole in holes], float(tool.value))


def ring_positions(loop: IncidentLoop, advance: Direction2[WorldXY]) -> list[tuple[float, float]]:
    """rho's at-rest probe ring: the loop's exact entry, then `LOOP_PROBE_COUNT_AT_REST` positions uniform in angle from *advance*.

    Built in degrees then converted, as rho's ``LOOP_PROBE_ANGLES_DEG`` is, so
    the positions are rho's to the bit.

    Args:
        loop: The authored loop.
        advance: Direction of travel into the loop.

    Returns:
        ``1 + LOOP_PROBE_COUNT_AT_REST`` positions, the entry first.
    """
    positions: list[tuple[float, float]] = [(loop.entry.x, loop.entry.y)]
    for index in range(LOOP_PROBE_COUNT_AT_REST):
        angle = math.radians(FULL_TURN_DEG * index / LOOP_PROBE_COUNT_AT_REST)
        cosine, sine = math.cos(angle), math.sin(angle)
        ux, uy = advance.x * cosine - advance.y * sine, advance.x * sine + advance.y * cosine
        positions.append((loop.centre.x + loop.radius.value * ux, loop.centre.y + loop.radius.value * uy))
    return positions


def straight_probe_positions(a: Point2[WorldXY], b: Point2[WorldXY], tool: ToolRadius) -> list[tuple[float, float]]:
    """Evenly spaced stations along a straight move, *a* first and *b* last.

    No decision reads these any more: containment is exact segment
    containment (`Admissibility.straight_in_domain`).

    Args:
        a: Start.
        b: End.
        tool: Tool radius (sets the spacing).

    Returns:
        The stations; the ends are exactly *a* and *b*.
    """
    length = math.hypot(b.x - a.x, b.y - a.y)
    count = max(MIN_LINK_PROBES, math.ceil(length / (LINK_PROBE_SPACING_TOOL_RADII * tool.value)) + 1)
    last = count - 1
    inner = [(a.x + (b.x - a.x) * index / last, a.y + (b.y - a.y) * index / last) for index in range(1, last)]
    return [(a.x, a.y), *inner, (b.x, b.y)]


@dataclass(frozen=True)
class Admissibility:
    """Exact in-process cap, containment and removal verdicts for one tool and cap.

    Attributes:
        tool: Tool radius.
        cap_ratio: The cap's exact chord surrogate.
        domain: The cutter-centre domain (the containment half; swappable).
        entry_rise_ratio: The declared entry-rise limit's exact chord
            surrogate (`entry_rise_ratio`), or None when this admissibility
            declares none; only the helical entry's pencil junctions read it.
        ring_order: Warm-start evaluation order of the at-rest ring's
            positions (cost only; excluded from equality and hashing).
            ``dataclasses.replace`` and ``copy.copy`` share this one
            `WarmStart`; harmless, since it only orders probes.
    """

    tool: ToolRadius
    cap_ratio: CapChordRatio
    domain: _coverage_2.CutterCentreDomain2
    entry_rise_ratio: EntryRiseChordRatio | None = None
    ring_order: WarmStart = field(default_factory=WarmStart, compare=False, repr=False)

    @classmethod
    def build(
        cls,
        tool: ToolRadius,
        cap_ratio: CapChordRatio,
        domain: _coverage_2.CutterCentreDomain2,
        entry_rise_ratio: EntryRiseChordRatio | None = None,
    ) -> Self:
        """Bundle the verdict inputs, refusing a chord ratio no angle in ``(0, pi]`` has.

        `CapChordRatio` and `EntryRiseChordRatio` are bare `NewType`s, so the
        ranges are checked here and a hand-made ratio cannot bypass what
        `cap_ratio` / `entry_rise_ratio` guarantee.

        Args:
            tool: Tool radius (positive by its type).
            cap_ratio: From `cap_ratio`.
            domain: From `centre_domain`.
            entry_rise_ratio: From `entry_rise_ratio`, injected once here beside
                the cap; None declares no entry-rise limit, and a pencil
                junction asked of it raises `MissingEntryRiseLimitError`.

        Returns:
            The admissibility.

        Raises:
            InvalidCapChordRatioError: *cap_ratio* is not in ``(0, 4]`` (or NaN).
            InvalidEntryRiseChordRatioError: *entry_rise_ratio* is not in ``(0, 4]`` (or NaN).
        """
        if not (cap_ratio > 0.0 and cap_ratio <= MAX_CAP_CHORD_RATIO):
            raise InvalidCapChordRatioError(f"a cap chord ratio lies in (0, {MAX_CAP_CHORD_RATIO}]; got {cap_ratio!r}.")
        if entry_rise_ratio is not None and not (entry_rise_ratio > 0.0 and entry_rise_ratio <= MAX_CAP_CHORD_RATIO):
            raise InvalidEntryRiseChordRatioError(f"an entry-rise chord ratio lies in (0, {MAX_CAP_CHORD_RATIO}]; got {entry_rise_ratio!r}.")
        return cls(tool=tool, cap_ratio=cap_ratio, domain=domain, entry_rise_ratio=entry_rise_ratio)

    def loop_within_cap(self, stock: Stock, loop: IncidentLoop, clockwise: bool) -> bool:
        """No station of the full loop exceeds the cap, in process, on the authored radius.

        Args:
            stock: Stock before the loop.
            loop: The authored loop; station 0 is its exact entry.
            clockwise: Turn direction.

        Returns:
            True iff every sampled station is under the cap.
        """
        return self._loop_stations_within_cap(stock, loop, clockwise, loop_parameters())

    def loop_within_cap_past_entry(self, stock: Stock, loop: IncidentLoop, clockwise: bool) -> bool:
        """`loop_within_cap` at its stations ``k >= 1``: station 0 is left to the at-rest ring's entry probe.

        Station 0 (the empty prefix) is exactly `engagement_at` at the loop's
        entry on the same stock -- the first position `loop_within_cap_at_rest`
        probes. So ``ring and loop_within_cap_past_entry`` decides exactly what
        ``ring and loop_within_cap`` decides, in either order, and the shared
        query runs once. Valid ONLY AND-ed with the ring on the same *stock*
        (`LoopConditions.ordered`); alone it is not the loop's decision.

        Args:
            stock: Stock before the loop (the ring's stock).
            loop: The authored loop.
            clockwise: Turn direction.

        Returns:
            True iff every sampled station past the entry is under the cap.
        """
        return self._loop_stations_within_cap(stock, loop, clockwise, loop_parameters()[1:])

    def loop_entry_within_rise(self, stock: Stock, loop: IncidentLoop, clockwise: bool) -> bool:
        """The loop's gate junction station (`JUNCTION_STATION`) is under the declared entry-rise limit, in process.

        The same station the gate's step reads at a junction into this loop,
        decided by the same exact native predicate as the cap, against the
        entry-rise surrogate: it bounds the loop's own engagement there,
        ``|phi_1| <= theta_e``. The gate's step cell is
        ``|phi_1(next) - phi_exit(prev)|``, so it is bounded by
        ``max(theta_e, phi_exit(prev))``, not certified ``<= theta_e`` by this
        decision. SAMPLED: one station, nothing between. Applied only to the
        helical entry's pencil members (USER 2026-10-16 Q4).

        Args:
            stock: Stock before the loop.
            loop: The authored loop.
            clockwise: Turn direction.

        Returns:
            True iff the junction station is under the limit.

        Raises:
            MissingEntryRiseLimitError: This admissibility declares no limit.
        """
        if self.entry_rise_ratio is None:
            raise MissingEntryRiseLimitError("a pencil junction needs a declared entry-rise limit; build the admissibility with entry_rise_ratio.")
        return not self._loop_stations_exceed(stock, loop, clockwise, [loop_parameters()[JUNCTION_STATION]], float(self.entry_rise_ratio))

    def _loop_stations_within_cap(self, stock: Stock, loop: IncidentLoop, clockwise: bool, parameters: list[float]) -> bool:
        return not self._loop_stations_exceed(stock, loop, clockwise, parameters, float(self.cap_ratio))

    def _loop_stations_exceed(self, stock: Stock, loop: IncidentLoop, clockwise: bool, parameters: list[float], ratio: float) -> bool:
        return _stock_2.engagement_circle_exceeded(
            stock.raw,
            float(loop.centre.x),
            float(loop.centre.y),
            float(loop.radius.value),
            float(loop.entry.x),
            float(loop.entry.y),
            clockwise,
            parameters,
            float(self.tool.value),
            ratio,
            GAP_CLOSE_RATIO,
        )

    def loop_within_cap_at_rest(self, stock: Stock, loop: IncidentLoop, advance: Direction2[WorldXY]) -> bool:
        """No position of rho's at-rest probe ring exceeds the cap (the frozen rho generator's loop rule).

        The positions are the loop's exact entry and `LOOP_PROBE_COUNT_AT_REST`
        points of the authored circle uniform in angle from *advance*, each
        decided against *stock* as it stands -- without the loop's own swept
        prefix and without the bridge into it. Material at a position only
        shrinks as the prefix is removed, so an at-rest verdict is never more
        permissive than the in-process one at the same position: this is the
        conservative half of the pair, and the rule the reference column
        (`rho_regulated`) advances by.

        EVALUATION ORDER (Task 1.2 lever L6). The verdict is "no position
        exceeds", a conjunction of the same exact per-position
        ``engagement_at`` verdicts over the fixed `ring_positions` set, so the
        positions are asked in `ring_order`'s warm-start order (the position
        that refused last first, then the natural order): the same positions,
        each asked the same way, at most once; an admitted ring has asked them
        all. VERDICTS are order-free. Native INPUT checks (finiteness, tool,
        ratio, gap) are position-independent. Position-dependent INVARIANT
        errors (``RimSpanNormalizationError``, engagement_2.cpp:527;
        ``ExactOneRootCoordinateMismatchError`` /
        ``InvalidMixedRadicalRootError``, audit_exact_station_2.cpp:374,394)
        may surface on a different probe or call than under HEAD's natural
        order, just as HEAD's own short-circuit already leaves later positions
        unasked; making errors order-free would require asking every position.

        Args:
            stock: Stock before the bridge and the loop.
            loop: The authored loop.
            advance: Direction of travel into the loop (origin centre to its centre).

        Returns:
            True iff every probed position is under the cap.
        """
        positions = ring_positions(loop, advance)
        raw, tool, ratio = stock.raw, float(self.tool.value), float(self.cap_ratio)

        def exceeds(index: int) -> bool:
            x, y = positions[index]
            return bool(_stock_2.engagement_at(raw, x, y, tool, ratio, GAP_CLOSE_RATIO)[2])

        return self.ring_order.first_refusal(len(positions), exceeds) is None

    def straight_within_cap(self, stock: Stock, a: Point2[WorldXY], b: Point2[WorldXY]) -> bool:
        """No station of the straight move exceeds the cap, in process.

        A zero-length move (``a == b`` exactly) is no motion and is admissible.

        Args:
            stock: Stock before the move.
            a: Start.
            b: End.

        Returns:
            True iff every sampled station is under the cap.
        """
        if a == b:
            return True
        return not _stock_2.engagement_line_exceeded(
            stock.raw,
            float(a.x),
            float(a.y),
            float(b.x),
            float(b.y),
            motion_parameters(),
            float(self.tool.value),
            float(self.cap_ratio),
            GAP_CLOSE_RATIO,
        )

    def arc_within_cap(self, stock: Stock, centre: Point2[WorldXY], start: Point2[WorldXY], end: Point2[WorldXY], clockwise: bool) -> bool:
        """No station of the arc exceeds the cap, in process (Ruling 1.1-I2).

        The native query the survey asks of an arc: the sweep from *start* to
        *end* about *centre* in the travel direction, at the gate's parameters.
        Callers pass the endpoints of the operation they will emit.

        Args:
            stock: Stock before the arc.
            centre: Arc centre.
            start: Where the travel begins.
            end: Where it ends.
            clockwise: Travel direction.

        Returns:
            True iff every sampled station is under the cap.
        """
        return not _stock_2.engagement_arc_exceeded(
            stock.raw,
            float(centre.x),
            float(centre.y),
            float(start.x),
            float(start.y),
            float(end.x),
            float(end.y),
            clockwise,
            motion_parameters(),
            float(self.tool.value),
            float(self.cap_ratio),
            GAP_CLOSE_RATIO,
        )

    def straight_in_domain(self, a: Point2[WorldXY], b: Point2[WorldXY]) -> bool:
        """The whole straight move lies in the cutter-centre domain (exact, not sampled).

        Decided by `CutterCentreDomain2.contains_segment`: the start is a legal
        centre and the segment keeps squared distance >= r^2 from every pocket
        and island edge, over exact rationals of the injected doubles. A chord
        whose ends are both legal can still cut a reflex corner's r-disk; this
        refuses it.

        Args:
            a: Start.
            b: End.

        Returns:
            True iff every point of the segment is a legal tool centre.
        """
        return bool(self.domain.contains_segment(float(a.x), float(a.y), float(b.x), float(b.y)))

    def removal(self, stock: Stock, op: ToolpathOperation) -> SweepRemoval:
        """What *op* would remove from *stock*: the gate's own exact bracket.

        Args:
            stock: Stock before the operation.
            op: The operation.

        Returns:
            REMOVES, WITHIN_BAND or NOTHING.
        """
        return _operation_removal(stock, op, float(self.tool.value))
