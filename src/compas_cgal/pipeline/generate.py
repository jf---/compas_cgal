"""The narrow-regime pipeline as typed stages: guide, entry, walk, next-chain choice, exit.

`generate_narrow` machines a pocket with ONE plunge into material where the
guide allows it: it walks a skeleton chain with the narrow ladder, then reaches
the next chain at depth through stock earlier chains already cleared (a link
from `compas_cgal.pipeline.linking`), and plunges again only where no
admissible link exists -- recording every such decision in the ledger.

The stages are separate functions so later tasks replace one without touching
the others: `enter` (plunge + root loop; 1.2 replaces it with a helix),
`walk` (one chain from a given station, callable per component), `choose_next`
(the next chain and its link), `link` (emit a chosen link). The generator and
the gate share ONE stock model: every emitted cut-plane motion depletes the
stock through `compas_cgal.engagement._subtract_operation`, the gate's replay
dispatch, and every removal verdict is the gate's own bracket.

Labels follow the gate's rule (a) honestly (Ruling 1.1-C17): a cut-plane motion
is emitted as `LINK` only when its exact removal bracket is NOTHING; a motion
that removes, or may remove (WITHIN_BAND), is a `CUT`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field

from compas.geometry import Polygon

from compas_cgal.engagement import _subtract_operation
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.direction import Direction2
from compas_cgal.pipeline.emit import arc_operation
from compas_cgal.pipeline.emit import line_operation
from compas_cgal.pipeline.emit import loop_operation
from compas_cgal.pipeline.emit import tessellate
from compas_cgal.pipeline.errors import CoincidentStationsError
from compas_cgal.pipeline.errors import EmptyNarrowPathError
from compas_cgal.pipeline.errors import InvalidClearanceHeightError
from compas_cgal.pipeline.errors import MixedTurnChainError
from compas_cgal.pipeline.guide import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.pipeline.ladder import CandidateCondition
from compas_cgal.pipeline.linking import ArcPiece
from compas_cgal.pipeline.linking import BeltLink
from compas_cgal.pipeline.linking import LinePiece
from compas_cgal.pipeline.linking import Link
from compas_cgal.pipeline.linking import LinkPiece
from compas_cgal.pipeline.linking import NoLink
from compas_cgal.pipeline.linking import SameCircleArc
from compas_cgal.pipeline.linking import StraightLink
from compas_cgal.pipeline.linking import link_between
from compas_cgal.pipeline.narrow import Advance
from compas_cgal.pipeline.narrow import NarrowCandidate
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.narrow import advance_ladder
from compas_cgal.pipeline.narrow import advance_window
from compas_cgal.pipeline.narrow import orientations
from compas_cgal.pipeline.result import ForcedDecision
from compas_cgal.pipeline.result import ForcedDecisionLedger
from compas_cgal.pipeline.result import ForcedEntry
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.pipeline.result import LinkKind
from compas_cgal.pipeline.result import LinkRecord
from compas_cgal.pipeline.result import OmittedLoop
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.units import Clearance
from compas_cgal.units import ClearanceZ
from compas_cgal.units import CutZ
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

# The plunge's own sweep as the gate decides it: a disk, i.e. a circle sweep of
# guide radius zero (benchmarks.depletion.PLUNGE_GUIDE_RADIUS; the same
# predicate, restated because pipeline/ imports no benchmarks module).
PLUNGE_GUIDE_RADIUS = 0.0

# Most skeleton chains the guide may emit; more raises GuideTruncatedError.
DEFAULT_MAX_PASSES = 1000


@dataclass
class PathState:
    """The mutable state of one generation: stream, stock, records.

    Attributes:
        stock: The depleting stock, the gate's model of it.
        admissibility: Exact verdicts for this tool and cap.
        cut: Cut-plane height.
        clearance: Rapid-plane height.
        operations: The emitted stream.
        records: Forced decisions, in the order taken.
        omitted: Loops skipped as removing nothing.
        links: Cut-plane links emitted between loops.
        path_index: Path index of the current chain walk.
    """

    stock: Stock
    admissibility: Admissibility
    cut: CutZ
    clearance: ClearanceZ
    operations: list[ToolpathOperation] = field(default_factory=list)
    records: list[ForcedDecision] = field(default_factory=list)
    omitted: list[OmittedLoop] = field(default_factory=list)
    links: list[LinkRecord] = field(default_factory=list)
    path_index: int = 0

    @property
    def tool_radius(self) -> float:
        return float(self.admissibility.tool.value)

    def emit_motion(self, op: ToolpathOperation) -> SweepRemoval:
        """Label *op* by its exact removal, append it, and deplete the stock through the gate's dispatch."""
        removal = self.admissibility.removal(self.stock, op)
        op.operation = OperationType.LINK if removal is SweepRemoval.NOTHING else OperationType.CUT
        self.operations.append(op)
        _subtract_operation(self.stock, op, self.tool_radius)
        return removal


def chain_link(origin: IncidentLoop, target: IncidentLoop, clockwise: bool, target_clockwise: bool) -> Link:
    """The link between two stations of ONE chain, whose entries sit on the pair's external tangent.

    Precondition: a chain turns one way throughout. `link_between` refuses
    loops turning opposite ways (`NoLink`) before it reads ``tangent_entries``,
    and that is the only way it refuses tangent entries, so a `NoLink` here is
    exactly a breach of the precondition.

    Args:
        origin: The loop the tool is on.
        target: The loop advanced to.
        clockwise: *origin*'s turn direction.
        target_clockwise: *target*'s turn direction.

    Returns:
        The link.

    Raises:
        MixedTurnChainError: The two stations turn opposite ways.
    """
    link = link_between(origin, target, clockwise, target_clockwise, tangent_entries=True)
    if isinstance(link, NoLink):
        raise MixedTurnChainError(f"consecutive chain stations turn opposite ways (clockwise {clockwise} then {target_clockwise}): a chain's turn direction must be uniform.")
    return link


@dataclass(frozen=True)
class LoopConditions:
    """The narrow ladder's injected conditions after the radius-step bound.

    Attributes:
        admissibility: The verdicts.
    """

    admissibility: Admissibility

    def bridge_admissible(self, stock: Stock, origin: NarrowCandidate, candidate: NarrowCandidate, /) -> bool:
        """The bridge from *origin* stays in the centre domain and under the cap, in process.

        Raises:
            MixedTurnChainError: *origin* and *candidate* turn opposite ways.
        """
        link = chain_link(origin.loop, candidate.loop, origin.station.clockwise, candidate.station.clockwise)
        return all(_piece_admissible(self.admissibility, stock, piece, origin.station.clockwise) for piece in link.pieces)

    def loop_admissible_at_rest(self, stock: Stock, origin: NarrowCandidate, candidate: NarrowCandidate, /) -> bool:
        """rho's rule: the candidate's at-rest probe ring, phased on the advance direction, stays under the cap.

        The NAMED EXCEPTION to plan R1(a) (User ruling (2026-10-16), decision
        1 = A): the one generator decision not measured in process. The ring
        reads each position against the stock as it stands, without the loop's
        own swept prefix. It is kept because it closes the engagement-step gaps
        at the bridge-to-loop junction that an in-process loop check alone
        leaves open (docs/pipeline.md). It is AND-ed with the in-process bridge
        and loop, so it only ever refuses more.

        Raises:
            CoincidentStationsError: The two centres coincide, so there is no
                advance direction to phase the ring on (Ruling 1.1-R1-loud).
        """
        if origin.loop.centre == candidate.loop.centre:
            raise CoincidentStationsError(f"stations {origin.index} and {candidate.index} share the centre {origin.loop.centre!r}: no advance direction exists.")
        dx, dy = candidate.loop.centre.x - origin.loop.centre.x, candidate.loop.centre.y - origin.loop.centre.y
        length = math.hypot(dx, dy)
        return self.admissibility.loop_within_cap_at_rest(stock, candidate.loop, Direction2[WorldXY].build(dx / length, dy / length))

    def loop_admissible(self, stock: Stock, origin: NarrowCandidate, candidate: NarrowCandidate, /) -> bool:
        """The candidate's full loop stays under the cap, in process (*origin* unused)."""
        del origin
        return self.admissibility.loop_within_cap(stock, candidate.loop, candidate.station.clockwise)

    def loop_admissible_past_entry(self, stock: Stock, origin: NarrowCandidate, candidate: NarrowCandidate, /) -> bool:
        """`loop_admissible` without station 0, which the at-rest ring's entry probe already asked (*origin* unused).

        Station 0 of the in-process loop is the same exact query as the ring's
        first position (`engagement_at` at the candidate's entry, same stock), so
        AND-ed with the ring this decides exactly what `loop_admissible` does.
        Only `ordered`, which always carries the ring, uses it.
        """
        del origin
        return self.admissibility.loop_within_cap_past_entry(stock, candidate.loop, candidate.station.clockwise)

    @property
    def ordered(self) -> tuple[CandidateCondition[NarrowCandidate], ...]:
        """The conditions in evaluation order (after the structural radius-step bound).

        Cheapest refusal first. The conditions are read-only and every one
        reads the same unchanged stock, and the ladder takes their
        conjunction, so the order sets only the cost, never a decision
        (`test_every_ladder_condition_reads_one_unchanged_stock` pins the
        one-unchanged-stock premise). The at-rest ring (short-circuiting point
        probes) runs before the bridge (in-process stations). The in-process
        loop leaves its station 0 to the ring's entry probe
        (`loop_admissible_past_entry`): the conjunction is unchanged.
        """
        return (self.loop_admissible_at_rest, self.bridge_admissible, self.loop_admissible_past_entry)


def _piece_admissible(admissibility: Admissibility, stock: Stock, piece: LinkPiece, clockwise: bool) -> bool:
    """A link piece's admissibility, decided like every emitted motion (Ruling 1.1-I2).

    Straights: centre domain, then cap in process. Arcs lie on a gouge-free
    loop's circle (in the domain by construction) and are decided against the
    cap in process at the endpoints and parameters of the operation that would
    be emitted, which is what the gate's survey reads of the motion.

    The stock is NOT what the survey reads: every piece of a link is judged
    against *stock* as it stands before the link, whereas the survey replays
    each piece against the stock its predecessors already depleted. Earlier
    pieces only remove material, so the pre-link stock engages at least as
    much: the check is conservative -- it can refuse more, never admit more.
    """
    if isinstance(piece, LinePiece):
        return admissibility.straight_in_domain(piece.start, piece.end) and admissibility.straight_within_cap(stock, piece.start, piece.end)
    geometry = arc_operation(piece.centre, piece.start, piece.end, clockwise, CutZ.build(0.0), 0, OperationType.LINK).geometry
    start, end = geometry.point_at(0.0), geometry.point_at(1.0)
    emitted_start = Point2[WorldXY].build(float(start[0]), float(start[1]))
    emitted_end = Point2[WorldXY].build(float(end[0]), float(end[1]))
    return admissibility.arc_within_cap(stock, piece.centre, emitted_start, emitted_end, clockwise)


def _plunge_into_material(stock: Stock, entry: Point2[WorldXY], tool_radius: float) -> bool:
    """The gate's plunge predicate: the plunge disk removes material."""
    return stock.circle_sweep_removal(entry.x, entry.y, PLUNGE_GUIDE_RADIUS, tool_radius) is SweepRemoval.REMOVES


def skip_loop_removing_nothing(path: PathState, chain: OrientedChain, index: int) -> bool:
    """The NOTHING-loop skip (one named rule, Ruling 1.1-M7): a loop whose exact annulus removes nothing is not emitted.

    Args:
        path: The generation state.
        chain: The chain being walked.
        index: The station, in walk order.

    Returns:
        True iff the loop was skipped (and recorded as `OmittedLoop`).
    """
    loop = chain.loops[index]
    assert loop is not None  # the walk visits floor-clearing stations only
    verdict = path.stock.circle_sweep_removal(loop.centre.x, loop.centre.y, loop.radius.value, path.tool_radius)
    if verdict is SweepRemoval.NOTHING:
        path.omitted.append(OmittedLoop(path_index=path.path_index, chain=chain.origin, station=chain.guide_index(index), verdict=verdict))
        return True
    return False


def emit_loop(path: PathState, chain: OrientedChain, index: int) -> None:
    """Emit the station's loop unless it removes nothing."""
    if skip_loop_removing_nothing(path, chain, index):
        return
    loop = chain.loops[index]
    assert loop is not None
    path.emit_motion(loop_operation(loop, path.cut, path.path_index, chain.stations[index].clockwise))


def _piece_operation(piece: LinkPiece, clockwise: bool, cut: CutZ, path_index: int) -> ToolpathOperation:
    if isinstance(piece, ArcPiece):
        return arc_operation(piece.centre, piece.start, piece.end, clockwise, cut, path_index, OperationType.LINK)
    return line_operation(piece.start, piece.end, cut.value, cut.value, OperationType.LINK, path_index)


def emit_link(path: PathState, link: Link, clockwise: bool) -> None:
    """Emit a link's pieces in travel order, each labelled by its own exact removal, and record it."""
    first = len(path.operations)
    for piece in link.pieces:
        path.emit_motion(_piece_operation(piece, clockwise, path.cut, path.path_index))
    if len(path.operations) > first:
        path.links.append(LinkRecord(path_index=path.path_index, kind=_LINK_KINDS[type(link)], first_operation=first, operation_count=len(path.operations) - first))


_LINK_KINDS: dict[type, LinkKind] = {SameCircleArc: LinkKind.SAME_CIRCLE_ARC, BeltLink: LinkKind.BELT, StraightLink: LinkKind.STRAIGHT}


def enter(path: PathState, chain: OrientedChain, station: int, reason: ForcedEntryReason, from_point: Point2[WorldXY] | None) -> None:
    """The entry stage: retract and rapid if the tool is down, plunge onto the station's exact entry, cut its root loop.

    One stage, so a later entry (a helix) replaces it whole: the plunge and
    the root loop together. The walk that follows starts at this station with
    ``emit_start_loop=False``.

    Args:
        path: The generation state.
        chain: The chain entered.
        station: The station entered at, in walk order.
        reason: Why the entry is forced.
        from_point: Where the tool is at the cut plane, or None before the first cut.
    """
    loop = chain.loops[station]
    assert loop is not None
    entry = loop.entry
    if from_point is not None:
        path.operations.append(line_operation(from_point, from_point, path.cut.value, path.clearance.value, OperationType.RETRACT, path.path_index))
        path.operations.append(line_operation(from_point, entry, path.clearance.value, path.clearance.value, OperationType.LINK, path.path_index))
    into_material = _plunge_into_material(path.stock, entry, path.tool_radius)
    path.operations.append(line_operation(entry, entry, path.clearance.value, path.cut.value, OperationType.PLUNGE, path.path_index))
    path.stock.subtract_disk(entry.x, entry.y, path.tool_radius)
    path.records.append(ForcedEntry(path_index=path.path_index, entry=entry, reason=reason, into_material=into_material))
    emit_loop(path, chain, station)


@dataclass(frozen=True)
class WalkEnd:
    """Where a chain walk left the tool.

    Attributes:
        loop: The last loop the walk reached (emitted or omitted); the tool is at its entry.
        clockwise: That loop's turn direction.
        resume: The station to re-enter the chain at after a split, in walk order, or None.
    """

    loop: IncidentLoop
    clockwise: bool
    resume: int | None


def walk(path: PathState, chain: OrientedChain, start_station: int, emit_start_loop: bool, window: int, conditions: LoopConditions) -> WalkEnd:
    """Walk one oriented chain from *start_station* to its last floor-clearing station, or to a split.

    Args:
        path: The generation state.
        chain: The chain.
        start_station: The station the tool is at, in walk order.
        emit_start_loop: Whether to emit the start station's loop (False when
            an earlier stage already cut it).
        window: The advance window, in stations.
        conditions: The ladder's injected conditions.

    Returns:
        Where the walk ended.

    Raises:
        MixedTurnChainError: Two stations the walk joins turn opposite ways
            (the chain's turn direction must be uniform).
    """
    index = start_station
    if emit_start_loop:
        emit_loop(path, chain, index)
    while index < chain.last:
        advance = advance_ladder(path.stock, chain, index, min(index + window, chain.last), conditions.ordered)
        if not isinstance(advance, Advance):
            path.records.append(advance.gap)
            current = chain.loops[index]
            assert current is not None
            return WalkEnd(loop=current, clockwise=chain.stations[index].clockwise, resume=advance.resume)
        if advance.forced is not None:
            path.records.append(advance.forced)
        origin, target = chain.loops[index], chain.loops[advance.index]
        assert origin is not None and target is not None
        link = chain_link(origin, target, chain.stations[index].clockwise, chain.stations[advance.index].clockwise)
        emit_link(path, link, chain.stations[index].clockwise)
        index = advance.index
        emit_loop(path, chain, index)
    last = chain.loops[index]
    assert last is not None
    return WalkEnd(loop=last, clockwise=chain.stations[index].clockwise, resume=None)


@dataclass(frozen=True)
class NextChain:
    """The next-chain choice.

    Attributes:
        pending: Index of the chosen chain among the pending ones.
        chain: The chosen orientation.
        link: The admissible link into it, or None when it must be entered from above.
    """

    pending: int
    chain: OrientedChain
    link: Link | None


def _link_admissible(path: PathState, link: Link, chain: OrientedChain) -> bool:
    """*link* and the root loop it arrives at are admissible: a pure conjunction, cheapest refusal first.

    The arriving loop is cut by `walk_chain(..., emit_start_loop=True)` with no
    ladder in front of it, so this is its only cap decision. Pieces and loop
    are judged on the PRE-link stock, the stock the ladder judges a bridge and
    its target loop on: it holds at least the material of the post-link stock,
    and engagement only grows with material, so the verdict is conservative --
    it refuses more, never admits more (sequential folding of pieces and loop
    is Task 1.6 work).

    Args:
        path: The generation state.
        link: The link into *chain*.
        chain: The orientation the link arrives at.

    Returns:
        True iff every piece and the arriving root loop are admissible.
    """
    clockwise = chain.entry_station.clockwise
    lines = [piece for piece in link.pieces if isinstance(piece, LinePiece)]
    if not all(path.admissibility.straight_in_domain(piece.start, piece.end) for piece in lines):
        return False
    if not all(_piece_admissible(path.admissibility, path.stock, piece, clockwise) for piece in link.pieces):
        return False
    arriving = chain.loops[chain.first]
    assert arriving is not None  # `orientations` starts every chain on a floor-clearing station
    return path.admissibility.loop_within_cap(path.stock, arriving, clockwise)


def _link_length(link: Link) -> float:
    """The ordering length of a link: its straight part (arcs travel through cleared stock)."""
    if isinstance(link, BeltLink):
        return link.tangent_length
    return 0.0


def choose_next(path: PathState, current: WalkEnd, pending: list[tuple[OrientedChain, ...]]) -> NextChain:
    """The next chain: an admissibly linked one first, then the shortest link, then guide order.

    The choice minimises ``(not linked, link length, pending position)``, the
    first offered orientation winning ties. Length and position are link
    geometry; only ``linked`` needs the in-process checks. So the links are
    checked lazily in ``(length, position, offer order)`` order -- a stable
    sort -- and the first admissible one IS the minimum; a link that cannot
    win is never checked. With none admissible, the minimum of ``(length or 0,
    position)`` over every offer is entered from above, as before.

    Args:
        path: The generation state.
        current: Where the tool is: the loop it is on and that loop's turn direction.
        pending: The chains not yet walked, each in its offered orientations.

    Returns:
        The choice.
    """
    offers: list[tuple[int, OrientedChain, Link | NoLink]] = []
    for position, offered in enumerate(pending):
        for chain in offered:
            first = chain.loops[chain.first]
            assert first is not None
            offers.append((position, chain, link_between(current.loop, first, current.clockwise, chain.entry_station.clockwise)))
    assert offers  # pending is non-empty
    linkable = sorted(((position, chain, link) for position, chain, link in offers if not isinstance(link, NoLink)), key=lambda offer: (_link_length(offer[2]), offer[0]))
    for position, chain, link in linkable:
        if _link_admissible(path, link, chain):
            return NextChain(pending=position, chain=chain, link=link)
    position, chain, _link = min(offers, key=lambda offer: (0.0 if isinstance(offer[2], NoLink) else _link_length(offer[2]), offer[0]))
    return NextChain(pending=position, chain=chain, link=None)


def generate_narrow(
    polygon: Polygon,
    holes: list[Polygon],
    tool: ToolRadius,
    cap: Radian,
    clearance_z: ClearanceZ,
    cut: CutZ | None = None,
    climb: bool = True,
    max_passes: int = DEFAULT_MAX_PASSES,
) -> PipelineResult:
    """Machine a pocket in the narrow regime with one plunge into material where the guide allows.

    Args:
        polygon: Outer pocket boundary in world XY.
        holes: Islands strictly inside *polygon*.
        tool: Tool radius.
        cap: Engagement cap in ``(0, pi]``.
        clearance_z: Rapid plane, strictly above the cut plane.
        cut: Cut plane; ``z = 0`` by default.
        climb: Climb milling (clockwise loops).
        max_passes: Most skeleton chains the guide may emit.

    Returns:
        The path, its ledger, and what the rules left out.

    Raises:
        InvalidClearanceHeightError: *clearance_z* is not above the cut plane.
        EmptyNarrowPathError: No guide station clears the degeneracy floor.
        GuideTruncatedError: *max_passes* truncated the guide.
        EmptyPipelineGuideError: The guide emitted no station.
        InvalidEngagementCapError: *cap* is outside ``(0, pi]``.
    """
    cut_plane = CutZ.build(0.0) if cut is None else cut
    if not clearance_z.value > cut_plane.value:
        raise InvalidClearanceHeightError(f"clearance plane {clearance_z.value!r} must be above the cut plane {cut_plane.value!r}.")
    diameter = 2.0 * tool.value
    guide = Guide.build(
        polygon,
        holes,
        tool,
        pipeline_spacing(diameter, cap),
        Clearance.build(RADIAL_CLEARANCE_FRACTION * diameter),
        climb=climb,
        max_passes=max_passes,
        clearance_z=clearance_z,
    )
    pending = [offered for offered in (orientations(chain, tool) for chain in guide.chains) if offered]
    if not pending:
        raise EmptyNarrowPathError(f"the guide produced {len(guide.chains)} chain(s), but no station clears the degeneracy floor at tool radius {tool.value!r}.")
    admissibility = Admissibility.build(tool, cap_ratio(cap), centre_domain(polygon, holes, tool))
    path = PathState(stock=Stock(polygon, holes=holes), admissibility=admissibility, cut=cut_plane, clearance=clearance_z)
    conditions = LoopConditions(admissibility)
    window = advance_window(GUIDE_STEP_TOOL_DIAMETERS)
    first = pending.pop(0)[0]
    path.records.extend(first.declined)
    enter(path, first, first.first, ForcedEntryReason.FIRST_CUT, None)
    end = walk_chain(path, first, first.first, window, conditions, emit_start_loop=False)
    while pending:
        path.path_index += 1
        choice = choose_next(path, end, pending)
        pending.pop(choice.pending)
        path.records.extend(choice.chain.declined)
        if choice.link is None:
            enter(path, choice.chain, choice.chain.first, ForcedEntryReason.NO_ADMISSIBLE_LINK, end.loop.entry)
            end = walk_chain(path, choice.chain, choice.chain.first, window, conditions, emit_start_loop=False)
            continue
        emit_link(path, choice.link, choice.chain.entry_station.clockwise)
        end = walk_chain(path, choice.chain, choice.chain.first, window, conditions, emit_start_loop=True)
    exit_point = end.loop.entry
    path.operations.append(line_operation(exit_point, exit_point, cut_plane.value, clearance_z.value, OperationType.RETRACT, path.path_index))
    return PipelineResult(
        operations=path.operations,
        polyline=tessellate(path.operations),
        ledger=ForcedDecisionLedger.build(path.records),
        omitted=tuple(path.omitted),
        links=tuple(path.links),
    )


def walk_chain(path: PathState, chain: OrientedChain, start: int, window: int, conditions: LoopConditions, *, emit_start_loop: bool) -> WalkEnd:
    """Walk *chain* from *start*, re-entering after every split, until its last station.

    Args:
        path: The generation state.
        chain: The chain.
        start: The station the tool is at, in walk order.
        window: The advance window, in stations.
        conditions: The ladder's injected conditions.
        emit_start_loop: False when the entry stage already cut the start
            station's loop; True when a link arrived at it.

    Returns:
        Where the walk ended.
    """
    end = walk(path, chain, start, emit_start_loop, window, conditions)
    while end.resume is not None:
        enter(path, chain, end.resume, ForcedEntryReason.RADIUS_STEP_SPLIT, end.loop.entry)
        end = walk(path, chain, end.resume, False, window, conditions)
    return end
