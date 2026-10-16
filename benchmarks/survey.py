"""One replay of a toolpath that records everything the quality groups reduce.

Four groups of metrics ask four different questions of the same path, and every
one of them needs the state of the material at the moment a motion cuts it. That
state exists once, during a depletion replay, and reconstructing it four times
would be four times the cost and four chances to drift apart. So the replay
happens HERE, once, and produces a `PathSurvey` of per-motion evidence that
`benchmarks.quality` reduces without touching a stock again.

The replay itself is `benchmarks.depletion.replay_events` -- the primitive that
owns the cut-plane model -- not a second copy of it. What this module adds is the
observation: engagement sampled along each motion, the exact effect of each
motion on the stock, and the geometric quantities (curvature, tangents, swept
area) that the speed and longevity groups read.

ENGAGEMENT IS MEASURED IN PROCESS. A station inside a motion meets the stock the
motion started from MINUS the motion's own swept prefix up to that station
(`_stock_2.engagement_line_prefixes` / `engagement_arc_prefixes`): a straight
slot into solid stock reads a half turn, not the full rim the static stock would
report. Stations are taken in travel order, so a clockwise loop is walked
clockwise. The generators' own probe rings still decide against the static
stock; that is theirs to change, not the gate's.

WHAT IS EXACT AND WHAT IS SAMPLED, per field, is stated on `MotionQuality`. The
short version: effects on the stock are exact set operations; everything read
along a motion is sampled at `samples_per_motion` positions, each position's
verdict exact, the reduction over the motion not. A sampled field can demonstrate
a defect and can never certify its absence.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import angle_vectors
from compas.tolerance import TOL

from benchmarks.depletion import CutMotion
from benchmarks.depletion import PlungeEvent
from benchmarks.depletion import replay_events
from benchmarks.errors import InconsistentRemovalVerdictError
from benchmarks.errors import InvalidMotionSampleCountError
from benchmarks.errors import UnmeasurableOperationLengthError
from benchmarks.errors import UnreplayableOperationError
from benchmarks.errors import UnsampleableMotionError
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.held_path_snapshot import snapshot_toolpath
from benchmarks.spec import PocketSpec
from benchmarks.units import OperationIndex
from compas_cgal import _coverage_2
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.engagement import _infer_cut_height
from compas_cgal.engagement import _operation_removal
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import classify_operation_replay
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.stock import _polygon_to_ccw_vertices
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Cutter positions probed per cut motion. Chosen against the instrument the
# generator regulates itself with -- a 32-position ring, `LOOP_PROBE_COUNT` in
# `compas_cgal.engagement_toolpath` -- and not against a cost budget: 45 samples
# sit 8 degrees apart against the generator's 11.25, and gcd(45, 32) = 1, so the
# seam is the only position the two rings share. A gate sampling at 32, 16, or 8
# would be asking the generator to confirm its own verdict. The two instruments
# also measure different things: the gate reads each station against the stock
# minus the motion's own swept prefix, the generator's ring reads the static
# stock it is about to cut.
QUALITY_SAMPLES_PER_MOTION = 45

# The engagement query is asked WITHOUT gap-closure pessimism, matching
# `benchmarks.exceedance.NO_GAP_CLOSURE`: this is a question about the material
# actually engaged at a position, not about what a certifier must assume could
# happen before the next station.
NO_GAP_CLOSURE = 0.0


class MotionKind(str, enum.Enum):
    """The geometric primitive a cut motion follows.

    Attributes:
        LOOP: A closed machining circle.
        ARC: An open circular arc.
        LINE: A straight cut, bridge, or cut-height link.
    """

    LOOP = "loop"
    ARC = "arc"
    LINE = "line"


@dataclass(frozen=True)
class EngagementSample:
    """One probed cutter position on a motion.

    Attributes:
        distance: Arc length from the motion's start to this position.
        position: Cutter-centre position in the world XY frame.
        engagement_deg: The largest engaged run of the cutter rim, in degrees.
            REPORTING ONLY: this value never decides cap exceedance.
        cap_exceeded: The exact cap predicate's verdict at this position.
        inside_centre_domain: Whether the cutter centre lies where a tool of this
            radius may legally be. False is a demonstrated gouge.

    Each retained ``cap_exceeded`` verdict is exact. Absence of exceedance across
    this finite station set remains sampled-negative evidence, not a certificate
    for the continuous motion.
    """

    distance: float
    position: Point2[WorldXY]
    engagement_deg: float
    cap_exceeded: bool
    inside_centre_domain: bool


@dataclass(frozen=True)
class MotionQuality:
    """Everything one cut motion did, measured once.

    Attributes:
        index: Position of the operation in its toolpath.
        operation: The operation's declared type.
        kind: The geometric primitive it follows.
        length: The primitive's exact length.
        loop_radius: Guide radius of a closed machining circle, or None when the
            motion is not one.
        curvature: Reciprocal of the guide radius; zero for a straight motion.
        swept_area: Area the tool sweeps, in closed form for the primitive.
            EXACT arithmetic, not a sample and not a boolean measurement.
        start: Cutter-centre position at the motion's start.
        end: Cutter-centre position at the motion's end.
        start_tangent: Unit travel direction at the start.
        end_tangent: Unit travel direction at the end.
        samples: The probed positions, in travel order.
        cap_exceeded: Whether the exact cap predicate fired at any sample.
            SAMPLED: a lower bound on true exceedance.
        slot_exceeded: Whether a STRAIGHT motion exceeded the cap at any
            sample: the straight-motion subset of *cap_exceeded*, judged by the
            same cap as a loop (gate rule c). SAMPLED.
        removes_material: Whether the motion may remove material: False only
            when its exact removal verdict is NOTHING, i.e. even a certified
            region CONTAINING the true swept area misses the stock's interior.
            A WITHIN_BAND motion counts as removing, because the band it may
            cut cannot be ruled out. EXACT.
        removal: The exact three-valued removal verdict (`SweepRemoval`) the
            replay decided for this motion, or None for a record assembled
            without one (synthetic fixtures). When present it determines
            *removes_material*.
        exit_sample: The station where the motion hands over to the next one,
            read for the junction step (gate rule f) and for nothing else. For a
            line or an arc it is the last of *samples*. For a closed loop it is
            an extra station at the seam with the FULL turn as prefix (Ruling
            0.3-f-loop): the seam is sampled once, at the start, so the last of
            *samples* sits a station spacing before the loop ends. It does not
            enter *cap_exceeded* or any per-motion statistic. None only on a
            record assembled without a survey (synthetic fixtures). The
            junction reads it only for a loop, and refuses a loop without it
            (`MissingLoopSeamStationError`); a line or arc junction reads the
            last of *samples* by definition.
    """

    index: int
    operation: OperationType
    kind: MotionKind
    length: float
    loop_radius: Optional[float]
    curvature: float
    swept_area: float
    start: Tuple[float, float]
    end: Tuple[float, float]
    start_tangent: Tuple[float, float]
    end_tangent: Tuple[float, float]
    samples: Tuple[EngagementSample, ...]
    cap_exceeded: bool
    slot_exceeded: bool
    removes_material: bool
    removal: Optional[SweepRemoval] = None
    exit_sample: Optional[EngagementSample] = None

    def __post_init__(self) -> None:
        """Refuse a record whose two removal fields disagree.

        Raises:
            InconsistentRemovalVerdictError: *removes_material* is not
                ``removal is not SweepRemoval.NOTHING``.
        """
        if self.removal is not None and self.removes_material != (self.removal is not SweepRemoval.NOTHING):
            raise InconsistentRemovalVerdictError(f"motion {self.index}: removes_material={self.removes_material!r} contradicts the verdict {self.removal.name}.")

    @property
    def peak_engagement_deg(self) -> float:
        """Largest sampled engaged run, in degrees; ``0.0`` for no samples."""
        return max((sample.engagement_deg for sample in self.samples), default=0.0)

    @property
    def gouges(self) -> bool:
        """Whether any sampled cutter centre lies outside the legal centre domain."""
        return any(not sample.inside_centre_domain for sample in self.samples)

    @property
    def is_engaged(self) -> bool:
        """Whether the motion touched material at any sampled position."""
        return any(sample.engagement_deg > 0.0 for sample in self.samples)


@dataclass(frozen=True)
class RapidMotion:
    """One motion that removes nothing.

    Attributes:
        index: Position of the operation in its toolpath.
        operation: The operation's declared type.
        length: The primitive's exact length.
        horizontal_at_cut_plane: Whether the move travels in XY while sitting on
            the cutting plane. In a 2.5D cut-plane model this is the only shape a
            rapid can take that would drive the cutter through material, so it is
            what `rapid_safety` tests. See `benchmarks.quality` on why that
            reduction is narrow but not vacuous.
    """

    index: int
    operation: OperationType
    kind: MotionKind
    length: float
    horizontal_at_cut_plane: bool


@dataclass(frozen=True)
class PathSurvey:
    """Every motion's findings, plus the state the path left the stock in.

    Attributes:
        spec: The instance the path was generated for.
        source_snapshot: Immutable structure of the exact operation stream that
            produced this survey.
        motions: One record per cut-plane motion, in toolpath order.
        rapids: One record per non-removing motion, in toolpath order.
        plunges: Downward bores in the path.
        retracts: Upward clearance moves in the path.
        plunge_indices: Source operation indices of the downward bores.
        material_plunge_indices: The plunges whose disk met the stock's
            interior when it was bored, a subset of *plunge_indices*. EXACT
            (`benchmarks.depletion.PlungeEvent`).
        retract_indices: Source operation indices of the upward clearance moves.
        final_stock: The stock after the whole path has been replayed.
        total_length: Analytic length of every operation.
        cut_length: Analytic length of the cut-plane motions.
        air_length: Analytic length of the motions that remove nothing.
        plunge_swept_area: Area bored by the plunges, in closed form.
    """

    spec: PocketSpec
    source_snapshot: tuple[HeldOperationSnapshot, ...]
    motions: Tuple[MotionQuality, ...]
    rapids: Tuple[RapidMotion, ...]
    plunges: int
    retracts: int
    plunge_indices: tuple[OperationIndex, ...]
    retract_indices: tuple[OperationIndex, ...]
    material_plunge_indices: tuple[OperationIndex, ...]
    final_stock: Stock
    total_length: float
    cut_length: float
    air_length: float
    plunge_swept_area: float

    @property
    def swept_area(self) -> float:
        """Total area the tool sweeps, cut motions and plunges together."""
        return sum(motion.swept_area for motion in self.motions) + self.plunge_swept_area

    @property
    def within_band_motions(self) -> int:
        """Motions whose exact removal verdict is WITHIN_BAND.

        Their stock change is undecidable on the replay's certified model: the
        model removes nothing, the true sweep may cut material inside the
        rectangle's side band. Counted as removing in the redundancy reduction; reported
        here so the undecided share is visible rather than folded in.
        """
        return sum(1 for motion in self.motions if motion.removal is SweepRemoval.WITHIN_BAND)


def survey_path(spec: PocketSpec, result: ToolpathResult, *, samples_per_motion: int = QUALITY_SAMPLES_PER_MOTION) -> PathSurvey:
    """Replay *result* once and record what every motion did.

    Args:
        spec: The pocket, tool, and engagement cap.
        result: The generated toolpath.
        samples_per_motion: Cutter positions probed along each cut motion.

    Returns:
        The per-motion findings and the depleted stock.

    Raises:
        InvalidMotionSampleCountError: *samples_per_motion* is below one.
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnsampleableMotionError: A cut motion carries no cutter-centre path.
        UnmeasurableOperationLengthError: An operation's length is undefined.
        InvalidPolygonError: The pocket boundary or a hole is degenerate.
        ReachableDomainConstructionError: The kernel could not build the legal
            centre domain for this pocket and tool.
    """
    if samples_per_motion < 1:
        raise InvalidMotionSampleCountError(f"samples_per_motion must be at least 1, got {samples_per_motion!r}.")

    _validate_cut_plane_curves(result.operations)

    # The cap crosses into exact-land the way the kernel's own cap does: the
    # transcendental intent (an angle) is converted once, here, to the exact
    # rational surrogate 4*sin^2(theta/2), and the predicate then decides
    # against that surrogate exactly. Slotting is judged by the same cap.
    cap_ratio = _cap_chord_ratio(spec.tea_cap_rad)
    holes = [_polygon_to_ccw_vertices(hole) for hole in spec.holes]
    centre_domain = _coverage_2.CutterCentreDomain2.build(
        _polygon_to_ccw_vertices(spec.polygon),
        holes,
        spec.tool_radius,
    )

    stock = Stock(spec.polygon, list(spec.holes))
    motions: List[MotionQuality] = []
    material_plunges: list[OperationIndex] = []
    for event in replay_events(spec, result, stock):
        if isinstance(event, PlungeEvent):
            if event.into_material:
                material_plunges.append(OperationIndex(event.index))
            continue
        motions.append(_measure_motion(event, spec.tool_radius, cap_ratio, centre_domain, samples_per_motion))

    rapids, plunge_indices, retract_indices, plunge_area = _classify_non_cutting(result, spec.tool_radius)
    cut_length = sum(motion.length for motion in motions)
    air_length = sum(rapid.length for rapid in rapids)
    return PathSurvey(
        spec=spec,
        source_snapshot=snapshot_toolpath(result),
        motions=tuple(motions),
        rapids=tuple(rapids),
        plunges=len(plunge_indices),
        retracts=len(retract_indices),
        plunge_indices=plunge_indices,
        retract_indices=retract_indices,
        material_plunge_indices=tuple(material_plunges),
        final_stock=stock,
        total_length=sum(_primitive_length(index, operation) for index, operation in enumerate(result.operations)),
        cut_length=cut_length,
        air_length=air_length,
        plunge_swept_area=plunge_area,
    )


def _measure_motion(
    motion: CutMotion,
    tool_radius: float,
    cap_ratio: float,
    centre_domain: "_coverage_2.CutterCentreDomain2",
    samples_per_motion: int,
) -> MotionQuality:
    """Measure one cut motion against every criterion the groups reduce.

    Args:
        motion: The motion and the stock it is about to cut.
        tool_radius: Tool radius.
        cap_ratio: Exact squared-chord surrogate of the engagement cap.
        centre_domain: Where a cutter of this radius may legally be centred.
        samples_per_motion: Cutter positions to probe.

    Returns:
        The motion's findings.

    Raises:
        UnsampleableMotionError: The motion carries no cutter-centre path.
        UnmeasurableOperationLengthError: The motion's length is undefined.
    """
    geometry = motion.operation.geometry
    straight = isinstance(geometry, Line)
    positions = _motion_travel_positions(motion, samples_per_motion)
    verdicts = _motion_prefix_samples(motion, samples_per_motion, tool_radius, cap_ratio)
    samples: List[EngagementSample] = []
    for (distance, x, y), (_total_tea, max_run_tea, exceeded) in zip(positions, verdicts, strict=True):
        samples.append(
            EngagementSample(
                distance=distance,
                position=Point2[WorldXY].build(x, y),
                engagement_deg=math.degrees(max_run_tea),
                cap_exceeded=exceeded,
                inside_centre_domain=centre_domain.contains(x, y),
            )
        )

    removal = _operation_removal(motion.stock, motion.operation, tool_radius)
    cap_exceeded = any(sample.cap_exceeded for sample in samples)
    exit_sample = _loop_seam_exit(motion, tool_radius, cap_ratio, centre_domain) if isinstance(geometry, Circle) else samples[-1]
    start, end = _endpoints(motion.operation)
    start_tangent, end_tangent = _tangents(motion.operation)
    return MotionQuality(
        index=motion.index,
        operation=motion.operation.operation,
        kind=_motion_kind(geometry),
        length=_primitive_length(motion.index, motion.operation),
        loop_radius=float(geometry.radius) if isinstance(geometry, Circle) else None,
        curvature=0.0 if straight else 1.0 / float(geometry.radius),
        swept_area=_swept_area(geometry, tool_radius),
        start=start,
        end=end,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
        samples=tuple(samples),
        cap_exceeded=cap_exceeded,
        slot_exceeded=straight and cap_exceeded,
        removes_material=removal is not SweepRemoval.NOTHING,
        removal=removal,
        exit_sample=exit_sample,
    )


# The prefix parameter of a closed loop's own end: the full turn.
FULL_TURN_PREFIX = 1.0


def _loop_seam_exit(motion: CutMotion, tool_radius: float, cap_ratio: float, centre_domain: "_coverage_2.CutterCentreDomain2") -> EngagementSample:
    """The junction station of a closed loop: its seam, after the whole turn.

    One native prefix query at the full-turn parameter, against the stock the
    loop started from (Ruling 0.3-f-loop). The position is the loop's entry
    (`_circle_entry`), where the loop both begins and ends.

    Args:
        motion: The loop and the stock it is about to cut.
        tool_radius: Tool radius.
        cap_ratio: Exact squared-chord surrogate of the cap.
        centre_domain: Where a cutter of this radius may legally be centred.

    Returns:
        The seam station with the full-turn prefix.
    """
    geometry = motion.operation.geometry
    assert isinstance(geometry, Circle)  # guaranteed by the caller
    cx, cy = _xy(geometry.frame.point)
    sx, sy = _circle_entry(motion.operation)
    ((_total, max_run, exceeded),) = _stock_2.engagement_circle_prefixes(
        motion.stock.raw, cx, cy, float(geometry.radius), sx, sy, motion.operation.clockwise, [FULL_TURN_PREFIX], tool_radius, cap_ratio, NO_GAP_CLOSURE
    )
    return EngagementSample(
        distance=float(geometry.circumference),
        position=Point2[WorldXY].build(sx, sy),
        engagement_deg=math.degrees(max_run),
        cap_exceeded=exceeded,
        inside_centre_domain=centre_domain.contains(sx, sy),
    )


def _validate_cut_plane_curves(operations: Sequence[ToolpathOperation]) -> None:
    """Validate every engaged circular primitive before replay can mutate stock.

    Lines independently anchor the cut height because their endpoints encode
    the plunge and clearance structure. A line-free circular stream instead
    establishes its common plane from its first engaged curve.

    Args:
        operations: The complete toolpath operation stream.

    Raises:
        UnreplayableOperationError: An engaged arc or circle is tilted or does
            not lie on the independently established cutting plane.
    """
    curves = [(index, operation) for index, operation in enumerate(operations) if isinstance(operation.geometry, (Arc, Circle))]
    if not curves:
        return

    cut_z = _infer_cut_height(list(operations))
    for index, operation in curves:
        try:
            category = classify_operation_replay(operation, Millimetre(cut_z))
        except (CutPlaneRampError, OffPlaneReplayCurveError) as error:
            raise UnreplayableOperationError(f"Operation {index} ({operation.operation.value}) {error}.") from error
        if category == "motion":
            _require_world_xy_curve(index, operation.geometry)


def _require_world_xy_curve(index: int, geometry: object) -> None:
    """Refuse circular geometry that the world-XY depletion cannot represent.

    Args:
        index: Position of the operation, for the error message.
        geometry: The motion primitive.

    Raises:
        UnreplayableOperationError: An arc or circle is tilted away from world XY.
    """
    if not isinstance(geometry, (Arc, Circle)):
        return

    axis_angle = angle_vectors(geometry.frame.zaxis, [0.0, 0.0, 1.0])
    axis_is_world_z = TOL.is_angle_zero(axis_angle) or TOL.is_angles_close(axis_angle, math.pi)
    if not axis_is_world_z:
        raise UnreplayableOperationError(
            f"Operation {index} carries {type(geometry).__name__} geometry outside the inferred world-XY cut plane; projecting it to XY would misrepresent material removal."
        )


def _sample_parameters(geometry: object, count: int) -> List[float]:
    """Prefix parameters of the probes: a closed circle wraps, so its seam is sampled once.

    Args:
        geometry: The motion primitive.
        count: Positions to place.

    Returns:
        Non-decreasing parameters in [0, 1].
    """
    stations = count if isinstance(geometry, Circle) else count + 1
    return [i / count for i in range(stations)]


def _stationary_point(geometry: object) -> Optional[Tuple[float, float]]:
    """The point of a zero-length cut-plane line, or None for any other motion.

    Zero length is decided by EXACT equality of the injected XY endpoints: a
    double enters exact-land as itself, so equal doubles are equal exact points
    and unequal ones are not; there is no tolerance (Ruling 0.2-zerolen).

    Args:
        geometry: The motion's primitive.

    Returns:
        ``(x, y)`` when *geometry* is a `Line` whose XY endpoints coincide.
    """
    if not isinstance(geometry, Line):
        return None
    start = (float(geometry.start[0]), float(geometry.start[1]))
    end = (float(geometry.end[0]), float(geometry.end[1]))
    return start if start == end else None


def _motion_travel_positions(motion: CutMotion, count: int) -> Sequence[Tuple[float, float, float]]:
    """Probe positions along one motion, in travel order, as ``(arc distance, x, y)``.

    A closed circle wraps, so its seam is sampled once; an open arc and a segment
    keep both endpoints, where a motion's engagement is typically extreme. This
    is `benchmarks.exceedance`'s station convention, carrying the arc distance as
    well so the engagement gradient and the time-at-engagement histogram can be
    built from the same probes. Circles and arcs are walked in the operation's
    turn direction (``ToolpathOperation.clockwise``, the direction the replay
    sweeps), not in their frame's parametrisation. For an arc that matters
    beyond order: the replay sweeps from ``point_at(0)`` to ``point_at(1)`` in
    the flag's direction, which for a clockwise-flagged arc (the snapshot
    convention keeps ``start_angle <= end_angle``) is the other way round the
    circle from ``point_at``, so arc positions come from `_arc_travel`.

    The positions are REPORTING doubles (and the centre-domain membership
    queries); the engagement at each is decided by `_motion_prefix_samples` on
    the same parameters. A zero-length line is ONE position: the tool sits
    there and sweeps nothing (Ruling 0.2-zerolen).

    Args:
        motion: The cut motion to sample.
        count: Positions to place.

    Returns:
        The probes, in travel order.

    Raises:
        UnsampleableMotionError: The geometry is neither a line, an arc, nor a
            circle.
    """
    geometry = motion.operation.geometry
    stationary = _stationary_point(geometry)
    if stationary is not None:
        return [(0.0, stationary[0], stationary[1])]
    parameters = _sample_parameters(geometry, count)
    if isinstance(geometry, Circle):
        span = float(geometry.circumference)
        # The frame's parametrisation runs counter-clockwise in world XY iff its
        # normal is +Z; `_require_world_xy_curve` has already pinned it to +-Z.
        if motion.operation.entry_xy is not None:
            # An authored entry: station t is the entry turned by the travelled
            # angle, so station 0 IS the entry (placement for reporting only).
            cx, cy = _xy(geometry.frame.point)
            sx, sy = _circle_entry(motion.operation)
            turn = -math.tau if motion.operation.clockwise else math.tau
            rotations = [(math.cos(turn * t), math.sin(turn * t)) for t in parameters]
            return [(span * t, cx + cos * (sx - cx) - sin * (sy - cy), cy + sin * (sx - cx) + cos * (sy - cy)) for t, (cos, sin) in zip(parameters, rotations, strict=True)]
        frame_runs_ccw = float(geometry.frame.zaxis[2]) > 0.0
        travel_reversed = frame_runs_ccw == motion.operation.clockwise
        return [(span * t, *_xy(geometry.point_at((1.0 - t) % 1.0 if travel_reversed else t))) for t in parameters]
    if isinstance(geometry, Arc):
        cx, cy, sx, sy, turn = _arc_travel(geometry, motion.operation.clockwise)
        radius = math.hypot(sx - cx, sy - cy)
        positions = []
        for t in parameters:
            angle = turn * t
            cosine, sine = math.cos(angle), math.sin(angle)
            x = cx + cosine * (sx - cx) - sine * (sy - cy)
            y = cy + sine * (sx - cx) + cosine * (sy - cy)
            positions.append((radius * abs(turn) * t, x, y))
        return positions
    if isinstance(geometry, Line):
        x0, y0 = float(geometry.start[0]), float(geometry.start[1])
        x1, y1 = float(geometry.end[0]), float(geometry.end[1])
        span = float(geometry.length)
        return [(span * t, x0 + (x1 - x0) * t, y0 + (y1 - y0) * t) for t in parameters]
    raise UnsampleableMotionError(f"Operation {motion.index} ({motion.operation.operation.value}) has geometry {type(geometry).__name__!r}, which carries no cutter-centre path.")


def _arc_travel(geometry: Arc, clockwise: bool) -> Tuple[float, float, float, float, float]:
    """The travel of an arc as the replay sweeps it: centre, start, signed turn.

    The replay (`compas_cgal.engagement._subtract_operation`) and the native
    prefix query sweep from ``point_at(0)`` to ``point_at(1)`` in the direction
    *clockwise* names. The turn mirrors the native query's intent sweep
    (``engagement_arc_prefixes``): the atan2 difference of the endpoints in the
    travel direction, lifted into ``(0, 2 pi]``. REPORTING and placement intent
    only; the engagement at each station is decided natively.

    Args:
        geometry: The arc.
        clockwise: The operation's turn direction.

    Returns:
        ``(cx, cy, sx, sy, turn)``: the centre, the start, and the signed travel
        angle in radians (negative for clockwise travel).
    """
    cx, cy = _xy(geometry.frame.point)
    sx, sy = _xy(geometry.point_at(0.0))
    ex, ey = _xy(geometry.point_at(1.0))
    start_angle = math.atan2(sy - cy, sx - cx)
    end_angle = math.atan2(ey - cy, ex - cx)
    sweep = start_angle - end_angle if clockwise else end_angle - start_angle
    if sweep <= 0.0:
        sweep += math.tau
    return cx, cy, sx, sy, -sweep if clockwise else sweep


def _motion_prefix_samples(motion: CutMotion, count: int, tool_radius: float, cap_ratio: float) -> Sequence[Tuple[float, float, bool]]:
    """In-process engagement at each probe: stock minus the motion's own swept prefix.

    One native call per motion. The motion's start, end, centre and turn
    direction are exactly those `compas_cgal.engagement._subtract_operation`
    sweeps. The prefix is the TRUE swept region: for a loop, the replay removes
    it exactly; for a non-Pythagorean line, the replay's model keeps a side band
    up to 1e-4 r wide that the true prefix covers, and later motions may read
    that band as engagement (a documented model limit, not physics). The stock
    the motion starts from is read, never mutated.

    Args:
        motion: The motion and the stock it is about to cut.
        count: Positions to place, as in `_motion_travel_positions`.
        tool_radius: Tool radius.
        cap_ratio: Exact squared-chord surrogate of the threshold to decide.

    Returns:
        ``(total_tea, max_run_tea, cap_exceeded)`` per probe, in travel order.

    Raises:
        UnsampleableMotionError: The geometry is neither a line, an arc, nor a
            circle.
    """
    geometry = motion.operation.geometry
    raw = motion.stock.raw
    stationary = _stationary_point(geometry)
    if stationary is not None:
        # One station with an empty prefix: the stock the motion starts from,
        # read at the point the tool sits on. The native prefix query refuses a
        # zero-length segment, and should: there is no prefix to subtract.
        return [_stock_2.engagement_at(raw, stationary[0], stationary[1], tool_radius, cap_ratio, NO_GAP_CLOSURE)]
    parameters = _sample_parameters(geometry, count)
    if isinstance(geometry, Line):
        x0, y0 = float(geometry.start[0]), float(geometry.start[1])
        x1, y1 = float(geometry.end[0]), float(geometry.end[1])
        return _stock_2.engagement_line_prefixes(raw, x0, y0, x1, y1, parameters, tool_radius, cap_ratio, NO_GAP_CLOSURE)
    if isinstance(geometry, Circle):
        # Station 0 is the entry itself (`_circle_entry`) -- where the tool is,
        # and where the previous motion ended; stations p > 0 lie on the
        # AUTHORED radius, the circle the replay's annulus is built on.
        cx, cy = _xy(geometry.frame.point)
        sx, sy = _circle_entry(motion.operation)
        return _stock_2.engagement_circle_prefixes(raw, cx, cy, float(geometry.radius), sx, sy, motion.operation.clockwise, parameters, tool_radius, cap_ratio, NO_GAP_CLOSURE)
    if isinstance(geometry, Arc):
        cx, cy = _xy(geometry.frame.point)
        sx, sy = _xy(geometry.point_at(0.0))
        ex, ey = _xy(geometry.point_at(1.0))
        return _stock_2.engagement_arc_prefixes(raw, cx, cy, sx, sy, ex, ey, motion.operation.clockwise, parameters, tool_radius, cap_ratio, NO_GAP_CLOSURE)
    raise UnsampleableMotionError(f"Operation {motion.index} ({motion.operation.operation.value}) has geometry {type(geometry).__name__!r}, which carries no cutter-centre path.")


def _swept_area(geometry: object, tool_radius: float) -> float:
    """Closed-form area swept by a tool of *tool_radius* along *geometry*.

    A segment of length L sweeps a capsule, ``2 r L + pi r^2``. A guide circle of
    radius rho sweeps an annulus of area ``4 pi rho r`` when ``rho > r``, and a
    filled disk of area ``pi (rho + r)^2`` when it does not -- which is the same
    boundary `benchmarks.quality.degenerate_loops` tests, arrived at from the
    area rather than from the hole. An arc of swept angle alpha sweeps the
    corresponding annular sector plus one full end disk.

    Args:
        geometry: The motion's primitive.
        tool_radius: Tool radius.

    Returns:
        The swept area.

    Raises:
        UnmeasurableOperationLengthError: The primitive is not a line, arc, or
            circle.
    """
    if isinstance(geometry, Circle):
        rho = float(geometry.radius)
        if rho > tool_radius:
            return 4.0 * math.pi * rho * tool_radius
        return math.pi * (rho + tool_radius) ** 2
    if isinstance(geometry, Arc):
        rho = float(geometry.radius)
        alpha = float(geometry.angle)
        inner = max(0.0, rho - tool_radius)
        return 0.5 * alpha * ((rho + tool_radius) ** 2 - inner**2) + math.pi * tool_radius**2
    if isinstance(geometry, Line):
        return 2.0 * tool_radius * float(geometry.length) + math.pi * tool_radius**2
    raise UnmeasurableOperationLengthError(f"Geometry {type(geometry).__name__!r} has no closed-form swept area.")


def _motion_kind(geometry: object) -> MotionKind:
    """Classify a primitive.

    Args:
        geometry: The motion's primitive.

    Returns:
        Its kind.

    Raises:
        UnmeasurableOperationLengthError: The primitive is not a line, arc, or
            circle.
    """
    if isinstance(geometry, Circle):
        return MotionKind.LOOP
    if isinstance(geometry, Arc):
        return MotionKind.ARC
    if isinstance(geometry, Line):
        return MotionKind.LINE
    raise UnmeasurableOperationLengthError(f"Geometry {type(geometry).__name__!r} is not a toolpath primitive.")


def _circle_entry(operation: ToolpathOperation) -> Tuple[float, float]:
    """Where a closed circle starts and ends: its exact entry when authored, else ``point_at(0)``.

    An `compas_cgal.incident_loop.IncidentLoop` carries its exact entry in
    ``entry_xy``; its ``point_at(0)`` is only a double rounding of it (Ruling
    1.1-C14-A). Every other circle is read at ``point_at(0)``, unchanged.

    Args:
        operation: A circle operation.

    Returns:
        The entry.
    """
    if operation.entry_xy is not None:
        return float(operation.entry_xy[0]), float(operation.entry_xy[1])
    geometry = operation.geometry
    assert isinstance(geometry, Circle)  # guaranteed by the caller
    return _xy(geometry.point_at(0.0))


def _endpoints(operation: ToolpathOperation) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Cutter-centre positions where a motion starts and ends.

    A closed circle starts and ends at the same point -- its entry
    (`_circle_entry`) -- which is what makes the continuity test meaningful
    across a loop.

    Args:
        operation: The motion.

    Returns:
        ``(start, end)``.
    """
    geometry = operation.geometry
    if isinstance(geometry, Circle):
        entry = _circle_entry(operation)
        return entry, entry
    if isinstance(geometry, Arc):
        return _xy(geometry.point_at(0.0)), _xy(geometry.point_at(1.0))
    line = geometry
    assert isinstance(line, Line)  # guaranteed by _motion_kind, which ran first
    return (float(line.start[0]), float(line.start[1])), (float(line.end[0]), float(line.end[1]))


def _tangents(operation: ToolpathOperation) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Unit travel directions where a motion starts and ends.

    Derived from the geometry rather than read off `ToolpathOperation`'s optional
    tangent fields, which not every generator populates: a metric that silently
    skipped the motions whose tangents happen to be None would under-count
    exactly the discontinuities it exists to find.

    Args:
        operation: The motion.

    Returns:
        ``(start_tangent, end_tangent)``, each a unit vector, or ``(0, 0)`` where
        the primitive is degenerate and has no direction.
    """
    geometry = operation.geometry
    if isinstance(geometry, Line):
        direction = _unit(float(geometry.end[0]) - float(geometry.start[0]), float(geometry.end[1]) - float(geometry.start[1]))
        return direction, direction
    centre = geometry.frame.point
    start, end = _endpoints(operation)
    turn = -1.0 if operation.clockwise else 1.0
    return (
        _perpendicular(start[0] - float(centre[0]), start[1] - float(centre[1]), turn),
        _perpendicular(end[0] - float(centre[0]), end[1] - float(centre[1]), turn),
    )


def _perpendicular(dx: float, dy: float, turn: float) -> Tuple[float, float]:
    """The unit travel direction at a point on a circle, given its turn sense.

    Args:
        dx: X offset from the centre to the point.
        dy: Y offset from the centre to the point.
        turn: ``+1`` counterclockwise, ``-1`` clockwise.

    Returns:
        The unit tangent, or ``(0, 0)`` at the centre.
    """
    return _unit(-turn * dy, turn * dx)


def _unit(dx: float, dy: float) -> Tuple[float, float]:
    """Normalise a vector, returning ``(0, 0)`` when it has no length.

    A zero return is a fact about the primitive -- a zero-length move has no
    direction -- and `benchmarks.quality` counts those separately rather than
    treating them as a direction of (1, 0).

    Args:
        dx: X component.
        dy: Y component.

    Returns:
        The unit vector, or ``(0.0, 0.0)``.
    """
    norm = math.hypot(dx, dy)
    if norm == 0.0:
        return (0.0, 0.0)
    return (dx / norm, dy / norm)


def _xy(point: Sequence[float]) -> Tuple[float, float]:
    """Drop a compas point's z, which the cut-plane model fixes."""
    return (float(point[0]), float(point[1]))


def _classify_non_cutting(
    result: ToolpathResult,
    tool_radius: float,
) -> tuple[List[RapidMotion], tuple[OperationIndex, ...], tuple[OperationIndex, ...], float]:
    """Record the motions the depletion replay does not yield.

    `replay_events` skips rapids and reports a plunge only by its foot and its
    into-material verdict, which is right for a depletion but loses the speed
    and longevity groups' inputs. This walk
    reclassifies with the SAME `_replay_kind`, so the two can never disagree
    about what counts as cutting.

    Args:
        result: The generated toolpath.
        tool_radius: Tool radius, for the plunges' bored area.

    Returns:
        Non-cutting observations, plunge indices, retract indices, and plunge
        swept area.

    Raises:
        UnreplayableOperationError: An operation lies outside the cut-plane model.
        UnmeasurableOperationLengthError: An operation's length is undefined.
    """
    cut_z = _infer_cut_height(result.operations)
    rapids: List[RapidMotion] = []
    plunge_indices: list[OperationIndex] = []
    retract_indices: list[OperationIndex] = []
    for index, operation in enumerate(result.operations):
        try:
            category = classify_operation_replay(operation, Millimetre(cut_z))
        except (CutPlaneRampError, OffPlaneReplayCurveError) as error:
            raise UnreplayableOperationError(f"Operation {index} ({operation.operation.value}) {error}.") from error
        if category == "motion":
            continue
        if category == "plunge":
            plunge_indices.append(OperationIndex(index))
            continue
        if category == "retract":
            retract_indices.append(OperationIndex(index))
        rapids.append(
            RapidMotion(
                index=index,
                operation=operation.operation,
                kind=_motion_kind(operation.geometry),
                length=_primitive_length(index, operation),
                horizontal_at_cut_plane=_is_horizontal_at(operation, cut_z),
            )
        )
    typed_plunges = tuple(plunge_indices)
    typed_retracts = tuple(retract_indices)
    return rapids, typed_plunges, typed_retracts, len(typed_plunges) * math.pi * tool_radius**2


def _is_horizontal_at(operation: ToolpathOperation, cut_z: float) -> bool:
    """Whether a motion travels in XY while sitting on the cutting plane.

    Args:
        operation: The motion.
        cut_z: The inferred cutting-plane height.

    Returns:
        True for an XY move at the cutting plane.
    """
    geometry = operation.geometry
    if not isinstance(geometry, Line):
        return float(geometry.frame.point[2]) <= cut_z
    if math.hypot(float(geometry.end[0]) - float(geometry.start[0]), float(geometry.end[1]) - float(geometry.start[1])) == 0.0:
        return False
    return max(float(geometry.start[2]), float(geometry.end[2])) <= cut_z


def _primitive_length(index: int, operation: ToolpathOperation) -> float:
    """Exact length of one operation's primitive.

    compas exposes `length` as a property on `Line` and `Arc` but as an
    unimplemented `Curve` method on `Circle`, whose length is `circumference`;
    reading `.length` uniformly raises `NotImplementedError` on every circle.

    Args:
        index: Position of the operation, for the error message.
        operation: The operation to measure.

    Returns:
        The primitive's length.

    Raises:
        UnmeasurableOperationLengthError: The primitive is not a line, arc, or
            circle.
    """
    geometry = operation.geometry
    if isinstance(geometry, Circle):
        return float(geometry.circumference)
    if isinstance(geometry, (Arc, Line)):
        return float(geometry.length)
    raise UnmeasurableOperationLengthError(f"Operation {index} ({operation.operation.value}) carries geometry {type(geometry).__name__!r}, whose length is not defined.")
