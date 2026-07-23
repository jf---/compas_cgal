"""Toolpath engagement audit: replay a toolpath against a depleting stock model.

An audit takes an emitted `ToolpathResult` and re-runs every motion in order
against an exact `Stock` that is depleted as the tool cuts, measuring the
tool-engagement angle (TEA) each motion actually sees *before* it removes its
own material. The result is an `EngagementReport` -- a per-operation record plus
report-level aggregates -- that states, truthfully and reproducibly, how hard
each move loads the tool and whether any move exceeds the engagement cap.

The audit is a diagnostic, not a regulator: it never rewrites the toolpath, so a
generator with no engagement control is reported as violating rather than
silently flattered. Every geometric decision is delegated to the exact-kernel
backend (`_stock_2.engagement_at`, `certify_segment_tea`); this module owns only
the replay bookkeeping, the cut-plane classification, and -- for circular
motions, which the C++ certifier does not cover -- a fixed-density mirror of the
C++ guarded-station method. Per the boundary doctrine (docs/exactness.md) the
sole transcendental-to-rational conversion on the Python side is isolated in
`_cap_chord_ratio`; the analytic TEA-growth guard is a refinement bound only,
never a geometric truth.
"""

import math
from dataclasses import dataclass

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon
from compas.tolerance import TOL

from compas_cgal import _stock_2  # type: ignore
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# Arc-length station step for circular motions, as a fraction of the FULL-circle
# circumference: 0.05 -> 20 stations per full turn, matched to a partial arc's
# extent so density is scale- and sweep-independent.
AUDIT_ARC_STEP_FRACTION = 0.05

# Operation kinds that engage material at the cutting plane and are therefore
# measured + subtracted. RETRACT and PLUNGE never contribute a TEA measurement
# (rapid up-move / plunge-rate-governed bore); links above the cut plane are
# rapid travel and are skipped by height, not by kind.
AUDIT_ENGAGED = frozenset(
    {
        OperationType.CUT,
        OperationType.LEAD_IN,
        OperationType.LEAD_OUT,
        OperationType.LINK,
    }
)

# Explicit integer safety factor on the analytic TEA-growth bound, mirroring
# TEA_GUARD_SAFETY_FACTOR in src/engagement_2.cpp. Too large a guard only forces
# a conservative "uncertified" verdict; it can never certify a violating motion.
TEA_GUARD_SAFETY_FACTOR = 2


class InvalidEngagementCapError(ValueError):
    """The engagement cap lies outside the admissible open-to-half-turn range ``(0, pi]``."""


class InvalidToolDiameterError(ValueError):
    """The tool diameter is not strictly positive (a real cutter has a positive radius)."""


class UnexpectedToolpathGeometryError(ValueError):
    """A motion's geometry is outside the audit's cut-plane replay model (e.g. a 3D ramp)."""


@dataclass(frozen=True)
class OperationEngagement:
    """Engagement record for a single replayed toolpath operation.

    Attributes:
        op_index: Position of the operation in the source ``ToolpathResult``.
        operation: The operation kind (:class:`OperationType`).
        max_tea: Largest contiguous engaged-run angle (radians) observed over the
            operation's sampled stations; ``0.0`` for unmeasured moves. A
            REPORTING quantity, never a decision input.
        cap_certified: ``True`` iff no cutter center on the motion can exceed the
            cap (exact station verdicts plus the analytic guard), or the move was
            not measured. ``False`` marks a certified violation or an unclosable
            margin.
        stations: Number of stations examined; ``0`` for unmeasured moves.
    """

    op_index: int
    operation: OperationType
    max_tea: float
    cap_certified: bool
    stations: int


@dataclass(frozen=True)
class EngagementReport:
    """Aggregate engagement audit over a whole toolpath.

    Attributes:
        tool_diameter: Tool diameter the audit ran with.
        tea_cap: Engagement-angle cap (radians) the audit certified against.
        operations: Per-operation records, one per source operation and in order.
        max_tea: Largest ``max_tea`` over all operations (radians).
        cap_violations: Number of operations whose cap could not be certified.
        engaged_ops: Number of operations that actually engaged material
            (``max_tea > 0``).
    """

    tool_diameter: float
    tea_cap: float
    operations: list[OperationEngagement]
    max_tea: float
    cap_violations: int
    engaged_ops: int


def _tea_growth_bound(d: float, r: float) -> float:
    """Conservative bound on how far a run's TEA can grow over center travel *d*.

    Four-line Python mirror of ``tea_growth_bound`` in ``src/engagement_2.cpp``
    (the factor-1 analytic lemma): endpoint drift ``4*asin(min(1, d/2r))`` across
    a run's two ends, plus newborn-contact ``2*acos(max(-1, 1 - d/r))`` for a
    feature first biting the rim. Monotone non-decreasing in *d*, saturating at
    ``d = 2r``. Evaluated in doubles purely as a REFINEMENT bound -- it selects
    how much to shrink the cap that the exact station predicate then tests, and
    its safe failure direction (too large -> conservative "uncertified") never
    manufactures a false certificate (docs/exactness.md, "Analytic bounds are
    not precision handling").

    Args:
        d: Euclidean center-travel distance (upper bound suffices).
        r: Tool radius.

    Returns:
        Upper bound on the run's angular growth (radians).
    """
    a = 4.0 * math.asin(min(1.0, d / (2.0 * r)))
    b = 2.0 * math.acos(max(-1.0, 1.0 - d / r))
    return a + b


def _tea_guard(d: float, r: float) -> float:
    """Safety-scaled TEA-growth guard: ``TEA_GUARD_SAFETY_FACTOR * _tea_growth_bound``."""
    return TEA_GUARD_SAFETY_FACTOR * _tea_growth_bound(d, r)


def _cap_chord_ratio(cap_radians: float) -> float:
    """Convert an engagement-angle cap (radians) to the exact squared-chord surrogate.

    BOUNDARY (docs/exactness.md, boundary doctrine): the transcendental cap
    crosses into the exact kernel ONLY as this dimensionless rational surrogate.
    This is the single Python-side transcendental->rational conversion in the
    audit; with ``0 < cap <= pi`` the ratio lies in ``(0, 4]``, exactly the valid
    domain of ``engagement_at``. The sub-ulp gap between this rational and the
    angle the caller typed is documented API semantics, never an in-core
    correction constant.

    Args:
        cap_radians: Engagement-angle cap in ``(0, pi]``.

    Returns:
        ``4 * sin(cap/2)**2`` in ``(0, 4]``.
    """
    return 4.0 * math.sin(cap_radians / 2.0) ** 2


def _subtract_operation(stock: Stock, op: ToolpathOperation, tool_radius: float) -> None:
    """Remove one operation's swept material from *stock* (geometry -> sweep dispatch).

    The single mapping from a :class:`ToolpathOperation`'s geometry to an exact
    stock boolean, reused by the audit and exposed for tools/tests: ``Line`` ->
    capsule of the XY segment; ``Circle`` -> full arc sweep (start == end);
    ``Arc`` -> arc sweep between its endpoints with the operation's turn
    direction. Plunges are handled by the replay (disk), not here.

    Args:
        stock: The depleting stock to subtract from.
        op: The operation whose swept area is removed.
        tool_radius: Tool radius (sweep half-width).

    Raises:
        UnexpectedToolpathGeometryError: If the geometry is not a supported
            ``Line``/``Arc``/``Circle`` primitive.
    """
    g = op.geometry
    if isinstance(g, Line):
        stock.subtract_capsule(
            float(g.start[0]),
            float(g.start[1]),
            float(g.end[0]),
            float(g.end[1]),
            tool_radius,
        )
    elif isinstance(g, Circle):
        c = g.frame.point
        s = g.point_at(0.0)
        stock.subtract_arc_sweep(
            float(c[0]),
            float(c[1]),
            float(s[0]),
            float(s[1]),
            float(s[0]),
            float(s[1]),
            op.clockwise,
            tool_radius,
        )
    elif isinstance(g, Arc):
        c = g.frame.point
        s = g.point_at(0.0)
        e = g.point_at(1.0)
        stock.subtract_arc_sweep(
            float(c[0]),
            float(c[1]),
            float(s[0]),
            float(s[1]),
            float(e[0]),
            float(e[1]),
            op.clockwise,
            tool_radius,
        )
    else:
        raise UnexpectedToolpathGeometryError(f"Cannot subtract operation with geometry type {type(g).__name__!r}.")


def _certify_arc_engagement(stock: Stock, geometry: Arc | Circle, tool_radius: float, tea_cap: float) -> tuple[float, bool, int]:
    """Certify the engagement cap along one circular cut motion by station sampling.

    Mirrors ``certify_segment_tea`` (src/engagement_2.cpp), which the C++
    certifier only implements for line segments, to the circular case: sample
    stations at a fixed arc-length step (``AUDIT_ARC_STEP_FRACTION`` of the full
    turn, matched to the sweep), measure each exactly with ``engagement_at``
    against a GUARDED cap ``tea_cap - guard(half the station spacing)``, and
    report the motion uncertified if any station exceeds the guarded cap or the
    spacing is too coarse to admit a positive guard. Every interior center lies
    within half a station spacing of the nearer station, and arc length upper-
    bounds that Euclidean distance, so the uniform guard is conservative. Unlike
    the C++ certifier the station density is FIXED, not adaptively refined.

    Args:
        stock: The current (frozen for this measurement) stock.
        geometry: The ``Arc`` or ``Circle`` traced by the tool center.
        tool_radius: Tool radius.
        tea_cap: Engagement-angle cap in ``(0, pi]``.

    Returns:
        ``(max_tea, cap_certified, stations)``: the largest engaged-run angle
        over the stations (radians, reporting), the cap verdict, and the station
        count.
    """
    r_arc = float(geometry.radius)
    full_turn = 2.0 * math.pi * r_arc
    step = AUDIT_ARC_STEP_FRACTION * full_turn
    # Arc length from radius and swept angle directly (compas exposes `length` as a
    # property on Arc but a method on Circle). A closed circle wraps (n stations
    # over n intervals, no duplicate seam); an open arc keeps both endpoints (n+1
    # stations over n intervals).
    if isinstance(geometry, Circle):
        arc_len = full_turn
        n_intervals = max(1, math.ceil(arc_len / step))
        params = [i / n_intervals for i in range(n_intervals)]
    else:
        arc_len = r_arc * abs(float(geometry.angle))
        n_intervals = max(1, math.ceil(arc_len / step))
        params = [i / n_intervals for i in range(n_intervals + 1)]
    spacing = arc_len / n_intervals  # uniform arc-length station spacing
    n_stations = len(params)

    # Uniform spacing -> a single guard for every interior center on the motion.
    cap_guarded = tea_cap - _tea_guard(0.5 * spacing, tool_radius)
    if cap_guarded > 0.0:
        ratio = _cap_chord_ratio(cap_guarded)
        cap_certified = True
    else:
        # No positive guarded cap exists at this spacing: unclosable at the fixed
        # density -> uncertified. Measure with the full-cap ratio only to report
        # max_tea; its per-station exceeded flag is intentionally not consulted.
        ratio = _cap_chord_ratio(tea_cap)
        cap_certified = False

    raw = stock.raw
    max_tea = 0.0
    for t in params:
        p = geometry.point_at(t)
        _total_tea, max_run_tea, exceeded = _stock_2.engagement_at(raw, float(p[0]), float(p[1]), tool_radius, ratio)
        if max_run_tea > max_tea:
            max_tea = max_run_tea
        if cap_certified and exceeded:
            cap_certified = False
    return max_tea, cap_certified, n_stations


def _unmeasured(op_index: int, operation: OperationType) -> OperationEngagement:
    """Engagement record for a move that engages no material (retract/plunge/clearance link)."""
    return OperationEngagement(op_index=op_index, operation=operation, max_tea=0.0, cap_certified=True, stations=0)


def _infer_cut_height(operations: list[ToolpathOperation]) -> float:
    """Infer the single cutting-plane z as the minimum motion height.

    A call machines one depth (the toolpath generator's documented single-depth
    contract), so every clearance, plunge, and retract move sits at or above the
    cutting plane and the global minimum endpoint z is that plane. Used only to
    classify horizontal links as cut-height (engaged) versus clearance-height
    (rapid, skipped).

    Args:
        operations: The operations to scan.

    Returns:
        The cutting-plane z (``0.0`` for an empty operation list).
    """
    heights: list[float] = []
    for op in operations:
        g = op.geometry
        if isinstance(g, Line):
            heights.append(float(g.start[2]))
            heights.append(float(g.end[2]))
        else:
            heights.append(float(g.frame.point[2]))
    return min(heights) if heights else 0.0


def _replay_line(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float, cut_z: float) -> OperationEngagement:
    """Replay a single linear operation: classify by height, measure + subtract if engaged."""
    operation = op.operation
    g = op.geometry
    assert isinstance(g, Line)
    z0, z1 = float(g.start[2]), float(g.end[2])
    xy_len = math.hypot(float(g.end[0]) - float(g.start[0]), float(g.end[1]) - float(g.start[1]))

    if abs(z0 - z1) > TOL.absolute:
        # Differing-z line: a plunge (down) or retract-shaped (up). The generator
        # emits both as pure-vertical lines; a differing-z line that also travels
        # in XY is a 3D ramp the cut-plane audit does not model.
        if xy_len > TOL.absolute:
            raise UnexpectedToolpathGeometryError(
                f"Operation {op_index} ({operation.value}) is a differing-z line with nonzero XY travel "
                f"(len={xy_len:.3e}); ramped 3D cutting is outside the engagement audit's cut-plane model."
            )
        if z1 < z0:
            # Plunge: full-immersion bore. Remove the tool disk at the end point;
            # plunge feed is governed by plunge-rate limits, not TEA, so the audit
            # deliberately records no engagement for it.
            stock.subtract_disk(float(g.end[0]), float(g.end[1]), tool_radius)
        # Upward (retract-shaped) moves remove no new material.
        return _unmeasured(op_index, operation)

    # Horizontal line at height z0.
    if z0 > cut_z + TOL.absolute:
        # Above the cutting plane (a link across an already-cleared corridor at the
        # clearance height): rapid travel, no material interaction.
        return _unmeasured(op_index, operation)

    if operation not in AUDIT_ENGAGED:
        return _unmeasured(op_index, operation)

    # Cut-height engaged linear motion. Certify against the material this move
    # actually sees, THEN deplete -- measuring before subtraction is what makes a
    # slotting cut report full immersion instead of its own post-cut void.
    max_tea, cap_certified, stations = _stock_2.certify_segment_tea(
        stock.raw,
        float(g.start[0]),
        float(g.start[1]),
        float(g.end[0]),
        float(g.end[1]),
        tool_radius,
        tea_cap,
    )
    _subtract_operation(stock, op, tool_radius)
    return OperationEngagement(op_index, operation, max_tea, cap_certified, stations)


def _replay_arc(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float) -> OperationEngagement:
    """Replay a single circular operation: measure (station sampling) + subtract if engaged."""
    operation = op.operation
    if operation not in AUDIT_ENGAGED:
        return _unmeasured(op_index, operation)
    geometry = op.geometry
    assert isinstance(geometry, (Arc, Circle))  # guaranteed by _replay_operation's Line dispatch
    max_tea, cap_certified, stations = _certify_arc_engagement(stock, geometry, tool_radius, tea_cap)
    _subtract_operation(stock, op, tool_radius)
    return OperationEngagement(op_index, operation, max_tea, cap_certified, stations)


def _replay_operation(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float, cut_z: float) -> OperationEngagement:
    """Replay one operation against the depleting stock, returning its engagement record."""
    if op.operation == OperationType.RETRACT:
        # Rapid clearance-plane up-move: no material interaction regardless of geometry.
        return _unmeasured(op_index, op.operation)
    if isinstance(op.geometry, Line):
        return _replay_line(stock, op_index, op, tool_radius, tea_cap, cut_z)
    # Arc / Circle: planar cut motions, always at the cutting plane by construction.
    return _replay_arc(stock, op_index, op, tool_radius, tea_cap)


def audit_toolpath_engagement(
    polygon: Polygon,
    result: ToolpathResult,
    tool_diameter: float,
    tea_cap: float,
    holes: list[Polygon] | None = None,
) -> EngagementReport:
    """Replay a toolpath against a depleting stock and report its engagement.

    Builds an exact `Stock` from *polygon* (and optional *holes*), then replays
    ``result.operations`` in order: each engaged motion at the cutting plane is
    measured for tool-engagement angle against the CURRENT stock and then removes
    its own swept material, so every measurement reflects the material that move
    truly sees. Retracts, plunges, and clearance-height links interact with no
    material and are recorded as unmeasured. The audit never modifies the
    toolpath; an unregulated generator is reported as violating.

    Args:
        polygon: Outer stock boundary in the world XY plane.
        result: The toolpath to audit (typed operations from the generator).
        tool_diameter: Tool diameter; the tool radius is half of this.
        tea_cap: Engagement-angle cap in radians, ``0 < tea_cap <= pi`` (a single
            engaged run subtends at most a half turn before the ``> pi`` case is an
            exact orientation verdict, so a larger cap is meaningless).
        holes: Optional island polygons strictly inside *polygon*.

    Returns:
        EngagementReport: Per-operation engagement records and report-level
        aggregates (peak TEA, cap-violation count, engaged-operation count).

    Raises:
        InvalidEngagementCapError: If *tea_cap* is not in ``(0, pi]``.
        InvalidToolDiameterError: If *tool_diameter* is not strictly positive.
        UnexpectedToolpathGeometryError: If an operation's geometry is outside the
            cut-plane replay model.
        InvalidPolygonError: If *polygon* or a hole is non-planar or degenerate.
    """
    if not (tea_cap > 0.0 and tea_cap <= math.pi):
        raise InvalidEngagementCapError(f"tea_cap must be in (0, pi]; got {tea_cap!r}.")
    if not tool_diameter > 0.0:
        raise InvalidToolDiameterError(f"tool_diameter must be strictly positive; got {tool_diameter!r}.")

    tool_radius = 0.5 * tool_diameter
    stock = Stock(polygon, holes=holes)
    cut_z = _infer_cut_height(result.operations)

    operations = [_replay_operation(stock, i, op, tool_radius, tea_cap, cut_z) for i, op in enumerate(result.operations)]

    max_tea = max((e.max_tea for e in operations), default=0.0)
    cap_violations = sum(1 for e in operations if not e.cap_certified)
    engaged_ops = sum(1 for e in operations if e.max_tea > 0.0)

    return EngagementReport(
        tool_diameter=tool_diameter,
        tea_cap=tea_cap,
        operations=operations,
        max_tea=max_tea,
        cap_violations=cap_violations,
        engaged_ops=engaged_ops,
    )
