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
from typing import Iterable
from typing import Never
from typing import NoReturn

from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.helix import Helix
from compas_cgal.replay_classification import AUDIT_ENGAGED as AUDIT_ENGAGED
from compas_cgal.replay_classification import CutPlaneRampError
from compas_cgal.replay_classification import OffPlaneHelixError
from compas_cgal.replay_classification import OffPlaneReplayCurveError
from compas_cgal.replay_classification import classify_operation_replay
from compas_cgal.replay_classification import infer_operation_cut_height
from compas_cgal.replay_classification import minimum_cut_height
from compas_cgal.stock import Stock
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import Millimetre

# Arc-length station step for circular motions, as a fraction of the FULL-circle
# circumference: 0.05 -> 20 stations per full turn, matched to a partial arc's
# extent so density is scale- and sweep-independent.
AUDIT_ARC_STEP_FRACTION = 0.05

# Operation kinds that engage material at the cutting plane and are therefore
# measured + subtracted. RETRACT and PLUNGE never contribute a TEA measurement
# (rapid up-move / plunge-rate-governed bore); links above the cut plane are
# rapid travel and are skipped by height, not by kind.
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


class UnhandledReplayCategoryError(RuntimeError):
    """A replay category reached a dispatch that has no branch for it."""


def _unhandled_category(category: Never) -> NoReturn:
    """Fail on a replay category no dispatch branch handles.

    Typed ``Never`` like ``typing.assert_never``: every category dispatch is an
    exhaustive ``==`` chain over `ReplayCategory` ending here, so ``mypy
    --strict`` rejects a chain missing a member; at runtime it raises a named
    error, never a silent default.

    Raises:
        UnhandledReplayCategoryError: Always.
    """
    raise UnhandledReplayCategoryError(f"no dispatch branch for replay category {category!r}")


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
    stock boolean, reused by the audit and exposed for tools/tests, each removed
    by LOCAL arrangement surgery (docs/local_depletion.md):

    - ``Line`` -> `Stock.subtract_line_sweep_local`: the exact capsule when the
      direction is Pythagorean, else exact end disks + a certified under-rectangle.
    - ``Circle`` -> `Stock.subtract_circle_sweep_local` with the AUTHORED
      ``g.radius``: the exact annulus, independent of where the loop starts.
    - ``Arc`` -> nothing when the authored radius's full-turn annulus is
      already swept (Ruling 1.1-C2, B-minimal); else the arc sweep between its
      endpoints with the operation's turn direction (a certified
      under-covering disk chain).
    - ``Helix`` -> `Stock.subtract_circle_sweep_local` on the footprint: the
      bore, the exact full-turn sweep of the footprint loop. It is the disk
      of radius ``rho_h + r`` when ``rho_h <= r`` (enforced by
      `compas_cgal.pipeline.pencil.helix_footprint`); for ``rho_h > r`` the
      replay removes the true annulus.

    Plunges are handled by the replay (disk), not here.

    Args:
        stock: The depleting stock to subtract from.
        op: The operation whose swept area is removed.
        tool_radius: Tool radius (sweep half-width).

    Raises:
        UnexpectedToolpathGeometryError: If the geometry is not a supported
            ``Line``/``Arc``/``Circle``/``Helix`` primitive.
    """
    g = op.geometry
    if isinstance(g, Helix):
        footprint = g.footprint
        stock.subtract_circle_sweep_local(float(footprint.centre.x), float(footprint.centre.y), float(footprint.radius.value), tool_radius)
    elif isinstance(g, Line):
        stock.subtract_line_sweep_local(float(g.start[0]), float(g.start[1]), float(g.end[0]), float(g.end[1]), tool_radius)
    elif isinstance(g, Circle):
        c = g.frame.point
        stock.subtract_circle_sweep_local(float(c[0]), float(c[1]), float(g.radius), tool_radius)
    elif isinstance(g, Arc):
        c = g.frame.point
        if stock.circle_sweep_removal(float(c[0]), float(c[1]), float(g.radius), tool_radius) is SweepRemoval.NOTHING:
            # The full-turn annulus of the authored radius is already swept, so
            # the arc removes nothing and the stock stays exactly as it is
            # (Ruling 1.1-C2, B-minimal).
            return
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


def _operation_removal(stock: Stock, op: ToolpathOperation, tool_radius: float) -> SweepRemoval:
    """What `_subtract_operation` would remove from *stock*, decided without removing it.

    The exact three-valued bracket (`SweepRemoval`) between the region the replay
    subtracts and a certified region containing the true swept area:

    - ``Line`` -> `Stock.line_sweep_removal`: two-valued for an exact capsule,
      three-valued for a canonical quad.
    - ``Circle`` -> `Stock.circle_sweep_removal` on the authored radius: the
      annulus is exact, so REMOVES or NOTHING.
    - ``Arc`` -> NOTHING first, decided on the full-turn annulus of the
      authored radius, which contains the arc's true sweep (Ruling 1.1-C2,
      B-minimal: an air arc on an already swept annulus never reads REMOVES on
      an ulp sliver of the replay's surrogate chain); otherwise the chain the
      replay removes decides REMOVES (one clone and an exact equality); else
      WITHIN_BAND.
    - ``Helix`` -> `Stock.circle_sweep_removal` on the footprint: the bore is
      exact, so REMOVES or NOTHING.

    *stock* is never mutated.

    Args:
        stock: The stock the operation is about to cut.
        op: The operation.
        tool_radius: Tool radius (sweep half-width).

    Returns:
        The verdict.

    Raises:
        UnexpectedToolpathGeometryError: If the geometry is not a supported
            ``Line``/``Arc``/``Circle``/``Helix`` primitive.
    """
    g = op.geometry
    if isinstance(g, Helix):
        footprint = g.footprint
        return stock.circle_sweep_removal(float(footprint.centre.x), float(footprint.centre.y), float(footprint.radius.value), tool_radius)
    if isinstance(g, Line):
        return stock.line_sweep_removal(float(g.start[0]), float(g.start[1]), float(g.end[0]), float(g.end[1]), tool_radius)
    if isinstance(g, Circle):
        c = g.frame.point
        return stock.circle_sweep_removal(float(c[0]), float(c[1]), float(g.radius), tool_radius)
    if isinstance(g, Arc):
        c = g.frame.point
        if stock.circle_sweep_removal(float(c[0]), float(c[1]), float(g.radius), tool_radius) is SweepRemoval.NOTHING:
            return SweepRemoval.NOTHING
        probe = stock.clone()
        _subtract_operation(probe, op, tool_radius)
        if not probe.exactly_equals(stock):
            return SweepRemoval.REMOVES
        return SweepRemoval.WITHIN_BAND
    raise UnexpectedToolpathGeometryError(f"Cannot decide the removal of geometry type {type(g).__name__!r}.")


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

    MERGE-CLOSURE WIRING (mirror of ``certify_segment_tea``). The station cap
    decision is measured with ``gap_close_ratio = 4*sin^2(min(gamma_guard, pi)/2)``,
    where ``gamma_guard == _tea_guard(half spacing) == 2*GROWTH(hs)`` -- the SAME
    angle subtracted from the cap. This pre-absorbs every void gap that could close
    within one station step, so a run-merge (two engaged runs fusing across a thin
    closing rim void -- an O(1) jump in max-run the O(sqrt d) growth lemma cannot
    bridge) is seen at the station instead of slipping between stations and
    false-passing. SAFE FAILURE DIRECTION: only ``GROWTH(hs)`` is strictly required
    to bound the closable span; ``gamma_guard`` uses the factor-2 ``_tea_guard`` and
    thus OVER-closes, which only enlarges the pessimistic runs -> a stricter test ->
    at most an extra (conservative) refusal, never a false certificate. The reported
    ``max_tea`` is the TRUE (unmerged) engaged-run measure regardless of the closure
    (the deciding/reporting split lives inside ``engagement_at``).

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
    # gamma_guard == 2*GROWTH(hs): TEA_GUARD_SAFETY_FACTOR == 2, so _tea_guard already
    # IS 2*GROWTH -- the merge-closure angle is the SAME quantity subtracted from the cap.
    gamma_guard = _tea_guard(0.5 * spacing, tool_radius)
    cap_guarded = tea_cap - gamma_guard
    if cap_guarded > 0.0:
        ratio = _cap_chord_ratio(cap_guarded)
        # Pre-absorb every void gap of span <= gamma_guard so a run-merge completing
        # within a station step is never invisible to the cap decision. gamma_guard
        # bounds a gap's closable span (its endpoints are run endpoints drifting under
        # the same growth moduli), so each station's PESSIMISTIC max-run already
        # accounts for any merge before the next station -- the growth lemma then
        # bridges stations with no merge term. Mirrors certify_segment_tea's
        # gap_close_ratio = 4*sin^2(min(gamma_guard, pi)/2) (src/engagement_2.cpp).
        gap_close_ratio = _cap_chord_ratio(min(gamma_guard, math.pi))
        cap_certified = True
    else:
        # No positive guarded cap exists at this spacing: unclosable at the fixed
        # density -> uncertified. Measure with the full-cap ratio only to report
        # max_tea; its per-station exceeded flag is intentionally not consulted. No
        # pessimism (mirrors the C++ guard>=cap branch passing gap_close_ratio 0.0).
        ratio = _cap_chord_ratio(tea_cap)
        gap_close_ratio = 0.0
        cap_certified = False

    raw = stock.raw
    max_tea = 0.0
    for t in params:
        p = geometry.point_at(t)
        _total_tea, max_run_tea, exceeded = _stock_2.engagement_at(raw, float(p[0]), float(p[1]), tool_radius, ratio, gap_close_ratio)
        if max_run_tea > max_tea:
            max_tea = max_run_tea
        if cap_certified and exceeded:
            cap_certified = False
    return max_tea, cap_certified, n_stations


def _unmeasured(op_index: int, operation: OperationType) -> OperationEngagement:
    """Engagement record for a move that is not cap-sampled.

    Rapids and retracts engage no material. A plunge and a helix bore their
    disk but are not cap-sampled by design (E8): a bore into virgin stock is
    full immersion by construction, and the motion after it is measured
    against the bored stock.
    """
    return OperationEngagement(op_index=op_index, operation=operation, max_tea=0.0, cap_certified=True, stations=0)


def _minimum_cut_height(heights: Iterable[float]) -> float:
    """Return the lowest motion height, or zero when no heights exist."""
    return float(minimum_cut_height(Millimetre(height) for height in heights))


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
    return float(infer_operation_cut_height(operations))


def _replay_line(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float) -> OperationEngagement:
    """Replay one engaged cut-plane linear motion: certify against the current stock, then subtract."""
    g = op.geometry
    assert isinstance(g, Line)  # guaranteed by _replay_motion's dispatch
    # Certify against the material this move actually sees, THEN deplete --
    # measuring before subtraction is what makes a slotting cut report full
    # immersion instead of its own post-cut void.
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
    return OperationEngagement(op_index, op.operation, max_tea, cap_certified, stations)


def _replay_arc(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float) -> OperationEngagement:
    """Replay one engaged cut-plane circular motion: measure (station sampling), then subtract."""
    geometry = op.geometry
    assert isinstance(geometry, (Arc, Circle))  # guaranteed by _replay_motion's dispatch
    max_tea, cap_certified, stations = _certify_arc_engagement(stock, geometry, tool_radius, tea_cap)
    _subtract_operation(stock, op, tool_radius)
    return OperationEngagement(op_index, op.operation, max_tea, cap_certified, stations)


def _replay_motion(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float) -> OperationEngagement:
    """Replay one ``"motion"`` operation by its geometry.

    Raises:
        UnexpectedToolpathGeometryError: The motion is not a ``Line``, ``Arc``
            or ``Circle``.
    """
    if isinstance(op.geometry, Line):
        return _replay_line(stock, op_index, op, tool_radius, tea_cap)
    if isinstance(op.geometry, (Arc, Circle)):
        return _replay_arc(stock, op_index, op, tool_radius, tea_cap)
    raise UnexpectedToolpathGeometryError(f"Operation {op_index} ({op.operation.value}) is a motion with geometry type {type(op.geometry).__name__!r}.")


def _replay_operation(stock: Stock, op_index: int, op: ToolpathOperation, tool_radius: float, tea_cap: float, cut_z: float) -> OperationEngagement:
    """Replay one operation against the depleting stock, returning its engagement record.

    An exhaustive dispatch over `ReplayCategory`: ``"motion"`` is measured and
    subtracted; ``"plunge"`` bores the tool disk at its foot and ``"helix"``
    its exact bore (`_subtract_operation`), both unmeasured by design (E8);
    ``"rapid"`` and ``"retract"`` remove nothing.

    Raises:
        UnexpectedToolpathGeometryError: The operation lies outside the
            cut-plane model.
        UnhandledReplayCategoryError: The category has no branch here.
    """
    try:
        category = classify_operation_replay(op, Millimetre(cut_z))
    except (CutPlaneRampError, OffPlaneReplayCurveError, OffPlaneHelixError) as error:
        raise UnexpectedToolpathGeometryError(f"Operation {op_index} ({op.operation.value}) {error}.") from error
    if category == "motion":
        return _replay_motion(stock, op_index, op, tool_radius, tea_cap)
    if category == "plunge":
        g = op.geometry
        assert isinstance(g, Line)  # only a line classifies as a plunge
        stock.subtract_disk(float(g.end[0]), float(g.end[1]), tool_radius)
        return _unmeasured(op_index, op.operation)
    if category == "helix":
        assert isinstance(op.geometry, Helix)  # only a helix classifies as a helix
        _subtract_operation(stock, op, tool_radius)
        return _unmeasured(op_index, op.operation)
    if category == "rapid" or category == "retract":
        return _unmeasured(op_index, op.operation)
    _unhandled_category(category)


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
    truly sees. Retracts and clearance-height links interact with no material;
    plunges and helices bore their disk without being cap-sampled (E8). All
    four are recorded as unmeasured. The audit never modifies the
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
        UnhandledReplayCategoryError: If a replay category has no dispatch branch.
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
