"""Properties every machining metric must satisfy, for every path.

WHY THIS FILE EXISTS. Three defects were found in `benchmarks.quality` in a
single session, and all three failed the same way: a metric was validated
against the generators that existed when it was written, so those generators'
HABITS were encoded as physics. Each defect surfaced only when a path of a new
KIND arrived.

| assumption baked in | the path that broke it | what was reported |
| --- | --- | --- |
| the generator retracts between chains | one that links at cutting depth | a 3.988 load step taken during a long transit |
| the loop list is one chain | a path with five | a step across a rapid, on every pocket |
| the angle is the textbook one | none -- it was a misreading | a "correction" that would call a plunge zero immersion |

Example-based tests cannot find this class, because the example is written by
the same person holding the assumption. What finds it is a property stated over
ALL paths. Three families are asserted here:

* INVARIANCE UNDER RIGID MOTION. A pocket machined in a different place, or
  turned on the table, is the same machining problem. Any metric that moves is
  reading a coordinate.
* COVARIANCE UNDER SCALE. Doubling the pocket, the tool and the path doubles
  every length and changes no ratio. This is the family that catches a unit or
  convention slip, which is what the angle misreading was.
* THE RUN SEMANTICS. `_loop_runs` decides which loop pairs may be compared, and
  its two boundaries are asserted here in both directions: merging chains may
  only ever raise the reported step, inserting a retract may only ever lower it.

Kept separate from `test_quality.py`, which pins specific measured values. That
file says what the numbers ARE; this one says what they must always OBEY.
"""

from __future__ import annotations

from typing import List
from typing import Sequence
from typing import Tuple

import pytest
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import numpy as np

from benchmarks.quality import _cut
from benchmarks.quality import _elementary
from benchmarks.quality import _program
from benchmarks.quality import _speed
from benchmarks.quality import radial_immersion
from benchmarks.survey import survey_path
from benchmarks.spec import PocketSpec
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

# The reference instance every property is stated against, at unit scale: a
# 12x8 pocket and a 2.0 cutter, matching `test_quality.py` so a failure here and
# a failure there describe the same geometry.
BASE_WIDTH = 12.0
BASE_HEIGHT = 8.0
BASE_TOOL_DIAMETER = 2.0
BASE_CAP_DEG = 120.0

# Few probes per motion: these properties compare two measurements of the same
# path, so sampling density cancels and only its consistency matters.
SAMPLES = 4

# Ratios that must be identical at any scale, read off `PathQuality` by path.
DIMENSIONLESS = (
    ("cut", "immersion_steady_fraction"),
    ("cut", "immersion_at_design_fraction"),
    ("cut", "loop_radius_cv"),
    ("cut", "max_loop_radius_step"),
    ("cut", "max_engagement_deg"),
    ("cut", "engagement_p95_deg"),
    ("cut", "max_chip_thickness_ratio"),
    ("speed", "air_fraction"),
    ("program", "arc_length_fraction"),
    ("program", "block_length_cv"),
)

# Quantities that must scale exactly with the geometry.
LENGTHS = (
    ("speed", "cutting_length"),
    ("speed", "air_length"),
    ("program", "min_block_length"),
    ("program", "median_block_length"),
)

# Counts that must not change at all under any of these transforms.
COUNTS = (
    ("cut", "immersion_excursions"),
    ("cut", "cap_exceedances"),
    ("cut", "slotting_motions"),
    ("speed", "tangent_breaks"),
    ("program", "block_count"),
    ("program", "cut_blocks"),
    ("elementary", "degenerate_loops"),
    ("elementary", "redundant_operations"),
)

# A tolerance that absorbs floating-point rework of the same geometry through a
# different coordinate frame, and nothing larger. Both measurements run the same
# code on the same shape, so agreement should be near machine precision; 1e-9
# leaves four orders of headroom over that without admitting a real difference.
REWORK_TOL = 1e-9


@st.composite
def chained_paths(draw: st.DrawFn) -> List[Tuple[int, float, float, float]]:
    """A plausible trochoidal path as ``(chain, x, y, radius)`` stations.

    Shaped like what the generators actually emit -- a handful of chains, each a
    run of machining circles joined by bridges -- rather than arbitrary
    geometry, because a property is only useful if the paths it quantifies over
    are ones the metrics will really meet.

    Args:
        draw: Hypothesis' draw function.

    Returns:
        Stations in toolpath order.
    """
    chains = draw(st.integers(min_value=1, max_value=3))
    stations: List[Tuple[int, float, float, float]] = []
    for chain in range(chains):
        count = draw(st.integers(min_value=2, max_value=4))
        y = draw(st.floats(min_value=-2.0, max_value=2.0, allow_nan=False, allow_infinity=False))
        x0 = draw(st.floats(min_value=-4.0, max_value=0.0, allow_nan=False, allow_infinity=False))
        for step in range(count):
            radius = draw(st.floats(min_value=0.3, max_value=2.5, allow_nan=False, allow_infinity=False))
            stations.append((chain, x0 + 0.9 * step, y, radius))
    return stations


def _spec(scale: float, quarter_turns: int, shift: Tuple[float, float]) -> PocketSpec:
    """The reference pocket, scaled then rigidly moved.

    Args:
        scale: Uniform scale factor applied to the pocket and the tool.
        quarter_turns: Whole 90-degree turns about the origin.
        shift: Translation applied after rotation.

    Returns:
        The transformed instance.
    """
    half_w, half_h = 0.5 * BASE_WIDTH * scale, 0.5 * BASE_HEIGHT * scale
    corners = ((-half_w, -half_h), (half_w, -half_h), (half_w, half_h), (-half_w, half_h))
    points = [_place(x, y, quarter_turns, shift) for x, y in corners]
    return PocketSpec.build(
        name=f"invariant_{scale}_{quarter_turns}",
        family="analytic",
        polygon=Polygon([[x, y, 0.0] for x, y in points]),
        tool_diameter=BASE_TOOL_DIAMETER * scale,
        tea_cap_deg=BASE_CAP_DEG,
    )


def _place(x: float, y: float, quarter_turns: int, shift: Tuple[float, float]) -> Tuple[float, float]:
    """Rotate ``(x, y)`` by whole quarter turns, then translate.

    QUARTER TURNS RATHER THAN AN ARBITRARY ANGLE, AND NOT FOR CONVENIENCE. A
    quarter turn is a coordinate swap and a sign flip, so it is exact in binary
    and maps an axis-aligned pocket to an axis-aligned pocket. Any other angle
    makes the pocket OBLIQUE, and an oblique pocket cannot be measured here in
    reasonable time: `_coverage_2.ReachableDomain2(...).center_domain()` takes
    5 ms on an axis-aligned rectangle with integer vertices and 17.5 SECONDS on
    an oblique quadrilateral with equally plain integer vertices, rising past a
    minute once the vertices carry decimals.

    So this file cannot assert rotation invariance in general, and says so
    rather than quietly testing the easy case. What it does assert -- that a
    quarter turn changes nothing -- still catches a metric that reads an axis,
    which is the defect class the property exists for. See
    `docs/machining_metric_validity.md` for the measurement.

    Args:
        x: Abscissa.
        y: Ordinate.
        quarter_turns: Number of 90-degree turns, counter-clockwise.
        shift: Translation applied after the rotation.

    Returns:
        The transformed point.
    """
    for _ in range(quarter_turns % 4):
        x, y = -y, x
    return (x + shift[0], y + shift[1])


def _frame(x: float, y: float, quarter_turns: int) -> Frame:
    """A circle frame whose AXES are turned with the pocket, not just its origin.

    `Circle.point_at(0)` is ``centre + radius * frame.xaxis``, and that point is
    where the survey reads a loop's entry and therefore its tangents. Turning
    only the origin leaves the seam on the +x side while everything else turns,
    which is not a rotation of the path -- it is a different path, with the loop
    entered at a different place. An earlier draft did exactly that and the
    rotation property duly failed on `tangent_breaks`, correctly reporting that
    the two paths were not the same.

    Args:
        x: Centre abscissa.
        y: Centre ordinate.
        quarter_turns: Whole 90-degree turns, matching `_place`.

    Returns:
        The frame.
    """
    axis_x = _place(1.0, 0.0, quarter_turns, (0.0, 0.0))
    axis_y = _place(0.0, 1.0, quarter_turns, (0.0, 0.0))
    return Frame([x, y, 0.0], [axis_x[0], axis_x[1], 0.0], [axis_y[0], axis_y[1], 0.0])


def _result(
    stations: Sequence[Tuple[int, float, float, float]],
    *,
    scale: float = 1.0,
    quarter_turns: int = 0,
    shift: Tuple[float, float] = (0.0, 0.0),
    merge_chains: bool = False,
    retract_between_chains: bool = True,
) -> ToolpathResult:
    """Build a toolpath from *stations* under a transform and a linking policy.

    Args:
        stations: ``(chain, x, y, radius)`` in toolpath order.
        scale: Uniform scale applied to every coordinate and radius.
        quarter_turns: Whole 90-degree turns about the origin.
        shift: Translation applied after rotation.
        merge_chains: Label every operation with one ``path_index``, collapsing
            the chain boundaries without moving any geometry.
        retract_between_chains: Emit a retract and a plunge at each chain change
            rather than a cut-depth link.

    Returns:
        The toolpath.
    """
    operations: List[ToolpathOperation] = []
    previous_chain: int | None = None
    previous_point: Tuple[float, float] | None = None
    for chain, raw_x, raw_y, raw_radius in stations:
        x, y = _place(raw_x * scale, raw_y * scale, quarter_turns, shift)
        radius = raw_radius * scale
        index = 0 if merge_chains else chain
        if previous_chain is None:
            operations.append(_op(Line([x, y, 4.0 * scale], [x, y, 0.0]), OperationType.PLUNGE, index))
        elif chain != previous_chain:
            if retract_between_chains:
                assert previous_point is not None
                operations.append(_op(Line([*previous_point, 0.0], [*previous_point, 4.0 * scale]), OperationType.RETRACT, index))
                operations.append(_op(Line([x, y, 4.0 * scale], [x, y, 0.0]), OperationType.PLUNGE, index))
            else:
                assert previous_point is not None
                operations.append(_op(Line([*previous_point, 0.0], [x, y, 0.0]), OperationType.LINK, index))
        else:
            assert previous_point is not None
            operations.append(_op(Line([*previous_point, 0.0], [x, y, 0.0]), OperationType.LINK, index))
        operations.append(_op(Circle(radius, frame=_frame(x, y, quarter_turns)), OperationType.CUT, index))
        previous_chain, previous_point = chain, (x, y)
    return ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float))


def _op(geometry: object, operation: OperationType, path_index: int) -> ToolpathOperation:
    """One operation on a named chain.

    Args:
        geometry: The primitive.
        operation: Its role.
        path_index: The chain it belongs to.

    Returns:
        The operation.
    """
    return ToolpathOperation(geometry=geometry, operation=operation, path_index=path_index)


def _field(quality: _Groups, group: str, name: str) -> float:
    """Read one metric off a `PathQuality` by group and field name.

    Args:
        quality: The measurement.
        group: Group attribute, such as ``"cut"``.
        name: Field within the group.

    Returns:
        The value.
    """
    return float(getattr(getattr(quality, group), name))


class _Groups:
    """The four groups a property needs, measured without the coverage grid.

    `measure_quality` also runs a coverage grid, which at the resolution the
    residue guard demands is some seventeen thousand exact point-location
    queries and dominates the cost by two orders of magnitude. None of the
    properties here assert a coverage-derived field -- they compare a path
    against a transformed copy of itself, and coverage is the one part that
    would only restate the transform. So the groups are built straight off the
    survey and the coverage inputs are passed as zero.

    That makes `uncut_fraction`, `recut_fraction` and `wall_scallop_height`
    MEANINGLESS on this object, and no property may read them. The fields the
    properties do read are listed in `DIMENSIONLESS`, `LENGTHS` and `COUNTS`,
    none of which touch the grid.

    Attributes:
        elementary: Validity, with the two coverage fields invalid.
        cut: Cut mechanics.
        speed: Machine cost.
        program: Program feasibility.
    """

    def __init__(self, spec: PocketSpec, result: ToolpathResult) -> None:
        """Measure *result* on *spec*, skipping coverage.

        Args:
            spec: The instance.
            result: The toolpath.
        """
        survey = survey_path(spec, result, samples_per_motion=SAMPLES)
        chain_of = {index: operation.path_index for index, operation in enumerate(result.operations)}
        self.elementary = _elementary(spec, survey, 0.0, 0.0)
        self.cut = _cut(spec, survey, 0.0, chain_of)
        self.speed = _speed(survey)
        self.program = _program(survey)


def _measure(spec: PocketSpec, result: ToolpathResult) -> _Groups:
    """Measure a path. Any refusal is a FAILURE, never a skipped example.

    Deliberately has no error handling. An earlier draft caught `BenchmarkError`
    and skipped, and it hid a real problem immediately: the grid sat exactly on
    the coverage guard's boundary, every example was discarded, and the property
    reported green while asserting nothing at all. A property test that quietly
    discards its inputs is worse than no test, because it also reports success.

    Args:
        spec: The instance.
        result: The toolpath.

    Returns:
        The groups.
    """
    return _Groups(spec, result)


PROPERTY_SETTINGS = settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture])


@PROPERTY_SETTINGS
@given(stations=chained_paths(), shift=st.tuples(st.floats(-50.0, 50.0), st.floats(-50.0, 50.0)))
def test_moving_the_pocket_across_the_table_changes_no_metric(stations: List[Tuple[int, float, float, float]], shift: Tuple[float, float]) -> None:
    """A pocket machined somewhere else is the same machining problem.

    Any metric that moves under translation is reading an absolute coordinate,
    which nothing about cutting metal depends on.
    """
    here = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations))
    there = _measure(_spec(1.0, 0, shift), _result(stations, shift=shift))
    for group, name in DIMENSIONLESS + LENGTHS:
        assert _field(here, group, name) == pytest.approx(_field(there, group, name), rel=REWORK_TOL, abs=REWORK_TOL), f"{group}.{name} moved with the pocket"
    for group, name in COUNTS:
        assert _field(here, group, name) == _field(there, group, name), f"{group}.{name} moved with the pocket"


@PROPERTY_SETTINGS
@given(stations=chained_paths())
def test_turning_the_pocket_a_quarter_turn_changes_no_metric(stations: List[Tuple[int, float, float, float]]) -> None:
    """Machining is not a function of which way the stock sits on the table.

    A quarter turn is used rather than an arbitrary angle so the bounding box
    keeps its longer side, and the coverage grid therefore resolves the two
    measurements identically. A metric that survives translation but not
    rotation is reading an axis.
    """
    upright = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations))
    turned = _measure(_spec(1.0, 1, (0.0, 0.0)), _result(stations, quarter_turns=1))
    for group, name in DIMENSIONLESS + LENGTHS:
        assert _field(upright, group, name) == pytest.approx(_field(turned, group, name), rel=1e-6, abs=1e-6), f"{group}.{name} turned with the pocket"
    for group, name in COUNTS:
        assert _field(upright, group, name) == _field(turned, group, name), f"{group}.{name} turned with the pocket"


@PROPERTY_SETTINGS
@given(stations=chained_paths(), scale=st.floats(min_value=0.25, max_value=4.0))
def test_a_bigger_pocket_and_a_bigger_cutter_is_the_same_cut(stations: List[Tuple[int, float, float, float]], scale: float) -> None:
    """Every length scales, every ratio holds, every count is unchanged.

    THE FAMILY THAT CATCHES A CONVENTION SLIP. A metric that mixes an angle with
    a length, or a radius with a diameter, generally survives translation and
    rotation and fails here -- because only a scale change separates a ratio
    that is genuinely dimensionless from one that merely looked it.
    """
    unit = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations))
    scaled = _measure(_spec(scale, 0, (0.0, 0.0)), _result(stations, scale=scale))
    for group, name in DIMENSIONLESS:
        assert _field(unit, group, name) == pytest.approx(_field(scaled, group, name), rel=1e-6, abs=1e-6), f"{group}.{name} is not dimensionless"
    for group, name in LENGTHS:
        assert scale * _field(unit, group, name) == pytest.approx(_field(scaled, group, name), rel=1e-6, abs=1e-6), f"{group}.{name} did not scale"
    for group, name in COUNTS:
        assert _field(unit, group, name) == _field(scaled, group, name), f"{group}.{name} changed with scale"


@PROPERTY_SETTINGS
@given(stations=chained_paths())
def test_merging_chains_can_only_raise_the_reported_radius_step(stations: List[Tuple[int, float, float, float]]) -> None:
    """One boundary of `_loop_runs`, asserted as an inequality over all paths.

    Relabelling every operation onto a single `path_index` moves no geometry; it
    only removes boundaries, so strictly more loop pairs become comparable and
    the reported maximum can rise or hold but never fall. A rule that splits on
    the wrong thing breaks this in one direction or the other.
    """
    split = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations, retract_between_chains=False))
    merged = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations, retract_between_chains=False, merge_chains=True))
    assert _field(merged, "cut", "max_loop_radius_step") >= _field(split, "cut", "max_loop_radius_step") - REWORK_TOL


@PROPERTY_SETTINGS
@given(stations=chained_paths())
def test_retracting_between_chains_can_only_lower_the_reported_radius_step(stations: List[Tuple[int, float, float, float]]) -> None:
    """The other boundary, in the opposite direction.

    Linking chains at cutting depth and retracting between them differ only in
    whether the tool lifts. Lifting adds a boundary, so it can only remove
    comparable pairs. Together with the merge property this pins the run rule
    from both sides without naming a single expected value.
    """
    linked = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations, retract_between_chains=False))
    lifted = _measure(_spec(1.0, 0, (0.0, 0.0)), _result(stations, retract_between_chains=True))
    assert _field(lifted, "cut", "max_loop_radius_step") <= _field(linked, "cut", "max_loop_radius_step") + REWORK_TOL


@given(rim=st.floats(min_value=0.0, max_value=360.0))
@settings(max_examples=200, deadline=None)
def test_radial_immersion_is_a_monotone_map_onto_the_unit_interval(rim: float) -> None:
    """The angle conversion, quantified rather than sampled at anchors.

    `test_quality.py` pins this function at the angles a text agrees on and
    against the kernel. This adds the shape between them: bounded, and never
    decreasing, so no rim arc anywhere in the range can report LESS immersion
    than a smaller one -- the specific way the superseded reading failed, by
    turning back down past a slot and calling a plunge light.
    """
    value = radial_immersion(rim)
    assert 0.0 <= value <= 1.0
    assert value >= radial_immersion(max(0.0, rim - 1.0)) - REWORK_TOL
