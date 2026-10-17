"""The survey observes a HELIX entry: ends, tangents, length, its measured exit, no cap sampling (Task 1.2).

A helix is not cap-sampled (E8): a bore into virgin stock is full immersion by
construction, scored by its ramp. So it is never a `MotionQuality`; it is a
`HelixObservation`, and the motion after it is sampled as usual. Its exit is
MEASURED (amendment M6): the footprint circle's seam query with the full turn
as prefix, against the stock the helix started from, so the helix -> first loop
junction is measurable rather than assumed. HELIX operations come from the real
emitter (`compas_cgal.pipeline.emit.helix_operation`).
"""

from __future__ import annotations

import dataclasses
import math

import pytest
from compas.geometry import Line
from compas.tolerance import TOL

from benchmarks.depletion import CutMotion
from benchmarks.errors import UnmeasurableOperationLengthError
from benchmarks.errors import UnreplayableOperationError
from benchmarks.errors import UnsampleableMotionError
from benchmarks.helix_observation import HelixObservation
from benchmarks.survey import HelixHasNoRapidVerdictError
from benchmarks.survey import PathSurvey
from benchmarks.survey import UnknownToolpathGeometryError
from benchmarks.survey import _endpoints
from benchmarks.survey import _is_horizontal_at
from benchmarks.survey import _measure_motion
from benchmarks.survey import _primitive_length
from benchmarks.survey import _tangents
from benchmarks.survey import survey_path
from compas_cgal.helix import Helix
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from tests.benchmarks.helix_fixtures import CUT_Z
from tests.benchmarks.helix_fixtures import FOOTPRINT_RADIUS
from tests.benchmarks.helix_fixtures import POCKET
from tests.benchmarks.helix_fixtures import TOOL_RADIUS
from tests.benchmarks.helix_fixtures import TOP_Z
from tests.benchmarks.helix_fixtures import emitted
from tests.benchmarks.helix_fixtures import first_member
from tests.benchmarks.helix_fixtures import make_helix
from tests.benchmarks.helix_fixtures import result

SAMPLES = 8
# A footprint wider than the tool: the bore is an annulus, not a disk (I1).
WIDE_FOOTPRINT_RADIUS = 1.5 * TOOL_RADIUS


def _survey_entry(clockwise: bool = True) -> tuple[Helix, PathSurvey]:
    helix = make_helix(clockwise)
    return helix, survey_path(POCKET, result([emitted(helix), first_member(helix)]), samples_per_motion=SAMPLES)


@pytest.mark.parametrize("clockwise", [True, False])
def test_a_helix_starts_and_ends_at_its_exact_entry_with_the_circle_tangent(clockwise: bool) -> None:
    helix = make_helix(clockwise)
    operation = emitted(helix)
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert _endpoints(operation) == (entry, entry)
    # P = (rho_h, 0) about the origin: travel is -y clockwise, +y counter-clockwise.
    expected = (0.0, -1.0) if clockwise else (0.0, 1.0)
    assert _tangents(operation) == (expected, expected)


def test_a_helix_length_is_the_helix_length() -> None:
    helix = make_helix()
    assert _primitive_length(0, emitted(helix)) == float(helix.length())


def test_the_helix_is_excluded_from_cap_sampling_like_a_plunge() -> None:
    _helix, survey = _survey_entry()
    assert [motion.index for motion in survey.motions] == [1]
    assert len(survey.motions[0].samples) == SAMPLES
    (observed,) = survey.helices
    assert isinstance(observed, HelixObservation)
    assert observed.index == 0
    assert observed.into_material
    assert survey.rapids == () and survey.plunge_indices == () and survey.material_plunge_indices == ()


def test_the_helix_observation_carries_its_ends_descent_and_length() -> None:
    helix, survey = _survey_entry()
    (observed,) = survey.helices
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert observed.start == entry and observed.end == entry
    assert (observed.start_tangent, observed.end_tangent) == ((0.0, -1.0), (0.0, -1.0))
    assert (observed.top_z, observed.bottom_z) == (TOP_Z, CUT_Z)
    assert observed.turns == helix.turns and observed.pitch == helix.pitch
    assert observed.length == float(helix.length())
    assert survey.total_length == float(helix.length()) + survey.motions[0].length
    assert TOL.is_close(observed.swept_area, math.pi * (FOOTPRINT_RADIUS + TOOL_RADIUS) ** 2)
    assert TOL.is_close(survey.swept_area, survey.motions[0].swept_area + observed.swept_area)


def test_a_wide_footprint_helix_sweeps_an_annulus_like_its_replay() -> None:
    helix = make_helix(radius=WIDE_FOOTPRINT_RADIUS)
    survey = survey_path(POCKET, result([emitted(helix)]), samples_per_motion=SAMPLES)
    (observed,) = survey.helices
    assert TOL.is_close(observed.swept_area, 4.0 * math.pi * WIDE_FOOTPRINT_RADIUS * TOOL_RADIUS)


def test_the_helix_exit_is_measured_on_the_post_bore_stock() -> None:
    helix, survey = _survey_entry()
    (observed,) = survey.helices
    exit_sample = observed.exit_sample
    entry = (helix.entry_xy.x, helix.entry_xy.y)
    assert (exit_sample.position.x, exit_sample.position.y) == entry
    assert exit_sample.distance == float(helix.length())
    # Measured, not hard-coded: the bored disk (rho_h + r about c_h) holds the
    # tool at P, internally tangent, so no rim run is engaged.
    assert exit_sample.engagement_deg == 0.0
    assert not exit_sample.cap_exceeded
    assert exit_sample.inside_centre_domain


def test_the_helix_to_first_loop_junction_is_measurable() -> None:
    _helix, survey = _survey_entry()
    (observed,) = survey.helices
    (member,) = survey.motions
    assert observed.end == member.start
    step = abs(member.samples[1].engagement_deg - observed.exit_sample.engagement_deg)
    assert step > 0.0, "the first member cuts past the bore, so the junction step is a real measurement"


def test_a_helix_off_the_cut_plane_is_unreplayable() -> None:
    helix = make_helix()
    below = ToolpathOperation(geometry=Line([5.0, 3.0, CUT_Z - 1.0], [5.0, 2.0, CUT_Z - 1.0]), operation=OperationType.CUT, path_index=0)
    with pytest.raises(UnreplayableOperationError):
        survey_path(POCKET, result([emitted(helix), below]), samples_per_motion=SAMPLES)


def test_a_helix_is_never_measured_as_a_cut_motion() -> None:
    motion = CutMotion(index=0, operation=emitted(make_helix()), stock=Stock(POCKET.polygon, []))
    with pytest.raises(UnsampleableMotionError):
        _measure_motion(motion, TOOL_RADIUS, 1.0, None, SAMPLES)  # type: ignore[arg-type]


def test_a_helix_has_no_rapid_verdict_by_name() -> None:
    with pytest.raises(HelixHasNoRapidVerdictError):
        _is_horizontal_at(emitted(make_helix()), CUT_Z)


def test_geometry_outside_the_union_is_refused_by_name() -> None:
    stranger = dataclasses.replace(emitted(make_helix()), geometry=object())
    with pytest.raises(UnknownToolpathGeometryError):
        _endpoints(stranger)
    with pytest.raises(UnknownToolpathGeometryError):
        _tangents(stranger)
    with pytest.raises(UnknownToolpathGeometryError):
        _is_horizontal_at(stranger, CUT_Z)
    with pytest.raises(UnmeasurableOperationLengthError):
        _primitive_length(0, stranger)
