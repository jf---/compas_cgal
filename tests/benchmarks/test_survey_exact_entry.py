"""The survey reads a loop's exact entry when the operation carries one (Ruling 1.1-C14-A).

An `IncidentLoop`'s ``Circle.point_at(0)`` is a double rounding of its exact
entry, so a plunge onto the exact entry would be misplaced under rule (e) if the
survey read only ``point_at(0)``. A circle without an entry is read exactly as
before.
"""

from __future__ import annotations

import dataclasses

import pytest
from compas.geometry import Circle
from compas.geometry import Frame

from benchmarks.quality import measure_quality
from benchmarks.survey import survey_path
from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.emit import loop_operation
from compas_cgal.units import CutZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY
from tests.benchmarks.test_quality import FAST_SAMPLES
from tests.benchmarks.test_quality import SYNTHETIC
from tests.benchmarks.test_quality import SYNTHETIC_GRID
from tests.benchmarks.test_quality import _circle
from tests.benchmarks.test_quality import _plunge
from tests.benchmarks.test_quality import _result

# Phases at which the authored entry and point_at(0) differ (P2: not ulp corners).
OFF_AXIS_PHASES = (0.3, 1.1, 3.7, 4.4)


def _loop(phase: float) -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(2.0), Radian(phase))


@pytest.mark.parametrize("phase", OFF_AXIS_PHASES)
def test_a_plunge_onto_an_incident_loops_exact_entry_is_not_misplaced(phase: float) -> None:
    loop = _loop(phase)
    operation = loop_operation(loop, CutZ.build(0.0), 0, clockwise=True)
    rounded = operation.geometry.point_at(0.0)
    assert (float(rounded[0]), float(rounded[1])) != (loop.entry.x, loop.entry.y), "fixture precondition (C14): point_at(0) is not the entry"
    quality = measure_quality(SYNTHETIC, _result([_plunge(loop.entry.x, loop.entry.y), operation]), samples_per_motion=FAST_SAMPLES, grid=SYNTHETIC_GRID)
    assert quality.elementary.misplaced_plunges == 0


@pytest.mark.parametrize("phase", OFF_AXIS_PHASES)
def test_the_loop_starts_ends_and_seams_at_its_exact_entry(phase: float) -> None:
    loop = _loop(phase)
    result = _result([_plunge(loop.entry.x, loop.entry.y), loop_operation(loop, CutZ.build(0.0), 0, clockwise=True)])
    (motion,) = survey_path(SYNTHETIC, result, samples_per_motion=FAST_SAMPLES).motions
    entry = (loop.entry.x, loop.entry.y)
    assert motion.start == entry and motion.end == entry
    assert (motion.samples[0].position.x, motion.samples[0].position.y) == entry
    assert motion.exit_sample is not None and (motion.exit_sample.position.x, motion.exit_sample.position.y) == entry


def test_a_frozen_circle_without_an_entry_reads_point_at_0_unchanged() -> None:
    circle = _circle(0.0, 0.0, 2.0)
    assert circle.entry_xy is None
    tilted = dataclasses.replace(circle, geometry=Circle(2.0, frame=Frame([0.0, 0.0, 0.0], [0.6, 0.8, 0.0], [-0.8, 0.6, 0.0])))
    start = tilted.geometry.point_at(0.0)
    (motion,) = survey_path(SYNTHETIC, _result([_plunge(float(start[0]), float(start[1])), tilted]), samples_per_motion=FAST_SAMPLES).motions
    point = (float(start[0]), float(start[1]))
    assert motion.start == point and motion.end == point
    assert (motion.samples[0].position.x, motion.samples[0].position.y) == point
