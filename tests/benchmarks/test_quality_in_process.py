"""The cut and material conversions read the gate's IN-PROCESS engagement.

The survey reports the engaged rim arc of each station against the stock minus
the motion's own swept prefix. For a cut into a straight wall that arc runs from
the entry point to the exit point, so it IS the textbook engagement angle; the
static reading, which these conversions were first derived for, was twice it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from compas.geometry import Line
from compas.geometry import Polygon

from benchmarks.models import MaterialModel
from benchmarks.quality import chip_thickness_ratio_from_in_process_rim
from benchmarks.quality import in_process_radial_immersion
from benchmarks.quality import in_process_textbook_engagement_deg
from benchmarks.quality import material_outcome
from benchmarks.spec import PocketSpec
from benchmarks.survey import survey_path
from compas_cgal import _stock_2
from compas_cgal.engagement import _cap_chord_ratio
from compas_cgal.stock import Stock
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult

TOOL_RADIUS = 1.0
# degrees; reporting comparison of the kernel's acos/atan2 span against the closed form.
REPORTED_DEG_TOLERANCE = 1e-6
# mm; the closed-form immersion against a depth the fixture sets exactly.
IMMERSION_DEPTH_TOLERANCE = 1e-9
# degrees; the replay subtracts a line as a disk chain that under-covers the true
# capsule by at most 1e-4 r (`Stock.subtract_capsule`), so a side cut along the
# chain's wall can see an extra sliver of rim of acos(1 - 1e-4) ~ sqrt(2e-4) rad.
DISK_CHAIN_SCALLOP_DEG = math.degrees(math.acos(1.0 - 1e-4))


def test_the_in_process_rim_arc_of_a_side_cut_is_the_textbook_angle() -> None:
    """Pinned against the kernel: a side cut into a straight wall at radial depth a_e.

    The cutter travels along +x below a half-plane. Its in-process engaged arc
    is ``acos(1 - a_e / r)``, which is the textbook entry-to-exit angle, so a slot
    (a_e = D) reads a half turn and half immersion (a_e = r) a quarter turn.
    """
    material = Polygon([[-20.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, 20.0, 0.0], [-20.0, 20.0, 0.0]])
    stock = Stock(material, [])
    ratio = _cap_chord_ratio(math.pi)
    # (radial depth into the material, the in-process rim arc the kernel must report)
    for depth, expected_deg in ((0.0, 0.0), (0.5, 60.0), (1.0, 90.0), (1.5, 120.0), (2.0, 180.0)):
        centre_y = depth - TOOL_RADIUS
        [(_total, run, _exceeded)] = _stock_2.engagement_line_prefixes(stock.raw, -2.0, centre_y, 2.0, centre_y, [0.5], TOOL_RADIUS, ratio, 0.0)
        measured = math.degrees(run)
        assert measured == pytest.approx(expected_deg, abs=REPORTED_DEG_TOLERANCE), f"depth {depth}: kernel reports {measured} deg"
        assert in_process_textbook_engagement_deg(measured) == pytest.approx(expected_deg)
        assert 2.0 * TOOL_RADIUS * in_process_radial_immersion(measured) == pytest.approx(depth, abs=IMMERSION_DEPTH_TOLERANCE)


def test_in_process_conversions_anchor_slot_and_half_immersion() -> None:
    assert in_process_radial_immersion(180.0) == pytest.approx(1.0)
    assert in_process_radial_immersion(90.0) == pytest.approx(0.5)
    assert chip_thickness_ratio_from_in_process_rim(90.0) == pytest.approx(1.0)
    assert chip_thickness_ratio_from_in_process_rim(30.0) == pytest.approx(0.5)
    # A motion's first station is the static reading and can exceed a half turn
    # (a loop entering solid stock); it is a slot, not more than one.
    assert in_process_textbook_engagement_deg(360.0) == pytest.approx(180.0)
    assert in_process_radial_immersion(360.0) == pytest.approx(1.0)


def test_material_outcome_reads_a_half_immersion_side_cut_at_full_chip() -> None:
    """A side cut at a_e = r reads 90 deg in process: the chip is at full thickness.

    Fed through the static conversion the same reading would be taken as a
    quarter of a static rim, and the light end of the chip distribution would
    report sin(45 deg) of the feed.
    """
    spec = PocketSpec.build(
        name="in_process_side_cut",
        family="analytic",
        polygon=Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]]),
        tool_diameter=2.0 * TOOL_RADIUS,
        tea_cap_deg=120.0,
    )
    operations = [
        ToolpathOperation(geometry=Line([-5.0, 0.0, 4.0], [-5.0, 0.0, 0.0]), operation=OperationType.PLUNGE, path_index=0),
        ToolpathOperation(geometry=Line([-5.0, 0.0, 0.0], [5.0, 0.0, 0.0]), operation=OperationType.CUT, path_index=0),
        ToolpathOperation(geometry=Line([-4.0, 1.0, 0.0], [4.0, 1.0, 0.0]), operation=OperationType.CUT, path_index=0),
    ]
    survey = survey_path(spec, ToolpathResult(operations=operations, polyline=np.zeros((0, 3), dtype=float)))
    side_cut = survey.motions[1]
    # The slot runs past both ends of the side cut, so every station after the
    # first is a half-immersion side cut: 90 deg, plus at most the disk-chain sliver.
    after_entry = [sample.engagement_deg for sample in side_cut.samples[1:]]
    assert all(90.0 - REPORTED_DEG_TOLERANCE <= value <= 90.0 + DISK_CHAIN_SCALLOP_DEG for value in after_entry)

    material = MaterialModel.build()
    outcome = material_outcome(survey, material)

    assert outcome.low_chip_thickness_mm == pytest.approx(material.feed_per_tooth_mm)
    assert outcome.max_chip_thickness_mm == pytest.approx(material.feed_per_tooth_mm)
