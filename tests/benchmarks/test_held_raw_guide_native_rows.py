"""The raw guide read straight from the native rows equals the compas-operation read (lever 2, step 2)."""

from __future__ import annotations

import pytest

from benchmarks.held_figure5_raw_guide import FIGURE5_MAX_GUIDE_RUNS
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import Figure5RawGuideRun
from benchmarks.held_figure5_raw_guide import Figure5RawGuideStation
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import InvalidFigure5RawGuideError
from benchmarks.held_figure5_raw_guide import _guide_station_centres
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.engagement_toolpath import GUIDE_STEP_TOOL_DIAMETERS
from compas_cgal.engagement_toolpath import _guide_chains
from compas_cgal.engagement_toolpath import _GuideStation
from compas_cgal.toolpath import RADIAL_CLEARANCE_FRACTION
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# The three generation workloads' pockets (figure5_refined, crossed skis, Monstera).
WORKLOAD_CASES = ("figure5", "figure8_crossed_skis", "figure8_monstera")


def _oracle_source_runs(case: HeldReferenceCase) -> list[list[_GuideStation]]:
    """Verbatim pre-step-2 read: the guide through compas operations."""
    spec = case.pocket_spec()
    return _guide_chains(
        spec.polygon,
        spec.tool_diameter,
        GUIDE_STEP_TOOL_DIAMETERS * spec.tool_diameter,
        RADIAL_CLEARANCE_FRACTION * spec.tool_diameter,
        False,
        FIGURE5_MAX_GUIDE_RUNS,
        list(spec.holes),
    )


def _oracle_station(run_id: int, ordinal_on_run: int, station: _GuideStation) -> Figure5RawGuideStation:
    """Verbatim pre-step-2 ``Figure5RawGuideStation.build``."""
    if run_id < 0 or ordinal_on_run < 0:
        raise InvalidFigure5RawGuideError("Figure 5 guide identities must be non-negative.")
    return Figure5RawGuideStation(
        GuideRunId(run_id),
        GuideRunStationOrdinal(ordinal_on_run),
        Point2[WorldXY].build(station.cx, station.cy),
    )


def _oracle_guide(case: HeldReferenceCase) -> Figure5RawGuide:
    """Verbatim pre-step-2 ``build_held_reference_raw_guide`` body."""
    runs = tuple(
        Figure5RawGuideRun.build(
            run_index,
            tuple(_oracle_station(run_index, ordinal, station) for ordinal, station in enumerate(source_run)),
        )
        for run_index, source_run in enumerate(_oracle_source_runs(case))
    )
    return Figure5RawGuide.build_reference(case, runs)


@pytest.mark.parametrize("case_name", WORKLOAD_CASES)
def test_native_station_centres_equal_operation_station_centres(case_name: str) -> None:
    case = load_held_reference_case(case_name)

    native = _guide_station_centres(case)

    assert native == [[(station.cx, station.cy) for station in run] for run in _oracle_source_runs(case)]


@pytest.mark.parametrize("case_name", WORKLOAD_CASES)
def test_raw_guide_is_byte_identical_to_operation_read(case_name: str) -> None:
    case = load_held_reference_case(case_name)

    guide = build_held_reference_raw_guide(case)
    oracle = _oracle_guide(case)

    assert repr(guide.runs) == repr(oracle.runs)
    assert guide == oracle


def test_station_rejects_negative_identity() -> None:
    with pytest.raises(InvalidFigure5RawGuideError, match="non-negative"):
        Figure5RawGuideStation.build(-1, 0, Point2[WorldXY].build(0.0, 0.0))
