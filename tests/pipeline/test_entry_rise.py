"""The declared entry ramp and its exact junction surrogate (Task 1.2, USER 2026-10-16 Q4).

`EntryRampLength` is declared user intent (default one tool diameter). The
rise limit at the gate's station 1 is ``theta_e = cap * s_e / L`` with
``s_e = 2 pi r / N`` (the station-1 arc on the smallest legal loop), i.e.
``cap * pi / 45`` at the default; its chord surrogate is computed once as a
double and injected exactly in `Admissibility.build`, beside the cap's. The
junction decision is the native exact predicate at that one station.
"""

from __future__ import annotations

import math

import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.helix_units import EntryRampLength
from compas_cgal.helix_units import InvalidEntryRampLengthError
from compas_cgal.incident_loop import IncidentLoop
from benchmarks.quality_observations import JUNCTION_STATION as GATE_JUNCTION_STATION
from benchmarks.survey import QUALITY_SAMPLES_PER_MOTION
from compas_cgal.pipeline.admissibility import GAP_CLOSE_RATIO
from compas_cgal.pipeline.admissibility import GATE_SAMPLES_PER_MOTION
from compas_cgal.pipeline.admissibility import JUNCTION_STATION
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import EntryRiseChordRatio
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.admissibility import centre_domain
from compas_cgal.pipeline.admissibility import entry_rise_ratio
from compas_cgal.pipeline.errors import InvalidEntryRiseChordRatioError
from compas_cgal.pipeline.errors import MissingEntryRiseLimitError
from compas_cgal.stock import Stock
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])
TOOL = ToolRadius.build(1.0)
GATE_CAPS_DEG = (120.0, 40.0)


def test_the_junction_station_is_the_gates_own() -> None:
    """The restated constant equals the gate's junction station (pipeline/ imports no benchmarks module)."""
    assert JUNCTION_STATION == GATE_JUNCTION_STATION


def test_a_ramp_shorter_than_the_station_one_arc_is_refused_by_name() -> None:
    """theta_e = cap * 2 pi / (45 * L / r) > pi when L / r < 2 cap / 45: no chord surrogate exists."""
    with pytest.raises(InvalidEntryRiseChordRatioError):
        entry_rise_ratio(Radian(math.pi), EntryRampLength.build(0.1))


@pytest.mark.parametrize("value", [0.0, -1.0, math.nan, math.inf, True])
def test_an_entry_ramp_length_that_is_not_finite_positive_is_refused_by_name(value: float) -> None:
    with pytest.raises(InvalidEntryRampLengthError):
        EntryRampLength.build(value)


@pytest.mark.parametrize("cap_deg", GATE_CAPS_DEG)
def test_the_default_rise_surrogate_is_the_chord_ratio_of_cap_pi_over_45(cap_deg: float) -> None:
    cap = math.radians(cap_deg)
    assert entry_rise_ratio(Radian(cap), EntryRampLength.default()) == _stock_2.cap_chord_ratio(cap * math.pi / GATE_SAMPLES_PER_MOTION)


def test_admissibility_injects_the_rise_surrogate_beside_the_cap() -> None:
    ratio = entry_rise_ratio(Radian(math.radians(120.0)), EntryRampLength.default())
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(math.radians(120.0))), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=ratio)
    assert admissibility.entry_rise_ratio == ratio


@pytest.mark.parametrize("ratio", [0.0, -0.5, 4.5, math.nan])
def test_a_rise_surrogate_no_angle_in_0_pi_has_is_refused_by_name(ratio: float) -> None:
    with pytest.raises(InvalidEntryRiseChordRatioError):
        Admissibility.build(TOOL, cap_ratio(Radian(1.0)), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=EntryRiseChordRatio(ratio))


def _loop() -> IncidentLoop:
    return IncidentLoop.build(Point2[WorldXY].build(0.0, 0.0), GuideRadius.build(2.0), Radian(0.0))


def test_the_junction_without_a_declared_rise_limit_raises_by_name() -> None:
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(1.0)), centre_domain(RECT_12X8, [], TOOL))
    with pytest.raises(MissingEntryRiseLimitError):
        admissibility.loop_entry_within_rise(Stock(RECT_12X8), _loop(), True)


def test_the_junction_asks_only_the_gate_station_one_at_the_rise_surrogate(monkeypatch: pytest.MonkeyPatch) -> None:
    ratio = entry_rise_ratio(Radian(math.radians(120.0)), EntryRampLength.default())
    admissibility = Admissibility.build(TOOL, cap_ratio(Radian(math.radians(120.0))), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=ratio)
    asked: list[tuple[list[float], float]] = []
    native = _stock_2.engagement_circle_exceeded

    def recorded(*args: object) -> bool:
        asked.append((list(args[7]), float(args[9])))  # type: ignore[call-overload, arg-type]
        return bool(native(*args))  # type: ignore[arg-type]

    monkeypatch.setattr(_stock_2, "engagement_circle_exceeded", recorded)
    admissibility.loop_entry_within_rise(Stock(RECT_12X8), _loop(), True)
    assert asked == [([1 / GATE_SAMPLES_PER_MOTION], ratio)]


# Candidate centres beside a pre-cut loop at the origin: some junction stations fall in air, some in material.
CANDIDATE_CENTRES_X = (-3.0, -2.0, -1.0, 0.0, 1.0, 2.5)


@pytest.mark.parametrize("cap_deg", GATE_CAPS_DEG)
def test_the_junction_verdict_is_the_gate_samples_exceeded_flag_at_its_junction_station(cap_deg: float) -> None:
    """The decision equals the reporting query the survey samples, at the gate's own station ``JUNCTION_STATION / QUALITY_SAMPLES_PER_MOTION``; both verdicts occur."""
    cap = Radian(math.radians(cap_deg))
    ratio = entry_rise_ratio(cap, EntryRampLength.default())
    admissibility = Admissibility.build(TOOL, cap_ratio(cap), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=ratio)
    stock = Stock(RECT_12X8)
    stock.subtract_circle_sweep_local(0.0, 0.0, 1.5, TOOL.value)
    verdicts = set()
    for centre_x in CANDIDATE_CENTRES_X:
        loop = IncidentLoop.build(Point2[WorldXY].build(centre_x, 0.0), GuideRadius.build(1.5), Radian(0.0))
        ((_station, _run, exceeded),) = _stock_2.engagement_circle_prefixes(
            stock.raw,
            loop.centre.x,
            loop.centre.y,
            loop.radius.value,
            loop.entry.x,
            loop.entry.y,
            True,
            [GATE_JUNCTION_STATION / QUALITY_SAMPLES_PER_MOTION],
            TOOL.value,
            ratio,
            GAP_CLOSE_RATIO,
        )
        verdict = admissibility.loop_entry_within_rise(stock, loop, True)
        assert verdict is (not exceeded)
        verdicts.add(verdict)
    assert verdicts == {True, False}


def test_the_junction_refuses_full_immersion_at_the_rise_limit_and_admits_it_at_a_half_turn() -> None:
    """Virgin stock: the loop's junction station is deep in material, so the 8.4 deg rise limit refuses and a half-turn limit admits."""
    loop, stock = _loop(), Stock(RECT_12X8)
    tight = Admissibility.build(
        TOOL, cap_ratio(Radian(math.pi)), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=entry_rise_ratio(Radian(math.radians(120.0)), EntryRampLength.default())
    )
    loose = Admissibility.build(TOOL, cap_ratio(Radian(math.pi)), centre_domain(RECT_12X8, [], TOOL), entry_rise_ratio=EntryRiseChordRatio(_stock_2.cap_chord_ratio(math.pi)))
    assert not tight.loop_entry_within_rise(stock, loop, True)
    assert loose.loop_entry_within_rise(stock, loop, True)
