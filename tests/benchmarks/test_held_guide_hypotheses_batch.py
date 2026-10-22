"""The whole-guide hypothesis build: the per-station hypotheses, without keeping the station analyses alive.

``build_guide_distance_admissible_hypotheses`` must equal the concatenation of
``build_distance_admissible_hypotheses`` over every run and station, run the
native batch once, and leave the guide's memo empty when it was empty, so the
analyses (about 95 MB on Monstera) are released once the hypotheses exist.
"""

from __future__ import annotations

import pytest

import benchmarks.held_figure5_raw_guide as raw_guide_module
from benchmarks.held_figure5_raw_guide import Figure5RawGuide
from benchmarks.held_figure5_raw_guide import build_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import build_guide_distance_admissible_hypotheses
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import load_held_reference_case


def _guide(case_name: str) -> Figure5RawGuide:
    case = load_held_reference_case(case_name)
    return build_figure5_raw_guide(case) if case_name == "figure5" else build_held_reference_raw_guide(case)


def _counted_batch(monkeypatch: pytest.MonkeyPatch) -> list[Figure5RawGuide]:
    calls: list[Figure5RawGuide] = []
    wired = raw_guide_module._native_station_boundary_analyses

    def counted(guide: Figure5RawGuide) -> object:
        calls.append(guide)
        return wired(guide)

    monkeypatch.setattr(raw_guide_module, "_native_station_boundary_analyses", counted)
    return calls


@pytest.mark.parametrize("case_name", ["figure5", "figure8_crossed_skis"])
def test_the_whole_guide_build_equals_the_per_station_builds(case_name: str) -> None:
    batch = build_guide_distance_admissible_hypotheses(_guide(case_name))
    reference_guide = _guide(case_name)
    per_station = tuple(hypothesis for run in reference_guide.runs for station in run.stations for hypothesis in build_distance_admissible_hypotheses(reference_guide, station))
    assert batch == per_station
    assert repr(batch) == repr(per_station)


def test_an_empty_memo_stays_empty_and_the_batch_runs_once(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _counted_batch(monkeypatch)
    guide = _guide("figure5")
    assert build_guide_distance_admissible_hypotheses(guide)
    assert calls == [guide]
    assert guide.derived.station_boundary_analyses is None


def test_a_filled_memo_is_read_and_kept(monkeypatch: pytest.MonkeyPatch) -> None:
    guide = _guide("figure5")
    analyses = guide.station_boundary_analyses
    calls = _counted_batch(monkeypatch)
    assert build_guide_distance_admissible_hypotheses(guide)
    assert calls == []
    assert guide.derived.station_boundary_analyses is analyses
