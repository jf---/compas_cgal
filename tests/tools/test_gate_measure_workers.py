"""The gate measurement's process pool is bounded by the core count (Ruling B-I1)."""

from __future__ import annotations

from typing import Callable
from typing import Iterable
from typing import Iterator

import pytest

import tools.gate_measure as gate_measure
from tools.gate_measure import GateCase
from tools.gate_measure import UnknownCoreCountError
from tools.gate_measure import gate_cases
from tools.gate_measure import measure_gate_vectors

SMALL_CORE_COUNT = 2
CASES_ASKED = 5


class RecordingPool:
    """A stand-in for `ProcessPoolExecutor` that records its worker count and maps serially."""

    asked: list[int] = []

    def __init__(self, max_workers: int) -> None:
        RecordingPool.asked.append(max_workers)

    def __enter__(self) -> RecordingPool:
        return self

    def __exit__(self, *_exc: object) -> None:
        return None

    def map(self, fn: Callable[[object], object], items: Iterable[object]) -> Iterator[object]:
        return iter([fn(item) for item in items])


def _fake_measure(case: GateCase) -> tuple[str, dict[str, str]]:
    cap_id, _cap, generator, pocket = case
    return gate_measure.case_key(cap_id, generator, pocket), {}


@pytest.fixture
def recording_pool(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    RecordingPool.asked = []
    monkeypatch.setattr(gate_measure, "ProcessPoolExecutor", RecordingPool)
    monkeypatch.setattr(gate_measure, "measure_case", _fake_measure)
    return RecordingPool.asked


def test_the_worker_count_never_exceeds_the_core_count(recording_pool: list[int], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gate_measure.os, "cpu_count", lambda: SMALL_CORE_COUNT)
    measure_gate_vectors(gate_cases()[:CASES_ASKED])
    assert recording_pool == [SMALL_CORE_COUNT]


def test_fewer_cases_than_cores_ask_only_for_the_cases(recording_pool: list[int], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gate_measure.os, "cpu_count", lambda: CASES_ASKED + 1)
    measure_gate_vectors(gate_cases()[:CASES_ASKED])
    assert recording_pool == [CASES_ASKED]


def test_an_unknown_core_count_is_refused(recording_pool: list[int], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(gate_measure.os, "cpu_count", lambda: None)
    with pytest.raises(UnknownCoreCountError):
        measure_gate_vectors(gate_cases()[:CASES_ASKED])
    assert recording_pool == []
