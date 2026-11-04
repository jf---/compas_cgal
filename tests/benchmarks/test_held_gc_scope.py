"""The scoped cyclic-GC pause: disabled inside, the prior state restored on every exit, nesting-safe."""

from __future__ import annotations

import gc
import weakref
from collections.abc import Iterator

import pytest

from benchmarks import held_omega_toolpath
from benchmarks import held_toolpath_benchmark
from benchmarks.held_gc_scope import paused_cyclic_gc


class _Sentinel(Exception):
    """Raised inside a scope to exercise the exceptional exit."""


@pytest.fixture
def gc_enabled() -> Iterator[None]:
    gc.enable()
    yield
    gc.enable()


@pytest.fixture
def gc_disabled() -> Iterator[None]:
    gc.disable()
    yield
    gc.enable()


@pytest.mark.usefixtures("gc_enabled")
def test_scope_disables_inside_and_reenables_after_normal_exit() -> None:
    with paused_cyclic_gc():
        assert not gc.isenabled()
    assert gc.isenabled()


@pytest.mark.usefixtures("gc_enabled")
def test_scope_reenables_after_an_exception_and_propagates_it() -> None:
    with pytest.raises(_Sentinel), paused_cyclic_gc():
        assert not gc.isenabled()
        raise _Sentinel
    assert gc.isenabled()


@pytest.mark.usefixtures("gc_enabled")
def test_nested_scope_does_not_reenable_early() -> None:
    with paused_cyclic_gc():
        with paused_cyclic_gc():
            assert not gc.isenabled()
        assert not gc.isenabled()
        with pytest.raises(_Sentinel), paused_cyclic_gc():
            raise _Sentinel
        assert not gc.isenabled()
    assert gc.isenabled()


@pytest.mark.usefixtures("gc_disabled")
def test_prior_disabled_state_is_preserved_without_a_collection(monkeypatch: pytest.MonkeyPatch) -> None:
    collections: list[int] = []

    def _collect(*args: int) -> int:
        collections.append(1)
        return 0

    monkeypatch.setattr(gc, "collect", _collect)
    with paused_cyclic_gc():
        assert not gc.isenabled()
    assert not gc.isenabled()
    assert collections == []


@pytest.mark.usefixtures("gc_enabled")
def test_outermost_exit_collects_once_after_reenabling(monkeypatch: pytest.MonkeyPatch) -> None:
    states: list[bool] = []

    def _collect(*args: int) -> int:
        states.append(gc.isenabled())
        return 0

    monkeypatch.setattr(gc, "collect", _collect)
    with paused_cyclic_gc(), paused_cyclic_gc():
        pass
    assert states == [True]


@pytest.mark.usefixtures("gc_enabled")
def test_outermost_exit_collects_the_young_generation_only(monkeypatch: pytest.MonkeyPatch) -> None:
    generations: list[tuple[int, ...]] = []

    def _collect(*args: int) -> int:
        generations.append(args)
        return 0

    monkeypatch.setattr(gc, "collect", _collect)
    with paused_cyclic_gc():
        pass
    assert generations == [(0,)]


class _Node:
    """A cyclic-garbage cell: ``node.other`` closes the cycle."""

    other: _Node | None = None


@pytest.mark.usefixtures("gc_enabled")
def test_a_cycle_created_inside_the_scope_is_reclaimed_on_exit() -> None:
    gc.collect()
    with paused_cyclic_gc():
        first, second = _Node(), _Node()
        first.other, second.other = second, first
        alive = weakref.ref(first)
        del first, second
        assert alive() is not None
    assert alive() is None


@pytest.mark.usefixtures("gc_enabled")
def test_generate_omega_runs_with_the_collector_paused(monkeypatch: pytest.MonkeyPatch) -> None:
    observed: list[bool] = []

    def _load(name: str) -> None:
        observed.append(gc.isenabled())
        raise _Sentinel

    monkeypatch.setattr(held_omega_toolpath, "load_held_reference_case", _load)
    with pytest.raises(_Sentinel):
        held_omega_toolpath.generate_omega("figure5_refined")
    assert observed == [False]
    assert gc.isenabled()


@pytest.mark.usefixtures("gc_enabled")
@pytest.mark.parametrize("workload", ["figure5_refined", "figure8_monstera_boundary_qualified_corners"])
def test_polygon_generators_run_with_the_collector_paused(monkeypatch: pytest.MonkeyPatch, workload: str) -> None:
    observed: list[bool] = []

    def _load(name: str) -> None:
        observed.append(gc.isenabled())
        raise _Sentinel

    monkeypatch.setattr(held_toolpath_benchmark, "load_held_reference_case", _load)
    with pytest.raises(_Sentinel):
        held_toolpath_benchmark.generate(workload)
    assert observed == [False]
    assert gc.isenabled()
