"""The ledger's floor-declined runs are the walked orientation's, not the first offered one's.

The floor is decided per orientation (the authored entry phase differs), so the
two orientations of a chain may decline different runs. The ledger must report
the runs of the orientation the generator actually walked.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable

import pytest

from benchmarks.gate import gate_pocket
from benchmarks.gate import pipeline
from compas_cgal.pipeline import generate
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.result import DeclineCause
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.units import ToolRadius

# Sentinel station indices far past any gate chain, one per walk direction, so
# each orientation's declined runs are distinguishable in the ledger.
SENTINEL_FORWARD = 1_000_000
SENTINEL_REVERSED = 2_000_000
# The L pocket has several chains, so both the first chain and `choose_next` are exercised.
POCKET = "L_shape"
CAP_DEG = 120.0


def _sentinel(chain: OrientedChain) -> tuple[DeclinedRun, ...]:
    station = SENTINEL_REVERSED if chain.reversed_walk else SENTINEL_FORWARD
    return (DeclinedRun(chain=chain.origin, first_station=station, last_station=station),)


def test_the_ledger_reports_the_walked_orientations_declined_runs(monkeypatch: pytest.MonkeyPatch) -> None:
    real_orientations: Callable[[GuideChain, ToolRadius], tuple[OrientedChain, ...]] = generate.orientations
    real_walk_chain = generate.walk_chain
    walked: list[OrientedChain] = []

    def tagged(chain: GuideChain, tool: ToolRadius) -> tuple[OrientedChain, ...]:
        return tuple(dataclasses.replace(oriented, declined=_sentinel(oriented)) for oriented in real_orientations(chain, tool))

    def spied(path: generate.PathState, chain: OrientedChain, start: int, *args: object, **kwargs: object) -> generate.WalkEnd:
        if start == chain.first:
            walked.append(chain)
        return real_walk_chain(path, chain, start, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(generate, "orientations", tagged)
    monkeypatch.setattr(generate, "walk_chain", spied)
    result = pipeline(gate_pocket(POCKET, tea_cap_deg=CAP_DEG))
    assert isinstance(result, PipelineResult)
    assert len(walked) > 1 and any(chain.reversed_walk for chain in walked), "the fixture must walk a reversed orientation"
    floor_runs = tuple(run for run in result.ledger.declined if run.cause is DeclineCause.FLOOR)
    assert floor_runs == tuple(run for chain in walked for run in _sentinel(chain))
