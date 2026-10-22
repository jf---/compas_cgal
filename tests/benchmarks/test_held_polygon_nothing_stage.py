"""The exact NOTHING-circle pass runs inside polygon generation as its own timed stage, closed over the seam."""

from __future__ import annotations

import math
from functools import cache

import pytest

from benchmarks import held_toolpath_benchmark as tb
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_nothing_circles import NothingCirclePass
from compas_cgal.motion import EngagementCap

pytestmark = pytest.mark.xdist_group("held_polygon_figure5_nothing")

# Figure 5 polygon circles the pass proves NOTHING (indices into the refined
# circles; NR's acceptance table, nr-nothing-circles-report.md, equal to the
# full exact replay's NOTHING set).
FIGURE5_POLYGON_NOTHING = (980, 987, 995, 1017, 1036, 1045)
# Index of the emitted circle duplicated in the defect-injection test: any
# interior index works, the duplicate's successor pair is an emitted pair.
DUPLICATED = 10


@cache
def _figure5() -> tuple[tb.GeneratedToolpath, tuple[tuple[Figure5CounterclockwiseCircle, ...], NothingCirclePass[Figure5CounterclockwiseCircle]]]:
    calls: list[tuple[tuple[Figure5CounterclockwiseCircle, ...], NothingCirclePass[Figure5CounterclockwiseCircle]]] = []
    stage = tb.remove_polygon_nothing_circles

    def recording(circles: tuple[Figure5CounterclockwiseCircle, ...], *args: object, **kwargs: object) -> NothingCirclePass[Figure5CounterclockwiseCircle]:
        result = stage(circles, *args, **kwargs)  # type: ignore[arg-type]
        calls.append((tuple(circles), result))
        return result

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(tb, "remove_polygon_nothing_circles", recording)
        path = tb.generate("figure5_refined")
    assert len(calls) == 1, "generation runs the pass exactly once"
    return path, calls[0]


def _cap(path: tb.GeneratedToolpath) -> EngagementCap:
    return EngagementCap.build(math.radians(float(path.case.tea_cap)))


def test_generation_times_the_nothing_stage_between_refine_and_connectors() -> None:
    stages = list(_figure5()[0].stage_seconds)
    assert stages.index("refine") < stages.index("nothing") < stages.index("connectors") < stages.index("connector_certificate")


def test_figure5_emits_the_refined_circles_minus_the_proven_nothing_set() -> None:
    path, (refined, nothing) = _figure5()
    assert nothing.removed == FIGURE5_POLYGON_NOTHING
    assert nothing.kept_for_successor_cap == ()
    removed = set(nothing.removed)
    assert path.circles == nothing.circles == tuple(c for k, c in enumerate(refined) if k not in removed)
    assert path.connectors.seam is Figure5SeamClosure.CONNECTOR
    assert len(path.transitions) == len(path.circles)


def test_the_emitted_path_has_no_nothing_circle_left() -> None:
    """The pass is idempotent: each kept circle was not covered by the union of all earlier input circles, which equals the kept union."""
    path, _ = _figure5()
    again = tb.remove_polygon_nothing_circles(path.circles, path.case.tool_radius, _cap(path))
    assert again.removed == ()
    assert again.circles == path.circles


def test_an_exact_duplicate_of_an_emitted_circle_is_removed() -> None:
    path, _ = _figure5()
    circles = path.circles
    injected = (*circles[: DUPLICATED + 1], circles[DUPLICATED], *circles[DUPLICATED + 1 :])
    result = tb.remove_polygon_nothing_circles(injected, path.case.tool_radius, _cap(path))
    assert result.removed == (DUPLICATED + 1,)
    assert result.circles == circles


def test_the_emitted_connectors_are_rebuilt_on_the_pruned_circles() -> None:
    path, _ = _figure5()
    assert path.connectors == Figure5BoundaryConnectors.build(path.boundary, path.circles)
