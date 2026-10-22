"""Ω* station ownership bridges junction arcs, as the polygon ring rule does.

A junction arc is the offset of a G1 junction of the declared boundary: a
near-zero piece that hypotheses almost never decide onto. A station whose
hypotheses straddle one lands on the pieces either side of it and must still be
one neighbourhood with one owner. Pieces separated by a non-junction piece (a
neck) stay separate neighbourhoods.
"""

from __future__ import annotations

from types import SimpleNamespace

from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_omega_path import _owner_indices
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

LINE = PieceKind.OFFSET_LINE
ARC = PieceKind.OFFSET_ARC
JUNCTION = PieceKind.JUNCTION_ARC


def _hypothesis(station: int, depth: float, x: float) -> SimpleNamespace:
    """A station hypothesis whose middle point lies ``depth`` from its footpoint."""
    return SimpleNamespace(
        run_id=GuideRunId(0),
        station_ordinal=GuideRunStationOrdinal(station),
        middle_point=Point2[WorldXY].build(x, depth),
        boundary_footpoint=Point2[WorldXY].build(x, 0.0),
        contact_point=Point2[WorldXY].build(x, 0.0),
    )


def test_a_station_straddling_a_junction_arc_has_one_owner() -> None:
    near, far = _hypothesis(0, 1.0, 0.0), _hypothesis(0, 2.0, 1.0)
    owners = _owner_indices((near, far), (0, 2), (LINE, JUNCTION, ARC, LINE), near.contact_point)  # type: ignore[arg-type]
    assert owners == (0,)


def test_a_junction_arc_before_the_first_piece_bridges_cyclically() -> None:
    near, far = _hypothesis(0, 1.0, 0.0), _hypothesis(0, 2.0, 1.0)
    owners = _owner_indices((near, far), (3, 1), (JUNCTION, LINE, ARC, LINE), near.contact_point)  # type: ignore[arg-type]
    assert owners == (0,)


def test_a_hypothesis_on_the_junction_arc_joins_its_neighbours() -> None:
    near, on_junction, far = _hypothesis(0, 1.0, 0.0), _hypothesis(0, 1.5, 0.5), _hypothesis(0, 2.0, 1.0)
    owners = _owner_indices((near, on_junction, far), (0, 1, 2), (LINE, JUNCTION, ARC, LINE), near.contact_point)  # type: ignore[arg-type]
    assert owners == (0,)


def test_pieces_separated_by_a_non_junction_piece_stay_two_neighbourhoods() -> None:
    near, far = _hypothesis(0, 1.0, 0.0), _hypothesis(0, 2.0, 1.0)
    owners = _owner_indices((near, far), (0, 2), (LINE, ARC, LINE, LINE), near.contact_point)  # type: ignore[arg-type]
    assert owners == (0, 1)
