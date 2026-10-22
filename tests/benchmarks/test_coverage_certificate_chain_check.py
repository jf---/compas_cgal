"""The ∂Ω* chain check's remaining failure branches (exact_connectors_2.cpp, ``verify_boundary_chains``).

- **A piece ends where two pieces start.** This is also how a pinch point of Ω*
  fails: two boundary cycles through one vertex make that vertex the start of
  two pieces. Loud and correct: the cap-circle theorem assumes simple chains.
- **A junction reverses the travel.** Two collinear pieces travelled in
  opposite directions: the turn is zero and the tangents point apart.
"""

from __future__ import annotations

import dataclasses
from fractions import Fraction

import pytest

from compas_cgal.coverage_certificate.errors import BoundaryChainError
from compas_cgal.coverage_certificate.exact_target import OffsetLinePiece
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.units import ExactMillimetre
from tests.benchmarks.test_coverage_certificate_exact_connectors import INTERIOR
from tests.benchmarks.test_coverage_certificate_exact_connectors import TOOL
from tests.benchmarks.test_coverage_certificate_exact_connectors import UNUSED_PATH
from tests.benchmarks.test_coverage_certificate_exact_connectors import _rational_point
from tests.benchmarks.test_coverage_certificate_exact_connectors import _sweeps
from tests.benchmarks.test_coverage_certificate_exact_connectors import _whole_boundary
from tests.benchmarks.test_coverage_certificate_exact_connectors import square

F = Fraction


def _run(target: object) -> None:
    run_certificate(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=target, exact_connectors=_whole_boundary(target))  # type: ignore[arg-type]


def test_a_piece_ending_where_two_pieces_start_is_a_named_error() -> None:
    _, _, target = square()
    branched = dataclasses.replace(target, pieces=(*target.pieces, target.pieces[1]))
    with pytest.raises(BoundaryChainError, match="ends where two pieces start"):
        _run(branched)


def test_a_junction_that_reverses_the_travel_is_a_named_error() -> None:
    # The square's bottom piece y = 1 from (1, 1) to (19, 1) replaced by
    # (1, 1) -> (10, 1) -> (5, 1) -> (19, 1): a closed chain whose first junction
    # turns back on itself.
    _, _, target = square()
    index, bottom = next((i, p) for i, p in enumerate(target.pieces) if isinstance(p, OffsetLinePiece) and p.normal == (F(0), F(1)))
    ten, five = _rational_point(F(10), F(1)), _rational_point(F(5), F(1))
    back = OffsetLinePiece.build((F(0), F(-1)), ExactMillimetre(F(-1)), ten, five)
    pieces = (
        *target.pieces[:index],
        OffsetLinePiece.build(bottom.normal, bottom.level, bottom.source, ten),
        back,
        OffsetLinePiece.build(bottom.normal, bottom.level, five, bottom.target),
        *target.pieces[index + 1 :],
    )
    with pytest.raises(BoundaryChainError, match="reverses the travel"):
        _run(dataclasses.replace(target, pieces=pieces))
