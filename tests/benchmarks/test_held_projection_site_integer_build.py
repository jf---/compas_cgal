"""``ProjectionBoundarySite.build`` decides its range on integers and keeps the old sites and errors (lever 2, step 4)."""

from __future__ import annotations

from collections.abc import Callable
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks import held_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import InvalidFigure5RawGuideError
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_figure5_raw_guide import ProjectionBoundaryVertexId


def _oracle_build(
    segment_id: ProjectionBoundarySegmentId,
    parameter: Fraction,
    segment_count: int,
) -> ProjectionBoundarySite:
    """Verbatim pre-step-4 ``ProjectionBoundarySite.build`` body."""
    if segment_count <= 0 or not 0 <= int(segment_id) < segment_count:
        raise InvalidFigure5RawGuideError("Projected boundary site has no segment in the closed ring.")
    if not Fraction(0) <= parameter <= Fraction(1):
        raise InvalidFigure5RawGuideError("Projected boundary side parameter must lie in the closed unit interval.")
    canonical_segment = int(segment_id)
    canonical_parameter = parameter
    if parameter == 1:
        canonical_segment = (canonical_segment + 1) % segment_count
        canonical_parameter = Fraction(0)
    vertex_id = ProjectionBoundaryVertexId(canonical_segment) if canonical_parameter == 0 else None
    return ProjectionBoundarySite(
        ProjectionBoundarySegmentId(canonical_segment),
        canonical_parameter,
        vertex_id,
    )


def _outcome(call: Callable[[], ProjectionBoundarySite]) -> tuple[object, ...]:
    try:
        site = call()
    except InvalidFigure5RawGuideError as error:
        return ("error", str(error))
    return ("site", repr(site), type(site.parameter), hash(site))


_PARAMETERS = st.one_of(
    st.fractions(min_value=Fraction(-2), max_value=Fraction(3)),
    st.fractions(min_value=Fraction(0), max_value=Fraction(1)),
    st.sampled_from((Fraction(0), Fraction(1), Fraction(-1, 2**80), Fraction(2**80 + 1, 2**80), Fraction(2**80 - 1, 2**80))),
    st.integers(min_value=-2, max_value=3),
    st.floats(min_value=0.0, max_value=1.0).map(Fraction),
)


@settings(max_examples=1500, deadline=None)
@given(segment=st.integers(min_value=-3, max_value=12), parameter=_PARAMETERS, segment_count=st.integers(min_value=-1, max_value=12))
def test_integer_build_matches_oracle(segment: int, parameter: Fraction, segment_count: int) -> None:
    segment_id = ProjectionBoundarySegmentId(segment)

    assert _outcome(lambda: ProjectionBoundarySite.build(segment_id, parameter, segment_count)) == _outcome(lambda: _oracle_build(segment_id, parameter, segment_count))


def test_build_constructs_no_fraction(monkeypatch: pytest.MonkeyPatch) -> None:
    """The range, the end-of-side and the vertex tests are integer tests on the given rational."""

    def _refuse(*_args: object) -> Fraction:
        raise AssertionError("ProjectionBoundarySite.build constructed a Fraction")

    monkeypatch.setattr(held_figure5_raw_guide, "Fraction", _refuse)

    inside = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(1), Fraction(1, 3), 4)
    at_end = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(3), Fraction(1), 4)

    assert inside == ProjectionBoundarySite(ProjectionBoundarySegmentId(1), Fraction(1, 3), None)
    assert at_end == ProjectionBoundarySite(ProjectionBoundarySegmentId(0), Fraction(0), ProjectionBoundaryVertexId(0))
