"""Contact qualification preserves evidence and fails on missing run coverage."""

from fractions import Fraction

import pytest

from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySegmentId
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_reference_contacts import MissingReferenceContactRunError
from benchmarks.held_reference_contacts import qualify_reference_contacts
from compas_cgal import _circle_geometry_2
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


BOUNDARY = tuple(Point2[WorldXY].build(x, y) for x, y in ((0, 0), (4, 0), (4, 4), (0, 4)))


def _hypothesis(run: int, contact: tuple[float, float]) -> ProjectionAdmissibleBoundaryHypothesis:
    site = ProjectionBoundarySite.build(ProjectionBoundarySegmentId(0), Fraction(1, 2), 4)
    point = Point2[WorldXY].build(*contact)
    return ProjectionAdmissibleBoundaryHypothesis(
        GuideRunId(run), GuideRunStationOrdinal(0), site, (site,), Millimetre(0), point, point, point, Point2[WorldXY].build(contact[0], contact[1] + 1), GuideRadius.build(1)
    )


def test_rejected_alternative_is_retained_as_evidence_beside_valid_contact() -> None:
    valid = _hypothesis(0, (2, 0.125))
    invalid = _hypothesis(0, (2, -1))
    sites, rejected = qualify_reference_contacts((invalid, valid), BOUNDARY, Millimetre(0.25))
    assert tuple(sites) == (valid,)
    assert sites[valid].segment_id == 0
    assert sites[valid].parameter == Fraction(1, 2)
    assert rejected == (invalid,)


def test_a_run_without_any_qualified_contact_fails_loudly() -> None:
    with pytest.raises(MissingReferenceContactRunError, match="7"):
        qualify_reference_contacts((_hypothesis(0, (2, 0)), _hypothesis(7, (2, -1))), BOUNDARY, Millimetre(0.25))


def test_ambiguous_contacts_are_not_discarded_as_out_of_bound() -> None:
    with pytest.raises(_circle_geometry_2.AmbiguousBoundaryProjectionError):
        qualify_reference_contacts((_hypothesis(0, (2, 2)),), BOUNDARY, Millimetre(3))


def test_empty_hypotheses_cannot_pass_contact_qualification() -> None:
    with pytest.raises(MissingReferenceContactRunError):
        qualify_reference_contacts((), BOUNDARY, Millimetre(0.25))
