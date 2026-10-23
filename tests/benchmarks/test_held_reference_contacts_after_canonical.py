"""Offset contacts are projected after canonicalisation, on the canonical hypotheses only.

User ruling 2026-10-22: ``build_held_reference_path`` (``qualify_contacts=False``)
projects only the contacts whose sites the path reads, those of the canonical
hypotheses. The sites it reads are unchanged; a contact that is never read can
no longer raise a bound or tie error. ``qualify_contacts=True`` filters
hypotheses by their projection before ownership, so it keeps projecting every
hypothesis.
"""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

import benchmarks.held_figure5_path as hfp
import benchmarks.held_reference_contacts as reference_contacts
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_offset_contact_projection import project_bounded_offset_contacts
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from compas_cgal import _circle_geometry_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

CASE = "figure8_upper"
# Far beyond any admissible distance gap of the reference cases (mm): the
# displaced contact's nearest offset site is outside the evidence bound.
FAR_OFFSET_MM = 1000.0


@pytest.fixture(scope="module")
def reference():
    case = load_held_reference_case(CASE)
    guide = build_held_reference_raw_guide(case)
    components = figure7_inward_offset(case).components
    hypotheses = hfp.build_guide_distance_admissible_hypotheses(guide)
    return guide, components, hypotheses


def _owner_indices(guide, components, hypotheses) -> tuple[int, ...]:
    return hfp._boundary_neighborhood_owner_indices(hypotheses, len(components[0]), guide.case.start_marker)


def _canonical_contacts(guide, components, hypotheses) -> tuple[Point2[WorldXY], ...]:
    owners = _owner_indices(guide, components, hypotheses)
    owned = tuple(hypotheses[index] for index in owners)
    return tuple(owned[index].contact_point for index, _ in hfp._canonical_indices(owned))


def _path_with_every_contact_projected(guide, components, hypotheses) -> hfp.HypothesisFigure5Path:
    """The pre-ruling construction: every hypothesis projected, sites passed parallel."""
    gap = guide.site_budget.admissible_distance_gap
    sites = project_bounded_offset_contacts(tuple(hypothesis.contact_point for hypothesis in hypotheses), components[0], gap)
    return hfp.build_hypothesis_figure5_path(
        inward_components=components,
        hypotheses=hypotheses,
        reached_run_ids=frozenset(run.run_id for run in guide.runs),
        published_start_center=guide.case.start_marker,
        published_start_radius=guide.case.start_marker_radius,
        boundary_evidence_bound=gap,
        start_evidence_bound=Millimetre(float(guide.site_budget.reconstruction_bound) + float(guide.site_budget.projection_bound)),
        tool_radius=guide.case.tool_radius,
        cap=EngagementCap.build(math.radians(float(guide.case.tea_cap))),
        offset_sites=sites,
        preserve_source_runs=True,
    )


def _spy(monkeypatch, module, name):
    calls = []
    original = getattr(module, name)

    def spy(contacts, *args, **kwargs):
        result = original(contacts, *args, **kwargs)
        calls.append((tuple(contacts), result))
        return result

    monkeypatch.setattr(module, name, spy)
    return calls


def test_only_canonical_contacts_are_projected_and_their_sites_are_unchanged(reference, monkeypatch) -> None:
    guide, components, hypotheses = reference
    canonical_contacts = _canonical_contacts(guide, components, hypotheses)
    assert len(canonical_contacts) < len(hypotheses), "fixture must have non-canonical hypotheses"
    every = project_bounded_offset_contacts(tuple(h.contact_point for h in hypotheses), components[0], guide.site_budget.admissible_distance_gap)
    site_of = dict(zip((h.contact_point for h in hypotheses), every, strict=True))
    before = _path_with_every_contact_projected(guide, components, hypotheses)

    calls = _spy(monkeypatch, hfp, "project_bounded_offset_contacts")
    after = hfp.build_held_reference_path(guide, components)

    assert len(calls) == 1
    contacts, sites = calls[0]
    assert contacts == canonical_contacts
    assert sites == tuple(site_of[contact] for contact in contacts)
    assert after == before


def test_unread_contact_beyond_the_bound_no_longer_raises(reference, monkeypatch) -> None:
    guide, components, hypotheses = reference
    owners = set(_owner_indices(guide, components, hypotheses))
    unread = next(index for index in range(len(hypotheses)) if index not in owners)
    contact = hypotheses[unread].contact_point
    far = Point2[WorldXY].build(float(contact.x) + FAR_OFFSET_MM, float(contact.y) + FAR_OFFSET_MM)
    displaced = hypotheses[:unread] + (replace(hypotheses[unread], contact_point=far),) + hypotheses[unread + 1 :]
    with pytest.raises(_circle_geometry_2.BoundaryProjectionDistanceError):
        _path_with_every_contact_projected(guide, components, displaced)
    baseline = hfp.build_held_reference_path(guide, components)

    monkeypatch.setattr(hfp, "build_guide_distance_admissible_hypotheses", lambda _guide: displaced)
    result = hfp.build_held_reference_path(guide, components)

    assert result == baseline


def test_qualified_contacts_still_project_every_hypothesis(reference, monkeypatch) -> None:
    guide, components, hypotheses = reference
    bounded = _spy(monkeypatch, hfp, "project_bounded_offset_contacts")
    qualified = _spy(monkeypatch, reference_contacts, "project_offset_contacts")

    hfp.build_held_reference_path(guide, components, qualify_contacts=True)

    assert bounded == []
    assert len(qualified) == 1
    assert qualified[0][0] == tuple(hypothesis.contact_point for hypothesis in hypotheses)
