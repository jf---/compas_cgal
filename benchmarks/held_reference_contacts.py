"""Qualify draft offset contacts without erasing rejected input evidence."""

from __future__ import annotations

from dataclasses import dataclass

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_raw_guide import ProjectionAdmissibleBoundaryHypothesis
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_offset_contact_projection import project_offset_contacts
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


class MissingReferenceContactRunError(BenchmarkError):
    """The qualified offset contacts cannot represent every input guide run."""


@dataclass(frozen=True)
class QualifiedReferenceContacts:
    """Bounded contacts with their offset sites (parallel tuples) and the rejected evidence."""

    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...]
    sites: tuple[ProjectionBoundarySite, ...]
    rejected: tuple[ProjectionAdmissibleBoundaryHypothesis, ...]


def qualify_reference_contact_sites(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> QualifiedReferenceContacts:
    """Separate bounded contacts from distance rejections, in input order.

    Projection, tie and distance decisions are exact, in one clearance batch
    (``project_offset_contacts``); site parameters are the correctly rounded
    exact parameters. A rejected hypothesis is preserved verbatim; ambiguous
    projections and malformed geometry raise. Empty input or loss of an
    entire source run also raises. Equal input hypotheses stay separate
    entries; ``qualify_reference_contacts`` merges them.
    """
    if not hypotheses:
        raise MissingReferenceContactRunError("Contact qualification requires input guide hypotheses.")
    projected = project_offset_contacts(tuple(hypothesis.contact_point for hypothesis in hypotheses), boundary, evidence_bound)
    qualified = tuple((hypothesis, site) for hypothesis, site in zip(hypotheses, projected, strict=True) if site is not None)
    missing = {hypothesis.run_id for hypothesis in hypotheses} - {hypothesis.run_id for hypothesis, _ in qualified}
    if missing:
        raise MissingReferenceContactRunError(f"No bounded offset contact for guide runs {sorted(map(int, missing))}.")
    return QualifiedReferenceContacts(
        tuple(hypothesis for hypothesis, _ in qualified),
        tuple(site for _, site in qualified),
        tuple(hypothesis for hypothesis, site in zip(hypotheses, projected, strict=True) if site is None),
    )


def qualify_reference_contacts(
    hypotheses: tuple[ProjectionAdmissibleBoundaryHypothesis, ...],
    boundary: tuple[Point2[WorldXY], ...],
    evidence_bound: Millimetre,
) -> tuple[dict[ProjectionAdmissibleBoundaryHypothesis, ProjectionBoundarySite], tuple[ProjectionAdmissibleBoundaryHypothesis, ...]]:
    """``qualify_reference_contact_sites`` keyed by hypothesis.

    Raises:
        MissingReferenceContactRunError: empty input, or a run with no
            bounded contact.
        _circle_geometry_2.AmbiguousBoundaryProjectionError: a bounded
            contact has distinct nearest points at equal distance.
    """
    qualified = qualify_reference_contact_sites(hypotheses, boundary, evidence_bound)
    return dict(zip(qualified.hypotheses, qualified.sites, strict=True)), qualified.rejected
