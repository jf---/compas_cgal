"""Exact connector-engagement certificate for Ω* paths: ∂Ω* lines and arcs, plus the circle junctions.

An Ω* connector runs along exact ∂Ω* sub-pieces between exact sites
(:mod:`benchmarks.held_omega_connectors`). Its circles are binary64 and pass
their site within the centre's rounding, not through it, so the motion model
adds one JUNCTION line at each end (``docs/held_connector_engagement.md``,
Junctions):

- ``J_k = c_k + rho_k·u_k`` is a RATIONAL point exactly on circle ``k``:
  ``u_k`` is an exactly unit Pythagorean direction proposed from the binary64
  direction of the site (it only chooses where on the circle the tool leaves;
  any point of the circle is a valid model, and being on it is exact by
  construction);
- connector ``i`` is ``J_i -> q_i`` (line), the ∂Ω* sub-pieces ``q_i -> q_{i+1}``
  (line and arc sweeps on their rational supports), then ``q_{i+1} -> J_{i+1}``
  (line). Both junction lines join a rational point to a one-root site, so
  they live in the site's extension. A fan (two circles on one site) is the
  two junction lines.

Circle ``i + 1`` is not cut before its entry junction, so the priors stay
``k <= i``; the seam ends at ``J_0`` and circle 0 is not re-cut.
"""

from __future__ import annotations

import math
from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_centre_domain_supports import PieceKind
from benchmarks.held_connector_engagement_arcs import ArcSweep
from benchmarks.held_connector_engagement_exact import ExactConnector
from benchmarks.held_connector_engagement_exact import ExactConnectorAssessment
from benchmarks.held_connector_engagement_exact import ExactConnectorCertificate
from benchmarks.held_connector_engagement_exact import ExactConnectorPath
from benchmarks.held_connector_engagement_exact import ExactSweep
from benchmarks.held_connector_engagement_lines import LineSweep
from benchmarks.held_connector_engagement_lines import RationalLine
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_connector_engagement_range import points_apart
from benchmarks.held_declared_boundary_emission import choose_direction
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_connectors import OmegaConnectorPiece
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.motion import EngagementCap

_ARC_KINDS = {PieceKind.OFFSET_ARC: SweepKind.OFFSET_ARC, PieceKind.JUNCTION_ARC: SweepKind.JUNCTION_ARC}


def circle_junction_point(centre: tuple[float, float], radius: float, site: OneRootPoint2) -> OneRootPoint2:
    """The rational point ``c + rho·u`` exactly on the binary64 circle ``(centre, radius)``, towards ``site``.

    ``u`` is an exactly unit rational direction chosen from the binary64
    direction ``site - centre`` (a proposal: any ``u`` gives a point exactly on
    the circle; this one keeps the junction line as short as the grid allows).
    """
    sx, sy = site.as_floats()
    ux, uy = choose_direction(sx - centre[0], sy - centre[1]).vector
    rho = Fraction(radius)
    return OneRootPoint2.rational(Fraction(centre[0]) + rho * ux, Fraction(centre[1]) + rho * uy)


def sweep_of(piece: OmegaConnectorPiece) -> ExactSweep:
    """The exact sweep of one ∂Ω* sub-piece on its rational support."""
    support = piece.piece.support
    if isinstance(support, OffsetLineSupport):
        return LineSweep.build(piece.start, piece.end, RationalLine.build(support.anchor, support.direction), SweepKind.OFFSET_LINE)
    return ArcSweep.build(support.centre, support.radius, support.sense, piece.start, piece.end, _ARC_KINDS[piece.piece.kind])


def _junction(first: OneRootPoint2, second: OneRootPoint2) -> tuple[LineSweep, ...]:
    return (LineSweep.build(first, second, kind=SweepKind.CIRCLE_JUNCTION),) if points_apart(first, second) else ()


def omega_connector_path(circles: tuple[OmegaCircle, ...], connectors: OmegaConnectors) -> ExactConnectorPath:
    """The exact connector path of an Ω* traversal, junction lines included."""
    prior = tuple(PriorCircle.from_binary64(circle.center, float(circle.radius.value)) for circle in circles)
    junctions = tuple(circle_junction_point((float(circle.center.x), float(circle.center.y)), float(circle.radius.value), circle.site.point()) for circle in circles)
    exact = []
    for i, connector in enumerate(connectors.connectors):
        following = junctions[(i + 1) % len(circles)]
        sweeps = (*_junction(junctions[i], connector.start), *(sweep_of(piece) for piece in connector.pieces), *_junction(connector.end, following))
        exact.append(ExactConnector.build(sweeps))
    return ExactConnectorPath.build(prior, tuple(exact), junctions)


def _tool_and_ratio(path: OmegaGeneratedToolpath) -> tuple[Fraction, float]:
    return Fraction(float(path.case.tool_radius.value)), EngagementCap.build(math.radians(float(path.case.tea_cap))).chord_ratio


def assess_omega_connectors(path: OmegaGeneratedToolpath) -> ExactConnectorAssessment:
    """Tier 1 / 1b verdicts of every connector of an Ω* path (never raises on a refusal)."""
    tool, ratio = _tool_and_ratio(path)
    return ExactConnectorAssessment.build(omega_connector_path(path.circles, path.connectors), tool, ratio)


def require_exact_omega_connector_certificate(path: OmegaGeneratedToolpath) -> ExactConnectorCertificate:
    """Post-gate on Ω*: every connector, seam and junctions included, certified by tier 1 or 1b.

    Raises:
        UncertifiedExactConnectorError: a connector is left without a
            certificate; every refused sweep is named with its kind.
    """
    tool, ratio = _tool_and_ratio(path)
    return ExactConnectorCertificate.build(omega_connector_path(path.circles, path.connectors), tool, ratio)
