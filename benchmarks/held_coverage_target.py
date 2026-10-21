"""Adapter: the exact Held centre domain Omega* as the coverage certificate's exact target T = Omega* (+) D_r.

The certificate (``compas_cgal.coverage_certificate``) takes the declared
boundary, every dOmega* piece with its one-root end vertices, and the
convex-vertex disks; this module reads them off :class:`CentreDomain` and its
:class:`DeclaredBoundary` without any rounding (every value is a Fraction or a
one-root number on both sides).
"""

from __future__ import annotations

from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.coverage_certificate.exact_target import CentreDomainPiece
from compas_cgal.coverage_certificate.exact_target import DeclaredCircleArc
from compas_cgal.coverage_certificate.exact_target import DeclaredPrimitive
from compas_cgal.coverage_certificate.exact_target import DeclaredSegment
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.exact_target import ExactXY
from compas_cgal.coverage_certificate.exact_target import OffsetCirclePiece
from compas_cgal.coverage_certificate.exact_target import OffsetLinePiece
from compas_cgal.coverage_certificate.exact_target import OneRootPoint
from compas_cgal.coverage_certificate.target import OneRootMillimetre
from compas_cgal.coverage_certificate.target import VertexDisk
from compas_cgal.units import ExactMillimetre
from compas_cgal.units import WorldXY


def _xy(point: ExactPoint2[WorldXY]) -> ExactXY:
    return ExactMillimetre(Fraction(point.x)), ExactMillimetre(Fraction(point.y))


def _mm(value: OneRootNumber) -> OneRootMillimetre:
    return OneRootMillimetre(value.rational, value.coefficient)


def _one_root(point: OneRootPoint2) -> OneRootPoint:
    return OneRootPoint.build(_mm(point.x), _mm(point.y), point.radicand)


def _primitive(primitive: DeclaredLine | DeclaredArc) -> DeclaredPrimitive:
    if isinstance(primitive, DeclaredLine):
        return DeclaredSegment.build(_xy(primitive.start), primitive.direction.vector, ExactMillimetre(Fraction(primitive.length)))
    return DeclaredCircleArc.build(_xy(primitive.start), _xy(primitive.end), _xy(primitive.centre), ExactMillimetre(Fraction(primitive.radius)), primitive.sense)


def exact_coverage_target(boundary: DeclaredBoundary, domain: CentreDomain) -> ExactCentreDomainTarget:
    """T = Omega* (+) D_r of ``domain``, built on the declared ``boundary`` it was constructed from.

    Raises:
        InvalidExactTargetError: an exact part fails the facade's validation.
    """
    pieces: list[CentreDomainPiece] = []
    for piece in domain.pieces:
        support = piece.support
        source, target = _one_root(piece.source), _one_root(piece.target)
        if isinstance(support, OffsetLineSupport):
            pieces.append(OffsetLinePiece.build(support.normal, ExactMillimetre(Fraction(support.level)), source, target))
        else:
            pieces.append(
                OffsetCirclePiece.build(
                    _xy(support.centre), ExactMillimetre(Fraction(support.radius)), ExactMillimetre(Fraction(support.boundary_radius)), support.sense, source, target
                )
            )
    disks = [VertexDisk.build(_mm(disk.centre.x), _mm(disk.centre.y), disk.centre.radicand) for disk in domain.vertex_disks]
    return ExactCentreDomainTarget.build([_primitive(p) for p in boundary.primitives], pieces, disks, components=len(domain.components))
