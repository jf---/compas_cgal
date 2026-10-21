"""Exact ∂Ω* connectors: the CCW sub-chain of ∂Ω* between two sites, segments and arcs.

User ruling 2026-10-17: connectors are exact ∂Ω* chains (not chords). A
connector from site ``a`` to site ``b`` runs CCW along Ω*'s one loop: the rest
of ``a``'s piece, every piece in between, and the head of ``b``'s piece. Its
sub-pieces carry exact one-root endpoints and their rational supports, so its
clearance from ∂P is exactly r by construction of Ω*, and is re-decided by the
declared-boundary audit (:mod:`benchmarks.held_omega_gouge_audit`).

The loop closes with a seam connector from the last circle back to the first
(:class:`benchmarks.held_figure5_boundary_path.Figure5SeamClosure`), unless
their sites coincide exactly. Two circles on one exact site form a fan: their
connector makes no progress.

Progress and length are REPORTED binary64 (``to_double`` of exact quantities,
``atan2`` for arc sweeps) and never decide.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Self

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_clearance import ArcElement
from benchmarks.held_declared_clearance import Element
from benchmarks.held_declared_clearance import SegmentElement
from benchmarks.held_exact_centre_domain import CentreDomainPiece
from benchmarks.held_figure5_boundary_path import Figure5SeamClosure
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import Millimetre


class OmegaConnectorError(BenchmarkError):
    """Base error of exact ∂Ω* connectors."""


class EmptyOmegaPathError(OmegaConnectorError):
    """No circles to connect."""


@dataclass(frozen=True)
class OmegaConnectorPiece:
    """The part of ∂Ω* piece ``index`` from ``start`` to ``end`` (exact, in travel order)."""

    index: int
    piece: CentreDomainPiece
    start: OneRootPoint2
    end: OneRootPoint2

    @property
    def full(self) -> bool:
        """A whole closed circle piece travelled from its source back to it."""
        return self.start == self.end

    def element(self) -> Element:
        """The exact clearance element of this sub-piece."""
        support = self.piece.support
        if isinstance(support, OffsetLineSupport):
            return SegmentElement(self.start, self.end, support.direction, support.normal, support.level)
        return ArcElement(support.centre, support.radius, support.sense, self.start, self.end)

    def reported_length(self) -> float:
        """Binary64 length (reporting only)."""
        support = self.piece.support
        if isinstance(support, OffsetLineSupport):
            return abs(float(support.parameter(self.end)) - float(support.parameter(self.start)))
        return float(support.radius) * _reported_sweep(support, self.start, self.end)


def _reported_sweep(support: OffsetCircleSupport, start: OneRootPoint2, end: OneRootPoint2) -> float:
    cx, cy = float(support.centre.x), float(support.centre.y)
    sx, sy = start.as_floats()
    ex, ey = end.as_floats()
    sweep = support.sense * (math.atan2(ey - cy, ex - cx) - math.atan2(sy - cy, sx - cx))
    sweep %= math.tau
    return math.tau if start == end else sweep


@dataclass(frozen=True)
class OmegaConnector:
    """CCW motion along ∂Ω* from one circle's site to the next; ``within_fan`` when they share the site."""

    pieces: tuple[OmegaConnectorPiece, ...]
    start: OneRootPoint2
    end: OneRootPoint2
    within_fan: bool

    @property
    def boundary_progress(self) -> Millimetre:
        """Reported binary64 length along ∂Ω*."""
        return Millimetre(sum(piece.reported_length() for piece in self.pieces))


def _chain(pieces: tuple[CentreDomainPiece, ...], start: OmegaSite, end: OmegaSite) -> tuple[OmegaConnectorPiece, ...]:
    """Sub-pieces CCW from ``start`` to ``end``.

    Forward on one piece it is the one sub-piece between them. Otherwise it is
    the rest of ``start``'s piece (the whole piece from its source, a full
    circle for a closed one), every piece strictly between (all others when
    the chain wraps onto ``start``'s piece), and the head of ``end``'s piece
    up to ``end`` (nothing when ``end`` is that piece's source).
    """
    head, tail = start.point(), end.point()
    if start.index == end.index and start.key < end.key:
        return (OmegaConnectorPiece(start.index, start.piece, head, tail),)
    count = len(pieces)
    chain = [OmegaConnectorPiece(start.index, start.piece, head, start.piece.target)]
    between = range(1, count) if start.index == end.index else range(1, (end.index - start.index) % count)
    for offset in between:
        index = (start.index + offset) % count
        chain.append(OmegaConnectorPiece(index, pieces[index], pieces[index].source, pieces[index].target))
    if not end.at_vertex:
        chain.append(OmegaConnectorPiece(end.index, end.piece, end.piece.source, tail))
    return tuple(chain)


def connect(pieces: tuple[CentreDomainPiece, ...], start: OmegaCircle, end: OmegaCircle) -> OmegaConnector:
    """The connector from ``start``'s site to ``end``'s: a fan when they share it, a full lap when ``start is end``."""
    first, second = start.site, end.site
    if start is not end and first.same_position(second):
        point = first.point()
        return OmegaConnector((), point, point, True)
    return OmegaConnector(_chain(pieces, first, second), first.point(), second.point(), False)


@dataclass(frozen=True)
class OmegaConnectors:
    """The connectors of one closed traversal; ``connectors[i]`` runs ``circles[i]`` to ``circles[i + 1]`` (the last closes the seam)."""

    connectors: tuple[OmegaConnector, ...]
    seam: Figure5SeamClosure
    seam_progress: Millimetre

    @classmethod
    def build(cls, pieces: tuple[CentreDomainPiece, ...], circles: tuple[OmegaCircle, ...]) -> Self:
        """Connect consecutive circles CCW along ∂Ω* and close the seam.

        Raises:
            EmptyOmegaPathError: If there are no circles.
        """
        if not circles:
            raise EmptyOmegaPathError("An Ω* path needs at least one circle.")
        opened = tuple(connect(pieces, first, second) for first, second in zip(circles, circles[1:]))
        last, first = circles[-1], circles[0]
        if last is not first and last.site.same_position(first.site):
            return cls(opened, Figure5SeamClosure.COINCIDENT, Millimetre(0.0))
        closing = connect(pieces, last, first)
        return cls((*opened, closing), Figure5SeamClosure.CONNECTOR, closing.boundary_progress)
