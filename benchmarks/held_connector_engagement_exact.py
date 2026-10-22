"""Tiers 1 and 1b of the connector-engagement certificate on exact connectors (line and arc sweeps).

The binary64 path (:mod:`benchmarks.held_connector_engagement`) reads
straight pieces between binary64 contacts on one integer grid. This module
reads exact connectors: chains of :class:`LineSweep` / :class:`ArcSweep` with
rational or one-root endpoints, between exact JUNCTIONS (the point at which
the tool is on circle ``k`` when it leaves or enters it). Binary64 inputs are
injected exactly (:meth:`ExactConnectorPath.from_figure5`); nothing is rounded.

The decision is the same lemma as the binary64 path, over each sweep's
realised squared-distance range (:class:`RangeEngagementLemma`): tier 1 with
the circle just cut (``k = i``); tier 1b per sub-sweep with any circle cut
before (``k <= i``; every circle for the seam), splitting lines at exact
midpoints and arcs at rational points of their circle, to
``MAX_BISECTION_DEPTH``. A refusal means "not proved", never "violated", and
it is named per sweep (:class:`ExactRefusal`).
"""

from __future__ import annotations

import enum
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

import numpy as np
import numpy.typing as npt

from benchmarks.held_connector_engagement import MAX_BISECTION_DEPTH
from benchmarks.held_connector_engagement import ConnectorTier
from benchmarks.held_connector_engagement import InvalidConnectorEngagementInputError
from benchmarks.held_connector_engagement import UncertifiedConnectorEngagementError
from benchmarks.held_connector_engagement import prior_circle_indices
from benchmarks.held_connector_engagement import reach_bracket
from benchmarks.held_connector_engagement_arcs import ArcSweep
from benchmarks.held_connector_engagement_lines import LineSweep
from benchmarks.held_connector_engagement_range import DiscontinuousConnectorError
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from benchmarks.held_connector_engagement_range import SweepKind
from benchmarks.held_dyadic_interval import FRACTION_BITS
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_figure5_boundary_path import Figure5BoundaryConnectors
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_one_root import OneRootPoint2
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

ExactSweep = LineSweep | ArcSweep
FloatArray = npt.NDArray[np.float64]


class RefusalReason(enum.Enum):
    """Why a sweep is left without a certificate."""

    NOT_PROVED_AT_DEPTH = "not proved at the bisection depth"
    SPLIT_UNRESOLVED = "no exact split point strictly inside"


@dataclass(frozen=True)
class ExactRefusal:
    """Sweep ``sweep`` of connector ``connector`` (a ``kind``) has no certificate, for ``reason``."""

    connector: int
    sweep: int
    kind: SweepKind
    reason: RefusalReason


class UncertifiedExactConnectorError(UncertifiedConnectorEngagementError):
    """At least one exact connector is left without a tier-1 or tier-1b certificate; every refused sweep is named."""

    def __init__(self, refusals: tuple[ExactRefusal, ...]) -> None:
        super().__init__(tuple(sorted({refusal.connector for refusal in refusals})))
        self.refusals = refusals


@dataclass(frozen=True)
class ExactConnector:
    """One connector: exact sweeps, each starting exactly where the previous ends; ``()`` is a stationary fan."""

    sweeps: tuple[ExactSweep, ...]

    @classmethod
    def build(cls, sweeps: tuple[ExactSweep, ...]) -> Self:
        """Validate continuity.

        Raises:
            DiscontinuousConnectorError: a sweep does not start exactly at the previous end.
        """
        for k, (first, second) in enumerate(zip(sweeps, sweeps[1:])):
            if first.end != second.start:
                raise DiscontinuousConnectorError(f"sweep {k + 1} starts at {second.start!r}, not at the previous end {first.end!r}")
        return cls(tuple(sweeps))


@dataclass(frozen=True)
class ExactConnectorPath:
    """Circles in cut order, their junctions, and the connectors between them.

    ``connectors[i]`` runs from ``junctions[i]`` (on circle ``i``) to
    ``junctions[i + 1]``; with ``len(connectors) == len(circles)`` the last
    one is the seam back to ``junctions[0]``.
    """

    circles: tuple[PriorCircle, ...]
    connectors: tuple[ExactConnector, ...]
    junctions: tuple[OneRootPoint2, ...]

    @classmethod
    def build(cls, circles: tuple[PriorCircle, ...], connectors: tuple[ExactConnector, ...], junctions: tuple[OneRootPoint2, ...]) -> Self:
        """Validate counts and that every connector joins its circles' junctions exactly.

        Raises:
            InvalidConnectorEngagementInputError: no circle, a count mismatch,
                or a connector that does not run between its junctions.
        """
        count = len(circles)
        if count == 0 or len(junctions) != count:
            raise InvalidConnectorEngagementInputError(f"{count} circles need {count} junctions, got {len(junctions)}")
        if len(connectors) not in (count, count - 1):
            raise InvalidConnectorEngagementInputError(f"{count} circles take {count} (seam) or {count - 1} connectors, not {len(connectors)}")
        for k, connector in enumerate(connectors):
            start, end = junctions[k], junctions[(k + 1) % count]
            if connector.sweeps:
                joined = connector.sweeps[0].start == start and connector.sweeps[-1].end == end
            else:
                joined = start == end
            if not joined:
                raise InvalidConnectorEngagementInputError(f"connector {k} does not run between the junctions of circles {k} and {(k + 1) % count}")
        return cls(tuple(circles), tuple(connectors), tuple(junctions))

    @property
    def closed(self) -> bool:
        """The last connector is a seam."""
        return len(self.connectors) == len(self.circles)

    @classmethod
    def from_figure5(cls, circles: tuple[Figure5CounterclockwiseCircle, ...], connectors: Figure5BoundaryConnectors) -> Self:
        """Inject a binary64 path exactly: contacts are the junctions, sample polylines the line sweeps."""
        prior = tuple(PriorCircle.from_binary64(circle.center, float(circle.radius.value)) for circle in circles)
        junctions = tuple(_injected(circle.contact_point) for circle in circles)
        exact = []
        for transition in connectors.transitions:
            vertices = [_injected(point) for point in transition.samples]
            exact.append(ExactConnector.build(tuple(LineSweep.build(a, b) for a, b in zip(vertices, vertices[1:]) if a != b)))
        return cls.build(prior, tuple(exact), junctions)


def _injected(point: Point2[WorldXY]) -> OneRootPoint2:
    return OneRootPoint2.rational(Fraction(float(point.x)), Fraction(float(point.y)))


@dataclass(frozen=True)
class ExactSubSweepCertificate:
    """The part of sweep ``sweep`` from ``start`` to ``end`` is certified by circle ``circle``."""

    sweep: int
    start: OneRootPoint2
    end: OneRootPoint2
    circle: int


@dataclass(frozen=True)
class ExactConnectorVerdict:
    """One connector's verdict: certified sub-sweeps covering every sweep, or the refused sweeps."""

    connector: int
    seam: bool
    tier: ConnectorTier
    sweeps: int
    sub_sweeps: tuple[ExactSubSweepCertificate, ...]
    refusals: tuple[ExactRefusal, ...]

    @property
    def certified(self) -> bool:
        return self.tier is not ConnectorTier.REFUSED


class _CandidateIndex:
    """Prefilter of the circles that could certify a sub-sweep, libm-free.

    A circle certifies a sweep only if the outer test passes, which needs
    ``D <= rho² + 2 rho r < (rho + 2 r)²`` at every point, the sweep's start
    included. A circle is dropped only when its centre is PROVABLY farther
    than ``rho + 2 r`` from the start: the start's interval enclosure
    (outward to binary64) against the outward-rounded reach
    (:func:`benchmarks.held_connector_engagement.reach_bracket`, IEEE
    ``+ - ×`` only). The filter chooses which exact questions to ask: it never
    drops a circle that could certify, so it loses no certificate and creates none.
    """

    def __init__(self, circles: tuple[PriorCircle, ...], tool_radius: Fraction) -> None:
        self._centres: FloatArray = np.array([(float(c.centre.x), float(c.centre.y)) for c in circles], dtype=np.float64).reshape(-1, 2)
        reach = np.array([float(c.guide_radius + 2 * tool_radius) for c in circles], dtype=np.float64)
        self._reach_low, self._reach_high = np.nextafter(reach, -np.inf), np.nextafter(reach, np.inf)

    @staticmethod
    def _outward(enclosure: Interval) -> tuple[float, float]:
        low = float(np.nextafter(float(Fraction(enclosure.lo, 1 << FRACTION_BITS)), -np.inf))
        high = float(np.nextafter(float(Fraction(enclosure.hi, 1 << FRACTION_BITS)), np.inf))
        return low, high

    def near(self, sweep: ExactSweep, last: int) -> list[int]:
        """Candidate circles ``k <= last``, most recent first."""
        stop = last + 1
        box = sweep.start_box
        beyond, _ = reach_bracket(self._centres[:stop], self._reach_low[:stop], self._reach_high[:stop], self._outward(box.x), self._outward(box.y))
        return [int(k) for k in np.nonzero(~beyond)[0][::-1]]


def _bisect(
    lemma: RangeEngagementLemma, circles: tuple[PriorCircle, ...], index: _CandidateIndex, last: int, position: int, sweep: ExactSweep
) -> tuple[tuple[ExactSubSweepCertificate, ...], RefusalReason | None]:
    """Tier 1b on one sweep: any circle ``k <= last`` per sub-sweep, exact splits, depth-first."""
    certified: list[ExactSubSweepCertificate] = []
    stack: list[tuple[ExactSweep, int]] = [(sweep, 0)]
    while stack:
        part, depth = stack.pop()
        winner = next((k for k in index.near(part, last) if lemma.certifies(part, circles[k])), None)
        if winner is not None:
            certified.append(ExactSubSweepCertificate(position, part.start, part.end, winner))
            continue
        if depth >= MAX_BISECTION_DEPTH:
            return (), RefusalReason.NOT_PROVED_AT_DEPTH
        halves = part.split()
        if halves is None:
            return (), RefusalReason.SPLIT_UNRESOLVED
        stack.append((halves[1], depth + 1))
        stack.append((halves[0], depth + 1))
    return tuple(certified), None


@dataclass(frozen=True)
class ExactConnectorAssessment:
    """Tier 1 and 1b verdicts for every exact connector, refusals included; never raises on a refusal."""

    verdicts: tuple[ExactConnectorVerdict, ...]
    seconds: float

    @classmethod
    def build(cls, path: ExactConnectorPath, tool_radius: Fraction, chord_ratio: float) -> Self:
        """Decide every connector: tier 1 with the circle just cut, else tier 1b per sweep over every prior circle.

        Raises:
            InvalidConnectorEngagementInputError: bad tool radius or chord surrogate.
        """
        started = time.perf_counter()
        lemma = RangeEngagementLemma.build(tool_radius, chord_ratio)
        index = _CandidateIndex(path.circles, lemma.tool_radius)
        seam_index = len(path.circles) - 1 if path.closed else None
        verdicts = []
        for k, connector in enumerate(path.connectors):
            seam = k == seam_index
            sweeps = connector.sweeps
            if not sweeps:
                verdicts.append(ExactConnectorVerdict(k, seam, ConnectorTier.STATIONARY, 0, (), ()))
                continue
            own = path.circles[k]
            first_try = [lemma.certifies(sweep, own) for sweep in sweeps]
            if all(first_try):
                subs = tuple(ExactSubSweepCertificate(p, sweep.start, sweep.end, k) for p, sweep in enumerate(sweeps))
                verdicts.append(ExactConnectorVerdict(k, seam, ConnectorTier.TIER_1, len(sweeps), subs, ()))
                continue
            last = prior_circle_indices(k, len(path.circles))[0]
            sub_sweeps: list[ExactSubSweepCertificate] = []
            refusals: list[ExactRefusal] = []
            for p, (sweep, by_own) in enumerate(zip(sweeps, first_try)):
                if by_own:
                    sub_sweeps.append(ExactSubSweepCertificate(p, sweep.start, sweep.end, k))
                    continue
                found, reason = _bisect(lemma, path.circles, index, last, p, sweep)
                if reason is not None:
                    refusals.append(ExactRefusal(k, p, sweep.kind, reason))
                sub_sweeps.extend(found)
            if refusals:
                verdicts.append(ExactConnectorVerdict(k, seam, ConnectorTier.REFUSED, len(sweeps), (), tuple(refusals)))
            else:
                verdicts.append(ExactConnectorVerdict(k, seam, ConnectorTier.TIER_1B, len(sweeps), tuple(sub_sweeps), ()))
        return cls(tuple(verdicts), time.perf_counter() - started)

    @property
    def refused(self) -> tuple[int, ...]:
        return tuple(verdict.connector for verdict in self.verdicts if not verdict.certified)

    @property
    def refusals(self) -> tuple[ExactRefusal, ...]:
        return tuple(refusal for verdict in self.verdicts for refusal in verdict.refusals)

    def count(self, tier: ConnectorTier) -> int:
        return sum(1 for verdict in self.verdicts if verdict.tier is tier)


@dataclass(frozen=True)
class ExactConnectorCertificate:
    """Every exact connector of a path certified within the cap, continuously, by tier 1 or 1b."""

    verdicts: tuple[ExactConnectorVerdict, ...]

    @classmethod
    def build(cls, path: ExactConnectorPath, tool_radius: Fraction, chord_ratio: float) -> Self:
        """Post-gate over all connectors, the seam included.

        Raises:
            InvalidConnectorEngagementInputError: bad tool radius or chord surrogate.
            UncertifiedExactConnectorError: a connector is left without a
                certificate; every refused sweep is named.
        """
        assessment = ExactConnectorAssessment.build(path, tool_radius, chord_ratio)
        if assessment.refusals:
            raise UncertifiedExactConnectorError(assessment.refusals)
        return cls(assessment.verdicts)
