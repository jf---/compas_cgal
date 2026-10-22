"""Binary64 machining circles placed on exact ∂Ω* sites, contained in Ω* exactly.

User ruling 2026-10-17: the certified object is the exact path; circles stay
binary64 and must lie in Ω*, decided exactly. A circle for a site ``q`` with
guide radius ``R`` is centred on the inward normal: ``c* = q + R·n`` (rational on
a line, one-root on a circle). Its binary64 representation is then made safe by
two exact steps, never by a tolerance:

1. **Local support.** ``c`` is the binary64 nearest ``c*``; the radius is the
   largest binary64 ``<= min(R, bound)``, where ``bound`` is the exact clearance
   of ``c`` from the piece's whole support (``n̂·c - level`` on a line,
   ``√|c - C|² ∓ ρ'`` on a circle). This is directed rounding of an exact
   quantity: the circle is on the interior side of its own support exactly,
   and touches it to within the centre's rounding.
2. **Global.** ``dist(circle, ∂P) >= r`` is decided by a
   :class:`benchmarks.held_declared_clearance.CircleClearanceOracle` (the
   native oracle in generation, the Python oracle as its reference) and one
   circle point is located in Ω* (``CentreDomain.locate``); together they decide
   ``circle ⊂ Ω*`` exactly. When it fails (a convex vertex of Ω* or a neck
   nearby), the circle is shrunk along the same normal: the exact tangent
   circles ``q + t·n`` are nested in ``t``, so containment is monotone in ``t``
   and a bisection over the binary64 ``t`` finds a passing one. Every accepted
   circle is verified exactly; the bisection only searches.
   :func:`contain_all` decides a whole set of circles in one oracle batch and
   searches only the ones that fail.

The circle touches ∂Ω* within the centre's rounding of the exact site, not at
it: the circle-to-connector junction carries that sub-ulp gap, reported, never
deciding.
"""

from __future__ import annotations

import math
import struct
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from benchmarks.errors import BenchmarkError
from benchmarks.held_centre_domain_native import OUTSIDE_PIECE
from benchmarks.held_centre_domain_native import PROBE_AT_CIRCLE_CENTRE
from benchmarks.held_centre_domain_native import native_place_binary64
from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import CircleClearanceOracle
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_dyadic_interval import Interval
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_exact_centre_domain import PointLocation
from benchmarks.held_figure5_raw_guide import GuideRunId
from benchmarks.held_figure5_raw_guide import GuideRunStationOrdinal
from benchmarks.held_figure5_raw_guide import ProjectionBoundarySite
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_boundary import round_down
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_sites import ProbeAtCircleCentreError
from benchmarks.held_omega_sites import SiteOutsidePieceError
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_standard_placement import PaperCircleCandidate
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

# Bits of the integer square-root bracket used to approximate a one-root
# centre coordinate before rounding it to binary64: far below one ulp of any
# coordinate in the corpus (|x| < 2^10 mm), so the double is the nearest or
# adjacent one. The approximation only proposes; the exact steps decide.
APPROXIMATION_BITS = 128
# Halvings of the guide radius tried before a site is declared to admit no
# contained circle: 2^-64 R is far below the ulp of any corpus coordinate
# (|x| < 2^10 mm, ulp >= 2^-42), where a centre off the support is not
# representable at all.
MAX_SHRINK_HALVINGS = 64


class OmegaCircleError(BenchmarkError):
    """Base error of Ω* circle placement."""


class NoContainedCircleError(OmegaCircleError):
    """No binary64 circle on the site's inward normal lies in Ω* (the site is a convex vertex, or a neck closes on it)."""


def approximate(value: OneRootNumber) -> float:
    """A binary64 within an ulp of a one-root number (proposal only, never a decision)."""
    if value.is_rational:
        return float(value.rational)
    scale = 1 << APPROXIMATION_BITS
    numerator, denominator = value.radicand.numerator, value.radicand.denominator
    root = Fraction(math.isqrt(numerator * denominator * scale * scale), denominator * scale)
    return float(value.rational + value.coefficient * root)


def _largest_double(estimate: float, holds: Callable[[Fraction], bool]) -> float:
    """The largest binary64 ``x`` with ``holds(x)``, for a predicate true below a threshold (exact steps from an estimate)."""
    value = estimate
    while not holds(Fraction(value)):
        value = math.nextafter(value, -math.inf)
    while holds(Fraction(upper := math.nextafter(value, math.inf))):
        value = upper
    return value


def local_radius_bound(site: OmegaSite, centre: ExactPoint2[WorldXY], wanted: float) -> float:
    """``wanted`` when its circle about ``centre`` stays on the interior side of the site's whole support, else the largest binary64 radius that does."""
    support = site.piece.support
    holds: Callable[[Fraction], bool]
    if isinstance(support, OffsetLineSupport):
        clearance = support.normal[0] * centre.x + support.normal[1] * centre.y - support.level
        estimate = float(clearance)

        def holds(radius: Fraction) -> bool:
            return radius <= clearance

    else:
        dx, dy = centre.minus(support.centre)
        squared = dx * dx + dy * dy
        rho = support.radius
        if support.sense == -1:
            # Ω* is outside the circle: R <= √D - ρ'  ⇔  R + ρ' <= 0 or (R + ρ')² <= D.
            estimate = math.sqrt(float(squared)) - float(rho)

            def holds(radius: Fraction) -> bool:
                return radius + rho <= 0 or (radius + rho) ** 2 <= squared

        else:
            # Ω* is inside the circle: R <= ρ' - √D  ⇔  ρ' - R >= 0 and D <= (ρ' - R)².
            estimate = float(rho) - math.sqrt(float(squared))

            def holds(radius: Fraction) -> bool:
                return rho - radius >= 0 and squared <= (rho - radius) ** 2

    if holds(Fraction(wanted)):
        return wanted
    return _largest_double(min(estimate, wanted), holds)


@dataclass(frozen=True)
class OmegaCircle:
    """One binary64 machining circle based at an exact ∂Ω* site.

    ``contact_point`` is the circle's binary64 point towards the site (the
    engagement model's ``q``; reporting geometry). ``limited`` records that the
    global containment shrank the guide radius.
    """

    run_id: GuideRunId
    station_ordinal: GuideRunStationOrdinal
    source_boundary_site: ProjectionBoundarySite
    site: OmegaSite
    center: Point2[WorldXY]
    contact_point: Point2[WorldXY]
    radius: GuideRadius
    limited: bool = False
    clockwise: Literal[False] = False

    def paper_candidate(self) -> PaperCircleCandidate:
        """The standard-placement view ``(c, ρ, q)``."""
        return PaperCircleCandidate.build(center=self.center, guide_radius=self.radius, contact_point=self.contact_point)

    def element(self) -> CircleElement:
        """The exact clearance element (the binary64 values injected exactly)."""
        return CircleElement(ExactPoint2[WorldXY].from_binary64(float(self.center.x), float(self.center.y)), Fraction(float(self.radius.value)))


@dataclass(frozen=True)
class _Ray:
    """The exact ray ``q + t·n`` of tangent-circle centres at a site, and the binary64 view of ``n``.

    Built without materialising a one-root site when the site is interior: on
    a line ``q`` and ``n̂`` are rational; on a circle with probe direction ``w``
    the centre is ``C + (ρ' - σt)·w/√|w|²`` (one root).
    """

    site: OmegaSite
    normal_view: tuple[float, float]

    @classmethod
    def build(cls, site: OmegaSite) -> _Ray:
        support = site.piece.support
        if isinstance(support, OffsetLineSupport):
            return cls(site, (float(support.normal[0]), float(support.normal[1])))
        if site.probe is None:
            normal = site.inward_normal()
            return cls(site, (approximate(normal[0]), approximate(normal[1])))
        wx, wy = site.probe.minus(support.centre)
        length = math.sqrt(float(wx * wx + wy * wy))
        return cls(site, (-support.sense * float(wx) / length, -support.sense * float(wy) / length))

    def centre(self, distance: float) -> ExactPoint2[WorldXY]:
        """The binary64 point nearest ``q + t·n``."""
        site, step = self.site, Fraction(distance)
        support = site.piece.support
        if site.probe is None:
            point, normal = site.point(), site.inward_normal()
            return ExactPoint2[WorldXY].from_binary64(approximate(point.x + step * normal[0]), approximate(point.y + step * normal[1]))
        if isinstance(support, OffsetLineSupport):
            along = site.key.along
            assert along is not None
            foot = support.anchor.offset(support.direction, along).offset(support.normal, step)
            return ExactPoint2[WorldXY].from_binary64(float(foot.x), float(foot.y))
        # Proposal from enclosures: C + (ρ' - σt)·w/|w| (the exact steps decide).
        frame = site.frame
        wx = Interval.of(site.probe.x) - frame.origin_box[0]
        wy = Interval.of(site.probe.y) - frame.origin_box[1]
        length = (wx.square() + wy.square()).sqrt()
        if length.lo <= 0:
            # A probe within the enclosure resolution of the centre: propose from the exact point.
            point, normal = site.point(), site.inward_normal()
            return ExactPoint2[WorldXY].from_binary64(approximate(point.x + step * normal[0]), approximate(point.y + step * normal[1]))
        scale = Interval.of(support.radius - support.sense * step) / length
        return ExactPoint2[WorldXY].from_binary64((frame.origin_box[0] + scale * wx).midpoint(), (frame.origin_box[1] + scale * wy).midpoint())


def _proposal(ray: _Ray, distance: float) -> tuple[ExactPoint2[WorldXY], float] | None:
    """Centre nearest ``q + t·n`` and its locally safe radius ``min(t, bound)``; ``None`` when no positive radius fits."""
    centre = ray.centre(distance)
    radius = local_radius_bound(ray.site, centre, distance)
    return (centre, radius) if radius > 0 else None


def _circle(ray: _Ray, centre: ExactPoint2[WorldXY], radius: float, limited: bool, provenance: tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]) -> OmegaCircle:
    site = ray.site
    cx, cy = float(centre.x), float(centre.y)
    contact = Point2[WorldXY].build(cx - radius * ray.normal_view[0], cy - radius * ray.normal_view[1])
    run_id, station, source = provenance
    return OmegaCircle(run_id, station, source, site, Point2[WorldXY].build(cx, cy), contact, GuideRadius.build(radius), limited)


def place_local(site: OmegaSite, guide_radius: float, provenance: tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]) -> OmegaCircle:
    """Circle on the site's normal, safe against its own support only (global containment is :func:`contain`).

    Raises:
        NoContainedCircleError: If no positive radius fits the support at the rounded centre.
    """
    ray = _Ray.build(site)
    distance, limited = guide_radius, False
    support = site.piece.support
    if isinstance(support, OffsetCircleSupport) and support.sense == 1 and support.radius < Fraction(guide_radius):
        # Inside a convex arc of radius ρ' the tangent circles q + t·n stay
        # inside its circle only for t <= ρ' (t = ρ' is the circle itself).
        distance, limited = round_down(support.radius), True
    proposal = _proposal(ray, distance)
    if proposal is None:
        raise NoContainedCircleError(f"No positive circle fits the support of piece {site.index} at radius {distance}.")
    centre, radius = proposal
    return _circle(ray, centre, radius, limited, provenance)


def place_local_batch(
    loop: OmegaLoop,
    contacts: Sequence[DecidedContact],
    guide_radii: Sequence[float],
    provenances: Sequence[tuple[GuideRunId, GuideRunStationOrdinal, ProjectionBoundarySite]],
) -> tuple[OmegaCircle, ...]:
    """:func:`place_local` of every contact's :meth:`OmegaSite.decided` site, decided in one native call.

    The native batch (``CentreDomain2.place_binary64``) makes the same exact
    decisions: the site's placement and range, the correctly rounded centre
    of ``q + t·n``, and the largest binary64 radius against the support. The
    circles are the ones :func:`place_local` builds, field for field
    (contract-tested).

    Raises:
        SiteOutsidePieceError: If a contact projects outside its piece.
        ProbeAtCircleCentreError: If a circle contact is the circle's centre.
        NoContainedCircleError: If no positive radius fits a support at the rounded centre.
    """
    placements = native_place_binary64(
        loop.domain.native, [contact.query for contact in contacts], [contact.piece.index for contact in contacts], [contact.is_vertex for contact in contacts], list(guide_radii)
    )
    circles: list[OmegaCircle] = []
    for contact, placement, (run_id, station, source) in zip(contacts, placements, provenances, strict=True):
        index = contact.piece.index
        if placement.status == OUTSIDE_PIECE:
            raise SiteOutsidePieceError(f"Contact {contact.query} projects outside piece {index}.")
        if placement.status == PROBE_AT_CIRCLE_CENTRE:
            raise ProbeAtCircleCentreError(f"Contact {contact.query} is the centre of piece {index}.")
        radius = placement.radius
        if not radius > 0:
            raise NoContainedCircleError(f"No positive circle fits the support of piece {index}.")
        site = OmegaSite.from_decision(loop.frames[index], contact.query, placement.at_source, placement.half)
        cx, cy = placement.centre
        normal_x, normal_y = placement.normal_view
        contact_point = Point2[WorldXY].build(cx - radius * normal_x, cy - radius * normal_y)
        circles.append(OmegaCircle(run_id, station, source, site, Point2[WorldXY].build(cx, cy), contact_point, GuideRadius.build(radius), placement.limited))
    return tuple(circles)


def _decided_contained(verdict: ClearanceVerdict, domain: CentreDomain, centre: ExactPoint2[WorldXY], radius: float) -> bool:
    if verdict.clearance is Clearance.GOUGE:
        return False
    # dist >= r everywhere, so the circle is wholly inside Ω* or wholly outside P: one point decides.
    return domain.locate(ExactPoint2[WorldXY].build(centre.x + Fraction(radius), centre.y)) is not PointLocation.EXTERIOR


def _contained(clearance: CircleClearanceOracle, domain: CentreDomain, centre: ExactPoint2[WorldXY], radius: float) -> bool:
    (verdict,) = clearance.circles([CircleElement(centre, Fraction(radius))])
    return _decided_contained(verdict, domain, centre, radius)


def _exact_centre(circle: OmegaCircle) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].from_binary64(float(circle.center.x), float(circle.center.y))


def _bits(value: float) -> int:
    return int(struct.unpack("<q", struct.pack("<d", value))[0])


def _from_bits(bits: int) -> float:
    return float(struct.unpack("<d", struct.pack("<q", bits))[0])


def contain(circle: OmegaCircle, clearance: CircleClearanceOracle, domain: CentreDomain) -> OmegaCircle:
    """The circle itself when it lies in Ω* exactly, else the largest found contained circle on its normal.

    Raises:
        NoContainedCircleError: If no binary64 circle on the normal is contained.
    """
    if _contained(clearance, domain, _exact_centre(circle), float(circle.radius.value)):
        return circle
    return _shrunk(circle, clearance, domain)


def contain_all(circles: Sequence[OmegaCircle], clearance: CircleClearanceOracle, domain: CentreDomain) -> tuple[OmegaCircle, ...]:
    """:func:`contain` of every circle, in order, with the first containment decision of all of them in one oracle batch.

    Raises:
        NoContainedCircleError: If a circle admits no contained circle on its normal.
    """
    centres = [_exact_centre(circle) for circle in circles]
    radii = [float(circle.radius.value) for circle in circles]
    verdicts = clearance.circles([CircleElement(centre, Fraction(radius)) for centre, radius in zip(centres, radii, strict=True)])
    return tuple(
        circle if _decided_contained(verdict, domain, centre, radius) else _shrunk(circle, clearance, domain)
        for circle, verdict, centre, radius in zip(circles, verdicts, centres, radii, strict=True)
    )


def _shrunk(circle: OmegaCircle, clearance: CircleClearanceOracle, domain: CentreDomain) -> OmegaCircle:
    """The largest found contained circle on the normal of a circle that is not contained.

    Raises:
        NoContainedCircleError: If no binary64 circle on the normal is contained.
    """
    site = circle.site
    ray = _Ray.build(site)
    provenance = (circle.run_id, circle.station_ordinal, circle.source_boundary_site)
    # Halve t from the guide radius until a circle passes, then bisect over the
    # binary64 bit patterns between the passing t and the failing one above it
    # (positive doubles order as their bit patterns). Halving first because
    # tiny t round the centre onto the support, where no radius fits.
    high = _bits(float(circle.radius.value))
    best: tuple[ExactPoint2[WorldXY], float] | None = None
    low = high
    for _ in range(MAX_SHRINK_HALVINGS):
        low = _bits(_from_bits(low) / 2)
        proposal = _proposal(ray, _from_bits(low))
        if proposal is not None and _contained(clearance, domain, *proposal):
            best = proposal
            break
        high = low
    if best is not None:
        while high - low > 1:
            middle = (low + high) // 2
            proposal = _proposal(ray, _from_bits(middle))
            if proposal is not None and _contained(clearance, domain, *proposal):
                low, best = middle, proposal
            else:
                high = middle
    if best is None:
        raise NoContainedCircleError(f"No binary64 circle on the inward normal of piece {site.index} lies in Ω*.")
    return _circle(ray, best[0], best[1], True, provenance)
