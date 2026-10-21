"""Members of the machinable target T = Omega (+) D_r: the centre domain, its boundary sites, and one-root vertex disks."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from compas_cgal import _coverage_certificate_2
from compas_cgal.coverage_certificate.errors import InvalidCentreBoundarySiteError
from compas_cgal.coverage_certificate.errors import InvalidCentreDomainError
from compas_cgal.coverage_certificate.errors import InvalidVertexDiskError
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY

Loop = tuple[Point2[WorldXY], ...]


def loop_coordinates(loop: Loop) -> list[tuple[float, float]]:
    """The binary64 wire form of a loop; each double is injected exactly natively."""
    return [(float(point.x), float(point.y)) for point in loop]


@dataclass(frozen=True)
class CentreDomain:
    """Validated centre domain Omega in world XY mm: outer loop CCW, holes CW.

    Loops are open (the last vertex connects back to the first). Construct
    through ``build``; the raw constructor does not validate.
    """

    outer: Loop
    holes: tuple[Loop, ...]

    @classmethod
    def build(cls, outer: Sequence[Point2[WorldXY]], holes: Sequence[Sequence[Point2[WorldXY]]] = ()) -> CentreDomain:
        """Validate the polygon invariants with exact predicates.

        Raises:
            InvalidCentreDomainError: a loop is not simple, the outer loop is not
                CCW, a hole is not CW or not strictly inside, loops touch, a loop
                has fewer than three vertices or a zero-length edge.
        """
        domain = cls(tuple(outer), tuple(tuple(hole) for hole in holes))
        try:
            _coverage_certificate_2.validate_centre_domain(loop_coordinates(domain.outer), [loop_coordinates(h) for h in domain.holes])
        except _coverage_certificate_2.InvalidCentreDomainError as error:
            raise InvalidCentreDomainError(str(error)) from error
        return domain

    @property
    def loops(self) -> tuple[Loop, ...]:
        """Loop 0 is the outer loop; loops 1.. are the holes."""
        return (self.outer, *self.holes)


@dataclass(frozen=True)
class CentreBoundarySite:
    """An exact point on loop ``loop``, segment ``segment``: v_s + parameter (v_{s+1} - v_s)."""

    loop: int
    segment: int
    parameter: Fraction

    @classmethod
    def build(cls, loop: int, segment: int, parameter: Fraction) -> CentreBoundarySite:
        """Build an exact boundary site; the loop and segment are checked against a domain later.

        Raises:
            InvalidCentreBoundarySiteError: a negative index, a non-Fraction
                parameter, or a parameter outside [0, 1).
        """
        if loop < 0 or segment < 0:
            raise InvalidCentreBoundarySiteError(f"site indices must be non-negative, got loop {loop}, segment {segment}")
        if not isinstance(parameter, Fraction):
            raise InvalidCentreBoundarySiteError(f"site parameter must be an exact Fraction, got {type(parameter).__name__}")
        if not Fraction(0) <= parameter < Fraction(1):
            raise InvalidCentreBoundarySiteError(f"site parameter must lie in [0, 1), got {parameter}")
        return cls(loop, segment, parameter)

    def check_in(self, domain: CentreDomain) -> None:
        """Check that ``domain`` has this site's loop and segment.

        Raises:
            InvalidCentreBoundarySiteError: the domain has no such loop or segment.
        """
        if self.loop >= len(domain.loops):
            raise InvalidCentreBoundarySiteError(f"domain has {len(domain.loops)} loop(s), site names loop {self.loop}")
        if self.segment >= len(domain.loops[self.loop]):
            raise InvalidCentreBoundarySiteError(f"loop {self.loop} has {len(domain.loops[self.loop])} segment(s), site names segment {self.segment}")


@dataclass(frozen=True)
class OneRootMillimetre:
    """A world-XY mm coordinate rational + coefficient * sqrt(radicand); the radicand lives on the point."""

    rational: Fraction
    coefficient: Fraction


@dataclass(frozen=True)
class VertexDisk:
    """Target member D_r(w) about a convex vertex w of the exact centre domain Omega*.

    The centre is a one-root point: both coordinates lie in Q(sqrt(radicand)).
    Its radius is the tool radius of the certificate it is passed to.
    """

    x: OneRootMillimetre
    y: OneRootMillimetre
    radicand: Fraction

    @classmethod
    def build(cls, x: OneRootMillimetre, y: OneRootMillimetre, radicand: Fraction) -> VertexDisk:
        """Build a vertex disk about a one-root centre.

        Raises:
            InvalidVertexDiskError: a part is not a Fraction or the radicand is negative.
        """
        for value in (x.rational, x.coefficient, y.rational, y.coefficient, radicand):
            if not isinstance(value, Fraction):
                raise InvalidVertexDiskError(f"vertex-disk parts must be exact Fractions, got {type(value).__name__}")
        if radicand < 0:
            raise InvalidVertexDiskError(f"vertex-disk radicand must be non-negative, got {radicand}")
        return cls(x, y, radicand)

    def wire(self) -> list[tuple[str, str]]:
        """(x0, x1, y0, y1, radicand) as decimal numerator/denominator text."""
        parts = (self.x.rational, self.x.coefficient, self.y.rational, self.y.coefficient, self.radicand)
        return [(str(part.numerator), str(part.denominator)) for part in parts]
