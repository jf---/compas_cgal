"""Contract tests: exact dOmega* connectors (line and arc sub-pieces) in the coverage certificate.

An exact connector runs along dOmega* itself, so its outer flank lies on the
declared boundary, its support circles coincide with declared circles and its
end disks coincide with the target's vertex disks: every coincident support the
exact-target certificate used to reject with ``CoincidentTargetSupportError``.
They are now decided exactly (side probes, ``side_probe_2.h``) and arc
sub-pieces sweep annular sectors (``exact_connectors_2.h``).

Every witness is re-checked independently of the certificate:

- in T: the exact squared distance to Omega* (native centre-domain module) is at
  most r^2;
- outside the connectors: the connectors cover dOmega*, so they sweep exactly
  dOmega* (+) D_r, and a witness must lie farther than r from dOmega* (the
  native exact nearest-boundary query);
- outside every cut annulus: Fraction arithmetic.
"""

from __future__ import annotations

import dataclasses
import hashlib
import math
from collections.abc import Iterator
from collections.abc import Sequence
from fractions import Fraction
from functools import cache

import pytest

from benchmarks.held_centre_domain_native import native_nearest
from benchmarks.held_coverage_target import exact_coverage_target
from benchmarks.held_declared_boundary import DeclaredArc
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_exact_centre_domain import CentreDomain as ExactCentreDomain
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.errors import BoundaryChainError
from compas_cgal.coverage_certificate.errors import ConnectorBoundaryCoverageError
from compas_cgal.coverage_certificate.errors import CoverageCertificateNativeError
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorPiece
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorSweep
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.exact_target import OffsetCirclePiece
from compas_cgal.coverage_certificate.exact_target import OffsetLinePiece
from compas_cgal.coverage_certificate.exact_target import OneRootPoint
from compas_cgal.coverage_certificate.native_certificate import CertificateReport
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import OneRootMillimetre
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.witness import ExposedPieceKind
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

F = Fraction
R = F(1)
TOOL = ToolRadius.build(1.0)
Cut = tuple[float, float, float]
XY = tuple[Fraction, Fraction]


def _point(x: Fraction | int, y: Fraction | int) -> ExactPoint2[WorldXY]:
    return ExactPoint2[WorldXY].build(F(x), F(y))


class _Pocket:
    """Sequential builder of a declared CCW cycle."""

    def __init__(self, start: ExactPoint2[WorldXY]) -> None:
        self.cursor = start
        self.primitives: list[DeclaredLine | DeclaredArc] = []

    def line(self, quarter: int, length: Fraction | int) -> _Pocket:
        line = DeclaredLine.build(self.cursor, PythagoreanDirection.build(F(0), quarter), F(length), len(self.primitives))
        self.primitives.append(line)
        self.cursor = line.end
        return self

    def arc(self, end: ExactPoint2[WorldXY], t: Fraction, quarter: int, sense: int) -> _Pocket:
        arc = DeclaredArc.build(self.cursor, end, PythagoreanDirection.build(t, quarter), 1 if sense == 1 else -1, len(self.primitives))
        self.primitives.append(arc)
        self.cursor = arc.end
        return self

    def close(self, name: str) -> DeclaredBoundary:
        return DeclaredBoundary.build(name, tuple(self.primitives), hashlib.sha256(name.encode()).hexdigest())


Fixture = tuple[DeclaredBoundary, ExactCentreDomain, ExactCentreDomainTarget]


@cache
def square() -> Fixture:
    boundary = _Pocket(_point(0, 0)).line(0, 20).line(1, 20).line(2, 20).line(3, 20).close("square")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def convex_arc() -> Fixture:
    """20 x 20 square whose bottom edge is a CCW arc about (10, 7.5), rho = 12.5: one-root kinks on an offset arc."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(-1, 2), 0, 1).line(1, 20).line(2, 20).line(3, 20).close("convex_arc")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def concave_arc() -> Fixture:
    """20 x 20 square whose bottom edge is a CW arc about (10, -7.5), rho = 12.5, bulging into P."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(1, 2), 0, -1).line(1, 20).line(2, 20).line(3, 20).close("concave_arc")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def l_shape() -> Fixture:
    """L-shape with a reflex corner at (10, 10): a junction r-circle arc sub-piece (R = r, a pie sector)."""
    boundary = _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, 10).line(1, 10).line(2, 10).line(3, 20).close("l_shape")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def lens() -> Fixture:
    """Convex arcs top and bottom: both line pieces run between one-root kinks, so both caps of a line are vertex circles."""
    boundary = _Pocket(_point(0, 0)).arc(_point(20, 0), F(-1, 2), 0, 1).line(1, 20).arc(_point(0, 20), F(-1, 2), 2, 1).line(3, 20).close("lens")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


@cache
def disk() -> Fixture:
    """Disk pocket of radius 10 about (10, 10): dOmega* is ONE closed circle piece (source == target)."""
    boundary = _Pocket(_point(10, 0)).arc(_point(10, 20), F(0), 0, 1).arc(_point(10, 0), F(0), 2, 1).close("disk")
    domain = ExactCentreDomain.build(boundary, R)
    return boundary, domain, exact_coverage_target(boundary, domain)


FIXTURES = (square, convex_arc, concave_arc, l_shape)


@pytest.fixture(autouse=True, scope="module")
def _release_native_domains() -> Iterator[None]:
    """Drop the cached native Omega* handles before interpreter shutdown."""
    yield
    for fixture in (*FIXTURES, lens, disk):
        fixture.cache_clear()


# The polygon path domain only carries polygon connectors' sites; unused here.
UNUSED_PATH = CentreDomain.build([Point2[WorldXY].build(x, y) for x, y in ((5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0))])


def _sweeps(cuts: Sequence[Cut]) -> list[FullCircleSweep]:
    return [FullCircleSweep.build(Point2[WorldXY].build(x, y), GuideRadius.build(big_r)) for x, y, big_r in cuts]


def _grid(x0: float, y0: float, step: float, nx: int, ny: int, big_r: float = 1.0) -> list[Cut]:
    return [(x0 + step * i, y0 + step * j, big_r) for i in range(nx) for j in range(ny)]


def _whole_boundary(target: ExactCentreDomainTarget, split: bool = True) -> list[ExactConnectorSweep]:
    """One connector per dOmega* piece (or one for all), each sub-piece the whole piece."""
    pieces = [ExactConnectorPiece.build(target, i, piece.source, piece.target) for i, piece in enumerate(target.pieces)]
    if split:
        return [ExactConnectorSweep.build([piece]) for piece in pieces]
    return [ExactConnectorSweep.build(pieces)]


def _doubled_arcs(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    """One connector per piece, and a second, identical one on every arc piece (overlapping sub-pieces)."""
    connectors = []
    for i, piece in enumerate(target.pieces):
        connectors.append(ExactConnectorSweep.build([ExactConnectorPiece.build(target, i, piece.source, piece.target)]))
        if isinstance(piece, OffsetCirclePiece):
            connectors.append(ExactConnectorSweep.build([ExactConnectorPiece.build(target, i, piece.source, piece.target)]))
    return connectors


def _both_modes(target: ExactCentreDomainTarget, cuts: Sequence[Cut], connectors: Sequence[ExactConnectorSweep]) -> CertificateReport:
    lazy = run_certificate(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target, exact_connectors=connectors, filter_mode=FilterMode.LAZY, audit=True)
    exact = run_certificate(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target, exact_connectors=connectors, filter_mode=FilterMode.EXACT_ONLY, audit=True)
    assert lazy.canonical_text() == exact.canonical_text()
    return lazy


def _in_target(domain: ExactCentreDomain, point: XY) -> bool:
    return domain.squared_distance_to(_point(*point)).compare(R * R) <= 0


def _beyond_connectors(domain: ExactCentreDomain, point: XY) -> bool:
    return native_nearest(domain.native, _point(*point)).squared_distance.compare(R * R) > 0


def _cut(point: XY, cuts: Sequence[Cut]) -> bool:
    for x, y, big_r in cuts:
        d2 = (point[0] - F(x)) ** 2 + (point[1] - F(y)) ** 2
        if (max(F(big_r) - R, F(0))) ** 2 <= d2 <= (F(big_r) + R) ** 2:
            return True
    return False


def _check_witnesses(domain: ExactCentreDomain, report: CertificateReport, cuts: Sequence[Cut]) -> None:
    assert report.exposed
    for piece in report.exposed:
        point = (F(piece.uncovered[0]), F(piece.uncovered[1]))
        assert _in_target(domain, point), f"witness {piece} is not in T"
        assert _beyond_connectors(domain, point), f"witness {piece} is swept by a connector"
        assert not _cut(point, cuts), f"witness {piece} is cut"


# Interior grids: disk cuts of reach 2 every 2.5 mm, centres at least 2 mm from
# the walls, so only the connectors can cover the strip along the walls.
INTERIOR = {
    "square": _grid(2.0, 2.0, 2.5, 7, 7),
    "convex_arc": _grid(2.0, -2.0, 2.5, 7, 9),
    "concave_arc": _grid(2.0, 7.0, 2.5, 7, 5),
    "l_shape": _grid(2.0, 2.0, 2.5, 7, 3) + _grid(2.0, 9.5, 2.5, 3, 4),
}


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda f: f.__name__)
def test_exact_connectors_on_the_whole_boundary_decide_with_the_interior_grid(fixture: object) -> None:
    # The connector flanks and support circles lie on the declared boundary and
    # the end disks on the vertex circles: every one is a coincident target
    # support, decided by a side probe instead of raising.
    assert callable(fixture)
    boundary, domain, target = fixture()
    cuts = INTERIOR[boundary.case_name]
    report = _both_modes(target, cuts, _whole_boundary(target))
    assert report.counters["coincident_target_decisions"] > 0
    assert report.counters["connector_pieces"] > 0
    if report.exposed:
        _check_witnesses(domain, report, cuts)


def test_the_square_is_certified_covered_by_exact_connectors_and_the_interior_grid() -> None:
    _, _, target = square()
    certificate = CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=target, exact_connectors=_whole_boundary(target))
    assert certificate.exposed_edge_count == 0


def test_without_the_connectors_the_wall_strip_is_proven_uncovered() -> None:
    _, domain, target = square()
    cuts = INTERIOR["square"]
    report = run_certificate(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target)
    assert report.exposed
    for piece in report.exposed:
        point = (F(piece.uncovered[0]), F(piece.uncovered[1]))
        assert _in_target(domain, point)
        assert not _cut(point, cuts)


def test_arc_sub_pieces_are_checked_as_connector_arcs_and_a_gap_is_proven() -> None:
    # Remove the interior cuts near the convex arc: the arc connector's inner
    # support circle borders the gap and must be exposed, as CONNECTOR_ARC.
    _, domain, target = convex_arc()
    cuts = [c for c in INTERIOR["convex_arc"] if c[1] > 3.0]
    report = _both_modes(target, cuts, _whole_boundary(target))
    _check_witnesses(domain, report, cuts)
    assert ExposedPieceKind.CONNECTOR_ARC in {piece.piece_kind for piece in report.exposed}
    with pytest.raises(UncoveredPocketError):
        CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target, exact_connectors=_whole_boundary(target))


def test_one_connector_for_the_whole_loop_decides_like_one_per_piece() -> None:
    _, _, target = convex_arc()
    cuts = [c for c in INTERIOR["convex_arc"] if c[1] > 3.0]
    split = _both_modes(target, cuts, _whole_boundary(target))
    joined = _both_modes(target, cuts, _whole_boundary(target, split=False))
    assert [(p.piece_kind, p.uncovered) for p in split.exposed] == [(p.piece_kind, p.uncovered) for p in joined.exposed]


def test_overlapping_connectors_on_one_arc_report_each_edge_once() -> None:
    _, domain, target = convex_arc()
    cuts = [c for c in INTERIOR["convex_arc"] if c[1] > 3.0]
    single = _both_modes(target, cuts, _whole_boundary(target))
    doubled = _both_modes(target, cuts, _doubled_arcs(target))
    _check_witnesses(domain, doubled, cuts)
    assert len(doubled.exposed) == len(single.exposed)


def test_connectors_that_leave_part_of_the_boundary_unswept_are_a_named_error() -> None:
    _, _, target = square()
    connectors = _whole_boundary(target)[:-1]
    with pytest.raises(ConnectorBoundaryCoverageError):
        CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=target, exact_connectors=connectors)


# --- (b) fix round 1: the cover check, the chain check and the check paths ---


def _rational_point(x: Fraction, y: Fraction) -> OneRootPoint:
    return OneRootPoint.build(OneRootMillimetre(x, F(0)), OneRootMillimetre(y, F(0)), F(0))


def _approx(point: OneRootPoint) -> tuple[float, float]:
    root = math.sqrt(point.radicand)
    return float(point.x.rational) + float(point.x.coefficient) * root, float(point.y.rational) + float(point.y.coefficient) * root


def _on_line(piece: OffsetLinePiece, along: Fraction) -> OneRootPoint:
    """The rational point level n + along u of a line piece (u = (n_y, -n_x))."""
    nx, ny = piece.normal
    return _rational_point(piece.level * nx + along * ny, piece.level * ny - along * nx)


def _line_along(piece: OffsetLinePiece, point: OneRootPoint) -> float:
    x, y = _approx(point)
    return float(piece.normal[1]) * x - float(piece.normal[0]) * y


def _on_circle(piece: OffsetCirclePiece, angle: float) -> OneRootPoint:
    """A rational point of the piece's circle near `angle` (stereographic s = tan(angle / 2), rationalised)."""
    s = Fraction(math.tan(0.5 * angle)).limit_denominator(4096)
    cx, cy = piece.centre
    return _rational_point(cx + piece.radius * (1 - s * s) / (1 + s * s), cy + piece.radius * 2 * s / (1 + s * s))


def _angle(piece: OffsetCirclePiece, point: OneRootPoint) -> float:
    x, y = _approx(point)
    return math.atan2(y - float(piece.centre[1]), x - float(piece.centre[0]))


def _arc_point(piece: OffsetCirclePiece, fraction: float) -> OneRootPoint:
    """A rational point `fraction` of the way along an open arc piece, in its sense."""
    a0, a1 = _angle(piece, piece.source), _angle(piece, piece.target)
    sweep = (a1 - a0) % (2 * math.pi) if piece.sense > 0 else -((a0 - a1) % (2 * math.pi))
    return _on_circle(piece, a0 + fraction * sweep)


def _interior_point(piece: OffsetLinePiece | OffsetCirclePiece, fraction: float) -> OneRootPoint:
    if isinstance(piece, OffsetCirclePiece):
        return _arc_point(piece, fraction)
    t0, t1 = _line_along(piece, piece.source), _line_along(piece, piece.target)
    return _on_line(piece, Fraction(t0 + fraction * (t1 - t0)).limit_denominator(4096))


def _sub(target: ExactCentreDomainTarget, index: int, start: OneRootPoint, end: OneRootPoint, *, full: bool = False) -> ExactConnectorPiece:
    return ExactConnectorPiece.build(target, index, start, end, full=full)


def _split_boundary(target: ExactCentreDomainTarget, cut: float = 0.5) -> list[ExactConnectorSweep]:
    """Every piece split at a rational interior site into two sub-pieces, each its own connector, listed in reverse order."""
    connectors = []
    for i, piece in enumerate(target.pieces):
        site = _interior_point(piece, cut)
        connectors.append(ExactConnectorSweep.build([_sub(target, i, site, piece.target)]))
        connectors.append(ExactConnectorSweep.build([_sub(target, i, piece.source, site)]))
    return connectors[::-1]


def _exposures(report: CertificateReport) -> list[tuple[ExposedPieceKind, tuple[object, object]]]:
    return [(piece.piece_kind, piece.uncovered) for piece in report.exposed]


def _replace_piece(target: ExactCentreDomainTarget, index: int, connectors: list[ExactConnectorSweep], *pieces: ExactConnectorPiece) -> list[ExactConnectorSweep]:
    out = list(connectors)
    out[index] = ExactConnectorSweep.build(list(pieces))
    return out


GAP_CUTS = {
    # The convex arc's strip uncut: CONNECTOR_ARC exposures.
    "convex_arc": [c for c in INTERIOR["convex_arc"] if c[1] > 3.0],
    # The square's bottom strip uncut (y in (2, 2.5)): exact line-flank exposures.
    "square": [c for c in INTERIOR["square"] if c[1] > 3.0],
}


@pytest.mark.parametrize("fixture", (square, convex_arc), ids=lambda f: f.__name__)
def test_real_sub_piece_splits_at_interior_sites_decide_like_whole_pieces(fixture: object) -> None:
    # Sub-pieces that end at interior sites, overlap nothing, and arrive unsorted
    # sweep the same set: one composite region per piece either way.
    assert callable(fixture)
    boundary, domain, target = fixture()
    cuts = GAP_CUTS[boundary.case_name]
    whole = _both_modes(target, cuts, _whole_boundary(target))
    split = _both_modes(target, cuts, _split_boundary(target))
    _check_witnesses(domain, split, cuts)
    assert _exposures(split) == _exposures(whole)


def test_overlapping_sub_pieces_out_of_order_decide_like_whole_pieces() -> None:
    _, _, target = convex_arc()
    cuts = GAP_CUTS["convex_arc"]
    connectors = _split_boundary(target, 0.7) + _split_boundary(target, 0.3)
    assert _exposures(_both_modes(target, cuts, connectors)) == _exposures(_both_modes(target, cuts, _whole_boundary(target)))


@pytest.mark.parametrize(("fixture", "kind"), ((square, OffsetLinePiece), (convex_arc, OffsetCirclePiece)), ids=("line", "arc"))
def test_a_gap_strictly_inside_a_piece_is_a_named_error(fixture: object, kind: type) -> None:
    assert callable(fixture)
    _, _, target = fixture()
    index = next(i for i, piece in enumerate(target.pieces) if isinstance(piece, kind))
    piece = target.pieces[index]
    connectors = _replace_piece(
        target,
        index,
        _whole_boundary(target),
        _sub(target, index, piece.source, _interior_point(piece, 0.4)),
        _sub(target, index, _interior_point(piece, 0.6), piece.target),
    )
    with pytest.raises(ConnectorBoundaryCoverageError):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=target, exact_connectors=connectors)


DISK_CUTS = [c for c in _grid(2.5, 2.5, 2.5, 7, 7) if (c[0] - 10.0) ** 2 + (c[1] - 10.0) ** 2 <= 42.0]


def _disk_piece() -> tuple[ExactCentreDomainTarget, OffsetCirclePiece]:
    _, _, target = disk()
    assert len(target.pieces) == 1
    piece = target.pieces[0]
    assert isinstance(piece, OffsetCirclePiece)
    assert piece.source == piece.target
    return target, piece


def _passes_the_cover_check(target: ExactCentreDomainTarget, connectors: list[ExactConnectorSweep]) -> None:
    # The disk has no offset-line piece, so a run stops at the interior anchor,
    # which needs one; reaching that limit means the chain and cover checks,
    # which run first, accepted the connectors.
    with pytest.raises(CoverageCertificateNativeError, match="no offset-line piece to anchor"):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(DISK_CUTS), exact_target=target, exact_connectors=connectors)


def test_a_closed_circle_piece_is_covered_by_a_full_turn_and_by_a_split_wrapping_its_source() -> None:
    target, piece = _disk_piece()
    top, bottom = _on_circle(piece, 0.5 * math.pi), _on_circle(piece, -0.5 * math.pi)
    _passes_the_cover_check(target, [ExactConnectorSweep.build([_sub(target, 0, piece.source, piece.source, full=True)])])
    # Without `full`, source -> source is a zero-length sub-piece: it covers nothing.
    with pytest.raises(ConnectorBoundaryCoverageError):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(DISK_CUTS), exact_target=target, exact_connectors=_whole_boundary(target))
    # The second sub-piece runs from the bottom through the source back to the top.
    _passes_the_cover_check(target, [ExactConnectorSweep.build([_sub(target, 0, top, bottom)]), ExactConnectorSweep.build([_sub(target, 0, bottom, top)])])


def test_a_gap_on_a_closed_circle_piece_across_its_source_is_a_named_error() -> None:
    target, piece = _disk_piece()
    connectors = [
        ExactConnectorSweep.build([_sub(target, 0, _on_circle(piece, 0.5 * math.pi), _on_circle(piece, -0.5 * math.pi))]),
        ExactConnectorSweep.build([_sub(target, 0, _on_circle(piece, -0.3 * math.pi), _on_circle(piece, 0.3 * math.pi))]),
    ]
    with pytest.raises(ConnectorBoundaryCoverageError):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(DISK_CUTS), exact_target=target, exact_connectors=connectors)


def _invalid(fixture: object, make: object, match: str) -> None:
    assert callable(fixture) and callable(make)
    _, _, target = fixture()
    connectors = make(target)
    with pytest.raises(CoverageCertificateNativeError, match=match):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=target, exact_connectors=connectors)


def _first(target: ExactCentreDomainTarget, kind: type) -> int:
    return next(i for i, piece in enumerate(target.pieces) if isinstance(piece, kind))


def _with(target: ExactCentreDomainTarget, index: int, start: OneRootPoint, end: OneRootPoint, *, full: bool = False) -> list[ExactConnectorSweep]:
    return _replace_piece(target, index, _whole_boundary(target), _sub(target, index, start, end, full=full))


def _line_end_off(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetLinePiece)
    x, y = target.pieces[i].target.x.rational, target.pieces[i].target.y.rational
    return _with(target, i, target.pieces[i].source, _rational_point(x + target.pieces[i].normal[0], y + target.pieces[i].normal[1]))


def _line_backwards(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetLinePiece)
    return _with(target, i, target.pieces[i].target, target.pieces[i].source)


def _line_beyond(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetLinePiece)
    piece = target.pieces[i]
    return _with(target, i, piece.source, _on_line(piece, Fraction(_line_along(piece, piece.target) + 0.5).limit_denominator(64)))


def _line_full(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetLinePiece)
    return _with(target, i, target.pieces[i].source, target.pieces[i].source, full=True)


def _arc_end_off(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetCirclePiece)
    return _with(target, i, target.pieces[i].source, _rational_point(F(10), F(-7)))


def _arc_backwards(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetCirclePiece)
    return _with(target, i, _arc_point(target.pieces[i], 0.6), _arc_point(target.pieces[i], 0.4))


def _arc_beyond(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetCirclePiece)
    return _with(target, i, target.pieces[i].source, _arc_point(target.pieces[i], 1.1))


def _arc_full_on_open(target: ExactCentreDomainTarget) -> list[ExactConnectorSweep]:
    i = _first(target, OffsetCirclePiece)
    return _with(target, i, target.pieces[i].source, target.pieces[i].source, full=True)


@pytest.mark.parametrize(
    ("fixture", "make", "match"),
    (
        (square, _line_end_off, "not on its offset line"),
        (square, _line_backwards, "runs backwards or leaves its piece"),
        (square, _line_beyond, "runs backwards or leaves its piece"),
        (square, _line_full, "cannot be a full turn"),
        (convex_arc, _arc_end_off, "not on its offset circle"),
        (convex_arc, _arc_backwards, "arc sub-piece runs backwards"),
        (convex_arc, _arc_beyond, "arc sub-piece leaves its piece"),
        (convex_arc, _arc_full_on_open, "full-turn sub-piece must start and end"),
    ),
    ids=lambda v: getattr(v, "__name__", v) if callable(v) else None,
)
def test_each_sub_piece_validation_error_is_raised(fixture: object, make: object, match: str) -> None:
    _invalid(fixture, make, match)


def test_an_open_boundary_chain_is_a_named_error() -> None:
    _, _, target = square()
    opened = dataclasses.replace(target, pieces=target.pieces[:-1])
    with pytest.raises(BoundaryChainError, match="ends where no piece starts"):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(INTERIOR["square"]), exact_target=opened, exact_connectors=_whole_boundary(opened))


def test_a_reflex_junction_is_a_named_error() -> None:
    # The L-shape's reflex corner is rounded by a junction arc in Omega*; replace
    # the line-arc-line run by two lines meeting at the corner (9, 9): a right
    # turn, which an erosion never has.
    _, _, target = l_shape()
    before, arc, after = target.pieces[0], target.pieces[1], target.pieces[2]
    assert isinstance(before, OffsetLinePiece) and isinstance(arc, OffsetCirclePiece) and isinstance(after, OffsetLinePiece)
    corner = _rational_point(F(9), F(9))
    reflex = (
        OffsetLinePiece.build(before.normal, before.level, before.source, corner),
        OffsetLinePiece.build(after.normal, after.level, corner, after.target),
        *target.pieces[3:],
    )
    bent = dataclasses.replace(target, pieces=reflex)
    with pytest.raises(BoundaryChainError, match="is reflex"):
        run_certificate(UNUSED_PATH, TOOL, _sweeps(INTERIOR["l_shape"]), exact_target=bent, exact_connectors=_whole_boundary(bent))


def test_an_exposed_exact_line_flank_is_reported_as_a_connector_flank() -> None:
    # The square's bottom strip uncut: the bottom connector's inner flank y = 2
    # borders the gap and must be exposed, as CAPSULE_FLANK of an exact connector.
    _, domain, target = square()
    cuts = GAP_CUTS["square"]
    report = _both_modes(target, cuts, _whole_boundary(target))
    _check_witnesses(domain, report, cuts)
    flanks = [p for p in report.exposed if p.piece_kind is ExposedPieceKind.CAPSULE_FLANK]
    assert flanks
    assert all(len(cuts) <= p.region_index < len(cuts) + len(target.pieces) for p in flanks)
    assert any(2 < F(p.uncovered[1]) < F(5, 2) for p in flanks)


def test_a_gap_at_a_convex_corner_is_found_through_flank_edges() -> None:
    # The cap-circle theorem's configuration: the cuts near the convex vertex
    # (1, 1) removed. No cap circle is checked; the uncovered corner must still be
    # proven through the flank and cut edges around it.
    _, domain, target = square()
    cuts = [c for c in INTERIOR["square"] if not (c[0] < 5.0 and c[1] < 5.0)]
    report = _both_modes(target, cuts, _whole_boundary(target))
    _check_witnesses(domain, report, cuts)
    assert any(F(p.uncovered[0]) < 5 and F(p.uncovered[1]) < 5 for p in report.exposed)
    with pytest.raises(UncoveredPocketError):
        CoverageCertificate.build(UNUSED_PATH, TOOL, _sweeps(cuts), exact_target=target, exact_connectors=_whole_boundary(target))


LENS_GRID = _grid(2.0, -2.0, 2.5, 7, 11)
# Exposure digests: sha256 over (kind, exact witness) in report order, for
# uncovered configurations of every fixture shape. The certificate is exact and
# deterministic, so any change to the edge arrangement (a dropped cut, a
# misordered event) changes a digest; LAZY == EXACT_ONLY cannot see a wrong rule
# that both modes share (the fast paths), this can.
EXPOSURE_CASES = {
    "square_bottom_gap": (square, GAP_CUTS["square"]),
    "square_corner_gap": (square, [c for c in INTERIOR["square"] if not (c[0] < 5.0 and c[1] < 5.0)]),
    "convex_arc_gap": (convex_arc, GAP_CUTS["convex_arc"]),
    "lens_left_gap": (lens, [c for c in LENS_GRID if c[0] > 3.0]),
}
EXPOSURE_DIGESTS = {
    "convex_arc_gap": "2377c32f42bd6e80",  # 10 exposures
    "lens_left_gap": "dea0e1e1c902f878",  # 20 exposures
    "square_bottom_gap": "8d0091bf931cafc7",  # 10 exposures
    "square_corner_gap": "2a532d7562788eda",  # 6 exposures
}


def _digest(report: CertificateReport) -> str:
    text = repr([(p.piece_kind.name, str(p.uncovered[0]), str(p.uncovered[1])) for p in report.exposed])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


@pytest.mark.parametrize("case", sorted(EXPOSURE_CASES))
def test_exposures_of_the_uncovered_configurations_are_pinned(case: str) -> None:
    fixture, cuts = EXPOSURE_CASES[case]
    _, domain, target = fixture()
    report = _both_modes(target, cuts, _whole_boundary(target))
    _check_witnesses(domain, report, cuts)
    assert _digest(report) == EXPOSURE_DIGESTS[case]
