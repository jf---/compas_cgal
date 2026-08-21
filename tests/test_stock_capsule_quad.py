"""Exact contract of the quad capsule depletion (`Stock2.subtract_capsule_quad`).

The quad path removes ``disk(A, r) U disk(B, r) U rect(A, B, h)`` in place of the
disk chain's hundreds of arcs. Two properties define it, and both are decided
here EXACTLY, over `fractions.Fraction`, with no tolerance anywhere in the
verdict:

SAFETY (the direction that matters)
    Every point the quad path removes lies within *radius* of the segment. The
    removed region is a SUBSET of the true swept capsule, so the model never
    reports material gone that the tool did not in fact reach.

COVERAGE (the budget)
    Every point within ``radius * (1 - CHAIN_SLACK_FRACTION)`` of the segment IS
    removed. That pins the under-coverage to the same budget the disk chain
    documents, rather than leaving it unstated.

Sample points are GENERATED in doubles -- any point will do -- but which
property a point must satisfy is decided by its EXACT rational distance to the
segment, never by the double that generated it.
"""

import math
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _stock_2

# The construction-density budget `stock_2.cpp` spends on the quad's half-width.
# NOT a decision tolerance: it is the documented under-coverage of the removed
# region, mirrored here as the exact rational the backend forms from the same
# double. Kept in lockstep with CHAIN_SLACK_FRACTION in src/stock_2.cpp.
CHAIN_SLACK_FRACTION = 1e-4

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)

# Grid resolution of the breadth sweep over the stock square. 81 x 81 = 6561
# point locations per segment; the discriminating samples are the targeted ones
# below, this sweep is here to catch a gross error anywhere else in the square.
GRID_SAMPLES_PER_AXIS = 81

# Targeted sampling of the two surfaces the properties are stated on: stations
# along the segment, and rays around each end cap.
BAND_STATIONS = 41
CAP_DIRECTIONS = 36

# Non-vacuity floor. A property test that never reaches its assertion passes for
# the wrong reason, so each sweep asserts it actually saw material removed / band
# points. The scarcest segment in SEGMENTS yields 446 removed and 422 in-band
# samples; this floor sits under both and would trip long before a sampler or
# geometry change hollowed a case out.
MIN_DISCRIMINATING_SAMPLES = 400

# (x0, y0, x1, y1, radius) -- axis-aligned, diagonal, shallow, steep, very
# short, long, boundary-crossing, and both travel directions of each.
SEGMENTS = [
    (2.0, 5.0, 8.0, 5.0, 0.5),  # horizontal
    (8.0, 5.0, 2.0, 5.0, 0.5),  # horizontal, reversed
    (5.0, 2.0, 5.0, 8.0, 0.75),  # vertical
    (5.0, 8.0, 5.0, 2.0, 0.75),  # vertical, reversed
    (2.0, 2.0, 8.0, 8.0, 1.0),  # 45 degrees
    (8.0, 2.0, 2.0, 8.0, 1.0),  # -45 degrees
    (1.0, 4.0, 9.0, 4.6, 0.6),  # shallow
    (9.0, 4.0, 1.0, 4.6, 0.6),  # shallow, reversed
    (4.0, 1.0, 4.6, 9.0, 0.6),  # steep
    (4.6, 9.0, 4.0, 1.0, 0.6),  # steep, reversed
    (5.0, 5.0, 5.001, 5.0, 0.4),  # very short, horizontal
    (5.0, 5.0, 5.0, 5.001, 0.4),  # very short, vertical
    (5.0, 5.0, 5.0007, 5.0003, 0.4),  # very short, oblique
    (5.0, 5.0, 5.0 + 1e-9, 5.0 + 1e-9, 0.4),  # near-degenerate, huge r/len ratio
    (0.5, 0.5, 9.5, 9.5, 1.5),  # long diagonal
    (-2.0, 5.0, 12.0, 5.0, 1.25),  # crosses both boundaries
    (5.0, -2.0, 5.0, 12.0, 1.25),  # crosses both boundaries, vertical
    (-1.0, -1.0, 4.0, 4.0, 0.9),  # enters through a corner
    (7.0, 11.0, 11.0, 7.0, 2.0),  # mostly outside, clips the corner
    (0.0, 0.0, 10.0, 3.0, 0.8),  # starts on a boundary corner
    (3.0, 0.0, 7.0, 0.0, 0.7),  # runs along a boundary edge
    (2.0, 3.0, 3.0, 2.0, 0.25),  # short anti-diagonal, small tool
    (1.3, 8.7, 8.9, 1.1, 2.0),  # long, large tool
]


def exact_squared_distance_to_segment(px, py, x0, y0, x1, y1):
    """Exact squared distance from a point to a segment, over the rationals.

    Every input is a binary64 value and therefore IS a rational, so
    `fractions.Fraction` carries the whole computation without rounding: the
    projection parameter, its clamp to the segment, and the squared distance are
    all exact. This is the oracle the two properties are decided against.

    Args:
        px: X coordinate of the query point.
        py: Y coordinate of the query point.
        x0: X coordinate of the segment start.
        y0: Y coordinate of the segment start.
        x1: X coordinate of the segment end.
        y1: Y coordinate of the segment end.

    Returns:
        The squared distance as an exact `Fraction`.
    """
    p_x, p_y = Fraction(px), Fraction(py)
    a_x, a_y = Fraction(x0), Fraction(y0)
    d_x, d_y = Fraction(x1) - a_x, Fraction(y1) - a_y
    length_squared = d_x * d_x + d_y * d_y
    if length_squared == 0:
        return (p_x - a_x) ** 2 + (p_y - a_y) ** 2
    t = ((p_x - a_x) * d_x + (p_y - a_y) * d_y) / length_squared
    t = min(max(t, Fraction(0)), Fraction(1))
    foot_x, foot_y = a_x + t * d_x, a_y + t * d_y
    return (p_x - foot_x) ** 2 + (p_y - foot_y) ** 2


def slack_half_width(radius):
    """The exact half-width the quad is certified to cover: ``r * (1 - f)``.

    Formed exactly as the backend forms it -- the product of the injected radius
    with ``1 - CHAIN_SLACK_FRACTION`` as rationals -- so the test threshold and
    the certified bound are the same number, not two roundings of one intent.

    Args:
        radius: Tool radius as a binary64 value.

    Returns:
        The exact `Fraction` lower bound on the removed rectangle's half-width.
    """
    return Fraction(radius) * (1 - Fraction(CHAIN_SLACK_FRACTION))


def sample_points(x0, y0, x1, y1, radius):
    """Generate query points: a grid over the square plus targeted extremals.

    The targeted points sit on the two surfaces the properties are stated on --
    the band at ``+/- (1 - f) * r`` alongside the segment, the capsule wall at
    ``+/- r``, and rays around both end caps -- because a blind grid almost never
    lands where under- and over-coverage would actually show.

    Args:
        x0: X coordinate of the segment start.
        y0: Y coordinate of the segment start.
        x1: X coordinate of the segment end.
        y1: Y coordinate of the segment end.
        radius: Tool radius.

    Returns:
        A list of ``(x, y)`` double pairs.
    """
    axis = np.linspace(0.0, 10.0, GRID_SAMPLES_PER_AXIS)
    points = [(float(x), float(y)) for x in axis for y in axis]

    length = math.hypot(x1 - x0, y1 - y0)
    if length > 0.0:
        ux, uy = (x1 - x0) / length, (y1 - y0) / length
        nx, ny = -uy, ux
        offsets = [f * float(slack_half_width(radius)) for f in (-1.0, -0.5, 0.0, 0.5, 1.0)]
        offsets += [f * radius for f in (-1.0, -0.999, 0.999, 1.0)]
        for i in range(BAND_STATIONS):
            t = i / (BAND_STATIONS - 1)
            cx, cy = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
            points += [(cx + s * nx, cy + s * ny) for s in offsets]

    radii = [float(slack_half_width(radius)), radius, 0.999 * radius, 1.001 * radius]
    for i in range(CAP_DIRECTIONS):
        angle = 2.0 * math.pi * i / CAP_DIRECTIONS
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        for rho in radii:
            points.append((x0 + rho * cos_a, y0 + rho * sin_a))
            points.append((x1 + rho * cos_a, y1 + rho * sin_a))
    return points


@pytest.mark.parametrize("segment", SEGMENTS)
def test_capsule_quad_removes_nothing_outside_the_true_capsule(segment):
    """SAFETY: every removed point lies within *radius* of the segment.

    Zero violations permitted -- this is the direction in which an error would
    hand the caller cleared material the tool never reached.
    """
    x0, y0, x1, y1, radius = segment
    before = _stock_2.Stock2(SQUARE, [])
    after = _stock_2.Stock2(SQUARE, [])
    after.subtract_capsule_quad(x0, y0, x1, y1, radius)

    radius_squared = Fraction(radius) ** 2
    removed = 0
    violations = []
    for px, py in sample_points(*segment):
        if after.contains(px, py) or not before.contains(px, py):
            continue
        removed += 1
        if exact_squared_distance_to_segment(px, py, x0, y0, x1, y1) > radius_squared:
            violations.append((px, py))
    assert violations == []
    assert removed >= MIN_DISCRIMINATING_SAMPLES


@pytest.mark.parametrize("segment", SEGMENTS)
def test_capsule_quad_covers_the_documented_slack_band(segment):
    """COVERAGE: everything within ``r * (1 - f)`` of the segment is removed.

    Pins the under-coverage to the documented budget: material may survive only
    in the outermost `CHAIN_SLACK_FRACTION` of the tool radius.
    """
    x0, y0, x1, y1, radius = segment
    before = _stock_2.Stock2(SQUARE, [])
    after = _stock_2.Stock2(SQUARE, [])
    after.subtract_capsule_quad(x0, y0, x1, y1, radius)

    band_squared = slack_half_width(radius) ** 2
    in_band = 0
    survivors = []
    for px, py in sample_points(*segment):
        if not before.contains(px, py):
            continue
        if exact_squared_distance_to_segment(px, py, x0, y0, x1, y1) > band_squared:
            continue
        in_band += 1
        if after.contains(px, py):
            survivors.append((px, py))
    assert survivors == []
    assert in_band >= MIN_DISCRIMINATING_SAMPLES


@pytest.mark.parametrize("segment", SEGMENTS)
def test_capsule_quad_keeps_the_boolean_set_canonical(segment):
    """The removal leaves a `General_polygon_set_2` that still satisfies its own
    representation invariant -- orthogonal to point-set equality, and silently
    degrading every later query if it fails."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule_quad(*segment)
    assert stock.representation_is_valid()


def test_capsule_quad_arrangement_stays_orders_below_the_chain():
    """The point of the quad: six curves, not hundreds of arcs.

    Observed on this 10x10 square, one capsule from (2,5) to (8,5):

        radius 0.5  ->  quad 10 vertices, chain 1206 vertices  (120x)
        radius 2.0  ->  quad 10 vertices, chain  306 vertices  ( 30x)

    The asserted floors are well under both, so the test measures the property
    (bounded, radius-independent growth per capsule) rather than the numbers.
    """
    for radius in (0.5, 2.0):
        quad = _stock_2.Stock2(SQUARE, [])
        quad.subtract_capsule_quad(2.0, 5.0, 8.0, 5.0, radius)
        chain = _stock_2.Stock2(SQUARE, [])
        chain.subtract_capsule(2.0, 5.0, 8.0, 5.0, radius)

        quad_vertices = quad.arrangement_stats()[0]
        chain_vertices = chain.arrangement_stats()[0]
        assert quad_vertices <= 12
        assert chain_vertices >= 20 * quad_vertices


def test_capsule_quad_growth_is_independent_of_segment_length():
    """A long capsule costs the arrangement no more than a short one.

    The chain's vertex count grows with the segment length (one disk per spacing
    step); the quad's does not, which is what removes the chain-length term from
    the depletion cost.
    """
    short = _stock_2.Stock2(SQUARE, [])
    short.subtract_capsule_quad(4.0, 5.0, 4.5, 5.0, 0.5)
    long = _stock_2.Stock2(SQUARE, [])
    long.subtract_capsule_quad(1.0, 5.0, 9.0, 5.0, 0.5)
    assert short.arrangement_stats()[0] == long.arrangement_stats()[0]


def test_capsule_quad_degenerate_segment_is_exactly_the_disk():
    """A zero-length motion sweeps the plain disk, decided by exact equality."""
    quad = _stock_2.Stock2(SQUARE, [])
    quad.subtract_capsule_quad(5.0, 5.0, 5.0, 5.0, 1.0)
    disk = _stock_2.Stock2(SQUARE, [])
    disk.subtract_disk(5.0, 5.0, 1.0)
    assert quad.exactly_equals(disk)


def test_capsule_quad_matches_the_chain_where_both_are_certain():
    """Both paths under-cover, so neither region contains the other -- but both
    must agree on every point the certified band forces in and on every point
    outside the capsule."""
    quad = _stock_2.Stock2(SQUARE, [])
    chain = _stock_2.Stock2(SQUARE, [])
    quad.subtract_capsule_quad(2.0, 3.0, 7.0, 6.0, 1.0)
    chain.subtract_capsule(2.0, 3.0, 7.0, 6.0, 1.0)

    band_squared = slack_half_width(1.0) ** 2
    radius_squared = Fraction(1.0) ** 2
    inside_band_survivors = []
    outside_capsule_disagreements = []
    for px, py in sample_points(2.0, 3.0, 7.0, 6.0, 1.0):
        distance_squared = exact_squared_distance_to_segment(px, py, 2.0, 3.0, 7.0, 6.0)
        if distance_squared <= band_squared:
            if quad.contains(px, py) or chain.contains(px, py):
                inside_band_survivors.append((px, py))
        elif distance_squared > radius_squared:
            if quad.contains(px, py) != chain.contains(px, py):
                outside_capsule_disagreements.append((px, py))
    assert inside_band_survivors == []
    assert outside_capsule_disagreements == []


def test_capsule_quad_rejects_non_positive_radius():
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(ValueError, match="radius should be positive"):
        stock.subtract_capsule_quad(2.0, 5.0, 8.0, 5.0, 0.0)
    with pytest.raises(ValueError, match="radius should be positive"):
        stock.subtract_capsule_quad(2.0, 5.0, 8.0, 5.0, -1.0)


@pytest.mark.parametrize(
    "arguments",
    [
        (math.nan, 5.0, 8.0, 5.0, 0.5),
        (2.0, math.inf, 8.0, 5.0, 0.5),
        (2.0, 5.0, -math.inf, 5.0, 0.5),
        (2.0, 5.0, 8.0, math.nan, 0.5),
        (2.0, 5.0, 8.0, 5.0, math.inf),
    ],
)
def test_capsule_quad_rejects_non_finite_input(arguments):
    """A non-finite double has no exact rational to inject, so it is refused at
    the boundary -- by its OWN named exception, a ValueError subclass."""
    stock = _stock_2.Stock2(SQUARE, [])
    with pytest.raises(_stock_2.NonFiniteCapsuleInputError):
        stock.subtract_capsule_quad(*arguments)
    assert issubclass(_stock_2.NonFiniteCapsuleInputError, ValueError)


def test_capsule_quad_certificate_error_is_not_an_argument_fault():
    """The half-width certificate failing is a broken internal invariant, so it
    must never be catchable as malformed input."""
    assert issubclass(_stock_2.CapsuleQuadCertificateError, RuntimeError)
    assert not issubclass(_stock_2.CapsuleQuadCertificateError, ValueError)


def test_capsule_quad_refuses_an_unrepresentable_radius_to_length_ratio():
    """A radius:length ratio that overflows the half-width scale fails loudly.

    The endpoints differ (so this is not the disk case) but by a subnormal step,
    while the radius is astronomically large: the double scale factor overflows
    to infinity and there is no rational to inject. The construction refuses
    rather than removing a region it cannot certify.
    """
    stock = _stock_2.Stock2(SQUARE, [])
    assert 5e-324 != 0.0  # the endpoints really are two different doubles
    with pytest.raises(_stock_2.CapsuleQuadCertificateError, match="representable range"):
        stock.subtract_capsule_quad(0.0, 0.0, 5e-324, 0.0, 1e300)
