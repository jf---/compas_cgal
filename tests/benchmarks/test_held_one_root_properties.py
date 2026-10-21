"""Hypothesis properties of :mod:`benchmarks.held_one_root` against an independent certified oracle.

The oracle never calls the module under test. It keeps its own raw (unfolded)
``a + b√c`` triples and evaluates sums ``A + B√α + C√β + D√(αβ)`` by rigorous
rational interval enclosures: ``m = isqrt(p·q·4^K)`` gives
``m ≤ √(p/q)·q·2^K < m + 1``, so ``√(p/q) ∈ [m, m + 1] / (q·2^K)``, and interval
arithmetic over :class:`~fractions.Fraction` adds no rounding. That enclosure is
the proven error bound. A draw whose enclosure straddles zero (exact ties such
as ``√8`` vs ``2√2``, and nonzero near-ties below ``2^-K``) is never discarded:
it is decided by exact squaring in this file, grouping by ``√α`` (the module
groups by ``√β``). Angular order on a circle is checked against the diamond
angle, a strictly increasing reparametrisation of the polar angle onto
``[0, 4)`` that is rational in the coordinates, so the swept-angle order is the
true angle order with no transcendental evaluation.

Inputs come in two families: generic binary64 values injected exactly
(``Fraction(float)``), and small integers with radicands ``s·k²/m²`` that make
equal values appear under different radicands, collinear points across roots,
coincident and antipodal points on a circle.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal
from typing import cast

import pytest
from hypothesis import event
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_boundary import Sense
from benchmarks.held_one_root import CrossRootArithmeticError
from benchmarks.held_one_root import OneRootNumber
from benchmarks.held_one_root import OneRootPoint2
from benchmarks.held_one_root import compare_on_arc
from benchmarks.held_one_root import orientation
from compas_cgal.units import WorldXY

# Bits of each √ enclosure: its absolute width is 2^-256 / q (q the radicand's
# denominator), far below every generic-double gap yet above the 1e-300 offsets
# drawn by the near-tie strategies, so both oracle routes are exercised.
ENCLOSURE_BITS = 256
# Refinement cap for the angle oracle (bits). Termination is guaranteed for
# distinct angles; the cap only turns a would-be infinite loop into a loud error.
MAX_ENCLOSURE_BITS = 1 << 16
# Bits of the lower bound L <= √R² used to place generic chords |h| < L inside a circle.
CHORD_BOUND_BITS = 64
EXAMPLES = 300

Interval = tuple[Fraction, Fraction]
Route = Literal["enclosure", "exact squaring"]


class OracleFormError(Exception):
    """Two oracle forms over different radicand pairs were combined."""


class OracleRefinementExhaustedError(Exception):
    """The angle oracle hit MAX_ENCLOSURE_BITS without separating two values."""


# ---------------------------------------------------------------------------
# Oracle: rigorous enclosures + exact squaring fallback
# ---------------------------------------------------------------------------


def _sgn(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def _sqrt_enclosure(radicand: Fraction, bits: int) -> Interval:
    """Rigorous ``[lo, hi] ∋ √radicand`` (``radicand >= 0``)."""
    p, q = radicand.numerator, radicand.denominator
    root = math.isqrt((p * q) << (2 * bits))
    scale = q << bits
    return Fraction(root, scale), Fraction(root + 1, scale)


def _scaled(factor: Fraction, interval: Interval) -> Interval:
    lo, hi = factor * interval[0], factor * interval[1]
    return (lo, hi) if factor >= 0 else (hi, lo)


def _sum(*intervals: Interval) -> Interval:
    return sum((i[0] for i in intervals), Fraction(0)), sum((i[1] for i in intervals), Fraction(0))


@dataclass(frozen=True)
class Radical4:
    """The oracle's form ``A + B√α + C√β + D√(αβ)``, all rational, ``α, β >= 0``."""

    a: Fraction
    b: Fraction
    c: Fraction
    d: Fraction
    alpha: Fraction
    beta: Fraction

    def enclosure(self, bits: int) -> Interval:
        return _sum(
            (self.a, self.a),
            _scaled(self.b, _sqrt_enclosure(self.alpha, bits)),
            _scaled(self.c, _sqrt_enclosure(self.beta, bits)),
            _scaled(self.d, _sqrt_enclosure(self.alpha * self.beta, bits)),
        )

    def minus(self, other: Radical4) -> Radical4:
        if (self.alpha, self.beta) != (other.alpha, other.beta):
            raise OracleFormError(f"(√{self.alpha}, √{self.beta}) vs (√{other.alpha}, √{other.beta})")
        return Radical4(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d, self.alpha, self.beta)

    def plus(self, other: Radical4) -> Radical4:
        return self.minus(Radical4(-other.a, -other.b, -other.c, -other.d, other.alpha, other.beta))


def _exact_sign2(p: Fraction, q: Fraction, r: Fraction) -> int:
    """Exact ``sign(p + q√r)``, ``r >= 0``, by squaring."""
    sign_p, sign_q = _sgn(p), (_sgn(q) if r != 0 else 0)
    if sign_q == 0:
        return sign_p
    if sign_p == 0:
        return sign_q
    if sign_p == sign_q:
        return sign_p
    return sign_p * _sgn(p * p - q * q * r)


def _exact_sign4(f: Radical4) -> int:
    """Exact sign grouping by ``√α``: ``(A + C√β) + √α·(B + D√β)``."""
    sign_u = _exact_sign2(f.a, f.c, f.beta)
    sign_w = _exact_sign2(f.b, f.d, f.beta) if f.alpha != 0 else 0
    if sign_w == 0:
        return sign_u
    if sign_u == 0:
        return sign_w
    if sign_u == sign_w:
        return sign_u
    # (A + C√β)² − α(B + D√β)² = (A² + C²β − αB² − αD²β) + 2(AC − αBD)√β
    rational = f.a * f.a + f.c * f.c * f.beta - f.alpha * (f.b * f.b + f.d * f.d * f.beta)
    return sign_u * _exact_sign2(rational, 2 * (f.a * f.c - f.alpha * f.b * f.d), f.beta)


@dataclass(frozen=True)
class Decision:
    sign: int
    route: Route


def _decide(f: Radical4, bits: int = ENCLOSURE_BITS) -> Decision:
    lo, hi = f.enclosure(bits)
    if lo > 0:
        return Decision(1, "enclosure")
    if hi < 0:
        return Decision(-1, "enclosure")
    return Decision(_exact_sign4(f), "exact squaring")


def _oracle_sign(f: Radical4, bits: int = ENCLOSURE_BITS) -> int:
    decision = _decide(f, bits)
    event(f"oracle route: {decision.route}")
    return decision.sign


@dataclass(frozen=True)
class Surd:
    """A raw drawn ``a + b√c`` (unfolded; the module folds its own copy)."""

    a: Fraction
    b: Fraction
    c: Fraction

    def number(self) -> OneRootNumber:
        return OneRootNumber.build(self.a, self.b, self.c)

    def shifted(self, offset: Fraction) -> Surd:
        return Surd(self.a + offset, self.b, self.c)

    def rescaled(self, k: int, m: int) -> Surd:
        """Same value under radicand ``c·k²/m²`` (``√(c·k²/m²) = (k/m)√c``)."""
        return Surd(self.a, self.b * m / k, self.c * k * k / (m * m))

    def form(self) -> Radical4:
        return Radical4(self.a, self.b, Fraction(0), Fraction(0), self.c, Fraction(0))


def _difference(x: Surd, y: Surd) -> Radical4:
    return Radical4(x.a - y.a, x.b, -y.b, Fraction(0), x.c, y.c)


def _product(x: Surd, y: Surd) -> Radical4:
    return Radical4(x.a * y.a, x.b * y.a, x.a * y.b, x.b * y.b, x.c, y.c)


def _minus_folded(f: Radical4, folded: OneRootNumber) -> Radical4:
    """``f − folded`` where ``folded``'s radicand is 0 or ``f.alpha``."""
    if folded.radicand not in (0, f.alpha):
        raise OracleFormError(f"√{folded.radicand} is not in ℚ(√{f.alpha})")
    return Radical4(f.a - folded.rational, f.b - folded.coefficient, f.c, f.d, f.alpha, f.beta)


@dataclass(frozen=True)
class RawPoint:
    """A drawn point; both coordinates carry the same raw radicand."""

    x: Surd
    y: Surd

    def point(self) -> OneRootPoint2:
        return OneRootPoint2.build(self.x.number(), self.y.number())

    def rescaled(self, k: int, m: int) -> RawPoint:
        return RawPoint(self.x.rescaled(k, m), self.y.rescaled(k, m))

    def reflected(self, centre: tuple[Fraction, Fraction]) -> RawPoint:
        return RawPoint(Surd(2 * centre[0] - self.x.a, -self.x.b, self.x.c), Surd(2 * centre[1] - self.y.a, -self.y.b, self.y.c))


def _cross(centre: tuple[Fraction, Fraction], first: RawPoint, second: RawPoint) -> Radical4:
    fx, fy = first.x.shifted(-centre[0]), first.y.shifted(-centre[1])
    sx, sy = second.x.shifted(-centre[0]), second.y.shifted(-centre[1])
    return _product(fx, sy).minus(_product(fy, sx))


def _dot(centre: tuple[Fraction, Fraction], first: RawPoint, second: RawPoint) -> Radical4:
    fx, fy = first.x.shifted(-centre[0]), first.y.shifted(-centre[1])
    sx, sy = second.x.shifted(-centre[0]), second.y.shifted(-centre[1])
    return _product(fx, sx).plus(_product(fy, sy))


def _diamond(dx: Surd, dy: Surd, bits: int) -> Interval | None:
    """Enclosure of the diamond angle of ``(dx, dy) != 0`` in ``[0, 4)``; ``None`` if ``bits`` is too coarse."""
    sign_x, sign_y = _oracle_sign(dx.form(), bits), _oracle_sign(dy.form(), bits)
    x, y = dx.form().enclosure(bits), dy.form().enclosure(bits)
    neg_x, neg_y = (-x[1], -x[0]), (-y[1], -y[0])
    if sign_x > 0 and sign_y >= 0:
        quadrant, numerator, denominator = 0, y, _sum(x, y)
    elif sign_x <= 0 and sign_y > 0:
        quadrant, numerator, denominator = 1, neg_x, _sum(neg_x, y)
    elif sign_x < 0 and sign_y <= 0:
        quadrant, numerator, denominator = 2, neg_y, _sum(neg_x, neg_y)
    else:
        quadrant, numerator, denominator = 3, x, _sum(x, neg_y)
    if denominator[0] <= 0:
        return None
    # The true numerator is >= 0 and the true ratio is < 1: clamping is sound.
    lo = max(Fraction(0), numerator[0]) / denominator[1]
    hi = min(Fraction(1), numerator[1] / denominator[0])
    return quadrant + lo, quadrant + hi


def _swept(centre: tuple[Fraction, Fraction], sense: int, start: RawPoint, point: RawPoint, bits: int) -> Interval | None:
    """Enclosure of the swept diamond angle from ``start`` to ``point`` in ``sense``, in ``[0, 4)``."""
    side = sense * _oracle_sign(_cross(centre, start, point), bits)
    if side == 0:  # exactly on the start ray or its opposite
        return (Fraction(0), Fraction(0)) if _oracle_sign(_dot(centre, start, point), bits) > 0 else (Fraction(2), Fraction(2))
    begin = _diamond(start.x.shifted(-centre[0]), start.y.shifted(-centre[1]), bits)
    end = _diamond(point.x.shifted(-centre[0]), point.y.shifted(-centre[1]), bits)
    if begin is None or end is None:
        return None
    lo, hi = end[0] - begin[1], end[1] - begin[0]
    if sense < 0:
        lo, hi = -hi, -lo
    turns = math.floor(lo / 4)
    lo, hi = lo - 4 * turns, hi - 4 * turns
    if lo <= 0 or hi >= 4:  # touches the start ray, which ``side != 0`` excludes: refine
        return None
    return lo, hi


def _arc_order(centre: tuple[Fraction, Fraction], sense: int, start: RawPoint, first: RawPoint, second: RawPoint) -> int:
    bits = ENCLOSURE_BITS
    while bits <= MAX_ENCLOSURE_BITS:
        a, b = _swept(centre, sense, start, first, bits), _swept(centre, sense, start, second, bits)
        if a is not None and b is not None:
            if a[1] < b[0]:
                return -1
            if b[1] < a[0]:
                return 1
            if _oracle_sign(_cross(centre, first, second), bits) == 0 and _oracle_sign(_dot(centre, first, second), bits) > 0:
                event("arc oracle: same direction (exact)")
                return 0
        event("arc oracle: refined")
        bits *= 2
    raise OracleRefinementExhaustedError(f"{first} vs {second} from {start}")


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

_GENERIC = st.floats(min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False, allow_subnormal=False).map(Fraction)
_GENERIC_RADICAND = st.floats(min_value=1e-3, max_value=1e3, allow_nan=False, allow_infinity=False).map(Fraction)
_SMALL = st.integers(-6, 6).map(Fraction)
_SQUAREFREE = st.sampled_from([2, 3, 5, 6, 7, 10])
_SCALE = st.integers(1, 4)
_TINY = st.floats(min_value=1e-300, max_value=1e-250).map(Fraction)


@st.composite
def _small_radicands(draw: st.DrawFn) -> Fraction:
    """``s·k²/m²`` (equal surds under different radicands), squares (fold to rationals) and 0."""
    if draw(st.booleans()):
        return Fraction(draw(_SQUAREFREE) * draw(_SCALE) ** 2, draw(_SCALE) ** 2)
    return draw(st.sampled_from([Fraction(0), Fraction(1), Fraction(4), Fraction(9, 4)]))


_FAMILIES = {"generic": (_GENERIC, _GENERIC_RADICAND), "small": (_SMALL, _small_radicands())}
_FAMILY = st.sampled_from(sorted(_FAMILIES))


@st.composite
def surds(draw: st.DrawFn, radicand: Fraction | None = None) -> Surd:
    rational, radicands = _FAMILIES[draw(_FAMILY)]
    return Surd(draw(rational), draw(rational), draw(radicands) if radicand is None else radicand)


@st.composite
def comparable_pairs(draw: st.DrawFn) -> tuple[Surd, Surd]:
    """Independent, same-root, equal-under-rescaled-radicand and sub-enclosure near-tie pairs."""
    x = draw(surds())
    kind = draw(st.sampled_from(["independent", "same root", "equal rescaled", "near tie rational", "near tie radicand"]))
    event(f"pair: {kind}")
    if kind == "independent":
        return x, draw(surds())
    if kind == "same root":
        return x, draw(surds(radicand=x.c))
    rescaled = x.rescaled(draw(_SCALE), draw(_SCALE))
    if kind == "equal rescaled":
        return x, rescaled
    if kind == "near tie rational":
        return x, rescaled.shifted(draw(_TINY) * draw(st.sampled_from([1, -1])))
    return x, Surd(x.a, x.b, x.c + draw(_TINY))


@st.composite
def raw_points(draw: st.DrawFn, family: str) -> RawPoint:
    rational, radicands = _FAMILIES[family]
    radicand = draw(radicands)
    return RawPoint(Surd(draw(rational), draw(rational), radicand), Surd(draw(rational), draw(rational), radicand))


@st.composite
def orientation_triples(draw: st.DrawFn) -> tuple[tuple[Fraction, Fraction], RawPoint, RawPoint]:
    """Independent points, points collinear with the centre across roots, coincident and rescaled points."""
    family = draw(_FAMILY)
    rational, radicands = _FAMILIES[family]
    centre = (draw(rational), draw(rational))
    first = draw(raw_points(family))
    kind = draw(st.sampled_from(["independent", "collinear across roots", "same point", "rescaled point", "reflected point"]))
    event(f"triple: {kind}")
    if kind == "independent":
        return centre, first, draw(raw_points(family))
    if kind == "same point":
        return centre, first, first
    if kind == "rescaled point":
        return centre, first, first.rescaled(draw(_SCALE), draw(_SCALE))
    if kind == "reflected point":
        return centre, first, first.reflected(centre)
    u, v = draw(rational), draw(rational)
    t1, t2 = draw(surds()), draw(surds())

    def along(t: Surd) -> RawPoint:
        return RawPoint(Surd(centre[0] + u * t.a, u * t.b, t.c), Surd(centre[1] + v * t.a, v * t.b, t.c))

    return centre, along(t1), along(t2)


@dataclass(frozen=True)
class Circle:
    centre: tuple[Fraction, Fraction]
    radius_squared: Fraction
    family: str


@st.composite
def circles(draw: st.DrawFn) -> Circle:
    family = draw(_FAMILY)
    if family == "generic":
        centre = (draw(_GENERIC), draw(_GENERIC))
        radius_squared = Fraction(draw(st.floats(min_value=1e-2, max_value=1e6, allow_nan=False, allow_infinity=False)))
    else:
        centre = (draw(_SMALL), draw(_SMALL))
        radius_squared = Fraction(draw(st.integers(1, 50)))
    return Circle(centre, radius_squared, family)


@st.composite
def circle_points(draw: st.DrawFn, circle: Circle) -> RawPoint:
    """A point of ``circle`` on a rational vertical or horizontal chord: ``(h, ±√(R² − h²))``, radicand rescaled."""
    if circle.family == "generic":
        bound = _sqrt_enclosure(circle.radius_squared, CHORD_BOUND_BITS)[0]
        offset = Fraction(draw(st.floats(min_value=-1, max_value=1, exclude_min=True, exclude_max=True))) * bound
    else:
        limit = math.isqrt(int(circle.radius_squared) - 1)
        offset = Fraction(draw(st.integers(-limit, limit)))
    k, m = draw(_SCALE), draw(_SCALE)
    radicand = (circle.radius_squared - offset * offset) * k * k / (m * m)
    root = Surd(Fraction(0), Fraction(draw(st.sampled_from([1, -1])) * m, k), radicand)
    chord = Surd(offset, Fraction(0), radicand)
    cx, cy = circle.centre
    if draw(st.booleans()):
        return RawPoint(chord.shifted(cx), root.shifted(cy))
    return RawPoint(root.shifted(cx), chord.shifted(cy))


@st.composite
def arc_cases(draw: st.DrawFn) -> tuple[Circle, Sense, RawPoint, RawPoint, RawPoint]:
    """Start and two points of one circle, with ties: start itself, its antipode, coincident and rescaled copies."""
    circle = draw(circles())
    sense = cast(Sense, draw(st.sampled_from([1, -1])))
    start = draw(circle_points(circle))
    first_kind = draw(st.sampled_from(["fresh", "start", "antipode of start"]))
    second_kind = draw(st.sampled_from(["fresh", "start", "antipode of start", "first", "first rescaled", "antipode of first"]))
    event(f"arc: {first_kind} / {second_kind}")
    choices = {"start": start, "antipode of start": start.reflected(circle.centre)}
    first = draw(circle_points(circle)) if first_kind == "fresh" else choices[first_kind]
    choices |= {"first": first, "first rescaled": first.rescaled(draw(_SCALE), draw(_SCALE)), "antipode of first": first.reflected(circle.centre)}
    second = draw(circle_points(circle)) if second_kind == "fresh" else choices[second_kind]
    return circle, sense, start, first, second


def _centre(centre: tuple[Fraction, Fraction]) -> ExactPoint2[WorldXY]:
    return ExactPoint2.build(centre[0], centre[1])


# ---------------------------------------------------------------------------
# Oracle self-checks
# ---------------------------------------------------------------------------


def test_oracle_decides_folded_ties_by_exact_squaring() -> None:
    sqrt8_minus_2sqrt2 = _difference(Surd(Fraction(0), Fraction(1), Fraction(8)), Surd(Fraction(0), Fraction(2), Fraction(2)))
    assert _decide(sqrt8_minus_2sqrt2) == Decision(0, "exact squaring")
    near = _difference(Surd(Fraction(10) ** -280, Fraction(1), Fraction(8)), Surd(Fraction(0), Fraction(2), Fraction(2)))
    assert _decide(near) == Decision(1, "exact squaring")
    assert _decide(_difference(Surd(Fraction(3), Fraction(1), Fraction(2)), Surd(Fraction(2), Fraction(1), Fraction(5)))) == Decision(1, "enclosure")


@settings(max_examples=EXAMPLES, deadline=None)
@given(x=surds(), y=surds(), z=surds())
def test_oracle_exact_fallback_agrees_with_certified_enclosure(x: Surd, y: Surd, z: Surd) -> None:
    f = _product(x, y).plus(Radical4(z.a, z.b, Fraction(0), Fraction(0), x.c, y.c))
    decision = _decide(f)
    event(f"oracle route: {decision.route}")
    if decision.route == "enclosure":
        assert _exact_sign4(f) == decision.sign


# ---------------------------------------------------------------------------
# One-root numbers
# ---------------------------------------------------------------------------


@settings(max_examples=EXAMPLES, deadline=None)
@given(x=surds())
def test_build_folds_without_changing_value(x: Surd) -> None:
    folded = x.number()
    assert folded.radicand >= 0 and (folded.coefficient == 0) == (folded.radicand == 0)
    assert _oracle_sign(_difference(x, Surd(folded.rational, folded.coefficient, folded.radicand))) == 0


@settings(max_examples=EXAMPLES, deadline=None)
@given(x=surds())
def test_sign_matches_oracle(x: Surd) -> None:
    assert x.number().sign() == _oracle_sign(x.form())


@settings(max_examples=EXAMPLES, deadline=None)
@given(pair=comparable_pairs())
def test_compare_matches_oracle(pair: tuple[Surd, Surd]) -> None:
    x, y = pair
    expected = _oracle_sign(_difference(x, y))
    left, right = x.number(), y.number()
    assert left.compare(right) == expected
    assert right.compare(left) == -expected
    assert (left == right, left < right, left > right) == (expected == 0, expected < 0, expected > 0)


@settings(max_examples=EXAMPLES, deadline=None)
@given(x=surds(), k=_SCALE, m=_SCALE)
def test_equal_values_under_different_radicands_compare_equal(x: Surd, k: int, m: int) -> None:
    assert x.number().compare(x.rescaled(k, m).number()) == 0
    assert x.number() == x.rescaled(k, m).number()


@settings(max_examples=EXAMPLES, deadline=None)
@given(x=surds())
def test_same_root_arithmetic_matches_oracle(x: Surd) -> None:
    y = Surd(x.a - x.b, x.a * x.b, x.c)
    left, right = x.number(), y.number()
    assert _oracle_sign(_minus_folded(x.form().plus(y.form()), left + right)) == 0
    assert _oracle_sign(_minus_folded(x.form().minus(y.form()), left - right)) == 0
    product = _product(x, y)
    as_one_root = Radical4(product.a + product.d * x.c, product.b + product.c, Fraction(0), Fraction(0), x.c, Fraction(0))
    assert _oracle_sign(_minus_folded(as_one_root, left * right)) == 0


@settings(max_examples=EXAMPLES, deadline=None)
@given(pair=comparable_pairs())
def test_cross_root_arithmetic_is_refused(pair: tuple[Surd, Surd]) -> None:
    left, right = pair[0].number(), pair[1].number()
    if left.radicand != 0 and right.radicand != 0 and left.radicand != right.radicand:
        with pytest.raises(CrossRootArithmeticError):
            left + right
        with pytest.raises(CrossRootArithmeticError):
            left * right


# ---------------------------------------------------------------------------
# Points: orientation and order on a circle
# ---------------------------------------------------------------------------


@settings(max_examples=EXAMPLES, deadline=None)
@given(triple=orientation_triples())
def test_orientation_matches_oracle(triple: tuple[tuple[Fraction, Fraction], RawPoint, RawPoint]) -> None:
    centre, first, second = triple
    expected = _oracle_sign(_cross(centre, first, second))
    assert orientation(_centre(centre), first.point(), second.point()) == expected
    assert orientation(_centre(centre), second.point(), first.point()) == -expected


@settings(max_examples=EXAMPLES, deadline=None)
@given(case=arc_cases())
def test_compare_on_arc_matches_angle_order(case: tuple[Circle, Sense, RawPoint, RawPoint, RawPoint]) -> None:
    circle, sense, start, first, second = case
    expected = _arc_order(circle.centre, sense, start, first, second)
    centre, start_point, first_point, second_point = _centre(circle.centre), start.point(), first.point(), second.point()
    assert compare_on_arc(centre, sense, start_point, first_point, second_point) == expected
    assert compare_on_arc(centre, sense, start_point, second_point, first_point) == -expected
