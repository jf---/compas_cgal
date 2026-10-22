"""Pins on the exact-target certificate's vertex-circle cuts: every swept circle is cut where it crosses a one-root vertex circle.

∂T ⊆ declared supports ∪ vertex circles, so the certificate cuts every swept
circle at its crossings with each one-root vertex circle (``events.add_vertex``)
and decides each resulting sub-arc on its own. A swept circle with k > 0
distinct cut points has exactly k sub-arcs; with none it is one edge.

The fixture is the convex-arc pocket (r = 1): its kink at (0, 0) gives the
Ω* vertex w = (1, 15/2 - √205/2), where the offset line x = 1 meets the
offset arc of radius 23/2 about (10, 15/2). One cut of guide radius R < r has a
single boundary circle C of radius R + r. The generator keeps C strictly away
from every declared support (x = 0, x = 20, y = 20 and the declared arc of
radius 25/2), from the other kink's vertex circle and from the rational vertex
circles at the square corners, all checked exactly here. C then meets only
D_r(w)'s boundary, and the count of its sub-arcs is 2 when C crosses that
circle (R < |c - w| < R + 2) and 1 otherwise; |c - w|² = a + b√205 with
b = cy - 15/2 != 0, so a tangency is impossible for rational c.
"""

from __future__ import annotations

from collections.abc import Iterator
from fractions import Fraction

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.coverage_certificate.circle_core import FilterMode
from compas_cgal.coverage_certificate.native_certificate import CertificateReport
from compas_cgal.coverage_certificate.native_certificate import run_certificate
from tests.benchmarks.test_coverage_certificate_exact_target import TOOL
from tests.benchmarks.test_coverage_certificate_exact_target import UNUSED_PATH
from tests.benchmarks.test_coverage_certificate_exact_target import _sweeps
from tests.benchmarks.test_coverage_certificate_exact_target import convex_arc

F = Fraction
R = F(1)
# The kink vertex w = (1, 15/2 - √RADICAND / 2) of Omega* in the convex-arc pocket.
RADICAND = F(205)
ARC_CENTRE = (F(10), F(15, 2))
ARC_RADIUS = F(25, 2)


@pytest.fixture(autouse=True, scope="module")
def _release_native_domains() -> Iterator[None]:
    """Drop the cached native Omega* handle before interpreter shutdown."""
    yield
    convex_arc.cache_clear()


def _sign_plus_root(a: Fraction, b: Fraction, radicand: Fraction) -> int:
    """Exact sign of a + b·√radicand (radicand > 0)."""
    sign_a, sign_b = (a > 0) - (a < 0), (b > 0) - (b < 0)
    if sign_b == 0 or sign_a == sign_b:
        return sign_a
    if sign_a == 0:
        return sign_b
    gap = a * a - b * b * radicand
    return sign_a * ((gap > 0) - (gap < 0))


def _squared_distance_to_kink(x: Fraction, y: Fraction) -> tuple[Fraction, Fraction]:
    """|c - w|² = a + b·√205 with w = (1, 15/2 - √205/2)."""
    dy = y - F(15, 2)
    return (x - 1) ** 2 + dy * dy + RADICAND / 4, dy


def _crosses_kink_circle(x: Fraction, y: Fraction, guide: Fraction) -> bool:
    """R < |c - w| < R + 2 r for the boundary circle of radius R + r against the vertex circle of radius r."""
    a, b = _squared_distance_to_kink(x, y)
    inner, outer = guide * guide, (guide + 2 * R) ** 2
    return _sign_plus_root(a - inner, b, RADICAND) > 0 and _sign_plus_root(a - outer, b, RADICAND) < 0


def _clear_of_other_supports(x: Fraction, y: Fraction, guide: Fraction) -> bool:
    """C (radius R + r about c) misses every declared support and every other vertex circle, exactly."""
    reach = guide + R
    if not (x - reach > 0 and x + reach < 20 and y + reach < 20):
        return False
    dx, dy = x - ARC_CENTRE[0], y - ARC_CENTRE[1]
    if not (ARC_RADIUS - reach > 0 and dx * dx + dy * dy < (ARC_RADIUS - reach) ** 2):
        return False
    # The other kink's vertex circle is centred at x = 19; the rational corner ones at y = 19.
    return x + reach + R < 19 and y + reach + R < 19


def _run(x: float, y: float, guide: float, mode: FilterMode) -> CertificateReport:
    _, _, target = convex_arc()
    return run_certificate(UNUSED_PATH, TOOL, _sweeps([(x, y, guide)]), exact_target=target, filter_mode=mode, audit=True)


def test_the_fixture_has_the_kink_vertex_circle() -> None:
    _, _, target = convex_arc()
    assert sum(1 for disk in target.vertex_disks if disk.radicand == RADICAND / 4) == 2


@pytest.mark.parametrize("mode", [FilterMode.LAZY, FilterMode.EXACT_ONLY], ids=lambda m: m.name)
def test_a_swept_circle_crossing_the_kink_vertex_circle_is_cut_into_two_sub_arcs(mode: FilterMode) -> None:
    """c = (1.472, 0.506), R = 1/16: C crosses D_r(w)'s boundary twice and nothing else."""
    x, y, guide = 1.472, 0.506, 0.0625
    assert _clear_of_other_supports(F(x), F(y), F(guide))
    assert _crosses_kink_circle(F(x), F(y), F(guide))
    report = _run(x, y, guide, mode)
    assert report.counters["vertex_circle_pairs"] == 1
    assert report.counters["sub_arcs"] == 2
    assert report.edge_count == 2


@pytest.mark.parametrize("mode", [FilterMode.LAZY, FilterMode.EXACT_ONLY], ids=lambda m: m.name)
def test_a_swept_circle_in_the_kink_box_but_off_its_circle_stays_one_edge(mode: FilterMode) -> None:
    """Control: c = (2.9, 2.4), R = 1/16 overlaps the vertex circle's box only; no cut, one edge."""
    x, y, guide = 2.9, 2.4, 0.0625
    assert _clear_of_other_supports(F(x), F(y), F(guide))
    assert not _crosses_kink_circle(F(x), F(y), F(guide))
    report = _run(x, y, guide, mode)
    assert report.counters["vertex_circle_pairs"] == 1
    assert report.counters["sub_arcs"] == 1
    assert report.edge_count == 1


_GUIDE = st.sampled_from([F(1, 64), F(1, 16), F(1, 4), F(1, 2), F(3, 4)])
_X = st.integers(min_value=1040, max_value=4000).map(lambda k: F(k, 1024))
_Y = st.integers(min_value=-1024, max_value=3072).map(lambda k: F(k, 1024))


@settings(max_examples=150, deadline=None)
@given(_X, _Y, _GUIDE)
def test_sub_arcs_of_a_lone_circle_equal_its_exact_vertex_circle_crossings(x: Fraction, y: Fraction, guide: Fraction) -> None:
    """A lone swept circle near the kink: sub-arcs == 2 iff it crosses the one-root vertex circle, else 1."""
    if not _clear_of_other_supports(x, y, guide):
        return
    report = _run(float(x), float(y), float(guide), FilterMode.LAZY)
    expected = 2 if _crosses_kink_circle(x, y, guide) else 1
    assert report.counters["sub_arcs"] == expected
    assert report.edge_count == expected
