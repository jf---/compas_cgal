"""Exact rigid motions for congruence invariants.

A rotation built from a Pythagorean triple (a, b, c) with a^2 + b^2 = c^2 has
cos = a/c and sin = b/c, both rational. Every rotated coordinate therefore stays
in the rationals the exact kernel injects without loss, so "congruent inputs must
produce identical verdicts" is an EXACTLY checkable invariant rather than an
approximate one. Irrational rotations would smear the test with representation
noise and could not distinguish a real bug from a rounding artifact.
"""

from __future__ import annotations

from typing import Callable
from typing import Tuple

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# (a, b, c) with a^2 + b^2 = c^2. Chosen to span the first quadrant coarsely:
# ~36.87 deg, ~22.62 deg, ~67.38 deg, ~16.26 deg. Each gives exactly rational
# cos = a/c and sin = b/c.
PYTHAGOREAN_ROTATIONS: Tuple[Tuple[int, int, int], ...] = (
    (4, 3, 5),
    (12, 5, 13),
    (5, 12, 13),
    (24, 7, 25),
)


def rotate_point(x: float, y: float, triple: tuple[int, int, int]) -> tuple[float, float]:
    """Rotate a point by the exact rational rotation given by a Pythagorean triple.

    Args:
        x: Point x coordinate.
        y: Point y coordinate.
        triple: ``(a, b, c)`` with ``a**2 + b**2 == c**2``.

    Returns:
        The rotated ``(x, y)``.
    """
    a, b, c = triple
    cos_t, sin_t = a / c, b / c
    return (cos_t * x - sin_t * y, sin_t * x + cos_t * y)


def rotate_spec(spec: PocketSpec, triple: tuple[int, int, int]) -> PocketSpec:
    """Return *spec* rotated about the origin by an exact rational rotation.

    Args:
        spec: The instance to rotate.
        triple: ``(a, b, c)`` with ``a**2 + b**2 == c**2``.

    Returns:
        A new spec whose geometry is congruent to the input.
    """
    a, b, c = triple
    return _remap(spec, lambda px, py: rotate_point(px, py, (a, b, c)), f"_rot{a}_{b}_{c}")


def translate_spec(spec: PocketSpec, dx: float, dy: float) -> PocketSpec:
    """Return *spec* translated by ``(dx, dy)``.

    Args:
        spec: The instance to translate.
        dx: Shift along x.
        dy: Shift along y.

    Returns:
        A new spec whose geometry is congruent to the input.
    """
    return _remap(spec, lambda px, py: (px + dx, py + dy), f"_t{dx:g}_{dy:g}")


def _remap(spec: PocketSpec, fn: Callable[[float, float], tuple[float, float]], suffix: str) -> PocketSpec:
    """Apply a coordinate map to every ring of *spec*.

    Args:
        spec: The instance to transform.
        fn: Callable mapping ``(x, y)`` to ``(x, y)``.
        suffix: Appended to the instance name.

    Returns:
        The transformed spec.
    """

    def ring(polygon: Polygon) -> Polygon:
        return Polygon([[*fn(p[0], p[1]), 0.0] for p in polygon.points])

    return PocketSpec.build(
        name=spec.name + suffix,
        family=spec.family,
        polygon=ring(spec.polygon),
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        holes=tuple(ring(h) for h in spec.holes),
        params=dict(spec.params),
    )
