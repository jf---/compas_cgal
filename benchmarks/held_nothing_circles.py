"""Exact removal of NOTHING circles: circles whose swept region removes no material.

A full machining circle of guide radius rho about c, swept by a tool of radius
r, sweeps the closed annulus A = {p : max(rho - r, 0) <= |p - c| <= rho + r}.
A circle removes NOTHING when A lies in the union of the regions swept before
it. This module proves a sufficient form of that statement exactly: A lies in
the union of the swept annuli of EARLIER CIRCLES (connectors ignored). Leaving
connectors out only shrinks the union, so every circle proven here removes
nothing on the true stock as well; circles covered only with a connector's help
are kept.

Decision (exact over the injected doubles): a local exact stock
(``_stock_2.Stock2``, Epeck circle-segment traits) is seeded with an axis square
containing A, every selected earlier annulus is subtracted exactly
(``subtract_circle_sweep``), and ``Stock2.circle_sweep_removal`` decides A
against the rest (a regularized interior intersection). ``NOTHING`` is the
proof.

Selection is not decision. Which earlier circles enter the union is chosen in
binary64 (``_overlapping_earlier``), and a binary64 point screen
(``_screen``) rejects circles before the exact test: a sample point of A
clearly outside every selected annulus means material is left. Leaving a sweep
out, or rejecting a circle, can only keep a circle, never remove one. The
screen tests membership with the outward relative slack
``SCREEN_RELATIVE_SLACK``, which dwarfs the binary64 error of its squared
distances (see the constant), so a circle whose sample points are all exactly
covered (an exact duplicate, whose rim samples sit on an earlier rim) passes it.

Removing circle k makes circle k+1's predecessor the last kept circle; that
pair's engagement is re-decided by the exact standard-model predicate
(``exact_predecessor_engagement(...).within(cap)``) before k is removed, and a
NOTHING circle whose removal would put its successor over the cap is kept.
Connectors are built after this pass from the surviving circles, so the
connector from the last kept circle to k+1 is formed by the ordinary connector
builder.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Generic
from typing import TypeVar

import numpy as np

from benchmarks.errors import BenchmarkError
from benchmarks.held_standard_placement import PaperCircleCandidate
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import exact_predecessor_engagement
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius

C = TypeVar("C")

# Relative outward slack of the binary64 screen on squared distances. The
# operands are coordinates and radii below 2^10 mm; one squared distance or
# squared radius carries a few ulps of relative error (about 1e-15) plus the
# sample point's own rounding (a few ulps of a coordinate below 2^11 mm, so below
# 1e-12 relative on any squared distance above 2^-20 mm^2). 1e-9 exceeds that by
# three orders of magnitude; the screen only rejects, so a larger slack costs
# exact work, never soundness.
SCREEN_RELATIVE_SLACK = 1e-9

# Screen sample grids (directions, radii) over the annulus, coarse first: each
# grid runs only on circles the previous one passed. Measured on Figure 5
# (2,010 circles, build/nr/dense.py): 128 x 9 leaves 21 circles, 512 x 17 leaves
# 14 for the exact test, 6 of them NOTHING.
SCREEN_GRIDS = ((48, 5), (128, 9), (512, 17))

# Half-width factor of the local stock's seed square about the annulus's outer
# radius. The square's exact containment of the annulus is checked, not assumed
# (``_seed_square``); the factor only keeps that check from ever failing.
SEED_SQUARE_FACTOR = 1.25


class NothingCircleError(BenchmarkError):
    """Base error of the NOTHING-circle pass."""


class SeedSquareContainmentError(NothingCircleError):
    """The local stock's binary64 seed square does not contain the annulus exactly."""


@dataclass(frozen=True)
class _Annulus:
    """A swept full circle, exact: centre (x, y), guide radius, tool radius."""

    x: Fraction
    y: Fraction
    guide: Fraction
    tool: Fraction

    @classmethod
    def of(cls, circle: PaperCircleCandidate, tool_radius: ToolRadius) -> _Annulus:
        return cls(Fraction(float(circle.center.x)), Fraction(float(circle.center.y)), Fraction(float(circle.guide_radius.value)), Fraction(float(tool_radius.value)))

    def contains(self, px: Fraction, py: Fraction) -> bool:
        """Closed membership of the exact point (px, py)."""
        squared = (px - self.x) ** 2 + (py - self.y) ** 2
        outer = self.guide + self.tool
        if squared > outer * outer:
            return False
        if self.guide <= self.tool:
            return True
        inner = self.guide - self.tool
        return squared >= inner * inner


def _arrays(circles: Sequence[PaperCircleCandidate]) -> np.ndarray:
    return np.array([(float(c.center.x), float(c.center.y), float(c.guide_radius.value)) for c in circles], dtype=np.float64).reshape(-1, 3)


def _overlapping_earlier(geometry: np.ndarray, index: int, tool: float) -> np.ndarray:
    """Indices below ``index`` whose annulus may meet ``index``'s (binary64 selection, not a decision)."""
    x, y, rho = geometry[index]
    earlier = geometry[:index]
    distance = np.hypot(earlier[:, 0] - x, earlier[:, 1] - y)
    inner = np.maximum(rho - tool, 0.0)
    earlier_inner = np.maximum(earlier[:, 2] - tool, 0.0)
    reach = (distance <= earlier[:, 2] + rho + 2.0 * tool) & (distance >= earlier_inner - rho - tool) & (distance >= inner - earlier[:, 2] - tool)
    return np.nonzero(reach)[0]


def _screen(geometry: np.ndarray, index: int, selected: np.ndarray, tool: float) -> bool:
    """Binary64 rejection: False when some sample point of the annulus is clearly outside every selected annulus."""
    if selected.size == 0:
        return False
    x, y, rho = geometry[index]
    other = geometry[selected]
    outer = (other[:, 2] + tool) ** 2 * (1.0 + SCREEN_RELATIVE_SLACK)
    inner = np.where(other[:, 2] > tool, (other[:, 2] - tool) ** 2, 0.0) * (1.0 - SCREEN_RELATIVE_SLACK)
    for directions, radii in SCREEN_GRIDS:
        angle = np.linspace(0.0, 2.0 * np.pi, directions, endpoint=False)
        radius = np.linspace(max(rho - tool, 0.0), rho + tool, radii)
        px = (x + radius[:, None] * np.cos(angle)[None, :]).reshape(-1, 1)
        py = (y + radius[:, None] * np.sin(angle)[None, :]).reshape(-1, 1)
        squared = (px - other[None, :, 0]) ** 2 + (py - other[None, :, 1]) ** 2
        if not ((squared <= outer[None, :]) & (squared >= inner[None, :])).any(axis=1).all():
            return False
    return True


def _seed_square(target: _Annulus) -> np.ndarray:
    x, y = float(target.x), float(target.y)
    half = SEED_SQUARE_FACTOR * float(target.guide + target.tool)
    corners = ((x - half, y - half), (x + half, y - half), (x + half, y + half), (x - half, y + half))
    outer = target.guide + target.tool
    xs = sorted({Fraction(cx) for cx, _ in corners})
    ys = sorted({Fraction(cy) for _, cy in corners})
    if xs[0] > target.x - outer or xs[-1] < target.x + outer or ys[0] > target.y - outer or ys[-1] < target.y + outer:
        raise SeedSquareContainmentError(f"seed square of half-width {half!r} does not contain the annulus about ({x!r}, {y!r}).")
    return np.array([(cx, cy, 0.0) for cx, cy in corners], dtype=np.float64)


def _union_covers(target: _Annulus, earlier: Sequence[_Annulus]) -> bool:
    stock = _stock_2.Stock2(_seed_square(target), [])
    for annulus in earlier:
        stock.subtract_circle_sweep(float(annulus.x), float(annulus.y), float(annulus.guide), float(annulus.tool))
    verdict = stock.circle_sweep_removal(float(target.x), float(target.y), float(target.guide), float(target.tool))
    return verdict == _stock_2.SweepRemovalVerdict.NOTHING


def union_covers_annulus(target: PaperCircleCandidate, earlier: Sequence[PaperCircleCandidate], tool_radius: ToolRadius) -> bool:
    """Step 2 alone, exact: ``target``'s swept annulus lies in the union of ``earlier``'s swept annuli."""
    return _union_covers(_Annulus.of(target, tool_radius), [_Annulus.of(circle, tool_radius) for circle in earlier])


def removes_nothing(target: PaperCircleCandidate, earlier: Sequence[PaperCircleCandidate], tool_radius: ToolRadius) -> bool:
    """Exact sufficient test: True proves ``target``'s swept annulus lies in the union of ``earlier``'s.

    False means either material is left or the proof was not attempted (the
    selection and the screen only ever reject).
    """
    geometry = _arrays((*earlier, target))
    return _decide(geometry, len(earlier), [_Annulus.of(circle, tool_radius) for circle in (*earlier, target)], float(tool_radius.value))


def _decide(geometry: np.ndarray, index: int, annuli: Sequence[_Annulus], tool: float) -> bool:
    selected = _overlapping_earlier(geometry, index, tool)
    if not _screen(geometry, index, selected, tool):
        return False
    return _union_covers(annuli[index], [annuli[j] for j in selected])


@dataclass(frozen=True)
class NothingCirclePass(Generic[C]):
    """The surviving circles and the indices (into the input) the pass removed or had to keep.

    ``kept_for_successor_cap``: proven NOTHING, kept because the successor
    would exceed the cap after the new predecessor (or the model has no
    engagement for that pair).
    """

    circles: tuple[C, ...]
    removed: tuple[int, ...]
    kept_for_successor_cap: tuple[int, ...]
    proven_nothing: int
    union_tests: int


def _successor_within_cap(predecessor: PaperCircleCandidate, successor: PaperCircleCandidate, tool_radius: ToolRadius, cap: EngagementCap) -> bool:
    try:
        return exact_predecessor_engagement(predecessor, successor, tool_radius).within(cap)
    except StandardPlacementFragmentationError:
        return False


def remove_nothing_circles(circles: Sequence[C], paper: Callable[[C], PaperCircleCandidate], tool_radius: ToolRadius, cap: EngagementCap) -> NothingCirclePass[C]:
    """Remove, in path order, every circle proven NOTHING whose successor stays within ``cap`` after the new predecessor.

    Circle 0 (the entry) is never removed. Earlier circles are all the input
    circles before k: a removed circle's annulus lies in the union of the kept
    ones before it, so that union equals the union of all earlier circles.
    """
    candidates = [paper(circle) for circle in circles]
    if len(candidates) < 2:
        return NothingCirclePass(tuple(circles), (), (), 0, 0)
    geometry = _arrays(candidates)
    annuli = [_Annulus.of(circle, tool_radius) for circle in candidates]
    tool = float(tool_radius.value)
    kept = [0]
    removed: list[int] = []
    refused: list[int] = []
    union_tests = 0
    for index in range(1, len(candidates)):
        selected = _overlapping_earlier(geometry, index, tool)
        nothing = False
        if _screen(geometry, index, selected, tool):
            union_tests += 1
            nothing = _union_covers(annuli[index], [annuli[j] for j in selected])
        if not nothing:
            kept.append(index)
            continue
        if index + 1 < len(candidates) and not _successor_within_cap(candidates[kept[-1]], candidates[index + 1], tool_radius, cap):
            refused.append(index)
            kept.append(index)
            continue
        removed.append(index)
    return NothingCirclePass(tuple(circles[k] for k in kept), tuple(removed), tuple(refused), len(removed) + len(refused), union_tests)
