"""Exact window geometry of the last-16 actual-stock certificate: the clip proof and exact sweep omission.

A circle's certificate is decided on a WINDOW stock: Omega minus a set of
prior sweeps, minus the closed annulus ``{R <= |x - c| <= outer}`` about the
circle's guide centre ``c`` (the CLIP). Two exact lemmas make the window stock
interchangeable with the unclipped, unselected stock for the decision:

**Read set.** ``engagement_circle_exceeded`` classifies each loop station ``s``
by walking the cutter circle ``circle(s, r)`` through the stock arrangement and
reading ``face->contained()`` on its sub-curves
(``src/audit_exact_station_2.cpp``, ``exact_engaged_arcs``); the swept-prefix
clip that follows is analytic (start disk and annular sector) and reads no
stock. Station 0 is the authored entry point; every other station is
``c + rho * u`` with ``u`` an exact rational unit vector rotated by the exact
chart (``src/engagement_2.cpp``, ``circle_prefix_motion``), so it lies EXACTLY
on the circle of radius ``rho`` about ``c``. The decision therefore depends on
the stock only through its trace on the union of the closed disks
``D(s, r)``.

**Clip.** The clip removes only points with ``|x - c| >= R``. If every
station's closed tool disk lies in the OPEN disk ``|x - c| < R`` the clip and
the stock agree on an open neighbourhood of every read point, so the regularized
difference leaves the read trace unchanged. ``prove_clip_contains_stations``
proves exactly that, over the rationals the doubles ARE: ``rho + r < R`` for
the chart stations and ``|e - c|^2 < (R - r)^2`` with ``R > r`` for the entry
``e``. ``R`` itself is chosen in binary64 (a cost choice: the window's size);
only the exact proof admits it.

**Omission.** A prior sweep may be left out of the window stock only if its
removed region cannot meet the open window disk, so that Omega minus it and
Omega agree on the window. ``prove_sweep_misses_window`` decides that exactly
against a superset of the removed region: the closed annulus
``[max(rho_j - r, 0), rho_j + r]`` for a circle (``full_turn_annulus_bounds``,
removed exactly), the closed true capsule for a connector segment (the replay
removes ``Stock_UNDER``, a subset of it). A failed proof means "include": the
default is inclusion, which can only cost time.

Every comparison is on ``fractions.Fraction`` values injected exactly from the
binary64 inputs. No tolerance, margin or float decision is in this module's
decision path.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


class WindowGeometryError(Exception):
    """Base class of the exact window-geometry failures."""


class WindowClipContainmentError(WindowGeometryError):
    """A station's closed tool disk is not inside the open clip disk: the clip could remove material the decision reads."""


class NonPositiveToolRadiusError(WindowGeometryError):
    """The tool radius is not positive."""


@dataclass(frozen=True)
class ExactDisk:
    """The open window disk ``|x - centre| < radius``; coordinates are the exact rationals of binary64 inputs."""

    cx: Fraction
    cy: Fraction
    radius: Fraction

    @classmethod
    def build(cls, cx: float, cy: float, radius: float) -> ExactDisk:
        return cls(Fraction(cx), Fraction(cy), Fraction(radius))


@dataclass(frozen=True)
class ExactCircleReach:
    """A circle sweep's removed region bound: the closed annulus ``inner <= |x - centre| <= outer``."""

    cx: Fraction
    cy: Fraction
    inner: Fraction
    outer: Fraction

    @classmethod
    def build(cls, cx: float, cy: float, guide_radius: float, tool_radius: float) -> ExactCircleReach:
        """Exact bounds of ``Stock2.subtract_circle_sweep(cx, cy, guide_radius, tool_radius)``.

        Raises:
            NonPositiveToolRadiusError: ``tool_radius <= 0``.
        """
        rho, r = Fraction(guide_radius), Fraction(tool_radius)
        if r <= 0:
            raise NonPositiveToolRadiusError(f"tool radius {tool_radius!r} must be positive")
        return cls(Fraction(cx), Fraction(cy), max(rho - r, Fraction(0)), rho + r)


@dataclass(frozen=True)
class ExactSegmentReach:
    """A connector segment's removed region bound: the closed capsule of radius ``radius`` about it."""

    x0: Fraction
    y0: Fraction
    x1: Fraction
    y1: Fraction
    radius: Fraction

    @classmethod
    def build(cls, x0: float, y0: float, x1: float, y1: float, tool_radius: float) -> ExactSegmentReach:
        """Exact capsule bounding ``Stock2.subtract_line_sweep`` (whose UNDER region it contains).

        Raises:
            NonPositiveToolRadiusError: ``tool_radius <= 0``.
        """
        r = Fraction(tool_radius)
        if r <= 0:
            raise NonPositiveToolRadiusError(f"tool radius {tool_radius!r} must be positive")
        return cls(Fraction(x0), Fraction(y0), Fraction(x1), Fraction(y1), r)


def _boxes_disjoint(window: ExactDisk, lo_x: Fraction, lo_y: Fraction, hi_x: Fraction, hi_y: Fraction) -> bool:
    """The closed box ``[lo, hi]`` misses the open bounding box of ``window``."""
    return hi_x <= window.cx - window.radius or lo_x >= window.cx + window.radius or hi_y <= window.cy - window.radius or lo_y >= window.cy + window.radius


def circle_misses_window(window: ExactDisk, reach: ExactCircleReach) -> bool:
    """Exact: the closed annulus ``reach`` does not meet the open disk ``window``.

    Disjoint iff the window lies beyond the outer circle (``d >= R + outer``) or
    inside the hole (``d + R <= inner``), both decided on squared rationals.
    """
    if _boxes_disjoint(window, reach.cx - reach.outer, reach.cy - reach.outer, reach.cx + reach.outer, reach.cy + reach.outer):
        return True
    dx, dy = reach.cx - window.cx, reach.cy - window.cy
    squared = dx * dx + dy * dy
    beyond = window.radius + reach.outer
    if squared >= beyond * beyond:
        return True
    hole = reach.inner - window.radius
    return hole >= 0 and squared <= hole * hole


def segment_misses_window(window: ExactDisk, reach: ExactSegmentReach) -> bool:
    """Exact: the closed capsule ``reach`` does not meet the open disk ``window`` (``dist(c, segment) >= R + r``)."""
    r = reach.radius
    if _boxes_disjoint(window, min(reach.x0, reach.x1) - r, min(reach.y0, reach.y1) - r, max(reach.x0, reach.x1) + r, max(reach.y0, reach.y1) + r):
        return True
    ux, uy = reach.x1 - reach.x0, reach.y1 - reach.y0
    wx, wy = window.cx - reach.x0, window.cy - reach.y0
    length = ux * ux + uy * uy
    t = Fraction(0) if length == 0 else min(max((wx * ux + wy * uy) / length, Fraction(0)), Fraction(1))
    nx, ny = wx - t * ux, wy - t * uy
    reach_sum = window.radius + r
    return nx * nx + ny * ny >= reach_sum * reach_sum


def prove_clip_contains_stations(window: ExactDisk, guide_radius: float, entry_x: float, entry_y: float, tool_radius: float) -> None:
    """Prove exactly that every station's closed tool disk lies in the open clip disk ``window``.

    Chart stations lie exactly on the guide circle: ``rho + r < R``. The entry
    station ``e``: ``R > r`` and ``|e - c|^2 < (R - r)^2``.

    Raises:
        NonPositiveToolRadiusError: ``tool_radius <= 0``.
        WindowClipContainmentError: either inequality fails.
    """
    r = Fraction(tool_radius)
    if r <= 0:
        raise NonPositiveToolRadiusError(f"tool radius {tool_radius!r} must be positive")
    if not Fraction(guide_radius) + r < window.radius:
        raise WindowClipContainmentError(f"chart stations: rho + r = {float(Fraction(guide_radius) + r)!r} is not below the clip radius {float(window.radius)!r}")
    slack = window.radius - r
    dx, dy = Fraction(entry_x) - window.cx, Fraction(entry_y) - window.cy
    if not (slack > 0 and dx * dx + dy * dy < slack * slack):
        raise WindowClipContainmentError(f"entry station ({entry_x!r}, {entry_y!r}): its tool disk is not inside the clip disk of radius {float(window.radius)!r}")
