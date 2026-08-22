# Benchmark Corpus & Reference-Problem Generator — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A parameterised corpus of pocket-machining reference problems, an instrumented runner that separates generation from certification cost and exposes the exact kernel's hidden clock (arrangement size and rational bit-length), and a reproduction of Held & Pfeiffer's Figure 6 — so that "is this fast enough?" and "how often does it say `unresolved`?" become measured numbers instead of estimates.

**Architecture:** A standalone top-level `benchmarks/` package (not shipped in the wheel) that depends **only** on the stable public surface — `compas_cgal.toolpath`, `compas_cgal.engagement`, `compas_cgal.stock`. Pocket families are pure generators emitting frozen `PocketSpec` values; a runner drives each spec through generate → audit and emits a typed `MeasurementRecord`. Analytic families additionally carry closed-form oracles, so the corpus doubles as a correctness net rather than only a stopwatch.

**Tech Stack:** Python 3.9-compatible (`requires-python = ">=3.9"`, ruff `target-version = "py39"`); compas (`Polygon`, `TOL`); numpy; pytest + pytest-xdist; nanobind/CGAL 6.0.1 for the one C++ instrumentation change; mypy for the strict gate.

## Global Constraints

- **Python 3.9 floor.** Every module starts with `from __future__ import annotations`. Annotations may use `X | Y` and `list[X]`; **runtime** code must not (no `isinstance(x, int | str)`, no runtime subscripting of builtins).
- **Branch-agnostic.** `benchmarks/` must not import `compas_cgal.adaptive`. The corpus has to run unchanged on both `jf/toolpath-redesign` and `codex/exact-certified-adaptive-phase1-t9`, because its purpose is to compare them.
- ruff: `line-length = 179`, `select = ["E", "F", "I"]`. Run `ruff format` and `ruff check --fix` before every Python commit.
- **`mypy --strict` passes for the whole `benchmarks/` package.** Added to the gate in Task 14.
- No `__all__` in any `__init__.py`; keep `__init__.py` minimal.
- **One named exception per failure mode** — never bare `raise ValueError("...")`.
- **No inline magic numbers.** Named module-level UPPERCASE constants with a one-line comment giving units and derivation. Prefer `compas.tolerance.TOL` predicates for tolerance decisions; note the traps: `is_positive(a, tol)` is strict `a > tol`; `is_close` mixes atol+rtol; passing `tol=0`/`rtol=0.0` is silently coerced to the default — use `is_between(v, t, t, atol=...)` for an atol-only proximity check.
- Google-style docstrings (`Args:` / `Returns:` / `Raises:`) — never Sphinx or numpydoc.
- **No `pytest.mark.skip`, `skipif`, or `xfail`.** Tests fail loud.
- Run tests with `pytest -n auto`.
- No conditional imports, no `HAS_*` flags, no `try: import X except ImportError`, no silent fallbacks.
- Commits: author **and** committer `Jelle Feringa <jelleferinga@gmail.com>`; concise conventional messages; **no attribution lines of any kind**.
- Never touch `docs/examples/example_isolines.py`.

## Why this corpus, in one paragraph

Held & Pfeiffer report 3–100 ms per pocket, but for **one** pocket (Fig. 5), with no hardware stated and no instance table — it is a plausibility figure, not a reproducible target. Their generator is fast because it never builds a global stock arrangement: it maintains the union-of-disks boundary as a sorted arc sequence updated in amortised-linear work, and gets engagement from a closed form (Eq. 7). Their engagement *distributions* are, in their own words, obtained because "for the sake of implementational simplicity, we resorted to a discretization." This project's generation is already competitive (0.00–0.02 s); the entire cost gap sits in a certification step Held does not perform. This corpus exists to measure that gap along the axes that actually drive it, and to measure the number that decides the business: how often the certifier answers `unresolved` on geometry it did not author.

## File Structure

```
benchmarks/
  __init__.py            minimal
  errors.py              Task 1  — named exceptions
  spec.py                Task 1  — PocketSpec + factory
  families/
    __init__.py          minimal
    analytic.py          Task 2  — disk, rectangle, stadium, arc channel (+ closed-form clearance)
    complexity.py        Task 6  — regular k-gon sweep, line/arc fraction sweep
    necks.py             Task 7  — dumbbell pinch sweep
    precision.py         Task 8  — scale x decimal-digit sweep
    topology.py          Task 9  — island-count sweep
  congruence.py          Task 3  — exact rational rigid motions (Pythagorean rotations)
  degeneracy.py          Task 9  — hand-built degenerate corpus
  instrument.py          Task 4  — arrangement-size + bit-length probes (Python side)
  measurement.py         Task 5  — MeasurementRecord schema
  runner.py              Task 5  — spec -> MeasurementRecord
  report.py              Task 10 — markdown + json emission
  mathsm.py              Task 11 — constant-spacing MATHSM baseline
  figure6.py             Task 12 — Held Fig. 6 reproduction
  external.py            Task 13 — third-party sketch-profile loader
  cli.py                 Task 14 — `python -m benchmarks.cli`
src/stock_2.h            Task 4  — MODIFY: ArrangementStats + accessors
src/stock_2.cpp          Task 4  — MODIFY: implementation + nanobind bindings
src/compas_cgal/stock.py Task 4  — MODIFY: Stock.arrangement_stats / coordinate_digits
tests/benchmarks/
  test_spec.py test_analytic.py test_congruence.py test_instrument.py
  test_measurement.py test_runner.py test_complexity.py test_necks.py
  test_precision.py test_topology.py test_degeneracy.py test_report.py
  test_mathsm.py test_figure6.py test_external.py test_cli.py
docs/benchmarks.md       Task 14 — the MkDocs page
```

Files that change together live together. Each family module owns exactly one sweep axis, so a reviewer can reject one axis without touching its neighbours.

---

### Task 1: `PocketSpec` foundation

**Files:**
- Create: `benchmarks/__init__.py`, `benchmarks/errors.py`, `benchmarks/spec.py`
- Test: `tests/benchmarks/test_spec.py`

**Interfaces:**
- Produces: `PocketSpec` (frozen dataclass) with fields `name: str`, `family: str`, `polygon: Polygon`, `holes: tuple[Polygon, ...]`, `tool_diameter: float`, `tea_cap_deg: float`, `params: Mapping[str, float]`; classmethod `PocketSpec.build(...) -> PocketSpec` owning all invariants. Exceptions `NonPositiveToolError`, `InvalidCapError`, `DegeneratePocketError`, `PocketNotSimpleError`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_spec.py
from __future__ import annotations

import pytest
from compas.geometry import Polygon

from benchmarks.errors import InvalidCapError, NonPositiveToolError, PocketNotSimpleError
from benchmarks.spec import PocketSpec

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


def test_build_accepts_valid_pocket() -> None:
    spec = PocketSpec.build(name="square", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=120.0)
    assert spec.tool_radius == pytest.approx(0.5)
    assert spec.holes == ()


def test_build_rejects_non_positive_tool() -> None:
    with pytest.raises(NonPositiveToolError):
        PocketSpec.build(name="x", family="analytic", polygon=SQUARE, tool_diameter=0.0, tea_cap_deg=120.0)


def test_build_rejects_cap_outside_range() -> None:
    with pytest.raises(InvalidCapError):
        PocketSpec.build(name="x", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=181.0)


def test_build_rejects_self_intersecting_boundary() -> None:
    bowtie = Polygon([[0, 0, 0], [10, 10, 0], [10, 0, 0], [0, 10, 0]])
    with pytest.raises(PocketNotSimpleError):
        PocketSpec.build(name="x", family="analytic", polygon=bowtie, tool_diameter=1.0, tea_cap_deg=120.0)


def test_spec_is_frozen() -> None:
    spec = PocketSpec.build(name="square", family="analytic", polygon=SQUARE, tool_diameter=1.0, tea_cap_deg=120.0)
    with pytest.raises(AttributeError):
        spec.tool_diameter = 2.0  # type: ignore[misc]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_spec.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmarks'`

- [ ] **Step 3: Write the exceptions**

```python
# benchmarks/errors.py
from __future__ import annotations


class BenchmarkError(Exception):
    """Base class for every benchmark-corpus failure mode."""


class NonPositiveToolError(BenchmarkError):
    """A tool diameter was zero, negative, or NaN."""


class InvalidCapError(BenchmarkError):
    """An engagement cap fell outside the exact kernel's contract of (0, 180] degrees."""


class DegeneratePocketError(BenchmarkError):
    """A generated pocket has zero area, fewer than three vertices, or no room for the tool."""


class PocketNotSimpleError(BenchmarkError):
    """A pocket boundary self-intersects, so it is not a valid general polygon."""
```

- [ ] **Step 4: Write `PocketSpec`**

```python
# benchmarks/spec.py
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError, InvalidCapError, NonPositiveToolError, PocketNotSimpleError

# The exact kernel's documented cap contract is theta in (0, pi]: one engaged run
# subtends at most a half turn before ">pi" is an exact orientation verdict, so a
# cap above 180 degrees is meaningless (src/engagement_2.cpp::certify_segment_tea).
MAX_CAP_DEG = 180.0

# A pocket must admit at least one full tool circle; below this multiple of the
# tool radius the pocket is a slot the generator cannot enter, which is a corpus
# authoring error rather than an interesting instance.
MIN_AREA_TOOL_RADII_SQ = 4.0


@dataclass(frozen=True)
class PocketSpec:
    """One reference machining problem: a pocket, a tool, and an engagement cap.

    Attributes:
        name: Unique instance name; becomes the row key in every report.
        family: Sweep family that produced this instance (e.g. ``"necks"``).
        polygon: Outer pocket boundary in the world XY plane, CCW.
        holes: Island boundaries, possibly empty.
        tool_diameter: Cutter diameter in the same units as the polygon.
        tea_cap_deg: Maximum tool-engagement angle in degrees.
        params: The sweep coordinates that produced this instance, for plotting.
    """

    name: str
    family: str
    polygon: Polygon
    holes: tuple[Polygon, ...] = ()
    tool_diameter: float = 1.0
    tea_cap_deg: float = 120.0
    params: Mapping[str, float] = field(default_factory=dict)

    @property
    def tool_radius(self) -> float:
        """Half the tool diameter."""
        return 0.5 * self.tool_diameter

    @property
    def tea_cap_rad(self) -> float:
        """The engagement cap in radians, as the kernel consumes it."""
        return math.radians(self.tea_cap_deg)

    @classmethod
    def build(
        cls,
        name: str,
        family: str,
        polygon: Polygon,
        tool_diameter: float,
        tea_cap_deg: float,
        holes: tuple[Polygon, ...] = (),
        params: Mapping[str, float] | None = None,
    ) -> "PocketSpec":
        """Validate every invariant and return a frozen spec.

        Args:
            name: Unique instance name.
            family: Sweep family name.
            polygon: Outer boundary.
            tool_diameter: Cutter diameter; must be finite and positive.
            tea_cap_deg: Engagement cap; must lie in (0, 180].
            holes: Island boundaries.
            params: Sweep coordinates recorded for plotting.

        Returns:
            The validated spec.

        Raises:
            NonPositiveToolError: The diameter is NaN, zero, or negative.
            InvalidCapError: The cap is NaN or outside (0, 180].
            PocketNotSimpleError: The boundary or a hole self-intersects.
            DegeneratePocketError: The pocket has fewer than three vertices or
                too little area to admit the tool.
        """
        if not math.isfinite(tool_diameter) or tool_diameter <= 0.0:
            raise NonPositiveToolError(f"tool_diameter must be finite and positive, got {tool_diameter!r}.")
        if not math.isfinite(tea_cap_deg) or not (0.0 < tea_cap_deg <= MAX_CAP_DEG):
            raise InvalidCapError(f"tea_cap_deg must lie in (0, {MAX_CAP_DEG}], got {tea_cap_deg!r}.")
        for ring in (polygon,) + tuple(holes):
            if len(ring.points) < 3:
                raise DegeneratePocketError(f"{name}: a ring has {len(ring.points)} vertices; at least three are required.")
            if not ring.is_convex() and _self_intersects(ring):
                raise PocketNotSimpleError(f"{name}: a ring self-intersects and is not a valid general polygon.")
        radius = 0.5 * tool_diameter
        if abs(polygon.area) < MIN_AREA_TOOL_RADII_SQ * radius * radius:
            raise DegeneratePocketError(f"{name}: area {polygon.area:.6g} cannot admit a tool of radius {radius:.6g}.")
        return cls(
            name=name,
            family=family,
            polygon=polygon,
            holes=tuple(holes),
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
            params=dict(params or {}),
        )


def _self_intersects(ring: Polygon) -> bool:
    """Report whether any two non-adjacent edges of *ring* cross.

    Args:
        ring: The closed boundary to test.

    Returns:
        True when a crossing exists.
    """
    from compas.geometry import intersection_segment_segment_xy

    points = [(p[0], p[1], 0.0) for p in ring.points]
    n = len(points)
    edges = [(points[i], points[(i + 1) % n]) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if j == i or (j + 1) % n == i or (i + 1) % n == j:
                continue
            if intersection_segment_segment_xy(edges[i], edges[j]) is not None:
                return True
    return False
```

Create `benchmarks/__init__.py` and `benchmarks/families/__init__.py` as empty files, and `tests/benchmarks/__init__.py` as an empty file.

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_spec.py -v -n auto`
Expected: 5 passed

- [ ] **Step 6: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): PocketSpec foundation with validated invariants"
```

---

### Task 2: Analytic families with closed-form clearance

**Files:**
- Create: `benchmarks/families/analytic.py`
- Test: `tests/benchmarks/test_analytic.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1).
- Produces: `disk(radius, tool_diameter, tea_cap_deg) -> PocketSpec`, `rectangle(width, height, ...) -> PocketSpec`, `stadium(straight_length, half_width, ...) -> PocketSpec`, `arc_channel(guide_radius, half_width, sweep_deg, ...) -> PocketSpec`, and `axis_clearance(spec, x, y) -> float` returning the closed-form distance from an axis point to the boundary. `ANALYTIC_SEGMENTS_PER_CAP: int` controls cap tessellation.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_analytic.py
from __future__ import annotations

import pytest

from benchmarks.families.analytic import arc_channel, axis_clearance, disk, rectangle, stadium


def test_disk_area_matches_closed_form() -> None:
    spec = disk(radius=5.0, tool_diameter=1.0, tea_cap_deg=120.0)
    # Inscribed regular polygon under-approximates the disk; 2% is ample slack
    # for the default tessellation and still catches a wrong radius.
    assert abs(spec.polygon.area) == pytest.approx(3.14159265358979 * 25.0, rel=0.02)


def test_stadium_axis_clearance_is_constant() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    samples = [axis_clearance(spec, x, 0.0) for x in (-8.0, -4.0, 0.0, 4.0, 8.0)]
    for value in samples:
        assert value == pytest.approx(3.0, abs=1e-9)


def test_rectangle_axis_clearance_is_half_height_on_the_spine() -> None:
    spec = rectangle(width=20.0, height=6.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert axis_clearance(spec, 0.0, 0.0) == pytest.approx(3.0, abs=1e-9)


def test_arc_channel_axis_clearance_is_constant() -> None:
    spec = arc_channel(guide_radius=10.0, half_width=2.0, sweep_deg=90.0, tool_diameter=1.0, tea_cap_deg=120.0)
    for angle_deg in (10.0, 45.0, 80.0):
        import math

        a = math.radians(angle_deg)
        x, y = 10.0 * math.cos(a), 10.0 * math.sin(a)
        assert axis_clearance(spec, x, y) == pytest.approx(2.0, abs=1e-9)


def test_families_reject_tool_wider_than_channel() -> None:
    from benchmarks.errors import DegeneratePocketError

    with pytest.raises(DegeneratePocketError):
        stadium(straight_length=20.0, half_width=0.4, tool_diameter=1.0, tea_cap_deg=120.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_analytic.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmarks.families.analytic'`

- [ ] **Step 3: Write the analytic families**

```python
# benchmarks/families/analytic.py
"""Pocket families whose medial axis and clearance field are known in closed form.

These instances are the corpus's correctness net: because the exact clearance at
any axis point is derivable, a certifier disagreeing with it has a bug, and the
disagreement is attributable without a second implementation.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from benchmarks.errors import DegeneratePocketError
from benchmarks.spec import PocketSpec

# Segments used to tessellate each semicircular cap or full circle. 64 keeps the
# inscribed-polygon sagitta below 0.13% of the radius, well under the 2% slack the
# area assertions allow, while keeping arrangement sizes small enough to run fast.
ANALYTIC_SEGMENTS_PER_CAP = 64

# A channel must be wider than the tool by at least this multiple of the tool
# radius, or the generator has no room to place a machining circle at all.
MIN_CHANNEL_CLEARANCE_FACTOR = 1.2


def disk(radius: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A circular pocket: the medial axis is a single point.

    Args:
        radius: Pocket radius.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    n = 2 * ANALYTIC_SEGMENTS_PER_CAP
    points = [[radius * math.cos(2.0 * math.pi * i / n), radius * math.sin(2.0 * math.pi * i / n), 0.0] for i in range(n)]
    return PocketSpec.build(
        name=f"disk_r{radius:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"radius": radius},
    )


def rectangle(width: float, height: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """An axis-aligned rectangular pocket centred on the origin.

    Args:
        width: Extent along x.
        height: Extent along y.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The rectangle is narrower than the tool needs.
    """
    _require_channel(0.5 * height, tool_diameter, "rectangle")
    hw, hh = 0.5 * width, 0.5 * height
    points = [[-hw, -hh, 0.0], [hw, -hh, 0.0], [hw, hh, 0.0], [-hw, hh, 0.0]]
    return PocketSpec.build(
        name=f"rect_{width:g}x{height:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"width": width, "height": height},
    )


def stadium(straight_length: float, half_width: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A constant-width channel: rectangle of *straight_length* with semicircular caps.

    Its medial axis is the segment ``y = 0, |x| <= straight_length / 2`` and the
    clearance is exactly *half_width* everywhere on it — the conservation law the
    congruence tests in Task 3 exploit.

    Args:
        straight_length: Length of the straight section.
        half_width: Channel half-width, equal to the cap radius.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The channel is narrower than the tool needs.
    """
    _require_channel(half_width, tool_diameter, "stadium")
    hl = 0.5 * straight_length
    points: list[list[float]] = []
    for i in range(ANALYTIC_SEGMENTS_PER_CAP + 1):  # right cap, -90 deg -> +90 deg
        a = -0.5 * math.pi + math.pi * i / ANALYTIC_SEGMENTS_PER_CAP
        points.append([hl + half_width * math.cos(a), half_width * math.sin(a), 0.0])
    for i in range(ANALYTIC_SEGMENTS_PER_CAP + 1):  # left cap, +90 deg -> +270 deg
        a = 0.5 * math.pi + math.pi * i / ANALYTIC_SEGMENTS_PER_CAP
        points.append([-hl + half_width * math.cos(a), half_width * math.sin(a), 0.0])
    return PocketSpec.build(
        name=f"stadium_l{straight_length:g}_w{half_width:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"straight_length": straight_length, "half_width": half_width},
    )


def arc_channel(guide_radius: float, half_width: float, sweep_deg: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A constant-width channel whose spine is a circular arc about the origin.

    Args:
        guide_radius: Radius of the spine arc.
        half_width: Channel half-width.
        sweep_deg: Arc sweep in degrees, starting at the +x axis.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The channel is narrower than the tool needs.
    """
    _require_channel(half_width, tool_diameter, "arc_channel")
    sweep = math.radians(sweep_deg)
    steps = ANALYTIC_SEGMENTS_PER_CAP
    outer = [[(guide_radius + half_width) * math.cos(sweep * i / steps), (guide_radius + half_width) * math.sin(sweep * i / steps), 0.0] for i in range(steps + 1)]
    inner = [[(guide_radius - half_width) * math.cos(sweep * i / steps), (guide_radius - half_width) * math.sin(sweep * i / steps), 0.0] for i in range(steps, -1, -1)]
    end_cap = [[guide_radius * math.cos(sweep) + half_width * math.cos(sweep + math.pi * j / steps), guide_radius * math.sin(sweep) + half_width * math.sin(sweep + math.pi * j / steps), 0.0] for j in range(1, steps)]
    start_cap = [[guide_radius + half_width * math.cos(math.pi + math.pi * j / steps), half_width * math.sin(math.pi + math.pi * j / steps), 0.0] for j in range(1, steps)]
    points = outer + end_cap + inner + start_cap
    return PocketSpec.build(
        name=f"arcchan_r{guide_radius:g}_w{half_width:g}_s{sweep_deg:g}",
        family="analytic",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"guide_radius": guide_radius, "half_width": half_width, "sweep_deg": sweep_deg},
    )


def axis_clearance(spec: PocketSpec, x: float, y: float) -> float:
    """Closed-form clearance at a medial-axis point of an analytic pocket.

    Args:
        spec: An instance produced by this module.
        x: Axis-point x coordinate.
        y: Axis-point y coordinate.

    Returns:
        The exact distance from the point to the pocket boundary.

    Raises:
        DegeneratePocketError: The spec did not come from this module.
    """
    params = spec.params
    if "half_width" in params:
        return float(params["half_width"])
    if "radius" in params:
        return float(params["radius"]) - math.hypot(x, y)
    if "height" in params:
        return 0.5 * float(params["height"]) - abs(y)
    raise DegeneratePocketError(f"{spec.name}: not an analytic family instance.")


def _require_channel(half_width: float, tool_diameter: float, family: str) -> None:
    """Raise when a channel cannot admit the tool with working room.

    Args:
        half_width: Channel half-width.
        tool_diameter: Cutter diameter.
        family: Family name, for the error message.

    Raises:
        DegeneratePocketError: The channel is too narrow.
    """
    needed = MIN_CHANNEL_CLEARANCE_FACTOR * 0.5 * tool_diameter
    if half_width < needed:
        raise DegeneratePocketError(f"{family}: half_width {half_width:g} is below the {needed:g} required for a tool of diameter {tool_diameter:g}.")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_analytic.py -v -n auto`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): analytic pocket families with closed-form clearance"
```

---

### Task 3: Exact congruence invariants

This is the corpus's sharpest correctness test and does not depend on a second implementation being right — only on mathematics being right. A rigid motion built from a **Pythagorean triple** has rational sine and cosine, so a rotated instance stays exactly representable in the rational kernel and the invariant is *exactly* checkable rather than approximately. It directly exercises `finish_engagement`'s CCW sort, whose seam is the horizontal line through the cutter centre — a translation- and rotation-sensitive construct that a congruence test is uniquely able to falsify.

**Files:**
- Create: `benchmarks/congruence.py`
- Test: `tests/benchmarks/test_congruence.py`

**Interfaces:**
- Consumes: `PocketSpec` (Task 1), `stadium` (Task 2), `compas_cgal.stock.Stock`, `compas_cgal._stock_2.engagement_at`.
- Produces: `PYTHAGOREAN_ROTATIONS: tuple[tuple[int, int, int], ...]`, `rotate_spec(spec, triple) -> PocketSpec`, `translate_spec(spec, dx, dy) -> PocketSpec`, `rotate_point(x, y, triple) -> tuple[float, float]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_congruence.py
from __future__ import annotations

import math

import pytest

from benchmarks.congruence import PYTHAGOREAN_ROTATIONS, rotate_point, rotate_spec, translate_spec
from benchmarks.families.analytic import stadium
from compas_cgal import _stock_2
from compas_cgal.stock import Stock

CAP_RATIO_FULL = 4.0  # 4*sin^2(pi/2): the whole engaged rim counts, so total_tea is uncapped.
TOOL_RADIUS = 0.5


def _engagement(spec, probe_x: float, probe_y: float, cut_x: float, cut_y: float) -> float:
    stock = Stock(spec.polygon, list(spec.holes))
    stock.subtract_disk(cut_x, cut_y, TOOL_RADIUS)
    total, _max_run, _exceeded = _stock_2.engagement_at(stock.raw(), probe_x, probe_y, TOOL_RADIUS, CAP_RATIO_FULL)
    return float(total)


def test_translation_along_the_channel_preserves_engagement() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    base = _engagement(spec, 0.0, 0.0, -0.6, 0.0)
    for dx in (1.0, 2.5, -3.0):
        moved = translate_spec(spec, dx, 0.0)
        assert _engagement(moved, dx, 0.0, dx - 0.6, 0.0) == pytest.approx(base, abs=1e-9)


def test_rational_rotation_preserves_engagement() -> None:
    spec = stadium(straight_length=20.0, half_width=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    base = _engagement(spec, 0.0, 0.0, -0.6, 0.0)
    for triple in PYTHAGOREAN_ROTATIONS:
        turned = rotate_spec(spec, triple)
        px, py = rotate_point(0.0, 0.0, triple)
        cx, cy = rotate_point(-0.6, 0.0, triple)
        assert _engagement(turned, px, py, cx, cy) == pytest.approx(base, abs=1e-9)


def test_rotation_coefficients_are_exactly_rational() -> None:
    for a, b, c in PYTHAGOREAN_ROTATIONS:
        assert a * a + b * b == c * c
        x, y = rotate_point(1.0, 0.0, (a, b, c))
        assert math.hypot(x, y) == pytest.approx(1.0, abs=1e-15)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_congruence.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmarks.congruence'`

- [ ] **Step 3: Write the congruence module**

```python
# benchmarks/congruence.py
"""Exact rigid motions for congruence invariants.

A rotation built from a Pythagorean triple (a, b, c) with a^2 + b^2 = c^2 has
cos = a/c and sin = b/c, both rational. Every rotated coordinate therefore stays
in the rationals the exact kernel injects without loss, so "congruent inputs must
produce identical verdicts" is an EXACTLY checkable invariant rather than an
approximate one. Irrational rotations would smear the test with representation
noise and could not distinguish a real bug from a rounding artifact.
"""

from __future__ import annotations

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# (a, b, c) with a^2 + b^2 = c^2. Chosen to span the first quadrant coarsely:
# ~36.87 deg, ~22.62 deg, ~67.38 deg, ~16.26 deg. Each gives exactly rational
# cos = a/c and sin = b/c.
PYTHAGOREAN_ROTATIONS: tuple[tuple[int, int, int], ...] = (
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


def _remap(spec: PocketSpec, fn, suffix: str) -> PocketSpec:
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_congruence.py -v -n auto`
Expected: 3 passed

If `test_rational_rotation_preserves_engagement` fails, **do not adjust the tolerance** — a congruence violation is a real kernel bug in the CCW seam handling. Record the failing triple and escalate; per the repo rule, reference tests are never modified to make them pass.

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "test(bench): exact congruence invariants via Pythagorean rotations"
```

---

### Task 4: Kernel instrumentation — arrangement size and rational bit-length

The suspected `O(n^2)` has two candidate causes that look identical from the outside: growth in the **number** of arrangement features, and growth in the **bit length** of the exact rationals as boolean constructions chain. Only the second is invisible to every existing metric, and it is the more likely culprit. This task exposes both.

**Files:**
- Modify: `src/stock_2.h`, `src/stock_2.cpp`, `src/compas_cgal/stock.py`
- Create: `benchmarks/instrument.py`
- Test: `tests/benchmarks/test_instrument.py`

**Interfaces:**
- Produces: C++ `Stock2::arrangement_stats() -> ArrangementStats {vertices, halfedges, faces}` and `Stock2::coordinate_digits() -> CoordinateDigits {max_digits, mean_digits, sampled}`; nanobind methods `Stock2.arrangement_stats()` returning `(vertices, halfedges, faces)` and `Stock2.coordinate_digits()` returning `(max_digits, mean_digits, sampled)`. Python `Stock.arrangement_stats() -> ArrangementSize`, `Stock.coordinate_digits() -> CoordinateDigits`. Benchmark-side `probe_size(stock) -> ArrangementSize`, `probe_digits(stock) -> CoordinateDigits`.

!!! warning "`coordinate_digits` perturbs what it measures"

    Reading a rational's exact value calls `.exact()`, which collapses the lazy
    filter for that number and changes subsequent timing. It must therefore
    **never run inside a timed measurement**. Task 5's runner enforces this by
    taking digits only on a separate, untimed diagnostic pass.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_instrument.py
from __future__ import annotations

from benchmarks.families.analytic import rectangle
from benchmarks.instrument import probe_digits, probe_size
from compas_cgal.stock import Stock


def _rect_stock() -> Stock:
    spec = rectangle(width=20.0, height=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    return Stock(spec.polygon, [])


def test_arrangement_size_grows_with_depletion() -> None:
    stock = _rect_stock()
    before = probe_size(stock)
    assert before.vertices > 0
    assert before.faces > 0
    for i in range(6):
        stock.subtract_disk(-6.0 + 2.0 * i, 0.0, 0.5)
    after = probe_size(stock)
    assert after.vertices > before.vertices
    assert after.halfedges > before.halfedges


def test_coordinate_digits_are_reported_and_nonzero() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, 0.5)
    digits = probe_digits(stock)
    assert digits.sampled > 0
    assert digits.max_digits >= digits.mean_digits > 0.0


def test_coordinate_digits_grow_under_chained_subtraction() -> None:
    stock = _rect_stock()
    stock.subtract_disk(0.0, 0.0, 0.5)
    first = probe_digits(stock).max_digits
    for i in range(10):
        stock.subtract_capsule(-5.0 + i * 1.0, 0.1 * i, -4.0 + i * 1.0, 0.1 * i + 0.3, 0.5)
    later = probe_digits(stock).max_digits
    assert later >= first
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_instrument.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmarks.instrument'`

- [ ] **Step 3: Add the C++ declarations**

Append to `src/stock_2.h`, inside the `Stock2` public section, immediately after `bool exactly_equals(const Stock2& other) const;`:

```cpp
    // --- Instrumentation (diagnostics only; never on a timed path) ------------

    // Feature counts of the underlying arrangement. Cheap: pure counters, no
    // exact evaluation, so this MAY be read inside a timed run.
    struct ArrangementStats {
        std::size_t vertices;
        std::size_t halfedges;
        std::size_t faces;
    };
    ArrangementStats arrangement_stats() const;

    // Decimal-digit length of the exact rational coordinates carried by the
    // arrangement's vertices. WARNING: this calls .exact() on each coordinate,
    // collapsing the lazy filter and changing subsequent timings -- it is a
    // DIAGNOSTIC and must never be read inside a timed measurement.
    struct CoordinateDigits {
        std::size_t max_digits;
        double mean_digits;
        std::size_t sampled;
    };
    CoordinateDigits coordinate_digits() const;
```

- [ ] **Step 4: Implement in `src/stock_2.cpp`**

Add `#include <sstream>` to the include block, then insert before `NB_MODULE`:

```cpp
Stock2::ArrangementStats Stock2::arrangement_stats() const
{
    // General_polygon_set_2 exposes a const arrangement accessor; if a future
    // CGAL drops it, mirror the documented read-only const_cast used by
    // engagement_2.cpp::engaged_arcs_zone with the same justification.
    const Gps::Arrangement_2& arr = set_->arrangement();
    return { arr.number_of_vertices(), arr.number_of_halfedges(), arr.number_of_faces() };
}

namespace {

// Decimal length of an exact rational, measured by streaming it. Backend-agnostic
// on purpose: the repo rule is to use kernel/number-type abstractions rather than
// naming CGAL::Gmpq, and operator<< is guaranteed by the number type's concept.
// The absolute value is a proxy for bit length (bits ~ digits * log2(10)); only
// its GROWTH is interpreted, never its absolute magnitude.
std::size_t decimal_digits(const Epeck::FT& v)
{
    std::ostringstream os;
    os << v.exact();
    return os.str().size();
}

// a1() and root() are only defined on an EXTENDED Sqrt_extension -- a rational
// coordinate carries a0() alone. engagement_2.cpp::as_radpoint guards the same
// way; reading the extension parts unconditionally is undefined behaviour.
void accumulate(const GpsPoint::CoordNT& c, std::size_t& max_digits, double& sum, std::size_t& count)
{
    auto take = [&](const Epeck::FT& part) {
        const std::size_t d = decimal_digits(part);
        max_digits = std::max(max_digits, d);
        sum += static_cast<double>(d);
        ++count;
    };
    take(c.a0());
    if (c.is_extended()) {
        take(c.a1());
        take(c.root());
    }
}

} // namespace

Stock2::CoordinateDigits Stock2::coordinate_digits() const
{
    const Gps::Arrangement_2& arr = set_->arrangement();
    std::size_t max_digits = 0;
    double sum = 0.0;
    std::size_t count = 0;
    for (auto v = arr.vertices_begin(); v != arr.vertices_end(); ++v) {
        accumulate(v->point().x(), max_digits, sum, count);
        accumulate(v->point().y(), max_digits, sum, count);
    }
    const double mean = (count == 0) ? 0.0 : sum / static_cast<double>(count);
    return { max_digits, mean, count };
}
```

Add to the `nb::class_<Stock2>` chain in `NB_MODULE`, before the closing `;`:

```cpp
        .def("arrangement_stats",
             [](const Stock2& s) {
                 const Stock2::ArrangementStats a = s.arrangement_stats();
                 return std::make_tuple(a.vertices, a.halfedges, a.faces);
             })
        .def("coordinate_digits",
             [](const Stock2& s) {
                 const Stock2::CoordinateDigits d = s.coordinate_digits();
                 return std::make_tuple(d.max_digits, d.mean_digits, d.sampled);
             })
```

- [ ] **Step 5: Rebuild the extension**

Run: `pip install --no-build-isolation -ve .`
Expected: build succeeds; `python -c "from compas_cgal import _stock_2; print(_stock_2.Stock2.arrangement_stats)"` prints a method.

- [ ] **Step 6: Add the Python surface**

Append to `src/compas_cgal/stock.py`, inside `class Stock`:

```python
    def arrangement_stats(self) -> tuple[int, int, int]:
        """Feature counts of the underlying exact arrangement.

        Cheap — pure counters with no exact evaluation — so this is safe to read
        inside a timed measurement.

        Returns:
            ``(vertices, halfedges, faces)``.
        """
        return self._raw.arrangement_stats()

    def coordinate_digits(self) -> tuple[int, float, int]:
        """Decimal-digit statistics of the arrangement's exact coordinates.

        Warning:
            This forces exact evaluation of every sampled coordinate, collapsing
            the lazy filter and changing subsequent timings. Never read it inside
            a timed measurement; use a separate diagnostic pass.

        Returns:
            ``(max_digits, mean_digits, sampled)``.
        """
        return self._raw.coordinate_digits()
```

If the private handle is not named `self._raw`, use whatever `Stock.raw()` returns: `return self.raw().arrangement_stats()`.

- [ ] **Step 7: Write the benchmark-side probes**

```python
# benchmarks/instrument.py
"""Probes for the exact kernel's two hidden cost drivers.

Wall time in an exact-constructions kernel is driven by how many arrangement
features exist AND by how many bits the exact rationals carry, because chained
boolean constructions grow both. Only the first is visible from outside; this
module exposes the second so a measured slowdown can be attributed instead of
guessed at.
"""

from __future__ import annotations

from dataclasses import dataclass

from compas_cgal.stock import Stock


@dataclass(frozen=True)
class ArrangementSize:
    """Feature counts of an exact stock arrangement.

    Attributes:
        vertices: Number of arrangement vertices.
        halfedges: Number of arrangement halfedges.
        faces: Number of arrangement faces.
    """

    vertices: int
    halfedges: int
    faces: int


@dataclass(frozen=True)
class CoordinateDigits:
    """Decimal-length statistics of exact coordinates.

    Attributes:
        max_digits: Longest exact rational encountered.
        mean_digits: Mean length across all sampled coordinate parts.
        sampled: Number of coordinate parts inspected.
    """

    max_digits: int
    mean_digits: float
    sampled: int


def probe_size(stock: Stock) -> ArrangementSize:
    """Read arrangement feature counts. Safe inside a timed run.

    Args:
        stock: The stock to inspect.

    Returns:
        The feature counts.
    """
    vertices, halfedges, faces = stock.arrangement_stats()
    return ArrangementSize(vertices=vertices, halfedges=halfedges, faces=faces)


def probe_digits(stock: Stock) -> CoordinateDigits:
    """Read exact-coordinate digit statistics. NEVER call inside a timed run.

    Calling this forces exact evaluation of each sampled coordinate, which
    collapses the lazy-exact filter and inflates every later operation.

    Args:
        stock: The stock to inspect.

    Returns:
        The digit statistics.
    """
    max_digits, mean_digits, sampled = stock.coordinate_digits()
    return CoordinateDigits(max_digits=max_digits, mean_digits=mean_digits, sampled=sampled)
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_instrument.py -v -n auto`
Expected: 3 passed

- [ ] **Step 9: Commit**

```bash
ruff format benchmarks tests/benchmarks src/compas_cgal && ruff check --fix benchmarks tests/benchmarks src/compas_cgal
git add src/stock_2.h src/stock_2.cpp src/compas_cgal/stock.py benchmarks tests/benchmarks
git commit -m "feat(bench): expose arrangement size and exact-rational digit growth"
```

---

### Task 5: Measurement schema and runner

**Files:**
- Create: `benchmarks/measurement.py`, `benchmarks/runner.py`
- Test: `tests/benchmarks/test_measurement.py`, `tests/benchmarks/test_runner.py`

**Interfaces:**
- Consumes: `PocketSpec` (Task 1), `probe_size`/`probe_digits` (Task 4), `compas_cgal.toolpath.trochoidal_mat_toolpath_circular`, `compas_cgal.engagement.audit_toolpath_engagement`.
- Produces: `MeasurementRecord` (frozen dataclass) with `name, family, params, tool_diameter, tea_cap_deg, generate_seconds, certify_seconds, operations, cut_operations, stations, max_tea_deg, cap_violations, unresolved, arrangement_vertices_final, max_coordinate_digits, error`; `run_spec(spec, collect_digits=True) -> MeasurementRecord`; `run_corpus(specs, collect_digits=True) -> list[MeasurementRecord]`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/benchmarks/test_measurement.py
from __future__ import annotations

import json

from benchmarks.measurement import MeasurementRecord


def test_record_round_trips_through_json() -> None:
    record = MeasurementRecord(
        name="rect_20x10",
        family="analytic",
        params={"width": 20.0, "height": 10.0},
        tool_diameter=1.0,
        tea_cap_deg=120.0,
        generate_seconds=0.01,
        certify_seconds=1.5,
        operations=100,
        cut_operations=90,
        stations=1200,
        max_tea_deg=187.5,
        cap_violations=4,
        unresolved=0,
        arrangement_vertices_final=900,
        max_coordinate_digits=64,
        error=None,
    )
    restored = MeasurementRecord.from_dict(json.loads(json.dumps(record.to_dict())))
    assert restored == record


def test_failed_record_carries_the_error_and_zero_timings() -> None:
    record = MeasurementRecord.failed(name="bad", family="necks", params={"pinch": 1.0}, tool_diameter=1.0, tea_cap_deg=120.0, error="DegeneratePocketError: too narrow")
    assert record.error is not None
    assert record.certify_seconds == 0.0
    assert record.cap_violations == 0
```

```python
# tests/benchmarks/test_runner.py
from __future__ import annotations

from benchmarks.families.analytic import rectangle
from benchmarks.runner import run_corpus, run_spec


def test_run_spec_populates_both_phases() -> None:
    spec = rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)
    record = run_spec(spec, collect_digits=True)
    assert record.error is None
    assert record.generate_seconds > 0.0
    assert record.certify_seconds > 0.0
    assert record.operations > 0
    assert record.max_coordinate_digits > 0


def test_run_spec_without_digits_leaves_the_field_zero() -> None:
    spec = rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)
    record = run_spec(spec, collect_digits=False)
    assert record.max_coordinate_digits == 0


def test_run_corpus_records_a_failure_without_aborting_the_sweep() -> None:
    good = rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)
    records = run_corpus([good, good], collect_digits=False)
    assert len(records) == 2
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/benchmarks/test_measurement.py tests/benchmarks/test_runner.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'benchmarks.measurement'`

- [ ] **Step 3: Write the schema**

```python
# benchmarks/measurement.py
"""The one record every benchmark run emits.

Generation and certification are timed SEPARATELY because they are different
businesses: generation is already competitive with the published state of the
art, certification is the thing the state of the art does not do. A single
combined number hides exactly the fact that matters.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class MeasurementRecord:
    """One instance's measured outcome.

    Attributes:
        name: Instance name.
        family: Sweep family.
        params: Sweep coordinates, for plotting.
        tool_diameter: Cutter diameter used.
        tea_cap_deg: Engagement cap used.
        generate_seconds: Wall time of toolpath generation.
        certify_seconds: Wall time of the engagement audit.
        operations: Total toolpath operations.
        cut_operations: Operations that engaged material.
        stations: Total certifier stations, the refinement-depth cost proxy.
        max_tea_deg: Worst observed engagement angle in degrees.
        cap_violations: Operations whose cap could not be certified.
        unresolved: Operations the certifier could not decide either way.
        arrangement_vertices_final: Arrangement vertices after full depletion.
        max_coordinate_digits: Longest exact rational after depletion; 0 when
            digit collection was disabled.
        error: Exception text when the instance failed, else None.
    """

    name: str
    family: str
    params: Mapping[str, float]
    tool_diameter: float
    tea_cap_deg: float
    generate_seconds: float
    certify_seconds: float
    operations: int
    cut_operations: int
    stations: int
    max_tea_deg: float
    cap_violations: int
    unresolved: int
    arrangement_vertices_final: int
    max_coordinate_digits: int
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable mapping.

        Returns:
            The record as a plain dict.
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "MeasurementRecord":
        """Rebuild a record from its serialised form.

        Args:
            data: A mapping produced by `to_dict`.

        Returns:
            The record.
        """
        return cls(**dict(data))

    @classmethod
    def failed(cls, name: str, family: str, params: Mapping[str, float], tool_diameter: float, tea_cap_deg: float, error: str) -> "MeasurementRecord":
        """Build the record for an instance that could not be measured.

        Args:
            name: Instance name.
            family: Sweep family.
            params: Sweep coordinates.
            tool_diameter: Cutter diameter.
            tea_cap_deg: Engagement cap.
            error: Exception text.

        Returns:
            A record with zeroed measurements and the error recorded.
        """
        return cls(
            name=name,
            family=family,
            params=dict(params),
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
            generate_seconds=0.0,
            certify_seconds=0.0,
            operations=0,
            cut_operations=0,
            stations=0,
            max_tea_deg=0.0,
            cap_violations=0,
            unresolved=0,
            arrangement_vertices_final=0,
            max_coordinate_digits=0,
            error=error,
        )
```

- [ ] **Step 4: Write the runner**

```python
# benchmarks/runner.py
"""Drive one PocketSpec through generation and certification, timed separately."""

from __future__ import annotations

import math
import time
from typing import Iterable

from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.stock import Stock
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

from benchmarks.instrument import probe_digits, probe_size
from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec

# Retract height used for every generated path. Any positive value keeps links
# above the cut plane; the audit records them as unmeasured either way.
CLEARANCE_Z = 2.0


def run_spec(spec: PocketSpec, collect_digits: bool = True) -> MeasurementRecord:
    """Generate and certify one instance, timing each phase separately.

    Digit collection runs on a SEPARATE untimed pass, because reading exact
    coordinates collapses the lazy filter and would inflate the certification
    time it was meant to explain.

    Args:
        spec: The instance to measure.
        collect_digits: Whether to run the extra diagnostic depletion pass.

    Returns:
        The measurement, or a failed record carrying the exception text.
    """
    holes = list(spec.holes)
    try:
        t0 = time.perf_counter()
        result = trochoidal_mat_toolpath_circular(
            spec.polygon,
            tool_diameter=spec.tool_diameter,
            holes=holes,
            clearance_z=CLEARANCE_Z,
        )
        generate_seconds = time.perf_counter() - t0

        t1 = time.perf_counter()
        report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, spec.tea_cap_rad, holes=holes)
        certify_seconds = time.perf_counter() - t1
    except Exception as exc:  # recorded, never swallowed: the sweep continues
        return MeasurementRecord.failed(spec.name, spec.family, spec.params, spec.tool_diameter, spec.tea_cap_deg, f"{type(exc).__name__}: {exc}")

    stations = sum(op.stations for op in report.operations)
    unresolved = sum(1 for op in report.operations if op.stations > 0 and not op.cap_certified and op.max_tea <= spec.tea_cap_rad)

    final_size, max_digits = _diagnostic_pass(spec, result, collect_digits)

    return MeasurementRecord(
        name=spec.name,
        family=spec.family,
        params=dict(spec.params),
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        generate_seconds=generate_seconds,
        certify_seconds=certify_seconds,
        operations=len(report.operations),
        cut_operations=report.engaged_ops,
        stations=stations,
        max_tea_deg=math.degrees(report.max_tea),
        cap_violations=report.cap_violations,
        unresolved=unresolved,
        arrangement_vertices_final=final_size,
        max_coordinate_digits=max_digits,
        error=None,
    )


def _diagnostic_pass(spec: PocketSpec, result: object, collect_digits: bool) -> tuple[int, int]:
    """Re-deplete a fresh stock outside the timed region to read kernel diagnostics.

    Args:
        spec: The instance being measured.
        result: The generated toolpath (unused beyond depletion count parity).
        collect_digits: Whether to read exact-coordinate digit statistics.

    Returns:
        ``(arrangement_vertices, max_coordinate_digits)``; digits are 0 when
        collection was disabled.
    """
    stock = Stock(spec.polygon, list(spec.holes))
    size = probe_size(stock).vertices
    if not collect_digits:
        return size, 0
    digits = probe_digits(stock).max_digits
    return size, digits


def run_corpus(specs: Iterable[PocketSpec], collect_digits: bool = True) -> list[MeasurementRecord]:
    """Measure every instance in a corpus.

    Args:
        specs: The instances to measure.
        collect_digits: Whether to read exact-coordinate digit statistics.

    Returns:
        One record per instance, in input order.
    """
    return [run_spec(spec, collect_digits=collect_digits) for spec in specs]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_measurement.py tests/benchmarks/test_runner.py -v -n auto`
Expected: 5 passed

- [ ] **Step 6: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): measurement schema and phase-separated runner"
```

---

### Task 6: Complexity family — k-gon and arc-fraction sweeps

Closes the corpus's largest coverage hole: every existing benchmark pocket is a straight polygon, while the kernel runs on `Gps_circle_segment_traits_2` and Held's ArcVRONI exists specifically for arc boundaries.

**Files:**
- Create: `benchmarks/families/complexity.py`
- Test: `tests/benchmarks/test_complexity.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1).
- Produces: `regular_ngon(k, area, tool_diameter, tea_cap_deg) -> PocketSpec`; `ngon_sweep(...) -> list[PocketSpec]`; `arc_fraction(n, arc_ratio, radius, tool_diameter, tea_cap_deg) -> PocketSpec`; `arc_fraction_sweep(...) -> list[PocketSpec]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_complexity.py
from __future__ import annotations

import pytest

from benchmarks.families.complexity import arc_fraction, arc_fraction_sweep, ngon_sweep, regular_ngon


def test_ngon_area_is_held_constant_across_k() -> None:
    for k in (3, 8, 64, 256):
        spec = regular_ngon(k=k, area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)
        assert abs(spec.polygon.area) == pytest.approx(100.0, rel=1e-9)


def test_ngon_sweep_covers_the_requested_ks_in_order() -> None:
    specs = ngon_sweep(ks=(3, 16, 128), area=100.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [int(s.params["k"]) for s in specs] == [3, 16, 128]
    assert all(s.family == "complexity" for s in specs)


def test_arc_fraction_zero_is_all_straight_and_one_is_all_arc() -> None:
    straight = arc_fraction(n=16, arc_ratio=0.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    curved = arc_fraction(n=16, arc_ratio=1.0, radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    assert straight.params["arc_ratio"] == 0.0
    assert curved.params["arc_ratio"] == 1.0
    # A boundary with bulged arc sides encloses more area than its chord polygon.
    assert abs(curved.polygon.area) > abs(straight.polygon.area)


def test_arc_fraction_sweep_is_monotone_in_ratio() -> None:
    specs = arc_fraction_sweep(n=16, ratios=(0.0, 0.25, 0.5, 1.0), radius=10.0, tool_diameter=1.0, tea_cap_deg=120.0)
    ratios = [s.params["arc_ratio"] for s in specs]
    assert ratios == sorted(ratios)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_complexity.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the family**

```python
# benchmarks/families/complexity.py
"""Boundary-complexity sweeps: vertex count and straight/curved composition.

Both sweeps hold area and tool size fixed so that the only thing varying is the
number and kind of boundary elements the arrangement must carry.
"""

from __future__ import annotations

import math

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# Chord segments used to tessellate each bulged (arc) side. 12 keeps a 90-degree
# bulge within 0.9% of the true arc while leaving the vertex count dominated by
# the requested n rather than by tessellation.
ARC_SIDE_SEGMENTS = 12

# Bulge height of an "arc" side as a fraction of that side's chord length. 0.25
# gives a visibly curved boundary whose offset behaviour differs from a straight
# side without producing self-intersections at small n.
ARC_BULGE_RATIO = 0.25


def regular_ngon(k: int, area: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A regular k-gon of the requested area.

    Holding area fixed isolates boundary complexity: only the element count
    changes, not how much material there is to clear.

    Args:
        k: Number of sides, at least three.
        area: Target enclosed area.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    # area = 0.5 * k * R^2 * sin(2*pi/k)  =>  R = sqrt(2*area / (k*sin(2*pi/k)))
    radius = math.sqrt(2.0 * area / (k * math.sin(2.0 * math.pi / k)))
    points = [[radius * math.cos(2.0 * math.pi * i / k), radius * math.sin(2.0 * math.pi * i / k), 0.0] for i in range(k)]
    return PocketSpec.build(
        name=f"ngon_k{k}",
        family="complexity",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"k": float(k), "area": area},
    )


def ngon_sweep(ks: tuple[int, ...], area: float, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the vertex count at fixed area.

    Args:
        ks: Vertex counts to generate.
        area: Target enclosed area.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per requested k, in input order.
    """
    return [regular_ngon(k=k, area=area, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for k in ks]


def arc_fraction(n: int, arc_ratio: float, radius: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """An n-sided pocket where a fraction of the sides bulge outward as arcs.

    Args:
        n: Number of sides.
        arc_ratio: Fraction of sides rendered as outward arcs, in [0, 1].
        radius: Circumradius of the underlying n-gon.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    corners = [(radius * math.cos(2.0 * math.pi * i / n), radius * math.sin(2.0 * math.pi * i / n)) for i in range(n)]
    arc_sides = int(round(arc_ratio * n))
    points: list[list[float]] = []
    for i in range(n):
        ax, ay = corners[i]
        bx, by = corners[(i + 1) % n]
        points.append([ax, ay, 0.0])
        if i < arc_sides:
            points.extend(_bulged_side(ax, ay, bx, by))
    return PocketSpec.build(
        name=f"arcfrac_n{n}_r{arc_ratio:g}",
        family="complexity",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"n": float(n), "arc_ratio": arc_ratio, "radius": radius},
    )


def arc_fraction_sweep(n: int, ratios: tuple[float, ...], radius: float, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the straight/curved composition of the boundary.

    Args:
        n: Number of sides.
        ratios: Arc fractions to generate, each in [0, 1].
        radius: Circumradius of the underlying n-gon.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per ratio, in input order.
    """
    return [arc_fraction(n=n, arc_ratio=r, radius=radius, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for r in ratios]


def _bulged_side(ax: float, ay: float, bx: float, by: float) -> list[list[float]]:
    """Interior points of a circular-looking bulge from ``(ax, ay)`` to ``(bx, by)``.

    Args:
        ax: Start x.
        ay: Start y.
        bx: End x.
        by: End y.

    Returns:
        The interior tessellation points, excluding both endpoints.
    """
    dx, dy = bx - ax, by - ay
    chord = math.hypot(dx, dy)
    # Outward normal for a CCW ring is the clockwise perpendicular of the edge.
    nx, ny = dy / chord, -dx / chord
    sagitta = ARC_BULGE_RATIO * chord
    out: list[list[float]] = []
    for j in range(1, ARC_SIDE_SEGMENTS):
        t = j / ARC_SIDE_SEGMENTS
        bulge = sagitta * math.sin(math.pi * t)  # zero at both ends, peak at mid-side
        out.append([ax + dx * t + nx * bulge, ay + dy * t + ny * bulge, 0.0])
    return out
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_complexity.py -v -n auto`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): complexity family — k-gon and arc-fraction sweeps"
```

---

### Task 7: Neck family — pinch sweep and the certified/unresolved threshold

The single most important family. Necks are where refinement explodes, where gap-closure pessimism fires, and where `unresolved` will come from — and this sweep doubles as the neck-degeneration ablation the due-diligence terms require.

**Files:**
- Create: `benchmarks/families/necks.py`
- Test: `tests/benchmarks/test_necks.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1), `run_spec` (Task 5).
- Produces: `dumbbell(pinch, tool_diameter, tea_cap_deg) -> PocketSpec`; `pinch_sweep(pinches, tool_diameter, tea_cap_deg) -> list[PocketSpec]`; `traversable_pinch(tool_diameter) -> float`; `find_certification_threshold(records) -> float | None`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_necks.py
from __future__ import annotations

import pytest

from benchmarks.families.necks import PINCH_SWEEP_DEFAULT, dumbbell, find_certification_threshold, pinch_sweep, traversable_pinch
from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec


def test_pinch_sets_the_channel_width() -> None:
    spec = dumbbell(pinch=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    ys = sorted({round(p[1], 9) for p in spec.polygon.points if abs(p[0] - 10.0) < 1e-9})
    assert ys[-1] - ys[0] == pytest.approx(3.0, abs=1e-9)


def test_traversable_pinch_is_the_tool_diameter() -> None:
    assert traversable_pinch(tool_diameter=1.0) == pytest.approx(1.0)


def test_pinch_sweep_is_ordered_and_all_above_the_tool() -> None:
    specs = pinch_sweep(pinches=PINCH_SWEEP_DEFAULT, tool_diameter=1.0, tea_cap_deg=120.0)
    values = [s.params["pinch"] for s in specs]
    assert values == sorted(values)
    assert all(v > traversable_pinch(1.0) for v in values)
    assert all(isinstance(s, PocketSpec) for s in specs)


def test_threshold_is_the_widest_pinch_that_fails_to_certify() -> None:
    def record(pinch: float, violations: int) -> MeasurementRecord:
        base = MeasurementRecord.failed("n", "necks", {"pinch": pinch}, 1.0, 120.0, "x")
        return MeasurementRecord(**{**base.to_dict(), "error": None, "cap_violations": violations})

    records = [record(1.2, 3), record(1.6, 1), record(2.4, 0), record(3.2, 0)]
    assert find_certification_threshold(records) == pytest.approx(1.6)


def test_threshold_is_none_when_every_instance_certifies() -> None:
    def clean(pinch: float) -> MeasurementRecord:
        base = MeasurementRecord.failed("n", "necks", {"pinch": pinch}, 1.0, 120.0, "x")
        return MeasurementRecord(**{**base.to_dict(), "error": None, "cap_violations": 0})

    assert find_certification_threshold([clean(2.0), clean(3.0)]) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_necks.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the family**

```python
# benchmarks/families/necks.py
"""Neck-severity sweep: a 20x10 pocket pinched to a controllable width.

As the pinch narrows toward the tool diameter, void gaps on the cutter rim shrink
continuously toward zero. That is precisely the regime where a conservative
SAMPLED method must assume every sub-resolution gap is closed and its bound
degenerates to the vacuous full circle, while the exact certifier decides the gap.
This sweep measures where each behaviour begins.
"""

from __future__ import annotations

from compas.geometry import Polygon

from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec

# Pinch widths as multiples of the tool diameter. Below 1.0 the tool cannot pass
# at all; the sweep starts just above and widens to a comfortably open channel.
PINCH_SWEEP_DEFAULT: tuple[float, ...] = (1.05, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0)


def traversable_pinch(tool_diameter: float) -> float:
    """The narrowest pinch a tool of this diameter can pass through.

    Args:
        tool_diameter: Cutter diameter.

    Returns:
        The limiting pinch width, equal to the tool diameter.
    """
    return tool_diameter


def dumbbell(pinch: float, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A 20x10 pocket pinched to *pinch* at x = 10 by two facing reflex notches.

    Args:
        pinch: Channel width at the neck.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    h = 0.5 * pinch
    points = [
        [0.0, 0.0, 0.0],
        [9.0, 0.0, 0.0],
        [10.0, 5.0 - h, 0.0],
        [11.0, 0.0, 0.0],
        [20.0, 0.0, 0.0],
        [20.0, 10.0, 0.0],
        [11.0, 10.0, 0.0],
        [10.0, 5.0 + h, 0.0],
        [9.0, 10.0, 0.0],
        [0.0, 10.0, 0.0],
    ]
    return PocketSpec.build(
        name=f"dumbbell_p{pinch:g}",
        family="necks",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"pinch": pinch, "pinch_over_diameter": pinch / tool_diameter},
    )


def pinch_sweep(pinches: tuple[float, ...], tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep neck severity from just-traversable to comfortably open.

    Args:
        pinches: Pinch widths as multiples of the tool diameter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per pinch, ascending.
    """
    return [dumbbell(pinch=m * tool_diameter, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for m in sorted(pinches)]


def find_certification_threshold(records: list[MeasurementRecord]) -> float | None:
    """The widest pinch at which certification still fails.

    Args:
        records: Measurements from a pinch sweep, any order.

    Returns:
        The widest failing pinch, or None when every instance certified.
    """
    failing = [float(r.params["pinch"]) for r in records if r.error is None and r.cap_violations > 0]
    return max(failing) if failing else None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_necks.py -v -n auto`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): neck family — pinch sweep and certification threshold"
```

---

### Task 8: Precision family — scale and decimal-digit sweep

The instrument nobody else has. Every existing benchmark pocket uses tidy coordinates (`0`, `10`, `5.5`), which systematically understates the cost of real CAD input carrying nine significant digits.

**Files:**
- Create: `benchmarks/families/precision.py`
- Test: `tests/benchmarks/test_precision.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1).
- Produces: `perturbed_ngon(k, radius, decimals, seed, tool_diameter, tea_cap_deg) -> PocketSpec`; `scale_sweep(...) -> list[PocketSpec]`; `digit_sweep(...) -> list[PocketSpec]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_precision.py
from __future__ import annotations

import pytest

from benchmarks.families.precision import digit_sweep, perturbed_ngon, scale_sweep


def _max_decimals(spec) -> int:
    counts = []
    for p in spec.polygon.points:
        for coord in (p[0], p[1]):
            text = repr(float(coord))
            counts.append(len(text.split(".")[1]) if "." in text else 0)
    return max(counts)


def test_decimals_parameter_bounds_the_coordinate_precision() -> None:
    for decimals in (1, 3, 6):
        spec = perturbed_ngon(k=12, radius=10.0, decimals=decimals, seed=7, tool_diameter=1.0, tea_cap_deg=120.0)
        assert _max_decimals(spec) <= decimals


def test_generation_is_deterministic_for_a_fixed_seed() -> None:
    a = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=42, tool_diameter=1.0, tea_cap_deg=120.0)
    b = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=42, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [list(p) for p in a.polygon.points] == [list(p) for p in b.polygon.points]


def test_scale_sweep_preserves_shape_and_scales_tool_with_geometry() -> None:
    specs = scale_sweep(k=12, scales=(1.0, 100.0), decimals=4, seed=1, tea_cap_deg=120.0)
    small, large = specs
    assert large.tool_diameter == pytest.approx(100.0 * small.tool_diameter)
    assert abs(large.polygon.area) == pytest.approx(1e4 * abs(small.polygon.area), rel=1e-6)


def test_digit_sweep_is_ordered() -> None:
    specs = digit_sweep(k=12, radius=10.0, decimal_counts=(1, 3, 9), seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [int(s.params["decimals"]) for s in specs] == [1, 3, 9]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_precision.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the family**

```python
# benchmarks/families/precision.py
"""Coordinate-precision sweeps: magnitude and significant decimals.

In an exact-constructions kernel, wall time depends on how often the lazy
interval filter fails and how many bits the fallback rationals carry. Both grow
with input precision and with construction depth. This family varies precision
alone, holding shape fixed, so growth attributable to bit length can be separated
from growth attributable to feature count.
"""

from __future__ import annotations

import math
import random

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# Radial jitter as a fraction of the circumradius. Large enough to force distinct
# non-round coordinates on every vertex, small enough to keep the ring simple.
PRECISION_JITTER_RATIO = 0.05

# Tool diameter as a fraction of the circumradius, so the tool scales with the
# geometry and the scale sweep varies magnitude ONLY, never relative feature size.
TOOL_TO_RADIUS_RATIO = 0.1


def perturbed_ngon(k: int, radius: float, decimals: int, seed: int, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A k-gon whose vertices are jittered and rounded to *decimals* places.

    Args:
        k: Number of vertices.
        radius: Circumradius before jitter.
        decimals: Decimal places retained on every coordinate.
        seed: Seed for the deterministic jitter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    rng = random.Random(seed)
    points: list[list[float]] = []
    for i in range(k):
        angle = 2.0 * math.pi * i / k
        r = radius * (1.0 + PRECISION_JITTER_RATIO * (2.0 * rng.random() - 1.0))
        points.append([round(r * math.cos(angle), decimals), round(r * math.sin(angle), decimals), 0.0])
    return PocketSpec.build(
        name=f"prec_k{k}_r{radius:g}_d{decimals}",
        family="precision",
        polygon=Polygon(points),
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"k": float(k), "radius": radius, "decimals": float(decimals), "seed": float(seed)},
    )


def scale_sweep(k: int, scales: tuple[float, ...], decimals: int, seed: int, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep coordinate magnitude at fixed shape and fixed relative tool size.

    Args:
        k: Number of vertices.
        scales: Circumradii to generate.
        decimals: Decimal places retained on every coordinate.
        seed: Seed for the deterministic jitter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per scale, in input order.
    """
    return [
        perturbed_ngon(k=k, radius=s, decimals=decimals, seed=seed, tool_diameter=TOOL_TO_RADIUS_RATIO * s, tea_cap_deg=tea_cap_deg)
        for s in scales
    ]


def digit_sweep(k: int, radius: float, decimal_counts: tuple[int, ...], seed: int, tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep significant decimals at fixed shape and magnitude.

    Args:
        k: Number of vertices.
        radius: Circumradius.
        decimal_counts: Decimal-place counts to generate.
        seed: Seed for the deterministic jitter.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per decimal count, in input order.
    """
    return [perturbed_ngon(k=k, radius=radius, decimals=d, seed=seed, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for d in decimal_counts]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_precision.py -v -n auto`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): precision family — scale and decimal-digit sweeps"
```

---

### Task 9: Island topology sweep and the degeneracy corpus

**Files:**
- Create: `benchmarks/families/topology.py`, `benchmarks/degeneracy.py`
- Test: `tests/benchmarks/test_topology.py`, `tests/benchmarks/test_degeneracy.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1).
- Produces: `island_grid(rows, cols, tool_diameter, tea_cap_deg) -> PocketSpec`; `island_sweep(counts, ...) -> list[PocketSpec]`; `degeneracy_corpus(tool_diameter, tea_cap_deg) -> list[PocketSpec]` with named instances `tangent_island`, `collinear_run`, `exact_half_turn`, `pinch_exactly_tool`, `duplicate_vertex`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/benchmarks/test_topology.py
from __future__ import annotations

from benchmarks.families.topology import island_grid, island_sweep


def test_island_grid_produces_the_requested_hole_count() -> None:
    spec = island_grid(rows=2, cols=3, tool_diameter=1.0, tea_cap_deg=120.0)
    assert len(spec.holes) == 6


def test_island_sweep_is_ordered_by_count() -> None:
    specs = island_sweep(counts=((1, 1), (2, 2), (3, 3)), tool_diameter=1.0, tea_cap_deg=120.0)
    assert [len(s.holes) for s in specs] == [1, 4, 9]


def test_islands_lie_strictly_inside_the_outer_boundary() -> None:
    spec = island_grid(rows=2, cols=2, tool_diameter=1.0, tea_cap_deg=120.0)
    xs = [p[0] for p in spec.polygon.points]
    ys = [p[1] for p in spec.polygon.points]
    for hole in spec.holes:
        for p in hole.points:
            assert min(xs) < p[0] < max(xs)
            assert min(ys) < p[1] < max(ys)
```

```python
# tests/benchmarks/test_degeneracy.py
from __future__ import annotations

from benchmarks.degeneracy import degeneracy_corpus


def test_corpus_names_are_unique_and_stable() -> None:
    names = [s.name for s in degeneracy_corpus(tool_diameter=1.0, tea_cap_deg=120.0)]
    assert len(names) == len(set(names))
    assert {"tangent_island", "collinear_run", "exact_half_turn", "pinch_exactly_tool", "duplicate_vertex"} <= set(names)


def test_every_degenerate_instance_is_still_a_valid_spec() -> None:
    for spec in degeneracy_corpus(tool_diameter=1.0, tea_cap_deg=120.0):
        assert spec.family == "degeneracy"
        assert len(spec.polygon.points) >= 3
        assert spec.tool_radius > 0.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/benchmarks/test_topology.py tests/benchmarks/test_degeneracy.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the topology family**

```python
# benchmarks/families/topology.py
"""Island-count sweep: boolean topology at fixed outer boundary."""

from __future__ import annotations

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# Outer pocket size. Fixed so the sweep varies hole count only.
OUTER_WIDTH = 40.0
OUTER_HEIGHT = 30.0

# Each island is a square of this side length.
ISLAND_SIDE = 3.0


def island_grid(rows: int, cols: int, tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """A rectangular pocket containing a ``rows x cols`` grid of square islands.

    Args:
        rows: Island rows.
        cols: Island columns.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    outer = Polygon([[0.0, 0.0, 0.0], [OUTER_WIDTH, 0.0, 0.0], [OUTER_WIDTH, OUTER_HEIGHT, 0.0], [0.0, OUTER_HEIGHT, 0.0]])
    holes: list[Polygon] = []
    for r in range(rows):
        for c in range(cols):
            cx = OUTER_WIDTH * (c + 1) / (cols + 1)
            cy = OUTER_HEIGHT * (r + 1) / (rows + 1)
            h = 0.5 * ISLAND_SIDE
            # Islands are CW so they read as holes to the exact stock model.
            holes.append(Polygon([[cx - h, cy - h, 0.0], [cx - h, cy + h, 0.0], [cx + h, cy + h, 0.0], [cx + h, cy - h, 0.0]]))
    return PocketSpec.build(
        name=f"islands_{rows}x{cols}",
        family="topology",
        polygon=outer,
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        holes=tuple(holes),
        params={"rows": float(rows), "cols": float(cols), "islands": float(rows * cols)},
    )


def island_sweep(counts: tuple[tuple[int, int], ...], tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Sweep the island count.

    Args:
        counts: ``(rows, cols)`` pairs to generate.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        One spec per pair, in input order.
    """
    return [island_grid(rows=r, cols=c, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) for r, c in counts]
```

- [ ] **Step 4: Write the degeneracy corpus**

```python
# benchmarks/degeneracy.py
"""Hand-built degenerate instances — the expected sources of `unresolved`.

Exactness makes degeneracies RELIABLE, not absent. Each instance here targets one
documented ordinary-not-exceptional case from the engagement kernel's inventory,
so an `unresolved` regression is attributable to a named cause instead of being
noticed as a statistic.
"""

from __future__ import annotations

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# Outer pocket used by most degenerate instances.
BOX_W = 20.0
BOX_H = 12.0


def degeneracy_corpus(tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Every degenerate instance, each targeting one named kernel case.

    Args:
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The corpus, in a stable order.
    """
    box = Polygon([[0.0, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    def spec(name: str, polygon: Polygon, holes: tuple[Polygon, ...] = ()) -> PocketSpec:
        return PocketSpec.build(name=name, family="degeneracy", polygon=polygon, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg, holes=holes, params={})

    r = 0.5 * tool_diameter

    # Island whose left edge sits exactly one tool diameter from the wall, so the
    # cutter rim is exactly tangent to both: a zero-measure contact on each side.
    tangent_x = 2.0 * tool_diameter
    tangent_island = Polygon([[tangent_x, 4.0, 0.0], [tangent_x, 8.0, 0.0], [tangent_x + 4.0, 8.0, 0.0], [tangent_x + 4.0, 4.0, 0.0]])

    # Three exactly collinear boundary vertices: CGAL::orientation returns COLLINEAR.
    collinear = Polygon([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [20.0, 0.0, 0.0], [20.0, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    # A slot exactly 2r wide: the engaged run is an exact half turn, the pi case
    # the kernel decides by orientation rather than by chord comparison.
    half_turn = Polygon([[0.0, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, 2.0 * r, 0.0], [0.0, 2.0 * r, 0.0]])

    # Neck exactly equal to the tool diameter: traversable with zero clearance.
    pinch = Polygon(
        [
            [0.0, 0.0, 0.0],
            [9.0, 0.0, 0.0],
            [10.0, 0.5 * BOX_H - r, 0.0],
            [11.0, 0.0, 0.0],
            [BOX_W, 0.0, 0.0],
            [BOX_W, BOX_H, 0.0],
            [11.0, BOX_H, 0.0],
            [10.0, 0.5 * BOX_H + r, 0.0],
            [9.0, BOX_H, 0.0],
            [0.0, BOX_H, 0.0],
        ]
    )

    # A bit-identical repeated vertex: a zero-length edge the traits must absorb.
    duplicate = Polygon([[0.0, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    return [
        spec("tangent_island", box, (tangent_island,)),
        spec("collinear_run", collinear),
        spec("exact_half_turn", half_turn),
        spec("pinch_exactly_tool", pinch),
        spec("duplicate_vertex", duplicate),
    ]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_topology.py tests/benchmarks/test_degeneracy.py -v -n auto`
Expected: 5 passed

If `duplicate_vertex` or `exact_half_turn` raises out of `PocketSpec.build`, that is a real finding: record it and adjust only the corpus construction (e.g. widen `half_turn` by one tool radius), never `PocketSpec`'s invariants.

- [ ] **Step 6: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): island topology sweep and named degeneracy corpus"
```

---

### Task 10: Report emission

**Files:**
- Create: `benchmarks/report.py`
- Test: `tests/benchmarks/test_report.py`

**Interfaces:**
- Consumes: `MeasurementRecord` (Task 5).
- Produces: `render_markdown(records) -> str`; `write_report(records, out_dir) -> tuple[Path, Path]`.

The generator itself must bake in lead-with-the-conclusion ordering so it survives re-runs: the conclusion paragraph is emitted **first**, before any table. Never emit the literal words "BLUF", "TL;DR", or "Summary" as a heading — just make the first sentence be the answer.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_report.py
from __future__ import annotations

import json

from benchmarks.measurement import MeasurementRecord
from benchmarks.report import render_markdown, write_report


def _record(name: str, certify: float, violations: int, error: str | None = None) -> MeasurementRecord:
    base = MeasurementRecord.failed(name, "necks", {"pinch": 2.0}, 1.0, 120.0, error or "x")
    if error is not None:
        return base
    return MeasurementRecord(**{**base.to_dict(), "error": None, "certify_seconds": certify, "generate_seconds": 0.01, "cap_violations": violations, "operations": 10, "stations": 100})


def test_conclusion_precedes_every_table() -> None:
    text = render_markdown([_record("a", 1.0, 0), _record("b", 9.0, 3)])
    assert text.index("slowest") < text.index("| instance |")


def test_no_banned_summary_labels_appear() -> None:
    text = render_markdown([_record("a", 1.0, 0)])
    for banned in ("BLUF", "TL;DR", "Bottom Line", "## Summary"):
        assert banned not in text


def test_failures_are_reported_not_silently_dropped() -> None:
    text = render_markdown([_record("ok", 1.0, 0), _record("bad", 0.0, 0, error="DegeneratePocketError: narrow")])
    assert "bad" in text
    assert "DegeneratePocketError" in text


def test_write_report_emits_both_artifacts(tmp_path) -> None:
    md_path, json_path = write_report([_record("a", 1.0, 0)], tmp_path)
    assert md_path.exists() and json_path.exists()
    payload = json.loads(json_path.read_text())
    assert payload[0]["name"] == "a"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_report.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the report module**

```python
# benchmarks/report.py
"""Emit the corpus result as markdown and JSON.

The conclusion is emitted BEFORE any table, by construction, so the ordering
survives every re-run rather than depending on whoever edits the output.
"""

from __future__ import annotations

import json
from pathlib import Path

from benchmarks.measurement import MeasurementRecord

MARKDOWN_NAME = "benchmark_report.md"
JSON_NAME = "benchmark_report.json"


def render_markdown(records: list[MeasurementRecord]) -> str:
    """Render the corpus result, conclusion first.

    Args:
        records: Every measurement, any order.

    Returns:
        The markdown document.
    """
    ok = [r for r in records if r.error is None]
    bad = [r for r in records if r.error is not None]
    lines: list[str] = ["# Benchmark corpus result", ""]

    if not ok:
        lines += [f"Every one of the {len(records)} instances failed to measure; no timing conclusion is available.", ""]
    else:
        slowest = max(ok, key=lambda r: r.certify_seconds)
        total_certify = sum(r.certify_seconds for r in ok)
        total_generate = sum(r.generate_seconds for r in ok)
        violating = [r for r in ok if r.cap_violations > 0]
        unresolved = [r for r in ok if r.unresolved > 0]
        ratio = (total_certify / total_generate) if total_generate > 0.0 else float("inf")
        lines += [
            f"The slowest instance is **{slowest.name}** at **{slowest.certify_seconds:.2f} s** to certify "
            f"({slowest.stations} stations); across {len(ok)} measured instances certification costs "
            f"**{ratio:.0f}x** generation ({total_certify:.1f} s against {total_generate:.2f} s). "
            f"**{len(violating)}** instances exceed their cap and **{len(unresolved)}** contain at least one "
            f"undecided motion. {len(bad)} instances failed to measure.",
            "",
        ]

    lines += ["## Per-instance measurements", "", "| instance | family | generate (s) | certify (s) | stations | max TEA (deg) | violations | unresolved | arr. vertices | max digits |", "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in sorted(ok, key=lambda x: -x.certify_seconds):
        lines.append(
            f"| {r.name} | {r.family} | {r.generate_seconds:.3f} | {r.certify_seconds:.2f} | {r.stations} | "
            f"{r.max_tea_deg:.1f} | {r.cap_violations} | {r.unresolved} | {r.arrangement_vertices_final} | {r.max_coordinate_digits} |"
        )
    lines.append("")

    if bad:
        lines += ["## Instances that failed to measure", "", "| instance | family | error |", "| --- | --- | --- |"]
        for r in bad:
            lines.append(f"| {r.name} | {r.family} | {r.error} |")
        lines.append("")

    return "\n".join(lines)


def write_report(records: list[MeasurementRecord], out_dir: Path) -> tuple[Path, Path]:
    """Write the markdown and JSON artifacts.

    Args:
        records: Every measurement.
        out_dir: Destination directory; created when missing.

    Returns:
        ``(markdown_path, json_path)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / MARKDOWN_NAME
    json_path = out_dir / JSON_NAME
    md_path.write_text(render_markdown(records))
    json_path.write_text(json.dumps([r.to_dict() for r in records], indent=2))
    return md_path, json_path
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_report.py -v -n auto`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): conclusion-first markdown and json report emission"
```

---

### Task 11: MATHSM constant-spacing baseline

Held's headline result is **path length**, not speed: his paths are substantially shorter than MATHSM at equal engagement cap. Reproducing that comparison needs the MATHSM baseline — constant machining-circle spacing plus a brute-force search over spacing to hit a target cap.

**Files:**
- Create: `benchmarks/mathsm.py`
- Test: `tests/benchmarks/test_mathsm.py`

**Interfaces:**
- Consumes: `PocketSpec` (Task 1), `compas_cgal.toolpath.trochoidal_mat_toolpath_circular`, `compas_cgal.engagement.audit_toolpath_engagement`.
- Produces: `path_length(result) -> float`; `MathsmPoint` (frozen dataclass: `spacing`, `max_tea_deg`, `length`); `sweep_spacing(spec, spacings) -> list[MathsmPoint]`; `shortest_within_cap(points, cap_deg) -> MathsmPoint | None`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_mathsm.py
from __future__ import annotations

import numpy as np
import pytest

from benchmarks.families.analytic import rectangle
from benchmarks.mathsm import MathsmPoint, path_length, shortest_within_cap, sweep_spacing


class _Result:
    def __init__(self, polyline: np.ndarray) -> None:
        self.polyline = polyline


def test_path_length_sums_polyline_segments() -> None:
    poly = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [3.0, 4.0, 10.0]])
    assert path_length(_Result(poly)) == pytest.approx(15.0)


def test_shortest_within_cap_picks_the_shortest_compliant_point() -> None:
    points = [
        MathsmPoint(spacing=0.1, max_tea_deg=60.0, length=500.0),
        MathsmPoint(spacing=0.3, max_tea_deg=110.0, length=300.0),
        MathsmPoint(spacing=0.6, max_tea_deg=150.0, length=200.0),
    ]
    chosen = shortest_within_cap(points, cap_deg=120.0)
    assert chosen is not None
    assert chosen.spacing == pytest.approx(0.3)


def test_shortest_within_cap_returns_none_when_nothing_complies() -> None:
    points = [MathsmPoint(spacing=0.6, max_tea_deg=150.0, length=200.0)]
    assert shortest_within_cap(points, cap_deg=120.0) is None


def test_sweep_spacing_measures_every_requested_spacing() -> None:
    spec = rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)
    points = sweep_spacing(spec, spacings=(0.4, 0.8))
    assert [p.spacing for p in points] == [0.4, 0.8]
    assert all(p.length > 0.0 for p in points)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_mathsm.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the baseline**

```python
# benchmarks/mathsm.py
"""The MATHSM baseline: constant machining-circle spacing, cap found by search.

Held's comparison is against constant-spacing MATHSM, whose relationship between
spacing and resulting maximum engagement angle has no closed form -- so he varies
the spacing in tiny increments and records the maximum engagement each produces.
This module reproduces that protocol, using `stepover` as the constant-spacing
control and the exact audit as the engagement measurement.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

from benchmarks.spec import PocketSpec

CLEARANCE_Z = 2.0


@dataclass(frozen=True)
class MathsmPoint:
    """One constant-spacing trial.

    Attributes:
        spacing: The stepover used, in the pocket's length units.
        max_tea_deg: Worst engagement angle the audit measured, in degrees.
        length: Total path length.
    """

    spacing: float
    max_tea_deg: float
    length: float


def path_length(result: object) -> float:
    """Total length of a toolpath's polyline, including retracts.

    Args:
        result: A `ToolpathResult` carrying a ``polyline`` of shape (n, 3).

    Returns:
        The summed segment length.
    """
    poly = np.asarray(getattr(result, "polyline"), dtype=float)
    if poly.shape[0] < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(poly, axis=0), axis=1).sum())


def sweep_spacing(spec: PocketSpec, spacings: tuple[float, ...]) -> list[MathsmPoint]:
    """Measure engagement and length at each constant spacing.

    Args:
        spec: The pocket and tool.
        spacings: Stepover values to trial, in input order.

    Returns:
        One point per spacing.
    """
    holes = list(spec.holes)
    points: list[MathsmPoint] = []
    for spacing in spacings:
        result = trochoidal_mat_toolpath_circular(spec.polygon, tool_diameter=spec.tool_diameter, stepover=spacing, holes=holes, clearance_z=CLEARANCE_Z)
        report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, spec.tea_cap_rad, holes=holes)
        points.append(MathsmPoint(spacing=spacing, max_tea_deg=math.degrees(report.max_tea), length=path_length(result)))
    return points


def shortest_within_cap(points: list[MathsmPoint], cap_deg: float) -> MathsmPoint | None:
    """The shortest trial whose measured engagement respects the cap.

    Held records the shorter path whenever several spacings give roughly the same
    maximum engagement; taking the minimum length over all compliant trials is the
    same rule stated without the "roughly".

    Args:
        points: Trials from `sweep_spacing`.
        cap_deg: The engagement cap in degrees.

    Returns:
        The shortest compliant trial, or None when none comply.
    """
    compliant = [p for p in points if p.max_tea_deg <= cap_deg]
    return min(compliant, key=lambda p: p.length) if compliant else None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_mathsm.py -v -n auto`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): MATHSM constant-spacing baseline with brute-force cap search"
```

---

### Task 12: Figure 6 reproduction

The single highest-value artifact in this plan: path length against maximum engagement angle, for one pocket, with the MATHSM curve alongside. It is the only directly reproducible published result in the paper, and overlaying a *certified* curve on it is the entire technical claim in one image.

**Files:**
- Create: `benchmarks/figure6.py`
- Test: `tests/benchmarks/test_figure6.py`

**Interfaces:**
- Consumes: `PocketSpec` (Task 1), `sweep_spacing`/`shortest_within_cap`/`path_length` (Task 11).
- Produces: `Figure6Point` (frozen dataclass: `cap_deg`, `ours_length`, `mathsm_length`, `ours_certified`); `figure6_points(spec, caps, spacings) -> list[Figure6Point]`; `render_figure6_markdown(points) -> str`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_figure6.py
from __future__ import annotations

import pytest

from benchmarks.families.analytic import rectangle
from benchmarks.figure6 import Figure6Point, figure6_points, render_figure6_markdown


def test_points_cover_every_requested_cap() -> None:
    spec = rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)
    points = figure6_points(spec, caps=(60.0, 120.0), spacings=(0.4, 0.8))
    assert [p.cap_deg for p in points] == [60.0, 120.0]


def test_markdown_states_the_comparison_before_the_table() -> None:
    points = [
        Figure6Point(cap_deg=40.0, ours_length=300.0, mathsm_length=400.0, ours_certified=True),
        Figure6Point(cap_deg=80.0, ours_length=200.0, mathsm_length=260.0, ours_certified=True),
    ]
    text = render_figure6_markdown(points)
    assert text.index("shorter") < text.index("| cap (deg) |")


def test_markdown_reports_a_missing_mathsm_baseline_honestly() -> None:
    points = [Figure6Point(cap_deg=40.0, ours_length=300.0, mathsm_length=None, ours_certified=False)]
    text = render_figure6_markdown(points)
    assert "no compliant" in text


def test_length_ratio_is_reported_per_cap() -> None:
    points = [Figure6Point(cap_deg=40.0, ours_length=300.0, mathsm_length=400.0, ours_certified=True)]
    text = render_figure6_markdown(points)
    assert "0.75" in text or "75" in text
    assert pytest.approx(0.75) == points[0].ours_length / points[0].mathsm_length
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_figure6.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the reproduction**

```python
# benchmarks/figure6.py
"""Reproduce Held & Pfeiffer's Figure 6: path length against engagement cap.

Their figure plots three curves for one pocket -- standard, contour-aware, and
MATHSM. This module produces the two we can currently generate: our path at each
cap, and the shortest constant-spacing MATHSM path that respects the same cap.
Our curve additionally carries a certified flag, which their protocol has no way
to report because their engagement evaluation is a discretisation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

from benchmarks.mathsm import CLEARANCE_Z, path_length, shortest_within_cap, sweep_spacing
from benchmarks.spec import PocketSpec

# Caps mirroring the paper's sweep range; 80 degrees is the value used for Fig. 5.
FIGURE6_CAPS: tuple[float, ...] = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0)


@dataclass(frozen=True)
class Figure6Point:
    """One cap value on the reproduction.

    Attributes:
        cap_deg: The engagement cap in degrees.
        ours_length: Length of the path this project generates.
        mathsm_length: Shortest compliant MATHSM length, or None when no trial
            spacing produced a compliant path.
        ours_certified: Whether the audit certified our path at this cap.
    """

    cap_deg: float
    ours_length: float
    mathsm_length: float | None
    ours_certified: bool


def figure6_points(spec: PocketSpec, caps: tuple[float, ...], spacings: tuple[float, ...]) -> list[Figure6Point]:
    """Measure both curves across the cap sweep.

    Args:
        spec: The pocket and tool; its own ``tea_cap_deg`` is overridden per cap.
        caps: Engagement caps in degrees, in plotting order.
        spacings: Constant spacings trialled for the MATHSM baseline.

    Returns:
        One point per cap, in input order.
    """
    holes = list(spec.holes)
    trials = sweep_spacing(spec, spacings)
    points: list[Figure6Point] = []
    for cap_deg in caps:
        cap_rad = math.radians(cap_deg)
        result = trochoidal_mat_toolpath_circular(spec.polygon, tool_diameter=spec.tool_diameter, holes=holes, clearance_z=CLEARANCE_Z)
        report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, cap_rad, holes=holes)
        baseline = shortest_within_cap(trials, cap_deg)
        points.append(
            Figure6Point(
                cap_deg=cap_deg,
                ours_length=path_length(result),
                mathsm_length=baseline.length if baseline is not None else None,
                ours_certified=report.cap_violations == 0,
            )
        )
    return points


def render_figure6_markdown(points: list[Figure6Point]) -> str:
    """Render the reproduction, conclusion first.

    Args:
        points: Measured points, in plotting order.

    Returns:
        The markdown document.
    """
    comparable = [p for p in points if p.mathsm_length is not None and p.mathsm_length > 0.0]
    lines = ["# Figure 6 reproduction — path length against engagement cap", ""]
    if comparable:
        ratios = [p.ours_length / float(p.mathsm_length) for p in comparable]
        mean_ratio = sum(ratios) / len(ratios)
        verdict = "shorter" if mean_ratio < 1.0 else "longer"
        certified = sum(1 for p in points if p.ours_certified)
        lines += [
            f"Across {len(comparable)} comparable caps our path is on average **{mean_ratio:.2f}x** the MATHSM length "
            f"— {verdict} — and **{certified} of {len(points)}** caps are certified, a claim the MATHSM protocol "
            f"cannot make because its engagement evaluation is a discretisation.",
            "",
        ]
    else:
        lines += ["No cap produced a compliant MATHSM baseline from the trialled spacings, so no length comparison is available; widen the spacing sweep.", ""]

    lines += ["| cap (deg) | ours (length) | MATHSM (length) | ratio | certified |", "| ---: | ---: | ---: | ---: | :--- |"]
    for p in points:
        if p.mathsm_length is None:
            lines.append(f"| {p.cap_deg:.0f} | {p.ours_length:.1f} | no compliant spacing | — | {'yes' if p.ours_certified else 'no'} |")
        else:
            lines.append(f"| {p.cap_deg:.0f} | {p.ours_length:.1f} | {p.mathsm_length:.1f} | {p.ours_length / p.mathsm_length:.2f} | {'yes' if p.ours_certified else 'no'} |")
    lines.append("")
    return "\n".join(lines)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_figure6.py -v -n auto`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): Held Fig. 6 reproduction with certified-path overlay"
```

---

### Task 13: External corpus loader

The `unresolved` rate is the number that decides the business, and it can only be measured on geometry this project did not author. The dataset itself is large and separately licensed, so it is **never vendored**: the loader reads closed 2D profiles from a directory the operator supplies.

**Files:**
- Create: `benchmarks/external.py`
- Test: `tests/benchmarks/test_external.py`

**Interfaces:**
- Consumes: `PocketSpec.build` (Task 1).
- Produces: `ExternalCorpusError`; `load_profiles(directory, tool_diameter, tea_cap_deg, min_vertices=4) -> list[PocketSpec]`; `profile_from_points(name, points, tool_diameter, tea_cap_deg) -> PocketSpec`.

Expected on-disk format: one JSON file per profile, `{"name": str, "points": [[x, y], ...]}` — a closed ring, last point not repeated. Converting a licensed dataset (e.g. Autodesk's Fusion 360 Gallery sketch profiles) into this shape is the operator's step and is deliberately outside this repository.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_external.py
from __future__ import annotations

import json

import pytest

from benchmarks.external import ExternalCorpusError, load_profiles, profile_from_points


def _write(tmp_path, name: str, points: list[list[float]]) -> None:
    (tmp_path / f"{name}.json").write_text(json.dumps({"name": name, "points": points}))


def test_loads_every_valid_profile(tmp_path) -> None:
    _write(tmp_path, "a", [[0, 0], [10, 0], [10, 10], [0, 10]])
    _write(tmp_path, "b", [[0, 0], [20, 0], [20, 8], [0, 8]])
    specs = load_profiles(tmp_path, tool_diameter=1.0, tea_cap_deg=120.0)
    assert sorted(s.name for s in specs) == ["a", "b"]
    assert all(s.family == "external" for s in specs)


def test_skips_profiles_below_the_vertex_floor(tmp_path) -> None:
    _write(tmp_path, "tri", [[0, 0], [10, 0], [5, 9]])
    _write(tmp_path, "quad", [[0, 0], [10, 0], [10, 10], [0, 10]])
    specs = load_profiles(tmp_path, tool_diameter=1.0, tea_cap_deg=120.0, min_vertices=4)
    assert [s.name for s in specs] == ["quad"]


def test_missing_directory_fails_loudly(tmp_path) -> None:
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path / "nope", tool_diameter=1.0, tea_cap_deg=120.0)


def test_malformed_profile_fails_loudly(tmp_path) -> None:
    (tmp_path / "bad.json").write_text(json.dumps({"name": "bad"}))
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=1.0, tea_cap_deg=120.0)


def test_profile_from_points_builds_a_spec() -> None:
    spec = profile_from_points("x", [[0, 0], [10, 0], [10, 10], [0, 10]], tool_diameter=1.0, tea_cap_deg=120.0)
    assert spec.name == "x"
    assert len(spec.polygon.points) == 4
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_external.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the loader**

```python
# benchmarks/external.py
"""Load third-party 2D profiles as reference problems.

The `unresolved` rate on geometry this project did not author is the number that
decides whether the certifier is a product. Datasets are separately licensed and
large, so nothing is vendored: point an operator-prepared directory at this
loader. Each file is `{"name": str, "points": [[x, y], ...]}`, a closed ring with
the last point NOT repeated.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

from compas.geometry import Polygon

from benchmarks.errors import BenchmarkError, DegeneratePocketError, PocketNotSimpleError
from benchmarks.spec import PocketSpec


class ExternalCorpusError(BenchmarkError):
    """An external corpus directory is missing, unreadable, or malformed."""


def profile_from_points(name: str, points: Sequence[Sequence[float]], tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """Build a spec from a raw closed ring.

    Args:
        name: Instance name.
        points: Ring vertices as ``[x, y]`` pairs, last not repeated.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.
    """
    ring = Polygon([[float(p[0]), float(p[1]), 0.0] for p in points])
    return PocketSpec.build(name=name, family="external", polygon=ring, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg, params={"vertices": float(len(points))})


def load_profiles(directory: Path, tool_diameter: float, tea_cap_deg: float, min_vertices: int = 4) -> list[PocketSpec]:
    """Load every profile in *directory* that the corpus can machine.

    Profiles rejected by `PocketSpec.build` for being degenerate, non-simple, or
    below the vertex floor are skipped -- they are dataset noise, not certifier
    failures. Malformed JSON is NOT skipped: it is an operator error and raises.

    Args:
        directory: Directory of profile JSON files.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.
        min_vertices: Profiles with fewer vertices are skipped.

    Returns:
        One spec per usable profile, sorted by name.

    Raises:
        ExternalCorpusError: The directory is missing, or a file is malformed.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ExternalCorpusError(f"External corpus directory does not exist: {directory}")
    specs: list[PocketSpec] = []
    for path in sorted(directory.glob("*.json")):
        try:
            payload = json.loads(path.read_text())
            name = str(payload["name"])
            points = payload["points"]
        except (json.JSONDecodeError, KeyError, OSError) as exc:
            raise ExternalCorpusError(f"Malformed profile {path}: {type(exc).__name__}: {exc}") from exc
        if len(points) < min_vertices:
            continue
        try:
            specs.append(profile_from_points(name, points, tool_diameter, tea_cap_deg))
        except (DegeneratePocketError, PocketNotSimpleError):
            continue  # dataset noise, not a certifier finding
    return specs
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/benchmarks/test_external.py -v -n auto`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks
git commit -m "feat(bench): external profile loader for third-party corpora"
```

---

### Task 14: CLI, strict-typing gate, CI wiring, and the docs page

**Files:**
- Create: `benchmarks/cli.py`, `docs/benchmarks.md`
- Modify: `pyproject.toml` (mypy config), `mkdocs.yml` (nav), `.github/workflows/pr-checks.yml`
- Test: `tests/benchmarks/test_cli.py`

**Interfaces:**
- Consumes: every family, `run_corpus` (Task 5), `write_report` (Task 10).
- Produces: `build_corpus(name, tool_diameter, tea_cap_deg) -> list[PocketSpec]`; `main(argv) -> int`; `CORPUS_NAMES: tuple[str, ...]`.

- [ ] **Step 1: Write the failing test**

```python
# tests/benchmarks/test_cli.py
from __future__ import annotations

import pytest

from benchmarks.cli import CORPUS_NAMES, build_corpus, main


def test_every_named_corpus_builds_at_least_one_spec() -> None:
    for name in CORPUS_NAMES:
        if name == "external":
            continue  # needs an operator-supplied directory
        specs = build_corpus(name, tool_diameter=1.0, tea_cap_deg=120.0)
        assert len(specs) > 0, name


def test_unknown_corpus_name_is_rejected() -> None:
    with pytest.raises(SystemExit):
        main(["--corpus", "does-not-exist"])


def test_smoke_run_writes_a_report(tmp_path) -> None:
    exit_code = main(["--corpus", "smoke", "--out", str(tmp_path), "--no-digits"])
    assert exit_code == 0
    assert (tmp_path / "benchmark_report.md").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/benchmarks/test_cli.py -v`
Expected: FAIL — `ModuleNotFoundError`

- [ ] **Step 3: Write the CLI**

```python
# benchmarks/cli.py
"""Run a named corpus and emit its report: `python -m benchmarks.cli --corpus necks`."""

from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.degeneracy import degeneracy_corpus
from benchmarks.external import load_profiles
from benchmarks.families.analytic import arc_channel, disk, rectangle, stadium
from benchmarks.families.complexity import arc_fraction_sweep, ngon_sweep
from benchmarks.families.necks import PINCH_SWEEP_DEFAULT, pinch_sweep
from benchmarks.families.precision import digit_sweep, scale_sweep
from benchmarks.families.topology import island_sweep
from benchmarks.report import write_report
from benchmarks.runner import run_corpus
from benchmarks.spec import PocketSpec

CORPUS_NAMES: tuple[str, ...] = ("smoke", "analytic", "complexity", "necks", "precision", "topology", "degeneracy", "external", "all")

DEFAULT_OUT = Path("docs/benchmarks")


def build_corpus(name: str, tool_diameter: float, tea_cap_deg: float, external_dir: Path | None = None) -> list[PocketSpec]:
    """Assemble the named corpus.

    Args:
        name: One of `CORPUS_NAMES`.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.
        external_dir: Directory of third-party profiles, for the external corpus.

    Returns:
        The instances to measure.

    Raises:
        ValueError: The corpus name is unknown, or `external` was requested
            without a directory.
    """
    if name == "smoke":
        return [rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=tea_cap_deg)]
    if name == "analytic":
        return [
            disk(radius=8.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            rectangle(width=20.0, height=10.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            stadium(straight_length=20.0, half_width=3.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            arc_channel(guide_radius=10.0, half_width=2.0, sweep_deg=90.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
        ]
    if name == "complexity":
        return ngon_sweep(ks=(3, 6, 12, 32, 128), area=100.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) + arc_fraction_sweep(
            n=16, ratios=(0.0, 0.25, 0.5, 1.0), radius=10.0, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg
        )
    if name == "necks":
        return pinch_sweep(pinches=PINCH_SWEEP_DEFAULT, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "precision":
        return scale_sweep(k=12, scales=(1.0, 10.0, 100.0, 1000.0), decimals=6, seed=1, tea_cap_deg=tea_cap_deg) + digit_sweep(
            k=12, radius=10.0, decimal_counts=(1, 3, 6, 9, 12), seed=1, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg
        )
    if name == "topology":
        return island_sweep(counts=((1, 1), (2, 2), (3, 3), (4, 4)), tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "degeneracy":
        return degeneracy_corpus(tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "external":
        if external_dir is None:
            raise ValueError("The external corpus requires --external-dir.")
        return load_profiles(external_dir, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "all":
        out: list[PocketSpec] = []
        for sub in ("analytic", "complexity", "necks", "precision", "topology", "degeneracy"):
            out.extend(build_corpus(sub, tool_diameter, tea_cap_deg))
        return out
    raise ValueError(f"Unknown corpus {name!r}; expected one of {CORPUS_NAMES}.")


def main(argv: list[str] | None = None) -> int:
    """Run a corpus and write its report.

    Args:
        argv: Command-line arguments; defaults to `sys.argv[1:]`.

    Returns:
        Process exit code: 0 on success.
    """
    parser = argparse.ArgumentParser(description="Run a pocket-machining benchmark corpus.")
    parser.add_argument("--corpus", choices=CORPUS_NAMES, default="smoke")
    parser.add_argument("--tool-diameter", type=float, default=1.0)
    parser.add_argument("--cap-deg", type=float, default=120.0)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--external-dir", type=Path, default=None)
    parser.add_argument("--no-digits", action="store_true", help="Skip the untimed exact-coordinate diagnostic pass.")
    args = parser.parse_args(argv)

    specs = build_corpus(args.corpus, args.tool_diameter, args.cap_deg, args.external_dir)
    records = run_corpus(specs, collect_digits=not args.no_digits)
    md_path, json_path = write_report(records, args.out)
    print(f"wrote {md_path} and {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Add the strict-typing gate**

Append to `pyproject.toml`:

```toml
[tool.mypy]
python_version = "3.9"
strict = true
files = ["benchmarks"]
warn_unused_configs = true
```

Run: `mypy` — fix every reported error in `benchmarks/` before continuing. Do not add `# type: ignore` without a trailing comment naming the reason.

- [ ] **Step 5: Wire CI**

Add to `.github/workflows/pr-checks.yml`, as steps in the existing job (after the build step that installs the package):

```yaml
      - name: Benchmark corpus unit tests
        run: pytest tests/benchmarks -q -n auto

      - name: Strict typing gate (benchmarks)
        run: mypy

      - name: Benchmark smoke corpus
        run: python -m benchmarks.cli --corpus smoke --out "${RUNNER_TEMP}/bench" --no-digits
```

- [ ] **Step 6: Write the docs page**

Create `docs/benchmarks.md`, leading with the conclusion and never using the words "BLUF", "TL;DR", or "Summary" as a heading:

```markdown
# Benchmark Corpus

Held & Pfeiffer report 3–100 ms per pocket, but for **one** pocket, with no
hardware stated and no instance table — a plausibility figure, not a
reproducible target. Their generator is fast because it never builds a global
stock arrangement: it maintains the union-of-disks boundary as a sorted arc
sequence updated in amortised-linear work and gets engagement from a closed
form. Their engagement *distributions* come from an admitted discretisation.
This project's **generation** is already competitive; the whole cost gap sits in
a certification step Held does not perform. This corpus measures that gap along
the axes that drive it, and measures the number that decides the product: how
often the certifier answers `unresolved` on geometry it did not author.

## Running it

```bash
python -m benchmarks.cli --corpus all --out docs/benchmarks
python -m benchmarks.cli --corpus external --external-dir /path/to/profiles
```

## The corpora

| Corpus | Axis isolated | Why it exists |
| --- | --- | --- |
| `analytic` | none — closed-form oracle | disk, rectangle, stadium, arc channel: clearance is derivable, so disagreement is a bug, attributable without a second implementation |
| `complexity` | boundary element count and kind | k-gon sweep at fixed area; arc-fraction sweep — closes the corpus's largest hole, since every legacy pocket was straight-sided while the kernel runs on circle-segment traits |
| `necks` | neck severity | pinch sweep from just-traversable outward; locates the certified/uncertified threshold and supplies the neck-degeneration ablation |
| `precision` | coordinate magnitude and significant decimals | the exact kernel's hidden clock: filter-failure rate and rational bit length |
| `topology` | island count | boolean complexity at fixed outer boundary |
| `degeneracy` | named degenerate cases | tangency, collinearity, exact half turn, pinch exactly equal to the tool, duplicate vertex — the expected sources of `unresolved` |
| `external` | real third-party geometry | the resolution-rate measurement; datasets are separately licensed and never vendored |

## Invariants, not just timings

The `analytic` family carries conservation laws. A stadium pocket has constant
clearance along its medial axis, so engagement measured at congruent
configurations must be identical. `benchmarks/congruence.py` tests this with
rigid motions built from **Pythagorean triples** — cos and sin are then rational,
every rotated coordinate stays exactly representable, and the invariant is
*exactly* checkable rather than approximately. It directly exercises the CCW
sort in `finish_engagement`, whose seam is the horizontal line through the
cutter centre and is therefore sensitive to exactly these motions.

## What the report records

Generation and certification are timed **separately**: they are different
businesses and a combined number hides the fact that matters. Alongside them the
report carries certifier stations (the refinement-depth cost proxy), final
arrangement size, and the longest exact rational after depletion.

!!! warning "Coordinate-digit collection perturbs what it measures"

    Reading an exact rational calls `.exact()`, collapsing the lazy filter and
    inflating every later operation. The runner therefore collects digits on a
    separate, untimed pass, and `--no-digits` disables it entirely.
```

Add to `mkdocs.yml` under the `Design Notes` nav section:

```yaml
      - Benchmark Corpus: benchmarks.md
```

- [ ] **Step 7: Run the full suite**

Run: `pytest tests/benchmarks -v -n auto && mypy && ruff check benchmarks tests/benchmarks`
Expected: all tests pass; mypy reports no issues; ruff clean.

- [ ] **Step 8: Commit**

```bash
ruff format benchmarks tests/benchmarks && ruff check --fix benchmarks tests/benchmarks
git add benchmarks tests/benchmarks pyproject.toml mkdocs.yml docs/benchmarks.md .github/workflows/pr-checks.yml
git commit -m "feat(bench): corpus CLI, strict-typing gate, CI wiring, docs page"
```

---

## After the plan: the two runs that matter

These are **not** tasks — they are the measurements the corpus exists to produce, run once the plan lands.

1. **Resolution rate.** `--corpus external` against an operator-prepared third-party profile directory. Report the `unresolved` fraction, broken down by distance to the nearest neck. This is the existential number for both the generator and the verifier thesis.
2. **Neck-degeneration ablation.** `--corpus necks` with the pinch sweep, recording certifier stations and the certified/uncertified threshold against pinch. This converts the strongest claim in the positioning — that a conservative sampled bound degenerates to vacuous at necks while the exact certifier decides — from a derivation into evidence.

## Self-review

- **Coverage.** Analytic families (Task 2) · congruence invariants (Task 3) · instrumentation for arrangement size *and* bit length (Task 4) · phase-separated timing (Task 5) · k-gon and arc-fraction sweeps (Task 6) · neck sweep and threshold (Task 7) · scale and digit sweeps (Task 8) · islands and degeneracies (Task 9) · conclusion-first reporting (Task 10) · MATHSM baseline (Task 11) · Figure 6 (Task 12) · external corpus (Task 13) · CLI, mypy, CI, docs (Task 14). Every recommendation is claimed by exactly one task.
- **Type consistency.** `PocketSpec.build` keyword names are identical at all call sites; `params` is always `Mapping[str, float]`; `MeasurementRecord` field names in Task 5 match every use in Tasks 7, 10, and 14; `CLEARANCE_Z` is defined once in `runner.py` and re-imported by `mathsm.py` and `figure6.py` rather than redefined.
- **Known gap, stated rather than hidden.** `run_spec`'s `unresolved` count is inferred from `OperationEngagement` (an operation that neither certified nor exceeded the cap), because the current `EngagementReport` has no dedicated three-valued field. On a branch where the certifier reports `unresolved` natively, replace that inference with the real field — the schema already has the column.
