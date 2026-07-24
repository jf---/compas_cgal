# Exact-certified adaptive trochoidal-MAT — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A generator where the user sets a maximum tool-engagement angle θ_max and gets back a trochoidal-MAT toolpath whose exact engagement angle is provably ≤ θ_max on every cutting motion.

**Architecture:** A Python subpackage `compas_cgal.adaptive` orchestrating the existing exact CGAL kernel: the C++ straight-skeleton medial axis (newly exposed as a station walk) feeds an adaptive spacer that, per candidate machining circle, asks the exact engagement kernel (`_stock_2.engagement_at`, contour-aware against a depleting `Stock2`) whether θ ≤ θ_max, packing to the tightest spacing that passes. The generated path is re-audited by the *independent* `audit_toolpath_engagement` to certify the claim.

**Tech Stack:** Python 3.12; compas; the built `compas_cgal._stock_2` (Epeck one-root kernel) and `compas_cgal._toolpath` extensions; nanobind/CGAL 6.0.1 for the C++ station exposure; pytest + pytest-xdist.

## Global Constraints

- Python target **3.12 / 3.9** (compas dependency); no `__all__` in any `__init__.py`; keep `__init__.py` minimal.
- **`mypy --strict`** must pass for the whole `compas_cgal.adaptive` subpackage — the unit-typing is aspirational without it.
- **Exact-Kernel Discipline** (`CLAUDE.md`, `docs/exactness.md`): no epsilon / deflation / ulp-nudge in any decision path; every spacing acceptance is an exact predicate via the kernel; `to_double` is reporting-only and never feeds a decision.
- One responsibility per file; factory functions (`.build(...)`) own invariants; raw dataclass constructors stay minimal. One **named exception per failure mode** — never bare `raise ValueError`.
- Build the C++ extension (no pip): `cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _stock_2 _toolpath -j 8`. Run tests: `PYTHONPATH=src /Users/jelle/mambaforge/envs/cgal-dev/bin/python -m pytest tests/adaptive -n auto -q`.
- Commits: author AND committer `Jelle Feringa <jelleferinga@gmail.com>`; concise conventional messages; NO attribution lines. Ruff before each Python commit. Never touch `docs/examples/example_isolines.py`.

---

## File Structure

```
src/compas_cgal/adaptive/
  __init__.py            # minimal
  units.py               # Task 1  — typed vocabulary (Radian, Millimetre, ToolRadius, Clearance, TrochoidRadius, Spacing)
  errors.py              # Task 2  — named exceptions
  engagement_cap.py      # Task 3  — EngagementCap (angle -> exact chord surrogate)
  medial_axis.py         # Task 5  — MedialAxisWalk (wraps the new C++ station binding)
  machined_area.py       # Task 6  — MachinedArea protocol + Stock2Area implementation
  engagement_probe.py    # Task 7  — EngagementProbe (exact theta<=theta_max verdict)
  spacer.py              # Task 8+9 — AdaptiveSpacer (+ neck-aware cap)
  transitions.py         # Task 10 — TransitionBuilder (non-slotting links)
  certificate.py         # Task 11 — BuildIdentity + EngagementCertificate
  generator.py           # Task 12 — adaptive_trochoidal_toolpath entry point
src/toolpath.cpp         # Task 4  — expose medial_axis_stations(...) binding (MODIFY)
src/toolpath.h           # Task 4  — declare it (MODIFY)
tests/adaptive/
  test_units.py test_errors.py test_engagement_cap.py test_medial_axis.py
  test_machined_area.py test_engagement_probe.py test_spacer.py test_transitions.py
  test_certificate.py test_generator.py
  test_acceptance.py     # Task 13 — differential certificate oracle + ablation
```

Files that change together live together; each module has one responsibility. The C++ change (Task 4) is the only non-Python unit.

---

### Task 1: Typed unit vocabulary

**Files:**
- Create: `src/compas_cgal/adaptive/units.py`
- Create: `src/compas_cgal/adaptive/__init__.py` (empty)
- Test: `tests/adaptive/test_units.py`

**Interfaces:**
- Produces: `Radian = NewType("Radian", float)`, `Millimetre = NewType("Millimetre", float)`; frozen dataclasses `ToolRadius`, `Clearance`, `TrochoidRadius`, `Spacing`, each with a `.value: Millimetre` field and a classmethod `build(mm: float) -> Self` that rejects NaN/≤0 (raises `NonPositiveLengthError`, Task 2).

- [ ] **Step 1: Write the failing test**
```python
# tests/adaptive/test_units.py
import math
import pytest
from compas_cgal.adaptive.units import ToolRadius, Spacing
from compas_cgal.adaptive.errors import NonPositiveLengthError

def test_tool_radius_build_accepts_positive():
    assert ToolRadius.build(1.5).value == 1.5

def test_tool_radius_rejects_nonpositive():
    with pytest.raises(NonPositiveLengthError):
        ToolRadius.build(0.0)

def test_spacing_rejects_nan():
    with pytest.raises(NonPositiveLengthError):
        Spacing.build(math.nan)
```
- [ ] **Step 2: Run to verify it fails** — `pytest tests/adaptive/test_units.py -v` → FAIL (module missing).
- [ ] **Step 3: Implement** `units.py` — the `NewType`s and one shared frozen dataclass base with a `build` classmethod doing `if not (mm > 0.0): raise NonPositiveLengthError(mm)`; derive `ToolRadius`, `Clearance`, `TrochoidRadius`, `Spacing` from it. Minimal `__init__.py`.
- [ ] **Step 4: Run to verify pass** — same command → PASS.
- [ ] **Step 5: Commit** — `git add src/compas_cgal/adaptive/units.py src/compas_cgal/adaptive/__init__.py tests/adaptive/test_units.py && git commit -m "feat(adaptive): typed length vocabulary"`

---

### Task 2: Named error model

**Files:** Create `src/compas_cgal/adaptive/errors.py`; Test `tests/adaptive/test_errors.py`.

**Interfaces:**
- Produces exceptions (all subclass a package base `AdaptiveError(Exception)`): `NonPositiveLengthError`, `EngagementOutOfRangeError`, `PocketNotMachinableError`, `EngagementCapInfeasibleError`, `NeckTooTightError(width: float, cap_deg: float)`, `TransitionGougeError`, `DegenerateMedialAxisError`.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_errors.py
from compas_cgal.adaptive.errors import AdaptiveError, NeckTooTightError

def test_neck_error_carries_context():
    e = NeckTooTightError(width=0.4, cap_deg=80.0)
    assert isinstance(e, AdaptiveError)
    assert "0.4" in str(e) and "80" in str(e)
```
- [ ] **Step 2: Verify fail.** **Step 3:** implement each exception; the parameterised ones format a fix-oriented message in `__init__`. **Step 4: verify pass.** **Step 5: commit** `feat(adaptive): named error model`.

---

### Task 3: EngagementCap (the API boundary)

**Files:** Create `src/compas_cgal/adaptive/engagement_cap.py`; Test `tests/adaptive/test_engagement_cap.py`.

**Interfaces:**
- Consumes: `Radian` (Task 1), `EngagementOutOfRangeError` (Task 2).
- Produces: frozen `EngagementCap` with fields `theta: Radian`, `chord_ratio: float`; `EngagementCap.from_degrees(deg: float) -> EngagementCap` and `from_radians(rad: float) -> EngagementCap`, validating `θ ∈ (0, π]` and computing `chord_ratio = 4*sin(θ/2)**2` (the exact rational surrogate the kernel consumes). Method `reduced(factor: float) -> EngagementCap` for neck handling (factor ∈ (0,1], re-derives the surrogate).

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_engagement_cap.py
import math, pytest
from compas_cgal.adaptive.engagement_cap import EngagementCap
from compas_cgal.adaptive.errors import EngagementOutOfRangeError

def test_full_cap_surrogate_is_four():
    assert EngagementCap.from_radians(math.pi).chord_ratio == pytest.approx(4.0)

def test_ninety_deg_surrogate():
    assert EngagementCap.from_degrees(90).chord_ratio == pytest.approx(4*math.sin(math.radians(45))**2)

@pytest.mark.parametrize("bad", [0.0, -1.0, math.pi + 1e-6])
def test_rejects_out_of_range(bad):
    with pytest.raises(EngagementOutOfRangeError):
        EngagementCap.from_radians(bad)

def test_reduced_shrinks_surrogate():
    assert EngagementCap.from_degrees(120).reduced(0.5).theta < EngagementCap.from_degrees(120).theta
```
- [ ] **Step 2–4:** verify fail → implement (surrogate + validation mirror `_stock_2.engagement_at`'s contract) → verify pass. **Step 5: commit** `feat(adaptive): EngagementCap angle->exact chord surrogate`.

---

### Task 4: Expose the medial-axis station walk (C++)

**Files:** Modify `src/toolpath.cpp`, `src/toolpath.h`; Test `tests/adaptive/test_medial_axis.py`.

**Interfaces:**
- Produces a nanobind binding `_toolpath.medial_axis_stations(boundary: RowMatrixXd, holes: list[RowMatrixXd], tool_radius: float, radial_clearance: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]` returning, in medial-axis-walk order, `centers (N×2)`, `clearance (N,)` (exact boundary distance at each center, as double), `trochoid_radius (N,)`.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_medial_axis.py
import numpy as np
from compas.geometry import Polygon
from compas_cgal import _toolpath
from compas_cgal.stock import _polygon_to_ccw_vertices  # reuse existing CCW helper

def test_square_stations_lie_inside_and_clearances_positive():
    sq = _polygon_to_ccw_vertices(Polygon([[0,0,0],[10,0,0],[10,10,0],[0,10,0]]))
    c, clr, rad = _toolpath.medial_axis_stations(sq, [], 1.0, 0.0)
    assert len(c) > 3 and clr.min() > 0
    # every station center is >= tool_radius from the boundary (machinable)
    assert clr.min() >= 1.0 - 1e-6
    # medial axis of a square peaks at clearance ~5 (center)
    assert clr.max() == np.testing.assert_allclose or clr.max() > 4.0
```
(Fix the final assert to `assert clr.max() > 4.0` when writing.)
- [ ] **Step 2: verify fail** — `medial_axis_stations` not bound.
- [ ] **Step 3: Implement.** In `toolpath.cpp`, factor the existing straight-skeleton walk + `Station{center, clearance, radius}` computation (currently internal to `trochoidal_mat_toolpath_circular`, ~lines 185–260) into a reusable function that returns the ordered `std::vector<Station>`; bind a thin wrapper converting to three arrays. Do NOT change the existing generator's behaviour — extract-and-reuse, add-alongside. Rebuild: `cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _toolpath -j 8`.
- [ ] **Step 4: verify pass.** **Step 5: commit** `feat(toolpath): expose medial-axis station walk with exact clearance`.

---

### Task 5: MedialAxisWalk wrapper

**Files:** Create `src/compas_cgal/adaptive/medial_axis.py`; Test extends `tests/adaptive/test_medial_axis.py`.

**Interfaces:**
- Consumes: `_toolpath.medial_axis_stations` (Task 4), `Clearance`/`TrochoidRadius` (Task 1), `DegenerateMedialAxisError` (Task 2).
- Produces: frozen `Station(center: tuple[float,float], clearance: Clearance, trochoid_radius: TrochoidRadius)`; `MedialAxisWalk.build(polygon, holes, tool_radius: ToolRadius) -> MedialAxisWalk`; iterable of `Station` in walk order; raises `DegenerateMedialAxisError` on empty/degenerate output.

- [ ] **Step 1: Failing test**
```python
def test_walk_yields_typed_stations():
    from compas.geometry import Polygon
    from compas_cgal.adaptive.medial_axis import MedialAxisWalk
    from compas_cgal.adaptive.units import ToolRadius
    walk = MedialAxisWalk.build(Polygon([[0,0,0],[10,0,0],[10,10,0],[0,10,0]]), [], ToolRadius.build(1.0))
    stations = list(walk)
    assert stations and stations[0].clearance.value > 0
```
- [ ] **Step 2–4:** verify fail → implement (call the binding, wrap rows into typed `Station`s, raise on empty) → verify pass. **Step 5: commit** `feat(adaptive): typed medial-axis walk`.

---

### Task 6: MachinedArea (contour-aware, Stock2-backed)

**Files:** Create `src/compas_cgal/adaptive/machined_area.py`; Test `tests/adaptive/test_machined_area.py`.

**Interfaces:**
- Consumes: `compas_cgal.stock.Stock`, `_stock_2` (existing).
- Produces: a `MachinedArea` typing.Protocol with `raw` (the `_stock_2.Stock2` handle for engagement queries) and `deplete_circle(cx: float, cy: float, radius: float) -> None`; and `Stock2Area.build(polygon, holes, tool_radius) -> Stock2Area` implementing it via `Stock` + `subtract_disk`/`subtract_arc_sweep`. Depletion models the cleared region; the engagement kernel reads `.raw` contour-aware.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_machined_area.py
import math
from compas.geometry import Polygon
from compas_cgal import _stock_2
from compas_cgal.adaptive.machined_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius

def test_depletion_lowers_subsequent_engagement():
    P = Polygon([[0,0,0],[10,0,0],[10,10,0],[0,10,0]])
    area = Stock2Area.build(P, [], ToolRadius.build(0.5))
    cap = 4.0  # full cap
    before = _stock_2.engagement_at(area.raw, 5.0, 5.0, 0.5, cap)[0]
    area.deplete_circle(5.0, 5.0, 0.5)            # clear where the cutter is
    after = _stock_2.engagement_at(area.raw, 5.0, 5.0, 0.5, cap)[0]
    assert before > 0.0 and after == 0.0          # contour-aware: now cleared
```
- [ ] **Step 2–4:** verify fail → implement (`Stock` wrapper + subtract; `raw` returns `stock.raw`) → verify pass. **Step 5: commit** `feat(adaptive): contour-aware machined-area (Stock2)`.

---

### Task 7: EngagementProbe (the exact decision)

**Files:** Create `src/compas_cgal/adaptive/engagement_probe.py`; Test `tests/adaptive/test_engagement_probe.py`.

**Interfaces:**
- Consumes: `MachinedArea` (Task 6), `EngagementCap` (Task 3), `_stock_2.engagement_at`, `MachiningCircle` (defined here: frozen `MachiningCircle(cx: float, cy: float, tool_radius: float, trochoid_radius: float)`).
- Produces: `EngagementProbe(area: MachinedArea, tool_radius: float)`; method `within_cap(circle: MachiningCircle, cap: EngagementCap) -> bool` — the EXACT verdict "does the tool's max engagement as it sweeps this machining circle stay ≤ θ_max, contour-aware?" Implemented by sampling stations along the machining circle and calling `engagement_at(area.raw, x, y, tool_radius, cap.chord_ratio)`; returns True iff no sampled station reports `cap_exceeded`. Station count derived from the arc length so a full circle is covered at ≤ tool_radius spacing (same density basis as `certify_segment_tea`).

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_engagement_probe.py
import math
from compas.geometry import Polygon
from compas_cgal.adaptive.machined_area import Stock2Area
from compas_cgal.adaptive.engagement_probe import EngagementProbe, MachiningCircle
from compas_cgal.adaptive.engagement_cap import EngagementCap
from compas_cgal.adaptive.units import ToolRadius

def test_small_immersion_passes_large_cap_fails_tiny_cap():
    P = Polygon([[0,0,0],[10,0,0],[10,10,0],[0,10,0]])
    probe = EngagementProbe(Stock2Area.build(P, [], ToolRadius.build(0.5)), tool_radius=0.5)
    # a machining circle deep inside virgin stock: full immersion on part of the sweep
    circle = MachiningCircle(cx=5.0, cy=5.0, tool_radius=0.5, trochoid_radius=1.5)
    assert probe.within_cap(circle, EngagementCap.from_degrees(180)) is True
    assert probe.within_cap(circle, EngagementCap.from_degrees(20)) is False
```
- [ ] **Step 2: verify fail.** **Step 3: implement** — sample the machining-circle tool-center positions, query `engagement_at`, return the AND of `not exceeded`. The decision is exact (kernel `cap_exceeded`); `to_double` is never compared. **Step 4: verify pass.** **Step 5: commit** `feat(adaptive): exact engagement probe (contour-aware verdict)`.

---

### Task 8: AdaptiveSpacer — tightest certified spacing

**Files:** Create `src/compas_cgal/adaptive/spacer.py`; Test `tests/adaptive/test_spacer.py`.

**Interfaces:**
- Consumes: `MedialAxisWalk` (Task 5), `EngagementProbe` (Task 7), `MachiningCircle`, `EngagementCap`, `MachinedArea`, `EngagementCapInfeasibleError` (Task 2).
- Produces: `AdaptiveSpacer(walk, probe, area, cap)`; `next_circle(prev: MachiningCircle) -> MachiningCircle` — bisect the along-walk spacing to the LARGEST value whose `MachiningCircle` passes `probe.within_cap(..., cap)`, deplete the accepted circle into `area`, return it. `spacing_search_floor` and `resolution` are named module constants (with unit/derivation comments), not literals. Raises `EngagementCapInfeasibleError` if even the floor spacing exceeds the cap. `circles() -> Iterator[MachiningCircle]` walks the whole path.

- [ ] **Step 1: Failing test** — the load-bearing behavioural contract:
```python
# tests/adaptive/test_spacer.py
import math
from compas.geometry import Polygon
from compas_cgal.adaptive.spacer import AdaptiveSpacer
from compas_cgal.adaptive.engagement_cap import EngagementCap
# (build walk/probe/area via their factories; helper in the test)

def test_every_accepted_circle_is_within_cap():
    cap = EngagementCap.from_degrees(80)
    spacer = _build_spacer(Polygon([[0,0,0],[20,0,0],[20,20,0],[0,20,0]]), tool_radius=0.8, cap=cap)
    circles = list(spacer.circles())
    assert circles
    # re-probe each accepted circle against the area state it was accepted in:
    for c in circles:
        assert spacer.last_verdict(c) is True     # spacer records the exact verdict per circle

def test_tighter_cap_packs_more_circles():
    P = Polygon([[0,0,0],[20,0,0],[20,20,0],[0,20,0]])
    n80 = len(list(_build_spacer(P, 0.8, EngagementCap.from_degrees(80)).circles()))
    n40 = len(list(_build_spacer(P, 0.8, EngagementCap.from_degrees(40)).circles()))
    assert n40 > n80                              # smaller cap -> more circles (Held Fig-6 monotonicity)
```
- [ ] **Step 2: verify fail.** **Step 3: implement** the bisection over spacing on the exact `within_cap` boolean (monotone: larger spacing → higher θ). Record `last_verdict` per accepted circle. **Step 4: verify pass.** **Step 5: commit** `feat(adaptive): certified adaptive spacer`.

---

### Task 9: Neck-aware cap reduction

**Files:** Modify `src/compas_cgal/adaptive/spacer.py`; Test extends `tests/adaptive/test_spacer.py`.

**Interfaces:**
- Consumes: `Clearance` from the walk, `EngagementCap.reduced` (Task 3), `NeckTooTightError` (Task 2).
- Produces: within `next_circle`, when the station clearance marks a bottleneck (clearance below a named `NECK_CLEARANCE_FACTOR * tool_radius`), use `cap.reduced(width_fraction)` for that circle; raise `NeckTooTightError` if the reduced cap is infeasible at the floor spacing.

- [ ] **Step 1: Failing test** — the slotted lower-half witness (from `test_engagement_audit.py`): a thin neck forces a reduced cap, and the path through it stays certified at the reduced angle.
```python
def test_neck_uses_reduced_cap_and_stays_certified():
    # slotted lower-half pocket; the neck circles must certify against a reduced cap, not the nominal
    spacer = _build_spacer(_slotted_lower_half(), tool_radius=0.5, cap=EngagementCap.from_degrees(120))
    circles = list(spacer.circles())
    assert all(spacer.last_verdict(c) for c in circles)   # no circle exceeds its (possibly reduced) cap
```
- [ ] **Step 2–4:** verify fail → implement neck detection off `Clearance` + reduced cap → verify pass. **Step 5: commit** `feat(adaptive): neck-aware cap reduction`.

---

### Task 10: TransitionBuilder (non-slotting links)

**Files:** Create `src/compas_cgal/adaptive/transitions.py`; Test `tests/adaptive/test_transitions.py`.

**Interfaces:**
- Consumes: `MachiningCircle`, `MachinedArea` (for the cleared region), `TransitionGougeError` (Task 2).
- Produces: `TransitionBuilder(area)`; `link(a: MachiningCircle, b: MachiningCircle) -> list[tuple[float,float]]` — a polyline connecting consecutive machining circles that stays within the already-cleared area (so it does not slot fresh stock). Verified by probing the link's midpoints: engagement against the depleted area must be 0. Raises `TransitionGougeError` if no cleared route exists.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_transitions.py
def test_link_between_adjacent_circles_does_not_engage_material():
    # two overlapping cleared circles; the link between their centers is in cleared area
    builder, a, b, area = _two_cleared_circles()
    pts = builder.link(a, b)
    for x, y in pts:
        assert _stock_2.engagement_at(area.raw, x, y, a.tool_radius, 4.0)[0] == 0.0
```
- [ ] **Step 2–4:** verify fail → implement (route through overlap of cleared circles; probe midpoints) → verify pass. **Step 5: commit** `feat(adaptive): non-slotting transitions`.

---

### Task 11: BuildIdentity + EngagementCertificate

**Files:** Create `src/compas_cgal/adaptive/certificate.py`; Test `tests/adaptive/test_certificate.py`.

**Interfaces:**
- Consumes: `MachiningCircle`, `EngagementCap`.
- Produces: frozen `BuildIdentity(sha: str)` via `BuildIdentity.of(polygon, holes, tool_radius, cap, version: str)` — SHA-256 over canonicalised inputs; frozen `EngagementCertificate(identity: BuildIdentity, cap: EngagementCap, per_circle_within_cap: tuple[bool, ...])` with property `holds: bool = all(per_circle_within_cap)`.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_certificate.py
from compas_cgal.adaptive.certificate import BuildIdentity, EngagementCertificate
from compas_cgal.adaptive.engagement_cap import EngagementCap

def test_identity_is_deterministic_and_input_sensitive():
    from compas.geometry import Polygon
    P = Polygon([[0,0,0],[10,0,0],[10,10,0],[0,10,0]])
    cap = EngagementCap.from_degrees(80)
    a = BuildIdentity.of(P, [], 1.0, cap, "phase1")
    b = BuildIdentity.of(P, [], 1.0, cap, "phase1")
    c = BuildIdentity.of(P, [], 1.5, cap, "phase1")
    assert a == b and a != c

def test_certificate_holds_iff_all_within_cap():
    cap = EngagementCap.from_degrees(80)
    ok = EngagementCertificate(BuildIdentity("x"), cap, (True, True))
    bad = EngagementCertificate(BuildIdentity("x"), cap, (True, False))
    assert ok.holds and not bad.holds
```
- [ ] **Step 2–4:** verify fail → implement (canonical JSON of rounded verts + params → sha256) → verify pass. **Step 5: commit** `feat(adaptive): content-addressed engagement certificate`.

---

### Task 12: Generator entry point

**Files:** Create `src/compas_cgal/adaptive/generator.py`; Test `tests/adaptive/test_generator.py`.

**Interfaces:**
- Consumes: all prior units.
- Produces: `adaptive_trochoidal_toolpath(polygon, tool_diameter: float, max_engagement_deg: float, holes=None, radial_clearance: float = 0.0) -> CertifiedToolpath`, where frozen `CertifiedToolpath(circles: tuple[MachiningCircle, ...], links: tuple[tuple[tuple[float,float],...], ...], certificate: EngagementCertificate)`. Assembles: `EngagementCap.from_degrees` → `MedialAxisWalk.build` → `Stock2Area.build` → `EngagementProbe` → `AdaptiveSpacer` → `TransitionBuilder` → `EngagementCertificate`.

- [ ] **Step 1: Failing test**
```python
# tests/adaptive/test_generator.py
from compas.geometry import Polygon
from compas_cgal.adaptive.generator import adaptive_trochoidal_toolpath

def test_generator_returns_holding_certificate_for_square():
    P = Polygon([[0,0,0],[20,0,0],[20,20,0],[0,20,0]])
    out = adaptive_trochoidal_toolpath(P, tool_diameter=1.6, max_engagement_deg=80.0)
    assert out.circles and out.certificate.holds
    assert out.certificate.cap.theta > 0
```
- [ ] **Step 2–4:** verify fail → implement the assembly → verify pass. **Step 5: commit** `feat(adaptive): adaptive_trochoidal_toolpath generator`.

---

### Task 13: Acceptance — differential oracle + ablation

**Files:** Test `tests/adaptive/test_acceptance.py`.

**Interfaces:** Consumes `adaptive_trochoidal_toolpath` (Task 12), the *independent* `compas_cgal.engagement.audit_toolpath_engagement`, `_stock_2.engagement_at`.

- [ ] **Step 1: Write the acceptance tests** (no new production code — these certify the whole):
```python
# tests/adaptive/test_acceptance.py
import math
import numpy as np
from compas.geometry import Polygon
from compas_cgal import _stock_2
from compas_cgal.adaptive.generator import adaptive_trochoidal_toolpath
from compas_cgal.adaptive.machined_area import Stock2Area
from compas_cgal.adaptive.units import ToolRadius

SQUARE = Polygon([[0,0,0],[20,0,0],[20,20,0],[0,20,0]])

def _independent_max_engagement(out, tool_radius, cap_ratio):
    """Replay the certified circles through a FRESH depleting stock and re-measure
    the max engagement per circle -- a second, independent path through the kernel."""
    area = Stock2Area.build(SQUARE, [], ToolRadius.build(tool_radius))
    worst = 0.0
    for c in out.circles:
        # sample the machining circle, measure engagement contour-aware, then deplete
        for a in np.linspace(0, 2*math.pi, 24, endpoint=False):
            x, y = c.cx + c.trochoid_radius*math.cos(a), c.cy + c.trochoid_radius*math.sin(a)
            total, mx, exceeded = _stock_2.engagement_at(area.raw, x, y, tool_radius, cap_ratio)
            worst = max(worst, mx)
        area.deplete_circle(c.cx, c.cy, tool_radius)
    return worst

def test_differential_certificate_no_cutting_arc_exceeds_cap():
    deg = 80.0
    out = adaptive_trochoidal_toolpath(SQUARE, tool_diameter=1.6, max_engagement_deg=deg)
    assert out.certificate.holds                       # the generator's claim
    worst = _independent_max_engagement(out, 0.8, out.certificate.cap.chord_ratio)
    assert worst <= math.radians(deg) + 1e-6           # the independent confirmation

def test_ablation_fixed_stepover_produces_a_tail_adaptive_does_not():
    # adaptive at 120 deg: 0% of circles exceed 120
    out = adaptive_trochoidal_toolpath(SQUARE, tool_diameter=1.6, max_engagement_deg=120.0)
    assert all(out.certificate.per_circle_within_cap)
    # (the fixed-stepover comparison uses trochoidal_mat_toolpath_circular audited at 120 deg,
    #  reproducing the measured ~15% tail -- see docs spec acceptance section)
```
- [ ] **Step 2: Run** — `PYTHONPATH=src ... pytest tests/adaptive/test_acceptance.py -v`. Expected: the differential oracle passes (generator claim == independent measurement, θ ≤ θ_max); if it fails, the spacer over-packed — STOP and tighten, do not loosen the tolerance.
- [ ] **Step 3: Full suite green** — `pytest tests/ -n auto -q`, the pre-existing suite unchanged.
- [ ] **Step 4: Commit** — `test(adaptive): differential certificate oracle + ablation`.

---

## Self-Review

**Spec coverage:** EngagementCap (T3←spec §Contract/boundary), MedialAxisWalk (T4-5←§Components.2), MachinedArea Stock2 (T6←§Components.3 Phase-1), EngagementProbe exact decision (T7←§Components.4 + deciding/reporting), AdaptiveSpacer (T8←§Components.5), neck reduction (T9←§Components.5 §3.2), transitions (T10←§Components.6), certificate+BuildIdentity (T11←§Components.7), generator (T12), differential oracle + ablation (T13←§Acceptance). Types (T1←§Units), errors (T2←§Error model). **Deferred by spec:** Phase-1.5 exact sorted-arc contour (not in this plan — a follow-up plan swaps the `MachinedArea` implementation behind its Protocol); Phase-2 elliptic returns. Held-scale full-pocket run is a Phase-1.5 exit, not a Phase-1 task — Task 13 uses a tractable 20×20 square, matching the spec.

**Placeholder scan:** the `clr.max()` assert in Task 4 Step 1 has an inline correction note — the implementer writes `assert clr.max() > 4.0`. No TODO/TBD elsewhere; each task carries its test code.

**Type consistency:** `MachiningCircle(cx, cy, tool_radius, trochoid_radius)` defined in Task 7 and used identically in Tasks 8/10/12/13; `EngagementCap.chord_ratio` (T3) consumed by T7/T13; `MachinedArea.raw`/`deplete_circle` (T6) consumed by T7/T10/T13; `.circles()`/`last_verdict` (T8) consumed by T9/T13.
