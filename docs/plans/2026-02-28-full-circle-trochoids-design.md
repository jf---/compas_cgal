# Full-Circle Trochoidal Toolpaths — Design & Plan

**Branch:** `jf/toolpath-redesign`
**Problem:** Current `trochoid_chain` generates half-circle arcs with alternating tangent sides. The machine must decelerate at every arc-to-line junction. Full 360° circles keep the CNC controller in smooth continuous motion.

**Reference:** CNC Cookbook arc-to-corner pattern — full circles translate along the medial axis with tangent line transitions between them.

---

## Design Decisions

| # | Decision | Choice |
|---|----------|--------|
| 1 | Circle completeness | Full 360° circles (not half-arcs) |
| 2 | Transitions between circles | Tangent line when circles don't overlap; seamless tangent point when they do (pitch-dependent) |
| 3 | Tangent side selection | Always same side per chain (no alternation), determined by milling direction |
| 4 | `merge_circle_pairs()` | Remove entirely — circles are primary output, not a post-hoc merge |
| 5 | Implementation scope | Rewrite `trochoid_chain()` in place — single codepath for both polyline and circular outputs |
| 6 | Winding parameter | Rename `winding_ccw` → `climb_milling` (bool) for intent-revealing semantics |

---

## CGAL API Audit

Research confirms: CGAL's linear kernel (`Epick`) has **no** built-in tangent computation, arc parametric evaluation, or `Circular_arc_2` (requires incompatible Circular Kernel). Our `TrochoidArc` struct and `external_tangents()` are necessarily hand-rolled.

However, we're underusing several CGAL idioms. This plan incorporates them:

### Use `Circle_2` orientation instead of separate `clockwise` bool

`Circle_2` constructor accepts `CGAL::CLOCKWISE` / `CGAL::COUNTERCLOCKWISE` as third argument. Currently we always construct with default orientation and carry a separate `bool clockwise`. Instead:

```cpp
// Before: redundant orientation
Circle_2 circle(center, r*r);  // always default CCW
bool clockwise;                 // separate field

// After: orientation lives in the Circle_2
Circle_2 circle(center, r*r, CGAL::CLOCKWISE);     // for CW
circle.orientation()  // → CGAL::CLOCKWISE
circle.opposite()     // → returns CCW copy
```

**Impact:** Remove `clockwise` field from `TrochoidArc`. Derive winding from `circle.orientation()`. `reversed()` uses `circle.opposite()`. Cleaner, more idiomatic, one source of truth.

### Use `Circle_2::is_degenerate()` for line detection

```cpp
// Before
bool is_line() const { return circle.squared_radius() == K::FT(0); }

// After
bool is_line() const { return circle.is_degenerate(); }
```

### Use `CGAL::has_smaller_distance_to_point()` for nearest-neighbor ordering

```cpp
// Before: compare doubles
double dist = CGAL::to_double(CGAL::squared_distance(cur_end, paths[i].front().start));
if (dist < best_dist) { best_dist = dist; best_idx = i; }

// After: exact predicate, no to_double conversion
if (best_idx < 0 || CGAL::has_smaller_distance_to_point(cur_end, paths[i].front().start, paths[best_idx].front().start)) {
    best_idx = i;
}
```

### Use `CGAL::compare_squared_distance()` for tolerance checks

```cpp
// Before: convert to double, compare
if (CGAL::to_double(CGAL::squared_distance(a, b)) > tol_sq) ...

// After: exact comparison against threshold
if (CGAL::compare_squared_distance(a, b, K::FT(tol_sq)) == CGAL::LARGER) ...
```

### Things that must stay hand-rolled (CGAL has no alternative)

| Feature | Why |
|---------|-----|
| `external_tangents()` | No tangent-line function in CGAL |
| `sweep()` via `atan2` | `CGAL::angle()` returns enum, not numeric |
| `radius()` via `sqrt` | `Circle_2` only stores `squared_radius` |
| `tessellate_chain()` cos/sin | No parametric circle evaluation |
| `TrochoidArc` struct | `Circular_arc_2` requires different kernel |

---

## Architecture Change

### Before (half-arc zigzag)

```
circle[0]  circle[1]  circle[2]  circle[3]
   ○──────────○──────────○──────────○
        ╲        ╱        ╲        ╱
    arc(180°) arc(180°) arc(180°) arc(180°)
    LEFT      RIGHT     LEFT      RIGHT
```

Each iteration alternates tangent sides. Arcs are ~180° on each intermediate circle. The pattern zigzags. `merge_circle_pairs()` can sometimes detect that two consecutive half-arcs complete 360°, but this is incidental.

### After (full-circle + tangent advance)

```
circle[0]     circle[1]     circle[2]
   ⊙─────────────⊙─────────────⊙
   │  tangent     │  tangent     │
   │  line        │  line        │
   ○ (360°)       ○ (360°)       ○ (360°)
   SAME SIDE      SAME SIDE      SAME SIDE
```

Each circle is a full 360° revolution. Tool exits via tangent, advances to next circle's tangent entry, completes another full revolution. When pitch is small enough that circles overlap, the tangent exit of one circle meets the tangent entry of the next — zero-length transition, seamless handoff.

### Key: the tangent connection geometry

For consecutive circles `c[i]` and `c[i+1]`:
1. Compute external tangent on the `climb_milling` side (always the same side)
2. Full circle on `c[i]`: starts at tangent departure point, sweeps 360° back to same point
3. Tangent line from `c[i]` departure to `c[i+1]` arrival (may be zero-length if overlapping)
4. Full circle on `c[i+1]`: starts at tangent arrival point, sweeps 360°

The departure/arrival points are where the external tangent touches each circle. Because we always pick the same side, these points are geometrically consistent.

---

## What Changes

### C++ (`src/toolpath.cpp`)

**Refactor: `TrochoidArc` struct**

- Remove `clockwise` field. Winding direction lives in `Circle_2::orientation()`.
- `is_line()` → `circle.is_degenerate()`
- `reversed()` → uses `circle.opposite()` instead of `!clockwise`
- Add `is_clockwise()` helper: `circle.orientation() == CGAL::CLOCKWISE`
- Tangent methods use `circle.orientation()` instead of reading `clockwise` field
- `make_arc(center, r, start, end, cw)` constructs `Circle_2(center, r*r, cw ? CGAL::CLOCKWISE : CGAL::COUNTERCLOCKWISE)`
- `make_line(s, e)` constructs degenerate `Circle_2(s)` (degenerate, orientation irrelevant)

**New: `TrochoidArc::make_circle()` factory**

```cpp
static TrochoidArc make_circle(const Point_2& center, double r, const Point_2& on_circle, Orientation ori)
```

Returns a `TrochoidArc` where `start == end == on_circle`, circle carries `ori`.

**Modify: `TrochoidArc::sweep()` (lines 73-95)**

Current: returns `0.0` when `start == end`.
New: returns `2π` when `start == end` and `!is_degenerate()` (full circle case).

**Modify: `trochoid_chain()` (lines 212-265)**

New logic:
- Always pick same tangent side (determined by `climb_milling`)
- Emit full-circle `TrochoidArc` centered on `circles[i]` via `make_circle()`
- Emit tangent line segment to next circle (skip if zero-length / overlapping)
- For the last circle, emit one more full circle

Signature change: `bool winding_ccw` → `bool climb_milling`

**Modify: `completes_circle()`**

Remove entirely — no longer needed (circles are primary, not merged pairs).

**Modify: nearest-neighbor ordering in `pmp_trochoidal_mat_toolpath_circular`**

Replace double comparison with `CGAL::has_smaller_distance_to_point()`.

**Remove: `merge_circle_pairs()` (lines 654-680)**

Dead code.

**Remove: `merge_circles` parameter** from `pmp_trochoidal_mat_toolpath_circular` and its binding.

**Modify: serialization (type detection)**

Current: detects circles via `start == end` after merge. New: same check, but circles come from `trochoid_chain` directly. Type code: `0` = line, `1` = arc (not used for cut ops anymore, only lead-in/out), `2` = circle.

### Python (`src/compas_cgal/toolpath.py`)

- Remove `merge_circles` parameter from `trochoidal_mat_toolpath_circular()`
- No other changes — `_row_to_circle()` already handles type=2

### Tests (`tests/test_toolpath.py`)

**Modify: `test_circular_merge_circles`** → **`test_full_circle_output`**
Verify cut operations contain `Circle` geometry objects.

**New: `test_full_circle_sweep`**
Verify circles have sweep ≈ 2π.

**All other tests pass unchanged** — geometric contract is maintained.

### nanobind binding

Remove `merge_circles` kwarg.

---

## Task Sequence

### Task 1: Verify baseline
Run `pytest tests/test_toolpath.py -n auto -v` — all 10 tests pass.

### Task 2: Refactor `TrochoidArc` for CGAL idioms
- Remove `clockwise` field, use `Circle_2::orientation()` instead
- `is_line()` → `circle.is_degenerate()`
- `reversed()` → `circle.opposite()`
- Add `is_clockwise()` helper
- Update `make_arc()` to pass orientation to `Circle_2` constructor
- Update `sweep()`: return `2π` when `start == end && !is_degenerate()`
- Add `make_circle()` factory
- Remove `completes_circle()` method
- Build, verify all existing tests still pass (behavioral equivalence)

### Task 3: Rewrite `trochoid_chain()` for full circles
- Always pick same tangent side (no `i % 2` alternation)
- Emit full circle per intermediate circle via `make_circle()`
- Emit tangent line (skipped when zero-length)
- Rename parameter `winding_ccw` → `climb_milling`
- Update caller `mat_edge_chains()`

Build, run tests. Polyline tests should pass (tessellated full circles ⊂ polygon). Circular tests may need merge test update.

### Task 4: CGAL-idiomatic nearest-neighbor + remove dead code
- Replace double-based nearest-neighbor with `CGAL::has_smaller_distance_to_point()`
- Replace `CGAL::to_double(CGAL::squared_distance(...)) > tol` patterns with `CGAL::compare_squared_distance()`
- Delete `merge_circle_pairs()`
- Remove `merge_circles` param from C++ function, binding, and Python
- Remove `completes_circle()` if not already done

Build, run all tests.

### Task 5: Update tests for full-circle contract
- `test_circular_merge_circles` → `test_full_circle_output`: verify `Circle` geometry in cut ops
- New `test_full_circle_sweep`: verify sweep ≈ 2π
- Run full test suite

### Task 6: Visual verification
Run `docs/examples/example_toolpath_trochoidal_mat.py` — verify smooth full-circle geometry in viewer.

---

## Risks

| Risk | Mitigation |
|------|-----------|
| `sweep()` returning 2π for `start==end` may affect degenerate arcs | Only `make_circle()` produces `start==end` non-degenerate arcs; `make_arc()` never does |
| Removing `clockwise` field touches many call sites | Mechanical refactor — `is_clockwise()` helper keeps callers clean |
| `has_smaller_distance_to_point()` may not be in older CGAL | We use CGAL 6.0.1, confirmed available |
| Tangent alignment test may regress | Full circles exit tangentially by definition — should improve |
| Removing `merge_circles` is API break | Parameter only exists in this branch, not released |

## Build & Test

```bash
pip install --no-build-isolation -ve .
pytest tests/test_toolpath.py -n auto -v
ruff check --fix src/compas_cgal/toolpath.py tests/test_toolpath.py && ruff format src/compas_cgal/toolpath.py tests/test_toolpath.py
```
