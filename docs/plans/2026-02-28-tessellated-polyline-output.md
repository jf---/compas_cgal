# Tessellated Polyline Output for Circular Toolpath

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a C++-tessellated polyline to `trochoidal_mat_toolpath_circular` so traversal order is authoritative, not heuristically reconstructed in Python.

**Architecture:** The C++ function already builds a `std::vector<ToolpathPrimitive>` with correct CW/CCW winding. A new `tessellate_operations` function walks this vector, sampling arcs at `samples_per_radian` resolution and lines at their endpoints, with exact `Point_2` junction dedup. The polyline is returned as a 6th matrix in the tuple. Python wraps it in a `ToolpathResult` dataclass alongside the existing `list[ToolpathOperation]`.

**Tech Stack:** C++20, CGAL Epick, nanobind, Eigen, numpy

---

### Task 1: C++ tessellate_operations function

**Files:**
- Modify: `src/toolpath.cpp:245` (after existing `tessellate_chain`)

**Step 1: Write `tessellate_operations`**

Insert after `tessellate_chain` (line 282), inside the anonymous namespace:

```cpp
compas::RowMatrixXd
tessellate_operations(const std::vector<ToolpathPrimitive>& ops, double samples_per_radian)
{
    if (ops.empty()) return compas::RowMatrixXd(0, 3);

    std::vector<std::array<double, 3>> pts;
    pts.reserve(ops.size() * 32);

    for (std::size_t oi = 0; oi < ops.size(); ++oi) {
        const auto& op = ops[oi];
        const bool first = pts.empty();

        if (op.arc.is_line()) {
            if (first) {
                pts.push_back({CGAL::to_double(op.arc.start.x()),
                               CGAL::to_double(op.arc.start.y()),
                               op.z_start});
            }
            pts.push_back({CGAL::to_double(op.arc.end.x()),
                           CGAL::to_double(op.arc.end.y()),
                           op.z_end});
        } else {
            const double r = op.arc.radius();
            const double cx = CGAL::to_double(op.arc.circle.center().x());
            const double cy = CGAL::to_double(op.arc.circle.center().y());
            const double sx = CGAL::to_double(op.arc.start.x()) - cx;
            const double sy = CGAL::to_double(op.arc.start.y()) - cy;
            const double start_angle = std::atan2(sy, sx);

            const double sw = op.arc.sweep();
            const double signed_sweep = op.arc.is_clockwise() ? -sw : sw;
            const int n = std::max(2, static_cast<int>(std::ceil(sw * samples_per_radian)));

            // Skip first sample (junction duplicate) unless first operation
            const int i0 = first ? 0 : 1;
            for (int i = i0; i < n - 1; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(n - 1);
                const double theta = start_angle + signed_sweep * t;
                const double z = op.z_start + t * (op.z_end - op.z_start);
                pts.push_back({cx + r * std::cos(theta),
                               cy + r * std::sin(theta), z});
            }
            // Exact endpoint — avoids cos/sin drift
            pts.push_back({CGAL::to_double(op.arc.end.x()),
                           CGAL::to_double(op.arc.end.y()),
                           op.z_end});
        }
    }

    const int m = static_cast<int>(pts.size());
    compas::RowMatrixXd result(m, 3);
    for (int i = 0; i < m; ++i) {
        result(i, 0) = pts[i][0];
        result(i, 1) = pts[i][1];
        result(i, 2) = pts[i][2];
    }
    return result;
}
```

Note: `ToolpathPrimitive` is defined later (line 626) in a second anonymous namespace block.
`tessellate_operations` must be placed AFTER `ToolpathPrimitive` is defined — insert it
just before `pmp_trochoidal_mat_toolpath_circular` (line 650), inside the same namespace
block or after the closing `}` of the second anonymous namespace.

**Step 2: No test yet — this is internal. Verified via integration in Task 3.**

---

### Task 2: Wire tessellation into `pmp_trochoidal_mat_toolpath_circular`

**Files:**
- Modify: `src/toolpath.cpp:650-882` (function signature + return)
- Modify: `src/toolpath.h:58-80` (declaration)

**Step 1: Add `samples_per_radian` parameter and 6-tuple return**

In `toolpath.h`, change the declaration:

```cpp
std::tuple<compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd,
           compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd>
pmp_trochoidal_mat_toolpath_circular(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    double lead_in,
    double lead_out,
    bool link_paths,
    bool optimize_order,
    double cut_z,
    double clearance_z,
    bool has_clearance_z,
    bool retract_at_end,
    double samples_per_radian
);
```

In `toolpath.cpp`, apply matching changes to the function definition (line 650).

**Step 2: Update empty returns**

Change both early-return sites (lines 679-682, 701-705) to return 6-tuple:

```cpp
compas::RowMatrixXd empty_meta(0, 4);
compas::RowMatrixXd empty3(0, 3);
compas::RowMatrixXd empty1(0, 1);
return std::make_tuple(empty_meta, empty3, empty3, empty3, empty1, empty3);
```

**Step 3: Tessellate and return**

Replace the final `return` (line 882) with:

```cpp
auto polyline = tessellate_operations(operations, samples_per_radian);
return std::make_tuple(meta, starts, ends, centers_out, radii_out, polyline);
```

**Step 4: Update nanobind binding**

In the `NB_MODULE` block (line 945-972), add the new parameter and update docstring:

```cpp
"- polyline (Mx3 float): tessellated 3D point sequence\n",
```

Add after `"retract_at_end"_a = true`:

```cpp
"samples_per_radian"_a = 10.0
```

**Step 5: Build**

Run: `pip install --no-build-isolation -ve .`
Expected: clean build

---

### Task 3: Python `ToolpathResult` and updated return

**Files:**
- Modify: `src/compas_cgal/toolpath.py`
- Test: `tests/test_toolpath.py`

**Step 1: Write the failing test**

Add to `tests/test_toolpath.py`:

```python
def test_circular_returns_toolpath_result():
    """trochoidal_mat_toolpath_circular returns ToolpathResult with operations + polyline."""
    from compas_cgal.toolpath import ToolpathResult

    result = trochoidal_mat_toolpath_circular(
        SQUARE, tool_diameter=2.0, pitch=1.0,
        max_trochoid_radius=float("inf"), max_passes=20,
    )
    assert isinstance(result, ToolpathResult)
    assert len(result.operations) > 0
    assert isinstance(result.operations[0], ToolpathOperation)
    assert result.polyline.ndim == 2
    assert result.polyline.shape[1] == 3
    assert result.polyline.shape[0] >= 2
    assert np.isfinite(result.polyline).all()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_toolpath.py::test_circular_returns_toolpath_result -v`
Expected: FAIL (ToolpathResult doesn't exist yet)

**Step 3: Implement `ToolpathResult` and update function**

In `src/compas_cgal/toolpath.py`:

Add to `__all__`:
```python
__all__ = [
    "ToolpathOperation",
    "ToolpathResult",
    ...
]
```

Add dataclass after `ToolpathOperation`:
```python
@dataclass
class ToolpathResult:
    """Toolpath output: typed operations for G-code + tessellated polyline for visualization."""

    operations: list[ToolpathOperation]
    polyline: np.ndarray
```

Update function signature — add `samples_per_radian: float = 10.0` parameter.

Update the C++ call to pass the new parameter:
```python
meta, starts, ends, centers, radii, polyline = _toolpath.trochoidal_mat_toolpath_circular(
    V, ..., bool(retract_at_end), float(samples_per_radian),
)
```

Parse polyline and return `ToolpathResult`:
```python
polyline = np.asarray(polyline, dtype=np.float64, order="C")
return ToolpathResult(
    operations=_matrices_to_operations(meta, starts, ends, centers, radii),
    polyline=polyline,
)
```

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_toolpath.py::test_circular_returns_toolpath_result -v`
Expected: PASS

**Step 5: Commit**

```
git add src/toolpath.h src/toolpath.cpp src/compas_cgal/toolpath.py tests/test_toolpath.py
git commit -m "feat: C++ tessellated polyline output for circular toolpath"
```

---

### Task 4: Fix existing tests for new return type

**Files:**
- Modify: `tests/test_toolpath.py` (all tests calling `trochoidal_mat_toolpath_circular`)

**Step 1: Update all call sites**

Every test that does `ops = trochoidal_mat_toolpath_circular(...)` must change to
`result = trochoidal_mat_toolpath_circular(...)` and use `result.operations`.

Tests affected:
- `test_circular_returns_toolpath_operations` → `result.operations`
- `test_circular_with_leads_and_links` → `result.operations`
- `test_circular_with_clearance_z` → `result.operations`
- `test_circular_continuity` → `result.operations`
- `test_full_circle_output` → `result.operations`
- `test_full_circle_sweep` → `result.operations`
- `test_circular_exit_tangent_alignment` → `result.operations`
- `test_tangent_continuity_known_polygons` → `result.operations`
- `test_tangent_continuity_random_polygons` → `result.operations`

**Step 2: Run full suite**

Run: `pytest tests/test_toolpath.py -v`
Expected: 17 passed (+ new test = 18)

**Step 3: Commit**

```
git commit -m "fix: update tests for ToolpathResult return type"
```

---

### Task 5: Polyline continuity test

**Files:**
- Modify: `tests/test_toolpath.py`

**Step 1: Write the failing test**

```python
def test_polyline_continuity():
    """Tessellated polyline has no positional jumps between consecutive points."""
    result = trochoidal_mat_toolpath_circular(
        IRREGULAR, tool_diameter=0.5, pitch=0.4,
        lead_in=0.15, lead_out=0.15, link_paths=True,
        optimize_order=True, cut_z=-0.2, clearance_z=2.0,
    )
    pts = result.polyline
    diffs = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    max_step = diffs.max()
    # No single step should exceed pitch + tool_diameter (generous bound)
    assert max_step < 3.0, f"Max step {max_step:.4f} exceeds bound"
    # No zero-length steps (dedup should eliminate exact duplicates)
    assert (diffs > 1e-14).all(), "Found duplicate consecutive points"


@given(polygon=_simple_polygons())
@settings(max_examples=50, deadline=None)
def test_polyline_continuity_random(polygon):
    """Property: tessellated polyline has no jumps for any simple polygon."""
    result = trochoidal_mat_toolpath_circular(
        polygon, tool_diameter=0.5, pitch=0.4,
        lead_in=0.15, lead_out=0.15, link_paths=True,
        optimize_order=True, cut_z=-0.2, clearance_z=2.0,
    )
    assume(result.polyline.shape[0] > 2)
    diffs = np.linalg.norm(np.diff(result.polyline, axis=0), axis=1)
    assert diffs.max() < 3.0
    assert (diffs > 1e-14).all()
```

**Step 2: Run to verify**

Run: `pytest tests/test_toolpath.py::test_polyline_continuity tests/test_toolpath.py::test_polyline_continuity_random -v`
Expected: PASS (the C++ tessellation should produce clean polylines)

**Step 3: Commit**

```
git commit -m "test: polyline continuity via hypothesis"
```

---

### Task 6: Update example to use C++-tessellated polyline

**Files:**
- Modify: `docs/examples/example_toolpath_trochoidal_mat.py`

**Step 1: Replace `_pts()` + heuristic reversal with `result.polyline`**

Change the toolpath call:
```python
result = trochoidal_mat_toolpath_circular(...)
operations = result.operations
```

Replace the entire `_pts` / `_dist_sq` / `path_points` assembly block with:
```python
path_points = result.polyline.tolist()
```

Remove `_pts`, `_dist_sq`, `MIN_PTS`, `PTS_PER_UNIT` — all dead code.

Update `groups` loop to use `operations` instead of `result`:
```python
for op in operations:
    groups[op.path_index].append(op)
```

**Step 2: Verify visually** (manual — run the example)

Run: `python docs/examples/example_toolpath_trochoidal_mat.py`
Expected: slider moves tool smoothly along the toolpath

**Step 3: Commit**

```
git commit -m "example: use C++ tessellated polyline, remove heuristic reversal"
```

---

### Task 7: Update design note

**Files:**
- Modify: `docs/toolpath_tangent_continuity.md`

Add a section noting that the winding representation gap is now bypassed for
visualization/verification by the C++-tessellated polyline, which uses the
correct `signed_sweep` from the `TrochoidArc` directly.

**Commit:**

```
git commit -m "docs: note C++ polyline bypasses winding repr gap"
```
