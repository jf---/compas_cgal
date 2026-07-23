# SP1 — Exact Stock Model + TEA Kernel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Exact in-process 2D stock model (CGAL circle-segment arrangement via the Boolean layer) with a tool-engagement-angle (TEA) kernel, a Monte-Carlo differential oracle, an engagement audit of the existing trochoid generator, and a measured wall-clock exit gate.

**Architecture:** `Stock2` holds remaining material as `General_polygon_set_2<Gps_circle_segment_traits_2<Epeck>>` (the regularized-boolean layer of the 2D Arrangements package owns face-containment bookkeeping); sweep regions are subtracted per committed motion — capsules and disks exact, arc sweeps as certified under-covering capsule chains (tool-offset radii of rational-r² arcs have irrational squares). TEA at a station = angular extent of the cutter circle inside material, computed exactly via GPS intersection + supporting-curve filtering; the hard cap is decided by exact chord predicates, band statistics are reported in doubles. Everything is test-driven from Python against the nanobind bindings.

**Tech Stack:** C++20, CGAL 6.0.1 (`Boolean_set_operations_2`, `Arr_circle_segment_traits_2`, Epeck), nanobind (NB_STATIC — one module for all SP1 C++), Eigen, scikit-build-core; Python 3.12, numpy, hypothesis, pytest-xdist.

## Global Constraints

- Repo discipline (CLAUDE.md): exact predicates for geometric decisions, kernel value types, no epsilon decisions; construction guards allowed only as named constants with unit + derivation comments.
- Conservatism rule (spec §3): every approximation in the stock model must UNDER-cover removal (over-estimate material). TEA and residual are then over-estimated, keeping both certificates sound.
- Never modify existing passing tests to make new code pass; the 45 toolpath tests must stay green.
- Python: ruff clean, Google docstrings, no `__all__` in `__init__.py`, fail loud (named exceptions), no fallbacks/conditional imports.
- pytest always with `-n auto`.
- **Environment (all commands run from repo root `/Users/jelle/Code/CADCAM/compas_cgal_prs`):**
  - `PY=/Users/jelle/mambaforge/envs/cgal-dev/bin/python`
  - Test: `PYTHONPATH=src $PY -m pytest tests/test_stock.py -n auto -q`
  - Rebuild after C++ edits: `cmake -S . -B build/cp312-abi3-macosx_15_0_arm64 && cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _stock_2 -j 8`
  - First build of the new module: after the cmake step, `ln -sf "$PWD/build/cp312-abi3-macosx_15_0_arm64/_stock_2.abi3.so" src/compas_cgal/_stock_2.abi3.so` (matches the existing symlink pattern; `.so` is gitignored).
- Commits: concise messages, no attribution lines, author is the configured repo user.

---

### Task 1: `_stock_2` module scaffold — `Stock2` init + `contains`

**Files:**
- Create: `src/stock_2.h`
- Create: `src/stock_2.cpp`
- Create: `src/engagement_2.h` (declarations only this task)
- Create: `src/engagement_2.cpp` (registration stub this task)
- Modify: `CMakeLists.txt` (one line in the module list)
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: `src/compas.h` (`compas::RowMatrixXd`), vendored CGAL.
- Produces (Python, module `compas_cgal._stock_2`):
  - `Stock2(boundary: float64[N,3], holes: list[float64[M,3]])` — CCW outer, CCW holes (converted internally), raises `ValueError` on non-simple input.
  - `Stock2.contains(x: float, y: float) -> bool` — point strictly inside remaining material.
  - `Stock2.is_empty() -> bool`
- Produces (C++, `stock_2.h`, used by every later task):
  ```cpp
  typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
  typedef CGAL::Gps_circle_segment_traits_2<Epeck> GpsTraits;
  typedef CGAL::General_polygon_set_2<GpsTraits> Gps;
  typedef GpsTraits::Polygon_2 GpsPolygon;             // circle-segment general polygon
  typedef GpsTraits::Polygon_with_holes_2 GpsPolygonWithHoles;
  typedef GpsTraits::Point_2 GpsPoint;                  // one-root coordinates
  typedef GpsTraits::X_monotone_curve_2 GpsXCurve;

  class Stock2 {
  public:
      Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
             const std::vector<compas::RowMatrixXd>& holes);
      bool contains(double x, double y) const;
      bool is_empty() const;
      const Gps& set() const { return set_; }           // engagement kernel reads this
  private:
      Gps set_;
  };
  ```

- [ ] **Step 1: Write the failing tests**

Create `tests/test_stock.py`:

```python
import numpy as np
import pytest

from compas_cgal import _stock_2

SQUARE = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=np.float64)
ISLAND = np.array([[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]], dtype=np.float64)


def test_stock_init_square():
    stock = _stock_2.Stock2(SQUARE, [])
    assert stock.contains(5.0, 5.0)
    assert not stock.contains(-1.0, 5.0)
    assert not stock.contains(11.0, 5.0)
    assert not stock.is_empty()


def test_stock_init_with_island():
    stock = _stock_2.Stock2(SQUARE, [ISLAND])
    assert stock.contains(2.0, 2.0)
    assert not stock.contains(5.0, 5.0)  # inside island = void


def test_stock_rejects_non_simple():
    bowtie = np.array([[0, 0, 0], [6, 0, 0], [0, 4, 0], [6, 4, 0]], dtype=np.float64)
    with pytest.raises(ValueError, match="simple"):
        _stock_2.Stock2(bowtie, [])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTHONPATH=src $PY -m pytest tests/test_stock.py -n auto -q`
Expected: FAIL/ERROR with `ImportError: cannot import name '_stock_2'`.

- [ ] **Step 3: Write the C++ skeleton**

`src/stock_2.h`:

```cpp
#pragma once

#include "compas.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
typedef CGAL::Gps_circle_segment_traits_2<Epeck> GpsTraits;
typedef CGAL::General_polygon_set_2<GpsTraits> Gps;
typedef GpsTraits::Polygon_2 GpsPolygon;
typedef GpsTraits::Polygon_with_holes_2 GpsPolygonWithHoles;
typedef GpsTraits::Point_2 GpsPoint;
typedef GpsTraits::X_monotone_curve_2 GpsXCurve;
typedef Epeck::Point_2 EPoint;
typedef Epeck::Circle_2 ECircle;
typedef Epeck::Segment_2 ESegment;

class Stock2 {
public:
    Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes);

    bool contains(double x, double y) const;
    bool is_empty() const;

    const Gps& set() const { return set_; }
    Gps& set() { return set_; }

private:
    Gps set_;
};
```

`src/stock_2.cpp` (module macro lives here; engagement registration is a
declared hook so both sources share one nanobind module — NB_STATIC forbids
cross-module type sharing):

```cpp
#include "stock_2.h"
#include "engagement_2.h"

#include <stdexcept>

namespace {

// Convert an Nx3 double matrix (rationals by construction) to a linear
// circle-segment general polygon, validating simplicity via Polygon_2.
GpsPolygon data_to_gps_polygon(Eigen::Ref<const compas::RowMatrixXd> vertices)
{
    if (vertices.rows() < 3) {
        throw std::invalid_argument("Expected at least three polygon vertices.");
    }
    CGAL::Polygon_2<Epeck> simple_check;
    for (int i = 0; i < vertices.rows(); ++i) {
        simple_check.push_back(EPoint(vertices(i, 0), vertices(i, 1)));
    }
    if (!simple_check.is_simple()) {
        throw std::invalid_argument("Polygon boundary must be simple (no self-intersections).");
    }
    if (simple_check.is_clockwise_oriented()) {
        simple_check.reverse_orientation();
    }

    GpsPolygon polygon;
    const std::size_t n = simple_check.size();
    for (std::size_t i = 0; i < n; ++i) {
        const EPoint& a = simple_check[i];
        const EPoint& b = simple_check[(i + 1) % n];
        polygon.push_back(GpsXCurve(ESegment(a, b)));
    }
    return polygon;
}

} // namespace

Stock2::Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
               const std::vector<compas::RowMatrixXd>& holes)
{
    set_.insert(data_to_gps_polygon(boundary));
    for (const auto& hole : holes) {
        Gps hole_set;
        hole_set.insert(data_to_gps_polygon(hole));
        set_.difference(hole_set);
    }
}

bool Stock2::contains(double x, double y) const
{
    return set_.oriented_side(GpsPoint(Epeck::FT(x), Epeck::FT(y))) == CGAL::ON_POSITIVE_SIDE;
}

bool Stock2::is_empty() const
{
    return set_.is_empty();
}

NB_MODULE(_stock_2, m)
{
    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty);

    register_engagement(m);
}
```

`src/engagement_2.h`:

```cpp
#pragma once

#include "compas.h"

class Stock2;

void register_engagement(nanobind::module_& m);
```

`src/engagement_2.cpp`:

```cpp
#include "engagement_2.h"
#include "stock_2.h"

void register_engagement(nanobind::module_& m)
{
    (void)m;  // engagement bindings land in Task 4
}
```

Note on `GpsPoint` construction: `Gps_circle_segment_traits_2` points take
one-root coordinates constructible from `Epeck::FT`. If the two-argument
constructor does not compile against the vendored CGAL, use
`GpsPoint(GpsTraits::CoordNT(Epeck::FT(x)), GpsTraits::CoordNT(Epeck::FT(y)))`.

- [ ] **Step 4: Register the module and build**

In `CMakeLists.txt`, after `add_nanobind_module(_toolpath src/toolpath.cpp)` add:

```cmake
add_nanobind_module(_stock_2 "src/stock_2.cpp;src/engagement_2.cpp")
```

Run:
```bash
cmake -S . -B build/cp312-abi3-macosx_15_0_arm64
cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _stock_2 -j 8
ln -sf "$PWD/build/cp312-abi3-macosx_15_0_arm64/_stock_2.abi3.so" src/compas_cgal/_stock_2.abi3.so
```
Expected: clean build, symlink present.

- [ ] **Step 5: Run tests to verify they pass**

Run: `PYTHONPATH=src $PY -m pytest tests/test_stock.py -n auto -q`
Expected: 3 passed.

- [ ] **Step 6: Commit**

```bash
git add src/stock_2.h src/stock_2.cpp src/engagement_2.h src/engagement_2.cpp CMakeLists.txt tests/test_stock.py
git commit -m "feat(stock): Stock2 exact circle-segment stock set, init + contains"
```

---

### Task 2: Exact sweep subtraction — `subtract_capsule`, `subtract_disk`

**Files:**
- Modify: `src/stock_2.h`, `src/stock_2.cpp`
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: Task 1 `Stock2`, `GpsPolygon`, `GpsXCurve`.
- Produces:
  - `Stock2.subtract_capsule(x0, y0, x1, y1, radius: float) -> None` — removes the exact capsule (rectangle + two half-disks) swept by a tool of `radius` along the segment. Degenerate (coincident endpoints) delegates to `subtract_disk`.
  - `Stock2.subtract_disk(cx, cy, radius: float) -> None`
  - C++ helpers in anonymous namespace: `GpsPolygon capsule_polygon(const EPoint&, const EPoint&, const Epeck::FT& r)`, `GpsPolygon disk_polygon(const EPoint&, const Epeck::FT& r)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_stock.py`:

```python
def test_subtract_disk_clears_point():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)
    assert not stock.contains(5.9, 5.0)   # inside disk
    assert stock.contains(6.1, 5.0)       # outside disk
    assert stock.contains(2.0, 2.0)


def test_subtract_capsule_geometry():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    assert not stock.contains(5.0, 5.0)       # on the axis
    assert not stock.contains(5.0, 5.4)       # inside half-width
    assert stock.contains(5.0, 5.6)           # outside half-width
    assert not stock.contains(8.4, 5.0)       # inside end cap
    assert stock.contains(8.6, 5.0)           # outside end cap


def test_overlapping_subtractions_are_regularized():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(2.0, 5.0, 8.0, 5.0, 0.5)
    stock.subtract_capsule(2.0, 5.2, 8.0, 5.2, 0.5)  # heavy overlap
    stock.subtract_disk(5.0, 5.0, 0.5)               # fully inside cleared
    assert not stock.contains(5.0, 5.3)
    assert stock.contains(5.0, 6.0)


def test_subtract_degenerate_capsule_is_disk():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(5.0, 5.0, 5.0, 5.0, 1.0)
    assert not stock.contains(5.5, 5.0)
    assert stock.contains(6.5, 5.0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTHONPATH=src $PY -m pytest tests/test_stock.py -n auto -q`
Expected: 4 new FAIL with `AttributeError: ... no attribute 'subtract_disk'`.

- [ ] **Step 3: Implement**

Add to the anonymous namespace of `src/stock_2.cpp`:

```cpp
// Full disk as a two-arc general polygon (counterclockwise).
GpsPolygon disk_polygon(const EPoint& center, const Epeck::FT& radius)
{
    const Epeck::FT r_sq = radius * radius;
    const ECircle circle(center, r_sq, CGAL::COUNTERCLOCKWISE);
    // Split at the two x-extreme points (rational: center ± (r, 0) is NOT
    // rational when r is FT-exact? radius comes from a double => rational, so
    // center.x() ± radius is rational).
    const GpsPoint p_min(center.x() - radius, center.y());
    const GpsPoint p_max(center.x() + radius, center.y());
    GpsPolygon polygon;
    polygon.push_back(GpsXCurve(circle, p_min, p_max, CGAL::COUNTERCLOCKWISE));  // lower half
    polygon.push_back(GpsXCurve(circle, p_max, p_min, CGAL::COUNTERCLOCKWISE));  // upper half
    return polygon;
}

// Exact capsule: two straight edges offset by r plus two half-disk caps.
// All coordinates rational: doubles in, rational unit offsets require the
// segment direction normalized — avoided by constructing offset points as
// center ± r * (dy, -dx) / len. len is irrational in general, so instead the
// capsule is built as the regularized union of the axis-aligned construction:
//   capsule = disk(p0, r) ∪ disk(p1, r) ∪ oriented rectangle
// The rectangle needs the unit normal — irrational. Resolution: build the
// rectangle from the two tangency chords of the disks: its corners are
// p ± r*n with n the unit normal. To stay rational, the rectangle corners are
// computed with the SAME one-root machinery the traits uses: corner
// coordinates are (x ± r*dy/len, y ∓ r*dx/len) — one-root numbers over
// len² = dx²+dy² (rational). GpsPoint accepts one-root coordinates, and
// SEGMENT curves between one-root points are NOT supported by the traits
// (segment endpoints must be rational for Curve_2 construction from
// ESegment).
//
// Therefore the exact capsule uses the supported construction:
//   X_monotone_curve_2(Line_2, one-root source, one-root target)
// Gps_circle_segment_traits_2 supports segments on a rational LINE with
// one-root endpoints (this is exactly how it represents clipped segments).
// The supporting line of each capsule side is rational:
//   side line: dy*x - dx*y + (c ± r*len) = 0  -> NOT rational (r*len term).
// Since the side line itself is irrational, the fully exact oriented capsule
// is not representable either. CONSERVATIVE RESOLUTION (spec §3 rule):
// under-cover the capsule by the union of disks along the segment at spacing
// s with radius r' = r - sagitta_bound(s), plus the exact end disks — the
// same certified chain used for arc sweeps (Task 3), applied with an
// infinite-radius (straight) guide. Under-covering keeps all certificates
// sound; the slack is bounded by CHAIN_SAGITTA_FRACTION * radial slack.
GpsPolygon capsule_polygon(const EPoint& p0, const EPoint& p1, const Epeck::FT& r);
// (implemented in Task 3 as chain_sweep; Task 2 implements subtract_disk and
//  a subtract_capsule that forwards to the chain with straight-guide spacing)
```

**Correction incorporated (implementer: read this):** the analysis above shows
the *oriented* exact capsule is also outside the traits (irrational side
lines), which the design doc's conservatism rule already covers. Task 2
therefore implements `subtract_disk` exactly, and `subtract_capsule` as a
certified under-covering **disk chain**:

```cpp
// Named constants — units and derivation.
// Fraction of the tool radius allowed as chain under-coverage slack; the
// generator's radial_clearance margin (1e-3 * D = 2e-3 * r) dominates it.
constexpr double CHAIN_SLACK_FRACTION = 1e-4;

void Stock2::subtract_disk(double cx, double cy, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    Gps region;
    region.insert(disk_polygon(EPoint(cx, cy), Epeck::FT(radius)));
    set_.difference(region);
}

void Stock2::subtract_capsule(double x0, double y0, double x1, double y1, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double len_sq = dx * dx + dy * dy;
    if (len_sq == 0.0) {
        subtract_disk(x0, y0, radius);
        return;
    }
    // Disk chain spacing s with per-disk radius r: the chain covers the true
    // capsule shrunk by delta = r - sqrt(r^2 - (s/2)^2) <= s^2/(4r).
    // Choose s so delta <= CHAIN_SLACK_FRACTION * r  =>  s = 2r*sqrt(fraction).
    const double len = std::sqrt(len_sq);
    const double spacing = 2.0 * radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const int n = std::max(1, static_cast<int>(std::ceil(len / spacing)));
    Gps region;
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        region.join(/* accumulate */);
        // implementer: region.insert on first, region.join(Gps(disk)) after;
        // or collect polygons and use CGAL::join range overload.
        Gps disk_set;
        disk_set.insert(disk_polygon(EPoint(x0 + t * dx, y0 + t * dy), Epeck::FT(radius)));
        region.join(disk_set);
    }
    set_.difference(region);
}
```

Remove the placeholder `region.join(/* accumulate */);` line — the two
statements after it are the implementation. Bind both methods in `NB_MODULE`:

```cpp
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a)
```

Test-geometry note: with `CHAIN_SLACK_FRACTION = 1e-4` the under-coverage is
≤ `1e-4 * r` = 5e-5 for r = 0.5 — invisible at the 0.1-unit resolution the
tests probe, so the capsule assertions hold as written.

- [ ] **Step 4: Rebuild and run tests**

Run the rebuild command, then: `PYTHONPATH=src $PY -m pytest tests/test_stock.py -n auto -q`
Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add src/stock_2.cpp src/stock_2.h tests/test_stock.py
git commit -m "feat(stock): exact disk + certified under-covering capsule subtraction"
```

---

### Task 3: Arc sweeps — certified under-covering chains

**Files:**
- Modify: `src/stock_2.h`, `src/stock_2.cpp`
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: Task 2 `disk_polygon`, `subtract_disk`, `CHAIN_SLACK_FRACTION`.
- Produces:
  - `Stock2.subtract_arc_sweep(cx, cy, sx, sy, ex, ey, cw: bool, tool_radius: float) -> None` — removes (an under-cover of) the region swept by the tool disk along the circular arc from `(sx,sy)` to `(ex,ey)` about `(cx,cy)`; full circle when start == end.
  - Internal: `void subtract_point_chain(const std::vector<std::pair<double,double>>&, double radius)` shared by capsule and arc paths.

- [ ] **Step 1: Write the failing tests**

```python
import math


def test_subtract_arc_sweep_full_circle():
    stock = _stock_2.Stock2(SQUARE, [])
    # tool r=0.4 swept around full circle of radius 2 centered at (5,5)
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    assert not stock.contains(7.0, 5.0)        # on guide circle
    assert not stock.contains(5.0, 7.0)        # opposite side of circle
    assert not stock.contains(7.3, 5.0)        # within tool band (outer)
    assert not stock.contains(6.7, 5.0)        # within tool band (inner)
    assert stock.contains(5.0, 5.0)            # annulus center survives
    assert stock.contains(7.6, 5.0)            # outside band + slack


def test_subtract_arc_sweep_quarter_arc():
    stock = _stock_2.Stock2(SQUARE, [])
    # CCW quarter arc from (7,5) to (5,7) about (5,5), tool r=0.4
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 5.0, 7.0, False, 0.4)
    mid = (5.0 + 2.0 * math.cos(math.pi / 4), 5.0 + 2.0 * math.sin(math.pi / 4))
    assert not stock.contains(*mid)
    assert stock.contains(5.0, 3.0)            # untouched opposite quadrant
    assert stock.contains(3.0, 5.0)


def test_arc_sweep_under_covers_never_over():
    """Conservatism: nothing farther than tool_radius from the guide arc is removed."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_arc_sweep(5.0, 5.0, 7.0, 5.0, 7.0, 5.0, True, 0.4)
    rng = np.random.default_rng(3)
    for _ in range(400):
        x, y = rng.uniform(0.2, 9.8, 2)
        d_guide = abs(math.hypot(x - 5.0, y - 5.0) - 2.0)
        if d_guide > 0.4:               # strictly outside the true sweep
            assert stock.contains(x, y), (x, y, d_guide)
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: `AttributeError: ... no attribute 'subtract_arc_sweep'`.

- [ ] **Step 3: Implement**

```cpp
// Subtract the union of exact tool disks centered at the given points.
void Stock2::subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                                  double radius)
{
    Gps region;
    for (const auto& [x, y] : centers) {
        Gps disk_set;
        disk_set.insert(disk_polygon(EPoint(x, y), Epeck::FT(radius)));
        region.join(disk_set);
    }
    set_.difference(region);
}

void Stock2::subtract_arc_sweep(double cx, double cy, double sx, double sy,
                                double ex, double ey, bool cw, double tool_radius)
{
    if (tool_radius <= 0.0) throw std::invalid_argument("tool_radius should be positive.");
    const double rx = sx - cx, ry = sy - cy;
    const double guide_r = std::hypot(rx, ry);
    if (guide_r == 0.0) { subtract_disk(cx, cy, tool_radius); return; }

    double a0 = std::atan2(ry, rx);
    double a1 = std::atan2(ey - cy, ex - cx);
    double sweep = cw ? a0 - a1 : a1 - a0;
    if (sweep <= 0.0) sweep += 2.0 * std::numbers::pi;
    const bool full = (sx == ex && sy == ey);
    if (full) sweep = 2.0 * std::numbers::pi;

    // Chain spacing along the guide: disks of tool_radius at arc-length step s
    // under-cover the true sweep by delta <= s^2/(4*tool_radius) (chord
    // sagitta bound; straight-chord bound is conservative for the arc chain
    // because consecutive centers are closer along the chord than the arc).
    const double spacing = 2.0 * tool_radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const double arc_len = guide_r * sweep;
    const int n = std::max(4, static_cast<int>(std::ceil(arc_len / spacing)));

    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double a = cw ? a0 - sweep * t : a0 + sweep * t;
        centers.emplace_back(cx + guide_r * std::cos(a), cy + guide_r * std::sin(a));
    }
    subtract_point_chain(centers, tool_radius);
}
```

Refactor `subtract_capsule` (Task 2) to build its center list and call
`subtract_point_chain` — one chain implementation, two callers. Bind:

```cpp
        .def("subtract_arc_sweep", &Stock2::subtract_arc_sweep,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
```

- [ ] **Step 4: Rebuild, run tests**

Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/stock_2.cpp src/stock_2.h tests/test_stock.py
git commit -m "feat(stock): certified under-covering arc sweep chains"
```

---

### Task 4: TEA at a station — exact engagement query

> **AMENDMENT (2026-07-23, supersedes the decision-core details below):** the
> original text below contains two exact-kernel violations — a
> `CAP_BOUND_DEFLATION` deflated-double cap comparison and a `1e-12`
> double-angle merge gap. Both are FORBIDDEN by the repo CLAUDE.md
> "Exact-Kernel Discipline" section. The binding design is:
>
> 1. **API contract:** the C++ signature is
>    `engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio)` where
>    `cap_chord_ratio = 4·sin²(cap/2)`, dimensionless, in `(0, 4]`, computed
>    by the CALLER (tests: `4.0 * math.sin(cap/2.0)**2`) and injected exactly
>    (`Epeck::FT(double)`). Caps are contractually `0 < cap ≤ π`. The
>    certificate is an exact statement about the exact rational threshold
>    `T = FT(cap_chord_ratio) · FT(r)²`; the sub-ulp gap to the transcendental
>    angle is documented API semantics at the parameter, not a correction
>    constant. NO deflation factor exists anywhere.
> 2. **Run merging is exact:** abutting sub-arcs share their endpoint one-root
>    point exactly (same arrangement vertex); merge maximal runs by exact
>    CoordNT point equality (`GpsPoint operator==`), including the
>    wrap-around pair. Delete `MERGE_GAP`/1e-12 entirely. Angular sorting for
>    run assembly uses exact point comparisons (half-plane wrt the rational
>    horizontal line through the center via `CGAL::compare(y, c_y)`, then
>    `compare_x` within each half) — `atan2` doubles are used ONLY for the
>    reported angles, never for adjacency or ordering decisions.
> 3. **Cap decision is an exact predicate.** For a maximal run with endpoints
>    p, q (one-root, possibly different roots α, β) on the cutter circle:
>    - angle(run) > π ⟺ exact orientation test of (center, start, end)
>      honoring the run's winding; any >π run violates (cap ≤ π). Equality
>      (collinear/antipodal) means angle = π exactly: violation iff
>      `cap_chord_ratio < 4` exactly.
>    - angle ≤ π: violation ⟺ squared chord |pq|² > T, decided exactly.
>    Both reduce to the sign of `A + B√α + C√β + D√(αβ)` with RATIONAL
>    A,B,C,D extracted via `Sqrt_extension::a0()/a1()/root()` from the point
>    coordinates. Implement one exact primitive
>    `CGAL::Sign sign_mixed_radical(A, B, C, D, alpha, beta)` by repeated
>    squaring over FT only (group `u = A + B√α`, `w = C + D√α`; the sign of
>    `u + √β·w` follows from sign(u), sign(w), and the one-root comparison
>    `u²  vs  β·w²`, which is again a 2-term `A' + B'√α` sign — bottoming out
>    in rational signs). Cross-root Sqrt_extension ARITHMETIC is never
>    performed (documented UB); only rational coefficient arithmetic.
>    Unit-test the primitive directly against high-precision references
>    (`Fraction`-based cases in Python through a test-only binding
>    `_sign_mixed_radical(a, b, c, d, alpha, beta)` exposed on the module).
> 4. **Reporting stays double** (`atan2` for total/max angles) and never
>    feeds decisions.
> 5. **Test changes vs the text below:** every `engagement_at(...)` call
>    passes `4.0 * math.sin(cap / 2.0) ** 2` instead of `cap`; the
>    `math.radians(181)` case is replaced by `cap = math.pi` (caps are ≤ π;
>    a ≈π engaged run under chain under-coverage measures strictly below π,
>    so `not exceeded` holds); add one test asserting `engagement_at` raises
>    `ValueError` for `cap_chord_ratio <= 0` or `> 4`.
>
> Where the text below conflicts with this amendment, the amendment governs.

**Files:**
- Modify: `src/engagement_2.h`, `src/engagement_2.cpp`
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: `Stock2::set()` (`const Gps&`), `disk_polygon`.
- Produces (Python):
  - `_stock_2.engagement_at(stock: Stock2, cx, cy, tool_radius, cap_radians) -> tuple[float, float, bool]` — `(total_tea, max_run_tea, cap_exceeded)`. Angles in radians reported as doubles; `cap_exceeded` decided EXACTLY: any maximal engaged run spanning > π with `cap ≤ π` violates outright; runs ≤ π compared by squared-chord predicate against a conservatively rounded-down bound.
- Produces (C++, `engagement_2.h`):
  ```cpp
  struct EngagementSample { double total_tea; double max_run_tea; bool cap_exceeded; };
  EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                                 double tool_radius, double cap_radians);
  ```

- [ ] **Step 1: Write the failing tests**

```python
def test_engagement_full_material():
    stock = _stock_2.Stock2(SQUARE, [])
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.pi)
    assert total == pytest.approx(2.0 * math.pi, abs=1e-9)
    assert exceeded  # 2π run > π cap, decided exactly


def test_engagement_zero_in_cleared():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.pi)
    assert total == pytest.approx(0.0, abs=1e-9)
    assert not exceeded


def test_engagement_half_plane():
    """Cutter centered on the boundary of a cleared half: TEA ≈ π."""
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 7.0, 10.0, 7.0, 2.0)  # clears the band y in [5, 9]
    total, max_run, exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.pi)
    assert total == pytest.approx(math.pi, abs=0.02)   # chain slack ≪ tolerance
    assert max_run == pytest.approx(math.pi, abs=0.02)


def test_engagement_cap_exact_flag():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(0.0, 7.0, 10.0, 7.0, 2.0)
    # cap 100° < measured ~180° -> exceeded; cap 175°+ margin -> also exceeded;
    # cap just above π/2? measured run π > cap -> exceeded
    _, _, exceeded_tight = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.radians(100))
    assert exceeded_tight
    _, _, exceeded_loose = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.radians(179.9))
    assert exceeded_loose  # run is ~180.0°, slack pushes it above 179.9°? NO —
    # chain slack makes material SMALLER never larger; measured run ≤ true π.
    # Use a strictly-below cap for the positive case and re-check:
```

Replace the last three lines of `test_engagement_cap_exact_flag` with:

```python
    _, _, not_exceeded = _stock_2.engagement_at(stock, 5.0, 5.0, 0.5, math.radians(181))
    assert not not_exceeded
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: `AttributeError: module '_stock_2' has no attribute 'engagement_at'`.

- [ ] **Step 3: Implement in `src/engagement_2.cpp`**

Algorithm (GPS-intersection route — documented, no zone dependency):

```cpp
#include "engagement_2.h"
#include "stock_2.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace {

// Conservative rounding inflation for the double-computed chord-bound
// constant 4*r^2*sin^2(cap/2): one part in 1e12 covers the few ulps of
// std::sin/std::pow error, biased DOWN so violations are never missed.
constexpr double CAP_BOUND_DEFLATION = 1.0 - 1e-12;

struct AngularArc { double a0; double a1; };  // CCW from a0 to a1, radians

} // namespace

EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_radians)
{
    // 1. region = stock ∩ tool disk  (regularized, exact)
    Gps disk_set;
    disk_set.insert(/* disk_polygon(EPoint(cx, cy), Epeck::FT(tool_radius)) —
        move disk_polygon to stock_2.h as a free declaration so both TUs share it */);
    Gps region = stock.set();
    region.intersection(disk_set);

    EngagementSample out{0.0, 0.0, false};
    if (region.is_empty()) return out;

    // 2. Engaged arcs = boundary curves of the intersection supported by the
    //    cutter circle (rim-in-material pieces).
    std::vector<GpsPolygonWithHoles> pwhs;
    region.polygons_with_holes(std::back_inserter(pwhs));

    const Epeck::FT r_sq = Epeck::FT(tool_radius) * Epeck::FT(tool_radius);
    const EPoint center(cx, cy);

    std::vector<AngularArc> arcs;
    auto harvest = [&](const GpsPolygon& poly) {
        for (auto cit = poly.curves_begin(); cit != poly.curves_end(); ++cit) {
            const GpsXCurve& xc = *cit;
            if (!xc.is_circular()) continue;
            if (xc.supporting_circle().center() != center) continue;
            if (xc.supporting_circle().squared_radius() != r_sq) continue;
            const double sx = CGAL::to_double(xc.source().x());
            const double sy = CGAL::to_double(xc.source().y());
            const double tx = CGAL::to_double(xc.target().x());
            const double ty = CGAL::to_double(xc.target().y());
            double a0 = std::atan2(sy - cy, sx - cx);
            double a1 = std::atan2(ty - cy, tx - cx);
            // orientation: boundary of the intersection is CCW on outer
            // boundaries; on hole boundaries the rim piece is CW — normalize
            // every piece to its CCW angular interval:
            if (xc.orientation() == CGAL::CLOCKWISE) std::swap(a0, a1);
            arcs.push_back({a0, a1});
        }
    };
    for (const auto& pwh : pwhs) {
        harvest(pwh.outer_boundary());
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) harvest(*hit);
    }
    if (arcs.empty()) {
        // No rim pieces: the rim is either fully in material (region touches
        // the full circle) or fully clear. Decide by one rim point.
        const bool rim_in = stock.contains(cx + tool_radius, cy);
        if (rim_in) {
            out.total_tea = 2.0 * std::numbers::pi;
            out.max_run_tea = 2.0 * std::numbers::pi;
            out.cap_exceeded = (cap_radians < 2.0 * std::numbers::pi * CAP_BOUND_DEFLATION);
        }
        return out;
    }

    // 3. Merge angularly-adjacent pieces into maximal runs (endpoints of
    //    x-monotone splits coincide exactly; compare in doubles with the
    //    knowledge that true gaps are >> ulp because distinct crossings are
    //    distinct arrangement vertices).
    auto span = [](const AngularArc& a) {
        double s = a.a1 - a.a0;
        if (s <= 0.0) s += 2.0 * std::numbers::pi;
        return s;
    };
    std::sort(arcs.begin(), arcs.end(),
              [](const AngularArc& u, const AngularArc& v) { return u.a0 < v.a0; });
    std::vector<AngularArc> runs;
    for (const auto& a : arcs) {
        if (!runs.empty()) {
            double gap = a.a0 - runs.back().a1;
            if (std::abs(gap) < 1e-12 /* MERGE_GAP: split-point reunification,
                    pure double-identity check, see header comment */) {
                runs.back().a1 = a.a1;
                continue;
            }
        }
        runs.push_back(a);
    }
    // wrap-around merge
    if (runs.size() > 1) {
        double gap = (runs.front().a0 + 2.0 * std::numbers::pi) - runs.back().a1;
        if (std::abs(gap) < 1e-12) {
            runs.back().a1 = runs.front().a1 + 2.0 * std::numbers::pi;
            runs.erase(runs.begin());
        }
    }

    // 4. Totals + exact-leaning cap decision.
    const double cap_bound = cap_radians * CAP_BOUND_DEFLATION;
    for (const auto& r : runs) {
        const double s = span(r);
        out.total_tea += s;
        out.max_run_tea = std::max(out.max_run_tea, s);
        if (s > cap_bound) out.cap_exceeded = true;
    }
    return out;
}

void register_engagement(nanobind::module_& m)
{
    m.def("engagement_at",
          [](const Stock2& stock, double cx, double cy, double tool_radius, double cap) {
              EngagementSample s = engagement_at(stock, cx, cy, tool_radius, cap);
              return std::make_tuple(s.total_tea, s.max_run_tea, s.cap_exceeded);
          },
          "stock"_a, "cx"_a, "cy"_a, "tool_radius"_a, "cap_radians"_a);
}
```

Implementer notes (read before coding):
- Move `disk_polygon` from the anonymous namespace of `stock_2.cpp` into
  `stock_2.h` as a declared free function so `engagement_2.cpp` can use it.
- `GpsXCurve` API: `is_circular()`, `supporting_circle()`, `source()`,
  `target()`, `orientation()` — all provided by
  `Arr_circle_segment_traits_2::X_monotone_curve_2`.
- The rational supporting-circle equality test (`center`, `squared_radius`)
  is exact — the cutter disk's circle data is rational.
- The plan's cap decision uses double angles with a deflated bound rather than
  the chord predicate; the chord-predicate upgrade is noted as a follow-up in
  the SP2 spec if the deflated-double margin (1e-12 relative) is ever the
  deciding factor — at TEA scales (radians ~1) it is 12 orders below process
  noise. Document this in a comment verbatim.

- [ ] **Step 4: Rebuild, run tests**

Expected: all stock tests pass (14).

- [ ] **Step 5: Commit**

```bash
git add src/engagement_2.cpp src/engagement_2.h src/stock_2.h src/stock_2.cpp tests/test_stock.py
git commit -m "feat(engagement): exact station TEA via GPS intersection rim harvest"
```

---

### Task 5: TEA certificate along a motion

> **AMENDMENT (2026-07-23, supersedes conflicting details below):** the
> certificate must respect the Exact-Kernel Discipline. Redesign:
>
> - **Stations decide exactly.** Every station measurement is the Task 4
>   (amended) exact cap predicate — but evaluated against a GUARDED
>   threshold: at refinement level with half-spacing `d`, the station test
>   uses `cap_guarded = cap − guard(d, r)` converted to its chord ratio by
>   the caller-side rule (C++ receives the guarded ratio per level from its
>   own conversion — see below).
> - **The guard is an analytic lemma bound with proof-level slack** (per the
>   CLAUDE.md analytic-bounds clause):
>   `guard(d, r) = 2·(4·asin(min(1, d/(2r))) + 2·acos(max(−1, 1 − d/r)))` —
>   the Task-5 endpoint-drift + newborn-contact bound at an explicit integer
>   safety factor 2, stated in the derivation comment with the safe failure
>   direction (a too-large guard causes more refinement or a false violation,
>   never a false certificate). Double evaluation of the guard is sanctioned
>   by the slack factor; the STATION test itself remains the exact predicate.
> - `certify_segment_tea(stock, x0, y0, x1, y1, tool_radius, cap_radians)`
>   keeps the ergonomic radians signature at the BINDING layer only: the
>   lambda validates `0 < cap ≤ π`, and per recursion level computes
>   `guarded_ratio = 4·sin²((cap − guard)/2)` as a double, injected exactly —
>   this is threshold SELECTION at the declared boundary (documented at the
>   parameter), after which all station decisions are exact. `guard ≥ cap`
>   forces refinement (station test unsatisfiable at that spacing) until the
>   spacing floor, then reports uncertified.
> - Recursion floor and depth cap stay as below (`STATION_FLOOR_FRACTION`,
>   depth 24); `max_tea` in the result stays a reported double.
> - Tests: same three scenarios as below, with `cap` passed in radians to
>   `certify_segment_tea` (unchanged), and expectations unchanged — plus one
>   new test: a corridor whose width puts stations exactly at the guarded
>   boundary must return `cap_certified` False rather than a false pass
>   (construct by using a cap barely above the corridor's station TEA).
>
> Where the text below conflicts, this amendment governs.

**Files:**
- Modify: `src/engagement_2.h`, `src/engagement_2.cpp`
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: Task 4 `engagement_at`.
- Produces (Python):
  - `_stock_2.certify_segment_tea(stock, x0, y0, x1, y1, tool_radius, cap_radians) -> tuple[float, bool, int]` — `(max_tea_observed, cap_certified, stations_used)`. `cap_certified=True` means: no point of the motion exceeds the cap, established by station measurements + the conservative growth bound; `False` means a measured or unclosable-margin violation.
- Produces (C++): `CertifiedTea certify_segment_tea(const Stock2&, double, double, double, double, double, double);` with `struct CertifiedTea { double max_tea; bool cap_certified; int stations; };`

- [ ] **Step 1: Write the failing tests**

```python
def test_certify_segment_in_open_field():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_capsule(1.0, 5.0, 9.0, 5.0, 1.5)  # wide cleared corridor
    # re-cut along the corridor axis: rim never touches material
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 2.0, 5.0, 8.0, 5.0, 0.5, math.radians(120))
    assert ok
    assert max_tea == pytest.approx(0.0, abs=1e-9)


def test_certify_segment_flags_slotting():
    stock = _stock_2.Stock2(SQUARE, [])
    # virgin slotting move: TEA is π (~180°) > 120° cap
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 2.0, 5.0, 8.0, 5.0, 0.5, math.radians(120))
    assert not ok
    assert max_tea > math.radians(170)


def test_certify_stations_refine_near_feature():
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 2.0)  # cleared pocket in the middle
    # move crossing from cleared into material: certificate must refine
    max_tea, ok, n = _stock_2.certify_segment_tea(stock, 5.0, 5.0, 9.0, 5.0, 0.5, math.radians(120))
    assert not ok            # exits the cleared disk into full slotting
    assert n > 2             # refinement actually happened
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: `AttributeError ... certify_segment_tea`.

- [ ] **Step 3: Implement**

```cpp
namespace {

// Conservative TEA growth bound over center travel d with tool radius r,
// stock frozen. Two effects: (a) existing engagement-interval endpoints slide
// along the rim by at most the rim-arc subtended by travel d: 2*asin(min(1, d/(2r)))
// per endpoint, two endpoints -> 4*asin(min(1, d/(2r))); (b) a NEW contact
// born between stations can grow to at most the angular width of the chord
// cut by a feature at depth d into the rim: 2*acos(max(-1, 1 - d/r)).
// GROWTH = (a) + (b). Validated empirically against dense sampling in
// test_growth_bound_sound (Task 7); the certificate never relies on
// smoothness beyond this bound.
double tea_growth_bound(double d, double r)
{
    const double a = 4.0 * std::asin(std::min(1.0, d / (2.0 * r)));
    const double b = 2.0 * std::acos(std::max(-1.0, 1.0 - d / r));
    return a + b;
}

// Station-spacing floor as a fraction of tool radius: below this the growth
// bound is < 1e-3 rad and further refinement is numerically meaningless.
constexpr double STATION_FLOOR_FRACTION = 1e-3;

void certify_recursive(const Stock2& stock, double x0, double y0, double x1, double y1,
                       double r, double cap, CertifiedTea& acc, int depth)
{
    const double d = std::hypot(x1 - x0, y1 - y0);
    const EngagementSample s0 = engagement_at(stock, x0, y0, r, cap);
    // measure only the far endpoint; near endpoint measured by caller level
    const EngagementSample s1 = engagement_at(stock, x1, y1, r, cap);
    acc.stations += 1;
    acc.max_tea = std::max({acc.max_tea, s0.max_run_tea, s1.max_run_tea});
    if (s0.cap_exceeded || s1.cap_exceeded) { acc.cap_certified = false; return; }

    const double margin = cap - std::max(s0.max_run_tea, s1.max_run_tea);
    if (tea_growth_bound(0.5 * d, r) <= margin) return;      // certified between
    if (d < STATION_FLOOR_FRACTION * r || depth > 24) {
        acc.cap_certified = false;                            // unclosable margin
        return;
    }
    const double mx = 0.5 * (x0 + x1), my = 0.5 * (y0 + y1);
    certify_recursive(stock, x0, y0, mx, my, r, cap, acc, depth + 1);
    if (!acc.cap_certified) return;
    certify_recursive(stock, mx, my, x1, y1, r, cap, acc, depth + 1);
}

} // namespace

CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians)
{
    CertifiedTea acc{0.0, true, 0};
    certify_recursive(stock, x0, y0, x1, y1, tool_radius, cap_radians, acc, 0);
    return acc;
}
```

Bind in `register_engagement` (tuple return, same style as Task 4).
Note the double-measurement of shared endpoints across recursion levels is
accepted for SP1 clarity; memoization is an SP2 optimization if the gate needs it.

- [ ] **Step 4: Rebuild, run tests**

Expected: 17 passed.

- [ ] **Step 5: Commit**

```bash
git add src/engagement_2.cpp src/engagement_2.h tests/test_stock.py
git commit -m "feat(engagement): station certificate with conservative growth bound"
```

---

### Task 6: Python API — `stock.py`, `engagement.py`, audit replay

**Files:**
- Create: `src/compas_cgal/stock.py`
- Create: `src/compas_cgal/engagement.py`
- Test: `tests/test_engagement_audit.py`

**Interfaces:**
- Consumes: `_stock_2` bindings (Tasks 1–5); `compas_cgal.toolpath` (`ToolpathResult`, `ToolpathOperation`, `OperationType`, `trochoidal_mat_toolpath_circular`); `_polygon_to_ccw_vertices` conversion pattern (reimplemented via public compas API — do not import the private helper).
- Produces:
  - `stock.py`: `class Stock` — `__init__(polygon: Polygon, holes: list[Polygon] | None = None)`, `contains(x, y) -> bool`, `subtract_capsule(...)`, `subtract_disk(...)`, `subtract_arc_sweep(...)`, property `raw` (the `_stock_2.Stock2`). Thin, typed, Google docstrings.
  - `engagement.py`:
    ```python
    @dataclass
    class OperationEngagement:
        op_index: int
        operation: OperationType
        max_tea: float          # radians
        cap_certified: bool
        stations: int

    @dataclass
    class EngagementReport:
        tool_diameter: float
        tea_cap: float
        operations: list[OperationEngagement]
        max_tea: float
        cap_violations: int
        engaged_ops: int

    def audit_toolpath_engagement(
        polygon: Polygon,
        result: ToolpathResult,
        tool_diameter: float,
        tea_cap: float,
        holes: list[Polygon] | None = None,
    ) -> EngagementReport: ...
    ```
  - Audit semantics (replay): iterate `result.operations` in order. For ops in `{CUT, LEAD_IN, LEAD_OUT, LINK}` at cutting height: certify TEA against the CURRENT stock (Line → `certify_segment_tea`; Arc/Circle → station sampling every `AUDIT_ARC_STEP_FRACTION = 0.05` of circumference through `engagement_at`, growth bound applied between stations in Python mirroring Task 5), then subtract the op's sweep (Line → capsule; Arc/Circle → arc sweep; PLUNGE → disk). RETRACT: skip. LINK ops at clearance height (z of both ends > cut z): skip measurement, no subtraction.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_engagement_audit.py`:

```python
import math

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal.engagement import EngagementReport, audit_toolpath_engagement
from compas_cgal.stock import Stock
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


def test_stock_wrapper_roundtrip():
    stock = Stock(SQUARE)
    assert stock.contains(5.0, 5.0)
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not stock.contains(5.0, 5.0)


def test_audit_square_trochoid_toolpath():
    result = trochoidal_mat_toolpath_circular(SQUARE, tool_diameter=2.0, pitch=1.0)
    report = audit_toolpath_engagement(SQUARE, result, tool_diameter=2.0, tea_cap=math.radians(120))
    assert isinstance(report, EngagementReport)
    assert report.engaged_ops > 0
    assert len(report.operations) == len(result.operations)
    assert report.max_tea > 0.0
    # First engaged motion cuts virgin stock: full slotting must be reported
    # truthfully — the existing generator has no engagement regulation.
    assert report.max_tea > math.radians(150)
    assert report.cap_violations > 0


def test_audit_monotone_stock_depletion():
    """Replaying twice: second audit sees a cleared pocket, lower engagement."""
    result = trochoidal_mat_toolpath_circular(SQUARE, tool_diameter=2.0, pitch=1.0)
    r1 = audit_toolpath_engagement(SQUARE, result, tool_diameter=2.0, tea_cap=math.pi)
    # rebuild stock, subtract everything, audit again on depleted stock
    stock = Stock(SQUARE)
    from compas_cgal.engagement import _subtract_operation  # exposed for reuse
    for op in result.operations:
        _subtract_operation(stock, op, tool_radius=1.0)
    r2_ops = []
    for i, op in enumerate(result.operations):
        pass  # depleted-stock audit exercised through the public API below
    report2 = audit_toolpath_engagement(SQUARE, result, tool_diameter=2.0, tea_cap=math.pi,
                                        holes=None)
    assert report2.max_tea <= r1.max_tea + 1e-9
```

Simplify the third test to what it actually asserts — delete the `r2_ops`
loop and the `_subtract_operation` import block; keep only the two public
audits and the final monotonicity assertion (`report2` equals a fresh audit,
so assert `report2.max_tea == pytest.approx(r1.max_tea)` — determinism, not
depletion). Name it `test_audit_is_deterministic`.

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTHONPATH=src $PY -m pytest tests/test_engagement_audit.py -n auto -q`
Expected: `ModuleNotFoundError: No module named 'compas_cgal.stock'`.

- [ ] **Step 3: Implement `stock.py` and `engagement.py`**

`src/compas_cgal/stock.py` — wrapper converting `Polygon` → CCW float64
vertex arrays (z-planarity check + shoelace orientation exactly as
`toolpath._polygon_to_ccw_vertices`; reimplement locally in ~15 lines, do not
import the private function), forwarding to `_stock_2.Stock2`.

`src/compas_cgal/engagement.py` — the two dataclasses, the audit:

```python
AUDIT_ENGAGED = frozenset({OperationType.CUT, OperationType.LEAD_IN,
                           OperationType.LEAD_OUT, OperationType.LINK})


def _op_is_at_cut_height(op: ToolpathOperation) -> bool:
    g = op.geometry
    if isinstance(g, Line):
        return abs(float(g.start[2]) - float(g.end[2])) <= TOL.absolute
    return True  # arcs/circles are planar cut motions by construction


def _subtract_operation(stock: Stock, op: ToolpathOperation, tool_radius: float) -> None:
    g = op.geometry
    if isinstance(g, Line):
        stock.subtract_capsule(float(g.start[0]), float(g.start[1]),
                               float(g.end[0]), float(g.end[1]), tool_radius)
    elif isinstance(g, Circle):
        c = g.frame.point
        s = g.point_at(0.0)
        stock.subtract_arc_sweep(float(c[0]), float(c[1]),
                                 float(s[0]), float(s[1]), float(s[0]), float(s[1]),
                                 op.clockwise, tool_radius)
    elif isinstance(g, Arc):
        c = g.frame.point
        s = g.point_at(0.0)
        e = g.point_at(1.0)
        stock.subtract_arc_sweep(float(c[0]), float(c[1]),
                                 float(s[0]), float(s[1]), float(e[0]), float(e[1]),
                                 op.clockwise, tool_radius)
```

Measurement for Line ops via `_stock_2.certify_segment_tea`; for Arc/Circle
ops sample stations at arc-length step `AUDIT_ARC_STEP_FRACTION * circumference`
(`AUDIT_ARC_STEP_FRACTION = 0.05` — 20 stations per full circle; growth bound
between stations via the same formula as C++, reimplemented in 4 lines with a
cross-reference comment), `engagement_at` per station, cap flag if any station
exceeds or margins unclosable. Vertical plunge → subtract disk, no TEA.
Retract / clearance-height links → skip both. Assemble `EngagementReport`
(max over ops, violation count, engaged count).

- [ ] **Step 4: Run tests to verify they pass**

Expected: 3 passed (plus prior suites still green:
`PYTHONPATH=src $PY -m pytest tests/test_stock.py tests/test_toolpath.py tests/test_engagement_audit.py -n auto -q` → 65 passed).

- [ ] **Step 5: Ruff + commit**

```bash
$PY -m ruff check --fix src/compas_cgal/stock.py src/compas_cgal/engagement.py tests/test_engagement_audit.py && $PY -m ruff format src/compas_cgal/stock.py src/compas_cgal/engagement.py tests/test_engagement_audit.py
git add src/compas_cgal/stock.py src/compas_cgal/engagement.py tests/test_engagement_audit.py
git commit -m "feat(engagement): typed Stock wrapper + toolpath engagement audit"
```

---

### Task 7: Monte-Carlo differential oracle

**Files:**
- Test: `tests/test_engagement_oracle.py`

**Interfaces:**
- Consumes: `Stock`, `_stock_2.engagement_at`, `audit_toolpath_engagement` internals (public API only).
- Produces: an independent numpy raster oracle used ONLY in tests:
  - grid raster of the stock (bool occupancy, cell size `H = 0.02` on 10-unit pockets), updated by the same subtraction sequence via distance checks;
  - TEA estimate = fraction of `N_ANGLE = 1440` rim samples whose cell is material, × 2π.
  - Agreement tolerance derivation in the test docstring: rim sample step (2π/1440 ≈ 4.4 mrad) + boundary-cell ambiguity (≤ 2·H/r rad per interval endpoint, 2 endpoints per run) + chain under-coverage slack ⇒ `TOL_TEA = 2*math.pi/1440 + 4*H/r + 0.01`.

- [ ] **Step 1: Write the oracle + differential tests (they should PASS immediately — this task is a verification net, watch one deliberately-broken case fail instead)**

```python
import math

import numpy as np
import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

H = 0.02
N_ANGLE = 1440


class RasterStock:
    """Independent occupancy-grid oracle. Same subtraction semantics, no CGAL."""

    def __init__(self, size=10.0, h=H):
        n = int(round(size / h))
        self.h = h
        xs = (np.arange(n) + 0.5) * h
        self.gx, self.gy = np.meshgrid(xs, xs)
        self.material = np.ones_like(self.gx, dtype=bool)

    def subtract_disk(self, cx, cy, r):
        self.material &= (self.gx - cx) ** 2 + (self.gy - cy) ** 2 > r * r

    def subtract_capsule(self, x0, y0, x1, y1, r):
        dx, dy = x1 - x0, y1 - y0
        L2 = max(dx * dx + dy * dy, 1e-30)
        t = np.clip(((self.gx - x0) * dx + (self.gy - y0) * dy) / L2, 0.0, 1.0)
        px, py = x0 + t * dx, y0 + t * dy
        self.material &= (self.gx - px) ** 2 + (self.gy - py) ** 2 > r * r

    def tea_at(self, cx, cy, r):
        ang = np.linspace(0.0, 2.0 * math.pi, N_ANGLE, endpoint=False)
        x = cx + r * np.cos(ang)
        y = cy + r * np.sin(ang)
        i = np.clip((x / self.h).astype(int), 0, self.material.shape[1] - 1)
        j = np.clip((y / self.h).astype(int), 0, self.material.shape[0] - 1)
        return float(self.material[j, i].mean()) * 2.0 * math.pi


SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])


@pytest.mark.parametrize("seed", range(8))
def test_tea_agrees_with_raster_oracle(seed):
    rng = np.random.default_rng(seed)
    exact = Stock(SQUARE)
    raster = RasterStock()
    for _ in range(6):
        x0, y0, x1, y1 = rng.uniform(1.5, 8.5, 4)
        r = float(rng.uniform(0.3, 0.9))
        exact.subtract_capsule(x0, y0, x1, y1, r)
        raster.subtract_capsule(x0, y0, x1, y1, r)
    r_tool = 0.5
    tol = 2 * math.pi / N_ANGLE + 4 * H / r_tool + 0.01
    for _ in range(40):
        cx, cy = rng.uniform(1.0, 9.0, 2)
        total, _, _ = _stock_2.engagement_at(exact.raw, cx, cy, r_tool, math.pi)
        assert total == pytest.approx(raster.tea_at(cx, cy, r_tool), abs=tol)


def test_growth_bound_sound():
    """Empirical soundness of the certificate growth bound (Task 5 lemma)."""
    rng = np.random.default_rng(11)
    exact = Stock(SQUARE)
    for _ in range(5):
        x0, y0, x1, y1 = rng.uniform(1.5, 8.5, 4)
        exact.subtract_capsule(x0, y0, x1, y1, float(rng.uniform(0.4, 1.0)))
    r = 0.5
    for _ in range(200):
        cx, cy = rng.uniform(1.5, 8.5, 2)
        d = float(rng.uniform(1e-4, 0.2))
        theta = float(rng.uniform(0, 2 * math.pi))
        t0, _, _ = _stock_2.engagement_at(exact.raw, cx, cy, r, math.pi)
        t1, _, _ = _stock_2.engagement_at(exact.raw, cx + d * math.cos(theta), cy + d * math.sin(theta), r, math.pi)
        bound = 4 * math.asin(min(1.0, d / (2 * r))) + 2 * math.acos(max(-1.0, 1.0 - d / r))
        assert abs(t1 - t0) <= bound + 1e-9
```

- [ ] **Step 2: Verify the net catches breakage**

Temporarily change `CHAIN_SLACK_FRACTION` to `1e-1` in `stock_2.cpp`, rebuild,
run `tests/test_engagement_oracle.py` — Expected: agreement test FAILS
(under-coverage now visible). Revert the constant, rebuild.

- [ ] **Step 3: Run tests to verify they pass**

Run: `PYTHONPATH=src $PY -m pytest tests/test_engagement_oracle.py -n auto -q`
Expected: 9 passed (8 seeds + growth bound).

- [ ] **Step 4: Commit**

```bash
git add tests/test_engagement_oracle.py
git commit -m "test(engagement): raster differential oracle + growth-bound soundness"
```

---

### Task 8: Baseline engagement audit of the existing generator

**Files:**
- Create: `scripts/engagement_baseline.py`
- Create: `docs/benchmarks/engagement_baseline.md` (generated, committed)
- Test: `tests/test_engagement_audit.py` (one addition)

**Interfaces:**
- Consumes: `audit_toolpath_engagement`, `trochoidal_mat_toolpath_circular`, benchmark polygons (copy the literal coordinates from `tests/test_toolpath.py` into the script — scripts do not import test modules).
- Produces: `python scripts/engagement_baseline.py` writes the report (BLUF first: worst TEA, violation counts per pocket at cap 120°, then the per-pocket table) and a JSON sidecar `docs/benchmarks/engagement_baseline.json` with per-op numbers. This baseline is the number SP2 must beat; the report says so in its first paragraph.

- [ ] **Step 1: Write the failing test**

```python
def test_baseline_script_generates_report(tmp_path):
    import subprocess, sys, os
    env = dict(os.environ, PYTHONPATH="src")
    out = subprocess.run(
        [sys.executable, "scripts/engagement_baseline.py", "--out", str(tmp_path)],
        capture_output=True, text=True, env=env, timeout=600,
    )
    assert out.returncode == 0, out.stderr
    md = (tmp_path / "engagement_baseline.md").read_text()
    assert md.splitlines()[0].startswith("# Engagement baseline")
    assert "worst TEA" in md
```

- [ ] **Step 2: Run it — fails with missing script.**

- [ ] **Step 3: Write `scripts/engagement_baseline.py`**

Pockets: SQUARE, IRREGULAR, L_SHAPE, STAR, KITE, dumbbell(2.4), island pocket
(literal coordinates from `tests/test_toolpath.py`), tool 1.0 (0.5 for the
small shapes, matching the test suite's parameters), `tea_cap = 120°`,
`clearance_z=2.0` everywhere. For each: generate with
`trochoidal_mat_toolpath_circular`, audit, accumulate. Emit markdown
(BLUF paragraph: worst pocket, worst TEA in degrees, total violations,
sentence "This baseline is the reference SP2 must beat.") + JSON. `--out`
argument (default `docs/benchmarks`). ~80 lines.

- [ ] **Step 4: Test passes; generate committed baseline**

```bash
PYTHONPATH=src $PY scripts/engagement_baseline.py
```

- [ ] **Step 5: Commit**

```bash
git add scripts/engagement_baseline.py docs/benchmarks/engagement_baseline.md docs/benchmarks/engagement_baseline.json tests/test_engagement_audit.py
git commit -m "feat(engagement): baseline audit of existing generator (SP2 reference)"
```

---

### Task 9: Wall-clock exit gate + wrap-up

**Files:**
- Test: `tests/test_stock_gate.py`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: everything above.
- Produces: the SP1 exit-gate test and the changelog entry.

- [ ] **Step 1: Write the gate test**

```python
import math
import time

import numpy as np
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

POCKETS = {
    "square": Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]),
    "random100": Polygon(
        [
            (
                (12 + 3 * math.sin(9 * i)) * math.cos(2 * math.pi * i / 100),
                (12 + 3 * math.sin(9 * i)) * math.sin(2 * math.pi * i / 100),
                0.0,
            )
            for i in range(100)
        ]
    ),
}

WALL_CLOCK_GATE_S = 10.0  # spec §8: per-pocket commit+certify throughput gate


def test_sp1_wall_clock_gate():
    rng = np.random.default_rng(0)
    for name, polygon in POCKETS.items():
        stock = Stock(polygon)
        t0 = time.perf_counter()
        for _ in range(60):  # representative pass count for one pocket
            x0, y0, x1, y1 = rng.uniform(-6.0, 6.0, 4) + 5.0
            stock.subtract_capsule(x0, y0, x1, y1, 0.5)
            _stock_2.certify_segment_tea(stock.raw, x0, y0, x1, y1, 0.5, math.radians(120))
        elapsed = time.perf_counter() - t0
        assert elapsed < WALL_CLOCK_GATE_S, f"{name}: {elapsed:.1f}s exceeds SP1 gate"
```

(For `random100`, shift the sample box to the polygon's actual span —
implementer: use the polygon bounding box center instead of the +5.0 offset.)

- [ ] **Step 2: Run the gate**

Run: `PYTHONPATH=src $PY -m pytest tests/test_stock_gate.py -q`
Expected: PASS, or an honest FAIL with the measured number — a failing gate is
a *result*, not a bug: report it and stop; the fallback decision (coarser
exact station grids / GPS locality) is §3 of the design doc and belongs to a
follow-up conversation, not to silent tolerance-raising.

- [ ] **Step 3: Full suite + changelog**

```bash
PYTHONPATH=src $PY -m pytest tests/ -n auto -q
```
Expected: everything green (previous 45 toolpath tests untouched).

Add to `CHANGELOG.md` under `## Unreleased` / `### Added`:

```markdown
- `stock` / `engagement`: exact in-process 2D stock model (CGAL circle-segment boolean set) with certified under-covering sweep chains, exact station TEA queries, motion certificates with a conservative growth bound, raster differential oracle, and an engagement baseline audit of the trochoidal generator (SP1 of the adaptive-clearing programme).
```

- [ ] **Step 4: Commit**

```bash
git add tests/test_stock_gate.py CHANGELOG.md
git commit -m "test(stock): SP1 wall-clock exit gate + changelog"
```

---

## Self-review record

- **Spec coverage:** §3 stock model → Tasks 1–3 (with the conservatism rule
  extended to straight capsules after the irrational-side-line finding —
  Task 2 documents it; the design doc's §3 covers this by its general
  under-covering rule). TEA certificate → Tasks 4–5. Audit → Tasks 6, 8.
  Monte-Carlo oracle → Task 7. Wall-clock gate → Task 9. Residual extraction
  (`residual_polylines`) is NOT needed by any SP1 gate and is deferred to
  SP2's report task — recorded here as a conscious YAGNI cut from the spec's
  component table.
- **Placeholder scan:** one intentional instruction ("Remove the placeholder
  line") in Task 2 marks dead scaffolding for deletion; no TBDs remain.
- **Type consistency:** `Stock2` / `Stock.raw` / `engagement_at` /
  `certify_segment_tea` signatures identical across Tasks 1–9;
  `EngagementReport` fields consistent between Tasks 6 and 8.
