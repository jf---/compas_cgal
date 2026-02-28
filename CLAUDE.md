# compas_cgal

C++ CGAL bindings for compas, exposed to Python via nanobind + Eigen.

## CGAL Idiomatic Code — MANDATORY

This project uses `Exact_predicates_inexact_constructions_kernel` (Epick). All C++ code MUST follow these rules:

### Exact predicates: NEVER use epsilon for geometric decisions
- Use `==` / `!=` on `Point_2` for coincidence (exact in Epick)
- Use `CGAL::orientation()` for left/right/collinear — NEVER manual cross products with epsilon
- Use `CGAL::compare_xy()` for deterministic ordering — NEVER nested epsilon comparison
- Use `CGAL::squared_distance()` for distance comparisons — compare squared values, avoid sqrt

### Kernel types: NEVER decompose into raw doubles for geometry
- Use `p1 - p0` for `Vector_2`, NEVER `Vector_2(p1.x() - p0.x(), p1.y() - p0.y())`
- Use `CGAL::barycenter(p0, w0, p1, w1)` for interpolation, NEVER manual coordinate lerp
- Use `Vector_2::perpendicular(CGAL::CLOCKWISE)` for normals, NEVER manual rotation
- Use `Segment_2`, `Circle_2`, `Direction_2` as value types — NEVER store geometry as scalar tuples
- Track positions as `Point_2` + separate z double — NEVER `double x, y, z` triples

### Construction guards: epsilons ARE allowed for inherently inexact operations
- Tangent normalization (`length > 1e-12` before dividing) — OK
- Degenerate edge filtering (`length < 1e-12`) — OK
- These guard against division-by-zero in double arithmetic, not geometric predicates

### Utilities
- Use `approx_distance(Point_2, Point_2)` and `approx_length(Vector_2)` utilities from toolpath.cpp
- Use `std::numbers::pi` (C++20), NEVER hand-rolled PI constants

### Anti-patterns to REJECT in review
- `CGAL::to_double(x) < 1e-12` for geometric decisions (use exact predicate)
- `double sx = CGAL::to_double(p.x()); double sy = CGAL::to_double(p.y());` when Point_2 arithmetic works
- Storing geometry as `(sx, sy, ex, ey, cx, cy, radius)` raw doubles — use kernel types
- `std::sqrt(CGAL::to_double(CGAL::squared_distance(a, b)))` — use `approx_distance`

## Build & Test

```bash
pip install --no-build-isolation -ve .
pytest tests/ -v
```

## Stack
- C++20, CGAL 6.0.1 (Epick kernel), nanobind, Eigen, scikit-build-core
- Python 3.9+ (compas dependency), pixi for package management
- ruff for Python formatting
