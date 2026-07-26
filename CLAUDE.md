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

## Exact-Kernel Discipline (Epeck / circle-segment traits) — MANDATORY

Full binding doctrine: `docs/exactness.md` (MkDocs page, repo-anchored,
14-question review checklist). `docs/dev/exact-cgal-idioms.md` is a generic
primer. The rules below are the enforced distillation.

The stock/engagement modules run on `Exact_predicates_exact_constructions_kernel`
with `Gps_circle_segment_traits_2` (one-root point coordinates,
`CGAL::Sqrt_extension`). CGAL's own doctrine (Developer Manual, Robustness):
"imprecise calculations can cause wrong or, much worse, mutually contradictory
decisions." In an exact-constructions kernel there is NO legitimate epsilon,
deflation factor, ulp-nudge, or `nextafter` in any decision path. C-style
floating-point precision handling is a defect, always.

### The deciding / reporting split
- A comparison that changes control flow, output topology, or a certificate is
  a DECISION: it must be an exact predicate on exact quantities.
- A number produced for humans (angles, lengths, areas in reports and
  statistics) is REPORTING: `to_double` is fine there, and reported values
  must never feed back into decisions.

### Boundary doctrine for doubles
- Doubles enter exact-land ONCE, at declared API boundaries, by exact
  injection — every `double` IS a rational; there is no snapping and no
  tolerance at the seam.
- Transcendental user intent (an angle cap, a feed angle) is NOT enforced
  transcendentally. The API contract is an exact rational surrogate (e.g. a
  squared-chord threshold `4·sin²(cap/2)` computed by the caller as a double
  and injected exactly). Certificates are exact statements about the exact
  surrogate; the sub-ulp gap between the surrogate and the transcendental
  intent is an API semantic, documented at the parameter — never an in-core
  correction constant.

### One-root numbers (`Sqrt_extension` / CoordNT)
- Arithmetic (`+ − × /`) requires operands in the SAME extension — the
  documented precondition is `a.root()==0 or b.root()==0 or a.root()==b.root()`;
  cross-root arithmetic is undefined behavior. Never write it.
- Cross-root COMPARISON of arrangement/traits points is exact and supported —
  use `==`, `compare_x`, `compare_xy` on points instead of double-gap checks.
- Abutting sub-curves share endpoints EXACTLY (same arrangement vertex).
  Reconstructing adjacency with double-angle tolerances is a design error:
  merge by exact endpoint equality.
- When a needed exact comparison is not a stock predicate (e.g. sign of
  `A + B√α + C√β + D√(αβ)` with rational A,B,C,D — squared distances and
  orientations of points from different circles reduce to this), implement it
  exactly over the rational coefficients (`a0()/a1()/root()` + repeated
  squaring). Do not fall back to `to_double`.

### Analytic bounds are not precision handling
- A mathematical lemma bound (Lipschitz/growth bounds over transcendental
  functions) may be evaluated in doubles ONLY to drive refinement, and only
  when the lemma's constants carry proof-level slack (integer-factor margins
  stated in the derivation comment) that dwarfs floating-point evaluation by
  orders of magnitude. The safe failure direction (more refinement / flag as
  violation) must be stated. This clause never applies to geometric truth —
  point membership, adjacency, thresholds — which is predicate territory.

### Numeric comparison is the exact-kernel idiom
- Comparisons are made at the NUMBER-TYPE level with CGAL's comparison
  machinery: `CGAL::sign`, `CGAL::compare`, `Comparison_result`
  (SMALLER/EQUAL/LARGER) — never `operator<` on `to_double` values, never
  re-deriving what `Real_embeddable_traits` already provides.
- `Sqrt_extension` is RealEmbeddable: `CGAL::sign` and same-root
  `CGAL::compare` (plus same-root `+ − ×`) are exact and supported — build
  derived quantities INSIDE one extension and compare there. A mixed-radical
  sign (`u + √β·w` with u, w ∈ ℚ(√α)) decomposes into `sign(u)`, `sign(w)`,
  and `compare(u², β·w²)` — three supported exact calls, not a hand-rolled
  bignum routine.

### Kernel primitives first
- Per the Developer Manual: "use kernel primitives whenever possible."
  Reach for `CGAL::orientation`, `compare_squared_distance`,
  `has_smaller_distance_to_point`, traits functors — before writing any
  coordinate arithmetic of your own.

### Exact algebraic code
- Read `docs/exactness.md`, especially **Algebraic kernels: identify roots
  without minimal polynomials**, before changing event partitions or
  algebraic certificates.
- Root identity is primitive square-free polynomial + ordered-real-root
  ordinal. Coalesce equal projections with
  `Polynomial_traits_d::Gcd_up_to_constant_factor` and `Compare_1`; NEVER
  trial-factor bignum coefficients to manufacture a minimal polynomial.
- Use `Algebraic_structure_traits`, `Fraction_traits`,
  `Polynomial_traits_d`, and algebraic-kernel functors before backend
  representation arithmetic. Explicit backend construction is allowed at
  the coefficient boundary; representation access is allowed for canonical
  bytes/build attestation.
- Zero-set canonicalization does not imply sign preservation. Signed
  predicates first detect a vanishing source factor, then retain the original
  unit/content sign and factor multiplicity parity.

## Build & Test

```bash
pixi run pytest -- tests/ -n auto
pixi run lint
```

## Stack
- C++20, CGAL 6.0.1 (Epick kernel), nanobind, Eigen, scikit-build-core
- Python 3.9+ (compas dependency), pixi for package management
- ruff for Python formatting
