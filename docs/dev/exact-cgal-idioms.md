# Idiomatic Exact CGAL

> Canonical exact-kernel doctrine for this repository. The repo `CLAUDE.md`
> carries the enforced distillation; this file is the full reference.
> The deepest idiom: **exactness is a property of the complete dataflow, not
> of the kernel typedef.**

## The central rule

Keep geometric decisions inside the kernel. Keep approximation at the system
boundary.

Exact CGAL code is not merely code compiled with EPECK. It is code whose
control flow, topology, and combinatorics never depend on values that have
been prematurely converted to floating point.

CGAL itself distinguishes:

* Predicates: decisions such as orientation, ordering, sidedness, incidence
  and intersection existence.
* Constructions: producing new coordinates or geometric objects.

Predicates return discrete results; constructions produce numerical data.
This distinction is the foundation of CGAL's robustness model.

## 1. Choose the weakest kernel that is actually sufficient

```cpp
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
```

Use Epick when:

* only predicates must be reliable;
* constructed coordinates are outputs, not subsequent decision inputs;
* the package explicitly supports it.

```cpp
using K = CGAL::Exact_predicates_exact_constructions_kernel;
```

Use Epeck when constructed geometry is fed into later predicates,
intersections, arrangements, Boolean operations, offsets, repeated clipping,
or topology-building.

Epick gives exact predicates but generally inexact constructions. Epeck gives
both exact predicates and exact constructions, using filtering and lazy exact
evaluation to avoid paying the full exact-arithmetic cost in ordinary cases.

A useful test is:

Can an approximate constructed point later change an orientation, ordering,
equality, incidence, or classification decision?

If yes, use exact constructions through that pipeline.

## 2. Kernel exactness does not repair inexact input semantics

This is legal:

```cpp
K::Point_2 p(0.1, 0.2);
```

With Epeck, CGAL represents the supplied double values exactly. But it
represents the exact binary floating-point values passed into the
constructor — not necessarily the mathematical decimals 1/10 and 1/5.

That distinction matters.

For measured data, this is often exactly what you want: preserve the
observations as received.

For symbolic or decimal-defined geometry, construct the intended rational
values explicitly:

```cpp
using K  = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
const FT one_tenth = FT(1) / FT(10);
const FT one_fifth = FT(1) / FT(5);
const K::Point_2 p(one_tenth, one_fifth);
```

Do not write:

```cpp
const FT x = 0.1;  // exact representation of the binary double 0.1
```

when your model says "exactly one tenth."

## 3. Do not leave the kernel to make a decision

Bad:

```cpp
double d = CGAL::to_double(CGAL::squared_distance(p, q));
if (d < tolerance * tolerance) {
    merge_vertices();
}
```

Now topology depends on an approximation.

Better:

```cpp
const K::FT tolerance_sq = tolerance * tolerance;
if (CGAL::squared_distance(p, q) < tolerance_sq) {
    merge_vertices();
}
```

Better still, when a kernel predicate expresses the actual question:

```cpp
if (CGAL::compare_squared_distance(p, q, tolerance_sq)
    == CGAL::SMALLER) {
    merge_vertices();
}
```

Use `CGAL::to_double()` only for:

* rendering;
* logging;
* user-interface display;
* export to an explicitly approximate format;
* constructing a conservative external bound.

CGAL documents `to_double()` as an approximation with no general accuracy
guarantee. For lazy exact numbers, repeated calls may even improve after an
intervening exact evaluation. `to_interval()` is preferable when you need a
certified enclosing double interval.

## 4. Prefer predicates over coordinate arithmetic

Bad:

```cpp
const auto ux = b.x() - a.x();
const auto uy = b.y() - a.y();
const auto vx = c.x() - a.x();
const auto vy = c.y() - a.y();
if (ux * vy - uy * vx > 0) {
    // ...
}
```

Mathematically this can be exact with Epeck::FT, but it is less idiomatic and
less generic.

Prefer:

```cpp
switch (CGAL::orientation(a, b, c)) {
case CGAL::LEFT_TURN:
    break;
case CGAL::RIGHT_TURN:
    break;
case CGAL::COLLINEAR:
    break;
}
```

Why?

* It states the geometric intent.
* It permits CGAL's filtered predicate machinery.
* It handles representation details.
* It generalizes better across kernels and traits.
* It preserves degenerate cases explicitly.

Likewise prefer `CGAL::compare_x`, `CGAL::compare_y`, `CGAL::collinear`,
`CGAL::do_intersect`, `CGAL::bounded_side_2`,
`CGAL::compare_distance_to_point` over manually reproducing their arithmetic.

## 5. Avoid square roots unless the problem genuinely requires them

Prefer squared quantities:

```cpp
const K::FT d2 = CGAL::squared_distance(p, q);
```

rather than:

```cpp
const auto d = CGAL::sqrt(CGAL::squared_distance(p, q));
```

Ordinary Epeck is based primarily on exact rational arithmetic. A generic
square root of a rational is usually irrational and cannot be represented by
a rational field type.

If the algorithm truly needs exact values involving square roots, use a
kernel or traits class whose number type models FieldWithSqrt, such as:

```cpp
using K =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
```

CGAL provides separate predefined kernels for square roots, kth roots, and
more general algebraic roots because "exact" does not mean "closed under
every mathematical operation."

## 6. Treat degeneracy as an ordinary input case

Exact arithmetic exposes degeneracies reliably; it does not make them
disappear.

Idiomatic code handles all documented outcomes:

```cpp
const auto result = CGAL::intersection(segment, line);
if (!result) {
    // Empty intersection.
} else if (const auto* point =
               std::get_if<K::Point_2>(&*result)) {
    // Point intersection.
} else if (const auto* overlap =
               std::get_if<K::Segment_2>(&*result)) {
    // Collinear overlap.
}
```

Do not assume that two segments intersect only at a point. Modern CGAL
intersection APIs return an optional variant precisely because several
geometrically valid result types may occur.

The same applies to: collinear points; coplanar faces; repeated vertices;
zero-length segments; tangent intersections; coincident curves; non-unique
nearest points; lower-dimensional Boolean results.

## 7. Never use epsilon to imitate a predicate

Suspicious:

```cpp
if (std::abs(CGAL::to_double(value)) < 1e-12) {
    // Treat as zero.
}
```

There are two separate questions:

1. Is the mathematical value zero?
2. Is the value close enough for an application-specific tolerance model?

For the first:

```cpp
if (CGAL::is_zero(value)) {
    // Exactly zero.
}
```

For the second, compare exactly against an explicitly modelled tolerance:

```cpp
if (CGAL::abs(value) <= tolerance) {
    // Within the application's tolerance.
}
```

A tolerance is a domain decision — machining tolerance, measurement
uncertainty, resolution — not a numerical patch for unreliable arithmetic.

## 8. Do not mix kernels casually

This is a common source of accidental loss of exactness:

```cpp
Epeck::Point_3 exact_point = ...;
Epick::Point_3 approximate_point(
    CGAL::to_double(exact_point.x()),
    CGAL::to_double(exact_point.y()),
    CGAL::to_double(exact_point.z()));
```

That conversion is a deliberate approximation boundary. It should look like
one and preferably live in a named adapter:

```cpp
Render_point to_render_point(const Exact_point& p);
```

Avoid interleaving exact and approximate types throughout the algorithm.

A useful architecture is:

```
parse/model input
      ↓
exact geometric domain
      ↓
topology and classification
      ↓
explicit approximation adapter
      ↓
rendering / visualization / machine output
```

## 9. Parameterize algorithms by kernel or traits

Library-quality CGAL code generally should not hard-code Epeck in every
function:

```cpp
template<class Kernel>
typename Kernel::Orientation
classify_turn(
    const typename Kernel::Point_2& a,
    const typename Kernel::Point_2& b,
    const typename Kernel::Point_2& c)
{
    return CGAL::orientation(a, b, c);
}
```

Or infer the kernel from the object. For package development, specify the
actual concept requirements: Kernel, FieldNumberType, FieldWithSqrt, exact
predicates only, exact constructions, Cartesian versus homogeneous
compatibility, package-specific traits concepts.

CGAL's developer manual explicitly says packages should state their
arithmetic/kernel requirements and, where appropriate, work with both
Cartesian and homogeneous representations rather than depending on
representation internals.

## 10. Use kernel types, not assumptions about their implementation

Prefer:

```cpp
using FT       = typename K::FT;
using Point_3  = typename K::Point_3;
using Vector_3 = typename K::Vector_3;
```

Do not assume:

```cpp
using FT = CGAL::Gmpq;
```

Epeck::FT is commonly a lazy exact wrapper around an exact rational backend.
Its exact implementation can depend on the CGAL build and available
libraries. The abstraction is intentional.

Also avoid reaching into lazy-number internals with `.exact()` unless
profiling has demonstrated a specific need. Forcing exact evaluation
systematically defeats laziness.

## 11. Respect the package's traits contract

The correct kernel is package-dependent.

Examples:

* Convex hull over original input points generally needs exact predicates,
  not exact constructions.
* Arrangement insertion of intersecting segments creates new vertices and
  generally requires exact constructions.
* A triangulation may work with Epick, but exact dual/Voronoi constructions
  may require Epeck.
* Circular arcs and conics may require algebraic kernels rather than merely
  rational EPECK.

CGAL's package documentation should be considered part of the type contract,
not optional commentary.

## 12. Separate exact topology from approximate metric policy

This is particularly important for CAD/CAM.

An exact kernel can answer:

* Does this intersection exist?
* Which side is the point on?
* Are these entities exactly incident?
* What is the exact ordering?
* Is the Boolean result topologically valid under the mathematical model?

It cannot decide:

* Whether two measured surfaces should be treated as coincident.
* Whether a 3 μm gap is acceptable.
* Whether a sliver is manufacturable.
* Whether two nominally equal CAD edges should be sewn.

Those require an explicit uncertainty or tolerance policy.

A strong architecture therefore distinguishes:

```cpp
struct Geometric_model {
    Exact_geometry geometry;
};
struct Manufacturing_policy {
    Exact_FT linear_tolerance;
    Exact_FT angular_tolerance;
    Exact_FT minimum_feature_size;
};
```

Exact arithmetic makes that policy deterministic. It does not replace it.

## A compact review checklist

Proper exact-kernel code should pass these questions:

1. Is the kernel appropriate for every construction whose result is reused?
2. Does any `to_double()` result affect control flow or topology?
3. Are decimal inputs represented according to their actual semantics?
4. Are kernel predicates used instead of hand-coded determinants and
   epsilons?
5. Are squared quantities used instead of unnecessary roots?
6. Does the selected number type support every required algebraic operation?
7. Are all degeneracies and variant intersection results handled?
8. Are tolerance decisions explicit domain policy rather than numerical
   repair?
9. Are exact-to-approximate conversions isolated at named boundaries?
10. Does the code use kernel and traits abstractions rather than backend
    types?
11. Does it satisfy the specific CGAL package's traits requirements?
12. Has it been tested on degenerate and nearly degenerate data?

You can use Epeck and still write fundamentally inexact code by converting
one coordinate to double before a branch. Conversely, carefully written
Epick code can be fully robust when no nontrivial construction feeds a later
predicate.
