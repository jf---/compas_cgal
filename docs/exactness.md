# Exact-Kernel Discipline

**Exactness is a property of the complete dataflow, not of the kernel
typedef.** You can instantiate `Exact_predicates_exact_constructions_kernel`
and still write fundamentally inexact code by converting one coordinate to
`double` before a branch; conversely, carefully written Epick code is fully
robust when no nontrivial construction ever feeds a later predicate. This
page is the discipline this repository holds itself to: geometric decisions
stay inside the kernel, approximation lives only at named system boundaries,
and every rule below is anchored to code that exists in this codebase.

!!! warning "Why this page exists"

    During SP1 of the adaptive-clearing programme, a generated design placed
    a `CAP_BOUND_DEFLATION = 1.0 - 1e-12` fudge factor and a `1e-12`
    double-angle merge gap inside the engagement kernel's decision core —
    C-style floating-point precision handling inside an exact-constructions
    pipeline. The exact formulation existed all along (see the
    [case study](#case-study-the-deflation-constant-that-wasnt-needed)).
    Nothing in an exact kernel ever needs an epsilon to make a geometric
    decision; if you find yourself typing one, the design is wrong, not the
    arithmetic.

## The central rule

CGAL's robustness model rests on one distinction:

- **Predicates** decide: orientation, ordering, sidedness, incidence,
  intersection existence. They return discrete results
  (`Orientation`, `Comparison_result`, `Sign`, `bool`).
- **Constructions** produce: new coordinates, new geometric objects.

Exact code keeps control flow, topology, and combinatorics downstream of
predicates only. The Developer Manual's phrasing is blunt: imprecise
calculations cause "wrong or, much worse, mutually contradictory decisions."

This repository adds a project-level restatement, the
**deciding / reporting split**:

- A comparison that changes control flow, output topology, or a certificate
  is a *decision* — exact predicate territory, no exceptions.
- A number produced for humans — angles, lengths, areas in reports and
  statistics — is *reporting* — `CGAL::to_double` is fine there, and a
  reported value must never feed back into a decision.

```cpp
// DECIDING — exact, in-kernel (engagement kernel):
if (CGAL::compare_squared_distance(p, q, tolerance_sq) == CGAL::SMALLER) { ... }

// REPORTING — approximate, at the boundary (EngagementReport statistics):
const double tea = std::atan2(CGAL::to_double(sy), CGAL::to_double(sx));
```

## Choosing a kernel

```mermaid
flowchart TD
    A[New module] --> B{Do constructed coordinates\nfeed later predicates?}
    B -- no --> C["Epick\nExact_predicates_inexact_constructions_kernel\n(toolpath.cpp: trochoid generator)"]
    B -- yes --> D{Closed under the\noperations you need?}
    D -- "rational ops only" --> E["Epeck\nExact_predicates_exact_constructions_kernel\n(stock_2.cpp: boolean stock model)"]
    D -- "square roots required" --> F["Epeck_with_sqrt /\nalgebraic kernels"]
    D -- "circles & segments" --> G["Epeck + Arr_circle_segment_traits_2\n(one-root coordinates)\n(engagement_2.cpp: TEA kernel)"]
```

The test is always the same question: **can an approximate constructed point
later change an orientation, ordering, equality, incidence, or
classification decision?** If yes, exact constructions must carry the whole
pipeline to that decision.

Both regimes coexist in this repository by design:

| Module | Kernel | Why it is sufficient |
| --- | --- | --- |
| `toolpath.cpp` (trochoid generator) | Epick | decisions are predicates on input-derived values; constructed tangent points are *outputs* (polylines, arcs) |
| `stock_2.cpp` (in-process stock) | Epeck + `Gps_circle_segment_traits_2` | boolean results are re-queried by later predicates — constructions feed decisions |
| `engagement_2.cpp` (TEA kernel) | same as stock | engagement intervals have one-root endpoints; the cap certificate is a predicate on constructed points |

!!! note "The package's traits contract is part of the type contract"

    `Gps_circle_segment_traits_2` requires *rational* circles (rational
    center, rational squared radius) and represents intersection points as
    one-root numbers. That contract is why tool-offset sweeps of arcs are
    not exactly representable here (their offset radii have irrational
    squares) — and why the stock model uses certified under-covering disk
    chains instead of pretending. Read the package documentation as
    normative, not as commentary.

## Input semantics: doubles are exact, decimals are not

`Epeck::FT(0.1)` represents the *binary double* `0.1` exactly — not the
mathematical decimal 1/10. Both facts are useful; confuse them and you have
a semantic bug no kernel can fix.

- **Measured or computed data** (everything crossing the nanobind boundary
  in this repo): exact injection of the double is precisely correct — the
  observation is preserved bit-for-bit, and there is *no snapping and no
  tolerance at the seam*.
- **Symbolic decimal geometry**: construct the intended rational explicitly.

```cpp
// Measured/computed input (our nanobind boundary) — exact injection:
Stock2::Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary, ...)
// every double IS a rational; Epeck::FT(vertices(i, 0)) is exact.

// Symbolic decimal — say what you mean:
const FT one_tenth = FT(1) / FT(10);   // the mathematical 1/10
const FT not_this  = FT(0.1);          // the binary double 0.1
```

## The boundary doctrine for transcendental intent

Some user intent is inherently transcendental — an engagement-angle cap, for
instance. The discipline: **do not chase the transcendental with epsilons;
change the contract to an exact surrogate.**

The TEA kernel's cap crosses the API as a *chord ratio*:

```python
# Python boundary (compas_cgal.engagement): the ONLY conversion, documented.
cap_chord_ratio = 4.0 * math.sin(tea_cap / 2.0) ** 2   # dimensionless, (0, 4]
```

```cpp
// C++ core: T is an exact rational threshold. The certificate is an exact
// statement about T. No deflation factor exists anywhere downstream.
const Epeck::FT T = Epeck::FT(cap_chord_ratio) * FT(r) * FT(r);
```

The sub-ulp gap between the rational surrogate and the transcendental angle
the user typed is an *API semantic documented at the parameter* — not a
correction constant inside the decision core. After the boundary: exact
only.

```mermaid
flowchart LR
    A[parse / model input\ndoubles, Polygons] --> B[exact geometric domain\nEpeck, one-root points]
    B --> C[topology & classification\npredicates, certificates]
    C --> D[explicit approximation adapter\nto_double, atan2 — named, one-way]
    D --> E[rendering / reports /\nG-code / visualization]
```

Exact-to-approximate conversions are isolated in named adapters at the right
edge of that diagram. Interleaving `to_double` through the middle is how
Epeck code becomes inexact code with extra steps.

## Prefer predicates over coordinate arithmetic

Hand-rolled determinants can be exact over `Epeck::FT` — and are still the
wrong idiom. `CGAL::orientation(a, b, c)` states intent, engages the
filtered-predicate machinery, survives kernel changes, and forces the
degenerate case (`COLLINEAR`) into view instead of letting `> 0` silently
absorb it.

```cpp
switch (CGAL::orientation(a, b, c)) {
    case CGAL::LEFT_TURN:  ...; break;
    case CGAL::RIGHT_TURN: ...; break;
    case CGAL::COLLINEAR:  ...; break;   // degeneracy is an ordinary case
}
```

Reach for `CGAL::compare_x`, `compare_xy`, `compare_squared_distance`,
`has_smaller_distance_to_point`, `do_intersect`, `bounded_side_2` before
writing any coordinate arithmetic of your own. The repo's Epick module
already lives by this (see `toolpath.cpp`'s exact-predicate winding and
repositioning gates); the exact modules have no excuse at all.

## Numeric comparison is the exact-kernel idiom

Comparisons happen at the *number-type* level with CGAL's machinery:
`CGAL::sign`, `CGAL::compare`, `Comparison_result` — never `operator<` on
`to_double` values, never re-deriving what `Real_embeddable_traits`
provides.

Two questions that look similar and are not:

```cpp
if (CGAL::is_zero(value))          { ... }  // is the mathematical value zero?
if (CGAL::abs(value) <= tolerance) { ... }  // within an explicit DOMAIN tolerance?
```

A tolerance is a *domain decision* — machining tolerance, measurement
uncertainty, feature resolution — modelled as an exact `FT` and compared
exactly. It is never a patch for arithmetic you didn't trust.

### One-root numbers (`Sqrt_extension` / CoordNT)

The circle-segment traits represents point coordinates as
`a0 + a1·√root` (`CGAL::Sqrt_extension`). The closure rules are documented
preconditions, not folklore:

- **Arithmetic** (`+ − × /`) requires operands in the *same* extension —
  the precondition is `a.root()==0 or b.root()==0 or a.root()==b.root()`.
  Cross-root arithmetic is undefined behavior. Never write it.
- **Comparison** across roots *is* exact and supported for the traits'
  points (`==`, `compare_x`, `compare_xy`) — the sweep line depends on it.
  Adjacency of sub-curves is therefore decided by **exact endpoint
  equality** (abutting sub-arcs share the same arrangement vertex), never
  by double-angle gap thresholds.
- When a needed comparison is not a stock predicate, *decompose it into
  supported exact calls* rather than hand-rolling bignum arithmetic. The
  TEA kernel's mixed-radical sign — `sign(A + B√α + C√β + D√(αβ))` with
  rational coefficients, which both the chord-vs-threshold predicate and
  the >π orientation test reduce to — decomposes as:

```cpp
// u = A + B√α and w = C + D√α are SAME-root values: legal arithmetic.
// sign(u + √β·w) from three supported exact calls:
const CGAL::Sign su = CGAL::sign(u);
const CGAL::Sign sw = CGAL::sign(w);
if (sw == CGAL::ZERO) return su;
if (su == CGAL::ZERO) return sw;
if (su == sw)         return su;
// opposite signs: the dominant magnitude wins
switch (CGAL::compare(u * u, w * w * CoordNT(beta))) {   // same-root products
    case CGAL::LARGER:  return su;
    case CGAL::SMALLER: return sw;
    case CGAL::EQUAL:   return CGAL::ZERO;
}
```

## Algebraic kernels: identify roots without minimal polynomials

CGAL's `Algebraic_kernel_d_1` model represents a real root with a square-free
polynomial plus an isolating interval. It deliberately does not require a
minimal polynomial: computing one is usually expensive and does not improve
the root operations the kernel needs.

The repository therefore uses this canonical identity for a real event root:

1. denominator-clear and primitive-normalize every source polynomial;
2. square-free-factorize it with the algebraic-kernel or polynomial-traits
   functor;
3. solve the factors and group equal roots with `Compare_1`;
4. for one equality group, fold
   `Polynomial_traits_d::Gcd_up_to_constant_factor` over its source factors;
5. primitive-normalize that nonconstant square-free GCD with positive leading
   coefficient;
6. solve the GCD and locate the equal root with `Compare_1`; its ordinal is
   among all ordered real roots of that GCD.

For a root seen through two different projections, the GCD removes
irrelevant factors shared by neither source. For a root supported by only one
reducible square-free projection, the complete projection remains the
representative. This is intentional: trial-dividing arbitrary-precision
coefficients in search of a minimal polynomial is neither part of the
identity contract nor an acceptable exact-arithmetic algorithm.

Multiplicity is evidence about a source projection, not the number of times
the same source evidence was submitted. The aggregate scalar therefore takes
the maximum multiplicity while the event incidences preserve the distinct
source projections.

An isolating interval is a witness for locating and verifying the root, not
the root's identity. Refinement may change the interval without changing the
ID, but every accepted interval must still isolate the represented root.

## Use CGAL concepts before number backends

Exact code should be written against the algebraic concepts that state the
operation:

- `CGAL::is_square`, dispatched through
  `Algebraic_structure_traits<NT>::Is_square`, for exact square tests;
- `Fraction_traits<NT>::Decompose` and `Compose` for rational
  numerator/denominator access, including the recursive polynomial
  specialization;
- `Polynomial_traits_d<P>::Canonicalize`,
  `Gcd_up_to_constant_factor`, and the square-free operations for polynomial
  normalization;
- algebraic-kernel `Solve_1`, `Compare_1`, `Sign_at_1`/`Sign_at_2`, and
  `Bound_between_1` for roots, signs, ordering, and separating witnesses.

In the pinned CGAL 6.0.1 algebraic kernel,
`Algebraic_real_1::is_rational()` means that a rational representation is
already known; it is not a rationality decision procedure. `Solve_1` can
therefore return an implicit quadratic representation even when that
quadratic has rational roots. For the bounded degree-two convenience path,
the implementation derives the finite candidates with `CGAL::is_square` and
accepts one only after exact `Compare_1` equality. The algebraic root remains
the certified value; the rational is only a proven convenience rendering.

Backend-specific representation access, handwritten GCD/LCM, manual Newton
square roots, backend bit scans, and rational-root trial division bypass
those contracts. Explicit backend construction is allowed at the declared
coefficient boundary, and representation access is allowed for canonical
byte serialization and build attestation.

!!! warning "Materialize multiprecision expressions before locals die"

    CORE and Boost multiprecision arithmetic uses lazy expression templates.
    An `auto`-deduced function or lambda return can therefore retain references
    to local operands instead of returning a number. The expression later
    dereferences dead storage, often crashing only in optimized arithmetic.
    Give such returns the concrete exact number type (`CORE::BigRat`,
    `CORE::Expr`, or the kernel `FT`) so materialization happens before the
    operands leave scope.

!!! warning "Zero-set normalization is not sign normalization"

    Multiplying by a negative unit preserves a polynomial's zero set while
    changing its sign. Replacing an even power \(q^{2k}\) by the square-free
    factor \(q\) also preserves the zero set but changes the sign across
    \(q=0\). Use the algebraic kernel's documented sign functor on the
    original polynomial. If an implementation pre-factors a polynomial,
    first return zero when any source factor vanishes; otherwise multiply the
    original unit/content sign by factor signs according to multiplicity
    parity. Never reuse a zero-set-only canonical form as a signed predicate.

## Avoid square roots; prefer squared quantities

Epeck's field type is rational: a generic square root leaves the field.
Compare `squared_distance` against squared thresholds; keep lengths squared
until the reporting adapter. If the algorithm genuinely needs `√` *inside
decisions*, that is a kernel-selection fact
(`Exact_predicates_exact_constructions_kernel_with_sqrt`, algebraic
kernels), not a license to `std::sqrt(to_double(...))`. "Exact" is not
"closed under every operation" — CGAL ships distinct kernels for sqrt, kth
roots, and algebraic numbers precisely because of this.

In this repository the rule has teeth in both directions:
`approx_length`/`approx_distance` in the Epick module are *constructions*
feeding outputs (allowed); the Epeck stock model never takes a square root
at all — the trochoid guide radius appears only as `guide_r²` (rational) or
as reporting doubles.

## Degeneracy is an ordinary input case

Exact arithmetic makes degeneracies *reliable*, not absent. Idiomatic code
enumerates every documented outcome — CGAL's intersection APIs return
optional variants because a segment-segment intersection *is* sometimes a
segment:

```cpp
const auto result = CGAL::intersection(segment, line);
if (!result) { /* empty */ }
else if (const auto* pt = std::get_if<K::Point_2>(&*result))   { /* point */ }
else if (const auto* seg = std::get_if<K::Segment_2>(&*result)) { /* overlap */ }
```

The engagement kernel's inventory of ordinary-not-exceptional cases: an
empty intersection region (rim fully clear), a crossing-free rim fully in
material (2π run), tangent touches (zero-measure arcs), the exact-π run
(collinear center/start/end), coincident stations, zero-radius guides.
Every one has a branch and a test.

!!! tip "Test on degenerate and nearly degenerate data"

    The property suites in `tests/test_stock.py` exist because exactness
    claims are cheap and degenerate inputs are not. A hypothesis strategy
    plus a handful of constructed degeneracies (bit-identical endpoints,
    antipodal runs, tangent circles) is the difference between "uses Epeck"
    and "is exact."

## Analytic bounds are not precision handling

One clause separates legitimate double arithmetic from smuggled epsilons. A
*mathematical lemma bound* — e.g. the TEA growth bound that closes a
certificate between exact stations — may be evaluated in doubles **only**
when all three hold:

1. its constants carry proof-level slack (integer safety factors stated in
   the derivation comment) that dwarfs floating-point evaluation error by
   orders of magnitude;
2. the safe failure direction is stated (more refinement / flag violation —
   never a false certificate);
3. it drives *refinement*, never geometric truth — membership, adjacency,
   thresholds remain predicate territory.

In the certifier this shows up as *threshold selection*: the guard shrinks
the cap that the exact station predicate then tests. The lemma picks which
exact question to ask; it never answers one.

## Exact topology, approximate metric policy

For CAD/CAM the last separation is the one that keeps the system honest. An
exact kernel answers: does this intersection exist, which side, exactly
incident, exact ordering, topologically valid boolean. It cannot answer:
should these measured surfaces be treated as coincident, is a 3 µm gap
acceptable, is this sliver manufacturable. Those belong to an explicit
policy object with exact-`FT` fields:

```cpp
struct Manufacturing_policy {          // SP2 direction of travel
    Epeck::FT linear_tolerance;        // domain decisions, exactly modelled
    Epeck::FT angular_surrogate;       // e.g. a chord ratio, not an angle
    Epeck::FT minimum_feature_size;
};
```

Exact arithmetic makes the policy *deterministic*; it does not replace it.
`radial_clearance` and the TEA cap are early instances — deliberate domain
margins, exactly represented, never numerical repair.

## Case study: the deflation constant that wasn't needed

The incident that produced this page, in three acts:

1. **The defect.** The engagement kernel's first design compared *reported*
   double angles against `cap · (1.0 − 1e-12)` and merged runs whose double
   angular gap was below `1e-12`. Two decisions — a certificate and output
   topology — depended on floating point, inside an Epeck pipeline.
2. **The observation** (review, one sentence): *it's an exact kernel.*
   Both decisions had exact formulations already latent in the data: run
   adjacency is endpoint identity in the arrangement (exact `==` on
   one-root points), and the cap comparison is a chord-length predicate
   against a rational threshold once the cap crosses the boundary as
   `4·sin²(cap/2)`.
3. **The lesson.** Neither fix required new mathematics — only refusing the
   reflex of "handle precision" and asking instead *what is the exact
   question here, and which supported predicate answers it?* That question
   is this page.

## Review checklist

Run every exact-kernel change through these fourteen questions:

1. Is the kernel appropriate for every construction whose result is reused?
2. Does any `to_double()` result affect control flow or topology?
3. Are decimal inputs represented according to their actual semantics?
4. Are kernel predicates used instead of hand-coded determinants and
   epsilons?
5. Are squared quantities used instead of unnecessary roots?
6. Does the selected number type support every required algebraic
   operation?
7. Are all degeneracies and variant intersection results handled?
8. Are tolerance decisions explicit domain policy rather than numerical
   repair?
9. Are exact-to-approximate conversions isolated at named boundaries?
10. Does the code use kernel and traits abstractions (`K::FT`, not
    `CGAL::Gmpq`; no un-profiled `.exact()`) rather than backend types?
11. Does it satisfy the specific CGAL package's traits requirements?
12. Has it been tested on degenerate and nearly degenerate data?
13. Are algebraic roots represented by square-free polynomial plus isolation,
    without minimal-polynomial factoring?
14. Does every normalization used by a sign predicate preserve the original
    polynomial's sign, not merely its zero set?

## References

- [CGAL Developer Manual — Robustness][cgal-robustness] — exact computation
  and "use kernel primitives whenever possible."
- [CGAL Algebraic Kernel manual][cgal-ak] — square-free polynomial,
  isolating-interval representation and why minimal polynomials are avoided.
- [CGAL Polynomial traits][cgal-polynomial-traits] and
  [Algebraic foundations][cgal-algebra] — canonicalization, GCD, square-free
  factorization, fraction decomposition, and algebraic structure concepts.
- [CGAL `Sqrt_extension` documentation][cgal-sqrt-extension] — accessor API
  and same-extension arithmetic preconditions.
- 2D Arrangements & 2D Regularized Boolean Set-Operations — the
  circle-segment traits contract this repository's stock model lives under.
- Repo enforcement: the distilled rules in `CLAUDE.md` (binding for every
  change) and the SP1 design documents under `docs/plans/`.

[cgal-robustness]: https://doc.cgal.org/6.0.1/Manual/devman_robustness.html
[cgal-ak]: https://doc.cgal.org/6.0.1/Algebraic_kernel_d/index.html
[cgal-polynomial-traits]: https://doc.cgal.org/6.0.1/Polynomial/classPolynomialTraits__d.html
[cgal-algebra]: https://doc.cgal.org/6.0.1/Algebraic_foundations/index.html
[cgal-sqrt-extension]: https://doc.cgal.org/6.0.1/Number_types/classCGAL_1_1Sqrt__extension.html
