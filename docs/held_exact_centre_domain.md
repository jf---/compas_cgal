# Held exact centre domain

!!! note "Ω\* paths implemented beside the polygon path (2026-10-22)"
    Ω\* is built natively and queried exactly (see
    [Implementation](#implementation)). The Held paths on Ω\* run beside the
    polygon-Ω paths, behind `PathDomain`, and the polygon path stays the
    default. Containment and the gouge audits on these paths use the native
    declared-clearance oracle by default
    ([Integration status](#integration-status-2026-10-22),
    [Held native clearance](held_native_clearance.md)).
    `figure7_inward_offset` stays the Held lane's domain until the
    [acceptance](#acceptance) holds. Removing it needs a separate decision.

**Conclusion.** Ω\* = {x ∈ P : dist(x, ∂P) ≥ r} can be built exactly in the
circle-segment machinery `_stock_2` already uses
(`Gps_circle_segment_traits_2<Epeck>`, one-root points). It cannot be built from
the primitives as they are stored today. Doing it needs one new boundary
contract, the **Pythagorean (offset-rational) boundary**:

- every line has a rational unit direction;
- every arc has a rational centre and a rational radius;
- every number is derived from binary64 inputs by field operations only.

Under that contract every curve of Ω\* is a rational line or a rational circle:

- offset lines;
- offset circles of radius ρ ± r;
- r-circles at junctions;
- ∂P itself.

Every vertex of Ω\* is then a one-root point, so every Gps operation stays inside
the traits' documented one-root contract.

Measured on the four reference cases:

- **Closure.** The re-emitted boundary closes exactly.
- **Displacement.** It moves the published primitives by at most
  2.3·10⁻¹¹ mm. For comparison, the reconstruction certificate against the PDF
  is 0.074–0.102 mm, so the change is about 10⁹ times smaller.
- **Number size.** Its largest number is a 1062-bit squared offset radius
  (Monstera).
- **Cost.** A binary64 cost proxy of the same Gps workload runs in 0.85 s on
  Monstera.

The budget is **2 s for Ω\* on Monstera**. The baseline generation time is
449.8 s, so this is under 0.5 %.

What changes for consumers:

- **Connectors become exact ∂Ω\* sub-chains** (segments and arcs), not chords.
- **Circles stay binary64**, but each one is certified inside Ω\* by an exact
  predicate.
- **The coverage certificate needs one new target member:** a disk with a
  one-root centre, at the convex vertices of Ω\*.

## The decision being made

Notation:

- **P**: the pocket. Its boundary ∂P is the analytic `ReferenceBoundary` (lines
  and arcs, CCW, holes allowed).
- **r**: the tool radius, a binary64 value (1.0 in every case).
- **Ω\***: the tool-centre domain, Ω\* = P ⊖ D_r = P ∖ (∂P ⊕ D_r). It is the
  closure of the set of centres whose tool disk lies inside P.
- **T**: the coverage target, T = Ω\* ⊕ D_r, as approved for the
  [coverage certificate](held_coverage_certificate.md).

Ω\* decides the following:

- contact projection (`benchmarks/held_figure5_path.py:521-532`,
  `benchmarks/held_reference_contacts.py:40`);
- boundary-order ownership and selection;
- the connectors;
- circle containment;
- the base of T.

Today's Ω is a float straight-skeleton offset of the polygonal projection
`case.projection.points` (`benchmarks/held_reference_figures.py:493-510`). The
[gouge audit](#evidence-the-old-domain) shows what that costs:

- **Gouges.** At concave arcs Ω overshoots the true erosion by exactly the chord
  sagitta, so every published path gouges.
- **Unreachable slivers.** At convex arcs Ω falls short, so slivers are left
  that no path can reach.

Shrinking Ω or inflating a margin is excluded. Equivalence with the old Ω is
not a goal, because the old Ω is wrong. Acceptance is defined
[below](#acceptance).

## Evidence: the old domain

Generator: `build/omega-audit/audit.py`. Report:
`.superpowers/sdd/2026-10-13-product-pipeline/omega-gouge-audit-report.md`.

| case | max overshoot (gouge) mm | = max concave chord sagitta | convex sliver area mm² | gouging circles / connectors (> 10⁻⁹ mm) |
|---|---:|---|---:|---|
| figure5 | 0.021009 | yes, to 10⁻¹⁴ | 0.721 | 214 / 138 |
| figure8_upper | 0.028547 | yes | 2.142 | 578 + 193 / 424 + 171 (lane + boundary) |
| figure8_crossed_skis | 0.038065 | yes | 1.028 | 804 / 809 |
| figure8_monstera | 0.046181 | yes | 1.680 | 581 / 319 |

## What the published data gives

Generator: `build/held-exact-omega/boundary_number_census.py`. Output:
`census.json`. Every value is injected exactly as a `Fraction`, and every
decision below is an exact comparison.

| | figure5 | upper | crossed skis | monstera |
|---|---:|---:|---:|---:|
| primitives (lines / arcs) | 31 (5 / 26) | 70 (4 / 66) | 58 (2 / 56) | 317 (103 / 214) |
| convex / concave arcs | 16 / 10 | 39 / 27 | 31 / 25 | 126 / 88 |
| lines with rational length \|b − a\| | 1 | 0 | 0 | 6 |
| arcs with rational radius \|s − c\| | 1 | 0 | 0 | 7 |
| **arcs whose end lies on their start circle** | **0** | **0** | **0** | **0** |
| max end-radius gap (float) mm | 1.8e-15 | 4.0e-15 | 5.3e-15 | 4.9e-15 |
| junctions exactly G1 | 2 | 0 | 0 | 10 |
| junction turns ≤ 1e-12 rad (float) | 18 | 44 | 28 | 116 |
| junction turns 1e-6 … 1° (float) | 13 | 24 | 26 | 176 |
| corners > 1° | 0 | 2 | 4 | 19 (one at −112.7°) |
| arcs with ρ < r (of which convex) | 0 | 2 (2) | 8 (0) | 12 (2) |
| certified reconstruction deviation mm | 0.0738 | 0.0813 | 0.0965 | 0.1023 |

Four consequences follow.

1. **The stored boundary is not a closed exact curve.** No arc's stored `end`
   lies on the circle through its `start` about its `centre`. The audit had to
   choose a closing convention (c + ρ·unit(e − c)), and that convention leaves
   a gap of about 10⁻¹⁵ mm to the next primitive's `start`. An exact Ω\* of the
   data "as stored" is therefore not even well defined. A contract decision is
   needed whatever machinery is chosen.
2. **Radii and line lengths are irrational.** ρ = √|s − c|² is irrational in
   354 of 362 arcs, and |b − a| is irrational in 107 of 114 lines (all four
   cases together).
3. **The boundary is not G1.** Only 12 of 476 junctions are exactly tangent.
   Most of the rest are real kinks of 10⁻⁶ rad to 1° that the biarc
   reconstruction produces. An exact construction must handle every junction
   as a possible kink, so preserving G1 buys nothing.
4. **Some arcs are tighter than the tool.** Upper and Monstera each have two
   convex arcs with ρ < r. There Ω\* has no offset piece, and the material is
   truly unreachable.

## Number fields: as stored versus Pythagorean

| object | as stored | Pythagorean contract |
|---|---|---|
| r-offset of a line through a, b | n·x = C + r√α, α = \|b − a\|² ∉ ℚ²: a line in ℚ(√α) | n̂·x = C + r with n̂ rational: **rational line** |
| r-offset of an arc (c, ρ) | radius ρ ± r with (ρ ± r)² = ρ² + r² ± 2r√(ρ²) ∈ ℚ(√ρ²): **not a rational circle** | radius ρ ± r rational: **rational circle** |
| r-circle at a junction J | J rational, r² rational: rational circle | the same |
| offset line ∩ offset line | ℚ(√α_i, √α_j): two roots | rational |
| offset line ∩ offset circle | root of a quadratic over ℚ(√α, √β): depth 3 | one-root ℚ(√δ) |
| offset circle ∩ offset circle | the same | one-root |
| contact of a centre q (q − r·n(q)) | needs √ per point | **rational affine image of q**, in the same extension |

As stored, each primitive lives in its own quadratic extension, and every
intersection lands in a tower of depth two or three. That is outside the
`Sqrt_extension` contract: arithmetic across different roots is undefined
behaviour (CONVENTIONS.md, [Exactness](exactness.md)).

Under the Pythagorean contract every support curve is rational. Every
intersection is a one-root point, built and compared by
`Gps_circle_segment_traits_2` itself. Comparisons between points with different
roots are exact traits predicates, and no hand-written arithmetic mixes roots.

## Machinery options

| option | can it represent Ω\*? | numbers | cost / risk on Monstera | verdict |
|---|---|---|---|---|
| `Gps_circle_segment_traits_2<Epeck>` on the **stored** data | no: offset circles and offset lines are not rational | — | — | impossible |
| `Gps_circle_segment_traits_2` over `Epeck_with_sqrt` (FT = `CORE::Expr`) | yes: the curves are rational over the FT | nested radicals in every predicate | every stadium flank is tangent to its end disk *identically*, so 317+ identically-zero decisions fall to CORE's separation bound. CORE stalls for minutes on generic doubles in exactly this case (repo memory, *CORE kernel stalls on identically-zero decisions*) | rejected |
| `Circular_kernel_2` (`Exact_circular_kernel_2`) | only rational circles and lines | same limitation as the first row | — | rejected (no gain over Gps) |
| `Minkowski_sum_2` `offset_polygon_2` / `inset_polygon_2` (conic traits, CORE) | **polygons only**; arcs would have to be polygonised again, which is the defect | algebraic, CORE | CORE cost and identically-zero risk | rejected |
| `approximated_offset_2` / `approximated_inset_2` | approximation with an ε bound | rational | cheap | rejected: an approximate certificate is margin by another name |
| `Segment_Delaunay_graph_2` / `src/segment_site_*` MAT, offset traced through Voronoi cells (Held's VRONI method) | segment and point sites only; no arc sites in CGAL or in the repo MAT | — | the repo MAT is bounded slices only (see [segment-site MAT](segment_site_mat.md)) | rejected now; it is the right long-term structure if arc sites ever exist |
| **`Gps_circle_segment_traits_2<Epeck>` on a Pythagorean boundary** | **yes** | rational supports, one-root points, ≤ 1062-bit rationals (measured) | proxy 0.85 s, budget 2 s | **chosen** |

## Boundary contract: the Pythagorean boundary

A new typed stage, `PythagoreanBoundary.build(reference_boundary)`, is derived
deterministically from the recorded `ReferenceBoundary`. It sits beside the
existing stage, which is left unchanged. All of its data is binary64, injected
exactly.

- **`PythagoreanDirection(t, quarter)`.** The unit vector
  u = i^quarter · ((1 − t²), 2t) / (1 + t²), with t a binary64, |t| ≤ tan(π/8).
  Then |u| = 1 exactly, and every rational unit vector has this form. The float
  `atan2` / `tan` only **chooses** t. No exact value is ever derived through a
  float square root.
- **`PythagoreanLine(anchor, direction, length)`.** The segment
  anchor + λu, λ ∈ [0, L], with L a binary64. Its interior normal n̂ = u^⊥ is
  rational, and its r-offset is the rational line n̂·x = n̂·anchor + r.
- **`PythagoreanArc(start, end, start_tangent, sense)`.**
  - Inputs: start S and end E are exact junction points, the start tangent is a
    `PythagoreanDirection`, and N is the unit normal towards the centre.
  - Centre and radius: s = |E − S|² / (2 N·(E − S)), c = S + sN, ρ = s. All
    three are rational.
  - Consequences: |E − c| = ρ holds exactly, and the end tangent (E − c)^⊥/ρ is
    rational too.
  - Degenerate chord: N·(E − S) ≤ 0 raises `DegenerateArcChordError`.
- **Junctions.**
  - An arc ends exactly on its published binary64 junction.
  - A line starts at its predecessor's end and ends at anchor + Lu.
  - Lines only ever accumulate: Monstera's longest run of consecutive lines is 3
    (census), so number size stays bounded.
  - The cycle closes exactly, by construction (`OpenPythagoreanCycleError`
    otherwise).
- **Tangent grid.** t is rounded to the grid 2⁻⁵⁶ℤ. Its resolution,
  2.8·10⁻¹⁷ rad, is below the float direction error, and it caps every
  direction denominator at 113 bits. Without the grid, near-axis directions
  reach 211 bits (measured).
- **G1.** G1 is not imposed. Junctions that were G1 to rounding become exact
  kinks of about 10⁻¹⁶ rad. On ∂Ω\*, reflex kinks become junction arcs of about
  10⁻¹⁶ mm and convex kinks become vertices. Both are exact and ordinary.

Generators: `build/held-exact-omega/pythagorean_reemission.py` (t unrounded) and
`build/held-exact-omega/pythagorean_reemission_tgrid56.py` (t on 2⁻⁵⁶ℤ). Outputs:
`pythagorean_reemission{,_tgrid56}.json`. The numbers below are for the 2⁻⁵⁶ grid.

| case | closes exactly | max line-end move mm | max line rotation rad | max arc centre move mm | max radius change mm | exact G1 / left / right junctions | max bits: junction · line normal · arc centre · arc radius · (ρ ± r)² |
|---|---|---:|---:|---:|---:|---|---|
| figure5 | yes | 7.1e-15 | 3.2e-16 | 3.6e-13 | 3.6e-13 | 1 / 15 / 15 | 163 · 109 · 431 · 322 · 643 |
| upper | yes | 7.3e-15 | 5.3e-16 | 4.4e-12 | 4.4e-12 | 0 / 37 / 33 | 163 · 111 · 435 · 321 · 641 |
| crossed skis | yes | 0 | 4.6e-17 | 2.3e-11 | 2.3e-11 | 0 / 25 / 33 | 167 · 113 · 440 · 323 · 644 |
| monstera | yes | 1.4e-14 | 8.0e-16 | 1.2e-11 | 1.2e-11 | 6 / 140 / 171 | 380 · 113 · 855 · 532 · 1062 |

**Status of the change.** The re-emitted boundary is a *new declared boundary*.
It is not a tolerance applied to the old one.

- **Its own certificate.** Its distance to the PDF source must be certified
  again, by the same reconstruction certificate run on the re-emitted
  primitives. The measured displacement is 2.3·10⁻¹¹ mm against bounds of
  0.074–0.102 mm, so the re-run is expected to pass, but it must be run. The
  stage records the re-certified bound in its `BuildIdentity`.
- **The authority moves.** "Exact" from now on means exact with respect to the
  Pythagorean boundary.

!!! note "Implemented (2026-10-17): the declared boundary"
    The contract above is the declared Held boundary, built beside the stored
    primitives (consumers switch later): `benchmarks/held_declared_boundary.py`
    (types, exact closure, provenance, content hash),
    `held_declared_boundary_emission.py` (re-emission) and
    `held_declared_boundary_certificate.py` (re-certification). Generator:
    `pixi run held-boundary-census` → `benchmarks/results/held_{boundary_census,declared_boundary}.json`.

    - **The boundary moves far less than 2.3·10⁻¹¹ mm.** That figure is the
      arc *centre* and *radius* move, not a boundary move. A rigorous rational
      upper bound on the two-sided Hausdorff move, per primitive, is at most
      2.6e-14 / 1.2e-14 / 1.1e-14 / 3.5e-14 mm (figure5 / upper / crossed skis /
      Monstera).
    - **The bound is composed (B + δ), conditional on one open lemma.** The
      loader re-runs the publisher reconstruction, which gives B; δ is a proved
      exact rational move bound; B + δ is checked exactly against the
      reconstruction limit (stroke / 4). B enters as `nextafter(B, +inf)` of
      the recorded double (declared API semantic: covers a record rounded to
      nearest from below). **Open item `OpenLemma.BIARC_POINT_SET_BOUND`:**
      that B is a two-sided point-set bound for *biarc* spans is read from the
      float-input parameter correspondence, not proved; circle spans are
      proved. `DeclaredBoundaryCertificate.conditional_on` carries the lemma
      and its identity hash includes it. It is not a certificate for
      biarc-derived spans until the lemma is proved or the per-primitive
      certificate is ported to rationals.
    - **The junction split is not stable.** Exact G1 counts match the prototype
      (1 / 0 / 0 / 6). The left/right split differs from the prototype because
      it reduced the quarter index mod 4 before subtracting, which put some t
      on other grid points. The ~1e-16 rad kinks then flip sign.
    - **t is dyadic, not always binary64.** t is k·2⁻⁵⁶ with |k| ≤ 2⁵⁵. For
      |t| ≥ 1/8 such a value need not be a binary64.

## Construction of Ω\*

Native code, in a new `src/held_centre_domain_2.{h,cpp}` inside `_stock_2`. It
holds only the construction (one responsibility per file).

1. **P.** Build P as a Gps polygon with holes, from rational segments and
   rational-circle arcs with rational endpoints.
2. **The members of ∂P ⊕ D_r.** Each one is regular-closed with rational
   supports.
   - *Line band*: the rectangle between the two r-offsets of the support
     (n̂·x = C ± r), closed by the two normal segments through the endpoints.
     Each normal line has a rational direction and passes through a rational
     point, so it is rational. The end caps are the junction disks below.
   - *Arc band*: the annular sector max(ρ − r, 0) ≤ |x − c| ≤ ρ + r between
     the endpoint radii. When ρ ≤ r it degenerates to a pie sector (this is
     structural, not a tolerance).
   - *Junction disk*: D_r(J) at every junction J, rational.
3. **The union.** U = join of all members, as one aggregated
   (divide-and-conquer) `join(begin, end)`, never sequential differences.
   Ω\* = P ∖ U.
4. **Components.** Faces of the result become `CentreDomainComponent`s, each an
   outer CCB with holes. `EmptyCentreDomainError` is raised if there are none.
   More than one component is a hard error, `NeckSplitsCentreDomainError`
   (user ruling, 2026-10-17): a neck no wider than 2r that splits Ω\*. Both
   counts are checked before anything queries the native handle. A neck or
   dead-end slot no wider than 2r that only removes a lobe leaves one
   component; its material is reported as a `NECK` stretch of the
   [unreachable attribution](#implementation).
5. **Degeneracies are constructed, not accidental.** Shared end disks, band ends
   lying on junction-disk diameters, and collinear or concentric neighbours are
   exact one-root equalities. They are decided by the traits, with no CORE.
   Checklist item 17 (zero radicands folded away before an extension is
   constructed) applies to every native predicate written for steps 2 and 6.
6. **Provenance.** Every ∂Ω\* curve lies on a support that is the offset of
   exactly one thing:
   - a line (`OFFSET_LINE i`);
   - an arc (`OFFSET_ARC i`);
   - a junction (`JUNCTION_ARC k`).

   It is recovered by an exact support-key lookup (line coefficients, or centre
   plus squared radius). Coincident supports, such as collinear lines or
   co-circular arcs, carry the full owner set. Gps does not keep curve data, and
   this lookup is exact and needs no change to the traits.

The traits-ownership rules for `Gps` copies apply unchanged
([Exactness](exactness.md), checklist 18).

## ∂Ω\* parameterisation, sites and contacts

- **`CentreDomainPiece`.** A piece is `(component, ccb, index, kind, owner,
  support, source, target)`. Source and target are one-root points, in CCW
  order along the CCB.
- **`CentreDomainSite(component, piece, position)`.** It replaces
  `ProjectionBoundarySite`, whose cyclic modulus is today the vertex count of
  the Ω polygon.
  - Order is lexicographic in (piece, position).
  - On an `OFFSET_LINE`, position is the exact parameter u·(q − anchor). It is
    one-root when q is.
  - On an arc piece, positions are compared by exact orientation about the
    support centre, as in the arc-prefix idiom of [Exactness](exactness.md).
  - Rational witness positions use the stereographic chart of
    [the coverage certificate](held_coverage_certificate.md).
  - Arc length is reporting only.
- **Contact map.** The boundary contact of a centre q is a rational affine image
  of q, so it adds no root:
  - `OFFSET_LINE i`: q − r n̂_i.
  - `OFFSET_ARC i`: c_i + ρ_i (q − c_i)/(ρ_i ± r).
  - `JUNCTION_ARC k`: J_k.
- **Projection of a hypothesis contact.** A binary64 point p is projected to its
  nearest ∂Ω\* piece by exact comparison of squared distances. To an offset
  circle that distance is (|p − c| − ρ′)². Its sign comparisons reduce to
  one-root comparisons by squaring, never to `to_double`. The AABB tree over the
  pieces only prunes, as in [held native clearance](#integration).

## Implementation

- **Native construction**: `src/held_centre_domain_2.{h,cpp}` in its own
  CORE-free extension `_held_centre_domain_2`, not inside `_stock_2`, so a
  rebuild never swaps a library other processes have loaded. It is one aggregated
  `join` of the 2N members (junction disk and band per primitive) and one
  `difference` from P, on `Gps_circle_segment_traits_2<Epeck>`. Boundary cycles
  are read from the arrangement's halfedges, so the travel direction is the
  halfedge's, not the stored curve's. Consecutive halfedges on one support are
  merged into pieces.
- **Input guard** (`benchmarks/held_centre_domain_supports.py`,
  `validate_simple_counter_clockwise`): the native construction does not yet
  check its precondition, and a clockwise cycle there returns garbage
  silently. Python therefore refuses, exactly and before the native call, a
  cycle that is not simple (`NonSimpleDeclaredBoundaryError`: crossings,
  touches and collinear or co-circular overlaps away from a shared junction)
  or not counter-clockwise (`ClockwiseDeclaredBoundaryError`: the turn at the
  lexicographically lowest point). Adjacent primitives meet a second time
  only at a rational point (the second intersection through a rational point
  of a rational line or circle), so the check needs no square root there.
  This is a stopgap until the native build checks `CGAL::is_valid_polygon`.
- **Matching to supports** (`benchmarks/held_exact_centre_domain.py`): every
  piece is matched to the rational supports by an exact key (normalised line
  coefficients, or centre plus R²). Its owners are every support on that key
  it travels along: a coincident support of the opposite sense (an offset arc
  of ρ = 2r about a reflex junction is the junction's r-circle, travelled the
  other way) is not an owner. A piece that matches no support, or that no
  candidate travels along, raises `UnmatchedCentreDomainPieceError`.
- **Queries**: `locate` uses the arrangement's `oriented_side`. Nearest point,
  `squared_distance_to` and `project_contact` / `decide_contacts` are in
  `src/held_centre_domain_query_2.cpp`, which uses CGAL's static-filter idiom
  per predicate:
    - certified intervals decide first, with vertex coordinates enclosed
      tightly at build time by an exact bracket of √r;
    - an undecided predicate is re-decided on the exact rationals;
    - a certified box per piece only prunes.

  Contact projection reports a vertex on the piece that starts there. It
  raises `ContactBeyondGapError` or `AmbiguousContactProjectionError`. The
  pure-`Fraction` module `benchmarks/held_centre_domain_distance.py` is its
  differential oracle, compared in both directions (native unique ⇔ oracle
  unique, distinct nearest points deduplicated exactly).
- **Batch projection** (`CentreDomain.decide_contacts`): one native call
  decides every contact and returns `ContactDecisions`, which carries the
  number of predicates the interval stage left to the exact stage. The count
  is read right after the call, so no caller reads the extension's module
  global. Each decision is a `DecidedContact`, whose one-root site is built
  lazily by `site()` in pure `Fraction` arithmetic and with no native call: a
  vertex is the piece's source, and an interior point is the orthogonal
  projection onto the piece's support. `project_contacts` returns the same
  lazy contacts.
- **Vertex disks**: a vertex disk sits at every convex vertex of a cycle, which
  is a strict left turn or a π-turn cusp (cross = 0, dot < 0). A π-turn is
  always a convex spike of Ω\*, because Ω\*'s complement is a union of open
  r-disks, so it has no inward cusp. `VertexDisk.owners` are the owners of the
  two meeting pieces, because a neck or a tight-arc crossing has no single
  junction. T itself is not built: its vertex disks have one-root centres.
- **Unreachable attribution, by unreached boundary.**
    - *Rule.* A boundary point b is reached iff its inward centre
      b + r·n̂ lies in Ω\*. Per primitive, the reached set is the closed
      contact-map image of the ∂Ω\* pieces it owns, plus the points vertex
      disks touch (a foot or a junction at distance exactly r from a vertex),
      plus the whole arc for an arc with ρ = r whose centre the native
      `locate` puts in Ω\*. The complement intervals are the unreachable
      boundary (`UnreachableInterval`, exact one-root ends, open or closed).
    - *Completeness.* For x ∈ P ∖ T with nearest boundary point b, the
      inward centre of b is not in Ω\*, else x would lie in its disk. So every
      part of P ∖ T touches a reported interval, including dead-end slots no
      wider than 2r, which a structural rule misses.
    - *Cause.* Intervals connected through unreached junctions form one
      stretch (`UnreachableAttribution`). A stretch is local when the Ω\*
      points that reach its two ends share one vertex. A local stretch with a
      convex arc of ρ < r is `TIGHT_CONVEX_ARC`; otherwise a local stretch with
      a convex kink is `CONVEX_KINK_CUSP`. Every other stretch is `NECK`.
    - *Exact ordering.* Order along a line and along an arc of sweep ≤ π uses
      a scalar key: u·(p − S), or the chord projection (E − S)·(p − S), which
      increases strictly along such an arc. Arc membership is one chord-side
      sign. `compare_on_arc` is needed only on major arcs.

Measured by `pixi run held-centre-domain-timing` →
`benchmarks/results/held_centre_domain_timing.json`. Each figure is the median
of 3 samples, each in a fresh process. The machine was heavily loaded (1-min
load average about 347), so wall-clock figures are inflated; CPU is shown
beside them. The table predates the input guard and the boundary attribution,
which add exact Python work to every build. Refresh it by re-running the task.

| case | pieces | Ω\* build wall / CPU s (native) | hypotheses | projection wall / CPU s | exact predicates |
|---|---:|---|---:|---|---:|
| figure5 | 46 | 0.103 / 0.079 (0.065) | 31,742 | 0.575 / 0.276 | 6,243 |
| upper | 104 | 0.278 / 0.188 (0.194) | 62,730 | 1.114 / 0.599 | 25,732 |
| crossed skis | 92 | 0.168 / 0.138 (0.119) | 32,204 | 0.585 / 0.391 | 12,283 |
| monstera | 487 | 1.262 / 0.805 (0.898) | 120,272 | 1.222 / 0.908 | 19,342 |

The exact predicates are intrinsic, not filter slack. The hypothesis contacts
come from a float offset of the same lines, so many lie within ulps of ∂Ω\*,
and there candidate distances tie inside any double enclosure.

## Integration

Every consumer gets a new path *beside* the old one. The old path is removed
only after acceptance, and only with a separate decision.

| consumer | today | new |
|---|---|---|
| `figure7_inward_offset` (`benchmarks/held_reference_figures.py:493`) | float straight skeleton of `projection.points` → `InwardOffsetEvidence` polygon | `exact_centre_domain(case) -> ExactCentreDomain` in a new `benchmarks/held_exact_centre_domain.py`, backed by native `CentreDomain2` |
| contact projection (`held_figure5_path.py:521-532`, `held_reference_contacts.py:40`, `_circle_geometry_2.project_boundary_contact`) | nearest Ω polygon edge, `(side, Fraction(parameter))` | `CentreDomain2.project_contact(p, gap) -> CentreDomainSite`, with the same native distance-gap rejection semantics |
| canonical offset point and circle placement (`held_figure5_boundary_path.py:180-196`) | contact on the polygon edge; centre shifted in floats | exact contact; centre and guide radius stay binary64; **containment is decided exactly** (next row) |
| containment (`benchmarks/held_containment_report.py::containment_report`, the gate's "containment rejections") | circle against the Ω polygon | `circle_in_centre_domain(c, R)` ⇔ annulus max(R − r, 0) ≤ \|x − c\| ≤ R + r lies in P. This is decided against the Pythagorean primitives directly (line: rational; arc: one-root after squaring), with no Ω\* construction needed |
| ownership and order (`_boundary_neighborhood_owners(..., len(component), ...)`, `_site_key`, `held_figure5_path.py:362-367`) | cyclic polygon vertex index | cyclic `(piece, position)` order per CCB |
| connectors (`_transition`, `held_figure5_boundary_path.py:210-250`; `_connectors` in `held_toolpath_benchmark.py:173`, `held_gouge_workloads.py:54`) | contact → Ω vertices → contact polyline; progress from float side lengths | `ExactConnector`: the CCW sub-chain of ∂Ω\* pieces between two sites, made of segments and **arcs**; progress is reported as an exact-order chain with a float length |
| refinement and spacing (`refine_figure5_engagement`, `refine_boundary_ordered_path`, `select_contour_bound_path`, `build_coverage_preserving_path`, `_corner_radius_limit`) | take the Ω polygon | take `ExactCentreDomain`; each decision they make on the polygon is ported to a piece predicate, and engagement *metrics* may read the reported float view |
| coverage target T ([certificate](held_coverage_certificate.md)) | Ω ∪ ⋃ Cap(e_k), Ω a binary64 polygon | see [Coverage target](#coverage-target-and-the-exact-path) |
| native AABB clearance (`src/held_boundary_clearance_2*`, not yet wired) | stations against the **projection ring** | the same prune-then-decide design against Pythagorean sites. Arc sites are a new site kind: the distance to an arc in its wedge is \|\|p − c\| − ρ\|, compared by exact squaring; the wedge is decided by orientation |

## Integration status (2026-10-22)

The Ω\* consumers in the table above are built as the **Ω\* path**, beside
the polygon path:

- **Selection.** `PathDomain` in `benchmarks/held_omega_toolpath.py`, and
  `--domain omega` on the benchmark CLI. The polygon path is the default and
  is unchanged.
- **Modules.**
  - `held_omega_sites.py`: exact ∂Ω\* sites and their order.
  - `held_omega_circles.py`: binary64 circles on the inward normal, contained
    exactly.
  - `held_omega_connectors.py`: exact ∂Ω\* sub-chains, segments and arcs.
  - `held_omega_path.py` and `held_omega_refinement.py`: initial path and
    refinement.
- **Clearance.** Every clearance decision on these paths is one exact oracle
  against the declared primitives. This covers circle containment in Ω\*
  during generation and both gouge audits (`held_omega_gouge_audit.py`).
  - The oracle runs natively by default (`ClearanceBackend.NATIVE`,
    `src/held_declared_clearance_2/`).
  - The Python oracle `held_declared_clearance.py` stays as the differential
    reference (`ClearanceBackend.PYTHON`).
  - The two make identical decisions on all five workloads, on the polygon
    path and on the Ω\* path.
  - Details, the owner-by-construction TOUCH lemmas, the evidence and the
    timings: [Held native clearance](held_native_clearance.md).
- **Owners.** `CentreDomainPiece.owners`, `VertexDisk.owners` and
  `UnreachableAttribution.owners` are `OwnerSet`s. They iterate and print in
  canonical order (kind, then index), so Ω\* output and its digests do not
  depend on `PYTHONHASHSEED`.

### Native site and placement

`CentreDomain2.place_binary64(xs, ys, pieces, vertices, guide_radii)`
(`src/held_centre_domain_placement_2.{h,cpp}`, Python
`benchmarks.held_omega_circles.place_local_batch`) places, in one native call,
the tangent circle of every decided contact. It makes the decisions that
`OmegaSite.decided` + `place_local` make in Python, and returns the same
circle field for field (contract test
`tests/benchmarks/test_held_omega_native_placement.py`: the batch against
`place_local` on line, concave-arc and convex-arc pockets and on Hypothesis
probes, plus the round-to-nearest-even and sub-ulp near-tie witnesses).

- **Site.** The site's placement (source vertex or interior) and the
  outside-the-piece verdict are exact signs: the line parameter `u·q` against
  `u·s`, `u·t` (one-root ends), or on a circle the half turn of `w = q − C`
  from `s − C` and the before-target cross sign. Each is filtered by a
  certified `Interval_nt` enclosure and decided exactly when uncertain.
- **Support data.** Every line piece of ∂Ω\* is an r-offset of a declared
  line, so it carries that line's exact unit normal `n̂` (matched by exact
  slope, oriented left of travel) and `level = n̂·x` — the declared form, whose
  rationals are small. The arrangement's own `(a, b, c)` is used only to find it.
- **Centre.** The correctly rounded (nearest, ties to even) binary64 of each
  exact coordinate of `c* = q + t·n`: `p + (t − h)·n̂` on a line,
  `C + (ρ' − σt)·w/√|w|²` on a circle, `(1 − k)·s + k·C` at a vertex. A
  candidate double is accepted when the exact value lies between its rational
  midpoints with both neighbours; `nextafter` only enumerates candidates. This
  removes the Python arc proposal's dependence on the probe's enclosure.
- **Radius.** The guide radius when `t ≤ bound`, else the largest binary64
  `≤ bound`, with `bound = n̂·c − level` on a line and `√D ∓ ρ'` on a circle
  (rational predicates `(t + ρ')² ≤ D`, `D ≤ (ρ' − t)²`); the rounding-down is
  exact. Inside a convex arc of radius `ρ' < t`, `t` is first capped at the
  largest double `≤ ρ'` (`limited`).
- **Reporting.** The normal view (`n̂` rounded on a line,
  `−σ·fl(w)/√fl(|w|²)` on a circle) is returned so the contact point is the
  same binary64 expression as before; it decides nothing.
- **Statuses.** `OUTSIDE_PIECE` and `PROBE_AT_CIRCLE_CENTRE` are raised by
  the Python wrapper as `SiteOutsidePieceError` / `ProbeAtCircleCentreError`;
  a radius ≤ 0 as `NoContainedCircleError`. A line piece with no rational unit
  normal raises `NoRationalUnitNormalError` (native).

The exact number type is `Epeck::FT::ET` (boost `cpp_rational`; GMP is
disabled in this build), whose every operation reduces by gcd. So the circle
branch does its per-placement arithmetic on unreduced integer numerator /
denominator pairs (CGAL `Fraction_traits`) and compares by
cross-multiplication; the radius bound's search is seeded without
cancellation, `(D − ρ'²)/(√D + ρ')`.

**Out of scope, but on the same defect.** The raw guide and the
`ProjectionAdmissibleBoundaryHypothesis` stations still read the projection
ring. Their admissible gap 2·(reconstruction + projection) absorbs the
projection error as *evidence*: they propose hypotheses, and do not decide
machining. They can stay until the projection is retired. Every decision that
cuts metal moves to Ω\* in the table above.

## Coverage target and the exact path

T = Ω\* ⊕ D_r, decomposed into members for the certificate's local theorem:

- **Ω\* itself.**
- **The outward band of every ∂Ω\* piece.** These are the regions between the
  piece and the primitive it offsets, and they are rational: the outer flank is
  the original line or circle, and the ends are normals.
- **D_r(w) at every convex vertex w of Ω\* that is not G1.** These are convex
  boundary kinks (140 in Monstera after re-emission) and neck crossings. w is a
  one-root point, so **D_r(w) is a circle with an irrational centre**. That is a
  new target piece type for the certificate, with one-root witness points
  instead of rational ones.

**Consequence.** Points of ∂D_r(w) in the normal fan at w are reachable only
from w itself. A tool path whose centres are all binary64 therefore cannot cover
T exactly, because w is irrational. The uncovered depth is about
|x − w|²/(2r), around 10⁻³⁰ mm, but the certificate correctly decides it as
uncovered. The same holds for tangency to a Pythagorean arc at a binary64 point.

So the certified object is the **exact path**:

- **Circles** stay binary64 and must lie in Ω\*, decided exactly.
- **Connectors** are exact ∂Ω\* chains. They contain every vertex w and run on
  every band's inner flank, so T ⊆ U is decided exactly with coincident
  supports.
- **Emission to the machine is a separate, reported quantisation.** Connector
  arcs and segments are rounded to G-code doubles, and the stage records the
  measured maximum displacement. It feeds no decision.

**Implemented (2026-10-21).** `CoverageCertificate.build(..., exact_target=...)`
certifies T exactly, with these members:

- the bands as above: line bands, annular sectors for offset and junction arcs;
- the one-root vertex disks;
- Ω\* membership as p ∈ P and dist(p, ∂P) ≥ r on the declared primitives.

The declared supports and the vertex circles cut every swept piece. Instead of
sweeping target pieces, one interior point of T is checked for cover: int T is
connected because Ω\* is one component.

Exact ∂Ω\* connectors are not modelled yet. They coincide with the cut
supports by design and raise `CoincidentTargetSupportError`.

Details and the Monstera measurement:
[Held coverage certificate, the exact target](held_coverage_certificate.md#the-exact-target).

## Acceptance

Acceptance is *not* equivalence with the old Ω.

1. **Exact gouge audit: 0 gouges.** The method of `build/omega-audit/audit.py` is
   ported to Pythagorean primitives and to exact connector pieces (arc pieces
   need radial-range candidates against arcs, which the audit already has for
   circles). The requirements:
   - every circle and every connector piece of every Held workload has
     clearance ≥ r exactly;
   - no escalation is left undecided.

   The emitted binary64 path is audited separately. It reports its maximum
   rounding-scale gouge, which must be ≤ the measured emission displacement.
   That number is reported, not gated by a tolerance.
2. **No unreachable sliver except true unreachability.**
   - For every primitive piece that Ω\* reaches, exact clearance holds:
     dist(∂P piece, Ω\*) = r.
   - The material left out of T must be exactly the set P ∖ (Ω\* ⊕ D_r). Every
     part of it touches an unreached stretch of ∂P, and each stretch is
     attributed to one of: a convex kink cusp, a convex arc with ρ < r, or a
     neck (see [Implementation](#implementation), "Unreachable attribution").
   - The analytic sliver areas of the old Ω (0.72 / 2.14 / 1.03 / 1.68 mm²)
     must disappear, except for these attributed parts.
3. **Figure counts are reported, old beside new.** For each case: circles,
   connectors, planar length, max engagement, containment rejections, Ω\*
   components, pieces by kind, convex vertices. They are reported with
   old-versus-new figures shown before any commit, and nothing is gated on
   being equal.
4. **Coverage.** The certificate decides T ⊆ U exactly on the exact path for all
   four cases.
5. **Time.** The budget below holds.

## Time budget and cost model

The budget is **Ω\* construction ≤ 2.0 s on Monstera** (M1 Max, one process,
`_stock_2` release build). Contact projection for all Monstera hypotheses adds
≤ 1.0 s. If construction exceeds 5 s, implementation stops: profile, report,
and do not change the representation without a decision.

The cost model, with N primitives and J junctions:

- **Members.** About 2N + J (Monstera: 317 bands plus 317 junction disks, about
  1,600 curves before x-monotone splitting).
- **The union.** One aggregated join is O((n + k) log n). The intersection
  count k is linear: each band meets its two neighbours and two disks, plus
  O(#necks) far interactions.
- **Predicates.** Most are filtered by Epeck's interval stage. The exact
  fallbacks are the *constructed* degeneracies, about 10–20 per junction, on
  rationals of up to about 1,100 bits plus one root; at 20–100 µs each that is
  0.1–0.6 s.
- **Measured proxy.** Generator: `build/held-exact-omega/gps_cost_proxy.py`, which
  wrote `gps_cost_proxy.json`. It runs the existing `Stock2` with binary64 data:
  the projection polygon, minus a junction disk and a line sweep per primitive.

  | case | Gps ops | time | final arrangement |
  |---|---:|---:|---|
  | figure5 | 62 | 0.031 s | 40 vertices |
  | upper | 140 | 0.078 s | 100 vertices |
  | crossed skis | 116 | 0.058 s | 68 vertices |
  | monstera | 634 | 0.85 s | 454 vertices |

  The proxy errs in both directions. It uses sequential differences, which
  overstate the work an aggregated join does. It uses 53-bit numbers, which
  understate the per-predicate cost of the real ~1,000-bit data.
- **Context.** Monstera's whole generation takes 449.8 s (`docs/held_performance.md`,
  `benchmarks/results/held_toolpath_baseline.json`). The budget is set by the
  rule that generation time is a gate, not by this headroom.

## Tests

Every fixture is Pythagorean by construction, using integers and 3-4-5
triples, so expected values are exact rationals or one-root numbers. In every
fixture r = 1 unless the case says otherwise.

**Contract, `tests/benchmarks/test_held_exact_centre_domain.py`** (reference
tests; the design's original list, with these names):

1. `test_concave_arc_pocket_old_domain_gouges_new_does_not`
   - Pocket: (0,0) → an arc about (10, −7.5), ρ = 12.5, clockwise → (20,0) →
     (20,20) → (0,20).
   - Old: the old Ω overshoots by the chord sagitta.
   - New: an `OFFSET_ARC` of radius 13.5 with apex (10, 6), and the one-root
     corner (19, −7.5 + √101.25).
2. `test_convex_arc_pocket_old_domain_leaves_a_sliver_new_reaches_it`.
3. `test_reflex_corner_has_a_junction_arc`
   - The L-shape gets a `JUNCTION_ARC` about (10,10) from (10, 9) to (9, 10).
   - (9.4, 9.2) lies on ∂Ω\*.
4. Necks:
   - `test_neck_of_width_exactly_2r_is_a_hard_error`: two components, so
     `NeckSplitsCentreDomainError`;
   - `test_neck_wider_than_2r_keeps_one_component` at width 2 + 2⁻²⁰.
5. `test_convex_arc_tighter_than_tool_has_no_offset_piece`: ρ = 0.5 gives no
   offset piece and is attributed `TIGHT_CONVEX_ARC`.
6. Over the four cases:
   - `test_centre_domain_boundary_closes_exactly[case]`;
   - `test_centre_domain_is_built_on_the_recertified_boundary[case]`.
7. `test_generic_double_witness` (Hypothesis): generic binary64 pockets build
   within budget, so identically-zero stalls are caught.
8. One-root number and site order:
   - `test_one_root_numbers_fold_and_compare_exactly`;
   - `test_cross_root_arithmetic_and_mixed_points_are_refused`;
   - `test_tool_radius_must_be_a_positive_binary64`;
   - `test_junction_arc_sites_are_ordered_exactly_across_roots`;
   - `test_contact_map_adds_no_root`.

**Queries, `tests/benchmarks/test_held_centre_domain_queries.py`.**

- Projection:
  - `test_reflex_corner_projects_onto_the_junction_arc_at_a_one_root_point`;
  - `test_projection_beyond_the_gap_raises`;
  - `test_a_vertex_site_is_reported_on_the_following_piece`;
  - `test_equidistant_distinct_contacts_are_ambiguous`;
  - `test_concave_offset_arc_projection_lands_on_the_offset_circle`;
  - `test_arc_with_radius_equal_to_the_tool_is_reached_not_unreachable`.
- `test_native_construction_rejects_a_non_unit_line_normal`.
- Differential oracle, in both directions with exact dedupe:
  - `test_native_nearest_matches_the_python_oracle` (Hypothesis) over the L,
    concave, convex, full-circle, stadium and rich fixtures.
- Tests that kill the defects the review found surviving:
  - **gap ties.** `test_contact_exactly_at_the_gap_is_projected_not_beyond`
    covers a rational tie and a tie the exact stage must decide.
    `test_batch_decisions_return_their_exact_fallback_count` checks the count.
  - **at-centre filter.** `test_query_at_a_full_circle_centre_is_ambiguous`.
  - **`piece.full` path.** `test_full_circle_piece_projects_radially`.
  - **antipode branch.** `test_half_turn_piece_decides_both_sides_of_its_diameter`
    uses the stadium cap offsets, which span exactly π.
  - **arc box extremes.** `test_queries_beyond_every_arc_extreme_match_the_oracle`,
    on the fixtures, figure5 and crossed skis.
  - **tight enclosures.**
    `test_queries_off_irrational_vertex_directions_by_less_than_an_ulp_match_the_oracle`
    queries 2⁻⁷⁰ off one-root vertex directions.
  - **batch versus single.** `test_batch_projection_equals_single_projection`
    checks lazy sites against single projections.

  Native mutants cannot be built in the Python lane, so the C++ slot must
  confirm that each of these native defect classes is killed.

**Attribution and input,
`tests/benchmarks/test_held_centre_domain_unreachable.py`.**

- The review's 2r slot probe:
  `test_dead_end_slot_exactly_2r_wide_is_attributed_to_a_neck` gives exact
  intervals, walls and cap. `test_dead_end_slot_wider_than_2r_is_reached`
  covers 2.5r.
- Exact cusp and tight-arc intervals:
  - `test_square_corner_cusp_intervals_are_exact`;
  - `test_tight_convex_arc_part_spans_the_arc_and_the_line_tails`.
- The ρ = r exemption: `test_equal_radius_cap_is_exempt_only_when_its_centre_is_in_the_domain`.
- The lemma: `test_every_unreached_interval_is_disjoint_from_the_reached_contacts`.
- π-turn cusps:
  - `test_pi_turn_cusp_of_the_domain_gets_a_vertex_disk`;
  - `test_vertex_disk_reaches_an_isolated_boundary_point`, a wall reached at a
    single point.
- Multi-owner sense: `test_coincident_support_of_opposite_sense_does_not_block_the_match`.
- Co-circular owners: `test_co_circular_arcs_share_one_piece_and_are_both_reached`.
- Input guards:
  - `test_clockwise_cycle_is_a_named_error`;
  - `test_crossing_cycle_is_a_named_error`;
  - `test_self_touching_cycle_is_a_named_error`;
  - `test_tool_that_fits_nowhere_is_an_empty_domain_error`.

**Acceptance.** The exact gouge audit of the published paths is
`benchmarks/acceptance/test_held_gouge_audit_figure5_slow.py` and
`pixi run held-gouge-audit`. On the Ω\* paths, acceptance item 1 is measured
by the differential audit in
[Held native clearance](held_native_clearance.md#evidence): 0 gouging circles
and 0 gouging connectors on all five workloads. Item 2 waits for the
coverage certificate to model exact ∂Ω\* connectors.

## Rulings (2026-10-17)

The four questions the design left open have been ruled on:

1. **The boundary contract.** Yes. The Pythagorean re-emission
   (`benchmarks/held_declared_boundary.py`) *is* the boundary, and "exact"
   means exact with respect to it.
2. **The exact path.** Yes. The certified object is the exact path: exact ∂Ω\*
   connectors and binary64 circles. Machine emission is a reported
   quantisation.
3. **Multiple components.** No. A neck no wider than 2r that splits Ω\* is a
   hard error (`NeckSplitsCentreDomainError`). A neck that only removes a
   lobe is attributed as `NECK`.
4. **The coverage certificate.** It gains the one-root-centre `VERTEX_DISK`
   member. The certificate API accepts `vertex_disks` today, but refuses to
   certify while their one-root circles are not cut into its arrangement.

## Reproducibility

Every number on this page names its generator. Run each pixi task from the
worktree root.

- **Ω\* timing table.** `pixi run held-centre-domain-timing`
  (`tools/held_centre_domain_timing.py`) →
  `benchmarks/results/held_centre_domain_timing.json`.
- **Data census.** `pixi run held-boundary-census` runs both:
  - `tools/held_boundary_census.py`, promoted from
    `build/held-exact-omega/boundary_number_census.py`;
  - `tools/held_declared_boundary_census.py`, promoted from
    `build/held-exact-omega/pythagorean_reemission_tgrid56.py`. It gives the
    re-emission closure, displacement and bit sizes.
- **Exact gouges of the published paths.** `pixi run held-gouge-audit`
  (`benchmarks/held_gouge_audit_table.py`). This supersedes the gouge part of
  `build/omega-audit/audit.py`.
- **Native clearance differential and Ω\* gouge counts.**
  `pixi run held-native-clearance-results -- --measure`
  (`tools/held_native_clearance_results.py`) →
  `benchmarks/results/held_native_clearance.{json,md}`.
- **Monstera generation time.** `docs/held_performance.md` /
  `benchmarks/results/held_toolpath_baseline.json`.
- **Not reproducible from the tree.** Two design-time figures came from
  generators that were never promoted:
  - the Gps cost-proxy table, from `build/held-exact-omega/gps_cost_proxy.py`;
  - the old domain's analytic sliver areas (0.72 / 2.14 / 1.03 / 1.68 mm²),
    from `build/omega-audit/chord_area.py`.

  They are kept as the record of the decision. The cost proxy is superseded
  by the measured timing table above. Promote these generators to `tools/`
  before they back any acceptance claim.
