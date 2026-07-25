# Task 3 single-arrangement reachable-domain design

**Status:** architecture approved; written spec pending review
**Date:** 2026-07-25
**Scope:** replace the rejected, uncommitted Task 3 spike
**Governing plan:** `2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`

## Decision

Build `D`, `C_r = D ⊖ B_r`, `M_r = C_r ⊕ B_r`, and `D \ M_r` once per
reachable-domain request. Construct `C_r` and its certificate from the same
provenance-carrying CGAL arrangement. Certification must never replay the
geometry or recover provenance by scanning every source curve.

Independent full replay remains a Task 13 artifact-verification obligation. It
does not run inside each Task 3 constructor.

## Why the spike is rejected

The uncommitted spike violates the project performance contract before a
profiler is needed:

- one successful owner executes two complete geometry and certificate passes;
- one owner performs `4h + 20` exact set differences and `4h + 22` deep
  `ReachSet` operand copies for `h` holes;
- `C_r` is enumerated six times and the source arrangement is constructed
  twice;
- provenance recovery costs `(E + B)P` exact source-piece tests for `E`
  arrangement edges, `B` selected boundary pieces, and `P` source pieces;
- ring and cycle normalization enumerate and copy every rotation;
- read-only region accessors deep-copy complete exact arrangements;
- coverage repeats sweep construction, residual construction, and
  equality-against-self checks.

The focused Task 3 tests currently request 30 successful owners. The spike
turns those into 60 full geometry passes before accessor clones and coverage
replays.

Running the same deterministic algorithm twice is not independent evidence. It
only proves that two invocations agreed.

## Performance axiom

For one reachable-domain request:

1. canonical input normalization runs once;
2. exactly one provenance arrangement is constructed;
3. `D`, `C_r`, `M_r`, and `D \ M_r` are each constructed once;
4. each arrangement vertex, halfedge, and face is visited a constant number of
   times independent of the number of source curves;
5. certificate ordering is linear plus sorting;
6. read-only region access and alias-safe clones do not copy geometry.

Let:

- `n` be the number of input boundary vertices;
- `p = 4n` be the number of x-monotone offset-line and vertex-circle source
  pieces before intersections;
- `k` be the total arrangement feature count after intersections;
- `z` be the total size of propagated provenance labels;
- `q` be the number of canonical records;
- `b` be the total serialized certificate byte count.

Outside CGAL's unavoidable arrangement construction, the reachable-domain
overhead is bounded by:

`O(n log n + k + z + b log q)`.

There is no `O(kn)` post-construction geometric rematching and no quadratic
rotation enumeration.

For one coverage transition, the exact sweep, accumulated union, and next
residual are each constructed once. A ledger clone shares immutable geometry.

## Exact set model

The input design is one exact polygon with holes:

`D = outer \ union(holes)`.

For every input edge, construct the closed radius-`r` boundary neighborhood
from:

- one exact rectangular edge strip;
- one exact radius-`r` disk at each input vertex.

Their union is the radius-`r` neighborhood of the polygon boundary. The open
arrangement cells selected for `C_r` are inside `D` and outside every boundary
neighborhood primitive. Taking the regularized closure of their union retains
exact tangency:

`C_r = closure(D \ boundary_neighborhood_r)`.

This is equivalent to closed-disk containment for polygonal `D`. The decision
uses exact arrangement topology; it does not sample offsets or introduce a
tolerance.

`M_r` is one batched union of `C_r` with the exact sweep primitives of every
union-boundary curve of `C_r`. `D \ M_r` is one exact difference.

## Canonical identities

Canonical input normalization uses a linear minimal-rotation algorithm for each
ring. Hole rings are then sorted by their canonical records. Reversing,
rotating, or reordering equivalent input rings does not change the result.

Stable IDs never contain handles or iteration order:

- outer ring: `outer`;
- hole ring: its ordinal in the canonically sorted hole sequence;
- source curve: ring ID, edge or vertex ordinal, construction role, and exact
  binary64 radius record;
- source piece: source-curve ID plus its x-monotone piece ordinal;
- boundary primitive: canonical ring ID, feature ordinal, and primitive role.

The certificate already binds the canonical input recipe. Source records
reference its ring IDs instead of embedding a complete ring record in every
source record. This prevents quadratic certificate growth.

The public Python digest remains the sole CCAN content digest. Native records
are deterministic structural bytes, not a second hashing or algebraic-number
serialization scheme.

## Provenance-carrying arrangement

Use CGAL 6.0.1 `Arr_curve_data_traits_2` over the locked
`Gps_circle_segment_traits_2` and exact-with-sqrt kernel.

Each inserted x-monotone curve carries:

- zero or more stable source-piece IDs;
- zero or more closed-boundary primitive IDs.

The traits adapter propagates these labels through curve splitting. Its overlap
merge functor produces sorted unique unions. Curve data never contains
arrangement handles.

Insert, in one aggregate arrangement build:

- the outer and hole boundaries;
- both exact edge-offset sides;
- the construction-only ends of every edge strip;
- both x-monotone halves of every vertex-radius circle.

Design boundaries and construction-only strip ends carry primitive membership
labels but no source-curve identity. Offset sides and vertex circles carry both.

If a selected-cell boundary halfedge lacks source-piece provenance, construction
fails. Such an edge would mean the primitive decomposition or cell-selection
logic is wrong.

## Exact face classification

Crossing a closed primitive boundary toggles membership in that primitive.
Orientation is therefore irrelevant and overlapping boundaries can carry
multiple primitive IDs.

Starting from the unbounded face with every membership bit false, traverse a
dual spanning forest depth-first. Maintain one mutable membership state while
entering and leaving tree edges:

- whether the outer polygon is active;
- the number of active hole polygons;
- the number of active forbidden boundary-neighborhood primitives.

A face is selected exactly when:

- the outer polygon is active;
- no hole polygon is active;
- no forbidden primitive is active.

Each halfedge label is toggled a constant number of times. The traversal does
not copy a complete primitive bitset per face.

Each visited face stores its outer, active-hole-count, and
active-forbidden-count tuple. Reaching it through a non-tree dual edge must
predict the same tuple. A contradiction raises
`ReachableArrangementTopologyError`.

Selected-face adjacency determines connected components. Zero selected
components or more than one selected component raises
`PocketNotMachinableError` under the Phase 1 single-entry contract.

## Geometry and certificate from one DCEL

The classified arrangement is the common source of geometry and proof data.

### `C_r`

Halfedges with a selected face on the left and an unselected face on the right
form the regularized boundary of `C_r`. Extract their cycles once, validate one
outer cycle plus any hole cycles, and construct one exact polygon with holes.

Internal edges between selected cells remain available for cell and adjacency
records but do not enter the union boundary.

### Vertex records

Traverse arrangement vertices once. For each vertex:

1. gather and sort incident source-piece records;
2. group vertices with the same incidence multiset;
3. sort each group by exact `compare_xy`;
4. store each vertex's lexicographic ordinal directly.

Later cell serialization performs map lookup. It never scans an incidence group
or the source list.

### Selected-cell records

Emit one record for every selected arrangement face, not merely one record per
connected polygon-set component. Each boundary cycle contains:

- stable vertex records;
- propagated source-piece records;
- exact curve direction and endpoint roles.

Cycle normalization uses a linear minimal-rotation algorithm. Inner cycles and
cell records are sorted canonically.

### Component records

Component records bind:

- the sorted selected-cell IDs;
- the sorted adjacency edges between selected cells;
- the regularized outer and hole boundary-cycle IDs.

The Phase 1 component count is exactly one.

### Certificate semantics

The rejected spike's `exact_reconstruction` flag is removed. It conflates
structural closure with independent replay.

Task 3 records:

- exact cell-selection closure;
- complete source provenance;
- complete selected-cell and component structure;
- exact `M_r ⊆ D`;
- the construction recipes for `D`, `C_r`, `M_r`, and `D \ M_r`.

Task 13 independently reconstructs the domain from canonical inputs and compares
its content identities during artifact replay.

## Batched exact dilation

Enumerate the extracted `C_r` union boundary once. Each line or circular arc
emits exact sweep polygons:

- line: its exact edge strip plus endpoint disks;
- circular arc: its exact annular sector plus endpoint disks.

Collect all sweep polygons before modifying a set. Use CGAL's range `join`,
which performs a divide-and-conquer overlay, to union:

- the single `C_r` polygon with holes;
- every sweep polygon.

There is no growing sequence of `result.join(part)` global overlays and no
temporary `ReachSet` per boundary curve.

Perform one exact subset decision `M_r ⊆ D`. A failure raises
`ReachableMaterialContainmentError`. Construct `D \ M_r` once afterward.

## Immutable exact-region ownership

`ExactRegion2` stores `shared_ptr<const ReachSet>` plus its role and recipe.
There is no public mutator.

Read-only accessors and `clone()` share the immutable set. A coverage transition
creates new exact output sets and commits them atomically; it never mutates a
shared set.

This makes Python property reads and transaction forks alias-safe without
deep-copying arrangements.

## Incremental coverage

Initial state:

- `A_0` is the exact precleared-entry disk;
- `R_0 = M_r \ A_0`.

For sweep `W_i`:

- construct `W_i` once;
- construct `A_{i+1} = A_i ∪ W_i` once;
- construct `R_{i+1} = R_i \ W_i` once.

The induction is exact:

`R_i = M_r \ A_i`

implies

`R_i \ W_i = M_r \ (A_i ∪ W_i)`.

The operation lineage and parent residual identity certify the transition.
Neither the sweep nor the residual is replayed inside `Coverage2`.

`CoverageLedger.add_sweep` may fork a trial owner, but the fork shares all
immutable geometry. Only the two next-state sets are allocated. A failed exact
operation leaves the original ledger unchanged.

Final coverage remains the exact predicate `R_N.is_empty()`. The separate
unreachable digest still binds `D \ M_r`.

## Error model

Failures are loud and preserve the prior owner state:

- `InvalidReachableDomainInputError`: malformed rings, invalid hole topology,
  or invalid radius;
- `ReachableArrangementTopologyError`: inconsistent parity, missing propagated
  provenance, or invalid selected-cycle extraction;
- `PocketNotMachinableError`: empty or disconnected `C_r`;
- `ReachableMaterialContainmentError`: exact `M_r ⊆ D` failure;
- `InvalidCoverageGeometryError`: nonfinite, degenerate, or invalid-radius
  sweep geometry;
- `CoverageTransitionError`: exact union, difference, or atomic state
  transition failure;
- existing Python certificate and sweep errors validate the native records at
  their ownership boundaries.

No fallback kernel, sampled boundary, rounded disk chain, skip, or xfail is
permitted.

## File responsibilities

Keep each file below the project size thresholds and give it one job:

- `exact_region_2.h/.cpp`: immutable exact-region ownership and basic set
  predicates;
- `exact_sweep_2.h/.cpp`: exact disk, capsule, arc-sweep, and full-circle
  primitives;
- `reachable_arrangement_2.h/.cpp`: primitive generation, propagated metadata,
  face classification, and selected topology;
- `reachable_certificate_2.h/.cpp`: canonical vertex, cell, component, and
  source records;
- `reachable_domain_2.h/.cpp`: public owner and one-pass build orchestration;
- `coverage_2.h/.cpp`: incremental exact coverage state;
- Python modules remain split between reachable-domain ownership and coverage
  lineage/certificates.

## Verification

TDD preserves and strengthens the existing RED fixtures:

1. rectangle, acute, reflex, narrow-neck, island, and disconnected cases;
2. exact tangency and immediate-across-boundary rejection;
3. insertion-order, rotation, and reversal invariance;
4. source-provenance propagation through intersections, splits, and overlaps;
5. one selected-cell record per selected arrangement face;
6. exact selected-cell adjacency and one connected component;
7. exact `M_r ⊆ D` and exact nonempty convex-corner `D \ M_r`;
8. capsule and annulus endpoint, side, radius, and seam membership;
9. coverage order and duplicate semantics;
10. alias-safe region and ledger clones;
11. rounded guide-radius counterexample.

A native algorithm gate records structural work counts and asserts:

- one reachable-domain build pass;
- one provenance arrangement;
- one `C_r` extraction;
- one batched `M_r` union;
- one subset decision;
- one residual difference;
- no replay builds;
- no post-arrangement source-rematching predicates.

These are operation-count assertions, not wall-clock assertions. They are the
falsifiable performance axiom.

Every pytest invocation uses `-n auto`. After one editable native build:

- each isolated diagnostic has a hard 180-second process-group watchdog;
- the complete focused Task 3 pair has a hard 180-second watchdog;
- a timeout stops broader testing and returns to the violated structural stage;
- no profiler runs unless all structural counts pass and a bounded fixture is
  still slow.

If profiling becomes necessary, it tests one written hypothesis about an
unavoidable CGAL stage. It is not used to discover the architecture.

## Non-goals

- changing the locked exact-with-sqrt kernel or CGAL version;
- approximating offsets, sweeps, residuals, or membership;
- test-fixture caching that conceals construction cost;
- independent artifact replay before Task 13;
- modifying accepted Tasks 0 through 2;
- optimizing CGAL internals before the application-level operation bound holds.
