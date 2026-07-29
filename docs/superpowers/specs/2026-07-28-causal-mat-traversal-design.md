# Causal MAT traversal and covered-generation design

## Scope

Task 13 turns the exact MAT, classified neck inventory, finite candidate
lattice, qualified entry, and atomic Task 12 evaluator into one complete
covered-generation loop. It owns:

- one well-founded directed cursor for every MAT edge;
- deterministic component and branch traversal from stable graph identities;
- causal neck-scope assignment from exact separating-cut side history;
- the distinct first-circle launch from authenticated precleared stock;
- joint link/circle search and commit after launch;
- terminal traversal plus exact empty-reachable-residual acceptance; and
- fresh replay of the completed operation stream.

Task 13 returns an internal `GenerationResult`. Task 14 remains responsible for
the public `CertifiedToolpath`, build provenance, and artifact identity.

The current production MAT backend is the exact polygonal L-pocket backend.
Task 13 must complete that supported backend without presenting the result as
arbitrary-pocket coverage. The design nevertheless preserves cycles,
multi-partition necks, and multiple MAT components as explicit state. Any
topology that cannot be assigned by the exact projection fails with a named
error; it is never approximated from coordinates.

## Interface audit

The Task 12 boundary is locally complete but intentionally insufficient for
global generation:

1. native `MatExactNeckEvidence2` retains the complete ordered
   `edge_partitions()` and exact four-variant neck location;
2. the current 20-field numeric MAT projection flattens those partitions into
   one cut-edge union;
3. the two production L-pocket plateau necks each have **three**, not two,
   certified partitions;
4. the reproduced Task 12 candidate edge belongs to the large partition of
   both necks, so cut-union membership cannot identify a crossing;
5. `GenerationState` correctly owns the last accepted physical boundary, but
   its single cursor is not the global per-edge traversal ledger;
6. the first circle has no direct link and must be certified inside the
   precleared entry before ordinary link/circle transactions can begin; and
7. canonical MAT edge direction does not in general agree with the direction
   in which a connected traversal discovers that edge.

These are architectural boundaries, not missing conditionals. Task 13 adds
separate lifetimes and one exact topology projection rather than expanding
`GenerationState` into a god object or parsing native certificate bytes in
Python.

## Additive exact neck-topology projection

The frozen 20-field numeric projection remains unchanged. The retained native
`SegmentSiteMedialAxis` owner gains an additive read-only neck-topology
projection derived from the same `MatNeckEvidenceV1` objects that produced the
canonical evidence bytes.

For every classified neck it exposes:

- the exact location tag;
- canonical locus edge IDs;
- canonical locus node IDs;
- the strict/endpoint algebraic parameter-root identity when one exists; and
- the complete canonical tuple of separating-cut edge partitions.

Python maps the four native location variants to one closed typed union:

- `StrictEdgeNeckLocus`;
- `ClearanceEndpointNeckLocus`;
- `SharedVertexNeckLocus`; and
- `PlateauNeckLocus`.

`NeckSide` binds one neck evidence digest, one canonical partition, and a
SHA-256 side identity. `ClassifiedNeck` retains the existing cut union and also
owns its exact `locus` and `sides`. Construction requires the union of all
projected side edges to equal the existing numeric cut union. This
cross-check keeps the established projection useful while preventing it from
being mistaken for the deciding topology.

No Python decoder interprets `evidence_bytes`. Native construction, canonical
evidence serialization, classification, and traversal topology all consume
the same exact evidence object.

## Global traversal state

Static and evolving facts have separate lifetimes.
`TraversalSampleIndex` binds every refinement-dependent exact parameter,
cursor, edge, and ordinal once. `TraversalGraph` projects exact component,
branch, node, edge, and endpoint-cursor authority. `MatTraversalAuthority`
derives both plus the deterministic route and complete causal neck-frontier
schedule from the exact owners, policies, and authenticated entry exactly
once. It binds the MAT sampling policy separately from its realized sample
index so coincident outputs from distinct policies cannot alias.

`MatTraversalState` is immutable and content-addressed. It retains that static
authority and only the state that evolves on the graph lifetime:

- one `DirectedEdgeCursor` per exact MAT edge;
- the active edge;
- explicit visited `(edge, node)` incidences;
- the current exact side of every classified neck; and
- at most one pending causal neck transit.

`DirectedEdgeCursor` owns:

- component, edge, and branch identities;
- entry and exit node identities;
- one native or proof-carrying derived cursor;
- accepted-candidate count; and
- terminal state.

The entry/exit node pair carries direction structurally. There is no
coordinate-derived adjacency and no boolean whose meaning changes with call
site. A cursor can advance only through a candidate whose
`AdvanceTraversalDecision` starts at its exact identity. Every accepted
candidate advances exactly one edge cursor; rejected candidates change none.

`GenerationState` remains the physical accepted boundary: stock, coverage,
phase, operation prefix, neck passages, and the cursor named by its final
circle. `MatTraversalState` remains the global graph boundary. A
`TraversalCommit` cross-binds the parent and child digests of both. This split
keeps stock/coverage replay independent from graph-frontier evolution while
making an atomic generator step authenticate both.

## Directed spans

Connected traversal may discover an edge from either exact endpoint.
`MiddleCurveSpan` therefore supports monotonically increasing or decreasing
native sample ordinals. Direction is inferred from the ordered cursor pair and
retained by `DerivedCandidateCursor`; no second candidate type or reversed
geometry copy is introduced.

Candidate progress remains positive in the selected direction. The span's
terminal limit is the endpoint reached in that direction. Candidate identity
already binds cursor-before, cursor-limit, exact progress, and cursor-after,
so forward and reverse candidates cannot collide.

Fresh replay enumerates the same directed window and requires a unique
candidate match. It never guesses direction from circle orientation or
reporting coordinates.

## Route and side causality

The route planner operates only on exact node and edge identities:

1. order MAT components with `TraversalPolicy.order_components()`;
2. seed the component containing the entry-launched candidate;
3. discover incident unfinished edges through shared node IDs;
4. order available branches with `TraversalPolicy.order_branches()`;
5. choose edge direction from the shared discovery node;
6. retain both visited incidences explicitly so a cycle-closing edge cannot be
   dropped; and
7. move to the next component only after the current component has no
   unfinished edge.

A neck side is a complete separating-cut partition, not one edge-membership
test. For each planned transition, the traversal authority carries the
previous exact side through any neck locus and resolves the next unique side.
Different source and target side identities create one `CausalNeckTransit`.
The transit binds neck owner, source side, target side, and the canonical
orientation:

- `FORWARD` when source-side identity sorts before target-side identity;
- `REVERSE` otherwise.

This ordering is exact and deterministic for two or more partitions. It does
not pretend a three-way plateau is binary. The independent forward/reverse
`NeckPassage` histories remain conservative per owner and orientation.

An operation receives `NoNeckScope` unless the global route proves that it is
the accepted motion carrying one causal transit. A transit supplies the
current passage's `NeckCapDecision`; the passage advances only when that
candidate commits.

If one motion would cross multiple independently active necks, or if a locus
does not resolve unique source and target sides, generation raises a named
topology error. The current operation schema owns one exact neck decision; it
must not silently choose one restriction, merge caps without provenance, or
infer a side from geometry.

## Entry bootstrap

`InitialCandidateEvaluator` owns the physically distinct launch transaction.
It enumerates entry-compatible candidates across directed native spans and
accepts only candidates whose exact phase equals the qualified entry center.
For each trial it:

1. derives scope from the seeded global traversal state;
2. proves the complete cutter sweep lies inside the authenticated precleared
   disk with `PreclearedEntry.certify_first_circle()`;
3. containment-checks the same circle in the design domain;
4. event-certifies it against frozen post-entry stock;
5. depletes stock;
6. appends the exact coverage sweep;
7. advances exactly one traversal cursor and, when causal, one passage; and
8. builds `GenerationState` beginning with approach, plunge, and the circle.

The immutable `InitialCandidateTransaction` binds the entry, candidate,
circle witness, traversal parent/child, and resulting generation-state digest.
Commit independently reproduces byte-identical evidence.

There is no zero-length synthetic link and no exemption from the motion
oracle.

## Continuation and branch switches

Task 12 remains the sole deciding engine for a direct link plus following
circle. Its internal trial accepts an explicit authoritative
`TraversalCursorState`. The existing same-edge consumer supplies
`GenerationState.traversal`; Task 13's branch-switch consumer supplies the
active cursor authenticated by `MatTraversalState`.

Both consumers call the same containment, event certification, depletion,
coverage, passage, witness, and transaction construction path. There is no
second acceptance implementation and no alternate proof order.

After independent candidate commit, `TraversalCommit` verifies that:

- the physical transaction and global traversal share the same candidate;
- exactly one directed cursor advanced;
- every other cursor and visited-side record is byte-identical;
- any causal transit equals the candidate scope and cap transition; and
- the resulting `GenerationState` and `MatTraversalState` digests match the
  committed child.

## Finite search and failure semantics

Candidate families are materialized once per exact span and already use the
invariant furthest-progress, largest-radius, canonical-identity order.
Generation evaluates them in that order and may stop at the first accepted
transaction: no later candidate can outrank it.

Expected proved infeasibility can continue the finite search:

- gouge containment failure;
- exact cap exceedance; and
- a candidate-specific degenerate direct link.

Unsupported or corrupted authority fails immediately:

- unresolved event partition;
- depletion construction or center-limit failure;
- foreign candidate/cursor/policy identity;
- ambiguous or overlapping causal neck topology; and
- any witness or replay inconsistency.

If the finite family is exhausted, the generator raises a named aggregate
failure:

- `EngagementCapInfeasibleError` for a non-neck family proved cap-infeasible;
- `NeckTooTightError` for a causal neck family proved cap-infeasible;
- `NoFeasibleCandidateError` for a mixed exact infeasibility set.

The exception records canonical counts and the active cursor identity in its
message; it does not carry mutable rejected trials or manufacture a success
result.

## Terminal coverage and replay

An edge is terminal only after finite-lattice exhaustion has proved that its
last accepted positive-radius candidate is the terminal candidate in the
selected direction. A component is terminal only when every owned edge cursor
is terminal and every incidence required by its cycles has been visited.

After all components terminate:

1. `CoverageLedger.require_complete()` proves exact reachable residual
   emptiness;
2. the unreachable `D \ M_r` digest remains separately bound;
3. fresh replay rebuilds the input, MAT, neck topology, directed traversal,
   candidate families, stock, and coverage;
4. replay causally re-derives every submitted neck scope;
5. every lateral motion is independently certified before depletion; and
6. every cursor is terminal and the fresh reachable residual is empty.

Traversal exhaustion with nonempty exact residual raises
`IncompletePocketCoverageError`. Coverage cannot terminalize a dropped branch,
and terminal cursors cannot excuse missing coverage.

`GenerationResult` binds the terminal physical state, terminal global
traversal state, ordered traversal commits, coverage certificate, and fresh
replay evidence. It is internal and content-addressed.

## Structural performance contract

Before timing, Task 13 enforces these counts:

- one MAT and neck-inventory construction per generation or fresh replay;
- one native sample index, graph projection, route, and causal-frontier build
  per traversal authority;
- one exact side index per neck;
- constant-time accepted-candidate native cursor-limit lookup;
- one candidate-family materialization per attempted span;
- one isolated stock/coverage fork per candidate trial;
- no repeated full MAT build, cut reconstruction, or `O(k n)` coordinate
  rematching;
- no evaluation of lower-ranked candidates after the first accepted member;
- one independent re-evaluation for each committed winner; and
- one complete fresh replay after terminal coverage.

Exact event certification may still dominate runtime. Profiling begins only
after the complete fixture is GREEN and structural counters confirm these
bounds.

## Acceptance

The committed tractable fixture must prove:

- qualified approach and plunge;
- one entry-launched certified full circle;
- at least one certified direct link;
- at least one later certified full circle;
- ordered operation/witness bijection with no holes;
- at least one causally derived neck scope or an explicit fixture reason that
  no neck is crossed;
- every MAT edge/component cursor terminal;
- exact empty reachable residual with unreachable residual separately bound;
- fresh replay with zero cap violations and byte-identical causal scope; and
- mutation rejection for empty, entry-only, first-circle-only,
  dropped-final-branch, relabelled side, reversed transit, stale traversal
  parent, and certify-after-deplete histories.

The Held–Pfeiffer comparison changes only after these gates pass. Design
strength alone is not end-to-end evidence or a performance result.
