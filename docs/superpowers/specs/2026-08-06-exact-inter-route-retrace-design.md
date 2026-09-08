# Exact inter-route retrace design

**Date:** 2026-08-06

**Status:** approved 2026-08-07

**Scope:** exact-certified adaptive trochoidal MAT Task 13F route-2 boundary

## Problem

Task 13F now reconstructs and commits two exact route families from the
authenticated launch root:

| Route | Family | Accepted physical result |
| --- | ---: | --- |
| `0` | 16 circle candidates | trial 4 link/circle transaction |
| `1` | 36 zero-guide candidates | rank-1 advancing segment to `(5, 1)` |
| `2` | 56 circle candidates | all 56 direct links from `(5, 1)` gouge |

Route `2` is not circle-infeasible. Thirty circles are individually contained.
The failure occurs before their circle proofs because the global route is a
deterministic depth-first edge-discovery ledger, not a continuous cutter walk.
Route `1` exits at a different exact MAT node from route `2`'s entry.

Reversing the accepted route-1 `AdvanceSegmentOperation` is exactly contained.
From its start phase, 10 route-2 links and six complete link/circle pairs are
contained. The missing capability is therefore a physical return through
accepted motion lineage between two graph routes.

The stage must cross this one proved boundary without widening the candidate
lattice, weakening containment, manufacturing traversal progress, or claiming
terminal coverage.

## Decision

Add a distinct exact cut-depth retrace path:

```text
terminal nonincident route
    -> RouteRetraceDecision from the immediately preceding accepted advance
    -> RetraceSegmentOperation with exactly reversed source motion
    -> RouteRetraceTransaction with one swept-prefix witness
    -> RouteRetraceCommit atomically publishes physical return + route activation
    -> materialize the next route family from the restored physical phase
```

The retrace changes physical phase, stock/coverage lineage, and the operation
stream. It does not consume a MAT candidate, advance a MAT cursor, or advance a
neck passage. The activated global route remains a separate authority.

This is deliberately a bounded first retrace contract. It admits the Task 13F
source proved at the current boundary: the immediately preceding, no-neck
`AdvanceSegmentOperation` owned by one `ZeroGuideLinkTransaction`. A future
boundary may justify a nonempty multi-operation retrace plan. This stage does
not guess that more general grammar before evidence requires it.

## Rejected alternatives

### Candidate-lattice widening

All 56 current links leave the reachable design domain. Adding more progress,
radius, or phase samples cannot certify a disconnected physical transition and
would hide the actual chronology defect.

### Traversal-only branch activation

Changing `active_route_index` does not move the cutter. Treating graph
activation as physical transport would let coverage and replay omit a required
motion while still appearing globally advanced.

### Reconstructing transport from MAT coordinates

Re-solving a return path from node coordinates would create new geometry that
was never accepted, lose operation provenance, and require a separate path
selection theorem. Exact accepted operation bytes already provide the required
motion.

### Retract, rapid traverse, and re-plunge

That is a valid later process strategy, but it adds clearance-height,
collision, approach, plunge, and machine-kinematic contracts and abandons the
current continuous cut-depth methodology. It is not a fallback for this stage.

### Generic segment event certification

The generic exact event oracle remained active beyond 10 minutes on this
bounded reverse. The existing swept-prefix theorem certified the same motion
in `0.040943 s` with two exact strata. Routing retrace through the generic
oracle is both the wrong motion-class theorem and a measured performance cliff.

## Exact trigger

Retrace is required only after an accepted candidate terminalizes its active
route and another route remains. Let `completed` be that terminal route step
and `next` the next canonical route step.

```text
completed.exit_node_id != next.entry_node_id
```

is the exact branch-switch predicate. It uses stable MAT node identities, not
distance, tolerance, or reporting coordinates.

The logical transition order is:

1. retain the terminal pre-activation `MatTraversalState`;
2. derive exactly one `activate_next()` state;
3. compare the completed exit and activated entry identities;
4. if they are equal, publish the deterministic activation without motion;
5. if they differ, derive and commit retrace before materializing candidates;
6. expose the activated route to candidate materialization only after the
   retrace physical child exists.

A nonincident transition with no admissible source fails loud before candidate
family construction. It is not translated into candidate exhaustion.

## Causal retrace decision

`RouteRetraceDecision` is a frozen content-addressed record binding:

- completed and activated route indices;
- completed exit-node and activated entry-node identities;
- terminal pre-activation traversal-state digest;
- activated traversal-state digest;
- source `TraversalCommit` digest;
- source `ZeroGuideLinkTransaction` digest;
- source operation index;
- SHA-256 digest of the source operation's canonical bytes; and
- decision schema/strategy identity.

Route-boundary derivation receives the available preimages and proves:

- the source commit is the final physical/global commit in the prefix;
- its global child is the supplied terminal traversal state;
- `activate_next()` reproduces the supplied activated state byte for byte;
- activated route index equals completed route index plus one;
- the completed and activated route steps are nonincident by exact node ID;
- the source transaction is exactly `ZeroGuideLinkTransaction`;
- its only witnessed operation is the final parent operation;
- that operation is exactly `AdvanceSegmentOperation`;
- its traversal decision terminalizes the completed route;
- its scope is `NoNeckScope` with one `FullCapDecision`; and
- its endpoint equals the current physical phase.

Only after these checks does `RouteRetraceDecision.build(...)` seal the typed
indices, node IDs, and digests. No caller may provide a free motion, target
point, or source digest without the corresponding owned preimage.

## Operation grammar

`RetraceSegmentOperation` is a separate exact lateral operation with:

- `ExactSegmentMotion`;
- authenticated `CutZ`;
- `FullCapDecision`; and
- `RouteRetraceDecision`.

Its factory derives the motion; callers do not choose it:

```text
motion.start = source AdvanceSegmentOperation.motion.end
motion.end   = source AdvanceSegmentOperation.motion.start
```

The source and retrace cut depths must be identical. The retrace cap is the
same complete user-cap decision admitted by the no-neck source. The operation
has no `HoldTraversalDecision`, `AdvanceTraversalDecision`, neck scope, circle
fields, or candidate identity.

| Operation | Physical phase | MAT cursor | Passage state |
| --- | --- | --- | --- |
| `LinkSegmentOperation` | changes | holds | holds |
| `CutFullCircleOperation` | returns to its phase | advances | may advance |
| `AdvanceSegmentOperation` | changes | advances | may advance |
| `RetraceSegmentOperation` | changes by exact reversal | holds globally | holds |

The bounded lateral grammar becomes:

```text
entry circle
(
    (hold link, advancing circle)
  | advancing segment
  | route retrace
)*
```

Additional chronology rules make this more restrictive than the regular
expression:

- retrace cannot be the entry lateral operation;
- it must immediately follow the source advancing segment named by its
  decision;
- the source index must equal the retrace index minus one;
- it cannot follow a hold link or full circle in this stage;
- two retraces cannot be adjacent; and
- it does not satisfy the required hold-link predecessor of a later circle.

The operation uses a distinct canonical tag and joins `CanonicalOperation`.

## Generation-state chronology

`GenerationState` treats retrace as a lateral stock/coverage motion while
preserving the last accepted traversal advance.

Construction proves:

- current phase equals the retrace start;
- the immediately preceding operation is the exact source named by the
  decision;
- source and retrace motions are exact reverses;
- source operation digest and ordinal match;
- cut depth and full-cap decision match;
- no hold link is pending;
- resulting phase equals the retrace endpoint;
- local `TraversalCursorState` remains the source advance's terminal result;
- neck passages remain byte-identical; and
- stock and coverage each append the retrace motion at the same ordinal.

Even if the geometric stock boundary and covered union are unchanged because
the corridor was already swept, both lineages append an authenticated witness.
There is no no-op shortcut: chronology is part of the artifact.

## Physical trial and transaction

One `RouteRetraceTransaction` binds:

- authoritative parent-state digest;
- exact `RouteRetraceDecision`;
- one `ReplayLateralWitness` for its `RetraceSegmentOperation`; and
- authoritative child-state digest.

Evaluation uses independent stock and coverage forks in this fixed order:

1. validate parent state, decision, source preimage, and cut-depth authority;
2. construct the exact reverse operation;
3. certify the cutter capsule is a subset of the reachable design domain;
4. build `MotionCertifier` over the unchanged pre-motion stock;
5. call `certify_swept_prefix_segment(...)` exactly once;
6. require the sealed strategy, theorem, and two-stratum witness identities;
7. deplete stock with the existing exact depletion policy;
8. add the identical sweep to coverage;
9. build the child with unchanged local traversal and passages; and
10. build the immutable transaction.

The theorem is applicable because the source operation already removed the
same cutter sweep in the opposite direction. The native audit still verifies
the actual current stock; lineage alone never substitutes for physical proof.

Commit repeats derivation and evaluation from the authoritative parent and
requires byte-identical transaction and child identities. Evaluation failure
leaves parent stock, coverage, phase, operations, traversal, and passages
unchanged.

The existing advancing-segment public function remains. Shared physical trial
steps may be factored into one private swept-prefix segment implementation only
when characterization tests prove the existing advancing witness remains byte
identical. There must be one deciding proof call graph, not copied retrace
logic.

## Cross-axis commit and continuation

`RouteRetraceCommit` atomically binds:

- physical parent and child digests;
- terminal pre-activation `MatTraversalState`;
- activated `MatTraversalState`;
- source `TraversalCommit` digest; and
- exact `RouteRetraceTransaction`.

Its constructor independently re-derives `activate_next()`, the nonincident
predicate, the decision, the one-operation physical suffix, and the unchanged
cursor/passage semantics.

`GenerationContinuation.commits` becomes the closed union
`TraversalCommit | RouteRetraceCommit`. Ordering is authoritative:

- a `TraversalCommit` advances exactly one active cursor;
- an incident completed route activates deterministically without a physical
  commit;
- a nonincident completed route must be followed immediately by exactly one
  `RouteRetraceCommit` before any new `TraversalCommit`; and
- `RouteRetraceCommit` changes the physical digest and publishes exactly one
  route activation, while advancing no MAT cursor.

No second retrace ledger or compatibility field is added. Canonical
continuation bytes preserve the single ordered commit stream.

## Generator integration

Candidate dispatch remains unchanged within one active route. Route-boundary
orchestration gains one exact transition step:

```text
accepted candidate commit
    -> active route terminal?
       -> no: continue candidate search
       -> yes, no next route: activate terminal state
       -> yes, incident next route: activate without motion
       -> yes, nonincident next route: commit exact retrace + activation
    -> materialize next candidate family
```

For Task 13F, route `0 -> 1` remains an incident activation. Route `1 -> 2`
creates exactly one retrace derived from the committed zero-guide advance.
Route-2 candidates are first materialized against the restored phase. The
previous 56-gouge run remains a diagnostic boundary contract; production does
not redispatch those known-disconnected links before retracing.

The stage succeeds when full continuation commits the retrace and advances
beyond the pinned route-2 exhaustion cursor. It does not require terminal
coverage if a later exact boundary appears.

## Fresh replay

Fresh operation replay admits `RetraceSegmentOperation` but never reconstructs
it as a MAT candidate.

Replay must:

- rebuild the source advancing segment from fresh MAT zero-guide proof bytes;
- verify the retrace decision names the immediately preceding source ordinal
  and canonical digest;
- derive the exact reverse motion again;
- skip retrace during circle/zero-guide candidate enumeration;
- preserve candidate indices across the noncandidate operation;
- rerun containment, swept-prefix certification, depletion, and coverage;
- reproduce the live retrace witness byte for byte; and
- reject any operation stream in which retrace is missing, duplicated,
  reordered, forward-directed, or relabelled.

Live `GenerationContinuation` validation independently reproduces route
activation and the cross-axis commit. The later terminal certificate stage
still owns fresh reconstruction of the complete global route and exact empty
reachable residual. This stage does not manufacture a successful
`ReplayCertificate` before those gates exist.

## Error model

The retrace path fails loud with named errors:

- `InvalidRouteRetraceDecisionError` for malformed or cross-wired causal
  identity;
- `InvalidRetraceSegmentOperationError` for foreign geometry, depth, cap, or
  source reversal;
- `UnsupportedRouteRetraceError` when a nonincident transition has no admitted
  source operation;
- `InvalidRouteRetraceTransactionError` for malformed physical evidence;
- `StaleRouteRetraceTransactionError` when authoritative parent state changed;
- `InvalidRouteRetraceCommitError` for physical/global activation mismatch;
  and
- `ReplayRouteRetraceError` for fresh grammar, lineage, or witness divergence.

Containment, cap, unresolved-theorem, stock, and coverage errors retain their
existing named types. None is caught as ordinary candidate infeasibility.

## File responsibilities

Implementation preserves focused ownership:

- `operation.py`: canonical decision and operation grammar only;
- `retrace_segment_trial.py`: shared ordered swept-prefix physical trial;
- `retrace_transaction.py`: physical evaluate/commit and transaction identity;
- `generation_state.py`: operation chronology and physical-state invariants;
- `generator.py`: exact route trigger, cross-axis commit, and ordered
  continuation integration;
- `replay_trace.py`: closed lateral witness grammar;
- `replay.py`: small dispatch points into retrace replay, not another retrace
  implementation;
- `errors.py`: named public failure modes; and
- new focused retrace test modules: operation/state, transaction/generator,
  and fresh replay.

Replay remains orchestration; the deciding retrace proof stays in its focused
trial module. Every new test has a contextual Google-style docstring explaining
the geometric or causal contract it protects.

## TDD gates

### Decision and operation

- exact nonincident route identities derive one stable decision;
- input reversal and fresh reconstruction preserve decision bytes;
- incident routes produce no retrace decision;
- wrong route index, node, traversal digest, source commit, transaction,
  operation ordinal, or operation digest fails;
- forward, shortened, lengthened, offset, and zero-length motions fail; and
- neck-owning, circle, hold-link, nonfinal, and foreign-depth sources fail.

### State and transaction

- the real route-1 reverse certifies with exactly one two-stratum
  swept-prefix audit;
- containment and certification precede depletion and coverage;
- child phase returns to the source start;
- local traversal and passages remain byte-identical;
- stock and coverage lineages each append one matching motion;
- parent state remains unchanged after every failed trial;
- stale parent, changed policy, wrong cap, wrong witness theorem, wrong stock
  boundary, and cross-wired lineage fail; and
- independent commit reproduces transaction and child bytes.

### Generator and continuation

- route `0 -> 1` remains an incident zero-motion activation;
- route `1 -> 2` emits exactly one `RouteRetraceCommit` before materialization;
- no route-2 candidate is evaluated against phase `(5, 1)`;
- route-2 family identity and invariant ordering remain unchanged;
- at least one route-2 pair commits from the restored phase;
- missing, duplicate, delayed, or reordered retrace commits fail;
- retrace cannot advance the active route-2 cursor or a neck passage; and
- continuation advances beyond cursor
  `def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0`.

### Fresh replay

- fresh MAT proof reconstructs the source advance before retrace validation;
- replay reproduces the live retrace witness byte for byte;
- retrace is never counted as a circle or zero-guide candidate;
- altered source ordinal/digest, direction, route activation, cap, theorem,
  stock owner, operation order, or omitted retrace fails; and
- the expected later nonterminal or residual gate remains explicit until the
  complete path exists.

### Publication

- focused RED/GREEN tests with `-n auto --testmon`;
- affected adaptive tests;
- complete adaptive suite with `-n auto`;
- Ruff and formatting;
- strict mypy;
- strict MkDocs;
- diff check;
- independent Critical/Important review; and
- author/committer, remote SHA, and clean-worktree verification.

## Performance contract

One Task 13F retrace transition performs:

- one exact node-identity comparison;
- one adjacent source-operation lookup and digest verification;
- one exact segment-containment proof;
- one native swept-prefix audit with exactly two strata;
- one stock depletion;
- one coverage update; and
- one independent commit replay of the same bounded trial.

It performs no generic event partition, coordinate rematch, candidate-family
probe before retrace, MAT rebuild, graph search, or retract/re-plunge planning.
Structural counts are acceptance gates; wall time is recorded as bounded
regression evidence rather than a flaky pass/fail threshold.

Matched Held–Pfeiffer performance comparison remains downstream of terminal
generation. The maturity statement remains:

```text
stronger exact bounded proof contract;
weaker end-to-end and measured-performance evidence
```

## Documentation contract

At the coherent implementation checkpoint, update:

- `docs/segment_site_mat.md`, including the Held–Pfeiffer comparison;
- `docs/continuous_engagement.md`;
- the causal traversal and zero-guide continuation plans; and
- this design's implementation status and measured structural evidence.

Document the exact trigger, causal source, operation grammar, transaction,
fresh replay, performance counts, current limitation, and next boundary while
the evidence is live. Tri-dexel removal/thermal replay and matched performance
remain later tasks; neither substitutes for the exact retrace certificate.

## Acceptance

This stage is complete only when the production generator:

1. reaches the established route-1 terminal child unchanged;
2. detects the route-2 branch switch by exact node identity;
3. derives exactly one retrace from the immediately preceding accepted
   zero-guide operation;
4. independently commits its containment, swept-prefix, stock, and coverage
   evidence;
5. publishes the activated route and physical return under one ordered
   continuation identity;
6. materializes and commits route-2 work from the restored phase;
7. reproduces the retrace under fresh operation replay; and
8. stops fail-closed at the next exact boundary or returns only after the
   already specified terminal traversal and residual gates are genuinely met.
