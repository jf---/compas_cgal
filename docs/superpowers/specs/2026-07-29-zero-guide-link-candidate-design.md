# Exact zero-guide link candidate

**Date:** 2026-07-29  
**Status:** approved direction  
**Scope:** exact-certified adaptive trochoidal MAT Task 13F / terminal traversal

## Problem

The radius-1 Task 13F fixture now selects and independently commits its fourth
circle candidate. That commit terminalizes route 0. Route 1 is a straight
width-2 arm whose exact MAT squared-clearance function is identically the
radius-1 tool square. Its MATHSM guide radius is therefore zero over the whole
edge.

The current candidate grammar admits only positive-radius MATHSM circles. It
correctly materializes no circle on route 1 and fails loud with an empty
family. That is not a terminal route: moving the cutter centre along the arm
is required material removal.

The missing operation is an exact link-only MAT advance. It must not be
represented by a zero-radius circle, inferred from reporting floats, or
replaced by traversal-only terminalization.

## Decision

Add a second, structurally distinct candidate and transaction variant for an
exact constant tool-clearance MAT run.

The existing path remains:

```text
positive guide cell
    -> MiddleCurveCandidate
    -> hold LinkSegmentOperation
    -> CutFullCircleOperation owns traversal advance
    -> CandidateTransaction with two lateral witnesses
```

The new path is:

```text
exact zero-guide run
    -> ZeroGuideLinkCandidate
    -> AdvanceSegmentOperation owns traversal advance
    -> ZeroGuideLinkTransaction with one lateral witness
```

Both paths use the same exact segment containment, TEA, depletion, coverage,
cap, neck, traversal, state-identity, and fresh-replay authorities. They differ
only where their physical motion and chronology genuinely differ.

## Rejected alternatives

### Optional circle fields

Generalizing `CandidateTransaction` so that its circle and circle witness are
optional would give one type two chronologies. Consumers would need presence
tests to learn whether a link holds or advances traversal. That weakens the
current pair invariant and is rejected.

### Traversal-only exhaustion

Marking the zero-guide edge terminal without a lateral operation would certify
graph visitation while omitting required removal and coverage. Traversal
terminality is not machining completeness and cannot authorize this shortcut.

### Zero-radius circle

The exact circle grammar requires a nonzero phase vector. A degenerate circle
would have no meaningful full-circle sweep, orientation, event partition, or
MATHSM radius identity. It is rejected rather than special-cased.

### Tolerance classification

Reporting `guide_radius == 0.0`, near-zero tests, or policy comparison cannot
authorize the new path. The decision comes only from the exact MAT
squared-clearance polynomial and exact tool-radius square.

## Native zero-guide theorem

### Predicate

For one clipped MAT edge with squared-clearance polynomial
\(C_e(t)\) and exact tool-radius square \(r^2\), a zero-guide run exists iff

```text
C_e(t) - r²
```

is the zero polynomial over the complete certified edge interval.

This is an identity test on exact rational coefficients, not an endpoint or
sample test. It proves zero guide radius throughout the edge, including
between reporting stations.

### Native record

A new native zero-guide module emits one canonical record per proved edge. The
record binds:

- schema and strategy version;
- the complete MAT-certificate digest;
- centre-domain digest;
- exact edge identity;
- defining generator-site identities;
- lower and upper exact parameter-root identities;
- exact squared-clearance coefficients;
- exact tool-radius square; and
- the identity verdict.

The inventory is emitted only after replaying the submitted MAT certificate
against the retained exact owner. Inventory construction visits every
certificate edge exactly once and emits records only for zero-polynomial
differences. Nonconstant and unequal constant profiles are ordinary nonmembers,
not errors. The explicit record verifier fails with distinct named native
exceptions for duplicate or missing edge ownership and for any forged
nonconstant, unequal, or otherwise mismatched record.

`SegmentSiteMatBundle2` constructs the complete canonical inventory once.
The nanobind owner exposes immutable `(edge_id, certificate_bytes)` pairs
without adding another field to the fixed 20-position numeric projection.

### Python ownership

`MedialAxis` projects the native inventory into an immutable
`MatZeroGuideInventory`. Each `MatZeroGuideRun` retains:

- its owned `MatEdgeId`;
- its native certificate bytes; and
- the MAT-certificate digest to which those bytes belong.

Construction validates uniqueness, exact edge ownership, and byte equality
with the retained native owner. `run_by_edge_id` is graph-lifetime authority;
candidate materialization performs an indexed lookup and never reclassifies
the edge.

## Candidate model

`ZeroGuideLinkCandidate` is a frozen content-addressed record with:

- `MatZeroGuideRun`;
- `CandidatePolicy`;
- exact spatial progress and contributing spatial refinement levels;
- exact target `Point2[WorldXY]`;
- native cursor-limit identity;
- neck scope and effective-cap decision;
- `AdvanceTraversalDecision`; and
- SHA-256 identity over the complete canonical record.

It deliberately has no guide radius, radius levels, phase index, circle
orientation, material side, or generator-site choice.

Spatial cells reuse the existing exact middle-curve spatial lattice and
analytic span interpolation. A zero-guide candidate is emitted only when the
span edge owns a native `MatZeroGuideRun`. The terminal native endpoint remains
an ordinary terminal candidate because a nonzero segment can reach it.

The common invariant ordering remains furthest progress, then squared guide
radius, then canonical identity. The zero-guide variant contributes exact
squared radius zero. Exact edge classification prevents positive-guide and
zero-guide variants from competing on the same span.

## Operation grammar

`AdvanceSegmentOperation` is a new exact operation type with:

- `ExactSegmentMotion`;
- `CutZ`;
- neck scope;
- effective-cap decision; and
- `AdvanceTraversalDecision`.

It is not a second calling convention for `LinkSegmentOperation`.

| Operation | Traversal semantics | Required successor |
| --- | --- | --- |
| `LinkSegmentOperation` | hold | immediately paired full circle |
| `CutFullCircleOperation` | advance | none |
| `AdvanceSegmentOperation` | advance | none |

The operation is canonicalized under a distinct tag. Operation-stream
continuity advances the accepted phase to the segment endpoint.

## Transaction and state chronology

`ZeroGuideLinkTransaction` binds:

- authoritative parent-state digest;
- exact zero-guide candidate;
- one `ReplayLateralWitness` for its `AdvanceSegmentOperation`;
- resulting traversal cursor;
- optional advanced neck passage;
- authoritative child-state digest.

Its constructor proves:

- candidate, operation, scope, cap, and traversal identity agree;
- segment endpoint equals the candidate target;
- the motion witness observes the pre-depletion stock lineage;
- depletion and coverage share the same motion and parent lineages;
- the resulting phase point is the segment endpoint; and
- passage advancement agrees with the candidate cap decision.

`CandidateTransaction` remains the exact two-witness link/circle record.
`AcceptedCandidateTransaction` is a closed union used by orchestration.

`GenerationState` admits this lateral chronology:

```text
entry circle
((hold link, advancing circle) | advancing segment)*
```

It rejects a dangling hold link, a hold link followed by an advancing segment,
two different operations claiming one traversal advance, or passage state
advanced by a hold link.

The state may end at either an accepted full circle or an advancing segment.
Its `phase_point` and traversal cursor must match that final advancing
operation.

## Evaluation and commit

`CandidateEvaluator` gains two explicit methods:

- `evaluate_zero_guide_from_cursor(...)`; and
- `commit_zero_guide(...)`.

The methods use the existing short-lived proof authority but return the new
transaction type. The current circle methods and canonical bytes do not
change.

Evaluation order for one zero-guide trial is:

1. validate the physical parent and authenticated traversal cursor;
2. validate the candidate's internally bound run and cursor identities;
3. derive the effective cap and causal passage result;
4. construct the segment from `state.phase_point` to candidate target;
5. prove containment;
6. certify TEA against the pre-depletion stock;
7. deplete exact stock;
8. add exact coverage;
9. build the child `GenerationState`; and
10. build the content-addressed transaction.

Commit independently repeats those steps and requires byte-identical
transaction and child identities.

## Generator dispatch

`TraversalCandidate` and `AcceptedCandidateTransaction` are closed type
unions.

For every active span:

- if its edge owns a `MatZeroGuideRun`, materialize only zero-guide link
  candidates;
- otherwise materialize only existing positive-guide circle candidates.

Before dispatch, orchestration validates the run owner against the traversal
axis rebuilt from the same `InputIdentity`. The choice is a graph-lifetime
exact lookup. It is not a caught exception, empty-family retry, or
reporting-value fallback.

The combined family is still validated for unique identities and invariant
order. Selection continues only after the same named local failures:
containment, cap exceedance, or degenerate segment. Any native proof gap aborts
immediately. The first accepted transaction is committed through its matching
typed commit method.

## Fresh replay

Fresh replay rebuilds the MAT and its zero-guide inventory from
`InputIdentity`.

It reconstructs:

- circle candidates from `CutFullCircleOperation`;
- zero-guide candidates from `AdvanceSegmentOperation`; and
- hold links only as the immediately preceding member of a circle pair.

An advancing segment must match exactly one freshly enumerated zero-guide
candidate. Its native zero-guide certificate bytes must match the rebuilt MAT
owner. Pairing maps each operation index to the candidate that owns its cap and
traversal decision.

Fresh physical replay uses the same segment containment, motion certification,
depletion, and coverage sequence as live evaluation. The final terminal gate
continues to require complete traversal and exact empty reachable residual.

## Error model

The new path fails loud with named errors:

- `InvalidZeroGuideCertificateError`;
- `UncertifiedZeroGuideEdgeError`;
- `InvalidZeroGuideCandidateError`;
- `InvalidAdvanceSegmentOperationError`;
- `InvalidZeroGuideTransactionError`;
- `StaleZeroGuideTransactionError`; and
- `ReplayZeroGuideCandidateError`.

No new error is translated to candidate-local infeasibility unless it proves
containment failure, cap exceedance, or degenerate motion under the existing
selection policy.

## TDD gates

### Native proof

- radius-1 L emits exactly its constant tool-clearance line runs;
- radius-0.5 L emits no false zero-guide record for positive-guide lines;
- record order and bytes repeat and survive input reversal;
- MAT-certificate, edge, site, endpoint, coefficient, tool-radius, and verdict
  mutations fail;
- endpoint-only equality and near-equality do not produce a run.

### Candidate

- forward and reverse spatial lattices are deterministic and identity-distinct;
- only owned zero-guide edges enumerate;
- target, progress, cursor, run, scope, and cap cross-wiring fail;
- no radius, phase, or circle field exists;
- terminal endpoint and derived interior cursor retain exact lineage.

### Transaction

- one real segment is contained, TEA-certified, depleted, and covered in order;
- failure leaves parent stock, coverage, traversal, and passage unchanged;
- stale parent, foreign run, altered cap, wrong endpoint, wrong traversal,
  certify-after-deplete, and cross-wired witness mutations fail;
- independent commit reproduces bytes and child digest.

### Generator and replay

- Task 13F route 0 remains the same fourth-circle winner;
- route 1 materializes a nonempty zero-guide family;
- invariant order stops at its first accepted link-only winner;
- no lower-ranked dispatch occurs;
- fresh MAT reconstruction uniquely recovers the advancing segment candidate;
- changing the MAT certificate or relabelling a hold link as an advancing
  segment fails;
- full continuation advances beyond the current `attempts=0` boundary.

### Regression and publication

- affected `--testmon`;
- complete adaptive suite;
- Ruff;
- strict mypy;
- strict MkDocs;
- diff check;
- source/replay mutation slice; and
- bounded structural count/timing for Task 13F after a complete path exists.

## Performance contract

Zero-guide classification runs once during native MAT construction and is
linear in edge count times the already bounded clearance-polynomial degree.
Materialization uses one immutable edge lookup and the existing spatial
lattice. It performs no repeated MAT build, coordinate rematch, trial fallback,
or radius/phase enumeration.

Runtime comparison with Held–Pfeiffer remains prohibited until Task 13F
reaches a complete terminal path and the bounded planner workload is measured.

## Documentation contract

At each coherent implementation checkpoint, update:

- `docs/segment_site_mat.md`;
- `docs/continuous_engagement.md`;
- the Task 13F / causal traversal plan; and
- the Held–Pfeiffer comparison table.

The maturity statement remains:

```text
stronger exact bounded proof contract;
weaker end-to-end and measured-performance evidence
```

until terminal traversal, exact residual emptiness, fresh replay, and matched
performance gates are all complete.
