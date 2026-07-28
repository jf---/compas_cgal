# Atomic candidate transaction design

## Scope

Task 12 evaluates one finite-lattice `MiddleCurveCandidate` as a joint
link-and-circle proposal. It either returns immutable, content-addressed
acceptance evidence or raises the exact named proof failure. Evaluation never
changes authoritative stock, coverage, traversal, passage, or operation state.

Task 13 remains responsible for choosing the active MAT branch and deriving
causal neck scope from global separating-cut traversal. Task 12 validates the
scope and transition supplied by that traversal authority; it does not infer
them from geometry.

## State ownership

`GenerationState` is an immutable authoritative snapshot. It owns independent
stock and coverage snapshots plus value state that changes at the same accepted
candidate boundary:

- current phase point;
- current component, edge, branch, and cursor identity;
- canonical oriented-neck passage states;
- accepted operation prefix.

Factories clone mutable native owners, validate cursor and passage uniqueness,
and compute a canonical state digest. Read access returns immutable evidence or
fresh forks, never the owned mutable native objects.

State digesting binds the stock boundary and depletion lineage, coverage
certificate, phase point, traversal cursor, passage states, and operation
prefix. A transaction evaluated against one digest cannot commit against
another.

## Evaluation

`CandidateEvaluator` is short-lived and owns only invariant evaluation inputs:
reachable-domain containment authority, tool radius, user cap, neck policy,
depletion policy, cut depth, and material side.

For each candidate it:

1. validates that the candidate advances the state's exact active cursor;
2. independently reconstructs the full or oriented-neck effective cap;
3. constructs the direct segment from the current phase point to the candidate
   circle phase;
4. forks stock and coverage;
5. containment-checks and event-certifies the segment against pre-link stock;
6. depletes the segment and appends its exact coverage sweep;
7. containment-checks and event-certifies the circle against post-link stock;
8. depletes the circle and appends its exact coverage sweep;
9. advances the copied traversal and passage value state;
10. returns a frozen `CandidateTransaction` binding every operation, witness,
    parent digest, and resulting state digest.

The certify-before-deplete order is structural. Each motion witness must bind
the exact stock lineage visible immediately before its own depletion.

## Commit

The transaction contains no mutable stock or coverage payload. Explicit commit
first requires an exact parent-state digest match, then re-evaluates the
selected candidate from the authoritative state. The recomputed transaction
must be byte-identical to the submitted transaction before the newly evaluated
forks become the next `GenerationState`.

This adds one bounded evaluation for the single winner, not for rejected
candidates. It is preferred to embedding mutable native owners in a nominally
frozen transaction, because the latter would make immutability cosmetic and
permit post-acceptance proof tampering. Runtime measurement follows after the
complete path exists.

## Rejection and errors

Success returns `CandidateTransaction`. Existing named failures propagate for
gouge containment, cap exceedance, unresolved event partitions, depletion
proof failure, and coverage proof failure. Task 12 adds distinct named errors
for:

- structurally invalid generation state;
- candidate traversal or passage mismatch;
- malformed transaction evidence;
- stale parent state;
- empty or cross-parent winner selection.

There is no rejection enum, result wrapper, rollback path, fallback, or partial
commit.

## Deterministic winner

Winner selection is a pure function over accepted transactions. All
transactions must bind one parent state and one `CandidatePolicy`. Selection
reuses the policy's invariant order: furthest exact progress, then largest
exact guide radius, then canonical candidate identity. Input order cannot
change the winner.

## Verification

Contract tests use the real L-pocket MAT candidate lattice and exact native
stock consumers. They prove:

- link pass followed by circle failure leaves every authoritative digest and
  lineage unchanged;
- link failure, unresolved event, and gouge failure are loud and atomic;
- successful evaluation is deterministic and still non-mutating;
- explicit commit advances stock, coverage, cursor, phase, passage, and
  operation state together;
- a stale-parent transaction cannot commit;
- deplete-before-certify and passage-before-commit mutations are rejected;
- winner selection is invariant under input permutation.

Every test carries a docstring explaining the geometric and proof-boundary
reason for the contract.
