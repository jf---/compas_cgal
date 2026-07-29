# Exact Continuous Engagement Certification

Continuous tool-engagement certification decides whether every cutter station
along one declared segment or full circle respects an effective tool
engagement angle (TEA) cap. It is a proof boundary, not a dense station audit:
acceptance follows from a complete exact event partition and sign-invariant
cells.

!!! warning "Current maturity"

    The native segment and full-circle oracles and the typed
    `MotionCertifier` consumer are implemented for the bounded Phase-1
    contracts. Both motion families return digest-bound event traces, and the
    Python boundary distinguishes certified, cap-exceeded, and unresolved
    outcomes. Task 11A now independently reconstructs a two-circle no-neck
    L-pocket prefix and its derived-cursor link. It validates the complete
    link/circle relation before stock mutation, then certifies the first circle,
    link, and following nonuniform circle against their frozen pre-cut states
    before each depletion and coverage mutation. Every returned containment,
    motion, depletion, and sweep proof is now retained in a canonical,
    cross-validated per-operation bundle; two independent rebuilds produce the
    same trace bytes. Replay next fails at the explicit
    nonterminal-traversal boundary because other MAT edges remain untouched.
    Fresh exact neck inventory plus submitted oriented-owner/state/cap replay
    are also gated, including the production `80 degrees` second-passage
    rejection. Task 13A now retains the exact neck loci and all three certified
    sides of each production plateau neck; global traversal has not yet used
    that side history to derive scope independently. Complete
    traversal/coverage, replay/artifact certification, arbitrary-pocket
    evidence, and matched Held–Pfeiffer performance remain incomplete.

## Authority hierarchy

The records have deliberately different jobs:

| Layer | Circle authority | Segment authority | Consumer role |
| --- | --- | --- | --- |
| Reconstructed partition | `EventPartitionCertificate2` | `SegmentEventPartition2`, including its nested `EventPartitionCertificate2` | Proves that every declared projection, algebraic root, cell, fibre, overlap, and seam can be reconstructed |
| Deciding data | Exact uniform-sweep disposition, dedicated violation witness, or `FullCircleCellAuthority2` with reconstructed stationary strata | Ordered segment strata, active branches, and pair orientation/cap dispositions | Authorizes the exact verdict |
| Trace | `EventTrace2` v2 | `EventTrace2` v2 | Orders canonical events and binds the deciding record by SHA-256 |
| Typed witness | `MotionWitness` | `MotionWitness` | Canonically binds operation ordinal/kind, exact motion, caps, stock lineage, strategy, and trace digest under one content identity |

A generic event certificate is topology authority, not generally decision
authority for either motion family. Segment cell classification consumes the
larger `SegmentEventPartition2`, especially the ordered active-branch inventory
and each material run's pair-cap disposition. A nonuniform circle likewise
reconstructs one stationary stratum per exact parameter cell and records those
dispositions in a separate `FullCircleCellAuthority2`.

`EventTrace2` v2 therefore carries two views:

- `partition` is the generic certificate used by common trace consumers;
- `decision_authority_bytes` and `decision_authority_digest` identify the
  complete reconstructed record that authorized the verdict.

For uniform circles and the dedicated line/vertex violation path, the
authority bytes equal `partition.canonical_bytes`. For nonuniform cell
classification and for segments, they deliberately differ from the nested
generic certificate. The trace canonical record includes both the generic
partition digest and the decision-authority digest. A
`MotionWitness.event_trace_digest` is consequently transitive to the full
deciding record.

!!! danger "A digest beside a proof is not proof ownership"

    Recording a full segment or circle-cell digest as metadata without
    including it in the trace identity would allow the same trace digest to be
    relabelled onto a different deciding partition. The authority digest must
    be inside the canonical trace record.

The planning `InputIdentity` separately binds
`motion-certificate-schema-v1` and the native
`event-exact-motion-oracle-v3` component. Changing the witness schema or
oracle implementation family therefore changes the planning root before
replay begins.

## Exact segment partition

`SegmentEventSource2` lifts the supplied binary64 endpoints, tool radius, and
cap surrogate exactly once. The motion parameter is the closed rational
interval `[0, 1]`. The native producer derives and reconstructs:

- two rim charts for each reachable boundary support;
- support tangencies and support overlaps;
- trimmed-vertex passages and endpoint-order resultants;
- pair-orientation boundaries that distinguish zero from `pi`;
- pair-cap equalities on the correct CCW branch;
- every isolated algebraic root, open cell, and zero-dimensional fibre; and
- original-equation, orientation, trim, incidence, and multiplicity checks at
  each fibre.

Within an open cell, the ordered active rim intersections cannot change.
One exact reference point determines whether the cyclic reference interval is
material or clear. The oracle then pairs the alternating intersections into
material runs. Every run must have one reconstructed pair disposition:

- `below-cap` or `equal-cap` is safe;
- `above-cap` proves cap exceedance;
- a missing/unknown disposition is unresolved.

A rim with no intersections is not automatically clear. Its exact reference
classification distinguishes a fully clear rim from a fully material rim.
The latter proves cap exceedance.

At event fibres, all algebraic roots must be evaluated and their original
geometry, orientation, and trims rechecked. An unresolved overlap or failed
recheck makes the whole motion unresolved unless another cell already proves
cap exceedance. A proof of physical exceedance takes precedence over proof
incompleteness.

## Exact full-circle partition

A full circle is defined by its center, nonzero phase vector, and orientation;
no rounded guide radius replaces the phase vector. Four exact center-circle
charts own the complete cycle and its seams. Canonical event order follows
oriented motion order and then event identity, so CW/CCW results do not depend
on chart insertion order.

Uniform empty/material cases still produce replayable four-seam partitions.
Nonuniform cases reconstruct the complete line/circle pullback event set,
including tangent, overlap, endpoint-order, cap, and seam fibres. Unsupported
or unreconstructed degeneracy is unresolved, never sampled into acceptance.

### Nonuniform full-circle cell authority

The event partition proves where the engagement combinatorics may change. It
does not prove the disposition of the open cells between those events.
`event-exact-motion-oracle-v3` closes that second proof obligation for
same-support material runs.

Each global rational cell witness is assigned to its unique owned quarter
chart. The rational Pythagorean parametrization

```text
u(t) = ((1 - t²) / (1 + t²), 2t / (1 + t²))
```

then turns the exact binary64 circle center and phase vector into an exact
rational cutter station. The station uses the same branch construction,
cyclic material-run pairing, orientation test, and pair-cap classifier as the
segment oracle. This is one shared stationary theorem, not a second circle-only
interpretation of engagement.

For every partial material run, the authority additionally requires both
boundary branches to lie on one line or circle support and requires the
chart-specific cap projection that makes its sign invariant over the complete
cell. An equal-cap run requires the corresponding identically-equal overlap
record. Cross-support pairs and unsupported overlap families remain
`unresolved`; an exact witness is never promoted beyond the theorem admitted
by the authority.

`FullCircleCellAuthority2` commits to:

- the reconstructed generic partition;
- the canonical stock-feature inventory;
- exact center, phase, tool-radius, and cap inputs;
- every rational station source and complete stationary cell stratum;
- each reference-side and material-run disposition; and
- the exact cap-projection identities used by each decided partial run.

Its SHA-256 digest is therefore deliberately distinct from the generic
partition digest and is included inside `EventTrace2`. The strategy is
versioned `full-circle-cell-strata-exact-v1`.

!!! warning "A verified event partition is not an engagement verdict"

    Reconstructing every root, cell, fibre, and seam proves that the parameter
    decomposition is complete. It does not classify those strata as clear,
    below-cap, or above-cap. Partition verification and cell-disposition
    authority must remain separate; otherwise a structurally valid partition
    can be mistaken for a motion certificate.

The whole-motion cell authority uses fail-closed verdict precedence:
an unresolved cell dominates locally material or above-cap cells. A separate
dedicated line/vertex witness may still prove existential cap exceedance
without completing every cell; that is a different proof object with a
different contract.

The first implementation accidentally rebuilt the partition three times,
reconstructed every stationary cell twice, and re-extracted the stock boundary
once per cell. Hoisting verifiedness and boundary records across the cell loop
reduced a warm Apple M1 Max Release median from 252 ms to 128 ms for the
10-cell concentric certificate (10 measured repetitions), and from 572 ms to
301 ms for a 43-cell cross-support audit (5 repetitions): 1.96x and 1.90x,
respectively. These are bounded development timings for one oracle call, not a
Held–Pfeiffer planner comparison.

### Exact clear sweeps in nonempty stock

The first Task 11A replay circle exposed a proof case that global stock
emptiness cannot represent. Its complete cutter sweep lies inside the
qualified entry void, while material correctly remains elsewhere in the
pocket. The boundary pullback partition reconstructed successfully, but the
old uniform classifier recognized `clear` only when the stock had no boundary
records. It therefore returned `unresolved` for a physically and exactly
zero-engagement motion.

The retained uniform path in `event-exact-motion-oracle-v3` derives a
sqrt-capable exact-region
view of the current stock. Every linear/circular support and every one-root
endpoint is lifted without a double conversion:

```text
x = a0 + a1 sqrt(root)
```

It then constructs the same exact full-circle sweep used by containment and
tests the regularized intersection

```text
full_circle_sweep(c, phase, r) intersect current_stock = empty.
```

An empty result certifies a uniformly clear rim even when remote material
remains. A nonempty result proves nothing about the cap and continues through
the event-partition path; it is not relabelled as material or exceeded. The
uniform-circle strategy is versioned
`full-circle-uniform-event-exact-v2`.

The regression uses only dyadic geometry. A phase vector `(3/32, 4/32)` has
guide radius `5/32`; with tool radius `1/2`, its complete sweep fits exactly
inside an entry disk of radius
`2(5/32) + 1/2 = 13/16` centred at the phase point. The square stock remains
nonempty, so the fixture kills any return to an empty-stock shortcut while
also gating exact boundary tangency.

### First state-dependent continuation boundary

Fresh replay now reaches beyond the uniformly clear entry recut. From the
first accepted radius-`1/8` circle it reconstructs the exact intermediate
`DerivedCandidateCursor`, enumerates a terminal radius-`1/8`,
phase-index-`3` circle owned by the point generator at `(2, 2)`, and derives
the direct link between their phase points. The link is exactly contained and
event-certified against the post-first-circle stock, then depleted and added
to coverage.

The following nonuniform circle now receives a distinct cell-decision
authority, certifies against the frozen post-link stock, and only then depletes
stock and adds its exact coverage sweep. The chronology regression requires
all three actions. Replay subsequently raises `ReplayTraversalError` because
other MAT edges remain nonterminal; it still emits no partial certificate.
Before that rejection, `FreshReplayTrace` retains the entry depletion,
first-circle-in-entry proof, initial coverage state, and each circle/link proof
bundle plus terminal snapshots. The result closes the first post-link circle
disposition and its evidence handoff without claiming empty residual or
complete-path certification.

Fresh oriented-neck replay is now a separate consumer of the same exact oracle.
It rebuilds the neck inventory, seeds independent forward/reverse passage
state, resolves the recorded owner, reconstructs the policy cap from the
current state, and advances only after one unique circle candidate matches. A
real first-passage circle reaches the certifier at `90 degrees`, not the
`120 degrees` user cap.

!!! warning "Validated scope is not yet derived scope"

    The candidate lattice currently accepts scope as an identity input. The
    same exact motion/traversal can be rebuilt with no-neck/full-cap identity
    or with oriented-neck/restricted-cap identity. Replay validates the latter
    internally but cannot infer which separator was crossed from motion
    geometry or one edge cursor. Task 13 must carry global certified-side
    history and assign the causal scope before final replay certification.

!!! warning "A legal passage does not certify its circle"

    The next radius-`1/8`, phase-index-`3` circle advances the same neck passage
    legally, but exceeds its reconstructed `80 degrees` cap against post-link
    stock. That negative is an integration gate. An equal-cap fixture separately
    proves the two state transitions. Task 12 now rejects that complete joint
    trial; Task 13 must continue the feasible-candidate search instead of
    weakening the oracle or conflating passage legality with motion acceptance.

### Task 12 atomic consumer

The event-exact oracle now has a production candidate-level consumer.
`CandidateEvaluator` evaluates the direct phase link and its proposed full
circle on private `Stock2Area` and `CoverageLedger` forks. For both motions it
preserves the only valid causal order:

```text
containment -> MotionCertifier.certify -> deplete -> add_sweep
```

The second certification therefore observes the link-depleted trial stock,
while the parent `GenerationState` remains unchanged. A named containment,
cap-exceeded, unresolved-event, depletion, or coverage failure destroys the
trial by ordinary scope exit; no partial link can escape.

Successful evaluation returns immutable `CandidateTransaction` evidence, not
mutable stock. Commit rejects a stale parent digest, independently repeats the
winning evaluation, and requires identical canonical bytes before returning
the new state. The state factory separately replays every oriented passage
transition from `UNVISITED`, pairs every historical link and circle by
scope/cap/cursor, and requires all lateral `CutZ` values to equal the qualified
entry plane. The evaluator binds candidate, depletion, cap, and cut-direction
policy before extending that state. This closes planar-consumer gaps that TEA
alone cannot observe: unearned neck state, mixed policies, or a geometrically
valid motion emitted at a foreign depth.

The focused L-pocket gate proves both positive and negative composition. One
real no-neck candidate commits its link and circle together. A real oriented
second-passage link certifies at the reconstructed `80 degrees` cap, but its
circle exceeds that cap and the complete joint trial is discarded. A separate
larger lattice circle passes its direct link and fails exact containment.

Task 13 remains the sole authority for assigning causal neck scope, choosing
the active global MAT branch, and deciding which candidate failures constitute
normal search exhaustion. Task 12 validates a supplied finite candidate and
state transition; it does not infer global traversal facts.

!!! danger "Exact-number backend definitions are an ABI contract"

    `_continuous_tea_2` compiles with its target-local CORE/Boost backend
    definitions. Linking a sweep object compiled under a different exact
    backend can give the same `ReachKernel` spelling different effective
    number types across a C++ boundary. The oracle therefore compiles the
    authoritative `exact_sweep_2.cpp` source inside its own target. Reuse the
    deciding source, but never share compiled exact-kernel objects until every
    consumer inherits one verified interface definition set.

## Exact reach pruning before pair resultants

Pair-event generation is quadratic in the number of scheduled boundary
features. If `n` features survive and the fixed event grammar emits `K`
resultants per unordered pair, the scheduled count is:

```text
K n (n - 1) / 2
```

Before constructing those resultants, the segment producer computes the exact
rational minimum squared distance from each circular support center to the
closed motion segment. A circular support of radius `R` is removed only when
that minimum distance `d_min` proves:

```text
d_min > R + r
```

where `r` is the tool radius. The implementation uses the squared exact
equivalent and retains every line support. It is conservative for nested
circle cases: ambiguity retains work; only proven non-reach removes it. A
structural assertion verifies that removing features strictly reduces the
scheduled pair-resultant count.

During the Task-6 recovery, the same 13-test segment-oracle command exceeded
176 seconds and was interrupted with two cases unfinished. After tangent-root
normalization and exact support pruning, it completed in 39.42 seconds. That
is a greater-than-4.4x observed lower-bound improvement for the combined
repair, not an isolated benchmark or a Held–Pfeiffer comparison. The
structural count identifies reach pruning as the primary algorithmic lever;
the matched benchmark still has to measure it independently.

!!! note "Fixture noise can masquerade as topology"

    A historical capsule fixture expected a merge/split event generated from
    supports that could not participate in the relevant physical rim
    topology. Exact reach pruning removed that algebraic noise and exposed the
    real event: a tangent birth. The replacement fixture uses a concave
    line-only pocket and proves two genuine equal-cardinality `2 -> 2`
    endpoint-order fibres.

## Event identity under reversal

A boundary `vertex_id` identifies a physical vertex, not one traversal event.
The same motion can cross the tool-radius locus of that vertex twice—once
entering and once leaving. Keying events only by `vertex_id` collapses a split
and a merge.

Reversal tests therefore preserve fibre order. They compare forward fibres
with reversed reverse-motion fibres and require:

- the same physical vertex and incident supports;
- equal active cardinality on both sides;
- exact incidence-permutation rechecks; and
- opposite split/merge disposition at the same physical crossing.

## Typed failure model

`MotionCertifier.certify(...)` accepts only lateral operations:

| Native outcome | Python result | Meaning |
| --- | --- | --- |
| `certified` | immutable `MotionWitness` | Complete exact proof that the effective cap is respected |
| `cap_exceeded` | `EngagementCapExceededError` | Physical infeasibility proved; candidate evaluation may reject it normally |
| `unresolved` | `UnresolvedMotionEventError` | Proof completeness failed; never relabel as physical infeasibility |
| tuple/trace contradiction | `InvalidMotionCertificateError` | Native consumer contract is corrupt |

Before native dispatch, the certifier proves the exact cap surrogate satisfies
`effective_cap <= user_cap` and rejects kind/motion mismatch. Before returning
a witness, it independently verifies the trace digest, decision-authority
digest, authority inclusion, verdict agreement, and strategy identity.
`MotionWitness.__post_init__` repeats the kind/motion and cap-order invariants,
so direct construction or `dataclasses.replace` cannot bypass them.

### Canonical motion and replay handoff

`MotionWitness.canonical_bytes` is the single serialization authority for a
returned engagement proof. Its `motion-witness-v1` record contains every
witness field, using explicit versioned CUT/LINK and certified-verdict tags
rather than Python enum representation. `MotionWitness.digest` is SHA-256 of
that record. Changing only the operation index, effective cap, oracle strategy,
stock lineage, event-trace digest, or event-cell count changes the content
identity.

Replay wraps that witness in `ReplayLateralWitness` only after the same
operation has also produced exact containment, depletion, and coverage
evidence. The wrapper requires:

- the recorded operation and recomputed cap decision to agree;
- operation ordinal and semantic kind to agree with `MotionWitness`;
- all four proof records to own the same exact motion;
- containment, depletion, and coverage to use the same typed tool radius; and
- motion-witness user/effective cap bytes to equal the recomputed decision.

`FreshReplayTrace` orders those wrappers behind the one entry-depletion
witness. It independently recomputes stock and coverage parent lineages and
binds pristine, post-entry, and terminal state identities. Repeating a real
circle-link-circle replay twice emits byte-identical canonical traces;
swapping the two valid lateral bundles or borrowing a valid circle sweep for
the link raises `InvalidReplayTraceError`.

!!! warning "Canonical trace does not imply completion"

    The trace records authentic evidence available before the terminal gate.
    It has no terminal verdict. Only `ReplayCertificate` may bind terminal MAT
    traversal and exact empty reachable residual, and that certificate remains
    pending.

## Verification

The focused stage gates are:

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest --testmon -n auto -q \
  tests/adaptive/test_segment_event_substrate.py \
  tests/adaptive/test_segment_event_proof_contracts.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_replay.py \
  tests/adaptive/test_identity.py
pixi run -e docs docs
```

The comparative claim and matched performance boundary remain in
[Exact Segment-Site Medial Axis](segment_site_mat.md#relation-to-held-and-pfeiffer-2025).
