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
    rejection. Task 13A retains the exact neck loci and all three certified
    sides of each production plateau neck. Task 13B now carries those sides
    through one content-addressed route over every MAT edge and derives a
    unique pending oriented neck scope from causal side history. Task 13C now
    consumes increasing and decreasing native spans, retains direction through
    derived cursors, and freshly reconstructs one reverse P–S candidate without
    using cutter rotation as topology. Task 13D now composes the qualified
    entry, pristine exact stock, empty coverage, first full-circle oracle, and
    global cursor into one independently replayed launch transaction. The
    Task 13E now cross-binds physical and global continuation, materializes
    each directed forward span once, searches one invariant-ordered finite
    family, and stops at the first accepted transaction. Expected
    gouge/cap/degenerate-link rejection may continue; unresolved native
    authority aborts immediately. A launch-rooted `GenerationContinuation`
    records the resulting prefix, but the adopted real L family currently
    exposes an unresolved segment-partition case before global completion.
    A shorter radius-1 Task 13F family then exposed a one-root algebraic trim
    boundary after its first continuation link. Source v4 admits that
    coordinate family without binary64 collapse and separates physical roots
    from norm conjugates; rational-coordinate source v3 remains independently
    reconstructible. The real 16-cell ordered family now proves its first
    three trials are gouges and independently commits trial 4. At that
    circle's phase seam, incomparable exact active sets prove a mixed
    transition even though tangent incidence evidence remains inactive on
    both adjacent cells. The accepted transaction binds the independently
    reproduced native trace digest, and the first MAT route becomes terminal.
    Full continuation next reaches a constant-clearance arm whose clearance
    equals the tool radius everywhere. Its guide radius is therefore exactly
    zero: the arm requires a certified link-only cutting advance, not a fake
    zero-radius circle or a traversal-only skip.
    Traversal/coverage closure, fresh terminal replay, arbitrary-pocket
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
`event-exact-motion-oracle-v5` component. Changing the witness schema or
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
`event-exact-motion-oracle-v5` closes that second proof obligation for
same-support and noncoincident cross-support material runs over rational and
shared-radical one-root trim vertices.

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

For a same-support partial material run, the authority requires the existing
chart-specific support-cap projection. An equal-cap same-support run requires
the corresponding identically-equal overlap record.

### One-root vertex passages

Source v4 encodes every coordinate with one grammar:

```text
q = a + b sqrt(alpha)
```

where `a`, `b`, and the positive radicand `alpha` are exact rationals.
Rational coordinates canonicalize to `b = alpha = 0`. A vertex may use one
shared nonzero radicand across `x` and `y`; distinct nonzero coordinate
radicands raise `UnsupportedAlgebraicVertexProjectionError`.

On one owned quarter chart, a cutter-boundary passage through such a vertex
has the exact form

```text
F(t) = A(t) + B(t) sqrt(alpha).
```

The integer polynomial

```text
N(t) = den(alpha) A(t)^2 - num(alpha) B(t)^2
```

discovers every candidate partition root. It also contains roots of the
conjugate equation `A - B sqrt(alpha) = 0`. The algebraic kernel therefore
retains every root of `N`, but endpoint events are attached only after exact
re-evaluation of the original `F = 0` equation. Conjugate-only roots remain
valid partition boundaries without acquiring physical topology.

The active sheet cannot use `sign(N)`, because the norm loses the sign of the
original radical equation. Each one-root projection consequently binds jointly
normalized `A`, `B`, and `alpha`; cell classification evaluates
`sign(A + B sqrt(alpha))` exactly from coefficient signs and squared
magnitudes. `one_root_physical_root_ids(...)` and
`one_root_conjugate_root_ids(...)` recompute this provenance even after
coincident roots have been reduced to a common governing factor.

### Exact active-pair demand

For a cross-support run, the final partition contributes the eliminant of the
two boundary pullbacks with the oriented-rim determinant

```text
o = x₁ y₂ - y₁ x₂.
```

For caps below `pi`, it also contributes the reduced chord predicate

```text
g = (2 - C) d₁ d₂ - 2 (x₁ x₂ + y₁ y₂),
```

where `(xᵢ / dᵢ, yᵢ / dᵢ)` is the rational rim parametrization and `C` is the
exact cap chord ratio.

A circle-support pullback contains one exact rim-chart denominator factor
`1 + u²`. That factor has no real roots, but its complex roots are projective
chart poles. Leaving it in both pullbacks makes the pair resultant
identically zero through a nonphysical shared component. The producer now
proves exact divisibility and saturates each circle pullback by that factor
before elimination. The pair predicate likewise omits the known-positive real
factor `d₁ d₂ = (1 + u²)(1 + v²)`. Both transformations preserve every real
sign and root used by the engagement proof.

The algebraic v4 producer does not construct every unordered feature/chart
product. It:

1. builds the physical topology partition without cross-support pair factors;
2. evaluates one exact rational station in each topology cell;
3. records every cross-support branch pair active there, together with its
   owned center and rim charts;
4. sorts and deduplicates that canonical request set; and
5. rebuilds with only the requested orientation and sub-`pi` cap factors; and
6. unions any requests exposed by the refined witnesses and repeats until the
   canonical request set reaches a fixed point.

Pair roots cannot create a physical boundary branch, but a first-pass witness
can lie exactly on a pair-order coincidence that was absent from the topology
partition. Refinement then moves the child witnesses off that coincidence and
can expose a request hidden at the parent witness. Monotone closure retains
every earlier request and terminates because the feature/chart/rim request
universe is finite. The closed request set is part of the v4 source payload
and therefore reconstructs byte-identically. Completeness remains
independently fail-closed:
`FullCircleCellAuthority2` requires the matching factor for every active
material run and returns unresolved if one is absent.

These pair roots refine the motion partition but are not physical boundary
events. Their projection inputs deliberately carry an empty event list;
tangent, endpoint-order, overlap, and seam fibres remain the only topology
authority. A cell decision binds every required pair projection identity.
A zero or missing pair eliminant leaves only the cell that needs it unresolved;
an irrelevant zero does not poison the entire motion. Cross-support equal-cap
runs and unsupported overlap families remain unresolved.

`FullCircleCellAuthority2` commits to:

- the reconstructed generic partition;
- the canonical stock-feature inventory;
- exact center, phase, tool-radius, and cap inputs;
- every rational station source and complete stationary cell stratum;
- each reference-side and material-run disposition; and
- the exact cap-projection identities used by each decided partial run.

Its SHA-256 digest is therefore deliberately distinct from the generic
partition digest and is included inside `EventTrace2`. The strategy is
versioned `full-circle-cell-strata-exact-v3`. Rational-coordinate audits retain
the independently replayed `full-circle-boundary-pullbacks-v3` source.
Algebraic audits use `full-circle-boundary-pullbacks-v4`, whose seventh source
field binds the closed pair-request set.

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

The first same-support implementation accidentally rebuilt the partition
three times, reconstructed every stationary cell twice, and re-extracted the
stock boundary once per cell. Hoisting verifiedness and boundary records
reduced its warm Apple M1 Max Release median from 252 ms to 128 ms for a
10-cell concentric certificate, and from 572 ms to 301 ms for a 43-cell
cross-support audit.

The earlier complete-product theorem had a different cost. On the same machine and Release
worktree on 2026-07-29, five warm calls gave a 1,737.242 ms median for the
13-cell/576-projection concentric certificate and 2,018.540 ms for the
54-cell/452-projection cross-support audit.

The post-capsule positive control contains `309` boundary records, of which
`290` are admitted line/circle pair features. The former product would build

```text
290 * 289 / 2 * 4 center charts * 2 * 2 rim charts
    = 670,480
```

orientation resultant chains before any sub-`pi` cap factors. The first
integration run remained unfinished after `100 s`. Exact topology demand
reduced that workload to `233` canonical active-pair requests. After
circle-pullback saturation, the isolated pytest positive control proved the
expected cap exceedance in `30.46 s`. This is a structural reduction of more
than three orders of magnitude in requested pair combinations, but it is
still a slow development fixture and not Held–Pfeiffer runtime parity.

!!! warning "One-root scope is exact but deliberately narrow"

    Source v4's coordinate grammar can represent rationals and one shared
    radical extension per vertex, but the production oracle retains rational
    inputs on source v3 until a separate demand-path migration is validated.
    Neither source tag reconstructs the other. V4 does not silently combine
    unrelated radicals or approximate an unsupported algebraic field. A
    vertex with distinct nonzero `x`/`y` radicands remains an explicit
    unsupported proof boundary.

### Exact clear sweeps in nonempty stock

The first Task 11A replay circle exposed a proof case that global stock
emptiness cannot represent. Its complete cutter sweep lies inside the
qualified entry void, while material correctly remains elsewhere in the
pocket. The boundary pullback partition reconstructed successfully, but the
old uniform classifier recognized `clear` only when the stock had no boundary
records. It therefore returned `unresolved` for a physically and exactly
zero-engagement motion.

The retained uniform path in `event-exact-motion-oracle-v5` derives a
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

!!! warning "Derived scope is not yet integrated scope"

    The candidate lattice currently accepts scope as an identity input. The
    same exact motion/traversal can be rebuilt with no-neck/full-cap identity
    or with oriented-neck/restricted-cap identity. Replay validates the latter
    internally but cannot infer which separator was crossed from motion
    geometry or one edge cursor. Task 13B now carries complete certified-side
    history and derives that causal scope from the global route. The later
    generator/replay stages must consume `MatTraversalState.neck_scope`;
    accepting a separately supplied scope would leave two authorities.

!!! warning "A legal passage does not certify its circle"

    The next radius-`1/8`, phase-index-`3` circle advances the same neck passage
    legally, but exceeds its reconstructed `80 degrees` cap against post-link
    stock. That negative is an integration gate. An equal-cap fixture separately
    proves the two state transitions. Task 12 now rejects that complete joint
    trial; the later Task 13 candidate loop must continue the feasible search
    instead of weakening the oracle or conflating passage legality with motion
    acceptance.

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

Task 13B is now the sole authority for assigning causal neck scope and choosing
the active global MAT branch. Task 13C lets that authority consume either
native ordinal direction without changing positive MATHSM progress or cutter
rotation. Later Task 13 stages decide which candidate failures constitute
normal search exhaustion. Task 12 validates a supplied finite candidate and
state transition; it does not infer global traversal facts.

### Task 13D entry-circle bootstrap

The first physical circle is not a degenerate instance of the Task 12
link/circle transaction. Its exact phase already equals the qualified bore
center, so a synthetic link would have zero progress and cannot enter the
closed motion grammar. `InitialCandidateEvaluator` instead rebuilds pristine
stock from `InputIdentity`, applies `PreclearedEntry` once, creates an empty
`CoverageLedger`, and evaluates only the full circle.

The launch preserves the same deciding proof order as every later circle:

```text
entry-sweep containment
    -> design-domain containment
    -> MotionCertifier.certify(post-entry stock)
    -> exact circle depletion
    -> exact coverage sweep
    -> functional traversal child
```

The immutable `InitialCandidateTransaction` binds the complete input digest,
global traversal parent and child, pristine-stock boundary, entry-depletion
witness, empty initial coverage certificate, circle-in-entry certificate,
per-operation replay witness, and resulting physical-state digest. Evaluation
operates on private stock/coverage owners and returns evidence only. Commit
reconstructs the complete launch and requires byte-identical transaction
evidence before returning the separate `GenerationState` and
`MatTraversalState` children.

The adopted L fixture exercises a real reverse P–S route. Its explicit launch
lattice admits a guide radius of `1/32 mm`; with the `1/2 mm` tool, the complete
sweep fits exactly in a `9/16 mm` qualified entry disk centred at the phase:

```text
2(1/32 mm) + 1/2 mm = 9/16 mm
```

That equality is an exact containment case, not a tolerance. Shrinking the
entry to `11/20 mm` leaves the same circle inside the design pocket but rejects
it as a launch because the cutter sweep is no longer proved to lie in the
qualified void.

!!! note "A seeded entry establishes a neck side; it does not cross one"

    `MatTraversalState.seed()` has no previously occupied side and therefore
    no pending `CausalNeckTransit`. The entry circle must use `NoNeckScope`,
    advance exactly one MAT cursor, and leave every oriented `NeckPassage`
    `UNVISITED`. A causal passage first becomes possible after a terminal edge
    activates a later route branch with distinct authenticated source and
    target sides. Inventing a causal launch fixture would encode history that
    never occurred.

The focused mutation gate rejects a nearby but unequal phase, undersized
entry, zero-length synthetic link, foreign MAT sampling, foreign tool,
foreign bore evidence, foreign design-domain proof, foreign cut depth, stale
global parent, certification after self-depletion, cross-wired circle witness,
and a child with two advanced cursors. Two independent trials and the commit
replay produce identical canonical bytes.

!!! danger "Exact-number backend definitions are an ABI contract"

    `_continuous_tea_2` compiles with its target-local CORE/Boost backend
    definitions. Linking a sweep object compiled under a different exact
    backend can give the same `ReachKernel` spelling different effective
    number types across a C++ boundary. The oracle therefore compiles the
    authoritative `exact_sweep_2.cpp` source inside its own target. Reuse the
    deciding source, but never share compiled exact-kernel objects until every
    consumer inherits one verified interface definition set.

### Task 13E global continuation and finite search

Task 13E keeps Task 12 as the only deciding link-and-circle engine.
`CandidateEvaluator.evaluate_from_cursor(...)` accepts the physical cursor
authenticated by the global route; the existing same-edge API delegates to
that same implementation. Containment, motion certification, depletion,
coverage, passage advance, and witness construction therefore retain one
proof order across same-edge and branch-switch consumers.

`TraversalCommit` closes the square after independent replay:

```text
physical parent  --Task 12 transaction-->  physical child
       |                                        |
       | same candidate                         | same candidate
       v                                        v
global MAT parent  --------------------->  global MAT child
```

It binds both parent/child digests, requires exactly one active global cursor
to change, preserves every other cursor byte-for-byte, and verifies causal
neck scope against the pending transit. A stale physical parent, stale global
parent, valid child from another candidate, or manufactured neck scope cannot
acquire one commit identity.

For one active route cursor, `materialize_active_candidate_family(...)`:

1. derives the directed native limits from `TraversalPolicy.forward_window`;
2. derives `FullCapDecision` or the pending transit's exact
   `NeckCapDecision` from current physical passage state;
3. enumerates each exact span once with the shared candidate grammar;
4. combines all cells in furthest-progress, largest-radius, canonical-identity
   order; and
5. submits that immutable tuple to the single Task 12 trial engine.

The search may continue only after `GougeContainmentError`,
`EngagementCapExceededError`, or a candidate-specific
`DegenerateSegmentMotionError`. It stops at the first accepted transaction and
does not evaluate lower-ranked cells. Exhaustion reports fixed canonical
counts plus the active cursor identity as
`EngagementCapInfeasibleError`, `NeckTooTightError`, or
`NoFeasibleCandidateError`; rejected mutable trial owners are not retained.

`GenerationContinuation` roots the post-launch chain in the independently
replayed `InitialCandidateTransaction`, normalizes deterministic zero-cut
route activation between commits, and binds the current physical state,
global state, and ordered `TraversalCommit` tuple. It deliberately does not
assert exact residual emptiness. The terminal-seal stage must still call the
coverage gate and rebuild the complete run through fresh replay before a
`GenerationResult` exists.

!!! warning "Three cursor endpoints carry different proof"

    A candidate may stop inside a native span, exactly on an intermediate
    native window limit, or at the finite family's terminal cell. The first
    retains a `DerivedCandidateCursor`; the second must retain the owned
    `MatSample` so the next exact span starts from the same algebraic parameter;
    the third retains either the native route endpoint or an
    `ExhaustedCandidateCursor` at the last positive-radius cell. Exact finite
    exhaustion before a tool-radius MAT leaf is terminal traversal evidence,
    not a claim that the cutter reached the native endpoint.

!!! danger "Unresolved substrate is not candidate infeasibility"

    The first real ordered continuation family on the adopted scaled L fixture
    reaches a native `IncompleteSegmentPartitionError` on a high-ranked link.
    `MotionCertifier` translates the closed native proof gap to
    `UnresolvedMotionEventError`, retains the native exception as its cause,
    and the search aborts before later accepted cells. Searching around that
    event would silently convert missing proof into feasibility.

Focused gates cover a real exact branch-switch commit, stale/cross-wired
parents, invariant ordering, first-feasible stopping, all three exhaustion
classes, causal neck-cap derivation, one materialization per forward span,
launch-rooted continuation, intermediate native cursors, and immediate
propagation of the real unresolved substrate. Complete terminal coverage is
not part of this stage.

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
| incomplete/unsupported native event substrate | `UnresolvedMotionEventError` with native cause | Backend authority is incomplete; generation aborts immediately |
| malformed native input or trace-verification failure | `InvalidMotionCertificateError` with native cause | Typed/native proof contract is corrupt |
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
  tests/adaptive/test_identity.py \
  tests/adaptive/test_transaction.py \
  tests/adaptive/test_traversal.py \
  tests/adaptive/test_acceptance.py \
  tests/adaptive/test_generator.py
pixi run -e docs docs
```

At the Task 13E stage gate on 2026-07-29, affected selection passed `113`
tests and the complete adaptive suite passed `518` tests in `223.88 s`. Ruff,
strict mypy over `27` source files, strict MkDocs, and diff checks passed.
These are correctness gates, not planner-performance measurements.

At the cross-support v4 prerequisite gate on 2026-07-29, the complete circle
oracle file passed `17` tests, the affected consumer/replay slice passed `97`,
and the complete adaptive suite passed `521` tests in `230.33 s`. Ruff, strict
mypy over `27` source files, and strict MkDocs also passed.

At the source-v4 implementation checkpoint on 2026-07-29, the exact
physical/conjugate fixture passed with `6` physical and `14` conjugate-only
one-root roots, and all four one-root mutations failed reconstruction. The
`290`-feature post-capsule positive control separately passed with the named
cap-exceedance exception in `30.46 s`. The complete circle/motion files passed
`30` tests in `40.25 s`, and the complete adaptive suite passed `522` tests in
`242.54 s`. Ruff, strict mypy over `27` source files, strict MkDocs, affected
`--testmon`, and diff checks passed.

At the Task 13F route-0 integration checkpoint on 2026-07-29, the real
16-cell family dispatched exactly four trials: three named gouges followed by
the accepted winner. Independent commit replay reproduced the physical child,
and a separately reconstructed eight-cell native circle trace matched the
transaction's event-trace digest. Its phase-seam fibre retained exact inactive
endpoint incidences while certifying mixed directions. Full continuation then
failed loud at the next route with an empty circle family because the complete
MAT arm has exact zero guide radius. Task 13F is not terminal; the required
link-only advance remains an implementation stage.

The focused seam/Task 13F slice passed `4` tests in `28.22 s`; affected
`--testmon` selection passed `5` in `29.37 s`; and the complete adaptive suite
passed `523` in `245.18 s`. Ruff, strict mypy over `27` source files, strict
MkDocs, and `git diff --check` also passed.

The comparative claim and matched performance boundary remain in
[Exact Segment-Site Medial Axis](segment_site_mat.md#relation-to-held-and-pfeiffer-2025).
