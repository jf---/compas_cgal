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
    outcomes. Candidate selection, certify-before-deplete orchestration,
    fresh artifact replay, arbitrary-pocket evidence, and matched
    Held–Pfeiffer performance remain later tasks.

## Authority hierarchy

The records have deliberately different jobs:

| Layer | Circle authority | Segment authority | Consumer role |
| --- | --- | --- | --- |
| Reconstructed partition | `EventPartitionCertificate2` | `SegmentEventPartition2`, including its nested `EventPartitionCertificate2` | Proves that every declared projection, algebraic root, cell, fibre, overlap, and seam can be reconstructed |
| Deciding data | Generic cells/fibres and uniform/full-circle dispositions | Ordered segment strata, active branches, and pair orientation/cap dispositions | Authorizes the exact verdict |
| Trace | `EventTrace2` v2 | `EventTrace2` v2 | Orders canonical events and binds the deciding record by SHA-256 |
| Typed witness | `MotionWitness` | `MotionWitness` | Binds operation ordinal/kind, exact motion, caps, stock lineage, strategy, and trace digest |

A generic event certificate is sufficient authority for a circle trace. It is
not sufficient for a segment verdict: segment cell classification consumes
the larger `SegmentEventPartition2`, especially the ordered active-branch
inventory and each material run's pair-cap disposition.

`EventTrace2` v2 therefore carries two views:

- `partition` is the generic certificate used by common trace consumers;
- `decision_authority_bytes` and `decision_authority_digest` identify the
  complete reconstructed record that authorized the verdict.

For circles, the authority bytes equal `partition.canonical_bytes`. For
segments, they equal the full segment partition bytes and deliberately differ
from the nested generic certificate. The trace canonical record includes both
the generic partition digest and the decision-authority digest. A
`MotionWitness.event_trace_digest` is consequently transitive to the full
deciding record.

!!! danger "A digest beside a proof is not proof ownership"

    Recording the full segment digest as metadata without including it in the
    trace identity would allow the same trace digest to be relabelled onto a
    different deciding partition. The authority digest must be inside the
    canonical trace record.

The planning `InputIdentity` separately binds
`motion-certificate-schema-v1` and the native
`event-exact-motion-oracle-v1` component. Changing the witness schema or
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
  tests/adaptive/test_identity.py
pixi run -e docs docs
```

The comparative claim and matched performance boundary remain in
[Exact Segment-Site Medial Axis](segment_site_mat.md#relation-to-held-and-pfeiffer-2025).
