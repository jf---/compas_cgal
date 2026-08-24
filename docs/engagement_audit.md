# Truthful Engagement Audit

## Maturity

Stage 1 Tasks 1–3 define the authenticated input and operation-record boundary
and the exact partial-arc surrogate used by stock depletion. The audit replay,
native motion certification, report aggregation, regulated-generator evidence,
and release gate remain incomplete. This stage therefore proves admissible
input classification and exact-on-surrogate depletion, not toolpath compliance.

## Authority boundary

`EngagementAuditInput.build(...)` is the only mutable COMPAS ingress. It
canonically encodes every audit-relevant `ToolpathOperation` field and invokes
the `_stock_2` Epeck classifier. The completed input retains only:

- immutable records carrying opaque Epeck-classified motion values;
- the ordered source-operation digest;
- canonical design, cut-plane, tool, and cap values; and
- complete native, Python, component, and lockfile build identity.

Caller-owned `ToolpathOperation` objects are never retained. Mutating their
geometry, role, orientation, or path metadata after construction cannot change
the audit request.

The ordered `operation_index` is assigned at this ingress boundary. It is the
authoritative stream ordinal. `ToolpathOperation.path_index` remains generator
metadata in the authenticated source bytes and is never used as replay order.

## Exact classification

All branch-producing geometric decisions execute in
`audit_classification_2.cpp` after exact binary64 injection into Epeck:

| Primitive | Native decisions |
| --- | --- |
| line | exact XY coincidence, endpoint Z equality, cut/clearance plane, vertical direction, ramp rejection |
| circle | exact orthonormal world-XY frame, positive orientation and radius, cut/clearance plane |
| arc | circle decisions plus one exact-injected authored binary64 sweep observation, full-turn bound, orientation consistency |

Operation labels can contradict a native geometric classification and cause a
named failure. They never prove that a motion is non-engaging. Vertical retract
and clearance-plane transport are the non-mutating native alternatives. A
native-proved plunge is retained separately because it must deplete a disk
before later motion certification. Its endpoint exists only inside the opaque
native value; Python retains no parallel geometry. Every cut-plane lateral
line, circle, or arc becomes an opaque Epeck motion value for a later native
certifier.

The six nanobind motion classes have no Python constructor and expose no
reconstructive coordinate getters. An arc uses the versioned
`audit-arc-quarter-chart-binary64-v2` seam. Ingress observes the authored
binary64 subtraction `end_angle - start_angle` once, rejects a nonfinite
result, and exact-injects that observation. Epeck then owns sweep sign,
orientation, full-turn extent, normalization, and quarter-chart selection.
After those exact decisions, `tan(local_angle / 2)` is evaluated directly from
the original authored angle and the exact-selected, binary64-exact integer turn
and chart indices; no exact local angle is converted back to binary64. That one
tan result is exact-injected as the rational parameter of a frozen Pythagorean
quarter chart. Exact seam angles map structurally to parameter zero and never
call `tan`.

The native arc retains its canonical start/end chart coordinates, ordered
trimmed chart intervals, exact guide radius, cut plane, authored-angle identity,
strategy versions, and a 32-byte SHA-256 motion digest. There is no Python
`point_at`, `atan2`, phase subtraction, chart-coordinate getter, or angle
reconstruction path. The native strategy identifier remains bound into audit
input identity. Task 4 will additionally bind the opaque motion digest into the
native replay request and ordered authenticated-operation digest.

Python snapshots curve radius as a finite, millimetre-bearing observation. It
does not decide positivity. `CGAL::sign` after exact injection is the sole
positive-radius authority and its named native rejection crosses the adapter.

Each authenticated lateral, plunge, and non-engaging carrier has its own
versioned canonical encoding and SHA-256 digest. The encoding binds stream
ordinal, source-operation digest, and the closed native classification tag;
source identity binds the opaque geometry itself.

Task 2 accepts only frames whose binary64 axes satisfy exact orthonormal
world-XY predicates. Scaled, skewed, tilted, and inexactly normalized rotated
frames fail closed. Expanding that domain requires a proved native
normalization contract; no tolerance or COMPAS fallback is permitted.

## Failure model

The native boundary exposes distinct exceptions for nonfinite input, invalid
cut/clearance planes, unsupported geometry, off-plane motion, contradictory
roles, and contradictory orientation. Python maps those exception types—not
message text—to the engagement-audit error model.
Empty streams, malformed source records, unsupported primitives, mixed-Z
ramps, and undeclared depths fail before stock construction or mutation.

## Exact partial-arc depletion

`exact_circle_chart_atlas_2.h` is the single frozen atlas consumed by both the
Epeck guide evaluator and the continuous-TEA polynomial charts. Its four CCW
records fix chart identifiers, numerator/denominator coefficient triples,
domain `[0,1]`, and parameter-zero seam ownership. Full-circle depletion and
partial-arc depletion therefore evaluate the same exact rational maps as the
later certifier substrate; duplicated coefficient tables are prohibited.

`Stock2.subtract_exact_arc(...)` accepts only the opaque `AuditArcMotion2`.
For each declared trimmed interval, native code dyadically refines its exact
parameter domain until every consecutive squared chord is at most the exact
policy bound. Cardinal centers use the owning chart's parameter-zero
representation exactly once. Every partial-arc trace is non-cyclic, including
a full-turn `AuditArcMotion2`. A full turn retains its explicit duplicate
terminal anchor and parameter, so the closing chord is the final consecutive
segment; it is never represented by an implicit cyclic last-to-first edge.

The structural validator regenerates the complete ordered parameter sequence
from the motion and chord policy and compares every exact chart/parameter pair.
Reordered, missing, duplicated, reversed, or complementary sequences fail even
when all their points remain incident on the same guide circle. It also checks
exact incidence, anchors, direction, seam ownership, center-count limit, and
the exact policy relation `0 < chord_bound < tool_radius`.

After construction, stock depletion subtracts full-radius exact disks from a
cloned polygon set, validates the immutable trace, then swaps the clone into the
stock. Every construction, policy, limit, or trace failure occurs before that
swap. The trace canonical bytes bind the native motion digest, tool radius,
chord bound, center-count limit, non-cyclic state, strategy version, and ordered
exact parameters; its witness digest is SHA-256 over those bytes. Exact
rationals use the repository's authoritative CCAN encoder after
`CGAL::Fraction_traits` decomposition—never decimal text or `to_double`.

Every disk center is exactly incident on the declared rational-chart guide, so
the disk union is a subset of that surrogate sweep. The strict chord/tool
relation proves adjacent disks overlap. It does not prove a quantitative bound
on retained material between disks, and it does not identify the rational
surrogate with an ideal transcendental COMPAS arc.

## Evidence

The Task 2 TDD sequence first produced eight missing-native-API failures. The
implemented native and Python boundary is covered by exact-nextafter plane
tests, vertical and clearance classification, ramp and role contradiction,
tilted/skewed/scaled-frame rejection, arc-direction rejection, unforgeable
native values, coherent one-shot capture, immutable source mutation,
content-address mutation, exact cap-surrogate identity, and strict consumer
typing. The dedicated
`types-audit` task applies `mypy --strict` to the package and its consumer
contract.

Task 3 RED first failed because the arc carrier had no canonical intervals,
motion digest, shared atlas, exact arc depletion trace, witness identity, or
atomic public subtraction. The adversarial RED round then required canonical
seam ownership, foreign-sequence rejection, exact operand-level sweep and
quadrant decisions, authoritative CCAN bytes, and shared-atlas consumers. The
review repair additionally covers ordinary non-seam full turns in both
directions, one authored sweep observation, explicit terminal closure,
foreign-motion trace identity, finite public inputs, and named factory-bypass
failures. The
focused native gate covers awkward and exact-seam angles, negative and
multi-turn starts, CW/CCW minor/major/full traversals, nextafter seam values,
rational rotation/translation/scale, digest mutations, every malformed trace
shape, nonpositive and incoherent policies, and bounded allocation. The focused
Python gate additionally covers opaque motion identity, source guards, the
linked public `Stock2` transaction, and exact stock equality after every
rejection.

The earlier Python-owned geometric classifier was rejected during review and
is not an accepted implementation or evidence source. The production path has
one classifier: `_stock_2` with Epeck predicates.

## Remaining work

Stage 1 must still implement measure-before-deplete replay, the native
three-way segment/circle/arc certification adapters, truthful report
aggregation, consumer migration, and the bounded evidence matrix. Before replay
can become authoritative, Task 4 must make the native certifier consume this
same opaque rational-chart motion and complete the atomic certify-and-deplete
transaction for segment, circle, and arc. The legacy trigonometric
`subtract_arc_sweep` remains a separate non-authoritative API and is unreachable
from the audit arc sources. Until those steps pass, no audit result can claim
complete engagement compliance.
