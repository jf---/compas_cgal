# Truthful Engagement Audit

## Maturity

Stage 1 Tasks 1 and 2 define the authenticated input and operation-record
boundary. The audit replay, native motion certification, report aggregation,
regulated-generator evidence, and release gate remain incomplete. This stage
therefore proves admissible input classification, not toolpath compliance.

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
| arc | circle decisions plus exact injected sweep sign, full-turn bound, orientation consistency |

Operation labels can contradict a native geometric classification and cause a
named failure. They never prove that a motion is non-engaging. Vertical retract
and clearance-plane transport are the non-mutating native alternatives. A
native-proved plunge is retained separately because it must deplete a disk
before later motion certification. Its endpoint exists only inside the opaque
native value; Python retains no parallel geometry. Every cut-plane lateral
line, circle, or arc becomes an opaque Epeck motion value for a later native
certifier.

The six nanobind motion classes have no Python constructor and expose no
reconstructive coordinate getters. For arcs, `audit-arc-phase-binary64-v1`
computes the start phase with `sin`/`cos` once inside the native boundary,
exact-injects that binary64 surrogate, and retains the Epeck vector. There is no
Python `point_at`, `atan2`, phase subtraction, or angle reconstruction path.
The native module exposes that strategy identifier and every audit-input digest
binds it, so changing the transcendental seam changes request identity.
Circle and arc values also retain their exact-injected guide radius internally,
so later native depletion and certification never reconstruct it from a phase
vector through an unavailable square root.

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

The earlier Python-owned geometric classifier was rejected during review and
is not an accepted implementation or evidence source. The production path has
one classifier: `_stock_2` with Epeck predicates.

## Remaining work

Stage 1 must still implement measure-before-deplete replay, the native
three-way segment/circle/arc certification adapters, truthful report
aggregation, consumer migration, and the bounded evidence matrix. Before replay
can become authoritative, Task 3 must replace the current independently
injected arc phase/radius pair with one exact rational-chart surrogate shared
by certification and exact-on-surrogate depletion. The legacy trigonometric
`subtract_arc_sweep` is not proof-bearing depletion. Until those steps pass, no
audit result can claim complete engagement compliance.
