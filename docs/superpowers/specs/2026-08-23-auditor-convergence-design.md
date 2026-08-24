# Exact Auditor Convergence Design

**Date:** 2026-08-23

**Status:** approved 2026-08-23

**Scope:** reconcile the exact-certified adaptive branch family and reach one
truthful, replayable, CI-enforced engagement-audit milestone

## Decision

Converge the existing exact MAT/zero-guide spine; do not rewrite it and do not
merge the divergent branch histories wholesale.

The programme has four separately reviewable deliverables:

1. a content-equivalence reconciliation of the two branch frontiers;
2. an audit contract in which unmeasured can never mean certified;
3. a successful fresh replay plus verified motion-theorem boundary; and
4. CI-enforced, content-addressed evidence from the authoritative path.

Generator quality and commercial positioning remain downstream decisions. The
quality instrument is the nearer-term product candidate, but it is not called
saleable until all four deliverables are green and a fresh third-party corpus
run has been committed.

This design decomposes into four implementation plans. A plan may begin only
after its predecessor's acceptance commit has passed focused tests, strict
typing, documentation, and independent review.

## Audited starting state

The external review assessed commit `e0ae116`. This design is based on a fresh
read of both live frontiers:

| Role | Branch | Audited tip | State |
| --- | --- | --- | --- |
| integration base | `codex/exact-certified-adaptive-phase1-t9-zero-guide` | `073a0f7` | clean, five commits ahead of origin |
| certifier source | `jf/toolpath-redesign` | `73d5372` | one commit ahead of origin; unrelated dirty checkout content excluded |
| common ancestor | both | `1860167` | 254 commits on integration side, 36 on certifier side |

The integration base retains the exact MAT, zero-guide candidate,
transaction, traversal, replay substrate, benchmark corpus, and later
generator experiments. The certifier source retains the release-build
shared-root guard, corrected full-turn reporting theorem, false-certificate
corpus, and incomplete adaptive arc certifier.

Current defects that define this programme are:

- unmeasured operations are represented by `cap_certified=True`;
- the integration base uses fixed-density Python arc certification;
- `replay_certificate` has no value-returning path;
- the swept-prefix certifier certifies from two booleans without a geometric
  acceptance test;
- the strict adaptive mypy gate reports six errors;
- package metadata promises Python 3.9 while adaptive modules import
  `typing.Self`;
- CI does not execute the adaptive proof suite and contains dead or invalid
  tasks;
- production component identities do not bind the actual native source tree;
- generator override decisions disappear at public result boundaries; and
- benchmark and comparison claims are not derived from committed measurement
  artifacts.

## Goal and exit claim

The programme exits only when the following claim is true:

> For one authenticated planar toolpath, the audit either proves every
> engagement-relevant motion compliant, identifies a proved cap exceedance, or
> fails loud with an unresolved/unsupported reason. It never converts missing
> measurement into compliance. The same operation stream reconstructs a fresh,
> content-addressed terminal replay certificate, and CI reruns the exact
> contracts that justify the claim.

The exit claim does not assert:

- that the generator produces useful paths on every corpus pocket;
- that exact certification meets Held and Pfeiffer's runtime;
- that between-position engagement is bounded by sampled reporting;
- that unresolved rates are commercially acceptable before measurement; or
- that the certifier, generator, or quality instrument is ready for sale.

## Rejected approaches

### Wholesale branch rebase

Replaying all 36 `jf/toolpath-redesign` commits over the 254-commit integration
line would combine independently implemented benchmark foundations, build
changes, legacy generator repairs, and native certifier work. Conflict
resolution would silently choose between semantically different copies. Commit
ancestry is not evidence of content equivalence, so this route is rejected.

### Fresh rewrite

The exact MAT, zero-guide distinction, operation grammar, transaction
chronology, and mutation corpus are already strong. Rewriting them would
discard verified proof boundaries while leaving the demonstrated audit and CI
defects untouched.

### CI-first cleanup

Connecting a red or semantically false contract to CI only automates the wrong
answer. CI wiring follows audit and replay truth; it does not substitute for
them.

### Product-track split before convergence

Forking a generator product and an auditor product now would duplicate the
same unresolved motion, identity, and evidence boundaries. The shared proof
substrate converges first. Product branches may split after the evidence gate.

## Programme invariants

### Exactness

- Exact-kernel decisions use CGAL predicates or exact number-type comparisons.
- Reporting doubles never feed back into certification, topology, identity, or
  operation classification.
- Transcendental caps cross one declared seam as exact rational surrogates.
- Sampled or raster checks are falsifiers and reporting instruments, never
  acceptance authorities.
- An incomplete theorem or unsupported geometry returns unresolved or raises a
  named exception; there is no fallback certificate.

### Types and construction

- Python is 3.12 throughout metadata, Pixi, Ruff, mypy, CI, and wheel policy.
- Frame, length, angle, and squared-length information reuse the adaptive typed
  units; new audit APIs do not expose bare physical floats.
- Domain records use validated `build(...)` factories and named exceptions.
- Closed unions carry genuinely different chronologies or proof states.
- No enum, result wrapper, or class hierarchy is introduced for ordinary
  control flow.
- No `__all__` is introduced; the existing prohibited declaration is removed
  only in its own reviewed compatibility change.

### Architecture

- One production call graph owns each deciding question.
- The legacy diagnostic remains available during validation but no acceptance
  consumer may call it after the authoritative audit lands.
- Existing code is not deleted during add-and-validate work. Removal requires
  explicit user permission after independent equivalence evidence.
- Files split by responsibility: input identity, motion classification,
  native certification adaptation, replay, reporting, and artifact emission do
  not share one module.
- Generated measurements and figures derive from checked-in source data and a
  repository task; hand-transcribed numbers are not evidence.

### Git and publication

- All work occurs in isolated, agent-owned worktrees and non-main branches.
- The two source tips remain unchanged throughout reconciliation.
- Integration history is linear and consists of focused, reviewable commits.
- No source branch or worktree is deleted during this programme.
- Every publication verifies local ancestry, remote SHA, and clean status.
- In-flight CI for the target branch is cancelled before a later push.

### Implementation and verification

- Pixi exclusively owns dependency resolution, builds, tests, lint, typing,
  documentation, artifact generation, and release tooling.
- Every pytest invocation uses `-n auto`; changed-code gates additionally use
  `--testmon` after the explicit RED has been observed.
- Reference tests are immutable evidence. Repairs change production code or add
  new tests; they do not weaken existing assertions.
- Tests are never skipped or marked expected-failure. Deliberately red product
  criteria live in a separately invoked task.
- Python edits pass Ruff formatting and lint before commit.
- Missing dependencies fail at import; no conditional imports, capability
  flags, silent fallbacks, or swallowed exceptions are permitted.
- Every commit uses author and committer
  `Jelle Feringa <jelleferinga@gmail.com>` and contains one coherent deliverable.

## Authority graph

```text
Authenticated audit input
    |-- canonical design rings + holes
    |-- typed cut plane
    |-- typed tool radius + cap
    |-- canonical operation-stream digest
    `-- native BuildIdentity
                |
                v
Geometry-derived operation classification
    |-- proved non-engaging motion ------------> non-engaging record
    |-- supported lateral cut -----------------> native motion certifier
    `-- unsupported / contradictory geometry --> named failure
                                                    |
                                                    v
                                  certified | cap_exceeded | unresolved
                                                    |
                       +----------------------------+------------------+
                       |                                               |
                       v                                               v
            content-addressed audit report                  fresh replay
                                                                  |
                                                                  v
                                                    terminal ReplayCertificate
                       |                                               |
                       +----------------------------+------------------+
                                                    v
                                      schema-validated evidence artifact
```

No label supplied by the toolpath author can bypass geometry classification.
No reporting aggregate can turn `unresolved` into `certified`.

CGAL owns the deciding geometry and proof semantics. COMPAS geometry exists
only at the authenticated ingress adapter. One frozen, frame-and-unit-typed
capture feeds both canonical identity and the native call. Successful
classification returns opaque Epeck-backed nanobind values that later native
certifiers consume directly; Python exposes no reconstructive scalar getters.
Any ambiguous or contradictory adapter mapping fails before native execution.
The native arc-phase strategy version is part of audit-input identity.

## Stage 0: branch reconciliation and baseline ledger

### Purpose

Determine what must move from `jf/toolpath-redesign` without importing
duplicate or superseded implementations.

### Reconciliation ledger

Create a committed state document binding:

- source tip, integration tip, merge base, and source-tree status;
- every one of the 36 source commits;
- files and public/native symbols affected;
- disposition `equivalent`, `superseded`, `required`, `dependent`, or
  `unrelated`;
- the integration-side evidence for equivalent/superseded dispositions;
- the exact tests that authenticate each required change; and
- the intended destination plan and commit.

`equivalent` requires symbol-level comparison plus a passing contract test. A
larger file or similar comment is insufficient. `superseded` requires a
strictly stronger current contract with tests. `required` changes are reapplied
as focused commits from the integration base. Source commits are never merged
merely to preserve their SHAs.

The minimum source closure to adjudicate includes:

- exact seam validation and bounded disk-chain construction;
- the shared compiled TEA guard source;
- false-certificate witness liveness tests;
- swept-annulus segment certification;
- rotation-invariant interior bounds;
- sub-ulp rim reporting repair;
- release-build shared-root enforcement;
- corrected full-turn slack/precondition proof; and
- adaptive arc certification with its negative and positive controls.

### Baseline ledger

Record commands, exit codes, immutable tree SHA, and results for:

- focused engagement tests;
- focused adaptive replay/oracle tests;
- strict adaptive mypy;
- Ruff;
- MkDocs strict build;
- benchmark instrument tests; and
- the full suite, with aspirational quality gates reported separately.

Known-red product-quality tests remain an explicit `quality-gate` command.
They are not skipped, marked expected-failure, or included in the green
regression baseline.

### Stage 0 acceptance

- Both original tips remain byte-identical and clean relative to their starting
  status.
- Every source commit has exactly one ledger disposition.
- Every required patch has a named target stage and proving test.
- No code patch lands in Stage 0.
- The baseline ledger distinguishes environmental, regression, and deliberate
  product-gate failures.

## Stage 1: truthful engagement audit

### New authoritative package

Add a focused `compas_cgal.engagement_audit` package:

| File | Responsibility |
| --- | --- |
| `identity.py` | validate content-addressed Python/native component inputs as `BuildIdentity` |
| `operation_identity.py` | snapshot and canonically encode mutable COMPAS operation ingress |
| `input.py` | retain immutable classified operations and content-address the authoritative request |
| `classification.py` | translate the closed native classifier result into typed Python motion records |
| `records.py` | invariant-bearing measured and non-engaging operation records |
| `native.py` | dispatch opaque motions into the native transactional replay owner |
| `replay.py` | ordered measure-before-deplete orchestration |
| `report.py` | report construction, aggregates, certification assertion, canonical bytes |
| `errors.py` | one named exception per audit failure mode |

The existing `compas_cgal.engagement` module remains unchanged until the new
path is validated. Once benchmark and downstream consumers use the new path,
the old function is labelled legacy diagnostic in documentation. Deletion or
renaming is a later user decision.

### Input contract

`EngagementAuditInput.build(...)` consumes:

```python
EngagementAuditInput.build(
    *,
    design_boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    engagement_cap: EngagementCap,
    depletion_policy: DepletionPolicy,
    operations: tuple[ToolpathOperation, ...],
    build_identity: BuildIdentity,
) -> EngagementAuditInput
```

The factory binds canonical rings, frame/unit-bearing physical parameters, the
ordered operation-stream digest, audit schema version, motion-certifier
versions, depletion policy, and native build identity. The cap encoding binds
both the authored binary64 angle and its native chord-ratio surrogate. The
native policy boundary recomputes that binary64 surrogate once and requires
bit-identical equality before exact injection. The input encoding also binds
the arc-surrogate, audit decision-contract, and depletion strategy
versions plus the ordered authenticated-operation digests. Arc-operation
identity therefore includes the exact chart-motion digest, not only authored
angles and a strategy label. It also binds a native request digest over the
exact-injected stock rings, immutable native policy, and ordered native-motion
digests. Empty operation streams raise
`EmptyToolpathAuditError`. Multi-depth or unsupported 3D motion raises a named
geometry error before stock mutation.

Digest domains have no generic raw-byte constructor or cross-domain retagger.
External input and authenticated-operation expectations enter through distinct
named ingress functions; generated request, motion, policy, decision,
depletion, lineage, and result digests arise only from their domain's canonical-
hash authority. All six
native motion carriers are private-constructor invariant classes. The request
factory accepts only their opaque Python values and retains their ordered typed
digests; no raw digest tuple or generic motion union crosses Python.

Stage 1 defines `BuildIdentity` and requires complete component-version and
32-byte source/lock digests at the audit boundary. Stage 4 replaces the
caller-supplied construction site with the deterministic build-generated
manifest and proves source mutation changes the digest. The type and audit
schema do not change between those stages.

The cut plane is never inferred. Every operation endpoint is classified
against its exact declared `cut_z`/`clearance_z` contract at the input seam.
Binary64 coordinates are injected into Epeck exactly; tolerance does not decide
whether a motion is cutting. COMPAS objects are mutable one-shot ingress only.
The factory retains immutable classified snapshots and their ordered source
digest, never caller-owned `ToolpathOperation` instances. Replay therefore has
no mutation window to reauthenticate.

### Operation result domain

Replace the ambiguous boolean in the authoritative path with a closed union:

```python
NativeMotionVerdict = Literal["certified", "cap_exceeded", "unresolved"]

OperationAudit = (
    MeasuredOperationAudit
    | PlungeOperationAudit
    | NonEngagingOperationAudit
)
```

`MeasuredOperationAudit` binds the native verdict, reporting maximum TEA,
station/event counts, pre- and post-motion stock lineage, native decision
digest, depletion-witness digest, and authenticated-operation digest.
`PlungeOperationAudit` binds its disk-depletion witness and lineage transition
without inventing a lateral TEA verdict. `NonEngagingOperationAudit` binds a
geometry-derived reason limited to vertical retract or clearance-plane
transport. It has no certification verdict, does not mutate stock, and cannot
enter compliant-motion counts.

A native-proved plunge remains a separate authenticated input operation because
its terminal cut-plane disk mutates stock even though it needs no lateral TEA
measurement. The authenticated carrier stores only opaque native geometry;
Python neither copies the endpoint nor reconstructs a curve radius. Circle and
arc native values retain their exact-injected guide radius for downstream
native consumers. Retracts and clearance transports remain non-mutating native
values and become `NonEngagingOperationAudit` only during replay/reporting.

Ingress radius is a finite unit-bearing observation, not a positivity
claim. Exact native sign owns acceptance. Every authenticated input-operation
carrier has a distinct versioned digest binding its stream index, source
operation digest, closed native classification tag, and opaque native-motion
digest.

A cut-height lateral motion always reaches a native certifier. Unsupported
geometry raises; it never produces a non-engaging record. A `RETRACT` label on
cutting geometry is contradictory input and raises
`ContradictoryOperationRoleError`.

### Native transaction and lineage

Python never supplies an authoritative `Stock2`. Native request identity
recomputes and verifies the stock/policy/motion digest from actual opaque native
values; the input digest separately accepts the verified request digest, so no
circular native authentication claim is made. One opaque, nonconstructible
native replay owner constructs stock from the canonical input rings, accepts
the input digest as an external typed root, binds the ordered authenticated-
operation digests, and only then derives
`StockLineageDigest = SHA256(stock-lineage-seed-v1 {input_digest,
verified_native_request_digest})`. Input-seeded means domain-separated
derivation, never retagging the input digest bytes. Each lateral call must match
the next bound identity and consumes one opaque native motion. It performs one
atomic transaction: decide from immutable pre-motion stock, clone exactly once,
apply validated depletion in-place to the trial for every returned verdict,
validate decision and depletion witnesses, derive the complete result and next
state, then perform one no-throw state swap. No allocation, hashing, validation,
or cursor update remains after the swap. Any exception exposes neither a result
nor partial report. Plunges use the same chronology; retracts and clearance
transports preserve lineage.

The native request value owns its canonical stock/policy identity and ordered
typed motion digests; only replay can inspect that retained sequence. The Python
request factory accepts the six opaque motion classes, never a generic motion
union or raw digest tuple. `begin_audit_replay(...)` accepts external bytes only
for the named input and authenticated-operation digest domains, then rebuilds
and matches the actual rings, opaque policy, and native request before state
exists. Native request digests are generated-only.

Replay's one-clone trial applicators do not clone or swap internally. Each
returns a validated nonconstructible depletion witness binding motion, policy,
exact trace or plunge disk, pre/post canonical stock-state digests, and strategy
identity. The authority swap occurs only after the complete decision, witness,
lineage, result, and next-state values exist. Controlled failure tests at every
stage prove stock, cursor, lineage, and finalization remain unchanged.

The native result boundary is a closed union of distinct read-only values for
lateral, plunge, and non-engaging chronologies. Reporting-only TEA and work
counts cannot affect the exact verdict or native decision digest. Every result
binds its authenticated-operation digest and adjacent lineage, and the report
validates the complete chain through the terminal lineage. Native finalization
fails unless every bound operation was consumed exactly once.

### Report semantics

`EngagementAuditReport.build(...)` derives, rather than accepts:

- certified motion count;
- proved cap-exceedance count;
- unresolved motion count;
- plunge-depletion count;
- non-engaging motion count;
- maximum reported TEA over measured motions;
- ordered operation-audit digest; and
- complete report digest bound to `EngagementAuditInput`.

`report.require_certified()` succeeds only when the stream is nonempty,
contains at least one measured lateral motion, and every measured verdict is
`certified`. It raises distinct `CapExceededToolpathError` and
`UnresolvedEngagementAuditError`. Reporting APIs may display all verdicts but
cannot reinterpret them.

### Exact arc surrogate, certification, and depletion

The authoritative partial-arc semantic is an explicit exact rational-chart
surrogate constructed once at the native classification seam. It retains exact
center, guide radius, radius-aligned base phase, ordered quarter-chart
intervals, direction, cut plane, canonical bytes, and digest. A versioned
binary64 angle-to-chart seam normalizes authored angles into `[0, tau)`, chooses
quadrants by exact comparison against exact-injected seam values, maps exact
seams structurally to parameter zero, and exact-injects
`tan(local_angle / 2)` only for non-seam interval endpoints. Rounding that
collapses distinct requested endpoints fails loudly. Every reconstructed phase
therefore lies exactly on the same Epeck guide circle. The canonical motion
binds its exact chart parameters and strategy version; Python observes only its
digest, never reconstructive fields.

One shared exact quarter-chart evaluator serves full circles, partial arcs,
certification, and depletion; duplicate formula implementations are forbidden.
The four frozen rational Pythagorean quarter charts used by full-circle event
motion are reused, but a partial arc owns clipped interval endpoints and seam
ownership. Its structural certificate proves endpoint and sample incidence,
ordered direction, declared complement, seam ownership, exact chord bounds,
and finite center count. The resulting full-radius disk union is a proved
under-cover because every disk center lies on the declared surrogate. Exact
`chord_bound < tool_radius` proves adjacent disks overlap; no stronger retained-
sliver bound is claimed without a separate theorem. The full-circle oracle's
cyclic four-chart coverage cannot certify a trimmed arc unchanged: partial arcs
require trimmed-domain endpoint/seam ownership and non-cyclic event ordering.
The legacy trigonometric `subtract_arc_sweep` is never evidence for, or
reachable from, the authoritative path.

Commit `73d5372` is theorem, falsifier, and refinement-strategy input, not
transferable authority. Its false `cap_certified` result conflates an observed
violation with proof exhaustion, its guarded station threshold is stricter than
the authored cap, and its double trigonometric centers do not describe the Task
3 surrogate. `CAP_EXCEEDED` therefore requires a replay-validated exact station
witness against the unguarded policy chord surrogate. Guarded-cap failure or an
inconclusive swept upper bound can only refine or return `UNRESOLVED`.

The native arc adapter and depleter consume the same opaque
`AuditArcMotion2`. Exact existential witnesses take precedence over incomplete
coverage. `CERTIFIED` requires complete exact coverage and exact
between-station closure of every clipped interval; until a general trimmed
closure exists, a sound safe-superset proof such as complete full-circle
authority may certify the subset, otherwise a witness-free valid arc is
`UNRESOLVED`. Together the adapter and depleter must:

- own the sole circular-motion acceptance path;
- use adaptive refinement with an absolute spatial floor and finite depth;
- preserve exact station predicates and exact between-station closure;
- return all three native verdicts;
- deplete certified, cap-exceeded, and unresolved motions through the same
  exact-on-surrogate disk-chain contract;
- refuse annular, machined, and spiral rib counterexamples;
- certify non-vacuous clear positive controls; and
- use double analytic quantities only to schedule refinement, never to decide;
- report refinement work without using it as a decision input;
- validate every proof and depletion trace before atomically swapping stock;
  and
- bind native decision, depletion witness, and post-lineage digests.

The fixed-density Python mirror remains available only to legacy diagnostics
during validation. No benchmark or certification consumer may call it after
the native path becomes authoritative.

### Stage 1 acceptance

- Empty, multi-depth, mislabeled, off-plane, ramped, and unsupported inputs fail
  with their named errors before any audit success object exists.
- Every cut-height line, arc, and circle yields one native three-way verdict.
- Every mutating operation advances one continuous input-seeded stock lineage;
  returned exceeded and unresolved motions still deplete.
- An unmeasured cut cannot be constructed through public or raw dataclass
  construction.
- False-certificate mutation tests kill both segment and arc paths.
- Exact partial-arc incidence, chart order, seam ownership, center bounds, and
  no-legacy-depletion controls pass for CW/CCW minor, major, and full-turn arc
  surrogates.
- Positive controls prove the certifier does not pass by refusing everything.
- The 12x8 generator-circle acceptance measurement is recorded with exact
  command, input digest, certified fraction, unresolved fraction, and runtime.
- Benchmarks consume only `EngagementAuditReport`.
- Strict mypy and Ruff pass for the new package and its consumers.

## Stage 2: replay and motion-theorem closure

### Replay certificate

`_replay_fresh_state(...)` already produces a validated `FreshReplayTrace` and
checks terminal traversal plus complete coverage. `replay_certificate(...)`
must retain that trace and construct a certificate instead of discarding it.

The authoritative `ReplayCertificate` is expanded into a content-addressed
domain record binding:

- input identity digest and complete build identity;
- ordered operation digest and per-operation identity chain;
- rebuilt MAT certificate and sampling-policy identity;
- fresh replay-trace digest;
- terminal stock boundary and lineage digests;
- terminal coverage-certificate digest;
- motion-oracle component versions; and
- replay schema/strategy version.

`ReplayCertificate.build(...)` validates all lengths, exact types, canonical
record structure, chain cardinality, terminal trace, and cross-digest
relations. `replay_certificate(...)` returns only the factory-built value.

Acceptance tests must include a genuinely terminal complete-path positive
control. Existing rejection-only fixtures remain and mutation tests alter each
bound component independently. A successful replay is then repeated in a
fresh process and must reproduce canonical bytes and digest.

### Swept-prefix theorem

The current fast swept-prefix path is not accepted as production authority
solely from `start_clear && cap_is_exact_pi`.

The implementation stage first writes a proof note that states the quantified
geometric theorem, all stock and motion preconditions, and the exact conclusion.
The native certificate binds evidence for every precondition. Independent
tests construct non-vacuous stock geometries across translation, rotation,
scale, multiple components, thin ribs, grazing boundaries, and reversed
source motion.

Until that proof and test package passes independent review, production
retrace routes through the generic event-exact segment certifier and the fast
swept-prefix implementation is research-only. If generic certification is too
slow, that is a measured product-gate failure, not authority to restore the
unproved fast path.

The latent empty-partition contradiction is repaired in the same stage:
empty segment strata are unresolved, and a certified verdict may never coexist
with `whole_rim_disposition == "unresolved"`.

### Independent acceptance side

For every certified segment or circle verdict, run a disjoint falsifier that
checks the exact `cap_exceeded` predicate, not only reported total TEA. The
falsifier cannot certify; it can only reject a candidate certificate. Native
TEA tests cover the C++ authorities directly, and property tests generate
adversarial event-corpus variants rather than only fat star-shaped polygons.

### Stage 2 acceptance

- At least one complete toolpath returns a `ReplayCertificate` outside a
  rejection context.
- Fresh-process replay is byte-identical.
- Every certificate field is mutation-killed by a named error.
- Terminal traversal and complete coverage remain mandatory.
- Empty event partitions fail closed.
- Swept-prefix production use is backed by the reviewed theorem gate or remains
  disabled; there is no ambiguous middle state.
- Direct native TEA tests and independent accept-side falsifiers pass.

## Stage 3: generator qualification and typed overrides

### Override ledger

Warnings are reporting only. Each regulated generator returns an immutable
`EngagementOverrideLedger` whose entries bind:

- operation/candidate identity;
- exact refusing predicate or unresolved certificate digest;
- forced-entry, forced-loop, or forced-advance reason;
- selected alternative identity;
- pre-motion stock lineage; and
- canonical record digest.

The public generator result contains the complete ledger. Forced advances are
not dropped at the public boundary. A consumer claiming cap compliance must
call `require_no_overrides()` and then run the independent engagement audit.

### Search authority

Non-monotonic engagement forbids bisection as the admissibility authority. Add
one shared ladder-search module with typed spatial/radius inputs and explicit
finite enumeration. Migrate generators one at a time with characterization and
counterexample tests. The duplicate private bisection and downward scan remain
until the new authority is validated; their later removal requires explicit
permission.

The engagement generators adopt frame/unit-bearing types at their public and
shared-private boundaries. Point, tangent, radius, advance, and cap do not share
one scalar carrier. `_GuideStation`-shaped data is decomposed into typed point,
direction, and radius values.

### Quality gates

Repair tangent-continuity tests so their non-vacuity floor counts actual
comparisons. Structural Line/Arc alternation cannot exempt a transition; any
bridge exemption must prove its expected geometry. Early returns on absent arcs
become loud failed preconditions in fixtures that promise arcs.

The corpus runs the regulated generators as the primary subjects. The
unregulated generator remains a named baseline. Add a second cap that separates
generator behaviour. Unsupported true helical motion fails loud until a typed
3D depletion contract exists.

### Stage 3 acceptance

- Every forced loop and forced advance survives the public result boundary.
- Warning filters cannot change audit or compliance results.
- Ladder search and exhaustive small-oracle enumeration agree.
- No regulated consumer uses the monotonic bisection path.
- Tangent tests assert a comparison-count floor and validate every exemption
  geometrically.
- Corpus identities state generator, cap, input, build, and override-ledger
  digests.
- Product-quality failures remain honest failures in `quality-gate`.

## Stage 4: enforcement, build identity, and evidence

### Python and type-system enforcement

Raise all declared floors and targets to Python 3.12. Fix the six existing
strict-mypy errors by narrowing unions from their discriminating domain
preimages; do not silence them with casts or ignores. Remove obsolete ignores.

Strict gates cover:

- the complete adaptive package;
- the authoritative engagement-audit package;
- benchmark schemas and runners; and
- typed consumer-contract fixtures.

### Native build identity

Production constructs `BuildIdentity`; tests are not its only caller. The
build embeds a SHA-256 manifest over:

- every native source/header in the exact audit call graph;
- binding stubs and canonical-schema constants;
- CMake configuration affecting those targets;
- compiler identity and flags;
- CGAL, CORE/GMP, nanobind, Python ABI, and component versions; and
- the Pixi lock digest.

The manifest is generated deterministically by a Pixi-owned build task and
embedded by CMake. Runtime bindings expose the canonical manifest bytes and
digest. Python `BuildIdentity.build(...)` verifies the bytes and binds them into
audit and replay input identities. A source mutation without a digest change is
an acceptance-test failure.

### CI topology

CI uses Pixi tasks exclusively. Pull requests and branch pushes run without a
`main`-only trigger. Concurrency cancels superseded runs.

Separate jobs provide independently readable evidence:

1. native build and direct native TEA tests;
2. adaptive exact tests with `-n auto`;
3. authoritative audit tests with `-n auto`;
4. strict mypy and Ruff;
5. schema/artifact round trips;
6. benchmark-instrument validity tests; and
7. strict MkDocs.

Dead schema and mutation tasks are implemented or removed from the public task
surface. Every task must fail when it collects zero tests. The product
`quality-gate` is not hidden inside the regression baseline.

### Evidence artifacts

One repository command regenerates a complete evidence bundle containing:

- canonical input and build identities;
- schema-versioned raw measurement JSON;
- audit verdict counts and operation-level digests;
- replay-certificate digest;
- runtime distributions with warm/cold protocol;
- regulated and baseline generator identities;
- override-ledger counts;
- unresolved rate on the committed third-party corpus; and
- derived tables and figures.

JSON round-trips through `json.loads` and schema validation. Every figure reads
the committed raw data. Held digitisation data and its provenance are committed
before any comparison figure is regenerated. Documentation numbers are
generated from the same artifact; a stale headline fails a docs evidence test.

### Stage 4 acceptance

- All CI jobs run on the programme branch and are green.
- Zero-test, invalid CLI, and stale-schema false greens are impossible.
- Runtime audit/replay artifacts expose the actual native source-tree identity.
- One clean-checkout command reproduces byte-identical deterministic evidence
  and equivalent statistical summaries where wall-clock bytes are excluded.
- Documentation states exact maturity, known-red quality gates, unresolved
  rates, runtime gap, and unsupported claims without promotional inference.

## Disposition of remaining review findings

The programme includes findings that can change proof truth, public truth, or
release enforcement. Other maintainability findings are retained in a bounded
post-convergence ledger:

| Finding | Programme disposition |
| --- | --- |
| private engagement module used as shared library | Stage 1 new focused package; legacy deletion later |
| large adaptive authority modules | decompose only files changed by Stages 1-2, by responsibility |
| production MAT spike fixtures | post-convergence extraction unless reconciliation touches the file |
| oversized nanobind module body | extract only binding-leaked core invariants needed by Stage 2 |
| factory validation concentrated in `__post_init__` | repair only factories changed by this programme |
| prohibited root `__all__` | isolated compatibility change in Stage 4 |
| undocumented medial-axis public symbols | update the coherent stage page and touched public APIs |
| no certifier Hypothesis coverage | Stage 2 adversarial event-corpus strategies |
| missing circle certified cross-check | Stage 2 independent accept-side falsifier |
| missing Figure 6 data/generator | Stage 4 evidence bundle |

This prevents aesthetic cleanup from delaying the truth boundary while keeping
every review finding explicitly owned.

## SDD plan set and review gates

After this design is approved, create these plans:

1. `2026-08-23-auditor-convergence-p0-reconciliation.md`
2. `2026-08-23-auditor-convergence-p1-truthful-audit.md`
3. `2026-08-23-auditor-convergence-p2-replay-theorem.md`
4. `2026-08-23-auditor-convergence-p3-generator-qualification.md`
5. `2026-08-23-auditor-convergence-p4-enforcement-evidence.md`

Each plan uses TDD, 2-5 minute executable steps, exact file/symbol interfaces,
focused commits, and its own independent review gate. P1 cannot begin until P0
is approved; P2 cannot begin until P1 is green; P3 and P4 may begin only after
P2 because both consume the authoritative audit/replay contracts. P3 and P4
may then execute independently in separate worktrees before a final linear
rebase onto the convergence branch.

Maintain one cold-start ledger at
`.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`. It records current
plan/task, branch/worktree/SHA, dirty state, last verified commands, accepted
review commits, open blockers, and next exact command. It never substitutes for
durable MkDocs documentation or committed evidence.

## Final release decision

After P4, run a fresh whole-branch adversarial review and classify three theses
independently:

| Thesis | Required decision evidence |
| --- | --- |
| quality instrument | truthful third-party verdicts, unresolved rate, replay identity, reproducible artifact, CI |
| regulated generator | quality-gate results, override-free fraction, coverage/gouge evidence, runtime |
| exact certifier | native false-certificate resistance, acceptance controls, exact replay, runtime distribution |

Each thesis receives `release candidate`, `research prototype`, or `rejected at
current evidence`. Strength in one thesis cannot compensate for a failed gate
in another.
