# Exact Zero-Guide Link Candidate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prove constant tool-clearance MAT spans exactly, machine them with a
typed advancing segment, and reconstruct the same decision under fresh replay.

**Architecture:** The retained native MAT owner classifies zero-guide edges by
an exact polynomial identity and exposes canonical proof records. Python
projects those records into graph-lifetime authority, then uses separate
candidate, operation, and transaction variants so advancing segments cannot be
confused with hold links. Generator and replay dispatch over closed unions and
authenticate every record against the rebuilt MAT.

**Tech Stack:** C++20, CGAL 6.0.1 exact kernels, CORE exact rationals, nanobind,
Python 3.12, COMPAS frame-tagged geometry, CCAN, SHA-256, pytest-xdist,
pytest-testmon, Ruff, strict mypy, strict MkDocs.

## Global Constraints

- Exact zero-guide membership is the coefficient identity
  `C_e(t) - r² == 0`; no float comparison, epsilon, endpoint test, sample test,
  caught exception, or empty-family retry may decide it.
- The fixed 20-position numerical MAT projection remains unchanged. New exact
  evidence is an immutable property of `SegmentSiteMedialAxis`.
- `LinkSegmentOperation` remains a hold operation paired with a full circle.
  `AdvanceSegmentOperation` is the only link-only traversal advance.
- `CandidateTransaction` retains its two-witness chronology.
  `ZeroGuideLinkTransaction` has exactly one advancing-segment witness.
- Existing circle candidate/evaluation/replay canonical bytes do not change.
- Every new failure mode raises a named exception and every proof-boundary test
  has a Google-style contextual docstring.
- Every pytest command uses `-n auto`; affected tests use `--testmon`.
- Run Ruff before committing Python and run strict mypy for every type-system
  claim.
- Update `docs/segment_site_mat.md`, `docs/continuous_engagement.md`, the
  causal traversal plan, and the Held–Pfeiffer comparison at each coherent
  stage.
- Commit with author and committer
  `Jelle Feringa <jelleferinga@gmail.com>`, cancel in-progress branch CI before
  each push, and verify the remote SHA after each push.

---

### Task 1: Native Exact Zero-Guide Proof Inventory

**Files**

- Create: `src/segment_site_zero_guide.h`
- Create: `src/segment_site_zero_guide.cpp`
- Modify: `src/segment_site_mat_certificate.h`
- Modify: `src/segment_site_mat_certificate.cpp`
- Modify: `src/segment_site_mat_bundle.h`
- Modify: `src/segment_site_mat_bundle.cpp`
- Modify: `src/medial_axis_2.cpp`
- Modify: `src/compas_cgal/_medial_axis_2.pyi`
- Modify: `CMakeLists.txt`
- Test: `tests/adaptive/test_medial_axis.py`

**Interfaces**

- Consumes:
  `MatClearanceProfileGraph2`,
  `MatCertificateV1`,
  `MatToolRadiusMm2::squared_exact()`.
- Produces:
  `MatZeroGuideRecordV1`,
  `canonical_mat_clearance_profile_v1(...)`,
  `certified_mat_zero_guide_records_v1(...)`,
  `validate_mat_zero_guide_records_v1(...)`,
  `SegmentSiteMedialAxis.zero_guide_records`,
  `SegmentSiteMedialAxis.validate_zero_guide_records(...)`.

- [ ] **Step 1: Write the failing native-owner tests**

Add radius-1 and radius-0.5 owner tests with proof context:

```python
def test_native_owner_proves_complete_radius_one_zero_guide_inventory() -> None:
    """The two width-2 arms have clearance square exactly equal to 1."""
    owner = _owned_mat(tool_radius=1.0)

    assert len(owner.zero_guide_records) == 2
    assert owner.zero_guide_records == tuple(sorted(owner.zero_guide_records))
    assert owner.validate_zero_guide_records(
        owner.projection[19],
        owner.zero_guide_records,
    ) == tuple(edge_id for edge_id, _ in owner.zero_guide_records)


def test_native_owner_does_not_classify_positive_guide_lines_as_zero() -> None:
    """Radius 0.5 leaves positive guide radius on the same width-2 arms."""
    owner = _owned_mat(tool_radius=0.5)

    assert owner.zero_guide_records == ()
```

Require repeat builds and reversed input to yield the same sorted `(edge_id,
record_bytes)` pairs. Change one record byte, delete a record, duplicate an
edge ID with foreign bytes, and pass a foreign MAT certificate. Require the
specific native binding errors `MissingMatZeroGuideEdgeError`,
`DuplicateMatZeroGuideEdgeError`,
`MismatchedMatZeroGuideRecordError`, or
`InvalidMatCertificateReplayError`, never a generic construction error. The
radius-1 parabolic edges touch tool clearance only at an endpoint but remain
absent because their complete polynomial is nonconstant. A
`numpy.nextafter(1.0, 0.0)` radius also emits no records, proving binary64
near-equality is not zero-guide authority.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest tests/adaptive/test_medial_axis.py \
  -n auto --testmon -q
```

Expected: collection or attribute failures because the owner exposes no
zero-guide proof inventory.

- [ ] **Step 3: Define the exact native record**

Declare one value type with validated construction and canonical bytes:

```cpp
class MatZeroGuideRecordV1 {
public:
  static MatZeroGuideRecordV1 build(
      const std::string &mat_certificate_digest,
      const std::string &center_domain_digest,
      const MatClearanceEdgeProfile2 &profile,
      const CORE::BigRat &tool_radius_squared);

  const std::string &edge_id() const noexcept;
  const std::string &canonical_bytes() const noexcept;
};
```

`build` requires a constant coefficient vector of length one whose coefficient
equals `tool_radius_squared`. Canonical bytes bind schema/strategy, MAT digest,
center-domain digest, exact tool-radius square, the identity verdict, and the
authoritative `canonical_mat_clearance_profile_v1(profile)` record. Promote the
existing MAT-certificate profile encoder to that public function instead of
duplicating endpoint/root/site/coefficient encoding. The nested record already
binds edge ID, defining sites, both exact endpoint root IDs, and all clearance
coefficients. Materialize concrete exact-number values before returning; never
return a lazy CORE expression.

- [ ] **Step 4: Build and verify the complete inventory**

Implement classification as one pass over certificate profiles:

```cpp
for (const MatClearanceEdgeProfile2 &profile : exact.profiles()) {
  const RationalPolynomial &coefficients =
      profile.squared_clearance().coefficients();
  if (coefficients.size() == 1 &&
      coefficients.front() == tool_radius_squared) {
    records.push_back(MatZeroGuideRecordV1::build(
        mat_digest, center_domain_digest, profile, tool_radius_squared));
  }
}
```

Use the repository's exact polynomial and certificate accessors rather than
text parsing. Sort by exact `edge_id`, reject duplicate IDs, and have the
explicit verifier replay the MAT certificate, rebuild the full inventory, and
compare edge IDs and canonical bytes exactly. An absent positive-guide edge is
normal; an absent proved zero-guide record is an error.

- [ ] **Step 5: Retain and bind the records**

Store `MatToolRadiusMm2` plus the records in `SegmentSiteMatBundle2` at graph
lifetime and expose one read-only tuple:

```python
owner.zero_guide_records: tuple[tuple[bytes, bytes], ...]
```

Add:

```python
owner.validate_zero_guide_records(
    mat_certificate: bytes,
    records: tuple[tuple[bytes, bytes], ...],
) -> tuple[bytes, ...]
```

Translate only the named native proof errors. Do not add fields to
`projection_tuple`.

- [ ] **Step 6: Run GREEN and native regressions**

```bash
pixi run pytest tests/adaptive/test_medial_axis.py \
  -n auto --testmon -q
pixi run task9-mat-compile-gate
git diff --check
```

Expected: exact inventory tests pass; the existing fixed projection and native
Task-9 gates remain unchanged.

- [ ] **Step 7: Document, commit, and push**

Update the native proof inventory and exact polynomial theorem in
`docs/segment_site_mat.md`, then run:

```bash
pixi run docs
git add CMakeLists.txt src/segment_site_zero_guide.h \
  src/segment_site_zero_guide.cpp src/segment_site_mat_certificate.h \
  src/segment_site_mat_certificate.cpp src/segment_site_mat_bundle.h \
  src/segment_site_mat_bundle.cpp src/medial_axis_2.cpp \
  src/compas_cgal/_medial_axis_2.pyi tests/adaptive/test_medial_axis.py \
  docs/segment_site_mat.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(mat): prove zero-guide runs"
```

---

### Task 2: Typed Python Proof Projection

**Files**

- Modify: `src/compas_cgal/adaptive/medial_axis.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_medial_axis.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- Consumes:
  `SegmentSiteMedialAxis.zero_guide_records`,
  `MedialAxis.proof.mat_certificate`.
- Produces:
  `MatZeroGuideRun.build(...)`,
  `MatZeroGuideInventory.build(...)`,
  `MedialAxis.zero_guide_inventory`,
  `MedialAxis.zero_guide_run_by_edge_id`.

- [ ] **Step 1: Write RED projection and mutation tests**

Require immutable typed projection and exact owner cross-validation:

```python
def test_typed_axis_projects_owned_zero_guide_runs() -> None:
    """Python preserves native proof bytes instead of reclassifying floats."""
    axis = _typed_axis(tool_radius=1.0)

    assert tuple(run.edge_id for run in axis.zero_guide_inventory.runs) == tuple(
        edge_id for edge_id, _ in axis.native_owner.zero_guide_records
    )
    assert {
        edge_id: run.native_certificate
        for edge_id, run in axis.zero_guide_run_by_edge_id.items()
    } == dict(axis.native_owner.zero_guide_records)
```

Use `dataclasses.replace` to cross-wire edge ID, certificate bytes, MAT digest,
and native owner. Require `InvalidZeroGuideCertificateError`. Add strict-mypy
consumer assertions for immutable bytes and `MatEdgeId` keys.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest tests/adaptive/test_medial_axis.py \
  -n auto --testmon -q
pixi run types-adaptive
```

Expected: missing Python proof types and axis properties.

- [ ] **Step 3: Implement focused proof types**

Define these graph-lifetime proof types in `medial_axis.py` beside
`MatEdgeId`, before `MatProof`, so no module needs a duplicate identity type or
circular import:

```python
@dataclass(frozen=True)
class MatZeroGuideRun:
    edge_id: MatEdgeId
    mat_certificate_digest: IdentityDigest
    native_certificate: bytes

    @classmethod
    def build(
        cls,
        *,
        edge_id: MatEdgeId,
        mat_certificate: bytes,
        native_certificate: bytes,
    ) -> Self: ...


@dataclass(frozen=True)
class MatZeroGuideInventory:
    runs: tuple[MatZeroGuideRun, ...]
    mat_certificate_digest: IdentityDigest
```

Factories validate exact concrete types, nonempty canonical bytes, sorted
unique edge IDs, and one common MAT digest. The inventory owns only proof
records; `MedialAxis` validates edge existence and native byte equality.

- [ ] **Step 4: Project once during `MedialAxis.build`**

Extend `MatProof`, not the 20-field projection, with the immutable inventory:

```python
proof = MatProof(
    center_domain_digest=center_domain_digest,
    mat_certificate=mat_certificate,
    native_owner=native_owner,
    zero_guide_inventory=MatZeroGuideInventory.build(...),
)
```

Require `native_owner.validate_zero_guide_records(...)` to return the exact
same edge sequence before constructing runs. Add read-only convenience
properties on `MedialAxis`; do not cache a second mutable dictionary. Require
`wc -l src/compas_cgal/adaptive/medial_axis.py` to remain below 1,000 after the
minimal implementation.

- [ ] **Step 5: Run GREEN and type gates**

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_medial_axis.py \
  -n auto --testmon -q
pixi run types-adaptive
```

- [ ] **Step 6: Document, commit, and push**

Document native-to-Python ownership in `docs/segment_site_mat.md`, then:

```bash
git add src/compas_cgal/adaptive/medial_axis.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_medial_axis.py \
  tests/adaptive/typecheck/consumer_contract.py \
  docs/segment_site_mat.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): project zero-guide proof"
```

---

### Task 3: Finite-Lattice Zero-Guide Candidates

**Files**

- Modify: `src/compas_cgal/adaptive/candidates.py`
- Modify: `src/compas_cgal/adaptive/policy.py`
- Modify: `src/compas_cgal/adaptive/traversal.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_candidates.py`
- Modify: `tests/adaptive/test_policy.py`
- Modify: `tests/adaptive/test_traversal.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- Consumes:
  `MiddleCurveSpan`,
  `MatZeroGuideRun`,
  `_refined_spatial_values(...)`,
  `_middle_point(...)`,
  `_cursor_after(...)`.
- Produces:
  `ZeroGuideLinkCandidate.build(...)`,
  `enumerate_zero_guide_link_candidates(...)`,
  `TraversalCandidate`,
  a `DerivedCandidateCursor` accepting either exact candidate variant.

- [ ] **Step 1: Write RED candidate contracts**

Build forward and reverse spans on the two radius-1 zero-guide edges and
require:

```python
candidates = enumerate_zero_guide_link_candidates(
    span=span,
    policy=policy,
    neck_scope=scope,
    effective_cap_decision=cap,
    makes_cursor_terminal_at_limit=True,
)
assert candidates
assert all(candidate.zero_guide_run.edge_id == span.edge.identity for candidate in candidates)
assert candidates[0].spatial_progress == max(candidate.spatial_progress for candidate in candidates)
assert candidates[0].traversal_decision.makes_cursor_terminal
assert not hasattr(candidates[0], "guide_radius")
assert not hasattr(candidates[0], "phase_index")
assert not hasattr(candidates[0], "motion")
```

Require repeatable order, forward/reverse identity separation, and failures for
foreign run, positive-guide edge, wrong target, progress, levels, cursor limit,
scope, cap, and traversal identity. Candidate-record failures raise
`InvalidZeroGuideCandidateError`; asking the enumerator to classify an
unproved edge raises `UncertifiedZeroGuideEdgeError`.

Add a policy test proving that an exact zero squared radius is valid and ranks
below a positive radius at equal progress:

```python
zero_key = CandidateOrderKey.build(
    progress=ExactMillimetre(Fraction(2)),
    squared_radius=SquaredMillimetre(Fraction(0)),
    canonical_identity=b"zero",
)
zero_candidate = _Candidate("zero", zero_key)
circle_candidate = _Candidate(
    "circle",
    CandidateOrderKey.build(
        progress=ExactMillimetre(Fraction(2)),
        squared_radius=SquaredMillimetre(Fraction(1)),
        canonical_identity=b"circle",
    ),
)
ordered = policy.order_candidates(
    (zero_candidate, circle_candidate),
    key=lambda candidate: candidate.order_key,
)
assert ordered == (circle_candidate, zero_candidate)
```

- [ ] **Step 2: Run RED**

```bash
pixi run pytest tests/adaptive/test_candidates.py \
  tests/adaptive/test_policy.py \
  tests/adaptive/test_traversal.py \
  -n auto --testmon -q
```

Expected: missing zero-guide candidate type and closed union.

- [ ] **Step 3: Implement the candidate record**

Define the frozen content-addressed type beside the existing candidate domain
types so the closed union and derived cursor have no import cycle:

```python
ZeroGuideLinkCandidateId = NewType("ZeroGuideLinkCandidateId", bytes)


@dataclass(frozen=True)
class ZeroGuideLinkCandidate:
    identity: ZeroGuideLinkCandidateId
    zero_guide_run: MatZeroGuideRun
    policy: CandidatePolicy
    spatial_progress: ExactMillimetre
    spatial_levels: tuple[int, ...]
    target: Point2[WorldXY]
    cursor_limit_identity: CursorIdentity
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision
```

Its canonical bytes use a distinct `zero-guide-link-candidate-v1` tag and bind
every field. Its order key is:

```python
CandidateOrderKey.build(
    progress=self.spatial_progress,
    squared_radius=SquaredMillimetre(Fraction(0)),
    canonical_identity=bytes(self.identity),
)
```

Change `CandidateOrderKey` validation from positive to non-negative squared
radius. The candidate contains no circle-derived field.

- [ ] **Step 4: Enumerate only owned exact runs**

Reuse the existing spatial lattice and interpolation helpers without copying
their deciding logic. At function entry:

```python
run = span.axis.zero_guide_run_by_edge_id.get(span.edge.identity)
if run is None:
    raise UncertifiedZeroGuideEdgeError(
        "zero-guide enumeration requires native proof for the complete span edge."
    )
```

Build one candidate per refined spatial value, including the terminal endpoint.
There is no radius or phase loop.

- [ ] **Step 5: Extend traversal with a closed union**

Define:

```python
TraversalCandidate: TypeAlias = MiddleCurveCandidate | ZeroGuideLinkCandidate
```

Update `DerivedCandidateCursor` and cursor-lineage validation through explicit
exact-type branches. Both variants must prove progress, target/middle point,
cursor limit, and traversal decision against the same `MiddleCurveSpan`;
circle candidates additionally prove their generator-site ownership.
`ExhaustedCandidateCursor` remains circle-only because only the
positive-radius lattice can exhaust before a native endpoint. Neither
candidate may enter through duck typing. Preserve existing circle candidate
bytes and behavior. Require `wc -l src/compas_cgal/adaptive/candidates.py` to
remain below 1,200 after the minimal implementation; if it cannot, stop for a
separate module-boundary design instead of introducing a circular import.

- [ ] **Step 6: Run GREEN, lint, and types**

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_candidates.py \
  tests/adaptive/test_policy.py \
  tests/adaptive/test_traversal.py \
  -n auto --testmon -q
pixi run types-adaptive
```

- [ ] **Step 7: Document, commit, and push**

Update the finite-family section of `docs/segment_site_mat.md`, then:

```bash
git add src/compas_cgal/adaptive/candidates.py \
  src/compas_cgal/adaptive/policy.py \
  src/compas_cgal/adaptive/traversal.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_candidates.py tests/adaptive/test_policy.py \
  tests/adaptive/test_traversal.py \
  tests/adaptive/typecheck/consumer_contract.py \
  docs/segment_site_mat.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): enumerate zero-guide links"
```

---

### Task 4: Advancing Segment Operation and State Grammar

**Files**

- Modify: `src/compas_cgal/adaptive/operation.py`
- Modify: `src/compas_cgal/adaptive/generation_state.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_operation.py`
- Modify: `tests/adaptive/test_generation_state.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- Consumes:
  `ExactSegmentMotion`,
  `AdvanceTraversalDecision`,
  `EffectiveCapDecision`,
  `NeckScope`.
- Produces:
  `AdvanceSegmentOperation.build(...)`,
  extended `CanonicalOperation`.

- [ ] **Step 1: Write RED operation identity tests**

Require:

```python
operation = AdvanceSegmentOperation.build(
    motion=segment,
    cut_z=cut_z,
    neck_scope=scope,
    effective_cap_decision=cap,
    traversal_decision=advance,
)
assert require_canonical_record(operation.canonical_bytes) == (
    operation.canonical_bytes
)
assert b"advance-segment-v1" in operation.canonical_bytes
assert operation != LinkSegmentOperation.build(
    motion=segment,
    cut_z=cut_z,
    neck_scope=scope,
    effective_cap_decision=cap,
    traversal_decision=hold,
)
```

Reject a hold decision, mismatched scope/cap, foreign motion, subclass, and
changed endpoint. Every invalid construction raises
`InvalidAdvanceSegmentOperationError`, not the generic operation identity
error. Add mypy exhaustiveness for all five operation variants.

- [ ] **Step 2: Write RED chronology tests**

Require `GenerationState` to accept:

```text
entry circle, advance segment
entry circle, hold link, circle, advance segment
entry circle, advance segment, hold link, circle
```

Reject a dangling hold link, hold link followed by advancing segment,
advancing segment followed by a circle claiming the same cursor, and final
phase/traversal state not matching the last advancing operation.

- [ ] **Step 3: Run RED**

```bash
pixi run pytest tests/adaptive/test_operation.py \
  tests/adaptive/test_generation_state.py \
  -n auto --testmon -q
```

- [ ] **Step 4: Implement the distinct operation**

Add:

```python
@dataclass(frozen=True)
class AdvanceSegmentOperation:
    motion: ExactSegmentMotion
    cut_z: CutZ
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision
```

Use `advance-segment-v1`; do not add a mode field to `LinkSegmentOperation`.
Extend `CanonicalOperation` explicitly.

- [ ] **Step 5: Implement the closed state grammar**

Parse lateral operations with explicit cases:

```python
if type(operation) is LinkSegmentOperation:
    require_next_circle_pair(...)
elif type(operation) is CutFullCircleOperation:
    require_not_unpaired(...)
elif type(operation) is AdvanceSegmentOperation:
    advance_phase_to(operation.motion.end)
else:
    raise InvalidGenerationStateError(...)
```

Require exactly one traversal advance per accepted lateral chronology and
allow the state to end on either a circle or advancing segment.

- [ ] **Step 6: Run GREEN, lint, and types**

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_operation.py \
  tests/adaptive/test_generation_state.py \
  -n auto --testmon -q
pixi run types-adaptive
```

- [ ] **Step 7: Document, commit, and push**

Update the operation/replay grammar in `docs/continuous_engagement.md`, then:

```bash
git add src/compas_cgal/adaptive/operation.py \
  src/compas_cgal/adaptive/generation_state.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_operation.py \
  tests/adaptive/test_generation_state.py \
  tests/adaptive/typecheck/consumer_contract.py \
  docs/continuous_engagement.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): add advancing segment"
```

---

### Task 5: Atomic Zero-Guide Evaluation and Commit

**Files**

- Modify: `src/compas_cgal/adaptive/transaction.py`
- Modify: `src/compas_cgal/adaptive/replay_trace.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_transaction.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- Consumes:
  `ZeroGuideLinkCandidate`,
  `AdvanceSegmentOperation`,
  `ReplayLateralWitness`,
  existing short-lived proof authorities.
- Produces:
  `ZeroGuideLinkTransaction.build(...)`,
  `AcceptedCandidateTransaction`,
  `CandidateEvaluator.evaluate_zero_guide_from_cursor(...)`,
  `CandidateEvaluator.commit_zero_guide_from_cursor(...)`.

- [ ] **Step 1: Write RED atomicity and chronology tests**

Use a real radius-1 arm segment and require:

```python
transaction = evaluator.evaluate_zero_guide_from_cursor(
    state,
    traversal_before,
    candidate,
)
assert type(transaction.segment_witness.operation) is AdvanceSegmentOperation
assert transaction.segment_witness.motion_witness.stock_lineage_digest == (
    stock_lineage_digest(state.stock.lineage)
)
child = evaluator.commit_zero_guide_from_cursor(
    state,
    traversal_before,
    transaction,
)
assert child.digest == transaction.result_state_digest
assert child.operations[-1] == transaction.segment_witness.operation
```

Instrument containment, TEA certification, depletion, and coverage to prove
that exact order. Force each proof stage to fail and assert the immutable
parent stock, coverage, traversal cursor, operations, and passage remain
unchanged.

- [ ] **Step 2: Write RED mutation and stale-parent tests**

Cross-wire parent digest, candidate run, target, cap, scope, operation,
pre-depletion lineage, coverage lineage, traversal result, passage result, and
child digest. Require `InvalidZeroGuideTransactionError`. Commit against a
different parent or cursor and require `StaleZeroGuideTransactionError`.

- [ ] **Step 3: Run RED**

```bash
pixi run pytest tests/adaptive/test_transaction.py \
  -n auto --testmon -q
```

- [ ] **Step 4: Implement the one-witness transaction**

Define the second transaction beside `CandidateTransaction` so both variants
share the one authoritative stock-lineage and passage-validation helpers
without a circular import:

```python
@dataclass(frozen=True)
class ZeroGuideLinkTransaction:
    parent_state_digest: IdentityDigest
    candidate: ZeroGuideLinkCandidate
    segment_witness: ReplayLateralWitness
    traversal_after: TraversalCursorState
    passage_after: NeckPassage | None
    result_state_digest: IdentityDigest
```

Canonicalize under `zero-guide-link-transaction-v1`. Validate the operation
against every candidate field, require the segment endpoint to equal target,
prove certify-before-deplete lineage, and prove one matching coverage update.
Extend `ReplayLateralWitness` with the exact `AdvanceSegmentOperation` variant;
like a hold link it requires `OperationType.LINK` and a
`SegmentContainmentCertificate`.

- [ ] **Step 5: Add explicit evaluator methods**

Implement:

```python
def evaluate_zero_guide_from_cursor(
    self,
    state: GenerationState,
    traversal_before: TraversalCursorState,
    candidate: ZeroGuideLinkCandidate,
) -> ZeroGuideLinkTransaction: ...

def commit_zero_guide_from_cursor(
    self,
    state: GenerationState,
    traversal_before: TraversalCursorState,
    transaction: ZeroGuideLinkTransaction,
) -> GenerationState: ...
```

Evaluation constructs one segment from `state.phase_point` to
`candidate.target`, then performs containment, TEA certification,
depletion, coverage, traversal/passage advancement, child state construction,
and transaction construction. Commit independently repeats the same complete
path and compares transaction bytes plus child digest. Do not dispatch both
calling conventions inside one method body.

- [ ] **Step 6: Define the accepted closed union**

```python
AcceptedCandidateTransaction: TypeAlias = (
    CandidateTransaction | ZeroGuideLinkTransaction
)
```

Keep existing circle evaluator methods unchanged and add strict-mypy narrowing
tests. Require `wc -l src/compas_cgal/adaptive/transaction.py` to remain below
1,200 after the complete minimal implementation.

- [ ] **Step 7: Run GREEN, lint, and types**

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_transaction.py \
  -n auto --testmon -q
pixi run types-adaptive
```

- [ ] **Step 8: Document, commit, and push**

Update the atomicity section in `docs/continuous_engagement.md`, then:

```bash
git add src/compas_cgal/adaptive/transaction.py \
  src/compas_cgal/adaptive/replay_trace.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_transaction.py \
  tests/adaptive/typecheck/consumer_contract.py \
  docs/continuous_engagement.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): commit zero-guide links"
```

---

### Task 6: Generator Dispatch Through Task 13F Route 1

**Files**

- Modify: `src/compas_cgal/adaptive/generator.py`
- Modify: `tests/adaptive/test_generator.py`
- Modify: `src/compas_cgal/adaptive/errors.py`

**Interfaces**

- Consumes:
  `TraversalCandidate`,
  `AcceptedCandidateTransaction`,
  `MedialAxis.zero_guide_run_by_edge_id`.
- Produces:
  exact variant selection in `advance_active_candidate_family(...)`.

- [ ] **Step 1: Write the Task 13F RED route-1 test**

Extend the authenticated fixture from the already committed route-0 child:

```python
def test_task13f_route_one_accepts_exact_zero_guide_link(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The width-2 arm is machined by its exact MAT centerline segment."""
    fixture = Task13FFixture.build()
    physical_zero, traversal_zero, _ = advance_active_candidate_family(
        physical=fixture.physical,
        traversal=fixture.traversal,
        evaluator=fixture.evaluator,
    )
    traversal_one = traversal_zero.activate_next()
    attempts: list[TraversalCandidate] = []
    real_evaluate = CandidateEvaluator.evaluate_zero_guide_from_cursor

    def track_trial(
        self: CandidateEvaluator,
        state: GenerationState,
        traversal_before: TraversalCursorState,
        candidate: ZeroGuideLinkCandidate,
    ) -> ZeroGuideLinkTransaction:
        attempts.append(candidate)
        return real_evaluate(self, state, traversal_before, candidate)

    monkeypatch.setattr(
        CandidateEvaluator,
        "evaluate_zero_guide_from_cursor",
        track_trial,
    )

    physical_after, traversal_after, commit = advance_active_candidate_family(
        physical=physical_zero,
        traversal=traversal_one,
        evaluator=fixture.evaluator,
    )

    assert attempts
    assert all(type(candidate) is ZeroGuideLinkCandidate for candidate in attempts)
    assert type(commit.transaction) is ZeroGuideLinkTransaction
    assert physical_after.operations[-1].motion.end == physical_after.phase_point
    assert traversal_after.active_cursor != traversal_one.active_cursor
```

Monkeypatch `CandidateEvaluator.evaluate_zero_guide_from_cursor` to append each
candidate to `attempts` before calling the real method. This instruments the
existing calling convention instead of adding a callback API. Require the
first accepted winner to stop dispatch with no lower-ranked evaluation, and
require route 0 to remain the exact same fourth-circle winner.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest tests/adaptive/test_generator.py \
  -n auto --testmon -q
```

Expected: route 1 still exhausts with `attempts=0`.

- [ ] **Step 3: Implement exact one-time variant selection**

At the active span:

```python
run = traversal.authority.axis.zero_guide_run_by_edge_id.get(
    span.edge.identity
)
if run is None:
    family = enumerate_middle_curve_candidates(...)
else:
    family = enumerate_zero_guide_link_candidates(...)
```

Before enumeration, require the run to be byte-identical to the inventory
owned by `traversal.authority.axis`; retain the existing evaluator/MAT
authority cross-check. Extend `TraversalCommit`, family validation, evaluation,
global advance, physical commit, and suffix validation over
`AcceptedCandidateTransaction` through explicit exact-type branches. Preserve
invariant ordering and the existing named local-rejection set; zero-guide
proof failures abort.

- [ ] **Step 4: Run GREEN and structural counts**

```bash
pixi run pytest tests/adaptive/test_generator.py \
  -n auto --testmon -q
```

Record route number, family size, dispatch count, winner rank, operation count,
and elapsed time for route 1. This is structural measurement, not a
Held–Pfeiffer runtime comparison.

- [ ] **Step 5: Document, commit, and push**

Update Task 13F evidence in `docs/segment_site_mat.md` and check off the
corresponding causal traversal plan step:

```bash
git add src/compas_cgal/adaptive/generator.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_generator.py \
  docs/segment_site_mat.md \
  docs/superpowers/plans/2026-07-28-causal-mat-traversal.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): advance zero-guide route"
```

---

### Task 7: Fresh Replay of Advancing Segments

**Files**

- Modify: `src/compas_cgal/adaptive/replay.py`
- Modify: `src/compas_cgal/adaptive/replay_trace.py`
- Modify: `tests/adaptive/test_replay.py`
- Modify: `tests/adaptive/test_generator.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- Consumes:
  `AdvanceSegmentOperation`,
  `ZeroGuideLinkCandidate`,
  `ZeroGuideLinkTransaction`.
- Produces:
  unique operation-to-candidate reconstruction and fresh physical replay for
  both lateral variants.

- [ ] **Step 1: Write RED reconstruction tests**

From the Task 13F route-1 child, intercept the existing fresh-replay terminal
gate without weakening it:

```python
selected: list[ZeroGuideLinkCandidate] = []
traces: list[FreshReplayTrace] = []
real_match = replay_module._match_zero_guide_candidate
real_terminal = replay_module._require_terminal_replay

def tracked_match(**kwargs: object) -> ZeroGuideLinkCandidate:
    candidate = real_match(**kwargs)  # type: ignore[arg-type]
    selected.append(candidate)
    return candidate

def tracked_terminal(**kwargs: object) -> None:
    trace = kwargs["trace"]
    assert type(trace) is FreshReplayTrace
    traces.append(trace)
    real_terminal(**kwargs)  # type: ignore[arg-type]

monkeypatch.setattr(replay_module, "_match_zero_guide_candidate", tracked_match)
monkeypatch.setattr(replay_module, "_require_terminal_replay", tracked_terminal)
with pytest.raises(
    ReplayTraversalError,
    match="fresh MAT traversal remains nonterminal",
):
    _replay(task13f.identity, route_one_child.operations)

assert selected == [route_one_transaction.candidate]
assert traces[0].lateral_witnesses[-1].canonical_bytes == (
    route_one_transaction.segment_witness.canonical_bytes
)
```

The matcher must require the rebuilt MAT owner's zero-guide record bytes to
match the committed candidate before `_replay_fresh_state` begins physical
replay.

- [ ] **Step 2: Write RED grammar and mutation tests**

Relabel a hold link as `AdvanceSegmentOperation`, relabel an advancing segment
as a hold link, change its edge/run certificate, alter the endpoint, delete the
native record, and supply two matching candidates. Require
`ReplayZeroGuideCandidateError`; no mutation may fall through to circle
reconstruction.

- [ ] **Step 3: Run RED**

```bash
pixi run pytest tests/adaptive/test_replay.py \
  tests/adaptive/test_generator.py \
  -n auto --testmon -q
```

- [ ] **Step 4: Extend operation pairing**

Extend the existing `dict[int, TraversalCandidate]` owner mapping. Pair
`LinkSegmentOperation` only with its immediate circle and map both indices to
that circle candidate. Map each `AdvanceSegmentOperation` index to its one
zero-guide candidate. Reject every other grammar before candidate search; no
new wrapper type is needed for this simple mapping.

- [ ] **Step 5: Reconstruct uniquely from the fresh MAT**

For an advancing segment, enumerate only the freshly proved zero-guide family
on the authenticated directed window and require exactly one candidate whose
target, scope, cap, traversal decision, operation bytes, and run certificate
match. Extend `_replay_candidate_stream`, `_candidate_effective_cap`,
`_pair_lateral_operations`, and `_replay_fresh_state` over the closed
candidate union. Then run the same segment containment, TEA, depletion, and
coverage sequence as live evaluation using `OperationType.LINK`.

- [ ] **Step 6: Run GREEN, lint, and types**

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_replay.py \
  tests/adaptive/test_generator.py \
  -n auto --testmon -q
pixi run types-adaptive
```

- [ ] **Step 7: Document, commit, and push**

Update replay strength and remaining gaps in both developer pages:

```bash
git add src/compas_cgal/adaptive/replay.py \
  src/compas_cgal/adaptive/replay_trace.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_replay.py tests/adaptive/test_generator.py \
  tests/adaptive/typecheck/consumer_contract.py \
  docs/continuous_engagement.md docs/segment_site_mat.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "feat(adaptive): replay zero-guide links"
```

---

### Task 8: Continue to the Next Exact Boundary and Publish Evidence

**Files**

- Modify: `tests/adaptive/test_generator.py`
- Modify: `tests/adaptive/test_replay.py`
- Modify: `docs/segment_site_mat.md`
- Modify: `docs/continuous_engagement.md`
- Modify:
  `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`
- Modify:
  `docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md`

**Interfaces**

- Consumes:
  `generate_exact_adaptive_continuation(...)`,
  fresh replay terminal gate.
- Produces:
  either a complete Task 13F terminal certificate or a named, reproducible next
  exact boundary with all prior routes committed.

- [ ] **Step 1: Run the full Task 13F continuation**

```bash
pixi run pytest \
  tests/adaptive/test_generator.py::test_task13f_full_continuation \
  -n auto --testmon -q
```

If it fails, preserve the exact route/cursor/attempt/cap/exception evidence in
the test and documentation. Do not widen scope by guessing at the next
geometric theorem. Add this named test in Step 1; it does not pre-exist.

- [ ] **Step 2: Close terminal traversal and physical coverage when reached**

Require:

```python
continuation = generate_exact_adaptive_continuation(
    initial_evaluator=task13f.initial_evaluator,
    evaluator=task13f.evaluator,
    seeded_traversal=task13f.seeded_traversal,
    launch_transaction=task13f.launch_transaction,
)
continuation.traversal.require_terminal()
continuation.physical.coverage.require_complete()
certificate = replay_certificate(
    input_identity=task13f.identity,
    pocket=_task13f_pocket(),
    holes=(),
    cut_plane=task13f.identity.cut_plane,
    tool_radius=task13f.identity.tool_radius,
    user_cap=task13f.identity.user_cap,
    entry=task13f.identity.entry,
    operations=continuation.physical.operations,
    candidate_policy=task13f.identity.candidate_policy,
    neck_policy=task13f.identity.neck_policy,
    depletion_policy=task13f.identity.depletion_policy,
    traversal_policy=task13f.identity.traversal_policy,
    cut_direction_policy=task13f.identity.cut_direction_policy,
)
assert certificate.input_digest == task13f.identity.digest
```

Also require exact empty reachable residual under the existing terminal gate.
The current replay function intentionally has no success return after its
terminal gate. If generation reaches that boundary, preserve the exact failure
as the next stage and create a focused certificate-seal design; do not add a
return value inside this zero-guide plan.

- [ ] **Step 3: Measure the complete bounded path**

Only after the full path exists, record:

```text
native MAT builds
candidate families
candidate dispatches
exact segment certifications
exact circle certifications
fresh-replay certifications
wall-clock seconds
```

Compare against Held–Pfeiffer only on a matched disclosed workload. Until
then, retain “stronger proof contract; weaker measured-performance evidence.”

- [ ] **Step 4: Run the publication gate**

```bash
pixi run pytest tests/adaptive -n auto --testmon -q
pixi run pytest tests/adaptive -n auto -q
pixi run lint
pixi run types-adaptive
pixi run docs
git diff --check
```

- [ ] **Step 5: Update MkDocs and the insight record**

Update:

- implemented theorem and operation grammar;
- exact mutation evidence;
- Task 13F route and terminal status;
- structural counts and bounded runtime;
- current limitations;
- next exact boundary;
- Held–Pfeiffer `stronger | equivalent | weaker | incomplete` table.

- [ ] **Step 6: Commit and push the coherent checkpoint**

```bash
git add tests/adaptive/test_generator.py \
  tests/adaptive/test_replay.py \
  docs/segment_site_mat.md docs/continuous_engagement.md \
  docs/superpowers/plans/2026-07-28-causal-mat-traversal.md \
  docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md
git commit -m "docs(adaptive): record zero-guide stage"
```

Before pushing, list and cancel every `in_progress` run for
`codex/exact-certified-adaptive-phase1-t9` in `jf---/compas_cgal`. Push
without force, verify `git ls-remote` equals `git rev-parse HEAD`, and require
a clean worktree.
