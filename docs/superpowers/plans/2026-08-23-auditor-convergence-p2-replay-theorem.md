# Auditor Convergence P2 Replay and Theorem Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Return a mutation-resistant fresh replay certificate for a complete
toolpath and make every fast swept-prefix certification depend on a reviewed,
geometrically exercised theorem contract.

**Architecture:** Complete terminal generation without weakening traversal or
coverage, move replay-certificate ownership into a focused module, and bind the
existing `FreshReplayTrace` into a content-addressed success result. Treat the
swept-prefix fast path as research-only until its quantified theorem and exact
preconditions pass direct native and independent falsification gates.

**Tech Stack:** Python 3.12, C++20, CGAL exact kernels, CORE rationals, nanobind,
CCAN, SHA-256, Pixi, pytest-xdist, pytest-testmon, strict mypy, Ruff, MkDocs.

**Spec:** `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md`

## Global Constraints

- P1 acceptance commit is mandatory.
- No reference test is edited to accommodate changed production behaviour.
  A stale contract discovered by P0/P2 blocks at its evidence gate and is
  presented for an explicit user decision.
- Traversal terminality and exact residual emptiness are independent mandatory
  conditions.
- No sampled/raster check can certify; independent checks reject only.
- Empty event strata are unresolved.
- Generic event-exact certification is the production authority until the fast
  theorem gate is independently approved.
- Every pytest command uses `-n auto`; changed-code GREEN uses `--testmon`.
- No fallback, tolerance decision, swallowed exception, or silent retry.

---

### Task 1: Establish a complete terminal operation stream

**Files:**

- Create: `tests/adaptive/test_terminal_continuation.py`
- Modify: `src/compas_cgal/adaptive/generator.py`
- Modify: `src/compas_cgal/adaptive/generation_state.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `docs/segment_site_mat.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: `generate_exact_adaptive_continuation(...)`,
  `GenerationContinuation`, exact route-retrace commits.
- Produces: one terminal `GenerationContinuation` whose traversal has no active
  route and whose physical operation stream is complete enough for fresh replay.

- [ ] **Step 1: Reproduce the current Task 13F boundary**

```bash
pixi run pytest -- tests/adaptive/test_generator.py::test_task13f_full_continuation -n auto -q
```

Record the exact actual/expected behaviour in the progress ledger. Do not edit
the test. If it asserts an intentionally retired boundary, classify the
contract conflict with P0 evidence before proceeding.

- [ ] **Step 2: Write a new RED terminal-continuation acceptance test**

```python
def test_task13f_continuation_reaches_terminal_traversal(task13f: Task13FFixture) -> None:
    continuation = generate_exact_adaptive_continuation(
        initial_evaluator=task13f.initial_evaluator,
        evaluator=task13f.evaluator,
        seeded_traversal=task13f.seeded_traversal,
        launch_transaction=task13f.launch_transaction,
    )

    continuation.traversal.require_terminal()
    assert continuation.physical.operations[-1].canonical_bytes
    assert continuation.commits
```

Add non-vacuity assertions for route indices, circle/zero-guide/retrace commit
variants, and final physical phase.

- [ ] **Step 3: Run RED**

```bash
pixi run pytest -- tests/adaptive/test_terminal_continuation.py -n auto -q
```

Expected: the first unresolved exact event or terminal-seal defect is named.

- [ ] **Step 4: Repair one production boundary at a time**

For each failure, add one focused negative/positive test before changing code.
Narrow unions through exact commit/operation variants; do not use casts. Every
accepted operation must retain its motion, stock, coverage, traversal, and
canonical lineage.

- [ ] **Step 5: Run GREEN and regression gates**

```bash
pixi run ruff format src/compas_cgal/adaptive tests/adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/adaptive/test_terminal_continuation.py tests/adaptive/test_generator.py tests/adaptive/test_route_retrace_generator.py tests/adaptive/test_zero_guide_transaction.py -n auto --testmon -q
```

If the untouched reference test and the new acceptance test demand mutually
exclusive outcomes, stop with both command outputs and request the user's
contract decision.

- [ ] **Step 6: Document and commit**

```bash
git add src/compas_cgal/adaptive/generator.py src/compas_cgal/adaptive/generation_state.py src/compas_cgal/adaptive/errors.py tests/adaptive/test_terminal_continuation.py docs/segment_site_mat.md .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "feat(adaptive): reach terminal traversal"
```

### Task 2: Define and return the replay certificate

**Files:**

- Create: `src/compas_cgal/adaptive/replay_certificate.py`
- Modify: `src/compas_cgal/adaptive/replay.py`
- Modify: `src/compas_cgal/adaptive/replay_trace.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Create: `tests/adaptive/test_replay_certificate.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces:**

- Consumes: terminal `GenerationContinuation`, `FreshReplayTrace`,
  `BuildIdentity`, `InputIdentity`.
- Produces: `ReplayCertificate.build(...)`,
  `replay_certificate(...) -> ReplayCertificate`.

- [ ] **Step 1: Write RED success and mutation tests**

```python
def test_complete_program_returns_fresh_replay_certificate(task13f_terminal: TerminalFixture) -> None:
    certificate = replay_certificate(**task13f_terminal.replay_arguments)

    assert type(certificate) is ReplayCertificate
    assert certificate.digest == hashlib.sha256(certificate.canonical_bytes).digest()
    assert certificate.fresh_replay_trace_digest == task13f_terminal.expected_trace_digest


@pytest.mark.parametrize("field", ReplayCertificate.semantic_field_names())
def test_replay_certificate_rejects_each_mutated_semantic_field(
    terminal_certificate: ReplayCertificate,
    field: str,
) -> None:
    with pytest.raises(InvalidReplayCertificateError):
        _mutate_and_rebuild(terminal_certificate, field)
```

`semantic_field_names()` is a test helper over dataclass fields, not production
reflection. Mutations cover input/build/MAT/trace/stock/coverage/oracle/version
digests and operation-chain cardinality.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/adaptive/test_replay_certificate.py -n auto -q
```

Expected: missing module or unconditional final raise.

- [ ] **Step 3: Implement the focused certificate record**

```python
@dataclass(frozen=True)
class ReplayCertificate:
    input_digest: IdentityDigest
    build_identity_digest: IdentityDigest
    ordered_operation_digest: bytes
    operation_index_chain: tuple[IdentityDigest, ...]
    mat_certificate_digest: bytes
    fresh_replay_trace_digest: bytes
    terminal_stock_boundary_digest: bytes
    terminal_stock_lineage_digest: bytes
    terminal_coverage_certificate_digest: bytes
    component_versions: tuple[ComponentVersionBinding, ...]

    @classmethod
    def build(cls, *, input_identity: InputIdentity, build_identity: BuildIdentity,
              operation_index_chain: tuple[IdentityDigest, ...],
              mat_certificate_digest: bytes, trace: FreshReplayTrace,
              component_versions: tuple[ComponentVersionBinding, ...]) -> Self:
        operation_digest = _ordered_operation_digest(operation_index_chain)
        terminal = _validated_terminal_trace(trace)
        return cls(
            input_digest=input_identity.digest,
            build_identity_digest=build_identity.digest,
            ordered_operation_digest=operation_digest,
            operation_index_chain=operation_index_chain,
            mat_certificate_digest=_require_sha256(mat_certificate_digest),
            fresh_replay_trace_digest=trace.digest,
            terminal_stock_boundary_digest=terminal.stock_boundary_digest,
            terminal_stock_lineage_digest=terminal.stock_lineage_digest,
            terminal_coverage_certificate_digest=terminal.coverage_certificate_digest,
            component_versions=_validate_component_versions(component_versions),
        )
```

The factory derives ordered and terminal digests from owned preimages; callers
do not pass redundant aggregate bytes. Raw construction revalidates invariants.

- [ ] **Step 4: Retain `_replay_fresh_state` output and return**

Assign the returned trace, validate it once, and call
`ReplayCertificate.build(...)`. Remove the unconditional raise only after the
positive control reaches it. Do not add a partial certificate variant.

- [ ] **Step 5: Prove fresh-process identity**

Add a subprocess test that writes canonical input/operation bytes, reconstructs
all owners in a fresh Pixi process, and asserts identical certificate bytes and
digest.

- [ ] **Step 6: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal/adaptive tests/adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/adaptive/test_replay_certificate.py tests/adaptive/test_replay.py tests/adaptive/test_zero_guide_replay.py -n auto --testmon -q
git diff --check
git add src/compas_cgal/adaptive/replay_certificate.py src/compas_cgal/adaptive/replay.py src/compas_cgal/adaptive/replay_trace.py src/compas_cgal/adaptive/errors.py tests/adaptive/test_replay_certificate.py tests/adaptive/typecheck/consumer_contract.py
git commit -m "feat(adaptive): return replay certificate"
```

### Task 3: Close the empty-partition contradiction

**Files:**

- Modify: `src/continuous_tea_2/segment_oracle.cpp`
- Modify: `src/continuous_tea_2/event_trace.cpp`
- Create: `tests/adaptive/test_segment_oracle_empty_partition.py`
- Create: `tests/native/test_segment_oracle.cpp`

**Interfaces:**

- Consumes: `SegmentEventPartition2`, `ContinuousTeaVerdict`.
- Produces: invariant `CERTIFIED => nonempty strata && whole rim resolved`.

- [ ] **Step 1: Write RED empty-strata and contradictory-trace tests**

Use the native test seam to construct an empty verified partition and a trace
with certified/unresolved fields. Require `UNRESOLVED_DEGENERACY` for the first
and named rejection for the second.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/adaptive/test_segment_oracle.py -n auto -q
```

- [ ] **Step 3: Implement the fail-closed guards**

Include `partition.strata.empty()` in the unresolved disjunction. Tighten event
trace self-consistency so certified cannot coexist with unresolved whole-rim
disposition. Do not infer reachability; enforce the invariant structurally.

- [ ] **Step 4: Run GREEN and commit**

```bash
pixi run pytest -- tests/adaptive/test_segment_oracle.py tests/adaptive/test_event_substrate.py -n auto --testmon -q
git diff --check
git add src/continuous_tea_2/segment_oracle.cpp src/continuous_tea_2/event_trace.cpp tests/adaptive/test_segment_oracle.py tests/native/test_segment_oracle.cpp
git commit -m "fix(tea): reject empty segment strata"
```

### Task 4: Prove or disable the swept-prefix fast path

**Files:**

- Create: `docs/swept_prefix_theorem.md`
- Modify: `src/continuous_tea_2/segment_oracle.cpp`
- Modify: `src/continuous_tea_2/segment_oracle.h`
- Modify: `src/compas_cgal/adaptive/motion_certificate.py`
- Create: `tests/native/test_swept_prefix_theorem.cpp`
- Create: `tests/adaptive/test_swept_prefix_theorem.py`
- Modify: `mkdocs.yml`

**Interfaces:**

- Consumes: exact stock region, exact segment motion, exact pi-cap surrogate,
  source depletion lineage.
- Produces: reviewed `SweptPrefixSegmentTeaAudit2` or an explicit research-only
  unresolved result; never unconditional production certification.

- [ ] **Step 1: Write the quantified theorem document before code**

State domains, quantifiers, assumptions, and conclusion mathematically. Bind
each assumption to a native field. Explicitly distinguish a previously swept
capsule from a merely clear start disk.

- [ ] **Step 2: Write RED geometric tests**

Cover translation, Pythagorean rotation, binary64 scale, multiple stock
components, thin ribs, grazing contact, reversed motion, non-pi cap, uncleared
start, and a mutated source-lineage digest. Each positive case has exact contact
liveness; each negative case identifies the failed theorem precondition.

- [ ] **Step 3: Run RED against the two-boolean implementation**

```bash
pixi run pytest -- tests/adaptive/test_swept_prefix_theorem.py -n auto -q
```

Expected: at least the source-sweep/precondition mutations are accepted or not
represented.

- [ ] **Step 4: Encode every proved precondition or disable authority**

If the quantified proof survives review, bind exact prior-sweep containment and
lineage evidence into canonical bytes and native self-verification. Otherwise,
make the fast function return unresolved and route production retrace through
the generic event-exact segment certifier. The decision is one reviewed commit;
there is no runtime fallback between paths.

- [ ] **Step 5: Run theorem and retrace gates**

```bash
pixi run pytest -- tests/adaptive/test_swept_prefix_theorem.py tests/adaptive/test_route_retrace_transaction.py tests/adaptive/test_replay_certificate.py -n auto --testmon -q
pixi run -e docs docs
git diff --check
```

- [ ] **Step 6: Commit**

```bash
git add docs/swept_prefix_theorem.md mkdocs.yml src/continuous_tea_2/segment_oracle.cpp src/continuous_tea_2/segment_oracle.h src/compas_cgal/adaptive/motion_certificate.py tests/native/test_swept_prefix_theorem.cpp tests/adaptive/test_swept_prefix_theorem.py
git commit -m "fix(tea): bind swept-prefix theorem"
```

### Task 5: Add independent acceptance-side falsifiers and gate P2

**Files:**

- Create: `src/compas_cgal/adaptive/motion_acceptance_falsifier.py`
- Create: `tests/adaptive/test_motion_acceptance_falsifier.py`
- Modify: `docs/auditor_convergence.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: every certified segment/circle witness.
- Produces: rejection-only exact-cap cross-check; no certificate type.

- [ ] **Step 1: Write RED mutation and non-vacuity tests**

Inject an always-certified segment and circle verdict. Require the disjoint
exact `cap_exceeded` predicate to kill both. Require positive controls with
real stock contact to survive.

- [ ] **Step 2: Implement rejection-only dispatch**

The module receives a sealed certificate plus stock/motion preimages, evaluates
the independent exact cap flag at its mandated dyadic/event witnesses, and
raises `FalseMotionCertificateError` on counterevidence. It exposes no success
object and cannot be called as an acceptance authority.

- [ ] **Step 3: Run complete P2 gates**

```bash
pixi run ruff format src/compas_cgal/adaptive tests/adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/adaptive/test_replay_certificate.py tests/adaptive/test_motion_acceptance_falsifier.py tests/adaptive/test_segment_oracle.py tests/adaptive/test_circle_oracle.py tests/adaptive/test_swept_prefix_theorem.py -n auto -q
pixi run -e docs docs
git diff --check
```

- [ ] **Step 4: Document, update progress, and commit**

```bash
git add src/compas_cgal/adaptive/motion_acceptance_falsifier.py tests/adaptive/test_motion_acceptance_falsifier.py docs/auditor_convergence.md .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "test(tea): falsify accepted motions"
```

Record replay certificate digest, fresh-process equality, theorem disposition,
and exact P3/P4 start SHAs.
