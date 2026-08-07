# Exact Inter-Route Retrace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cross the Task 13F route-1-to-route-2 discontinuity with one exact, causally authenticated cut-depth retrace and reproduce it under fresh replay.

**Architecture:** The immediately preceding accepted no-neck `AdvanceSegmentOperation` is the sole retrace source. A distinct decision, operation, physical transaction, and cross-axis commit reverse that exact source while preserving MAT cursor and passage state; `GenerationContinuation` keeps traversal and retrace commits in one ordered closed union. Live evaluation and fresh replay share one swept-prefix physical proof implementation, while route activation remains independently authenticated global state.

**Tech Stack:** Python 3.12, COMPAS frame-tagged geometry, CGAL exact native stock/coverage owners, frozen typed dataclasses, CCAN, SHA-256, pytest-xdist, pytest-testmon, Ruff, strict mypy, strict MkDocs.

## Global Constraints

- Governing design: `docs/superpowers/specs/2026-08-06-exact-inter-route-retrace-design.md`, approved 2026-08-07.
- Governing master plan: `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`, bounded Task 13F continuation.
- Work only in `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-exact-certified-adaptive-phase1-t9-zero-guide` on `codex/exact-certified-adaptive-phase1-t9-zero-guide`; never mutate `main`, `master`, or the protected checkout.
- Exact branch-switch trigger: `completed.exit_node_id != activated.entry_node_id`. No coordinate comparison, distance, tolerance, epsilon, sampling, or reporting surrogate may decide it.
- The admitted source is exactly the final `AdvanceSegmentOperation` witnessed by the immediately preceding `ZeroGuideLinkTransaction` and `TraversalCommit`; it must be no-neck, full-cap, route-terminal, and physically current.
- Retrace motion is derived, never submitted: `start = source.motion.end` and `end = source.motion.start` at the identical `CutZ`.
- Retrace changes phase, stock lineage, coverage lineage, and operations only. It advances no local MAT cursor, global route cursor, neck passage, or candidate index.
- Evaluation order is containment, pre-motion stock certifier, one `certify_swept_prefix_segment` call, depletion, then coverage. No generic event certification or copied deciding proof path.
- `GenerationContinuation.commits` is one ordered `TraversalCommit | RouteRetraceCommit` stream. No parallel retrace ledger or compatibility field.
- Incident route activation remains zero-motion. Nonincident activation is unavailable to candidate materialization until the retrace physical child is committed.
- Fresh replay rebuilds the source zero-guide candidate first, skips retrace during candidate enumeration, preserves candidate indices, and reruns the shared physical proof.
- Exact failures stay loud through the seven named retrace errors. No fallback, conditional import, `HAS_*` branch, skip, skipif, xfail, or weakened reference assertion.
- Preserve the route-2 direct-path negative control: cursor `def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0`, 56 attempts, 56 gouges, and family digest `56d92fcf2089c1e771a9add79bd356d172904a69e6c7acb18ab02522f4ed093b`.
- Preserve existing byte identities through characterization: advancing witness `0d7f0ac6c7b58ebe114b5c5453ce12b014a990687bda16bfcffa199874762cc7`, zero-guide transaction `5668c64ae62f1aa1b52f7338db626626ec0636eef66c731934917b9b331a40b4`, and traversal commit `0a20db00bf444336ffd0138da9e4a938f28dd8bf7bc0baacb69ebdc8f1656c73`.
- Every new test has a contextual Google-style docstring explaining the geometry, causal authority, and failure boundary.
- Run every pytest command with `-n auto`. After code changes, use `--testmon` before the broader affected and full suites.
- Keep one responsibility per new file. Do not introduce result wrappers, enums, factories, or class hierarchies beyond the domain records and proof authorities required by the approved design.
- Use Pixi exclusively. Format Python with Ruff before each code commit and gate type claims with `mypy --strict`.
- Commit author and committer must both be `Jelle Feringa <jelleferinga@gmail.com>`. Never add attribution trailers, force-push, or leave in-flight branch CI running before a push.
- Documentation and insight capture are part of the coherent implementation checkpoint, not retrospective cleanup.

---

### Task 1: Seal the causal decision and canonical operation grammar

**Files:**

- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `src/compas_cgal/adaptive/operation.py`
- Modify: `tests/adaptive/task13f_fixture.py`
- Create: `tests/adaptive/test_route_retrace_operation.py`

**Interfaces:**

- Consumes: `AdvanceSegmentOperation`, `NoNeckScope`, `FullCapDecision`, `ExactSegmentMotion.build(start, end)`, `IdentityDigest`, and existing CCAN encoders.
- Produces: `RouteNodeId`, `RouteRetraceDecision.build(...) -> RouteRetraceDecision`, `RetraceSegmentOperation.build(*, source_operation, decision) -> RetraceSegmentOperation`, and the expanded `CanonicalOperation` union.

- [ ] **Step 1: Add the three operation-boundary errors**

Add these named public errors beside the existing operation and traversal errors:

~~~python
class InvalidRouteRetraceDecisionError(ArtifactIdentityError):
    """Raised when causal route-boundary identities are malformed or cross-wired."""


class InvalidRetraceSegmentOperationError(ArtifactIdentityError):
    """Raised when a retrace is not the exact admitted source reversal."""


class UnsupportedRouteRetraceError(RuntimeError):
    """Raised when a nonincident route boundary has no admitted retrace source."""
~~~

- [ ] **Step 2: Add one shared Task 13F established-prefix helper**

Add this tuple-returning test helper to `task13f_fixture.py` so every retrace
test uses the same accepted production prefix without inventing a fixture
result type:

~~~python
def task13f_route_one_terminal(
    fixture: Task13FFixture,
) -> tuple[
    GenerationState,
    MatTraversalState,
    tuple[TraversalCommit, TraversalCommit],
]:
    physical, traversal, route_zero = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=fixture.physical,
        traversal=fixture.traversal,
    )
    assert traversal.active_cursor.terminal
    physical, traversal, route_one = advance_active_candidate_family(
        evaluator=fixture.evaluator,
        physical=physical,
        traversal=traversal.activate_next(),
    )
    assert traversal.active_route_index == 1
    assert traversal.active_cursor.terminal
    return physical, traversal, (route_zero, route_one)


def task13f_retrace_decision(
    *,
    physical: GenerationState,
    terminal: MatTraversalState,
    source_commit: TraversalCommit,
) -> RouteRetraceDecision:
    """Seal fixture preimages for physical-unit tests before generator wiring."""
    activated = terminal.activate_next()
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    return RouteRetraceDecision.build(
        completed_route_index=terminal.active_route_index,
        activated_route_index=activated.active_route_index,
        completed_exit_node_id=RouteNodeId(
            bytes(terminal.active_cursor.route_step.exit_node_id),
        ),
        activated_entry_node_id=RouteNodeId(
            bytes(activated.active_cursor.route_step.entry_node_id),
        ),
        terminal_traversal_digest=terminal.digest,
        activated_traversal_digest=activated.digest,
        source_commit_digest=source_commit.digest,
        source_transaction_digest=source_commit.transaction.digest,
        source_operation_index=len(physical.operations) - 1,
        source_operation_digest=IdentityDigest(
            hashlib.sha256(source.canonical_bytes).digest(),
        ),
    )
~~~

- [ ] **Step 3: Write RED tests for local decision invariants**

Use a small `decision_fields` fixture with distinct 32-byte digests and node identities. Cover stable canonical bytes plus every locally decidable rejection:

~~~python
@pytest.fixture
def decision_fields() -> dict[str, object]:
    return {
        "completed_route_index": 1,
        "activated_route_index": 2,
        "completed_exit_node_id": RouteNodeId(b"completed-exit"),
        "activated_entry_node_id": RouteNodeId(b"activated-entry"),
        "terminal_traversal_digest": IdentityDigest(
            hashlib.sha256(b"terminal").digest(),
        ),
        "activated_traversal_digest": IdentityDigest(
            hashlib.sha256(b"activated").digest(),
        ),
        "source_commit_digest": IdentityDigest(
            hashlib.sha256(b"source-commit").digest(),
        ),
        "source_transaction_digest": IdentityDigest(
            hashlib.sha256(b"source-transaction").digest(),
        ),
        "source_operation_index": 5,
        "source_operation_digest": IdentityDigest(
            hashlib.sha256(b"source-operation").digest(),
        ),
    }


def test_route_retrace_decision_is_content_addressed(
    decision_fields: dict[str, object],
) -> None:
    """Bind a nonincident route switch and its accepted source by exact identity."""
    decision = RouteRetraceDecision.build(**decision_fields)
    rebuilt = RouteRetraceDecision.build(
        **dict(reversed(tuple(decision_fields.items()))),
    )

    assert decision.activated_route_index == decision.completed_route_index + 1
    assert decision.completed_exit_node_id != decision.activated_entry_node_id
    assert rebuilt.canonical_bytes == decision.canonical_bytes
    assert decision.digest == IdentityDigest(
        hashlib.sha256(decision.canonical_bytes).digest(),
    )


@pytest.mark.parametrize(
    ("field", "replacement"),
    (
        ("completed_route_index", -1),
        ("activated_route_index", 4),
        ("completed_exit_node_id", RouteNodeId(b"same-node")),
        ("terminal_traversal_digest", IdentityDigest(b"short")),
        ("source_operation_index", 1),
    ),
)
def test_route_retrace_decision_rejects_malformed_identity(
    decision_fields: dict[str, object],
    field: str,
    replacement: object,
) -> None:
    """Reject a decision that cannot name one adjacent authenticated boundary."""
    changed = dict(decision_fields)
    changed[field] = replacement
    if field == "completed_exit_node_id":
        changed["activated_entry_node_id"] = replacement

    with pytest.raises(InvalidRouteRetraceDecisionError):
        RouteRetraceDecision.build(**changed)
~~~

- [ ] **Step 4: Run the decision tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_operation.py -n auto -q
~~~

Expected: collection fails because `RouteNodeId` and `RouteRetraceDecision` do not exist.

- [ ] **Step 5: Add `RouteRetraceDecision` with exact typed fields**

Define `RouteNodeId = NewType("RouteNodeId", bytes)` and these fixed fields:

~~~python
ROUTE_RETRACE_STRATEGY_VERSION = b"reverse-final-zero-guide-v1"


@dataclass(frozen=True)
class RouteRetraceDecision:
    completed_route_index: int
    activated_route_index: int
    completed_exit_node_id: RouteNodeId
    activated_entry_node_id: RouteNodeId
    terminal_traversal_digest: IdentityDigest
    activated_traversal_digest: IdentityDigest
    source_commit_digest: IdentityDigest
    source_transaction_digest: IdentityDigest
    source_operation_index: int
    source_operation_digest: IdentityDigest
    strategy_identity: bytes = ROUTE_RETRACE_STRATEGY_VERSION

    @classmethod
    def build(
        cls,
        *,
        completed_route_index: int,
        activated_route_index: int,
        completed_exit_node_id: RouteNodeId,
        activated_entry_node_id: RouteNodeId,
        terminal_traversal_digest: IdentityDigest,
        activated_traversal_digest: IdentityDigest,
        source_commit_digest: IdentityDigest,
        source_transaction_digest: IdentityDigest,
        source_operation_index: int,
        source_operation_digest: IdentityDigest,
    ) -> Self:
        return cls(
            completed_route_index,
            activated_route_index,
            completed_exit_node_id,
            activated_entry_node_id,
            terminal_traversal_digest,
            activated_traversal_digest,
            source_commit_digest,
            source_transaction_digest,
            source_operation_index,
            source_operation_digest,
        )
~~~

In `__post_init__` require exact nonnegative integers, adjacent route indices,
distinct nonempty node bytes, five 32-byte SHA-256 identities, a source
operation index of at least 2, and the fixed strategy identity. Encode every
field under `route-retrace-decision-v1`; expose `digest` as SHA-256 of
`canonical_bytes`.

- [ ] **Step 6: Write RED operation derivation and tamper tests**

Build the source from the real Task 13F route-1 advancing operation. Assert
exact reverse derivation, same depth/cap, no traversal or neck fields, and
stable bytes. Parameterize locally decidable source replacements covering
neck-owning, non-full-cap, nonterminal, wrong source digest, circle, and
hold-link variants. Constructing forward, shortened, lengthened, offset, or
zero-length submitted geometry must be impossible because the factory accepts
no motion argument:

~~~python
def test_retrace_operation_derives_only_the_exact_source_reverse(
    task13f: Task13FFixture,
) -> None:
    """Reverse the accepted horizontal leaf without accepting caller geometry."""
    physical, terminal, commits = task13f_route_one_terminal(task13f)
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=commits[-1],
    )
    operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=decision,
    )

    assert operation.motion.start == source.motion.end
    assert operation.motion.end == source.motion.start
    assert operation.cut_z == source.cut_z
    assert operation.effective_cap_decision == source.effective_cap_decision
    assert not hasattr(operation, "traversal_decision")
    assert not hasattr(operation, "neck_scope")
~~~

The state/transaction tasks test source ordinal, physical finality, and
foreign cut-depth authority because those checks require the owned parent
preimage and cannot soundly be guessed by an operation-local factory.

- [ ] **Step 7: Run the operation tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_operation.py -n auto -q
~~~

Expected: decision tests pass and operation tests fail because `RetraceSegmentOperation` is absent.

- [ ] **Step 8: Implement the derived retrace operation**

Add:

~~~python
@dataclass(frozen=True, init=False)
class RetraceSegmentOperation:
    motion: ExactSegmentMotion
    cut_z: CutZ
    effective_cap_decision: FullCapDecision
    decision: RouteRetraceDecision

    @classmethod
    def build(
        cls,
        *,
        source_operation: AdvanceSegmentOperation,
        decision: RouteRetraceDecision,
    ) -> Self:
        if type(source_operation) is not AdvanceSegmentOperation:
            raise InvalidRetraceSegmentOperationError(
                "retrace source must be one exact advancing segment.",
            )
        if (
            type(source_operation.neck_scope) is not NoNeckScope
            or type(source_operation.effective_cap_decision) is not FullCapDecision
            or not source_operation.traversal_decision.makes_cursor_terminal
        ):
            raise InvalidRetraceSegmentOperationError(
                "retrace source must be no-neck, full-cap, and route-terminal.",
            )
        source_digest = IdentityDigest(
            hashlib.sha256(source_operation.canonical_bytes).digest(),
        )
        if decision.source_operation_digest != source_digest:
            raise InvalidRetraceSegmentOperationError(
                "retrace decision does not name its source operation.",
            )
        motion = ExactSegmentMotion.build(
            source_operation.motion.end,
            source_operation.motion.start,
        )
        operation = object.__new__(cls)
        object.__setattr__(operation, "motion", motion)
        object.__setattr__(operation, "cut_z", source_operation.cut_z)
        object.__setattr__(
            operation,
            "effective_cap_decision",
            source_operation.effective_cap_decision,
        )
        object.__setattr__(operation, "decision", decision)
        operation._validate()
        return operation
~~~

Make `_validate` enforce the closed types and fixed decision strategy.
`init=False` makes the source-derived factory the only normal construction
path; state and replay still independently revalidate the owned source
preimage. Encode `cut-z`, `decision`, `effective-cap-decision`, and `motion`
under the distinct `retrace-segment-v1` tag. Add the type to
`CanonicalOperation`.

- [ ] **Step 9: Run focused GREEN, format, lint, and strict types**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_route_retrace_operation.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
~~~

Expected: all commands pass.

- [ ] **Step 10: Commit and publish Task 1**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/errors.py src/compas_cgal/adaptive/operation.py tests/adaptive/task13f_fixture.py tests/adaptive/test_route_retrace_operation.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "feat(adaptive): add retrace operation"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 2: Factor one swept-prefix physical proof path

**Files:**

- Create: `src/compas_cgal/adaptive/retrace_segment_trial.py`
- Modify: `src/compas_cgal/adaptive/advancing_segment_trial.py`
- Modify: `src/compas_cgal/adaptive/replay_trace.py`
- Create: `tests/adaptive/test_retrace_segment_trial.py`

**Interfaces:**

- Consumes: `AdvanceSegmentOperation | RetraceSegmentOperation`, `GougeContainment`, `Stock2Area`, `CoverageLedger`, and `MotionCertifier.certify_swept_prefix_segment(...)`.
- Produces: `evaluate_retrace_segment_trial(...) -> ReplayLateralWitness` and private `_evaluate_swept_prefix_segment_trial(...) -> ReplayLateralWitness`. The existing `evaluate_advancing_segment_trial(...)` signature remains unchanged and delegates to the shared implementation.

- [ ] **Step 1: Characterize the established advancing witness before refactoring**

Use the real Task 13F route-1 accepted transaction:

~~~python
def test_shared_trial_preserves_established_advancing_identity(
    task13f: Task13FFixture,
) -> None:
    """Keep the accepted route-1 proof byte-identical while sharing its theorem."""
    physical, traversal, _ = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )
    traversal = traversal.activate_next()
    _, _, route_one_commit = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=physical,
        traversal=traversal,
    )

    witness = route_one_commit.transaction.segment_witness
    assert witness.digest.hex() == (
        "0d7f0ac6c7b58ebe114b5c5453ce12b014a990687bda16bfcffa199874762cc7"
    )
    assert route_one_commit.transaction.digest.hex() == (
        "5668c64ae62f1aa1b52f7338db626626ec0636eef66c731934917b9b331a40b4"
    )
    assert route_one_commit.digest.hex() == (
        "0a20db00bf444336ffd0138da9e4a938f28dd8bf7bc0baacb69ebdc8f1656c73"
    )
~~~

Run:

~~~bash
pixi run pytest tests/adaptive/test_retrace_segment_trial.py::test_shared_trial_preserves_established_advancing_identity -n auto -q
~~~

Expected: PASS on the unrefactored production path.

- [ ] **Step 2: Write RED retrace order and theorem-count tests**

Build the exact reverse arguments directly from the shared accepted prefix:

~~~python
def _retrace_trial_arguments(
    task13f: Task13FFixture,
) -> dict[str, object]:
    physical, terminal, commits = task13f_route_one_terminal(task13f)
    source_commit = commits[-1]
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=source_commit,
    )
    return {
        "containment_authority": GougeContainment.build(
            task13f.identity.reachable_domain,
        ),
        "stock": physical.fork_stock(),
        "coverage": physical.clone_coverage(),
        "operation_index": len(physical.operations),
        "operation": RetraceSegmentOperation.build(
            source_operation=source,
            decision=decision,
        ),
        "tool_radius": task13f.identity.tool_radius,
        "user_cap": task13f.identity.user_cap,
        "effective_cap": task13f.identity.user_cap,
        "depletion_policy": task13f.identity.depletion_policy,
    }
~~~

Instrument the exact authorities with wrappers that append `"containment"`,
`"swept-prefix"`, `"deplete"`, and `"coverage"`. Patch
`MotionCertifier.certify` to raise immediately, proving the generic event
oracle is unreachable:

~~~python
def test_retrace_uses_one_ordered_two_stratum_swept_prefix_proof(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Certify the reverse corridor once before either trial owner mutates."""
    order: list[str] = []
    calls = {"generic": 0, "swept_prefix": 0}
    real_containment = GougeContainment.certify_segment
    real_swept_prefix = MotionCertifier.certify_swept_prefix_segment
    real_deplete = Stock2Area.deplete
    real_coverage = CoverageLedger.add_sweep

    def tracked_containment(self: GougeContainment, *args: object) -> object:
        order.append("containment")
        return real_containment(self, *args)  # type: ignore[arg-type]

    def tracked_swept_prefix(
        self: MotionCertifier,
        **kwargs: object,
    ) -> SweptPrefixMotionWitness:
        order.append("swept-prefix")
        calls["swept_prefix"] += 1
        return real_swept_prefix(self, **kwargs)  # type: ignore[arg-type]

    def reject_generic(self: MotionCertifier, **kwargs: object) -> NoReturn:
        calls["generic"] += 1
        raise AssertionError("generic event certification is forbidden")

    def tracked_deplete(self: Stock2Area, *args: object) -> object:
        order.append("deplete")
        return real_deplete(self, *args)  # type: ignore[arg-type]

    def tracked_coverage(self: CoverageLedger, *args: object) -> object:
        order.append("coverage")
        return real_coverage(self, *args)  # type: ignore[arg-type]

    monkeypatch.setattr(GougeContainment, "certify_segment", tracked_containment)
    monkeypatch.setattr(
        MotionCertifier,
        "certify_swept_prefix_segment",
        tracked_swept_prefix,
    )
    monkeypatch.setattr(MotionCertifier, "certify", reject_generic)
    monkeypatch.setattr(Stock2Area, "deplete", tracked_deplete)
    monkeypatch.setattr(CoverageLedger, "add_sweep", tracked_coverage)

    witness = evaluate_retrace_segment_trial(
        **_retrace_trial_arguments(task13f),
    )  # type: ignore[arg-type]

    assert order == ["containment", "swept-prefix", "deplete", "coverage"]
    assert calls == {"generic": 0, "swept_prefix": 1}
    assert type(witness.motion_witness) is SweptPrefixMotionWitness
    assert witness.motion_witness.event_cell_count == SWEPT_PREFIX_MOTION_STRATA == 2
~~~

- [ ] **Step 3: Run retrace trial tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_retrace_segment_trial.py -n auto -q
~~~

Expected: characterization passes; retrace test fails because `evaluate_retrace_segment_trial` is absent.

- [ ] **Step 4: Create the single deciding implementation**

Add the current ordered algorithm to `retrace_segment_trial.py` under this
exact closed union:

~~~python
SweptPrefixSegmentOperation: TypeAlias = (
    AdvanceSegmentOperation | RetraceSegmentOperation
)


def _evaluate_swept_prefix_segment_trial(
    *,
    containment_authority: GougeContainment,
    stock: Stock2Area,
    coverage: CoverageLedger,
    operation_index: int,
    operation: SweptPrefixSegmentOperation,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    effective_cap: EngagementCap,
    depletion_policy: DepletionPolicy,
) -> ReplayLateralWitness:
    containment = containment_authority.certify_segment(
        operation.motion,
        tool_radius,
    )
    certifier = MotionCertifier.build(stock=stock, tool_radius=tool_radius)
    motion_witness = certifier.certify_swept_prefix_segment(
        operation_index=operation_index,
        motion=operation.motion,
        user_cap=user_cap,
        effective_cap=effective_cap,
    )
    depletion = stock.deplete(operation.motion, tool_radius, depletion_policy)
    sweep = coverage.add_sweep(operation.motion, tool_radius)
    return ReplayLateralWitness(
        operation_index=operation_index,
        operation=operation,
        effective_cap_decision=operation.effective_cap_decision,
        stock_boundary_digest=certifier.canonical_boundary_digest,
        containment_certificate=containment,
        motion_witness=motion_witness,
        depletion_witness=depletion,
        sweep_witness=sweep,
    )
~~~

`evaluate_retrace_segment_trial` delegates to this function with the same
keyword signature, except its operation type is
`RetraceSegmentOperation`. At this step, leave the established advancing body
unchanged.

- [ ] **Step 5: Extend `ReplayLateralWitness` without weakening ordinary motions**

Add `RetraceSegmentOperation` to the operation union. Require `SweptPrefixMotionWitness` with the fixed strategy, theorem, stock boundary, and exactly two strata for both advancing and retrace operations:

~~~python
if type(witness.operation) in (
    AdvanceSegmentOperation,
    RetraceSegmentOperation,
):
    if type(witness.motion_witness) is not SweptPrefixMotionWitness:
        raise InvalidReplayTraceError(
            "swept-prefix segment requires one exact swept-prefix witness.",
        )
    if (
        witness.motion_witness.strategy_identity
        != SWEPT_PREFIX_STRATEGY_VERSION
        or witness.motion_witness.theorem_identity
        != SWEPT_PREFIX_THEOREM_VERSION
        or witness.motion_witness.event_cell_count
        != SWEPT_PREFIX_MOTION_STRATA
        or witness.motion_witness.stock_boundary_digest
        != witness.stock_boundary_digest
    ):
        raise InvalidReplayTraceError(
            "swept-prefix segment contradicts its exact theorem.",
        )
~~~

Keep `MotionWitness` mandatory for links and circles.

- [ ] **Step 6: Validate the additive retrace path before switching advancing**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_retrace_segment_trial.py -n auto -q
~~~

Expected: the new retrace proof passes and the unchanged advancing
characterization still reports all three baseline digests.

- [ ] **Step 7: Switch advancing to the validated shared implementation**

Replace only the body of `evaluate_advancing_segment_trial` with a keyword
delegation:

~~~python
return _evaluate_swept_prefix_segment_trial(
    containment_authority=containment_authority,
    stock=stock,
    coverage=coverage,
    operation_index=operation_index,
    operation=operation,
    tool_radius=tool_radius,
    user_cap=user_cap,
    effective_cap=effective_cap,
    depletion_policy=depletion_policy,
)
~~~

The public signature and Google-style docstring remain unchanged. The
additive GREEN result from Step 6 is the authorization gate for this switch.

- [ ] **Step 8: Run identity regression and strict gates**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_retrace_segment_trial.py tests/adaptive/test_zero_guide_transaction.py tests/adaptive/test_zero_guide_replay.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
~~~

Expected: all tests pass and all three established digests remain exact.

- [ ] **Step 9: Commit and publish Task 2**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/retrace_segment_trial.py src/compas_cgal/adaptive/advancing_segment_trial.py src/compas_cgal/adaptive/replay_trace.py tests/adaptive/test_retrace_segment_trial.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "refactor(adaptive): share swept-prefix trial"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 3: Admit retrace into authoritative physical chronology

**Files:**

- Modify: `src/compas_cgal/adaptive/generation_state.py`
- Modify: `src/compas_cgal/adaptive/transaction.py`
- Create: `tests/adaptive/test_route_retrace_state.py`

**Interfaces:**

- Consumes: `RetraceSegmentOperation` and a trial-local stock/coverage pair already mutated by `evaluate_retrace_segment_trial`.
- Produces: `GenerationState.build(...)` validation for retrace chronology while retaining the prior `TraversalCursorState` and passage tuple byte for byte. `CandidateEvaluator._validate_physical_parent(...)` explicitly recognizes the new parent operation.

- [ ] **Step 1: Write RED state-mutation tests**

Build the real route-1 parent, evaluate its exact reverse on forks, and return
ordinary domain values:

~~~python
def _evaluated_retrace(
    task13f: Task13FFixture,
) -> tuple[
    GenerationState,
    AdvanceSegmentOperation,
    RetraceSegmentOperation,
    Stock2Area,
    CoverageLedger,
]:
    parent, terminal, commits = task13f_route_one_terminal(task13f)
    source = parent.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    operation = RetraceSegmentOperation.build(
        source_operation=source,
        decision=task13f_retrace_decision(
            physical=parent,
            terminal=terminal,
            source_commit=commits[-1],
        ),
    )
    stock = parent.fork_stock()
    coverage = parent.clone_coverage()
    evaluate_retrace_segment_trial(
        containment_authority=GougeContainment.build(
            task13f.identity.reachable_domain,
        ),
        stock=stock,
        coverage=coverage,
        operation_index=len(parent.operations),
        operation=operation,
        tool_radius=task13f.identity.tool_radius,
        user_cap=task13f.identity.user_cap,
        effective_cap=task13f.identity.user_cap,
        depletion_policy=task13f.identity.depletion_policy,
    )
    return parent, source, operation, stock, coverage


def test_retrace_state_changes_only_physical_chronology(
    task13f: Task13FFixture,
) -> None:
    """Return through removed stock without manufacturing MAT progress."""
    parent, source, operation, stock, coverage = _evaluated_retrace(task13f)
    child = GenerationState.build(
        stock=stock,
        coverage=coverage,
        tool_radius=parent.tool_radius,
        phase_point=operation.motion.end,
        traversal=parent.traversal,
        passages=parent.passages,
        operations=parent.operations + (operation,),
    )

    assert child.phase_point == source.motion.start
    assert child.traversal.canonical_bytes == parent.traversal.canonical_bytes
    assert child.passages == parent.passages
    assert len(child.fork_stock().lineage) == len(parent.fork_stock().lineage) + 1
    assert len(child.clone_coverage().lineage) == len(parent.clone_coverage().lineage) + 1
    assert child.operations[-1] == operation
~~~

Add focused corruptions for wrong predecessor type, source ordinal, source
digest, a deliberately forged non-reverse stored operation, depth, cap, start
phase, pending hold link, changed local cursor, changed passage, missing stock
lineage, and missing coverage lineage.

- [ ] **Step 2: Run state tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_state.py -n auto -q
~~~

Expected: child construction fails because retrace is outside the closed state grammar.

- [ ] **Step 3: Extend the closed operation and chronology checks**

Add `RetraceSegmentOperation` to both state operation unions. Enumerate lateral operations with their global ordinal:

~~~python
for operation_index, operation in enumerate(lateral, start=2):
    if type(operation) is RetraceSegmentOperation:
        if pending_link is not None:
            raise InvalidGenerationStateError(
                "route retrace cannot consume a pending hold link.",
            )
        source_index = operation.decision.source_operation_index
        if source_index != operation_index - 1:
            raise InvalidGenerationStateError(
                "route retrace must immediately follow its named source.",
            )
        source = self.operations[source_index]
        if type(source) is not AdvanceSegmentOperation:
            raise InvalidGenerationStateError(
                "route retrace source is not an advancing segment.",
            )
        source_digest = IdentityDigest(
            hashlib.sha256(source.canonical_bytes).digest(),
        )
        if (
            source_digest != operation.decision.source_operation_digest
            or operation.motion.start != current_phase
            or operation.motion.start != source.motion.end
            or operation.motion.end != source.motion.start
            or operation.cut_z != source.cut_z
            or operation.effective_cap_decision
            != source.effective_cap_decision
        ):
            raise InvalidGenerationStateError(
                "route retrace contradicts its source or physical phase.",
            )
        current_phase = operation.motion.end
        continue
~~~

Do not update `last_advance` in this branch. The existing final cursor check must therefore still authenticate the source advance. `_validate_passages` must ignore retrace explicitly because it owns no passage transition.

- [ ] **Step 4: Make candidate-parent validation variant-aware**

Branch before reading `neck_scope`:

~~~python
for raw_operation in state.operations[2:]:
    if type(raw_operation) is RetraceSegmentOperation:
        if raw_operation.effective_cap_decision != expected_full_cap:
            raise CandidateStateMismatchError(
                "candidate parent retrace contradicts evaluator full cap.",
            )
        continue
    if type(raw_operation) not in (
        LinkSegmentOperation,
        CutFullCircleOperation,
        AdvanceSegmentOperation,
    ):
        raise CandidateStateMismatchError(
            "candidate parent contains a foreign lateral operation.",
        )
    operation = cast(
        LinkSegmentOperation | CutFullCircleOperation | AdvanceSegmentOperation,
        raw_operation,
    )
~~~

This preserves one evaluator policy authority and avoids inventing a neck field on retrace.

- [ ] **Step 5: Run focused GREEN and parent-consumer regression**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_route_retrace_state.py tests/adaptive/test_generation_state.py tests/adaptive/test_transaction.py tests/adaptive/test_zero_guide_transaction.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
~~~

Expected: all commands pass.

- [ ] **Step 6: Commit and publish Task 3**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/generation_state.py src/compas_cgal/adaptive/transaction.py tests/adaptive/test_route_retrace_state.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "feat(adaptive): validate retrace state"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 4: Add atomic retrace evaluation and independent commit

**Files:**

- Create: `src/compas_cgal/adaptive/retrace_transaction.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Create: `tests/adaptive/test_route_retrace_transaction.py`

**Interfaces:**

- Consumes: `CandidateEvaluator` as the sole physical-policy authority, `GenerationState`, `RouteRetraceDecision`, and `evaluate_retrace_segment_trial(...)`.
- Produces: `RouteRetraceTransaction` and `RouteRetraceEvaluator.build(*, evaluator)` with `evaluate(state, decision) -> RouteRetraceTransaction` and `commit(state, transaction) -> GenerationState`.

- [ ] **Step 1: Add transaction errors**

~~~python
class InvalidRouteRetraceTransactionError(ValueError):
    """Raised when retrace physical evidence is malformed or cross-wired."""


class StaleRouteRetraceTransactionError(RuntimeError):
    """Raised when a retrace no longer names the authoritative parent."""
~~~

- [ ] **Step 2: Write RED transaction identity tests**

~~~python
def test_route_retrace_transaction_binds_one_physical_suffix(
    task13f: Task13FFixture,
) -> None:
    """Bind the exact reverse witness to one immutable parent and child."""
    parent, terminal, commits = task13f_route_one_terminal(task13f)
    evaluator = RouteRetraceEvaluator.build(evaluator=task13f.evaluator)
    decision = task13f_retrace_decision(
        physical=parent,
        terminal=terminal,
        source_commit=commits[-1],
    )
    transaction = evaluator.evaluate(parent, decision)
    child = evaluator.commit(parent, transaction)

    assert transaction.parent_state_digest == parent.digest
    assert transaction.decision == transaction.segment_witness.operation.decision
    assert transaction.result_state_digest == child.digest
    assert child.operations == (
        parent.operations + (transaction.segment_witness.operation,)
    )
~~~

Raw-construction tests must reject a foreign operation, wrong decision, wrong theorem/strategy/stratum count, wrong stock boundary, broken stock or coverage parent lineage, wrong operation index, and malformed parent/child digest.

- [ ] **Step 3: Write RED atomic evaluation tests**

Capture all parent identities before each trial and inject failures at containment, swept-prefix proof, depletion, coverage, and child construction:

~~~python
def test_failed_retrace_trial_leaves_parent_byte_identical(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep every authoritative parent owner immutable after proof failure."""
    parent, terminal, commits = task13f_route_one_terminal(task13f)
    evaluator = RouteRetraceEvaluator.build(evaluator=task13f.evaluator)
    decision = task13f_retrace_decision(
        physical=parent,
        terminal=terminal,
        source_commit=commits[-1],
    )
    before = (
        parent.digest,
        parent.stock_boundary_digest,
        parent.stock_lineage_digest,
        parent.clone_coverage().certificate.digest,
        parent.phase_point,
        parent.traversal.canonical_bytes,
        tuple(passage.canonical_bytes for passage in parent.passages),
        tuple(operation.canonical_bytes for operation in parent.operations),
    )

    def reject_swept_prefix(
        self: MotionCertifier,
        **kwargs: object,
    ) -> NoReturn:
        raise UnresolvedMotionEventError("injected retrace theorem failure")

    monkeypatch.setattr(
        MotionCertifier,
        "certify_swept_prefix_segment",
        reject_swept_prefix,
    )
    with pytest.raises(
        UnresolvedMotionEventError,
        match="injected retrace theorem failure",
    ):
        evaluator.evaluate(parent, decision)

    after = (
        parent.digest,
        parent.stock_boundary_digest,
        parent.stock_lineage_digest,
        parent.clone_coverage().certificate.digest,
        parent.phase_point,
        parent.traversal.canonical_bytes,
        tuple(passage.canonical_bytes for passage in parent.passages),
        tuple(operation.canonical_bytes for operation in parent.operations),
    )
    assert after == before
~~~

Repeat the same complete before/after assertion in four focused tests that
patch, respectively, `GougeContainment.certify_segment`,
`Stock2Area.deplete`, `CoverageLedger.add_sweep`, and
`GenerationState.build` to raise their existing named failure. Also test stale
parent, changed depletion policy, wrong full cap, and cross-wired source
lineage.

- [ ] **Step 4: Run transaction tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_transaction.py -n auto -q
~~~

Expected: collection fails because the retrace transaction module is absent.

- [ ] **Step 5: Implement the immutable transaction**

~~~python
@dataclass(frozen=True)
class RouteRetraceTransaction:
    parent_state_digest: IdentityDigest
    decision: RouteRetraceDecision
    segment_witness: ReplayLateralWitness
    result_state_digest: IdentityDigest

    @classmethod
    def build(
        cls,
        *,
        parent_state_digest: IdentityDigest,
        decision: RouteRetraceDecision,
        segment_witness: ReplayLateralWitness,
        result_state_digest: IdentityDigest,
    ) -> Self:
        return cls(
            parent_state_digest,
            decision,
            segment_witness,
            result_state_digest,
        )
~~~

Validate exact types; require a `RetraceSegmentOperation` whose decision matches; call `ReplayLateralWitness.validate`; require the fixed swept-prefix strategy/theorem/two strata; require certify-before-deplete lineage; and encode all four fields under `route-retrace-transaction-v1`.

- [ ] **Step 6: Implement the short-lived evaluator authority**

~~~python
@dataclass(frozen=True)
class RouteRetraceEvaluator:
    evaluator: CandidateEvaluator

    def __post_init__(self) -> None:
        if type(self.evaluator) is not CandidateEvaluator:
            raise InvalidRouteRetraceTransactionError(
                "retrace evaluator requires one exact candidate authority.",
            )

    @classmethod
    def build(cls, *, evaluator: CandidateEvaluator) -> Self:
        return cls(evaluator)

    def evaluate(
        self,
        state: GenerationState,
        decision: RouteRetraceDecision,
    ) -> RouteRetraceTransaction:
        transaction, _ = self._evaluate_trial(state, decision)
        return transaction

    def commit(
        self,
        state: GenerationState,
        transaction: RouteRetraceTransaction,
    ) -> GenerationState:
        if type(transaction) is not RouteRetraceTransaction:
            raise InvalidRouteRetraceTransactionError(
                "retrace commit requires one exact transaction.",
            )
        if transaction.parent_state_digest != state.digest:
            raise StaleRouteRetraceTransactionError(
                "retrace transaction parent is no longer authoritative.",
            )
        reproduced, child = self._evaluate_trial(state, transaction.decision)
        if reproduced.canonical_bytes != transaction.canonical_bytes:
            raise InvalidRouteRetraceTransactionError(
                "retrace transaction differs from independent commit replay.",
            )
        return child
~~~

`_evaluate_trial` must:

1. require exact `GenerationState` and `RouteRetraceDecision`;
2. call the existing candidate evaluator's physical-parent authority;
3. retrieve `state.operations[decision.source_operation_index]` and require it is the final operation;
4. build `RetraceSegmentOperation` from that owned source;
5. fork stock and clone coverage;
6. call `evaluate_retrace_segment_trial` with operation index `len(state.operations)` and user cap as both user/effective cap;
7. build the child with identical traversal/passages and one appended operation; and
8. seal `RouteRetraceTransaction` from the parent, witness, and child identities.

- [ ] **Step 7: Prove independent byte reproduction**

~~~python
def test_route_retrace_commit_reproduces_transaction_and_child_bytes(
    task13f: Task13FFixture,
) -> None:
    """Repeat the bounded physical proof before publishing its child."""
    parent, terminal, commits = task13f_route_one_terminal(task13f)
    evaluator = RouteRetraceEvaluator.build(evaluator=task13f.evaluator)
    decision = task13f_retrace_decision(
        physical=parent,
        terminal=terminal,
        source_commit=commits[-1],
    )
    transaction = evaluator.evaluate(parent, decision)
    child = evaluator.commit(parent, transaction)
    reproduced = evaluator.evaluate(parent, decision)

    assert reproduced.canonical_bytes == transaction.canonical_bytes
    assert child.digest == transaction.result_state_digest
~~~

- [ ] **Step 8: Run focused GREEN and strict gates**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_route_retrace_transaction.py tests/adaptive/test_route_retrace_state.py tests/adaptive/test_retrace_segment_trial.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
~~~

Expected: all commands pass.

- [ ] **Step 9: Commit and publish Task 4**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/retrace_transaction.py src/compas_cgal/adaptive/errors.py tests/adaptive/test_route_retrace_transaction.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "feat(adaptive): commit retrace transaction"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 5: Authenticate the route boundary and ordered cross-axis commit

**Files:**

- Modify: `src/compas_cgal/adaptive/generator.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/task13f_fixture.py`
- Create: `tests/adaptive/test_route_retrace_generator.py`

**Interfaces:**

- Consumes: `MatTraversalState.activate_next()`, `TraversalCommit`, `ZeroGuideLinkTransaction`, `RouteRetraceEvaluator`, and `GenerationState`.
- Produces: `_route_retrace_required(terminal, activated) -> bool`, `_derive_route_retrace_decision(...) -> RouteRetraceDecision`, `RouteRetraceCommit.build(...)`, `ContinuationCommit = TraversalCommit | RouteRetraceCommit`, and generator boundary orchestration.

- [ ] **Step 1: Add the cross-axis error**

~~~python
class InvalidRouteRetraceCommitError(ValueError):
    """Raised when physical retrace and global route activation disagree."""
~~~

- [ ] **Step 2: Add a reusable valid retrace-prefix test helper**

The route-1 helper already owns the established prefix. Add a second helper
that exercises the production decision, transaction, commit, and first real
route-2 acceptance, returning the existing domain artifact rather than a new
result type:

~~~python
def task13f_retrace_continuation(
    fixture: Task13FFixture,
) -> GenerationContinuation:
    physical, terminal, traversal_commits = task13f_route_one_terminal(fixture)
    activated = terminal.activate_next()
    source_commit = traversal_commits[-1]
    decision = _derive_route_retrace_decision(
        physical=physical,
        terminal=terminal,
        activated=activated,
        source_commit=source_commit,
    )
    evaluator = RouteRetraceEvaluator.build(evaluator=fixture.evaluator)
    transaction = evaluator.evaluate(physical, decision)
    physical_after = evaluator.commit(physical, transaction)
    retrace_commit = RouteRetraceCommit.build(
        physical_before=physical,
        traversal_before=terminal,
        source_commit=source_commit,
        transaction=transaction,
        physical_after=physical_after,
        traversal_after=activated,
    )
    physical_final, traversal_final, route_two_commit = (
        advance_active_candidate_family(
            evaluator=fixture.evaluator,
            physical=physical_after,
            traversal=activated,
        )
    )
    continuation = GenerationContinuation.build(
        launch_transaction=fixture.launch_transaction,
        physical=physical_final,
        traversal=traversal_final,
        commits=(
            *traversal_commits,
            retrace_commit,
            route_two_commit,
        ),
    )
    assert tuple(type(commit) for commit in continuation.commits) == (
        TraversalCommit,
        TraversalCommit,
        RouteRetraceCommit,
        TraversalCommit,
    )
    return continuation
~~~

Before deriving retrace, the helper also asserts both baseline traversal commit
digests from Global Constraints. It is a consumer-boundary fixture: if any
production cross-axis contract is wrong, helper construction fails loud.

- [ ] **Step 3: Write RED exact-trigger and causal-derivation tests**

~~~python
def test_task13f_route_trigger_distinguishes_incident_and_nonincident_edges(
    task13f: Task13FFixture,
) -> None:
    """Use stable MAT node identity, never coordinates, for route transport."""
    physical, route_one_terminal, commits = task13f_route_one_terminal(task13f)
    route_two_active = route_one_terminal.activate_next()

    assert not _route_retrace_required(
        commits[0].traversal_after,
        commits[0].traversal_after.activate_next(),
    )
    assert _route_retrace_required(route_one_terminal, route_two_active)
    assert physical.phase_point == Point2[WorldXY].build(5.0, 1.0)
~~~

For `_derive_route_retrace_decision`, mutate each supplied preimage independently: final source commit, global child, activated state, route index, node, source transaction variant, source witness, source operation ordinal/digest, terminal flag, neck scope, cap, and physical endpoint. Every case must raise `InvalidRouteRetraceDecisionError` or `UnsupportedRouteRetraceError` before candidate materialization.

- [ ] **Step 4: Run trigger tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_generator.py -n auto -q
~~~

Expected: collection fails because the route trigger and retrace commit do not exist.

- [ ] **Step 5: Implement the exact trigger and causal decision derivation**

~~~python
def _route_retrace_required(
    terminal: MatTraversalState,
    activated: MatTraversalState,
) -> bool:
    if (
        type(terminal) is not MatTraversalState
        or type(activated) is not MatTraversalState
        or terminal.active_route_index is None
        or not terminal.active_cursor.terminal
        or activated != terminal.activate_next()
        or activated.active_route_index is None
    ):
        raise InvalidRouteRetraceDecisionError(
            "route trigger requires one exact terminal-to-active transition.",
        )
    completed = terminal.active_cursor.route_step
    next_step = activated.active_cursor.route_step
    return completed.exit_node_id != next_step.entry_node_id
~~~

`_derive_route_retrace_decision` receives `physical`, `terminal`, `activated`, and `source_commit`. It must prove every causal item listed in the approved design before calling `RouteRetraceDecision.build`. In particular:

~~~python
source_transaction = source_commit.transaction
if type(source_transaction) is not ZeroGuideLinkTransaction:
    raise UnsupportedRouteRetraceError(
        "nonincident route requires a final zero-guide source.",
    )
source_index = len(physical.operations) - 1
source_operation = physical.operations[source_index]
if (
    source_commit.physical_child_digest != physical.digest
    or source_commit.traversal_after != terminal
    or source_transaction.result_state_digest != physical.digest
    or source_transaction.segment_witness.operation != source_operation
    or type(source_operation) is not AdvanceSegmentOperation
    or source_operation.motion.end != physical.phase_point
):
    raise InvalidRouteRetraceDecisionError(
        "route retrace source is not the final physical/global commit.",
    )
~~~

Then require no-neck, full-cap, terminal traversal, exact nonincidence, and build the decision from owned preimages only.

- [ ] **Step 6: Implement `RouteRetraceCommit`**

~~~python
@dataclass(frozen=True)
class RouteRetraceCommit:
    physical_parent_digest: IdentityDigest
    traversal_before: MatTraversalState
    transaction: RouteRetraceTransaction
    physical_child_digest: IdentityDigest
    traversal_after: MatTraversalState
    source_commit_digest: IdentityDigest

    @classmethod
    def build(
        cls,
        *,
        physical_before: GenerationState,
        traversal_before: MatTraversalState,
        source_commit: TraversalCommit,
        transaction: RouteRetraceTransaction,
        physical_after: GenerationState,
        traversal_after: MatTraversalState,
    ) -> Self:
        decision = _derive_route_retrace_decision(
            physical=physical_before,
            terminal=traversal_before,
            activated=traversal_after,
            source_commit=source_commit,
        )
        if transaction.decision != decision:
            raise InvalidRouteRetraceCommitError(
                "retrace transaction contradicts route activation.",
            )
        if (
            transaction.parent_state_digest != physical_before.digest
            or transaction.result_state_digest != physical_after.digest
            or physical_after.operations
            != physical_before.operations
            + (transaction.segment_witness.operation,)
            or physical_after.traversal != physical_before.traversal
            or physical_after.passages != physical_before.passages
        ):
            raise InvalidRouteRetraceCommitError(
                "retrace commit breaks physical or held-state lineage.",
            )
        return cls(
            physical_before.digest,
            traversal_before,
            transaction,
            physical_after.digest,
            traversal_after,
            source_commit.digest,
        )
~~~

`__post_init__` repeats all closed type/digest/local cross-links. Canonical bytes use `route-retrace-commit-v1`.

- [ ] **Step 7: Replace blind activation with an exact boundary state machine**

Retain automatic activation only for terminal completion and incident next routes:

~~~python
def _activate_completed_incident_routes(
    traversal: MatTraversalState,
) -> MatTraversalState:
    current = traversal
    while current.active_route_index is not None and current.active_cursor.terminal:
        activated = current.activate_next()
        if (
            activated.active_route_index is not None
            and _route_retrace_required(current, activated)
        ):
            return current
        current = activated
    return current
~~~

Add:

~~~python
ContinuationCommit: TypeAlias = TraversalCommit | RouteRetraceCommit
~~~

Update `GenerationContinuation.commits` and validation as one ordered state machine:

1. `TraversalCommit` must consume the expected physical/global parents and advance exactly one cursor.
2. Normalize only incident or final activations.
3. If normalization stops at a nonincident terminal state, the next union member must be exactly one `RouteRetraceCommit`.
4. `RouteRetraceCommit` must consume the current physical digest, terminal state, and immediately preceding `TraversalCommit.digest`; it publishes exactly `terminal.activate_next()`.
5. No `TraversalCommit` is admitted while retrace is pending; no `RouteRetraceCommit` is admitted otherwise.
6. Final `physical` and `traversal` must equal the fold result.

- [ ] **Step 8: Write RED continuation-order tests, then make them GREEN**

Construct valid prefixes and mutate the commit tuple:

~~~python
@pytest.mark.parametrize(
    "mutation",
    ("missing", "duplicate", "delayed", "reordered"),
)
def test_continuation_rejects_invalid_retrace_order(
    task13f: Task13FFixture,
    mutation: str,
) -> None:
    """Keep physical return adjacent to the branch activation it authorizes."""
    valid = task13f_retrace_continuation(task13f)
    route_zero, route_one, retrace, route_two = valid.commits
    variants: dict[str, tuple[ContinuationCommit, ...]] = {
        "missing": (route_zero, route_one, route_two),
        "duplicate": (
            route_zero,
            route_one,
            retrace,
            retrace,
            route_two,
        ),
        "delayed": (route_zero, route_one, route_two, retrace),
        "reordered": (route_zero, retrace, route_one, route_two),
    }
    commits = variants[mutation]

    with pytest.raises(InvalidRouteRetraceCommitError):
        GenerationContinuation.build(
            launch_transaction=valid.launch_transaction,
            physical=valid.physical,
            traversal=valid.traversal,
            commits=commits,
        )
~~~

Add assertions that route `0 -> 1` adds no retrace commit and route `1 -> 2` adds exactly one.

- [ ] **Step 9: Integrate retrace before route-2 materialization**

After each accepted traversal commit:

~~~python
commits.append(commit)
traversal = _activate_completed_incident_routes(traversal)
if traversal.active_route_index is not None and traversal.active_cursor.terminal:
    activated = traversal.activate_next()
    decision = _derive_route_retrace_decision(
        physical=physical,
        terminal=traversal,
        activated=activated,
        source_commit=commit,
    )
    retrace_transaction = retrace_evaluator.evaluate(physical, decision)
    physical_after = retrace_evaluator.commit(
        physical,
        retrace_transaction,
    )
    retrace_commit = RouteRetraceCommit.build(
        physical_before=physical,
        traversal_before=traversal,
        source_commit=commit,
        transaction=retrace_transaction,
        physical_after=physical_after,
        traversal_after=activated,
    )
    physical = physical_after
    traversal = activated
    commits.append(retrace_commit)
~~~

Build `RouteRetraceEvaluator` once beside the candidate evaluator, not once per boundary. Candidate materialization remains unchanged and runs only after this block.

- [ ] **Step 10: Prove route-2 work begins from the restored phase**

Wrap `materialize_active_candidate_family` and `advance_active_candidate_family` in the integration test. Record the first route-2 physical parent and the first accepted route-2 commit; stop with a test-local sentinel only on the second route-2 dispatch:

~~~python
assert route_two_parent.phase_point == source_operation.motion.start
assert route_two_parent.operations[-1] == retrace_commit.transaction.segment_witness.operation
assert route_two_family_digest == (
    "56d92fcf2089c1e771a9add79bd356d172904a69e6c7acb18ab02522f4ed093b"
)
assert type(first_route_two_commit.transaction) is CandidateTransaction
assert first_route_two_commit.traversal_before.active_route_index == 2
assert (
    first_route_two_commit.traversal_after.active_cursor.cursor.cursor_identity.hex()
    != "def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0"
)
~~~

Also assert no route-2 candidate evaluator call observed phase `Point2[WorldXY].build(5.0, 1.0)`.

- [ ] **Step 11: Run focused GREEN and structural regressions**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_route_retrace_generator.py tests/adaptive/test_generator.py tests/adaptive/test_route_retrace_transaction.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
wc -l src/compas_cgal/adaptive/generator.py src/compas_cgal/adaptive/retrace_transaction.py
~~~

Expected: tests and static gates pass; the line count is recorded for
morphology review, and the focused retrace transaction file—not generator
or replay—owns physical evaluation.

- [ ] **Step 12: Commit and publish Task 5**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/generator.py src/compas_cgal/adaptive/errors.py tests/adaptive/task13f_fixture.py tests/adaptive/test_route_retrace_generator.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "feat(adaptive): cross route boundary"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 6: Reconstruct retrace under fresh operation replay

**Files:**

- Create: `src/compas_cgal/adaptive/route_retrace_replay.py`
- Modify: `src/compas_cgal/adaptive/replay.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Create: `tests/adaptive/test_route_retrace_replay.py`

**Interfaces:**

- Consumes: `RetraceSegmentOperation`, fresh `AdvanceSegmentOperation` reconstruction, `evaluate_retrace_segment_trial(...)`, and existing replay candidate/pairing maps.
- Produces: retrace-aware grammar, continuity, candidate enumeration, pairing,
  fresh global route-boundary replay, and physical replay.
  `ReplayCertificate` remains terminal-only.

- [ ] **Step 1: Add the replay-specific error**

~~~python
class ReplayRouteRetraceError(ValueError):
    """Raised when fresh retrace grammar, lineage, or proof diverges."""
~~~

- [ ] **Step 2: Write RED grammar and candidate-index tests**

Create one partial operation prefix ending after the first accepted route-2 pair. Test:

~~~python
def _pocket() -> Polygon:
    return Polygon([[x, y, 0.0] for x, y in TASK13F_OUTER])


def _replay(
    identity: InputIdentity,
    operations: tuple[CanonicalOperation, ...],
) -> None:
    replay_certificate(
        input_identity=identity,
        pocket=_pocket(),
        holes=(),
        cut_plane=identity.cut_plane,
        tool_radius=identity.tool_radius,
        user_cap=identity.user_cap,
        entry=identity.entry,
        operations=operations,
        candidate_policy=identity.candidate_policy,
        neck_policy=identity.neck_policy,
        depletion_policy=identity.depletion_policy,
        traversal_policy=identity.traversal_policy,
        cut_direction_policy=identity.cut_direction_policy,
    )


@pytest.fixture
def task13f_retrace_prefix(
    task13f: Task13FFixture,
) -> tuple[CanonicalOperation, ...]:
    continuation = task13f_retrace_continuation(task13f)
    return continuation.physical.operations


def test_fresh_candidate_reconstruction_skips_retrace(
    task13f_retrace_prefix: tuple[CanonicalOperation, ...],
    task13f: Task13FFixture,
) -> None:
    """Keep the noncandidate physical return out of MAT candidate chronology."""
    axis = task13f.seeded_traversal.authority.axis
    candidates, _, inventory = _replay_candidate_stream(
        axis=axis,
        operations=task13f_retrace_prefix,
        user_cap=task13f.identity.user_cap,
        candidate_policy=task13f.identity.candidate_policy,
        neck_policy=task13f.identity.neck_policy,
        traversal_policy=task13f.identity.traversal_policy,
        cut_direction_policy=task13f.identity.cut_direction_policy,
    )
    candidate_operations = tuple(
        operation
        for operation in task13f_retrace_prefix[2:]
        if type(operation) in (
            CutFullCircleOperation,
            AdvanceSegmentOperation,
        )
    )

    assert len(candidates) == len(candidate_operations)
    assert all(
        type(candidate) in (MiddleCurveCandidate, ZeroGuideLinkCandidate)
        for candidate in candidates
    )
    assert type(inventory) is NeckInventory
~~~

Parameterize missing, duplicated, reordered, forward-directed, source-ordinal, source-digest, cap, and relabelled operation mutations. Require `ReplayRouteRetraceError` before physical mutation.

- [ ] **Step 3: Run replay grammar tests and confirm RED**

Run:

~~~bash
pixi run pytest tests/adaptive/test_route_retrace_replay.py -n auto -q
~~~

Expected: retrace is rejected as a foreign operation.

- [ ] **Step 4: Extend replay grammar and exact continuity**

In `_validate_grammar` admit retrace only when:

~~~python
if type(operation) is RetraceSegmentOperation:
    if lateral_index == 0:
        raise ReplayRouteRetraceError(
            "route retrace cannot be the entry lateral operation.",
        )
    source_index = operation.decision.source_operation_index
    operation_index = lateral_index + 2
    source = operations[source_index]
    if (
        source_index != operation_index - 1
        or type(source) is not AdvanceSegmentOperation
        or hashlib.sha256(source.canonical_bytes).digest()
        != bytes(operation.decision.source_operation_digest)
    ):
        raise ReplayRouteRetraceError(
            "route retrace does not immediately follow its named source.",
        )
    lateral_index += 1
    continue
~~~

This branch naturally rejects adjacency, a preceding link/circle, and using retrace as the hold-link predecessor of a circle. In `_validate_continuity` require its start equals current phase, exact reversed source endpoints, equal cut depth/cap, then advance current phase to retrace end.

- [ ] **Step 5: Skip retrace in candidate reconstruction and pairing**

In `_replay_candidate_stream`, retain the freshly built `NeckInventory` in the
return tuple, then `continue` on retrace without appending or changing
`current_by_edge` or passages. In `_pair_lateral_operations`, validate the
immediately preceding source operation and digest, require no pending link,
and `continue` without incrementing `candidate_index` or adding a pairing
entry.

Add an assertion around the recorded Task 13F prefix:

~~~python
paired = _pair_lateral_operations(
    operations=operations,
    candidates=candidates,
)
retrace_index = next(
    index
    for index, operation in enumerate(operations)
    if type(operation) is RetraceSegmentOperation
)
assert retrace_index not in paired
assert retrace_index - 1 in paired
assert retrace_index + 1 in paired
assert retrace_index + 2 in paired
~~~

- [ ] **Step 6: Replay the exact global route boundary from the fresh MAT**

Create `route_retrace_replay.py` with one orchestration helper; it makes no
physical decision:

~~~python
def _replay_route_boundaries(
    *,
    axis: MedialAxis,
    inventory: NeckInventory,
    traversal_policy: TraversalPolicy,
    operations: tuple[CanonicalOperation, ...],
    paired_candidates: dict[int, TraversalCandidate],
) -> None:
    first = paired_candidates.get(2)
    if type(first) is not MiddleCurveCandidate:
        raise ReplayRouteRetraceError(
            "fresh route replay requires the entry-circle candidate.",
        )
    sample_index = TraversalSampleIndex.build(axis)
    graph = TraversalGraph.from_axis(
        axis=axis,
        sample_index=sample_index,
    )
    edge_id = MatEdgeId(bytes(first.traversal_decision.edge_id))
    edge = next(
        graph_edge for graph_edge in graph.edges
        if graph_edge.edge_id == edge_id
    )
    cursor_before = first.traversal_decision.cursor_before
    if cursor_before == edge.source_cursor_id:
        entry_node_id = edge.source_node_id
    elif cursor_before == edge.target_cursor_id:
        entry_node_id = edge.target_node_id
    else:
        raise ReplayRouteRetraceError(
            "entry candidate does not begin at a fresh MAT endpoint.",
        )
    traversal = MatTraversalState.seed(
        axis=axis,
        inventory=inventory,
        policy=traversal_policy,
        entry_edge_id=edge_id,
        entry_node_id=entry_node_id,
    )
    pending: tuple[
        MatTraversalState,
        MatTraversalState,
        int,
        AdvanceSegmentOperation,
    ] | None = None
    for operation_index, operation in enumerate(operations[2:], start=2):
        if type(operation) is RetraceSegmentOperation:
            if pending is None:
                raise ReplayRouteRetraceError(
                    "fresh route replay found an unexpected retrace.",
                )
            terminal, activated, source_index, source = pending
            decision = operation.decision
            if (
                source_index != operation_index - 1
                or decision.completed_route_index
                != terminal.active_route_index
                or decision.activated_route_index
                != activated.active_route_index
                or decision.completed_exit_node_id
                != RouteNodeId(
                    bytes(terminal.active_cursor.route_step.exit_node_id),
                )
                or decision.activated_entry_node_id
                != RouteNodeId(
                    bytes(activated.active_cursor.route_step.entry_node_id),
                )
                or decision.terminal_traversal_digest != terminal.digest
                or decision.activated_traversal_digest != activated.digest
                or decision.source_operation_index != source_index
                or decision.source_operation_digest
                != IdentityDigest(
                    hashlib.sha256(source.canonical_bytes).digest(),
                )
            ):
                raise ReplayRouteRetraceError(
                    "fresh route activation contradicts its retrace decision.",
                )
            traversal = activated
            pending = None
            continue
        if pending is not None:
            raise ReplayRouteRetraceError(
                "nonincident route activation is missing its retrace.",
            )
        if type(operation) is LinkSegmentOperation:
            continue
        candidate = paired_candidates.get(operation_index)
        if candidate is None:
            continue
        traversal = traversal.advance(candidate)
        if not traversal.active_cursor.terminal:
            continue
        activated = traversal.activate_next()
        if activated.active_route_index is None:
            traversal = activated
            continue
        completed_step = traversal.active_cursor.route_step
        activated_step = activated.active_cursor.route_step
        if completed_step.exit_node_id == activated_step.entry_node_id:
            traversal = activated
            continue
        if type(operation) is not AdvanceSegmentOperation:
            raise ReplayRouteRetraceError(
                "fresh nonincident route has no admitted advancing source.",
            )
        pending = (
            traversal,
            activated,
            operation_index,
            operation,
        )
    if pending is not None:
        raise ReplayRouteRetraceError(
            "operation prefix ends before its required route retrace.",
        )
~~~

Update the sole `replay_certificate` caller explicitly:

~~~python
candidates, current_by_edge, inventory = _replay_candidate_stream(...)
paired_candidates = _pair_lateral_operations(
    operations=operations,
    candidates=candidates,
)
_replay_route_boundaries(
    axis=axis,
    inventory=inventory,
    traversal_policy=traversal_policy,
    operations=operations,
    paired_candidates=paired_candidates,
)
trace = _replay_fresh_state(...)
~~~

`replay.py` only imports and calls this helper. It therefore runs after
candidate pairing and before
`_replay_fresh_state`. It
freshly reproduces route indices, node identities, terminal and activated
state digests, source ordinal, and source operation digest. The live
`GenerationContinuation` fold from Task 5 independently authenticates the
decision's source commit and transaction digests, which operation-only replay
cannot infer without their preimages.

Add missing, duplicate, delayed, reordered, and altered-activation tests
against this helper. Each must raise `ReplayRouteRetraceError` before
`evaluate_retrace_segment_trial` is called.

- [ ] **Step 7: Replay the physical proof through the shared function**

Add one branch in `_replay_fresh_state` before the full-circle branch:

~~~python
if type(operation) is RetraceSegmentOperation:
    segment_witness = evaluate_retrace_segment_trial(
        containment_authority=containment,
        stock=stock,
        coverage=coverage,
        operation_index=operation_index,
        operation=operation,
        tool_radius=tool_radius,
        user_cap=user_cap,
        effective_cap=user_cap,
        depletion_policy=depletion_policy,
    )
    lateral_witnesses.append(segment_witness)
    current_stock_state = MotionCertifier.build(
        stock=stock,
        tool_radius=tool_radius,
    )
    continue
~~~

The operation factory and grammar have already proved the exact full-cap source reversal; `ReplayLateralWitness` independently proves the physical theorem and owner lineages.

- [ ] **Step 8: Capture the partial fresh trace without weakening terminal replay**

Use the same tracking pattern as the established zero-guide replay test: wrap
`_require_terminal_replay`, retain its `trace` argument, and call the real
terminal gate:

~~~python
def test_fresh_replay_reproduces_live_retrace_before_terminal_gate(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rebuild the source first, then reproduce retrace before expected stop."""
    continuation = task13f_retrace_continuation(task13f)
    live_retrace_commit = next(
        commit
        for commit in continuation.commits
        if type(commit) is RouteRetraceCommit
    )
    order: list[str] = []
    traces: list[FreshReplayTrace] = []
    real_match = replay_module._match_zero_guide_candidate
    real_retrace_trial = replay_module.evaluate_retrace_segment_trial
    real_terminal = replay_module._require_terminal_replay

    def tracked_match(**kwargs: object) -> ZeroGuideLinkCandidate:
        candidate = real_match(**kwargs)  # type: ignore[arg-type]
        order.append("source-zero-guide")
        return candidate

    def tracked_retrace_trial(**kwargs: object) -> ReplayLateralWitness:
        order.append("retrace-physical-proof")
        return real_retrace_trial(**kwargs)  # type: ignore[arg-type]

    def tracked_terminal(**kwargs: object) -> None:
        trace = kwargs["trace"]
        assert type(trace) is FreshReplayTrace
        traces.append(trace)
        real_terminal(**kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(
        replay_module,
        "_match_zero_guide_candidate",
        tracked_match,
    )
    monkeypatch.setattr(
        replay_module,
        "evaluate_retrace_segment_trial",
        tracked_retrace_trial,
    )
    monkeypatch.setattr(
        replay_module,
        "_require_terminal_replay",
        tracked_terminal,
    )
    with pytest.raises(ReplayTraversalError):
        _replay(
            task13f.identity,
            continuation.physical.operations,
        )

    assert len(traces) == 1
    assert order.index("source-zero-guide") < order.index(
        "retrace-physical-proof",
    )
    captured_trace = traces[0]
    live_witness = live_retrace_commit.transaction.segment_witness
    fresh_witness = next(
        witness
        for witness in captured_trace.lateral_witnesses
        if type(witness.operation) is RetraceSegmentOperation
    )
    assert fresh_witness.canonical_bytes == live_witness.canonical_bytes
    assert fresh_witness.operation_index == live_witness.operation_index
    assert type(fresh_witness.motion_witness) is SweptPrefixMotionWitness
    assert fresh_witness.motion_witness.event_cell_count == 2
~~~

The wrapper calls the real gate, so the partial prefix still raises the
existing `ReplayTraversalError` or residual coverage error. No
`ReplayCertificate` is manufactured.

- [ ] **Step 9: Run focused GREEN and complete replay regressions**

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_route_retrace_replay.py tests/adaptive/test_replay.py tests/adaptive/test_zero_guide_replay.py tests/adaptive/test_route_retrace_generator.py -n auto -q
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
~~~

Expected: all commands pass.

- [ ] **Step 10: Commit and publish Task 6**

Run:

~~~bash
git diff --check
git add src/compas_cgal/adaptive/route_retrace_replay.py src/compas_cgal/adaptive/replay.py src/compas_cgal/adaptive/errors.py tests/adaptive/test_route_retrace_replay.py
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "feat(adaptive): replay route retrace"
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git status --short
~~~

Expected: one focused commit, matching local/remote SHA, clean worktree.

---

### Task 7: Preserve the negative control, document the proof delta, and close the stage

**Files:**

- Modify: `tests/adaptive/test_generator.py`
- Modify: `tests/adaptive/test_route_retrace_generator.py`
- Modify: `docs/segment_site_mat.md`
- Modify: `docs/continuous_engagement.md`
- Modify: `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`
- Modify: `docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md`
- Modify: `docs/superpowers/specs/2026-08-06-exact-inter-route-retrace-design.md`

**Interfaces:**

- Consumes: the production generator, ordered continuation, fresh replay witness, exact Task 13F fixture, and all measured structural evidence from Tasks 1–6.
- Produces: a retained direct-path negative test, implementation-status evidence, Held–Pfeiffer comparison, current limitation, and a clean published stage checkpoint.

- [ ] **Step 1: Refactor the old 56-gouge test into an explicit direct-path negative**

Do not call the production generator for this negative. Commit route 0 and route 1 through the real production candidate path, manually activate route 2 in the test without retrace, then dispatch the unchanged route-2 family:

~~~python
def test_task13f_direct_nonincident_route_two_remains_56_gouges(
    task13f: Task13FFixture,
) -> None:
    """Preserve the exact failure when physical return is deliberately omitted."""
    physical, route_one_terminal, _ = task13f_route_one_terminal(task13f)
    route_two = route_one_terminal.activate_next()

    with pytest.raises(
        NoFeasibleCandidateError,
        match=(
            "finite candidate family exhausted.*attempts=56; "
            "cap=0; gouge=56; degenerate-link=0"
        ),
    ) as raised:
        advance_active_candidate_family(
            evaluator=task13f.evaluator,
            physical=physical,
            traversal=route_two,
        )

    family = materialize_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=physical,
        traversal=route_two,
    )
    assert route_two.active_cursor.cursor_identity.hex() == (
        "def1bf1471e2df355ba488378ffad7b9e20116ac81909d9f9822a5a62b6abbb0"
    )
    assert len(family) == 56
    assert hashlib.sha256(
        b"".join(bytes(candidate.identity) for candidate in family),
    ).hexdigest() == (
        "56d92fcf2089c1e771a9add79bd356d172904a69e6c7acb18ab02522f4ed093b"
    )
    assert "gouge=56" in str(raised.value)
~~~

Run:

~~~bash
pixi run pytest --testmon tests/adaptive/test_generator.py::test_task13f_direct_nonincident_route_two_remains_56_gouges -n auto -q
~~~

Expected: PASS with the original cursor, counts, digest, and error disposition.

- [ ] **Step 2: Record structural counts and bounded wall time**

Run the real route-1-to-route-2 transition once for evaluation and once for
independent commit. Assert per-layer counters, not elapsed time:

~~~python
assert one_decision_derivation_counts == {
    "node_identity_comparison": 1,
    "adjacent_source_lookup": 1,
    "source_digest_verification": 1,
}
assert evaluation_counts == {
    "segment_containment": 1,
    "swept_prefix": 1,
    "swept_prefix_strata": 2,
    "stock_depletion": 1,
    "coverage_update": 1,
    "generic_event_certification": 0,
    "mat_rebuild": 0,
    "candidate_probe_before_retrace": 0,
}
assert cumulative_after_independent_physical_commit == {
    "segment_containment": 2,
    "swept_prefix": 2,
    "swept_prefix_strata": 4,
    "stock_depletion": 2,
    "coverage_update": 2,
    "generic_event_certification": 0,
    "mat_rebuild": 0,
    "candidate_probe_before_retrace": 0,
}
~~~

The `RouteRetraceCommit.build` test separately proves a second complete
decision derivation from preimages; do not merge that count into physical
trial counts. Record `time.perf_counter()` around the bounded transition and
report the observed value in docs without a pass/fail threshold. Keep the
established comparison that the generic oracle exceeded 10 minutes while the
prior swept-prefix probe took 0.040943 s; label both as fixture observations,
not universal performance claims.

- [ ] **Step 3: Update the primary MkDocs explanation while evidence is live**

In `docs/segment_site_mat.md` add:

1. the exact node-identity trigger;
2. the final zero-guide source proof and exact reverse derivation;
3. the decision → operation → transaction → commit pipeline;
4. one ordered continuation timeline;
5. the evaluation/commit structural counts;
6. the fresh-replay separation between candidate reconstruction and physical proof;
7. the retained 56-gouge direct-path negative;
8. the newly observed fail-closed boundary after the first route-2 acceptance; and
9. the initial comparison table:

~~~markdown
| Capability | Held–Pfeiffer (2025) | Exact-certified adaptive MAT | Current evidence |
| --- | --- | --- | --- |
| Engagement control | Full tool-engagement-angle control | Exact continuous cap certificates per accepted motion | Stronger bounded proof contract |
| Route construction | Published trochoidal pocket strategy | Exact MAT topology, causal neck state, finite candidate identity | Different architecture; bounded Task 13F evidence |
| Inter-route transport | Tool-path strategy described at algorithm level | Exact source-authenticated cut-depth retrace | Live + fresh replay at one proved boundary |
| Failure semantics | Algorithm/performance evaluation | Named fail-closed proof errors with retained negative controls | Stronger diagnostic evidence |
| End-to-end maturity | Published complete method and measurements | Route-2 progress; terminal coverage and matched benchmark pending | Held–Pfeiffer stronger today |
| Matched performance | Published evaluation | Structural count gate; matched terminal benchmark pending | No parity claim yet |
~~~

State the maturity conclusion exactly:

~~~text
stronger exact bounded proof contract;
weaker end-to-end and measured-performance evidence
~~~

- [ ] **Step 4: Update companion docs and continuation records**

In `docs/continuous_engagement.md` explain why swept-prefix is the correct theorem for retracing an already swept corridor and why generic event partition is excluded.

In the causal and zero-guide plans, append dated implementation evidence linking the retrace commit chain and explicitly preserving their earlier operation/transaction digests.

In the retrace design:

- change status to `implemented and boundedly verified 2026-08-07`;
- list the exact implementation commit SHAs;
- record live/fresh witness digest and structural counts;
- identify the exact next boundary;
- keep tri-dexel removal/thermal replay, kinematic SDF optimization, terminal coverage, and matched Held–Pfeiffer performance as later tasks.

- [ ] **Step 5: Capture the non-obvious implementation insights**

Use `insight-capture` without waiting for acknowledgement. Capture at minimum:

- global DFS route order is not a cutter trajectory;
- source operation bytes are reusable physical provenance, not coordinates to rematch;
- retrace is a physical mutation but a global/local traversal hold;
- a single ordered commit union is what makes cross-axis chronology replayable;
- repeated stock and coverage lineage remains material evidence even when the swept set is geometrically unchanged; and
- preserving the direct-path negative prevents the new route from weakening containment.

Do not duplicate facts already present in root `CLAUDE.md`; put durable geometric/proof knowledge in the relevant MkDocs page and concise repository navigation guidance in the local policy only if the skill directs it.

- [ ] **Step 6: Run the coherent stage gates**

Run:

~~~bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run affected
pixi run pytest tests/adaptive -n auto -q
pixi run -e docs docs
git diff --check
~~~

Expected: every command passes; no skip, skipif, xfail, warning-suppression, fallback, or changed reference value appears.

- [ ] **Step 7: Run an independent Critical/Important review**

Review the complete implementation against the approved spec and report only:

- Critical: unsound exact proof, mutable-authority leak, replay divergence, or missing acceptance gate;
- Important: contract ambiguity, untested consumer boundary, duplicated deciding call graph, or unsupported claim.

Resolve every Critical and Important finding through RED/GREEN tests, rerun Step 6, and retain the review disposition in the handoff.

- [ ] **Step 8: Commit the documentation checkpoint**

Run:

~~~bash
git diff --check
git add tests/adaptive/test_generator.py tests/adaptive/test_route_retrace_generator.py docs/segment_site_mat.md docs/continuous_engagement.md docs/superpowers/plans/2026-07-28-causal-mat-traversal.md docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md docs/superpowers/specs/2026-08-06-exact-inter-route-retrace-design.md
GIT_AUTHOR_NAME="Jelle Feringa" GIT_AUTHOR_EMAIL="jelleferinga@gmail.com" GIT_COMMITTER_NAME="Jelle Feringa" GIT_COMMITTER_EMAIL="jelleferinga@gmail.com" git commit -m "docs(adaptive): close retrace stage"
~~~

Expected: one focused test/documentation checkpoint by Jelle.

- [ ] **Step 9: Cancel stale CI, push, and verify durable publication**

Run:

~~~bash
git remote -v
for run_id in $(gh run list --repo jf---/compas_cgal --branch codex/exact-certified-adaptive-phase1-t9-zero-guide --limit 5 --json databaseId,status --jq '.[] | select(.status == "in_progress") | .databaseId'); do gh run cancel "$run_id" --repo jf---/compas_cgal; done
git push origin codex/exact-certified-adaptive-phase1-t9-zero-guide
test "$(git rev-parse HEAD)" = "$(git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide | cut -f1)"
git show -s --format='%an <%ae>%n%cn <%ce>%n%H%n%s' HEAD
git status --short
~~~

Expected: explicit `jf---/compas_cgal` target, Jelle author/committer, matching local/remote SHA, and clean worktree.

## Completion Evidence

The stage is complete only when the final handoff includes:

- exact local branch, worktree, and remote SHA;
- focused commit list and clean status;
- route `0 -> 1` incident activation with no retrace;
- route `1 -> 2` exactly one `RouteRetraceCommit`;
- first accepted route-2 candidate beyond the pinned exhaustion cursor;
- live and fresh retrace witness byte equality;
- unchanged local cursor and passage bytes across retrace;
- stock/coverage lineage append counts and shared two-stratum proof counts;
- retained 56-gouge direct-path negative identity;
- full test, Ruff, strict mypy, strict MkDocs, and diff-check results;
- independent review disposition;
- current fail-closed boundary and honest Held–Pfeiffer maturity statement; and
- next ordered tasks: terminal exact generation, matched Held–Pfeiffer performance campaign, tri-dexel thermal/removal integration, then kinematic-SDF optimization.
