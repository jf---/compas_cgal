# Atomic Candidate Transaction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evaluate each link-and-circle candidate on isolated exact state and commit only byte-identical, non-stale acceptance evidence.

**Architecture:** `GenerationState` owns one immutable authoritative snapshot and safe native forks. `CandidateEvaluator` constructs and certifies an isolated link/circle pair, returns a content-addressed `CandidateTransaction`, and independently repeats only the winning evaluation during functional commit.

**Tech Stack:** Python 3.12, COMPAS typed geometry, CGAL exact native owners, frozen dataclasses, CCAN canonical encoding, SHA-256, pytest-xdist, pytest-testmon, Ruff, mypy strict, MkDocs strict.

## Global Constraints

- Exact proof failures stay loud; there is no sampling, fallback, rollback, result wrapper, skip, or xfail.
- The evaluator validates submitted neck scope but never derives causal scope; Task 13 owns that proof.
- Authoritative stock, coverage, traversal, passage, and operation state remain unchanged until explicit commit.
- Every public factory validates exact runtime types and every immutable artifact owns canonical bytes plus SHA-256 identity.
- Every test has a docstring explaining its geometric and proof-boundary context.
- Run every pytest command with `-n auto`; use `--testmon` after code changes.

---

### Task 1: Canonical state primitives

**Files:**

- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `src/compas_cgal/adaptive/neck.py`
- Create: `src/compas_cgal/adaptive/generation_state.py`
- Test: `tests/adaptive/test_transaction.py`

**Interfaces:**

- Consumes: `Stock2Area.fork()`, `CoverageLedger.clone()`,
  `MotionCertifier.build(*, stock: Stock2Area, tool_radius: ToolRadius)`,
  `CanonicalOperation`, `NeckPassage`, and `AdvanceTraversalDecision`.
- Produces: `TraversalCursorState`, `GenerationState`, and named Task 12 integrity failures.

- [x] **Step 1: Write RED state-isolation and digest tests**

```python
def test_generation_state_clones_native_owners_and_binds_every_lineage() -> None:
    """Keep caller aliases and every accepted lineage outside authority."""
    fixture = _transaction_fixture()
    state = fixture.state
    digest = state.digest

    fixture.source_stock.deplete(
        fixture.second.motion,
        fixture.identity.tool_radius,
        fixture.identity.depletion_policy,
    )

    assert state.digest == digest
    assert state.stock_lineage != fixture.source_stock.lineage
```

Add exact failures for a foreign tool radius, stock/coverage chronology
divergence, duplicate passage scopes, phase/prefix mismatch, and traversal
cursor mismatch.

- [x] **Step 2: Run RED**

Run:

```bash
pixi run pytest tests/adaptive/test_transaction.py -n auto -q
```

Expected: collection fails because `generation_state` and its errors do not
exist.

- [x] **Step 3: Add named errors and canonical neck passage bytes**

Add:

```python
class InvalidGenerationStateError(ValueError):
    """An authoritative generation snapshot violates cross-owner chronology."""


class CandidateStateMismatchError(ValueError):
    """A candidate transition does not begin at the authoritative state."""


class InvalidCandidateTransactionError(ValueError):
    """Accepted candidate evidence is malformed, cross-wired, or reordered."""


class StaleCandidateTransactionError(RuntimeError):
    """A transaction parent digest no longer names authoritative state."""


class CandidateSelectionError(ValueError):
    """Accepted transactions cannot form one deterministic winner set."""
```

Give `ClassifiedNeck` and `NeckPassage` canonical CCAN records that bind owner,
evidence, class, cut, orientation, and passage state.

- [x] **Step 4: Implement immutable state**

Create frozen `TraversalCursorState` with `component_id: ComponentId`,
`edge_id: EdgeId`, `branch_id: BranchId`, `cursor: CursorIdentity`, and
`terminal: bool`. Its exact interfaces are:

- `TraversalCursorState.before(decision: AdvanceTraversalDecision) -> Self`;
- `advance(decision: AdvanceTraversalDecision) -> Self`;
- `canonical_bytes -> bytes`.

Create frozen `GenerationState` with private `Stock2Area` and
`CoverageLedger` owners plus `tool_radius: ToolRadius`,
`phase_point: Point2[WorldXY]`, `traversal: TraversalCursorState`,
`passages: tuple[NeckPassage, ...]`,
`operations: tuple[CanonicalOperation, ...]`, and an init-false
`_digest: IdentityDigest`. Its exact interfaces are:

- `GenerationState.build(*, stock: Stock2Area, coverage: CoverageLedger,
  tool_radius: ToolRadius, phase_point: Point2[WorldXY],
  traversal: TraversalCursorState, passages: tuple[NeckPassage, ...],
  operations: tuple[CanonicalOperation, ...]) -> Self`;
- `fork_stock() -> Stock2Area`;
- `clone_coverage() -> CoverageLedger`;
- `passage(scope: OrientedNeckScope) -> NeckPassage`;
- `digest -> IdentityDigest`.

`__post_init__` clones both native owners, proves entry/lateral stock lineage
matches coverage and the operation prefix, validates the current phase and
cursor, canonicalizes passage lookup, and derives the state digest.

- [x] **Step 5: Run focused state tests**

Run:

```bash
pixi run pytest tests/adaptive/test_transaction.py -n auto --testmon -q
```

Expected: state tests pass; evaluator tests remain RED.

### Task 2: RED evaluator and atomicity contracts

**Files:**

- Create: `src/compas_cgal/adaptive/transaction.py`
- Test: `tests/adaptive/test_transaction.py`

**Interfaces:**

- Consumes: `GenerationState`, `MiddleCurveCandidate`,
  `ReplayLateralWitness`, `GougeContainment`, `MotionCertifier`,
  `Stock2Area.deplete()`, and `CoverageLedger.add_sweep()`.
- Produces: `CandidateTransaction`, `CandidateEvaluator.evaluate()`,
  `CandidateEvaluator.commit()`, and `select_candidate_transaction()`.

- [x] **Step 1: Add real L-pocket fixture**

Build the exact L-pocket reachable domain, qualified entry, candidate lattice,
first accepted circle, stock/coverage lineage, traversal cursor, and neck
passages. Return a frozen test fixture whose first and second candidates are
fixed by exact progress, guide radius, phase, generator site, and cursor
identity.

- [x] **Step 2: Add RED atomicity tests**

Add:

```python
def test_circle_rejection_after_certified_link_is_atomic() -> None:
    """Discard a passing link when its joint circle exceeds the neck cap."""
    fixture = _transaction_fixture(neck_caps=(90.0, 80.0, 70.0))
    before = fixture.state.digest

    with pytest.raises(EngagementCapExceededError, match="exceeds"):
        fixture.evaluator.evaluate(fixture.state, fixture.second)

    assert fixture.state.digest == before
```

Also add link-cap failure, unresolved-circle propagation, real gouge failure,
deterministic accepted evidence, and proof that evaluation changes neither
stock nor coverage lineage.

- [x] **Step 3: Add RED chronology-mutation tests**

Construct a valid transaction, then use `dataclasses.replace` to:

- give the link motion witness the post-link stock lineage;
- swap link and circle proof bundles;
- cross-wire a containment certificate;
- alter operation indices;
- advance the neck passage in state before commit.

Each mutation must raise its exact Task 12 integrity error.

- [x] **Step 4: Add RED commit and selection tests**

Prove a successful commit returns a new state while preserving the parent,
changes stock and coverage together, appends link then circle, advances phase
and cursor, advances an oriented passage only at commit, rejects stale parent
state, and selects the same winner for every permutation.

- [x] **Step 5: Run RED**

Run:

```bash
pixi run pytest tests/adaptive/test_transaction.py -n auto -q
```

Expected: evaluator/transaction imports or assertions fail while state tests
stay GREEN.

### Task 3: Exact isolated evaluation

**Files:**

- Create: `src/compas_cgal/adaptive/transaction.py`
- Modify: `src/compas_cgal/adaptive/generation_state.py`
- Test: `tests/adaptive/test_transaction.py`

**Interfaces:**

- `CandidateTransaction.build(*, parent_state_digest: IdentityDigest,
  candidate: MiddleCurveCandidate, link_witness: ReplayLateralWitness,
  circle_witness: ReplayLateralWitness,
  traversal_after: TraversalCursorState,
  passage_after: NeckPassage | None,
  result_state_digest: IdentityDigest) -> CandidateTransaction`
- `CandidateEvaluator.build(*, reachable_domain: ReachableDomain,
  tool_radius: ToolRadius, user_cap: EngagementCap,
  candidate_policy: CandidatePolicy, neck_policy: NeckPolicy,
  depletion_policy: DepletionPolicy,
  cut_direction_policy: CutDirectionPolicy, cut_z: CutZ,
  material_side: MaterialSide) -> CandidateEvaluator`
- `CandidateEvaluator.evaluate(state, candidate) -> CandidateTransaction`
- `CandidateEvaluator.commit(state, transaction) -> GenerationState`
- `select_candidate_transaction(transactions) -> CandidateTransaction`

- [x] **Step 1: Implement immutable transaction evidence**

Create frozen `CandidateTransaction` with
`parent_state_digest: IdentityDigest`, `candidate: MiddleCurveCandidate`,
`link_witness: ReplayLateralWitness`,
`circle_witness: ReplayLateralWitness`,
`traversal_after: TraversalCursorState`,
`passage_after: NeckPassage | None`, and
`result_state_digest: IdentityDigest`. In addition to the exact `build`
signature above, expose `canonical_bytes -> bytes` and
`digest -> IdentityDigest`.

Validation binds link/circle operation order, motion kinds and indices,
candidate scope/cap/traversal, direct phase continuity, depletion parent
lineage, coverage parent lineage, certify-before-deplete lineage, and result
state identity.

- [x] **Step 2: Implement one isolated trial**

Private
`_evaluate_trial(state: GenerationState, candidate: MiddleCurveCandidate)`
returns exactly:

```python
tuple[CandidateTransaction, GenerationState]
```

It validates state/candidate traversal and passage, derives the effective cap,
builds link and circle operations, forks stock/coverage, and performs for each
motion:

```text
containment -> event certification -> depletion -> coverage
```

It then builds the next state and transaction from the complete witnesses.

- [x] **Step 3: Implement public evaluation and commit**

`evaluate()` returns the first element of `_evaluate_trial()` and discards the
isolated next-state owner. `commit()` rejects a parent-digest mismatch before
native work, reruns `_evaluate_trial()`, requires byte-identical transaction
evidence, and returns the recomputed next state. It never mutates its input.

- [x] **Step 4: Implement pure winner selection**

Require a nonempty tuple, one parent digest, one candidate policy, and unique
transaction digests. Call `CandidatePolicy.order_candidates()` with each
transaction's `candidate.order_key` and return the first result.

- [x] **Step 5: Run GREEN**

Run:

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_transaction.py -n auto --testmon -q
```

Expected: all Task 12 contracts pass.

### Task 4: Documentation, integration gates, and publication

**Files:**

- Modify: `docs/segment_site_mat.md`
- Modify: `docs/continuous_engagement.md`
- Modify: `docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`
- Modify: `docs/superpowers/plans/2026-07-28-atomic-candidate-transaction.md`

**Interfaces:**

- Consumes: verified Task 12 tests and exact transaction maturity.
- Produces: contemporaneous developer explanation, Held–Pfeiffer comparison
  update, task status, insight log, and published commit.

- [x] **Step 1: Document verified contracts**

Explain joint link/circle atomicity, isolated losing candidates, independent
winner replay, stale-parent rejection, exact chronology, and the remaining
Task 13 causal-scope/traversal boundary. Update comparison-table maturity
without claiming a complete generator or matched performance result.

- [x] **Step 2: Mark this plan accurately**

Change only completed checklist boxes to `[x]`; leave any incomplete
acceptance gate open.

- [x] **Step 3: Run complete bounded verification**

Run:

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run affected
pixi run pytest tests/adaptive -n auto -q
pixi run docs
git diff --check
```

- [x] **Step 4: Commit**

Stage only Task 12 code, tests, and documentation. Commit as Jelle Feringa:

```bash
git commit -m "feat(adaptive): add atomic candidate commit"
```

- [x] **Step 5: Cancel superseded CI, push, and verify**

List branch runs with explicit repository, cancel every `in_progress` run,
push `codex/exact-certified-adaptive-phase1-t9`, and require
`git ls-remote` to equal local `HEAD`.

Task 12 implementation was committed as
`a9fef135f72b4f3aded20ccba60510e90b2ab464`; the branch had no in-progress CI
run, push succeeded without force, and the remote branch resolved to that exact
SHA before this publication note.
