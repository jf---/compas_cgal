# Held Figure 5 2D Path Characterization Implementation Plan

> **status: stopped on 2026-09-05** - Tasks 1-10, Task 11 Steps 1-2, Task 11A,
> and Task 11B are implemented and preserved. The operator stopped Task 11
> Step 3 after 412.91 seconds and retired Steps 3-7 because characterization had
> become a detour from reproducing the paper figures. They are historical open
> work, not the active queue. Continue only with
> `2026-09-05-held-paper-figure-reproduction.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a durable, operation-attributed 2D characterization of the
current default generator on the fixed Held-Pfeiffer Figure 5 pocket and refuse
postprocessor qualification unless every declared Phase 1 criterion closes.

**Architecture:** Generate Figure 5 once, snapshot its complete operation stream,
then run the unchanged guarded engagement replay and one canonical quality
survey. Canonical quality reducers produce both existing `PathQuality` values
and typed operation attribution. A validated characterization may yield an
immutable candidate for later target-specific postprocessor qualification, but
never G-code or machine-release authority.

**Tech Stack:** Python 3.9-compatible typing on the Python 3.12 Pixi environment,
COMPAS geometry, existing CGAL stock/engagement bindings, pytest-xdist,
pytest-testmon, mypy strict, Ruff, and MkDocs.

**Spec:**
`docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md`

## Global Constraints

- Work only in `codex/held-reference-corpus`; never mutate `main` or `master`.
- Use exactly `load_held_reference_case("figure5")` and
  `benchmarks.runner.generate_toolpath`; never tune the reference geometry or
  generator during Phase 1.
- Preserve `stock_2`, `engagement_2`, native `audit_*`,
  `audit_toolpath_engagement`, stock depletion, and all twelve current quality
  thresholds as protected judges.
- Call the exercised engagement evidence guarded replay certification; do not
  imply that Phase 1 runs or qualifies the separate native audit protocol.
- Build one standard-density quality survey. Do not add a second survey,
  alternate decision path, silent fallback, conditional import, skipped test,
  or expected-failure marker.
- Reuse `Millimetre`, `Radian`, `Point2`, `Point3`, `Direction3`, `WorldXY`, and
  `WorldXYZ` from `compas_cgal.adaptive.units`; add only missing observation
  scalar types.
- Every new invariant-bearing record is frozen, disables direct initialization,
  and exposes a validating `.build(...)` factory with a named failure mode.
- Keep snapshot, criterion reduction, characterization, qualification, rendering,
  orchestration, adapter, and CLI responsibilities in separate files.
- Do not add provenance/identity machinery or freeze review inputs. Prove path
  binding by immutable structural values and behavior.
- Use Pixi exclusively. Every pytest invocation includes `-n auto`; after Python
  changes run an affected `--testmon` gate and Ruff format/check before commit.
- Preserve the exact known-red membership: twelve quality cells, one
  deterministic benchmark translation-invariance cell, and four adaptive
  cells. Generated property examples may extend discovery but may not be the
  sole authority for declared-red membership. Never edit a reference assertion
  to make it pass.
- Leave the user-owned untracked `tmp/` tree untouched and unstaged.
- Set both Git author and committer to
  `Jelle Feringa <jelleferinga@gmail.com>` for every commit.

---

### Task 1: Observation units and named failures

**Files:**

- Create: `benchmarks/units.py`
- Modify: `benchmarks/errors.py`
- Create: `tests/benchmarks/test_units.py`
- Create: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: canonical geometry units from `compas_cgal.adaptive.units`.
- Produces:
  - `Degrees = NewType("Degrees", float)`
  - `Seconds = NewType("Seconds", float)`
  - `UnitFraction = NewType("UnitFraction", float)`
  - `MotionCount = NewType("MotionCount", int)`
  - `ToolRadiusMultiple = NewType("ToolRadiusMultiple", float)`
  - `OperationIndex = NewType("OperationIndex", int)`
  - `seconds_value(value, *, name) -> Seconds`
  - `degrees_value(value, *, name) -> Degrees`
  - `closed_unit_fraction(value, *, name) -> UnitFraction`
  - `motion_count(value, *, name) -> MotionCount`
  - `tool_radius_multiple(value, *, name) -> ToolRadiusMultiple`
  - `operation_index(value, *, operation_count) -> OperationIndex`
  - the eight named spec failures in `benchmarks.errors`

- [x] **Step 1: Write RED runtime tests for each observation-unit validator**

Create tests that accept finite values, reject booleans, NaN, infinities,
negative seconds/counts, and fractions outside `[0, 1]`:

```python
def test_closed_unit_fraction_rejects_value_above_one() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        closed_unit_fraction(1.01, name="uncut fraction")


def test_motion_count_rejects_bool() -> None:
    with pytest.raises(InvalidHeldPathEvidenceError):
        motion_count(True, name="gouging motions")
```

- [x] **Step 2: Write the strict type contract**

The contract must prove the six scalar domains remain distinct and geometry
continues to use the existing unit/frame vocabulary:

```python
seconds = assert_type(seconds_value(1.0, name="audit"), Seconds)
degrees = assert_type(degrees_value(80.0, name="cap"), Degrees)
index = assert_type(operation_index(3, operation_count=4), OperationIndex)
point = assert_type(Point2[WorldXY].build(1.0, 2.0), Point2[WorldXY])
assert_type(Point3[WorldXYZ].build(1.0, 2.0, 0.0), Point3[WorldXYZ])
```

- [x] **Step 3: Run RED tests**

Run:

```bash
pixi run pytest -- tests/benchmarks/test_units.py -n auto --testmon -q
pixi run mypy --strict --warn-unused-ignores tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: collection/type checking fails because `benchmarks.units` and the new
failures do not exist.

- [x] **Step 4: Implement the minimal observation-unit module and failures**

Use Python 3.9-compatible `NewType` declarations and small functions returning
the typed values. Append these independent errors to `benchmarks/errors.py`:

```python
class InvalidHeldOperationSnapshotError(BenchmarkError):
    """A snapshotted operation contains malformed geometric or motion data."""


class InvalidHeldPathEvidenceError(BenchmarkError):
    """Typed values or operation coverage in Held path evidence are invalid."""


class InvalidHeldPathReportContextError(BenchmarkError):
    """Held report invocation metadata is empty, non-UTC, or malformed."""


class ContradictoryEngagementEvidenceError(BenchmarkError):
    """Guarded replay and sampled exact-predicate evidence contradict."""


class ContradictoryPathQualityEvidenceError(BenchmarkError):
    """Path-quality attribution does not reduce to its aggregate value."""


class MutatedHeldToolpathError(BenchmarkError):
    """A path-replay consumer changed the generated operation stream."""


class HeldPathNotEligibleForPostQualificationError(BenchmarkError):
    """A complete characterization has at least one open Phase 1 criterion."""


class UnexpectedHeldPathCaseError(BenchmarkError):
    """Characterization received a Held case other than Figure 5."""
```

- [x] **Step 5: Extend the `types-benchmarks` task and run GREEN gates**

Add `benchmarks/units.py` and the new type-contract file to the explicit mypy
file list, then run:

```bash
pixi run pytest -- tests/benchmarks/test_units.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/units.py benchmarks/errors.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/units.py benchmarks/errors.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: all named tests and strict typing pass.

- [x] **Step 6: Commit**

```bash
git add benchmarks/units.py benchmarks/errors.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): add path evidence units'
```

---

### Task 2: Behavior-complete immutable operation snapshots

**Files:**

- Create: `benchmarks/held_path_snapshot.py`
- Create: `tests/benchmarks/test_held_path_snapshot.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: `ToolpathResult`, `ToolpathOperation`, canonical frame/unit types,
  `OperationIndex`, and `InvalidHeldOperationSnapshotError`.
- Produces:
  - frozen/init-disabled `HeldLineSnapshot`, `HeldArcSnapshot`, and
    `HeldCircleSnapshot`
  - `HeldOperationSnapshot = HeldLineSnapshot | HeldArcSnapshot |
    HeldCircleSnapshot`
  - `snapshot_toolpath(result: ToolpathResult) -> tuple[HeldOperationSnapshot, ...]`
  - `assert_toolpath_matches_snapshot(result, snapshot) -> None`

- [x] **Step 1: Write RED factory tests for the closed primitive union**

Construct one line, arc, and circle `ToolpathOperation`. Assert each snapshot
retains ordinal, operation role, `path_index`, clockwise travel, optional
tangents, 3D geometry, curve frame axes, radius, and arc angles. Include:

```python
def test_circle_snapshot_retains_frame_phase() -> None:
    circle = Circle(2.0, frame=Frame([3.0, 4.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]))
    snapshot = snapshot_toolpath(_result([_operation(circle, clockwise=True)]))[0]
    assert isinstance(snapshot, HeldCircleSnapshot)
    assert snapshot.centre == Point3[WorldXYZ].build(3.0, 4.0, 0.0)
    assert snapshot.xaxis == Direction3[WorldXYZ].build(0.0, 1.0, 0.0)
```

- [x] **Step 2: Write RED mutation and malformed-input tests**

Parameterize one-field changes covering operation order/count, role,
`path_index`, clockwise, endpoints, centre, Z, frame axes, radius, arc angles,
and tangents. Every mutation must raise `MutatedHeldToolpathError`. Factories
must reject non-finite geometry, non-positive radius, non-unit or non-orthogonal
axes, invalid angles, unsupported geometry, and malformed tangents through
`InvalidHeldOperationSnapshotError`. Direct construction must raise `TypeError`.

- [x] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py -n auto --testmon -q
```

Expected: import fails because the snapshot module does not exist.

- [x] **Step 4: Implement structural snapshot factories**

Use three records rather than optional primitive fields. Build canonical
`Point3[WorldXYZ]`, `Direction3[WorldXYZ]`, `Millimetre`, and `Radian` values;
retain no COMPAS geometry or NumPy array. Structural comparison is exactly:

```python
def assert_toolpath_matches_snapshot(
    result: ToolpathResult,
    snapshot: tuple[HeldOperationSnapshot, ...],
) -> None:
    observed = snapshot_toolpath(result)
    if observed != snapshot:
        raise MutatedHeldToolpathError("Generated toolpath differs from the characterized operation snapshot.")
```

Do not import or extend the existing engagement-audit identity subsystem.

- [x] **Step 5: Run GREEN runtime, type, and formatting gates**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_snapshot.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_path_snapshot.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

- [x] **Step 6: Commit**

```bash
git add benchmarks/held_path_snapshot.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): snapshot Held operations'
```

---

### Task 3: Canonical survey retains exact sample evidence

**Files:**

- Modify: `benchmarks/survey.py:96-110, 284-345`
- Create: `tests/benchmarks/test_survey.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: `Point2[WorldXY]` and the existing `_stock_2.engagement_at` result.
- Produces: `EngagementSample.position: Point2[WorldXY]` and
  `EngagementSample.cap_exceeded: bool`; `MotionQuality.cap_exceeded` becomes
  `any(sample.cap_exceeded for sample in samples)`.

- [x] **Step 1: Write RED tests for retained sample positions and predicates**

Tests must establish that a circle has 45 seam-unique samples, while open lines
and arcs have 46 endpoint-inclusive samples at the standard 45 intervals. Assert
the first/last typed positions and:

```python
assert motion.cap_exceeded is any(sample.cap_exceeded for sample in motion.samples)
```

Monkeypatch only the reporting angle in a focused unit test to prove the Boolean
comes directly from the exact predicate result and is never reconstructed by
comparing `engagement_deg` with the cap.

- [x] **Step 2: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_survey.py -n auto --testmon -q
```

Expected: the new `EngagementSample` fields are absent.

- [x] **Step 3: Extend `EngagementSample` and `_measure_motion` additively**

For each `(distance, x, y)` station, retain the exact returned `exceeded` value
and `Point2[WorldXY].build(x, y)`. Update the docstring so
`engagement_deg` is explicitly reporting-only, per-position predicates are exact,
and absence across the finite station set remains sampled-negative evidence.
Do not change `QUALITY_SAMPLES_PER_MOTION` or `_motion_samples`.

- [x] **Step 4: Prove existing survey/quality behavior is unchanged**

```bash
pixi run pytest -- tests/benchmarks/test_survey.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/survey.py tests/benchmarks/test_survey.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/survey.py tests/benchmarks/test_survey.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: selected tests pass; no existing metric changes.

- [x] **Step 5: Commit**

```bash
git add benchmarks/survey.py tests/benchmarks/test_survey.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): retain sampled path evidence'
```

---

### Task 4: Add canonical typed quality observations beside the current judge

**Files:**

- Create: `benchmarks/quality_observations.py`
- Create: `tests/benchmarks/test_quality_observations.py`
- Create: `tests/benchmarks/typecheck/quality_observations_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: `PocketSpec`, `tuple[HeldOperationSnapshot, ...]`, `PathSurvey`, and
  `CoverageEstimate`.
- Produces:
  - closed `CriterionName`, `EvidenceKind`, and `ObservationOutcome` literals
  - init-disabled `FractionCriterion`, `CountCriterion`, `DegreesCriterion`, and
    `ToolRadiusMultipleCriterion`
  - `OperationPair`, `MeasuredStep`, `PathQualityAttribution`, and
    `PathQualityAssessment`
  - `assess_path_quality(spec, snapshot, survey, coverage) -> PathQualityAssessment`

- [x] **Step 1: Write RED criterion-factory tests**

Pin all twelve criterion names, thresholds, evidence kinds, and outcome
vocabulary. Required values are exactly `0`, the case cap, and `2` tool radii;
never infer them from current measurements. A wrong evidence kind, unit domain,
non-finite value, negative count, invalid index, or inconsistent outcome must
raise `InvalidHeldPathEvidenceError`.

- [x] **Step 2: Write RED simple-attribution tests**

Using synthetic `PathSurvey` values, require:

- uncut stock remains a spatial aggregate only;
- gouging, unsafe rapid, zero-length, degenerate-loop, redundant, cap-exceeded,
  and slotting counts equal their exact unique operation-index tuples;
- cap-exceeded indices come from `sample.cap_exceeded`, not reporting degrees;
- count observations reject attribution cardinality disagreement.

- [x] **Step 3: Write RED junction and step-attribution tests**

Pin the existing adjacency semantics:

- continuity compares only consecutive source operation indices and uses
  `CONTINUITY_TOOL_RADIUS_FRACTION * tool_radius`;
- reversal is also counted and attributed as a tangent break, while reversal
  excludes curvature and a non-reversal tangent break excludes curvature;
- engagement steps compare peak sampled engagement only for adjacent cut motions;
- loop-radius runs split at both rapid indices and generator `path_index` changes;
- a maximum stores the winning adjacent operation pair, with deterministic first
  occurrence on ties;
- failure-pair lists contain every pair over the unchanged threshold, not merely
  one pair reproducing the maximum.

- [x] **Step 4: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_quality_observations.py -n auto --testmon -q
pixi run mypy --strict --warn-unused-ignores tests/benchmarks/typecheck/quality_observations_contract.py
```

Expected: import fails because the canonical observation module is absent.

- [x] **Step 5: Implement the canonical findings and typed assessment**

Move no existing function yet. Reproduce the current deciding expressions in
new detailed reducers and construct the assessment from their outputs. The
evidence mapping is fixed:

```python
EVIDENCE_BY_CRITERION: dict[CriterionName, EvidenceKind] = {
    "uncut fraction": "sampled_diagnostic",
    "gouging motions": "sampled_diagnostic",
    "unsafe rapids": "tolerance_diagnostic",
    "continuity breaks": "tolerance_diagnostic",
    "zero-length motions": "derived_geometry",
    "degenerate loops": "derived_geometry",
    "redundant operations": "exact_depletion_replay",
    "cap exceedances": "sampled_exact_predicate",
    "slotting motions": "sampled_exact_predicate",
    "max engagement step (deg)": "sampled_diagnostic",
    "max loop radius step (tool radii)": "derived_geometry",
    "tangent breaks": "tolerance_diagnostic",
}
```

Use named existing tolerance constants; introduce no numeric literal at a
decision call site.

- [x] **Step 6: Run GREEN focused gates**

```bash
pixi run pytest -- tests/benchmarks/test_quality_observations.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/typecheck/quality_observations_contract.py
pixi run ruff check benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/typecheck/quality_observations_contract.py
```

- [x] **Step 7: Commit additive observations**

```bash
git add benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/typecheck/quality_observations_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): add attributed quality observations'
```

---

### Task 5: Prove old/new quality parity and stop for approval

**Files:**

- Modify: `tests/benchmarks/test_quality.py:82-137, 685-842`
- Modify: `tests/benchmarks/test_quality_invariants.py`
- Create: `docs/superpowers/state/held-figure5-quality-parity.md`

**Interfaces:**

- Consumes: untouched `measure_quality`, untouched private reducers,
  `assess_path_quality`, the same generated `ToolpathResult`, survey, coverage,
  and operation snapshot.
- Produces: reviewable proof that the new attributed observations preserve every
  old aggregate, threshold, test ID, and known-red cell.

- [x] **Step 1: Add parity assertions without changing the old verdict path**

Keep `_violations` and the final gate assertion untouched. Add a helper that
compares every gated `PathQuality` field with its typed observation, then add it
to synthetic machinery and invariant cases. Counts must reduce from every
attributed operation/pair; maximum observations must name the exact winning pair
or `None` when the maximum is zero.

- [x] **Step 2: Add parity to all twelve `QUALITY_GATE_CASES`**

In `test_the_generated_path_is_worth_running`, generate once and run the
unchanged quality path. Test spies capture its single `PathSurvey` and
`CoverageEstimate`; feed those exact objects plus the operation snapshot to the
additive assessment. Assert parity, then leave the existing
`_violations`/`assert not violations` lines unchanged. Per cell, prove exactly
one survey and one coverage measurement. The twelve tests must remain red for
the same reasons after the parity assertion passes.

- [x] **Step 3: Run synthetic and invariant parity GREEN**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
```

Expected: all selected tests pass.

- [x] **Step 4: make the exact known-red oracle deterministic**

```bash
zsh -c 'pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q; red_status=$?; if [[ $red_status -ne 1 ]]; then exit 2; fi; pixi run red-manifest'
```

Observed on 2026-09-02: exactly twelve
`test_the_generated_path_is_worth_running` failures retained their IDs and
criterion messages, but the full suite reported sixteen rather than seventeen
reds because the generated translation property did not discover its known
counterexample. A focused invocation on the unchanged code later found two
counterexamples, including `cap_exceedances` changing from `1` to `2` under a
pure translation. The defect remains; its discovery is nondeterministic.

Task 5A.1 must pin a minimized example before this step can be checked. The
required result remains exactly seventeen declared reds overall: twelve quality,
one deterministic benchmark translation-invariance cell, and four adaptive.
Exit status `1` is required for the focused expected-red invocation; collection,
configuration, interruption, and infrastructure statuses are blockers, not reds.

- [x] **Step 5: Record parity evidence**

Write `held-figure5-quality-parity.md` with:

- command and date;
- every old/new field compared;
- attribution cardinality and extremum-pair results;
- exact twelve quality-red membership and seventeen-red repository membership;
- explicit statement that current production still uses the old reducers; and
- the removal decision requested from Jelle.

- [x] **Step 6: Run focused hygiene and commit**

```bash
pixi run types-benchmarks
pixi run ruff format tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py
pixi run ruff check benchmarks tests/benchmarks
git diff --check
git add tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py docs/superpowers/state/held-figure5-quality-parity.md
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'test(benchmarks): prove quality observation parity'
```

- [x] **Step 7: STOP for independent review before removal approval**

The committed parity evidence was reviewed by six independent lenses. The
review converged on the Task 5A findings below. Do not request removal approval,
start Task 6, remove an old reducer, or redirect `measure_quality` until every
Task 5A checkpoint passes.

---

### Task 5A.1: Make declared-red membership deterministic

**Files:**

- Modify: `tests/benchmarks/test_quality_invariants.py:502-518`
- Modify: `docs/superpowers/state/held-figure5-quality-parity.md`

**Interfaces:**

- Consumes: the existing translation-invariance property and its protected
  quality judge.
- Produces: one pinned dyadic `@example` under the existing test ID; generated
  examples remain additional discovery evidence.

- [x] **Step 1: Pin the minimized translation witness**

Retain the generated property and insert this minimized example between
`@PROPERTY_SETTINGS` and the existing `@given` decorator:

```python
@PROPERTY_SETTINGS
@example(
    stations=[
        (0, 0.0, 0.0, 2.125),
        (0, 0.875, 0.0, 1.0),
        (1, 0.0, 0.0, 0.375),
        (1, 0.875, 0.0, 0.375),
    ],
    shift=(0.0, 2.0),
)
```

Do not mark the property as expected failure. The red manifest continues to own
the unchanged node ID and the defect remains visible as a normal assertion
failure.

- [x] **Step 2: Prove the red is independent of generated examples**

Run the existing test with fixed seeds that previously missed the defect:

```bash
pixi run pytest -- tests/benchmarks/test_quality_invariants.py::test_moving_the_pocket_across_the_table_changes_no_metric -n auto --hypothesis-seed=1 -q
pixi run pytest -- tests/benchmarks/test_quality_invariants.py::test_moving_the_pocket_across_the_table_changes_no_metric -n auto --hypothesis-seed=2 -q
```

Expected: both exit `1` on the pinned translation witness with a machining
metric changing under translation. Collection, configuration, interruption,
or process termination is not an expected red.

- [x] **Step 3: Record the corrected oracle interpretation and commit**

Update the parity state: the 16-red full run was a probabilistic miss, the
translation defect remains reproducible, and exact manifest membership now has
a deterministic witness.

```bash
pixi run ruff format tests/benchmarks/test_quality_invariants.py
pixi run ruff check tests/benchmarks/test_quality_invariants.py
git diff --check
git add tests/benchmarks/test_quality_invariants.py docs/superpowers/state/held-figure5-quality-parity.md
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'test(quality): pin translation defect'
```

---

### Task 5A.2: Fail closed at the 2D motion boundary

**Files:**

- Modify: `benchmarks/held_path_snapshot.py:97-141, 199-253, 312-387`
- Modify: `benchmarks/survey.py:237-289, 535-590`
- Modify: `tests/benchmarks/test_held_path_snapshot.py`
- Modify: `tests/benchmarks/test_survey.py`

**Interfaces:**

- Consumes: `HeldOperationSnapshot`, COMPAS line/arc/circle geometry, the
  inferred cut plane, and the shared `TOL` predicates.
- Produces: snapshots whose curve ranges and retained tangents are internally
  consistent with the source motion; the cut-plane survey refuses silent XY
  projection of tilted or vertically displaced curves. The structural snapshot
  may still record unsupported 3D input so the refusal remains observable.

- [x] **Step 1: Write RED malformed-motion tests**

Add direct-factory and real-operation tests named
`test_snapshot_rejects_descending_arc_range`,
`test_snapshot_retains_zero_sweep_for_zero_length_diagnosis`,
`test_snapshot_rejects_line_tangent_that_disagrees_with_geometry`,
`test_snapshot_rejects_arc_tangent_that_disagrees_with_travel`,
`test_snapshot_rejects_circle_tangent_that_disagrees_with_clockwise`,
`test_survey_rejects_tilted_cut_circle`, and
`test_survey_rejects_curve_outside_the_inferred_cut_plane`. Use the existing
`_operation`, `_result`, and `_survey_motion` fixtures; assert the named
snapshot or replay exception in every refusal test and structural equality in
the zero-sweep retention test.

Descending arcs are malformed because COMPAS reports negative length. A zero
sweep remains admissible so the protected zero-length quality criterion can
diagnose it. In-plane rotation and translation remain valid.

- [x] **Step 2: Run the focused RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_survey.py -k 'descending or zero_sweep or tangent_that_disagrees or tilted or inferred_cut_plane' -n auto --testmon --testmon-noselect -q
```

Expected: descending ranges, contradictory unit tangents, and non-planar curves
are accepted when they must fail.

- [x] **Step 3: Implement the narrow geometric validation**

Keep geometry in typed world coordinates. Enforce `end_angle >= start_angle`;
do not reject equality. Derive each primitive's expected travel tangent from its
geometry and clockwise flag, then compare it with retained tangents through the
shared COMPAS angular tolerance. At the cut-plane replay boundary require arc
and circle axes to lie in world XY and their centres to lie on the inferred cut
height. Raise `InvalidHeldOperationSnapshotError` for malformed snapshots and
the existing `UnreplayableOperationError` when survey input leaves the cut-plane
model. Introduce no numeric tolerance literal.

- [x] **Step 4: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_survey.py -n auto --testmon --testmon-noselect -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_snapshot.py benchmarks/survey.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_survey.py
pixi run ruff check benchmarks/held_path_snapshot.py benchmarks/survey.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_survey.py
git diff --check
git add benchmarks/held_path_snapshot.py benchmarks/survey.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_survey.py
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'fix(benchmarks): validate 2d path snapshots'
```

---

### Task 5A.3: Make quality evidence factories bypass-safe

**Files:**

- Create: `src/compas_cgal/replay_classification.py`
- Create: `tests/test_replay_classification.py`
- Modify: `src/compas_cgal/engagement.py`
- Modify: `benchmarks/depletion.py`
- Modify: `benchmarks/survey.py:129-229, 237-289, 550-588`
- Modify: `benchmarks/quality_observations.py:300-675`
- Modify: `tests/benchmarks/test_quality_observations.py`
- Modify: `tests/benchmarks/test_quality.py:745-880`
- Modify: `tests/benchmarks/test_quality_invariants.py:371-480`
- Modify: `tests/benchmarks/typecheck/quality_observations_contract.py`

**Interfaces:**

- Consumes: one `PocketSpec`, complete operation snapshot, survey produced from
  that same specification and operation stream, and unchanged criterion
  thresholds.
- Produces:
  - one dependency-neutral, four-way replay classification shared by production
    engagement, depletion, survey, and immutable evidence validation;
  - preserved consumer exceptions with the neutral classification failure as
    their retained cause;
  - operation-level plunge/retract observations sufficient to prove a complete
    operation partition;
  - a retained closed unit discriminator on `MeasuredStep[StepUnitT]`;
  - unconditional validation of `operation_count`;
  - fail-closed survey/spec/snapshot binding;
  - criterion outcomes and failure-pair attribution that cannot contradict.

- [x] **Step 1: Write RED input-binding tests**

Require `InvalidHeldPathEvidenceError` for a foreign `PocketSpec`, omitted or
reordered operation, mismatched operation role, wrong primitive kind, and
incomplete motion/rapid/plunge/retract partition. The survey must retain the
operation indices currently represented only by plunge/retract counts, and the
validated partition must equal `range(len(snapshot))` exactly.

Task 6's public `reduce_quality_evidence` must compute coverage directly from
the validated survey's `final_stock`; no public production boundary may accept
an independently supplied `CoverageEstimate`. The Task 5 transition reducer may
continue accepting captured coverage only until Task 6 is explicitly approved.

- [x] **Step 2: Write RED unit and count factory tests**

Add runtime and strict-type contracts for:

```python
with pytest.raises(InvalidHeldPathEvidenceError):
    MeasuredStep.build(value=ToolRadiusMultiple(1.0), pair=pair, unit=cast(Any, "seconds"))

with pytest.raises(InvalidHeldPathEvidenceError):
    PathQualityAttribution.build(operation_count=True, **empty_attribution)

with pytest.raises(InvalidHeldPathEvidenceError):
    PathQualityAttribution.build(operation_count=-1, **empty_attribution)
```

Retain `unit: Literal["degrees", "tool_radius_multiple"]` in each constructed
`MeasuredStep`. Use typed overloads to correlate `Degrees` with `"degrees"` and
`ToolRadiusMultiple` with `"tool_radius_multiple"`; validate the retained unit
again when installing a step into its engagement or loop-radius slot. Require an
exact non-Boolean integer operation count; zero remains valid only for a wholly
empty standalone attribution and is ineligible for Figure 5.

- [x] **Step 3: Write RED contradiction and exact-source tests**

Mutation tests must reject:

- a violated maximum with no failure pair;
- a satisfied maximum with a non-empty failure list;
- an attributed maximum absent from its failure list;
- an in-bounds but non-adjacent engagement pair;
- an in-bounds pair crossing a rapid or path-chain boundary; and
- replacement of any one of the nine count criteria's source tuples with a
  same-cardinality wrong tuple.

Raise `ContradictoryPathQualityEvidenceError` when valid component records do not
form a coherent assessment. Independently reconstruct exact expected source
tuples from raw survey observations in every synthetic, invariant, and twelve
gate invocation; aggregate equality and cardinality alone are insufficient.

- [x] **Step 4: Run RED, implement the smallest validators, and run GREEN**

The replay-classification correction must additionally cover strict tolerance
boundaries, retract precedence, invalid ramps, planar-curve height, and all four
real consumers. The neutral classifier depends only on toolpath roles, existing
millimetre types, COMPAS tolerance, and geometry facts; consumer-specific stock
effects and exception translation remain at their existing boundaries.

```bash
pixi run pytest -- tests/benchmarks/test_quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'binding or partition or operation_count or unit or contradiction or exact_sources' -n auto --testmon --testmon-noselect -q
```

Expected before implementation: every newly added negative contract reaches the
factory and is accepted when it must fail. Implement only the validation needed
by those contracts, then rerun the same command expecting all selected tests to
pass.

- [x] **Step 5: Run focused parity, types, Ruff, and commit**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_quality_observations.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon --testmon-noselect -q
pixi run types-benchmarks
pixi run pytest -- tests/test_replay_classification.py tests/test_engagement_audit.py tests/benchmarks/test_depletion.py tests/benchmarks/test_survey.py tests/benchmarks/test_quality_observations.py -n auto --testmon --testmon-noselect -q
pixi run ruff format benchmarks/survey.py benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/typecheck/quality_observations_contract.py
pixi run ruff check benchmarks/survey.py benchmarks/quality_observations.py tests/benchmarks
git diff --check
git add benchmarks/survey.py benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/typecheck/quality_observations_contract.py
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'fix(benchmarks): bind quality evidence'
```

---

### Task 5A.4: Establish a semantic transition gate and truthful checkpoint

**Files:**

- Create: `tests/benchmarks/test_quality_transition.py`
- Modify: `tests/benchmarks/test_quality.py`
- Modify: `tools/red_manifest.py`
- Modify: `tests/tools/test_red_manifest.py`
- Modify: `docs/red_manifest.json`
- Create: `docs/red_manifest-adaptive.json`
- Modify: `docs/superpowers/state/held-figure5-quality-parity.md`
- Modify: `docs/superpowers/plans/2026-09-02-held-figure5-2d-path-characterization.md`

**Interfaces:**

- Consumes: all twelve unchanged `QUALITY_GATE_CASES`, the protected old judge,
  and the validated attributed assessment.
- Produces: a green transition oracle that authenticates each structured
  criterion vector before the deliberate product-gate assertion; exact known-red
  membership remains a separate repository-state gate.

- [x] **Step 1: Write the permanent twelve-case green transition oracle**

Parameterize the unchanged gate cases without executing the final deliberate
`assert not violations`. For every case, generate once and require:

- one survey and one coverage evaluation;
- exact old/new values and thresholds for all twelve criteria;
- exact structured `ObservationOutcome` values;
- exact operation and pair attribution, including first-tie maxima; and
- the expected ordered criterion-message vector later consumed by the red test.

This test must fail if an arbitrary earlier assertion replaces the intended
product-gate failure, even when the red test ID and count remain unchanged.
Extract one test-owned pre-verdict evaluator and call that same evaluator from
both the green transition test and the protected red wrapper. Validate the red
wrapper behavior from its JUnit failure messages; do not freeze source shape or
use an AST/grep implementation-structure test.

- [x] **Step 2: Prove the transition and known-red oracles separately**

```bash
pixi run pytest -- tests/benchmarks/test_quality_transition.py -n auto --testmon --testmon-noselect -q
zsh -c 'pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q --junitxml=build/task5a4-quality.xml; red_status=$?; if [[ $red_status -ne 1 ]]; then exit 2; fi'
pixi run python - <<'PY'
from pathlib import Path

from benchmarks.gate import GATE_CAP_CASES, GATE_GENERATOR_NAMES, GATE_POCKET_NAMES
from tests.benchmarks.test_quality import _expected_quality_gate_violations
from tools.red_manifest import parse_junit_failures

expected = {
    f"tests.benchmarks.test_quality::test_the_generated_path_is_worth_running[{cap_id}-{generator}-{pocket}]": _expected_quality_gate_violations(cap, generator, pocket)
    for cap_id, cap in GATE_CAP_CASES
    for generator in GATE_GENERATOR_NAMES
    for pocket in GATE_POCKET_NAMES
}
failures = {failure.identity: failure.message for failure in parse_junit_failures(Path("build/task5a4-quality.xml"))}
assert failures.keys() == expected.keys()
for identity, violations in expected.items():
    message = failures[identity]
    marker = f"E       NOT MACHINABLE -- {len(violations)} criteria failed:\n"
    assert message.count(marker) == 1
    reported = tuple(line.removeprefix("E         ") for line in message.split(marker, 1)[1].splitlines() if line.startswith("E         ["))
    assert reported == violations
    assert "assert not violations" in message
print("12 exact quality reds reached the reviewed final violation vectors")
PY
pixi run red-manifest
```

Expected: the transition module is green; all twelve quality cells fail only at
the final product assertion; the repository manifest reports the deterministic
seventeen-red membership.
The quality JUnit must contain the exact twelve IDs and their expected final
criterion-message vectors; errors, duplicates, extra/missing failures, or an
earlier assertion block the task. Tighten the adaptive manifest entry to the
four exact node IDs rather than a module-wide count regex.

- [x] **Step 3: Resolve the prior native-process instability**

Run the previously implicated adaptive files three times under the acceptance
xdist configuration with Python fault handling enabled:

```bash
for run in 1 2 3; do
  report="build/task5a4-adaptive-${run}.xml"
  PYTHONFAULTHANDLER=1 pixi run pytest -- tests/adaptive/test_generator.py tests/adaptive/test_route_retrace_generator.py -n auto --dist=loadgroup -q --junitxml="$report"
  pytest_status=$?
  if [[ $pytest_status -ne 1 ]]; then exit 2; fi
  pixi run python -m tools.red_manifest "$report" --manifest docs/red_manifest-adaptive.json
done
```

Expected: each invocation completes with only the four manifest-owned adaptive
assertion failures. A segfault, truncated JUnit report, collection error, or
other process termination blocks Task 6 and starts a focused native-lifetime
diagnosis; it may not be classified as infrastructure without a concrete
external failure signature.
The loop must continue after expected pytest exit status `1`; each JUnit report
must prove the exact same four-member failure set independently.

- [x] **Step 4: Reconcile plan and durable evidence**

Correct the parity page so it attributes the Task 4 first-tie contracts to the
completed full run rather than the narrower focused parity command. Record the
deterministic translation witness, exact-source mutation coverage, semantic
transition result, native stability result, and unchanged production routing.

Mark Task 5 Step 4 and Tasks 5A.1-5A.4 complete only after their commands produce
the stated evidence. Change the plan header to `awaiting Task 6 removal
approval`; do not mark Task 6 active.

- [x] **Step 5: Run hygiene and commit the checkpoint**

```bash
pixi run types-benchmarks
pixi run ruff format tests/benchmarks/test_quality.py tests/benchmarks/test_quality_transition.py tools/red_manifest.py tests/tools/test_red_manifest.py
pixi run ruff check benchmarks tests/benchmarks tools/red_manifest.py tests/tools/test_red_manifest.py
git diff --check
git add tests/benchmarks/test_quality.py tests/benchmarks/test_quality_transition.py tools/red_manifest.py tests/tools/test_red_manifest.py docs/red_manifest.json docs/red_manifest-adaptive.json docs/superpowers/state/held-figure5-quality-parity.md docs/superpowers/plans/2026-09-02-held-figure5-2d-path-characterization.md
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'test(benchmarks): harden quality transition'
```

- [ ] **Step 6: STOP and request explicit Task 6 approval**

Present the four Task 5A commits, deterministic seventeen-red reconciliation,
semantic twelve-case transition result, and native stability evidence. Task 6
remains a separate user decision because it redirects production and removes
superseded deciding bodies.

---

### Task 6: Converge `measure_quality` to the canonical reducers

**Preconditions:** Every Task 5A checkbox is complete; Task 5 Step 4 reports the
deterministic seventeen-red set; the permanent twelve-case semantic transition
oracle is green; the implicated native tests complete stably under xdist; and
Jelle explicitly approved Task 6 after reviewing that evidence. Approval alone
cannot override a failed technical precondition.

**Files:**

- Modify: `benchmarks/quality.py:475-616, 1074-1285`
- Modify: `benchmarks/quality_observations.py`
- Modify: `tests/benchmarks/test_quality.py`
- Modify: `tests/benchmarks/test_quality_invariants.py`
- Modify: `tests/benchmarks/typecheck/quality_observations_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: canonical snapshot, survey, coverage, and attributed assessment.
- Produces:
  - `QualityEvidence` containing `PathQuality`, `PathQualityAssessment`, and
    `CoverageEstimate`
  - `reduce_quality_evidence(spec, snapshot, survey, *,
    grid: int = COVERAGE_GRID_SAMPLES) -> QualityEvidence`
  - unchanged `measure_quality(spec, result, *, samples_per_motion: int =
    QUALITY_SAMPLES_PER_MOTION, grid: int = COVERAGE_GRID_SAMPLES) -> PathQuality`
    compatibility API delegating to the single survey/reducer path

- [x] **Step 1: Write RED additive-reducer and compatibility tests**

Instrument `survey_path`, `measure_coverage`, and `assess_path_quality`. Require
one call each, coverage computed from the validated survey rather than accepted
as an argument, identical public `PathQuality` output, and the existing named
exceptions for zero-length paths and invalid grids. Keep separate tests for the
new additive reducer and the still-unmodified `measure_quality` compatibility
route.

- [x] **Step 2: Run RED focused tests**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'single_survey or compatibility' -n auto --testmon -q
```

Expected: the new reducer API does not exist and the compatibility path does not
delegate.

- [x] **Step 3: Add `QualityEvidence` beside the current production path**

`reduce_quality_evidence` measures coverage exactly once, evaluates the
canonical assessment once, and builds all five existing `PathQuality` groups.
It first applies the Task 5A survey/spec/snapshot validator and computes coverage
internally from `survey.final_stock`; no independently supplied coverage record
enters the public production boundary. For the twelve gated fields, take values
only from the assessment. Continue to compute report-only fields through their
existing reducers. Do not redirect `measure_quality` and do not delete an old
body in this step.

- [x] **Step 4: Validate and commit the additive path**

Run the new reducer tests directly while the old production path remains active:

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_quality_transition.py -k 'quality_evidence or transition' -n auto --testmon --testmon-noselect -q
pixi run types-benchmarks
pixi run ruff format benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_quality_transition.py
pixi run ruff check benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks
git diff --check
git add benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_quality_transition.py tests/benchmarks/typecheck/quality_observations_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): add canonical quality evidence'
```

Expected: the additive path reproduces every structured criterion and source;
`measure_quality` still uses the old protected path.

- [x] **Step 5: Redirect production with all old bodies retained**

Make `measure_quality` delegate to `reduce_quality_evidence`. Do not remove the
old reducers or the temporary transition checks. Run the compatibility tests,
permanent twelve-case semantic transition oracle, non-gate quality/invariant
set, all twelve gate cells, `red-manifest`, strict types, and Ruff. The routed
path must preserve both the structured criterion vectors and the deterministic
seventeen-red repository set before deletion begins.

- [x] **Step 6: Remove only superseded decision bodies**

Remove the old continuity, junction, engagement-step, and loop-radius-step
decision functions after all callers use canonical findings. Remove duplicate
gated count expressions from `_elementary`, `_cut`, and `_speed`. Retain
report-only curvature/reversal data by deriving it from canonical junction
findings, and retain all non-gated metric functions.

Change the test-local `_violations` helper to accept only
`PathQualityAssessment` and project its failed criteria. In each gate cell, build
one `QualityEvidence` from one snapshot, survey, and coverage evaluation; pass
`evidence.path_quality` to reporting and `evidence.assessment` to `_violations`.
Compatibility tests separately prove
`measure_quality(spec, result) == evidence.path_quality`. Preserve the
parametrization, test IDs, literal threshold contract tests, and final
`assert not violations`; after the approved convergence there must be one survey,
one coverage evaluation, and one threshold-decision call graph.

- [x] **Step 7: Retire temporary spies and retain permanent transition contracts**

Remove the temporary parity spies from Task 5 only after
`test_quality_transition.py` exercises the routed production path. Keep that
permanent green module, literal threshold assertions, criterion
name/value/evidence tests, unchanged `QUALITY_GATE_CASES`, unchanged test IDs,
and unchanged final gate assertion.

- [x] **Step 8: Rerun the same quality and type gates after deletion**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_qualityfigures.py tests/benchmarks/test_quality_transition.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
zsh -c 'pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q; red_status=$?; if [[ $red_status -ne 1 ]]; then exit 2; fi; pixi run red-manifest'
pixi run types-benchmarks
pixi run ruff format benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py
pixi run ruff check benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks
```

Expected: non-gate/non-translation tests pass; exactly the same twelve quality
cells remain red; the repository red manifest still reports exactly seventeen
declared reds.

- [x] **Step 9: Commit the validated authority transition**

```bash
git add benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/typecheck/quality_observations_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'refactor(benchmarks): unify quality decisions'
```

---

### Task 7: Validated Figure 5 engagement and path evidence

**Files:**

- Create: `benchmarks/held_path_evidence.py`
- Create: `tests/benchmarks/test_held_path_evidence.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: fixed `HeldReferenceCase`, operation snapshot, `EngagementReport`,
  `PathSurvey`, `QualityEvidence`, and four typed phase durations.
- Produces: init-disabled `EngagementExceedanceWitness`,
  `EngagementDispositionCounts`, and:

```python
HeldFigure5Characterization.build(
    case: HeldReferenceCase,
    snapshot: tuple[HeldOperationSnapshot, ...],
    audit: EngagementReport,
    survey: PathSurvey,
    quality: QualityEvidence,
    generation_seconds: Seconds,
    audit_seconds: Seconds,
    survey_seconds: Seconds,
    reduction_seconds: Seconds,
) -> HeldFigure5Characterization
```

- [x] **Step 1: Write RED engagement-partition tests**

Cover these falsifiers independently:

- `stations == 0` with `cap_certified=True` remains TEA-audit-excluded;
- five exact witness rows on one operation produce one demonstrated operation;
- certified/witnessed overlap raises `ContradictoryEngagementEvidenceError`;
- missing, duplicate, reordered, or wrong-kind audit rows fail;
- audit `stations > 0` indices differ from survey-motion indices;
- a material-removing plunge remains TEA-audit-excluded;
- a repeated cleared cut remains TEA-audited while sampled material contact is
  zero; and
- no `max_tea` or `engagement_deg` reporting comparison changes disposition.

- [x] **Step 2: Write RED cross-consumer and factory tests**

Require exact Figure 5 case, complete source-operation partition, in-bounds
witness indices, finite typed world-XY witness coordinates, non-negative finite
timings, all twelve criteria, and `PathQuality` values identical to the assessment
values. Explicitly permit witness coordinates outside the pocket and legal
cutter-centre domain because those positions are valid gouge evidence. Each
inconsistency must raise its named error. Direct construction must raise
`TypeError`.

Require `tea_audited | excluded` to equal every source-operation index with an
empty intersection, `len(excluded) == len(survey.rapids) + survey.plunges`, and
`audit.cap_violations == len(demonstrated_exceeded | unresolved)`. Derive and
store sampled material-contact count only from `MotionQuality.is_engaged`.

Mutate `audit.operations`, survey stock, and source arrays after construction;
the characterization must not change. It retains only immutable derived
operation records, index tuples, counts, typed timings, `PathQuality`, assessment,
coverage scalars, and witnesses—not raw `EngagementReport`, `PathSurvey`, `Stock`,
`ToolpathResult`, COMPAS geometry, NumPy arrays, or mutable collections.

- [x] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_evidence.py -n auto --testmon -q
```

- [x] **Step 4: Implement construction from authoritative index sets**

The factory derives, never accepts, disposition counts:

```python
tea_audited = frozenset(row.op_index for row in audit.operations if row.stations > 0)
certified = frozenset(row.op_index for row in audit.operations if row.stations > 0 and row.cap_certified)
witnessed = frozenset(witness.operation_index for witness in witnesses)
if certified & witnessed:
    raise ContradictoryEngagementEvidenceError("A certified operation has a sampled exact-predicate exceedance witness.")
unresolved = tea_audited - certified - witnessed
excluded = frozenset(range(len(snapshot))) - tea_audited
```

Witness rows come only from retained survey samples where
`sample.cap_exceeded` is true. Validate index coverage before constructing the
record.

- [x] **Step 5: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_evidence.py tests/test_engagement_audit.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_evidence.py tests/benchmarks/test_held_path_evidence.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_path_evidence.py tests/benchmarks/test_held_path_evidence.py tests/benchmarks/typecheck/held_path_characterization_contract.py
git add benchmarks/held_path_evidence.py tests/benchmarks/test_held_path_evidence.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): characterize Held path evidence'
```

---

### Task 8: Bypass-safe post-qualification candidate

**Files:**

- Create: `benchmarks/held_post_qualification.py`
- Create: `tests/benchmarks/test_held_post_qualification.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: complete `HeldFigure5Characterization`.
- Produces:
  - `post_qualification_failures(characterization) -> tuple[str, ...]`
  - frozen/init-disabled `HeldPostQualificationCandidate.build(...)`
  - `require_post_qualification_candidate(characterization) ->
    HeldPostQualificationCandidate`

- [x] **Step 1: Write RED refusal tests for every open gate**

Starting from complete constructible synthetic characterizations, vary exactly
one post-entry condition: one witnessed exceedance, one unresolved motion, or
each of the twelve open quality criteria. Every case must raise
`HeldPathNotEligibleForPostQualificationError` and return no candidate. The
fully closed fixture must produce a candidate. Wrong-case, mutated-snapshot, and
incomplete-classification inputs stay at the Task 7/9 construction boundaries;
do not bypass the validated characterization type to recreate impossible states
at the candidate boundary.
Assert `post_qualification_failures` returns every open condition in stable order
and an empty tuple only for the closed fixture.

- [x] **Step 2: Write RED authority tests**

Require direct candidate construction to fail. Successful `.build(...)` must
retain the closed characterization, exact immutable snapshot, normalized tool
diameter, 80-degree cap, and Figure 5 case context. The functional helper must
delegate to `.build(...)` and return the candidate, never `None` or a free
snapshot.

- [x] **Step 3: Implement one validating construction path**

Implement one public pure `post_qualification_failures` collector. `.build(...)`
consumes it and owns the only conversion from an empty failure tuple into a
candidate; the functional helper calls the factory. Do not create a second
Boolean gate. Failure messages list every open criterion plus
certified/demonstrated/unresolved counts.

- [x] **Step 4: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_post_qualification.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_post_qualification.py tests/benchmarks/test_held_post_qualification.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_post_qualification.py tests/benchmarks/test_held_post_qualification.py tests/benchmarks/typecheck/held_path_characterization_contract.py
git add benchmarks/held_post_qualification.py tests/benchmarks/test_held_post_qualification.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): gate Held post qualification'
```

---

### Task 9: Generate-once Figure 5 orchestration

**Files:**

- Create: `benchmarks/held_consumer_adapters.py`
- Create: `benchmarks/held_path_characterize.py`
- Create: `tests/benchmarks/test_held_path_characterize.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227`

**Interfaces:**

- Consumes: the four exact Protocol contracts in the spec and the fixed case
  loader.
- Produces:
  - `audit_figure5_engagement(spec, result) -> EngagementReport`
  - `_monotonic_seconds() -> Seconds`
  - `characterize_figure5(generator, engagement_auditor, path_surveyor,
    quality_evidence_reducer, *, phase_observer: Callable[[CharacterizationPhase],
    None], clock: Clock = _monotonic_seconds) -> HeldFigure5Characterization`

The orchestration-local ports are exact:

```python
class Generator(Protocol):
    def __call__(self, spec: PocketSpec) -> ToolpathResult: ...


class EngagementAuditor(Protocol):
    def __call__(self, spec: PocketSpec, result: ToolpathResult) -> EngagementReport: ...


class PathSurveyor(Protocol):
    def __call__(self, spec: PocketSpec, result: ToolpathResult) -> PathSurvey: ...


class QualityEvidenceReducer(Protocol):
    def __call__(
        self,
        spec: PocketSpec,
        snapshot: tuple[HeldOperationSnapshot, ...],
        survey: PathSurvey,
    ) -> QualityEvidence: ...
```

- [x] **Step 1: Write RED adapter tests**

Monkeypatch `audit_toolpath_engagement` and assert the adapter passes exactly
`spec.polygon`, the same result object, `spec.tool_diameter`,
`spec.tea_cap_rad`, and `list(spec.holes)`, then returns the same report. No
exception is caught or rewritten.

- [x] **Step 2: Write RED orchestration order and identity tests**

Use injected spies and a deterministic typed clock. Require this exact sequence:
load Figure 5, generate once, snapshot, assert unchanged, audit, assert unchanged,
survey, assert unchanged, reduce the same survey/snapshot, build evidence. Both
replay consumers receive the object-identical generated result. Record generation,
audit, survey, and coverage/reduction durations independently. The required phase
observer receives `generation`, `guarded_audit`, `survey`, and
`quality_reduction` immediately before their corresponding calls.

- [x] **Step 3: Write RED mutation and exception tests**

Have each replay consumer mutate one operation before returning; the following
post-call assertion must raise `MutatedHeldToolpathError`. Verify unexpected
generator/auditor/survey/reducer exceptions propagate unchanged and no partial
characterization exists.

- [x] **Step 4: Implement the adapter and orchestration**

Keep the adapter as a value-only signature translation. Keep orchestration free
of rendering, file writes, command-line parsing, and post-candidate construction.

- [x] **Step 5: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_characterize.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_consumer_adapters.py benchmarks/held_path_characterize.py tests/benchmarks/test_held_path_characterize.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_consumer_adapters.py benchmarks/held_path_characterize.py tests/benchmarks/test_held_path_characterize.py tests/benchmarks/typecheck/held_path_characterization_contract.py
git add benchmarks/held_consumer_adapters.py benchmarks/held_path_characterize.py tests/benchmarks/test_held_path_characterize.py tests/benchmarks/typecheck/held_path_characterization_contract.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): orchestrate Figure 5 evidence'
```

---

### Task 10: Deterministic evidence report and thin CLI

**Files:**

- Create: `benchmarks/held_path_report.py`
- Create: `tools/held_path_characterization.py`
- Create: `tests/benchmarks/test_held_path_report.py`
- Create: `tests/tools/test_held_path_characterization.py`
- Create: `tests/tools/typecheck/held_path_characterization_contract.py`
- Create: empty package markers under `tests/`, `tests/benchmarks/`,
  `tests/benchmarks/typecheck/`, `tests/tools/`, and `tests/tools/typecheck/`
- Modify: `benchmarks/held_path_evidence.py`
- Modify: `tests/benchmarks/test_held_path_evidence.py`
- Modify: `tests/benchmarks/typecheck/held_path_characterization_contract.py`
- Modify: `pyproject.toml:227-233`

The additive source-count dependency repair is checkpointed separately as
`fix(benchmarks): retain Held source counts`; the renderer consumes those
validated scalar projections and never reloads the case. The empty package
markers only distinguish the two same-named strict type-contract modules.

**Interfaces:**

- Consumes: complete characterization, `post_qualification_failures`, and
  injected aware UTC clock.
- Produces:
  - frozen/init-disabled `HeldPathReportContext.build(generated_at_utc: datetime,
    pixi_command: str, generator_policy_name: str) -> HeldPathReportContext`
  - `render_held_figure5_2d_path(characterization, context) -> str`
  - `DEFAULT_REPORT_PATH = Path("docs/benchmarks/held_figure5_2d_path.md")`
  - `_utc_now() -> datetime`
  - `_print_phase(phase: CharacterizationPhase) -> None`
  - `write_held_figure5_report(*, pixi_command: str, phase_observer:
    Callable[[CharacterizationPhase], None], path: Path = DEFAULT_REPORT_PATH)
    -> HeldFigure5Characterization`
  - `main() -> None`

- [x] **Step 1: Write RED report-context and rendering tests**

Reject naive/non-UTC datetimes, empty command, and empty policy name. With a
fixed context, assert deterministic output begins with the historical/non-release
banner and includes normalized-scale warning, case/tool/cap/primitive/projection
counts, four timings, TEA disposition counts, sampled-contact count, every
witness row, all twelve values/evidence kinds/outcomes, report-only
`PathQuality`, and the post-entry refusal. Pin Markdown escaping and operation
row ordering.

Pin the sentence that the exercised guarded replay is not evidence of native
audit compatibility. The rendered post-entry verdict must come from
`post_qualification_failures`, never copied threshold logic or an exception caught
from candidate construction.

The outcome vocabulary is exact: sampled evidence uses
`no_failure_observed`/`failure_observed`; exact depletion and derived geometry
use `criterion_satisfied`/`criterion_violated`; tolerance diagnostics use
`within_declared_tolerance`/`outside_declared_tolerance`.

- [x] **Step 2: Write RED CLI boundary tests**

Monkeypatch the module-local production characterization, renderer, UTC clock,
and destination seams. Assert exactly one write, exact content, the fixed invocation
`pixi run held-figure5-characterize`, and successful return for a complete but
ineligible characterization. Require stdout to name the report path, say
`CHARACTERIZATION COMPLETED`, state the explicit post-qualification verdict, and
say that exit zero means evidence completion rather than path eligibility. Assert
consumer/render/write errors propagate.

Assert the normal CLI passes `_print_phase` into `characterize_figure5` and emits
the four closed phase names in order. The explicit live oracle replaces this with
its ledger-writing observer; neither caller may omit or silently ignore progress.

- [x] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py -n auto --testmon -q
```

- [x] **Step 4: Implement pure rendering and thin CLI wiring**

The renderer rejects non-exact characterization and context objects before
projection and performs no geometry, timing, independent gate decision, or write. It
renders the canonical `post_qualification_failures` result. The tool
alone owns `pathlib.Path`, production dependency wiring, required observer
injection, UTC capture, and
`write_text(markdown, encoding="utf-8")`. It does not call the post-candidate gate;
the report records the characterization verdict without turning the expected
current failure into a process failure.

Every timing label states seconds. Every report-only `PathQuality` row carries
an explicit semantic unit, including count, fraction, ratio,
benchmark-normalized millimetres, inverse/squared length, and angular units.
Dynamic cells escape both Markdown table delimiters and HTML-sensitive text;
only renderer-inserted newline markers remain literal `<br>` elements.

- [x] **Step 5: Add the Pixi CLI task and strict type coverage**

Add:

```toml
held-figure5-characterize = { cmd = "python -m tools.held_path_characterization", depends-on = ["_editable-rebuild"], description = "Characterize the default generator on Held Figure 5" }
```

Extend `types-benchmarks` with every new/modified benchmark, tool, and typecheck
file from Tasks 1-10.

- [x] **Step 6: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py
git add benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py tests/__init__.py tests/benchmarks/__init__.py tests/benchmarks/typecheck/__init__.py tests/tools/__init__.py tests/tools/typecheck/__init__.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'feat(benchmarks): report Held Figure 5 evidence'
```

---

### Task 11: Explicit live oracle and durable Phase 1 evidence

**Files:**

- Create: `tests/benchmarks/held_figure5_path_live_oracle.py`
- Create: `docs/benchmarks/held_figure5_2d_path.md` through the production CLI
- Create: `docs/superpowers/state/held-figure5-live-run.md` during the live run
- Modify: `pyproject.toml:227-235`
- Modify: `docs/held_pfeiffer_reference_pockets.md`
- Modify: `docs/benchmarks.md`
- Modify: `docs/machining_quality.md`
- Modify: `docs/segment_site_mat.md`
- Modify: `mkdocs.yml:165-172`
- Modify: this plan's status/checklist after evidence review

**Interfaces:**

- Consumes: production Figure 5 orchestration and report writer.
- Produces: `held-figure5-characterize-live` Pixi task, the committed real report,
  documentation links/claim boundaries, and Phase 1 closeout evidence.

- [x] **Step 1: Write the explicit live-oracle assertions**

The explicitly named file must load strict Figure 5, call production
characterization exactly once, require exactly 2,289 source operations, validate
the complete
TEA-audited/excluded partition and all twelve observations, write the report,
and assert the current report refuses postprocessor qualification whenever any
criterion is open or any TEA-audited operation is unresolved. Pass
`pixi run held-figure5-characterize-live` into the report writer and require the
rendered invocation to match. Print the report path, `CHARACTERIZATION COMPLETED`,
the explicit refusal/eligibility verdict, and the exit-zero claim boundary.

- [x] **Step 2: Add the explicit serial live task**

Use the established editable-build environment and include `-n auto`:

```toml
held-figure5-characterize-live = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests/benchmarks/held_figure5_path_live_oracle.py -n auto -q -s''', depends-on = ["_editable-rebuild"], description = "Run the Held Figure 5 2D path oracle" }
```

One case is logically serial even when xdist is enabled. Do not add a product
timeout or weaker fallback.

#### Task 11A: Remove the measured centre-domain construction blocker

**Files:**

- Create: `src/cutter_centre_domain_2.h`
- Create: `src/cutter_centre_domain_2.cpp`
- Modify: `src/reachable_arrangement_2.h`
- Modify: `src/reachable_arrangement_2.cpp`
- Modify: `src/coverage_bindings_2.cpp`
- Modify: `src/compas_cgal/_coverage_2.pyi`
- Modify: `benchmarks/survey.py`
- Modify: `tests/adaptive/test_reachable_domain.py`
- Modify: `tests/benchmarks/test_survey.py`
- Modify: `CMakeLists.txt`
- Modify: this plan and its design spec

**Measured evidence:** Figure 5 generation completed in 0.063150 seconds and
produced 2,289 operations. The live run entered `survey` and timed out after
1839.18 seconds without a report. Bounded probes isolated the stall before
replay: constructing `ReachableDomain2(...).center_domain()` from the 65-vertex
Figure 5 boundary exceeded 50 seconds, as did surveys restricted to 3 and 25
source operations. One first-motion engagement query took 0.000640 seconds.
Downstream replay cost remains unmeasured and is outside this repair.

**Interface:** Add one lightweight exact native cutter-centre-domain predicate
whose public `build(...)` factory owns canonical validated reach input and whose
instances expose only `contains(x, y)`. Direct Python construction is disabled.
It requires a finite query to lie inside or on the canonical outer polygon,
outside every canonical hole interior, and at exact squared distance greater
than or equal to the exact squared tool radius from every outer and hole
boundary segment. Exact tangency is accepted; the next representable point on
the illegal side is rejected. Reuse the canonical reach-input validator and its
named input errors, and preserve exact polygon-with-holes relationship
validation through one shared helper. Keep `ReachableDomain2` unchanged and
route only `benchmarks.survey.survey_path` through the new predicate. Global
nonempty/connected erosion certification remains owned by `ReachableDomain2`;
the specialized point predicate neither rebuilds that arrangement nor claims
constructor-error equivalence.

- [x] **Step 11A.1: Write and observe focused RED tests**

Cover native/Python equivalence with the existing full centre domain on small
fast convex, concave, and holed fixtures; interior, exterior, exact tangency,
and the adjacent binary64 value on the illegal side; existing named invalid
input and polygon-with-holes errors; the global-certification ownership
boundary; and the survey boundary. Run the focused Python tests with
`-n auto --testmon` and observe the missing native API failure before production
implementation.

- [x] **Step 11A.2: Implement the minimum exact predicate and survey routing**

Use canonical reach input, exact injected binary64 coordinates, kernel polygon
membership, and exact squared-distance comparisons only. Do not add a tolerance,
snapping, approximation, alternate result, generator change, replay
optimization, survey-density change, threshold change, or report change.

- [x] **Step 11A.3: Prove GREEN and bounded Figure 5 construction**

Run the focused native/Python and survey gates, configured strict benchmark
types, Ruff on touched Python, and diff hygiene. Measure construction of the new
predicate on the real 65-vertex Figure 5 boundary under a bounded command. This
gate proves only removal of the constructor stall; it makes no downstream
runtime claim.

- [x] **Step 11A.4: Commit the reviewed repair and return to Step 3**

Record RED/GREEN evidence and the bounded measurement in the durable SDD report,
then commit only Task 11A files with subject
`perf(benchmarks): bound centre-domain survey`. Independent review precedes the
unchanged live-oracle retry owned by Task 11 Step 3.

#### Task 11B: Close the measured coverage-reduction blockers

**Measured evidence:** The second live run completed generation, guarded audit,
and survey, then failed after 780.83 seconds with `CoarseCoverageGridError`.
Figure 5 spans 66.75915853538243 by 45.76553230727174 at tool radius 1.0.
The unchanged default long-axis count 200 gives a 200 by 137 grid with
0.3340549803450492 cells. Under the existing shorter-axis rounding, 667 by 457
is still too coarse at 0.10014339673363619; the exact minimum is 668 by 458 at
0.09993886008290782, or 305,944 samples. A separate bounded probe found that
the current full reachable-material owner did not construct within 218 seconds.
The first focused material-only implementation passed its exact gates but also
exceeded the mandatory 600-second Figure 5 bound. A subsequent 10-second stack
sample placed all 8,424 main-thread samples in the single material-dilation
`reach_join_parts(parts, {center})` overlay; exact-number arithmetic dominated,
while the preceding center construction completed within 34.1 seconds.

**Interface:** Keep `COVERAGE_GRID_SAMPLES = 200`, all public reducer defaults,
scientific judges, and evidence semantics unchanged. Add
`minimum_coverage_grid(spec) -> int` in `benchmarks.coverage`, sharing the
existing aspect-ratio grid calculation, and a Figure 5 quality adapter that
passes `max(COVERAGE_GRID_SAMPLES, minimum_coverage_grid(spec))` exactly once to
the canonical reducer. Wire only the production Held report tool to it.

Add one factory-only exact reachable-material owner beside `ReachableDomain2`.
Its `build(...)` factory owns canonicalization, canonical validation, and
polygon-with-holes validation; constructs the forbidden boundary band by
divide-and-conquer union of the existing exact segment-capsule parts; exact-
differences that band from the design; applies the existing nonempty and
one-component entry rule; performs exactly one existing exact reachable-
material subset-of-design containment decision with its named
`ReachableMaterialContainmentError`; and reuses `build_reachable_material_once`.
It omits only the provenance arrangement, residual, and certificate products.
Route only `benchmarks.coverage.measure_coverage` through it.

The measured revisions below supersede that initial material-owner design.
Keep global `ExactRegion2` on its original exact `oriented_side` implementation.
Only the coverage-private predicate owns cached design and center point
locators; each is bound to its final `ReachSet`, never copied or moved after
binding, and mutex-protected for shared const queries. At the excluded live-
oracle boundary, record terminal `failed` state on any report-writer
`Exception`, including observed phases, phase-entry elapsed time, no report
completion, qualification not evaluated, and exact exception type/message,
then bare re-raise.

- [x] **Step 11B.1: Amend plan/spec and observe focused RED**

Prove the Figure 5 minimum 668, the 667/668 grid boundary, a shorter-axis
rounding counterexample, one-call adapter behavior and exception identity,
production wiring, exact material parity on small convex/concave/holed fixtures,
factory-only and named input/topology errors, zero material subset decisions,
sweeps, unions, arrangements, residuals, or certificates in the private
predicate, exactly one design locator and one center locator, unchanged legacy
`ReachableDomain2` material subset decision and
`ReachableMaterialContainmentError`, and terminal live-ledger failure with
identical bare re-raise.

Historical only, not final acceptance: the rejected eager material-owner
candidate proved exactly one material subset decision and one globally cached
`ExactRegion2` locator. Revision 4 superseded both choices; they are not
satisfied or claimed by the final predicate.

- [x] **Step 11B.2: Implement only the approved exact paths**

No generator, survey, replay, public default, threshold, reference test,
approximation, fallback, identity, or report-semantic change is permitted.

The measured fix revision changes only dilation operands. First add a native
RED audit/equality gate over mixed linear/circular outer and hole boundaries.
It must prove `boundary_curves == body_operands == vertex_disks`, total parts
`2C`, one material batch union, and exact equality with the historical `3C`
construction. Then emit one unchanged sweep body per boundary curve and one
full endpoint disk per distinct cycle vertex, retaining `center` in the one
range union. Do not remove unique disks, inward sweep halves, or change set
semantics.

That exact 2C revision passed its focused equality/count gates but again
exceeded the 600-second Figure 5 bound. Revision 2 preserves the interleaved 2C
sweep operand order and CGAL range tree, range-unions only those sweep operands,
then performs one binary exact join with the already-completed center set at the
root. Before production, extend the native RED/equality gate to prove one sweep
range union, one center root join, no center operand below root, and exact
equality with the flat historical union. Do not reorder operands, introduce a
custom tree, or remove one-sided sweep geometry.

Revision 2 also exceeded the 600-second Figure 5 bound. Replacing Boost's exact
number backend with GMP is rejected here: the exact kernel is shared across the
connected native target, and enabling GMP/MPFR therefore changes repository-wide
ABI and binary-wheel dependency policy rather than this isolated coverage path.

The user-approved revision 3 is a bounded, falsifiable experiment only. Add a
candidate beside the historical 3C and current 2C paths that labels every closed
center and sweep boundary cycle and invokes CGAL Minkowski_sum_2's
`Union_of_curve_cycles_2` once, producing the same exact reachable set through
one arrangement. Before routing the material-only factory through the candidate,
the native gate must prove historical 3C = current 2C = candidate on a mixed
line/arc outer-and-hole fixture; equal component and hole counts; center subset
of candidate and candidate subset of a containing design; complete unique
component/cycle labels; and exactly one candidate arrangement. The class is an
auxiliary API and its headers are GPL-3.0-or-later or commercially licensed, so
API fragility and package-license compatibility remain explicit adoption risks.
No old path is removed and there is no fallback or approximation. A real Figure
5 factory run has a hard 120-second bound; timeout, assertion, or exact mismatch
rejects the candidate immediately and forbids a production commit.

Revision 3 passed its focused exact and structural gates but produced no factory
result before the external 120-second Figure 5 bound. It is rejected; no query
projection, live oracle, or production commit follows from this experiment.

Revision 4 supersedes the rejected uncommitted material-region owner only for
coverage. Add a factory-only `ReachableMaterialPredicate2` that constructs the
same validated exact center set once, enforces its nonempty/one-component rule,
owns one cached exact center locator, and answers
`q in C or exact_distance(q, boundary(C)) <= r`. For a nonempty closed erosion
`C = D eroded by B(r)`, this is exactly membership in `C + B(r)`, while the
erosion definition proves `C + B(r)` is a subset of `D`. This theorem replaces
the eager material construction and global subset decision only on the private
coverage path. Legacy `ReachableDomain2` retains its material set, explicit
subset decision, and `ReachableMaterialContainmentError` unchanged.

The predicate scans every outer and hole x-monotone center-boundary curve with
exact line squared distance or exact circular-arc endpoint/radial-projection
tests. No tolerance, fallback, spatial index, or material sweep/union/arrangement
is permitted. Native parity with legacy material on convex, concave, holed,
narrow-bay, and mixed-arc fixtures plus exact tangency/adjacent-binary64 and
radius-branch tests precedes a hard-120-second stratified projection and full
305,944-query Figure 5 batch. Only `benchmarks.coverage` uses this predicate.

The first revision-4 ordering passed those exact gates but failed the bounded
projection: construction took 1.155735 seconds and 8,192 stratified queries
took 7.352074 seconds, projecting 275.731305 seconds total. The measured cause
is that every point outside `C`, including the 30.7251 percent outside `D`,
scans the full center boundary. Because `C + B(r)` is already proven a subset
of `D`, preserve the validated exact design as a second immutable
`ExactRegion2` and decide in the exact order `q not in D -> false`, `q in C ->
true`, then boundary distance. A read-only simulation reduced boundary-scan
candidates from 3,021 to 504 of 8,192 and projected 36.443539 seconds for the
full query batch. The audit must prove exactly one design locator and one center
locator. No batch API or spatial index is permitted unless this measured
ordering fails the unchanged 120-second gates.

The design-first ordering passed. Its identical 8,192-point stratified run took
0.917544 seconds after 1.174590 seconds construction, projecting 35.441826
seconds total. The complete 668 by 458 batch then evaluated all 305,944 points
in 33.795099 seconds after 1.167958 seconds construction, 34.963057 seconds
total, and classified 211,804 points inside. No batch API or spatial index was
needed.

- [x] **Step 11B.3: Prove GREEN and bounded Figure 5 material construction**

Run focused native/Python gates, configured strict types, Ruff, and diff hygiene.
Measure the real Figure 5 predicate construction, a stratified exact membership
batch, and the complete 305,944-query batch under the revision-4 external
bounds. Preserve the exact verdict count across locator-hardening changes.

- [x] **Step 11B.4: Correct evidence, report, review, and commit**

Correct the current live ledger to the failed second run while preserving its
start, phases, and 778.885679-second last phase timestamp; record pytest failure
at 780.83 seconds, exact `CoarseCoverageGridError`, no report, and qualification
not evaluated. Preserve the prior timeout in Git history. Complete final Task
11B independent review and checkpoint commit, then return to Step 3.

- [ ] **Step 3: Run the live oracle under the operator budget**

**Retired on 2026-09-05:** do not run this step under the active paper-figure
mission. The interrupted run and reason are recorded in
`docs/superpowers/state/held-figure5-live-run.md`.

Run:

```bash
pixi run held-figure5-characterize-live
```

Operator rule: stop after 30 minutes if it has not completed. On budget
exhaustion, write no new characterization report. The live oracle writes its UTC
run start and updates the active phase before each consumer in
`docs/superpowers/state/held-figure5-live-run.md`; record elapsed time and timeout
there, mark any pre-existing report as prior/stale evidence, leave Phase 1
blocked, and do not begin generator repair. On success, require the report UTC
instant to postdate the recorded run start. Report publication must first fully
write the sibling `.pending.md` file, then atomically `Path.replace` the report;
ordinary write or replacement exceptions preserve any prior report and may
leave the pending file as evidence.

- [ ] **Step 4: Inspect the generated report and obtain independent review**

**Retired on 2026-09-05:** no generated report exists from the interrupted run.

Verify the banner, normalized units, separately measured consumer timings, three-way
engagement partition, witnesses, all twelve criteria, sampled/certified claim
boundaries, and refusal verdict against the live objects. Independent review
must reject any unsupported native-audit, Fanuc, Held-superiority, or
machine-release claim. Geometry claims must remain explicitly projection-only.
The report must expose raw coverage `nx`, `ny`, reachable sample count, uncut
sample count, and cell area; decide any report correction only from this live
inspection.

- [ ] **Step 5: Integrate durable documentation**

**Retired on 2026-09-05:** this documentation integration depended on Step 3.

Link the report from the Held reference page and benchmark index; add it to
MkDocs navigation. Update `docs/machining_quality.md` with the single canonical
reducer/evidence vocabulary. Update `docs/segment_site_mat.md` with exact maturity:
guarded/sampled Figure 5 2D characterization only, with native audit,
target-specific Fanuc translation, machine setup, and Held superiority still
unproven.

- [ ] **Step 6: Run final repository gates**

**Retired on 2026-09-05:** these are characterization closeout gates, not the
paper-figure reproduction gate.

```bash
pixi run affected
pixi run types-benchmarks
pixi run ruff format benchmarks tools tests/benchmarks tests/tools
pixi run ruff check benchmarks tools tests
pixi run red-manifest
pixi run docs
pixi run plan-headers
git diff --check
```

Expected: affected tests introduce no new failure, strict types/Ruff/docs pass,
and the red manifest reports exactly the inherited seventeen reds.

- [ ] **Step 7: Close the plan and commit Phase 1 evidence**

**Retired on 2026-09-05:** preserve the implemented checkpoints; do not present
the uncompleted live characterization as complete.

Mark completed tasks and the plan status truthfully. Commit only the named
production, test, report, documentation, configuration, spec, and plan files;
do not stage `tmp/`:

```bash
git add benchmarks/errors.py benchmarks/held_consumer_adapters.py benchmarks/held_path_characterize.py benchmarks/held_path_evidence.py benchmarks/held_path_report.py benchmarks/held_path_snapshot.py benchmarks/held_post_qualification.py benchmarks/quality.py benchmarks/quality_observations.py benchmarks/survey.py benchmarks/units.py
git add tests/benchmarks/held_figure5_path_live_oracle.py tests/benchmarks/test_held_path_characterize.py tests/benchmarks/test_held_path_evidence.py tests/benchmarks/test_held_path_report.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/test_held_post_qualification.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_quality_observations.py tests/benchmarks/test_survey.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py tests/benchmarks/typecheck/quality_observations_contract.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py tools/held_path_characterization.py
git add docs/benchmarks/held_figure5_2d_path.md docs/benchmarks.md docs/held_pfeiffer_reference_pockets.md docs/machining_quality.md docs/segment_site_mat.md docs/superpowers/plans/2026-09-02-held-figure5-2d-path-characterization.md docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md docs/superpowers/state/held-figure5-live-run.md docs/superpowers/state/held-figure5-quality-parity.md mkdocs.yml pyproject.toml
git diff --cached --name-only
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'docs(held): record Figure 5 path evidence'
```

Phase 1 is complete only after the live report exists and independent review is
clean. Its measured failure signature is then the sole input to a new Phase 2
generator-policy design; this plan does not repair the generator or emit G-code.
