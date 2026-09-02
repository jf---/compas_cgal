# Held Figure 5 2D Path Characterization Implementation Plan

> **status: ready** - the governing design is approved; execution starts at
> Task 1 and must stop at the Task 5 parity checkpoint for explicit approval
> before superseded quality reducers are removed.

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
- Preserve the exact known-red membership: twelve quality cells, the benchmark
  translation-invariance cell, and four adaptive cells. Never edit a reference
  assertion to make it pass.
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

- [ ] **Step 1: Write RED runtime tests for each observation-unit validator**

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

- [ ] **Step 2: Write the strict type contract**

The contract must prove the six scalar domains remain distinct and geometry
continues to use the existing unit/frame vocabulary:

```python
seconds = assert_type(seconds_value(1.0, name="audit"), Seconds)
degrees = assert_type(degrees_value(80.0, name="cap"), Degrees)
index = assert_type(operation_index(3, operation_count=4), OperationIndex)
point = assert_type(Point2[WorldXY].build(1.0, 2.0), Point2[WorldXY])
assert_type(Point3[WorldXYZ].build(1.0, 2.0, 0.0), Point3[WorldXYZ])
```

- [ ] **Step 3: Run RED tests**

Run:

```bash
pixi run pytest -- tests/benchmarks/test_units.py -n auto --testmon -q
pixi run mypy --strict --warn-unused-ignores tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: collection/type checking fails because `benchmarks.units` and the new
failures do not exist.

- [ ] **Step 4: Implement the minimal observation-unit module and failures**

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

- [ ] **Step 5: Extend the `types-benchmarks` task and run GREEN gates**

Add `benchmarks/units.py` and the new type-contract file to the explicit mypy
file list, then run:

```bash
pixi run pytest -- tests/benchmarks/test_units.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/units.py benchmarks/errors.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/units.py benchmarks/errors.py tests/benchmarks/test_units.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: all named tests and strict typing pass.

- [ ] **Step 6: Commit**

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

- [ ] **Step 1: Write RED factory tests for the closed primitive union**

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

- [ ] **Step 2: Write RED mutation and malformed-input tests**

Parameterize one-field changes covering operation order/count, role,
`path_index`, clockwise, endpoints, centre, Z, frame axes, radius, arc angles,
and tangents. Every mutation must raise `MutatedHeldToolpathError`. Factories
must reject non-finite geometry, non-positive radius, non-unit or non-orthogonal
axes, invalid angles, unsupported geometry, and malformed tangents through
`InvalidHeldOperationSnapshotError`. Direct construction must raise `TypeError`.

- [ ] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py -n auto --testmon -q
```

Expected: import fails because the snapshot module does not exist.

- [ ] **Step 4: Implement structural snapshot factories**

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

- [ ] **Step 5: Run GREEN runtime, type, and formatting gates**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_snapshot.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_snapshot.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_path_snapshot.py tests/benchmarks/test_held_path_snapshot.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

- [ ] **Step 6: Commit**

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

- [ ] **Step 1: Write RED tests for retained sample positions and predicates**

Tests must establish that a circle has 45 seam-unique samples, while open lines
and arcs have 46 endpoint-inclusive samples at the standard 45 intervals. Assert
the first/last typed positions and:

```python
assert motion.cap_exceeded is any(sample.cap_exceeded for sample in motion.samples)
```

Monkeypatch only the reporting angle in a focused unit test to prove the Boolean
comes directly from the exact predicate result and is never reconstructed by
comparing `engagement_deg` with the cap.

- [ ] **Step 2: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_survey.py -n auto --testmon -q
```

Expected: the new `EngagementSample` fields are absent.

- [ ] **Step 3: Extend `EngagementSample` and `_measure_motion` additively**

For each `(distance, x, y)` station, retain the exact returned `exceeded` value
and `Point2[WorldXY].build(x, y)`. Update the docstring so
`engagement_deg` is explicitly reporting-only, per-position predicates are exact,
and absence across the finite station set remains sampled-negative evidence.
Do not change `QUALITY_SAMPLES_PER_MOTION` or `_motion_samples`.

- [ ] **Step 4: Prove existing survey/quality behavior is unchanged**

```bash
pixi run pytest -- tests/benchmarks/test_survey.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/survey.py tests/benchmarks/test_survey.py tests/benchmarks/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/survey.py tests/benchmarks/test_survey.py tests/benchmarks/typecheck/held_path_characterization_contract.py
```

Expected: selected tests pass; no existing metric changes.

- [ ] **Step 5: Commit**

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

- [ ] **Step 1: Write RED criterion-factory tests**

Pin all twelve criterion names, thresholds, evidence kinds, and outcome
vocabulary. Required values are exactly `0`, the case cap, and `2` tool radii;
never infer them from current measurements. A wrong evidence kind, unit domain,
non-finite value, negative count, invalid index, or inconsistent outcome must
raise `InvalidHeldPathEvidenceError`.

- [ ] **Step 2: Write RED simple-attribution tests**

Using synthetic `PathSurvey` values, require:

- uncut stock remains a spatial aggregate only;
- gouging, unsafe rapid, zero-length, degenerate-loop, redundant, cap-exceeded,
  and slotting counts equal their exact unique operation-index tuples;
- cap-exceeded indices come from `sample.cap_exceeded`, not reporting degrees;
- count observations reject attribution cardinality disagreement.

- [ ] **Step 3: Write RED junction and step-attribution tests**

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

- [ ] **Step 4: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_quality_observations.py -n auto --testmon -q
pixi run mypy --strict --warn-unused-ignores tests/benchmarks/typecheck/quality_observations_contract.py
```

Expected: import fails because the canonical observation module is absent.

- [ ] **Step 5: Implement the canonical findings and typed assessment**

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

- [ ] **Step 6: Run GREEN focused gates**

```bash
pixi run pytest -- tests/benchmarks/test_quality_observations.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/typecheck/quality_observations_contract.py
pixi run ruff check benchmarks/quality_observations.py tests/benchmarks/test_quality_observations.py tests/benchmarks/typecheck/quality_observations_contract.py
```

- [ ] **Step 7: Commit additive observations**

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

- [ ] **Step 1: Add parity assertions without changing the old verdict path**

Keep `_violations` and the final gate assertion untouched. Add a helper that
compares every gated `PathQuality` field with its typed observation, then add it
to synthetic machinery and invariant cases. Counts must reduce from every
attributed operation/pair; maximum observations must name the exact winning pair
or `None` when the maximum is zero.

- [ ] **Step 2: Add parity to all twelve `QUALITY_GATE_CASES`**

In `test_the_generated_path_is_worth_running`, generate once and run the
unchanged quality path. Test spies capture its single `PathSurvey` and
`CoverageEstimate`; feed those exact objects plus the operation snapshot to the
additive assessment. Assert parity, then leave the existing
`_violations`/`assert not violations` lines unchanged. Per cell, prove exactly
one survey and one coverage measurement. The twelve tests must remain red for
the same reasons after the parity assertion passes.

- [ ] **Step 3: Run synthetic and invariant parity GREEN**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
```

Expected: all selected tests pass.

- [ ] **Step 4: Run all twelve gate cells and reconcile expected reds**

```bash
zsh -c 'pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q; red_status=$?; if [[ $red_status -ne 1 ]]; then exit 2; fi; pixi run red-manifest'
```

Expected: exactly twelve `test_the_generated_path_is_worth_running` failures
with unchanged IDs and criterion messages; `red-manifest` confirms the repository
still has exactly seventeen declared reds overall: twelve quality, one benchmark
translation-invariance, and four adaptive.
Exit status `1` is required for the focused expected-red invocation; collection,
configuration, interruption, and infrastructure statuses are blockers, not reds.

- [ ] **Step 5: Record parity evidence**

Write `held-figure5-quality-parity.md` with:

- command and date;
- every old/new field compared;
- attribution cardinality and extremum-pair results;
- exact twelve quality-red membership and seventeen-red repository membership;
- explicit statement that current production still uses the old reducers; and
- the removal decision requested from Jelle.

- [ ] **Step 6: Run focused hygiene and commit**

```bash
pixi run types-benchmarks
pixi run ruff format tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py
pixi run ruff check benchmarks tests/benchmarks
git diff --check
git add tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py docs/superpowers/state/held-figure5-quality-parity.md
GIT_AUTHOR_NAME='Jelle Feringa' GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' GIT_COMMITTER_NAME='Jelle Feringa' GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' git commit -m 'test(benchmarks): prove quality observation parity'
```

- [ ] **Step 7: STOP and request explicit removal approval**

Present the committed parity evidence and exact known-red reconciliation. Do not
start Task 6, remove an old reducer, or redirect `measure_quality` until Jelle
explicitly approves convergence.

---

### Task 6: Converge `measure_quality` to the canonical reducers

**Precondition:** Jelle explicitly approved removal after reviewing Task 5.

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

- [ ] **Step 1: Write RED single-survey and compatibility tests**

Instrument `survey_path`, `measure_coverage`, and `assess_path_quality`. Require
one call each, identical public `PathQuality` output, and the existing named
exceptions for zero-length paths and invalid grids.

- [ ] **Step 2: Run RED focused tests**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'single_survey or compatibility' -n auto --testmon -q
```

Expected: the new reducer API does not exist or the compatibility path does not
delegate.

- [ ] **Step 3: Add `QualityEvidence` and route `measure_quality` through it**

`reduce_quality_evidence` measures coverage exactly once, evaluates the
canonical assessment once, and builds all five existing `PathQuality` groups.
For the twelve gated fields, take values only from the assessment. Continue to
compute report-only fields through their existing reducers.

- [ ] **Step 4: Validate the routed path before deletion**

With old bodies present but unused, run the compatibility tests, the non-gate
quality/invariant set excluding the declared translation red, all twelve gate
cells, `red-manifest`, strict types, and Ruff. Record that the routed canonical
path preserves the exact seventeen-red repository set.

- [ ] **Step 5: Remove only superseded decision bodies**

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

- [ ] **Step 6: Replace transition parity hooks with permanent contract tests**

Remove the temporary parity spies from Task 5. Keep literal
threshold assertions, criterion name/value/evidence tests, unchanged
`QUALITY_GATE_CASES`, unchanged test IDs, and unchanged final gate assertion.

- [ ] **Step 7: Rerun the same quality and type gates after deletion**

```bash
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py tests/benchmarks/test_qualityfigures.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
zsh -c 'pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q; red_status=$?; if [[ $red_status -ne 1 ]]; then exit 2; fi; pixi run red-manifest'
pixi run types-benchmarks
pixi run ruff format benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py
pixi run ruff check benchmarks/quality.py benchmarks/quality_observations.py tests/benchmarks
```

Expected: non-gate/non-translation tests pass; exactly the same twelve quality
cells remain red; the repository red manifest still reports exactly seventeen
declared reds.

- [ ] **Step 8: Commit**

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

- [ ] **Step 1: Write RED engagement-partition tests**

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

- [ ] **Step 2: Write RED cross-consumer and factory tests**

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

- [ ] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_evidence.py -n auto --testmon -q
```

- [ ] **Step 4: Implement construction from authoritative index sets**

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

- [ ] **Step 5: Run GREEN gates and commit**

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

- [ ] **Step 1: Write RED refusal tests for every open gate**

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

- [ ] **Step 2: Write RED authority tests**

Require direct candidate construction to fail. Successful `.build(...)` must
retain the closed characterization, exact immutable snapshot, normalized tool
diameter, 80-degree cap, and Figure 5 case context. The functional helper must
delegate to `.build(...)` and return the candidate, never `None` or a free
snapshot.

- [ ] **Step 3: Implement one validating construction path**

Implement one public pure `post_qualification_failures` collector. `.build(...)`
consumes it and owns the only conversion from an empty failure tuple into a
candidate; the functional helper calls the factory. Do not create a second
Boolean gate. Failure messages list every open criterion plus
certified/demonstrated/unresolved counts.

- [ ] **Step 4: Run GREEN gates and commit**

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

- [ ] **Step 1: Write RED adapter tests**

Monkeypatch `audit_toolpath_engagement` and assert the adapter passes exactly
`spec.polygon`, the same result object, `spec.tool_diameter`,
`spec.tea_cap_rad`, and `list(spec.holes)`, then returns the same report. No
exception is caught or rewritten.

- [ ] **Step 2: Write RED orchestration order and identity tests**

Use injected spies and a deterministic typed clock. Require this exact sequence:
load Figure 5, generate once, snapshot, assert unchanged, audit, assert unchanged,
survey, assert unchanged, reduce the same survey/snapshot, build evidence. Both
replay consumers receive the object-identical generated result. Record generation,
audit, survey, and coverage/reduction durations independently. The required phase
observer receives `generation`, `guarded_audit`, `survey`, and
`quality_reduction` immediately before their corresponding calls.

- [ ] **Step 3: Write RED mutation and exception tests**

Have each replay consumer mutate one operation before returning; the following
post-call assertion must raise `MutatedHeldToolpathError`. Verify unexpected
generator/auditor/survey/reducer exceptions propagate unchanged and no partial
characterization exists.

- [ ] **Step 4: Implement the adapter and orchestration**

Keep the adapter as a value-only signature translation. Keep orchestration free
of rendering, file writes, command-line parsing, and post-candidate construction.

- [ ] **Step 5: Run GREEN gates and commit**

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
- Modify: `pyproject.toml:227-233`

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
  - `write_held_figure5_report(*, pixi_command: str, path: Path =
    DEFAULT_REPORT_PATH) -> HeldFigure5Characterization`
  - `main() -> None`

- [ ] **Step 1: Write RED report-context and rendering tests**

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

- [ ] **Step 2: Write RED CLI boundary tests**

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

- [ ] **Step 3: Run RED tests**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py -n auto --testmon -q
```

- [ ] **Step 4: Implement pure rendering and thin CLI wiring**

The renderer performs no geometry, timing, independent gate decision, or write. It
renders the canonical `post_qualification_failures` result. The tool
alone owns `pathlib.Path`, production dependency wiring, UTC capture, and
`write_text(markdown, encoding="utf-8")`. It does not call the post-candidate gate;
the report records the characterization verdict without turning the expected
current failure into a process failure.

- [ ] **Step 5: Add the Pixi CLI task and strict type coverage**

Add:

```toml
held-figure5-characterize = { cmd = "python -m tools.held_path_characterization", depends-on = ["_editable-rebuild"], description = "Characterize the default generator on Held Figure 5" }
```

Extend `types-benchmarks` with every new/modified benchmark, tool, and typecheck
file from Tasks 1-10.

- [ ] **Step 6: Run GREEN gates and commit**

```bash
pixi run pytest -- tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py -n auto --testmon -q
pixi run types-benchmarks
pixi run ruff format benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py
pixi run ruff check benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py
git add benchmarks/held_path_report.py tools/held_path_characterization.py tests/benchmarks/test_held_path_report.py tests/tools/test_held_path_characterization.py tests/tools/typecheck/held_path_characterization_contract.py pyproject.toml
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

- [ ] **Step 1: Write the explicit live-oracle assertions**

The explicitly named file must load strict Figure 5, call production
characterization exactly once, require exactly 2,289 source operations, validate
the complete
TEA-audited/excluded partition and all twelve observations, write the report,
and assert the current report refuses postprocessor qualification whenever any
criterion is open or any TEA-audited operation is unresolved. Pass
`pixi run held-figure5-characterize-live` into the report writer and require the
rendered invocation to match. Print the report path, `CHARACTERIZATION COMPLETED`,
the explicit refusal/eligibility verdict, and the exit-zero claim boundary.

- [ ] **Step 2: Add the explicit serial live task**

Use the established editable-build environment and include `-n auto`:

```toml
held-figure5-characterize-live = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests/benchmarks/held_figure5_path_live_oracle.py -n auto -q -s''', depends-on = ["_editable-rebuild"], description = "Run the Held Figure 5 2D path oracle" }
```

One case is logically serial even when xdist is enabled. Do not add a product
timeout or weaker fallback.

- [ ] **Step 3: Run the live oracle under the operator budget**

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
instant to postdate the recorded run start.

- [ ] **Step 4: Inspect the generated report and obtain independent review**

Verify the banner, normalized units, separately measured consumer timings, three-way
engagement partition, witnesses, all twelve criteria, sampled/certified claim
boundaries, and refusal verdict against the live objects. Independent review
must reject any unsupported native-audit, Fanuc, Held-superiority, or
machine-release claim.

- [ ] **Step 5: Integrate durable documentation**

Link the report from the Held reference page and benchmark index; add it to
MkDocs navigation. Update `docs/machining_quality.md` with the single canonical
reducer/evidence vocabulary. Update `docs/segment_site_mat.md` with exact maturity:
guarded/sampled Figure 5 2D characterization only, with native audit,
target-specific Fanuc translation, machine setup, and Held superiority still
unproven.

- [ ] **Step 6: Run final repository gates**

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
