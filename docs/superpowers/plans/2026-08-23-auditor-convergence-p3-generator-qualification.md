# Auditor Convergence P3 Generator Qualification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve every forced engagement decision structurally, replace
monotonic search assumptions with one finite ladder authority, and make the
regulated-generator quality gate measure what it claims.

**Architecture:** Add a content-addressed override ledger shared by generator
results, then introduce one typed exhaustive ladder-search path alongside the
existing algorithms and migrate consumers with characterization tests. Repair
quality-test non-vacuity and run the authoritative audit over regulated and
baseline generators under identities that name every cap and override.

**Tech Stack:** Python 3.12, COMPAS, adaptive typed units, CCAN, SHA-256, Pixi,
pytest-xdist, pytest-testmon, Hypothesis, Ruff, strict mypy, MkDocs.

**Spec:** `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md`

## Global Constraints

- P2 acceptance commit is mandatory.
- Warnings remain human reporting only; no consumer parses them.
- Forced loops and forced advances cross every public result boundary.
- Engagement admissibility is non-monotonic; no bisection may decide it.
- Existing algorithms remain until the new path is independently validated and
  the user authorizes removal.
- Public/shared boundaries carry frame and units, not bare physical floats.
- Every pytest command uses `-n auto`; no skipped or expected-failure tests.
- Product-quality failures remain failing assertions in `quality-gate`.

---

### Task 1: Define the content-addressed override ledger

**Files:**

- Create: `src/compas_cgal/engagement_override.py`
- Create: `tests/test_engagement_override.py`
- Modify: `src/compas_cgal/engagement_toolpath.py`
- Modify: `src/compas_cgal/engagement_radial_toolpath.py`
- Modify: `src/compas_cgal/engagement_rho_toolpath.py`
- Modify: `src/compas_cgal/engagement_ordered_toolpath.py`
- Modify: `src/compas_cgal/engagement_spiral_entry_toolpath.py`

**Interfaces:**

- Consumes: operation/candidate digests, exact refusal or unresolved witness,
  stock lineage.
- Produces: `EngagementOverrideRecord.build(...)`,
  `EngagementOverrideLedger.build(...)`, `require_no_overrides()`.

- [ ] **Step 1: Write RED invariant tests**

```python
def test_override_ledger_binds_forced_advance_preimage() -> None:
    record = _forced_advance_record()
    ledger = EngagementOverrideLedger.build((record,))

    assert record.motion_certificate_digest in ledger.canonical_bytes
    assert ledger.digest == hashlib.sha256(ledger.canonical_bytes).digest()
    with pytest.raises(EngagementOverridePresentError, match="forced_advance"):
        ledger.require_no_overrides()


def test_warning_filter_cannot_erase_override() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = _forced_generator_result()

    assert result.override_ledger.records
```

Mutation tests cover operation index/digest, reason, refusal witness, selected
alternative, stock lineage, duplicate operation ownership, and record order.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/test_engagement_override.py -n auto -q
```

- [ ] **Step 3: Implement records and ledger**

```python
OverrideReason: TypeAlias = Literal[
    "forced_entry",
    "forced_loop",
    "forced_advance",
]


@dataclass(frozen=True)
class EngagementOverrideRecord:
    operation_index: int
    operation_digest: bytes
    reason: OverrideReason
    refusal_witness_digest: bytes
    selected_alternative_digest: bytes
    pre_motion_stock_lineage: bytes
```

Factories enforce exact types, 32-byte digests, unique operation ownership, and
canonical order. The ledger derives counts and digest.

- [ ] **Step 4: Add ledgers to public generator results**

Add a result type for each generator that does not already own one. Each result
contains operations/polyline plus `override_ledger`; do not overload one method
body with two calling conventions. Preserve current warnings as summaries of
the ledger.

- [ ] **Step 5: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run pytest -- tests/test_engagement_override.py tests/test_engagement_toolpath.py tests/test_engagement_radial_toolpath.py tests/test_engagement_rho_toolpath.py tests/test_engagement_ordered_toolpath.py tests/test_engagement_spiral_entry_toolpath.py -n auto --testmon -q
git diff --check
git add src/compas_cgal/engagement_override.py src/compas_cgal/engagement_*toolpath.py tests/test_engagement_override.py
git commit -m "feat(engagement): retain forced decisions"
```

### Task 2: Add one typed finite ladder-search authority

**Files:**

- Create: `src/compas_cgal/engagement_search.py`
- Create: `tests/test_engagement_search.py`
- Modify: `src/compas_cgal/engagement_toolpath.py`
- Modify: `src/compas_cgal/engagement_radial_toolpath.py`
- Modify: `src/compas_cgal/engagement_rho_toolpath.py`
- Modify: `src/compas_cgal/engagement_ordered_toolpath.py`
- Modify: `src/compas_cgal/engagement_spiral_entry_toolpath.py`

**Interfaces:**

- Consumes: typed station sequence, exact stock/cap authority, bounded integer
  window.
- Produces: `select_largest_admissible_advance(...) -> AdmissibleAdvance`.

- [ ] **Step 1: Write RED exhaustive-oracle tests**

```python
@given(_nonmonotonic_station_case())
def test_ladder_search_equals_exhaustive_oracle(case: SearchCase) -> None:
    selected = select_largest_admissible_advance(**case.arguments)
    expected = exhaustive_admissible_advance(**case.arguments)

    assert selected == expected
```

Include the pinned case where bisection disagrees with downward enumeration,
empty admissible set, first/last station, and exact cap equality.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/test_engagement_search.py -n auto -q
```

- [ ] **Step 3: Implement typed search inputs/result**

`SearchStation` carries `Point2[WorldXY]`, `Spacing`, and stable ordinal.
`AdmissibleAdvance` carries selected ordinal, refusal witness digest when
forced, and exact candidate digest. The function enumerates the complete finite
window and ranks admissible candidates; it never assumes monotonicity.

- [ ] **Step 4: Migrate one consumer at a time**

For each generator, first pin current operations, selected ordinals, and
override ledger on representative fixtures. Switch its import/call to the new
authority, rerun its focused suite, then proceed to the next generator. Do not
delete old private functions.

- [ ] **Step 5: Prove no regulated call site uses bisection**

```bash
rg -n '_largest_admissible_advance|bisect' src/compas_cgal/engagement_*toolpath.py
```

Expected: old definitions remain only as non-authoritative validation-era code;
regulated public generator call graphs import `select_largest_admissible_advance`.

- [ ] **Step 6: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run pytest -- tests/test_engagement_search.py tests/test_engagement_toolpath.py tests/test_engagement_radial_toolpath.py tests/test_engagement_rho_toolpath.py tests/test_engagement_ordered_toolpath.py tests/test_engagement_spiral_entry_toolpath.py -n auto --testmon -q
git add src/compas_cgal/engagement_search.py src/compas_cgal/engagement_*toolpath.py tests/test_engagement_search.py
git commit -m "feat(engagement): use finite ladder search"
```

### Task 3: Type the engagement-generator geometry boundary

**Files:**

- Create: `src/compas_cgal/engagement_geometry.py`
- Create: `tests/test_engagement_geometry.py`
- Modify: `src/compas_cgal/engagement_*toolpath.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces:**

- Consumes: `Point2[WorldXY]`, `Vector2[WorldXY]`, `ToolRadius`, `GuideRadius`,
  `Spacing`, `EngagementCap`.
- Produces: `GuideStation.build(...)`, `RegulationPolicy.build(...)`.

- [ ] **Step 1: Write RED scalar-confusion tests**

Require factories to reject a point where a tangent belongs, tool radius where
guide radius belongs, non-unit/zero direction, nonfinite values, and the known
rho-generator rotated-radius-as-tangent construction.

- [ ] **Step 2: Implement focused typed records**

```python
@dataclass(frozen=True)
class GuideStation:
    center: Point2[WorldXY]
    tangent: Vector2[WorldXY]
    guide_radius: GuideRadius


@dataclass(frozen=True)
class RegulationPolicy:
    tool_radius: ToolRadius
    engagement_cap: EngagementCap
    advance_window: int
```

Factories own validation; the records do not combine per-run counters or
mutable stock state.

- [ ] **Step 3: Migrate shared/private boundaries**

Convert to raw COMPAS/double values only at the legacy geometry/native seam.
Use separate scalar-positional and sequence `@overload` factories for vector
quantities. No conditional compatibility path.

- [ ] **Step 4: Run strict typing and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/test_engagement_geometry.py tests/test_engagement_search.py -n auto --testmon -q
git add src/compas_cgal/engagement_geometry.py src/compas_cgal/engagement_*toolpath.py tests/test_engagement_geometry.py tests/adaptive/typecheck/consumer_contract.py
git commit -m "refactor(engagement): type guide geometry"
```

### Task 4: Repair quality-gate non-vacuity

**Files:**

- Create: `tests/benchmarks/test_quality_nonvacuity.py`
- Modify: `benchmarks/quality.py`

**Interfaces:**

- Consumes: regulated generator results and authoritative audit records.
- Produces: comparison-count floors and geometric bridge qualification.

- [ ] **Step 1: Add a RED adversarial tangent fixture without changing reference assertions**

Create a new test whose Line/Arc alternation carries deliberately wrong
tangents. Invoke the current continuity helper and assert it must reject. The
test should expose the 7.7% structural exemption.

- [ ] **Step 2: Add comparison-count accounting**

The checker returns or records `eligible`, `checked`, and
`geometrically_exempt` counts. A fixture promising engaged transitions raises
`VacuousQualityMetricError` when `checked` is below its declared floor.

- [ ] **Step 3: Make bridge exemption geometric**

Verify exact shared endpoint, operation continuity, and the documented tangent
relation using named COMPAS tolerance predicates with an explicitly justified
angular tolerance. Structural Line/Arc shape alone never exempts.

- [ ] **Step 4: Replace absent-arc soft return with fixture failure**

Add a new fixture precondition used by tests that promise arcs. Preserve the
reference test file's assertions; production/fixture construction must satisfy
them rather than silently returning.

- [ ] **Step 5: Run GREEN and commit**

```bash
pixi run pytest -- tests/benchmarks/test_quality_nonvacuity.py tests/benchmarks/test_quality.py tests/test_toolpath.py -n auto --testmon -q
pixi run lint
git add benchmarks/quality.py tests/benchmarks/test_quality_nonvacuity.py
git commit -m "test(quality): enforce nonvacuity"
```

If satisfying the new test requires editing a protected reference assertion,
stop and present the exact conflict instead of changing it.

### Task 5: Run regulated corpus and gate P3

**Files:**

- Modify: `benchmarks/runner.py`
- Modify: `benchmarks/models.py`
- Modify: `benchmarks/cli.py`
- Create: `tests/benchmarks/test_regulated_runner.py`
- Create: `tests/benchmarks/test_regulated_quality.py`
- Modify: `docs/machining_quality.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: authoritative audit, override ledgers, regulated and baseline
  generators.
- Produces: identity-bound corpus rows and separate `quality-gate` result.

- [ ] **Step 1: Add RED corpus-subject tests**

Require every measurement row to identify generator, cap, build, input,
override ledger, and audit report. Require at least one cap in `40..100` degrees
that produces distinct regulated-generator paths.

- [ ] **Step 2: Register regulated generators as primary subjects**

Keep the unregulated generator under an explicit baseline identifier. The CLI
accepts exact subject names and rejects unknown/duplicate subjects.

- [ ] **Step 3: Define separate Pixi tasks**

Add `benchmark-instrument` for green schema/metric/figure contracts and
`quality-gate` for product criteria. Neither task skips tests; the latter may
exit nonzero and is recorded honestly.

- [ ] **Step 4: Run P3 gates and record results**

```bash
pixi run benchmark-instrument
pixi run quality-gate
pixi run lint
pixi run types-adaptive
pixi run -e docs docs
git diff --check
```

Record both exit codes, all metric values, and override counts. Do not weaken a
threshold to obtain green.

- [ ] **Step 5: Commit documentation and progress**

```bash
git add benchmarks pyproject.toml tests/benchmarks docs/machining_quality.md .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "feat(bench): qualify regulated generators"
```
