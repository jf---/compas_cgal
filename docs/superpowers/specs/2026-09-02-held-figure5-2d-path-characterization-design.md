# Held Figure 5 Machinability Characterization Design

> **status: draft for review** - the architectural brainstorm is approved;
> implementation starts only after this written design is approved.

## Goal

Characterize the current default generator on the fixed Held-Pfeiffer Figure 5
reference pocket through the repository's existing engagement and machining
quality consumers. The result is a durable failure report with
operation-attributed engagement witnesses, aggregate machining criteria, and a
fail-closed release verdict. This phase does not repair generation and does not
emit G-code.

## Why this phase comes first

The existing Held qualification proves only that the default generator returns
a non-empty `ToolpathResult`. Figure 5 currently produces 2,289 operations, but
no committed evidence says whether that path covers the pocket, avoids gouging,
respects the 80-degree engagement cap, remains continuous, or avoids redundant
cuts.

Generator repair cannot be designed honestly before those failures are
measured. Likewise, a postprocessor must not turn an unaudited path into a file
that looks machine-ready. Phase 1 therefore creates one complete diagnostic
vertical slice and makes its release refusal explicit.

## Scope

Phase 1 uses exactly:

- `load_held_reference_case("figure5")` as the fixed input;
- `benchmarks.runner.generate_toolpath` as the only generator;
- the case's normalized 2 mm tool and 80-degree cap;
- `audit_toolpath_engagement` as the independent conservative engagement
  consumer;
- the exact sampled-position predicate through the existing exceedance
  measurement;
- `measure_quality` as the existing stock, coverage, motion, and program
  quality consumer; and
- a serial live command that writes one committed Markdown report.

The other three Held cases enter only after the Figure 5 consumer contract is
complete. No alternative engagement-controlled, radial, rho, ordered,
spiral-entry, or adaptive generator is compared in this phase.

## Non-goals

Phase 1 does not:

- change any Held reference geometry, normalization, projection, or cap;
- change generator behavior or select a replacement generator;
- change engagement predicates, certification, replay, or quality thresholds;
- claim that sampled absence of gouge, overload, or residual stock is a proof
  of global absence;
- emit Fanuc or any other NC program;
- estimate real machine cycle time;
- model controller look-ahead, machine kinematics, fixtures, holders, tool
  deflection, or collision; or
- claim parity or superiority over Held and Pfeiffer.

## Protected scientific boundary

Benchmark pressure must not alter the judge. The following remain unchanged:

- native `stock_2`, `engagement_2`, `continuous_tea_2`, and `audit_*` decision
  semantics;
- exact stock-depletion behavior;
- `compas_cgal.engagement` replay and certification behavior;
- `compas_cgal.stock`;
- native and Python engagement falsifier/replay tests;
- `benchmarks.quality.measure_quality` semantics;
- the twelve existing machinability requirements; and
- all four Held case documents and their reconstruction code.

The new code orchestrates these consumers and records their results. It does
not make a new geometric decision and introduces no new artifact-identity
mechanism.

## Evidence hierarchy

The report keeps four forms of evidence separate:

1. **Continuous engagement certification.** An operation is certified only
   when the protected audit closes its cap proof.
2. **Demonstrated engagement exceedance.** An operation is demonstrated over
   the cap only when at least one sampled cutter position returns the exact
   cap-exceeded predicate. This is a lower bound on violating operations.
3. **Unresolved engagement.** An engaging operation is unresolved when the
   continuous audit does not certify it and no sampled position demonstrates
   exceedance. Unresolved blocks release but is not reported as overload.
4. **Machining diagnostics.** Coverage and most gouge observations are sampled
   diagnostics. Redundant-cut detection is exact under the existing replay.
   The report retains those claim boundaries beside the values.

For every engaging operation, the three engagement dispositions are mutually
exclusive and exhaustive:

```text
certified
or demonstrated_exceeded
or unresolved
```

Non-engaging rapid, retract, and plunge operations are counted separately and
do not inflate the engagement denominator.

## Architecture

```text
HeldReferenceCase("figure5")
        |
        v
generate_toolpath(spec) -------------------- generation timing
        |
        +--> audit_toolpath_engagement ----- continuous audit timing
        |
        +--> exceedance_positions ---------- demonstrated witnesses
        |
        +--> measure_quality ---------------- quality timing
        |
        v
HeldFigure5Characterization.build(...)
        |
        +--> render Markdown evidence
        |
        +--> require_releasable(...)
                  |
                  +--> return normally only if every gate closes
                  +--> raise HeldPathNotReleasableError otherwise
```

The exact same generated `ToolpathResult` is supplied to every consumer. The
orchestrator never regenerates between consumers, so a disagreement cannot be
explained by comparing different paths. Each consumer owns fresh stock through
its existing interface.

## Components

### `benchmarks/held_machinability.py`

Owns orchestration-neutral domain records and the release decision:

- `Seconds = NewType("Seconds", float)`;
- `Degrees = NewType("Degrees", float)`;
- `OperationIndex = NewType("OperationIndex", int)`;
- `EngagementDispositionCounts`;
- `EngagementExceedanceWitness` with `Point2[WorldXY]` and measured degrees;
- `MachinabilityViolation` with criterion, measured value, required maximum,
  and evidence kind;
- `HeldFigure5Characterization.build(...)`; and
- `require_releasable(characterization) -> None`.

Factories validate finite unit values, non-negative counts, operation-index
bounds, mutually exclusive engagement dispositions, and exact case name. Raw
records store already typed values. Release refuses incomplete evidence as well
as measured failure.

The twelve requirements are represented in production without altering their
current values:

1. uncut fraction equals zero;
2. gouging motions equal zero;
3. unsafe rapids equal zero;
4. continuity breaks equal zero;
5. zero-length motions equal zero;
6. degenerate loops equal zero;
7. redundant operations equal zero;
8. sampled cap exceedances equal zero;
9. slotting motions equal zero;
10. maximum engagement step is at most the case cap;
11. maximum loop-radius step is at most two tool radii; and
12. tangent breaks equal zero.

The existing `test_quality.py` gate remains intact as an independent shadow
oracle. New contract tests prove the production evaluator applies the same
named requirements to synthetic `PathQuality` records without executing the
expensive generator matrix; Phase 1 does not edit the old reference assertions.

### `tools/held_machinability_qualification.py`

Owns the slow orchestration and report writing:

```python
def characterize_figure5(
    generator: Generator,
    engagement_auditor: EngagementAuditor,
    quality_measurer: QualityMeasurer,
    exceedance_locator: ExceedanceLocator,
    *,
    clock: Clock = time.perf_counter,
) -> HeldFigure5Characterization:
    ...
```

The injected interfaces make normal tests fast without adding alternate runtime
behavior. Production `main()` wires only the existing concrete consumers.
Unexpected consumer exceptions propagate; they are not converted into a
successful report or an empty result.

The command writes
`docs/benchmarks/held_figure5_machinability.md`. Characterization itself exits
successfully when complete even if the path is not releasable: its product is
truthful evidence. A separate `require_releasable` call is the downstream NC
release seam and raises `HeldPathNotReleasableError` for the expected current
failure.

### Explicit live oracle

`tests/benchmarks/held_figure5_machinability_live_oracle.py` runs only through a
named Pixi task. It is not part of ordinary pytest discovery. The task runs
serially because exact accumulated-stock depletion dominates runtime and
concurrent cases would provide worse evidence at higher cost.

The live oracle:

- loads Figure 5 through the strict case loader;
- generates exactly once;
- completes all three consumers;
- validates the engagement partition;
- validates the complete twelve-criterion evaluation;
- writes the report;
- confirms the report says `not releasable` when any criterion or unresolved
  operation remains; and
- records generation, engagement-audit, and quality timings separately.

There is no arbitrary timeout in the product contract. If the protected audit
cannot complete in the available execution window, Phase 1 reports a measured
performance blocker rather than switching to a weaker consumer.

## Release semantics

`require_releasable` returns normally only when all of these are true:

- the case is exactly Figure 5;
- all generated operations were classified;
- every engaging operation is continuously certified;
- demonstrated exceedance count is zero;
- unresolved count is zero;
- all twelve machinability requirements pass; and
- every required consumer completed on the same generated result.

Any failure raises `HeldPathNotReleasableError` containing the criterion names
and engagement disposition counts. It does not include a partial NC program.

## Report contract

The committed Markdown report contains:

- case, tool diameter, cap, primitive count, and projection count;
- total and engaging operation counts;
- generation, continuous-audit, and quality timings;
- certified, demonstrated-exceeded, and unresolved engagement counts;
- operation-indexed exceedance witnesses with world-XY millimetre positions;
- all twelve criterion values and pass/fail outcomes;
- the remaining report-only `PathQuality` quantities;
- the release verdict; and
- explicit sampled-versus-certified claim boundaries.

The report does not include G-code, cycle-time estimates, or comparative Held
performance claims.

## Error model

New named failures are independent:

- `InvalidHeldMachinabilityEvidenceError`: typed counts, units, or operation
  coverage are inconsistent;
- `IncompleteHeldMachinabilityEvidenceError`: a required consumer result is
  absent;
- `HeldPathNotReleasableError`: complete evidence exists but at least one
  engagement or machinability gate is open; and
- `UnexpectedHeldMachinabilityCaseError`: orchestration receives anything
  except the Figure 5 case.

Existing consumer failures propagate with their existing types. No broad
exception handler turns a broken audit into a characterization.

## Testing

Ordinary fast tests use injected deterministic consumers and synthetic typed
results to cover:

- exactly one generator invocation and object-identical result consumption;
- certified/exceeded/unresolved partitioning;
- non-engaging operation exclusion;
- witness operation-index and coordinate validation;
- all twelve production criteria;
- parity with every named requirement in the unchanged existing quality gate;
- refusal on one open criterion;
- refusal on one unresolved operation;
- refusal when evidence is incomplete;
- unexpected consumer exception propagation;
- Markdown escaping and deterministic row ordering; and
- strict type contracts for every public interface.

Repository gates are focused pytest with `-n auto`, `types-benchmarks`, Ruff,
strict MkDocs, and `git diff --check`. The explicit live task supplies the real
Figure 5 evidence.

## Phase boundary

Phase 1 ends when the real Figure 5 report is committed and independently
reviewed. Its measured failure signature becomes the sole input to the Phase 2
generator-policy design. No generator repair begins in the Phase 1 plan.

The later Fanuc phase remains downstream of `require_releasable`; it cannot
start emitting the Held path merely because the Phase 1 characterization
completed.
