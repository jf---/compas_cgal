# Held Figure 5 2D Path Characterization Design

> **status: approved** - approved for implementation planning on 2026-09-02;
> generator repair and G-code remain outside Phase 1.

## Goal

Characterize the current default generator on the fixed Held-Pfeiffer Figure 5
reference pocket through the repository's existing engagement and machining
quality consumers. The result is a durable failure report with
operation-attributed failure evidence, aggregate 2D machining diagnostics, and
a fail-closed verdict on whether this exact path may enter target-specific
postprocessor qualification. This phase does not repair generation, authorize
machining, qualify a controller, or emit G-code.

## Why this phase comes first

The existing Held qualification proves only that the default generator returns
a non-empty `ToolpathResult`. Figure 5 currently produces 2,289 operations, but
no committed evidence says whether that path covers the pocket, avoids gouging,
respects the 80-degree engagement cap, remains continuous, or avoids redundant
cuts.

Generator repair cannot be designed honestly before those failures are
measured. Likewise, a postprocessor must not turn an unaudited path into a file
that looks machine-ready. Phase 1 therefore creates one complete diagnostic
vertical slice and makes its post-qualification refusal explicit.

## Scope

Phase 1 uses exactly:

- `load_held_reference_case("figure5")` as the fixed input;
- `benchmarks.runner.generate_toolpath` as the only generator;
- the case's normalized 2 mm tool and 80-degree cap;
- `audit_toolpath_engagement` as the independent conservative engagement
  consumer;
- one canonical standard-density `PathSurvey` (45 intervals per motion),
  retaining the exact sampled-position predicate verdict and world position for
  each sample;
- the unchanged `measure_quality` reduction semantics, now fed by that canonical
  survey rather than reconstructing it; and
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
  deflection, axial entry suitability, centre-cutting capability, plunge load,
  stock thickness, stock-top or cut-depth Z, work coordinate systems, safe
  planes, or collision; or
- claim parity or superiority over Held and Pfeiffer.

All millimetres in this phase are benchmark-normalized coordinates reconstructed
from the paper, not an authorized manufacturing scale or setup. Passing Phase 1
is not a machine-release decision. Global absence of residual stock and gouging,
controller semantics, machine kinematics, fixtures, holders, and collisions all
require downstream evidence that this phase does not own.

## Protected scientific boundary

Benchmark pressure must not alter the judge. The following remain unchanged:

- native `stock_2`, `engagement_2`, `continuous_tea_2`, and `audit_*` decision
  semantics;
- exact stock-depletion behavior;
- `compas_cgal.engagement` replay and certification behavior;
- `compas_cgal.stock`;
- native and Python engagement falsifier/replay tests;
- `benchmarks.quality.measure_quality` semantics;
- the twelve existing path-quality requirements; and
- all four Held case documents and their reconstruction code.

The new code orchestrates these consumers and records their results. It does
not make a new geometric decision or add an unrelated provenance mechanism.

Phase 1 exercises `audit_toolpath_engagement`, not the separate native audit
replay protocol. Its segment certifier and guarded fixed-density arc/circle
certifier are conservative over their motions, but a Phase 1 result establishes
no native-audit compatibility claim. Native audit integration remains a later,
independently gated concern.

## Evidence hierarchy

The report keeps four forms of evidence separate:

1. **Guarded replay certification.** A measured cut operation is certified only
   when `audit_toolpath_engagement` closes its conservative cap proof. For arcs
   and circles this is the protected guarded fixed-density replay, not the
   native audit protocol.
2. **Demonstrated engagement exceedance.** An operation is demonstrated over
   the cap only when at least one sampled cutter position returns the exact
   cap-exceeded predicate. This is a lower bound on violating operations.
3. **Unresolved engagement.** A measured cut operation is unresolved when the
   guarded replay does not certify it and no sampled position demonstrates
   exceedance. Unresolved blocks postprocessor qualification but is not reported
   as overload.
4. **Machining diagnostics.** Coverage and most gouge observations are sampled
   diagnostics. Redundant-cut detection is exact under the existing replay.
   The report retains those claim boundaries beside the values.

Every gate observation carries one explicit evidence kind, represented by a
closed `Literal` rather than an enum:

- `guarded_replay_certificate`;
- `sampled_exact_predicate`;
- `exact_depletion_replay`;
- `sampled_diagnostic`;
- `derived_geometry`; or
- `tolerance_diagnostic`.

A passing sampled diagnostic says that this measurement found no failure. It
never upgrades into a proof of global absence and therefore never authorizes
machine release.

The criterion-to-evidence mapping is fixed in production: uncut stock and
gouging are sampled diagnostics; cap exceedance and slotting are sampled exact
predicates; redundant cutting is exact depletion replay; zero length,
degenerate loops, and loop-radius steps are derived geometry; unsafe rapids,
continuity, and tangent breaks are tolerance diagnostics; engagement steps are
sampled diagnostics. A missing or wrong kind is invalid evidence.

For every measured cut operation, the three engagement dispositions are
mutually exclusive and exhaustive:

```text
certified
or demonstrated_exceeded
or unresolved
```

Unmeasured rapid, retract, plunge, and clearance operations are counted
separately and do not inflate this denominator. The report also gives a separate
sampled material-contact count from `MotionQuality.is_engaged`; it does not call
every measured cut an actually engaging cut.

The partition is constructed from operation indices, not reporting doubles:

```text
tea_audited_lateral = audit operations with stations > 0
certified = tea_audited_lateral with cap_certified
witness_rows = canonical survey samples whose exact cap predicate fired
witnessed_operation_indices = unique operation indices from witness_rows
demonstrated_exceeded = witnessed_operation_indices
unresolved = tea_audited_lateral - certified - demonstrated_exceeded
```

If any witnessed operation is certified, construction raises
`ContradictoryEngagementEvidenceError`. The system does not choose the more
convenient consumer. Multiple witness rows on one operation still produce one
disposition member. Audit operation order, count, and operation kind must also
match the canonical operation snapshot before any disposition is accepted.
The `tea_audited_lateral` indices must equal the canonical `PathSurvey.motion`
indices, and the audit-excluded indices must partition every remaining source
operation. A plunge may remove a disk of material, but remains explicitly
TEA-audit-excluded because axial entry requires a different governing model.

## Architecture

```text
HeldReferenceCase("figure5")
        |
        v
generate_toolpath(spec) -------------------- generation timing
        |
        +--> snapshot complete operation stream
        |
        +--> audit_toolpath_engagement ----- guarded replay timing
        |
        +--> survey_path -------------------- one canonical survey timing
        |
        +--> reduce_quality_evidence -------- coverage + reduction timing
        |
        +--> HeldFigure5Characterization.build(...)
        |
        +--> render Markdown evidence
        |
        +--> require_post_qualification_candidate(...)
                  |
                  +--> return immutable HeldPostQualificationCandidate
                  +--> raise HeldPathNotEligibleForPostQualificationError
```

Generation is snapshotted once into immutable, unit- and frame-typed primitive
values. The same generated `ToolpathResult` is supplied to both path-replay
consumers, and its complete stream is compared behaviorally with the snapshot
before and after every call. Count, order, operation kind, coordinates, Z,
radius, angular extent, travel direction, and generator chain identity are all
covered. Consumer mutation therefore fails loudly instead of creating mixed
evidence. Each replay consumer owns fresh stock through its existing interface.

`survey_path` builds exactly one `PathSurvey`; `reduce_quality_evidence` feeds it
to one canonical set of criterion reducers. Each reducer returns its aggregate
and the operation index or adjacent-operation pair that produced it.
`PathQuality` and the Held assessment consume those same reducer outputs; there
is no second survey and no second implementation of a decision rule. The
existing `measure_quality(spec, result)` remains the compatibility entry point
and delegates to these two stages. Residual stock is spatial evidence and is
not falsely attributed to one motion.

## Components

### `benchmarks/units.py`

Owns the small shared benchmark observation-unit vocabulary. It defines
`Degrees`, `Seconds`, `UnitFraction`, `MotionCount`, `ToolRadiusMultiple`, and
`OperationIndex`. Geometry reuses the canonical `Millimetre`, `Radian`,
`Point2[WorldXY]`, `Point3[WorldXYZ]`, and `Direction3[WorldXYZ]` types from
`compas_cgal.adaptive.units`; no parallel point or length vocabulary is added.
The report contract, rather than a duplicate physical unit, records that these
millimetres are benchmark-normalized.

### `benchmarks/survey.py`

Extends `EngagementSample` additively with `position` and the exact
`cap_exceeded` Boolean already returned by `_stock_2.engagement_at`. No caller
reconstructs an exact predicate verdict by comparing reported floating-point
degrees. Existing sampling density and quality decisions remain unchanged.

### `benchmarks/held_path_snapshot.py`

Owns the immutable, behavior-complete `HeldOperationSnapshot` primitive records,
snapshot construction, and pre/post-consumer equality assertion. Snapshot
is a closed union of `HeldLineSnapshot`, `HeldArcSnapshot`, and
`HeldCircleSnapshot`. Its variants preserve 3D endpoints or centre; frame x/y
axes; radius; arc start/end angles; operation role; clockwise travel;
generator `path_index`; and optional start/end tangents. `.build(...)` rejects
non-finite coordinates, non-positive radii, non-unit/non-orthogonal frame axes,
invalid angle ranges, and malformed tangent vectors. Those fields use the shared
benchmark length units plus typed angles and directions, so a circle seam/frame
or Z mutation is observable. No mutable COMPAS geometry escapes into the later
candidate.

### `benchmarks/quality_observations.py`

Owns the canonical reducers for the twelve existing requirements, their typed
criterion records, and their operation attribution:

- unit-specific validated criterion types for fractions, counts, degrees, and
  tool-radius multiples, each carrying a closed criterion name and its required
  `EvidenceKind`; and
- `PathQualityAssessment` with the twelve typed observations and
  operation-indexed attributions.

Distinct criterion types, rather than one runtime-erased generic, prevent a
fraction from being constructed as degrees, a count, or a tool-radius multiple.
Factories validate finite unit values, closed unit fractions, non-negative
counts, operation-index bounds, closed criterion names, and required evidence
kinds. Raw records store already typed values.

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

Each reducer returns both the typed aggregate and its attribution.
`benchmarks.quality` delegates its corresponding `PathQuality` fields to these
reducers, so the new Held report and the old quality gate share one decision
call graph.

Migration is additive and has two explicit checkpoints. First, add detailed
finding reducers and compare them against the untouched private functions over
the complete synthetic suite and all twelve current red gate cells. After that
parity is reviewed and removal is explicitly approved, make `measure_quality`
delegate, remove the superseded private reducer bodies, and retain direct
criterion name/value tests plus the unchanged twelve-red matrix. Public
signatures, thresholds, parametrization, and expected red membership do not
change.

Motion attribution reports:

- every sampled engagement witness;
- gouging, zero-length, redundant, degenerate-loop, and slotting operation
  indices;
- unsafe rapid indices;
- both operation indices at continuity, tangent, engagement-step, and
  loop-radius-step failures; and
- residual-stock fraction as spatial aggregate evidence without an invented
  causal operation.

Attribution reducers consume `PocketSpec`, the immutable operation snapshot,
`PathSurvey`, and coverage evidence as required. Generator chain identity is
needed for loop-run boundaries and is never guessed from `PathSurvey` alone.
Contract tests require every aggregate count and extremum to be the exact value
used to construct the corresponding `PathQuality` field.

### `benchmarks/held_path_evidence.py`

Owns `EngagementDispositionCounts`, `EngagementExceedanceWitness` with
`Point2[WorldXY]`, and `HeldFigure5Characterization.build(...)`.
Its factory validates exact case name, complete operation classification,
operation bounds, finite witness coordinates, mutually exclusive dispositions,
and cross-consumer consistency. Witness positions are not required to lie inside
the pocket or legal cutter-centre domain: an outside-domain sample is valid gouge
evidence. A characterization exists only after every consumer returns complete
typed evidence.

The factory copies only immutable derived values into the characterization. It
retains no `EngagementReport`, `PathSurvey`, `Stock`, `ToolpathResult`, COMPAS
geometry object, NumPy array, or mutable collection supplied by a consumer.

### `benchmarks/held_post_qualification.py`

Owns frozen, `init=False` `HeldPostQualificationCandidate`, the pure
`post_qualification_failures(...) -> tuple[str, ...]` decision function, and
`require_post_qualification_candidate(...)`. The candidate's only public
construction path is `HeldPostQualificationCandidate.build(characterization)`,
which consumes that one decision function and raises
`HeldPathNotEligibleForPostQualificationError` if one is open. A successful
candidate owns both the closed characterization and its exact immutable
operation snapshot. The retained characterization carries the exact Figure 5
case name, normalized tool diameter, TEA cap, and quality evidence required by
the future post consumer; callers cannot supply a free-standing snapshot.
Measured failure still yields a complete characterization, but never a
candidate.
`require_post_qualification_candidate` is the functional spelling of that same
factory and delegates to it; it is not a second gate implementation.
The Markdown renderer consumes `post_qualification_failures(...)` to state the
same verdict without constructing or catching a candidate.

All new invariant-bearing public records in this design are frozen with direct
initialization disabled and expose `.build(...)` as their validated construction
path. Their factories own the named error model; raw mutable field bags cannot
bypass invariants.

### `benchmarks/held_path_characterize.py`

Owns only the slow orchestration:

```python
class Generator(Protocol):
    def __call__(self, spec: PocketSpec) -> ToolpathResult: ...

class EngagementAuditor(Protocol):
    def __call__(
        self, spec: PocketSpec, result: ToolpathResult
    ) -> EngagementReport: ...

class PathSurveyor(Protocol):
    def __call__(
        self, spec: PocketSpec, result: ToolpathResult
    ) -> PathSurvey: ...

class QualityEvidenceReducer(Protocol):
    def __call__(
        self,
        spec: PocketSpec,
        snapshot: tuple[HeldOperationSnapshot, ...],
        survey: PathSurvey,
    ) -> QualityEvidence: ...

def characterize_figure5(
    generator: Generator,
    engagement_auditor: EngagementAuditor,
    path_surveyor: PathSurveyor,
    quality_evidence_reducer: QualityEvidenceReducer,
    *,
    phase_observer: Callable[[CharacterizationPhase], None],
    clock: Clock = _monotonic_seconds,
) -> HeldFigure5Characterization:
    ...
```

The injected interfaces make normal tests fast without adding alternate runtime
behavior. The CLI wires only the existing concrete consumers.
Unexpected consumer exceptions propagate; they are not converted into a
successful report or an empty result.

`QualityEvidence` contains the existing `PathQuality`, the typed
`PathQualityAssessment`, and coverage evidence derived from the one survey.
`benchmarks/held_consumer_adapters.py` owns the only signature adapter:
`audit_figure5_engagement(spec, result)` maps `PocketSpec` into the existing
`audit_toolpath_engagement(polygon, result, tool_diameter, tea_cap, holes)` call.
It changes no values or decisions. `generate_toolpath`, `survey_path`, and
`reduce_quality_evidence` already satisfy their respective protocols directly.

`CharacterizationPhase` is the closed literal vocabulary `generation`,
`guarded_audit`, `survey`, and `quality_reduction`; the required observer is
called immediately before each stage so a supervised live run exposes where it
stopped. `Clock = Callable[[], Seconds]`; `_monotonic_seconds()` is the typed
wrapper around `time.perf_counter`. Characterization contains no optional
consumer fields: a consumer returns its complete typed evidence or raises. There
is no partially constructed characterization state.

### `benchmarks/held_path_report.py`

Owns frozen, validated `HeldPathReportContext` with an aware UTC generation
instant, exact Pixi invocation string, and generator policy name. It renders
deterministic Markdown from `(HeldFigure5Characterization,
HeldPathReportContext)`. It does not run geometry, decide criteria, read a clock,
or write files; it obtains the post-entry verdict from
`post_qualification_failures(...)` rather than duplicating its conditions.
It rejects non-exact characterization and context objects before projection,
labels every timing and report-only quantity with its semantic unit, and
escapes Markdown and HTML-sensitive dynamic table content while retaining
renderer-owned `<br>` newline markers.

### `tools/held_path_characterization.py`

Owns only CLI argument handling, production wiring, and the write to
`docs/benchmarks/held_figure5_2d_path.md`.

The report writer requires the actual Pixi invocation and a keyword-only
`Callable[[CharacterizationPhase], None]` observer as arguments. The normal
CLI supplies `pixi run held-figure5-characterize`; the explicit live oracle
supplies `pixi run held-figure5-characterize-live`. It constructs
`HeldPathReportContext` using its injected UTC wall clock. This wall clock is
distinct from the monotonic timing clock and is injected in tests as `UtcClock =
Callable[[], datetime]`; `_utc_now()` returns an aware UTC `datetime`.
The normal CLI passes `_print_phase` explicitly to the writer and orchestration
and emits all four phase
markers in order. The live oracle instead passes its ledger-writing observer;
there is no omitted or no-op production observer path.

Characterization itself exits successfully when complete even if the path is
not eligible: its product is truthful evidence. A separate
`require_post_qualification_candidate` call is the downstream seam and raises
`HeldPathNotEligibleForPostQualificationError` for the expected current failure.

On success, the returned candidate owns the immutable operation snapshot that
was characterized. A future postprocessor accepts that candidate, never an
unrelated or regenerated `ToolpathResult`. The target-specific Fanuc phase must
still prove a named control/option set, coordinate quantization, supported arc
forms, modal behavior, and semantic round trip before it may claim translation
eligibility.

### Explicit live oracle

`tests/benchmarks/held_figure5_path_live_oracle.py` runs only through a
named Pixi task. It is not part of ordinary pytest discovery. The task runs
serially because exact accumulated-stock depletion dominates runtime and
concurrent cases would provide worse evidence at higher cost.

The live oracle:

- loads Figure 5 through the strict case loader;
- generates exactly once;
- completes the guarded engagement replay, one canonical survey, and one
  quality-evidence reduction;
- validates the engagement partition;
- validates the complete twelve-criterion evaluation;
- writes the report;
- confirms the report says `not eligible for postprocessor qualification` when
  any criterion or unresolved operation remains; and
- records generation, guarded-audit, survey, and coverage-plus-reduction timings
  separately.

Before generation, it writes the run start to
`docs/superpowers/state/held-figure5-live-run.md`; the phase observer updates the
same ledger before each stage. A report is current evidence only when its UTC
generation instant postdates that run start. On operator-budget exhaustion, the
ledger's last phase remains durable and any older report is explicitly prior,
stale evidence rather than proof of this run.

There is no arbitrary timeout in the product contract. If the protected audit
cannot complete in the available execution window, Phase 1 reports a measured
performance blocker rather than switching to a weaker consumer.

The execution plan grants the single Figure 5 live run a 30-minute operator
budget. This is an attention/compute bound, not a geometry tolerance and not a
runtime claim. Exceeding it writes no characterization report and leaves Phase 1
blocked. Elapsed wall time and active phase go in the SDD progress ledger or
operator log because the current consumer exposes no trustworthy partial
operation count. Budget exhaustion does not authorize sampling reduction,
consumer removal, or a fallback result.

### Lightweight survey centre-domain predicate

The first live oracle generated 2,289 operations in 0.063150 seconds, completed
the guarded audit, then exceeded the 30-minute operator budget in `survey`
without producing a report. Bounded diagnosis isolated the stall before replay:
the survey constructed the full `ReachableDomain2` state from the 65-vertex
Figure 5 polygon even though it consumed only centre-domain point membership.
That constructor exceeded 50 seconds for the full boundary and for surveys
restricted to 3 and 25 source operations; one first-motion engagement query took
0.000640 seconds. Downstream replay cost is therefore unmeasured, not presumed
fast.

The survey instead builds one lightweight exact cutter-centre-domain predicate
through `CutterCentreDomain2.build(...)`. The public factory owns the same
canonical validated reach input; direct Python construction is disabled, and the
result answers only `contains(x, y)`. A finite binary64 query is injected exactly
once and is legal
exactly when all of these conditions hold:

- it lies inside or on the canonical outer polygon;
- it lies outside every canonical hole interior; and
- its exact squared distance to every outer and hole boundary segment is greater
  than or equal to the exact squared tool radius.

Exact tangency is accepted. The adjacent representable point on the illegal side
is rejected. The implementation uses exact kernel predicates and squared-distance
comparisons without tolerance, snapping, approximation, or fallback. Canonical
reach-input validation retains the existing named malformed-ring, non-finite,
and invalid-radius error boundary. One shared exact polygon-with-holes builder
also preserves the existing relationship validation for holes outside the outer
ring, intersecting rings, and otherwise invalid design topology.

This specialized predicate owns pointwise legal-centre semantics only. It does
not certify that the complete eroded centre domain is nonempty or connected;
those global construction results and their named failures remain exclusively
owned by `ReachableDomain2`. Avoiding that arrangement is the measured purpose
of this repair, so the predicate does not claim full constructor-error
equivalence.

`ReachableDomain2` remains unchanged for consumers that require reachable
material, residual, or certificate products. Only `benchmarks.survey.survey_path`
uses the lightweight predicate. The repair does not change generator behavior,
survey density, quality thresholds, engagement decisions, stock depletion,
evidence vocabulary, report semantics, or the unchanged Task 11 live-oracle
gate.

## Post-qualification entry semantics

`require_post_qualification_candidate` returns a
`HeldPostQualificationCandidate` owning the closed characterization and its
immutable operation snapshot only when all of these are true:

- the case is exactly Figure 5;
- every consumer observed the unchanged complete operation snapshot;
- all generated operations were classified;
- every TEA-audited lateral operation is guarded-replay certified;
- demonstrated exceedance count is zero;
- unresolved count is zero;
- all twelve 2D criteria close under their declared evidence semantics; and
- the guarded replay and single-survey quality consumer completed.

Any failure raises `HeldPathNotEligibleForPostQualificationError` containing the
criterion names and engagement disposition counts. It returns neither a
candidate nor a partial NC program.

Phase 1 deliberately defines no `require_machine_releasable` function. Sampled
negative evidence cannot implement that contract. The later machine-release
design must add globally sufficient stock/gouge evidence plus controller,
machine, setup, fixture, holder, and collision validation before it may use that
name.

The later target-specific postprocessor must accept only
`HeldPostQualificationCandidate`. Its provisional output remains outside the
ordinary machine-program export location and uses a non-`.nc` unreleased
artifact form until a distinct machine-release consumer succeeds. A header
comment alone is not quarantine.

## Report contract

The committed Markdown report contains:

- case, tool diameter, cap, primitive count, and projection count;
- total, TEA-audited lateral, TEA-audit-excluded, and sampled material-contact
  operation counts;
- generation, guarded-audit, survey, and coverage-plus-reduction timings;
- certified, demonstrated-exceeded, and unresolved engagement counts;
- operation-indexed exceedance witnesses with world-XY millimetre positions;
- all twelve criterion values, evidence kinds, and outcomes: sampled negatives
  use `no_failure_observed`, sampled positives use `failure_observed`, exact
  depletion/derived checks use `criterion_satisfied` or `criterion_violated`,
  and tolerance diagnostics use `within_declared_tolerance` or
  `outside_declared_tolerance`;
- the remaining report-only `PathQuality` quantities;
- the post-qualification entry verdict;
- explicit sampled-versus-certified claim boundaries.

The report does not include G-code, cycle-time estimates, or comparative Held
performance claims.

Tolerance outcomes mean only that the declared path-geometry predicate evaluated
inside or outside its named tolerance. They neither erase that tolerance nor
establish a global stock, controller, or machine proof.

The report begins with a prominent "historical 2D benchmark characterization -
not a manufacturing release" banner, the UTC generation time, the exact Pixi
command, and the generator policy name. Phase 2 must regenerate the live report
from the current worktree before using it; a committed older report is evidence
of its recorded run, never standing authorization.

## Error model

New named failures are independent:

- `InvalidHeldOperationSnapshotError`: a primitive contains non-finite,
  non-positive, non-unit, non-orthogonal, or otherwise malformed geometric
  fields;
- `InvalidHeldPathReportContextError`: report invocation metadata is empty,
  non-UTC, or otherwise malformed;
- `InvalidHeldPathEvidenceError`: typed counts, units, or operation
  coverage are inconsistent;
- `ContradictoryEngagementEvidenceError`: a certified operation also has an
  exact sampled exceedance witness, or consumer operation sequences disagree;
- `ContradictoryPathQualityEvidenceError`: operation attributions do not
  reduce to the corresponding aggregate quality values;
- `MutatedHeldToolpathError`: the generated operation stream differs from its
  canonical snapshot before or after a consumer call;
- `HeldPathNotEligibleForPostQualificationError`: complete evidence exists but
  at least one engagement or 2D criterion is open; and
- `UnexpectedHeldPathCaseError`: orchestration receives anything
  except the Figure 5 case.

Existing consumer failures propagate with their existing types. No broad
exception handler turns a broken audit into a characterization.

## Testing

Ordinary fast tests use injected deterministic consumers and synthetic typed
results to cover:

- exactly one generator invocation, object-identical consumer input, and
  snapshot-bound candidate output;
- certified/exceeded/unresolved partitioning;
- certified-plus-exceeded contradiction rejection;
- audit/source operation count, order, and kind rejection;
- equality of TEA-audited lateral and survey-motion indices;
- complete TEA-audit-excluded operation classification, including
  material-removing plunges;
- witness operation-index and coordinate validation;
- multiple witness rows reducing to one demonstrated operation;
- all twelve production criteria;
- parity with every named requirement in the unchanged existing quality gate;
- unit-safe fraction, count, degree, and tool-radius-multiple criteria;
- operation attribution reducing exactly to aggregate quality;
- pre- and post-consumer mutation rejection across every snapshotted primitive
  field;
- malformed snapshot rejection and disabled direct initialization for every
  invariant-bearing record;
- candidate retention of the closed characterization, operation snapshot,
  normalized tool, and cap;
- refusal on one open criterion;
- refusal on one unresolved operation;
- unexpected consumer exception propagation;
- Markdown escaping and deterministic row ordering; and
- strict type contracts for every public interface.

Repository gates are focused pytest with `-n auto`, `types-benchmarks`, Ruff,
strict MkDocs, and `git diff --check`. The explicit live task supplies the real
Figure 5 evidence.

## Phase boundary

Phase 1 ends only when the real Figure 5 report is committed and independently
reviewed. Its measured failure signature becomes the sole input to the Phase 2
generator-policy design. Budget exhaustion prevents Phase 1 completion and
therefore blocks Phase 2 generator work; its elapsed timing may guide bounded
audit-performance work and a resumed Phase 1 run, but it is not path evidence.
No generator repair begins in the Phase 1 plan.

The later Fanuc phase remains downstream of
`require_post_qualification_candidate`; it cannot claim translation eligibility
or emit a machine-loadable Held program merely because Phase 1 completed.
