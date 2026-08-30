# Wave-1 Task-7 review continuation kickoff prompt

Paste everything below the rule as the FIRST message of a fresh session.

---

I'm starting **Wave 1 continuation — Task-7 semantic review, then Tasks 8–11**
of the `compas_cgal` coherence programme (worktree
`/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence`, branch
`codex/sdd-coherence`). Commit `f3a6d2e8ac2f6625d7bb2bb420054b04f2365bb9`
is the last hard technical checkpoint. Its whole-branch review passed, but it
is **not** a final human-acceptance checkpoint: after that review, the user
rejected the presentation/inference around MC-013 and asked to see the actual
paths. Task 7 is therefore under semantic review; Tasks 8–11 have not started.

## What I want you to do

Use `plan-status` first and keep exactly one active step: resolve the MC-013
review boundary without rerunning or rewriting the authenticated one-shot
Figure-6 artifact. Use `superpowers:receiving-code-review` to evaluate the
rejection, and `superpowers:systematic-debugging` if any asserted path fact
fails to reproduce. Present root cause and the smallest technically sound
options before changing repository files, as required by the session-injected
AGENTS instructions.

Do **not** start Task 8 while Task 7 remains under review. Once the MC-013
ruling is explicit, implement only the approved correction with TDD, rerun the
Task-7 acceptance boundary, obtain independent review, and only then execute
Tasks 8, 9, 10, and 11 from
`docs/superpowers/plans/2026-08-28-coherence-wave1.md` in order. Use
`superpowers:subagent-driven-development` where tasks are independent enough
for fresh implementer/reviewer agents; agents use the same model and effort as
the parent.

## Read first, in dependency order

1. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/CLAUDE.md`
   — binding exact-kernel, documentation, Pixi, and development-stage policy.
2. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/docs/superpowers/plans/2026-08-28-coherence-wave1.md`
   — authoritative plan. Read Task 7 and Tasks 8–11 plus the plan header. The
   header says Tasks 2–7 landed; do not interpret that as human acceptance of
   MC-013 after the later review objection.
3. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/.superpowers/sdd/2026-08-28-coherence-wave1/progress.md`
   — execution/resume record and the last whole-review evidence. It is ignored
   working state, not a substitute for committed evidence.
4. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/docs/measurement_claims.md`
   — committed 14/14 ledger; MC-013 is row 013 and currently says `corrected`.
5. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/benchmarks/mathsm.py`
   — the disputed source wording: a fixed sweep is currently promoted into
   the broad statement `THE RELATION IS NOT MONOTONE`.
6. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/benchmarks/pathmetrics.py`
   and
   `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/benchmarks/figure6.py`
   — reference measurement/generation semantics. The audit excludes each
   chain's entry operation and reports a maximum per operation, not the exact
   cutter-centre station where the maximum occurred.
7. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/tools/measurement_claim_benchmark_schema.py`,
   `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/tools/measurement_claim_benchmark_validation.py`,
   `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/tools/measurement_claim_artifact_validation.py`,
   and
   `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/tools/measurement_claim_ledger.py`
   — the authenticated benchmark payload, artifact validator, renderer, and
   exact-two-artifact ledger consumer. Read before proposing any correction;
   the committed payload reason and disposition are content-addressed.
8. `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/docs/superpowers/state/wave1-gate-analysis.md`
   and
   `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence/docs/superpowers/state/wave1-interim-review.md`
   — accepted deviation rulings. The claims framework remains frozen through
   Task 11; do not generalize it or retire producer code during Wave 1.
9. `/Users/jelle/Code/CADCAM/compas_cgal_prs-visualizations/mc013-path-diagnostic/mc013-peaks-detail.png`
   and
   `/Users/jelle/Code/CADCAM/compas_cgal_prs-visualizations/mc013-path-diagnostic/mc013-paths.png`
   — review-only path diagnostics generated after the user objection. They are
   durable but outside Git and are **not** authenticated measurement artifacts.
10. `/Users/jelle/.codex/memories/MEMORY.md` — historical T9/reconciliation
    context only. Live files and Git state outrank it. No project-local
    `feedback_*.md` override exists; session-injected AGENTS instructions and
    the worktree `CLAUDE.md` are binding.

## Boundary contract (= your input)

The accepted Task-7 system consumes one generator artifact and one benchmark
artifact, renders their exact evidence, and requires an exact 14-row union:

```python
class MC013SelectedValuesPayload(TypedDict):
    fine_spacing: ToolDiameters
    fine_max_tea_after_entry: Degrees
    comparison_spacing: ToolDiameters
    comparison_max_tea_after_entry: Degrees
    angle_unit: Literal["degree"]
    spacing_unit: Literal["tool-diameter"]


class MC013ClaimPayload(TypedDict):
    claim_id: Literal["MC-013"]
    source: Literal["benchmarks/mathsm.py:47"]
    disposition: Literal["corrected"]
    source_commit: GitObjectId
    history_commit: GitObjectId
    reason: str
    missing_inputs: List[MissingBenchmarkInput]
    selected_values: MC013SelectedValuesPayload


def validate_figure6_payload(
    figure6_payload: object,
    figure6_markdown: bytes,
) -> tuple[str, Degrees, Degrees]: ...


def validate_benchmark_payload(
    payload: object,
    *,
    figure6_payload: object,
    figure6_markdown: bytes,
) -> BenchmarkClaimPayload: ...


def validate_ledger_evidence(
    ledger: pathlib.Path,
    artifact_directories: tuple[pathlib.Path, pathlib.Path],
) -> None: ...
```

The exact accepted artifact pair is:

```text
benchmarks/measurement_claim_results/2026-08-29-3f62304932dd-generator-033fc8b07e19
benchmarks/measurement_claim_results/2026-08-29-f65154fb2075-benchmark-9e0c635f36bd
```

The benchmark artifact is authenticated to source commit
`f65154fb20755f33065e919209adfc928ac6490e`; it contains exactly four files:
`benchmark-claims.json`, `figure6.json`, `figure6.md`, and `stamp.json`. It was
run once. Do not edit it, restamp it, move it, or rerun it merely to change an
adjudication. A semantic correction must preserve the raw artifact and make
any new authority explicit; never mutate content-addressed evidence in place.

### Epistemic state

| Claim or failure | Classification | Evidence | Content identity | Consequence |
| --- | --- | --- | --- | --- |
| The committed benchmark artifact is structurally valid and its raw values are authentic. | verified invariant | Task-7 final whole review; artifact validator; joint ledger validator | artifact source `f65154f`; payload SHA-256 `68230f5a43557428fd03c63b43465e7c30247dfac9295a22c89550d15e8a64ed` | Preserve artifact bytes and identities. |
| In that one `rect_20x12`, tool-diameter-2 run, spacing `0.025` reported `131.135610799756°` and spacing `0.1` reported `98.72798380871845°` after entry. | verified invariant for this artifact only | `figure6.json`; fresh exact rerun reproduced both floats | source `f65154f`; semantic command/config embedded in artifact | May be stated only as a one-run observation unless stronger evidence is added. |
| Both reported maxima came from path 1's upper-right branch: op 423 at circle radius `0.039176 mm`, and op 143 at radius `0.184047 mm`. | verified diagnostic observation, not artifact field | fresh live audit plus `mc013-peaks-detail.png` | live HEAD `f3a6d2e`; detail PNG SHA-256 `9c81f481e110d05686a309ab70b5dda05bf256f2e3366ab22e98ff8ffcdaee29` | Useful for review; regenerate and authenticate before promoting into a claim. |
| The exact cutter-centre station of either maximum is known. | false premise | `OperationEngagement` retains `op_index` and per-operation maximum, not the peak station | identity-independent API fact at live HEAD | Do not label a point as the peak station or infer contact geometry from one. Highlight the whole operation only. |
| `THE RELATION IS NOT MONOTONE` is established as a general spacing law. | unresolved/overbroad hypothesis | one fixed 12-value sweep plus local path diagnostic | only the `f65154f` configuration | Do not present it as a general theorem. Human acceptance remains unresolved. |
| `selected[0.025] > selected[0.1]` proves the raw artifact was generated as configured. | verified content-specific validator fact | `validate_figure6_payload` | benchmark artifact identity above | It can authenticate those two values without proving a general causal relationship. |
| Task 7 has an implementation failure. | not established | all committed Task-7 technical gates and independent reviews passed | HEAD `f3a6d2e` | Treat the open item as semantic review first, not an invitation to rewrite production geometry. |
| The suite's 11 red tests are regressions to fix during this continuation. | false oracle premise | `docs/red_manifest.json`; two-way red-manifest reconciliation | identity-independent manifest contract at HEAD | Six product reds, one tangency red, and four superseded retrace reds are deliberate. Any twelfth red is a defect. |
| Tasks 8–11 are incomplete. | confirmed implementation/documentation gap | unchecked plan steps and clean live tree | HEAD `f3a6d2e` | Execute only after Task-7 review closes. |

## Deliverables

First, close the Task-7 review boundary:

- A written technical ruling in chat distinguishing the authenticated one-run
  observation from the unsupported general inference. Do not ask the user to
  repeat the objection.
- If a repository correction is approved, the smallest coherent TDD change to
  the existing Task-7 authority. Candidate files are
  `benchmarks/mathsm.py`, `tools/measurement_claim_benchmark_schema.py`,
  `tools/measurement_claim_benchmark_validation.py`,
  `tools/measurement_claim_ledger.py`, their focused tests, and
  `docs/measurement_claims.md`; touch only files required by the chosen
  authority model. Do not mutate the existing artifact.
- Fresh Task-7 consumer-boundary evidence and independent review. Keep the
  plan header at `Landed: 2, 3, 4, 5, 6, 7` only if the revised state still
  satisfies zero-pending, authenticated 14/14 closure.

Then complete the authoritative plan:

- Task 8: add axis-parallel best-case qualifiers to living performance pages.
  The progress ledger's accepted ruling is semantic, not a stale six-file
  count: modify the pages that actually make performance claims;
  `docs/oblique_edge_cost.md` is already the source warning.
- Task 9: create `docs/radius_ladder_ablation.md` with the exact historical
  table and provenance warning; modify `mkdocs.yml` only after its contention
  check.
- Task 10: create `tests/benchmarks/test_gate_attribution.py`. Follow the
  progress-ledger ruling: compare canonical defining geometry, suppress only
  the expected warning, and close C3 only when cap 60 is wired into the real
  three-pocket by two-generator gate. A no-divergence result is a STOP finding.
- Task 11: run the full Wave-1 oracle, close the backlog items, and change the
  plan header to `landed` only when every exit gate is satisfied.

## Exit gates

From plan Task 11, verbatim:

> `pixi run red-manifest && pixi run plan-headers && pixi run -e docs docs && pixi run lint && pixi run mypy --strict tools` — all green/exit 0.

Also required before Task 8 starts:

- The benchmark artifact remains byte-identical and validates.
- The joint generator/benchmark ledger consumer validates the exact ordered
  MC-001 through MC-014 union with zero pending rows, or the plan/header is
  explicitly reopened if the ruling makes 14/14 impossible.
- Focused Task-7 tests, Ruff, Python-3.9 AST parsing, and strict mypy are green.
- Independent review reports zero Critical, Important, or Minor findings.

If any gate slips, first write
`docs/superpowers/state/wave1-task7-review-gate-analysis.md` naming the failed
criterion, measured evidence, and smallest options. Do not weaken a validator,
rewrite a reference test, rerun the artifact, add a fallback path, or move on
to Task 8 to hide the slip.

## Decoupling rule

This branch remains independent of the auditor-convergence fleet. Work only in
`/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence`; never run
Git from the main checkout or auditor worktree. Never mutate
`jf/toolpath-redesign@73d5372`,
`codex/exact-certified-adaptive-phase1-t9-zero-guide@073a0f7`, or
`codex/auditor-convergence-sdd@fa59120`.

Task-7 semantic review is decoupled from a new Figure-6 run: validate the
committed raw artifact, inspect its source wording and authenticated payload,
and use the existing review-only path image. If the diagnostic needs promotion,
add a separately authenticated derived diagnostic; do not regenerate the raw
measurement merely to obtain a preferred adjudication.

Tasks 8 and 9 are documentation-only once their cited evidence is verified.
Task 10 uses `benchmarks.gate` directly and does not depend on Task 8/9 text.
Task 11 is the only step that combines all gates.

## Load-bearing conventions

- Before substantive work, re-read the worktree-root `CLAUDE.md`; also read
  `.claude/CLAUDE.md` if it appears. Re-read after any branch/worktree switch.
- CGAL and exact-kernel authority outrank COMPAS convenience. Exact decisions
  use exact CGAL predicates; reporting doubles never feed deciding paths.
- Pixi exclusively. Run pytest with `-n auto`; after changes use
  `pytest --testmon`. Never bare `python`, `pytest`, pip, conda, poetry, or venv.
- TDD for every behavior change: observe RED for the intended reason, then the
  smallest GREEN implementation. Never edit reference tests to force green;
  never skip or xfail.
- No workaround without root-cause analysis first. CI/upstream-impacting fixes
  require options before implementation.
- No fallback behavior, conditional imports, `HAS_*` patterns, bare
  `except: pass`, or alternate acceptance paths.
- No magic tolerances. With COMPAS use `compas.tolerance.TOL` while respecting
  its strict-positive, mixed-relative, and falsy-default traps.
- Python 3.9-compatible in `tools/`; strict mypy is the type gate. No `__all__`
  in `__init__.py`; no Java-shaped result/factory/class hierarchies.
- One responsibility per file; 1,500 lines requires refactoring. The existing
  Task-6 framework is frozen through Task 11, not an invitation to generalize.
- Docs use MkDocs Markdown and Google-style docstrings, never rST.
- Commit by explicit pathspec on this non-main branch. Author and committer:
  `Jelle Feringa <jelleferinga@gmail.com>`. No Co-Authored-By or Codex text.
- Never mutate `main`/`master`; no push unless explicitly requested. Never
  remove a worktree or branch without explicit user consent. No force flags.

## Environment booby traps

- The Pixi environment and editable CGAL build already exist in this worktree.
  If lost, recover with `pixi install && pixi run _editable-rebuild` from this
  worktree; compilation takes minutes.
- Bare PATH Python is stale and may lack `_stock_2`. A resulting import failure
  is environment damage, not production evidence.
- The baseline deliberately exits through pytest's failure status because 11
  reds are expected. `pixi run red-manifest` owns the both-directions verdict;
  do not interpret raw pytest exit 1 as a new defect without reading JUnit.
- `pixi run -e docs docs` is the strict docs build. Plain `pixi run docs` uses
  the wrong environment for this gate.
- Native CGAL gates must be serialized. Do not run concurrent native audits
  against the shared editable build.
- The Task-7 benchmark artifact is a one-shot content-addressed result. Changing
  any of its four files invalidates its stamp and source identity.
- The path diagnostic files are outside Git. Their pixels are review aids, not
  authority. The detailed image is more legible than the saturated full-path
  overview at spacing `0.025`.
- `.superpowers/sdd/.../progress.md` is ignored but load-bearing resume state;
  do not confuse its presence with a tracked branch modification.
- `pyproject.toml` and `mkdocs.yml` are fleet-rebase surfaces. Preserve additive
  changes and run the plan's contention checks before editing.

## Suggested file structure

```text
benchmarks/mathsm.py                                  existing; modify only after MC-013 ruling
benchmarks/measurement_claim_results/
  2026-08-29-3f62304932dd-generator-033fc8b07e19/    existing; immutable
  2026-08-29-f65154fb2075-benchmark-9e0c635f36bd/    existing; immutable
tools/measurement_claim_benchmark_schema.py          existing; possible Task-7 correction
tools/measurement_claim_benchmark_validation.py      existing; possible Task-7 correction
tools/measurement_claim_artifact_validation.py       existing; preserve raw-artifact authority
tools/measurement_claim_ledger.py                    existing; exact joint consumer
tests/tools/test_measurement_claim_benchmark_result.py existing; focused Task-7 contract
tests/tools/test_measurement_claim_ledger.py          existing; live 14/14 boundary
docs/measurement_claims.md                            existing; possible Task-7 evidence correction
docs/benchmarks.md                                    existing; Task 8
docs/continuous_engagement.md                         existing; Task 8
docs/continuous_engagement_cost.md                    existing; Task 8
docs/exactness.md                                     existing; inspect semantically for Task 8
docs/segment_site_mat.md                              existing; Task 8
docs/oblique_edge_cost.md                             existing source warning; avoid redundancy
docs/radius_ladder_ablation.md                        create in Task 9
tests/benchmarks/test_gate_attribution.py             create in Task 10
docs/superpowers/state/backlog.md                     modify in Task 11 only
docs/superpowers/plans/2026-08-28-coherence-wave1.md modify checkboxes/header only
mkdocs.yml                                            modify in Task 9 only after contention check
```

## Risk + budget

No token or time budget is specified in the master plan. Remaining planned
scope is one Task-7 review/correction plus Tasks 8–11.

Known risks:

- MC-013's reason and disposition live inside a content-addressed artifact.
  A naive wording edit either breaks validation or launders a new judgment
  through old evidence. Preserve raw bytes and make adjudication authority
  explicit.
- The review-only path diagnostic identifies an operation but not a peak
  station. Treating the plotted marker as a station would create false evidence.
- Task 8's plan names six files, but the accepted progress ruling found actual
  performance claims in fewer pages and an existing warning in the source page.
  Semantic classification outranks blind sentence insertion.
- Task 9 records historical configurations that are not reconstructible from
  shipped knobs. The provenance warning is mandatory; do not manufacture a
  runnable command.
- Task 10 may prove the generators do not diverge at cap 60. That is a genuine
  STOP finding, not permission to weaken the assertion or choose a convenient
  cap.
- After Task 11, Wave 2 still owns absorption into the auditor frontier. This
  continuation does not merge, rebase, push, or mutate `main`.

## Start

Run this exact command first:

```bash
cd /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence && \
git status --short --branch && \
git rev-parse HEAD && \
pixi run pytest tests/tools/test_measurement_claim_benchmark_result.py \
  tests/tools/test_measurement_claim_ledger.py -q -n auto
```

Expected starting state: clean `codex/sdd-coherence` at `f3a6d2e`; focused
tests green. If either differs, stop and reconcile live state before applying
this prompt. Then print:

```text
Plan: docs/superpowers/plans/2026-08-28-coherence-wave1.md
Step: Task 7 semantic review — MC-013
Criterion: authenticated 14/14 ledger without a general claim stronger than its evidence
Evidence: one-shot artifact valid; path diagnostic review-only; human acceptance unresolved
Next: present the smallest technically sound MC-013 ruling; do not start Task 8
```

End of prompt.
