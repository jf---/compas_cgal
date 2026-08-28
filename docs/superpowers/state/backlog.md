# Backlog — the single live document

**This file is the one place open work lives.** Everything else under
`docs/superpowers/` is immutable history: plans carry a status header stating
what landed (verified by artifact audit, 2026-08-28 — checkbox state in plan
bodies was never maintained and is noise), specs record what was designed, and
the other `state/` documents are frozen kickoffs and ledgers. When an item here
closes, the closing commit removes it here.

## Authority and frozen refs

- **Canonical frontier:** `codex/auditor-convergence-sdd`. It strictly contains
  `codex/exact-certified-adaptive-phase1-t9-zero-guide` and is driven by the
  live five-phase programme `2026-08-23-auditor-convergence-p0…p4` with its own
  progress tracker (`.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`)
  and reconciliation ledger
  (`state/2026-08-23-auditor-reconciliation-ledger.md`). Do not stamp, edit, or
  commit into the programme's plans from outside its own review discipline.
- **Frozen sources — do not commit to either while the programme runs:**
  `jf/toolpath-redesign@73d5372` (certifier source; every one of its 36 unique
  commits carries a disposition in the ledger, 36/36) and
  `codex/exact-certified-adaptive-phase1-t9-zero-guide@073a0f7` (integration
  source, fully contained in the frontier).
- This file was authored on `codex/sdd-coherence`, a one-commit docs branch off
  the frontier tip, kept disjoint from every file the programme touches so it
  rebases or fast-forwards in with zero conflict. Absorbing it — and then
  retiring `codex/sdd-coherence` — is itself an item below.

## Decision needed from the user

| id | decision |
| --- | --- |
| D1 | **`2026-08-07-exact-inter-route-retrace` — finish or supersede.** The only genuinely in-flight legacy plan: `route_retrace_replay.py` + its test were never written, and the 4 standing red tests in `tests/adaptive/` are this plan's. It predates the pivot to the quality instrument. Finishing means writing the replay leg; superseding means stamping the plan and retiring its red tests' claim on the suite. |
| D2 | **Deadwood execution.** Verdicts in the ledger at the bottom of this file; branch/worktree deletion requires explicit consent and has not been performed. |
| D3 | **Main checkout dirt** (`compas_cgal_prs` on `jf/toolpath-redesign`): `docs/examples/example_isolines.py` (+32 lines, an abandoned experiment — `from asyncio import log` is an accidental import) and the empty stray file `0`. The reconciliation ledger explicitly excludes both and never touches them; they are the user's to discard or keep. |

## Open items

Origins: **R*n*** = `review-2026-08-22-t9-zero-guide.md` recommendation *n* ·
**A/B/C** = the benchmark-corpus plan's "Open after the plan landed" appendix.

### In execution — owned by the auditor-convergence programme (do not duplicate)

| id | item | where in the programme |
| --- | --- | --- |
| R1 | Reconcile the split frontier (`jf/toolpath-redesign` × codex tip) | P0 ledger, 36/36 dispositions recorded; consumption of `dependent`/`required` commits runs through P1 |
| R2 | C1 — verdict must not conflate violation with exhaustion (tri-valued, not Boolean) | P1 Task 4 (source audit finding: "Boolean conflates violation/exhaustion; proof/refinement input only") |
| R3 | C2 — a station-green segment can hide a cap violation | P1 Task 4: the annular-rib false-certificate control is disposed `required`; spiral-rib generalization `dependent` |
| R6 | C3 — replay/reproducibility instrument | P2 (`p2-replay-theorem`); native replay transaction landed through P1 Task 4C |
| R8 | C7 — swept-prefix theorem needs falsifiable footing | P1 Task 4: theorem + falsifier retained as refinement input, never certification authority |
| R4 | Wire the existing gates into CI | P4 (`p4-enforcement-evidence`); the programme's stated open blocker: identify the extreme full-suite tail before enforcement |

### Open — unowned

| id | item | notes |
| --- | --- | --- |
| R5 | Python floor: metadata says `>=3.9`, `adaptive/` imports `typing.Self` (3.11+) | the 6 `types-adaptive` errors are a tracked programme baseline; the floor is not |
| R7 | Commit one measured corpus run | turns the prose performance claims into checkable ones; none committed to date |
| A | Audit the 34 measurement-asserting comments across 12 files | 2 of 2 checked so far had a defect: one false (`MEAS` placeholder, fixed in `29050b0`), one correct but mislabelled (an unreachable "off" column). Re-run each claim; name the knob that produces every number |
| B | Write the two unrecorded ablations into docs | (1) peak/count provably opposed at a forced station — a constraint, not a knob; (2) the four-way ladder decomposition incl. the 213.6° slotting cut that justifies the gate. Both currently live only in a code comment (`29050b0`) and this line |
| C1 | Identify which layer decides the exact-tangency case | engagement is translation-variant when the rim lies exactly on a cleared boundary; the red property in `test_quality_invariants.py` is deliberate. Candidates: arrangement representation of a measure-zero contact, run extraction in `engagement_at`, boundary convention. `machining_metric_validity.md` case 4 |
| C2 | The oblique-edge cliff, and the parity qualifier | `center_domain()`: 5 ms axis-aligned → 17.5 s oblique integer vertices → >90 s generic rotation; mechanism explicitly unestablished (`oblique_edge_cost.md`). Every corpus pocket is axis-parallel, so **every performance figure including Held parity is best-case; `review.md` still states parity unqualified** |
| C3 | Second gate cap (40–100°) | at `GATE_CAP_DEG=120` both registered generators emit byte-identical paths, so the 3×2 gate cannot attribute a defect to either |
| C4 | Depletion model: score a true helical entry | `_replay_kind` refuses Z+XY motion; the workaround encoding silently degrades to a no-op rapid. Until fixed, the entry criterion is unsatisfiable — the 240° floor and the `V ≥ √3·r` requirement are proved in `loop_radius_degeneracy.md` |
| C5 | Arc motion vocabulary | inter-chain links cannot be tangent-continuous while loops are full circles; loops need distinct entry/exit tangency points (arcs) |
| S1 | Absorb `codex/sdd-coherence` into the frontier, then retire the branch and its worktree | one docs commit, disjoint files, rebases clean |

## Deadwood ledger (verdicts recorded; execution awaits explicit consent — D2)

Measured 2026-08-28 by `git cherry` patch-equivalence against the frontier tip
`fa59120`, not by branch age or name.

| branch | unique patches | verdict |
| --- | --- | --- |
| `codex/exact-certified-adaptive-phase1-t7-runtime` | 0 | delete |
| `codex/exact-certified-adaptive-phase1-t9` | 0 | delete |
| `codex/adaptive-clearing-sp1-completion` | 0 | delete |
| `jf/2D_Minkowski_Sums` | 0 | delete |
| `jf/adaptive-clearing-sp1`, `jf/geodesics-module`, `jf/isolines`, `jf/polylines-module`, `jf/standardize-cpp-params` | 0 | delete (long merged) |
| `worktree-agent-a00e2496`, `-a0ec673d`, `-a668f212` | 0 | delete (agent residue) |
| `codex/exact-certified-adaptive-phase1` | 5 | archive-tag then delete — early exact-TEA line, re-implemented rather than merged |
| `codex/…-t6` / `-t7` / `-t8` / `-t10` | 2–5 (shared) | archive-tag then delete — the shelved circle-oracle research line (`continuous_engagement_cost.md` is its verdict: "research result, not a component") |
| `perf/exact-rational-representation` | 3 | **keep branch** — the falsification record cited by `continuous_engagement_cost.md`; its worktree (clean) can go |
| `jf/update_reconstruction_example` | 1 | park — unrelated old example fix, candidate for a `main` PR |
| `jf/toolpath-redesign` | 36 | **frozen** — certifier source of the reconciliation ledger; untouchable until the programme closes it |
| `codex/exact-certified-adaptive-phase1-t9-zero-guide` | 0 vs frontier | keep until the programme completes, then retire (fully contained) |
| `main` | — | keep |

Worktrees: the 8 `/private/tmp/compas_cgal_prs-*` records point at directories
that no longer exist — `git worktree prune` removes only the stale records.
Live worktrees: main checkout (frozen source), `…-t9-zero-guide` (frozen
source), `…-auditor-convergence-sdd` (the programme's, ACTIVE and dirty — never
touch), `…-perf-exact-rational` (clean; removable, branch stays),
`…-sdd-coherence` (this pass; retire under S1).

## Conventions, so this does not regrow

1. **The commit that lands a task updates its plan's status header in the same
   commit.** Plan bodies are immutable once execution starts; the header is the
   truth, the body is the record of what was planned.
2. **One live backlog — this file.** A TODO written anywhere else is lost by
   construction; the three-place scatter this file replaced (unchecked boxes,
   a plan appendix, an unintegrated review) proved it.
3. **Status is derived from artifacts, never asserted from memory.** Checkbox
   counting called a fully-landed plan 0/80 done; `git cherry` and test runs
   are the measurement.
