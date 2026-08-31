# Coherence programme — the plan, and the single live backlog

**This file is the one place open work lives**, and it is now a *plan*, not an
inventory: it defines the end state, gives a checkable oracle for it, and orders
every open item into the wave that reaches it. Everything else under
`docs/superpowers/` is immutable history — plans carry artifact-verified status
headers, specs record designs, other `state/` files are frozen ledgers. When an
item closes, the closing commit removes it here.

## The end state — five invariants

**Maximally coherent** means all five hold at once. Each is stated so its
violation is detectable, not aspirational.

- **I1 — one history.** `main` is the only long-lived branch and contains every
  effort; what was deliberately not merged is pinned by an `archive/*` tag. No
  standing divergence, no frozen refs, one checkout, zero auxiliary worktrees.
- **I2 — one truth per claim.** Every numeric or performance claim in comments,
  docs, and memos is either generated from a committed measured artifact or
  names the exact command + configuration that reproduces it. A claim that
  cannot be re-earned is deleted. (The corpus runner exists to make this cheap.)
- **I3 — one meaning for red.** The full suite's expected-red set is enumerated
  in one manifest with a reason and a backlog link per entry; red outside the
  manifest is a defect *by definition*, green inside it is a finding. No skip,
  no xfail, ever — deliberate reds are the mechanism, the manifest is their
  accounting.
- **I4 — one live work document.** This file. Plans immutable once execution
  starts; status lives in headers, derived from artifacts; a lint enforces that
  every plan has a header.
- **I5 — one build story.** Docs and CLAUDE.md point only at checkouts that
  exist; env bootstrap is documented; lint, strict typing, baseline, and strict
  docs run in CI on `main`.

**Non-goal, stated so nobody "fixes" it:** the 12 quality-gate reds are the
*product* gate — red until the generators earn green. Coherence is their
accounting (I3), never their suppression.

## The convergence oracle

Run these when the waves complete; all must hold:

```
git branch                    → main only (plus archive/* tags)
git worktree list             → one line
pytest (full baseline)        → red set == the manifest, exactly, both directions
claim audit                   → every registered claim carries artifact or command
CI on main                    → lint + strict mypy gates + baseline + docs --strict green
backlog.md                    → "Open" empty except items the user explicitly parks
```

## The path — three waves

### Wave 1 — claims and instruments (landed locally)

The landed slice is on `codex/sdd-coherence`, the one sanctioned side branch
to be absorbed in Wave 2 (S1). C1 remains open; its deliberate red stays in
the canonical manifest until the exact-tangency investigation closes.

The completed claim audit covers only the frozen Task-5 extractor population;
it does not establish repository-wide I2. Local red-manifest and plan-header
instruments exist, while CI enforcement remains Wave 2/P4.

| id | item | notes |
| --- | --- | --- |
| C1 | *(instrument, optional in this wave)* identify which layer decides the exact-tangency case | translation-variance at a rim-on-boundary contact; the deliberate red property stays until this closes. `machining_metric_validity.md` case 4 |

### Wave 2 — programme close (the critical path; owned by the auditor-convergence fleet)

In execution inside the programme — do not duplicate: R1 (frontier
reconciliation; P0 ledger 36/36), R2 (C1-verdict tri-valuation; P1 T4), R3
(false-certificate control; P1 T4), R6 (replay theorem; P2), R8 (swept-prefix
falsifiability; P1 T4), R4 (CI enforcement; P4 — its stated blocker: the
full-suite tail).

**At programme close, the close itself must include:**

| id | item |
| --- | --- |
| W2.1 | Every `required`/`dependent` disposition of the 36 certifier commits consumed or explicitly discarded — then `jf/toolpath-redesign` archive-tagged and deleted; `codex/exact-certified-adaptive-phase1-t9-zero-guide` retired (fully contained, verified 0 missing) |
| W2.2 | The programme's own five plans stamped with status headers, same convention as the legacy eight |
| S2 | The 4 red retrace tests retired (their plan is superseded — D1); the manifest shrinks 17 → 13 |
| S1 | `codex/sdd-coherence` (this branch, grown by Wave 1) absorbed into the frontier; branch retired |
| R5 | Python floor raised to match `typing.Self` imports (≥3.11); the 6-error `types-adaptive` baseline cleared, not carried |

### Wave 3 — the singleton

| id | item |
| --- | --- |
| W3.1 | Fast-forward `main` to the closed frontier (it is 0 behind today — pure ff); delete `codex/auditor-convergence-sdd`; remove its worktree. One branch, one checkout |
| W3.2 | CI green on `main` with the full gate set + the red-manifest diff + the plan-header lint. Fix stale build pointers (CLAUDE.md worktree paths) in the same pass |
| W3.3 | Re-issue the due-diligence memo from committed measured artifacts — the current one predates the oblique-edge finding and states parity unqualified |
| W3.4 | Run the convergence oracle; park or close every remaining line of this file |

**Parked lane — capability work, deliberately not coherence-blocking:** C4
(depletion model: score a true helical entry; until then the entry criterion is
unsatisfiable — 240° floor and `V ≥ √3·r` proved in `loop_radius_degeneracy.md`)
and C5 (arc motion vocabulary for tangent-continuous inter-chain links). These
resume after Wave 3 or in parallel by explicit choice.

## Current red snapshot

`docs/red_manifest.json` is the sole current red authority. This table is its
human-readable snapshot.

| tests | count | reason | closes with |
| --- | ---: | --- | --- |
| `cap-120-default` quality product (2 generators × 3 pockets) | 6 | default-cap product gate: both registered generators emit identical streams and fail 5-7 criteria per pocket | generator work (parked lane C4/C5 feeds this) |
| `cap-40-attribution` quality product (2 generators × 3 pockets) | 6 | cap-40 attribution product gate: both registered generators fail current quality criteria on all three pockets | generator work (parked lane C4/C5 feeds this) |
| `test_quality_invariants.py::test_moving_the_pocket…` | 1 | deliberate: engagement is translation-variant at exact rim-on-boundary tangency (machining_metric_validity.md case 4) | backlog C1 |
| `tests/adaptive/{test_generator,test_route_retrace_generator}` | 4 | plan superseded (backlog D1); retirement scheduled through the auditor programme | backlog S2 |

Anything red beyond these 17 is a defect, full stop.

## Authority and frozen refs (until Wave 2 closes them)

- **Canonical frontier:** `codex/auditor-convergence-sdd`, driven by the
  five-phase programme with its own progress tracker
  (`.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`) and
  reconciliation ledger. Its plans are off-limits to outside edits.
- **Frozen sources — commit to neither:** `jf/toolpath-redesign@73d5372`
  (certifier source, 36/36 dispositions) and
  `codex/exact-certified-adaptive-phase1-t9-zero-guide@073a0f7` (integration
  source, fully contained; checkout removed 2026-08-28, branch intact).

## Decisions — all ruled 2026-08-28

| id | ruling | executed |
| --- | --- | --- |
| D1 | Retrace plan **superseded** by auditor P2 | plan stamped; test retirement is S2 |
| D2 | Deadwood: **archive-tag + delete, local only** | 17 branches deleted (zero-unique re-verified at deletion time; 5 pinned by `archive/*` first), 8 stale records pruned, perf worktree removed (branch kept), origin untouched |
| D3 | Main-checkout dirt **discarded** | certifier source bit-clean at `73d5372` |
| D4 | t9-zero-guide + sdd-coherence **checkouts removed** for disk (~2.7 GB) | branches intact at `073a0f7` / verified SHAs |

## Deadwood ledger (EXECUTED 2026-08-28 under D2 — the record of what went and why)

Measured by `git cherry` patch-equivalence against the frontier tip, not by
branch age or name. Deleted: 12 zero-unique branches
(`…-t7-runtime`, `…-t9`, `codex/adaptive-clearing-sp1-completion`,
`jf/2D_Minkowski_Sums`, `jf/adaptive-clearing-sp1`, `jf/geodesics-module`,
`jf/isolines`, `jf/polylines-module`, `jf/standardize-cpp-params`, three
`worktree-agent-*`) and, behind `archive/*` tags, the early exact-TEA line
(`codex/exact-certified-adaptive-phase1`) plus the shelved circle-oracle line
(`…-t6/-t7/-t8/-t10` — `continuous_engagement_cost.md` is its verdict). Kept:
`perf/exact-rational-representation` (falsification record),
`jf/update_reconstruction_example` (parked for a `main` PR), the frozen
sources, the frontier, `main`.

## Conventions, so this does not regrow

1. **The commit that lands a task updates its plan's status header in the same
   commit.** Bodies are immutable; the header is the truth.
2. **One live document — this file.** A TODO anywhere else is lost by
   construction; the three-place scatter this replaced proved it.
3. **Status is derived from artifacts, never asserted.** Checkbox counting
   called a fully-landed plan 0/80 done; `git cherry` and test runs are the
   measurement.
4. **A new deliberate red enters the manifest in the same commit that
   introduces it** — with its reason and the item that will close it.
