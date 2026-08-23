# Auditor Convergence P0 Reconciliation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a complete, evidence-backed equivalence ledger for the two
divergent branch frontiers and a reproducible baseline without changing
production code.

**Architecture:** Treat Git history as provenance, not as integration authority.
Inventory every source-side commit, prove equivalent/superseded dispositions at
symbol and contract-test level, and route only non-equivalent proof-bearing work
to later plans. Record all baseline commands against immutable SHAs and keep
deliberately red product gates separate from regression failures.

**Tech Stack:** Git, Pixi, pytest-xdist, pytest-testmon, Ruff, strict mypy,
MkDocs, Markdown.

**Spec:** `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md`

## Global Constraints

- Work only in `codex/auditor-convergence-sdd`; never checkout, rebase, commit,
  push, or mutate either source branch/worktree.
- Bind integration source `073a0f7`, certifier source `73d5372`, and merge base
  `1860167` in every ledger.
- Stage P0 changes documentation only. No production/test/configuration patch
  is allowed.
- A commit is `equivalent` only with symbol-level evidence and a passing
  contract test on the integration side.
- A commit is `superseded` only when the integration contract is strictly
  stronger and its test names are recorded.
- Every pytest command uses `-n auto`; Pixi owns every command.
- Do not skip, mark expected-failure, or edit a reference test.
- Commit as `Jelle Feringa <jelleferinga@gmail.com>` with concise messages.

---

### Task 1: Freeze the branch topology

**Files:**

- Create: `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: local refs `073a0f7`, `73d5372`, `1860167`.
- Produces: immutable `Source topology` section used by every later P0 task.

- [ ] **Step 1: Verify the isolated branch and both source refs**

Run each command from the convergence worktree:

```bash
git branch --show-current
git status --short --branch --untracked-files=all
git show-ref --verify refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide
git show-ref --verify refs/heads/jf/toolpath-redesign
git merge-base 073a0f7 73d5372
git rev-list --left-right --count 073a0f7...73d5372
```

Expected: convergence branch is active; source refs are `073a0f7` and
`73d5372`; merge base is `1860167`; divergence is `254 36`.

- [ ] **Step 2: Create the ledger header with exact evidence**

Create the state document with this complete schema:

```markdown
# Auditor reconciliation ledger

**Recorded:** 2026-08-23
**Integration source:** `073a0f7f833da440fe323a0a046ad31550579f82`
**Certifier source:** `73d53729564851d69eb1c28d6c2f1e68350cd20f`
**Merge base:** `1860167929e50f38fdae8d67ef77e9c967364f1c`
**Integration divergence:** 254 commits
**Certifier divergence:** 36 commits

## Source topology

## Commit dispositions

## Required dependency closure

## Baseline

## P0 acceptance
```

Under `Source topology`, record `git status --short --branch
--untracked-files=all` for all three worktrees and explicitly exclude the
unrelated `jf/toolpath-redesign` dirty files from every comparison.

- [ ] **Step 3: Create the cold-start progress ledger**

Create `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md` with:

```markdown
# Auditor convergence progress

**Plan:** P0 reconciliation
**Task:** 1 — freeze topology
**Branch:** `codex/auditor-convergence-sdd`
**Worktree:** `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-auditor-convergence-sdd`
**Integration source:** `073a0f7`
**Certifier source:** `73d5372`

## Accepted commits

- `993d4f8` — approved convergence design

## Last verified commands

## Open blockers

- none

## Next exact command

`git log --reverse --format='%H%x09%s' 1860167..73d5372`
```

- [ ] **Step 4: Validate documentation scope**

Run:

```bash
git diff --check
pixi run -e docs docs
```

Expected: no whitespace errors; strict MkDocs exits zero.

### Task 2: Inventory all 36 certifier-side commits

**Files:**

- Modify: `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: frozen source topology from Task 1.
- Produces: exactly one disposition row for every certifier-side commit.

- [ ] **Step 1: Generate the ordered source inventory**

Run:

```bash
git log --reverse --format='%H%x09%s' 1860167..73d5372
git diff --name-status 1860167..73d5372
```

Expected: exactly 36 commit rows.

- [ ] **Step 2: Add the disposition table**

Add one row per source commit with these columns:

```markdown
| Source commit | Subject | Symbols/files | Disposition | Integration evidence | Target plan |
| --- | --- | --- | --- | --- | --- |
```

Allowed dispositions are exactly `equivalent`, `superseded`, `required`,
`dependent`, and `unrelated`. No row may contain an undecided value.

- [ ] **Step 3: Inspect every commit rather than classifying by subject**

For each SHA, run:

```bash
git show --stat --oneline <source-sha>
git show --format=fuller --no-ext-diff <source-sha>
git log --all --oneline -- <each-touched-path>
```

Replace `<source-sha>` and `<each-touched-path>` with the literal SHA and paths
from the inventory before execution. Record exact integration symbols and tests
beside each disposition.

- [ ] **Step 4: Prove table cardinality**

Run:

```bash
git rev-list --count 1860167..73d5372
rg -c '^\| `[0-9a-f]{40}` ' docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md
```

Expected: both commands report `36`.

### Task 3: Measure the immutable baseline

**Files:**

- Modify: `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: clean convergence tree at the P0 planning commit.
- Produces: command/result matrix separating regression, product-gate, and
  environmental outcomes.

- [ ] **Step 1: Record the exact tested tree**

Run:

```bash
git rev-parse HEAD
git status --short --branch --untracked-files=all
```

Record the SHA before every command group. If tracked files change during a
test run, discard that result and rerun from a frozen tree.

- [ ] **Step 2: Run bounded proof and typing gates**

Run serially to avoid concurrent editable native rebuilds:

```bash
pixi run pytest -- tests/test_false_certificate.py tests/test_growth_bound.py tests/test_engagement_audit.py -n auto -q
pixi run pytest -- tests/adaptive/test_replay.py tests/adaptive/test_zero_guide_replay.py tests/adaptive/test_segment_oracle.py tests/adaptive/test_circle_oracle.py -n auto -q
pixi run types-adaptive
pixi run lint
pixi run -e docs docs
```

Record exit code, pass/fail counts, wall time, and complete failing node IDs.
The known six mypy errors are regression baseline, not permission to silence
them.

- [ ] **Step 3: Run benchmark-instrument and product gates separately**

Run:

```bash
pixi run pytest -- tests/benchmarks -n auto -q
pixi run pytest -- tests/benchmarks/test_quality.py -n auto -q
```

Classify schema/metric/figure contract failures as instrument regressions.
Classify `test_the_generated_path_is_worth_running` failures as product-gate
failures only when the assertion and measured values match the declared gate.

- [ ] **Step 4: Run the full baseline once**

Run:

```bash
pixi run baseline
```

Record every failure and reconcile it against Steps 2-3. A new or unexplained
failure blocks P0 acceptance.

### Task 4: Prove the required dependency closure

**Files:**

- Modify: `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: all 36 commit dispositions and baseline evidence.
- Produces: ordered patch/test closure for P1 and P2; no code transfer.

- [ ] **Step 1: Trace source dependencies for each required commit**

For each `required` row, inspect imports, native symbols, bindings, and test
fixtures with:

```bash
git show <source-sha>^:<path>
git show <source-sha>:<path>
rg -n '<symbol>' src tests
```

Record every prerequisite source commit as `dependent` unless the integration
side already proves it equivalent or superseded.

- [ ] **Step 2: Compare the nine load-bearing findings explicitly**

Create one subsection each for exact seam validation, bounded disk chain,
compiled growth guard, false-certificate liveness, swept-annulus segment
certificate, rotation-invariant bound, rim reporting, release shared-root
guard/full-turn theorem, and adaptive arc certification.

Each subsection records:

```markdown
- source symbols and commits
- integration symbols and commits
- behavioural difference
- proving negative control
- proving positive control
- destination P1/P2 task
```

- [ ] **Step 3: Validate no production diff exists in P0**

Run:

```bash
git diff --name-only 993d4f8...HEAD
```

Expected before the P0 commit: only the approved spec, five plan documents,
the reconciliation state document, and the SDD progress ledger.

### Task 5: Gate and commit P0

**Files:**

- Modify: `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: Tasks 1-4.
- Produces: reviewer-ready P0 acceptance commit and exact P1 start command.

- [ ] **Step 1: Run the P0 acceptance checks**

Run:

```bash
git diff --check
pixi run -e docs docs
rg -c '^\| `[0-9a-f]{40}` ' docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md
git status --short --branch --untracked-files=all
```

Expected: clean formatting/docs; 36 disposition rows; only planned
documentation is modified.

- [ ] **Step 2: Update the progress ledger**

Set P0 Task 5 complete, record every verified command/result, list no
unadjudicated source commits, and set the next exact command to the first RED
test in P1 Task 1.

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "docs(sdd): reconcile auditor branches"
```

- [ ] **Step 4: Reverify source immutability**

Run:

```bash
git show-ref --verify refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide
git show-ref --verify refs/heads/jf/toolpath-redesign
git status --short --branch --untracked-files=all
```

Expected: original refs remain `073a0f7` and `73d5372`; convergence worktree is
clean.
