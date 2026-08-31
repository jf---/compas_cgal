# Codebase Distillation Review Implementation Plan

> **status: approved for review execution** — design approved 2026-08-31; no review pass begun.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` to execute this plan task-by-task.
> `supervising-plan-execution` is required for the three long-running reading
> passes. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Read every scoped first-party code variant three complete times and
produce an evidence-backed, adversarially converged surgery ledger that
preserves every exceptional capability while proposing one authority per
domain truth.

**Architecture:** A compact review instrument produces an ordinary path/ref
manifest, a capability ledger, a findings ledger, and a final surgery ledger.
The primary reviewer performs three full readings with different questions;
semantic-family panels propose falsifiable challenges, and validated oracles
adjudicate all machine-checkable findings. Product code remains unchanged.

**Tech Stack:** Git worktrees and refs, Python 3.12, `pathlib`, Pixi, CMake,
pytest-xdist, pytest-testmon, Ruff, strict mypy, MkDocs, Markdown, TSV.

**Spec:**
`docs/superpowers/specs/2026-08-31-codebase-distillation-review-design.md`

## Global Constraints

- Execute in a new durable agent-owned worktree on a non-main branch created
  from the committed plan branch.
- Before reading any worktree or branch, read its root `CLAUDE.md` completely
  and also `.claude/CLAUDE.md` when present.
- Never mutate, checkout, rebase, merge, stash, clean, delete, or commit from a
  source worktree.
- Pause edits while reviewers read a worktree; invalidate and reread any
  variant that changes.
- Do not introduce application-level fingerprinting, receipt, identity, or
  review-pinning machinery. Git IDs are descriptive state only.
- Preserve all user-owned dirty and untracked content.
- Production code, product tests, build configuration, source branches, and
  source worktrees remain unchanged throughout this plan.
- Exact textual duplicates may share one row only after direct Git comparison;
  textually different variants remain separate through Pass 2.
- The preservation unit is a capability, never a file or branch.
- No disposition is assigned before Pass 3; `UNKNOWN` is the default when
  evidence is insufficient.
- Every reviewer finding follows the spec's finding contract and names a
  falsifier.
- Panel agreement never accepts a finding; oracles or Jelle adjudicate it.
- Native editable builds run sequentially, never concurrently.
- Every pytest command uses `-n auto`; affected Python gates use `--testmon`
  after a focused baseline.
- Pixi owns every Python, test, build, type, lint, docs, and tool command.
- Tests are never skipped, marked expected-failure, or weakened.
- Seeded defects live only in an isolated liveness worktree, are never
  committed, and are fully reversed before a reading begins.
- Every commit uses author and committer
  `Jelle Feringa <jelleferinga@gmail.com>`, has no attribution trailer, and
  carries one concise deliverable.
- The plan stops after Jelle receives the surgery ledger. It does not authorize
  surgery implementation.

---

## File structure

### Created by this plan

- `tools/distillation_review.py` — validate manifest, capability, finding, and
  surgery artifact structure and coverage; no value judgments.
- `tests/tools/test_distillation_review.py` — focused validator contracts.
- `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md` — one
  active step, evidence, drift, blockers, and next command.
- `docs/superpowers/state/2026-08-31-distillation-manifest.tsv` — scoped
  variants and three pass states.
- `docs/superpowers/state/2026-08-31-distillation-capabilities.md` — crown-jewel
  and capability ledger.
- `docs/superpowers/state/2026-08-31-distillation-findings.tsv` — proposed and
  adjudicated findings.
- `docs/superpowers/state/2026-08-31-distillation-surgery.md` — final authority
  graph and dispositions.
- `docs/codebase_distillation_review.md` — readable methodology, exact maturity,
  evidence, limitations, and handoff.

### Modified by this plan

- `pyproject.toml` — add `distillation-review` validation task.
- `mkdocs.yml` — add the developer page.

### Explicitly untouched

- all files under `src/`;
- all product and reference tests outside
  `tests/tools/test_distillation_review.py`;
- CMake and CI;
- benchmark results and figures; and
- all source branch/worktree content.

---

### Task 1: Establish the review execution boundary

**Files:**

- Create:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: this plan, its governing spec, local refs, live worktrees, and all
  applicable policy files.
- Produces: ratified execution worktree, source-worktree policy map, and one
  active-step ledger.

- [ ] **Step 1: Create the isolated execution worktree**

From the current planning worktree, run:

```bash
git status --short --branch
git worktree add \
  /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-codebase-distillation-review \
  -b codex/codebase-distillation-review \
  codex/codebase-distillation-review-plan
```

Expected: the planning worktree remains clean; the new worktree is on
`codex/codebase-distillation-review`.

- [ ] **Step 2: Read policy before any review command**

In the execution worktree, read:

```bash
cat CLAUDE.md
test ! -f .claude/CLAUDE.md || cat .claude/CLAUDE.md
```

Then list worktrees and branch tips:

```bash
git worktree list --porcelain
git branch -vv
git remote -v
```

For every live source worktree, read its root policy files from that worktree
before inspecting its code. Record any differing policy in the progress file.

- [ ] **Step 3: Create the progress ledger**

Use `apply_patch` to create:

```markdown
# Codebase distillation review progress

**Plan:** `docs/superpowers/plans/2026-08-31-codebase-distillation-review.md`
**Step:** Task 1 — establish review execution boundary
**Criterion:** every source worktree has a read policy and remains unmodified
**Evidence:** none yet
**Next:** ratify the branch and worktree scope

## Worktree policies

## Pass coverage

- Pass 1: not started
- Pass 2: not started
- Pass 3: not started

## Oracle liveness

## Drift

- none

## Blockers

- none

## Next exact command

`git for-each-ref --format='%(refname:short)' refs/heads`
```

- [ ] **Step 4: Verify the documentation-only boundary**

Run:

```bash
git status --short
git diff --check
```

Expected: only the progress file is untracked; no source worktree changed.

- [ ] **Step 5: Commit the boundary**

```bash
git add .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'start codebase distillation review'
```

### Task 2: Build the minimal artifact validator

**Files:**

- Create: `tools/distillation_review.py`
- Create: `tests/tools/test_distillation_review.py`
- Modify: `pyproject.toml`

**Interfaces:**

- Consumes:
  `validate_review(manifest: Path, capabilities: Path, findings: Path,
  surgery: Path | None) -> None`.
- Produces: `DistillationArtifactError`; CLI command
  `python -m tools.distillation_review validate`; Pixi task
  `distillation-review`.

- [ ] **Step 1: Write failing parser and coverage tests**

Create focused tests with these helpers and contracts:

```python
from pathlib import Path

import pytest

from tools.distillation_review import DistillationArtifactError
from tools.distillation_review import validate_review


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _manifest(
    tmp_path: Path,
    *,
    duplicate_variant: bool = False,
    pass_3: str = "complete",
) -> Path:
    header = (
        "VARIANT\tREF\tWORKTREE\tPATH\tFAMILY\tSCOPE\tSAME_AS\t"
        "PASS_1\tPASS_2\tPASS_3\tDRIFT\tNOTE\n"
    )
    row = (
        "frontier:src/stock_2.cpp\tfrontier\t/worktree/frontier\t"
        "src/stock_2.cpp\tstock\tincluded\t\tcomplete\tcomplete\t"
        f"{pass_3}\tclean\tCAP-0001\n"
    )
    return _write(
        tmp_path / "manifest.tsv",
        header + row + (row if duplicate_variant else ""),
    )


def _complete_manifest(tmp_path: Path) -> Path:
    return _manifest(tmp_path)


def _capabilities(tmp_path: Path) -> Path:
    return _write(tmp_path / "capabilities.md", "## CAP-0001 — Exact stock\n")


def _findings(tmp_path: Path) -> Path:
    header = (
        "ID\tPASS\tFAMILY\tVARIANT\tLOCATION\tQUOTE\tCAPABILITY\t"
        "ISSUE\tCONSUMER\tLOSS_IF_CHANGED\tPROPOSED_CONDENSATION\t"
        "ORACLE\tCOST\tFALSIFIER\tEVIDENCE\tSTATUS\n"
    )
    return _write(tmp_path / "findings.tsv", header)


def _surgery(
    tmp_path: Path,
    *,
    disposition: str = "KEEP",
    falsifier: str = "A smaller implementation proves consumer equivalence",
) -> Path:
    row = (
        f"| CAP-0001 | {disposition} | Exact stock | exact predicate | "
        f"Stock consumer | loss argument | focused test | {falsifier} | none |\n"
    )
    return _write(
        tmp_path / "surgery.md",
        "| Capability | Disposition | Retained owner | Valuable nucleus | "
        "Consumer contracts | Loss argument | Oracle | Falsifier | Residual risk |\n"
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n"
        + row,
    )


def test_manifest_requires_unique_variant_and_three_pass_states(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path, duplicate_variant=True)
    with pytest.raises(DistillationArtifactError, match="duplicate variant"):
        validate_review(manifest, _capabilities(tmp_path), _findings(tmp_path), None)


def test_complete_surgery_requires_three_reads_per_scoped_variant(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path, pass_3="not-started")
    surgery = _surgery(tmp_path)
    with pytest.raises(DistillationArtifactError, match="three complete readings"):
        validate_review(manifest, _capabilities(tmp_path), _findings(tmp_path), surgery)


def test_remove_requires_loss_oracle_falsifier_and_retained_owner(tmp_path: Path) -> None:
    surgery = _surgery(tmp_path, disposition="REMOVE", falsifier="")
    with pytest.raises(DistillationArtifactError, match="REMOVE.*falsifier"):
        validate_review(_complete_manifest(tmp_path), _capabilities(tmp_path), _findings(tmp_path), surgery)
```

Also test allowed pass states, finding statuses, capability references, and
`UNKNOWN` acceptance with named missing evidence.

- [ ] **Step 2: Run the tests to verify RED**

```bash
pixi run pytest -- tests/tools/test_distillation_review.py -n auto -q
```

Expected: import failure because `tools.distillation_review` does not exist.

- [ ] **Step 3: Implement the smallest validator**

Use `pathlib.Path`, `csv.DictReader`, Markdown heading extraction, and one named
exception. Keep the module below 500 lines. It must:

1. validate required TSV columns and allowed values;
2. reject duplicate ordinary variant IDs;
3. verify every capability referenced by a finding or surgery row exists;
4. verify every scoped row has three complete passes before surgery validation;
5. verify `CONDENSE`, `ABSORB`, and `REMOVE` carry the spec-required evidence;
6. accept `UNKNOWN` only with `MISSING_EVIDENCE` and `FALSIFIER`; and
7. print counts only after validation succeeds.

The validator must not rank, infer, or assign dispositions.

- [ ] **Step 4: Add the Pixi task**

Add under `[tool.pixi.tasks]`:

```toml
distillation-review = "python -m tools.distillation_review validate"
```

The CLI defaults to the four state paths from the spec and treats an absent
surgery ledger as an in-progress review.

- [ ] **Step 5: Run focused GREEN and affected gates**

```bash
pixi run pytest -- tests/tools/test_distillation_review.py -n auto -q
pixi run affected -- tests/tools/test_distillation_review.py
pixi run lint
```

Expected: focused and affected tests pass; Ruff passes.

- [ ] **Step 6: Commit the validator**

```bash
git add tools/distillation_review.py tests/tools/test_distillation_review.py pyproject.toml
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'validate distillation review artifacts'
```

### Task 3: Ratify the complete code universe

**Files:**

- Create: `docs/superpowers/state/2026-08-31-distillation-manifest.tsv`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: every local branch, active remote-only branch, live worktree,
  policy file, and the spec's inclusion/exclusion rules.
- Produces: one row per differing scoped variant with columns:
  `VARIANT, REF, WORKTREE, PATH, FAMILY, SCOPE, SAME_AS, PASS_1, PASS_2,
  PASS_3, DRIFT, NOTE`.

- [ ] **Step 1: Inventory branches, worktrees, and dirt**

Run from the execution worktree:

```bash
git for-each-ref --format='%(refname:short)' refs/heads
git for-each-ref --format='%(refname:short)' refs/remotes
git worktree list --porcelain
```

In every live worktree run, without changing directory outside that worktree
for Git operations:

```bash
git branch --show-current
git status --short --branch --untracked-files=all
git remote -v
```

Record user-owned dirty/untracked paths as scoped variants when they contain
first-party code; otherwise record the exclusion reason without reading their
contents unnecessarily.

- [ ] **Step 2: Enumerate first-party textual paths per ratified ref**

For each ratified ref, run:

```bash
review_ref=codex/held-reference-corpus
git ls-tree -r --name-only "$review_ref"
```

Repeat with `review_ref` set to each ratified ref. Include the extensions and
named configuration files in the spec. Exclude generated, vendored, binary,
result, and lock artifacts with one explicit reason per path class.

- [ ] **Step 3: Group only exact textual duplicates**

For the same path appearing on multiple refs, run pairwise:

```bash
review_ref_a=jf/toolpath-redesign
review_ref_b=codex/held-reference-corpus
review_path=src/stock_2.cpp
git diff --quiet "$review_ref_a" "$review_ref_b" -- "$review_path"
```

Repeat for every same-path ref pair. Exit zero permits `SAME_AS` to name the
ordinary retained variant label. Nonzero means both variants remain separate.
Never group semantically similar but textually different variants before Pass
3.

- [ ] **Step 4: Assign one primary semantic family**

Use exactly the eight families from the spec. Cross-family dependencies belong
in `NOTE`; they do not duplicate rows.

- [ ] **Step 5: Create the manifest**

Generate the TSV with all scoped rows at `not-started`, exclusions carrying a
reason, and no disposition field. Use ordinary labels such as
`held-reference-corpus:src/stock_2.cpp`.

- [ ] **Step 6: Prove manifest completeness**

Compare manifest path/ref counts against each `git ls-tree` inventory. Sample
at least one included and one excluded path from every class. Run:

```bash
pixi run distillation-review
git diff --check
```

Expected: validator succeeds in in-progress mode; no malformed rows.

- [ ] **Step 7: Update progress and commit**

Record branch/worktree counts, scoped variant count, duplicate groups,
exclusion counts, dirty variants, and the next exact command. Then commit:

```bash
git add docs/superpowers/state/2026-08-31-distillation-manifest.tsv \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'ratify distillation review scope'
```

### Task 4: Establish truthful baselines and prove oracle liveness

**Files:**

- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`
- Create: `docs/superpowers/state/2026-08-31-distillation-findings.tsv`

**Interfaces:**

- Consumes: manifest lineages, each lineage's policy-declared commands, known
  deliberate-red authority, and isolated liveness worktree.
- Produces: baseline matrix, live/dead oracle matrix, and initial missing-oracle
  findings.

- [ ] **Step 1: Record sequential branch-specific baselines**

Never run native editable builds concurrently. For the integrated frontier,
run sequentially:

```bash
pixi run baseline
pixi run lint
pixi run types-adaptive
pixi run types-audit
pixi run -e docs docs
pixi run red-manifest
pixi run task9-mat-compile-gate
pixi run audit-native
```

For each other executable lineage, run the commands declared by its own
`CLAUDE.md` and Pixi manifest. Record exit status, passed/failed counts,
deliberate-red classification, elapsed time, and any command that cannot run.
Do not repair failures.

- [ ] **Step 2: Create an isolated liveness worktree**

```bash
git worktree add \
  /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-distillation-oracle-liveness \
  -b codex/distillation-oracle-liveness \
  codex/codebase-distillation-review
```

Read that worktree's policies. Confirm it starts clean.

- [ ] **Step 3: Validate cheap oracle classes with reversible seeds**

Use `apply_patch` for one source defect at a time. For each seed:

1. run the named focused gate and observe the expected failure;
2. reverse the exact patch with `apply_patch`;
3. rerun the gate and observe the baseline result;
4. confirm `git diff --check` and `git status --short` are clean; and
5. record `live`, `dead`, or `misdirected` with evidence.

Validate at least:

- exact decision detection in the engagement kernel;
- frame/unit detection under strict adaptive typing;
- malformed motion detection at the public operation boundary;
- depletion error detection in the stock oracle;
- residual-stock detection in the coverage gate;
- cap-violation detection in the independent engagement audit; and
- planner/consumer disagreement in the quality gate.

Never alter reference tests. Never commit a seed.

- [ ] **Step 4: Create the findings ledger header**

Create the exact TSV header:

```text
ID	PASS	FAMILY	VARIANT	LOCATION	QUOTE	CAPABILITY	ISSUE	CONSUMER	LOSS_IF_CHANGED	PROPOSED_CONDENSATION	ORACLE	COST	FALSIFIER	EVIDENCE	STATUS
```

Add one `proposed` missing-oracle finding for each dead, misdirected, or absent
claim boundary. Do not add a condensation proposal before Pass 3.

- [ ] **Step 5: Remove the clean liveness worktree only after user approval**

Report that the source branch contains no commits and the worktree is clean.
Do not remove either automatically. Continue review with it present if approval
is not given.

- [ ] **Step 6: Validate and commit baseline evidence**

```bash
pixi run distillation-review
git diff --check
```

Update progress with baseline and oracle matrices, then commit:

```bash
git add docs/superpowers/state/2026-08-31-distillation-findings.tsv \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'record distillation review baselines'
```

### Task 5: Perform Pass 1 capability archaeology

**Files:**

- Create:
  `docs/superpowers/state/2026-08-31-distillation-capabilities.md`
- Modify:
  `docs/superpowers/state/2026-08-31-distillation-manifest.tsv`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: complete manifest, live oracle matrix, every scoped variant.
- Produces: capability ledger and `PASS_1=complete` for every scoped variant.

- [ ] **Step 1: Create the capability ledger schema**

Create one section per capability using exactly:

```markdown
## CAP-0001 — Name

- **Locations:**
- **Mathematical claim:**
- **Product claim:**
- **Invariants:**
- **Named failures:**
- **Consumers:**
- **Counterexamples:**
- **Oracles:**
- **Measured value:**
- **Uniqueness:**
- **Competing implementations:**
- **Replacement cost:**
- **Held relevance:**
- **Buchli relevance:**
- **Proposed nucleus:** not assessed before Pass 3
- **Surrounding machinery:** not assessed before Pass 2
- **Missing evidence:**
- **Disposition:** not assessed before Pass 3
```

- [ ] **Step 2: Read every family bottom-up**

Read all scoped variants in this order:

1. general bindings and compatibility;
2. stock/depletion/containment/coverage;
3. point/guarded/continuous engagement;
4. medial-axis/segment-site/reachability/neck;
5. motion vocabulary and legacy generation;
6. adaptive candidates/traversal/transactions/replay/generation;
7. benchmarks/quality/comparison/visualization;
8. build/bindings/typing/CI/scripts/tools/executable docs.

Within each family, start from pathological tests and fixtures, trace into the
deciding symbols, then trace to every consumer. Read every line of every scoped
variant; do not sample large files.

- [ ] **Step 3: Record capabilities without dispositions**

For each variant, add or reference every capability it carries. Preserve
branch-only ideas, counterexamples, exact distinctions, bounded algorithms,
performance results, and honest incomplete states. Mark `PASS_1=complete` only
after its full contents and consumers have been read.

- [ ] **Step 4: Run the historical salvage trigger**

When a current test, comment, plan, or consumer points to deleted unique code,
inspect only the named lineage between its documented branch points. Add any
recovered capability with an ordinary historical location and explicit current
consumer status. Do not broaden this into all-history review.

- [ ] **Step 5: Prove Pass 1 coverage**

Run:

```bash
pixi run distillation-review
```

Expected: every scoped row has `PASS_1=complete`; Pass 2/3 remain
`not-started`; every variant references at least one capability or carries an
explicit `no-unique-capability-observed` note without disposition.

- [ ] **Step 6: Independent Pass 1 omission review**

Dispatch one unbriefed reviewer against the manifest and capability ledger:

> Name capabilities, counterexamples, or replacement costs the ledger may have
> missed. Every finding must cite a scoped location and state what would be lost.

Verify each cheap finding; add verified omissions. Report unchecked expensive
findings explicitly.

- [ ] **Step 7: Commit Pass 1**

Update progress with exact coverage counts and next family order. Run strict
docs and commit:

```bash
pixi run -e docs docs
git diff --check
git add docs/superpowers/state/2026-08-31-distillation-capabilities.md \
  docs/superpowers/state/2026-08-31-distillation-manifest.tsv \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'map codebase capabilities'
```

### Task 6: Ratify semantic-family adversarial panels

**Files:**

- Modify:
  `docs/superpowers/state/2026-08-31-distillation-capabilities.md`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: Pass 1 capability map and family-specific failure surfaces.
- Produces: one ratification card per family with opponent, cost-bearer,
  resource, interface, domain, and null coverage.

- [ ] **Step 1: Enumerate and weight each family's failure surface**

Use 8–15 falsifiable failure classes per family. Weight each as coarse
probability `{0.1, 0.3, 0.6, 0.9}` times consequence `{1, 3, 10}`. Include
at minimum loss of exceptional capability, false equivalence, weakened
exactness, hidden consumer breakage, irrelevant performance evidence, and
product-claim overreach.

- [ ] **Step 2: Ratify Held and Buchli as standing global lenses**

Probe each in the third person, three independent samples, against checked
public anchors. Require calibrated commitments and no invented citation.
Replace a weak named prior with its de-named commitments. Never attribute
panel output to the real person.

- [ ] **Step 3: Construct each family panel**

Use `adversarial-panel-construction`. Every panel has:

- one documented opponent;
- one downstream cost-bearer;
- one resource/lifetime detector;
- one contract/interface detector;
- relevant mathematical/domain specialists; and
- one null reviewer.

Use 4–6 named or de-named lenses plus the null. No more than two seats may come
from one field. Record uncovered and single-point failure classes.

- [ ] **Step 4: Audit blind spots**

For every uncovered class, distinguish candidate-pool failure from panel
selection failure. Generate a new candidate only for a pool failure. Split any
class whose two best detectors rely on disjoint commitments.

- [ ] **Step 5: Commit ratification cards**

Append concise cards to the capability ledger, update progress, run strict
docs, and commit:

```bash
pixi run -e docs docs
git diff --check
git add docs/superpowers/state/2026-08-31-distillation-capabilities.md \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'ratify distillation review panels'
```

### Task 7: Perform Pass 2 authority and product-truth review

**Files:**

- Modify:
  `docs/superpowers/state/2026-08-31-distillation-capabilities.md`
- Modify: `docs/superpowers/state/2026-08-31-distillation-findings.tsv`
- Modify: `docs/superpowers/state/2026-08-31-distillation-manifest.tsv`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: Pass 1 map, family panels, baselines, live oracles, every scoped
  variant.
- Produces: authority/complexity map and `PASS_2=complete` for every scoped
  variant.

- [ ] **Step 1: Read every family top-down**

Read in reverse product order:

1. build/CI/tools/executable docs;
2. benchmark/quality/comparison/visualization;
3. adaptive generation/replay/transactions/traversal;
4. legacy generation and motion vocabulary;
5. MAT/reachability/neck;
6. engagement;
7. stock/depletion/coverage;
8. general bindings and compatibility.

Read every line again. Trace each output to its deciding kernels and every
consumer back to its authority.

- [ ] **Step 2: Apply the Held comparison trace**

For every link in the spec's Held path, record the authoritative symbol,
consumer evidence, limitation, and comparison maturity:
`stronger`, `equivalent`, `weaker`, or `incomplete`.

Do not infer G-code, controller, coverage, cap, or cycle-time capability when
the corresponding consumer boundary is absent.

- [ ] **Step 3: Apply the Buchli system trace**

Record where claims stop across contract, independent geometry, parsed G-code,
controller simulation, and machine execution. Separate a valid mathematical
claim from an incomplete system claim.

- [ ] **Step 4: Record authority and complexity**

Populate each capability with:

- current authority;
- competing authorities;
- essential mathematical complexity;
- accidental orchestration/ownership complexity;
- uncomposed exceptional kernels;
- candidate condensation seams; and
- missing consumer or oracle.

Do not assign final dispositions.

- [ ] **Step 5: Dispatch family panels blind**

Use the ratified panel, artifact, and finding contract. Do not brief reviewers
on owned failure classes. Dispatch in capacity-limited blind waves; do not
share earlier output with later reviewers. Include the null control.

Every returned finding must contain `LOCATION`, `QUOTE`, `ISSUE`, `ORACLE`,
`COST`, and `FALSIFIER`; reject incomplete submissions as not-yet-findings.

- [ ] **Step 6: Verify findings by class**

Batch C0 through existing gates. Verify C1 with focused characterization,
differential, metamorphic, or bounded benchmark checks. Record C2 with severity
if not run. Route C3 to Jelle without converting it to fact.

Return verification evidence only to the originating reviewer when a follow-up
is required.

- [ ] **Step 7: Prove Pass 2 coverage and commit**

Run:

```bash
pixi run distillation-review
pixi run -e docs docs
git diff --check
```

Expected: every scoped row has Pass 1 and Pass 2 complete; no final
dispositions. Update progress and commit:

```bash
git add docs/superpowers/state/2026-08-31-distillation-capabilities.md \
  docs/superpowers/state/2026-08-31-distillation-findings.tsv \
  docs/superpowers/state/2026-08-31-distillation-manifest.tsv \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'map codebase authority and complexity'
```

### Task 8: Perform Pass 3 condensation falsification

**Files:**

- Modify:
  `docs/superpowers/state/2026-08-31-distillation-capabilities.md`
- Modify: `docs/superpowers/state/2026-08-31-distillation-findings.tsv`
- Modify: `docs/superpowers/state/2026-08-31-distillation-manifest.tsv`
- Create: `docs/superpowers/state/2026-08-31-distillation-surgery.md`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: two complete readings, capability/authority map, verified findings,
  explicit C2/C3 queues, every scoped variant.
- Produces: third reading, falsified/retained reduction hypotheses, proposed
  authority graph, and disposition ledger.

- [ ] **Step 1: Create the surgery ledger structure**

Use:

```markdown
# Codebase distillation surgery proposal

## Executive claim boundary

## Held product path

## Buchli evidence ladder

## Proposed authority graph

## Capability dispositions

| Capability | Disposition | Retained owner | Valuable nucleus | Consumer contracts | Loss argument | Oracle | Falsifier | Residual risk |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

## C2 queue

## C3 decisions for Jelle

## Explicitly preserved research

## Proposed surgery slices

## Prohibited actions before approval
```

- [ ] **Step 2: Read by cross-branch semantic family**

Read every scoped variant a third time, grouping all variants of one capability
side by side. Use the eight families, but order each family from competing
implementations toward the proposed nucleus.

- [ ] **Step 3: State and attack one reduction hypothesis per candidate**

For every candidate condensation, absorption, quarantine, or removal, write:

> Capability X can be preserved by nucleus Y, contracts Z, and counterexamples
> Q; machinery M is not semantically required.

Then actively search hidden consumers, branch-only distinctions, exactness
loss, named failure loss, performance regression, ABI differences, and research
option loss. Record verified attacks as findings.

- [ ] **Step 4: Assign dispositions conservatively**

Use only the six spec dispositions. `REMOVE` requires no surviving unique
capability and a retained owner for every consumer outcome. `CONDENSE` requires
a named nucleus and preserved contracts. `ABSORB` requires matching semantics,
lifetime, failures, and performance. Unresolved exceptional work defaults to
`QUARANTINE` or `UNKNOWN`, never `REMOVE`.

- [ ] **Step 5: Mark Pass 3 only after the loss argument is complete**

Set `PASS_3=complete` per variant only when all carried capabilities have a
ledger disposition and every proposed reduction has an oracle/falsifier or is
explicitly `jelle-c3`.

- [ ] **Step 6: Validate complete coverage**

Run:

```bash
pixi run distillation-review
```

Expected: all scoped variants have three complete readings; every disposition
passes its evidence contract; C2 and C3 remain visible.

- [ ] **Step 7: Commit Pass 3**

Update progress, run strict docs and diff check, then commit:

```bash
pixi run -e docs docs
git diff --check
git add docs/superpowers/state/2026-08-31-distillation-capabilities.md \
  docs/superpowers/state/2026-08-31-distillation-findings.tsv \
  docs/superpowers/state/2026-08-31-distillation-manifest.tsv \
  docs/superpowers/state/2026-08-31-distillation-surgery.md \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'propose codebase distillation surgery'
```

### Task 9: Adversarially converge the surgery proposal

**Files:**

- Modify: `docs/superpowers/state/2026-08-31-distillation-surgery.md`
- Modify: `docs/superpowers/state/2026-08-31-distillation-findings.tsv`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: complete three-pass corpus, all dispositions, C2/C3 queues.
- Produces: zero-new-verified-disposition-change convergence result.

- [ ] **Step 1: Build the final synthesis failure surface**

Weight at minimum: false Held superiority, absent G-code outcome, hidden proof
defect, consumer breakage, discarded exceptional algorithm, weakened exactness,
machine/operator hazard, retained unconsumed machinery, duplicate authority,
irrelevant performance evidence, and dead oracle.

- [ ] **Step 2: Dispatch the final panel and null control**

Use `adversarial-convergence-review`. Give every reviewer the spec, manifest,
capability ledger, findings, and surgery proposal. Brief on no assigned classes.
Require the full finding contract.

- [ ] **Step 3: Route and verify every new finding**

Always run C0. Verify plausible or high-impact C1. Queue or run C2 by severity.
Present C3 to Jelle. Every verified defect in the review instrument leaves a
validator test; every verified capability omission updates the capability and
surgery ledgers.

- [ ] **Step 4: Return evidence only to originating reviewers**

When verification produces new information, send it only to the reviewer who
raised the finding. Do not broadcast it to the panel.

- [ ] **Step 5: Repeat fresh rounds**

Repeat until one fresh round produces zero new verified disposition changes,
all C0/C1 findings are adjudicated, and remaining C2/C3 items are explicit.

- [ ] **Step 6: Record convergence and commit**

Run:

```bash
pixi run pytest -- tests/tools/test_distillation_review.py -n auto -q
pixi run distillation-review
pixi run lint
pixi run -e docs docs
git diff --check
```

Record rounds, verified/rejected counts, C2/C3 remainder, null-only findings,
and blind spots. Commit:

```bash
git add docs/superpowers/state/2026-08-31-distillation-surgery.md \
  docs/superpowers/state/2026-08-31-distillation-findings.tsv \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md \
  tools/distillation_review.py tests/tools/test_distillation_review.py
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'converge codebase surgery proposal'
```

### Task 10: Publish the review and stop at the approval gate

**Files:**

- Create: `docs/codebase_distillation_review.md`
- Modify: `mkdocs.yml`
- Modify:
  `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`

**Interfaces:**

- Consumes: converged surgery ledger and complete validation evidence.
- Produces: readable engineering handoff and explicit Jelle approval gate.

- [ ] **Step 1: Write the developer page**

Document:

- reviewed scope and explicit exclusions;
- three-pass methodology and coverage count;
- crown jewels and their retained homes;
- proposed authority graph;
- Held comparison maturity by boundary;
- Buchli evidence maturity by boundary;
- proposed condensation slices;
- explicit quarantined research;
- C2/C3 remainder;
- limitations and blind spots; and
- statement that no product-code surgery is authorized.

- [ ] **Step 2: Add the page to MkDocs navigation**

Add one entry under the developer/engineering section:

```yaml
- Codebase distillation review: codebase_distillation_review.md
```

- [ ] **Step 3: Run final gates**

```bash
pixi run pytest -- tests/tools/test_distillation_review.py -n auto -q
pixi run distillation-review
pixi run lint
pixi run -e docs docs
git diff --check
```

Expected: all commands exit zero; the manifest reports every scoped variant
read three times.

- [ ] **Step 4: Run the mandated prohibited-mechanism vocabulary scan**

Apply the user-level scan from `AGENTS.md` to the complete diff. Any newly
introduced prohibited mechanism is removed; references required to describe
pre-existing code remain clearly observational and never become design
recommendations.

- [ ] **Step 5: Close progress at the approval boundary**

Set:

```markdown
**Step:** Task 10 — Jelle surgery approval
**Criterion:** Jelle approves or amends every C3 decision and surgery slice
**Evidence:** three-pass coverage and convergence gates complete
**Next:** wait for Jelle; do not write a surgery implementation plan
```

- [ ] **Step 6: Commit the publication**

```bash
git add docs/codebase_distillation_review.md mkdocs.yml \
  .superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md
GIT_AUTHOR_NAME='Jelle Feringa' \
GIT_AUTHOR_EMAIL='jelleferinga@gmail.com' \
GIT_COMMITTER_NAME='Jelle Feringa' \
GIT_COMMITTER_EMAIL='jelleferinga@gmail.com' \
git commit -m 'publish codebase distillation review'
```

- [ ] **Step 7: Stop**

Present the surgery ledger, C2 queue, C3 decisions, blind spots, and branch
status to Jelle. Do not implement, remove, merge, rebase, push, or clean
anything.

---

## Plan self-review map

| Spec requirement | Implementing task |
| --- | --- |
| Isolated, non-mutating review boundary | Task 1 |
| Compact mechanically checked artifacts | Task 2 |
| All active code variants and exclusions | Task 3 |
| Truthful baseline and live oracles | Task 4 |
| Complete capability archaeology | Task 5 |
| Failure-surface-specific panels | Task 6 |
| Held/Buchli authority and product review | Task 7 |
| Third complete reading and loss falsification | Task 8 |
| Adversarial convergence | Task 9 |
| Publication and explicit user approval | Task 10 |

The plan contains one instrumentation task and no product-code task. The
validator checks shape and coverage only; it never decides what code is
valuable.
