# Wave-1 Task-6 Historical Assertions Implementation Plan

> **status: landed** — Tasks 1–3 complete; canonical v2 Task-6 evidence accepted.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve every source number removed by `53135e0` in one authenticated
historical artifact, record the two Wave-1 deviation rulings, then land the
canonical v2 Task-6 measurement artifact and ledger acceptance.

**Architecture:** Store the exact raw Git correction diff plus a minimal
content-addressed manifest. A focused repository-boundary test recomputes the
diff and identity without adding another producer or measurement schema. The
existing v2 generator remains the sole accepted current-measurement path.

**Tech Stack:** Python 3.9-compatible stdlib, Git raw objects, JSON, SHA-256,
pytest, Ruff, strict mypy, Pixi.

**Spec:**
`docs/superpowers/specs/2026-08-29-wave1-task6-historical-assertions.md`

## Global Constraints

- Work only in the durable `codex/sdd-coherence` worktree.
- Do not mutate `main`, `master`, the auditor worktree, or frozen refs.
- Use `git --no-replace-objects` for every authenticated Git read.
- Historical assertions are provenance, never accepted measurements.
- Do not add a CLI, Pixi task, producer family, fallback, or alternate parser.
- Preserve the rejected v1 bundle; never place it in the accepted-result root.
- Use TDD, `pytest -n auto`, Ruff, and strict mypy.
- Commit with Jelle Feringa as author and committer and explicit pathspecs.

---

### Task 1: Historical assertion artifact contract

**Files:**
- Create: `tests/tools/test_measurement_claim_history.py`
- Create: `benchmarks/measurement_claim_history/2026-08-29-53135e04390e/manifest.json`
- Create: `benchmarks/measurement_claim_history/2026-08-29-53135e04390e/source-correction.patch`

**Interfaces:**
- Consumes: raw Git commits `b531d215e6ca06741d040b070ba43164f61abd58`
  and `53135e04390e84bf69aa74dc4d0c1ce6ca308eb4`.
- Produces: one byte-reproducible provenance artifact for MC-001 through
  MC-010. It is not imported by runtime code.

- [x] **Step 1: write the failing repository-boundary test**

The test must:

1. load the manifest with duplicate-key rejection;
2. require its exact key set and literal schema/status values;
3. authenticate the correction commit and its single parent from raw commit
   bytes using SHA-1 or SHA-256 according to object-ID length;
4. run the exact two-path `git --no-replace-objects diff --binary` command;
5. compare the committed patch byte-for-byte with stdout;
6. verify the patch SHA-256;
7. require ordered unique paths and exact MC-001 through MC-010 coverage;
8. require MC-001…008 map to the radial source and MC-009…010 to the advance
   source.

- [x] **Step 2: run the focused test and observe missing-artifact failure**

Run:

```bash
pixi run pytest tests/tools/test_measurement_claim_history.py -q -n auto
```

Expected: failure naming the absent manifest or patch.

- [x] **Step 3: create the exact patch and manifest**

Create `source-correction.patch` from the command in the spec without
normalization. Set `patch_sha256` to the SHA-256 of those exact bytes. Map all
ten claim IDs once in canonical order.

- [x] **Step 4: run focused and adjacent gates**

```bash
pixi run pytest tests/tools/test_measurement_claim_history.py -q -n auto
pixi run pytest tests/tools/test_measurement_claim_ledger.py -q -n auto
```

Expected: history test green; ledger suite retains only the intentional live
Task-6 acceptance red.

### Task 2: Record deviation rulings and preserve v1 rejection

**Files:**
- Create: `docs/superpowers/state/wave1-gate-analysis.md`
- Modify: `docs/superpowers/state/wave1-interim-review.md`
- Create: `benchmarks/measurement_claim_rejections/2026-08-29-53135e04390e-generator-8a73a14731d2/`
- Restore: `docs/measurement_claims.md` to its committed pending state.

**Interfaces:**
- Consumes: interim review findings and the rejected v1 bundle.
- Produces: explicit rulings, S1 absorption instructions, and durable rejected
  evidence outside the accepted root.

- [x] **Step 1: write both deviation rulings**

Record the exact fleet diagnostic test command, the observed twelfth red before
`681f045`, current focused green, why revert violates I3, and the S1
patch-equivalence rule. Record the framework guarantee, measured 7,757-line
cost, frozen ten-claim scope, and post-Task-11 retirement boundary.

- [x] **Step 2: move the rejected bundle without changing its bytes**

Move the untracked v1 directory from `benchmarks/measurement_claim_results/`
to `benchmarks/measurement_claim_rejections/`, add a rejection note containing
its input/result/payload/stamp hashes and semantic rejection reason, and verify
the payload/stamp hashes after the move.

- [x] **Step 3: restore the ledger with `apply_patch`**

Restore the status line and MC-001 through MC-010 evidence cells to the exact
committed Task-5 pending state. Do not use checkout, reset, or stash.

- [x] **Step 4: validate and commit the prerequisite evidence state**

Run the history contract, ledger structural tests, Ruff, and strict mypy for
the new test boundary. Commit only the historical artifact, rejected evidence,
gate analysis, interim review, spec, plan, and restored ledger state by explicit
pathspec.

### Task 3: Canonical v2 Task-6 closure

**Files:**
- Create: `benchmarks/measurement_claim_results/<canonical-v2-name>/`
- Modify: `docs/measurement_claims.md`
- Modify: `docs/superpowers/plans/2026-08-28-coherence-wave1.md`
- Modify: `.superpowers/sdd/2026-08-28-coherence-wave1/progress.md`

**Interfaces:**
- Consumes: clean execution commit, source correction `53135e0`, canonical v2
  generator, historical assertion artifact.
- Produces: exactly one accepted v2 artifact and 10/14 ledger acceptance.

- [x] **Step 1: generate once from the clean prerequisite commit**

```bash
pixi run measurement-claims-generator
```

Expected: exactly one directory below
`benchmarks/measurement_claim_results/`, with v2 envelope, payload, distinct
execution/correction identities, and `dirty: false`.

- [x] **Step 2: validate and render ledger evidence**

Use the canonical validator/renderer only. MC-001 must be `re-earned`;
MC-002 and MC-003 must be independently `corrected`; MC-004 through MC-010
must match their v2 adjudications. MC-011 through MC-014 remain pending.

- [x] **Step 3: run Task-6 acceptance gates**

```bash
pixi run pytest tests/tools -q -n auto
pixi run pytest --testmon tests -q -n auto -k "radial or rho or ordered or spiral"
pixi run red-manifest
pixi run plan-headers
pixi run ruff format --check tools tests/tools
pixi run ruff check tools tests/tools
pixi run mypy --strict tools
```

Expected: all green except only the reds already accepted by
`docs/red_manifest.json`; the ledger acceptance test is green.

- [x] **Step 4: update completion records and commit Task 6**

Tick Task-6 checkboxes and add `6` to the Wave-1 plan header only now. Record
the v2 artifact identities, historical artifact SHA-256, exact gate outputs,
and the post-Task-11 retirement boundary in the ignored progress ledger.
Commit the accepted artifact, ledger, original Wave-1 plan, and review records
with an explicit pathspec and terse Task-6 message.
