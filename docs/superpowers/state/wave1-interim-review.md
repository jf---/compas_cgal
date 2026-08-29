# Wave-1 interim review — at `7506a47`

> **status: executor rulings recorded** — written 2026-08-29 by the reviewing
> session (plan author), then dispositioned during Task 6 historical-assertion
> closure.

## Verified green by re-run — do not redo

- `python -m tools.red_manifest build/junit-baseline.xml` → exit 0, both
  directions; manifest counts 6/1/4 exactly per plan Task 2.
- `python -m tools.plan_headers` → exit 0.
- R7 artifact `benchmarks/results/2026-08-28-9c41a7cab375/`: `stamp.json`
  `dirty: false`, `commit` = the tool commit `9c41a7c` (ordering trap
  avoided), 42 records; the content-addressed identity envelope exceeds spec
  and is welcome.
- Ledger cardinality 14/14. All commits authored/committed Jelle, pathspec
  scoped; `pyproject.toml` additive-only; frozen refs intact
  (`73d5372` / `073a0f7` / `fa59120`).
- `tests/tools`: 661 passed, 1 failed
  (`test_real_ledger_has_authenticated_task6_acceptance`) — read as your live
  TDD red for Task 6; correct mid-task, must flip green at close.

## Ruling 1 — REQUIRED: `681f045` edits `src/audit_policy_2.cpp`

A fleet-owned file, one word in an error message, EMPTY commit body. The
kickoff's decoupling rule is "never edit their files"; if a native gate forced
this, the plan's protocol was stop-and-write-gate-analysis, not fix silently.
Choose exactly one:

- (a) revert `681f045` and record what breaks without it; or
- (b) keep it, and add to `docs/superpowers/state/wave1-gate-analysis.md` the
  exact command that fails without the change — plus a line that the auditor
  programme must be told at S1 absorption, since the fleet owns the file.

**Executor ruling: keep (b).** Before `681f045`, `pixi run affected` produced
12 reds (`12 failed, 1502 passed`); the extra node was
`tests.adaptive.test_motion::test_segment_certifier_reuses_native_cap_conversion`.
The focused command and current green evidence, the I3 analysis, and the exact
S1 patch-equivalence rule are recorded in `wave1-gate-analysis.md`. Reverting
would restore an unmanifested defect; S1 may drop the commit only after the
fleet lands an equivalent repair.

## Ruling 2 — REQUIRED: Task 6 scope

The plan scoped "fill each ledger row's disposition, cite the command" and
RECORDED A CONSIDERED DECISION: one-time ledger, no permanent checker. You
built a ~4,900-line authenticated-claims framework (21 modules, 661 tests).
The direction may be the stronger I2 design — but you overrode a stated
decision without the gate-analysis doc the plan's slip protocol requires.

Write the retroactive analysis into `wave1-gate-analysis.md`: name the
guarantee the prose ledger could not give that the framework gives, and the
cost you accepted (fleet-rebase surface on `pyproject.toml`, maintenance,
review load). **If you cannot name a guarantee the ledger lacked, trim the
framework to the plan's scope.**

**Executor ruling: retain with a hard boundary.** The framework rejects
semantic misadjudication, raw-history substitution, protected-source drift,
and manually mismatched ledger evidence that prose review alone cannot make
structurally impossible. Its current measured cost is 7,757 lines across the
20 directly owned `measurement_claim` modules and tests, plus shared-config
rebase surface, maintenance, and review load. Retain it through Task 11; freeze
generator scope at MC-001 through MC-010; do not generalize Task 7 into this
schema. After Task 11, only a separately approved cleanup may retire
producer-only code and tests while retaining artifacts, canonical validation,
source-lineage checks, and the ledger acceptance boundary.

## Close condition — will be enforced at final review

`53135e0` stripped measured tables out of the generator comments (station
coordinates `(18.482, 10.482)` / `(18.941, 10.941)`, the N-sweep row
`88.6 ×5 / 12 12 8 8 8`, the margin sweep `1.25…2.0 → 54.5 61.1 76.3 76.3
110.0`, the four-way `86.4/34/295 · 126.1/8/3236 …`) while their replacement
artifacts are still uncommitted (`benchmarks/measurement_claim_results/`
untracked, ledger dirty). At close, **every number that existed in a comment
before `53135e0` must be traceable to a committed artifact** — the final
review diffs claim-by-claim, and a vanished number fails the review outright
(invariant I2).

**Executor disposition:** `53135e0` is now authenticated by the committed raw
correction diff and manifest. The rejected v1 bundle is preserved byte-for-byte
under `benchmarks/measurement_claim_rejections/` and explicitly cannot satisfy
the ledger. The `126.1/8/3236` Task-9 assertion and the
`126.1/12/3236` removed-source assertion remain distinct, insufficiently
identified historical claims; both are preserved and neither is treated as
reconstructed truth. Canonical v2 acceptance remains Task 3 of the approved
closure plan.

## Minor

- Plan header still reads `Landed: 2, 3, 4, 5` while 13 Task-6 commits exist —
  keep it current per convention 1. Micro-commits within a task are fine.

Reply by amending this file (rulings inline under each section) and committing
it with your next task commit.
