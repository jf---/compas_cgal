# Wave-1 kickoff prompt

Paste everything below the rule as the FIRST message of a fresh session.

---

I'm starting **Wave 1 — claims and instruments** of the coherence programme for
`compas_cgal` (worktree
`/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence`, branch
`codex/sdd-coherence`). Commit `4425fe8` (the Wave-1 implementation plan) is
the last hard checkpoint.

## What I want you to do

Execute `docs/superpowers/plans/2026-08-28-coherence-wave1.md` task-by-task
using superpowers:subagent-driven-development (fresh subagent per task, review
between tasks; every subagent runs on **opus, never fable**). Tasks 1–11, in
order — Task 4 commits its tool BEFORE running it so the provenance stamp can
record a clean tree. Tick the plan's checkboxes as you go and update its
status header in the same commit that lands each task.

## Read first, in dependency order

1. `docs/superpowers/plans/2026-08-28-coherence-wave1.md` — the plan you are
   executing. Complete code for the three tools is in it; do not redesign.
2. `docs/superpowers/state/backlog.md` — the programme: five invariants
   (I1–I5), three waves, the red manifest v1 (11 entries), authority and
   frozen refs. Wave 1 implements I2/I3/I4.
3. `CLAUDE.md` at the worktree root — exact-kernel discipline, pixi-exclusive
   build rules. (Its "Where the pixi manifest lives" section names the OLD
   t9-zero-guide worktree — that checkout was removed 2026-08-28; YOUR
   worktree is the one above, and the manifest is on your branch.)
4. `docs/superpowers/state/2026-08-23-auditor-reconciliation-ledger.md` — why
   two branches are frozen and what the fleet owns. You never touch theirs.
5. Memory (auto-loaded): `measure-the-thing-not-a-proxy`,
   `shared-worktree-stage-by-pathspec`,
   `ask-agents-for-unreported-measurements`,
   `subagent-execution-is-the-default`, `pixi-lives-in-the-worktree` — all
   load-bearing here, especially the first (re-run a number before writing it
   anywhere) and the second (commit by `git commit -- <paths>`, never trust a
   clean tree to stay clean).

## Boundary contract (= your input)

You consume existing, verified interfaces — do not modify them:

```python
# benchmarks/gate.py
def gate_pocket(name: str, tool_diameter: float = GATE_TOOL_DIAMETER,
                tea_cap_deg: float = GATE_CAP_DEG) -> PocketSpec: ...
GATE_GENERATORS: Dict[str, Callable[[PocketSpec], ToolpathResult]]  # 2 entries
GATE_CAP_DEG = 120.0

# benchmarks/cli.py  (run as: python -m benchmarks.cli corpus …)
#   corpus --name {smoke,analytic,complexity,necks,precision,topology,
#                  degeneracy,external,all} --tool-diameter F --cap-deg F
#          --out DIR [--no-digits]
# "all" aggregates the six authored corpora; external is NOT in the aggregate.

# benchmarks/report.py
MARKDOWN_NAME = "benchmark_report.md"
JSON_NAME = "benchmark_report.json"
def write_report(records: List[MeasurementRecord], out_dir: Path) -> Tuple[Path, Path]: ...
```

The suite's current red set (verify at Task 1 Step 3): exactly the 11 entries
of `docs/red_manifest.json` as specified in plan Task 2 — 6 product-gate, 1
tangency property, 4 superseded retrace.

## Deliverables

Per the plan: `tools/{red_manifest,plan_headers,measured_run}.py` +
`tests/tools/` (Tasks 2–4) · `docs/red_manifest.json` (Task 2) · a committed
run under `benchmarks/results/<date>-<sha>/` with `stamp.json` (Task 4) ·
`docs/measurement_claims.md` ledger, skeleton then N/N dispositions
(Tasks 5–7) · qualifiers on exactly six living docs pages (Task 8) ·
`docs/radius_ladder_ablation.md` with the provenance admonition (Task 9) ·
`tests/benchmarks/test_gate_attribution.py` (Task 10) · backlog + plan-header
closeout (Task 11). Three new pixi tasks: `red-manifest`, `plan-headers`,
`measured-run` (incantations verbatim in the plan).

## Exit gates

From plan Task 11 — all must hold:

> `pixi run red-manifest && pixi run plan-headers && pixi run -e docs docs &&
> pixi run lint && pixi run mypy --strict tools` — all green/exit 0.

Plus: claims ledger has zero `pending` rows; the Task-4 artifact's
`stamp.json` has `"dirty": false`. If a gate slips, write
`docs/superpowers/state/wave1-gate-analysis.md` stating which gate, the
measured evidence, and the options — BEFORE any workaround. Two gate outcomes
are findings, not failures: plan Task 2 Step 6 (manifest corrected to observed
reality, stated in the commit) and Task 10 Step 2 (generators fail to diverge
at cap 60 — stop and report).

## Decoupling rule

You are fully independent of the auditor-convergence fleet. Their worktree
(`…/compas_cgal_prs-auditor-convergence-sdd`) is ACTIVE and DIRTY — never run
git there, never edit their files, never touch the frozen refs
`jf/toolpath-redesign@73d5372` and
`codex/exact-certified-adaptive-phase1-t9-zero-guide@073a0f7`. Your only
shared file is `pyproject.toml` (you add three task keys); a future rebase
resolves by key-union. If `mkdocs.yml` shows fleet commits after `fa59120`,
skip nav edits (plan Tasks 5/9 say how).

## Load-bearing conventions

- Pixi EXCLUSIVELY, from your worktree. Never bare `python`/`pytest` — the
  PATH python has no `_stock_2` and every hand-rolled invocation tests the
  wrong binary.
- Commit per task, by pathspec: `git add <paths> && git commit -F msg -- <paths>`.
  Author AND committer `Jelle Feringa <jelleferinga@gmail.com>`. No
  attribution lines, ever. Extremely concise messages.
- No skip/xfail; named exceptions; no magic tolerances; Google docstrings;
  mkdocs admonitions, never rST; `ruff format` before committing;
  `mypy --strict` clean on `tools/`.
- Plans are immutable except their status header and checkboxes.
- Re-measure before transcribing ANY number — two of two previously audited
  measurement comments were defective. A claim you didn't re-run is a claim
  you didn't verify.

## Environment booby traps

- **There is no env yet in your worktree** — Task 1 creates it:
  `pixi install && pixi run _editable-rebuild` (CGAL compiles; minutes).
- `pixi run pytest`/`baseline` set `SKBUILD_EDITABLE_SKIP` internally — that's
  why hand-rolled pytest silently tests a stale extension.
- If `pyproject.toml`/`pixi.lock` get edited, the editable hook file is
  deleted; any `pixi run <task>` regenerates it (a following bare-interpreter
  `ImportError: _stock_2` is env damage, not a code bug).
- Task 4 ordering trap: commit `tools/measured_run.py` FIRST, then run it —
  otherwise `stamp.json` records `"dirty": true` and the artifact test fails
  by design.
- `pixi run -e docs docs` is the strict docs build (plain `pixi run docs`
  fails on a missing extension in the default env).

## Suggested file structure

```
tools/                                  create: __init__.py, red_manifest.py,
                                        plan_headers.py, measured_run.py
tests/tools/                            create: __init__.py + 3 test files
tests/benchmarks/test_gate_attribution.py   create
docs/red_manifest.json                  create
docs/measurement_claims.md              create
docs/radius_ladder_ablation.md          create
benchmarks/results/<date>-<sha>/        created by Task 4's run, committed
docs/superpowers/state/backlog.md       modify (Task 11 closeout only)
docs/superpowers/plans/2026-08-28-coherence-wave1.md   modify (checkboxes+header)
pyproject.toml                          modify (three task keys)
mkdocs.yml                              modify ONLY if fleet hasn't touched it
```

## Risk + budget

No token/time budget is set in the plan; the sizing is eleven tasks, one
commit each. Known risks, from the plan's self-review: `pyproject.toml` is the
single rebase-conflict surface against the fleet (key-union resolves);
Task 4's corpus run and Task 1's CGAL build are the only long steps; Task 10
may surface a genuine finding (no divergence at cap 60) — that stops the task,
it does not weaken the assertion.

## Start

```bash
cd /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence && \
git status --short --branch && pixi install && pixi run _editable-rebuild
```

Then Task 1 Step 3 (the sanity gate), then dispatch Task 2.

End of prompt.
