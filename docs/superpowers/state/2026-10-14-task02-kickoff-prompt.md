# Kickoff prompt: Product Pipeline, Task 0.2 → Task 0.3 → Stage 1

I'm starting **Task 0.2 (exact, idempotent sweeps), then Task 0.3 (gate rules and the ratchet), then Stage 1** of the compas_cgal toolpath project. The work happens in the worktree `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation`, on branch `codex/held-task5-continuation`.

The last hard checkpoint is the Task 0.1 commit at HEAD (subject `gate: in-process engagement …`). It sits on top of the Task 0.0 commits `bdb7534e`, `08ae8664` and `0f7484a5`. This prompt supersedes `docs/superpowers/state/2026-10-13-stage0-kickoff-prompt.md` for everything that remains. Read that file only for history.

## 1. What I want you to do

Execute `docs/superpowers/plans/2026-10-13-product-pipeline.md` with **superpowers:subagent-driven-development**, resuming at Task 0.2. The plan is the authority, as amended on 2026-10-14 (slotting moved from 0.1 to 0.3). Use **plan-status** as the scope gate.

**Resume from the ledger.** The ledger is `.superpowers/sdd/2026-10-13-product-pipeline/progress.md`. It is gitignored through `.git/info/exclude` of the main repo. Tasks 0.0 and 0.1 are complete there. The ledger also holds every ruling, deferred minor finding and open question. The same directory holds the step-level drafts for the next two tasks: `task-0.2-expansion-draft.md` and `task-0.3-expansion-draft.md`. A drafter self-reported low effort on the 0.2 draft, and that self-report is unreliable, so treat both drafts as leads and re-measure what they claim.

**Fan-out (unchanged):**
- **One implementation workstream.** Tasks run in order.
- **Read-only expansion drafters.** At most two run in parallel, each drafting one upcoming task: Task 1.1 next, once 0.2 is active.
- **Per task:** one implementer, then one fresh reviewer, then you integrate and run the gate yourself.

**Effort:** `subagent_type: "opus-medium"`, `model: "opus"`, and **`effort: "medium"` passed explicitly on every dispatch.** The definition's frontmatter is not trusted on its own, and a subagent's self-reported effort is unreliable (memory `subagent-effort-pass-explicitly.md`). Each subagent's first output line states its model; a non-Opus model means stop it. Never spawn Fable. Raising effort needs my go, per task, with a reason.

**Commits:** per-task commits on `codex/held-task5-continuation` are approved for this plan. Author and committer are Jelle Feringa <jelleferinga@gmail.com>, with no attribution lines. Commit by pathspec. **Implementers never commit.** The orchestrator commits after review and the gate run, because the task commit must carry the re-pinned vectors.

**Standing approvals (2026-10-14):**
- **Granted:** re-pinning `_EXPECTED_QUALITY_GATE_MEASUREMENTS` in `tests/benchmarks/test_quality.py:1258-1271` for Task 0.2. Put the old vector, the new vector and the deltas in the commit message.
- **Still need my explicit go, each time:** the Task 0.3 retirements (R1e); the CI changes (show me the YAML diff first); any other reference-test edit; any edit to `docs/red_manifest.json`.

## 2. Read first, in dependency order

1. The plan, `docs/superpowers/plans/2026-10-13-product-pipeline.md`: the Spine, Rulings R1–R4, Global Constraints, Tasks 0.2 and 0.3 (0.3 as amended), and Task 1.1.
2. The ledger, `.superpowers/sdd/2026-10-13-product-pipeline/progress.md`: every `Ruling:` line, the `Open (rule at 0.2 activation …)` line, the 0.3 open conflicts, and the deferred minor findings.
3. `.superpowers/sdd/2026-10-13-product-pipeline/task-0.2-expansion-draft.md`. Section 0 lists its conflicts. Its probe is `build/task02-draft/probe.py`.
4. `.superpowers/sdd/2026-10-13-product-pipeline/task-0.3-expansion-draft.md`. Its §8 lists conflicts C1–C12. Its probes are in `build/task03/`.
5. `docs/superpowers/state/task-0.1-gate-analysis.md`: the 12-case table, old → new, plus why slotting can't move under 0.1.
6. `docs/superpowers/plans/2026-10-13-product-pipeline-traceability.md`: the rows for 0.2, 0.3 and 1.1.
7. The census records in `docs/superpowers/census/2026-10-13/out/`, by finding ID. The verify scripts are in `docs/superpowers/census/2026-10-13/verify/`; `verify_sweep.py` is 0.2's probe.
8. Repo `CONVENTIONS.md` and `docs/exactness.md`.
9. Memory, in `~/.config/dev/projects/-Users-jelle-Code-CADCAM-compas-cgal-prs/memory/`:
   - `product-pipeline-rulings-2026-10-13.md`
   - `toolpaths-over-formalities.md`
   - `measure-the-thing-not-a-proxy.md`
   - `subagent-execution-is-the-default.md`
   - `subagent-effort-pass-explicitly.md`
   - `shared-worktree-stage-by-pathspec.md`
   - `rules-beat-system-attribution.md`
   - `ds-store-breaks-native-dependency-digest.md`

## 3. Boundary contract (your input)

The gate consumes a `ToolpathResult` (`src/compas_cgal/toolpath.py`), unchanged:

```python
@dataclass
class ToolpathOperation:
    geometry: Line | Arc | Circle          # world XY, z = cut plane or clearance
    operation: OperationType               # CUT, LEAD_IN, LEAD_OUT, LINK, RETRACT, PLUNGE (HELIX in Task 1.2)
    path_index: int
    clockwise: bool = False
    start_tangent: np.ndarray | None = None
    end_tangent: np.ndarray | None = None

@dataclass
class ToolpathResult:
    operations: list[ToolpathOperation]
    polyline: np.ndarray                   # (N, 3) visualisation only
```

**What Task 0.1 delivered (now the gate's engagement model):**
- `benchmarks/survey.py` samples each motion in travel order. It measures against the stock minus that motion's own swept prefix, through the native `engagement_line_prefixes` / `engagement_arc_prefixes`.
  - The C++ is in `src/swept_prefix_engagement_2.{h,cpp}`, bound through `engagement_2_bindings.cpp`.
  - The query is non-mutating on a `const Stock2&`.
  - The arc extent is the exact enum `ArcPrefixExtent2`, decided by `CGAL::orientation`.
  - Arc stations come from the exact circle chart atlas.
- `benchmarks/quality.py` converts engagement to textbook angle, immersion and chip thickness from the in-process reading (θ_textbook = θ_rim). The static conversions remain, but no consumer uses them now.
- The generators still decide on **static** engagement under the same cap number; the gate uses in-process engagement. This is documented in `docs/machining_quality.md`.

**Epistemic state**

Content identity "live" means the Task 0.1 commit (HEAD at session start).

| Claim or failure | Classification | Evidence | Content identity | Consequence |
| --- | --- | --- | --- | --- |
| A straight slot into solid stock reads 180° in the gate | verified invariant | `verify/verify_direction.py` → `[0.0, 180.0, 180.0] … [180.0 …]` | live | Task 0.1 done |
| 12 gate vectors re-pinned: cap exceedances fell on 10 of 12 cases; max step fell on all 12 (e.g. 329.49 → 157.15); uncut, degenerate, redundant, slotting and tangent unchanged | verified invariant | `build/task01/gate_vectors_measured.json`; `state/task-0.1-gate-analysis.md` | live | the baseline 0.2 moves from |
| Suite reds equal the manifest (16 = 12 gate + 4 adaptive) | verified invariant | `build/red-manifest-task01c.log` | live | `red-manifest` is the per-task regression check |
| Slotting cannot move under in-process engagement at a 0.5 × cap threshold | verified invariant | `survey.py:100`; gate analysis | live | slotting moves in Task 0.3 rule (c) (plan amended) |
| A replayed loop is not idempotent: plunge 8/12, re-sweep 12/12 | verified invariant (model defect) | `verify/verify_sweep.py` | pre-0.1 (`0f7484a5`); 0.1 did not touch the sweeps | Task 0.2 |
| The authored radius alone gives re-sweep 0/12, plunge still 8/12, so exact entry incidence is also needed | unresolved hypothesis (drafter's number) | `build/task02-draft/probe.py` | pre-0.1 | re-measure before relying on it |
| Four 16 mm chain runs give 8,014 arrangement vertices; the quad gives 32; the removal test drops from 404 ms to 30 ms | unresolved hypothesis (drafter's numbers) | `build/task02-draft/probe.py` | pre-0.1 | re-measure; 8,014 is the plan's own figure |
| 0.3 probes: 1 centre-domain component per pocket; plunges into material 5/5/6 at 120°; plunge XY equals next start 12/12; no LINK at cut depth removes nothing; junction step 204/221/196 vs peak 329/338/325 | unresolved hypothesis (drafter's numbers, pre-0.1 survey) | `build/task03/` | pre-0.1 | re-measure under live; the step numbers will have moved |
| A motion's first engagement sample reads 0°: the tool sits in its own hole and the prefix is empty | verified invariant (orchestrator probe) | `build/task01/attr_probe.py` | live | Task 0.3 rule (f): never take junction values at sample 0 |
| The attribution test now pins the max-step pair (1, 2): peaks are loop 240°, link 180°, loop 180° | verified invariant | `build/task01/attr_probe.py` | live | re-check when rule (f) lands |
| Tests asserting the 12 reds, identical streams at 120°, entries over the cap, plunge per pass, Arc rejection | **false test/oracle premise** | plan Task 0.3 retire list | live | retire under R1e in 0.3, with my go |
| `tests/tools/test_red_manifest.py` hard-codes a mirror of the repository manifest's red list | verified invariant (coupling) | `_repository_red_ids` | live | every manifest edit must also update this mirror |
| TEA-growth lemma false for concave voids (`engagement_2.cpp:407-499`) | unresolved hypothesis | census g4 #4 | live | Stage 4.3 |
| Analytic Held lap = Held on Figure 5; fails at sharp convex corners | verified invariant / gap | `verify_analytic*.py` | corpus `a29ff8c6` | Stage 2 |
| The deleted ordered generator took rect entries from 5 to 1 | the 5→1 is relayed from docs | `git show 9fa6eea8^:src/compas_cgal/engagement_ordered_toolpath.py` | `9fa6eea8^` | Task 1.1: re-measure |

## 4. Rulings in force (from the ledger; do not re-litigate)

- **Frozen generators untouched** (0.1-C1, 0.2-C1). Their probe ring and depletion calls stay as they are. The "generator" side of every shared model means the pipeline, starting in Task 1.1.
- **Replay owns depletion.** Queries are non-mutating (0.1-C3).
- **Witnesses.** Interval refinement with exact comparisons, used to construct a witness, is not "bisection" (0.1-C4).
- **0.2 timing baseline (0.2-Q3).** Use the per-case durations from `build/gate-task01-2026-10-14.log`: a gate-only run with `--durations=0`, 267.5 s total, taken on the 0.1 tree. The pre-0.1 per-case times in `build/junit-baseline.xml` are reported alongside.
- **Open, to rule on at 0.2 activation with the code in view:**
  - **0.2-C2:** entry snapping on frozen output. The drafter's `IncidentLoop.build` moves the phase by up to about 0.9°, which changes the authored path, not just its incidence. The alternative is an honest rule-(e) failure.
  - **0.2-Q1:** non-Pythagorean lines leave the quad under-covering, which conflicts with "exact decisions only on exact models". The alternative is an under/over bracket at roughly 2× replay cost.
  - Rule each one and write a `Ruling:` line to the ledger. Stop and ask only if every path is a guess.
- **0.3 open questions** (draft §8):
  - C1 retirements: need my go.
  - C2 ratchet "improves": the lean is "no cell worse, at least one better".
  - C3 the gate path set.
  - C4/C5 CI: need my sign-off on the diff.
  - C6–C12: rule on these yourself.

## 5. Deliverables

- **Task 0.2:**
  - `Stock2::subtract_circle_sweep` with the authored radius, replacing the `hypot(point_at(0))` surrogate.
  - One capsule model, shared by the replay and (from 1.1) the pipeline.
  - Exact entry incidence.
  - `can_remove_circle` / `_loop_reaches_material` as the redundancy tests.
  - Local depletion, capsules included.
  - Re-pinned redundant counts.
- **Task 0.3:**
  - Rules (a)–(f).
  - The retirements (with my go).
  - `benchmarks/results/gate_vectors.json`.
  - `tools/gate_ratchet.py`.
  - The CI fixes (with my sign-off).
  - The per-stage exit test.
  - `docs/machining_quality.md` numbers read from the JSON.
- **Then Stage 1** (Tasks 1.1–1.7), as far as the session reaches.
- At the end of Stage 0, the mkdocs page `docs/machining_quality.md` is updated once.

## 6. Exit gates (quoted from the plan, as amended)

- **Task 0.2:**
  - "re-cutting a loop from any start and plunging at its own entry change nothing (12/12)";
  - "the arrangement after four 16 mm runs stays under 200 vertices, down from 8,014";
  - "the gate's redundant-operation counts are re-pinned";
  - "generation plus gate replay on each gate pocket takes no more than 2× the 2026-10-13 time" (the baseline per ruling 0.2-Q3).
- **Task 0.3:** "the old generators' vectors are re-pinned under the new rules; slotting moves on at least one pocket under rule (c) (moved here from Task 0.1, 2026-10-14); the ratchet passes on this commit and fails on a synthetic regression."
- **Task 1.1 (120°, all three pockets):** "plunges into material: 1 per pocket; redundant 0, degenerate 0; tangent breaks ≤ rho's 0/0/2; no criterion worse than rho_regulated's 2026-10-13 vector under the Stage 0 gate."

**If a gate slips:** stop. Write `docs/superpowers/state/<task>-gate-analysis.md` with what failed, the measured evidence and the smallest plan amendment, then bring it to me. The precedent is `task-0.1-gate-analysis.md`.

## 7. Decoupling rule

- Task 0.2 is checked on `verify/verify_sweep.py` (idempotence) and on the four-run vertex count first.
- Task 0.3 is checked on a synthetic regression vector for the ratchet.
- The full product gate is the integration check, run by you: `pixi run pytest tests/benchmarks/test_quality.py -k worth_running -n auto --durations=0`, about 4.5 min.
- To re-measure all 12 vectors without the pin assert, run `build/task01/gate_vectors.py` through the rebuild-bypass incantation. It writes `build/task01/gate_vectors_measured.json`. Use it for the re-pin.

## 8. Load-bearing conventions

- **The spine:**
  - every commit moves a gate number on a gate pocket;
  - integrate before inventing;
  - no bisection;
  - exact decisions only on exact models;
  - forced decisions are structured outputs.
- CGAL: exact predicates and kernel primitives, no epsilons in decisions.
- **Python:** 3.12, pathlib, no `__all__`, Google docstrings, named exceptions, `Cls.build` factories, fail loud.
- **Tests:** failing test first, no skip or xfail.
- **Frozen code:** `engagement_controlled_toolpath`, `radius_regulated_toolpath`, `rho_regulated_toolpath`, `spiral_entry_toolpath`.
- **Measure the thing yourself.** A subagent's number is a lead until you re-run it.
- **Deferred minors** go to the final whole-branch review. They are in the ledger: 0.1 M1, M2, M4–M7 and O1–O6.

## 9. Environment booby traps

- **pixi only, and only in this worktree.** Bare `python` lacks `_stock_2`.
- **Ad-hoc scripts** must bypass the rebuild:
  ```bash
  pixi run bash -c 'd="$(python -c "import sys; print(next(f.path for f in sys.meta_path if hasattr(f, \"known_wheel_files\") and \"compas_cgal._stock_2\" in f.known_wheel_files))")"; SKBUILD_EDITABLE_SKIP="$d" PYTHONPATH=. python <script>'
  ```
  This uses the last built extension. After C++ edits, rebuild first with any `pixi run pytest …`.
- **`pixi run red-manifest`** takes 14–22 min. Run it in the background with its log under `build/`, append `echo "exit=$?"`, and read the `red_manifest` lines at the end. The wrapper's exit code is not the verdict.
- **`.DS_Store` in `external/`** breaks the CMake configure. Run `find external -name .DS_Store -delete && <build>` as one command.
- **Network drops or host sleep** kill background agents mid-task (`API Error … ENOTFOUND` / "computer went to sleep"). Resume the same agent with `SendMessage` to its agent ID; its working-tree changes survive.
- **The IDE's clangd** flags new `.h` files with no compile-DB entry. Those errors are spurious; verify with a standalone compile under the project flags.
- **`pixi run -e docs docs`** has one existing dead-link warning, in `segment_site_mat.md`. It is not yours.
- **Module paths:** the corpus uses `compas_cgal.adaptive.motion/.units`; live uses `compas_cgal.motion/.units`.
- **Scratch** goes under `build/`, never `/tmp`.
- **zsh** does not word-split unquoted variables.

## 10. Suggested file structure

```
src/
  stock_2.cpp                     modify (0.2: authored-radius circle sweep, one capsule model, local depletion)
  stock_local_2.cpp               modify (0.2: capsule zone surgery)
  exact_circle_chart_2.cpp        modify (0.2: exact entry incidence)
  swept_prefix_engagement_2.*     exists (0.1)
src/compas_cgal/
  engagement.py                   modify (0.2: Circle replay uses the authored radius)
  replay_classification.py        modify (0.3 rules; 1.2 HELIX)
  pipeline/                       create in 1.1
benchmarks/
  survey.py, quality.py           exist (0.1); modify per 0.3 rules
  depletion.py, gate.py           modify (0.2, 0.3, 1.1)
  results/gate_vectors.json       create (0.3)
tools/gate_ratchet.py             create (0.3)
tests/benchmarks/test_quality.py  modify (0.2 re-pin; 0.3 pins move to JSON)
tests/tools/test_red_manifest.py  modify whenever docs/red_manifest.json changes
.github/workflows/benchmarks.yml  modify (0.3, with my sign-off)
```

## 11. Risk and budget

The plan sets no time budget. The risks, verbatim from their sources:
- Task 0.2's affordability is "conditional: affordable on a clean stock (3–13 ms per test), not once lines are replayed as chains (≈0.4 s per test after four wall runs)" (census g4).
- The frozen generators deplete with `subtract_capsule_quad` and the surrogate `subtract_arc_sweep_local`. If the replay changes model and they don't, generator and replay disagree on the baselines until those are retired (ruling 0.2-C1 accepts this).

## 12. Start

```bash
cd /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation && git log --oneline -5 && git status --short && tail -30 .superpowers/sdd/2026-10-13-product-pipeline/progress.md
```

Then rule on 0.2-C2 and 0.2-Q1 with the code open, and write them to the ledger. Re-measure the 0.2 draft's numbers: `verify/verify_sweep.py` and the four-run vertex count. Then dispatch the Task 0.2 implementer and the Task 1.1 drafter.

End of prompt.
