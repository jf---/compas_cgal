# Kickoff prompt: Product Pipeline, Task 0.0 → Stage 0 → Stage 1

I'm starting **Task 0.0 and Stage 0 (gate physics), then Stage 1 (narrow-regime pipeline)** of the compas_cgal toolpath project. The work happens in the worktree `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation` on branch `codex/held-task5-continuation`.

The last hard checkpoint is commit `db549852`. On top of it sit:
- an uncommitted 11-file diff from 2026-09-08/09;
- the new plan, traceability matrix and census evidence, all uncommitted and written 2026-10-13.

## 1. What I want you to do

Execute `docs/superpowers/plans/2026-10-13-product-pipeline.md` with **superpowers:subagent-driven-development**, starting at Task 0.0. The plan is the authority. Use **plan-status** as the scope gate.

**Fan-out, in shape:**
- **One implementation workstream.** Tasks run in order, because they share the gate files and the ratchet.
- **Expansion agents, run in parallel.** You expand the active task to step level, per the plan's "rolling wave" rule. Meanwhile one read-only agent per upcoming task (at most the next two) drafts that task's step-level expansion from the sources the plan names. You review and integrate the drafts. Agents never decide scope.
- **Per task:** one implementer agent, then one fresh reviewer agent, then you integrate and run the gate yourself.

**Effort: `opus-medium` for every subagent** (defined in `~/.config/dev/agents/opus-medium.md`, pinning `model: opus` and `effort: medium`).
- The project's failure mode was depth without integration, so agents stay at medium and the orchestrator owns integration.
- Raising one agent's effort needs my explicit go, per task, with a reason. Candidates are Tasks 0.1 and 0.2 (exact C++ kernel ports), and only after a medium attempt fails review.
- Never spawn Fable. Each subagent's first output line must state its model; if one isn't Opus, stop it and tell me.

**Commits** are gated by my global RULES.md. At session start, ask me once: "per-task commits on `codex/held-task5-continuation` for this plan: approve?" Commit per task only after a yes. Author and committer are Jelle Feringa, with no attribution lines.

## 2. Read first, in dependency order

1. `docs/superpowers/plans/2026-10-13-product-pipeline.md`: the whole plan. Pay most attention to the Spine, Rulings R1–R4, Where we are, Global Constraints, Task 0.0 and Stage 0.
2. `docs/superpowers/plans/2026-10-13-product-pipeline-traceability.md`: which census finding each task owns. Read the rows for the task you're on.
3. `docs/superpowers/census/2026-10-13/out/`: one record file per slice. Read the `## Slice summary` of `g4_cpp_kernels.md`, `v_validation_tests_tools_build.md`, `g1_python_generators_and_guide.md` and `g7_unmerged_commits.md`; these four drive Stage 0 and Task 1.1.
4. `docs/superpowers/census/2026-10-13/verify/`: the six orchestrator re-measurement scripts. Their outputs are quoted in the plan's "Where we are".
5. `CONVENTIONS.md` (repo) and `docs/exactness.md`: the CGAL exactness rules, with no epsilons in decisions.
6. Memory, in `~/.config/dev/projects/-Users-jelle-Code-CADCAM-compas-cgal-prs/memory/`:
   - `product-pipeline-rulings-2026-10-13.md` and `toolpaths-over-formalities.md` (the spine);
   - `measure-the-thing-not-a-proxy.md`: re-run agent measurements, never relay them;
   - `subagent-execution-is-the-default.md`;
   - `ds-store-breaks-native-dependency-digest.md`;
   - `shared-worktree-stage-by-pathspec.md`;
   - `rules-beat-system-attribution.md`.

## 3. Boundary contract (your input)

The gate consumes a `ToolpathResult` (`src/compas_cgal/toolpath.py`). Every generator, including the pipeline, emits:

```python
@dataclass
class ToolpathOperation:
    geometry: Line | Arc | Circle          # compas geometry, world XY, z = cut plane or clearance
    operation: OperationType               # CUT, LEAD_IN, LEAD_OUT, LINK, RETRACT, PLUNGE (HELIX arrives in Task 1.2)
    path_index: int
    clockwise: bool = False
    start_tangent: np.ndarray | None = None
    end_tangent: np.ndarray | None = None

@dataclass
class ToolpathResult:
    operations: list[ToolpathOperation]
    polyline: np.ndarray                   # (N, 3) visualisation only
```

The gate entry point is `tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running`. It is parametrised over `benchmarks.gate.GATE_GENERATORS` × pockets × caps. The replay is `benchmarks/depletion.py` plus `benchmarks/survey.py`, through `compas_cgal.engagement._subtract_operation` and `compas_cgal.replay_classification`.

**Epistemic state**

Content identity "live" means `db549852` plus the 11-file diff. "corpus" means `codex/held-reference-corpus` at `a29ff8c6`, which has 4 uncommitted files from another session; do not touch them.

| Claim or failure | Classification | Evidence | Content identity | Consequence |
| --- | --- | --- | --- | --- |
| All 12 product-gate cases fail, with vectors equal to the pin | verified invariant | gate run 2026-10-13 | live | baseline the ratchet starts from |
| Gate scores a straight slot into solid stock at 360° (physical 180°) | verified invariant (model defect) | `verify/verify_direction.py` | live | Task 0.1 changes the measure; every engagement number moves |
| Replayed loop not idempotent (8/12 plunge, 12/12 re-sweep change the stock) | verified invariant (model defect) | `verify/verify_sweep.py` | live | Task 0.2; until then, "removes material" and the redundant count cannot be trusted |
| Tests asserting the 12 reds, identical streams at 120°, entries exceeding the cap, plunge per pass, Arc rejection | **false test/oracle premise** (they pin bad behaviour) | plan Task 0.3 retire list; census `v` finding 4 | live | retire under R1e in Task 0.3; never make a generator satisfy them |
| `engagement_controlled` / `radius_regulated` still bisect the advance | confirmed implementation gap | `engagement_toolpath.py:601-654` | live | frozen baselines; the pipeline uses the ladder (Global Constraint) |
| Analytic Held lap matches Held on Figure 5 (270 loops, 0 ≤ r, 9.16 m) | verified invariant | `verify/verify_analytic.py` | corpus | Stage 2 port source |
| Analytic Held lap raises at sharp convex corners | confirmed implementation gap | `verify/verify_analytic_gate.py`, `verify_analytic_orient.py` | corpus | Task 2.2 |
| Deleted `engagement_ordered_toolpath` cut rect entries from 5 to 1 | verified (the code exists in history); the 5→1 number is relayed from docs | `git show 9fa6eea8^:src/compas_cgal/engagement_ordered_toolpath.py`; `docs/machining_metric_validity.md:136-141` | `9fa6eea8^` | Task 1.1 source; re-measure 5→1 yourself |
| Exact segment-site MAT and `adaptive/` are bound to one L fixture | verified invariant | `segment_site_mat_bundle.cpp:39`, `segment_site_catalog_graph.cpp:172-192` | live | not a guide; harvest parts only |
| TEA-growth lemma in `engagement_2.cpp:407-499` is false for concave voids | unresolved hypothesis: a probe shows a bound violation, no false certificate built | census g4 finding 4 | live | Stage 4.3; nothing consumes `cap_certified` meanwhile |
| The 240° spiral-entry floor holds only for the static model; physically 120° | unresolved hypothesis (analysis, not measured) | census g1 finding 1 | live | re-measure after Task 0.1 |
| `center_domain` about 3,500× slower on oblique edges | unresolved hypothesis (mechanism unknown) | `docs/oblique_edge_cost.md` | live | Task 1.3 cost bullet |
| Swept-prefix in-process engagement exists only on the corpus branch (`src/swept_prefix_engagement_2.cpp`, commit `7be75e0d`) | verified (code exists); its port conflicts with live's rewrite of `engagement_2.cpp` / `audit_exact_station_2.cpp` | census g7 finding 4 | corpus | Task 0.1 port; expect manual merge |

## 4. Deliverables (this session)

- **Task 0.0:** the 11-file diff committed after `red-manifest` passes; one docs commit with the plan, the traceability matrix, `docs/superpowers/census/2026-10-13/`, this kickoff, and status headers on the number-type stage 0–3 plans and the Held paper-reproduction plan.
- **Task 0.1:** in-process engagement in `benchmarks/survey.py` and in the generators' probe ring, ported from corpus `7be75e0d`.
- **Task 0.2:** exact idempotent sweeps (`subtract_circle_sweep` with the authored radius, one capsule model shared by generator and replay, exact entry incidence, local depletion).
- **Task 0.3:** the gate rules (a)–(f), the retired pins, `benchmarks/results/gate_vectors.json`, `tools/gate_ratchet.py`, CI fixes, and the per-stage exit test.
- **Then Stage 1** (Tasks 1.1–1.7), as far as the session reaches.

## 5. Exit gates (quoted from the plan)

- **Task 0.1:** "a straight slot into solid stock reads 180° ± one sample spacing; the 12 gate vectors are re-pinned; every engagement criterion (cap exceedances, max step, slotting) moves on at least one pocket, and the uncut, gouge, tangent and continuity criteria stay unchanged."
- **Task 0.2:**
  - "re-cutting a loop from any start and plunging at its own entry change nothing (12/12)";
  - "the arrangement after four 16 mm runs stays under 200 vertices, down from 8,014";
  - "the gate's redundant-operation counts are re-pinned";
  - "generation plus gate replay on each gate pocket takes no more than 2× the 2026-10-13 time."
- **Task 0.3:** "the old generators' vectors are re-pinned under the new rules; the ratchet passes on this commit and fails on a synthetic regression."
- **Task 1.1 (120°, all three pockets):** "plunges into material: 1 per pocket; redundant 0, degenerate 0; tangent breaks ≤ rho's 0/0/2; no criterion worse than rho_regulated's 2026-10-13 vector under the Stage 0 gate."

**If a gate slips:** stop. Write `docs/superpowers/state/<task>-gate-analysis.md` (what failed, the measured evidence, and the smallest plan amendment). Bring it to me before changing scope.

## 6. Decoupling rule

Each Stage 0 task is checked on its own synthetic probe before the full gate:
- Task 0.1: the straight-slot probe in `verify/verify_direction.py`;
- Task 0.2: the idempotence probe in `verify/verify_sweep.py`;
- Task 0.3: a synthetic regression vector for the ratchet.

The full product gate (`pixi run pytest tests/benchmarks/test_quality.py -k worth_running -n auto`, about 5 min) is the integration check, run by you, not by the implementer.

## 7. Load-bearing conventions

- The spine:
  - **every commit moves a gate number on a gate pocket**; setup ships inside the task that moves the number;
  - **integrate before inventing**, citing the source part named in the plan;
  - **no bisection**;
  - **exact decisions only on exact models**;
  - **forced decisions are structured outputs**.
- CGAL: exact predicates and kernel primitives, no epsilons in decisions (`CONVENTIONS.md`, `docs/exactness.md`).
- Python 3.12; pathlib; no `__all__`; Google docstrings; named exceptions; `Cls.build` factories; fail loud, no fallbacks.
- Failing test first. No skip or xfail. Editing reference tests happens only under R1e in Task 0.3.
- Frozen code, which nothing may change: `engagement_controlled_toolpath`, `radius_regulated_toolpath`, `rho_regulated_toolpath`, `spiral_entry_toolpath`. Removing them needs a ruling.
- Measure the thing yourself. An agent's number is a lead, not evidence, until you re-run it.

## 8. Environment booby traps

- **pixi only.** Run it in this worktree, whose env is live, and in the corpus worktree. The main checkout `/Users/jelle/Code/CADCAM/compas_cgal_prs` has no usable env, and bare `python` lacks `_stock_2`.
- **Ad-hoc scripts** must bypass the editable rebuild:
  ```bash
  pixi run bash -c 'd="$(python -c "import sys; print(next(f.path for f in sys.meta_path if hasattr(f, \"known_wheel_files\") and \"compas_cgal._stock_2\" in f.known_wheel_files))")"; SKBUILD_EDITABLE_SKIP="$d" PYTHONPATH=. python <script>'
  ```
- **`.DS_Store` in `external/`** breaks the native dependency digest at CMake configure (`NATIVE_DEPENDENCY_INTEGRITY_ERROR package=cgal`). Fix it in one command: `find external -name .DS_Store -delete && <build>`.
- **Module paths differ by branch.** On the corpus branch `EngagementCap` and the units live in `compas_cgal.adaptive.motion` / `compas_cgal.adaptive.units`; on live, in `compas_cgal.motion` / `compas_cgal.units`. Ports must re-point the imports.
- **Shell:** zsh does not word-split unquoted variables, so `F="a b"; cmd $F` passes ONE argument. Write paths out literally.
- **Scratch** goes under `build/` (gitignored), never `/tmp`.
- **Commands that run long** (the full gate, builds) get an explicit `timeout` of up to 600000 ms, or a log redirected under `build/`.

## 9. Suggested file structure

```
src/compas_cgal/
  engagement.py                  modify (0.2: Circle replay uses the authored radius)
  replay_classification.py       modify (0.3 rules; 1.2 HELIX)
  pipeline/                      create in 1.1 (__init__.py minimal)
    guide.py  narrow.py  linking.py  generate.py  result.py
src/
  swept_prefix_engagement_2.cpp  create in 0.1 (port from corpus 7be75e0d) + binding
  stock_2.cpp                    modify (0.2: one capsule model, local depletion)
benchmarks/
  survey.py  quality.py  quality_observations.py  depletion.py  gate.py   modify (0.1-0.3, 1.1)
  results/gate_vectors.json      create (0.3)
tools/gate_ratchet.py            create (0.3)
tests/benchmarks/test_quality.py modify (0.3: pins move to gate_vectors.json)
tests/benchmarks/test_quality_transition.py, test_gate_attribution.py (parts)   retire (0.3, R1e)
.github/workflows/benchmarks.yml modify (0.3: smoke step, ratchet)
docs/machining_quality.md        modify (end of Stage 0: the gate's new rules)
```

## 10. Risk and budget

The plan sets no time budget. Its risks, verbatim from its sources:
- the Task 0.1 port "conflicts heavily with live's rewrite of `engagement_2.cpp` and `audit_exact_station_2.cpp`" (census g7);
- Task 0.2's affordability is "conditional: affordable on a clean stock (3–13 ms per test), not once lines are replayed as chains (≈0.4 s per test after four wall runs)" (census g4).

## 11. Start

```bash
cd /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation && git status --short && pixi run lint
```

Then ask me the commit-approval question from section 1, run `pixi run red-manifest` (Task 0.0 step 1), and report the result before any commit.

End of prompt.
