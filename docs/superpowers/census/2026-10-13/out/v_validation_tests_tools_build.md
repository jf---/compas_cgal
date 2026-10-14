# Census slice v_validation_tests_tools_build (2026-10-13)

Structural slice: 304 files (benchmarks outside Held, tests, tools, CI, build). Read in full: the gate definition (`benchmarks/quality.py`, `quality_observations.py`, `survey.py`, `depletion.py`, `gate.py`, `coverage.py`, `pathmetrics.py`, `exceedance.py`, `spec.py`, `benchmarks/families/*`), `tests/benchmarks/test_quality.py`, `test_gate_attribution.py`, `test_quality_transition.py`, the CI workflows and `pyproject.toml`. Also read `src/compas_cgal/replay_classification.py`, which belongs to another slice, because the gate depends on it. Every other file was read structurally: docstring, public API or test names, generator imports, and any assertion that pins generator behaviour. No code was run.

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **The gate punishes keeping the tool down and does not count retracts. It pushes generators toward the opposite of the world-class bar.** (verified)
   - `replay_classification.py:18-25,99-103`: a `LINK`, `LEAD_IN` or `LEAD_OUT` at cut depth is classified `"motion"`.
   - `survey.py:363-382`: `removes_material` is `not probe.exactly_equals(stock)`.
   - `quality_observations.py:904`: every motion that removes nothing is counted as a redundant operation, and the required count is 0.
   - So a stay-down link through cleared stock fails the gate, and so does a lead-in arc in cleared space. A retract → clearance `LINK` → plunge costs nothing: `retract_count` and `reentry_count` are reported only (`quality.py:320-321`, `test_quality.py:154-176`).
   - **Physically wrong.** Commercial adaptive clearing keeps the tool down through cleared material because retracts and plunges cost time and add entry shocks.
   - **Fix:** split "redundant" into two things.
     - (a) A cut motion whose swept region is exactly disjoint from the stock. That is a stay-down link: legal, and counted as air-at-depth length.
     - (b) A loop that removes nothing. Keep that one as a defect.
   - Then gate the retract and plunge counts against the number of connected regions.
2. **The "slotting" criterion is not slotting.** (verified, `survey.py:80-89,284,348-352,381`)
   - Any straight cut motion that reaches `0.5 × cap` of rim arc at any sample counts as slotting.
   - At the 40° cap the threshold is a 20° rim arc, i.e. `a_e = r(1−cos 10°) ≈ 0.015 r`. At the 120° cap it is `0.134 r`, while loops may legally take `0.5 r`.
   - Real slotting is full immersion (360° rim arc in this kernel's convention, `test_quality.py:280`).
   - **Physically wrong for a world-class path.** Straight segments at controlled engagement are standard in Fusion Adaptive and VoluMill. This rule forces every cut onto a circle and makes corner cleanup (plan B2a, D2b) much harder than the physics requires.
   - **Fix:** judge straight motions against the same exact cap as every other motion.
3. **The gate cannot represent a helical or ramped entry.** (verified)
   - `replay_classification.py:93-96`: any line with a z change and XY travel raises `CutPlaneRampError`.
   - `:113-114`: any engaged curve off the cut plane raises `OffPlaneReplayCurveError`.
   - `tests/benchmarks/test_quality_observations.py:439-451` pins that a ramp must raise.
   - `tests/test_replay_classification.py:49-57,74-76` pins both of those raises.
   - The only entry the gate can express is a vertical plunge that bores a disk (`depletion.py:106-110`). The spiral-entry tests pin a 240° floor after a single plunge (`test_engagement_spiral_entry_toolpath.py:127-139`).
   - **So `cap exceedances = 0` and `max engagement step ≤ cap` cannot be met on any pocket by any planar generator.** This is not a generator defect.
   - Plan D1 is right that the vocabulary must grow. It understates the change: the replay classification, its pinning test, and the survey/snapshot partition validation (`quality_observations.py:773-833`) all have to change together.
4. **Green tests pin today's bad behaviour, so they fight a better generator.** (verified)
   - `tests/benchmarks/test_quality_transition.py:16-31`: a permanently green test asserts that the 12 gate cases fail with exactly the reviewed vectors.
   - `test_quality.py:1312-1316` asserts the same pin inside the product test, before the verdict. A partial improvement therefore fails as a pin mismatch and loses the diagnostic report.
   - `_expected_quality_gate_violations` (`test_quality.py:1280-1294`) always emits the uncut, degenerate, redundant, cap, step and tangent lines. A generator that fixes any one of them breaks the pin.
   - `test_gate_attribution.py:164-169` requires the two generators to emit identical streams at 120°. `:97-115` rejects `Arc` geometry outright, so a gated generator that emits G2 arcs fails.
   - `tests/test_engagement_toolpath.py:148-170` requires every chain-entry loop to exceed the cap.
   - `tests/test_engagement_radial_toolpath.py:574-588` requires more plunges than chains and a retract after every pass. `:408-431` requires the virgin-stock entry to be the maximal full-slot circle.
   - `tests/benchmarks/test_pathmetrics.py:131-141` pins a 360° raw maximum on the entry (unregulated generator, harmless).
5. **Uncut is measured over tool-reachable material, so convex corners require the cutter centre to reach the apex.** (verified, `coverage.py:11-13,125-152`; `test_quality.py:379-396`)
   - A G1 curve cannot pass through the apex of a convex cone without reversing, and reversals count as tangent breaks (`quality_observations.py:919-921`). The apex is therefore reachable only as a motion endpoint followed by a retract.
   - Together with rule 2, closing a corner inside one operation may be infeasible at the 40° cap. That is unverified (not measured).
   - **Physically:** tool-reachable is the right target for a whole job. Commercial adaptive leaves corner rest material to a rest or finishing operation.
   - **The gate conflates operation and job.** It has no notion of a second operation with its own criteria. Plan D2(a) points the right way, but the ruling should be "multi-operation programs, each gated", not "shrink the denominator".
6. **Degenerate loops are defined only for full circles.** (verified, `survey.py:372`; `quality_observations.py:896-903`)
   - An `Arc` with ρ ≤ r is never degenerate, so a generator can split a bore into arcs and dodge the rule.
   - A ρ ≤ r loop in virgin stock is already caught by the cap. In partly cleared corner stock, a small-radius circle under the cap is a legitimate cleanup move that the rule bans.
   - **Partly wrong physics.** Keep the "annulus has no hole" fact as a report, and let the cap decide.
7. **The world-class items "short path" and "few retracts" are not gated at all.** (verified, `test_quality.py:154-176`)
   - Path length, air, retracts, plunges and cycle time are reported only.
   - `benchmarks/models.py` `MachineModel` with `quality.machine_outcome` (`quality.py:732-763`) already computes a curvature-limited feed time, and nothing uses it as an objective.
8. **"No air-cutting loops" is only weakly enforced.** (verified)
   - `removes_material` is exact, so a loop that removes a sliver passes even if 99% of it is air.
   - `recut_fraction` is reported, not gated (`quality.py:184-191`). `immersion_at_design_fraction` (`quality.py:933-952`) is the right "is it cutting at the commissioned load" measure, and it is reported only.
9. **Validation holes the plan relies on.** (verified)
   - Continuity, tangent breaks and engagement steps are evaluated only between operations with consecutive indices (`quality_observations.py:698-699,914`). Nothing checks that a plunge lands where the next cut starts, or that rapids connect.
   - Nothing counts plunges into material. The plan's `TreeWalkResult.virgin_plunges` is self-reported by the generator.
   - The engagement-step criterion compares per-motion peaks (`:926-928`), not the junction values.
10. **The gate pockets are too narrow to show a world-class path.** (verified, `gate.py:64-137`)
    - Three axis-aligned, integer-coordinate polygons, one tool (D = 2), two caps.
    - No islands, necks, curved walls, acute or obtuse corners, open sides, or large pockets where ordering dominates.
    - Families for most of these already exist (`families/topology.py`, `necks.py`, `analytic.py` stadium/disk/arc_channel, `complexity.py` n-gon/arc_fraction), but they are scored only by the certifier-cost corpus. That corpus runs the unregulated `trochoidal_mat_toolpath_circular` (`runner.py:33-49`), never a gated generator.
    - Curved walls exist only as tessellations, 64 segments per cap (`analytic.py:20`). `complexity.py:22-23` records that at equal area, k = 3→64 grows operations from 47 to 767: the straight skeleton explodes on tessellated curves (unverified number, quoted from the docstring).
11. **CI never measures generation.** (verified)
    - `.github/workflows/benchmarks.yml` runs on push/PR to `main`/`master` only. It runs `pytest tests/benchmarks`, which contains the 12 deliberately red gate cases, so the job is red by construction. It does not use `red-manifest`.
    - Its smoke step calls `python -m benchmarks.cli --corpus smoke`, but the CLI needs the `corpus --name smoke` subcommand (`cli.py:284-292`), so the step errors.
    - The segment-site MAT, task3 and audit native gates are not in CI.
    - `build.yml` uses conda and Python 3.10, while `wheel.py-api = cp312`. The cibuildwheel `test-command` runs the whole of `tests/`, including the red gate (`pyproject.toml` `[tool.cibuildwheel]`).

### Generation capabilities in this slice that the gated generators do not use

`gate.py:38-39` imports only `engagement_controlled_toolpath` and `radius_regulated_toolpath`.

- `rho_regulated_toolpath`: no bores; tangent bridges at any radius difference (`tests/test_engagement_rho_toolpath.py:145-201`); `DeclinedRegion` reporting (`:222-263`). Tested, not gated.
- `spiral_entry_toolpath`: centre-plunge spiral ramp that reaches the 240° floor (`tests/test_engagement_spiral_entry_toolpath.py`). Not gated.
- The certified `compas_cgal.adaptive` stack:
  - exact segment-site MAT (`tests/adaptive/test_medial_axis.py`);
  - neck inventory (`test_neck.py`);
  - pre-cleared entry bore (`test_entry.py:117-259`);
  - middle-curve circle proposals along a MAT edge (`test_middle_curve_path.py`);
  - reachable domain (`test_reachable_domain.py`);
  - route-retrace transactions.
  - There is no bridge to `ToolpathResult`, so it cannot be gate-scored: grep for `ToolpathResult|survey_path|measure_quality` in `tests/adaptive` and `src/compas_cgal/adaptive` returns nothing. Two of its tests are manifest reds (`docs/red_manifest.json`).
- Held native kernels, tested at top level, are the contour-aware placement primitives:
  - boundary-normal circle (`tests/test_boundary_normal_circle.py`);
  - disk contour correction (`test_held_disk_contour.py`);
  - whole-orbit leading-rim upper bound (`test_held_contour_bound.py`);
  - local circle deletion that preserves material (`test_circle_redundancy.py`);
  - analytic curved-boundary medial contacts (`test_native_boundary_medial.py`);
  - exact line/arc boundary construction (`test_native_boundary_curve.py`).
- Fast exact depletion: quad capsule subtraction (`tests/test_stock_capsule_quad.py`) and local depletion (`test_stock_local_depletion.py`). Relevant to generator-internal stock cost in plan Task 3.
- Time-based objective: `MachineModel.curvature_limited_feed_mm_per_s` plus `machine_outcome` (cycle time as a lower bound).
- The pocket families listed in finding 10, plus `external.py`, a loader for third-party profiles.
- The MATHSM baseline (`benchmarks/mathsm.py`): constant spacing plus a brute-force cap search. The gate never compares a generator's length against it.

### Where this slice shows the draft plan is wrong or incomplete

1. **Design, "reposition above the stock".** This optimises toward finding 1, a gate distortion. Every non-adjacent tree edge becomes a retract and a plunge, which the gate does not count, and the world-class bar penalises them. Rule on D3 (stay-down links, gate rule 1a above) before Stage A. Otherwise Stage A's linking is rework.
2. **"Productive-only emission" with `exactly_equals`.** This removes zero-removal loops, but keeps sliver loops that are mostly air (finding 8). "0 redundant" is not "no air-cutting loops". The plan should report per-loop removed area against the swept annulus, or the share at design immersion, and say which it optimises.
3. **D1 (helix).** The gate today raises on ramps and off-plane arcs (finding 3). The ruling must cover:
   - `replay_classification.py`;
   - `test_quality_observations.py:439` and `tests/test_replay_classification.py:49-57,74-76` (the pins that ramps and off-plane curves raise);
   - the survey partition in `quality_observations.py:773-833`;
   - a z-aware swept area.
4. **B2a/D2.** Judging the wall run-in at `cap/2` inherits a rule with the wrong physics (finding 2). Recommend one ruling: straight motions use the full cap. With that, D2(b)'s corner exemption is unnecessary, and B2a is more likely to close the corners inside the gate. Whether it does is unverified.
5. **Task 4, adding `tree_walk` to the gate.** It must also update:
   - `test_gate_attribution.py`, which hard-codes two generators and rejects `Arc` (`:97-115`);
   - `test_quality_transition.py` and the pin table at `test_quality.py:1258-1271`;
   - `docs/red_manifest.json`, whose regexes name only the two generators.
6. **Generality is not shown.** The product gate stays at three pockets. The plan's review focus (island, dumbbell) is pinned only in the unit tests for `skeleton_guide`. Score the tree walk on `families.topology.island_grid`, `families.necks.dumbbell` and `families.analytic.stadium`/`disk` as well. On tessellated curves the straight-skeleton tree has one leaf per tessellation vertex (finding 10). That is the case for the exact segment- and arc-site MAT (`adaptive.medial_axis`, native boundary medial) as the guide.
7. **`virgin_plunges` is generator-reported.** Add a gate-side exact count: plunge disks that intersect the stock at the moment of the plunge. Also add a check that each plunge's XY equals the next motion's start (finding 9).

### Answers to the slice questions

**(a) Gate rules and their physics.**

| Rule | Where | Physically right? |
|---|---|---|
| Redundant = any cut-plane motion with an exactly unchanged stock, LINK and leads included | `quality_observations.py:904`, `survey.py:382`, `replay_classification.py:18-25` | Wrong for links and leads; right for loops |
| Slotting = straight motion over cap/2 rim arc | `survey.py:89,351` | Wrong (not slotting; bans controlled straight cuts) |
| Degenerate loop = full circle with ρ ≤ r, exact | `quality_observations.py:896-903` | The physics fact is right; gating it is redundant with the cap, and arcs are exempt |
| Uncut over tool-reachable material, grid ≤ 0.1 r, required exactly 0 | `coverage.py:41-49,125-152` | Right for a job; too strict for a single roughing operation with no rest operation |
| Entries count as cap exceedances | `quality.py:218-221`, `survey.py:347` | Right physics, but unfixable because the gate cannot express a helix or ramp |
| Tangent break = dot < 1 − 1e-4 (≈0.81°), reversals included | `quality.py:132-135`, `quality_observations.py:918-923` | Right, slightly strict; achievable with exact tangent bridges |
| Engagement step = peak-to-peak between index-adjacent motions, ≤ cap | `quality_observations.py:926-934` | Weak proxy: peak, not junction value |
| Loop radius step ≤ 2 r within a run; runs split at rapids and at `path_index` changes | `quality_observations.py:848-878`, `test_quality.py:134-145` | Harmless. `path_index` is generator-controlled, so it can be gamed |
| Continuity tolerance 1e-3 r between index-adjacent motions only; nothing across rapids or plunges | `quality.py:125-130` | Hole: plunge → cut teleports are unchecked |

**(b) Which generator each test exercises.** The record per file below gives the generator. The tests that pin bad behaviour are listed in finding 4.

**(c) The gate pockets.** They are not enough to show a world-class path (finding 10). Missing:
- islands;
- necks narrower than a trochoid;
- true curved walls (only tessellations exist);
- acute and obtuse convex corners;
- open-sided and slot features;
- large pockets where linking and ordering dominate;
- non-integer CAD coordinates.

**(d) CI and build.** See finding 11. Nothing in CI runs a generator against the gate on a work branch. Two more defects:
- dead pixi tasks: `schema` targets `tests/adaptive/test_schema.py` and `mutations-adaptive` targets `scripts/run-adaptive-mutations.py`, and both files are missing;
- `cli.py:26` advertises `pixi run measured-run`, which does not exist.

None of these blocks local generator work, since `pixi run baseline` and `red-manifest` work locally.

---

### .editorconfig | lines 26 | read structural
- **role:** build
- **used by:** editors
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### .github/ISSUE_TEMPLATE/bug_report.md | lines 29 | read structural
- **role:** doc
- **used by:** GitHub
- **state:** live (upstream compas_cgal template)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### .github/ISSUE_TEMPLATE/feature_request.md | lines 22 | read structural
- **role:** doc
- **used by:** GitHub
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### .github/workflows/benchmarks.yml | lines 52 | read 1-52
- **role:** build
- **used by:** GitHub Actions, push/PR to main/master
- **state:** live
- **generation:** n/a. It runs no generator against any quality bar that can turn green.
- **defects:**
  - `:45` `pytest tests/benchmarks` includes the 12 deliberate gate reds and does not use `red-manifest`, so the job is red by construction.
  - `:51` `python -m benchmarks.cli --corpus smoke ...`: the CLI requires the `corpus --name smoke` subcommand (`benchmarks/cli.py:284-292`), so argparse errors.
  - Runs only on main/master; work branches never run it.
- **unused capability / ideas:** run `pixi run red-manifest` instead of raw pytest; add `task9-mat-compile-gate` and `audit-native`.
- **claims:**
  - `:3-6` "existing build workflow ... runs no tests at all": contradicted. The `[tool.cibuildwheel] test-command` runs `pytest {project}/tests`.

### .github/workflows/build.yml | lines 89 | read 1-89
- **role:** build
- **used by:** GitHub Actions, push/PR to main
- **state:** live (upstream)
- **generation:** n/a
- **defects:**
  - conda + Python 3.10 (`:17-29`) against the project's pixi + 3.12 (`pyproject.toml` `wheel.py-api = "cp312"`).
  - cibuildwheel runs the whole test suite, including the red gate.
  - Windows is in the matrix, but pixi lists only osx-arm64 and linux-64.
- **unused capability / ideas:** none
- **claims:** none

### .github/workflows/docs.yml | lines 32 | read 1-32
- **role:** build
- **used by:** GitHub Actions, push to main
- **state:** live
- **generation:** n/a
- **defects:** `:26-30` uses pip, bypassing the pixi `docs` environment. Not generation-relevant.
- **unused capability / ideas:** none
- **claims:** none

### .github/workflows/pr-checks.yml | lines 21 | read 1-21
- **role:** build
- **used by:** GitHub Actions, PRs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### .github/workflows/release.yml | lines 96 | read 1-96
- **role:** build
- **used by:** tag push v*
- **state:** live
- **generation:** n/a
- **defects:** cibuildwheel test step inherits the red gate (see build.yml).
- **unused capability / ideas:** none
- **claims:** none

### .gitignore | lines 130 | read structural
- **role:** build
- **used by:** git
- **state:** live
- **generation:** n/a
- **defects:** none found (`build/` is ignored at `:9-11`, so census outputs stay untracked)
- **unused capability / ideas:** none
- **claims:** none

### CMakeLists.txt | lines 537 | read structural
- **role:** build
- **used by:** scikit-build-core, pixi gate tasks
- **state:** live
- **generation:** builds `_toolpath` (`src/toolpath.cpp`, the gated generators' skeleton guide), `_stock_2`, `_coverage_2`, `_medial_axis_2` and others.
  - Native gates are `EXCLUDE_FROM_ALL`: `task3_algorithm_gate`, the `exact_*` gates, `audit_native_gate`, `task9_mat_compile_gate`.
- **defects:** none found at structural depth. Note that `task9`/`audit` gates run only through local pixi tasks, not CI.
- **unused capability / ideas:** none
- **claims:** none

### LICENSE | lines 164 | read structural
- **role:** doc
- **used by:** packaging
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/__init__.py | lines 0 | read structural
- **role:** build
- **used by:** package marker
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/acceptance/test_held_coverage_replay_budget.py | lines 60 | read structural
- **role:** test
- **used by:** `pixi run held-coverage-corpus` (not under `testpaths`)
- **state:** parked (exact coverage replay set aside 2026-09-08 per the task description)
- **generation:** n/a; exercises the Held Figure 5 path.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "400 circles replayed in 5.4 s" is unverified.

### benchmarks/acceptance/test_held_full_coverage.py | lines 82 | read structural
- **role:** test
- **used by:** `held-coverage-corpus`
- **state:** parked
- **generation:** n/a; Held paths on four pockets.
- **defects:** none found
- **unused capability / ideas:** an exact full-motion residual gate, an alternative to grid uncut, if made affordable.
- **claims:** none checked

### benchmarks/acceptance/test_held_thinning_coverage.py | lines 47 | read structural
- **role:** test
- **used by:** `held-coverage-corpus`
- **state:** parked
- **generation:** n/a; Held contour-bound thinning.
- **defects:** none found
- **unused capability / ideas:** coverage-preserving circle thinning (an anti-air-loop tool)
- **claims:** none

### benchmarks/acceptance/test_held_toolpath_performance_slow.py | lines 18 | read structural
- **role:** test
- **used by:** `pixi run held-toolpath-gate`
- **state:** live (slow tier)
- **generation:** n/a; Held upper and Monstera generators against the baseline.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "about ten minutes" is unverified.

### benchmarks/cli.py | lines 333 | read structural
- **role:** tool
- **used by:** `pixi run figures`; benchmarks.yml (wrong args)
- **state:** live
- **generation:** n/a; the corpus uses `runner.generate_toolpath`, which is the unregulated generator.
- **defects:**
  - `:26` documents `pixi run measured-run`, which is not a task.
  - The CI calls it with `--corpus`, which it does not accept.
- **unused capability / ideas:** none
- **claims:** `pixi run measured-run` is contradicted: no such task in `pyproject.toml` `[tool.pixi.tasks]`.

### benchmarks/congruence.py | lines 98 | read structural
- **role:** validation
- **used by:** `tests/benchmarks/test_congruence.py`
- **state:** live
- **generation:** n/a. Exact Pythagorean rotations, for invariance tests.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/coverage.py | lines 265 | read 1-265
- **role:** validation
- **used by:** `quality_observations.reduce_quality_evidence`, `quality.py`, `tests/benchmarks/test_quality.py`, `test_coverage_grid.py`
- **state:** live
- **generation:** n/a. Gate rule: uncut is grid-sampled over `ReachableMaterialPredicate2`, i.e. tool-reachable material (`:125-152`). The cell must be ≤ 0.1 r (`:49,193-198`). The gate requires exactly 0.
- **defects:**
  - `:11-13,125` the denominator is tool-reachable. A convex-corner apex is reachable only as a motion endpoint, so with the tangent and slot rules the gate demands a corner operation it has no vocabulary for (finding 5).
  - `:244-265` `_distance_to_boundary` treats walls as segments only. That is fine for polygons, but curved boundaries do not exist in `PocketSpec`.
- **unused capability / ideas:** `wall_scallop_height` is reported. A per-region residue attribution (which corner) would let a generator target its rest material.
- **claims:**
  - "every number a LOWER BOUND on residue": verified by construction (sampling).
  - "denominator is the reachable material": verified (`:125-152`).

### benchmarks/data/held_pfeiffer_2025/README.md | lines 16 | read structural
- **role:** doc
- **used by:** Held reference cases
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/degeneracy.py | lines 170 | read structural
- **role:** validation
- **used by:** `tests/benchmarks/test_degeneracy.py`, corpus
- **state:** live
- **generation:** n/a. Records that a slot exactly 2r wide makes the unregulated generator emit zero operations.
- **defects:** none found
- **unused capability / ideas:** degenerate instances (cocircular square, half-turn arm) could stress a new guide.
- **claims:** "405 vs 525 operations" is unverified.

### benchmarks/depletion.py | lines 204 | read 1-204
- **role:** validation
- **used by:** `survey.py`, `exceedance.py`, `pathmetrics.py`, `toolpath_coverage.py`; the plan's Design says the generator should deplete exactly like `replay_cuts`.
- **state:** live
- **generation:** n/a. Replay model:
  - a plunge bores one tool disk at its foot (`:106-110`);
  - a cut-plane motion subtracts its sweep (`:111-112`);
  - rapids are skipped.
- **defects:**
  - `:151-163` delegates to `classify_line_replay`, which raises on any ramp. The gate cannot replay a helical or ramped entry (finding 3).
  - A plunge into virgin stock is never flagged as such.
- **unused capability / ideas:** the gate could count virgin plunges here exactly (`stock` before `subtract_disk`).
- **claims:** "retracts and clearance-height links are rapid travel; plunges bore a disk": verified `:104-112` and `replay_classification.py:93-103`.

### benchmarks/errors.py | lines 269 | read structural
- **role:** validation
- **used by:** all benchmark modules
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/exceedance.py | lines 185 | read 1-185
- **role:** validation
- **used by:** `pathmetrics.py`, `runner.py`/`report.py` (`truly_exceeding`), `test_exceedance.py`
- **state:** live
- **generation:** n/a. Twelve samples per motion, exact predicate, no gap closure.
- **defects:** none found. Its sample count (12) differs from the gate's (45, `survey.py:72`), so corpus "truly exceeding" and gate "cap exceedances" are not comparable numbers.
- **unused capability / ideas:** `exceedance_positions` attributes over-cap positions spatially. It is the tool plan B3 asks for.
- **claims:** "count unchanged from 12 through 60" (`:54-58`) is unverified.

### benchmarks/external.py | lines 206 | read structural
- **role:** validation
- **used by:** `cli.py corpus --name external`, `test_external.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** load real CAD profiles (e.g. a Fusion 360 Gallery conversion) as gate pockets, not just certifier-cost pockets.
- **claims:** none checked

### benchmarks/families/__init__.py | lines 0 | read structural
- **role:** build
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/families/analytic.py | lines 194 | read 1-194
- **role:** validation
- **used by:** `gate.py` (`rectangle`), corpus, `test_analytic.py`, `test_quality.py`
- **state:** live
- **generation:** n/a. Provides disk, rectangle, stadium and arc_channel pockets, with closed-form clearance.
- **defects:**
  - `:20` curved walls are 64-segment tessellations, so true arcs never reach a generator.
  - `:24` `MIN_CHANNEL_CLEARANCE_FACTOR = 1.2` admits channels with half-width 1.2r to 2r, where no trochoid with ρ > r fits. These are useful as hard cases, but not marked as such.
- **unused capability / ideas:** stadium, arc_channel and disk are ideal gate pockets for curved-wall constant-engagement behaviour (Held's middle curve). None is in the gate.
- **claims:** "64 segments keep sagitta below 0.13%": `1 − cos(π/64) ≈ 0.0012`, verified arithmetic.

### benchmarks/families/complexity.py | lines 222 | read 1-222
- **role:** validation
- **used by:** corpus, `test_complexity.py`
- **state:** live
- **generation:** n/a. Regular n-gon at fixed area; arc-fraction pockets.
- **defects:** none found
- **unused capability / ideas:** use the n-gon sweep (acute and obtuse convex corners) and arc_fraction (bulged walls with reflex corners) as gate pockets for corner and curvature behaviour.
- **claims:** `:22-25` "k 3→64 operations 47→767; arc ratio 0→1 operations 191→1191" (unregulated generator) is unverified. If true, it is evidence that straight-skeleton guides scale with tessellation and not with geometry.

### benchmarks/families/necks.py | lines 243 | read 1-243
- **role:** validation
- **used by:** corpus, `test_necks.py`
- **state:** live
- **generation:** n/a. A 20x10 dumbbell with a wedge neck (`:98-142`), pinch ≥ D.
- **defects:** none found
- **unused capability / ideas:** the dumbbell is exactly plan review-focus item 2 (a neck narrower than a trochoid). It should be a gate pocket, not just a unit fixture.
- **claims:** the measured table `:28-31` is unverified.

### benchmarks/families/precision.py | lines 201 | read 1-201
- **role:** validation
- **used by:** corpus, `test_precision.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** non-integer jittered n-gons, a test of CAD-like coordinates. The gate pockets are integer-only.
- **claims:** "exact-kernel cost driven by construction depth, not precision" is unverified (measured 2026-08-21 on this kernel; it may not hold for Epeck/CORE lanes).

### benchmarks/families/topology.py | lines 111 | read 1-111
- **role:** validation
- **used by:** corpus, `test_topology.py`
- **state:** live
- **generation:** n/a. A 40x30 pocket with a grid of 3x3 islands; channels ≥ 1.2D (`:47,70-75`).
- **defects:** none found
- **unused capability / ideas:** islands are absent from the product gate. `island_grid(1,1)` is the plan's review-focus item 1 (skeleton with a cycle).
- **claims:** the timing table is unverified.

### benchmarks/figure6.py | lines 493 | read structural
- **role:** figure
- **used by:** `cli.py figure6`, `test_figure6.py`, `toolpath_coverage`
- **state:** live
- **generation:** n/a. Compares `engagement_controlled_toolpath` with MATHSM-protocol `trochoidal_mat_toolpath_circular` on the reconstructed Figure 5 polygon.
- **defects:** none found. Its docstring itself notes that the controlled generator's bridges are unregulated.
- **unused capability / ideas:** none
- **claims:** "Figure 6 is the only directly reproducible published result" is unverified.

### benchmarks/figures.py | lines 314 | read structural
- **role:** figure
- **used by:** `pixi run figures`, `test_figures.py`
- **state:** live
- **generation:** n/a. Draws the controlled path against constant spacing (`figure6.controlled_path`, `mathsm.constant_spacing_path`).
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/gate.py | lines 206 | read 1-206
- **role:** validation
- **used by:** `test_quality.py`, `test_quality_transition.py`, `test_gate_attribution.py`, `qualityfigures.py`, `test_engagement_spiral_entry_toolpath.py`
- **state:** live
- **generation:** n/a. It fixes the product: `rect_12x8`, `rect_20x12` and `L_shape` (arm 3D, 12x10), tool 2.0, caps 120°/40°, clearance 2D. The generators are `engagement_controlled_toolpath` and `radius_regulated_toolpath` (`:38-39,203-206`).
- **defects:**
  - The corpus is too narrow (finding 10).
  - `:98,203` the `Generator` type is `PocketSpec → ToolpathResult`, which has no slot for multi-operation programs (rest operation).
  - `:59-62` there is only one tool size.
- **unused capability / ideas:**
  - Add `topology.island_grid`, `necks.dumbbell`, `analytic.stadium`/`disk`, `complexity.regular_ngon(3|5)` and one Held pocket.
  - Add `rho_regulated`/`spiral_entry` as reported (non-product) columns.
- **claims:**
  - "rect_12x8 and L_shape have the same area, 96": verified arithmetic (12·8; 12·6 + 6·4).
  - The `rect_20x12` corner defect (loop 0.0215, link 1.956 at 106.5°) is unverified.

### benchmarks/instrument.py | lines 78 | read structural
- **role:** validation
- **used by:** runner (digits probe)
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/marks.py | lines 376 | read structural
- **role:** figure
- **used by:** `plotting.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/mathsm.py | lines 151 | read structural
- **role:** validation
- **used by:** `figure6.py`, `figures.py`, `test_mathsm.py`, `toolpath_coverage` tests
- **state:** live
- **generation:** MATHSM protocol baseline. It sweeps the unregulated constant-spacing generator and picks the shortest path that meets the cap. Quality: adequate as a baseline.
- **defects:** none found
- **unused capability / ideas:** the gate should report the length ratio against this baseline on every pocket (Held's headline claim). Today only Figure 6 does.
- **claims:** "MC-013 reference sweep is non-monotone" is unverified here (cites `docs/measurement_claims.md`).

### benchmarks/measurement.py | lines 148 | read structural
- **role:** validation
- **used by:** runner, report, necks
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "generation already competitive with the published state of the art" is unverified, and concerns time, not path quality.

### benchmarks/models.py | lines 238 | read structural
- **role:** validation
- **used by:** `quality.py` calibrated outcomes, `test_quality.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `MachineModel` curvature-limited feed gives a cycle-time objective. The gated generators optimise neither time nor length.
- **claims:** "defaults are not measurements": stated, consistent.

### benchmarks/palette.py | lines 176 | read structural
- **role:** figure
- **used by:** plotting, figures
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/pathgeometry.py | lines 263 | read structural
- **role:** figure
- **used by:** plotting
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/pathmetrics.py | lines 231 | read 1-231
- **role:** validation
- **used by:** `gate.py` (`CLEARANCE_Z_TOOL_DIAMETERS`), `figure6`, `mathsm`, `figures`, `test_pathmetrics.py`
- **state:** live
- **generation:** n/a. It defines `entry_cut_indices` (first cut after a plunge) and reports the after-entry maximum.
- **defects:** `:56` measures at `REFERENCE_CAP_DEG = 180`, using a different instrument (the `audit_toolpath_engagement` certifier) from the gate's sampled predicate. Two engagement numbers coexist.
- **unused capability / ideas:** none
- **claims:**
  - "every generator that enters solid stock cuts a FULL SLOT on its first circle": verified for plunge entry by `test_pathmetrics.py:131-141` (unregulated generator).
  - `spiral_entry` contradicts the generality: its first loop is under 360° (`test_engagement_spiral_entry_toolpath.py:138`).
  - "after-entry 98.73–205.73" is unverified.

### benchmarks/performance_record.py | lines 80 | read structural
- **role:** validation
- **used by:** Held benchmarks
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/plotting.py | lines 802 | read structural
- **role:** figure
- **used by:** figures, qualityfigures, `test_plotting.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/quality.py | lines 1125 | read 1-1125
- **role:** validation
- **used by:** `quality_observations.py`, `test_quality.py`, `test_quality_invariants.py`, `qualityfigures.py`
- **state:** live
- **generation:** n/a. Defines the gate constants:
  - `DEGENERATE_LOOP_RATIO = 1.0` (`:87`);
  - `MARGINAL_LOOP_RATIO = 1.5` (`:96`, reported);
  - `CONTINUITY_TOOL_RADIUS_FRACTION = 1e-3` (`:130`);
  - `TANGENT_CONTINUITY_SLACK = 1e-4` (≈0.81°, `:135`);
  - `IMMERSION_STEADY_BAND_FRACTION = 0.10` (`:151`).
  - Chip, immersion and time physics: `:631-911`.
- **defects:**
  - `:543` `marginal_loops` counts circles only.
  - `:320-321` retracts and re-entries are reported, never gated (finding 7).
  - `:184-191,933-952` `recut_fraction` and `immersion_at_design_fraction` are reported only (finding 8).
  - `:507` `measure_quality` imports `quality_observations` lazily to break an import cycle (a design smell, not a generation issue).
- **unused capability / ideas:**
  - `immersion_at_design_fraction` is the best single "is this adaptive clearing" measure and could be a generator objective.
  - `machine_outcome` gives cycle time.
  - `low_chip_thickness_ratio` detects rubbing.
- **claims:**
  - Rim arc = 2× textbook angle: verified by the kernel pin `test_quality.py:262-285`.
  - `cap_exceedances` "includes the entry cut" (`:218-221`): verified (`quality_observations.py:905`).
  - `:261-265` and `test_quality.py:173-175`, "a middle-curve implementation ... does not exist yet": contradicted in part. `src/compas_cgal/adaptive/middle_curve_path.py` (63 lines, commit 621cfb8d, 2026-09-05) emits per-MAT-edge middle-curve circle proposals, and the Held Figure 5 reproduction places centres on a middle curve. There is no whole-pocket generator.

### benchmarks/quality_observations.py | lines 1080 | read 1-1080
- **role:** validation
- **used by:** `quality.measure_quality`, `test_quality.py`, `test_quality_transition.py`, `test_quality_observations.py`, plan Task 0b/4
- **state:** live (canonical reducer)
- **generation:** n/a. The 12 gate criteria and their exact thresholds:
  - uncut 0, counts 0, step ≤ cap, loop step ≤ 2r (`:142-148`);
  - attribution: `:881-953`.
- **defects:**
  - `:904` redundant includes `LINK`/`LEAD_IN`/`LEAD_OUT` (finding 1).
  - `:896-903` degenerate covers circles only (finding 6).
  - `:906` slot comes from the survey (finding 2).
  - `:698-699,914` adjacency is by consecutive operation index, so continuity is never checked across plunges or rapids (finding 9).
  - `:926-928` the engagement step is per-motion peak minus peak, not the junction values.
  - `:226-227,259-260,315-316` thresholds are hard-wired inside the record factories, so every gate ruling (D2, slot) also edits closed-vocabulary validators and `test_quality_observations.py:588-610`.
  - `:855,864-866` loop-step runs split on `path_index`, which the generator controls.
- **unused capability / ideas:** per-operation attribution (`PathQualityAttribution`) is what plan B3/B4 need. It already exists, so B3/B4 can read it instead of writing new attribution code.
- **claims:** "reproduces the twelve decisions" (`:3-5`): verified by `_assert_quality_parity`, `test_quality.py:989-1107`.

### benchmarks/qualityfigures.py | lines 1332 | read structural
- **role:** figure
- **used by:** docs figures, `test_qualityfigures.py`
- **state:** live
- **generation:** n/a. Uses `GATE_GENERATORS` and the unregulated baseline.
- **defects:** none found. At 1332 lines it is over the file-size bar.
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/report.py | lines 126 | read structural
- **role:** validation
- **used by:** cli corpus
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/results/2026-08-28-9c41a7cab375/benchmark_report.md | lines 52 | read structural
- **role:** doc
- **used by:** none (artifact)
- **state:** live (frozen artifact)
- **generation:** n/a. A corpus report on the unregulated generator.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** not re-checked

### benchmarks/runner.py | lines 163 | read structural
- **role:** validation
- **used by:** cli corpus, `test_runner.py`, `test_depletion.py`, `test_exceedance.py`, `test_held_path_report.py`
- **state:** live
- **generation:** `:33-49` `generate_toolpath` = `trochoidal_mat_toolpath_circular`, the UNREGULATED generator. Every corpus family is measured on it.
- **defects:** the corpus never measures a gated or candidate generator, so family evidence says nothing about generation quality.
- **unused capability / ideas:** parameterise the generator so families can score tree_walk/rho.
- **claims:** none

### benchmarks/spec.py | lines 133 | read 1-133
- **role:** validation
- **used by:** all benchmarks, gate
- **state:** live
- **generation:** n/a
- **defects:**
  - `:40-46` the pocket is `Polygon` only, with no arc boundaries, so curved-wall generation can never be gated on true curves.
  - `:23` `MIN_AREA_TOOL_RADII_SQ = 4.0` is weaker than "admits a trochoid" (fine).
- **unused capability / ideas:** none
- **claims:** cap contract (0, 180] matches kernel `certify_segment_tea`. Unverified here (not in slice).

### benchmarks/survey.py | lines 711 | read 1-711
- **role:** validation
- **used by:** `quality.py`, `quality_observations.py`, `test_quality.py`, `test_survey.py`, spiral-entry test, plan Task 0b
- **state:** live
- **generation:** n/a. One replay with:
  - 45 samples per motion (`:72`);
  - the exact cap and slot surrogates (`:283-284`);
  - an exact `removes_material` (`:363-382`);
  - closed-form swept areas (`:470-503`).
- **defects:**
  - `:80-89,348-352` slot = straight motion over 0.5·cap (finding 2).
  - `:372` `loop_radius` is set only for `Circle` (finding 6).
  - `:406-411` arcs and circles must lie on the cut plane, otherwise it raises (finding 3).
  - `:667` plunge swept area = π r² per plunge, regardless of what was already cleared.
- **unused capability / ideas:** none
- **claims:**
  - "45 vs the generator's 32-ring, gcd = 1": the constant is verified (`:72`). The generator ring of 32 is unverified (`engagement_toolpath.py` is not in this slice).
  - "slot threshold is ENGINEERING JUDGEMENT, not a standard" (`:86-88`): verified as self-described. The judgement is wrong (finding 2).

### benchmarks/toolpath_coverage.py | lines 79 | read structural
- **role:** validation
- **used by:** `test_toolpath_coverage.py`
- **state:** live
- **generation:** n/a. Exact residual from emitted planar operations.
- **defects:** none found
- **unused capability / ideas:** an exact (not grid) uncut check for small paths.
- **claims:** none

### benchmarks/units.py | lines 70 | read structural
- **role:** validation
- **used by:** `quality_observations`, survey
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### cmake/NativeDependencies.cmake | lines 210 | read structural
- **role:** build
- **used by:** CMakeLists
- **state:** live
- **generation:** n/a
- **defects:** none found. Per the auto-memory, `.DS_Store` in `external/` breaks the digest; not re-checked.
- **unused capability / ideas:** none
- **claims:** none

### cmake/VerifyNativeSource.cmake | lines 102 | read structural
- **role:** build
- **used by:** CMakeLists, `tests/build/test_native_source_integrity.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### environment.yml | lines 14 | read 1-14
- **role:** build
- **used by:** none in the pixi workflow (legacy conda)
- **state:** superseded-by pixi (`pyproject.toml [tool.pixi]`)
- **generation:** n/a
- **defects:** a conda/pip environment contradicts the pixi-only rule. Harmless if unused.
- **unused capability / ideas:** none
- **claims:** none

### mkdocs.yml | lines 180 | read structural
- **role:** build
- **used by:** `pixi run -e docs docs`, docs.yml
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### pyproject.toml | lines 321 | read 1-321
- **role:** build
- **used by:** pixi, scikit-build, pytest, mypy, ruff, cibuildwheel
- **state:** live
- **generation:** n/a. `testpaths = ["tests"]`, so `benchmarks/acceptance` is not collected. pixi tasks: baseline, red-manifest, affected, figures, the gates, and the Held tasks.
- **defects:**
  - Task `schema` targets the missing `tests/adaptive/test_schema.py`.
  - Task `mutations-adaptive` targets the missing `scripts/run-adaptive-mutations.py`.
  - A comment cites the missing `tests/benchmarks/test_package_contract.py`.
  - `[tool.cibuildwheel] test-command` runs all of `tests/`, including the reds.
  - `wheel.py-api = "cp312"` against `requires-python >= 3.9` and the 3.10 CI.
- **unused capability / ideas:** add a `gate` task, i.e. the product test plus an attributed report, scoring any registered generator.
- **claims:** the comment "tests/benchmarks/test_package_contract.py holds the part ..." is contradicted: the file does not exist.

### recipe/sha256.py | lines 10 | read 1-10
- **role:** build
- **used by:** conda recipe (manual)
- **state:** dead (hard-codes `compas_cgal-0.7.3.tar.gz`; current version 0.9.1)
- **generation:** n/a
- **defects:** stale version.
- **unused capability / ideas:** none
- **claims:** none

### requirements-dev.txt | lines 16 | read 1-16
- **role:** build
- **used by:** docs.yml, environment.yml
- **state:** superseded-by pixi
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none
- **claims:** none

### requirements-docs.txt | lines 15 | read 1-15
- **role:** build
- **used by:** docs.yml
- **state:** superseded-by pixi docs env
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none
- **claims:** none

### requirements-viz.txt | lines 1 | read 1-1
- **role:** build
- **used by:** none found
- **state:** dead
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### requirements.txt | lines 1 | read 1-1
- **role:** build
- **used by:** docs.yml
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### scripts/engagement_baseline.py | lines 292 | read structural
- **role:** tool
- **used by:** manual (SP2 reference)
- **state:** parked (docstring: the full run is deferred)
- **generation:** n/a. Audits the unregulated generator.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "~1 hour full run" is unverified.

### tasks.py | lines 26 | read 1-26
- **role:** build
- **used by:** invoke (compas-actions in build.yml)
- **state:** live (upstream)
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none
- **claims:** none

### tests/__init__.py | lines 0 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/task13f_fixture.py | lines 393 | read structural
- **role:** test
- **used by:** `tests/adaptive/test_generator.py`
- **state:** parked (Task 13F plan superseded, `red_manifest` backlog D1)
- **generation:** n/a; adaptive generator fixture.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_acceptance.py | lines 1412 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.bootstrap`/`entry`/`generator`: entry-circle launch, branch switch, finite search over neck caps.
- **defects:** none found
- **unused capability / ideas:** the exact entry-circle launch is an entry capability the gated generators lack.
- **claims:** none

### tests/adaptive/test_candidates.py | lines 863 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.candidates` (MATHSM circle proposals on the exact MAT).
- **defects:** none found
- **unused capability / ideas:** exact MAT circle proposals
- **claims:** none

### tests/adaptive/test_canonical.py | lines 542 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.canonical`.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_circle_oracle.py | lines 666 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exact circle engagement oracle.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_containment.py | lines 318 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.containment` (motion inside the reachable domain).
- **defects:** none found
- **unused capability / ideas:** exact gouge containment, where the gate only samples gouges.
- **claims:** none

### tests/adaptive/test_coverage.py | lines 564 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.coverage`.
- **defects:** none found
- **unused capability / ideas:** exact coverage
- **claims:** none

### tests/adaptive/test_entry.py | lines 283 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.entry`: a pre-cleared entry bore (`:117`) and a first full circle that must fit the entry void (`:231`).
- **defects:** none found
- **unused capability / ideas:** pre-drilled/pre-cleared entry contract (plan D1 alternative)
- **claims:** none

### tests/adaptive/test_event_substrate.py | lines 1111 | read structural
- **role:** test
- **used by:** pytest, pixi `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_exact_algebraic_identity.py | lines 200 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_exact_binary64_contract.py | lines 76 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_exact_depletion.py | lines 638 | read structural
- **role:** test
- **used by:** pytest, pixi `audit-native`
- **state:** live
- **generation:** n/a. Exact arc depletion.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_generation_state.py | lines 511 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_generator.py | lines 634 | read structural
- **role:** test
- **used by:** pytest
- **state:** live, with 2 manifest reds (`test_task13f_full_continuation`, `test_real_active_family_stops_at_unresolved_exact_event`; `docs/red_manifest.json`)
- **generation:** n/a. Global adaptive orchestration.
- **defects:** the generator does not complete a full continuation (manifest reds, verified). It has no `ToolpathResult` output, so it cannot be gate-scored.
- **unused capability / ideas:** the certified adaptive generator, if bridged to `ToolpathResult`.
- **claims:** none

### tests/adaptive/test_identity.py | lines 636 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_medial_axis.py | lines 523 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.medial_axis`, the exact segment-site MAT.
- **defects:** none found
- **unused capability / ideas:** an exact MAT guide, an alternative to the straight skeleton (draft plan Task 1)
- **claims:** none

### tests/adaptive/test_middle_curve_path.py | lines 106 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.middle_curve_path` (max-radius, boundary-phase circles per MAT edge).
- **defects:** none found
- **unused capability / ideas:** middle-curve placement (Held Figure 5), which would address the `loop_radius_cv` defect
- **claims:** contradicts the "does not exist yet" claim in `quality.py`/`test_quality.py`.

### tests/adaptive/test_motion.py | lines 157 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_motion_certificate.py | lines 558 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Continuous-motion cap certificate.
- **defects:** none found
- **unused capability / ideas:** certified (not sampled) cap, where the gate samples
- **claims:** none

### tests/adaptive/test_motion_counterexamples.py | lines 85 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_motion_oracle_cache.py | lines 375 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_motion_refutation.py | lines 513 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exact witness-based cap refutation.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_neck.py | lines 352 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.neck`: an exact neck inventory and legal passage steps.
- **defects:** none found
- **unused capability / ideas:** exact neck detection, which plan review-focus item 2 needs and the tree walk does not use
- **claims:** none

### tests/adaptive/test_operation.py | lines 1003 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_policy.py | lines 371 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_reachable_domain.py | lines 827 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.reachable_domain`.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_remaining_material.py | lines 70 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Residual-only exact replay.
- **defects:** none
- **unused capability / ideas:** remaining-material extraction (rest-operation input)
- **claims:** none

### tests/adaptive/test_replay.py | lines 1265 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_retrace_segment_trial.py | lines 361 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Advancing and retrace segments: a stay-down retrace through cleared stock.
- **defects:** none
- **unused capability / ideas:** exact retrace through cleared material, i.e. a certified stay-down link (plan D3)
- **claims:** none

### tests/adaptive/test_route_retrace_generator.py | lines 518 | read structural
- **role:** test
- **used by:** pytest
- **state:** live, with 2 manifest reds
- **generation:** n/a
- **defects:** superseded-plan reds (verified, manifest)
- **unused capability / ideas:** route retrace = linking across branches without a retract
- **claims:** none

### tests/adaptive/test_route_retrace_operation.py | lines 707 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** see above
- **claims:** none

### tests/adaptive/test_route_retrace_state.py | lines 744 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_route_retrace_transaction.py | lines 1018 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_segment_event_proof_contracts.py | lines 136 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_segment_event_substrate.py | lines 553 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_segment_oracle.py | lines 433 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Segment engagement oracle; a certified straight-move cap exists.
- **defects:** none
- **unused capability / ideas:** an exact straight-motion engagement certificate, which supports judging lines by the full cap (finding 2)
- **claims:** none

### tests/adaptive/test_station_generic_double_witness.py | lines 287 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. A generic-double witness with a wall-clock ceiling.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_stock_monotonicity.py | lines 93 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_transaction.py | lines 1524 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none found. At 1524 lines it is over the size bar.
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_traversal.py | lines 931 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `adaptive.traversal`: exact identity-driven global MAT traversal.
- **defects:** none found
- **unused capability / ideas:** global MAT traversal order, which competes with the plan's DFS tree walk
- **claims:** none

### tests/adaptive/test_trimmed_line_branches.py | lines 523 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_uncut_coverage.py | lines 85 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Continuous coverage of explicitly uncut design and reachable regions.
- **defects:** none
- **unused capability / ideas:** exact uncut-region accounting (D2(a) `DeclinedRegion` check)
- **claims:** none

### tests/adaptive/test_units.py | lines 141 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_zero_guide_replay.py | lines 400 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Zero-guide segment advances.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/test_zero_guide_transaction.py | lines 473 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/typecheck/auditor_contract.py | lines 79 | read structural
- **role:** test
- **used by:** pixi `types-audit`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/adaptive/typecheck/consumer_contract.py | lines 521 | read structural
- **role:** test
- **used by:** pixi `types-adaptive`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/__init__.py | lines 0 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/depth_bench.cpp | lines 49 | read structural
- **role:** test
- **used by:** none found (not in CMakeLists)
- **state:** dead
- **generation:** n/a
- **defects:** an orphan native bench
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/filter_bench.cpp | lines 60 | read structural
- **role:** test
- **used by:** none found (not in CMakeLists)
- **state:** dead
- **generation:** n/a
- **defects:** orphan
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/held_figure5_path_live_oracle.py | lines 320 | read structural
- **role:** test
- **used by:** pixi `held-figure5-characterize-live` (not named `test_*`)
- **state:** live (opt-in)
- **generation:** n/a. Held Figure 5 2D path.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/held_reference_figures_live_oracle.py | lines 28 | read structural
- **role:** test
- **used by:** opt-in (needs publisher PDF)
- **state:** live (opt-in)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/held_reference_live_reconstruction_oracle.py | lines 68 | read structural
- **role:** test
- **used by:** pixi `held-reference-reconstruct-live`
- **state:** live (opt-in)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_analytic.py | lines 43 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a; pocket families only.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_cli.py | lines 125 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** does not cover the `--corpus` form CI uses (the CI arg bug goes uncaught)
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_complexity.py | lines 73 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_congruence.py | lines 85 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_coverage_grid.py | lines 85 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Coverage grid on Held cases.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_degeneracy.py | lines 81 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_depletion.py | lines 98 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Runner (unregulated) generator.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_exceedance.py | lines 94 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Unregulated generator.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_external.py | lines 137 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_figure6.py | lines 264 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `engagement_controlled_toolpath` and `trochoidal_mat_toolpath_circular`.
- **defects:** `:101` "controlled path honours its cap away from the entry on this pocket" pins a current measurement on one pocket. It is not harmful to a better generator.
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_figures.py | lines 154 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_gate_attribution.py | lines 202 | read 1-202
- **role:** test
- **used by:** pytest
- **state:** live (green)
- **generation:** n/a. Exercises `GATE_GENERATORS` (engagement_controlled, radius_regulated).
- **defects:**
  - `:164-169` requires byte-identical streams for both generators at 120°, which pins radius regulation as inert at the default cap.
  - `:97-115,180-190` the C3 witness rejects `Arc`, so any gated generator emitting arcs (G2 connectors, Held-style arcs) fails here (finding 4).
  - `:161` exactly two generators are assumed.
- **unused capability / ideas:** none
- **claims:** the "default-cap streams identical" premise is verified as a pin (and by the identical pin rows at `test_quality.py:1259-1264`).

### tests/benchmarks/test_held_boundary_order_spacing.py | lines 65 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held boundary-order spacing.
- **defects:** none
- **unused capability / ideas:** Held contour-aware spacing
- **claims:** none

### tests/benchmarks/test_held_concentric_engagement.py | lines 79 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held standard placement.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_contour_bound_spacing.py | lines 87 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** contour-bound spacing
- **claims:** none

### tests/benchmarks/test_held_coverage_thinning.py | lines 119 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Circle deletion that preserves material.
- **defects:** none
- **unused capability / ideas:** remove air loops while preserving coverage, a direct tool for finding 8
- **claims:** none

### tests/benchmarks/test_held_curved_circle_benchmark.py | lines 96 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_curved_circle_performance.py | lines 22 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Timing regression.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_curved_circle_performance_plot.py | lines 30 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (modified in the working tree)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_exact_boundary_transition.py | lines 96 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exact boundary contacts.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_exact_motion_coverage.py | lines 69 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_boundary_path.py | lines 261 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held Figure 5 boundary path; plan Task 0b's adapter source.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_engagement_refinement.py | lines 130 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_exact_boundary.py | lines 37 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_path.py | lines 499 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held Figure 5 path.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_publisher.py | lines 59 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_raw_guide.py | lines 139 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure5_reference_guided.py | lines 93 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Publisher-guided start; the brief's "needs publisher-derived start data".
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_figure6_comparison.py | lines 170 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. `figure6.controlled_path` against Held.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_live_oracle_boundary.py | lines 55 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_motion_coverage.py | lines 43 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_native_boundary_draft.py | lines 61 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Native corner events.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_native_circle_geometry.py | lines 62 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_native_curve_import.py | lines 40 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Native line/arc boundaries.
- **defects:** none
- **unused capability / ideas:** a true curved boundary import, which `PocketSpec` lacks
- **claims:** none

### tests/benchmarks/test_held_native_curved_circle_budget.py | lines 32 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Wall-clock ceiling.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** "Monstera 2.6–3.1 s" is unverified.

### tests/benchmarks/test_held_path_characterize.py | lines 487 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_path_evidence.py | lines 484 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_path_report.py | lines 278 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Uses `runner.generate_toolpath` (unregulated).
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_path_snapshot.py | lines 636 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. `HeldOperationSnapshot`, the gate's source snapshot.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_post_qualification.py | lines 286 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_reference_cases.py | lines 561 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held reference pockets (Figure 5/8).
- **defects:** none
- **unused capability / ideas:** Held pockets as gate pockets (curved, real-world)
- **claims:** none

### tests/benchmarks/test_held_reference_contacts.py | lines 54 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_reference_figures.py | lines 369 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_reference_geometry.py | lines 910 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_reference_toolpaths.py | lines 30 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Figure 8 inputs reach the circle/connector generator.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_standard_placement.py | lines 156 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held standard placement.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_stock_replay.py | lines 51 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_toolpath_benchmark.py | lines 79 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_held_toolpath_performance.py | lines 22 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Timing regression.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_instrument.py | lines 62 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_mathsm.py | lines 115 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Unregulated constant-spacing generator.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_measurement.py | lines 84 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_necks.py | lines 94 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_pathmetrics.py | lines 158 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. `trochoidal_mat_toolpath_circular`.
- **defects:** `:131-141` pins a raw maximum of 360° from a plunge entry. It is fine for the unregulated generator and does not involve gated ones.
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_plotting.py | lines 408 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_precision.py | lines 99 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_quality.py | lines 1340 | read 1-1340
- **role:** test
- **used by:** pytest; `docs/red_manifest.json` (12 reds); `test_quality_transition.py` imports its helpers
- **state:** live
- **generation:** n/a. Exercises `GATE_GENERATORS` on the 3 pockets × 2 caps.
- **defects:**
  - `:1312-1316` asserts the pinned failure vector inside the product test, before the verdict. An improved generator fails on a pin mismatch, without the diagnostic report.
  - `:1280-1294` the pin builder always emits 6 violation lines (only slotting is conditional). A generator that fixes degenerate loops alone is still a mismatch.
  - `:154-176` path length, air, retracts and plunges are explicitly ungated.
  - `:160-162` "material_entries: one per chain is unavoidable without a helical entry" bakes per-chain plunging into the rationale (contradicted by the tree-walk design: one entry per connected region).
  - `:122-125` the slot rationale says a straight transfer "travels through material the loops already cleared". That holds only for the current generators' bridges, not for straight cutting moves in general (finding 2).
- **unused capability / ideas:**
  - The `_report` printout (`:1149-1242`) is the full diagnostic. Make it available for any generator without the pin.
  - At 120°, `engagement_controlled` and `radius_regulated` measure identically (`:1259-1264`). At 40°, `radius_regulated` has more redundant operations (311 vs 26 on rect_12x8, `:1265-1266`), meaning more air loops from the radius ladder.
- **claims:**
  - "THRESHOLDS DERIVED, NEVER READ OFF A MEASUREMENT" (`:93-95`): the loop step of 2r is derived (`:134-145`). Slot = cap/2 is self-described engineering judgement (`survey.py:86-88`), so it is not derived. Contradicted for that one rule.
  - Kernel rim-arc table: verified as a pinned test (`:262-285`).

### tests/benchmarks/test_quality_invariants.py | lines 659 | read structural
- **role:** test
- **used by:** pytest; 1 manifest red (translation)
- **state:** live
- **generation:** n/a. Hypothesis properties over synthetic stations.
- **defects:** none found. Its docstring explicitly records that earlier metrics "encoded generator HABITS as physics". The slot and redundant-link rules are two more of those habits.
- **unused capability / ideas:** extend the property families to "retracting between chains may not lower redundant count unless …" style run semantics.
- **claims:** none

### tests/benchmarks/test_quality_observations.py | lines 1417 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:**
  - `:439-451` pins that a ramp line raises `UnreplayableOperationError`, which blocks plan D1 until changed (finding 3).
  - `:588-610` the closed vocabulary pins the thresholds.
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_quality_transition.py | lines 31 | read 1-31
- **role:** test
- **used by:** pytest
- **state:** live (permanent green)
- **generation:** n/a. `GATE_GENERATORS`.
- **defects:** `:16-31` a green test asserts that the 12 gate cases produce exactly today's failure vectors. It turns red on ANY change to the gated generators, improvements included (finding 4).
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_qualityfigures.py | lines 219 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_report.py | lines 114 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_runner.py | lines 78 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_spec.py | lines 37 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_survey.py | lines 221 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** `:48-66` pins rejection of tilted or off-plane circles. That is consistent with the cut-plane model, but it blocks helices (finding 3).
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_toolpath_coverage.py | lines 66 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Controlled and constant-spacing paths.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_topology.py | lines 67 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/test_units.py | lines 130 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/__init__.py | lines 0 | read structural
- **role:** test
- **used by:** mypy
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/gate_contract.py | lines 10 | read structural
- **role:** test
- **used by:** mypy (`[tool.mypy] files`)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/held_path_characterization_contract.py | lines 186 | read structural
- **role:** test
- **used by:** pixi `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/held_reference_cases_contract.py | lines 54 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/held_reference_contract.py | lines 98 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/held_reference_figures_contract.py | lines 50 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/benchmarks/typecheck/quality_observations_contract.py | lines 67 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/build/test_native_source_integrity.py | lines 254 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Native dependency digests.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/build/test_station_lane_decoder_ratchet.py | lines 242 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Number-type stage-2 ratchet (parked lane).
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/conftest.py | lines 19 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held publisher-PDF fixture plumbing.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_arc_authority.py | lines 46 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_classification.py | lines 204 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_identity.py | lines 117 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_input.py | lines 535 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_native_classification.py | lines 309 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_native_replay.py | lines 300 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_native_request_identity.py | lines 378 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/engagement_audit/test_records.py | lines 241 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/audit_certification_fixtures_2.h | lines 430 | read structural
- **role:** test
- **used by:** `audit_native_gate` sources (included by 5 files)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/audit_replay_fixtures_2.h | lines 207 | read structural
- **role:** test
- **used by:** `audit_native_gate` (included by 4 files)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/exact_canonical_gate.cpp | lines 184 | read structural
- **role:** test
- **used by:** pixi `exact-gates` (CI benchmarks.yml)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/exact_one_root_gate.cpp | lines 167 | read structural
- **role:** test
- **used by:** `exact-gates` (CI)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/exact_segment_attestation_gate.cpp | lines 562 | read structural
- **role:** test
- **used by:** `exact-gates` (CI)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/exact_station_attestation_gate.cpp | lines 336 | read structural
- **role:** test
- **used by:** `exact-gates` (CI)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/exact_vocabulary_gate.cpp | lines 95 | read structural
- **role:** test
- **used by:** `exact-gates` (CI)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/sign_mixed_radical_gate.cpp | lines 441 | read structural
- **role:** test
- **used by:** `exact-gates` (CI)
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task3_algorithm_gate.cpp | lines 1215 | read structural
- **role:** test
- **used by:** pixi `task3-algorithm-gate` (local only)
- **state:** live (modified in the working tree: the `shared_ptr<ReachSet>` change)
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task3_certificate_gate.cpp | lines 759 | read structural
- **role:** test
- **used by:** `task3_algorithm_gate` target (CMakeLists `:233`)
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task3_certificate_gate.h | lines 3 | read structural
- **role:** test
- **used by:** task3 gate sources
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_catalog_delaunay_gate.cpp | lines 274 | read structural
- **role:** test
- **used by:** pixi `task9-mat-compile-gate` (local only)
- **state:** live
- **generation:** n/a. Segment-site MAT.
- **defects:** not in CI
- **unused capability / ideas:** the exact segment-site MAT, the guide alternative to the straight skeleton
- **claims:** none

### tests/native/task9_catalog_graph_gate.cpp | lines 1363 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a. Segment-site MAT graph.
- **defects:** not in CI
- **unused capability / ideas:** MAT graph (a tree or cycle structure for the walk)
- **claims:** none

### tests/native/task9_catalog_neck_gate.cpp | lines 258 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** exact necks
- **claims:** none

### tests/native/task9_catalog_sampling_gate.cpp | lines 185 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_catalog_voronoi_gate.cpp | lines 239 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_domain_identity_gate.cpp | lines 99 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_graph_csr_gate.cpp | lines 613 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_indexed_endpoint_binding_gate.cpp | lines 1286 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_mat_certificate_gate.cpp | lines 163 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_mat_compile_gate.cpp | lines 1833 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate` (entry point)
- **state:** live
- **generation:** n/a
- **defects:** not in CI. At 1833 lines it is over the size bar.
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_mat_numeric_table_gate.cpp | lines 168 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_mat_proposal_table_gate.cpp | lines 169 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a. MAT circle proposal tables.
- **defects:** not in CI
- **unused capability / ideas:** native placement proposals
- **claims:** none

### tests/native/task9_mat_sampling_gate.cpp | lines 354 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_neck_classification_gate.cpp | lines 266 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_neck_clearance_gate.cpp | lines 398 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_neck_evidence_gate.cpp | lines 681 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_neck_gate.cpp | lines 516 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_point_graph_gate.cpp | lines 1067 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_rational_sources_gate.cpp | lines 242 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_segment_limiter_gate.cpp | lines 625 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_segment_segment_gate.cpp | lines 4183 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a. Segment/segment bisector contracts.
- **defects:** not in CI. At 4183 lines it is far over the size bar.
- **unused capability / ideas:** none
- **claims:** none

### tests/native/task9_site_catalog_gate.cpp | lines 377 | read structural
- **role:** test
- **used by:** `task9-mat-compile-gate`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_arc_motion_2.cpp | lines 304 | read structural
- **role:** test
- **used by:** pixi `audit-native` (`audit_native_gate`)
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_2.cpp | lines 14 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_contract_2.cpp | lines 430 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_falsifiers_2.cpp | lines 315 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_identity_2.cpp | lines 458 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_refinement_2.cpp | lines 345 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_certification_verdicts_2.cpp | lines 459 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_identity_2.cpp | lines 483 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_replay_2.cpp | lines 12 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_replay_finalization_2.cpp | lines 255 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_replay_identity_lineage_2.cpp | lines 802 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_replay_transaction_state_2.cpp | lines 340 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_audit_replay_trial_depletion_2.cpp | lines 674 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/native/test_exact_arc_depletion_2.cpp | lines 456 | read structural
- **role:** test
- **used by:** `audit-native`
- **state:** live
- **generation:** n/a
- **defects:** not in CI
- **unused capability / ideas:** none
- **claims:** none

### tests/test_boundary_circle_contact.py | lines 140 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Native exact proposal contacts (Held placement).
- **defects:** none
- **unused capability / ideas:** contour-aware placement primitive
- **claims:** none

### tests/test_boundary_circle_engagement.py | lines 115 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held predecessor decisions on native circle proposals.
- **defects:** none
- **unused capability / ideas:** predecessor-engagement placement
- **claims:** none

### tests/test_boundary_normal_circle.py | lines 123 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Boundary-normal circle proposals at the true first medial contact.
- **defects:** none
- **unused capability / ideas:** boundary-normal placement, which reaches near-wall and corner material that skeleton-centred loops miss
- **claims:** none

### tests/test_circle_redundancy.py | lines 69 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Local circle deletion that preserves material.
- **defects:** none
- **unused capability / ideas:** exact post-pass removal of redundant loops (finding 8)
- **claims:** none

### tests/test_circle_sweep_stock.py | lines 47 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_engagement_audit.py | lines 436 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `trochoidal_mat_toolpath_circular`.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_engagement_oracle.py | lines 324 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Monte-Carlo raster oracle for station TEA.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_engagement_radial_toolpath.py | lines 652 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `radius_regulated_toolpath` and `engagement_controlled_toolpath`.
- **defects:**
  - `:574-588` requires plunges > chains and plunges == retracts. This pins multi-pass retract/plunge behaviour, so a stay-down multi-pass fails.
  - `:408-431` pins the virgin-stock entry as the maximal full-slot circle.
  - `:590-614` "clears everything" exempts all material within r of a wall, which is weaker than the gate's tool-reachable uncut.
- **unused capability / ideas:** none
- **claims:** none

### tests/test_engagement_rho_toolpath.py | lines 263 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `rho_regulated_toolpath`.
- **defects:** `:204-213` requires PLUNGE and RETRACT in the vocabulary (mild).
- **unused capability / ideas:** tangent bridges at any radius difference (`:145-161`), no bores (`:185`), overlapping annuli (`:195`), and `DeclinedRegion` maximal runs (`:222-263`). These are proven pieces that plan Stage A should reuse.
- **claims:** none

### tests/test_engagement_spiral_entry_toolpath.py | lines 144 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `spiral_entry_toolpath` on `gate_pocket("rect_12x8")`.
- **defects:** none. `:127-139` asserts the first loop is in [240°, 360°). That correctly pins a planar-entry floor.
- **unused capability / ideas:** the spiral ramp entry. It is not gated, and it is the best planar entry available.
- **claims:** "no non-degenerate first loop can engage below 240 degrees after one plunge": stated as a theorem in the test, unverified here.

### tests/test_engagement_toolpath.py | lines 364 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `engagement_controlled_toolpath` against the unregulated generator.
- **defects:**
  - `:148-170` requires the set of over-cap operations to EQUAL the set of chain-entry loops. A generator whose entry is under the cap fails (finding 4).
  - `:293-331` "clears everything" exempts the wall band.
- **unused capability / ideas:** none
- **claims:** none

### tests/test_exact_region_traits_ownership.py | lines 217 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (modified in the working tree)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_false_certificate.py | lines 225 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_growth_bound.py | lines 173 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_held_contour_bound.py | lines 69 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Whole-orbit leading-rim upper bound from cleared-contour clearance.
- **defects:** none
- **unused capability / ideas:** a cheap sound engagement bound for placement, which could replace sampled probes in the generators
- **claims:** none

### tests/test_held_disk_contour.py | lines 84 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Held §3.1 clockwise contact correction.
- **defects:** none
- **unused capability / ideas:** Held contact correction
- **claims:** none

### tests/test_native_boundary_curve.py | lines 82 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exact C0 line/arc boundary construction.
- **defects:** none
- **unused capability / ideas:** true curved pockets
- **claims:** none

### tests/test_native_boundary_medial.py | lines 280 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (modified in the working tree)
- **generation:** n/a. Analytic curved-boundary medial contacts.
- **defects:** none
- **unused capability / ideas:** an arc-site medial axis (a guide on curved pockets without tessellation blow-up)
- **claims:** none

### tests/test_replay_classification.py | lines 76 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Shared cut-plane replay decision.
- **defects:** `:49-57` pins `CutPlaneRampError` for any ramp with XY travel, and `:74-76` pins `OffPlaneReplayCurveError` for any engaged curve off the plane. Both block a helical entry (finding 3, plan D1). It does not pin LINK-at-depth explicitly (only `LINK` above the plane = rapid, `:65`).
- **unused capability / ideas:** none
- **claims:** none

### tests/test_stock.py | lines 969 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_stock_capsule_quad.py | lines 352 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Quad capsule depletion.
- **defects:** none
- **unused capability / ideas:** fast exact link depletion for generator-side stock (plan Task 3)
- **claims:** none

### tests/test_stock_clone_lifetime.py | lines 154 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_stock_local_depletion.py | lines 213 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** local depletion (generator stock cost)
- **claims:** none

### tests/test_stock_traits_ownership.py | lines 201 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_straight_skeleton_2_interior_straight_skeleton.py | lines 34 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (upstream)
- **generation:** n/a. Straight skeleton (the gated guide).
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_straight_skeleton_2_interior_straight_skeleton_offset_polygon.py | lines 35 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (upstream)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** offset polygons (contour-parallel finishing pass for wall residue)
- **claims:** none

### tests/test_straight_skeleton_2_interior_straight_skeleton_weighted_offset_polygon.py | lines 30 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (upstream)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_straight_skeleton_2_interior_straight_skeleton_with_holes.py | lines 40 | read structural
- **role:** test
- **used by:** pytest
- **state:** live (upstream)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/test_toolpath.py | lines 879 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a. Exercises `trochoidal_mat_toolpath` and `trochoidal_mat_toolpath_circular` (C++ unregulated).
- **defects:**
  - `:741-744` `max_passes` truncation only WARNS. The plan's `skeleton_guide` raises, which is consistent with the plan, but the existing generator silently continues apart from the warning.
  - Tangent tests use `|dot|` (known issue in auto-memory; not re-read).
- **unused capability / ideas:** Hypothesis no-gouge and tangent properties on random polygons (`:573-593,844-879`). This property harness should be pointed at the gated and new generators.
- **claims:** none

### tests/tools/__init__.py | lines 0 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/held_reference_live_pdf_oracle.py | lines 40 | read structural
- **role:** test
- **used by:** pixi `held-reference-live`
- **state:** live (opt-in)
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_held_path_characterization.py | lines 205 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_held_reference_case_generator.py | lines 70 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_held_reference_extractor.py | lines 328 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_held_reference_qualification.py | lines 145 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_lane_equivalence.py | lines 174 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_plan_headers.py | lines 101 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/test_red_manifest.py | lines 570 | read structural
- **role:** test
- **used by:** pytest
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/typecheck/__init__.py | lines 0 | read structural
- **role:** test
- **used by:** mypy
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/typecheck/held_path_characterization_contract.py | lines 60 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/typecheck/held_reference_extractor_contract.py | lines 46 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tests/tools/typecheck/held_reference_qualification_contract.py | lines 26 | read structural
- **role:** test
- **used by:** `types-benchmarks`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/__init__.py | lines 0 | read structural
- **role:** tool
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/held_figure5_publisher_check.py | lines 23 | read structural
- **role:** tool
- **used by:** pixi `held-figure5-publisher-check`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/held_path_characterization.py | lines 79 | read structural
- **role:** tool
- **used by:** pixi `held-figure5-characterize`
- **state:** live
- **generation:** n/a. Characterises the "default generator" on Held Figure 5.
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/held_reference_case_generator.py | lines 340 | read structural
- **role:** tool
- **used by:** pixi `held-reference-generate`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/held_reference_extractor.py | lines 810 | read structural
- **role:** tool
- **used by:** pixi `held-reference-extract`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/held_reference_qualification.py | lines 239 | read structural
- **role:** tool
- **used by:** pixi `held-reference-qualify`
- **state:** live
- **generation:** n/a. Qualifies the existing generator on the four Held pockets.
- **defects:** none
- **unused capability / ideas:** the Held pockets as a qualification corpus. It should feed the product gate.
- **claims:** none

### tools/lane_equivalence.py | lines 403 | read structural
- **role:** tool
- **used by:** `tests/tools/test_lane_equivalence.py`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/plan_headers.py | lines 73 | read structural
- **role:** tool
- **used by:** pixi `plan-headers`
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### tools/red_manifest.py | lines 259 | read structural
- **role:** tool
- **used by:** pixi `red-manifest`
- **state:** live
- **generation:** n/a. Compares suite reds with `docs/red_manifest.json`, in both directions.
- **defects:** none found. Note that the manifest regexes name only the two gated generators, so adding a generator to the gate requires a manifest edit (plan Task 4).
- **unused capability / ideas:** none
- **claims:** none
