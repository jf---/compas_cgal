# Product Pipeline: One Guide, Two Regimes, One Gate

> **status: active** — the one active plan, 2026-10-13. It replaces the same-day
> tree-walk draft that once lived at this path; the census showed that draft
> re-derived a deleted generator and steered by a broken engagement model.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One toolpath pipeline whose paths a machinist would run. The gate must read it as world class on the gate pockets (`rect_12x8`, `rect_20x12`, `L_shape`) and on Held's published pockets: capped engagement everywhere, one entry, the tool kept down, G1, no air loops, full coverage, and a short path.

**Architecture:** A typed pipeline in `src/compas_cgal/pipeline/`, built from parts that already exist:
- **Guide.** One medial-axis guide.
- **Placement, in two regimes by local wall-to-axis distance `d`:**
  - narrow (`d ≤ 3r`): stay-down trochoid chains centred on the axis, from the deleted ordered generator plus rho's station rules;
  - wide (`d > 3r`): Held's boundary lap, from the off-live analytic generator.
- **Entry:** a helical bore.
- **Corners:** wall run-ins, then declared rest material.
- **Gate:** engagement measured in process, with exact sweeps.

**Tech Stack:** C++20, CGAL 6.0.1, nanobind, Python 3.12, compas, pixi, pytest-xdist, ruff, mypy.

**Spec:** The gate is the spec (`docs/machining_quality.md`), as amended by rulings R1–R4 below. Every claim in this plan is backed by the census in `docs/superpowers/census/2026-10-13/`: 753 files and 57 unmerged commits, coverage-certified, 230 registered findings, and measurements re-run by the orchestrator. The census disposition of every finding is in `2026-10-13-product-pipeline-traceability.md`.

## The spine (user, 2026-10-13)

Sessions optimised for rigour, and nobody owned moving the gate's numbers. This plan enforces the opposite:

1. **One generator pipeline.** The four existing generators are frozen as baseline columns and are parts sources, not competitors.
2. **One guide:** the pocket's medial axis.
3. **One gate.**
4. **Every commit moves a gate number on a gate pocket.** Setup code ships inside the task whose commit moves the number. A CI ratchet (Task 0.3) rejects a commit that touches `src/compas_cgal/pipeline/` or the gate and improves no pinned vector, or that worsens any.
5. **Integrate before inventing.** Every task names the existing parts it integrates. New code is allowed only where the census found no part.
6. **The orchestrator integrates; agents implement one task each.**

## Rulings (user, 2026-10-13)

- **R1, gate physics:**
  - (a) Engagement is measured in process, against the stock minus the motion's own swept prefix. This applies in the gate *and* in every generator decision.
  - (b) The sweep model is exact and idempotent.
  - (c) A stay-down link through cleared stock is legal; retracts and plunges are gated.
  - (d) Straight motions are judged by the same cap as loops.
  - (e) Tests that pin today's bad behaviour are retired.
- **R2, pipeline:** two regimes on one medial-axis guide. Narrow regions get axis-centred stay-down trochoid chains; wide regions get Held's boundary lap. Held's pockets join the gate.
- **R3, entry:** a helical bore as its own operation. Its replay removes its disk at the cut plane, it is scored by ramp angle against the machine model, and the step from the helix to the first loop stays gated.
- **R4, corners:** first, run-ins along each wall offset to the corner apex, under the cap, in the same operation. What remains is declared as named rest material for a second operation, reported separately from material that was missed.

## Where we are (re-measured by the orchestrator on 2026-10-13)

| fact | evidence |
| --- | --- |
| All 12 product-gate cases fail; the gated generators have not changed behaviour since 2026-08-22 | gate run 2026-10-13, 4m52s; pin matches to the digit |
| The gate scores a straight slot into solid stock at **360°** in its interior; physically it is 180° | `census/2026-10-13/verify/verify_direction.py` |
| A replayed loop is not idempotent: a plunge at its own entry changes the stock 8/12 times, and re-cutting it from another start 12/12 | `census/2026-10-13/verify/verify_sweep.py` |
| The live Held Figure 5 reproduction has 2,013 loops (1,059 with ρ ≤ r, smallest 0.00046 mm) and 25.9 m of loops; Held's published path has 265–278 turns, none ≤ r, smallest 1.167 mm, 8.8 m | `census/2026-10-13/verify/verify_g6.py` |
| The off-live analytic Held lap on Figure 5 gives **270 loops, 0 ≤ r, smallest 1.1658 mm, 9.16 m**: Held's own path | `census/2026-10-13/verify/verify_analytic.py` (run in the corpus worktree) |
| The same lap fails at sharp convex corners (rect_12x8, L_shape, 80° and 120°), needs a CCW boundary, and works with rounded corners (61 loops, 29 ≤ r: its ρ = (d−r)/2 needs `d > 3r`) | `census/2026-10-13/verify/verify_analytic_gate.py`, `verify_analytic_orient.py` (corpus worktree) |
| A stay-down chain-ordering generator existed (685 lines). It took rect_12x8 from 5 entries to 1, then was deleted as "unreachable" | `git show 9fa6eea8^:src/compas_cgal/engagement_ordered_toolpath.py` |
| The exact segment-site medial axis and the certified `adaptive/` generator run only on one hard-coded L fixture; adaptive never completes a pocket | `segment_site_mat_bundle.cpp:39`, `segment_site_catalog_graph.cpp:172-192`, `generator.py:1362-1364`, `replay.py:1276` |
| Plans written: 164 tasks across 22 plans; zero generation tasks ever landed on a gated generator | census d2 |

## Pipeline stages (the target architecture)

```
Pocket ─▶ Guide ─▶ Regions ─▶ Placement ─▶ Order+Link ─▶ Entry ─▶ Corners ─▶ ToolpathResult(+rest, +forced) ─▶ Gate
          medial    narrow/     narrow: axis    stay-down,     helical   wall run-ins,
          axis      wide/       trochoids       in-process     bore      declared rest
                    unreachable wide: Held lap  certified
```

One file per stage in `src/compas_cgal/pipeline/`, plus `generate.py`. Typed stage outputs; each stage fails loudly on a contract breach; `mypy --strict` gates the package.

## Global Constraints

- pixi only: `pixi run pytest <paths> -n auto -q`; full suite `pixi run baseline` → `pixi run red-manifest`.
- Failing test first. No skip, skipif or xfail.
- **No bisection** anywhere in a search. Engagement is not monotone (rho and radial docstrings), so use a finite ladder in exhaustive order (auditor P3 Task 2, carried).
- Exact decisions run on exact models: no exact predicate over an under-covering or surrogate-radius sweep (census g4).
- Forced decisions are structured outputs, not warnings: forced advances, forced loops, plunges into material, declined rest material (auditor P3 Task 1, carried).
- New package code meets the engineering bar: units and frames in types, `Cls.build` factories, one named exception per failure mode, one responsibility per file.
- Commits: author and committer Jelle Feringa, short subject naming the gate delta, no attribution lines.
- **Rolling wave:** only the active task is expanded to step level (failing-test code, implementation steps, exact commands), from the sources it names, and reviewed before execution. Every task's acceptance is fixed here, up front, as a gate-vector delta.

## Review Focus (inputs no gate pocket exercises yet)

1. Pockets with islands: the medial axis has cycles, and the centre domain may split. The reachable kernels throw `PocketNotMachinableError` on more than one component (`reachable_arrangement_2.cpp:997-1006`). Task 1.7.
2. Necks narrower than a trochoid, and channels `W ≤ 4r`, which came out as 68 plunges on the old generators. Tasks 1.7 and 2.3.
3. Curved walls. Only tessellations exist in the families, and the straight skeleton explodes on them. Task 2.1 (Held pockets are curved).
4. Acute and obtuse convex corners, not only 90°. Task 1.3.
5. Non-integer CAD coordinates, where the exact predicates see generic input and CORE stall risk applies (memory: identically-zero decisions). Task 1.7.

---

## Stage 0: the gate measures physics (R1). Three tasks, each moving every engagement number.

### Task 0.0: Clean baseline (precondition; outside the ratchet's scope)

The live worktree carries an uncommitted 11-file diff from 2026-09-08/09: import paths, formatting, a regenerated stub, and the `shared_ptr<ReachSet>` gate.
1. Run `pixi run lint` and `pixi run red-manifest`. Commit the diff only if the suite reds equal the manifest.
2. In one docs commit, land this plan, its traceability matrix and `docs/superpowers/census/2026-10-13/`. In the same commit, stamp status headers on the number-type stage 0–3 plans (landed 0–2, parked 3) and park `2026-09-05-held-paper-figure-reproduction.md` behind this plan, so that `pixi run plan-headers` passes.

### Task 0.1: In-process engagement, in the gate and in the generators' decisions

- **Parts:**
  - `src/swept_prefix_engagement_2.cpp` (`engagement_line_prefixes` / `engagement_arc_prefixes`, corpus commit `7be75e0d`), ported onto live's rewritten `engagement_2.cpp` and `audit_exact_station_2.cpp`;
  - live `HeldDiskContour2.engagement_bound` (leading rim) as the cross-check.
- **Change:**
  - `benchmarks/survey.py` samples each motion against the stock minus its own swept prefix. Full circles are arc prefixes from 0 to 2π.
  - The generators' probe rings (rho, the restored ordered generator) use the same query.
  - The two arc-prefix decisions currently made in doubles (`swept_angle <= pi`, census g7 #4) become exact.
- **Acceptance:**
  - a straight slot into solid stock reads 180° ± one sample spacing;
  - the 12 gate vectors are re-pinned;
  - every engagement criterion (cap exceedances, max step, slotting) moves on at least one pocket, and the uncut, gouge, tangent and continuity criteria stay unchanged.
- **Moves:** cap exceedances, max engagement step, slotting.

### Task 0.2: Exact, idempotent sweeps in replay and generator

- **Parts:**
  - `Stock2::subtract_circle_sweep(c, ρ_authored, r)` (`stock_2.cpp:408-421`) for `Circle`, replacing the `hypot(point_at(0))` surrogate radius (`engagement.py:195-207`, `stock_2.cpp:673-697`);
  - one capsule model shared by generator and replay: the 6-curve quad (`stock_2.cpp:516-540` replaced), or an exact capsule where the direction is Pythagorean;
  - entry points authored exactly incident on the circle (`exact_circle_chart_2.cpp:38-65`);
  - `Stock2.can_remove_circle` and the radial `_loop_reaches_material` as the cheap exact redundancy tests, in place of clone + subtract + compare.
- **Cost:** the replay and the generator deplete *locally*, by zone surgery with the same exact sweep (`stock_local_2.cpp:290-324`). This includes the capsule local depletion that docs/local_depletion.md:147-149 left unattempted, and bounds the O(n²) exact depletion that CHANGELOG.md:13 deferred.
- **Acceptance:**
  - re-cutting a loop from any start and plunging at its own entry change nothing (12/12);
  - the arrangement after four 16 mm runs stays under 200 vertices, down from 8,014;
  - the gate's redundant-operation counts are re-pinned;
  - generation plus gate replay on each gate pocket takes no more than 2× the 2026-10-13 time.
- **Moves:** redundant operations (now truthful), uncut (sliver residue gone), runtime.

### Task 0.3: Gate rules and the ratchet

- **Rules:**
  - (a) Redundant means a cut or loop that removes nothing. A `LINK` at cut depth that removes nothing and stays under the in-process cap is a legal stay-down transit, reported as air-at-depth length.
  - (b) Plunges into material are counted by the gate (an exact plunge-disk ∩ stock test), and required to be ≤ the number of centre-domain components. Retracts are reported.
  - (c) Slotting means a straight motion over the cap: the same rule as loops.
  - (d) Degenerate loops are reported, not gated; the cap decides.
  - (e) A plunge's XY must equal the next motion's start.
  - (f) Engagement step is measured at junction values, not per-motion peaks.
- **Retire** (R1e):
  - `test_quality_transition.py`;
  - the pin-equality assert inside the product test (`test_quality.py:1312`);
  - the identical-streams and Arc-rejection tests (`test_gate_attribution.py:97-115,164-169`);
  - the entries-must-exceed-cap test (`test_engagement_toolpath.py:148-170`);
  - the plunge-per-pass tests (`test_engagement_radial_toolpath.py:408-431,574-588`).
- **Ratchet:** the pinned vectors move to `benchmarks/results/gate_vectors.json`. `tools/gate_ratchet.py` fails a commit that improves no vector or worsens any. CI runs it, after fixing `benchmarks.yml`'s broken `--corpus smoke` step (`cli.py:284-292`) and stopping cibuildwheel from running the red gate.
- **Results, not prose:** gate vectors live in one results layout (`benchmarks/results/gate_vectors.json`). Every number in a stage's docs page is read from it, never typed. A stage exits only when the ratchet passes and every acceptance row of that stage is green, checked by one exit-gate test per stage.
- **Acceptance:** the old generators' vectors are re-pinned under the new rules; the ratchet passes on this commit and fails on a synthetic regression.
- **Moves:** redundant and slotting on the old generators, as truthful re-pins.

## Stage 1: narrow regime, the pipeline's first form, on rect_12x8, rect_20x12 and L_shape

### Task 1.1: Restore stay-down chain ordering as `pipeline/` (narrow regime)

- **Parts:**
  - `engagement_ordered_toolpath.py` from `9fa6eea8^`: both orientations, links at depth admitted by the exact predicate, cap-checked bridges, `isolated_chains`;
  - rho's station rules: ρ > r floor, radius-step bound, external-tangent entries, furthest-first ladder (`engagement_rho_toolpath.py:235-537`);
  - Task 0.1's in-process query for loops, bridges and links.
- **Files:** `pipeline/guide.py` (reads the C++ guide for now), `pipeline/narrow.py`, `pipeline/linking.py`, `pipeline/generate.py`, `pipeline/result.py` (forced-decision ledger). `benchmarks/gate.py` registers `pipeline`.
- **Acceptance (120°, all three pockets):**
  - plunges into material: 1 per pocket;
  - redundant 0, degenerate 0;
  - tangent breaks ≤ rho's 0/0/2;
  - no criterion worse than rho_regulated's 2026-10-13 vector under the Stage 0 gate.
- **Moves:** entries 5→1, redundant, degenerate, tangent breaks.

### Task 1.2: Helical entry at the widest point (R3)

- **Parts:**
  - `OperationType.HELIX`, replayed as `subtract_circle_sweep(c, ρ_h, r)` (an exact disk when ρ_h ≤ r, `stock_2.cpp:155-163`);
  - loop-in-entry-disk containment (`containment_2.cpp:322-356`);
  - the pre-cleared entry contract types (`adaptive/entry.py:45-263`);
  - `MachineModel` for the ramp angle;
  - the root chosen as the maximum clearance over the guide, not over skeleton vertices. The gate L ties three vertices at 3.0, while the true widest point is (3.515, 3.515) (census g1 P-2).
- **Gate:**
  - `replay_classification.py` admits HELIX;
  - `quality_observations.py:773-833` partitions it;
  - the ramp angle ≤ the machine limit is a criterion;
  - the step from the helix to the first loop is gated.
- **Acceptance:** cap exceedances 0 and max engagement step ≤ cap on all three pockets at 120°.
- **Moves:** cap exceedances, max step.

### Task 1.3: Corners, first run-ins and then declared rest (R4)

- **Parts:**
  - the exact cutter-centre boundary `build_center_boundary_cycle`, `ReachableBoundaryCycle2.ccw_transition`, and the `is_stationary` apex (`reachable_arrangement_2.cpp:673-890`);
  - `ReachableDomain2.reachable_material` / `unreachable_residual` (`reachable_domain_2.cpp:119-190`);
  - `ReachableMaterialPredicate2`;
  - rho's `DeclinedRegion`, now *consumed*;
  - ETH Task 31's contract: clear with a non-trochoidal motion, or declare.
- **Gate:** uncut is split into *declined* (named rest material, a separate criterion) and *missed* (required 0).
- **Cost:** the acute and obtuse pockets have oblique edges. On those, `center_domain` was measured about 3,500× slower than on axis-aligned edges, mechanism unknown (docs/benchmarks.md:13-16, oblique_edge_cost.md:13-18). Establish the mechanism before these pockets enter the gate.
- **Acceptance (120° and 40°, including acute and obtuse corner pockets added to the gate):**
  - missed = 0 on all pockets;
  - run-ins under the cap;
  - generation on an oblique pocket no more than 3× the time on its axis-aligned twin.
- **Moves:** uncut.

### Task 1.4: Residue that isn't at the corners

- **Parts:**
  - `_coverage_2.remaining_material` and `ReachableMaterialPredicate2` to localise the residue (3.28 from any vertex on the old generators, `machining_quality.md:150-158`);
  - the annulus-overlap invariant (radius step ≤ 2r);
  - coverage brought into the accept rule (INTEGRATECODEX.md:157-161).
- **Acceptance:** missed uncut = 0 on all gate pockets at 120° and 40°.
- **Moves:** uncut.

### Task 1.5: Regulate engagement, don't merely limit it

- **Parts:**
  - re-measure, under the in-process model, the 141.8° saturation and the 171.8° loops seen at both 40° and 120° (`benchmarks.md:120-133`, `jf:review.md:470-476`);
  - a ladder over advance *and* radius;
  - the design's `tea_target ± band`, with ≥ 95% of length in band (adaptive-clearing design :34-35), and a lower bound against rubbing;
  - root cause of cap-40 saturation (ETH Task 37 Step 1).
- **Acceptance:**
  - every 40° gate case loses at least one criterion;
  - immersion at the design load ≥ 0.9 (reported today: 0.07);
  - no 120° regression.
- **Moves:** all 40° cases, immersion.

### Task 1.6: G1 and G2 everywhere, without air

- **Parts:**
  - bridge slope over the *chosen* advance rather than neighbour stations (rho `:298-305`; B4's two L breaks);
  - arc links and leads (backlog C5) with exact arc depletion (`Stock2::subtract_exact_arc`, `stock_2.h:215`);
  - an exact gouge certificate for every bridge and link via `_containment_2` (`containment_2.cpp:225-356`), not `certify_bridges`, whose decisions run on rounded doubles and can exit with a bridge uncertified;
  - in the C++ guide, retire three epsilon decisions:
    - the 1e-18 squared-distance threshold that chooses between a repositioning arc and a coincident point (`toolpath_tangent_continuity.md:83-90`);
    - the `abs(dot)` parallelism test, which cannot see a reversal (`:152-165`);
    - the hand-rolled perpendicular in `external_tangent` (ETH Task 17), replaced by kernel primitives. Confirm the old 1e-9 clamp is gone (redesign plan :262-265);
  - the partial-arc trochoid against full-circle air (recut 0.90);
  - `can_remove_circle` post-pass deletion of air loops.
- **Gate:** the tangent test stops exempting bridges (`toolpath_tangent_continuity.md:33-51`), and the recut fraction is reported against a target.
- **Acceptance:** tangent breaks 0 on all pockets; recut fraction ≤ 0.6 (reported today: 0.90).
- **Moves:** tangent breaks, recut.

### Task 1.7: Guide hygiene, islands and necks

- **Parts:**
  - jf's `toolpath.cpp` contract fixes:
    - `mat_scale ∈ (0,1]` (`d90e18b8`);
    - `clearance_z`, `radial_clearance` and min/max validation (`92cedda2`);
    - the discontinuity-preserving tessellator (`b2fa087e`/`34bddcc0`);
    - `b0094cb6` only together with `446daca5`, never alone: it breaks the guide call;
  - hole validation into the straight skeleton (ETH I9);
  - `max_passes` truncation fails loud and is recorded in the forced-decision ledger; today it only warns (CHANGELOG.md:19, adaptive spec :583);
  - one plunge-disk region per centre-domain component.
- **Neck semantics, fixed here and used everywhere:**
  - an *entry neck* is a split of the centre domain (width < 2r), which costs one entry per component;
  - a *regime neck* is `d ≤ 3r` between wide regions, which switches to axis trochoids;
  - the exact neck evidence and classification (`segment_site_neck_*`) stay parked with the exact segment-site medial axis.
- **Gate:** adds `families.topology.island_grid` and `families.necks.dumbbell`.
- **Acceptance:**
  - both new pockets generate with plunges = components, missed uncut 0, and no new failing criterion;
  - a neck ablation, with neck handling on and off on the dumbbell, shows the neck rule moves at least one criterion.
- **Moves:** the island and neck pockets.

## Stage 2: wide regime, Held's boundary lap

### Task 2.1: Port the analytic Held lap into `pipeline/wide.py`, and add Held's pockets to the gate

- **Parts (corpus):** `held_figure5_analytic.py`, from commits `757b948c`, `5b5ed867`, `4685653b`, `567da10d` and `274fd1a4` (land `8ca2c9f6` together with `274fd1a4`).
- **Merge work:**
  - `ReferenceBoundary` gains a required `logical_source_indices`;
  - the boundary is normalised to CCW;
  - r and the cap come from `PocketSpec`;
  - the motion stream is adapted to `ToolpathResult`: loop start at the wall contact q; tool-centre ρ; zero-length fan connectors and the wrap-around connector dropped;
  - circles are built natively, exact on the offset (`_coverage_2.boundary_circle_at_contact`), so contacts lie on the offset by construction. The gouge check then sees no tangency noise; today 1,010 "rejections" are 6.6e-15 mm rounding;
  - exact circular-arc input gets a declared construction policy: stored arc endpoints differ in radius by about 8.5e-16 (segment_site_mat.md:1706-1715);
  - the Held reference-boundary reconstruction items are closed before Held's pockets enter the gate, since they are gate inputs: the reconstruction node budget, the re-review, and the biarc split.
- **Gate:**
  - adds Figure 5 and Figure 8 (upper, crossed skis, Monstera), with the publisher turn evidence corrected to 278 (`a662a9c3`);
  - the minimum coverage grid via `held_consumer_adapters.py:35`;
  - Held's published loop count and length shown beside each.
- **Acceptance (Figure 5):** loops ≤ 278, length ≤ 9.4 m, none with ρ ≤ r, missed uncut 0, under the Stage 0 gate.
- **Moves:** the Held pockets, from unscored to scored.

### Task 2.2: Sharp convex corners in the lap

- **Parts:**
  - Held's convex-corner regularisation, radius slightly > r (held_motion_coverage.md:43-58);
  - the `is_stationary` apex and Task 1.3's run-ins;
  - the oversized-circle risk in `_locus_circle_from_state` (census g7 #1).
- **Acceptance:** the lap runs on rect_12x8, rect_20x12 and L_shape (it raises today) with no criterion worse than Stage 1's vector.
- **Moves:** gate pockets now runnable by the wide regime.

### Task 2.3: Region split and joint ordering

- **Parts:**
  - local `d` from the guide: Held lap where `d > 3r`, axis trochoids where `d ≤ 3r`;
  - stay-down transitions between regimes (Task 0.3a legal links), *measured* by the in-process query;
  - the swept-prefix retrace for depth-first transit (zero-guide plan :875-890), used as a measured motion here. Using it as a *certificate* waits on the proof note still owed from auditor P2 T4 (Task 4.2).
- **Acceptance:** Figure 8 upper and Monstera (narrow arms plus wide bays) improve on the better single regime in loops, length and criteria.
- **Moves:** the Figure 8 pockets.

### Task 2.4: Contour-aware placement and boundary-order spacing

- **Parts:**
  - Held §3.1 contour-aware placement (`held_figure5_contour_placement.py`, corpus `8ca2c9f6`) and `HeldDiskContour2` exposed contact;
  - boundary-order spacing (`held_boundary_order_spacing.py`: 2.4× fewer circles on Figure 8 upper);
  - concave-corner shared-contact fans (`held_figure5_engagement_refinement.py:175-193`);
  - MATHSM on the same locus (`held_figure6_mathsm.py`) as the Figure 6 comparison.
- **Acceptance:**
  - Figure 5 length below the standard lap's, with missed uncut 0. The contour-bound draft halved the length but lost coverage, and preserving coverage cost +82% (`held_figure5_coverage_preserving.json`);
  - Figure 6 reproduced like-for-like.
- **Moves:** Held pocket length and loops.

### Task 2.5: One guide

- **Parts:**
  - the CGAL segment-Delaunay medial-axis graph (`src/segment_site_guide.cpp`, corpus `07c87ce0`, `f90d1ce4`, `334c028b`; any polygon with holes, analytic parabolas);
  - the exact boundary-normal medial query (`native_boundary_medial_2.cpp:78-274`, `BoundaryNormalCircle2`) for per-point disks.
- The straight skeleton is retired as a guide once the L hub measures better: the true widest point and G1 parabolic arcs, against the skeleton's 90° kink.
- **Acceptance:** L_shape and Figure 8 vectors at least as good as with the skeleton.
- **Moves:** L_shape.

### Task 2.6: Generality on the external corpus

- **Parts:** `benchmarks/families/external.py` (a third-party profile loader), plus the island, neck, curved and n-gon families.
- **Change:** run the pipeline over the whole external corpus. The resolution rate (the share of pockets fully generated, with no raise) becomes a gate number.
- **Acceptance:** resolution rate is reported, then ratcheted. Every failure is a named exception attributed to a stage.
- **Moves:** resolution rate.

## Stage 3: on the machine

- **Task 3.1:** a fail-closed G-code post with a parser round-trip (reference-pockets plan :865-896, never started). It includes a `ManufacturingPolicy` (linear tolerance, angular surrogate, minimum feature size; exactness.md:432-449, proposed and never built) as the one place where user tolerances enter.
- **Task 3.2:** feeds and cycle time as the objective:
  - curvature-limited feeds (`MachineModel.curvature_limited_feed_mm_per_s`, `quality.machine_outcome`, `quality.py:732-763`);
  - chip-thinning feed compensation from the measured `h_ex/f_z`.
- **Task 3.3:** a machine-simulation gate on a named model: the gate's exact 2.5D stock replay, plus a mechanistic cutting-force model (specific cutting force per tooth). Its force KPIs are peak force, force-step, and the share of path at design load. It also measures controller starvation: the block rate at the programmed feed against the controller's look-ahead, using the block-length statistics the gate already reports.
- **Task 3.4:** ordering by cycle time. Revisit the abandoned dynamic-programming transition solver (redesign plan :515-520, :745-760) once cycle time is the objective.
- **Moves:** cycle time and force KPIs on every gate pocket, which become gated criteria here, by ruling.

## Stage 4: certificates, after the paths are good

- **Task 4.1:** certified advance and radius from exact roots (`exact_algebraic_1.cpp:1240-1277`). Refute with probes first, then certify the accepted rung (census g5).
- **Task 4.2:** the segment oracle's forward-semicircle mask, so straight cuts certify below 180° (`segment_oracle.cpp:505-513,598-602`).
- **Task 4.3:** remove the falsified growth lemma (`engagement_2.cpp:407-499`; the jf swept-annulus repair `2e2422b5`, `03866fcc`, `73d53729`). Until then, nothing consumes `cap_certified` (`benchmarks/pathmetrics.py:198`, `runner.py:86`, `held_consumer_adapters.py:20`).
- **Task 4.4:** total-contact cap, not only per run (`exact_cap_certificate_2.h:84-140`).
- **Task 4.5:** robustness and cost of the oracle the certificates stand on: the 13 findings listed under Task 4.5 in the traceability matrix. These are the CORE identically-zero stalls, the all-pairs resultants, rebuild-to-verify, the unused univariate `sign_mixed_radical`, the stubbed circle refuter, and the oracle's cost of 100–150× the probes. The swept-prefix proof note from auditor P2 T4 is part of Task 4.2.

## Frozen

- `engagement_controlled_toolpath`, `radius_regulated_toolpath`, `rho_regulated_toolpath` and `spiral_entry_toolpath` stay as baseline gate columns and parts sources. No behaviour changes, no removal without a ruling.
- The live Held reproduction (`held_figure5_path.py` lane, 2,013 loops) is superseded by Task 2.1 and stays as a comparison until then.

## Parking lot

Each item is parked with its reason; reopening one needs a ruling.
- **auditor-convergence P0–P4:** P3 Tasks 1–2 are carried as Global Constraints; the rest is governance.
- **number-type stages 3–6, the global type-system refactor:** no toolpath change.
- **exact full-motion coverage replay:** set aside 2026-09-08; Task 0.2 makes the cheap replay exact instead.
- **`adaptive/` as a generator:** it never completes a pocket. Its parts are harvested in Tasks 1.2, 1.3 and 2.3.
- **generalising the exact segment-site medial axis (8 steps, census g3 (e)):** only if Task 2.5's segment-Delaunay guide proves insufficient.
- **NRVO lifetime hazard on by-value `Gps` returns** (census g4 #8): unverified that it fires.
- **dead code in `engagement_2.cpp:29-36,86-345`; dead pixi tasks `schema`, `mutations-adaptive`, `measured-run`:** no toolpath change.
- **translation variance at exact rim tangency (C1):** a deliberate red, until the gate is rebuilt in Stage 0.
- **`jf/exact-geodesics`:** upstream, unrelated.
- **Licensing:** the CGAL packages the pipeline uses (straight skeleton, segment Delaunay graph, arrangements) carry GPL terms. A licence gate is needed before any distribution of a binary or of G-code tooling. It changes no toolpath; it becomes a ruling at Stage 3.
- **Docs and governance hygiene:** stale or false status headers, dangling references, deleted-file pointers, and claims that contradict measurements. Also the governance backlog items and packaging metadata (the requires-python floor). Examples: "dropping degenerate loops changes coverage 0.000%" against a measured 0.284%→0.351%, and "V ≥ √3 r proved", which is absent from its page. Each is fixed in the docs update of the stage whose page it sits on, never as standalone work. The ids are listed under this line in the traceability matrix.
