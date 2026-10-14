# Census slice d2_plans: the 22 SDD plans in docs/superpowers/plans/

Reader: Opus 5.5 (opus-5-5), 2026-10-13. All 22 assigned plans read line by line (32,239 lines in all).
The draft under review (`2026-10-13-product-gate-tree-walk.md`) was read only for its header, "Where we are",
"Design", "Global Constraints", "Review Focus", Task 0, Task 0b, "Stage B", "Stage C" and "Parking lot"
(lines 1-160 and 1374-1411), so it could be compared against these plans. It is not one of this slice's records.
"Live branch" means `codex/held-task5-continuation` at `db549852` plus its uncommitted 11-file diff.
No code was run. All evidence below comes from `grep`, `ls`, `git log` and `sed`.

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **Not one task in the 22 plans targeted the gated generators or the 12 product-gate criteria.** I counted
   164 tasks (table below). The 33 generation tasks all built either the exact-certified adaptive pipeline
   or the Held reproduction lane. The one plan that would have changed the gated code,
   `engagement_controlled_toolpath` (auditor P3 Tasks 1-3: ladder search, forced-decision ledger, typed
   geometry), never ran. That is consistent with the brief's "unchanged since 2026-08-22". It also means the
   only generator work that ever touched the gated paths has no surviving plan. *verified*: `ls`
   of `src/compas_cgal/` shows no `engagement_override.py`, `engagement_search.py` or
   `engagement_geometry.py`.
2. **The exact segment-site MAT and the certified `adaptive/` generator only accept one hard-coded L
   pocket.** That pocket is (0,0)-(6,0)-(6,2)-(2,2)-(2,6)-(0,6). The guard is
   `src/segment_site_catalog_graph.cpp:159-193` (`require_l_shape_input`). The production proposal graph
   calls only `canonical_l_shape_mat_proposal_graph` (`src/segment_site_mat_bundle.cpp:39`,
   `src/segment_site_catalog_sampling.cpp:277-283`). A rectangle MAT graph exists
   (`canonical_rectangle_mat_graph`, `src/segment_site_catalog_graph.cpp:905`), but nothing calls it. Even on
   its own fixture the adaptive generator has never finished a pocket:
   `tests/adaptive/test_generator.py::test_task13f_full_continuation` is a declared red
   (`docs/red_manifest.json`). **This contradicts the brief's suggestion** that the exact segment-site MAT
   or the certified generator could replace the straight skeleton in Stage A. Neither can run on
   `rect_12x8`, `rect_20x12` or the gate's `L_shape` (`benchmarks/gate.py:122-125`, whose arm lengths are
   tool-diameter multiples, not the fixture) without first building a general MAT.
3. **The gated generator still searches advances by bisection, under a monotonicity assumption that the
   repository already knows is false.** `src/compas_cgal/engagement_toolpath.py:601-655` states "Engagement
   is taken to be monotone ... a stated MODELLING ASSUMPTION". `engagement_rho_toolpath.py:94-101` documents
   that it is not monotone, which is why rho uses a ladder. The fix (auditor P3 Task 2) was planned and never
   executed. The draft's choice to reuse rho's ladder, not the gated walk, is therefore correct, and the
   draft should say so.
4. **Depth-first route order is not a cutter trajectory.** This was measured in the adaptive work and is
   directly relevant to the draft's tree walk. On the Task 13F fixture, route 2 failed outright after a
   depth-first activation: 56 of 56 candidate links gouged. Reversing the last accepted segment through its
   own swept corridor restored 6 valid link/circle pairs (zero-guide plan :1295-1306; causal plan :439-449).
   The fix built for it was an exact cut-depth retrace. Its swept-prefix theorem certifies that a translation
   starting clear engages at most π, in 0.04 s, against more than 10 minutes for the generic oracle
   (zero-guide plan :875-890). The draft repositions above the stock for every non-incident branch switch.
   That is sound, but each switch costs a retract and a plunge. The retrace is the existing, measured
   stay-down alternative and belongs in D3.
5. **The Held reproduction says straight-skeleton guides break Held's coverage assumptions, and it found
   the generic substitute.** Held needs true medial-axis incidence, with each circle's diameter equal to
   the q-m distance. The draft "uses straight-skeleton guides and later moves circles without preserving
   that incidence". Its coverage-thinning regression "fails at every scale" (paper plan :369-386). The
   working generic replacement is a native boundary-normal q-m circle construction for any simple polygon,
   with segment interiors and reflex point-site sectors. It ran 65 Figure 5 queries in 1.87 s (:393-411,
   `src/boundary_normal_circle_2.cpp`, `src/native_boundary_medial_2.cpp`). The draft's Stage A keeps the
   skeleton and inherits rho's 1-3.5% uncut. The draft should at least record this capability as the
   alternative to the skeleton that actually exists.
6. **Draft Task 0b as written will hit two measured failures.** First, the default coverage grid is too
   coarse for Figure 5. The second live characterization run died after 780.83 s with
   `CoarseCoverageGridError`; the minimum grid is 668x458 (characterization plan :1553-1567). The adapter
   that applies it exists (`benchmarks/held_consumer_adapters.py:35`,
   `max(COVERAGE_GRID_SAMPLES, minimum_coverage_grid(spec))`), but the draft does not name it. Second, the
   Held refined path has known residuals that will show up in the gate's gouging and uncut criteria: 1,010
   near-tangent containment rejections (excess about 6.7e-15 mm, paper plan :196-199) and an unreachable
   2.04 µm boundary remnant (:596-600). The draft's decision gate ("fewer failed criteria than rho") needs
   to know about both.
7. **Machine-readiness work was planned twice and never started.** The reference-pocket plan's "next
   phase" (:865-896: per-case evidence, generator repairs, fail-closed G-code post, parser round-trip,
   machine simulation) never began. No G-code emitter exists: `grep` finds only a comment in
   `src/compas_cgal/toolpath.py:67`. The characterization plan's "Phase 2 generator-policy design"
   (:1798-1800) was never written.
8. **The exact rest-material region already exists, which bears on the draft's D2(a).**
   `ReachableDomain2.unreachable_residual()` gives the exact D \ M_r (`src/coverage_bindings_2.cpp:390`),
   and `ReachableMaterialPredicate2` (`src/reachable_material_predicate_2.h:15`) is a fast exact membership
   test: 305,944 queries in 34 s on Figure 5 (characterization plan :1693-1698). D2(a) as written ("measure
   uncut only over what a trochoid with ρ > r reaches") needs a region that is narrower than M_r. M_r is
   the floor any such ruling must state, and the draft should use it rather than a new sampled region.

### Generation capabilities in this slice that the gated generators do NOT use

| Capability | Where | State | Use for the draft |
| --- | --- | --- | --- |
| Native boundary-normal q-m circle (true-MA incidence, generic simple polygon, reflex sectors) | `src/boundary_normal_circle_2.cpp`, `src/native_boundary_medial_2.cpp`; paper plan :393-411 | live (Held lane) | an alternative to the skeleton guide that actually exists |
| Contour-aware contact (`HeldDiskContour2`, filled outer-disk union) | `src/held_disk_contour_2.cpp`; paper plan :610-620 | live, contact query only | contour-aware placement (Held Fig. 5b) |
| Concave-corner approach/rotation/departure split | `benchmarks/held_coverage_path.py`, `held_figure5_engagement_refinement.py` (`--corner-approaches`); paper plan :632-641 | live, Held lane | corner motion for B2 |
| Boundary-order spacing (predecessor carried across families) | `benchmarks/held_boundary_order_spacing.py`; paper plan :653-662 | live | spacing |
| Swept-prefix exact theorem (clear-start translation engages at most π) | `src/continuous_tea_2/segment_oracle.cpp`; zero-guide plan :875-890 | live, adaptive only | certifying stay-down links and retraces (D3) |
| Exact cut-depth retrace of the last advancing segment | `src/compas_cgal/adaptive/retrace_*`; retrace plan T1-T5 | live, L only | the stay-down alternative to retract and plunge |
| Exact zero-guide (clearance == r) MAT classification and advancing-segment operation | `src/segment_site_zero_guide.cpp`, `adaptive/advancing_segment_trial.py` | live, L only | the model for constant-clearance arms (gate ρ ≤ r regions) |
| Exact reachable domain C_r, M_r, D \ M_r; `ReachableMaterialPredicate2` | `src/reachable_domain_2.cpp`, `src/reachable_material_predicate_2.h` | live | D2 rest-material ruling; exact uncut floor |
| Exact partial-arc depletion `Stock2::subtract_exact_arc` | `src/stock_2.h:215`; P1 T3 | live | arc link vocabulary (backlog C5) |
| Exact gouge containment for segments and circles | `src/containment_2.cpp`, `adaptive/containment.py` | live | plunge and link clearance |
| MATHSM constant-spacing baseline | `benchmarks/mathsm.py`; corpus plan T11 | live | Stage C comparison column |
| Transactional candidate evaluation (fork, certify, deplete, commit) | `adaptive/transaction.py` | live, L only | the pattern for "productive-only emission" |

### Where this slice shows the draft plan wrong or incomplete

- **Design, "Tree walk":** it assumes the skeleton is a sound guide. The paper plan records that skeleton
  guides break Held's coverage assumptions (:380-386), and a generic true-MA circle construction exists
  (finding 5). The draft does not mention either. It also does not cite the adaptive finding that
  depth-first route order needs transit at non-incident switches (finding 4). It handles that case by
  repositioning, but should name the cost and the retrace alternative.
- **Brief/draft suggestion "exact segment-site MAT" or "certified adaptive/ generator":** contradicted.
  Both are L-fixture only (finding 2).
- **Task 0b:** it does not name `minimum_coverage_grid` or `held_consumer_adapters`. A 600 s budget is at
  risk: survey plus coverage on Figure 5 were never measured inside 600 s, and the only full run took
  780 s including a guarded audit. It also ignores the known Held-path residuals (finding 6).
- **Stage B D1 (helix):** consistent with corpus plan :3241-3245 and backlog C4. The depletion model
  rejects ramped 3D moves (`benchmarks/depletion.py` `_replay_kind` raises `UnreplayableOperationError`).
  *verified*.
- **Stage B / Stage C miss backlog C5:** "Inter-chain links cannot be tangent-continuous while loops are
  full circles" (corpus plan :3246-3248). B4 only attributes the L_shape tangent breaks; it does not adopt
  the arc vocabulary or the exact arc depletion that already exists (P1 T3).
- **Stage B D2(a):** should be defined on the existing exact D \ M_r and `ReachableMaterialPredicate2`
  (finding 8), not on a new sampled measure.
- **Parking lot "auditor-convergence P0-P4":** P3 is the only generation work in the auditor programme.
  T1 is the forced-decision ledger; T2 is the ladder search over bisection. The draft must not lose the
  constraint behind T2: no bisection anywhere in the tree walk's advance search. It inherits rho's ladder,
  so this is satisfied, but it should be pinned as a Global Constraint.
- **Task 0 Step 4 (status headers):** headers in this slice are false, or point at deleted files, beyond
  the five the draft fixes (register F-d2_plans-19 to 22). That is governance only, so it is listed for
  completeness.
- **Missing altogether:** a G-code or postprocessor stage. The reference-pocket plan made fail-closed
  G-code post a gate (:865-896). The brief's bar ("a machinist would run") needs it. The draft has no stage
  for it.

### Effort by category (task counts)

Categories: **G** = generation (code that produces or changes toolpath motion: guide, placement, ordering,
linking, entry, corner, engagement control, and the generator-internal stock/transaction machinery);
**V** = validation or certification (oracles, replay, audit, coverage proofs, quality judges);
**M** = measurement (benchmarks, corpora, figures, performance, reference inputs); **Gov** = governance
(build, identity, manifests, ledgers, headers, docs-only); **K** = number-type kernel plumbing (byte-identical
by design, no toolpath effect).

| Plan | G | V | M | Gov | K | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 07-24 phase1 master | 6 | 8 | 1 | 5 | 0 | 20 |
| 07-25 task3 single-arrangement | 0 | 5 | 1 | 0 | 0 | 6 |
| 07-28 atomic transaction | 3 | 0 | 0 | 1 | 0 | 4 |
| 07-28 causal MAT traversal | 5 | 1 | 1 | 0 | 0 | 7 |
| 07-29 one-root full-circle source | 1 | 3 | 0 | 1 | 0 | 5 |
| 07-29 zero-guide link | 6 | 1 | 1 | 0 | 0 | 8 |
| 08-07 inter-route retrace | 5 | 1 | 0 | 1 | 0 | 7 |
| 08-19 benchmark corpus | 0 | 1 | 13 | 0 | 0 | 14 |
| 08-23 auditor P0 | 0 | 0 | 0 | 5 | 0 | 5 |
| 08-23 auditor P1 | 0 | 6 | 0 | 0 | 0 | 6 |
| 08-23 auditor P2 | 1 | 4 | 0 | 0 | 0 | 5 |
| 08-23 auditor P3 | 2 | 2 | 1 | 0 | 0 | 5 |
| 08-23 auditor P4 | 0 | 0 | 1 | 4 | 0 | 5 |
| 08-28 coherence wave 1 | 0 | 1 | 2 | 8 | 0 | 11 |
| 08-29 wave1 task 6 | 0 | 0 | 0 | 3 | 0 | 3 |
| 08-31 Held reference pockets | 0 | 0 | 5 | 0 | 0 | 5 |
| 09-02 Held Fig. 5 characterization | 0 | 10 | 7 | 0 | 0 | 17 |
| 09-05 Held paper reproduction | 4 | 0 | 4 | 0 | 0 | 8 |
| 09-08 number-type stages 0-3 | 0 | 0 | 0 | 0 | 23 | 23 |
| **total** | **33 (20%)** | **43 (26%)** | **37 (23%)** | **28 (17%)** | **23 (14%)** | **164** |

Of the 33 G tasks, 29 sit in the adaptive or exact-MAT line, which is L-fixture only and never produced a
complete path. Four sit in the Held lane, which depends on publisher-derived starts. Two (P3 T2 and T3)
targeted a gated generator and never ran. **Zero G tasks landed on a gated generator.** Read by landing
rather than by planning: of the 131 non-G tasks, every completed one hardened a judge, an oracle or the
governance around them.

### Generation goals set and never reached

| Goal | Plan | Evidence it was not reached |
| --- | --- | --- |
| Non-vacuous polygonal-pocket path with certified engagement, gouge-freedom and exact coverage | phase1 :14-16, completion :2713-2731 | `GenerationResult`, `certificate.py`, `schema.py`, `schemas/` and `scripts/run-adaptive-mutations.py` are all absent; Task 13F full continuation is a declared red |
| True segment-site MAT for arbitrary pockets | phase1 T9 :1216ff | L-only guard, `segment_site_catalog_graph.cpp:159-193` |
| Terminal seal: covered generation plus fresh replay | causal T6 :355-489; zero-guide T8 Step 2; P2 T1-T2 | no `GenerationResult`; no `replay_certificate.py` |
| Held Fig. 5 runtime targets (10 ms at 20°, 0.9 ms at 80°) | phase1 :2618-2628 | never measured; Held lane takes 20.6 s on Fig. 5 (paper plan :541-543) |
| Ladder search replacing bisection; forced-decision ledger in public results | P3 T1-T2 | `engagement_toolpath.py:601-655` still bisects; no ledger module |
| Generator repairs plus fail-closed G-code on the four Held pockets | reference-pockets :865-896 | "That phase is not started"; no G-code code |
| Phase 2 generator-policy design from the Fig. 5 failure signature | characterization :1798-1800 | Task 11 Steps 3-7 retired, no live report |
| Continuous, publisher-free Fig. 5 standard path from MA traversal | paper T5 :180-234 | all five boxes open |
| Fig. 6 three-curve protocol (standard, contour-aware, MATHSM) | paper T6 :236-256 | open; gallery lists the panels as blocked (`docs/held_toolpath_gallery.md:25,70`) |
| Every Fig. 5-8 panel generated; coverage-accepted drafts | paper T8 :332-770 | Fig. 5(b)-(d) and Fig. 7 blocked; coverage acceptance suspended (:523-533) |
| Zero-guide event so upper and Monstera curved-circle runs complete | paper :567-574 | open; stops at upper 76 and Monstera 211 |
| Contour-aware Fig. 5/8 drafts | paper :610-620 | contact query only |

### Plan items that should carry into the draft, and items the draft misses

Carry in, as rulings or tasks:
1. **P3 T2 as a constraint:** no bisection; a finite ladder in exhaustive order. The draft inherits it via
   rho; make it explicit.
2. **P3 T1 forced-decision ledger:** the draft's `TreeWalkResult.virgin_plunges` is half of this. Forced
   advances and forced loops at the cap should also be structured outputs, not warnings.
3. **Zero-guide and retrace findings for D3:** depth-first order needs transit, and a stay-down reverse
   retrace is certified by the swept-prefix theorem.
4. **Held true-MA q-m construction** as the alternative guide, at least as a measured Stage C comparison
   against the skeleton walk on the gate pockets.
5. **Corner approaches (Monstera)** as a B2 candidate next to the wall run-in.
6. **Exact D \ M_r and the reachable predicate** as the basis for D2(a).
7. **Reference-pocket next phase:** G-code post and machine-simulation gates, as a Stage D.
8. **Corpus open item C5:** arc link vocabulary, plus exact arc depletion (P1 T3) to score it.
9. **Task 0b plumbing:** `held_consumer_adapters` and the minimum coverage grid.

Missed by the draft: items 3, 4, 5, 6, 7, 8 and 9 (item 1 is satisfied implicitly, item 2 only partly).

---

### docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md | lines 2731 | read 1-2731

- **role:** doc (master plan: exact-certified adaptive trochoidal MATHSM on a segment-site MAT).
- **used by:** `tools/plan_headers.py` (header lint); cross-referenced by the atomic, causal and one-root plans.
- **state:** superseded-by the auditor-convergence programme (header :3-7), which is itself parked. The code it produced is live but has no product consumer: `compas_cgal.adaptive` is imported only by its own tests and one docstring (`engagement_toolpath.py:29` says it "deliberately does not use" it).
- **generation:** the idea is an engagement-certified, gouge-free, coverage-certified adaptive path. Placement is a finite candidate lattice of one-sided MATHSM circles on the true MAT; transactional candidate commit; neck-aware caps. **Quality: broken as a product.** It runs only on the hard-coded 6-vertex L fixture and never completed even that.
- **(a) generation goal:** "non-vacuous polygonal-pocket toolpath" with machine-checkable engagement, gouge-freedom and coverage (:14-16).
- **(b) landed (verified by artifact):** T0-T12. `src/continuous_tea_2/*`, `segment_site_*`, `containment_2*`, `reachable_domain_2*`, `coverage_2*` and `src/compas_cgal/adaptive/*` (17,794 lines) all exist. T13 landed partially as `generate_exact_adaptive_continuation` (`adaptive/generator.py:1336`), which "proves traversal exhaustion, not complete pocket coverage" (its own docstring).
- **(c) abandoned or superseded:** T11A-S (terminal seal), T14 (`certificate.py`, `schema.py`, `schemas/` absent), T14A (tri-dexel validation: no `test_tridexel_interchange.py`), T15 (no `scripts/run-adaptive-mutations.py`; `scripts/` holds only `engagement_baseline.py`), T16. The tracker role went to auditor P0-P4 (header).
- **(d) open generation items:** T13 terminal traversal and covered path (:2122-2263). Generalizing T9's MAT beyond the L fixture (never planned as a task; the plan assumed T9 was already general). The Held Fig. 5 fixture and runtime targets in T15 (:2611-2628).
- **defects:** :3-7 says "landed (this repository's leg)". *Contradicted*: T13-T16 are not landed, and T9's "true segment-site Voronoi MAT" is L-only (`src/segment_site_catalog_graph.cpp:159-193`). The effect on generated paths is that no path is ever produced outside one fixture.
- **unused capability / ideas:** neck classes with passage-state caps (T10). Exact containment of the annular sweep, not the filled outer disk, so annuli around islands are kept (:1679-1687). Zero-guide terminalisation at exact clearance = r (:1601-1606).
- **claims:** "109/121 promised artifacts exist" (:4): *unverified* (no recount). "T10 complete and gated" (:1512): *verified* (`adaptive/candidates.py`, `neck.py` exist). "Task 13F continuation integration pending" (:1848-1850): *verified* (declared red in `docs/red_manifest.json`). "All 518 adaptive tests pass" (:2174): *unverified*, historical.

### docs/superpowers/plans/2026-07-25-task3-single-arrangement-performance.md | lines 1489 | read 1-1489

- **role:** doc (exact reachable domain C_r and M_r, plus an incremental coverage ledger).
- **used by:** `tools/plan_headers.py`; master plan T3.
- **state:** live code, landed.
- **generation:** n/a (validation). The generation idea it records: exact `C_r = D ⊖ B_r`, `M_r = C_r ⊕ B_r` and `D \ M_r` (the unreachable residual) are built once, so exact rest-material regions are available.
- **(a) generation goal:** none.
- **(b) landed:** *verified*. `src/exact_region_2.*`, `exact_sweep_2.*`, `reachable_arrangement_2.*`, `reachable_certificate_2.*`, `coverage_bindings_2.cpp`, `canonical_rotation_2.h`, and the pixi task `task3-algorithm-gate` (`pyproject.toml:269`).
- **(c) abandoned:** none recorded. The post-Task-5 refactor gate (:1394-1489) had no recorded outcome.
- **(d) open generation items:** none.
- **defects:** none found that bear on generation.
- **unused capability / ideas:** `ReachableDomain2.unreachable_residual()` (`src/coverage_bindings_2.cpp:390`) is the exact convex-corner residual. The draft's D2 should be stated against it.
- **claims:** "24/24 promised artifacts exist; 0/46 boxes are noise" (:3-5): *unverified* count; the spot-checked files exist.

### docs/superpowers/plans/2026-07-28-atomic-candidate-transaction.md | lines 363 | read 1-363

- **role:** doc (generation infrastructure: fork, certify, deplete, commit of link/circle candidates).
- **used by:** `tools/plan_headers.py`; master T12.
- **state:** live (L only).
- **generation:** linking/placement transaction. Each candidate's link and circle are evaluated on forked exact stock and committed only on byte-identical replay. **Quality: adequate as a pattern**; it has no product consumer.
- **(a) goal:** "evaluate each link-and-circle candidate on isolated exact state" (:8).
- **(b) landed:** *verified*. `adaptive/generation_state.py`, `adaptive/transaction.py`, `tests/adaptive/test_transaction.py`. Commit cited at :360 (`a9fef135`).
- **(c) abandoned:** none.
- **(d) open generation items:** none in the plan. The causal neck scope was deferred to the causal traversal plan.
- **defects:** none found.
- **unused capability / ideas:** the draft's "productive-only emission" (decide with `Stock.exactly_equals` after `_subtract_operation`) is a degenerate case of this fork-then-commit pattern. Reusing the pattern's state ownership avoids a second stock owner.
- **claims:** "The only plan whose checkbox state matches reality" (:3-4): *verified* for this file (every box ticked, artifacts present).

### docs/superpowers/plans/2026-07-28-causal-mat-traversal.md | lines 541 | read 1-541

- **role:** doc (ordering and traversal over the exact MAT; neck-side causality; entry bootstrap).
- **used by:** `tools/plan_headers.py`; master T13; the zero-guide and retrace plans.
- **state:** live (L only), partially landed.
- **generation:** ordering. A directed global edge cursor and causal neck sides; the entry-circle launch; finite first-feasible search. **Quality: weak.** It runs only on the L fixture and stops before terminal coverage.
- **(a) goal:** "traverse every supported MAT edge, generate a terminal covered operation stream, and close it with fresh replay" (:8-10).
- **(b) landed:** Tasks 1-5 (`adaptive/neck_topology.py`, `traversal.py`, `traversal_graph.py`, `bootstrap.py`, `generator.py`). *verified* by `ls`.
- **(c) abandoned:** T6 Steps 2-5 (terminal seal), T7. The route-2 failure (:439-449) became the retrace plan, which was itself superseded.
- **(d) open generation items:** T6 terminal seal and `GenerationResult`; T7 performance baseline.
- **defects:** header :3-4 says "The 10 unticked boxes are unmaintained planning residue, not open work". *Contradicted*: T6 Steps 2-5 are open work (`class GenerationResult` is absent; `test_task13f_full_continuation` is red).
- **unused capability / ideas:** **depth-first edge-discovery order is not physical continuity** (:439-449). 56 of 56 route-2 links gouged from the route-1 endpoint, and reversing the route-1 motion restored 6 valid pairs. This bears directly on the draft's tree walk.
- **claims:** "affected 113, adaptive 518 pass" (:2174 in master, mirrored at :232): *unverified*, historical.

### docs/superpowers/plans/2026-07-29-one-root-full-circle-source.md | lines 346 | read 1-346

- **role:** doc (kernel: full-circle event oracle admitting one-root trim vertices).
- **used by:** `tools/plan_headers.py`; causal T6.
- **state:** live.
- **generation:** n/a (certifier).
- **(a) goal:** none (oracle completeness).
- **(b) landed:** T1-T3 and Task 4 Step 1. *verified* by `src/continuous_tea_2/circle_vertex_source.cpp`.
- **(c) abandoned:** Task 4 Steps 2-4 (cross-wiring mutations) and Task 5 Steps 2 and 5 are unticked. No evidence they ran.
- **(d) open generation items:** none.
- **defects:** header :3-4 says "unticked boxes are residue". *Unverified* for the Task 4 Step 2 mutations: I found no evidence either way.
- **unused capability / ideas:** none for generation.
- **claims:** "30 passed in 40.25 s; 522 passed in 242.54 s" (:257, :305): *unverified*, historical.

### docs/superpowers/plans/2026-07-29-zero-guide-link-candidate.md | lines 1407 | read 1-1407

- **role:** doc (generation: constant-clearance arms machined by an advancing centerline segment).
- **used by:** `tools/plan_headers.py`; causal T6.
- **state:** live (L only).
- **generation:** entry/advance through a channel exactly one tool diameter wide (clearance == r, so no trochoid fits). It emits a typed `AdvanceSegmentOperation` along the MAT centerline, certified by the swept-prefix theorem. **Quality: adequate in idea** (it is the honest answer to ρ ≤ r arms). It is L-fixture only.
- **(a) goal:** "machine [constant tool-clearance MAT spans] with a typed advancing segment" (:12-13).
- **(b) landed:** T1-T7 and T8 Steps 1, 4-6. *verified*: `src/segment_site_zero_guide.cpp`, `adaptive/advancing_segment_trial.py`, commits `8ef2981f` through `8f00952a` in `git log`.
- **(c) abandoned:** T8 Steps 2-3 (terminal seal, bounded path measurement).
- **(d) open generation items:** a complete Task 13F path. Inter-route transit, which was handed to the retrace plan.
- **defects:** none in the plan text. The header (:3-5, "landed") is true for T1-T7 only.
- **unused capability / ideas:** (1) **The swept-prefix theorem** (:875-890): a translation starting from a clear disk has at most the closed forward semicircle engaged, so it certifies exactly π; it is sealed by a native record. That is the cheap exact certificate for stay-down links through cleared stock (the draft's D3). (2) **Zero-guide classification by exact polynomial identity** (:26-30): the exact test for "this arm admits no trochoid", which the draft's convex-corner analysis (ρ ≤ r) could use. (3) The measured route-2 boundary (:1295-1306).
- **claims:** "generic event certification of that reverse remained active beyond 10 min ... swept-prefix certified the same motion in 0.040943 s" (:1306-1309; also causal :450-456): *unverified* (historical; not rerun).

### docs/superpowers/plans/2026-08-07-exact-inter-route-retrace.md | lines 2259 | read 1-2259

- **role:** doc (linking: exact cut-depth retrace across a non-incident route switch).
- **used by:** `tools/plan_headers.py`.
- **state:** superseded-by auditor P2 (header :3-9). P2 is parked. Tasks 1-5 code is live.
- **generation:** linking. A non-incident branch switch is crossed by reversing the last advancing segment at cut depth. The trigger is node-identity inequality, not coordinates. **Quality: adequate as a mechanism**: physically sound (it retraces cleared stock), and certified. It is L-only and narrow (the source must be a final no-neck zero-guide segment, :25).
- **(a) goal:** "cross the Task 13F route-1-to-route-2 discontinuity with one exact ... cut-depth retrace" (:13).
- **(b) landed:** T1-T5. *verified*: `adaptive/retrace_segment_trial.py`, `retrace_transaction.py`, `RouteRetraceCommit` in `generator.py:580`, commits `f5b045fe`, `988a2a51`, `6572b8c1`.
- **(c) abandoned:** T6 (`route_retrace_replay.py` is absent) and T7. The user ruling D1 on 2026-08-28 superseded it (from the backlog in `git show 2875cf8d^:docs/superpowers/state/backlog.md`). Its 4 red tests remain in `docs/red_manifest.json`.
- **(d) open generation items:** a general inter-route transit policy, beyond "reverse the final zero-guide segment".
- **defects:** header :5 says replay authority "is being built by the auditor programme's P2". *Contradicted*: P2 never ran (no `replay_certificate.py`).
- **unused capability / ideas:** the stay-down retrace is the direct precedent for the draft's D3 and an alternative to retract, link and plunge. The draft does not cite it.
- **claims:** "22/24 artifacts that did land remain in service" (:9): *unverified* count.

### docs/superpowers/plans/2026-08-19-benchmark-corpus.md | lines 3248 | read 1-3248

- **role:** doc (measurement corpus, MATHSM baseline, Fig. 6 reproduction).
- **used by:** `tools/plan_headers.py`; `benchmarks/*` (runner, families, mathsm, figure6, cli).
- **state:** live (landed).
- **generation:** n/a. It records generation ideas in its open items (:3162-3248).
- **(a) goal:** none (measurement: "is this fast enough? how often unresolved?").
- **(b) landed:** *verified* by `ls benchmarks`: `spec.py`, `families/`, `congruence.py`, `instrument.py`, `measurement.py`, `runner.py`, `report.py`, `mathsm.py`, `figure6.py`, `external.py`, `cli.py`.
- **(c) abandoned:** the "two runs that matter" (:3147-3152: external resolution rate and neck ablation). No doc under `docs/*.md` reports either (`grep "resolution rate|neck-degeneration ablation"` is empty).
- **(d) open generation items:** (1) a peak/count opposition at forced stations, as a docstring constraint (:3207-3215, landed as `docs/radius_ladder_ablation.md` per wave 1); (2) **the depletion model refuses Z+XY moves, so a true helical entry cannot be scored** (:3241-3245); (3) **arc vocabulary: inter-chain links cannot be tangent-continuous while loops are full circles** (:3246-3248).
- **defects:** header :7 says open items are "superseded by `docs/superpowers/state/backlog.md`, which is the single live backlog". *Contradicted*: that file was deleted in `2875cf8d` and `docs/superpowers/state/` does not exist on the live branch. The effect is that the C4 and C5 generation items have no live home except the draft's Stage B.
- **unused capability / ideas:** MATHSM spacing sweep (`benchmarks/mathsm.py`) as the comparison column; the neck pinch sweep for the draft's dumbbell review case.
- **claims:** "259 collected, 6 deliberate quality reds, 1 deliberate red" (:4-6): *contradicted* by the current manifest, which has 12 quality reds after the cap-40 attribution was added (`docs/red_manifest.json`). "Generation is already competitive (0.00-0.02 s)" (:35): *unverified*.

### docs/superpowers/plans/2026-08-23-auditor-convergence-p0-reconciliation.md | lines 367 | read 1-367

- **role:** doc (governance: branch reconciliation ledger).
- **used by:** `tools/plan_headers.py`, which exempts it by prefix (`DEFAULT_ALLOW_PREFIXES`, `tools/plan_headers.py:20`).
- **state:** landed (`449b45b0 docs(sdd): reconcile auditor branches`). Its ledger file was later deleted (`2875cf8d`).
- **generation:** n/a.
- **(a)-(d):** no generation goal; landed; ledger dropped; no open generation items.
- **defects:** no status header. Exempted by the lint, so there is no tool failure. The "nine load-bearing findings" (:276-288) were all certifier items.
- **unused capability / ideas:** none for generation.
- **claims:** "divergence 254 / 36" (:66-67): *unverified*.

### docs/superpowers/plans/2026-08-23-auditor-convergence-p1-truthful-audit.md | lines 1622 | read 1-1622

- **role:** doc (validation: truthful three-way engagement audit).
- **used by:** `tools/plan_headers.py`.
- **state:** partially landed (T1-T4); T5-T6 not.
- **generation:** n/a. It records one useful idea: exact partial-arc depletion on rational quarter charts (T3).
- **(a) goal:** none.
- **(b) landed:** *verified*: `src/compas_cgal/engagement_audit/` (`classification.py`, `identity.py`, `input.py`, `records.py`, `digests.py`, `decision_limits.py`); `src/audit_*_2.*`; `Stock2::subtract_exact_arc` (`src/stock_2.h:215`); commit `1e1dd528`.
- **(c) abandoned:** T5 (public `audit_toolpath_engagement` in `engagement_audit`; no `replay.py` or `report.py`) and T6 (the benchmark runner still imports the legacy `compas_cgal.engagement.audit_toolpath_engagement`, `benchmarks/runner.py:21`; no `docs/auditor_convergence.md`).
- **(d) open generation items:** none.
- **defects:** none for generation. The gate still uses the legacy sampled audit and survey, which this plan meant to replace. That is a validation concern.
- **unused capability / ideas:** `subtract_exact_arc` for scoring arc links (backlog C5).
- **claims:** none checked beyond the artifacts.

### docs/superpowers/plans/2026-08-23-auditor-convergence-p2-replay-theorem.md | lines 386 | read 1-386

- **role:** doc (validation, plus one generation task: terminal continuation).
- **used by:** `tools/plan_headers.py`; cited as the successor by the retrace header.
- **state:** never executed (parked by the draft).
- **generation:** T1, "reach terminal traversal" on the L fixture. **Quality: n/a** (never run).
- **(a) goal:** a terminal `GenerationContinuation` (T1).
- **(b) landed:** nothing. *verified*: no `test_terminal_continuation.py`, no `replay_certificate.py`, no `docs/swept_prefix_theorem.md`.
- **(c) abandoned:** all five tasks.
- **(d) open generation items:** T1 terminal traversal. T4 "prove or disable the swept-prefix fast path" (:274-325). The draft's D3 should inherit this theorem obligation if it uses swept-prefix certification.
- **defects:** none.
- **unused capability / ideas:** T4's requirement to separate "previously swept capsule" from "merely clear start disk" (:298). This is the correct precondition for any stay-down link certificate.
- **claims:** none checkable.

### docs/superpowers/plans/2026-08-23-auditor-convergence-p3-generator-qualification.md | lines 357 | read 1-357

- **role:** doc (the only auditor plan that touches the gated generators: engagement control and spacing).
- **used by:** `tools/plan_headers.py`.
- **state:** never executed.
- **generation:** T2 replaces bisection with a finite ladder search across all regulated generators; T1 makes forced engagement decisions structured outputs; T3 types the guide geometry. **Quality: n/a** (not run). The underlying defect is live.
- **(a) goal:** "replace monotonic search assumptions with one finite ladder authority, preserve every forced engagement decision structurally" (:8-11).
- **(b) landed:** nothing. *verified*: `engagement_override.py`, `engagement_search.py` and `engagement_geometry.py` are absent from `src/compas_cgal/`; `engagement_toolpath.py:601-655` still bisects.
- **(c) abandoned:** all five tasks. P3 T5's "a second cap in 40-100 separates the generators" was done another way: wave 1 T10 added cap 40.
- **(d) open generation items:** **T2 (ladder search in `engagement_controlled_toolpath`)** and T1 (override ledger). T3's "known rho-generator rotated-radius-as-tangent construction" (:216) names a defect in `engagement_rho_toolpath.py` that the draft reuses. *unverified* whether it still exists; the rho slice should check it.
- **defects:** none in the plan. Live defect it names: `src/compas_cgal/engagement_toolpath.py:601-655`. Bisection under an admitted-false monotonicity returns an admissible but not maximal advance, which means extra loops and a longer path on the gated generator.
- **unused capability / ideas:** T4's bridge-exemption concern does **not** apply to the gate's tangent-break criterion. That criterion has no exemption (`benchmarks/quality_observations.py:913-923` checks every adjacent motion pair). *verified*.
- **claims:** "the 7.7% structural exemption" (:271): applies to `tests/test_toolpath.py` helpers, not the gate. *unverified* there.

### docs/superpowers/plans/2026-08-23-auditor-convergence-p4-enforcement-evidence.md | lines 396 | read 1-396

- **role:** doc (governance: Python 3.12 floor, build identity, CI, evidence bundle).
- **used by:** `tools/plan_headers.py`.
- **state:** never executed.
- **generation:** n/a.
- **(a)-(d):** no generation goal. Nothing landed: `requires-python = ">=3.9"` (`pyproject.toml:9`), ruff `py39` (:130), no `tools/build_identity_manifest.py`, no `benchmarks/evidence.py`, and `build.yml` still uses conda. No generation items.
- **defects:** none for generation.
- **unused capability / ideas:** none.
- **claims:** none.

### docs/superpowers/plans/2026-08-28-coherence-wave1.md | lines 824 | read 1-824

- **role:** doc (governance: red manifest, plan headers, claims ledger, attribution test).
- **used by:** `tools/plan_headers.py`; `pixi run red-manifest` / `plan-headers` (`pyproject.toml:230-231`).
- **state:** landed. The measured-run and claims subsystem was later removed (`6b4293ae drop measurement-claim subsystem`; `tools/measured_run.py` is absent).
- **generation:** n/a. It records one generation constraint: "at a forced station every candidate is over-cap ... lowering the peak always costs a circle" (:711-714).
- **(a)-(d):** no goal. Landed: `tools/red_manifest.py`, `tools/plan_headers.py`, `docs/red_manifest.json`, `docs/radius_ladder_ablation.md`, `tests/benchmarks/test_gate_attribution.py`. Abandoned: the claims subsystem (removed). Open: none.
- **defects:** none for generation.
- **unused capability / ideas:** the cap-40 attribution (T10 execution correction, :731-741) is the reason the gate has 12 cells, not 6.
- **claims:** "red set == manifest" at close (:3-7): *unverified* now. The draft's Task 0 re-runs it.

### docs/superpowers/plans/2026-08-29-wave1-task6-historical-assertions.md | lines 181 | read 1-181

- **role:** doc (governance: provenance artifact for the deleted measurement claims).
- **used by:** `tools/plan_headers.py`.
- **state:** superseded-by `6b4293ae`. `benchmarks/measurement_claim_history/` is absent.
- **generation:** n/a.
- **(a)-(d):** none. It landed, then its artifacts were removed. No generation items.
- **defects:** header :3 says "landed" with no note that the artifacts were later removed. Governance only.
- **unused capability / ideas:** none.
- **claims:** none checkable.

### docs/superpowers/plans/2026-08-31-held-pfeiffer-reference-pockets.md | lines 896 | read 1-896

- **role:** doc (measurement inputs: four Held-Pfeiffer pockets reconstructed from publisher PDF vectors).
- **used by:** `tools/plan_headers.py`; `benchmarks/held_reference_cases.py` and the gallery tasks.
- **state:** live (complete).
- **generation:** n/a. Its "next phase" (:865-896) is a generation plan that never started.
- **(a) goal:** the inputs; the next phase targets generator repair and G-code on the four pockets.
- **(b) landed:** *verified*: `benchmarks/held_reference_geometry.py`, `tools/held_reference_extractor.py`, `benchmarks/data/held_pfeiffer_2025/`, `docs/benchmarks/held_reference_qualification.md` (all four cases report "passed" = non-empty path only).
- **(c) abandoned:** nothing inside the plan.
- **(d) open generation items:** :885-894, in order: per-case residual/gouge/engagement/degeneracy/redundancy/continuity evidence for unchanged paths; generator repairs; fail-closed machine-dialect G-code; parser round-trip and machine simulation; matched comparison. **None started.** No G-code code exists (only a comment, `src/compas_cgal/toolpath.py:67`).
- **defects:** none.
- **unused capability / ideas:** the four pockets are the right Stage C generalisation set (the draft lists them). The "protected oracle / repairable policy layer" split (:872-883) is the right governance frame for the draft.
- **claims:** "qualification passed" (`docs/benchmarks/held_reference_qualification.md`): *verified* as stated, and it means only that the path is non-empty.

### docs/superpowers/plans/2026-09-02-held-figure5-2d-path-characterization.md | lines 1800 | read 1-1800

- **role:** doc (validation: one canonical quality reducer; Fig. 5 characterization of the default generator).
- **used by:** `tools/plan_headers.py`; `benchmarks/quality_observations.py`, `held_path_*`, `held_consumer_adapters.py`.
- **state:** stopped on 2026-09-05 (header :3-8). Tasks 1-10, 5A.1-5A.4, 11A and 11B landed.
- **generation:** n/a. It declares generator repair out of scope (:1798-1800).
- **(a) goal:** a refusal gate for postprocessor qualification (validation, not generation).
- **(b) landed:** *verified*: `benchmarks/units.py`, `held_path_snapshot.py`, `quality_observations.py`, `held_path_evidence.py`, `held_post_qualification.py`, `held_path_characterize.py`, `held_path_report.py`, `held_consumer_adapters.py`, `src/cutter_centre_domain_2.cpp`, `src/reachable_material_predicate_2.*`, `minimum_coverage_grid` (`benchmarks/coverage.py:58`).
- **(c) abandoned:** Task 11 Steps 3-7 (the live Fig. 5 report): "characterization had become a detour" (:3-8).
- **(d) open generation items:** the Phase 2 generator-policy design (:1798-1800), never written.
- **defects:** none in the plan text.
- **unused capability / ideas:** (1) `minimum_coverage_grid` plus the Held adapter, which the draft's Task 0b needs. (2) Measured blockers for Task 0b's budget: survey and coverage on Fig. 5 ran 780.83 s before `CoarseCoverageGridError` (:1555-1561). (3) `ReachableMaterialPredicate2`: exact membership `q ∈ C or dist(q, ∂C) ≤ r`, 305,944 queries in 34 s (:1693-1698).
- **claims:** "Figure 5 generation completed in 0.063150 seconds and produced 2,289 operations" (:1498-1499): *unverified*.

### docs/superpowers/plans/2026-09-05-held-paper-figure-reproduction.md | lines 770 | read 1-770

- **role:** doc (generation in the Held lane: placement, spacing, corners, contour-aware contact; plus measurement).
- **used by:** `tools/plan_headers.py`; `docs/segment_site_mat.md:1361`; pixi tasks `held-toolpath-figures`, `held-toolpath-gate` (`pyproject.toml:220-222`).
- **state:** in progress per header; parked by the draft.
- **generation:** the most complete toolpath generation in the slice. Paper-order q/c/ρ placement against an 80° predecessor cap. Polygon refinement repaired 108 over-cap successors (2,013 circles, maximum 79.992°, :188-199). Exact connectors. Corner approaches (Monstera 4,367 circles, maximum 79.998°, :632-641). Boundary-order spacing (upper 4,183 → 1,684 circles, :653-662). Native q-m circles on curved boundaries. **Quality: adequate as a reproduction, weak as a product.** It needs a publisher-derived start. Entry, connector engagement and depleted-stock qualification are open (:197-199). It is slow: Fig. 5 20.6 s, Monstera 436 s (:541-543), against Held's published 3-100 ms.
- **(a) goal:** "repository-generated reproduction drafts for every figure and panel", and an independent continuous Fig. 5 path (:16-17, T5 gate :231-234).
- **(b) landed:** T1-T4 and many T8 checkpoints. *verified*: gallery `docs/held_toolpath_gallery.md`; `benchmarks/held_toolpath_benchmark.py`; `benchmarks/results/held_toolpath_baseline.json`, `held_curved_circle_baseline.json`.
- **(c) abandoned or rejected:** exact full-motion coverage replay "set aside" (:523-533, :550-565). Dense-source coverage preservation was rejected (+161% length, :589-594).
- **(d) open generation items:** T5 (publisher-free MA traversal; continuous path; standard-predecessor placement with no forced successors), T6 (contour-aware and MATHSM peers), T7 (stage accounting and one optimization), T8 (all panels; coverage acceptance; zero-guide event for upper and Monstera, :567-574; contour-aware drafts, :610-620). The local reachable-residual repair (:579-588) is open.
- **defects:** :380-386 itself records a defect: the drafts use straight-skeleton guides and lose Held's MA incidence, and the coverage-thinning regression fails at every scale. That is directly relevant to the draft.
- **unused capability / ideas:** native boundary-normal q-m circle (:393-411); `HeldDiskContour2` (:610-620); corner approaches (:632-641); boundary-order spacing (:653-662); conditioning-not-degeneracy lessons (:470-504).
- **claims:** "first Figure 5 query from 337 s to 8 ms" (:11): *unverified*. "The current exposed native MAT builder supports only its canonical L-shape fixture" (:387-389): *verified* (`src/segment_site_catalog_graph.cpp:159-193`).

### docs/superpowers/plans/2026-09-08-number-type-coherence-stage-0.md | lines 879 | read 1-879

- **role:** doc (kernel plumbing: `src/exact/` vocabulary and the binary64 and canonical "doors").
- **used by:** `tools/plan_headers.py`.
- **state:** landed (`ls src/exact` shows `rational.*`, `canonical.*`, `one_root.*`, `errors.h`). Parked thereafter by user ruling. No status header, which draft Task 0 fixes.
- **generation:** n/a (no toolpath effect by design).
- **(a)-(d):** no goal; landed; nothing abandoned (completion criterion 2 is formally unmet, :836-860); no generation items.
- **defects:** none for generation.
- **unused capability / ideas:** none.
- **claims:** "5 failed, 786 passed; 26 failed, 3183 passed, 1 error" (:843-844): *unverified*, historical.

### docs/superpowers/plans/2026-09-08-number-type-coherence-stage-1.md | lines 136 | read 1-136

- **role:** doc (kernel: one `sign_mixed_radical` definition).
- **used by:** `tools/plan_headers.py`.
- **state:** landed. *verified*: `src/exact/one_root.h:89`; `src/audit_exact_station_2.cpp:245` delegates to it; no definition left in `src/engagement_2.cpp`.
- **generation:** n/a.
- **(a)-(d):** none; landed; none; none.
- **defects:** none.
- **unused capability / ideas:** none.
- **claims:** none checkable.

### docs/superpowers/plans/2026-09-08-number-type-coherence-stage-2.md | lines 1676 | read 1-1676

- **role:** doc (kernel: station event source carries exact rationals).
- **used by:** `tools/plan_headers.py`; `docs/number_types.md:503`.
- **state:** landed (`003e144c`, `6d879769`; `tests/build/test_station_lane_decoder_ratchet.py`; `tests/native/exact_station_attestation_gate.cpp`).
- **generation:** n/a.
- **(a)-(d):** none; landed; none; none.
- **defects:** none.
- **unused capability / ideas:** none.
- **claims:** "integer 0.149 ms vs generic 0.208 ms, only 1.4×" (:1205-1215): *unverified*.

### docs/superpowers/plans/2026-09-08-number-type-coherence-stage-3.md | lines 2872 | read 1-2872

- **role:** doc (kernel: segment event source inversion).
- **used by:** `tools/plan_headers.py`.
- **state:** parked. Only Task 1 (the byte anchor) landed: `tests/native/exact_segment_attestation_gate.cpp` exists. `src/exact/canonical.h` still has `canonical_bytes` (2 hits), and `segment_source.cpp` still has `lift_binary64`/`exact_binary64` (3 hits). Matches the draft's proposed header.
- **generation:** n/a.
- **(a)-(d):** none; T1 only; T2-T9 parked (user ruling "does not change a toolpath"); none.
- **defects:** none.
- **unused capability / ideas:** none.
- **claims:** "the real variable is arc-crossing, not degenerate sweep" (:2131-2156): *unverified* here (a kernel claim).

---

## Findings register

```
F-d2_plans-1 | src/segment_site_catalog_graph.cpp:159-193; src/segment_site_mat_bundle.cpp:39 | exact segment-site MAT and certified adaptive generator accept only the hard-coded 6-vertex L fixture; cannot run on any gate pocket; rectangle MAT graph exists (:905) with no caller | open (evidence: code read) | y
F-d2_plans-2 | 2026-07-24 phase1 :14-16, :2713-2731 | phase-1 goal (non-vacuous certified covered path) never reached; T13 partial, T11A-S/T14/T14A/T15/T16 absent | open (evidence: no GenerationResult/certificate.py/schemas; test_task13f_full_continuation in docs/red_manifest.json) | y
F-d2_plans-3 | 2026-07-28 causal :439-449; zero-guide :1295-1306 | DFS route order is not a cutter trajectory; non-incident branch switches need transit; draft tree walk must budget for it | open (evidence: plan measurement; red test) | y
F-d2_plans-4 | zero-guide :875-890 | swept-prefix theorem (clear-start translation engages at most π) certifies stay-down retraces in ~0.04 s vs >10 min generic oracle; candidate certificate for the draft's D3 | open (P2 T4 proof obligation never discharged) | y
F-d2_plans-5 | retrace plan T1-T5 | stay-down reverse retrace through the swept corridor exists (L only) as the alternative to retract/link/plunge; not cited by the draft | open | y
F-d2_plans-6 | P3 T2; src/compas_cgal/engagement_toolpath.py:601-655 | gated generator bisects advances under an admitted-false monotonicity; planned ladder fix never executed; draft must pin "no bisection" as a constraint | open (evidence: code read) | y
F-d2_plans-7 | P3 T1 | forced engagement decisions are warnings, not structured outputs; draft's virgin_plunges covers entries only | open | y
F-d2_plans-8 | paper plan :380-386 | straight-skeleton guides break Held's MA-incidence coverage assumption; coverage-thinning regression failed at every scale | open | y
F-d2_plans-9 | paper plan :393-411; src/boundary_normal_circle_2.cpp | generic native boundary-normal q-m circle construction (true-MA incidence, reflex sectors) exists for arbitrary simple polygons; the skeleton alternative to evaluate in Stage C | open | y
F-d2_plans-10 | paper plan :632-641 | corner-approach motion (approach / shared-contact rotation / departure) solved Monstera's 111.59° concave-corner jump at an 80° cap; candidate for draft B2 | open (not used by gated generators) | y
F-d2_plans-11 | paper plan :610-620; src/held_disk_contour_2.cpp | contour-aware contact query exists; contour-aware placement never built | open | y
F-d2_plans-12 | paper plan :180-256 | Held T5 (publisher-free continuous Fig. 5 path) and T6 (three-curve Fig. 6) never reached | open | y
F-d2_plans-13 | paper plan :541-543 | Held-lane generation takes 20.6 s (Fig. 5) to 436 s (Monstera); any "generalise the Held pipeline" retarget in draft Task 0b inherits this | open | y
F-d2_plans-14 | characterization :1553-1567; benchmarks/held_consumer_adapters.py:35 | Fig. 5 needs coverage grid >= 668x458 or CoarseCoverageGridError; draft Task 0b must use the Held adapter; 600 s budget unproven | open | y
F-d2_plans-15 | paper plan :196-199, :596-600 | Held refined path has 1,010 near-tangent containment rejections and an unreachable 2.04 µm remnant; will hit gate gouging/uncut criteria in Task 0b | open | y
F-d2_plans-16 | reference-pockets :865-896 | next phase (generator repair, fail-closed G-code, parser round-trip, machine simulation) never started; no G-code emitter exists | open (evidence: grep) | y
F-d2_plans-17 | corpus :3241-3248 | helical entry unscoreable (depletion refuses Z+XY) and arc link vocabulary (C5) both open; draft covers D1 only, not C5 | open | y
F-d2_plans-18 | task3 plan; src/coverage_bindings_2.cpp:390; src/reachable_material_predicate_2.h:15 | exact D \ M_r and fast exact reachable-material membership exist; base the draft's D2(a) rest-material ruling on them | open | y
F-d2_plans-19 | phase1 header :3-7 | "landed" header false for T13-T16 and for T9 generality | open | n
F-d2_plans-20 | causal header :3-4 | "unticked boxes are residue, not open work" false (T6 terminal seal open) | open | n
F-d2_plans-21 | corpus header :7; retrace header :3-9 | point at docs/superpowers/state/backlog.md (deleted in 2875cf8d) and at P2 (never run) as the live owners | open | n
F-d2_plans-22 | wave1-task6 header :3 | "landed" omits that its artifacts were removed in 6b4293ae | open | n
F-d2_plans-23 | effort table (this file) | 164 tasks: G 33 / V 43 / M 37 / Gov 28 / K 23; zero G tasks on a gated generator; 29 of 33 G tasks on the L-only adaptive line | closed (evidence: this census) | y
F-d2_plans-24 | P3 T4; benchmarks/quality_observations.py:913-923 | gate tangent-break criterion has no structural bridge exemption; the non-vacuity concern does not apply to the gate | closed (evidence: code read) | y
F-d2_plans-25 | P1 T5-T6; benchmarks/runner.py:21 | truthful three-way audit never wired; gate and corpus still use legacy sampled audit and survey | open | n
F-d2_plans-26 | corpus :3147-3152 | external-corpus resolution rate and neck ablation runs never reported | open (evidence: grep docs) | n
F-d2_plans-27 | characterization :1798-1800; paper T7 | no Phase 2 generator-policy design and no stage-accounted performance work on any generator | open | y
F-d2_plans-28 | zero-guide :26-30 | exact zero-guide (clearance == r) classification is the exact test for "no trochoid fits"; relevant to the draft's ρ <= r convex-corner analysis but L-only | open | y
```
