# Census slice g4_cpp_kernels: the C++ kernels the generators decide with

Read in full: all 126 assigned files (24,042 lines). Probes were run with `pixi run python -` (stdin scripts, nothing written to the repo). Every number below comes from a file:line or from one of those probes, and the probes are named where they are used.

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **The exact removal test is exact about an inexact model, so "removes material" does not mean that any material gets cut.** (verified, probe P1/P2)
   - The gate replay and the draft plan both deplete through `compas_cgal.engagement._subtract_operation` (engagement.py:186-221):
     - a `Line` becomes `Stock2::subtract_capsule`, a disk chain spaced 0.02 r that under-covers by up to 1e-4 r (stock_2.cpp:516-540, :352);
     - a `Circle` becomes an annulus whose guide radius is `hypot(point_at(0) - c)` (stock_2.cpp:673-697), a double that is not the authored `Circle.radius`.
   - Measured on a 12×8 rectangle, r = 1 (P1):
     - a capsule (3,4)→(7,4) lying entirely inside an earlier capsule (2,4)→(8,4) still "changes the stock": True for the chain and also True for `subtract_capsule_quad`;
     - re-running the same nominal loop (same centre, ρ = 1.7) from a different start angle changes the stock in 30 of 39 trials.
   - Consequence: any op that only grazes model slivers counts as productive. The gate's "redundant operations" criterion under-counts air-cutting, and a generator steered by this test cannot tell when a walk has reached cleared stock.
2. **A plunge at a loop's own entry point is never exactly clear.** (verified, P2: 5 of 5 angles)
   - After a full loop centred at (10,6) with ρ = 1.7, `subtract_disk(entry, r)` still removes material at every start angle tried.
   - Cause: the tool disk at the entry is internally tangent to both annulus circles, and the entry point is ulps off the surrogate radius.
   - This contradicts the plan's `_enter` premise: "a full circle's sweep contains the tool disk at its own entry".
   - Every re-entry in the tree walk will be counted as a virgin plunge, so `virgin == 1` should fail on the rectangle and the L (prediction from the measured primitive; Task 3 was not run).
3. **Lines replayed as disk chains blow up the stock arrangement and the cost of every later exact decision.** (verified, P3/P4)
   - One 12 mm chain capsule takes the arrangement from 80 to 1,282 vertices; the quad version gives 86.
   - Four parallel 16 mm chain runs spaced 2r apart leave 8,014 vertices and 1,606 faces, i.e. hundreds of sliver islands between the runs.
   - On that stock one removal test (clone + subtract + `exactly_equals`) costs 0.378 s for a loop and 0.398 s for a bridge. On a clean 80-vertex stock it costs 0.003 s and 0.013 s.
   - The plan's bridges and the Stage B2a wall run-ins are lines.
4. **The analytic TEA-growth lemma behind `certify_segment_tea` and the gate's arc certifier is false for concave boundaries.** (verified, P5)
   - Mechanism (a), "a run end drifts by at most the arc of chord d" (engagement_2.cpp:356-360), is false near tangency, where drift grows like √d.
   - Probe: cutter (r = 1) inside a cleared hole of radius 1.01, moved by d = 0.001 r. TEA jumped from 35.7° to 59.5° (+23.8°). The claimed bound is GROWTH(d) = 5.24°, and the 2× guard is 10.5°.
   - No false certificate was constructed, so the certificate's soundness is unproven rather than shown broken.
5. **Engagement is capped per run, not as total contact.** (verified by code)
   - `run_exceeds_cap` is applied to each maximal run (exact_cap_certificate_2.h:84-140; audit_exact_station_2.cpp:410-423).
   - A station with two runs of 100° each passes a 120° cap with 200° of rim in material.
   - Two-sided contact, such as bridges through partly cleared corridors or walls at both sides, is admitted.
6. **The generators sample; continuous certification exists but is unusable at the cap.** (verified)
   - The gated generators decide with `engagement_at` at probe positions (engagement_toolpath.py:552).
   - `certify_segment_tea` (engagement_2.cpp:445-500) cannot certify any motion within about 5–7° of the cap: the station floor is 1e-3 r and the residual guard there is 0.090–0.128 rad (arithmetic from :407-438).
   - Measured cost is 511 stations and 0.163 s for one 2 mm segment at 92° under a 120° cap (P4).
   - A partial arc is certified only if its whole circle certifies (audit_certification_2.cpp:394-450). Otherwise it is UNRESOLVED.
7. **Dead code sits in the deciding module.** (verified)
   - Since `engagement_at` began delegating to `classify_audit_unguarded_station_exact` (engagement_2.cpp:532), the following have no caller: `finish_engagement`, `pessimistic_runs`, `engaged_arcs_zone`, `EngagementVisitor`, `Arc` and the `PL` alias (engagement_2.cpp:29-36, 86-345).
   - A doc comment at :86 has lost its first line.
8. **Latent lifetime hazard: by-value `Gps` returns rely on NRVO, which C++ does not guarantee.** (unverified that it ever fires)
   - Affected sites:
     - `exact_annulus_region`, `exact_disk_union` (stock_2.cpp:81-86, 114-122);
     - `build_exact_disk_union_region_2` (exact_disk_region_2.cpp:24-31);
     - `reach_join_parts` and `reach_full_circle_sweep` (exact_sweep_2.cpp:224-231, 261-268);
     - `disk_set` and `capsule_set` (containment_2.cpp:76-85).
   - The headers state the hazard themselves (exact_disk_region_2.h:18-23): a copied `Gps` reads the traits of a dying source. Clang and GCC apply NRVO at these sites in practice.

### Generation capabilities in this slice that the gated generators do not use

The gated generators use only three of these kernels: `engagement_at` (sampled), `Stock2::subtract_capsule_quad`, and `subtract_arc_sweep_local`.

| capability | file:line | what it would give generation |
| --- | --- | --- |
| Exact reachable material and **unreachable residual** (design \ (centre domain ⊕ r)) | reachable_domain_2.cpp:119-190 | The exact corner material no tool position can reach. This is the honest target for D2 and for "uncut". |
| Exact centre-domain (wall-offset) boundary cycle and CCW boundary transitions | reachable_arrangement_2.cpp:673-890; coverage_bindings_2.cpp:225-260 | The exact wall-offset path that B2a's run-ins need, plus a contour-parallel finishing / rest-material pass. |
| Exact boundary-normal medial disk at any boundary point (lines + arcs) | native_boundary_medial_2.cpp:78-274 | Held's contour-aware circle placement and corner circles, without the straight skeleton. |
| Exact gouge predicates for segments, loops, entry disk, and loop-inside-entry-disk | containment_2.cpp:225-356 | D1 helical entry: certify that the first loop lies inside the bored disk (`evaluate_exact_full_circle_in_disk`). |
| Coverage with a pre-cleared entry disk; exact remaining material over a CORE √-kernel (true capsules, true annuli) | coverage_2.cpp:172-259; remaining_material_2.cpp:11-85 | An exact removal and coverage model with no under-cover slivers. It is slow and was set aside 2026-09-08. |
| Exact circle-redundancy test `can_remove_circle` | stock_circle_redundancy_2.cpp:5-53 | Post-pass deletion of loops whose material neighbours already remove, conservative in the safe direction. |
| Exact whole-orbit engagement bound and exposed-contour contact (`HeldDiskContour2`) | held_contour_clearance_2.cpp:46-88; held_disk_contour_2.cpp:76-113 | Placement-time engagement bound with no sampling; used only by Held benchmarks. |
| Continuous TEA certificate for segments; full-circle event authority via `certify_audit_tea_exact` | engagement_2.cpp:560-578; audit_certification_2.cpp:329-377 | Exact loop admissibility in place of probe sampling. The circle authority lives in continuous_tea_2 (another slice). |
| Exact rational points on a circle (quarter-chart Pythagorean atlas) | exact_circle_chart_2.cpp:38-65 | Entry and bridge endpoints exactly incident on the guide circle, which is the cure for finding 2. |
| Local depletion (arrangement surgery) | stock_local_2.cpp:290-324 | Already used for loops (engagement_toolpath.py:795). The gate replay does not use it. |

### Where this slice shows the draft plan is wrong or incomplete

- **Design, "Productive-only emission … 0 redundant operations holds by construction": contradicted in substance.**
  - It holds only because the generator and the gate share one sliver-polluted test (P1).
  - The walk's `_first_productive` / `_next_step` will treat ulp and chain slivers as material. That defeats "emit nothing on a cleared cycle-closing edge" (Review Focus 5) and emits near-air loops that the gate then accepts.
- **Task 3 `_enter` docstring: contradicted** (P2, 5/5). The tool disk at a loop entry is not inside that loop's replayed sweep.
- **Global constraint "removal, plunge clearance … are exact decisions": incomplete.** Exactness of the predicate does not compensate for an under-covering, surrogate-radius model. One of these two is needed:
  - **(a) Fix the model.** Change the replay's `Circle` branch (engagement.py:195-207) to use the authored radius: `subtract_circle_sweep(c, g.radius, r)` (stock_2.cpp:408-421). Make the generator author loop entries and bridge endpoints exactly incident, for example axis-aligned phase on a dyadic authoring grid, or exact chart points.
  - **(b) Make the decision three-valued.** Keep an exact under-cover / over-cover stock bracket, so that "sliver-only" is its own verdict. No area epsilon.
- **Affordability (Task 3): conditional.**
  - It is affordable on a clean stock: 3–13 ms per test (P2).
  - It is not affordable once lines are replayed as chains: about 0.4 s per test after four wall runs (P3).
  - The fix goes in the gate replay and in the generator together: quad (6 curves), or exact capsules where the segment direction is Pythagorean.
  - A cheaper exact equivalent of `not clone-subtract-equals` is "sweep ∩ stock has non-empty interior": one boolean or a local zone query instead of clone + 3 overlays (unverified speed-up).
- **Stage B2 / D2 ignores existing kernels.**
  - The exact unreachable corner residual (`ReachableDomain2.unreachable_residual`) and the exact wall-offset cycle (`build_center_boundary_cycle`, `ccw_transition`) already exist.
  - B2a's "run in along each wall offset" can be built on them, and D2(a)'s "measure uncut only over what is reachable" should be measured against `reachable_material` (tool reach) rather than a generator-defined trochoid reach.
- **Stage B / D1 ignores existing kernels.**
  - An exact pre-cleared entry disk (coverage_2.cpp:190-259) and loop-in-entry-disk containment (containment_2.cpp:322-356) already exist.
  - A HELIX op could be replayed as `subtract_circle_sweep(c, ρ_h, r)` (an exact disk when ρ_h ≤ r, stock_2.cpp:155-163) and certified with those predicates.
- **Review Focus 2 (dumbbell): the plan does not account for a kernel limit.**
  - The reachable/coverage kernels throw `PocketNotMachinableError` whenever the centre domain has more than one component (reachable_arrangement_2.cpp:997-1006; reachable_material_predicate_2.cpp:277-284).
  - The plan's dumbbell (3 mm neck, r = 1) stays connected. A neck under 2r would disconnect it, and any coverage measurement through these kernels would then fail loud.

### Probes (all `pixi run python -`, stdin only)

- **P1:** 12×8 rectangle, r = 1.
  - Inner collinear chain capsule changes the stock: True.
  - Inner collinear quad capsule changes the stock: True.
  - Same nominal loop from another start angle: changed in 30/39 trials.
- **P2:** 20×12 rectangle; full loop at (10,6), ρ = 1.7.
  - Disk at its own entry removes material: True for angles 0.3, 1.1, 2.0, 2.9, 4.4.
  - A trochoid row of 38 loops plus chain bridges leaves an (80, 160, 3) arrangement.
  - Removal test: loop 0.003 s, 2 mm bridge 0.013 s.
  - clone 0.0001 s; `exactly_equals` 0.002 s.
- **P3:** four 16 mm chain runs at y = 1.5, 3.5, 8.5, 10.5 give an (8014, 19232, 1606) arrangement. Removal test: loop 0.378 s, bridge 0.398 s.
- **P4:** on the P2 stock:
  - `engagement_at` 0.176 ms/call.
  - `certify_segment_tea` over 2 mm: (91.9°, True, 511 stations), 0.163 s.
  - A 12 mm chain capsule takes the arrangement to (1282, 2564, 4) in 0.065 s; the quad takes it to (86, 172, 4) in 0.0008 s.
- **P5:** stock minus a disk of R = 1.01 at (6,4); cutter at offsets 0.0105 → 0.0115: max-run TEA 35.69° → 59.51°. GROWTH(0.001) = 5.24°; 2·GROWTH = 10.48°.

---

### src/audit_arc_motion_2.cpp | lines 370 | read 1-370
- **role:** validation (exact arc identity for the audit lane)
- **used by:** _stock_2 (audit classification); tests/engagement_audit; no generator
- **state:** live (audit lane)
- **generation:** n/a. It normalises authored binary64 angles to exact quarter-chart parameters (:59-92) and builds chart intervals (:106-162).
- **defects:**
  - :77-80: the chart parameter is `tan(local/2)` evaluated in double and then injected. The arc's exact endpoints are therefore a surrogate of the authored angle, as the comment at :276 admits.
  - A generator's arc whose endpoint sits on a stock vertex is replayed off that vertex. Effect: audit-only.
- **unused capability / ideas:** none for generation.
- **claims:** "audit-arc-quarter-chart-binary64-v2" surrogate semantics: verified (:368).

### src/audit_arc_motion_2.h | lines 154 | read 1-154
- **role:** validation
- **used by:** audit_* units, exact_depletion_2.h
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_authority_evidence_2.cpp | lines 240 | read 1-240
- **role:** validation (certified-coverage witness minting)
- **used by:** audit_certification_2.cpp
- **state:** live
- **generation:** n/a
- **defects:**
  - :162-166 / :195-199 / :231-235: the full event authority is re-run in order to mint, and is replayed again at audit_certification_replay_2.cpp:386-422. Cost only.
- **unused capability / ideas:** none
- **claims:** none

### src/audit_certification_2.cpp | lines 502 | read 1-502
- **role:** validation (exact per-motion TEA verdict)
- **used by:** audit_replay_mutating_2.cpp:20; tests
- **state:** live (audit lane; not on any generator path)
- **generation:** n/a
- **defects:**
  - :394-450: a partial arc certifies only through its full-circle superset. Otherwise bisection refinement (:189-269) can find a violation or end UNRESOLVED, but can never certify.
  - Effect: every near-cap partial arc a generator emits (corner arcs, bridges as arcs) audits UNRESOLVED.
- **unused capability / ideas:** `certify_audit_tea_exact(stock, AuditCircleMotion2,…)` (:329-377) is an exact continuous loop verdict a generator could decide with. It requires an axis-aligned exact frame (audit_classification_core_2.cpp:174-201).
- **claims:** none

### src/audit_certification_2.h | lines 508 | read 1-508
- **role:** validation
- **used by:** engagement_2.cpp, exact_cap_certificate_2.h, audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `classify_audit_unguarded_station_exact` (:161-166) also returns the exact rim crossings and runs, which could seed exact engagement-aware placement.
- **claims:** none

### src/audit_certification_internal_2.h | lines 284 | read 1-284
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** :125-127 "no raw bytes or forged public constructor can mint generated decision evidence": verified (private constructors, friend factory).

### src/audit_certification_replay_2.cpp | lines 723 | read 1-723
- **role:** validation (independent replay of decision witnesses)
- **used by:** audit_replay_mutating_2.cpp:34
- **state:** live
- **generation:** n/a
- **defects:**
  - :271-276: `context_matches` hashes the whole stock (`AuditStockStateIdentity2::build`) on every check.
  - With certification (:280-327) re-running every event authority, each audited op pays about 2× certification plus several whole-stock digests. Cost only (unverified magnitude).
- **unused capability / ideas:** none
- **claims:** none

### src/audit_classification_2.cpp | lines 91 | read 1-91
- **role:** validation (bindings)
- **used by:** src/compas_cgal/engagement_audit/classification.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_classification_2.h | lines 9 | read 1-9
- **role:** validation
- **used by:** stock_2_bindings.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_classification_core_2.cpp | lines 394 | read 1-394
- **role:** validation (exact motion classification)
- **used by:** engagement_audit/classification.py
- **state:** live
- **generation:** n/a, but it constrains generation.
- **defects:**
  - :174-201: requires an exactly orthonormal world-XY frame (`|xaxis|² == 1` exactly). Rotated COMPAS frames built from cos/sin are rejected (`AuditUnsupportedGeometryError`).
  - :268-271: mixed-Z XY ramps are unsupported, so a helical or ramped entry (the plan's D1) cannot enter this audit.
- **unused capability / ideas:** none
- **claims:** none

### src/audit_classification_core_2.h | lines 81 | read 1-81
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_decision_evidence_values_2.cpp | lines 203 | read 1-203
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_decision_strategy_2.cpp | lines 112 | read 1-112
- **role:** validation (decision limits)
- **used by:** audit units, engagement_audit/decision_limits.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_depletion_witness_2.cpp | lines 48 | read 1-48
- **role:** validation
- **used by:** stock_exact_depletion_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_depletion_witness_2.h | lines 57 | read 1-57
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_digest_2.h | lines 243 | read 1-243
- **role:** validation (typed SHA-256 domains)
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_exact_station_2.cpp | lines 442 | read 1-442
- **role:** kernel (THE exact engagement station classifier; `engagement_at` delegates here)
- **used by:** engagement_2.cpp:532 → `_stock_2.engagement_at` → engagement_toolpath.py:552 (gated generators), benchmarks/quality.py, survey.py
- **state:** live
- **generation:**
  - It zones the exact cutter circle through the stock arrangement (:70-125), assembles runs by exact endpoint equality (:127-177), and decides the cap per run on the exact rational chord surrogate (:410-423).
  - Quality: **adequate**. The predicate is exact, but it caps each run, not total contact.
- **defects:**
  - :410-423 (via exact_cap_certificate_2.h:84): the cap applies per maximal run. Two runs can sum far above the cap and still pass. Effect: two-sided contact is admitted by the generators.
  - :60-67: `found_overlap` drops rim pieces that coincide with a stock boundary circle. Where the overlapped boundary bounds material inside the cutter (an island of radius exactly r centred at the cutter, e.g. the uncut core of a loop with ρ = 2r), engagement reads 0 instead of 360°. Measure-zero in authored data. Effect: unverified in practice.
- **unused capability / ideas:**
  - `true_run_arcs()` and `boundary_intersections()` (:358-374) expose exact contact geometry. A placement step could solve for the exact advance that puts the run at the cap, rather than ladder-scanning stations.
- **claims:** :31-35 "trapezoidal point location supports the valid topology … without a floating-point decision": verified (uses `Stock2::point_location`).

### src/audit_identity_bindings_2.cpp | lines 260 | read 1-260
- **role:** validation (bindings)
- **used by:** engagement_audit/input.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_identity_bindings_2.h | lines 7 | read 1-7
- **role:** validation
- **used by:** stock_2_bindings.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_lineage_2.cpp | lines 90 | read 1-90
- **role:** validation
- **used by:** audit replay
- **state:** live (tests only)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_lineage_2.h | lines 43 | read 1-43
- **role:** validation
- **used by:** audit replay
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_motion_identity_2.cpp | lines 227 | read 1-227
- **role:** validation
- **used by:** audit classification
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_motion_identity_2.h | lines 157 | read 1-157
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_motion_result_2.cpp | lines 433 | read 1-433
- **role:** validation
- **used by:** audit replay
- **state:** live (tests only)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_motion_result_2.h | lines 199 | read 1-199
- **role:** validation
- **used by:** audit replay
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_policy_2.cpp | lines 177 | read 1-177
- **role:** validation (cap surrogate policy)
- **used by:** audit lane; engagement_2.cpp:504-507 (`cap_chord_ratio` delegates here)
- **state:** live
- **generation:** The surrogate `4 sin²(cap/2)` (:47-65) is the single boundary conversion every generator uses.
- **defects:**
  - :102-107: the depletion chord bound must be smaller than the tool radius. The chord bound sets sliver depth in the exact-depletion lane. Effect: audit-only.
- **unused capability / ideas:** none
- **claims:** :29-34 bit-identical surrogate check: verified.

### src/audit_policy_2.h | lines 92 | read 1-92
- **role:** validation
- **used by:** audit units, engagement_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_refinement_cell_2.cpp | lines 183 | read 1-183
- **role:** validation (partial-arc bisection cells)
- **used by:** audit_certification_2.cpp:189-269
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_2.cpp | lines 338 | read 1-338
- **role:** validation (transactional audit replay)
- **used by:** tests/engagement_audit/test_native_replay.py only (grep: no src/ or benchmarks/ consumer)
- **state:** parked (exact coverage replay set aside 2026-09-08)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_2.h | lines 198 | read 1-198
- **role:** validation
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_bindings_2.cpp | lines 289 | read 1-289
- **role:** validation (bindings)
- **used by:** tests only
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_bindings_2.h | lines 7 | read 1-7
- **role:** validation
- **used by:** stock_2_bindings.cpp
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_finalization_2.cpp | lines 37 | read 1-37
- **role:** validation
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_internal_2.h | lines 195 | read 1-195
- **role:** validation (test instrumentation hooks)
- **used by:** stock_2.cpp:281-282, 296-297, 671 (clone/swap/legacy-arc counters on the hot path)
- **state:** live
- **generation:** n/a
- **defects:** none found. The thread-local counters are cheap.
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_mutating_2.cpp | lines 226 | read 1-226
- **role:** validation
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:**
  - :20-51: each op runs a full exact certification, a trial depletion, a full certification replay and a depletion replay. This is the cost that set exact replay aside (unverified magnitude).
- **unused capability / ideas:** none
- **claims:** none

### src/audit_replay_non_engaging_2.cpp | lines 76 | read 1-76
- **role:** validation
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_request_identity_2.cpp | lines 67 | read 1-67
- **role:** validation
- **used by:** engagement_audit/input.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_request_identity_2.h | lines 43 | read 1-43
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_station_evidence_2.cpp | lines 282 | read 1-282
- **role:** validation
- **used by:** audit certification
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_station_parameters_2.cpp | lines 280 | read 1-280
- **role:** validation (exact station parameters on segment/circle/arc)
- **used by:** audit certification
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_stock_identity_2.cpp | lines 160 | read 1-160
- **role:** validation
- **used by:** engagement_audit/input.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_stock_identity_2.h | lines 51 | read 1-51
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_stock_state_identity_2.cpp | lines 50 | read 1-50
- **role:** validation (whole-stock digest)
- **used by:** audit certification, replay, and depletion witnesses
- **state:** live
- **generation:** n/a
- **defects:** :12-19: O(stock) boundary extraction per call, and it is called many times per audited op. Cost only (unverified magnitude).
- **unused capability / ideas:** none
- **claims:** none

### src/audit_stock_state_identity_2.h | lines 20 | read 1-20
- **role:** validation
- **used by:** audit units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_strategy_versions_2.cpp | lines 23 | read 1-23
- **role:** validation
- **used by:** audit identity
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_strategy_versions_2.h | lines 6 | read 1-6
- **role:** validation
- **used by:** audit identity
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/audit_unresolved_evidence_2.cpp | lines 583 | read 1-583
- **role:** validation
- **used by:** audit certification
- **state:** live
- **generation:** n/a
- **defects:** :203-214, :457-468, :504-515, :551-563: every unresolved path replays the event authority yet again. Cost only.
- **unused capability / ideas:** none
- **claims:** none

### src/boundary_circle_engagement_2.cpp | lines 115 | read 1-115
- **role:** engagement-control (Held standard predecessor model, exact CORE √-kernel)
- **used by:** no Python consumer found (grep `boundary_circle_engagement`: none). The header is reached through coverage_bindings_2.cpp.
- **state:** dead (no consumer found)
- **generation:**
  - It gives the exact engagement of a successor circle against a fully cleared predecessor outer disk, with the Held Eq. 4 continuity check (:84-89).
  - Quality: **adequate** as a placement predicate, but nothing calls it.
- **defects:**
  - :74: `CGAL::sqrt(distance_squared)` on CORE::Expr. An identically-zero decision downstream (for example `distance + outer == previous_outer` for concentric equal circles) is the documented CORE stall pattern (native_boundary_medial_2.cpp:11-23). Unverified here.
- **unused capability / ideas:** an exact engagement predicate between consecutive loops for Held-style spacing, without stock replay.
- **claims:** :97-100 "conservative full-slot bound … not an exact maximum": stated honestly.

### src/boundary_circle_engagement_2.h | lines 36 | read 1-36
- **role:** engagement-control
- **used by:** no consumer found
- **state:** dead
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** see .cpp
- **claims:** none

### src/boundary_contact_projection_2.h | lines 86 | read 1-86
- **role:** placement (project a contact onto a polygon draft)
- **used by:** no include found within this slice; consumer outside this slice not traced
- **state:** unverified consumer
- **generation:** exact selection of the closest boundary side; the parameter is returned as a double.
- **defects:**
  - :59-61: the parameter is clamped to [0,1] by exact comparison. This is fine.
  - :74: `maximum_distance` squared in exact arithmetic. This is fine.
  - none found
- **unused capability / ideas:** none
- **claims:** :39-40 "approximate view … not a native exact-contact identity": verified.

### src/canonical_encoding.cpp | lines 236 | read 1-236
- **role:** tool (canonical bytes)
- **used by:** audit, reachable, adaptive identity
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/canonical_encoding.h | lines 27 | read 1-27
- **role:** tool
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/canonical_rotation_2.h | lines 43 | read 1-43
- **role:** tool (Booth minimal rotation)
- **used by:** reachable_input_2.cpp, reachable_certificate_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas.h | lines 281 | read 1-281
- **role:** build (PCH and mesh converters, Epick)
- **used by:** types_std and other modules
- **state:** live
- **generation:** n/a
- **defects:** :229: `int j` against `size_t n` (sign compare). Effect: none on generation.
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/__init__.py | lines 1 | read 1-1
- **role:** validation
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/classification.py | lines 185 | read 1-185
- **role:** validation (Python ingress to native classification)
- **used by:** engagement_audit/input.py; tests; benchmarks/held_path_characterize.py (via input)
- **state:** live
- **generation:** n/a. It constrains generation: arcs and circles need exact axis frames (see audit_classification_core_2.cpp).
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/decision_limits.py | lines 89 | read 1-89
- **role:** validation
- **used by:** engagement_audit/input.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/digests.py | lines 9 | read 1-9
- **role:** validation
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/errors.py | lines 105 | read 1-105
- **role:** validation
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/identity.py | lines 105 | read 1-105
- **role:** validation (BuildIdentity)
- **used by:** engagement_audit/input.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/input.py | lines 252 | read 1-252
- **role:** validation (authoritative audit request)
- **used by:** benchmarks/held_path_characterize.py; tests/engagement_audit
- **state:** live (not gate-blocking)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/operation_identity.py | lines 406 | read 1-406
- **role:** validation (immutable op snapshots)
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/engagement_audit/records.py | lines 387 | read 1-387
- **role:** validation
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/compas_matrix.h | lines 18 | read 1-18
- **role:** build
- **used by:** all bindings
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/containment_2.cpp | lines 356 | read 1-356
- **role:** kernel (exact gouge / containment, CORE √-kernel)
- **used by:** src/compas_cgal/adaptive/containment.py (certified adaptive generator, not gated)
- **state:** live
- **generation:**
  - Exact subset tests:
    - a true capsule (√ offsets) inside the design (:225-256);
    - a true loop annulus inside the design (:258-292);
    - an entry disk (:294-320);
    - a loop inside an entry disk (:322-356).
  - Quality: **excellent** for gouge-freedom.
- **defects:**
  - :76-85: `disk_set` and `capsule_set` return `ReachSet` by value. NRVO-dependent, the same hazard exact_sweep_2.h:17-23 documents (latent, unverified).
- **unused capability / ideas:**
  - `evaluate_exact_full_circle_in_disk` is exactly the D1 certificate: the first loop lies inside the helix-bored disk.
  - The gated generators guarantee wall clearance only through `radial_clearance` offsets.
- **claims:** :240-243 "exact sweep subset is the authoritative gouge predicate": verified.

### src/containment_2.h | lines 66 | read 1-66
- **role:** kernel
- **used by:** containment bindings
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** see .cpp
- **claims:** none

### src/containment_2_bindings.cpp | lines 145 | read 1-145
- **role:** kernel (bindings, `_containment_2`)
- **used by:** adaptive/containment.py
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/coverage_2.cpp | lines 557 | read 1-557
- **role:** kernel (exact incremental coverage, CORE √-kernel)
- **used by:** src/compas_cgal/adaptive/coverage.py; benchmarks/coverage.py, survey.py, held_exact_motion_coverage.py
- **state:** live (adaptive / benchmarks), parked on the critical path
- **generation:**
  - Accumulated sweeps and residual are maintained with TRUE capsules and annuli (exact_sweep_2), so there are no under-cover slivers.
  - It supports a pre-cleared entry disk (:190-259).
  - Quality: **excellent** model, but slow (set aside 2026-09-08 per project memory).
- **defects:**
  - :451-458: every transition copies and re-overlays the full accumulated set and the full residual. O(total) per op. Cost only.
  - :441-497: traits-owner bookkeeping is resolved correctly through `ExactRegion2::build` candidates.
- **unused capability / ideas:**
  - The residual is the sliver-free answer to "does this op remove material".
  - `add_full_circle_sweep` takes a phase vector, so the loop radius equals the authored phase length exactly.
- **claims:** :49 "without inventing an entry cut": verified (`from_uncut`).

### src/coverage_2.h | lines 111 | read 1-111
- **role:** kernel
- **used by:** coverage bindings
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** see .cpp
- **claims:** none

### src/coverage_bindings_2.cpp | lines 506 | read 1-506
- **role:** kernel (bindings, `_coverage_2`)
- **used by:** adaptive/coverage.py, adaptive/reachable_domain.py, Held benchmarks
- **state:** live
- **generation:** It exposes the exact wall-offset cycle (`build_center_boundary_cycle` :245-260, `ccw_transition` :234-243), `NativeBoundary2.circle_on_piece` (:317-322) and `remaining_material` (:324-325).
- **defects:** :278: `UNUSED_VALIDATION_RADIUS_MM = 1.0` is a validation placeholder for a design-only region. Harmless.
- **unused capability / ideas:**
  - The wall-offset cycle and boundary-normal circles are unused by the gated generators.
  - B2a corner and wall run-ins and contour-aware placement should be built on them.
- **claims:** :135 "Approximate coordinates for reporting only": verified.

### src/cutter_centre_domain_2.cpp | lines 93 | read 1-93
- **role:** kernel (exact point-in-free-centre-domain)
- **used by:** benchmarks/survey.py, benchmarks/tools/held_toolpath_progress.py
- **state:** live (gate survey)
- **generation:** an exact test of whether a tool centre clears every wall by r. Quality: **excellent**, O(n) per query.
- **defects:** :37-41: a strict `<` means a centre exactly r from a wall is admitted. This is correct for a closed domain.
- **unused capability / ideas:** a generator could validate every station centre exactly rather than by `radial_clearance` margin.
- **claims:** none

### src/cutter_centre_domain_2.h | lines 21 | read 1-21
- **role:** kernel
- **used by:** coverage bindings
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/engagement_2.cpp | lines 592 | read 1-592
- **role:** engagement-control (the `engagement_at` and `certify_segment_tea` entry points)
- **used by:** `_stock_2.engagement_at` → engagement_toolpath.py:552,596; engagement_radial / rho / spiral generators; benchmarks/quality.py, survey.py, exceedance.py. `certify_segment_tea` → engagement.py:369 (gate audit), benchmarks/spec.py.
- **state:** live
- **generation:**
  - Station verdicts are exact. Segment certification uses guarded bisection.
  - Quality: **weak** as a generation tool. It cannot certify within about 5–7° of the cap, and it costs 511 stations for a 2 mm segment (P4).
- **defects:**
  - :352-412: lemma (a) "drift ≤ arc of chord d" is contradicted. Measured jump 23.8° over 0.001 r against 2·GROWTH = 10.5° (P5). The certificate's soundness against concave conforming boundaries (radius in (r, ~4r/3)) is unproven.
  - :434: STATION_FLOOR_FRACTION 1e-3 leaves a residual guard of 0.090–0.128 rad (arithmetic), so near-cap motions always come out uncertified.
  - :29-36, :86-345: dead code (`Arc`, `pessimistic_runs`, `finish_engagement`, `engaged_arcs_zone`, `EngagementVisitor`, `PL`). No caller: grep shows only self-references.
  - :86: orphaned doc comment (first line lost).
  - :110-114 (header) and :521-558: the cap applies per run, not to total contact.
- **unused capability / ideas:** `gap_close_ratio` pessimism (:521-538) is implemented but the generators pass 0.0 (engagement_toolpath.py:552).
- **claims:**
  - Header :116-122 "O(cutter-crossings), not O(stock)": partially contradicted. Each new stock epoch rebuilds the trapezoid RIC (stock_2.cpp:265-272) at O(n log n) on first query.
  - :366-371 "GROWTH monotone": verified.
  - :386-399 QED: contradicted at its premise (a).

### src/engagement_2.h | lines 199 | read 1-199
- **role:** engagement-control
- **used by:** engagement_2.cpp, bindings
- **state:** live
- **generation:** see .cpp
- **defects:** :110-114: `EngagementSample` documents a per-run decision with no total-contact verdict.
- **unused capability / ideas:** none
- **claims:** :167-181 "two stations both under the guarded cap certify the whole span": unverified-sound (see P5).

### src/engagement_2_bindings.cpp | lines 78 | read 1-78
- **role:** engagement-control (bindings)
- **used by:** `_stock_2`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/engagement_2_bindings.h | lines 5 | read 1-5
- **role:** build
- **used by:** stock_2_bindings.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_build_audit_2.h | lines 41 | read 1-41
- **role:** validation (counters)
- **used by:** reachable, coverage
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_cap_certificate_2.h | lines 154 | read 1-154
- **role:** kernel (the single cap predicate)
- **used by:** audit_exact_station_2.cpp (therefore every generator), engagement_2.cpp
- **state:** live
- **generation:** an exact orientation sign plus an exact chord comparison over one-root points via a mixed-radical sign. Quality: **excellent** as a predicate.
- **defects:** :84-140: per-run semantics only (see finding 5).
- **unused capability / ideas:** none
- **claims:** :14-21 "the loud path is kept" (shared-root precondition throws): verified (:56-68).

### src/exact_circle_chart_2.cpp | lines 71 | read 1-71
- **role:** motion-geometry (exact rational points on a circle)
- **used by:** exact_depletion_2.cpp, audit units
- **state:** live
- **generation:** exact rational parametrisation of a circle of authored phase vector. Quality: **excellent**.
- **defects:** none found
- **unused capability / ideas:** Generators could author loop entries and bridge endpoints as exact chart points, giving exact incidence. That cures finding 2 if the op schema carries exact rationals.
- **claims:** none

### src/exact_circle_chart_2.h | lines 41 | read 1-41
- **role:** motion-geometry
- **used by:** as above
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_circle_chart_atlas_2.h | lines 45 | read 1-45
- **role:** motion-geometry
- **used by:** exact_circle_chart_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_depletion_2.cpp | lines 829 | read 1-829
- **role:** kernel (exact under-cover disk-chain constructions with traces)
- **used by:** `Stock2::subtract_exact_*` → src/compas_cgal/adaptive/stock_area.py; audit depletion
- **state:** live (adaptive lane)
- **generation:**
  - Exact chain centres on segments and circles (Pythagorean chart), with an exact chord bound.
  - Quality: **adequate**. The model is still an under-covering chain, so it leaves slivers of depth up to the chord sagitta.
- **defects:**
  - :187-201, :597-607: each doubling re-walks every chord, O(n log n) exact work. Cost only.
  - :147-168: the arc refinement throws once `intervals >= center_count_limit`, even when `intervals+1` would fit. Off by one, conservative.
- **unused capability / ideas:** none for the gated generators.
- **claims:** none

### src/exact_depletion_2.h | lines 156 | read 1-156
- **role:** kernel
- **used by:** stock_2.h
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_disk_region_2.cpp | lines 31 | read 1-31
- **role:** kernel
- **used by:** stock_2.cpp, stock_exact_depletion_2.cpp
- **state:** live
- **generation:** union of exact disks
- **defects:** :24-31: the by-value `Gps` return is NRVO-dependent (latent, as the header :18-23 itself warns).
- **unused capability / ideas:** none
- **claims:** none

### src/exact_disk_region_2.h | lines 30 | read 1-30
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** :18-23 lifetime warning: verified as a description. The by-value function is still used at stock_exact_depletion_2.cpp:308, 335, 357, 379.

### src/exact_motion_2.h | lines 20 | read 1-20
- **role:** kernel (Epeck types)
- **used by:** all Epeck units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_region_2.cpp | lines 276 | read 1-276
- **role:** kernel (immutable exact regions over the CORE √-kernel)
- **used by:** reachable, coverage, containment
- **state:** live
- **generation:** n/a (storage)
- **defects:** none found. Traits ownership is resolved (:158-190).
- **unused capability / ideas:** none
- **claims:** :251-256 audit accessor: verified.

### src/exact_region_2.h | lines 123 | read 1-123
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_stock_region_2.cpp | lines 133 | read 1-133
- **role:** kernel (lift an Epeck stock into the √-kernel)
- **used by:** no consumer found within this slice (callers would be in continuous_tea_2 or adaptive)
- **state:** unverified consumer
- **generation:** A bridge from the generator's stock to exact coverage queries. It would allow "stock ∩ true sweep" decisions without slivers on the sweep side.
- **defects:** :119-133: `reach_join_parts` returns by value (NRVO-dependent).
- **unused capability / ideas:** use it for an exact removal test against TRUE sweeps.
- **claims:** none

### src/exact_stock_region_2.h | lines 13 | read 1-13
- **role:** kernel
- **used by:** see .cpp
- **state:** unverified consumer
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact_sweep_2.cpp | lines 269 | read 1-269
- **role:** kernel (TRUE swept regions: capsule with √ offset, arc sector, full-circle annulus)
- **used by:** containment, coverage, reachable, remaining_material
- **state:** live
- **generation:** exact sweeps with no under-cover. Quality: **excellent**.
- **defects:** :191-204: the EQUAL and SMALLER branches are identical. Harmless.
- **unused capability / ideas:** This is the sliver-free sweep the gate replay should use where cost allows. The authored loop radius goes in as the phase length.
- **claims:** :17-23 lifetime note: verified.

### src/exact_sweep_2.h | lines 39 | read 1-39
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/held_contour_clearance_2.cpp | lines 88 | read 1-88
- **role:** placement (Held exact whole-orbit engagement bound)
- **used by:** benchmarks/held_coverage_thinning.py, held_contour_bound_spacing.py; tests
- **state:** live (Held reproduction only)
- **generation:** an exact upper bound on the uncut leading-semicircle engagement from the clearance to the filled-disk contour. Quality: **excellent** for placement.
- **defects:** none found. Cross-root `CoordNT` compares are exact (`Sqrt_extension` ACDE = true; Sqrt_extension_type.h:528-600).
- **unused capability / ideas:** A Held-style spacing predicate the product generators could use in place of probe sampling.
- **claims:** :81-84: verified by reading.

### src/held_disk_contour_2.cpp | lines 113 | read 1-113
- **role:** placement (Held filled-disk contour, exposed-arc contact)
- **used by:** Held benchmarks
- **state:** live (Held only)
- **generation:** the exact contact b on the latest outer circle, moved CW to the exposed contour. Quality: **excellent**.
- **defects:** :64,72: `contour_.join` per append overlays the whole union. O(n) per circle. Cost only.
- **unused capability / ideas:** same as above.
- **claims:** none

### src/held_disk_contour_2.h | lines 51 | read 1-51
- **role:** placement
- **used by:** as above
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** :26-27 "placement model, distinct from Stock2's physical annular sweeps": verified.

### src/held_disk_contour_2_bindings.cpp | lines 66 | read 1-66
- **role:** placement (bindings)
- **used by:** Held benchmarks
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/native_boundary_curve_2.cpp | lines 117 | read 1-117
- **role:** guide (exact line+arc boundary import)
- **used by:** benchmarks/held_native_curve_import*.py
- **state:** live (Held only)
- **generation:** exact curved pocket boundaries. Gated generators accept only polygons.
- **defects:** none found
- **unused capability / ideas:** curved-pocket support for the product generators.
- **claims:** :45-46 "measured 2026-09-08, 25 ms per competitor line": unverified (not re-measured).

### src/native_boundary_curve_2.h | lines 60 | read 1-60
- **role:** guide
- **used by:** as above
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/native_boundary_medial_2.cpp | lines 274 | read 1-274
- **role:** placement (exact boundary-normal medial disk, Held contour-aware circles)
- **used by:** benchmarks/held_native_boundary_draft.py (via `circle_on_piece`)
- **state:** live (Held only)
- **generation:** For a boundary point p, the first maximal disk tangent at p, and the largest tool circle touching the wall at p (:260-273). Quality: **excellent** for corner and wall-aware placement.
- **defects:** :209-213, :250-253: throws when the one-sided normal is not admissible. Fail-loud, correct.
- **unused capability / ideas:** a corner and wall placement primitive the plan's D2 and B2a ignore.
- **claims:** :11-17 CORE timings: unverified (not re-measured).

### src/reachable_arrangement_2.cpp | lines 1013 | read 1-1013
- **role:** kernel (exact free-centre domain by arrangement parity)
- **used by:** reachable_domain_2.cpp, coverage bindings (`build_center_boundary_cycle`)
- **state:** live
- **generation:** the exact Minkowski erosion of the pocket by r, with source provenance, plus CCW boundary transitions (:773-890). Quality: **excellent**.
- **defects:** :997-1006: rejects a centre domain with more than one component. A pocket with a neck narrower than 2r cannot be measured or covered through this kernel.
- **unused capability / ideas:** the exact wall-offset path for wall run-ins, finishing, and rest-material passes.
- **claims:** none

### src/reachable_arrangement_2.h | lines 92 | read 1-92
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_boundary_sampling_2.cpp | lines 122 | read 1-122
- **role:** kernel (exact sample points on boundary pieces)
- **used by:** native_boundary_medial_2.cpp, coverage bindings
- **state:** live
- **generation:** n/a
- **defects:**
  - :15: `SLIVER_SWEEP_RAD = 1e-9` is a construction-regime switch (a double angle threshold). Membership is still checked exactly (:113-120). Acceptable under the construction-guard rule.
- **unused capability / ideas:** none
- **claims:** :55-58 timings: unverified.

### src/reachable_boundary_sampling_2.h | lines 10 | read 1-10
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_certificate_2.cpp | lines 925 | read 1-925
- **role:** validation (reachable-domain certificate)
- **used by:** reachable_domain_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_certificate_2.h | lines 27 | read 1-27
- **role:** validation
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_certificate_encoding_2.cpp | lines 186 | read 1-186
- **role:** validation
- **used by:** a core library outside this slice (CMakeLists.txt, segment_site lib)
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_certificate_encoding_2.h | lines 40 | read 1-40
- **role:** validation
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_domain_2.cpp | lines 221 | read 1-221
- **role:** kernel (exact reachable material and unreachable residual)
- **used by:** src/compas_cgal/adaptive/reachable_domain.py; benchmarks/qualityfigures.py
- **state:** live (adaptive / figures; not the gate generators)
- **generation:** exact `reachable_material` = centre domain ⊕ r (:60-80, :126-144) and `unreachable_residual` = design \ material (:142-144). Quality: **excellent**.
- **defects:** none found
- **unused capability / ideas:**
  - The exact corner material the tool cannot reach at all. D2(a) should be stated against it.
  - The gate's "uncut fraction" should be normalised by it (whether the gate already does was not traced here).
- **claims:** none

### src/reachable_domain_2.h | lines 46 | read 1-46
- **role:** kernel
- **used by:** coverage bindings
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_errors_2.h | lines 28 | read 1-28
- **role:** kernel
- **used by:** reachable units
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_input_2.cpp | lines 347 | read 1-347
- **role:** kernel (canonical exact input rings)
- **used by:** reachable, cutter-centre, material predicate
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_input_2.h | lines 61 | read 1-61
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** :12-21 redefines `compas::RowMatrixXd` (duplicates compas_matrix.h). Harmless.
- **unused capability / ideas:** none
- **claims:** none

### src/reachable_material_predicate_2.cpp | lines 348 | read 1-348
- **role:** kernel (point-in-reachable-material)
- **used by:** benchmarks/coverage.py (gate coverage)
- **state:** live
- **generation:** n/a (measurement)
- **defects:** :281-284: same single-component restriction.
- **unused capability / ideas:** none
- **claims:** :199-205 lifetime rationale: verified.

### src/reachable_material_predicate_2.h | lines 33 | read 1-33
- **role:** kernel
- **used by:** as above
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/remaining_material_2.cpp | lines 85 | read 1-85
- **role:** kernel (exact full-motion remaining material, √-kernel)
- **used by:** benchmarks/held_exact_motion_coverage.py, toolpath_coverage.py, acceptance/test_held_coverage_replay_budget.py
- **state:** parked (set aside 2026-09-08)
- **generation:** sliver-free coverage keyed on the authored loop radius (:57-59 uses `(ρ,0)`). Quality: **excellent**, slow.
- **defects:** :55-73: one full overlay per motion. O(motions × stock). Cost only.
- **unused capability / ideas:** the reference oracle for the sliver question in finding 1.
- **claims:** none

### src/remaining_material_2.h | lines 14 | read 1-14
- **role:** kernel
- **used by:** coverage bindings
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** :8 "Final coverage only; no chronological engagement": verified.

### src/stock_2.cpp | lines 1089 | read 1-1089
- **role:** kernel (exact depleting stock: the generators' and the gate's material model)
- **used by:** every generator via src/compas_cgal/stock.py; gate replay (benchmarks/depletion.py:109-112, survey.py:363-382)
- **state:** live
- **generation:**
  - Exact `Gps` stock, exact annulus loops, disk-chain and quad capsules, exact-depletion lanes, and clone/swap with resolved traits ownership.
  - Quality: **adequate**. Exact booleans, but under-covering line sweeps and surrogate-radius loop sweeps feed every "removes material" decision (P1–P3).
- **defects:**
  - :516-540: chain capsule. A 12 mm run adds about 1,200 vertices (P4). It leaves sliver islands between adjacent runs (P3: 1,606 faces) and makes every later boolean about 100× costlier.
  - :673-697: the loop radius is `hypot` of a sampled start point, not the authored radius. Same-loop re-sweeps leave ulp slivers (P1: 30/39).
  - :580-666: the quad still under-covers by up to 5e-5 r. Slivers remain (P1: True).
  - :331-335: `exactly_equals` costs 2 full differences. A removal test costs clone + 3 overlays, where the alternative is one interior-intersection test.
  - :81-86, :114-122: by-value `Gps` helpers are NRVO-dependent (latent).
- **unused capability / ideas:**
  - `subtract_circle_sweep(cx,cy,ρ,r)` (:408-421) takes the authored radius exactly. The replay should use it for `Circle`.
  - `subtract_annulus_exact` takes exact radii.
  - `can_remove_circle` (stock_circle_redundancy_2.cpp).
- **claims:**
  - :677-683 "the annulus is the swept region exactly, nothing left over": contradicted for re-sweeps of the same nominal circle (P1). The surrogate caveat is stated at :689-693.
  - :350-352 slack dominated by `radial_clearance`: verified (arithmetic).
  - :1022-1025 comment says `unique_ptr` but `set_` is a `shared_ptr`. Cosmetic.

### src/stock_2.h | lines 358 | read 1-358
- **role:** kernel
- **used by:** all stock consumers
- **state:** live
- **generation:** see .cpp
- **defects:** none found beyond .cpp
- **unused capability / ideas:** none
- **claims:** :101-118, :286-303 lifetime model: verified by reading `resolve_traits_owner` and `replace_set` (stock_2.cpp:188-324).

### src/stock_2_bindings.cpp | lines 537 | read 1-537
- **role:** kernel (bindings, `_stock_2`)
- **used by:** stock.py and all generators
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `subtract_circle_sweep` and `can_remove_circle` are bound but unused by the gated generators.
- **claims:** none

### src/stock_circle_redundancy_2.cpp | lines 53 | read 1-53
- **role:** ordering (exact loop-deletion test)
- **used by:** benchmarks/held_coverage_thinning.py; tests/test_circle_redundancy.py
- **state:** live (Held only)
- **generation:** decides whether removing loop k loses material, given loops k±1 and quad connectors. The test is conservative: empty means removal is safe. Quality: **excellent**.
- **defects:** none found
- **unused capability / ideas:** a post-pass that deletes redundant loops from any generator. It is a sound alternative to the plan's per-op "removes material" test.
- **claims:** :50-52: verified.

### src/stock_exact_depletion_2.cpp | lines 544 | read 1-544
- **role:** validation (audit trial depletion and witnesses)
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:** :308,335,357,379: `Gps removal = build_exact_disk_union_region_2(…)` by value (NRVO-dependent, latent).
- **unused capability / ideas:** none
- **claims:** none

### src/stock_exact_depletion_2.h | lines 49 | read 1-49
- **role:** validation
- **used by:** audit replay
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/stock_local_2.cpp | lines 324 | read 1-324
- **role:** kernel (local depletion by arrangement surgery)
- **used by:** `Stock2::subtract_arc_sweep_local` → engagement_toolpath.py:795, rho/radial/spiral generators
- **state:** live
- **generation:** Removes loop annuli without a global overlay. Quality: **adequate** (sound argument :25-58).
- **defects:**
  - :299-301: inserting a boundary circle that OVERLAPS an existing stock edge (same circle as an earlier sweep) is not covered by the I1–I4 argument. The orientation of the overlapped edge may not satisfy the `Gps` "contained side left" invariant. Unverified; `representation_is_valid` exists to gate it.
  - The gate replays with the global path, so generator stock and replay stock can differ structurally.
- **unused capability / ideas:** none
- **claims:** :60-64 "No epsilon … anywhere": verified.

### src/stock_local_2.h | lines 65 | read 1-65
- **role:** kernel
- **used by:** stock_2.cpp
- **state:** live
- **generation:** see .cpp
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/types_std.cpp | lines 11 | read 1-11
- **role:** build
- **used by:** `_types_std`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/types_std.h | lines 4 | read 1-4
- **role:** build
- **used by:** types_std.cpp
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none
