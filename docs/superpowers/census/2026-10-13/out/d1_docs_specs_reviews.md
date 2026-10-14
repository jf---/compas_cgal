## Slice summary

Slice: 78 assigned files (all read in full; `docs/public/favicon.ico` is a binary PNG, inspected with `file(1)`) plus 14 off-branch documents. The four `figures/fig6_*.{svg,pdf}` added by `f369b3df` are images, not documents, and were not read. 202 findings in the register at the end.

### Top findings, ranked by impact on toolpath quality (all open, all generation-relevant)

1. **No entry strategy exists that the gate or the auditor can accept.** Every cap exceedance and the 325-338° engagement step are chain-entry loops (machining_quality.md:160-168). A rim plunge followed by a loop gives 360°. A centre plunge plus concentric ramp has a proved floor of 240° and a measured worst step of 124.64° (loop_radius_degeneracy.md:104-134). A helical or ramped entry cannot be replayed (loop_radius_degeneracy.md:136-149), and the truthful auditor rejects ramps outright (engagement_audit.md:40,127-128). The certified spec says a plunge of one tool radius cannot launch a capped lateral move at all (adaptive design spec :59-61). A **precleared bore** of radius > r is already built as typed input (`adaptive/entry.py`, segment_site_mat.md:1025-1055) and needs no Z+XY vocabulary. The draft plan does not consider it. [F-34, 64, 65, 69, 111, 119, 133, 150, 171]
2. **Engagement is limited, not regulated.** Measured TEA saturates near 141.8° for any cap from 20° to 140° (benchmarks.md:120-133). Loops come out at 171.8° whether the cap is 40° or 120°: "maximising the loop rather than regulating it" (review.md:470-476). `radius_regulated_toolpath` still imports the *bisecting* advance search (`engagement_radial_toolpath.py:138`), although engagement is non-monotone (measurement_claims.md:18-33). [F-13, 35, 78, 131, 172, 184, 192]
3. **Chain order and linking.** Chains are ordered by greedy nearest neighbour from the longest chain, are never reversed, and each gets its own plunge (`src/toolpath.cpp:961-992`). A stay-down chain-ordering generator was built (`engagement_ordered_toolpath.py`, `7af7a24c`) that took material entries from 5 to 1. It admitted links with the exact predicate, tried both chain orientations and cap-checked bridges. It was deleted as "unreachable" in `9fa6eea8` (machining_metric_validity.md:136-141). The certified lane's DFS route was measured not to be a physical walk: direct links between routes gouge (continuous_engagement.md:1292-1327). An exact retrace through accepted lineage was built and then superseded (red_manifest). [F-43, 67, 123, 127, 167]
4. **Convex corners and narrow regions have no motion.** A station with `ρ ≤ r` cannot carry a trochoid. On the middle curve that already happens at `d ≤ 3r` (loop_radius_degeneracy.md:39-42). `declined_regions` is reported and consumed by nothing (ETH plan Task 31; verified by grep). A `W ≤ 4r` channel comes out as 68 plunges (machining_quality.md:305-313). [F-62, 73, 101, 168]
5. **Residue is not a corner problem.** The farthest uncut sample is 3.28 units from any vertex. Residue is 5× the unregulated generator's, and "coverage is not in the accept/reject rule" (machining_quality.md:150-158; INTEGRATECODEX.md:157-161). [F-71, 179]
6. **Bridges are not regulated by either gated generator.** The worst bridge measures 141.8° (radius_regulated_toolpath.md:28-32,198-200). Walking a corner chain from the spine out to its tip makes every bridge advance into virgin stock, which shows up as slotting (deleted ordered-generator docstring). [F-100]
7. **G1 breaks come from the bridge tilt.** The bridge leaves the guide tangent by `atan(dρ/ds)`, 35.26° on every corner chain. The external common tangent removes it (loop_radius_degeneracy.md:151-190). The tangent test exempts bridges as an "inherent ~90°" junction (toolpath_tangent_continuity.md:33-51,167-171) and checks only 7.7% of transitions. [F-37, 66, 146, 132]
8. **Air cutting.** Full 360° loops re-cut about 90% of their swept area (recut fraction 0.8974). The partial-arc trochoid was abandoned (full-circle design :94-103). There are 10 zero-engagement redundant loops. [F-72, 85]
9. **No lower engagement bound and no constant-engagement target.** Nothing addresses chip thinning or rubbing, and there is no feed post-processor. The design's `tea_target ± tea_band` (≥95% of length in band) was never built (adaptive-clearing design :34-35). [F-74, 89, 91]
10. **Smaller items.** `mat_scale > 1` still emits gouging circles (`src/toolpath.cpp:607-613` does not validate it). Hole validation into the straight skeleton is incomplete, which is undefined behaviour with islands. Engagement differs under translation at exact rim/boundary tangency. [F-156, 165, 68]

### Generation capabilities in this slice that the gated generators do not use

- From `rho_regulated` and spiral entry: the `ρ > r` floor, external-tangent bridges, ladder scan of the advance, and centre plunge plus concentric ramp.
- Typed entry and region definitions: a precleared bore as entry (`adaptive/entry.py`), the exact reachable material `M_r` and the unreachable corner set `D \ M_r` (`ReachableDomain2`, any polygon with holes, single `C_r` component), `CutterCentreDomain2`, and `ReachableMaterialPredicate2`.
- Cheap exact tests: an exact station-refutation ladder (9 stations, 5.39× faster rejection), and `Stock2.can_remove_circle`, a local exact redundancy test that is cheaper than a full subtract-and-compare.
- Certification and motion types:
  - the swept-prefix theorem, which certifies only a cap of exactly π and is "research-only" until auditor P2;
  - `certify_segment_tea` / `certify_arc_tea` (their growth bound is unsound for concave voids where ρ → r⁺);
  - zero-guide advancing segments, a distinct operation type that "owns its advance".
- The certified lane's candidate order: furthest progress, then largest radius (MATHSM lattice).
- The Held lane:
  - contour-aware qm-diameter circles (`BoundaryNormalCircle2`, `NativeBoundary2.circle_on_piece`), exposed-contact contour (`HeldDiskContour2`);
  - **boundary-order spacing**: 2.4× fewer circles and 58.7% shorter on Figure 8 upper;
  - convex-corner radius limits plus **concave-corner shared-contact fans**, Monstera "corner approaches", qualified contacts.
- Local zone engagement query, local annulus depletion, and the 6-curve capsule quad.
- `offset_polygon` and the weighted straight skeleton.
- The MATHSM constant-spacing benchmark sweep.
- Quality reducers with operation attribution (`quality_observations.py`).
- The **deleted** chain-ordering stay-down generator, recoverable with `git show 9fa6eea8^:src/compas_cgal/engagement_ordered_toolpath.py`.

### Where this slice shows the draft plan is wrong or incomplete

1. **It re-derives a generator that was already built and deleted.** `engagement_ordered_toolpath.py` (7af7a24c, 685 lines, removed 9fa6eea8) already orders and orients chains so each later chain starts in cleared stock. It uses the same exact predicate to choose links and orientation, cap-checks bridges, and reports `forced_entries`. Task 2-3 should start from it, or say why not. Its docstring also warns that walking a corner chain from the spine outward makes the **bridges advance into virgin material and shows up as slotting**. The plan's Design does not cap-check bridges, so B3 (slotting) may get worse under the tree walk.
2. **Bridges and links are not cap-checked.** The gated generators' bridges reach 141.8°. Productive-only emission removes loops that cut nothing, but it does not stop a bridge that cuts too much.
3. **D1 entry:**
   - The plan excludes the HELIX from the planar cap, but loop_radius_degeneracy.md:98-102 forbids excluding the entry from the *step* criterion. The step from the helix to the first loop must stay measured.
   - The truthful auditor rejects ramps.
   - A **precleared bore** (built) and the spiral entry (built; worst step 124.64° vs 120°) are cheaper options. The plan mentions only the 240° floor.
4. **D2 corners:**
   - Held's paper regularizes convex corners with radius slightly > r. The machinable target is the paper's precedent for D2(a).
   - `D \ M_r` gives an exact rest-material definition.
   - ETH Task 31 already specified the contract: "either the roughing programme clears these corners with a non-trochoidal motion, or it declares them for a second operation with a smaller tool". Report uncut material that was **declined and named** separately from material **silently missed**.
   - held_motion_coverage.md:251-258 proves that a straight corner approach cannot be certified in the predecessor-only model. B2a's straight wall run-ins need a stock-aware check. The swept-prefix theorem certifies only π, not the 60° slot threshold.
5. **Coverage is unaddressed.** The gated generators' residue is not corner-confined. The tree walk reuses rho's station rules, which leave 1-3.5% uncut. Nothing in Stage A adds coverage pressure.
6. **The cost of productive-only emission is unmeasured.** A full `_subtract_operation` plus `Stock.exactly_equals` per candidate is the O(n²) exact depletion that CHANGELOG.md:13 deferred. `can_remove_circle` and the zone query are local alternatives.
7. **The plan is right not to switch the guide.**
   - The exact segment-site MAT builds only one exact six-edge L fixture and one rectangle contract (`require_l_shape_input`, segment_site_mat.md:3733-3734,4358-4359).
   - The middle curve was measured *not* gentler than the straight skeleton (CoV 0.5719 vs 0.5752) and has more degenerate loops (24.5% vs 13.0%) (ETH plan :3215-3225).
   - The certified `adaptive/` generator has never finished a pocket, certifies advancing cuts only at π, and costs about 15 s per segment.
   - Its pieces can still be borrowed: the refutation ladder, the "advance vs hold" operation distinction needed for D3, and the precleared entry.
8. **Retract-and-plunge repositioning** was explicitly rejected by the retrace design as adding clearance, collision and kinematic contracts (retrace spec :78-82). The plan should state the air-time and retract-count cost it accepts. The gate needs the advance-versus-hold operation distinction before D3.
9. **Not covered at all:**
   - the cap-40 saturation root cause (station set vs regulator, ETH Task 37 Step 1);
   - narrow channels with `W ≤ 4r`;
   - lead-in/lead-out (none emitted);
   - arc links (backlog C5);
   - the partial-arc trochoid and air time;
   - a lower engagement bound;
   - `mat_scale` and hole validation in the Epick guide (Review Focus 1 islands meet ETH I9 undefined behaviour).

### Open generation-relevant findings ranked by impact, and draft-plan coverage

| rank | finding (register ids) | covered by draft plan? |
| --- | --- | --- |
| 1 | entry: 360° loop / 329° step; no replayable entry (34, 64, 65, 69, 75, 111, 119, 133, 150, 171) | D1, as a ruling only; precleared bore and step criterion not covered |
| 2 | chain order, per-chain plunge, redundant loops (43, 72, 167, 182) | Stage A |
| 3 | stay-down links / transit (67, 123, 127, 126) | D3, Stage C; deleted generator and retrace not cited |
| 4 | convex corners below the ρ>r floor; declined regions (45, 62, 101, 168, 52, 117) | B2/D2; Task 31 contract and M_r not cited |
| 5 | engagement not regulated: saturation, loop maximised, bisection (13, 35, 78, 131, 172, 184, 192) | Stage C "re-measure" only — **not covered** |
| 6 | bridges unregulated (100) and slotting from tip-ward walks | **not covered** (B3 attributes slotting after the fact) |
| 7 | coverage / non-corner residue (71, 179, 38) | **not covered** |
| 8 | G1: bridge tilt, exempted tangent test (37, 66, 146, 132) | via rho rules (66); test exemption **not covered** |
| 9 | narrow channels W ≤ 4r (73) | **not covered** (dumbbell checks only a second entry) |
| 10 | air cutting: full circles, partial-arc trochoid (85) | **not covered** |
| 11 | no constant-engagement band / lower bound / feeds (74, 89, 91) | **not covered** |
| 12 | contour-aware placement, boundary-order spacing, Held corner fans (53, 54, 58, 107, 112) | Stage C rejoin only — **not covered** for the product generator |
| 13 | G2 / arc links / leads (110, 145, 151) | **not covered** |
| 14 | forced loops/advances erased at the public result (130, 185) | partly (`virgin_plunges` counter) |
| 15 | guide hygiene: mat_scale, hole validation, Epick decisions (156, 165, 163, 177) | **not covered** (islands in Review Focus 1 only) |
| 16 | translation variance at exact tangency (68) | **not covered** (parking lot C1) |

# Records

### .superpowers/sdd/2026-08-23-auditor-convergence-p1-truthful-audit/task-3-report.md | lines 136 | read 1-136
- **role:** doc (SDD task report, validation)
- **used by:** no consumer found (agent scaffolding; docs commits 2875cf8d/aa1f28c7 dropped other such scaffolding, this one survived)
- **state:** parked (auditor-convergence P0-P4 is in the draft plan's parking lot)
- **generation:** n/a. Records the exact rational-chart arc surrogate for depletion (`construct_exact_arc_depletion`, `Stock2::subtract_exact_arc`). No generation behaviour.
- **defects:** none found for generation. :128-132 the surrogate proves only removed-set ⊆ swept set; "No quantitative retained-sliver bound is claimed" — so an audit replay with this depletion under-removes, which can over-report uncut/engagement relative to the real sweep.
- **unused capability / ideas:** an exact, content-addressed arc depletion (`subtract_exact_arc`) that the generators do not use (they deplete with `_subtract_operation`).
- **claims:** "Implementation ready ... No commit created" (:5-6) — contradicted by progress.md:20 (`27f97c7`), which is an ancestor of the live branch (verified `git merge-base --is-ancestor 27f97c7 HEAD`). `src/audit_arc_motion_2.cpp` exists on live (verified `ls src`).

### .superpowers/sdd/2026-08-23-auditor-convergence/progress.md | lines 215 | read 1-215
- **role:** doc (SDD progress ledger, validation)
- **used by:** no consumer found
- **state:** parked (auditor-convergence P0-P4 parked by draft plan "Parking lot")
- **generation:** n/a. Every entry is audit/replay identity work. :41-42 records "benchmark instrument — 259 passed; 6 product-gate failures" at the time, i.e. the product gate was red then and the auditor programme did not touch generation.
- **defects:** none found for generation.
- **unused capability / ideas:** native typed replay (six motion routes, one clone + one swap per mutating call, :180-191) — an exact replay transaction a generator could deplete through so its stock is the auditor's stock. The draft plan instead depletes with `_subtract_operation` (the gate replay's primitive), which is the right one for the gate; noted only.
- **claims:** :7 "Current HEAD 7cda8f6" — verified ancestor of live. :9 "Certifier source 73d5372" — that sha is NOT an ancestor of live (verified `git merge-base --is-ancestor 73d53729 HEAD` → NOT); the certifier content may have arrived under other shas, unverified. :153-156 machined/spiral arc-adapter full-circle closure exceeded 120 s and "no machined/spiral arc-adapter coverage is claimed" — unverified on live.

### .superpowers/sdd/2026-08-31-held-pfeiffer-reference-pockets/task-1-report.md | lines 187 | read 1-187
- **role:** doc (SDD task report; Held input reconstruction)
- **used by:** no consumer found (records `benchmarks/held_*` input work)
- **state:** live (Held reference pockets are on live: commits da1512c4, bcfc1a1b etc. in `git log`)
- **generation:** n/a for generation; it reconstructs Held's four pocket boundaries (lines+arcs) with certified deviation bounds (:155-161). This is the input side of the only pipeline that reproduces a published adaptive path.
- **defects:** :19-21 reconstruction and proof recursion have a depth ceiling but no shared node budget — open. :123-124 / :147-150 biarc child arcs can never merge across spans → more primitives than needed (cost, not path quality).
- **unused capability / ideas:** certified curve→arc reconstruction (G1 biarcs, polar-monotonicity certificate) could feed arc-bounded product pockets; the gated generators accept only polygons.
- **claims:** live census Figure 5 `5+44=49` etc. (:158-161) — verified as later than `docs/benchmarks/held_reference_qualification.md:7-10`, which still shows the earlier `31/70/58/317` counts (contradiction between the two docs; the qualification page is stale). :128 "independent certificate re-review remains pending" — unknown whether done.

### CHANGELOG.md | lines 264 | read 1-264
- **role:** doc
- **used by:** `docs/changelog.md` (snippet include)
- **state:** live, partly stale
- **generation:** records the trochoid generator's certified gouge-free station radii (:14), stepover bounding radial engagement per cycle (:15), holes (:16), climb/conventional (:17), `max_passes` warning (:19).
- **defects:** :21-28 lists `mesh_skeleton_with_mapping`, `HeatGeodesicSolver`, `heat_geodesic_distances`, `simplify_polyline`, `closest_points_on_polyline` under Unreleased/Added, but the mesh API was cut (commit 71d4c7c7) and none of these names exist in `src/compas_cgal` (verified grep: no match). :13 cites `docs/superpowers/state/sp1-gate-c-analysis.md` — directory does not exist on live (verified `ls`). :19 "max_passes truncation now warns" — the draft plan's Review Focus 4 requires it to fail loud instead; conflict.
- **unused capability / ideas:** :13 "full generation is O(n²) exact stock depletion and is deferred" — the draft plan's productive-only emission is exactly that O(n²) depletion inside generation; cost unmeasured.
- **claims:** :14 gouge-free radii from exact clearance — unverified here. :15 `stepover` bounds radial engagement per cycle — unverified. :32 `polygon_skeleton_clearance` exists — verified `src/compas_cgal/toolpath.py:242`. :13 `audit_toolpath_engagement` — verified `src/compas_cgal/engagement.py:413`; `scripts/engagement_baseline.py` verified present.

### CONVENTIONS.md | lines 194 | read 1-194
- **role:** doc (project rules)
- **used by:** agents
- **state:** live
- **generation:** n/a; rules. :161-176 "Single-Path Novel Architecture": one authoritative call graph; "separate generator algorithms remain isolated peers. They may not enter the exact-certified pipeline as a prefilter, fallback, or substitute proof". :150-155 names `docs/segment_site_mat.md` as the current stage page for exact-certified adaptive/MAT work. :193-194 Held figures must be real toolpaths.
- **defects:** :7 says the project uses Epick, while :43-45 says stock/engagement run on Epeck — both true for different modules, but the Epick-only rules at :9-35 (e.g. `approx_distance` from toolpath.cpp) do not apply to Epeck code. Mild. :39-40 "14-question review checklist" vs the user CONVENTIONS.md "12-question" — doc drift (unverified which is current; see exactness.md record).
- **unused capability / ideas:** none.
- **claims:** :150-155 segment_site_mat.md is the current stage page — contradicted in effect by the draft plan, which is the "one active plan" and is built on the straight skeleton, not the segment-site MAT.

### README.md | lines 23 | read 1-23
- **role:** doc
- **used by:** packaging
- **state:** live, stale (upstream compas_cgal readme; says nothing about toolpaths)
- **generation:** n/a
- **defects:** none for generation.
- **unused capability / ideas:** none.
- **claims:** "bindings for specific algorithms of CGAL" — n/a.

### docs/acknowledgements.md | lines 1 | read 1-1
- **role:** doc
- **used by:** mkdocs nav (unverified)
- **state:** live (empty heading only)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/api/compas_cgal.adaptive.motion_certificate.md | lines 3 | read 1-3
- **role:** doc (mkdocstrings stub)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a; exposes `compas_cgal.adaptive.motion_certificate`, the certified adaptive path's motion certificate module, which the gated generators do not use.
- **defects:** none found
- **unused capability / ideas:** the `adaptive/` certified generator (see segment_site_mat.md record).
- **claims:** none

### docs/api/compas_cgal.straight_skeleton_2.md | lines 1 | read 1-1
- **role:** doc (stub)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/api/compas_cgal.toolpath.md | lines 1 | read 1-1
- **role:** doc (stub)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/api/compas_cgal.types.md | lines 1 | read 1-1
- **role:** doc (stub)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/benchmarks.md | lines 163 | read 1-163
- **role:** doc (benchmark corpus)
- **used by:** mkdocs; documents `benchmarks/cli.py corpus`
- **state:** live
- **generation:** records the key generation finding: the engagement-controlled generator "regulates advance but not trochoid radius" (:130-133), so its measured max TEA saturates at ~141.8° for any requested cap 20-140° (:124-127), with 264/72/40/20 forced minimum advances at caps 20/40/60/80. Also: engagement control gives 0.19-0.43× the length of best compliant constant-spacing path on rect_20x12 (:114-118); constant spacing has no compliant result below 100° (:93-94).
- **defects:** generation-relevant finding only; doc itself fine. :13-16 `center_domain()` ~3,500× slower on oblique geometry, mechanism not established (see oblique_edge_cost.md).
- **unused capability / ideas:** the `necks` corpus ("the family that shows engagement control working", :45) and `analytic`/`degeneracy` families are not in the product gate (only rect_12x8, rect_20x12, L_shape). Pythagorean-triple congruence tests (:67-72) could be a generator invariance gate.
- **claims:** :124-127 saturation table — unverified (historical, not re-run). :127 forced minimum advances at cap 40 = 72 vs draft plan Stage C "88-215 forced minimum advances per pocket" at cap 40 — different pockets/generators, not contradictory. "0.31× at a 120° cap" — unverified.

### docs/benchmarks/held_reference_qualification.md | lines 12 | read 1-12
- **role:** doc (report)
- **used by:** linked from docs/benchmarks.md:11
- **state:** live, stale
- **generation:** n/a. "exercises the existing generator serially"; a pass means only "a non-empty path" (:12). No quality score.
- **defects:** primitive counts 31/70/58/317 (:7-10) are the pre-certificate-fix census (task-1-report:109-112), not the current 49/80/60/319 (task-1-report:158-161) — stale. Which generator ran is not named.
- **unused capability / ideas:** none
- **claims:** "does not ... establish Held parity" — consistent with the draft plan.

### docs/changelog.md | lines 0 | read 1-1
- **role:** doc (include of CHANGELOG.md)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found (one line, no trailing newline)
- **unused capability / ideas:** none
- **claims:** none

### docs/continuous_engagement.md | lines 1364 | read 1-1364
- **role:** doc (certified adaptive generator: exact continuous TEA oracle, candidate transactions, MAT traversal)
- **used by:** mkdocs; documents `src/compas_cgal/adaptive/*` (`generator.py:1336 generate_exact_adaptive_continuation`, `transaction.py`, `motion_certificate.py`, `motion_refutation.py`, `retrace_*`) and `continuous_tea_2`
- **state:** parked (the certified `adaptive/` generator is not gated and not in the draft plan; its doc says "Exact inter-route retrace, traversal/coverage closure, fresh terminal sealing, arbitrary-pocket evidence ... remain incomplete", :69-73)
- **generation:** describes the only generator in the repo that (a) certifies every link and circle against the exact pre-motion stock before committing it (`containment -> certify -> deplete -> add_sweep`, :461), (b) walks the exact MAT with a positive-radius MATHSM circle family plus a zero-guide advancing-segment family (:1186-1197), (c) searches candidates in furthest-progress, largest-radius order and stops at the first feasible (:587-603), and (d) cheaply rejects over-cap links with an exact station-refutation ladder (:767-864, 5.39× wall clock). Quality: **weak as a generator today** — it has never produced a complete pocket: it reaches route 2 of the L fixture and stops because every direct link from the finished leaf gouges (:1292-1319); `test_task13f_full_continuation` fails and `test_real_active_family_stops_at_unresolved_exact_event` newly fails (:1355-1364). Excellent as a certification substrate.
- **defects:** :1309-1311 the exact route is "a deterministic DFS edge-discovery ledger, not an Euler walk and not a physical-motion sequence", so consecutive routes are physically disconnected — the same defect class the draft plan's tree walk must handle. :874-879/:943-944 `IncompleteSegmentPartitionError` on 15/125 audits and a process-terminating `float_next ... inf` abort in `continuous_tea_2` event discovery — open. :1029-1031 MAT arms with clearance equal to tool radius have zero guide radius and need link-only advance (built later, Task 3-7).
- **unused capability / ideas:**
  - **swept-prefix theorem** (:1105-1128): a straight advancing cut whose start disk is clear can engage only the forward semicircle, so it certifies exactly a π (180°) cap in 0.04 s (:1329-1337). Directly relevant to draft-plan B2a (straight wall run-ins) and to stay-down links: the gated generators and the plan never use it. Limitation: it certifies only π, never 60° (the plan's slot threshold).
  - **exact station refutation ladder** (:767-864): 9 exact rational stations refute a motion over cap — a cheap exact "does this link exceed the cap" test the draft plan could use when placing bridges/run-ins.
  - **candidate ordering** furthest progress first, largest radius second (:593): an engagement-capped greedy step policy, the MATHSM/Held idea, never used by the gated generators.
  - **zero-guide edges** (:1199-1205): on MAT arms where clearance = tool radius, emit a straight advancing segment instead of a degenerate circle — the principled answer to the draft plan's "ρ ≈ 0.02 degenerate loops" and convex-corner problem (B2).
  - **"graph coverage order is not toolpath continuity"** (:1321-1327) and "No retract/re-plunge fallback ... is admitted" (:1342): the adaptive stage designed an exact inter-route retrace through accepted motion lineage; the draft plan chooses retract + clearance LINK + plunge instead (and defers stay-down to D3).
  - "A real first-passage circle reaches the certifier at 90 degrees, not the 120 degrees user cap" (:428-430): oriented-neck policy caps (second passage 80°) — a neck-aware cap schedule.
- **claims:** :1004-1006 "complete adaptive suite passed 518 tests" (2026-07-29) — historical, unverified now. :781 refutation 5.39× — unverified (measurement). :1349-1364 suite "4 failed, 882 passed" with two generator tests failing — unverified now; both tests exist (verified `tests/adaptive/test_generator.py:140,491`). :1339-1342 "next coherent stage is an exact inter-route retrace" — partially built: `src/compas_cgal/adaptive/retrace_segment_trial.py` and `retrace_transaction.py` exist on live (verified `ls`); whether route 2 now completes is unknown.

### docs/continuous_engagement_cost.md | lines 120 | read 1-120
- **role:** doc (performance verdict)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a directly; decides which certifier a generator can afford: sampled + guarded `engagement_2.cpp` at 129 ms per segment vs `continuous_tea_2` ~15,000 ms (:14-18), superlinear in contacted features (:29-32).
- **defects:** none in the doc. Generation consequence: any generator that certifies each candidate with `continuous_tea_2` (the `adaptive/` generator) is ~5 orders from Held (:3-4).
- **unused capability / ideas:** :93-97 the events along a one-parameter motion are univariate (degree ≤ 2 in the motion parameter); bivariate curve-pair analysis may be unnecessary — never pursued.
- **claims:** :16-17 129 ms vs 15 s — unverified (measurement, not re-run). :119-120 "Answer the univariate question before reviving it" — still open; no commit found addressing it (unverified by search of log subjects).

### docs/dev/arrangement-redesign/.gitignore | lines 6 | read 1-6
- **role:** build
- **used by:** git
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/dev/arrangement-redesign/compile.sh | lines 14 | read 1-14
- **role:** tool (prototype build)
- **used by:** no consumer found
- **state:** dead (references `/Users/jelle/mambaforge/envs/cgal-dev`, CGAL 5.6.1; conda was removed 2026-09-10 per project memory, so this script cannot run as written; uses `c++` outside pixi)
- **generation:** n/a
- **defects:** :6 hard-coded conda env path; violates pixi-only rule.
- **unused capability / ideas:** none
- **claims:** "CGAL 5.6.1" vs project CGAL 6.0.1 — the prototypes were validated on a different CGAL; unverified on 6.0.1.

### docs/dev/arrangement-redesign/findings.md | lines 127 | read 1-127
- **role:** doc (performance finding: engagement query)
- **used by:** mkdocs (dev)
- **state:** live (fix landed 2026-07-24 per :108)
- **generation:** n/a directly. Relevant to generation cost: engagement queries read locally from the Gps arrangement via `Arrangement_zone_2`, ~3.8-5.1× faster (:123-125); depletion via aggregated disk-chain union 21× per circle (:21).
- **defects:** :121 cites `docs/superpowers/state/engagement-zone-divergence.md` — directory absent on live (dangling). :67-69 "exact annular sector ... deferred" — depletion still uses disk chains.
- **unused capability / ideas:** the zone query answers "does this cutter rim touch material" locally; the draft plan's productive-only test (`_subtract_operation` then `Stock.exactly_equals`, a full Boolean + full equality per candidate loop) could instead ask the zone/arrangement whether the swept region meets any `contained()` face. Not considered by the plan.
- **claims:** :114 "0 cap_exceeded mismatches over 7200 comparisons" — unverified (historical). "LANDED" — verified in part: `src/engagement_2.cpp` exists; `engaged_arcs_zone` presence unverified.

### docs/dev/arrangement-redesign/gps_zone.cpp | lines 82 | read 1-82
- **role:** tool (prototype)
- **used by:** `compile.sh` only
- **state:** dead (prototype; production equivalent in `src/engagement_2.cpp` per findings.md:108)
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none beyond findings.md
- **claims:** none

### docs/dev/arrangement-redesign/harvest.cpp | lines 105 | read 1-105
- **role:** tool (prototype)
- **used by:** `compile.sh` only
- **state:** dead (prototype)
- **generation:** n/a
- **defects:** :83-87 sums spans with `to_double`/`atan2` (reporting only, acceptable for a prototype); :87 `span <= 0.0` then +2π mislabels a full circle as 2π only by accident of seam.
- **unused capability / ideas:** none
- **claims:** none

### docs/dev/arrangement-redesign/min.cpp | lines 21 | read 1-21
- **role:** tool (prototype)
- **used by:** `compile.sh`
- **state:** dead
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/dev/arrangement-redesign/pl_bench.cpp | lines 62 | read 1-62
- **role:** tool (prototype benchmark)
- **used by:** `compile.sh`
- **state:** dead
- **generation:** n/a
- **defects:** none relevant (landmarks PL crashes under dynamic insertion, :58-59)
- **unused capability / ideas:** none
- **claims:** none

### docs/dev/arrangement-redesign/stock_arr.cpp | lines 73 | read 1-73
- **role:** tool (prototype)
- **used by:** `compile.sh`
- **state:** dead (superseded by zone-in-Gps)
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none
- **claims:** none

### docs/dev/cgal-traits-aliasing-report.md | lines 520 | read 1-520
- **role:** doc (upstream CGAL defect report, CGAL/cgal#9638)
- **used by:** linked from `docs/exactness.md`
- **state:** live
- **generation:** n/a. Indirectly relevant: any generator that clones `Stock2` (the draft plan's productive-only test subtracts on a copy and compares with `Stock.exactly_equals`) depends on the clone not aliasing a freed traits object.
- **defects:** none in the doc. The CGAL defect itself is open upstream (unverified status of #9638).
- **unused capability / ideas:** :476-479 detection must be structural (assert arrangement traits pointer == owned traits).
- **claims:** :216 20/20 structural, :236 30/30 crashes — unverified (not re-run). Local mitigation commits exist on live: `57482731 fix: stock2 clone family shares one traits object`, `df6fb020 fix: exact region adopts its set, one traits owner per family`, `8d2b76cc fix: stock2 resolves its traits owner` (verified `git log`). :490-497 says a shared family traits does NOT close it; 8d2b76cc's subject suggests the later fix addressed that — unverified.

### docs/dev/exact-cgal-idioms.md | lines 437 | read 1-437
- **role:** doc (generic primer)
- **used by:** referenced by CONVENTIONS.md:39-41
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** :378-412 separate exact topology from a `Manufacturing_policy` (linear/angular tolerance, minimum feature size) — the product gate has no minimum-feature-size notion; the convex-corner uncut problem (draft plan B2) is exactly a "is this sliver manufacturable" policy question, which the plan frames as a gate ruling (D2a rest material).
- **claims:** :414-433 12-question checklist (CONVENTIONS.md:40 says exactness.md has a 14-question checklist; consistent — different docs).

### docs/engagement_audit.md | lines 227 | read 1-227
- **role:** doc (truthful engagement audit, validation)
- **used by:** mkdocs; documents `src/compas_cgal/engagement_audit/*` and `src/audit_*_2.cpp`
- **state:** parked (auditor convergence P0-P4 in draft plan parking lot)
- **generation:** n/a. Constrains generation if it ever becomes the gate: the native classifier rejects ramps — line primitive "ramp rejection" (:40) and "mixed-Z ramps ... fail before stock construction" (:127-128). Only vertical plunge/retract, clearance-plane transport, and cut-plane lateral line/circle/arc are admissible (:44-51).
- **defects:** :40/:127-128 ramp rejection is incompatible with the draft plan's D1 (helical bore as `OperationType.HELIX`) and with any ramped entry — a world-class requirement. If this auditor is ever promoted to authority, it will reject the entry the plan recommends. :89-92 only exact world-XY orthonormal frames accepted.
- **unused capability / ideas:** none for generation.
- **claims:** :217-227 "Stage 1 must still implement measure-before-deplete replay, native three-way certification ... no audit result can claim complete engagement compliance" — partly superseded: progress.md records Task 4B/4C native certification and replay committed (`7cda8f6` on live). Report aggregation, consumer migration and evidence matrix: unknown. The doc is stale relative to progress.md.

### docs/engagement_controlled_toolpath.md | lines 293 | read 1-293
- **role:** doc (gated generator `engagement_controlled_toolpath`)
- **used by:** mkdocs; documents `src/compas_cgal/engagement_toolpath.py`
- **state:** live (gated generator; behaviour unchanged since 2026-08-22 per brief)
- **generation:** engagement-regulated advance between full-circle trochoid loops on the C++ skeleton guide: integer station bisection, 33 exact `engagement_at` probes per candidate, exact annulus depletion for loops and a 6-curve capsule quad for bridges. Quality: **weak** against the world-class bar — guarantee is only at 33 evaluated positions (:16-23); radius is not regulated (benchmarks.md:130-133); every chain starts with a full-slot entry loop (:263-270); guide tangency G1 only on straight guide segments (:290-291), which is the source of the gate's 24/32/22 tangent breaks.
- **defects:**
  - :40-45, :283-286 bisection assumes engagement monotone in advance; project memory "Engagement is not monotone ... ladder search, never bisection" says it is non-monotone, and the generator still bisects (verified `src/compas_cgal/engagement_toolpath.py:42,616`). Effect: the largest admissible advance is not found; paths are longer than needed.
  - :263-270 "A helical or ramped entry would remove them. That is not implemented here." — the D1 problem, recorded 2026-08 and still open.
  - :124-129 between-position engagement measured 43.4° on a 40° cap — sampled guarantee is not a cap.
  - :180-184 guide radius stays a double surrogate; exact annulus about the true guide circle not representable.
  - :188-191 bridge depletion is an under-approximation (capsule side lines irrational).
- **unused capability / ideas:** :95-104 the trailing-lateral quadrant is where the cap breaks; probe placement must be uniform. :250-261 no local depletion for capsules (seed completeness).
- **claims:** :5-6 "holds every machining circle except the unavoidable chain-entry loops at or below the requested cap" — consistent with the draft plan's probe ("Every over-cap motion at 120° is a chain-entry loop"), verified only as reported by the plan; the doc is silent about the gate's degenerate loops (4/4/5) and tangent breaks (24/32/22) (draft plan table). :168/:235 "49 cuts / 441.8 length / 5 over cap (12x8)" — matches the plan's 5 cap exceedances on rect_12x8 (verified against plan table). :133 test exists (verified `tests/test_engagement_toolpath.py:148`). :287 "~27x slower" vs :7 "~86" — internal inconsistency (different stages of the doc).

### docs/exactness.md | lines 732 | read 1-732
- **role:** doc (exact-kernel doctrine, review checklist)
- **used by:** CONVENTIONS.md, mkdocs
- **state:** live
- **generation:** n/a, but :78-82 records the kernel split: `toolpath.cpp` (the guide/trochoid generator the gated generators build on) is Epick, "decisions are predicates on input-derived values". :86-92 tool-offset sweeps of arcs are not exactly representable in circle-segment traits, so stock uses "certified under-covering disk chains". :432-449 `Manufacturing_policy` with `minimum_feature_size` ("SP2 direction of travel") — never built.
- **defects:** :677 "eighteen questions" vs CONVENTIONS.md:39-40 "14-question review checklist" and user RULES "12-question" — doc drift. :86-92 "disk chains" for depletion is stale: engagement_controlled_toolpath.md:140-160 says full turns are now exact annuli.
- **unused capability / ideas:** `Manufacturing_policy.minimum_feature_size` (:440-444) is the principled form of the draft plan's D2(a) "rest material" ruling: corner material narrower than a feature size is policy, not a gate exemption.
- **claims:** :652-654 status table: Stock2, ExactRegion2, ReachableMaterialPredicateStorage2 fixed (`df6fb020` verified on live via `git log`). :577-583 `57482731` did not fix Stock2 — consistent with later `8d2b76cc` (verified in log).

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton.md | lines 15 | read 1-15
- **role:** doc (example page)
- **used by:** mkdocs
- **state:** live (upstream compas_cgal example)
- **generation:** n/a
- **defects:** :3 image `../assets/images/example_straight_skeleton_2_interior_straight_skeleton.png` — existence unverified (commit 808cdcca "drop orphaned assets").
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton.py | lines 38 | read 1-38
- **role:** doc (example script, needs compas_viewer)
- **used by:** example md include
- **state:** live
- **generation:** n/a. Shows the skeleton graph exposes `inner_bisector` / `bisector` edge attributes — the same distinction the draft plan's `GuideTree` needs (bisector edges vs contour edges).
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_offset_polygon.md | lines 15 | read 1-15
- **role:** doc (example page)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_offset_polygon.py | lines 68 | read 1-68
- **role:** doc (example script)
- **used by:** example md
- **state:** live
- **generation:** n/a. Exposes `offset_polygon` / `offset_polygon_with_holes` (CGAL straight-skeleton offsets). An offset (contour-parallel) capability exists in the package and no generator uses it — e.g. for a final wall/finishing contour or to define the convex-corner rest region.
- **defects:** none found
- **unused capability / ideas:** `offset_polygon` for rest-material region (draft-plan D2a) or a final contour pass.
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_weighted_offset_polygon.md | lines 15 | read 1-15
- **role:** doc (example page)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_weighted_offset_polygon.py | lines 36 | read 1-36
- **role:** doc (example script)
- **used by:** example md
- **state:** live
- **generation:** n/a. `weighted_offset_polygon` (per-edge weights) — a weighted skeleton could bias a guide away from short edges; unused.
- **defects:** none found
- **unused capability / ideas:** weighted straight skeleton (unused).
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_with_holes.md | lines 16 | read 1-16
- **role:** doc (example page)
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_straight_skeleton_2_interior_straight_skeleton_with_holes.py | lines 56 | read 1-56
- **role:** doc (example script)
- **used by:** example md
- **state:** live
- **generation:** n/a. Skeleton with holes produces a cyclic graph — the draft plan's Review Focus 1 (island).
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/examples/example_toolpath_trochoidal_mat.md | lines 17 | read 1-17
- **role:** doc (example page for the C++ trochoid generator)
- **used by:** mkdocs
- **state:** live
- **generation:** documents `trochoidal_mat_toolpath_circular`: arcs+lines, lead-in/out, links, retract/plunge, "Radius maximization from centerline toward pocket boundary" (:13).
- **defects:** :3 calls it "adaptive" — it is stepover-based, not engagement-adaptive (engagement_controlled_toolpath.md:27-33 measures 136.4° at a nominal 120°).
- **unused capability / ideas:** none
- **claims:** :3 "adaptive trochoidal" — contradicted (engagement_controlled_toolpath.md:27-33).

### docs/examples/example_toolpath_trochoidal_mat.py | lines 294 | read 1-294
- **role:** doc (example script, viewer)
- **used by:** example md
- **state:** live
- **generation:** calls `trochoidal_mat_toolpath_circular(..., link_paths=True, optimize_order=True, cut_z, clearance_z)` (:38-50). `optimize_order` is a greedy nearest-neighbour order starting from the longest chain, never reversing a chain (verified `src/toolpath.cpp:961-992`) — this is the ordering the gated generators inherit, and a cause of the draft plan's root cause 1.
- **defects:** :136-143 `1e-24` epsilon in a viewer helper (reporting only; acceptable). `match_skeleton_edge` matches skeleton edges to groups by mean distance — visualization only.
- **unused capability / ideas:** tool-cone slider along `result.polyline` (visual playback) — a debugging aid for the tree walk.
- **claims:** :41 comment "max radial engagement per cycle (0.4 x tool diameter)" — stepover is a proxy, measured to under-read engagement (engagement_controlled_toolpath.md:27-33).

### docs/held_motion_coverage.md | lines 765 | read 1-765
- **role:** doc (Held reproduction: coverage, contour/medial construction, thinning)
- **used by:** mkdocs; documents `benchmarks/held_*`, `src/boundary_normal_circle_2.cpp`, `src/native_boundary_curve_2.cpp`, `src/stock_circle_redundancy_2.cpp` (`Stock2::can_remove_circle`, verified `src/stock_circle_redundancy_2.cpp:5`)
- **state:** live (Held lane); exact coverage replay set aside 2026-09-08 (:3-26)
- **generation:** records Held's construction and where this repo diverges. Held: machinable target first (opening by r, convex arcs regularized with offset slightly > r, :43-51); machining centre c halfway between offset contact q and medial point m, qm is a diameter (:60-66); spacing Eq. 4 `|c_i - c_prev| + rho_i - rho_prev <= 2r` (:70-75); paired boundary/medial advancement (:77-83). Repo divergences (:88-103): target is the original polygon; guide is the straight skeleton (`_guide_chains`), not the true MAT; `_circle` moves centres without recomputing m; traversal is not paired. Quality of the Held lane generation: **adequate for reproduction, weak as a general generator** (needs publisher start data — brief).
- **defects:**
  - :251-258 at a rectangular convex corner the stationary contact (r,r) cannot be contained in any predecessor's outer disk (`sqrt(5)·rho > rho`): the predecessor-only engagement model cannot certify a corner approach. Directly relevant to draft-plan B2/B2a.
  - :92-93 straight-skeleton guide ≠ medial axis (max sampled skeleton-to-native-medial discrepancy 0.38436 mm on Figure 5, :167; reconstruction drift up to 5.5102 mm, :116-118).
  - :701-707 shared thinning regression red at all three scales (engagement-only selection drops a circle that alone covers material).
  - :3-15 exact coverage replay does not finish (1,200 circles > 4 min).
  - :714-720 Monstera interior gap 1.0282 mm from nearest motion for a 1 mm tool.
  - :105-109 `SegmentSiteMedialAxis.build` delegates to `canonical_l_shape_mat_graph` — only one six-vertex L fixture supported (at time of writing).
- **unused capability / ideas:**
  - **Machinable target** (Held §2.1, :43-58): define the pocket as the union of contained tool disks before generating; corners outside it are reported residual, not uncut. This is the paper-backed version of draft-plan D2(a); the plan treats it as a gate ruling only and does not cite Held's precedent.
  - **Contour-aware qm-diameter circle** `BoundaryNormalCircle2` / `NativeBoundary2.circle_on_piece` (:132-160, :332-365): native exact construction of the Held circle from a boundary contact. Gated generators and the draft plan use skeleton-station clearance circles instead.
  - **Stationary events at convex corners** (:191-194, :211-216): clearance = r gives zero guide radius; Held treats it as a transition event, not a circle. Draft plan's ρ>r floor is the same fact without the event semantics.
  - **`Stock2.can_remove_circle`** (:639-645): native local-residual redundancy test (candidate annulus minus neighbours' annuli and connector quads). The draft plan's productive-only emission re-derives this with a full `_subtract_operation` + `exactly_equals`; `can_remove_circle` is local and already exists.
  - **local repair, not global densification** (:612-624): rejected +161%/+82% length experiments; add motion only for reachable gaps.
  - boundary-order spacing, qualified contacts, corner approaches (Held lane, see gallery record).
- **claims:** :14-15 / :756-760 exact replay does not finish — consistent with commit `334b2f2f coverage: exact full-motion replay set aside` (verified log). :155 "27 new tests ... pass" — unverified now. :528-537 "Actual Figure 5 runtime blocker ... 337 s" — contradicted by the same page :367-371 and :495-497 (now 8 ms / 1.33 ms; fix commit `282152b1 ... 337s to 8ms` verified in log). The page keeps both the old blocker and the fix.

### docs/held_performance.md | lines 72 | read 1-72
- **role:** doc (Held lane performance gate)
- **used by:** mkdocs; `tests/benchmarks/test_held_toolpath_performance.py`, `benchmarks/results/held_toolpath_baseline.json`, pixi tasks `held-toolpath-*`
- **state:** live
- **generation:** Held lane produces Figure 5 (2,013 circles, 26,436.6 mm, 79.992°), Figure 8 upper boundary spacing (1,684 circles), crossed skis (1,616), Monstera (4,367; 449.8 s) (:12-18). "In every workload the initial placement is the cost" (:9-10). Engagement is the paper's standard predecessor model (:21-23), not exact stock replay.
- **defects:** :21-23 max engagement measured in Held's predecessor model only; never scored on the product gate (consistent with draft plan). :18 Monstera 449.8 s generation — far from Held's 3-100 ms.
- **unused capability / ideas:** boundary spacing vs lane spacing: 1,684 vs 4,183 circles and 18,236.8 vs 44,186.9 mm on Figure 8 upper (:15-16) — boundary-order spacing halves the path. The gated generators space along the skeleton, not the boundary.
- **claims:** baseline at `fa1f53a3` — verified commit exists in log (`fa1f53a3 test: anchor segment attestation bytes`). Numbers unverified (not re-run).

### docs/held_pfeiffer_reference_pockets.md | lines 73 | read 1-73
- **role:** doc (Held input corpus)
- **used by:** mkdocs; pixi `held-reference-qualify`
- **state:** live
- **generation:** n/a; inputs. Polygon projection is consumed by `PocketSpec` and "the existing generator" (:26-27, :50-51).
- **defects:** :53-58 primitive counts 31/70/58/317 — stale vs task-1-report:158-161 (49/80/60/319). :71-73 claim boundary is honest.
- **unused capability / ideas:** these four pockets plus the gate's three would give the product gate curved/necked/island cases; draft plan Stage C proposes running the product generator on them.
- **claims:** :51 "The current run passed all four" means non-empty path only (qualification page :12) — verified consistent.

### docs/held_toolpath_gallery.md | lines 77 | read 1-77
- **role:** doc (gallery)
- **used by:** mkdocs; pixi `held-toolpath-figures`
- **state:** live
- **generation:** Held lane outputs: Figure 5 corner-aware refinement; Figure 8 boundary-order spacing; Monstera with **qualified contacts and corner approaches** — "concave corners partitioned into approach, shared-contact rotation and departure" (:56-58). Blocked: contour-aware (b), intermediate contour (c), MATHSM (d) on Figure 5 (:24-27); Figure 6 contour/constant-spacing curves (:69-70); Figure 7 engagement maps (:72-77).
- **defects:** :66-69 Figure 6 "cap compliance unaudited".
- **unused capability / ideas:** **corner approaches** (approach / shared-contact rotation / departure at concave corners) — a built corner strategy the gated generators lack; **corner-aware engagement refinement**. Both exist only in the Held lane.
- **claims:** :24-27 contour-aware and MATHSM blocked — draft plan says they exist on `codex/held-reference-corpus` (not on live). Consistent.

### docs/index.md | lines 4 | read 1-4
- **role:** doc
- **used by:** mkdocs
- **state:** live (upstream landing text)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/installation.md | lines 46 | read 1-46
- **role:** doc / build
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** :34-40 `coverage_exact_core` and nanobind modules depend on `external_downloads` — unverified.

### docs/license.md | lines 2 | read 1-2
- **role:** doc
- **used by:** mkdocs
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### docs/licenses/cgal-adaptive-package-audit.md | lines 48 | read 1-48
- **role:** doc (license gate)
- **used by:** mkdocs
- **state:** live, status `CGAL_LICENSE_GATE_ERROR` (:3)
- **generation:** n/a for path quality, but it blocks distribution of every generator: 2D Arrangements and 2D Regularized Boolean Set Operations (used by `_stock_2`, which every engagement generator depletes through) are GPL-3.0-or-later or commercial (:19,:21), while the repo is LGPL-3.0 (:5-7). Segment Delaunay Graph / Voronoi adaptor (the exact segment-site MAT route) are also GPL (:22-24).
- **defects:** :41-45 "Tasks 1-16 are blocked until the repository records" a GPL decision or commercial entitlement — no such record found on live (unverified by search; no doc in this slice records it).
- **unused capability / ideas:** none
- **claims:** :28-34 Minkowski sums deliberately not selected.

### docs/local_depletion.md | lines 197 | read 1-197
- **role:** doc (stock depletion performance)
- **used by:** mkdocs; documents `Stock2::subtract_*_local`
- **state:** live
- **generation:** n/a for path shape; depletion cost for generators. Local annulus depletion 2.4× faster on L-shape, 1.3× slower on 12x8 (:9-15). Insert zone walk is 92% of the cost (:65-74) — "removes the rebuild, not the sweep".
- **defects:** :17-22 "nothing calls it yet" — contradicted: `engagement_toolpath.py` and siblings call `subtract_arc_sweep_local` (verified grep: 3 occurrences across `src/compas_cgal/engagement_*_toolpath.py`). :127-131 "`exactly_equals` alone is not a sufficient gate" for representation validity — the draft plan uses `Stock.exactly_equals` as its removal decision; that is a point-set question, so it is the right primitive for "removed material", but the plan should not use it as a stock-validity check.
- **unused capability / ideas:** none new.
- **claims:** :49-51 "rect 12x8 = 49 cuts / 441.8 / 5 truly exceeding; L-shape = 112 / 595.3 / 7" — consistent with engagement_controlled_toolpath.md:168; the draft plan's L_shape cap exceedances = 6 (not 7) — different counter (gate vs 12-position replay), not a contradiction.

### docs/loop_radius_degeneracy.md | lines 202 | read 1-202
- **role:** doc (generation theory: degeneracy, entry, bridge tilt)
- **used by:** mkdocs; cited by `engagement_rho_toolpath.py`, `engagement_spiral_entry_toolpath.py`, `benchmarks.quality.DEGENERATE_LOOP_RATIO`
- **state:** live
- **generation:** three load-bearing results. (1) `rho <= r` is degenerate; on the middle-curve construction `rho_max = (d - r)/2`, so any station with clearance `d <= 3r` cannot carry a trochoid (:23-42); MAT-centred threshold is `d <= 2r + c` (:44-51). (2) Every non-degenerate entry loop after a rim plunge peaks at 360° (:69-90); with a **centre plunge** the first loop peaks at `360° − 2·acos(rho/2r)`, infimum 240° (:104-125); a concentric ramp measured 243.20 → 118.56 → ... max step 124.64° vs 329.49° for a rim plunge (:127-131), implemented in `engagement_spiral_entry_toolpath`. (3) Bridge tilt = `atan(dρ/ds)` independent of advance: 35.26° on a rectangle's corner bisector, all 24 tangent breaks on rect_12x8 (:151-175); fix = entry at the **external common tangent** (:184-190), valid while `|dρ| <= advance`.
- **defects:** none in the doc. Generation defects it identifies: gated generators use the normal-placed entry (24 tangent breaks), rim plunge (329.49° step).
- **unused capability / ideas:**
  - centre plunge + concentric ramp (spiral entry) — built, restored (`d4311356`), not used by gated generators; costs two tangent breaks per ramp turn (:132-134).
  - :136-149 a helical entry "is a change to the depletion model, not to a generator": `benchmarks.depletion._replay_kind` raises `UnreplayableOperationError` on any Z+XY move. The draft plan's D1 agrees (needs a ruling).
  - :98-102 "**Do not** exclude the entry motion from the step criterion" — the draft plan's D1 excludes the HELIX from the *planar engagement cap*; it must still keep the helix→first-loop step in `max_engagement_step`, or it reproduces exactly the exclusion this page forbids. The plan is silent on the step.
  - external common tangent bridge — used by `rho_regulated_toolpath` (:199-200), not by the gated generators.
- **claims:** :24-26 `_mathsm_geometry` at candidates.py:134 and :1115 — line numbers unverified. :125 "first admissible rung 1.048 r measures 243.20°" — unverified. :173-175 "all 24 tangent breaks are the four corner chains' six junctions each, all at 35°" — consistent with draft plan table (24 tangent breaks on rect_12x8).

### docs/machining_metric_validity.md | lines 333 | read 1-333
- **role:** doc (gate metric validity)
- **used by:** mkdocs; `tests/benchmarks/test_quality.py`
- **state:** live
- **generation:** records that a **chain-ordering generator** "sequences its chains so the tool never lifts ... takes `material_entries` from 5 to 1 ... links chains at cutting depth" (:136-141). That generator was `engagement_ordered_toolpath.py` (commit `7af7a24c`), **deleted** in `9fa6eea8 refactor: drop two unreachable toolpath generators` (verified `git show --stat 9fa6eea8`: 685 lines removed). Its docstring (verified `git show 9fa6eea8^:src/compas_cgal/engagement_ordered_toolpath.py`): orders and orients chains so every chain after the first begins in cleared stock, takes the link at cutting depth only when every evaluated position is under cap by exact `engagement_at`, offers both orientations of every chain to the same predicate, cap-checks bridges too, and names refused chains in `forced_entries`.
- **defects:** :289-333 engagement is not translation-invariant at exact rim/cleared-boundary tangency (90.57° vs 0.0° under an exact shift); `test_moving_the_pocket_across_the_table_changes_no_metric` deliberately red; layer not established. This can flip a generator's accept/reject at exact tangencies (e.g. dyadic fixtures). :148-151 with chain-ordered stay-down links rect_20x12 failed the radius-step criterion at 3.988 under the rapid-only run rule.
- **unused capability / ideas:**
  - **the deleted stay-down chain-ordering generator** is the draft plan's Task 2-3 tree walk plus D3 stay-down links, already built once with exact link admission and orientation choice. The draft plan neither cites nor recovers it.
  - run boundary = rapid OR `path_index` change (:127-153): a tree walk that links at depth across chains must keep emitting distinct `path_index` per chain, or the radius-step metric merges runs.
  - :234-242 "A straight plunge is the defect. Helical ramp entry is under construction" — never landed on live.
- **claims:** :112-116 corrected loop-radius steps 1.368/1.376/1.376 — unverified. :172-173 corner chain radii `0.0155, 1.3840, 2.7524, 2.9980` — consistent with draft plan "ρ ≈ 0.02" corner entries. :275-287 at 120° the two gated generators' streams are structurally identical on every pocket — consistent with draft plan table ("engagement_controlled = radius_regulated").

### docs/machining_quality.md | lines 324 | read 1-324
- **role:** doc (the product gate's specification; the draft plan's "Spec")
- **used by:** mkdocs; `tests/benchmarks/test_quality.py`, `benchmarks/quality.py`
- **state:** live
- **generation:** defines the 12-criterion gate and measures both gated generators (identical streams at 120°, :89-90, :176-184). Findings about generation:
  - uncut residue is **not corner-confined**: farthest uncut sample 3.28 units from any polygon vertex on rect_12x8; resolution-stable 0.904/0.804/0.865% (:154-158). The unregulated constant-spacing generator leaves 0.152%, a fifth as much (:157-158).
  - 10 redundant operations are full-size loops (ρ = 2.998, 2.370, 1.998) with zero engagement: "the guide revisits stations it has already cleared" (:152-154).
  - degenerate loops are *not* redundant: removing them raises residue 0.284% → 0.351% (:195-197).
  - every cap exceedance is an entry cut (:160-168); worst non-entry 119.19° (rect_12x8), 116.36° (L).
  - slotting = 4 only on rect_20x12, a corner phenomenon; the straight cut leaving the degenerate loop "runs 1.956 into the corner at 140°" (:190-195).
  - W ≤ 4r channels: both generators emit 68 disk-sweep "circles"; "decline or switch strategy" (:305-313).
  - generator operates on the rising limb of chip thinning (50-110° rim arc), no lower engagement bound anywhere (:213-224).
- **defects:** (as spec) :39 uncut is measured against the tool-reachable region on a grid ("not one reachable sample survives", :288) — sampled, not exact. :290 slot threshold = half the cap is a judgement. :296-298 `material_entries` not gated ("one per chain is unavoidable without a helical entry") — stale after the deleted chain-ordering generator showed 5→1.
- **unused capability / ideas:**
  - **minimum-engagement / rubbing floor** (:213-224, :285): the gate reports `low_chip_thickness_ratio` and rubbing but no generator holds a lower bound or compensates feed. World-class adaptive clearing (constant engagement) needs both bounds.
  - **narrow-channel strategy switch** (:305-313): decline or switch (e.g. straight slotting pass at reduced feed) when `W ≤ 4r`.
  - feed-limited time model `v ≤ sqrt(a_max/κ)` exists (◈ machine) — no generator optimises curvature.
- **claims:** :96-133 table matches the draft plan's table for both gated generators digit-for-digit on the listed criteria (verified by comparison: uncut 0.008041→0.0080, degenerate 4/4/5, redundant 10/10/10, cap 5/9/6, slotting 0/4/0, step 329.5/337.5/325.4 vs plan 329/338/325, tangent 24/32/22). :301-303 "a middle-curve implementation ... does not exist yet" — contradicted in part: `src/compas_cgal/adaptive/middle_curve_path.py` and `candidates.py` middle-curve enumeration exist (verified `ls`), though they do not produce gate-scorable paths.

### docs/measurement_claims.md | lines 33 | read 1-33
- **role:** doc (ruling)
- **used by:** mkdocs; `benchmarks/mathsm.py` depends on it (:15-16)
- **state:** live
- **generation:** ruling: maximum engagement after entry is **non-monotone in spacing** (rect_20x12, 2 mm tool: spacing 0.025 → 131.14°, 0.1 → 98.73°) and also in guide radius (:18-33). "A selector may not bisect on spacing or stop at the first compliant sample" (:25-26).
- **defects:** generation consequence: `engagement_controlled_toolpath` bisects the advance (spacing) (verified `engagement_toolpath.py:42,616`) — contradicts this ruling; the doc names only `benchmarks.mathsm` as honouring it.
- **unused capability / ideas:** `benchmarks.mathsm` sweeps the complete spacing set and takes the shortest compliant trial — MATHSM constant-spacing baseline exists as a benchmark, not a generator.
- **claims:** :18-21 numbers — unverified (measurement).

### docs/number_types.md | lines 1323 | read 1-1323
- **role:** doc (number-type doctrine; stages 0-2 record)
- **used by:** mkdocs; CONVENTIONS.md
- **state:** parked (stages 3-6 in draft plan parking lot)
- **generation:** n/a. Lane L1 = Epick `double` in `toolpath.cpp` (the guide generator) (:26). Exact cost model: CORE exact-zero identity decisions cost 10⁵-10⁷× a generic sign (:86); integer fixtures hide it (:90-95).
- **defects:** :1285-1301 stage 2 does not deliver invariant 4 (filter under every branch); `exact-rational-v1` names two incompatible framings; :669-675 `AttestationByteDriftError` declared, never raised; :472-486 real pipeline never run under GMP; filter-failure rate of real workload unknown.
- **unused capability / ideas:** none for generation. Generation-relevant warning: any new exact generator decision (e.g. the draft plan's exact root choice by "largest exact clearance") must be tested on generic doubles, not only integer rectangles (:90-95, checklist item 9).
- **claims:** :701 `541b6ede` single mixed-radical sign — verified in log. :1040 `003e144c` station source carries rationals — verified in log. :1280 "4 failed, 787 passed" at close of stage 2 — unverified now.

### docs/oblique_edge_cost.md | lines 101 | read 1-101
- **role:** doc (performance finding)
- **used by:** mkdocs; cited by benchmarks.md, continuous_engagement.md, continuous_engagement_cost.md
- **state:** live
- **generation:** generation-relevant through cost: `ReachableDomain2.center_domain()` (used by `survey_path`, i.e. the gate's replay, :57-60) costs 5 ms axis-aligned vs 17.5 s oblique integer quad, > 90 s rotated (:25-31). Every gate pocket is axis-parallel (:74-80). Any generator that builds the centre domain (e.g. for plunge-clear or legality decisions) inherits the cliff.
- **defects:** mechanism not established (:62-70); hypothesis of compounding radicals contradicted by the 3-4-5 row.
- **unused capability / ideas:** :99-101 "decide whether the domain is needed at all" — a cheaper sufficient legality test.
- **claims:** :15-18 "The 71.9 ms/pocket figure" Held-parity claim — the cited figure is not found elsewhere in this slice; unverified.

### docs/plans/2026-02-28-full-circle-trochoids-design.md | lines 266 | read 1-266
- **role:** doc (historical design + plan for `toolpath.cpp`)
- **used by:** no consumer (restored design record, `d6e3bb0e`)
- **state:** superseded-by current `toolpath.cpp` full-circle generator (implemented; `trochoidal_mat_toolpath_circular` emits circles)
- **generation:** decided full 360° circles with same-side external-tangent transitions (:12-19, :105-126) replacing alternating half-arcs; nearest-neighbour ordering by `has_smaller_distance_to_point` (:56-67) — verified present at `src/toolpath.cpp:984`.
- **defects:** :257 "Tangent alignment test may regress — Full circles exit tangentially by definition — should improve" — the gate measures 24/32/22 tangent breaks on the generators built on this guide, so the claim did not hold end-to-end (bridges tilt, loop_radius_degeneracy.md:151-175). :263 uses `pip install` (pre-pixi).
- **unused capability / ideas:** "seamless tangent point when circles overlap" (:15, :116) — zero-length transition; the half-arc zigzag (:94-103) is the commercial-style trochoid (partial arcs, less air cutting) and was abandoned for full circles; full circles cut air on the trailing half (gate recut_fraction ~0.9).
- **claims:** :25 "CGAL's linear kernel has no built-in tangent computation" — plausible, unverified.

### docs/plans/2026-02-28-tessellated-polyline-output.md | lines 397 | read 1-397
- **role:** doc (historical plan)
- **used by:** no consumer (restored design record)
- **state:** superseded-by implemented `ToolpathResult.polyline` (example uses `result.polyline.tolist()`, verified `docs/examples/example_toolpath_trochoidal_mat.py:258`)
- **generation:** n/a (visualization output). C++ tessellation from `ToolpathPrimitive` with signed sweep.
- **defects:** :54-63 tessellation with `atan2`/`cos`/`sin` — output only, acceptable. :166/:263 `pip install` (pre-pixi).
- **unused capability / ideas:** none
- **claims:** none checkable beyond implementation presence.

### docs/plans/2026-02-28-toolpath-cgal-redesign.md | lines 806 | read 1-806
- **role:** doc (historical plan: rewrite `toolpath.cpp` on Epick kernel types)
- **used by:** no consumer (restored design record)
- **state:** superseded (implemented and later changed; `toolpath.cpp` now has exact-clearance station radii per CHANGELOG.md:14)
- **generation:** records the guide generator's original design: per skeleton halfedge, `trochoid_circles` with **linearly interpolated** radius `r0 + t(r1 - r0)` (:224-235) and `compute_radius = mat_scale·(d − r − radial_clearance)` capped by max/min trochoid radius (:444-454); external-tangent chain (:240-358); **eliminated a DP tangent solver** (`choose_arc_transition`, `Cost`, `TransitionDecision`, :515-520, :745-760).
- **defects:** :262-265 `external_tangents` clamps `|Δr| >= length` to `length − 1e-9` — an epsilon fudge where no external tangent exists (`|dρ| > advance`); current `toolpath.cpp` has no such literal (verified grep: only `OUTPUT_DEDUP_TOL = 1e-9` at :47), so presumably removed — unverified how the case is handled now. :224-235 linear radius interpolation gouges at reflex vertices (fixed per CHANGELOG.md:14, "up to 89% of tool radius").
- **unused capability / ideas:** the deleted **DP transition solver** chose tangent sides/transitions by a cost — an optimisation over the chain that no current generator has; abandoned for the fixed-side rule.
- **claims:** :31 "1439 → ~450 lines" — current `src/toolpath.cpp` length unverified here (memory says ~970).

### docs/plans/2026-07-23-adaptive-clearing-design.md | lines 189 | read 1-189
- **role:** doc (programme design: certified adaptive clearing SP1-SP3)
- **used by:** no consumer (restored design record); SP1 built (`Stock2`, `engagement_2`)
- **state:** SP1 built; **SP2 generator never built** (no `src/clearing_2.cpp`, no `compas_cgal/clearing.py`, no history on any branch — verified `git log --all -- src/clearing_2.cpp src/compas_cgal/clearing.py` empty); SP2b, SP3 never built
- **generation:** the most complete world-class generator design in the repository, never implemented:
  - **contour-parallel offset-loop passes** from straight-skeleton interior offsets per connected region (:98-101), radial step derived from `tea_target` (:104-105);
  - **propose → certify → adapt → commit** (:102-114); adapt order: shrink radial step → insert certified trochoidal relief (corners, slots, first channels) → split the pass;
  - **entry per connected region at the skeleton's maximum-clearance vertex** (:115-117) — the same root rule as the draft plan's tree walk;
  - **target band** `tea_target ± tea_band` with ≥ 95% of engaged length in band (:34-35) — constant engagement, not just a cap;
  - **SP2b spiralization** (Held/Spielberger-style) to remove loop seams (:122-124);
  - **SP3 chip-thinning feeds** from `EngagementReport` (:144);
  - sweeps inserted "slightly fat" by `radial_clearance` to avoid grazing tangencies (:74-77); prune edges interior to cleared regions (:78-80).
- **defects:** :74-77 "slightly fat" insertion is a margin hack the later exactness doctrine forbids. :83 landmarks point location — later measured to crash under dynamic insertion (pl_bench.cpp:58-59). :86-94 station certificate "exact at stations plus conservative refinement" was implemented (CHANGELOG.md:11-12) but the gated generators use 33 probes without the refinement bound.
- **unused capability / ideas:** all of SP2/SP2b/SP3 (above). The draft plan keeps the trochoid-on-skeleton family only; contour-parallel offset passes with trochoidal relief only where regulation fails is how Fusion Adaptive / VoluMill clear open regions with far less air cutting (recut ~0.9 on full-circle trochoids, machining_quality.md:104).
- **claims:** :40-42 "differentiator over VoluMill / iMachining / Fusion Adaptive, which ... publish no certificates" — unverified marketing claim.

### docs/plans/2026-07-24-engagement-zone-query.md | lines 67 | read 1-67
- **role:** doc (plan: local zone engagement query)
- **used by:** no consumer
- **state:** closed — landed (findings.md:108-127; `engaged_arcs_zone` present at `src/engagement_2.cpp:310`, verified grep)
- **generation:** n/a (engagement query speed for all generators)
- **defects:** :50-51 failure-path doc `docs/superpowers/state/engagement-zone-divergence.md` referenced; directory absent on live.
- **unused capability / ideas:** none
- **claims:** "Same exact result" — supported by findings.md:114-121 (0/7200 mismatches), unverified now.

### docs/plans/2026-07-23-sp1-stock-tea-kernel.md | lines 1508 | read 1-1508
- **role:** doc (SP1 implementation plan: Stock2, engagement_at, certify_segment_tea, audit, oracle)
- **used by:** no consumer (restored design record)
- **state:** built (superseded in detail): `Stock2`, `engagement_at`, `certify_segment_tea` (verified `src/engagement_2_bindings.cpp:60`), `tests/test_engagement_oracle.py` present; `tests/test_stock_gate.py` (Task 9 wall-clock gate) absent; `docs/benchmarks/engagement_baseline.md` (Task 8) absent on live (verified `ls docs/benchmarks`).
- **generation:** n/a directly; builds the substrate every engagement generator uses. Records the **station certificate with a growth bound** (`tea_growth_bound`, :969-1010; amended to a factor-2 guarded threshold, :886-919) — a continuous-between-stations guarantee for straight motions. The gated generators do not use it for their 33 probes (engagement_controlled_toolpath.md:16-23 says no continuous guarantee).
- **defects:** :599-654 original design had `CAP_BOUND_DEFLATION` and a `1e-12` merge gap (amended away; the case study in exactness.md:656-673). :1097 audit samples arcs at 20 stations with a Python growth bound. :380-391 capsule as disk chain (later replaced by the quad, engagement_controlled_toolpath.md:186-236).
- **unused capability / ideas:** `certify_segment_tea` (exact stations + guarded growth bound) could certify the draft plan's straight wall run-ins (B2a) and bridges continuously instead of sampling; `certify_arc_tea` mirror exists per CHANGELOG.md:11.
- **claims:** :1363 "This baseline is the number SP2 must beat" — SP2 never built (see adaptive-clearing-design record); the baseline report file is absent on live.

### docs/radius_ladder_ablation.md | lines 56 | read 1-56
- **role:** doc (historical ablation for `radius_regulated_toolpath`)
- **used by:** mkdocs
- **state:** live (historical)
- **generation:** ranking vs refinement halves of the radius ladder: refinement alone produced a **213.6° slotting cut** on 6x4 (:36, :40-42); both ungated worse than ranking alone (:37, :42-44). Forced stations: "lowering the peak always costs a circle" (:46-52).
- **defects:** :19-23 two incompatible historical assertions (8 vs 12 over-cap circles) left unreconciled.
- **unused capability / ideas:** none
- **claims:** not reconstructible from shipped knobs (:12-17) — by its own statement unverifiable.

### docs/radius_regulated_toolpath.md | lines 280 | read 1-280
- **role:** doc (gated generator `radius_regulated_toolpath`)
- **used by:** mkdocs; documents `src/compas_cgal/engagement_radial_toolpath.py`
- **state:** live (gated; identical stream to engagement_controlled at 120°, machining_quality.md:176-184)
- **generation:** adds a 40-rung radius ladder **scanned, never bisected** (:59-65) with two admissibility conditions: no probe over cap, and the loop still reaches uncut stock (`Stock.contains` at deepest reach, :89-117). Multiple passes per chain until each station's maximal circle is emitted (:202-218); one plunge/retract per pass (:261-263; 9 plunges over 5 chains on 20x12 at caps 40-100, :208). Quality: **weak** — cuts over-cap circle count 1.6-2.7× (:5-8) but the worst loop stays 125.8° at caps 40/60/80 (:119-131): corner-apex stations with radius 0.057 against a 1.0 tool have one rung (:156-171); bridges unregulated, worst 141.8° (:198-200).
- **defects:**
  - :28-32 / :198-200 **bridge cuts are not regulated by either gated generator** (worst 141.8°). The draft plan reuses rho's external-tangent bridges and emits a bridge only if it removes material, but does not say the bridge is cap-checked.
  - :208 one plunge per pass multiplies material entries at tight caps.
  - :239 generation 7784 ms on 12x8 at cap 40.
- **unused capability / ideas:** condition 2 "does this loop still cut anything" by exact point location (:93-94) — a cheap productive-loop test, related to the draft plan's productive-only emission (which uses a full subtract + equality instead).
- **claims:** :21-22 "worst accepted circle within 3.4° of the cap at every cap 40-100" — unverified. :240 "zero uncleared grid points" on 6x4 and 10x6 — but machining_quality.md:96 reports uncut 0.008041 on rect_12x8 for this generator at 120°; different pockets, not contradictory. :198 "the radius knob is spent" — consistent with draft plan focusing on order/entry/corners.

### docs/refutation_soundness_crosscheck.md | lines 382 | read 1-382
- **role:** doc (validation evidence for the adaptive generator's link refutation)
- **used by:** mkdocs; linked from continuous_engagement.md
- **state:** live (adaptive lane parked)
- **generation:** n/a directly. The adaptive generator's link search rejects over-cap links via an exact station probe (160× cheaper than the partition, :8-10).
- **defects:** :213-264 partition incompleteness throw site `src/continuous_tea_2/segment_fibre.cpp:958` (`k` roots vs `m` identities, split/merge events uncovered); :279-331 process abort is `cpp_int` → double overflow in `float_next` after coordinates grow to 139 digits on a 135-vertex arrangement; 6.3 GB RSS. Both open.
- **unused capability / ideas:** :271-277 two cheap next steps (payload in exception message; correspondence rule for `k != m, k > 1`) — not done (unverified).
- **claims:** :21 "0 REFUTED-and-certified over 125 audits" — unverified (measurement). :333-361 retraction of the unsoundness claim — consistent with continuous_engagement.md:935-941.

### docs/segment_site_mat.md | lines 4820 | read 1-4820
- **role:** doc (exact segment-site MAT + certified adaptive pipeline stage page; also Held Figure 5/8 lane records)
- **used by:** mkdocs; named stage page in CONVENTIONS.md:150-155; documents `src/segment_site_*`, `src/medial_axis_2.cpp`, `src/compas_cgal/adaptive/*`, `benchmarks/held_*`
- **state:** parked (adaptive lane) / live (Held lane sections)
- **generation:**
  - **Exact MAT scope is two fixtures.** Production graphs exist only for an exact rectangle (`canonical_rectangle_mat_graph`, rejects non-orthogonal 4-edge rings, :4757-4758) and the exact six-edge L ring `(0,0)-(6,0)-(6,2)-(2,2)-(2,6)-(0,6)` (:3516-3524, :3733-3734 "arbitrary six-edge rings and general pockets are not implied"; verified `src/segment_site_catalog_graph.cpp:1207-1220` calls `require_l_shape_input`). The Python 20-field tuple and planner are gated for the production L only (:4381-4382). "General degeneracy-removal traversal: pending" (:4359). **So the exact segment-site MAT cannot replace the straight skeleton for the gate's `L_shape` or arbitrary pockets today**; whether the rectangle path accepts `rect_12x8`/`rect_20x12` dimensions is unverified.
  - MATHSM candidate lattice (:581-619): circle centre `c = midpoint(m, q)` where q is the tool-contact point on the boundary footpoint ray, `rho_max = (|m-p|-r)/2`; ordered by furthest progress, then largest radius — the one-sided Held/MATHSM construction, built for the adaptive generator only.
  - **Held lane (live):** corner-aware refinement (:1619-1680): limits radii at convex corners, inserts **shared-contact fans at concave corners**, subdivides by physical boundary distance; Figure 5 2,013 circles, 0 successors over 80°. Native contour-aware critical point `HeldDiskContour2` (:1370-1420) — Held §3.1 exposed-contact construction — built; contour-aware *placement* (max engagement over candidate against evolving contour) not built (:1417-1420, :1453-1456). Boundary-order spacing (:1460-1501) −58.7% path length. Qualified contacts and **corner approaches** (:1586-1593).
  - Precleared entry `QualifiedBore`/`PreclearedEntry` (:1025-1055): a declared pre-drilled bore of radius > tool radius as a physical input; first circle must lie inside the bore void.
- **defects:** :1324-1355 Figure 5 hypothesis generator places circles in lanes then reorders them, so predecessor checks don't hold on the emitted path; 801-1,010 exact containment rejections (:1637, :1662-1663) are unresolved (near-tangent float-built circles). :1347 earlier 1,532-circle path had 108 over-cap successors. :2101-2105 smoothing a certified arc/line path to higher order does not preserve the engagement bound (open problem). :4090-4093 rectangle profiles/neck evidence pending.
- **unused capability / ideas:**
  - **pre-drilled / precleared entry** as typed input (built in `adaptive/entry.py`) — an alternative to draft-plan D1's helix that needs no new replay vocabulary for Z+XY motion: a vertical bore of radius > r, then the first loop inside it. Not mentioned by the plan.
  - **concave-corner shared-contact fans and convex-corner radius limits** (Held lane) — corner strategy built, unused by product generators.
  - `HeldDiskContour2` exposed-contact query — the contour-aware placement primitive.
  - zero-guide advancing segments on MAT arms of clearance = r (see continuous_engagement.md record).
  - :2094-2099 retain Held's structural efficiency (backwards contour scan, graph-linear bottleneck search).
- **claims:** :1484 "58.7% reduction" — consistent with held_performance.md:15-16 (44,186.9 → 18,236.8 mm; verified arithmetic 1-18236.8/44186.9 = 0.587). :2059-2066 hard target 1/10 of Held's planner clocks — no evidence of being met (unverified/likely open). :1730-1735 reconstruction census Figure 8 upper bound `0.08131110669521248` differs from task-1-report:159 `0.08300467787418957` — contradiction between docs (unresolved which is current). :3-9 "Monstera has two sampled residual candidates" — consistent with held_motion_coverage.md:714-720.

### docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md | lines 797 | read 1-797
- **role:** doc (approved Phase-1 design spec for the certified adaptive generator)
- **used by:** no consumer (design record); implemented in part by `src/compas_cgal/adaptive/*`
- **state:** parked (implementation stalled at inter-route transit / terminal coverage, continuous_engagement.md:1339-1342)
- **generation:** defines a world-class-intent generator: contour-aware MATHSM on the true segment-site MAT (:193-252), exact reachable material `M_r = (D ⊖ B_r) ⊕ B_r` with the unreachable corner residual `D \ M_r` reported separately and **not** counted as a coverage failure (:81-86, :260-289), finite candidate lattice, joint link+circle transactions (:514-546), neck policy with oriented passage caps (:548-577), terminal coverage `M_r \ ∪W_i = ∅` (:579-600), `max_passes` truncation forbidden (:583).
- **defects:** :59-61 "every nonzero displacement from [a tool-radius plunge] engages more than π" against frozen stock — a plunge alone cannot launch a capped lateral move; Phase 1 requires a precleared bore. :381-385 **"The existing center-distance-only growth lemma is not a Phase-1 certificate: equal-radius circle coincidence is reachable after a tool disk is removed ... The current overlap visitor also drops positive-length coincident arcs. Those facts invalidate a universal endpoint-drift bridge."** This contradicts the soundness of `certify_segment_tea`'s growth bound (sp1 plan record) and is the same mechanism as the translation-invariance failure (machining_metric_validity.md:289-333).
- **unused capability / ideas:**
  - `M_r` reachable-material definition (:81-86): exactly the principled answer to draft-plan D2(a); built as `ReachableDomain` (`src/compas_cgal/adaptive/reachable_domain.py`, verified `ls`) but not used by the product gate's uncut metric (which uses a grid over the "tool-reachable region", machining_quality.md:39).
  - Phase 1.5: generated/qualified helical entry; Held's sorted circular-arc machined contour (:777-778). Phase 2: **paper-style offset-contour transitions and elliptic returns** (:779) — never built.
  - ablation requirement (:755-759): fixed-policy comparator must show a nonzero certified-violation tail.
- **claims:** :3 "approved architecture" — consistent; :91-93 "No path is returned when any guarantee is unresolved" — consistent with adaptive generator never returning a complete path.

### docs/superpowers/specs/2026-07-25-task3-single-arrangement-performance-design.md | lines 430 | read 1-430
- **role:** doc (design: exact reachable domain `D, C_r, M_r, D\M_r` and incremental coverage)
- **used by:** no consumer (design record); implemented as `ReachableDomain2`/`Coverage2` (`src/coverage_2.cpp`, `src/reachable_*`; `reachable_domain.py` verified present)
- **state:** live (built; used by the gate's `survey_path` centre-domain per oblique_edge_cost.md:57-60, and by the adaptive lane)
- **generation:** generic exact construction for any polygon with holes (not fixture-bound): `C_r = closure(D \ boundary_neighborhood_r)`, `M_r` by one batched dilation, `D \ M_r` once (:82-106, :260-278); incremental coverage `R_{i+1} = R_i \ W_i` (:292-321). Restriction: exactly one connected `C_r` component, else `PocketNotMachinableError` (:193-195).
- **defects:** :193-195 single-component restriction means a pocket whose admissible-centre domain splits at a neck narrower than the tool (the draft plan's dumbbell, Review Focus 2) cannot be represented by this owner. :417-419 isolated processes run 2-3 min under containment.
- **unused capability / ideas:** `D \ M_r` is the exact "corner material no tool can reach"; with `M_r` the gate could score uncut over reachable material exactly instead of on a grid, and the draft plan's D2(a) "rest material" ruling has a ready exact definition here. Also `C_r` is the exact plunge-legality/centre-legality set.
- **claims:** :71-73 "no O(kn) post-construction rematching" — unverified.

### docs/superpowers/specs/2026-07-26-segment-event-substrate-design.md | lines 86 | read 1-86
- **role:** doc (design: exact segment event partition, Task 5 of adaptive lane)
- **used by:** no consumer (design record); implemented in `src/continuous_tea_2/segment_*`
- **state:** parked (adaptive lane)
- **generation:** n/a (certification substrate)
- **defects:** none in the doc; implementation incompleteness recorded elsewhere (refutation_soundness_crosscheck.md:213-264).
- **unused capability / ideas:** none for product generation.
- **claims:** :56-58 "no soft unresolved construction path" — consistent with continuous_engagement.md:705-706.

### docs/superpowers/specs/2026-07-28-atomic-candidate-transaction-design.md | lines 119 | read 1-119
- **role:** doc (design: atomic link+circle candidate transaction, Task 12)
- **used by:** no consumer (design record); implemented in `src/compas_cgal/adaptive/transaction.py`
- **state:** parked (adaptive lane; implemented per continuous_engagement.md:453-491)
- **generation:** joint link+circle trial on forks; certify-before-deplete; winner = furthest progress, largest radius, canonical identity (:95-101). A generation pattern the draft plan partly mirrors (productive-only emission tests on the replay's stock) but without the fork/commit or cap check on the link.
- **defects:** none found
- **unused capability / ideas:** fork-evaluate-commit candidate transactions with cap-checked links — the draft plan emits bridges/plunges without a cap check on the link.
- **claims:** consistent with continuous_engagement.md:453-491.

### docs/superpowers/specs/2026-07-28-causal-mat-traversal-design.md | lines 343 | read 1-343
- **role:** doc (design: Task 13 global MAT traversal, entry launch, terminal coverage)
- **used by:** no consumer (design record); implemented in `adaptive/traversal.py`, `traversal_graph.py`, `bootstrap.py`, `generator.py`
- **state:** parked (implemented through route 2 of the L fixture; terminal coverage open)
- **generation:** the route rule (:155-167) is **the draft plan's tree walk**: seed the component containing the entry, discover unfinished edges through shared node IDs, orient each edge from its discovery node, keep both incidences so cycle-closing edges are not dropped, finish a component before the next. The draft plan re-specifies the same on the straight skeleton, rooted at the widest vertex.
- **defects:** :20-25 production MAT backend is the exact L-pocket only. Measured consequence (continuous_engagement.md:1309-1319): the DFS discovery order is not a physical walk — after finishing a leaf, the direct link to the next route gouges; a transit through cleared lineage (retrace) is needed. **The draft plan inherits this exact problem** (each chain after the first "starts in stock an earlier walk already cleared", but the tool is at the end of the previous leaf) and answers it with retract → clearance LINK → plunge; it does not cite this evidence or the alternative (retrace / stay-down through cleared lineage, built as `adaptive/retrace_*.py`).
- **unused capability / ideas:** `TraversalPolicy.order_branches()` (:162) — branch ordering policy hook; finite-search failure taxonomy (:245-275) `EngagementCapInfeasibleError` / `NeckTooTightError` / `NoFeasibleCandidateError` — named exhaustion the draft plan's `virgin_plunges` counter could mirror.
- **claims:** :326-340 acceptance criteria — not met on live (continuous_engagement.md:1292-1342).

### docs/superpowers/specs/2026-07-29-one-root-full-circle-source-design.md | lines 254 | read 1-254
- **role:** doc (design: one-root full-circle event source v4)
- **used by:** no consumer (design record); implemented per continuous_engagement.md:194-232, :1014-1021
- **state:** parked (adaptive lane; implemented)
- **generation:** n/a (certification substrate). :182-194 notes that a **diagonal disk-chain depletion produced 309 boundary records** — depletion representation drives certification cost; the capsule quad later reduced bridge depletion to 10 vertices (engagement_controlled_toolpath.md:228).
- **defects:** :67-69 distinct nonzero radicands per vertex unsupported (`UnsupportedAlgebraicVertexProjectionError`).
- **unused capability / ideas:** none for product generation.
- **claims:** :231-233 233 active pairs — consistent with continuous_engagement.md:349-353.

### docs/superpowers/specs/2026-07-29-zero-guide-link-candidate-design.md | lines 384 | read 1-384
- **role:** doc (design: link-only advance on MAT edges with clearance ≡ tool radius)
- **used by:** no consumer (design record); implemented (`adaptive/advancing_segment_trial.py`, `candidates.py`; continuous_engagement.md:1058-1290)
- **state:** parked (adaptive lane; implemented through replay)
- **generation:** where the MATHSM guide radius is identically zero (exact polynomial identity `C_e(t) − r² ≡ 0`, :79-94), emit a straight advancing segment (`AdvanceSegmentOperation`) instead of a degenerate circle (:24-46); rejects zero-radius circles, traversal-only terminalization and tolerance classification (:52-77). This is the principled form of what the gated generators get wrong at corners ("degenerate loops" with ρ≈0.02) — but it applies only where clearance is *exactly* r along a whole edge (a slot of width exactly 2r), not to corner spokes where clearance merely tends to r at the apex.
- **defects:** applicability is narrow (constant-clearance edges only); its TEA proof is the swept-prefix π theorem, which certifies π, not a 120° or 60° cap.
- **unused capability / ideas:** a "link-only advance" operation type distinct from a hold link — the draft plan's corner run-ins (B2a) and stay-down links (D3) need exactly this distinction (an advancing cut that owns material removal vs. a reposition), and the gate currently cannot tell them apart (brief: a LINK at depth that removes nothing is counted redundant).
- **claims:** consistent with continuous_engagement.md Task 3-7 records.

### docs/superpowers/specs/2026-08-06-exact-inter-route-retrace-design.md | lines 491 | read 1-491
- **role:** doc (design: exact cut-depth retrace between MAT routes)
- **used by:** no consumer (design record); implemented as `src/compas_cgal/adaptive/retrace_segment_trial.py`, `retrace_transaction.py` (commits `988a2a51`, `098aa942`, `2452fb1f`, verified `git log` and ancestor of live)
- **state:** superseded — `docs/red_manifest.json:21-45` lists `test_task13f_full_continuation`, `test_real_active_family_stops_at_unresolved_exact_event` and two `test_route_retrace_generator` tests as expected-red with reason "plan superseded (backlog D1); retirement scheduled through the auditor programme" (verified)
- **generation:** **stay-down transit between disconnected routes by exactly reversing the previous accepted advancing cut** at cut depth (:36-45, :158-163), certified with the swept-prefix π theorem in 0.04 s (:84-89) — i.e. a stay-down link through cleared stock, the thing the draft plan defers to D3. It explicitly rejected "retract, rapid traverse and re-plunge" as a fallback because it "adds clearance-height, collision, approach, plunge and machine-kinematic contracts" (:78-82); the draft plan adopts exactly that retract/plunge strategy for Stage A.
- **defects:** :51-55 bounded to one immediately preceding no-neck advancing segment; general multi-operation retrace not designed. Implementation's generator tests are red/superseded.
- **unused capability / ideas:** retrace-through-accepted-lineage as a stay-down link with a sound, fast π certificate; operation-grammar distinction of retrace (phase changes, cursor holds) — needed for the gate to stop counting stay-down links as redundant operations (draft plan D3).
- **claims:** :86-89 generic oracle >10 min vs swept-prefix 0.040943 s — consistent with continuous_engagement.md:1329-1337; unverified now.

### docs/superpowers/specs/2026-08-23-auditor-convergence-design.md | lines 836 | read 1-836
- **role:** doc (approved programme design P0-P4: truthful audit, replay theorem, generator qualification, enforcement)
- **used by:** no consumer (design record); P1 partly implemented (`src/compas_cgal/engagement_audit/`, `src/audit_*_2.cpp`, ancestor commit `7cda8f6`)
- **state:** parked (draft plan parking lot "auditor-convergence P0-P4")
- **generation:** Stage 3 "generator qualification and typed overrides" (:641-695) is generation work that was never started: (a) an `EngagementOverrideLedger` so forced entries/loops/advances survive the public result (:643-657) — today they are warnings only (:645, :60); (b) **one shared ladder-search module; non-monotonic engagement forbids bisection as admissibility authority** (:659-666) — `engagement_controlled_toolpath` still bisects; (c) typed units at generator boundaries (:668-671); (d) tangent tests must count comparisons and validate every bridge exemption geometrically (:673-678); (e) "Unsupported true helical motion fails loud until a typed 3D depletion contract exists" (:682-683).
- **defects (recorded at 2026-08-23):** :50 unmeasured operations represented by `cap_certified=True` (legacy audit); :53 **swept-prefix certifier certifies from two booleans without a geometric acceptance test**; :599-615 the fast swept-prefix path "is research-only" until a reviewed proof note and adversarial test package pass; production retrace must route through the generic certifier meanwhile. **This removes the swept-prefix theorem as a usable certificate for the draft plan's run-ins or stay-down links** until P2 lands. :60 generator override decisions disappear at public result boundaries.
- **unused capability / ideas:** override ledger; shared ladder search; independent accept-side falsifier (:621-627); evidence bundle (:749-766).
- **claims:** :36-40 branch tips `073a0f7` / `73d5372` — `073a0f7` verified ancestor of live; `73d5372` not an ancestor (verified). 

### docs/superpowers/specs/2026-08-29-wave1-task6-historical-assertions.md | lines 111 | read 1-111
- **role:** doc (design: measurement-claim history artifact)
- **used by:** no consumer
- **state:** superseded-by removal of the measurement-claim subsystem (`6b4293ae refactor: drop measurement-claim subsystem`, `67d636f3`, verified log; measurement_claims.md:8-11)
- **generation:** n/a. Records that historical generator claims (MC-001..MC-010 from `engagement_radial_toolpath.py`/`engagement_toolpath.py`) were removed in `53135e0`; the 8-vs-12 over-cap conflict (:64-71) is the same as radius_ladder_ablation.md:19-23.
- **defects:** :92 references `docs/superpowers/state/wave1-gate-analysis.md` — directory absent on live (dangling).
- **unused capability / ideas:** none
- **claims:** :98-111 completion conditions — moot after subsystem removal.

### docs/superpowers/specs/2026-08-31-held-pfeiffer-reference-pockets-design.md | lines 309 | read 1-309
- **role:** doc (approved design: Held reference pocket corpus)
- **used by:** no consumer (design record); implemented (`benchmarks/held_reference_*`, held_pfeiffer_reference_pockets.md)
- **state:** live (implemented; task-1-report records completion)
- **generation:** n/a (inputs). Supplies a narrow-neck pocket (crossed skis, :57) and high-complexity Monstera (:58) at an 80° cap — the product gate's three axis-parallel polygons have neither curved walls nor necks.
- **defects:** :48 "The publisher PDF is read locally and is not committed" — contradicted: held_pfeiffer_reference_pockets.md:11-14 says the PDF is stored in the repository (`assets/papers/held-pfeiffer-2025.pdf`). :233-242 qualification is non-empty-path only.
- **unused capability / ideas:** use the four Held pockets (with their 80° cap) as additional product-gate pockets — the draft plan's Stage C proposes running the product generator on them, but not gating on them.
- **claims:** :30-31 "Replace the current rectangular stand-in" — partially: Figure 6 uses the Figure 5 projection (held_pfeiffer_reference_pockets.md:64-66); the product gate still uses rectangles.

### docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md | lines 814 | read 1-814
- **role:** doc (approved design: score a generated Figure 5 path on the 12 gate criteria + guarded audit; post-qualification seam)
- **used by:** no consumer (design record); implemented in part: `benchmarks/quality_observations.py`, `held_path_*`, `CutterCentreDomain2`, `ReachableMaterialPredicate2` (exactness.md:654 references the latter)
- **state:** open (Phase 1 report `docs/benchmarks/held_figure5_2d_path.md` is absent on live, verified `ls`; Phase 2 generator-policy design never written)
- **generation:** characterizes "the current default generator" = `benchmarks.runner.generate_toolpath` → `trochoidal_mat_toolpath_circular` (verified `benchmarks/runner.py:33-45`), the unregulated stepover generator, on Figure 5 at 80°. :805-810 "Its measured failure signature becomes the sole input to the Phase 2 generator-policy design. ... No generator repair begins in the Phase 1 plan." Phase 2 never happened.
- **defects:** :492-500 first live run stalled >30 min in `survey` building the full `ReachableDomain2`; :538-606 coverage reduction exceeded 600 s through three revisions; revision 4's two-locator predicate ran 305,944 queries in 33.8 s (:646-651). The gate's coverage needs a 668×458 grid on Figure 5 (:538-547) — the 200-sample default is too coarse for curved pockets.
- **unused capability / ideas:**
  - this is **the scoring adapter the draft plan's Task 0b needs** (a Held path scored on the gate's 12 criteria): `quality_observations.py` reducers with operation attribution, `PathQualityAssessment`, typed evidence kinds (:114-133, :246-308). The draft plan says "no adapter from a Held path to ToolpathResult exists" — this design built the reducer/attribution side; the remaining gap is converting a Held-lane circle path into a `ToolpathResult`. Unverified whether `held_path_*` modules exist on live.
  - post-qualification seam and a target-specific (Fanuc) postprocessor phase (:680-709, :812-814) — never built; nothing in the repo emits G-code.
- **claims:** :18-22 "Figure 5 currently produces 2,289 operations" — consistent with held_reference_qualification.md:7.

### docs/superpowers/specs/2026-09-08-number-type-coherence-design.md | lines 403 | read 1-403
- **role:** doc (approved design: exact number vocabulary, stages 0-6)
- **used by:** number_types.md; stages 0-2 implemented (`src/exact/`, commits `9161bc84`, `541b6ede`, `003e144c` verified in log)
- **state:** parked (stages 3-6 in draft plan parking lot)
- **generation:** n/a (no generator uses `continuous_tea_2`/`segment_site_*` lanes; the gated generators run on Epick `toolpath.cpp` + Epeck `Stock2`/`engagement_2`, both already filtered lanes per number_types.md:26-28)
- **defects:** :336-357 invariant 1 violated by dead second door `exact_binary64`/`lift_binary64` (zero callers) — scheduled for stage 3 (parked). :359-389 two incompatible framings under tag `exact-rational-v1`; `CanonicalRational::canonical_bytes()` latent defect.
- **unused capability / ideas:** none for generation.
- **claims:** :17 "design approved, not implemented" — stale; stages 0-2 landed (number_types.md:677-697).

### docs/public/favicon.ico | lines 16 | read 1-16
- **role:** build (mkdocs asset; binary 32×32 PNG, 4,856 bytes, inspected with `file(1)`, not text)
- **used by:** mkdocs theme (unverified)
- **state:** live
- **generation:** n/a
- **defects:** none found (a PNG stored under a `.ico` name)
- **unused capability / ideas:** none
- **claims:** none

### docs/toolpath_drawing.md | lines 140 | read 1-140
- **role:** doc (figure API `benchmarks.plotting`)
- **used by:** mkdocs; `benchmarks/figures.py`, pixi `figures`
- **state:** live
- **generation:** n/a. Notes "lead moves ... the current generators emit none at all" (:110-111, :119) — no lead-in/lead-out on any engagement generator.
- **defects:** none for generation.
- **unused capability / ideas:** `ColourBy.ENGAGEMENT` per-operation plan view and `annotate_traversals` — useful to inspect the tree-walk order.
- **claims:** :12-14 numbers measured at draw time — unverified.

### docs/toolpath_tangent_continuity.md | lines 198 | read 1-198
- **role:** doc (tangent-continuity model of `toolpath.cpp` and its test)
- **used by:** mkdocs; `tests/test_toolpath.py` tangent tests
- **state:** live, partly stale
- **generation:** consistent winding: full circle + short repositioning arc at nominal winding at varying-radius stations (:53-62) — verified current (`src/toolpath.cpp:394-401`). Repositioning-arc gate `compare_squared_distance` against `1e-18` (:83-90).
- **defects:**
  - :33-51 says bridge lines are "~90° ... inherent to trochoidal milling geometry — not a defect" and the test exempts bridge↔circle transitions (:167-171). The product gate counts exactly those as tangent breaks (24/32/22, machining_quality.md:126), and loop_radius_degeneracy.md:151-190 shows the break is the clearance-slope tilt, fixable by the external common tangent. The exemption hides a real G1 defect.
  - :83-90 `1e-18` threshold "treated as coincident" is an epsilon on a geometric decision (merge vs reposition) — violates the CONVENTIONS.md Epick rule against epsilon coincidence, though framed as an exact predicate on a threshold.
  - :152-165 test checks parallelism with `abs(dot)`, so a 180° reversal passes.
- **unused capability / ideas:** none
- **claims:** :35 "~90° angle is inherent ... not a defect" — contradicted by machining_quality.md:126 and loop_radius_degeneracy.md:171-190.

### docs/toolpath_trochoid_chain_design.md | lines 114 | read 1-114
- **role:** doc (circle stitching at varying radii in `toolpath.cpp`)
- **used by:** mkdocs
- **state:** superseded-by full circle + repositioning arc (toolpath_tangent_continuity.md:53-62; `src/toolpath.cpp:394-401`)
- **generation:** records four stitching options; claims the "current implementation" is the adaptive-winding long-way arc that flips winding when the nominal sweep < π (:74-97).
- **defects:** :74 "current implementation" — contradicted by `src/toolpath.cpp:394-401` (full circle + repositioning arc, nominal winding) and toolpath_tangent_continuity.md:77-81 ("~60% of arcs were flipped ... The fix splits this"). :95-97 winding flips mix climb/conventional — a known machining defect the later fix removed.
- **unused capability / ideas:** none
- **claims:** :112 "≤1.5% short for varying radii" — moot.


# Off-branch documents (not ledger-checked)

### codex/sdd-coherence:docs/superpowers/state/backlog.md | lines 158 | read 1-158
- **role:** doc (coherence programme backlog; "the one place open work lives")
- **used by:** none on live (`docs/superpowers/state/` absent on live)
- **state:** parked/off-branch (the draft plan, dated 2026-10-13, is now "the one active plan")
- **generation:** :34-36 the 12 quality-gate reds are the product gate. :94-98 **parked lane C4** "depletion model: score a true helical entry; until then the entry criterion is unsatisfiable — 240° floor and `V ≥ √3·r` proved in `loop_radius_degeneracy.md`" and **C5** "arc motion vocabulary for tangent-continuous inter-chain links" — both generation items, never built.
- **defects:** :96 cites `V ≥ √3·r` as proved in loop_radius_degeneracy.md — no such statement in the live page (verified grep `sqrt(3)|√3`: no match); contradicted/unverified. :116-123 names `codex/auditor-convergence-sdd` as canonical frontier and `jf/toolpath-redesign@73d5372` frozen — superseded by the live branch `codex/held-task5-continuation`.
- **unused capability / ideas:** C5 arc-link vocabulary (G1 inter-chain links) — the draft plan uses straight LINKs at clearance and straight bridges; tangent-continuous arc links are absent from its design.
- **claims:** :112 "Anything red beyond these 17 is a defect" — consistent with live `docs/red_manifest.json` (6+6+1+4 entries, verified).

### jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md | lines 194 | read 1-194
- **role:** doc (ETH-grade audit spec of `toolpath`, `stock_2`, `engagement_2`; commit `5c1d492a`, not an ancestor of live)
- **used by:** none on live (file absent on live)
- **state:** off-branch; findings checked against live below
- **generation:** C1 (:35-61) **`mat_scale > 1` silently emits gouging toolpaths** — trochoid circles are never certified, safe only by construction; 1.5 → +2.2485 penetration. I5 (:149) `link_paths=False` concatenates every op into one polyline, a through-material traverse. I7 (:151) gouge decisions in `toolpath.cpp` are on Epick rounded doubles. I8 (:152) `certify_bridges` can exit with an uncertified bridge. I9 (:153) hole validation incomplete (UB into straight skeleton). I1 (:145) Python arc certifier issues 0/126 positive certificates in the generator's regime. C2 (:63-129) `tea_growth_bound` is not an upper bound near concave voids (ρ→r⁺), safe on the shipped path only because cut radius = query radius.
- **defects (live status):** C1 **open** — `validate_toolpath_params` (`src/toolpath.cpp:607-613`) still does not validate `mat_scale`; `toolpath.py:282,379` unchanged (verified grep). C3 **open** — `pyproject.toml:9 requires-python = ">=3.9"` (verified). `_unmeasured` returns `cap_certified=True` — **open** (verified `src/compas_cgal/engagement.py:320-322`). I1 — `AUDIT_ARC_STEP_FRACTION = 0.05` unchanged (`engagement.py:48`). I6 shared-root guard — CHANGELOG/brief commit `148a49ad fix(engagement): enforce the shared-root precondition in release builds` is on jf/toolpath-redesign; whether on live unverified.
- **unused capability / ideas:** curvature-aware growth bound (:127-129); `mat_scale` contract.
- **claims:** :18 "exact-kernel core is publication-grade" — opinion.

### jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md | lines 3472 | read 1-3472
- **role:** doc (remediation plan for the ETH audit + machining-quality track Tasks 29-37; created in `5c1d492a`, extended in `c90a4fef`)
- **used by:** none on live (absent on live)
- **state:** off-branch; Tasks 0-9, 22, 23 recorded complete on jf/toolpath-redesign (:67), most not ancestors of live; Tasks 10-37 open
- **generation:** the most concrete generation backlog in the slice:
  - **Task 30** (:3065-3095) register `rho_regulated_toolpath`, `chain_ordered_toolpath`, `spiral_entry_toolpath` in the gate; recorded claims on rect_12x8 @120°: degenerate loops 4→0, tangent breaks 24→0 (rho); **material entries 5→1 (chain ordering)**; worst entry 243.20° vs 240° floor (spiral). Live status: **open** — `benchmarks/gate.py` registers only `engagement_controlled_toolpath` and `radius_regulated_toolpath` (verified grep :38-39, :176, :194); chain-ordered generator was later deleted (`9fa6eea8`).
  - **Task 31** (:3099-3135) consume `declined_regions`: a corner-clearing pass (non-trochoidal motion) or a declared second operation with a smaller tool; split `uncut_fraction` into declined-and-named vs silently-missed. Measured rho uncut regression 0.0080→0.0295 / 0.0028→0.0117 / 0.0069→0.0352. Live status: **open** — `declined_regions` is consumed nowhere outside its own generators/tests (verified grep).
  - **Task 33** (:3164-3194) missing KPI families: cutting force (calibrated), **exit impact** (exits at high engagement; uncalibrated), controller starvation (block length vs interpolation period; "an early radius sweep produced 862 cuts at one cap").
  - **Task 34** (:3198-3233) generalise the segment-site MAT; **warning (:3215-3225): measured head-to-head on the fixture, the middle curve is NOT gentler than the straight skeleton — guide-radius CoV 0.5719 vs 0.5752, degenerate loops 24.5% vs 13.0%** (middle curve degenerates at `d ≤ 3r`). The case for the true MAT is structural and untested.
  - **Task 36** (:3266-3347) helical / ramped / pre-bored entry; both helical and ramped break `_infer_cut_height`'s single-plane assumption and the audit's ramp refusal; property: no motion exceeds the cap, entry included.
  - **Task 37** (:3351-3434) break the ~141.8° saturation; Step 1 must decide whether the station set (guide) or the regulator is binding; success = length curve steepens (Held 40-60×, ours 6.8×).
  - Task 29 (:3022-3061) second gate cap — **closed** on live (gate has cap-120 and cap-40 cells, verified `docs/red_manifest.json:4-13`).
  - Task 1/2/4/13/14/16/17 generator hygiene (mat_scale, implied traverse, parameter contracts, bridge certification exit, hole validation, sentinel, kernel external tangent) — open on live for Task 1 at least (verified `src/toolpath.cpp:607-613`).
- **defects:** :117-134 the plan itself notes "Two Held-parity items still have no task" (later added as 36/37). Task 7's curvature-aware bound and Task 10's arc certifier landed only on jf/toolpath-redesign (brief: `73d53729 feat(engagement): adaptive arc certification` is that branch's head; not an ancestor of live).
- **unused capability / ideas:** pre-bored entry "start the chain where material is already gone" (:3300) — the draft plan's tree walk is exactly this for chains after the first; exit-impact KPI; declined-vs-missed uncut split.
- **claims:** :3029-3033 at cap 120 both gate generators identical (cutting length 411.566 / 1018.665 / 565.311) — consistent with machining_quality.md:122 (411.6 / 1018.7 / 565.3).

### jf/toolpath-redesign:INTEGRATECODEX.md | lines 241 | read 1-241
- **role:** doc (integration decision memo, commit `f369b3df`)
- **used by:** none on live
- **state:** superseded (the auditor-convergence reconciliation and later live branch replaced the "rebase jf onto codex" plan; backlog.md:116-123)
- **generation:** "critical path to Held parity" (:136-169), ranked: **1 break the ~141.8° saturation** (station set suspected); **2 helical or ramped entry** (clears 2 of 6 gate failures); **3 close the residue** — 0.28-0.80% vs unregulated 0.152%, not corner-confined, "coverage is not in the accept/reject rule, so nothing drives the frontier where the predicate refuses"; **4 kill degenerate and redundant motions**; 5 speed last (~86× unregulated; `engagement_at` 34%, depletion 22%). :77-83 the swept-prefix theorem forces the cap to be **exactly π** for every advancing cut; any other cap makes the first zero-guide advance raise and abort — "The whole adaptive fixture runs at π."
- **defects:** :85-95 three exactness defects (non-finite doubles in five `Stock2` entry points — `subtract_disk(c, inf)` empties stock; shared-root guard compiled out; disk-chain int cast saturating to a ~34 GB reserve) — fixed on jf/toolpath-redesign; live status unknown.
- **unused capability / ideas:** item 3's diagnosis — **coverage must enter the accept/reject rule** (drive the frontier where the predicate refuses) — is absent from the draft plan, which only adds productive-only emission (removes non-productive loops) and does not add coverage pressure.
- **claims:** :36 "three [generators], all working end to end" on codex — at the time; the gate shows they fail 6-7 criteria (machining_quality.md:139). :4-5 codex's cap "strictly stronger than Held & Pfeiffer's float bisection" — only per evaluated position (engagement_controlled_toolpath.md:276-278).

### jf/toolpath-redesign:MOTIVATION.md | lines 16 | read 1-16
- **role:** doc (positioning note, commit `f369b3df`)
- **used by:** none on live
- **state:** off-branch
- **generation:** names the hard parts a world-class generator must handle: run-merges (O(1) engagement jumps), necks/bottlenecks where θ_max must be *reduced*, concave transitions where MRR spikes (:7); G¹/G² smoothness, linking moves, "the ~50% non-cutting air problem" (:13).
- **defects:** none (prose).
- **unused capability / ideas:** bottleneck cap reduction (built only in the adaptive lane's neck policy); air-time reduction as a target.
- **claims:** :9 "3-5× the feed" for constant-engagement commercial tools — unverified, uncited.

### /Users/jelle/Code/CADCAM/compas_cgal_prs/MOTIVATION.md | lines 16 | read 1-16
- **role:** doc (main-checkout copy)
- **used by:** none
- **state:** byte-identical to `jf/toolpath-redesign:MOTIVATION.md` (verified `diff`: no output); see that record.
- **generation:** see jf/toolpath-redesign:MOTIVATION.md record.
- **defects:** none beyond that record.
- **unused capability / ideas:** see above.
- **claims:** see above.

### jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md | lines 567 | read 1-567
- **role:** doc (ETH-level assessment of the codex tip `e0ae116`, commit `f369b3df`)
- **used by:** none on live
- **state:** off-branch; feeds auditor-convergence design
- **generation:** generation findings, checked against live:
  - **`_largest_admissible_advance` exists twice with opposite algorithms** (:401-405): `engagement_toolpath.py` bisects under a monotonicity assumption; `engagement_rho_toolpath.py` scans because that assumption is false; **`engagement_radial_toolpath.py` imports the bisecting one** despite its "never a bisection" doctrine. Live: still true — `engagement_radial_toolpath.py:138` imports from `engagement_toolpath` (bisecting, `:601`); rho and spiral use rho's scan (`engagement_rho_toolpath.py:480`) (verified grep). **The draft plan reuses rho's station rules ("ladder scan") — correct; it should not reuse `engagement_toolpath`'s advance search.**
  - C9 (:318-337) forced loops/advances survive only as warnings; `RadialToolpathResult` drops `forced_advances`.
  - :364-372 Hypothesis strategy generates only fat star-shaped polygons — no necks, pinches, islands.
  - :398-400 the corpus exercises the *unregulated* generator; `engagement_rho_toolpath.py` "is wired into nothing".
  - :343-363 tangent-continuity test checks 7.7% of transitions (structural bridge exemption).
  - :486-493 the 329.49° step and 89.7% recut are the generator's behaviour.
- **defects:** C2 (:231-254) arc certifier guard exceeds a full turn past 2× tool radius; uncertified 47-52% on every corpus instance (live: `AUDIT_ARC_STEP_FRACTION` still 0.05, `engagement.py:48`). C1 (:189-229) vacuous certification (live: `_unmeasured` still `cap_certified=True`). C7 (:294-308) swept-prefix verdict from two booleans, no geometric test. C4/C5/C6 CI never ran; 6 mypy errors; Python floor false.
- **unused capability / ideas:** :502-506 a separate `pixi run gate` outside `baseline` — since implemented as the red manifest (live `docs/red_manifest.json`).
- **claims:** :61-62 "generation measured at 72.4 ms ... inside Held's stated 3-100 ms range" — only on axis-parallel rect_12x8 (oblique_edge_cost.md:13-18); :464 baseline 1288 passed / 10 failed — historical.

### jf/toolpath-redesign:review.md | lines 827 | read 1-827
- **role:** doc (technical due-diligence memo, dated 2026-07-29 with 2026-08-20/22 updates; commit `f369b3df`)
- **used by:** none on live; backlog W3.3 says it must be re-issued (predates the oblique-edge finding)
- **state:** off-branch, partly stale
- **generation:** :460-491 **loop engagement is 171.8° at both a 40° and a 120° cap, bit-identical; "engagement runs from 0° at a small radius to 114.5° at the radius the generator actually picks; it is maximising the loop rather than regulating it."** Real work wants 30-60° caps (:480-482), "exactly the settings the generator currently ignores". :343-357 the generator "plunges into pocket corners and cuts straight out of them": four degenerate loops ρ 0.02-0.06, and "a 1.96-long straight link slotting into a corner at 106.5°". :761 air cutting (TEA = 0 at cutting feed) as a verdict with direct cycle-time value.
- **defects:** :352-353 "Dropping all four degenerate loops changes coverage by **0.000%**: they remove nothing" — **contradicted** by machining_quality.md:195-197 (removing them raises residue 0.284% → 0.351%). :432-437 "generator parity ... 71.9 ms ... inside Held & Pfeiffer's published 3-100 ms range" — axis-parallel only (oblique_edge_cost.md:13-18), and contradicted for the oblique case. :529-534 `const_cast` in `engaged_arcs_zone` blocks parallel station queries.
- **unused capability / ideas:** verdict menu (:753-765): slotting detection, **air-cutting detection**, rapid-into-material, residual, chip thickness, multi-op stock continuity — air-cutting is the gate's `recut_fraction`/redundant-ops idea at motion granularity; cheap conservative screen + exact certifier on flagged set (:792-796) — the right architecture for a generator's inner loop too (refutation probe then certify).
- **claims:** :35-39 "certifier is independently saleable" — opinion; :567 resolution-rate on third-party geometry "never measured" — still unmeasured (no artifact in slice).

### /Users/jelle/Code/CADCAM/compas_cgal_prs/review.md | lines 827 | read 1-827
- **role:** doc (main-checkout copy)
- **used by:** none
- **state:** byte-identical to `jf/toolpath-redesign:review.md` (verified `diff`: no output); see that record.
- **generation:** see jf/toolpath-redesign:review.md record.
- **defects:** see that record.
- **unused capability / ideas:** see that record.
- **claims:** see that record.

### jf/toolpath-redesign:CONVENTIONS.md | lines 153 | read 1-153
- **role:** doc (project rules on jf/toolpath-redesign; changed in `c90a4fef` "pixi build rules")
- **used by:** agents in the main checkout
- **state:** stale — points the pixi manifest at `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-exact-certified-adaptive-phase1-t9-zero-guide` (:126-130), a checkout removed 2026-08-28 (backlog.md:132 D4); the live branch is `codex/held-task5-continuation`
- **generation:** n/a. Lacks the live CONVENTIONS.md's "Single-Path Novel Architecture" and "Development-Stage Documentation" sections and says "12-question review checklist" (:39) vs live "14" and exactness.md "eighteen".
- **defects:** :126-133 dangling worktree path (also present in the user's project CONVENTIONS.md shown to this session).
- **unused capability / ideas:** none
- **claims:** :39 checklist count — contradicted by live docs/exactness.md:677 (18).

### jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md | lines 3152 | read 1-3152
- **role:** doc (implementation plan for `benchmarks/`; commit `c90a4fef`)
- **used by:** none on live (plan absent); implemented on live as `benchmarks/` (docs/benchmarks.md describes the same corpora)
- **state:** implemented (superseded in detail)
- **generation:** two generation-relevant baselines: **MATHSM constant-spacing baseline** (Task 11, :2353-2522) — sweep `stepover` on the unregulated generator, take the shortest compliant trial (no monotonicity assumed); **Figure 6** (Task 12, :2526-2705). As written, `figure6_points` generates "ours" with `trochoidal_mat_toolpath_circular` without passing the cap (:2645) — so "our" path does not depend on the cap and only the audit changes; benchmarks.md later uses the engagement-controlled generator (benchmarks.md:92).
- **defects:** :1471-1595 the "arc-fraction" complexity family tessellates bulged sides into polylines, so no true arc boundary reaches the kernel (the plan claims it closes the "straight-sided" coverage hole, :1400). :664-688 Pythagorean rotations computed as `a / c` in binary64 are not exactly rational (4/5 is not dyadic), contradicting "every rotated coordinate stays exactly representable" (:3100-3102, also docs/benchmarks.md:67-72); consistent with oblique_edge_cost.md:76-79 noting those cases are slow.
- **unused capability / ideas:** necks pinch sweep (`dumbbell`, PINCH_SWEEP_DEFAULT 1.05-4.0 × D) and island grids exist as pocket families but are not product-gate pockets; external third-party profile loader (Task 13).
- **claims:** :3071 "This project's generation is already competitive" — on axis-parallel pockets only (oblique_edge_cost.md). :3145-3146 the resolution-rate and neck-ablation runs — never recorded (no artifact in this slice).

# Findings register

F-d1_docs_specs_reviews-1 | .superpowers/sdd/2026-08-23-auditor-convergence-p1-truthful-audit/task-3-report.md:134-136 | Task 4 owes schema v2, depletion policy in identity, authenticated digest, native certification, all-motion atomic replay | closed (evidence: 7cda8f6 "exact native motion decision" is ancestor of live; src/audit_replay_2.cpp present) | n
F-d1_docs_specs_reviews-2 | .superpowers/sdd/2026-08-23-auditor-convergence-p1-truthful-audit/task-3-report.md:128-132 | exact arc surrogate claims no retained-sliver bound: audit replay may leave slivers the real tool removes | open | n
F-d1_docs_specs_reviews-3 | .superpowers/sdd/2026-08-23-auditor-convergence/progress.md:211 | P4 must identify the extreme full-suite tail before CI enforcement | open | n
F-d1_docs_specs_reviews-4 | .superpowers/sdd/2026-08-23-auditor-convergence/progress.md:153-156 | full-circle closure on machined and spiral arc fixtures exceeds 120 s; no arc-adapter coverage for those | open | n
F-d1_docs_specs_reviews-5 | .superpowers/sdd/2026-08-23-auditor-convergence/progress.md:213-215 | next command: engagement_audit replay/report tests (Task 4D+) never recorded done | unknown | n
F-d1_docs_specs_reviews-6 | .superpowers/sdd/2026-08-31-held-pfeiffer-reference-pockets/task-1-report.md:19-21 | reconstruction/proof recursion has depth ceiling but no shared node budget; decide whether a work budget is needed | open | n
F-d1_docs_specs_reviews-7 | .superpowers/sdd/2026-08-31-held-pfeiffer-reference-pockets/task-1-report.md:128,171-173 | independent certificate re-review pending | unknown | n
F-d1_docs_specs_reviews-8 | .superpowers/sdd/2026-08-31-held-pfeiffer-reference-pockets/task-1-report.md:123-124 | biarc children merge-ineligible until source cubic split at a sound arc-length correspondence | open | n
F-d1_docs_specs_reviews-9 | CHANGELOG.md:21-28 | Unreleased lists mesh/geodesic/polyline APIs that were removed (71d4c7c7) | open | n
F-d1_docs_specs_reviews-10 | CHANGELOG.md:13 | dangling reference docs/superpowers/state/sp1-gate-c-analysis.md; O(n^2) exact depletion in generation "deferred" | open | y
F-d1_docs_specs_reviews-11 | CHANGELOG.md:19 | max_passes truncation only warns (draft plan wants fail-loud) | open | y
F-d1_docs_specs_reviews-12 | CONVENTIONS.md:150-155 | segment_site_mat.md named as current stage page; draft plan (2026-10-13) supersedes without updating CONVENTIONS.md | open | y
F-d1_docs_specs_reviews-13 | docs/benchmarks.md:120-133 | engagement_controlled regulates advance not radius; saturates ~141.8° for caps 20-140; must also regulate trochoid radius (Held takes smaller circles) | open (radius_regulated exists but is also failing the gate; draft-plan Stage C cap 40 deferred) | y
F-d1_docs_specs_reviews-14 | docs/benchmarks.md:13-16 | center_domain ~3,500x slower on oblique edges, mechanism not established | unknown (see oblique_edge_cost.md) | n
F-d1_docs_specs_reviews-15 | docs/benchmarks.md:29-37 | benchmarks/results/ layout no longer enforced after measured-run publisher removal | open | n
F-d1_docs_specs_reviews-16 | docs/benchmarks/held_reference_qualification.md:5-12 | stale primitive counts; generator unnamed; pass = non-empty only, never scored on the gate | open | y
F-d1_docs_specs_reviews-17 | docs/continuous_engagement.md:69-73,1339-1342 | certified adaptive generator incomplete: inter-route retrace, traversal/coverage closure, terminal sealing, arbitrary pockets, Held-matched performance | open (retrace_*.py exist; no gate result) | y
F-d1_docs_specs_reviews-18 | docs/continuous_engagement.md:1309-1327 | MAT route is a DFS discovery ledger, not a physical walk; direct links between routes gouge; need a transit through cleared lineage | open | y
F-d1_docs_specs_reviews-19 | docs/continuous_engagement.md:874-879,914-944 | continuous_tea_2 event discovery incomplete (IncompleteSegmentPartitionError 15/125) and can abort the process (float_next inf) | open | n
F-d1_docs_specs_reviews-20 | docs/continuous_engagement.md:916-933,1355-1364 | test_task13f_full_continuation stale counts; test_real_active_family_stops_at_unresolved_exact_event now fails (proof of violation instead of gap) | open | y
F-d1_docs_specs_reviews-21 | docs/continuous_engagement.md:355-363 | production oracle keeps rational inputs on source v3; v4 one-root migration not validated | open | n
F-d1_docs_specs_reviews-22 | docs/continuous_engagement.md:432-441 | generator/replay must consume MatTraversalState.neck_scope, never a supplied scope | unknown | y
F-d1_docs_specs_reviews-23 | docs/continuous_engagement.md:1105-1128 | IDEA built in adaptive/: swept-prefix theorem certifies straight advancing cuts at exactly pi; unused by gated generators and draft plan | open | y
F-d1_docs_specs_reviews-24 | docs/continuous_engagement.md:587-603 | IDEA built in adaptive/: furthest-progress, largest-radius first-feasible candidate search under exact cap; unused by gated generators | open | y
F-d1_docs_specs_reviews-25 | docs/continuous_engagement_cost.md:112-120 | ship sampled+guarded certifier; continuous_tea_2 is research until univariate formulation question is answered | open | y
F-d1_docs_specs_reviews-26 | docs/continuous_engagement_cost.md:88-91 | Shear_transformation and resultant degree are the only lever at the required magnitude | open | n
F-d1_docs_specs_reviews-27 | docs/dev/arrangement-redesign/compile.sh:6 | prototype build uses removed conda env and CGAL 5.6.1 | open | n
F-d1_docs_specs_reviews-28 | docs/dev/arrangement-redesign/findings.md:63-73 | persistent arrangement + Arr_observer redesign "validated feasible; not built" (superseded by zone-in-Gps fix) | closed (evidence: findings.md:108 LANDED note; superseded) | n
F-d1_docs_specs_reviews-29 | docs/dev/arrangement-redesign/findings.md:67-69 | exact annular-sector sweep (guide_r ± r one-root) deferred; depletion uses disk chains | open | y
F-d1_docs_specs_reviews-30 | docs/dev/cgal-traits-aliasing-report.md:1-9 | Gps copy aliases source traits (use-after-free); filed upstream as CGAL/cgal#9638 | closed locally (evidence: 8d2b76cc, df6fb020, 57482731 on live); upstream open/unknown | n
F-d1_docs_specs_reviews-31 | docs/engagement_audit.md:40,127-128 | truthful audit rejects ramps/mixed-Z motion: a helical or ramped entry (draft plan D1) is inadmissible to this auditor | open | y
F-d1_docs_specs_reviews-32 | docs/engagement_audit.md:217-227 | remaining: replay, native 3-way certification, report aggregation, consumer migration, evidence matrix | open (replay/certification partly closed by 7cda8f6 on live; aggregation/migration unknown) | n
F-d1_docs_specs_reviews-33 | docs/engagement_audit.md:169-173 | exact partial-arc depletion proves no retained-sliver bound | open | n
F-d1_docs_specs_reviews-34 | docs/engagement_controlled_toolpath.md:263-270 | chain entries are full slots; helical or ramped entry not implemented | open (draft plan D1) | y
F-d1_docs_specs_reviews-35 | docs/engagement_controlled_toolpath.md:283-286 | bisection assumes monotone engagement in advance (memory says non-monotone; ladder required) | open (engagement_toolpath.py:42 still bisects) | y
F-d1_docs_specs_reviews-36 | docs/engagement_controlled_toolpath.md:16-23,124-129 | cap guaranteed only at 33 evaluated positions; 43.4 deg measured on a 40 deg cap between positions | open | y
F-d1_docs_specs_reviews-37 | docs/engagement_controlled_toolpath.md:290-291 | guide tangency G1 only on straight guide segments, approximate through turns | open (gate tangent breaks 24/32/22) | y
F-d1_docs_specs_reviews-38 | docs/engagement_controlled_toolpath.md:292-293 | coverage measured on a grid, not certified | open | y
F-d1_docs_specs_reviews-39 | docs/engagement_controlled_toolpath.md:250-261 | no local capsule depletion (seed completeness unsolved) | open | n
F-d1_docs_specs_reviews-40 | docs/exactness.md:432-449 | Manufacturing_policy (linear tolerance, angular surrogate, minimum_feature_size) proposed as SP2 direction; never built | open | y
F-d1_docs_specs_reviews-41 | docs/exactness.md:677 vs CONVENTIONS.md:39-40 | checklist count drift (18 vs 14 vs 12) | open | n
F-d1_docs_specs_reviews-42 | docs/exactness.md:86-92 | says stock uses disk chains for arc sweeps; full turns are now exact annuli (stale) | open | n
F-d1_docs_specs_reviews-43 | src/toolpath.cpp:961-992 (via docs/examples/example_toolpath_trochoidal_mat.py:38-50) | chain ordering is greedy nearest-neighbour from the longest chain, chains never reversed, each with own plunge | open (draft plan Task 2 GuideTree addresses) | y
F-d1_docs_specs_reviews-44 | docs/held_motion_coverage.md:88-103 | restore true medial incidence, qm diameter relation and ordered paired traversal before any local repair layer | open | y
F-d1_docs_specs_reviews-45 | docs/held_motion_coverage.md:251-269 | predecessor-only engagement model cannot certify a convex-corner stationary event; need stock-aware connector/entry | open | y
F-d1_docs_specs_reviews-46 | docs/held_motion_coverage.md:701-707 | shared thinning regression red at three scales (engagement-only selection loses required circle) | open | y
F-d1_docs_specs_reviews-47 | docs/held_motion_coverage.md:3-26 | exact full-motion coverage replay does not finish; fix = one D&C union then single difference, capsules as one polygon, number-type refactor | open (set aside, 334b2f2f) | n
F-d1_docs_specs_reviews-48 | docs/held_motion_coverage.md:511-526 | conditioning audit of other Epeck_with_sqrt lanes; exact representation of emitted circles at coverage seam is a design decision | open | n
F-d1_docs_specs_reviews-49 | docs/held_motion_coverage.md:562-565 | reachable-material target must be reported alongside full-design residuals | open | y
F-d1_docs_specs_reviews-50 | docs/held_motion_coverage.md:714-736 | Monstera interior gap (1.03 mm from motion); selectors omit strongest source circle; local repair needed | open | y
F-d1_docs_specs_reviews-51 | docs/held_motion_coverage.md:528-543 | stale "runtime blocker 337 s" text contradicts fixed 8 ms on same page | open (doc) | n
F-d1_docs_specs_reviews-52 | docs/held_motion_coverage.md:43-58 | IDEA: Held machinable target (opening by r, convex arcs regularized) — never built as generator input | open | y
F-d1_docs_specs_reviews-53 | docs/held_motion_coverage.md:132-160,332-365 | IDEA built (native): contour-aware qm-diameter circle construction; not used by gated generators | open | y
F-d1_docs_specs_reviews-54 | docs/held_performance.md:15-16 | IDEA built in Held lane: boundary-order spacing gives 2.4x fewer circles than lane spacing; unused by gated generators | open | y
F-d1_docs_specs_reviews-55 | docs/held_performance.md:9-10,18 | initial placement dominates generation time (Monstera 446.5 s) | open | n
F-d1_docs_specs_reviews-56 | docs/held_pfeiffer_reference_pockets.md:53-58 | stale primitive counts vs certified census | open | n
F-d1_docs_specs_reviews-57 | docs/held_toolpath_gallery.md:24-27,69-77 | contour-aware placement, MATHSM, intermediate contour, Fig 6 curves, Fig 7 maps blocked on live (exist on codex/held-reference-corpus per plan) | open | y
F-d1_docs_specs_reviews-58 | docs/held_toolpath_gallery.md:56-58 | IDEA built (Held lane): corner approaches = approach, shared-contact rotation, departure at concave corners; unused by gated generators | open | y
F-d1_docs_specs_reviews-59 | docs/licenses/cgal-adaptive-package-audit.md:3,41-45 | CGAL_LICENSE_GATE_ERROR: LGPL repo instantiates GPL CGAL packages (Arrangements, Boolean ops, SDG); needs GPL decision or commercial entitlement | open | n
F-d1_docs_specs_reviews-60 | docs/local_depletion.md:17-22 | says local depletion is not wired in; generators call subtract_arc_sweep_local (stale) | open (doc) | n
F-d1_docs_specs_reviews-61 | docs/local_depletion.md:147-149 | sound local depletion for disk-chain/capsule needs custom inserting zone visitor; not attempted | open | n
F-d1_docs_specs_reviews-62 | docs/loop_radius_degeneracy.md:39-42 | stations with clearance <= 3r (middle curve) or <= 2r+c (MAT-centred) cannot carry a trochoid; need a non-trochoid motion there | open (draft plan B2/D2) | y
F-d1_docs_specs_reviews-63 | docs/loop_radius_degeneracy.md:104-134 | IDEA built: centre plunge + concentric ramp entry, 240 deg floor, step 124.64 vs 329.49; spiral_entry module unused by gated generators | open | y
F-d1_docs_specs_reviews-64 | docs/loop_radius_degeneracy.md:136-149 | helical entry requires a depletion-model change (_replay_kind rejects Z+XY moves) | open (draft plan D1) | y
F-d1_docs_specs_reviews-65 | docs/loop_radius_degeneracy.md:98-102 | do not exclude the entry motion from the engagement-step criterion | open (draft plan D1 must honour it) | y
F-d1_docs_specs_reviews-66 | docs/loop_radius_degeneracy.md:177-190 | IDEA: external common tangent bridge removes clearance-slope tilt; radius-step bound alone makes tangent breaks worse | closed in rho_regulated (evidence: draft plan table tangent breaks 0/0/2 for rho); open for gated generators | y
F-d1_docs_specs_reviews-67 | docs/machining_metric_validity.md:136-141 | IDEA built and abandoned: chain-ordering stay-down generator (engagement_ordered_toolpath.py, 7af7a24c) took material_entries 5->1; deleted 9fa6eea8 as unreachable; draft plan re-derives it without citing | open | y
F-d1_docs_specs_reviews-68 | docs/machining_metric_validity.md:289-333 | engagement not translation-invariant at exact rim/boundary tangency; deciding layer unknown; test deliberately red | open | y
F-d1_docs_specs_reviews-69 | docs/machining_metric_validity.md:234-242 | straight plunge is the defect; helical ramp entry "under construction" never landed | open | y
F-d1_docs_specs_reviews-70 | docs/machining_metric_validity.md:120-153 | run boundary must be rapid OR path_index change; stay-down cross-chain links need distinct path_index | closed (evidence: doc says both tested in tests/benchmarks/test_quality.py) | y
F-d1_docs_specs_reviews-71 | docs/machining_quality.md:150-158 | uncut residue not corner-confined (3.28 units from nearest vertex); constant-spacing leaves 1/5 as much | open | y
F-d1_docs_specs_reviews-72 | docs/machining_quality.md:150-154 | guide revisits cleared stations: 10 zero-engagement full-size loops | open (draft plan productive-only emission targets it) | y
F-d1_docs_specs_reviews-73 | docs/machining_quality.md:305-313 | W <= 4r channels emitted as 68 plunges; generator must decline or switch strategy | open | y
F-d1_docs_specs_reviews-74 | docs/machining_quality.md:213-224 | no lower engagement bound / chip-thinning compensation in any generator | open | y
F-d1_docs_specs_reviews-75 | docs/machining_quality.md:160-168 | helical or pre-drilled entry would clear cap_exceedances and engagement_step on its own | open (draft plan D1) | y
F-d1_docs_specs_reviews-76 | docs/machining_quality.md:190-195 | rect_20x12 slotting: straight cut out of degenerate corner loop runs into corner at 140 deg | open (draft plan B3) | y
F-d1_docs_specs_reviews-77 | docs/machining_quality.md:296-303 | material_entries and loop_radius_cv deliberately ungated; middle-curve CV reference "does not exist yet" | open | y
F-d1_docs_specs_reviews-78 | docs/measurement_claims.md:18-33 | engagement non-monotone in spacing and guide radius: never bisect advance/spacing | open (engagement_controlled_toolpath still bisects; draft plan reuses rho's ladder scan) | y
F-d1_docs_specs_reviews-79 | docs/number_types.md:1285-1301 | stage 2 leaves invariant 4 undelivered; two incompatible exact-rational-v1 framings; station lane not bound to Python | open (parked) | n
F-d1_docs_specs_reviews-80 | docs/number_types.md:471-486 | filter-failure rate of the real workload unmeasured; pipeline never run under GMP | open | n
F-d1_docs_specs_reviews-81 | docs/number_types.md:90-95,1322-1323 | exact decisions must be tested with generic-double fixtures and a time bound (integer fixtures hide root-bound refinement) | open | y
F-d1_docs_specs_reviews-82 | docs/number_types.md:669-675 | AttestationByteDriftError declared, nothing raises it | open | n
F-d1_docs_specs_reviews-83 | docs/oblique_edge_cost.md:86-101 | add oblique pocket to corpus; profile center_domain symbolicated; decide whether domain is needed | open | y
F-d1_docs_specs_reviews-84 | docs/oblique_edge_cost.md:13-18 | performance/parity claims hold only for axis-aligned pockets | open | n
F-d1_docs_specs_reviews-85 | docs/plans/2026-02-28-full-circle-trochoids-design.md:94-103 | IDEA abandoned: half-arc (partial-arc) trochoid replaced by full 360 deg circles; full circles recut ~90% of swept area | open (no generator uses partial-arc loops) | y
F-d1_docs_specs_reviews-86 | docs/plans/2026-02-28-toolpath-cgal-redesign.md:515-520,745-760 | IDEA abandoned: DP tangent/transition solver (cost-optimal transitions) removed for fixed-side rule | open (not revisited) | y
F-d1_docs_specs_reviews-87 | docs/plans/2026-02-28-toolpath-cgal-redesign.md:262-265 | external tangent clamp with 1e-9 when |dr| >= advance | unknown (literal absent from current toolpath.cpp) | y
F-d1_docs_specs_reviews-88 | docs/plans/2026-07-23-adaptive-clearing-design.md:96-124 | IDEA never built: SP2 contour-parallel offset-loop clearing with propose/certify/adapt and trochoidal relief only where regulation fails | open | y
F-d1_docs_specs_reviews-89 | docs/plans/2026-07-23-adaptive-clearing-design.md:34-35 | IDEA never built: engagement target band (>=95% of cut length within target±band), i.e. constant engagement | open | y
F-d1_docs_specs_reviews-90 | docs/plans/2026-07-23-adaptive-clearing-design.md:122-124 | IDEA never built: SP2b spiralization of adjacent loops (Held/Spielberger) | open | y
F-d1_docs_specs_reviews-91 | docs/plans/2026-07-23-adaptive-clearing-design.md:144,173-174 | IDEA never built: SP3 chip-thinning feed post-processor | open | y
F-d1_docs_specs_reviews-92 | docs/plans/2026-07-23-adaptive-clearing-design.md:115-117 | entry per connected region at max-clearance skeleton vertex with certified slot primitive | open (draft plan adopts root rule; certified slot primitive not used) | y
F-d1_docs_specs_reviews-93 | docs/plans/2026-07-23-adaptive-clearing-design.md:186-189 | future: rest-machining chains, G2/jerk-limited blending, true-MA guides (Segment_Delaunay_graph_2) | open | y
F-d1_docs_specs_reviews-94 | docs/plans/2026-07-24-engagement-zone-query.md:7-13 | land local zone engagement query | closed (evidence: src/engagement_2.cpp:310 engaged_arcs_zone) | n
F-d1_docs_specs_reviews-95 | docs/plans/2026-07-23-sp1-stock-tea-kernel.md:1409-1490 | SP1 wall-clock exit gate test (test_stock_gate.py) | open (file absent on live) | n
F-d1_docs_specs_reviews-96 | docs/plans/2026-07-23-sp1-stock-tea-kernel.md:1354-1405 | baseline engagement audit report docs/benchmarks/engagement_baseline.md (SP2 reference) | open (absent on live) | y
F-d1_docs_specs_reviews-97 | docs/plans/2026-07-23-sp1-stock-tea-kernel.md:884-1037 | IDEA built, unused by gated generators: certify_segment_tea continuous cap certificate via guarded growth bound | open | y
F-d1_docs_specs_reviews-98 | docs/plans/2026-07-23-sp1-stock-tea-kernel.md:1494-1503 | residual extraction (residual_polylines) deferred to SP2 | open | y
F-d1_docs_specs_reviews-99 | docs/radius_ladder_ablation.md:46-56 | at forced stations every candidate is over cap; lowering peak costs a circle (structural opposition) | open | y
F-d1_docs_specs_reviews-100 | docs/radius_regulated_toolpath.md:28-32,198-200 | bridge cuts unregulated in both gated generators; worst 141.8 deg | open (draft plan does not require cap-checked bridges) | y
F-d1_docs_specs_reviews-101 | docs/radius_regulated_toolpath.md:156-171 | corner-apex stations (radius 0.057 vs tool 1.0) have one rung; geometry pins max TEA 125.8 deg | open (draft plan B2/D2) | y
F-d1_docs_specs_reviews-102 | docs/radius_regulated_toolpath.md:202-218 | multiple passes per chain, one plunge per pass; MAX_RADIAL_SWEEPS_PER_CHAIN empirical budget raises | open | y
F-d1_docs_specs_reviews-103 | docs/refutation_soundness_crosscheck.md:271-277 | put parameters/identities sizes in IncompleteSegmentPartitionError; decide split/merge correspondence rule k != m | open | n
F-d1_docs_specs_reviews-104 | docs/refutation_soundness_crosscheck.md:279-331 | cpp_int->double overflow abort (float_next inf) on 139-digit coordinates; missing overflow guard; GMP effect untested | open | n
F-d1_docs_specs_reviews-105 | docs/refutation_soundness_crosscheck.md:167-181 | MATERIAL vs CAP_EXCEEDED split not exposed at Python boundary | open | n
F-d1_docs_specs_reviews-106 | docs/segment_site_mat.md:4358-4359,3733-3734 | exact segment-site MAT production graph exists only for one exact rectangle contract and one exact six-edge L; arbitrary pockets pending | open | y
F-d1_docs_specs_reviews-107 | docs/segment_site_mat.md:1417-1420,1453-1456 | contour-aware placement (max engagement over candidate vs evolving contour, Fig 4(d) overlap) not built on live | open | y
F-d1_docs_specs_reviews-108 | docs/segment_site_mat.md:1672-1678 | Figure 5 refined path: 1,010 exact containment rejections, publisher-derived start, no depleted-stock engagement, Task 5 acceptance open | open | y
F-d1_docs_specs_reviews-109 | docs/segment_site_mat.md:2059-2066 | Phase-1 hard target: planner at 1/10 of Held's published clocks (10 ms at 20 deg ... 0.3-0.6 ms) | open | n
F-d1_docs_specs_reviews-110 | docs/segment_site_mat.md:2101-2105 | higher-order smoothing of certified arc/line paths does not preserve engagement bound (G2 open problem) | open | y
F-d1_docs_specs_reviews-111 | docs/segment_site_mat.md:1025-1055 | IDEA built (adaptive/entry.py): qualified precleared bore as entry; unused by product generators and not considered by draft plan D1 | open | y
F-d1_docs_specs_reviews-112 | docs/segment_site_mat.md:1619-1680 | IDEA built (Held lane): convex-corner radius limits + concave-corner shared-contact fans; unused by product generators | open | y
F-d1_docs_specs_reviews-113 | docs/segment_site_mat.md:1706-1715 | exact circular-arc input needs a declared construction policy (stored arc endpoints unequal radii ~8.5e-16) | open | n
F-d1_docs_specs_reviews-114 | docs/segment_site_mat.md:4090-4093 | rectangle producer emits no clearance profiles / neck evidence | open | n
F-d1_docs_specs_reviews-115 | docs/segment_site_mat.md:4394-4396 | planned: tri-dexel replay (Task 14A), GPU CDF scoring | open | n
F-d1_docs_specs_reviews-116 | docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md:381-385 | center-distance growth lemma (certify_segment_tea guard) is not sound: coincident equal-radius arcs, overlap visitor drops positive-length coincident arcs | open (contradicts F-97 as a usable certificate; related F-68) | y
F-d1_docs_specs_reviews-117 | docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md:81-86,260-289 | IDEA built (ReachableDomain M_r): unreachable residual D\M_r reported separately, not a coverage failure; product gate does not use it | open | y
F-d1_docs_specs_reviews-118 | docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md:777-779 | IDEA never built: Phase 1.5 generated helical entry + Held sorted arc contour; Phase 2 offset-contour transitions and elliptic returns | open | y
F-d1_docs_specs_reviews-119 | docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md:59-61 | a tool-radius plunge cannot launch a capped lateral move against frozen stock (>pi) | open (draft plan D1 / plunge placement) | y
F-d1_docs_specs_reviews-120 | docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md:583 | max_passes truncation forbidden | open (CHANGELOG.md:19 warns instead; draft plan Review Focus 4 fails loud) | y
F-d1_docs_specs_reviews-121 | docs/superpowers/specs/2026-07-25-task3-single-arrangement-performance-design.md:193-195 | ReachableDomain requires one connected C_r (single-entry contract); multi-region pockets raise PocketNotMachinableError | open | y
F-d1_docs_specs_reviews-122 | docs/superpowers/specs/2026-07-25-task3-single-arrangement-performance-design.md:74-77,404-415 | post-Task5 measurement gate on hash indexing | unknown | n
F-d1_docs_specs_reviews-123 | docs/superpowers/specs/2026-07-28-causal-mat-traversal-design.md:155-167 | DFS route over MAT edges from entry component = draft plan tree walk; measured to need an inter-route transit (direct links gouge) | open (draft plan uses retract/plunge; retrace_*.py exists unscored) | y
F-d1_docs_specs_reviews-124 | docs/superpowers/specs/2026-07-28-causal-mat-traversal-design.md:277-301,324-340 | terminal coverage, fresh replay, GenerationResult, acceptance on tractable fixture | open | n
F-d1_docs_specs_reviews-125 | docs/superpowers/specs/2026-07-29-one-root-full-circle-source-design.md:67-69 | vertices with distinct nonzero x/y radicands unsupported in v4 source | open | n
F-d1_docs_specs_reviews-126 | docs/superpowers/specs/2026-07-29-zero-guide-link-candidate-design.md:163-182 | IDEA built: distinct AdvanceSegment (owns removal) vs hold Link (reposition) operation grammar; product gate lacks the distinction (penalises stay-down links) | open (draft plan D3) | y
F-d1_docs_specs_reviews-127 | docs/superpowers/specs/2026-08-06-exact-inter-route-retrace-design.md:36-45,78-82 | IDEA built then superseded: exact stay-down retrace through accepted lineage (swept-prefix certified); retract/re-plunge explicitly rejected; draft plan adopts retract/plunge | open (tests red, docs/red_manifest.json "plan superseded (backlog D1)") | y
F-d1_docs_specs_reviews-128 | docs/red_manifest.json:21-45 (cited from this spec's implementation) | four adaptive generator/retrace tests expected-red pending retirement via auditor programme (backlog S2) | open | n
F-d1_docs_specs_reviews-129 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:599-615 | swept-prefix fast path research-only until proof note + adversarial tests; certifies from two booleans | open (P2 parked) | y
F-d1_docs_specs_reviews-130 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:643-657 | generator override ledger: forced entries/loops/advances must survive public result; require_no_overrides() | open (P3 parked) | y
F-d1_docs_specs_reviews-131 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:659-666 | one shared ladder-search module; no regulated generator may use monotone bisection | open (engagement_toolpath.py still bisects) | y
F-d1_docs_specs_reviews-132 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:673-678 | tangent tests: comparison-count floor; bridge exemptions validated geometrically | open | y
F-d1_docs_specs_reviews-133 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:682-683 | true helical motion fails loud until a typed 3D depletion contract exists | open (draft plan D1 needs it) | y
F-d1_docs_specs_reviews-134 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:50,617-619 | unmeasured=certified in legacy audit; empty segment strata must be unresolved | unknown (P1 work on live; legacy path unchanged) | n
F-d1_docs_specs_reviews-135 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:701-747 | Python 3.12 floors, six strict-mypy errors, native build-identity manifest, CI topology | open (progress.md:39 six errors persist) | n
F-d1_docs_specs_reviews-136 | docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:570-597 | ReplayCertificate value-returning path with terminal complete-path positive control | open | n
F-d1_docs_specs_reviews-137 | docs/superpowers/specs/2026-08-29-wave1-task6-historical-assertions.md:96-111 | Task 6 completion (v2 measurement artifact, ledger) | closed (evidence: superseded by 6b4293ae subsystem removal) | n
F-d1_docs_specs_reviews-138 | docs/superpowers/specs/2026-08-31-held-pfeiffer-reference-pockets-design.md:30-36,233-242 | Held pockets supplied as inputs; generator quality/cap/coverage on them are later measured stages | open | y
F-d1_docs_specs_reviews-139 | docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md:803-810 | Phase 1 Figure 5 characterization report never committed; Phase 2 generator-policy design never started | open | y
F-d1_docs_specs_reviews-140 | docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md:538-547 | coverage grid must scale per pocket (668x458 for Figure 5); 200-sample default too coarse for curved pockets | closed for Figure 5 adapter (evidence: doc :543-546 minimum_coverage_grid); general gate unknown | y
F-d1_docs_specs_reviews-141 | docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md:680-709,812-814 | IDEA never built: post-qualification candidate seam + target-specific postprocessor (no G-code output exists) | open | y
F-d1_docs_specs_reviews-142 | docs/superpowers/specs/2026-09-02-held-figure5-2d-path-characterization-design.md:588-605 | Union_of_curve_cycles_2 single-arrangement experiment rejected (timeout + GPL license) | closed (evidence: doc :602-605 rejected; d40ef841 "docs(sdd): reject single-arrangement experiment" in log) | n
F-d1_docs_specs_reviews-143 | docs/superpowers/specs/2026-09-08-number-type-coherence-design.md:301-313 | stages 3-6 (segment source inversion, segment_site BigRat->Rational, coefficient discipline, decoder removal) | open (parked) | n
F-d1_docs_specs_reviews-144 | docs/superpowers/specs/2026-09-08-number-type-coherence-design.md:336-346 | dead second double->exact door exact_binary64/lift_binary64 | open (stage 3) | n
F-d1_docs_specs_reviews-145 | docs/toolpath_drawing.md:106-111,119 | gated generators emit no lead-in/lead-out moves | open | y
F-d1_docs_specs_reviews-146 | docs/toolpath_tangent_continuity.md:33-51,167-171 | tangent test exempts bridge<->circle junctions as "inherent ~90°"; gate counts them as breaks; external tangent removes them | open | y
F-d1_docs_specs_reviews-147 | docs/toolpath_tangent_continuity.md:83-90 | 1e-18 squared-distance threshold decides repositioning arc vs coincident (epsilon in a geometric decision) | open | y
F-d1_docs_specs_reviews-148 | docs/toolpath_tangent_continuity.md:152-165 | parallelism-only check via abs(dot) cannot detect direction reversal | open | n
F-d1_docs_specs_reviews-149 | docs/toolpath_trochoid_chain_design.md:74-97 | doc describes superseded adaptive-winding arc as current | open (doc stale) | n
F-d1_docs_specs_reviews-150 | codex/sdd-coherence:docs/superpowers/state/backlog.md:94-97 | C4: depletion model to score a true helical entry (entry criterion unsatisfiable until then) | open (draft plan D1) | y
F-d1_docs_specs_reviews-151 | codex/sdd-coherence:docs/superpowers/state/backlog.md:97 | C5: arc motion vocabulary for tangent-continuous inter-chain links | open (not in draft plan) | y
F-d1_docs_specs_reviews-152 | codex/sdd-coherence:docs/superpowers/state/backlog.md:65 | C1: identify which layer decides exact-tangency translation variance | open | y
F-d1_docs_specs_reviews-153 | codex/sdd-coherence:docs/superpowers/state/backlog.md:79-83 | W2.1/W2.2/S1/S2/R5: dispositions of 36 certifier commits, plan headers, absorb sdd-coherence, retire 4 retrace reds, Python floor + 6 mypy errors | open/unknown | n
F-d1_docs_specs_reviews-154 | codex/sdd-coherence:docs/superpowers/state/backlog.md:89-92 | W3: one branch, CI green on main, re-issue due-diligence memo (predates oblique-edge finding), convergence oracle | open | n
F-d1_docs_specs_reviews-155 | codex/sdd-coherence:docs/superpowers/state/backlog.md:96 | claim "V >= sqrt(3) r proved in loop_radius_degeneracy.md" absent from live page | open (contradicted) | y
F-d1_docs_specs_reviews-156 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:35-61 | C1 mat_scale>1 emits gouging circles; circles never certified; validate (0,1] | open (live src/toolpath.cpp:607-613 lacks mat_scale check) | y
F-d1_docs_specs_reviews-157 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:63-129 | C2 tea_growth_bound not an upper bound at concave voids; curvature-aware bound needed | unknown (fix commits on jf/toolpath-redesign not ancestors of live) | y
F-d1_docs_specs_reviews-158 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:131-137 | C3 requires-python >=3.9 but modules need >=3.10 | open (pyproject.toml:9) | n
F-d1_docs_specs_reviews-159 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:145 | I1 Python arc certifier cannot certify any circle in generator regime (0/126) | open (engagement.py:48 unchanged) | y
F-d1_docs_specs_reviews-160 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:146 | I2 max_tea can report 720 deg (duplicate rim sub-arcs) | unknown | n
F-d1_docs_specs_reviews-161 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:148 | I4 geometry params unvalidated (negative tool radius, clearance_z <= cut_z silently flat links) | unknown | y
F-d1_docs_specs_reviews-162 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:149 | I5 link_paths=False emits uncertified through-material traverse in polyline | unknown | y
F-d1_docs_specs_reviews-163 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:151 | I7 Epick gouge decisions compare rounded doubles (segment_clear) | unknown | y
F-d1_docs_specs_reviews-164 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:152 | I8 certify_bridges can exit with an uncertified bridge | unknown | y
F-d1_docs_specs_reviews-165 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:153 | I9 hole validation incomplete before straight skeleton (UB) | unknown | y
F-d1_docs_specs_reviews-166 | jf/toolpath-redesign:docs/superpowers/specs/2026-08-19-eth-audit-findings.md:175-176 | _unmeasured returns cap_certified=True; cap_violations under-counts | open (engagement.py:320-322) | n
F-d1_docs_specs_reviews-167 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3065-3095 | Task 30: gate rho/chain-ordered/spiral generators; chain ordering 5->1 entries | open (gate.py registers 2 generators; chain-ordered deleted 9fa6eea8) | y
F-d1_docs_specs_reviews-168 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3099-3135 | Task 31: consume declined_regions with a corner-clearing pass or declared second tool; split uncut into declined vs silently missed | open (no consumer, verified grep) | y
F-d1_docs_specs_reviews-169 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3164-3194 | Task 33: cutting-force, exit-impact, controller-starvation KPIs | open | y
F-d1_docs_specs_reviews-170 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3198-3233 | Task 34: generalise segment-site MAT; middle curve measured not gentler (CoV 0.5719 vs 0.5752) and more degenerate loops (24.5% vs 13.0%) | open | y
F-d1_docs_specs_reviews-171 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3266-3347 | Task 36: helical/ramped/pre-bored entry; audit must replay it | open (draft plan D1) | y
F-d1_docs_specs_reviews-172 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3351-3434 | Task 37: break ~141.8 deg saturation; decide whether station set or regulator binds; length curve must steepen | open (draft plan Stage C cap 40) | y
F-d1_docs_specs_reviews-173 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:3022-3061 | Task 29: second gate cap so generator defects are attributable | closed (evidence: docs/red_manifest.json cap-40-attribution cells) | y
F-d1_docs_specs_reviews-174 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:2988-3018 | Task 28: replay_certificate unconditional raise; RetraceSegmentOperation rejected by replay grammar | unknown | n
F-d1_docs_specs_reviews-175 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:2879-2927 | Task 25: gate the four prose-only numbers (136.4/119.5, 914.3, 224 vs 2.6 ms, 125.8) | open | n
F-d1_docs_specs_reviews-176 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:2931-2954 | Task 26: exercise island pockets on engagement generators (holes never tested there) | open | y
F-d1_docs_specs_reviews-177 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:1994-2097,2101-2199 | Tasks 13/14: bridge certification exit, hole-edge/nested-hole validation | unknown | y
F-d1_docs_specs_reviews-178 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-eth-audit-remediation.md:2393-2462 | Task 17: external_tangent via kernel primitives (hand-rolled perpendicular) | unknown | n
F-d1_docs_specs_reviews-179 | jf/toolpath-redesign:INTEGRATECODEX.md:157-161 | residue 5x worse than unregulated, not corner-confined; coverage absent from the accept/reject rule | open | y
F-d1_docs_specs_reviews-180 | jf/toolpath-redesign:INTEGRATECODEX.md:77-83 | swept-prefix certifies only pi; any other cap aborts the first zero-guide advance; adaptive fixture runs at pi | open | y
F-d1_docs_specs_reviews-181 | jf/toolpath-redesign:INTEGRATECODEX.md:85-95 | non-finite doubles empty Stock2; shared-root guard compiled out; disk-chain int cast -> 34 GB reserve | unknown (fixes on jf/toolpath-redesign, not ancestors of live) | n
F-d1_docs_specs_reviews-182 | jf/toolpath-redesign:INTEGRATECODEX.md:163-165 | 4-5 degenerate loops (rho 0.0215) and 10 zero-engagement redundant loops; both exactly decidable | open (draft plan Stage A targets) | y
F-d1_docs_specs_reviews-183 | jf/toolpath-redesign:MOTIVATION.md:7,13 | IDEA: reduce theta_max at bottlenecks; minimise ~50% non-cutting air time; G1/G2 for high feed | open (bottleneck caps only in parked adaptive lane) | y
F-d1_docs_specs_reviews-184 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:401-405 | radius_regulated imports the bisecting _largest_admissible_advance from engagement_toolpath (non-monotone engagement) | open (verified src/compas_cgal/engagement_radial_toolpath.py:138) | y
F-d1_docs_specs_reviews-185 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:318-337 | forced loops/advances erased at public result (warnings only); RadialToolpathResult drops forced_advances | open | y
F-d1_docs_specs_reviews-186 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:364-372 | Hypothesis pockets are fat stars only: no necks/pinches/islands; generators untested on hard topology | open | y
F-d1_docs_specs_reviews-187 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:398-400 | corpus exercises unregulated generator; rho generator wired into nothing | open (gate still lacks rho) | y
F-d1_docs_specs_reviews-188 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:231-254 | arc certifier fixed 20 stations/turn: guard > full turn past 2x tool radius; 47-52% uncertified | open (engagement.py:48) | n
F-d1_docs_specs_reviews-189 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:294-308 | swept-prefix certifies from two booleans with no geometric test | open (auditor P2 parked) | y
F-d1_docs_specs_reviews-190 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:378-387 | segment oracle lacks empty-strata unresolved guard (sibling circle oracle has it) | unknown | n
F-d1_docs_specs_reviews-191 | jf/toolpath-redesign:review-2026-08-22-t9-zero-guide.md:388-397 | no committed measured artifact; prose numbers stale (224 ms vs 72.4 ms) | open (engagement_controlled_toolpath.md:7 vs :233 still inconsistent) | n
F-d1_docs_specs_reviews-192 | jf/toolpath-redesign:review.md:460-491 | loops maximised not regulated: 171.8 deg at 40 and 120 caps; useful 30-60 deg range unreachable | open (draft plan Stage C) | y
F-d1_docs_specs_reviews-193 | jf/toolpath-redesign:review.md:351-353 | claim "dropping degenerate loops changes coverage 0.000%" contradicts machining_quality.md:195-197 (0.284% -> 0.351%) | open (contradiction) | y
F-d1_docs_specs_reviews-194 | jf/toolpath-redesign:review.md:567-570,782-785 | measure unresolved rate on third-party geometry; neck-degeneration ablation | open | n
F-d1_docs_specs_reviews-195 | jf/toolpath-redesign:review.md:529-534 | const_cast in engaged_arcs_zone blocks parallel engagement queries | open | n
F-d1_docs_specs_reviews-196 | jf/toolpath-redesign:review.md:524-527,599-611 | CGAL licence gate unresolved; condition precedent to any commercial use | open | n
F-d1_docs_specs_reviews-197 | jf/toolpath-redesign:CONVENTIONS.md:126-133 | CONVENTIONS.md points pixi users at a removed t9-zero-guide worktree | open | n
F-d1_docs_specs_reviews-198 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md:1613-1700,2048-2091 | neck (dumbbell) and island families exist but are not product-gate pockets | open | y
F-d1_docs_specs_reviews-199 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md:1525-1595 | arc-fraction family is tessellated; true arc boundaries never exercised by the corpus | open | n
F-d1_docs_specs_reviews-200 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md:664-688,3097-3104 | Pythagorean rotations via binary64 a/c are not exactly rational; congruence "exact" claim false | open | n
F-d1_docs_specs_reviews-201 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md:3141-3146 | resolution-rate on external corpus and neck-degeneration ablation never run/recorded | open | n
F-d1_docs_specs_reviews-202 | jf/toolpath-redesign:docs/superpowers/plans/2026-08-19-benchmark-corpus.md:2353-2522 | IDEA built (benchmark only): MATHSM constant-spacing sweep baseline; no MATHSM generator on live (MATHSM traversal only on codex/held-reference-corpus per draft plan) | open | y
