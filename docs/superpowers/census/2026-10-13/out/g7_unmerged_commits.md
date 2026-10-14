# g7 — commits not on `codex/held-task5-continuation`

Live branch = `codex/held-task5-continuation`. Merge bases: `jf/toolpath-redesign` forked at `1860167` (2026-07-24; live has 503 commits since); `codex/held-reference-corpus` and `codex/stock-clone-ownership` forked at `43f49b5` (2026-09-06; live has 97 commits since). `git cherry` marks 3 of the 57 commits patch-equivalent on live (`23a1de0a`, `6f404354`, `4c058e33`); the other 54 are not on live in any form `git cherry` can see. Every "missing on live" claim below was checked by grepping the live tree, and the anchors are cited.

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **A publisher-free Held standard-path generator exists only on `codex/held-reference-corpus`** — verified. Commits `757b948c`, `5b5ed867`, `4685653b`, `567da10d` and `274fd1a4` build `benchmarks/held_figure5_analytic.py`. It works from boundary geometry alone, as a `ReferenceBoundary` of lines and arcs. The start is a geometric "south-east support" point on the r-offset contour, and a test pins that it does not change when the boundary primitives are rotated in storage. The generator walks the cutter-centre offset contour once and emits **one continuous stream**: full CCW machining circles that span the wall to the MAT, joined by exact-endpoint line/arc offset connectors. It stops where the stream re-enters the initial clearance circle. On Figure 5 at 80° it pins 270 circles, with predecessor engagement in (80° − 1e-7 rad, 80°] (test). That contradicts the draft plan's premise that the Held pipeline "needs the publisher-derived start circle and raw guide" (plan:38). Limits:
   - Engagement comes from Held's predecessor-only model, not the gate's stock replay.
   - The successor search assumes golden-section unimodality inside 0.25·r cells. That is a heuristic, not a certificate.
   - Figure 5 is declared one smooth composite (`logical_source_indices` all 0, `274fd1a4`).
   - On sharp convex corners, code reading suggests the MAT candidate jumps to the opposite wall once the adjacent-edge candidate falls to `≤ r`, giving an oversized circle (`_locus_circle_from_state` selects `min(d > tool)`). **Unverified**: it has not been run on `rect_12x8`.
2. **Contour-aware placement (Held §3.1) exists only off live** — verified. `8ca2c9f6` adds `benchmarks/held_figure5_contour_placement.py`. It maintains the exposed union contour of the accepted clearance disks and moves the Fig. 4 entry point `b` onto it, which gives a shorter path with fewer circles than the standard path (asserted in a test). A shared policy-driven traversal (`274fd1a4`) runs the standard, contour-aware and MATHSM constant-spacing policies on **one** locus, so Figure 6 becomes like-for-like. Live's `benchmarks/mathsm.py` uses `trochoidal_mat_toolpath_circular` instead, a different guide (live `benchmarks/figure6.py` docstring).
3. **Generic segment-site Voronoi guide.** Live's exact segment-site MAT accepts only the L fixture: `require_l_shape_input` is still called at `src/segment_site_catalog_graph.cpp:1220` (verified). `07c87ce0`, `f90d1ce4` and `334c028b` add `SegmentSiteGuide2` (`src/segment_site_guide.cpp`). It builds the CGAL segment-Delaunay graph for any polygon with holes and emits a connected MAT graph with generator-site provenance. Parabolic arcs are evaluated analytically (`benchmarks/held_segment_site_locus.py`), and spurious vertex/incident-segment branches are dropped. It is double-sampled (point merge at 8 ε, domain filter on sampled points), so it is a usable guide, not a certificate. This is the "MAT instead of the straight skeleton" capability the brief asks about. Only a non-exact version exists, and only off live.
4. **In-motion (swept-prefix) engagement for lines and arcs at any cap** — verified off-live-only. `7be75e0d` adds `engagement_line_prefixes` / `engagement_arc_prefixes` (`src/swept_prefix_engagement_2.cpp`). They measure engagement at stations *inside* a motion against the stock minus that motion's own swept prefix, with exact rim clipping and rational witnesses, and then deplete the stock. Live's swept-prefix theorem certifies only a π cap on segments (`src/continuous_tea_2/segment_oracle.h:50-57`). This is the instrument needed to score stay-down links and entries. Weaknesses: arc sweep angle and `swept_angle <= pi` are doubles in a covered/not-covered decision, and stations are sampled.
5. **The live audit still certifies on the falsified growth lemma** — verified:
   - `src/engagement_2.cpp:407-425` (`tea_growth_bound`, still carrying its "QED"), `:445-499` (`certify_recursive`), and `src/compas_cgal/engagement.py:286` (the arc mirror).
   - jf found five false certificates (rib/spiral rib/sector, `eee1c35b`, `3c270cb5`) and replaced the interior proof with an exact swept-annulus bound (`2e2422b5`, `03866fcc`, `73d53729`).
   - Live keeps the witnesses only as non-certifying negative controls (`tests/test_false_certificate.py`, added in `7cda8f6d`).
   - Consumers that can report a false `cap_certified=True`: `benchmarks/pathmetrics.py:198`, `benchmarks/runner.py:86`, `benchmarks/held_consumer_adapters.py:20`. The product gate does not use them (it uses sampled exceedance).
   - The jf bound is sound but cannot certify a centreline pass down a tool-width corridor. jf's own doc (`docs/engagement_certificate.md`, "Open") reports a 3.14 rad disk bound against a 0.01 rad true segment maximum. So it is an audit repair, not a generation instrument.
6. **`src/toolpath.cpp` is untouched on live since the fork** (verified: `git diff 1860167..live -- src/toolpath.cpp` is empty), so every jf generator-contract fix is missing:
   - `mat_scale > 1` gouges (`d90e18b8`; live `src/toolpath.cpp:600-613`).
   - The tessellator drops op start points, splicing discontinuous moves into the polyline (`b2fa087e`/`34bddcc0`; live `:825-870`).
   - `clearance_z <= cut_z` silently degrades to flat linking (`92cedda2`; live `:996`).
   - Negative `radial_clearance` and inverted min/max radius are accepted (`92cedda2`).

   **Conflict:** `b0094cb6` raises on any uncertified implied traverse when `link_paths=False` and no `clearance_z` is given. That is exactly how the gated generator reads its guide (`src/compas_cgal/engagement_toolpath.py:438-452`), so porting it alone can break `engagement_controlled_toolpath` on non-convex pockets (unverified which gate pockets). Port it together with `446daca5`, or make the guide call pass a `clearance_z`.
7. **`d0f32e8d` and live's `57482731` + `8d2b76cc` fix the same stock-clone traits-lifetime bug in different ways** — verified:
   - Live: resolves and pins the Gps whose traits the arrangement reads (`src/stock_2.cpp:188-207`, `:279-292`, `:306-325`).
   - `d0f32e8d`: rebuilds a fresh, self-owning Gps from `polygons_with_holes` (no ancestor pinning; its doc measures 0.725→0.784 ms per clone).
   - Both are correct by their own tests. They conflict in `Stock2::clone` and in the same new test filename. Keep live's.
8. **Live's publisher Figure 5 evidence miscounts turns, 265 vs 278** — verified at `benchmarks/held_figure5_publisher.py:33,198`. `a662a9c3` decodes 17 five-cubic circles that live classifies as connectors. Every live publisher-relative residual (circle count, connector coverage) inherits the error. This affects comparison, not generation.
9. **Figure 5 reference boundary G1 seam.** `567da10d` repairs a reconstruction-induced line/cubic tangent break that created a radius-to-zero MAT branch, within a bounded publisher deviation (`0.0738 ≤ 0.07386` mm, test-pinned). Live has no `_figure5_g1_sources` (grep, verified), so live Figure 5 guides may carry that spurious branch. **Unverified** on live.
10. **jf's docs (`c90a4fef`, `f369b3df`) record measured, generation-critical diagnoses that the live plan does not cite:**
    - the ~141.8° cap saturation (20° cap yields the same measured max as 140°, with 264 forced advances);
    - "every cap exceedance is a chain entry", verified by set equality;
    - `rho_max = (d-r)/2`, so `d <= 3r` corners cannot carry a trochoid;
    - "middle curve is not gentler than the skeleton" (CoV 0.5719 vs 0.5752; degenerate loops 24.5% vs 13.0%).

    All are unverified by me (relayed from those docs).

### Generation capabilities that exist ONLY off the live branch

- Publisher-free analytic Held **standard** path (one continuous circle + offset-connector stream, geometric start, clearance-circle termination, explicit `motions` stream): `benchmarks/held_figure5_analytic.py` (corpus).
- Held §3.1 **contour-aware** placement (exposed union-contour maintenance, adjusted entry point `b`): `benchmarks/held_figure5_contour_placement.py`, plus `predecessor_entry_point` / `engagement_from_predecessor_entry` in `held_standard_placement.py` (corpus).
- **Placement-policy traversal** (`build_analytic_figure5_path_with_policy`, `Figure5TraversalInterval` with one-sided piece ownership) and **MATHSM constant centre-spacing on the same locus** (`benchmarks/held_figure6_mathsm.py`) (corpus).
- **Figure 6 sweep** over publisher abscissae with typed success/failure samples and sampled cap compliance (`benchmarks/held_figure6_paths.py`, `held_figure6_path_metrics.py`) (corpus).
- **Generic segment-site Voronoi guide** for arbitrary polygons with holes (`src/segment_site_guide.{h,cpp}`, `_medial_axis_2.SegmentSiteGuide`), its analytic parabola evaluation (`benchmarks/held_segment_site_locus.py`), and a boundary-ordered circle placement on it (`benchmarks/held_figure5_voronoi.py`) (corpus).
- **Swept-prefix in-motion engagement** for lines and arcs at any cap, plus a sequential depleting replay at a declared pitch with a density-doubling residual (`src/swept_prefix_engagement_2.cpp`, `benchmarks/held_figure5_swept_engagement.py`) (corpus).
- **Analytic offset-contour construction** with convex joins trimmed to the support intersection and concave joins as tool-radius arcs about the vertex (`_repair_convex_offset_joins`, `_insert_concave_offset_joins`) (corpus).
- **Adaptive arc certification** `certify_arc_tea` (bisection over the arc parameter on the swept-annulus bound) (jf `73d53729`; incomplete per its message).
- **Generator contracts** on `trochoidal_mat_toolpath(_circular)`: `mat_scale ∈ (0,1]`, `radial_clearance ≥ 0`, `min ≤ max` trochoid radius, `clearance_z > cut_z`, implied-traverse certification or lift, and a discontinuity-preserving tessellator (jf).

### Where this slice shows the draft plan (`2026-10-13-product-gate-tree-walk.md`) is wrong or incomplete

- **plan:38 is contradicted.** "`build_held_reference_path` needs the publisher-derived start circle and raw guide, so it cannot run on rect_12x8…" holds for live's pipeline only. The corpus analytic generator needs neither. Task 0b should also score `build_analytic_figure5_path` and `build_analytic_figure5_contour_path`. Their `motions` stream (full circles + line/arc connectors) is already the explicit motion list the planned `held_path_to_toolpath` adapter must reconstruct. Running them on the gate pockets first needs a line-only `ReferenceBoundary` from the gate polygons, and the convex-corner behaviour must be tested (finding 1).
- **The design only ever considers skeleton-chain ordering.** Held's ordering walks the offset contour once, with circles spanning wall to MAT. By construction that gives one stream and one entry, and no chain finishes inside the centre (root causes 1 and 2 of the plan's "where we are"). The plan's Stage A should at least name it as the alternative its Task 0b decision gate is testing.
- **Stage C, "rebase the 20 commits", underestimates the conflicts:**
  - `CMakeLists.txt` (live +81/−… since base);
  - `src/engagement_2.cpp` and `src/audit_exact_station_2.cpp` (live removed about 130 and 122 lines in the regions `7be75e0d` edits);
  - `benchmarks/held_standard_placement.py` (live 66 lines changed);
  - `ReferenceBoundary` gains a required `logical_source_indices` field (an API change for every live constructor);
  - `figure5.json` and `figure5_publisher_turns.json` data changes.

  `8ca2c9f6` imports names that are only defined in the next commit `274fd1a4` (verified with `git show 8ca2c9f6:benchmarks/held_figure5_analytic.py`), so the two must land together.
- **The plan's tree walk reuses guide stations read from `trochoidal_mat_toolpath_circular(link_paths=False)`.** The plan should state that jf's contract fixes to that function are pending and that `b0094cb6` alone would make the guide call raise (finding 6).
- **"Productive-only emission" decides removal with `Stock.exactly_equals`, which only says whether something was removed.** It cannot score engagement *during* a stay-down link or entry. The corpus swept-prefix query (finding 4) is the existing instrument for D1/D3, and the plan does not mention it.
- **D1 (entry):** the corpus models entry as a pre-cleared disk (`stock.subtract_disk(... guide_radius + r)` in `held_figure5_swept_engagement.py`), the same gap the plan names. No helix exists on any branch in this slice.

---

## jf/toolpath-redesign

### commit 33cbcb6d docs(engagement): mkdocs page for zone-query speedup + measured A/B | files 4 | read full
- **generation relevance:** n. This is a performance doc for `engagement_at` (zone query: 92.2 s → 19.1 s audit, 12.4× per call, relayed).
- **on live branch?:** partly. Mermaid fences are present (`mkdocs.yml:79`); `docs/engagement_zone_query.md` is missing; the BLUF-label edits are present in some form.
- **value for toolpaths:** none directly. It is a good, honest write-up of a 4.8× audit speedup.
- **conflicts:** `findings.md` was edited on live afterwards. Minor.
- **should it be brought over:** partly. Bring the page only, re-measured against the live kernel. Low priority (docs-only; slice d1).

### commit 23a1de0a feat(bench): PocketSpec foundation with validated invariants | files 6 | read full
- **generation relevance:** n (benchmark spec).
- **on live branch?:** equivalent present (patch-identical: live `a1c01017`; `benchmarks/spec.py`).
- **value for toolpaths:** none beyond the corpus.
- **conflicts:** none.
- **should it be brought over:** no. It is already there.

### commit 6f404354 feat(bench): analytic pocket families with closed-form clearance | files 2 | read full
- **generation relevance:** n (test pockets).
- **on live branch?:** equivalent present (live `87ea5871`).
- **value for toolpaths:** indirect, through closed-form test pockets.
- **conflicts:** none.
- **should it be brought over:** no.

### commit 4c058e33 test(bench): exact congruence invariants via Pythagorean rotations | files 2 | read full
- **generation relevance:** n.
- **on live branch?:** equivalent present (live `e4495b60`).
- **value for toolpaths:** none.
- **conflicts:** none.
- **should it be brought over:** no.

### commit ebf0dcbd feat(bench): expose arrangement size and exact-rational digit growth | files 5 | read full
- **generation relevance:** n (instrumentation).
- **on live branch?:** equivalent present: `src/compas_cgal/stock.py:339-365`, `src/stock_2.h:245-263`, `benchmarks/instrument.py`.
- **value for toolpaths:** none.
- **conflicts:** a textual conflict in `stock_2.cpp` if cherry-picked.
- **should it be brought over:** no.

### commit 83b36ab1 build: pixi manifest for jf/toolpath-redesign | files 4 | read full (pixi.lock not read)
- **generation relevance:** n.
- **on live branch?:** superseded by live's own `[tool.pixi]` tables (`pyproject.toml:175-233`).
- **value for toolpaths:** none.
- **conflicts:** `pyproject.toml` and `pixi.lock` throughout.
- **should it be brought over:** no.

### commit 8b1d431e fix(build): scope strict-config to pixi, dedup probe, widen lint gate | files 3 | read full
- **generation relevance:** n.
- **on live branch?:** missing.
  - Live has `strict-config = false` globally (`pyproject.toml:78`).
  - Live tasks still use the inline `dir="$(…)" && pytest` pattern (`pyproject.toml:224-232`), which this commit reports, as measured, does not short-circuit under pixi's shell (relayed, unverified).
- **value for toolpaths:** none. Build hygiene.
- **conflicts:** `pyproject.toml` tasks.
- **should it be brought over:** partly. Port `scripts/pytest_editable.py` and the strict-config override as a build task, outside generation scope.

### commit e333c46f build: scope lint gate to owned code, green it | files 3 | read full
- **generation relevance:** n.
- **on live branch?:** superseded (live lint task differs, `pyproject.toml:233`).
- **value for toolpaths:** none.
- **conflicts:** `pyproject.toml`.
- **should it be brought over:** no.

### commit 2af221cd docs(build): correct strtobool claim in strict-config override comment | files 1 | read full
- **generation relevance:** n.
- **on live branch?:** n/a; it depends on `8b1d431e`.
- **value for toolpaths:** none.
- **conflicts:** none.
- **should it be brought over:** only together with `8b1d431e`.

### commit d90e18b8 fix(toolpath): mat_scale bounded to (0,1] — gouge-free holds by construction only there | files 4 | read full
- **generation relevance:** y. A trochoid circle's gouge-freedom rests on `mat_scale ≤ 1` (the circles are never certified). jf's spec C1 measured +2.25 wall penetration at 1.5 (relayed).
- **on live branch?:** missing. Live `src/toolpath.cpp:600-613` validates neither `mat_scale` nor `radial_clearance`.
- **value for toolpaths:** an essential safety contract for the guide every live generator walks. Quality: adequate, a one-line seam check plus parametrised gouge tests.
- **conflicts:** none. `toolpath.cpp` is untouched on live.
- **should it be brought over:** yes.

### commit 2568453b fix(pixi): _editable-rebuild actually rebuilds | files 1 | read full
- **generation relevance:** n.
- **on live branch?:** superseded. Live sets `editable.rebuild = true` (`pyproject.toml:76`).
- **value for toolpaths:** none.
- **conflicts:** `pyproject.toml`.
- **should it be brought over:** no.

### commit b0094cb6 fix(toolpath): certify the implied traverse when link_paths is off | files 4 | read full
- **generation relevance:** y (linking safety).
- **on live branch?:** missing. In live `src/toolpath.cpp:1072-1073`, the `link_paths=false` branch moves `cur_xy` without any check.
- **value for toolpaths:** closes an uncertified through-material traverse in the polyline. Quality: adequate.
- **conflicts: semantic, serious.** `engagement_controlled_toolpath` and `radius_regulated_toolpath` read their guide by calling `trochoidal_mat_toolpath_circular(link_paths=False, no clearance_z)` (`src/compas_cgal/engagement_toolpath.py:438-452`). With this commit alone, any guide whose consecutive chains have a wall-crossing straight join would **raise**, and the gated generators would fail on such pockets. It has not been measured which gate pockets trigger this.
- **should it be brought over:** partly. Bring it only with `446daca5`, and make the gated guide call pass a `clearance_z` (or certify the join), with a gate run before and after.

### commit 446daca5 fix(toolpath): lift the unlinked traverse to clearance_z | files 4 | read full
- **generation relevance:** y. With `clearance_z`, the unlinked traverse now rides a retract/plunge pair.
- **on live branch?:** missing (same anchor, `src/toolpath.cpp:1072`).
- **value for toolpaths:** correct air move. It adds a retract and plunge per chain boundary (13-path pocket: 1+1 → 12+12, relayed). That is honest, but more air time.
- **conflicts:** as `b0094cb6`.
- **should it be brought over:** yes, together with `b0094cb6` and `b2fa087e`.

### commit b2fa087e fix(toolpath): tessellator must carry a discontinuous op stream | files 2 | read full
- **generation relevance:** y. The emitted polyline is what a post-processor would read. Live drops every non-first op's start point, so a retract/link/plunge gap becomes one diagonal through the wall (11 ramps in jf's measurement, relayed).
- **on live branch?:** missing (`src/toolpath.cpp:825-870`).
- **value for toolpaths:** essential for polyline consumers. Quality: adequate.
- **conflicts:** none.
- **should it be brought over:** yes, with `34bddcc0`.

### commit 34bddcc0 fix(toolpath): dedup tessellator junctions at OUTPUT_DEDUP_TOL | files 2 | read full
- **generation relevance:** y (polyline output).
- **on live branch?:** missing. `OUTPUT_DEDUP_TOL` exists (`src/toolpath.cpp:47`) but is not used in `tessellate_operations`.
- **value for toolpaths:** removes ~1e-16 segments. It uses an existing named reporting tolerance, which is acceptable for output (not a decision).
- **conflicts:** none.
- **should it be brought over:** yes, with `b2fa087e`.

### commit 2f983d8b fix(engagement): validate tool radius and coordinate finiteness at the exact seam | files 3 | read diff of src; tests skimmed
- **generation relevance:** n (input contracts).
- **on live branch?:** superseded by live `validate_engagement_input_binary64` and `validate_segment_certification_input_binary64` (`src/audit_exact_station_2.cpp:248-300`, called at `src/engagement_2.cpp:525`).
- **value for toolpaths:** none further.
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** no.

### commit 832765b0 fix(stock): seal the exact-kernel seam on every Stock2 entry point | files 5 | read src diff; tests skimmed
- **generation relevance:** n directly. Robustness: `subtract_disk(c, +inf)` empties the stock.
- **on live branch?:** missing for the legacy entry points. Live `src/stock_2.cpp:354-357` (`subtract_disk`), `:516-518` (`subtract_capsule`) and `:668-672` (`subtract_arc_sweep`) check only `radius <= 0`. The newer exact entry points do have finiteness errors (e.g. `:585-588`).
- **value for toolpaths:** low. Generators pass finite values.
- **conflicts:** `stock_2.cpp`/`.h`, heavily rewritten on live (traits owner, `set()` accessor).
- **should it be brought over:** partly. Re-implement the guards on live's code; do not cherry-pick.

### commit 5c489fbd fix(stock): bound the disk-chain count; finiteness alone does not | files 3 | read src diff; tests skimmed
- **generation relevance:** n. A runaway-allocation guard: an `int` cast that saturates at a long move or a tiny radius.
- **on live branch?:** missing. `static_cast<int>(std::ceil(len / spacing))` is unbounded at `src/stock_2.cpp:532` and `:711`.
- **value for toolpaths:** low. `n = 50·len/r` stays small for real moves (relayed derivation).
- **conflicts:** `stock_2.cpp`.
- **should it be brought over:** partly. Re-implement on live.

### commit 92cedda2 fix(toolpath): contracts for radial_clearance, radius bounds, clearance plane | files 3 | read full
- **generation relevance:** y.
  - Negative `radial_clearance` licenses cutting past the wall.
  - An inverted min/max trochoid radius silently discards the cap.
  - `clearance_z ≤ cut_z` silently degrades to flat linking.
- **on live branch?:** missing (`src/toolpath.cpp:600-613`, `:996`). Live has no `src/exact_boundary.h`.
- **value for toolpaths:** a real safety contract on the guide generator. Quality: adequate.
- **conflicts:** it introduces `exact_boundary.h` (also created by `832765b0`). Cherry-pick the order carefully.
- **should it be brought over:** yes.

### commit 728b6c82 refactor(engagement): single compiled source of truth for the TEA growth guard | files 4 | read full
- **generation relevance:** n (audit).
- **on live branch?:** missing. Live still has a Python mirror with its own `TEA_GUARD_SAFETY_FACTOR` (`src/compas_cgal/engagement.py:118-145`).
- **value for toolpaths:** none. It removes drift risk in the audit.
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** no. It is moot if the growth-lemma path is retired in favour of an exact audit.

### commit ebb627d3 test(engagement): falsification harness for the TEA growth bound (RED) | files 3 | read full test diff
- **generation relevance:** n (validation). It shows the lemma under-bounds growth for every concave void ρ > r, by up to 24× (relayed from the message table).
- **on live branch?:** superseded by a trimmed `tests/test_growth_bound.py` (173 lines vs jf's 953+; live `7cda8f6d`).
- **value for toolpaths:** none for generation. Strong evidence for the audit.
- **conflicts:** same filename on live.
- **should it be brought over:** no. Live already keeps a version.

### commit eee1c35b test(engagement): false TEA certificate FOUND -- rib of the tool's own radius | files 2 | read test names and asserts
- **generation relevance:** n.
- **on live branch?:** equivalent present as non-certifying negative controls (`tests/test_false_certificate.py:151-225`).
- **value for toolpaths:** none.
- **conflicts:** same file.
- **should it be brought over:** no.

### commit 3c270cb5 test(engagement): spiral rib -- the false certificate is generic, not a radius coincidence | files 1 | read test names and asserts
- **generation relevance:** n. Note the generation implication: any spiral or ramped contour pass leaves a strip whose local curvature radius equals r (relayed).
- **on live branch?:** equivalent present (`tests/test_false_certificate.py:182`).
- **value for toolpaths:** none.
- **conflicts:** same file.
- **should it be brought over:** no.

### commit c3bdef84 test(engagement): separate witness liveness from certifier verdict in every red | files 1 | read asserts
- **generation relevance:** n.
- **on live branch?:** equivalent spirit present (the live controls assert liveness only).
- **value for toolpaths:** none.
- **conflicts:** same file.
- **should it be brought over:** no.

### commit 2e2422b5 fix(engagement): exact swept-annulus guard replaces the false growth lemma | files 3 | read src diff full; tests not read
- **generation relevance:** n (audit). It is not usable to drive generation: it saturates whenever a material component wraps the station, and it cannot certify a corridor centreline (jf doc).
- **on live branch?:** missing. Live `certify_segment_tea` still rests on `tea_growth_bound` (`src/engagement_2.cpp:407-499`).
- **value for toolpaths:** none for generation. For the audit it is a sound but conservative bound. Weak against the exactness doctrine:
  - adjacency is decided by `Aabb::near` on doubles with a 1e-12 margin;
  - the angular bound uses doubles with 1e-9 slack.

  The fail-safe direction is documented, but these are double decisions inside a certificate.
- **conflicts:** `engagement_2.cpp`, heavily rewritten on live.
- **should it be brought over:** no. Retire the legacy audit for live's exact `continuous_tea_2` / `audit_certification_2` path instead. Until then, label `cap_certified` from `audit_toolpath_engagement` as unsound (the live `engagement_2.cpp:399` "QED" comment is contradicted by live's own `tests/test_false_certificate.py`).

### commit 03866fcc fix(engagement): rotation-invariant interior bound -- exact direction spans replace the axis-aligned box | files 2 | read src code lines; tests not read
- **generation relevance:** n.
- **on live branch?:** missing.
- **value for toolpaths:** none. It tightens the jf audit bound (exact station-enclosure test, union of direction spans).
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** no (as above).

### commit e446d035 test(engagement): move the stations==1 witness from an assertion to the docstring | files 1 | read assert lines only
- **generation relevance:** n.
- **on live branch?:** n/a.
- **value for toolpaths:** none.
- **conflicts:** test file.
- **should it be brought over:** no.

### commit 3db92d64 fix(engagement): sub-ulp rim sliver was reported as a full turn | files 1 | read code lines
- **generation relevance:** n (reporting).
- **on live branch?:** superseded. Live `rim_sub_arc_span` uses the chord half-angle `2·atan2(|b−a|, |a+b|)` (`src/engagement_2.cpp:80-85`) plus a full-turn ceiling (`:52`, `:553`).
- **value for toolpaths:** none.
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** no.

### commit f227dab2 docs: the TEA certificate -- the false certificate and the swept-annulus guard | files 3 | read partial (intro, status, references of docs/engagement_certificate.md)
- **generation relevance:** n (docs-only).
- **on live branch?:** missing (`docs/engagement_certificate.md` absent).
- **value for toolpaths:** records why the growth lemma is false, and the corridor-centreline blind spot of the replacement.
- **conflicts:** none.
- **should it be brought over:** partly. Bring the false-certificate record; slice d1.

### commit 387bb7ef docs: harvest 2pi misreport is repaired at 3db92d6, not open | files 1 | read full
- **generation relevance:** n.
- **on live branch?:** n/a (page missing).
- **value for toolpaths:** none.
- **conflicts:** none.
- **should it be brought over:** with `f227dab2` only.

### commit 5c1d492a docs: ETH audit findings + 21-task remediation plan | files 2 | read partial (findings spec full; plan headings + tasks 13–19 full)
- **generation relevance:** y in the findings that touch `toolpath.cpp`:
  - I7: `Boundary::segment_clear` takes gouge decisions on Epick doubles.
  - I8: `certify_bridges` can exit on an untested refined run.
  - I9: hole edges crossing a concave outer, and nested holes, reach `create_interior_straight_skeleton_2` (UB).
  - Note: `external_tangent` unpacks to doubles.

  All are in the guide the gated generators walk. Docs-only.
- **on live branch?:** missing. Live `toolpath.cpp` is unchanged since the fork, so I7/I8/I9 stand (verified file untouched; the findings themselves are relayed).
- **value for toolpaths:** a ready, test-first plan for guide hygiene.
- **conflicts:** none.
- **should it be brought over:** partly. Turn tasks 13, 14, 16, 17 and 19 into live tasks; slice d1 extracts the findings.

### commit 148a49ad fix(engagement): enforce the shared-root precondition in release builds | files 2 | read src diff
- **generation relevance:** n.
- **on live branch?:** equivalent present: the runtime root check in `src/audit_exact_station_2.cpp:335`.
- **value for toolpaths:** none further.
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** no.

### commit 48b9e08f fix(engagement): correct the slack derivation, guard the <=pi precondition | files 2 | read src diff
- **generation relevance:** n (reporting). It adds a guard that the cutter circle splits into exactly two x-monotone pieces, and throws on read-out collapse.
- **on live branch?:** partly superseded. Live uses a different span formula and a 1e-12 absolute slack (`src/engagement_2.cpp:52`). jf's analysis says an absolute slack is wrong at large coordinate/radius ratios (relayed).
- **value for toolpaths:** none.
- **conflicts:** `engagement_2.cpp`.
- **should it be brought over:** partly. Check live's 1e-12 constant against jf's `ulp(scale)/r` model (wedge measurement in `73d53729`'s comments); reporting only.

### commit c90a4fef docs: machining-quality track (tasks 29-35), pixi build rules | files 3 | read partial (CONVENTIONS.md delta and the remediation-plan delta in full; corpus plan headings only, live has a superset)
- **generation relevance:** y (docs-only). Tasks 31 (declined corner residue), 34 (generalise MAT beyond the L fixture), 36 (helical/ramped entry) and 37 (break the ~141.8° cap saturation) are the generation agenda. It records two instrument defects:
  - `max_loop_radius_step` compared radii across a retract;
  - the engagement angle is the rim arc, twice the textbook angle.
- **on live branch?:**
  - Second gate cap: present (`benchmarks/gate.py:55-59`).
  - Corpus plan: present (live superset, 3248 lines).
  - Tasks 31, 34, 36, 37: open on live. 34 is verified (`require_l_shape_input` still at `src/segment_site_catalog_graph.cpp:1220`).
- **value for toolpaths:** high as a diagnosis record, especially the Task 37 analysis: (1) the station set, (2) the acceptance rule and (3) non-monotonicity are the candidate causes, to be measured.
- **conflicts:** CONVENTIONS.md wording.
- **should it be brought over:** partly. Fold tasks 31/34/36/37 into the live plan's Stage B/C; slice d1.

### commit f369b3df docs: due-diligence reviews, integration decision, fig6 plots | files 8 | read partial (INTEGRATECODEX.md and MOTIVATION.md full; review.md risks 0b–5; t9 review verdict; SVG/PDF not read)
- **generation relevance:** y (docs-only). INTEGRATECODEX ranks the critical path to Held parity:
  1. break the ~141.8° saturation;
  2. helical entry;
  3. close the residue (0.28–0.80% vs 0.152%);
  4. kill degenerate and redundant loops.

  `review.md` 0b: "the test suite does not encode toolpath quality".
- **on live branch?:** missing (no INTEGRATECODEX/review files on live). The recommended rebase of jf onto codex was never done: every jf fix listed as "clean" is still missing.
- **value for toolpaths:** the ranking matches the brief's root causes and adds the saturation item, which the live plan parks in Stage C.
- **conflicts:** none.
- **should it be brought over:** partly. Bring the INTEGRATECODEX "clean cherry-pick" list as a task (the toolpath.cpp fixes); slice d1.

### commit 73d53729 feat(engagement): adaptive arc certification on the swept-annulus guard | files 5 | read partial (src + engagement.py diff ~520 changed lines; tests not read)
- **generation relevance:** n (audit). Per its message it is incomplete and unreviewed; the certified fraction of generator circles was never measured (0 of 126 before).
- **on live branch?:** missing. Live arc certification is still the fixed-density Python mirror on the false lemma (`src/compas_cgal/engagement.py:226-306`). That mirror refuses every circular cut at 120° before measuring (t9 review, relayed).
- **value for toolpaths:** none for generation.
  - The turn-direction cross-check against `point_at(0.5)` is a good contract.
  - Weak: decisions use double angles (`a0`, `a1`, cos/sin station placement).
- **conflicts:** `engagement_2.cpp`, `engagement.py`, the tests.
- **should it be brought over:** no. Route the arc audit through live's exact circle oracle (`continuous_tea_2/circle_oracle`) instead.

## codex/held-reference-corpus

### commit 07c87ce0 feat(toolpath): add Figure 5 Voronoi guide | files 6 | read full
- **generation relevance:** y (guide). `SegmentSiteGuide2::build` builds the CGAL segment-Delaunay/Voronoi graph of any polygon with holes. It samples lines and parabolas (sagitta-bounded refinement), keeps samples whose exact point is at least r from every edge, and returns nodes, edges and generator-site provenance.
- **on live branch?:** missing. Live's exact MAT is L-fixture-only (`src/segment_site_catalog_graph.cpp:1220`).
- **value for toolpaths:** the only general MAT guide in the repo.
  - Quality: adequate as a guide; weak against the exactness doctrine. Sampled points merge at 8 ε (`REPORTING_POINT_MERGE_EPSILON_FACTOR`), the in-domain filter runs on rounded samples, and edges straddling the domain boundary are dropped rather than clipped, so branch ends near r are lost.
  - It needs `f90d1ce4` (bogus branches).
- **conflicts:** `CMakeLists.txt` and `src/medial_axis_2.cpp` (live changed CMake since base).
- **should it be brought over:** yes, as the general-polygon MAT guide candidate for the plan's guide decision. Test it on the gate pockets against the straight skeleton.

### commit 1dc71b8e docs(plan): focus continuous Figure 5 path | files 1 | read full
- **generation relevance:** n (plan text).
- **on live branch?:** missing (live plan file diverged).
- **value for toolpaths:** none.
- **conflicts:** `docs/superpowers/plans/2026-09-05-held-paper-figure-reproduction.md`.
- **should it be brought over:** no.

### commit 266bfec6 feat(toolpath): render independent Figure 5 path | files 4 | read full
- **generation relevance:** y (placement + ordering prototype). Paper circles from the sampled Voronoi guide are snapped to the exact cutter-centre boundary cycle (`_coverage_2.build_center_boundary_cycle`) and ordered by boundary progress. Successors come from `select_next_standard_candidate`, and connectors run along the exact offset cycle.
- **on live branch?:** missing.
- **value for toolpaths:** superseded within the branch by the analytic path (`757b948c`).
  - Weak: many reporting-epsilon (64 ulp) equality filters in placement decisions.
  - Weak: no successor within cap triggers an "entry" restart, so it is not one stream.
- **conflicts:** none (new files).
- **should it be brought over:** partly. Keep it only as the polygon-input path that pairs with `07c87ce0`; it is the only Held-style placement that runs on arbitrary polygons.

### commit 9db59ac5 feat(benchmarks): measure Figure 6 planar paths | files 2 | read full
- **generation relevance:** n (measurement). It measures planar length (circles + connectors only) on Figure 6 axes.
- **on live branch?:** missing. Live `benchmarks/held_figure6_comparison.py` measures the live controlled generator, with compliance unaudited.
- **value for toolpaths:** a like-for-like length metric.
- **conflicts:** none.
- **should it be brought over:** yes, with `0a3cd8a8`.

### commit f90d1ce4 fix(toolpath): remove false medial branches | files 3 | read full
- **generation relevance:** y. It drops Voronoi edges between a vertex site and its own incident segment, which are not medial-axis points.
- **on live branch?:** missing.
- **value for toolpaths:** a required correctness fix for `07c87ce0`.
- **conflicts:** none.
- **should it be brought over:** yes, with `07c87ce0`.

### commit 334c028b feat(toolpath): evaluate continuous MAT edges | files 2 | read full
- **generation relevance:** y. It evaluates point/segment guide edges on the true parabola rather than the chord.
- **on live branch?:** missing.
- **value for toolpaths:** adequate (double evaluation; the focal-length = 0 check is exact equality).
- **conflicts:** none.
- **should it be brought over:** yes, with `07c87ce0`.

### commit 757b948c feat(toolpath): reproduce analytic Figure 5 path | files 5 | read full (png excluded)
- **generation relevance:** y (**core**). This is Held's standard algorithm:
  - walk the r-offset of the line/arc boundary;
  - at each progress, the first equidistant boundary contact on the inward normal gives the MAT point m;
  - the circle is `(d−r)/2` between q and m;
  - the successor is bracketed and bisected to exactly the cap.
- **on live branch?:** missing. Live's Figure 5 paths need the publisher-derived start and raw guide (plan:38).
- **value for toolpaths:** high. It is the first publisher-free continuous Held standard path. This commit alone still has a wrong locus at joins, a storage-seam start, and bracketing that can skip a first crossing; all are fixed in later commits.
- **conflicts:** none (new file).
- **should it be brought over:** yes, with `5b5ed867`, `4685653b`, `567da10d` and `274fd1a4` as one unit.

### commit 0c792fc7 docs(plan): reopen Figure 5 after review | files 1 | read full
- **generation relevance:** n. It records review defects: join gaps up to 0.001852 mm, a storage-seam start/termination, and full circles supplied only by the renderer.
- **on live branch?:** n/a.
- **value for toolpaths:** an honest reopening.
- **conflicts:** plan file.
- **should it be brought over:** no.

### commit 5b5ed867 fix(toolpath): make Figure 5 path continuous | files 8 | read full
- **generation relevance:** y. It:
  - trims convex offset joins to the support intersection;
  - inserts concave joins as tool-radius arcs about the boundary vertex;
  - starts at a geometric south-east support point, invariant to primitive rotation (test);
  - adds a cyclic traversal that terminates on re-entering the initial clearance disk;
  - emits explicit `motions` (full CCW circles interleaved with connectors).
- **on live branch?:** missing.
- **value for toolpaths:** high: one continuous executable planar stream.
  - Weak: it removes the nearest-boundary validation of the MAT candidate (`min(d > tool)` instead). Correct on smooth Figure 5; risky on sharp convex corners (summary finding 1, unverified).
- **conflicts:** none.
- **should it be brought over:** yes (unit as above).

### commit 22ffc356 docs(plan): record Figure 5 review gaps | files 1 | read full
- **generation relevance:** n. It records that fixed-step bracketing skipped a first engagement crossing at 68.75°, and that the terminal overshot.
- **on live branch?:** n/a.
- **value for toolpaths:** none.
- **conflicts:** plan file.
- **should it be brought over:** no.

### commit 4685653b fix(toolpath): resolve Figure 5 review gaps | files 8 | read full
- **generation relevance:** y. It:
  - makes concave joins first-class offset pieces with a locus (circles placed along them);
  - resolves locus ownership at exact joins by maximum support slack;
  - replaces fixed bracketing with per-cell golden-section maximum search then bisection;
  - trims the terminal connector at the first entry into the initial clearance circle;
  - adds sampled Hausdorff transition residuals.
- **on live branch?:** missing.
- **value for toolpaths:** high. Still weak: golden-section search assumes unimodal engagement inside a 0.25·r cell, which is not a guarantee (the memory note says "engagement is not monotone"). A narrow peak plus a second crossing in one cell is still possible.
- **conflicts:** none.
- **should it be brought over:** yes (unit).

### commit a662a9c3 fix(toolpath): decode publisher turns geometrically | files 6 | read full (json excluded)
- **generation relevance:** n (comparison evidence). Publisher turns are 278, not 265: 17 turns are five-cubic circles.
- **on live branch?:** missing. Live has `FIGURE5_PUBLISHER_TURN_COUNT = 265` (`benchmarks/held_figure5_publisher.py:33`) and the open-turn range 123–136 (`:198`).
- **value for toolpaths:** corrects every live publisher-relative count and residual.
- **conflicts:** `benchmarks/held_figure5_publisher.py` (live changed it 8 lines) and the evidence JSON.
- **should it be brought over:** yes.

### commit 0a3cd8a8 feat(toolpath): measure Figure 6 motion stream | files 2 | read full
- **generation relevance:** n (measurement). Planar length is summed over the explicit continuous `motions` stream, including the terminal connector.
- **on live branch?:** missing.
- **value for toolpaths:** a correct metric contract (it rejects a discontinuous stream).
- **conflicts:** none.
- **should it be brought over:** yes.

### commit 567da10d fix(toolpath): close Figure 5 review findings | files 12 | read full (json/png excluded)
- **generation relevance:** y.
  - Terminal when no further crossing: the path runs to the end of the lap (valid 120° sweeps).
  - Upstream G1 repair of the Figure 5 line/cubic seam removes a radius-to-zero MAT branch, with an outward-rounded perturbation bound (`0.07380 ≤ 0.07386` mm, test).
  - A tracked raw publisher operator stream (`figure5_publisher_source.json`) and ordered-transition RMS.
- **on live branch?:** missing. Live has no `_figure5_g1_sources`; live `held_reference_cases.py` changed 12 lines since base.
- **value for toolpaths:** the G1 repair matters for any Figure 5 guide, live's included (unverified on live).
- **conflicts:** `benchmarks/held_reference_cases.py`, `figure5.json` (reference data), `held_figure5_publisher.py`.
- **should it be brought over:** yes, after checking live's Figure 5 guides for the spurious branch.

### commit 26c04d3f test(toolpath): enforce overlay freshness | files 1 | read full
- **generation relevance:** n.
- **on live branch?:** n/a.
- **value for toolpaths:** none.
- **conflicts:** none.
- **should it be brought over:** with the unit only.

### commit 4bf94eed docs(plan): close Figure 5 toolpath | files 1 | read full
- **generation relevance:** n. It claims a five-perspective zero-finding review after 114 focused tests (relayed, unverified).
- **on live branch?:** n/a.
- **value for toolpaths:** none.
- **conflicts:** plan file.
- **should it be brought over:** no.

### commit 7be75e0d feat(toolpath): replay engagement within motions | files 12 | read full
- **generation relevance:** y (instrument for generation decisions). Swept-prefix station classification for lines (rear half-disk covered) and arcs (start disk ∪ angular sector covered). Clipping is exact on one-root points, with rational witnesses, and the motion is depleted afterwards. A sequential replay runs at a 0.5·r pitch with a density-doubling residual.
- **on live branch?:** missing. Live has only the π-cap swept-prefix theorem for segments (`src/continuous_tea_2/segment_oracle.h:50-57`).
- **value for toolpaths:** high as a measurement of in-motion engagement (entries, stay-down links, connectors). Weak:
  - the arc `covered` predicate branches on a double `swept_angle` (`<= pi`, `== 2pi`);
  - arc station centres come from cos/sin;
  - sampled stations are not a certificate (it says so).
- **conflicts: substantial.** `src/engagement_2.cpp`, `src/audit_exact_station_2.cpp` and `audit_certification_2.h` were rewritten on live in the same regions (live −130/−122 lines vs base). `CMakeLists.txt` and `_stock_2.pyi` also conflict.
- **should it be brought over:** yes, re-applied by hand onto live's `audit_exact_station_2` (the `classify_audit_filtered_station_exact` seam).

### commit 8ca2c9f6 feat(toolpath): add contour-aware Figure 5 path | files 3 | read full
- **generation relevance:** y (**core**). Held §3.1: it keeps the exposed arcs of the union of accepted clearance disks and moves the predecessor entry point b clockwise onto that contour. Candidate queries are pure; acceptance is transactional. It yields fewer circles and a shorter path than the standard path (test).
- **on live branch?:** missing.
- **value for toolpaths:** high: a better path at the same cap.
  - Weak: all-double circle intersections with a 64-ulp angular slack in `contains`.
  - O(n²) contour maintenance (stated).
  - Engagement is still the predecessor model, not the stock.
- **conflicts:** `benchmarks/held_standard_placement.py` (live changed it 66 lines since base). **Not importable at this commit**: it imports `build_analytic_figure5_path_with_policy`, `Figure5CircleAtProgress` and `Figure5TraversalIntervals`, which first exist in `274fd1a4` (verified with `git show 8ca2c9f6:benchmarks/held_figure5_analytic.py`).
- **should it be brought over:** yes, squashed with `274fd1a4`.

### commit 274fd1a4 feat(toolpath): add shared MATHSM traversal | files 8 | read full
- **generation relevance:** y. Placement becomes a policy over one traversal (`Figure5PlacementPolicy`, one-sided interval ownership at piece joins). MATHSM constant Euclidean centre spacing runs on the same locus (12 spacings 0.025–0.6 D; 1313 circles at 0.075 D, stable under density halving). The commit also excludes adjacent primitives of the same logical source as MAT second contacts.
- **on live branch?:** missing. Live MATHSM uses a different generator (`benchmarks/mathsm.py`).
- **value for toolpaths:** high for honest baselines.
  - Weak: Figure 5 declares all 31 primitives one logical source (`held_reference_cases.py`). That is a pocket-specific declaration, and it changed the 68.75° witness result (advance 0.5006 → 1.8438).
- **conflicts:** `ReferenceBoundary` gets a new field `logical_source_indices` (`held_reference_geometry.py`), an API change for live's many `ReferenceBoundary` users; `benchmarks/units.py`.
- **should it be brought over:** yes, squashed with `8ca2c9f6`. Audit every live `ReferenceBoundary` construction site.

### commit a29ff8c6 feat(toolpath): sweep standard Figure 6 paths | files 2 | read full
- **generation relevance:** n (measurement). It sweeps the standard path over the publisher Figure 6 caps (0–200 by 20). Out-of-domain and unresolved caps are typed failures; 180° fails with "without leaving its initial clearance disk". Compliance comes from the sampled swept replay.
- **on live branch?:** missing.
- **value for toolpaths:** a like-for-like Figure 6 curve for the standard path. Contour-aware and MATHSM sweeps are not wired here.
- **conflicts:** depends on live's `benchmarks/held_figure6_publisher.py` (present).
- **should it be brought over:** yes.

## codex/stock-clone-ownership

The 20 other commits listed under this branch are the `codex/held-reference-corpus` commits recorded above (same SHAs).

### commit d0f32e8d fix(stock): retain cloned arrangement traits | files 5 | read full
- **generation relevance:** n directly. Stock clone and transaction lifetime: a crash in point location after the source died.
- **on live branch?:** superseded by `57482731` + `8d2b76cc`, a different design:
  - Live resolves and pins the traits-owning Gps (`src/stock_2.cpp:188-207`, `:279-292`, `:306-325`).
  - It raises `StockTraitsUnownedError` when no candidate owns the traits.
  - Live's exact transactions copy and resolve through `replace_set` (`:734-772`).
- **value for toolpaths:** equal correctness.
  - `d0f32e8d` is simpler and self-owning: no introspection of CGAL's borrowed pointer, no ancestor arrangement kept alive. Its own bounded timing shows a slower clone (0.725→0.784 ms with 49 holes; 236.5→244.3 ms at 6,908 vertices; relayed).
  - Live's approach can keep a whole ancestor Gps alive per clone family (memory, not correctness).
  - `tools/rebuild_stock_isolated.py` is a recovery tool and is not generation-relevant.
- **conflicts:** `Stock2::clone` and the three `subtract_exact_*` transactions in `src/stock_2.cpp`; the same new file `tests/test_stock_clone_lifetime.py`, with different contents on both sides. Its tests import `engagement_arc_prefixes`, which exists only with `7be75e0d`.
- **should it be brought over:** no. Keep live's fix. Partly: port its `docs/exactness.md` paragraph (the borrowed-traits mechanism) and its three-arc regression scenario as an extra test against live's implementation.
