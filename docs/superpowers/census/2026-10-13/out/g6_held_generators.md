# Census slice g6_held_generators

Measurements marked `measured` were run read-only with `pixi run python -c ...` in this worktree on 2026-10-13 (no files written). Plan = `docs/superpowers/plans/2026-10-13-product-gate-tree-walk.md`.

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **The reproduction is far from Held's own path. It is 2.8x longer and has 7.6x as many circles.** Held's published Figure 5(a) stream has 265 turns, about 8,831 mm of turn length plus about 541 mm of connector chords (chord lower bound). Its smallest turn radius is 1.167 mm (> r = 1), and no turn has rho <= r (measured from `benchmarks/data/held_pfeiffer_2025/figure5_publisher_turns.json` via `held_figure5_publisher.load_figure5_publisher_path`). The repository's `figure5_refined` path has 2,013 circles and 26,436.6 mm (`benchmarks/results/held_toolpath_baseline.json`, verified). 1,059 of those circles have rho <= r, 495 have rho <= 0.1 r, and the smallest is 0.00046 mm (measured). The plan calls this "the most complete toolpath result", and it is. Measured against Held, though, it is not a good toolpath yet. **Contradicts** the brief's premise that these are near world class.
2. **The raw guide is not publisher-derived. It is the product's own straight-skeleton guide.** `held_figure5_raw_guide.py:417-438` calls `compas_cgal.engagement_toolpath._guide_chains(spec.polygon, ...)`. Only the **start marker** and the **evidence tolerances** are publisher-derived. **Contradicts** plan lines 33 and 1374ff ("needs the publisher-derived start circle and raw guide").
3. **A native, exact, publisher-free Held circle construction already exists and runs on all three gate pockets** (measured). `_coverage_2.build_center_boundary_cycle` plus `_coverage_2.boundary_circle_at_contact` give, at any offset contact q, the exact medial contact m, a centre at the midpoint of q and m, and rho = |qm|/2. Convex offset vertices come back as `is_stationary` with rho = 0. The L_shape reflex corner comes back as an offset arc primitive. On rect_12x8 (r = 1) rho ranges from 0 to 1.5. Used only by the diagnostic `held_native_boundary_draft.py:132-144` and `held_native_curved_circle_plot.py:136`, never by a generator.
4. **The gate's degenerate-loop rule structurally excludes Held placement.** `benchmarks/quality.py:87` sets `DEGENERATE_LOOP_RATIO = 1.0`, so every loop with tool-centre radius <= r counts (`quality_observations.py:900`). On Figure 5, 1,059 of 2,013 reproduced circles would count (measured). On rect_12x8, Held circles have rho <= 1.5 r, so a large share would count too (measured, see finding 3). Task 0b's decision gate would therefore be decided by this rule, not by path quality. It needs a ruling before Task 0b's numbers mean anything.
5. **Spacing decisions are made in a different order from emission.** `build_hypothesis_figure5_path` places with the standard model per *lane* (one guide run x one boundary side, `held_figure5_path.py:393-436`), then re-sorts the selection into boundary order (`:444-454`, `held_figure5_boundary_path.py:279-301`). The emitted predecessor is therefore often not the placement predecessor, and `refine_figure5_engagement` has to insert circles: 1,532 grow to 2,013 (measured). `held_boundary_order_spacing.py` fixes this by selecting in emission order (`:51-65`), which is why Figure 8 upper drops from 4,183 to 1,684 circles (baseline, verified).
6. **The contour-aware spacing is about half the length but loses coverage.** The `held_contour_bound_spacing.py` draft is recorded at 468 circles and 13,087 mm on Figure 5, against 26,437 mm. Preserving the dense source's coverage needs 1,558 circles and 23,823 mm (+82%) (`docs/assets/images/held_figure5_coverage_preserving.json`, verified). The CLI labels coverage-preserving deletion "the rejected dense-source preservation experiment" (`tools/held_contour_bound_toolpaths.py:37`).
7. **No entry, anywhere.** The first circle is "the declared cleared seed" (`held_contour_bound_spacing.py:57-58`) or simply the first sweep (`held_coverage_thinning.py:72-73`). This is the same open D1 as the plan.
8. **Convex corners are handled by shrinking rho toward 0** (`held_figure5_engagement_refinement.py:51-74`, applied to every circle). That is the "bore it, degenerate loops" variant the plan's B2 already rejects, so it does **not** solve D2.
9. **The predecessor-only 79.99 deg is not the whole-history bound.** The native whole-history leading-rim bound (`_stock_2.HeldDiskContour2.engagement_bound`) fails to certify 12 successors of the 2,013-circle path, with a maximum bound of 138.6 deg (measured). Both are models. The bound is conservative, so these 12 are unproven, not shown to exceed. Connector engagement is audited nowhere (`held_contour_bound_spacing.py:3-5`).
10. **The "containment rejections" in the baseline are tangency rounding, not gouges.** There are 1,010 of 2,013 (baseline), but the maximum protrusion is 6.6e-15 mm (measured). Every circle is tangent to the r-offset at q, so the exact predicate on binary64-rounded centres splits roughly 50/50. An adapter or gate gouge check on these floats will show the same noise.

### Generation capabilities here that the gated generators do NOT use (verified: no `src/compas_cgal/*.py` imports them, grep)

- Boundary-order single-lap traversal: one CCW climb lap; each loop tangent to the wall offset at q and reaching the medial axis; connectors along the wall offset; zero retracts and zero links by construction (`held_figure5_boundary_path.py:251-304`).
- Exact native Held circles from an offset contact (`_coverage_2.boundary_circle_at_contact`, `NativeBoundary2.circle_on_piece`). This includes curved boundaries and exact reflex-corner arcs in the tool-centre offset (`_coverage_2.build_center_boundary_cycle`, `cycle.ccw_transition`).
- Held Section 2.3 standard-model engagement, with CGAL-owned corrected chord (`held_standard_placement.py:161-222`).
- Whole-history filled-disk contour with an exact squared-chord cap surrogate (`_stock_2.HeldDiskContour2.engagement_bound/append/contact_toward`), used in `held_contour_bound_spacing.py:72-78`.
- An exact deletion predicate that preserves coverage (`_stock_2.Stock2.can_remove_circle`, `held_coverage_thinning.py:93`).
- Exact residual stock after full circles plus capsule connectors (`_coverage_2.remaining_material`, `held_exact_motion_coverage.py:57-63`).
- Reflex-corner phase partition ("corner approaches", `held_figure5_engagement_refinement.py:175-193`).

### Where the plan is wrong or incomplete

- **Plan "Where we are", bullet 1, and Task 0b.** "Needs the publisher-derived start circle and raw guide" is contradicted (finding 2). The real dependencies are the start marker and three tolerance bounds (see (b)).
- **Task 0b, Step 4.** The plan scores Held on `figure5`, but cannot score Held on the gate pockets because it assumes no Held generator runs there. The native construction already runs on rect_12x8, rect_20x12 and L_shape (finding 3). Task 0b should score a native Held lap on the **gate pockets**, not only on Figure 5.
- **Task 0b, Step 5.** The decision rule ("fails fewer criteria than rho_regulated") is biased by `DEGENERATE_LOOP_RATIO = 1.0` (finding 4) and by tangency-noise gouges (finding 10). Add a ruling item: is a loop with rho <= r inside Held's engagement model a defect?
- **Plan, line 31.** It quotes 79.99 deg without the 1,010 containment rejections and without the whole-history bound result (findings 9 and 10).
- **Design, "Reposition above the stock".** The Held lap needs no repositioning at all: connectors are cutting moves along the wall offset. The plan's retract → LINK → plunge design is a regression compared with an existing capability.
- **Stage B, D2.** It ignores that every offset-contact construction marks convex offset corners as `is_stationary` (rho = 0). The exact corner witness is available natively, so B2a's "wall run-in, retract at apex" can be built from `cycle.ccw_transition` into the stationary contact, with no new geometry.
- **Stage C, "Held reproduction rejoins".** It treats Held as a comparison column. It should be a candidate generator topology (see (d)).

### (a) Generation pipeline, input to connectors

1. Input: `HeldReferenceCase` → `pocket_spec()`, the 65-vertex polygon projection with r = 1 and cap 80 deg, both hard-coded (`held_reference_cases.py:391-394, 464-475`).
2. Guide: the product's straight-skeleton `_guide_chains`. Figure 5 has 126 runs and 7,075 stations (measured) (`held_figure5_raw_guide.py:417-438`).
3. Circle hypotheses: for each station m, project to every polygon side exactly (Fraction). Keep sides within the publisher evidence gap. Then rho = (d - r)/2, q = foot + r·n, c = q + rho·n (`held_figure5_raw_guide.py:333-376`). Each circle runs from the wall offset to the guide station.
4. Ownership and canonicalization: keep the minimum-distance hypothesis per cyclic side neighbourhood, plus the hypothesis nearest the published start (`held_figure5_path.py:136-196`).
5. Contact projection onto the r-offset (`offset_polygon`, `held_reference_figures.py:493-510`): either CGAL `project_boundary_contact` within the evidence gap, or "qualified contacts" that reject out-of-bound hypotheses and keep them as evidence (`held_figure5_path.py:515-532`, `held_reference_contacts.py`).
6. Placement: per lane, bisection over ordered candidates to the 80 deg standard-model cap (`held_standard_placement.py:247-310`).
7. Ordering: sort by offset boundary site and phase at the start marker. That gives one CCW lap (`held_figure5_boundary_path.py:279-301`).
8. Spacing repair:
   - Recursive boundary-arclength midpoint insertion until every emitted pair is <= cap (`held_figure5_engagement_refinement.py:196-237`).
   - Alternative: boundary-order selection (`held_boundary_order_spacing.py`).
   - Alternative: whole-history contour-bound selection (`held_contour_bound_spacing.py`).
9. Corners: convex corners get a radius limit so the circle stays tangent inside both sides (rho → 0). Reflex corners (`--corner-approaches`) get two circles at the vertex, one per incident-side normal.
10. Connectors: a polyline along the r-offset from q_i to q_{i+1}. Inside a fan the move has zero length (`held_figure5_boundary_path.py:209-248`).
11. Entry: none.

**Against the world-class bar:**

| Criterion | Status |
| --- | --- |
| Single entry | partial. One start, but a virgin 360 deg loop. |
| Stays down | yes |
| No air loops | yes by construction, but fans repeat circles at one q: 16 zero-length fan transitions on Figure 5 (measured) |
| Capped engagement | predecessor-model only; connectors unaudited; whole-history bound unproven for 12 |
| G1 | tangent at q between loop and connector. Kinks: 62 connectors cross offset vertices, plus fans with differing normals (measured) |
| Coverage | unproven; the exact replay is set aside |
| Short path | no, 2.8x Held |

Overall: adequate topology, weak execution.

### (b) Publisher-derived dependencies, and what an arbitrary polygon needs

| Dependency | Where | Generalisation |
| --- | --- | --- |
| Start marker centre and radius | `held_figure5_path.py:508-511`, `:169-177`. `held_figure5_boundary_path.py:268-272` checks the marker radius equals r. `:290-301` seeds the phase. | Generator-chosen start, e.g. the contact of the entry helix. Drop the radius check. Note: `held_figure5_reference_guided.py:212-223` shows the marker is the publisher stream's **terminal** point. |
| Evidence tolerances: `admissible_distance_gap = 2·(reconstruction + projection)` | `held_figure5_raw_guide.py:84-94` | Used as the hypothesis admissibility gap (`:297-305`), the vertex-site resolution bound (`held_figure5_boundary_path.py:144-171`), the start bound, and the contact-projection bound. On an exact polygon the right value is 0, but `:268` requires > 0 and float contacts are not on the offset. Replace steps 3-5 with native `boundary_circle_at_contact`, whose contacts lie on the offset by construction. Then no tolerance is needed. |
| Raw guide | not publisher-derived | `_guide_chains(spec.polygon, ...)` already works on any `PocketSpec`. Only the `HeldReferenceCase` typing (`held_figure5_raw_guide.py:154-166, 419`) blocks it. |
| Reference vectors (publisher turns, Figure 6 digitisation) | `held_figure5_reference_guided.py`, `held_figure6_*`, `tools/held_figure*` | Comparison only. Generation never reads them. |
| Single offset component, no holes | `held_figure5_boundary_path.py:264-265`, `held_reference_figures.py:254-255`, `pocket_spec holes=()` | Fine for rect_12x8, rect_20x12 and L_shape (measured: one offset cycle each). A neck narrower than 2r splits the lap: one lap per component, with a link between them. |
| r = 1 and cap 80 deg hard-wired in the case type | `held_reference_cases.py:391-394` | Take r and cap from `PocketSpec`. |

### (c) Corner approaches and qualified contacts

- **Qualified contacts** filter ingestion noise. CGAL projects each hypothesis contact onto the offset. Those beyond the evidence bound are rejected but kept, and losing a whole guide run raises (`held_reference_contacts.py:21-49`). They have nothing to do with corners.
- **Corner approaches** handle **reflex** offset vertices only (`_turn < 0` on a CCW ring, `held_figure5_engagement_refinement.py:134, 188`). Two circles are emitted at the vertex, one with each incident side's normal, so interpolated phases never tilt through the boundary (`:100-101, 175-193`).
- **Convex corners** are handled by `_corner_radius_limit` (`:51-74`): rho <= d/(1 - n·m), which goes to 0 at the vertex.
- **Neither solves Stage B D2.** Shrinking toward a zero-radius loop at the offset vertex is the measured "bore it" variant, and it produces degenerate loops (495 Figure 5 circles below 0.1 r, measured). The useful piece for D2 is the native `is_stationary` event at each convex offset vertex: an exact apex for B2a's wall run-in.

### (d) Generalise Held rather than the tree walk?

**Pros**
- Single stay-down lap, with no retract, plunge or link. It removes the plan's root causes 1 and 2 by construction, except the first loop.
- Climb CCW. Every loop kisses the wall.
- A native exact construction exists for arbitrary line/arc boundaries, including reflex arcs (measured on the gate pockets).
- The whole-history bound and the exact deletion predicate already exist.
- It reproduces Held's Figures 5 and 8 at 79.99 deg in Held's model (baseline).

**Cons**
- Small loops: rho <= r for about half the circles. The gate's degenerate-loop rule needs a ruling.
- Convex corners stay D2.
- Entry stays D1.
- Kinks at offset vertices.
- The current Python is 2.8x Held's length and 21 s / 450 s per pocket (baseline `generation_seconds`, verified).
- Coverage is unproven: the contour-bound variant loses coverage.
- No islands.
- Everything is typed to `HeldReferenceCase` inside `benchmarks/`.

**Verdict:** the plan should not commit to the tree walk before Task 0b runs a native Held lap on the gate pockets with a ruling on rho <= r. Held's own published path, with rho >= 1.167 and 265 turns, shows that a short boundary lap is achievable. The tree walk keeps the rho > r floor but keeps N chains and their repositioning.

### (e) Representation and adapter

- **Circle.** `Figure5CounterclockwiseCircle` has `center`, `contact_point` q and `radius: GuideRadius`, always CCW. `radius` is the **tool-centre (machining-circle) radius rho**, not the cleared disk. The swept outer radius is rho + r (`held_standard_placement.py:181, 193-194`). Containment checks rho against the r-offset (`tools/held_toolpath_progress.py:146-148`).
- **Loop start and end** is q (`render_path` `:110-115`).
- **Connector.** `Figure5BoundaryTransition.samples` is a polyline on the r-offset, with float `boundary_progress` and a `within_fan` flag.
- **Adapter needs:**
  - PLUNGE at q0. D1 decides whether this becomes a helix.
  - Per circle, a CUT `Circle(Frame(c, xaxis = (q - c)/rho, yaxis = +90 deg), rho)` with `clockwise=False`.
  - CUT `Line`s for each connector segment. Skip `within_fan` zero-length moves.
  - Drop the wrap-around transition. `build_figure5_boundary_path` emits N transitions, including the closing one (`:303`). The benchmark length uses N - 1 (`held_toolpath_benchmark.py:173-175`). The plan's 26,436.6 mm = full circumferences + N - 1 connectors (verified, `:149-152`).
  - One RETRACT.
  - Drop native `is_stationary` (rho = 0) events, or emit them as point visits.
  - Expect the gouge criterion to see tangency noise (finding 10).
- Note: `HeldFigure5Characterization` already reduces the 12 gate criteria on Figure 5, but for the **product** generator (`benchmarks.runner.generate_toolpath`, `tests/tools/test_held_path_characterization.py:102`), not the Held path.

## Records

### benchmarks/held_boundary_domain_audit.py | lines 141 | read 1-141
- **role:** validation (reporting audit of the source curves)
- **used by:** `__main__` only (docstring `:3`). No pixi task found.
- **state:** live
- **generation:** n/a. It reports convex arcs with radius <= r, which are corners the cutter cannot follow.
- **defects:** none found. The 0.001 deg cutoff is display-only (`:26-27`).
- **unused capability / ideas:** the list of convex joins is a corner inventory the generator could consume.
- **claims:** "reporting only" (`:4`) verified. Nothing feeds a generator.

### benchmarks/held_boundary_order_spacing.py | lines 90 | read 1-90
- **role:** spacing
- **used by:** `held_toolpath_benchmark.py:25`, `tools/held_reference_toolpaths.py:10`, `tools/held_contour_bound_toolpaths.py:9`
- **state:** live (Figure 8 workloads)
- **generation:** selects circles in emitted (boundary) order with a single predecessor, using a linear lookahead and keeping each source run's last representative; then refines. Quality: adequate. It fixes finding 5 (4,183 → 1,684 circles on Figure 8 upper, baseline).
- **defects:**
  - `:53-62`: a run-retention rule forces circles in that engagement does not need. Source-run provenance constrains geometry.
  - `:56`: the lookahead tests `circles[index + 1]` against `selected[-1]` only (predecessor model).
- **unused capability / ideas:** linear lookahead rather than bisection, as the "engagement is not monotone" memory recommends.
- **claims:** "not a continuous swept-coverage certificate" (`:3`) verified (no stock).

### benchmarks/held_consumer_adapters.py | lines 36 | read 1-36
- **role:** validation adapter
- **used by:** `tools/held_path_characterization.py` (Figure 5 characterization of the product generator)
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** this is the hook Task 0b would reuse to score a Held-path `ToolpathResult`.
- **claims:** none

### benchmarks/held_contour_bound_spacing.py | lines 110 | read 1-110
- **role:** spacing / engagement-control
- **used by:** `held_coverage_path.py:5`, `tools/held_contour_bound_toolpaths.py:10`
- **state:** live (experimental)
- **generation:** contour-aware selection. The whole-history native bound with the exact `cap.chord_ratio` surrogate decides each cap. Over-cap gaps are subdivided by midpoint. Quality: weak. The draft is half the length (13,087 mm recorded) but loses coverage (finding 6).
- **defects:**
  - `:57-58, 72`: the first circle is assumed cleared, so there is no entry.
  - `:83`: the skip decision looks ahead one circle and ignores coverage, which is the cause of the coverage loss.
  - Connectors are not in the contour, so their engagement is unaudited (`:3-5`).
- **unused capability / ideas:** `HeldDiskContour2` is the exact contour-aware engagement primitive Held Section 3 calls for. The product generators do not use it.
- **claims:** "CGAL's whole-orbit leading-rim bound owns every cap decision" (`:3`) verified (`:77-78, 91`).

### benchmarks/held_coverage_path.py | lines 60 | read 1-60
- **role:** spacing (coverage-preserving pipeline)
- **used by:** `tools/held_contour_bound_toolpaths.py:11` (behind `--coverage-preserving`)
- **state:** parked. The CLI calls it a "rejected ... experiment" (`tools/held_contour_bound_toolpaths.py:37`).
- **generation:** dense contour-bound refinement of all sources, then exact thinning. Quality: weak (+82% length on Figure 5, recorded).
- **defects:** none in the code. The approach itself was rejected.
- **unused capability / ideas:** the ordering "refine, then delete under an exact coverage predicate" is sound.
- **claims:** "Source completeness and containment remain separate qualifications" (`:34`) verified.

### benchmarks/held_coverage_thinning.py | lines 114 | read 1-114
- **role:** spacing / validation hybrid
- **used by:** `held_coverage_path.py:6`
- **state:** parked (see above)
- **generation:** deletes a circle only when `Stock2.can_remove_circle` proves the replacement keeps the required stock. Connector chains are kept verbatim. Quality: adequate as a deletion rule.
- **defects:** `:72-73`: the first circle is subtracted with no entry modelled.
- **unused capability / ideas:** `Stock2.can_remove_circle` is an exact "does this loop remove anything needed" primitive. This is exactly the plan's "productive-only emission", already native.
- **claims:** "Delete circles only when native predicates preserve required swept stock" (`:54`) verified (`:93`).

### benchmarks/held_curved_circle_benchmark.py | lines 316 | read 1-316
- **role:** tool (performance gate)
- **used by:** pixi `held-curved-circle-benchmark` / `-baseline` (`pyproject.toml:208-209`), `held_curved_circle_performance_plot.py`
- **state:** live. Curved-circle performance is on the plan's parking lot.
- **generation:** n/a. It times `NativeBoundary2.circle_on_piece`, the exact Held circle on curved boundaries.
- **defects:** none found
- **unused capability / ideas:** proves the native curved query is under 0.25 s per piece (`:46`), so a native Held generator is fast enough per contact.
- **claims:** ceilings (`:36-47`) unverified (measured 2026-09-08 per comment).

### benchmarks/held_curved_circle_performance_plot.py | lines 88 | read 1-88
- **role:** figure
- **used by:** `__main__`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_exact_motion_coverage.py | lines 132 | read 1-132
- **role:** validation
- **used by:** `held_toolpath_benchmark.py:26`, the tools' `--check-coverage` / `--exact-coverage` flags
- **state:** parked (exact replay set aside 2026-09-08)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `_coverage_2.remaining_material` is the exact residual over circles and capsules.
- **claims:** "No entry clearing assumed" (`:100-101`) verified.

### benchmarks/held_figure5_boundary_path.py | lines 304 | read 1-304
- **role:** ordering / linking
- **used by:** all Held path modules
- **state:** live
- **generation:** projects circles onto the r-offset sides, sorts them by boundary site, phases at the published start, and emits CCW offset connectors. Quality: adequate. This is the single-lap topology.
- **defects:**
  - `:268-272`: requires the published start radius to equal r. Publisher-derived.
  - `:290-301`: the start seed comes from the publisher marker.
  - `:303`: emits N transitions, including the wrap back to the start, which is not a motion.
  - `:215-222`: fans are zero-length motions.
  - `:229-239`: connectors kink at offset vertices (62 on Figure 5, measured).
- **unused capability / ideas:** the native `cycle.ccw_transition` (`tools/held_exact_boundary_progress.py:51`) is the exact version of `_transition`, with arcs at reflex corners.
- **claims:** "one bounded, continuous CCW traversal" (`:263`) verified for position; G1 is not continuous at kinks.

### benchmarks/held_figure5_engagement_refinement.py | lines 237 | read 1-237
- **role:** corner / spacing
- **used by:** `held_boundary_order_spacing`, `held_contour_bound_spacing`, `held_toolpath_benchmark`, tools
- **state:** live
- **generation:** convex-corner radius limit, reflex-corner fan bridges, recursive midpoint insertion. Quality: weak.
- **defects:**
  - `:51-74`: drives rho → 0 at convex corners, giving degenerate loops (finding 8).
  - `:172`: inserted circles take the **average** radius rather than the clearance-derived (d - r)/2, so they need not reach the medial axis. Coverage is unproven.
  - `:153-159`: if the midpoint lands exactly on a vertex, it moves to 1/3. An ad-hoc nudge.
  - `:225`: decisions use the float predecessor model, not an exact surrogate.
- **unused capability / ideas:** none beyond the corner fan (`:175-193`).
- **claims:** "Containment must be assessed independently" (`:6-7`) verified.

### benchmarks/held_figure5_path.py | lines 548 | read 1-548
- **role:** placement
- **used by:** `held_toolpath_benchmark`, `held_boundary_order_spacing`, tools
- **state:** live
- **generation:** hypothesis canonicalisation, lane formation, Section 2.4 placement per lane, emission in boundary order. Quality: weak, because of the placement-versus-emission order mismatch (finding 5).
- **defects:**
  - `:169-177`, `:508-511`: dependency on the publisher start.
  - `:393-436` versus `:444`: the order mismatch.
  - `:412-416`: `preserve_source_runs` forces stops by provenance, not geometry.
  - `:520-532`: contact projection needs the publisher evidence gap.
- **unused capability / ideas:** none
- **claims:** "Uses the case's publisher-derived start; does not claim contour awareness" (`:502`) verified.

### benchmarks/held_figure5_publisher.py | lines 504 | read 1-504
- **role:** validation (reference data decoder)
- **used by:** `held_figure5_reference_guided.py`, `tools/held_figure5_reference_guided.py`
- **state:** live
- **generation:** n/a. It holds Held's real path: 265 turns, rho >= 1.167, about 9.4 m (measured), the efficiency target.
- **defects:** none found
- **unused capability / ideas:** use it as the quality oracle: length, turn count, rho distribution.
- **claims:** 265 turns, 1,218 cubics, 211 lines (`:33-35`) verified by load.

### benchmarks/held_figure5_raw_guide.py | lines 438 | read 1-438
- **role:** guide / placement (hypotheses)
- **used by:** `held_figure5_path`, `held_toolpath_benchmark`, tools
- **state:** live
- **generation:** builds the guide from the product skeleton, then the circle hypotheses: rho = (d - r)/2, contact on the r-offset, centre halfway between q and the station. Quality: adequate formula, but on sampled skeleton stations.
- **defects:**
  - `:154-166, 412-419`: restricted to the canonical Held cases by type.
  - `:84-94`: the publisher tolerance budget.
  - `:353-364`: float construction of q and c, so contacts are off the offset and need tolerance-based projection later.
- **unused capability / ideas:** the native `boundary_circle_at_contact` computes the same circle exactly, from the exact medial contact, without stations.
- **claims:** module title "pre-engagement-thinning guide runs" verified. Publisher-derived: **contradicted** (`:421-430`).

### benchmarks/held_figure5_reference_guided.py | lines 260 | read 1-260
- **role:** validation (correspondence to the publisher path)
- **used by:** `tools/held_figure5_reference_guided.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** finds the publisher marker is the **terminal** endpoint (`:221-223`), while generation uses it as the start.
- **claims:** "no independent path ... claim" (`:36-38`) verified.

### benchmarks/held_figure6_comparison.py | lines 67 | read 1-67
- **role:** validation / figure
- **used by:** `tools/held_figure6_same_axes.py`
- **state:** live
- **generation:** n/a. It measures the **product** `controlled_path` (`benchmarks/figure6.py:166`), not Held.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "compliance unaudited" verified.

### benchmarks/held_figure6_publisher.py | lines 174 | read 1-174
- **role:** validation (digitised reference)
- **used by:** `tools/held_figure6_same_axes.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_motion_coverage.py | lines 101 | read 1-101
- **role:** validation (sampled)
- **used by:** `tools/held_contour_bound_toolpaths.py:21`, `held_coverage_thinning.py:16`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "not a continuous proof" (`:19, 97`) verified.

### benchmarks/held_native_boundary_draft.py | lines 265 | read 1-265
- **role:** placement (diagnostic)
- **used by:** pixi `held-native-boundary-draft` (`pyproject.toml:249`)
- **state:** live, diagnostic only
- **generation:** exact native Held circles at dyadic offset contacts, plus exact CCW transitions. Quality: n/a (not emitted), but it is the right core.
- **defects:** none found
- **unused capability / ideas:** **the generalisation path for (b).** Measured working on rect_12x8, rect_20x12 and L_shape.
- **claims:** "no machining qualification" (`:1`) verified.

### benchmarks/held_native_curve_import.py | lines 39 | read 1-39
- **role:** kernel adapter (curved boundary import)
- **used by:** `held_curved_circle_benchmark`, the native plots
- **state:** live
- **generation:** n/a. It enables Held on true line/arc boundaries.
- **defects:** none found
- **unused capability / ideas:** curved pockets without polygonisation.
- **claims:** "no curve is flattened" (`:22-24`) unverified (native side).

### benchmarks/held_native_curve_import_plot.py | lines 121 | read 1-121
- **role:** figure
- **used by:** `__main__`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_native_curved_circle_plot.py | lines 176 | read 1-176
- **role:** figure / placement diagnostic
- **used by:** `__main__`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `circle_on_piece` on curved boundaries, as above.
- **claims:** none

### benchmarks/held_path_characterize.py | lines 111 | read 1-111
- **role:** validation (12-criterion harness)
- **used by:** `tools/held_path_characterization.py`, tests
- **state:** live
- **generation:** n/a. Takes any `Generator(spec) -> ToolpathResult`.
- **defects:** none found
- **unused capability / ideas:** Task 0b can call `characterize_figure5` with a Held adapter instead of writing a new harness.
- **claims:** none

### benchmarks/held_path_evidence.py | lines 468 | read 1-468
- **role:** validation
- **used by:** characterize, report, post-qualification
- **state:** live
- **generation:** n/a
- **defects:** `:208-210`: the characterisation is restricted to Figure 5, so it cannot score gate pockets.
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_path_report.py | lines 271 | read 1-271
- **role:** validation (report)
- **used by:** `tools/held_path_characterization.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_path_snapshot.py | lines 486 | read 1-486
- **role:** validation (snapshot)
- **used by:** evidence and characterize
- **state:** live
- **generation:** n/a
- **defects:** `:130-132`: radius must be > 0, so native stationary rho = 0 events cannot be represented.
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_post_qualification.py | lines 90 | read 1-90
- **role:** validation (refusal gate)
- **used by:** report
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_reference_cases.py | lines 898 | read 1-898
- **role:** validation / input (case loader)
- **used by:** all Held modules
- **state:** live
- **generation:** the input. r = 1 and cap 80 are hard-coded (`:391-394`). `pocket_spec` gives tool diameter 2.0 and no holes (`:464-475`). The start marker comes from the publisher (`:398-400, 842-848`).
- **defects:** none for its purpose. Hard-wiring blocks reuse.
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_reference_certification.py | lines 1141 | read 1-1141
- **role:** validation (ingestion certification of cubic → arc)
- **used by:** `held_reference_geometry.py`
- **state:** live
- **generation:** n/a
- **defects:** none relevant to generation
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_reference_contacts.py | lines 49 | read 1-49
- **role:** placement (contact qualification)
- **used by:** `held_figure5_path.py:520`
- **state:** live (Monstera workload)
- **generation:** rejects contacts farther than the evidence bound from the offset. Quality: adequate as a filter, but only needed because contacts are float constructions.
- **defects:** depends on the publisher evidence bound.
- **unused capability / ideas:** none
- **claims:** "A rejected hypothesis is preserved verbatim" (`:30`) verified.

### benchmarks/held_reference_figures.py | lines 915 | read 1-915
- **role:** figure, plus the inward offset used by generation
- **used by:** `figure7_inward_offset` in all Held generators
- **state:** live
- **generation:** `figure7_inward_offset` (`:493-510`) is the r-offset via `offset_polygon` (straight skeleton), one component required.
- **defects:** `:254-255`: one component only.
- **unused capability / ideas:** the native `build_center_boundary_cycle` is the exact equivalent, with reflex arcs.
- **claims:** none

### benchmarks/held_reference_geometry.py | lines 1129 | read 1-1129
- **role:** validation (PDF cubic → certified arcs, projection)
- **used by:** cases, figures
- **state:** live
- **generation:** n/a
- **defects:** none relevant
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/held_standard_placement.py | lines 310 | read 1-310
- **role:** engagement-control / placement
- **used by:** all Held spacing modules
- **state:** live
- **generation:** Held Equation 4 overlap bound, Section 2.3 engagement, and Section 2.4 bisection to the cap within 0.001 rad. Quality: adequate (paper-faithful).
- **defects:**
  - `:261-274`: bisection assumes monotone engagement over candidates, which the memory says is false in general.
  - `:268`: the cap decision uses float acos, not the exact chord surrogate.
  - `:276-286`: if nothing is admissible, a "forced" placement of the first candidate can exceed the cap.
- **unused capability / ideas:** `maximum_predecessor_engagement` is a closed-form engagement model the product does not use.
- **claims:** "neither owns nor consults a global depleted-stock contour" (`:5-6`) verified.

### benchmarks/held_stock_replay.py | lines 35 | read 1-35
- **role:** validation
- **used by:** `held_motion_coverage`, `tools/held_stock_plot`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "without inventing a plunge" (`:29`) verified.

### benchmarks/held_toolpath_benchmark.py | lines 402 | read 1-402
- **role:** tool (generation regression gate)
- **used by:** pixi `held-toolpath-benchmark` / `-baseline` (`pyproject.toml:213-214`)
- **state:** live
- **generation:** wires the five workloads end to end.
- **defects:** `:162`: `containment_rejections` counts tangency rounding (finding 10), so it is a noisy metric.
- **unused capability / ideas:** none
- **claims:** baseline numbers in the plan (2,013 circles, 26.4 m, 79.99 deg; 1,684; 4,367) verified against `benchmarks/results/held_toolpath_baseline.json`.

### benchmarks/tools/__init__.py | lines 1 | read 1-1
- **role:** build
- **used by:** package
- **state:** live
- **generation:** n/a
- **defects:** none
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_circle_distribution.py | lines 75 | read 1-75
- **role:** figure
- **used by:** `tools/held_reference_toolpaths.py:20`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "Cap compliance does not imply efficient spacing" (`:68`) verified by finding 1.

### benchmarks/tools/held_contour_bound_plot.py | lines 79 | read 1-79
- **role:** figure
- **used by:** `tools/held_contour_bound_toolpaths.py`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_contour_bound_toolpaths.py | lines 110 | read 1-110
- **role:** tool (CLI generator)
- **used by:** `__main__`
- **state:** live
- **generation:** contour-bound and coverage-preserving drafts.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "rejected dense-source preservation experiment" (`:37`) verified as the recorded status.

### benchmarks/tools/held_contour_contact_plot.py | lines 124 | read 1-124
- **role:** figure
- **used by:** `held_toolpath_progress`, `held_reference_toolpaths`
- **state:** live
- **generation:** n/a. It exercises `HeldDiskContour2.contact_toward`, the Section 3.1 exposed-contact correction.
- **defects:** none found
- **unused capability / ideas:** `contact_toward` is the contour-aware contact correction, not used in placement.
- **claims:** "Placement and engagement are not yet contour-aware" (`:103`) verified.

### benchmarks/tools/held_coverage_preserving_plot.py | lines 73 | read 1-73
- **role:** figure
- **used by:** `tools/held_contour_bound_toolpaths.py`
- **state:** parked
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_exact_boundary_progress.py | lines 111 | read 1-111
- **role:** figure (native connector)
- **used by:** `__main__`
- **state:** live
- **generation:** n/a. It shows `cycle.ccw_transition`, the exact offset connectors with reflex arcs.
- **defects:** none found
- **unused capability / ideas:** the exact connector for a native Held lap.
- **claims:** "Still pending: analytic segment/arc MAT traversal, machining circles" (`:102`) is out of date: the native circle construction exists (`held_native_boundary_draft.py`).

### benchmarks/tools/held_figure5_reference_guided.py | lines 121 | read 1-121
- **role:** figure
- **used by:** `__main__`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_figure6_same_axes.py | lines 149 | read 1-149
- **role:** figure
- **used by:** pixi `_held-toolpath-figure6`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_reference_toolpaths.py | lines 105 | read 1-105
- **role:** tool (CLI generator)
- **used by:** pixi `_held-toolpath-*` (`pyproject.toml:216-219`)
- **state:** live
- **generation:** standard-model drafts with lane or boundary spacing, optional qualified contacts and corner approaches.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** the scope text "Publisher-start polygon guide" (`:95`) is accurate on the start; "guide" is the product skeleton.

### benchmarks/tools/held_stock_plot.py | lines 70 | read 1-70
- **role:** figure
- **used by:** `held_toolpath_progress`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### benchmarks/tools/held_toolpath_progress.py | lines 211 | read 1-211
- **role:** figure / tool
- **used by:** `held_toolpath_benchmark:48` (`_containment_report`), tools
- **state:** live
- **generation:** n/a. `render_path` is the de facto contract for the representation: the loop starts at q, runs CCW, and is full.
- **defects:** `:139-155`: containment counts tangency noise (finding 10).
- **unused capability / ideas:** none
- **claims:** "Predecessor-circle model only; entry, connectors, depleted-stock engagement and coverage require separate checks" (`:90`) verified.
