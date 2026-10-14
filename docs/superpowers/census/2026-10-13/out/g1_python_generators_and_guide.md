## Slice summary

Slice `g1_python_generators_and_guide`: the four Python generator lineages, the C++ straight-skeleton guide under them, and the kernels and stubs they stand on. 28 of 28 files read in full.

Two probes were run with `pixi run python`, and the scripts were deleted afterwards:
- **P1.** `polygon_skeleton_clearance`, `_straight_skeleton_2.create_interior_straight_skeleton` and `_guide_chains` on the gate L, `[[0,0],[12,0],[12,6],[6,6],[6,10],[0,10]]`, and on rect_12x8, with D = 2 and guide step = 0.05.
- **P2.** `_stock_2.engagement_at` on rect_12x8 after `subtract_disk(6,4,1)`, at ρ ∈ {0.5, 1.0, 1.048, 1.5}.

### Top findings, ranked by impact on toolpath quality

1. **The engagement model is static and blind to direction, and both the generators and the gate use it.** (verified)
   - `engagement_at(stock, cx, cy, r, cap, gap)` takes a position only, with no travel direction (`src/engagement_2.cpp:521`). It reports the largest contiguous engaged arc over the whole cutter rim.
   - The gate measures every arc station against the stock as it was *before* that arc (`engagement.py:396-397`). The generators probe against the pre-loop stock (`engagement_toolpath.py:550-555`).
   - So a loop's own sweep behind the tool still counts as material.
   - The spiral module's "240° floor" is exact *for this model*. P2 measured 240.00° at ρ = r and 243.20° at ρ = 1.048 r, which matches the closed form to the digit.
   - It is not a physical floor. Analysis, not measured: for the same ρ = r loop around a centre plunge, the leading half of the cutter meets material over [−90°, 30°], which is 120°.
   - The codebase already has direction-aware measures:
     - `HeldDiskContour2.engagement_bound` ("leading-rim", `_stock_2.pyi:19-24`);
     - `audit_swept_prefix_segment_tea_exact` (`_continuous_tea_2.pyi:603-611`).
   - The plan's Stage B D1 (a new HELIX op type, "a planar circle cannot avoid it") rests on the static model. The model should be settled first.
2. **The guide is consumed narrow end first, one plunge per chain, in greedy nearest-neighbour order starting from the longest chain.**
   - `toolpath.cpp:656-661` orients each chain narrow end first. `toolpath.cpp:963-993` orders chains by nearest neighbour from the longest, and never reverses a chain.
   - `engagement_toolpath.py:438-453` reads that stream with `optimize_order=True`. Its plunge, retract and link are at `:784`, `:812` and `:932`.
   - P1 on rect_12x8: 5 chains. Chain 0 is the spine at ρ 2.998, and four corner chains start at ρ 0.016. (verified)
3. **Bridges in the Python generators are certified for neither wall clearance nor engagement.**
   - C++ certifies every bridge and lead against the boundary (`certify_bridges` at `toolpath.cpp:427-481`, `certified_lead` at `:885-903`).
   - The Python generators instead emit their own entry-to-entry lines (`engagement_toolpath.py:802`, `engagement_rho_toolpath.py:679`) that skip guide stations.
   - These lines have no `segment_clear` check and no engagement check. The radial module says so at `engagement_radial_toolpath.py:59-60`.
   - Gouge risk near reflex corners: unverified.
4. **The advance-only lineage, which the gate tests, still uses bisection, which assumes monotone engagement.**
   - `_largest_admissible_advance` bisects (`engagement_toolpath.py:601-654`). `radius_regulated_toolpath` calls the same bisection (`engagement_radial_toolpath.py:138`, `:852`).
   - Its own module docstring says engagement is not monotone and bisection must not be used (`:31-51`).
   - Only rho and spiral scan, from the furthest candidate down (`engagement_rho_toolpath.py:480-537`).
5. **Entry placement causes the 24 rect_12x8 tangent breaks; rho's external-tangent entry fixes them only between adjacent constant-slope stations.**
   - The entry is placed along the guide normal (`engagement_toolpath.py:353-354`), so the bridge tilts by atan(dρ/ds). The rho docstring derives 35.26° (`engagement_rho_toolpath.py:47-61`). The arithmetic checks: atan(sin 45°) = 35.26°.
   - Rho's slope is a central difference over neighbouring *guide* stations (`:298-305`), but the chosen advance skips up to 40 stations.
   - Hypothesis, unverified: at clearance kinks the bridge is not the common tangent. P1 found a kink on L chain 5. This would explain rho's 2 L_shape tangent breaks, the plan's B4.
6. **The generators deplete their internal stock with different primitives from the gate replay.**
   - Generators: `subtract_capsule_quad` (under-covers, keeping everything outside `radius*(1-1e-4)` per `stock.py:166`) and `subtract_arc_sweep_local`.
   - Gate replay: `_subtract_operation`, i.e. `subtract_capsule` and `subtract_arc_sweep` (`engagement.py:168-223`).
   - So probes see slivers the replay never has, and the two judge "redundant" differently. The plan's `_Emitter` fixes this by using `_subtract_operation`.
7. **The straight skeleton is adequate for the gate's coverage and gouge-safety, but loses quality at the L's reflex hub.** See (b).
   - P1: all three interior SS vertices of the gate L have clearance 3.000. The widest guide station (ρ 2.492, clearance ≈ 3.494) lies on the reflex-bisector chain, not at any vertex.
   - The exact medial-axis vertex is at (3.515, 3.515) with clearance 6(2−√2) = 3.515. (analysis)
8. **Convex corners below the ρ > r floor have no motion in the current vocabulary** (rho docstring `:73-88`).
   - Exact building blocks for the plan's B2a wall run-in exist, unused: `ReachableBoundaryCycle2.ccw_transition`, `CutterCentreDomain2`, `ReachableDomain2.reachable_material/unreachable_residual` (`_coverage_2.pyi:56-136`), and the SS inward offsets (`straight_skeleton_2.cpp:146-187`).

### Generation capabilities in this slice that the gated generators do not use

- **C++ guide (`toolpath.cpp`):**
  - wall-certified external-tangent bridges (`external_tangent` :317-360, `certify_bridges` :427-481);
  - certified leads (`certified_lead` :885-903);
  - flat-link gouge check (:1059-1069).
  - The gated generators set `lead_in=0`, `link_paths=False` and rebuild bridges in Python.
- **Rho lineage (`engagement_rho_toolpath.py`)**, not gated:
  - degeneracy floor ρ > r;
  - radius-step bound;
  - external-tangent G1 entries;
  - furthest-first scan;
  - `DeclinedRegion` as returned data.
- **Spiral entry (`engagement_spiral_entry_toolpath.py`):** a centre plunge plus a concentric ramp, scored on the gate in the docstring only.
- **Radial lineage:** `_loop_reaches_material`, an exact point-membership "does this loop bite" test (`engagement_radial_toolpath.py:507-544`). It is far cheaper than the plan's clone-and-compare.
- **`Stock` / `_stock_2`:**
  - `subtract_annulus` (exact full loop);
  - the `*_local` variants;
  - `is_subset_of` / `exactly_equals` / `clone`;
  - `can_remove_circle`, a proof that a circle can be deleted without losing material, i.e. air-loop thinning (`_stock_2.pyi:342-357`);
  - `HeldDiskContour2.engagement_bound` / `contact_toward`, leading-rim engagement and contact correction (`:14-30`);
  - `subtract_circle_sweep`.
- **`_continuous_tea_2`:**
  - continuous exact engagement certificates for full circles and segments (`audit_full_circle_tea_event_exact` :537-546, `audit_segment_tea_event_exact` :594-602);
  - swept-prefix segment audit (:603-611).
- **`_coverage_2`:**
  - `ReachableDomain2` with exact reachable material and unreachable residual. This is the right denominator for "uncut" and for D2(a).
  - `Coverage2` residual tracking;
  - `remaining_material`;
  - `CutterCentreDomain2`;
  - boundary-cycle transitions (`ccw_transition`), for wall-following run-ins;
  - `NativeBoundary2.circle_on_piece`.
- **`_containment_2`:** exact containment for a segment, a full circle, an entry disk, and a circle within an entry disk (`:57-96`). This is the gouge certificate the Python bridges lack, and support for a helical or pre-drilled entry disk.
- **`_circle_geometry_2` / `boundary_normal_circle_2`:**
  - contour-aware circle placement in the style of Held & Pfeiffer: `BoundaryNormalCircle2.query` and `query_vertex` for reflex point-site cones;
  - `boundary_circle_engagement`;
  - `disk_contains_disk`, `swept_disk_intersection`, `corrected_engagement_cosine_squared`.
- **`_medial_axis_2`:** the exact segment-site medial axis with node and edge topology, per-sample clearance and guide radius, and neck evidence and classification. It is limited to fixtures, see (b).
- **`straight_skeleton_2`:** inward, outward and weighted offset polygons, for contour-parallel and finishing passes and the cutter-centre boundary.
- **`motion.py` / `units.py`:** a typed motion vocabulary (`ExactSegmentMotion`, `ExactCircleMotion`, `EngagementCap`, `EntryRadius`, `CutPlane`). It is used by `adaptive/` and the Held benchmarks, never by the gated generators.

### Where this slice shows the draft plan is wrong or incomplete

- **P-1. The plan's L is not the gate's L.**
  - Plan Tasks 2 and 3 use `[[0,0],[16,0],[16,6],[6,6],[6,14],[0,14]]`, area 144.
  - The gate L is arm 6, long arm 12, short arm 10, i.e. `[[0,0],[12,0],[12,6],[6,6],[6,10],[0,10]]`, area 96 (`benchmarks/gate.py:84-89`, `:121-124`).
  - Task 3's stop rule ("virgin==1 fails on l_shape … a fact about the gate's L") therefore tests the wrong pocket. (verified)
- **P-2. The root is chosen among SS vertices, not over the guide.**
  - On the gate L, ids 12 (9,3), 13 (3,3) and 15 (3,7) all have clearance 3.000 (P1). The lowest-id tie-break puts the root at the far end of an arm.
  - The widest station, ρ 2.492, is on chain 5, the reflex bisector between vertices 3 and 13.
  - "Exact clearance" is also a misnomer: under Epick `K::FT` is `double` (`toolpath.cpp:12`), so the tie is decided on doubles.
- **P-3. "Away from the root" does not mean "wide to narrow".**
  - On the L, the DFS walk 13→3 runs chain 5 from ρ 1.998 up to ρ 2.492. On the dumbbell, the second room is entered from its neck end.
  - Only the rectangle test checks that radius falls along each walk (plan Task 2 test).
- **P-4. The DFS order produces one retract, link and plunge per leaf.**
  - Ascending-id neighbour order visits contour corners (ids 0..n−1) before interior vertices (P1 ids: corners 0–5, interior 12, 13, 15). The first walk from the root therefore goes to a corner, and every leaf ends a walk.
  - rect: 5 plunges and 5 retracts. L: 8 chains, so up to 8.
  - Virgin plunges fall to 1, but the retract count does not move toward the commercial bar of one entry with the tool kept down. D3 (stay-down) is the real fix and is deferred to Stage C.
- **P-5. "Productive" is defined as `not exactly_equals`.**
  - A 360° loop that removes a sliver of any size counts as productive, so air-cutting loops survive under the gate's definition.
  - Cost: `_removes_material` clones the whole exact arrangement for every candidate. The `_first_productive` fallback scans to the chain end. Neither is measured.
  - Cheaper exact filter: radial `_loop_reaches_material`. Exact pruning afterwards: `Stock2.can_remove_circle`.
- **P-6. The `_next_step` fallback breaks rho's rules.**
  - The fallback is `_first_productive(origin+1 .. last)`. It ignores `_radius_step_admissible` and the advance window, so it can jump more than MAX_ADVANCE = 1 D.
  - The plan says "station rules … reused unchanged". Rho's fallback respects the radius step (`engagement_rho_toolpath.py:533-537`).
- **P-7. Bridges stay uncertified against the walls,** as in rho. The C++ `certify_bridges` and the exact `_containment_2` are not used. The plan doesn't mention it.
- **P-8. D1's premise is model-relative.** See finding 1. The 240° floor holds only for the static, direction-blind `engagement_at`. Before ruling HELIX into the replay vocabulary, rule on the engagement model: leading-rim, or swept-prefix as already implemented for segments.
- **P-9. Smaller defects:**
  - Task 1 reports truncation with a bare `std::invalid_argument`, i.e. `ValueError`. The engineering bar asks for a named exception.
  - Task 3 stamps `DeclinedRegion.path_index` with the walk index even for walks that emit nothing. This breaks the join contract documented at `engagement_rho_toolpath.py:395-400`.
  - Task 3 uses global `_subtract_operation` for every cut and every probe clone. The existing generators use the `*_local` variants (`engagement_toolpath.py:795`). The performance impact is not measured.
- **P-10. Corners: B2a can be built exactly from what exists.**
  - The wall run-in follows the cutter-centre boundary: `build_center_boundary_cycle` and `ReachableBoundaryCycle2.ccw_transition` (`_coverage_2.pyi:56-66`), or the SS inner offset at r (`straight_skeleton_2.cpp:146-160`).
  - D2(a) "rest material" has an exact definition ready: `ReachableDomain2.reachable_material()` / `unreachable_residual()`.
- **P-11. "Use the exact MAT instead" is not a drop-in.** `_medial_axis_2` builds only an exact rectangle and one hard-coded L fixture, `(0,0),(6,0),(6,2),(2,2),(2,6),(0,6)` (`src/segment_site_catalog_graph.cpp:125-193`, `:905`, `:1207`). It raises `UnsupportedCanonicalMatLShapeGraphError` on the gate L.
- **Checked and found sound:**
  - `_enter`'s fallback plunge is clear by construction at both call sites. In-walk, the bridge capsule contains the target entry disk; at walk start, the preceding station's loop removes nothing.
  - All Task 1–3 imports and call signatures match the real code: `_guide_chains`, `_Regulation.build`, `_loop_operation(station, cut_z, path_index)`, `_station_is_admissible`, `_declined_regions`, `_first_trochoidal_station`, `_line_operation`, `_tessellate`, `_holes_to_vertex_arrays`, `_polygon_to_ccw_vertices`, `RADIAL_CLEARANCE_FRACTION`, `GUIDE_STEPOVER_ADVANCE_MULTIPLE`. Lines 469-483 of `_guide_chains` are exactly the station loop.

### (a) Generator lineages compared

| lineage | guide | placement | ordering / linking | entry | corners | continuity | engagement control | short of world class |
|---|---|---|---|---|---|---|---|---|
| C++ `trochoidal_mat_toolpath_circular` (`toolpath.cpp:908-1177`) | SS bisectors, narrow end first (:656-661) | stations at min(pitch, stepover); refined by the proxy `spacing+growth ≤ stepover` (:245-270) | greedy NN from the longest chain, never reversed (:963-993); retract, link at clearance, plunge (:1046-1057) | plunge at the chain start (narrow end, ρ≈0.016 on rect) | bores the corners (stations down to ρ→0) | external-tangent G1 bridges, certified against the walls (:427-481); full circle plus reposition arc where radius varies (:396-397) | none, geometric proxy | no engagement control; full-circle + reposition arc overcuts; per-chain plunges |
| `engagement_controlled_toolpath` (gated) | the C++ circles read back (`engagement_toolpath.py:438-467`) | **bisection** of the advance on a D/40 grid (:601-654), 33 sampled probes (:141, :492-522) | C++ NN order; plunge/retract/link per chain (:784, :812, :932) | full-immersion loop at the narrow end; the over-cap verdict is computed and then emitted anyway (:787) | bores them (degenerate loops) | entry on the guide normal, so 24 breaks on rect (rho docstring :53-56) | exact predicate at sampled positions; bridges unregulated | monotonicity assumption; corner entries; air loops in the cleared centre; no G1 |
| `radius_regulated_toolpath` (gated) | same | same bisection advance (:852) plus a scanned radius ladder (:691-754) | repeated sweeps per chain, each with its own plunge/retract (:825, :869, :966) | same as above | same | same | ladder; reported-double decisions (:749, :619-621), disclosed | at 120° the vector equals engagement_controlled (plan table); more entries at tighter caps |
| `rho_regulated_toolpath` (not gated) | same | furthest-first scan (:480-537), radius-step bound (:465-477), floor ρ > r (:351-365) | same per-chain plunges (:643, :684, :817) | full immersion at the first station with ρ > r | declines them, returning `DeclinedRegion` (uncut 1–3.5%) | external-tangent entry (:235-275), G1 between adjacent constant-slope stations | same sampled predicate | per-chain plunges; leaves the corners; slope estimate assumes adjacent stations |
| `spiral_entry_toolpath` (proof artefact) | same | rho walk after the ramp | same | centre plunge plus concentric ramp chosen by the predicate (:205-305) | as rho | radial links, two breaks per turn (:272-277) | same | fails 6 criteria per its docstring (:8-11); its 240° floor is a property of the model |

### (b) Is the straight skeleton the right guide?

**On convex pockets it is identical to the medial axis** (`toolpath.cpp:1185-1186`, `toolpath.py:6-8`).
- Every station radius uses the true Euclidean clearance (`toolpath.cpp:217-221`), so it is gouge-safe on any polygon.
- The rectangle failures are therefore not caused by the guide. They come from orientation, ordering and entry.

**On the gate L (P1, verified),** the SS has interior vertices (9,3), (3,3) and (3,7), all at clearance 3, plus a reflex bisector (6,6)→(3,3). The medial axis instead has:
- a vertex at (3.515, 3.515) with clearance 3.515;
- parabolic arcs y = ((x−6)²+36)/12 from it to (6,3), and the mirror arc to (3,6);
- straight spines only from x = 6 and y = 6 outward. (analysis)

**What the difference costs:**
1. Loop radius at the hub is ρ ≤ 1.998 on the SS spine between (3,3) and (6,3), against up to 2.513 on the medial axis.
2. The reflex chain has clearance slope 1 toward (6,6). The C++ dedup collapses the contained disks, leaving chain 5 with 16 stations and a kink near (3.53, 3.53). P1: chain 5 first station ρ 2.492, last ρ 1.998.
3. The guide has a 90° kink at (3,3) between the spines, where the medial axis is G1. This matters only once walks are chained without retracting, which the plan does not do.
4. Root and entry selection over SS vertices misses the widest point.

Coverage is not lost: the hub loops reach about 2.998 from (3,3), and chain 5 reaches toward the reflex corner.

**Verdict:** the SS is adequate for Stage A on the gate. It is the wrong long-term guide for contour-aware placement. The exact medial axis exists only for fixtures (P-11). The general exact alternative already in tree is the boundary-normal medial query (`BoundaryNormalCircle2.query/query_vertex`), which is Held's contour-aware placement and is not fixture-limited, but holes are excluded (`boundary_normal_circle_2.h:75-77`).

### (c) What the spiral module proved, and what the radial ladder found

**Spiral:**
- With a centre plunge (void V = r), a concentric loop at ρ has static `max_run_tea` = 360° − 2·acos(ρ/2r). This increases with ρ, so for ρ > r the minimum is 240°; as ρ → 0 it tends to 180° (`engagement_spiral_entry_toolpath.py:41-80`).
- P2 matches the closed form exactly at all four radii. Verified for the static model.
- The engagement step after the entry is therefore at least 240° − cap. At a 120° cap that infimum is exactly 120°, approached only from above (`:61-68`). The 243.20° → 124.64° ramp numbers are unverified.
- Physically, by leading-half analysis, the ρ = r first loop engages 120°. So the "floor" says what the gate's measure can and cannot score, not what entering stock must cost.

**Radial:**
- Engagement is not monotone in loop radius. At one 20x12 station it falls to 69.2° at R = 4.148, then rises to 72.4°, 77.0° and 85.0° as R shrinks. Bisection would return rung 36 (3.198) where the scan returns rung 17 (4.148) (`engagement_radial_toolpath.py:31-45`). Unverified; quoted from the docstring, which says a test pins the state.
- The gentlest refused rung is interior: 300.9° → 253.7° → 360° (`:569-576`). Unverified.
- Contradiction: the same module still bisects the advance (`:852` → `engagement_toolpath.py:636-654`).

### (d) Unused capabilities

See the summary list above. All were verified by an import grep over `src`, `benchmarks` and `tests`:
- `_continuous_tea_2`, `_containment_2` and `_medial_axis_2` are imported only by `src/compas_cgal/adaptive/*` and tests.
- `_circle_geometry_2` only by the Held benchmarks and tests.
- `straight_skeleton_2` only by `benchmarks/held_reference_figures.py` and tests.
- `can_remove_circle` and `HeldDiskContour2` only by the Held benchmarks and tests.
- `engagement_rho_toolpath` only by spiral and tests.
- `engagement_spiral_entry_toolpath` only by its own test.

### (e) The plan's Task 1–3 code against real signatures

See P-1 to P-11 above.

---

### src/boundary_normal_circle_2.cpp | lines 254 | read 1-254
- **role:** placement (contour-aware circle), kernel.
- **used by:** `_circle_geometry_2` (`circle_geometry_2_bindings.cpp:146-151`) → `benchmarks/held_standard_placement.py`, `held_figure5_path.py`, `held_reference_contacts.py`, `held_native_boundary_draft.py`, and tests `test_boundary_normal_circle.py` and others. Not used by the gated generators.
- **state:** live (Held lane).
- **generation:**
  - For a boundary point p on a segment interior, or a reflex vertex with a direction in its normal cone, it finds exactly the first medial contact m along the inward normal (`:150-196`).
  - It then builds the machining circle with centre midpoint(q, m) and guide radius (clearance − r)/2, whose outer cutter disk is tangent to the wall at p and nested in the medial disk (`:223-235`).
  - Excellent as a placement primitive: exact, general simple polygons, and it handles reflex point sites. It is Held & Pfeiffer's contour-aware placement.
- **defects:**
  - `:219-221` raises `NoPositiveBoundaryCircleError` when clearance ≤ r, so convex-corner tips have no circle. This is the corner gap again, now as a hard error.
  - Holes are not supported (header `:75-77`).
  - Epeck_with_sqrt `CGAL::sqrt` at `:69`, `:96`, `:179` is the CORE path. Memory notes that CORE stalls on identically-zero decisions; the risk here is unmeasured.
- **unused capability / ideas:** a contour-aware placement source for a product generator: boundary-spaced circles near walls and medial circles in the interior. `at_contact` (`:99-148`) gives exact inverse queries for transitions.
- **claims:**
  - "complete outer cutter disk is nested in the boundary-tangent medial disk": verified, checked exactly at `:228-234`.
  - "no earlier competing contact is possible" (`:135-139`): unverified, argued in a comment.

### src/boundary_normal_circle_2.h | lines 94 | read 1-94
- **role:** placement, kernel (declarations).
- **used by:** `boundary_normal_circle_2.cpp`, `circle_geometry_2_bindings.cpp`, `NativeBoundary2` (friend, `:72`).
- **state:** live.
- **generation:** declares the proposal (p, m, q, centre, guide radius, clearance, competing features) and the query API. Adequate.
- **defects:** none found. The Epeck_with_sqrt kernel choice is noted under the .cpp.
- **unused capability / ideas:** `exact_contact` / `exact_center` / `exact_guide_radius` (`:48-50`) give native-side consumers exact geometry, for example to chain circles without double round-trips.
- **claims:** "Only reporting properties cross into doubles" (`:44`): verified, since the `*_mm()` accessors are the only double outputs.

### src/circle_geometry_2_bindings.cpp | lines 164 | read 1-164
- **role:** kernel (bindings), placement / engagement helpers.
- **used by:** `_circle_geometry_2` importers listed in (d), all in the Held lane or tests.
- **state:** live.
- **generation:**
  - Exposes `BoundaryNormalCircle2`, `boundary_circle_engagement` (the standard predecessor engagement bound) and `disk_contains_disk`, exact with tangency (`:49-57`).
  - Also `swept_disk_intersection` (`:72-89`) and `corrected_engagement_cosine_squared` (`:90-115`); both decide existence exactly and report in doubles.
  - Adequate.
- **defects:**
  - `:86-88` reporting uses `std::sqrt(to_double)`. That is reporting only, which is fine.
  - The local exception classes at `:24-31` duplicate the pattern used elsewhere. Effect: none.
- **unused capability / ideas:** `disk_contains_disk` is an exact test for whether a loop is contained in an earlier one, a cheap redundant-loop filter. The plan uses arrangement cloning instead.
- **claims:** "Report squared chord cosine after exact corrected-intersection checks": verified at `:99-113`.

### src/compas_cgal/__init__.py | lines 11 | read 1-11
- **role:** build.
- **used by:** package import; the compas plugin loader reads `__all_plugins__`.
- **state:** live.
- **generation:** n/a.
- **defects:** none found. Only `compas_cgal.toolpath` is registered as a plugin (`:9-11`).
- **unused capability / ideas:** none.
- **claims:** "Only modules that still exist may appear": verified, `toolpath.py` exists.

### src/compas_cgal/_circle_geometry_2.pyi | lines 82 | read 1-82
- **role:** kernel (stub).
- **used by:** type checkers; mirrors `circle_geometry_2_bindings.cpp`.
- **state:** live.
- **generation:** n/a. It is the API surface of contour-aware placement (`:39-54`) and of the predecessor engagement bound (`:76-82`).
- **defects:** `boundary_circle_engagement` is declared as "conditional on its cleared outer disk" (`:81`). That condition is a model assumption, and callers must establish it.
- **unused capability / ideas:** see the .cpp.
- **claims:** "numerical properties are reporting only" (`:16`): verified against `boundary_normal_circle_2.cpp:245-252`.

### src/compas_cgal/_containment_2.pyi | lines 97 | read 1-97
- **role:** kernel (stub): containment and gouge certification.
- **used by:** `src/compas_cgal/adaptive/containment.py`, `adaptive/errors.py`, `tests/adaptive/test_containment.py`, `tests/test_exact_region_traits_ownership.py`.
- **state:** live (adaptive lane only).
- **generation:** exact containment records for a segment, a full circle, an entry disk, and a full circle inside an entry disk (`:57-96`).
- **defects:** none found.
- **unused capability / ideas:**
  - It is exactly what the Python generators' uncertified bridges need, as a gouge certificate.
  - `evaluate_exact_entry_disk_containment` and `evaluate_exact_full_circle_in_disk` support an entry that opens a disk (helical or pre-drilled) and then runs circles inside it. That is relevant to D1.
- **claims:** none beyond the signatures.

### src/compas_cgal/_continuous_tea_2.pyi | lines 631 | read 1-631
- **role:** kernel (stub): continuous exact engagement certification.
- **used by:** `src/compas_cgal/adaptive/{identity,motion_certificate,motion_oracle_cache,motion_refutation}.py` and `tests/adaptive/*`.
- **state:** live (adaptive lane).
- **generation:**
  - Continuous exact verdicts (certified, cap exceeded or unresolved) per full circle (`:537-546`) and per segment (`:594-602`).
  - A swept-prefix segment audit (`:603-611`) that accounts for material the motion itself has already removed.
  - Exact probes at rational parameters (`:547-558`, `:612-622`).
- **defects:** n/a for generation.
- **unused capability / ideas:**
  - These would replace the 33-probe sampling in every Python generator (`engagement_toolpath.py:141`) with continuous verdicts.
  - The swept-prefix concept is the fix for finding 1's direction-blind measure, but it exists for segments only, with no circle variant.
  - Cost is not measured here. Memory: the TEA2 segment audit cost is CGAL's algebraic kernel.
- **claims:** none.

### src/compas_cgal/_coverage_2.pyi | lines 279 | read 1-279
- **role:** kernel (stub): reachable domain and coverage.
- **used by:**
  - `benchmarks/coverage.py`, `survey.py`, `toolpath_coverage.py`, `qualityfigures.py`;
  - the Held benchmarks and `src/compas_cgal/adaptive/{coverage,reachable_domain,errors}.py`;
  - tests.
- **state:** live.
- **generation:**
  - `ReachableDomain2`: design region, centre domain, reachable material, unreachable residual (`:122-136`).
  - `CutterCentreDomain2.contains` (`:138-146`).
  - `Coverage2`: exact residual after disk, segment and full-circle sweeps (`:192-234`).
  - `remaining_material` (`:83-96`).
  - Boundary cycle with `ccw_transition` along the cutter-centre boundary (`:56-66`).
  - `NativeBoundary2` for curved boundaries (`:270-279`).
- **defects:** none found.
- **unused capability / ideas:**
  - The exact basis for (i) D2(a) "rest material" (reachable vs unreachable), (ii) the plan's B2a wall run-in (walk the cutter-centre boundary cycle), and (iii) a coverage-driven stop rule instead of "walk every station".
  - None of it is used by any generator in this slice.
- **claims:** "Exact final residual, without accumulated sweeps or chronological history" (`:90`): unverified here.

### src/compas_cgal/_medial_axis_2.pyi | lines 107 | read 1-107
- **role:** kernel (stub): exact segment-site medial axis.
- **used by:** `src/compas_cgal/adaptive/medial_axis.py`, `tests/adaptive/test_medial_axis.py`, `test_neck.py`, `typecheck/consumer_contract.py`.
- **state:** live (adaptive lane), fixture-limited.
- **generation:** graph nodes and edges, per-sample centre, clearance and guide radius (`:9-30`, indices 9-14), neck evidence and classification (`:88-93`), zero-guide records.
- **defects:** the binding surfaces `UnsupportedCanonicalMatLShapeGraphError` (`:39`). The builder only accepts the canonical rectangle and the one 6x6 L fixture (`src/segment_site_catalog_graph.cpp:159-193`).
- **unused capability / ideas:** the topology, clearance and neck data are exactly what the plan's tree walk needs, plus neck detection for the dumbbell. Not usable on the gate L until generalised.
- **claims:** none.

### src/compas_cgal/_stock_2.pyi | lines 575 | read 1-575
- **role:** kernel (stub): exact stock, depletion and engagement.
- **used by:** every generator in this slice through `stock.py` and direct `_stock_2.engagement_at` (`engagement_toolpath.py:552`, `:596`), plus `engagement.py`, `motion.py` and the benchmarks.
- **state:** live.
- **generation:**
  - `Stock2` booleans (`:331-444`).
  - `engagement_at(stock, cx, cy, r, cap, gap)`, position only (`:551-558`).
  - `certify_segment_tea` (`:559-567`).
  - `cap_chord_ratio` (`:549`).
  - Held contour helpers (`:14-30`).
  - `can_remove_circle` (`:342-357`).
  - Exact depletion traces (`:385-412`).
- **defects:** `engagement_at` takes no direction, so engagement on the trailing side counts. That inflates loop and slot engagement relative to cutting physics (finding 1, analysis).
- **unused capability / ideas:**
  - `HeldDiskContour2.engagement_bound` is the leading-rim model.
  - `can_remove_circle` provides exact air-loop removal.
  - Also `subtract_annulus`, `subtract_circle_sweep`, `subtract_exact_full_circle` and the exact induction checks (`:525-546`).
- **claims:** "Report whole-orbit leading-rim angle bound … Does not certify full-rim entry, connectors, containment, or actual stock" (`:19-24`): unverified here.

### src/compas_cgal/engagement.py | lines 472 | read 1-472
- **role:** validation (audit replay), and the depletion mapping reused by generation.
- **used by:** `benchmarks/{depletion,exceedance,pathmetrics,runner,survey,held_*}.py`, tests, and the plan's Task 3 (`_subtract_operation`).
- **state:** live.
- **generation:** n/a, but it constrains generation:
  - `_subtract_operation` (`:168-223`) is the gate's depletion. Lines become `subtract_capsule`; Circles and Arcs become `subtract_arc_sweep`.
  - Arc engagement is fixed at 20 stations per turn, with a guard, against the stock before the operation (`:226-317`, `:396-397`).
- **defects:**
  - `:396-397`: measured before its own subtraction. For a full loop, the trailing portion the loop itself has already cut counts as material, so the reported loop engagement is too high. Magnitude unmeasured.
  - `:361`: a plunge removes a plain disk, so a helical entry is unrepresentable. This is a known constraint (the spiral docstring at `:22-34` references `benchmarks.depletion`).
- **unused capability / ideas:** none for generation.
- **claims:**
  - "never rewrites the toolpath" (`:10`): verified.
  - "Every interior center lies within half a station spacing" (`:235-236`): verified by construction, `:280`.

### src/compas_cgal/engagement_radial_toolpath.py | lines 1144 | read 1-1144
- **role:** engagement-control (radius ladder), placement.
- **used by:** `benchmarks/gate.py` (gated as `radius_regulated`) and `tests/test_engagement_radial_toolpath.py`.
- **state:** live (gated).
- **generation:**
  - A descending scanned radius ladder per station (`:691-754`).
  - Repeated sweeps per chain until every station has emitted its maximal circle (`:907-988`).
  - The advance comes from the advance-only module's **bisection** (`:852`).
  - Weak: one plunge and retract per sweep (`:825`, `:869`), each with its own LINK (`:966`); otherwise inherits every engagement_controlled defect.
- **defects:**
  - `:138`, `:852`: bisection advance, which contradicts the module's own doctrine (`:31-51`). Effect: it can choose a non-furthest admissible advance, giving a longer path.
  - `:749`, `:619-621`: decisions made from reported doubles (disclosed).
  - `:822-829`: reduced-radius loop entries along the guide normal (via `_GuideStation.entry`), so bridges between unequal radii are not tangent.
  - `:288`: `MAX_RADIAL_SWEEPS_PER_CHAIN` is an empirical budget.
- **unused capability / ideas:** `_loop_reaches_material` (`:507-544`) is an exact membership bite test; the plan should use it instead of clone-and-compare.
- **claims:**
  - Non-monotone radius numbers (`:33-45`): unverified.
  - "at a 40 deg cap it takes the forced count from 88 to 44" (`:18-19`): unverified.
  - "worst requirement … TWO passes" (`:283-284`): unverified.

### src/compas_cgal/engagement_rho_toolpath.py | lines 849 | read 1-849
- **role:** engagement-control, entry placement (G1), corner (decline).
- **used by:** `engagement_spiral_entry_toolpath.py`, `tests/test_engagement_rho_toolpath.py`, `tests/test_engagement_spiral_entry_toolpath.py`, and the plan's Task 3. Not in `benchmarks/gate.py`.
- **state:** live, not gated.
- **generation:**
  - Degeneracy floor ρ > r (`:351-365`), radius-step bound (`:465-477`), furthest-first scan (`:480-537`).
  - External-tangent entry for G1 bridges (`:235-275`, `:278-319`).
  - Declined regions returned as data (`:368-462`).
  - Adequate, and the best lineage here. It still plunges per chain (`:643`, `:684`, `:817`) and leaves the corners uncut.
- **defects:**
  - `:298-305`: the slope is a central difference over neighbouring guide stations, but bridges join stations up to `window` = 40 apart (`engagement_toolpath.py:295`). The bridge is the exact common tangent only where the clearance slope is constant along the span. At clearance kinks (P1: L chain 5) G1 is lost. Hypothesis: these are the L's 2 breaks.
  - `:264-269`: |slope| ≤ 1 is checked on a double. A guide chain with slope exactly 1 (the reflex bisector before dedup) sits on the boundary of the check.
  - `:679`: bridges are not certified against the walls.
  - `:680`: depletion uses `subtract_capsule_quad`, while the gate uses `subtract_capsule` (finding 6).
- **unused capability / ideas:** all of it is unused by the gate. `DeclinedRegion` should feed a corner pass (B2a).
- **claims:**
  - Variant table (`:80-84`): rho row consistent with the plan's measured table (0.0295, 0, 0, 0) at 120°; other rows unverified.
  - 24 breaks at 35.26° = atan(sin 45°) (`:53-56`): the arithmetic is verified; the count matches the plan's pinned 24.
  - "Entering virgin stock … full-immersion cut at EVERY non-degenerate radius" (`:65-71`): verified for the static model (P2 extends this to centre plunges).

### src/compas_cgal/engagement_spiral_entry_toolpath.py | lines 542 | read 1-542
- **role:** entry (proof artefact).
- **used by:** `tests/test_engagement_spiral_entry_toolpath.py` only.
- **state:** parked. It describes itself as "not the entry strategy" (`:3-11`).
- **generation:**
  - Centre plunge (`:244-245`), then concentric turns, each the largest rung the predicate admits, forced to the smallest rung otherwise (`:251-300`).
  - After the ramp, a rho walk (`:375-392`).
  - Weak as a generator: two tangent breaks per ramp turn (`:272-277`).
- **defects:**
  - `:245`: the plunge is subtracted as a zero-length `subtract_capsule_quad`, while the replay uses `subtract_disk`. The generator stock differs from the replay stock.
  - `:268`: on a forced turn it takes the minimum reachable rung, without ranking.
  - `:329-330`: function-local imports.
- **unused capability / ideas:**
  - The law-of-cosines bound (`:82-93`, V ≥ sqrt(ρ²+r²−2ρr·cos((360−θ)/2))) says how large a hole an entry must open for a given cap. It is the design equation for a helical bore's diameter under the static model.
- **claims:**
  - 240° floor (`:41-59`): verified by P2 for the static `engagement_at`. As a physical floor it is contradicted by leading-half analysis (120° at ρ = r).
  - Ramp 243.20 → … → 104.10 and largest step 124.64 (`:65-68`): 243.20 verified by P2; the rest unverified.
  - "fails six criteria where rho fails four" (`:9-10`): unverified.

### src/compas_cgal/engagement_toolpath.py | lines 958 | read 1-958
- **role:** engagement-control, placement, ordering, linking, entry (gated).
- **used by:** `benchmarks/gate.py`, `survey.py`, `figure6.py`, `held_figure5_raw_guide.py`, all three sibling generators, tests, and the plan's Tasks 1–3.
- **state:** live (gated, frozen since 2026-08-22 per the brief).
- **generation:**
  - Guide read back from the C++ circular stream (`:400-489`).
  - Advance by bisection (`:601-654`), with 32 ring probes plus the entry probe (`:141`, `:492-522`).
  - Forced minimum advance (`:648-653`).
  - Plunge, retract and clearance link per chain (`:784`, `:812`, `:932`).
  - Weak against the bar: per-chain virgin plunges at corner tips, loops that only air-cut in the cleared centre, tangent breaks, bores at corners.
- **defects:**
  - `:614-617`, `:636-654`: bisection under a monotonicity assumption that the sibling modules and project memory say is false.
  - `:353-354`: entry on the guide normal, causing tangent breaks wherever the radius varies.
  - `:438-453`: `optimize_order=True` with direction fixed narrow end first by C++.
  - `:787`: the entry verdict is computed and the full-immersion loop is emitted anyway.
  - `:795`, `:808`: local and quad depletion differ from the gate replay (finding 6).
  - `:802`: bridges are not checked against walls or engagement.
  - `:932`: a LINK is emitted even between chains whose entry and exit coincide.
- **unused capability / ideas:**
  - The `_Regulation.build` factory (`:212-318`) is good and is reused by the plan.
  - The exact rational cap surrogate is crossed once (`:357-381`).
- **claims:**
  - "43.4 deg on a circle every probe accepted" (`:845-846`): unverified.
  - "no residual material survives further than 0.29 mm from a wall" (`:103-106`): unverified.
  - "Reusing it keeps … the chain ordering in the one C++ implementation" (`:36-41`): verified, and it is the root cause of the corner entries.

### src/compas_cgal/errors.py | lines 23 | read 1-23
- **role:** build (error model for units and motion).
- **used by:** `units.py`, `motion.py`, `adaptive/generator.py`, `adaptive/operation.py`, `benchmarks/held_reference_cases.py`, tests.
- **state:** live.
- **generation:** n/a.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** "Errors raised only by the certification lane live in `compas_cgal.adaptive.errors`" (`:5-6`): unverified.

### src/compas_cgal/motion.py | lines 87 | read 1-87
- **role:** motion-geometry (typed motion vocabulary).
- **used by:** `src/compas_cgal/adaptive/*`, the `benchmarks/held_*` modules, `benchmarks/tools/held_*`, acceptance tests. Not used by any generator in this slice.
- **state:** live (adaptive and Held lanes).
- **generation:** `ExactSegmentMotion`, `ExactCircleMotion` (with `squared_radius` as a Fraction, `:57-61`) and `EngagementCap` (native chord ratio plus bytes, `:64-87`). Adequate as a vocabulary.
- **defects:**
  - Typing: `typing.Self` (`:5`) needs Python 3.11 or later. The project stack says Python 3.9+, but pixi pins 3.12. Effect: none under pixi.
  - The `build()` factories add nothing over `__post_init__` (`:29-31`, `:48-55`, `:85-87`).
- **unused capability / ideas:** the gated generators pass raw floats. A product generator should emit these typed motions so the exact audit lanes can consume them.
- **claims:** none.

### src/compas_cgal/replay_classification.py | lines 140 | read 1-140
- **role:** validation (shared replay classification).
- **used by:** `engagement.py`, `benchmarks/{depletion,quality_observations,survey,toolpath_coverage}.py`, tests.
- **state:** live.
- **generation:** n/a, but it constrains generation:
  - `AUDIT_ENGAGED` includes LINK (`:18-25`).
  - A LINK at cut height classifies as "motion" (`:100-104`), so stay-down links are replayed as cuts and become redundant operations when nothing is removed.
  - Any line with both Z change and XY travel raises (`:96-98`), so ramps and helices are unrepresentable.
- **defects:** `:96-98` and `:18-25` block the two world-class entry and linking moves (ramp/helix entry, stay-down link). This needs a gate ruling (plan D1/D3).
- **unused capability / ideas:** none.
- **claims:** none.

### src/compas_cgal/stock.py | lines 365 | read 1-365
- **role:** kernel (typed stock wrapper).
- **used by:** all four generators, `engagement.py`, `benchmarks/{coverage,depletion,exceedance,instrument,qualityfigures,survey}.py`, `adaptive/{bootstrap,replay,stock_area}.py`, tests.
- **state:** live.
- **generation:** exact material removal and membership. Adequate.
- **defects:**
  - `:21-68` duplicates `toolpath._polygon_to_ccw_vertices` (`toolpath.py:107-134`), by design per `:24-26`.
  - `:152-183`: `subtract_capsule_quad` under-covers by design (contract at `:162-167`). The generators use it while the gate uses `subtract_capsule` (finding 6).
  - `can_remove_circle` and `subtract_circle_sweep` exist on `Stock2` but not on this wrapper.
- **unused capability / ideas:** `subtract_annulus` (exact full loop) and `is_subset_of`. Wrap `can_remove_circle` to thin air loops.
- **claims:**
  - "All queries and removals are exact (no tolerance)" (`:76-78`): contradicted in part. The quad capsule and the partial-arc disk chains under-approximate by contract (`:155-167`, `:230-232`). They are exact booleans of approximating shapes.

### src/compas_cgal/straight_skeleton_2.py | lines 264 | read 1-264
- **role:** guide (upstream API: skeleton graph, offsets).
- **used by:** `benchmarks/held_reference_figures.py`, `tests/test_straight_skeleton_2_*.py`. No generator uses it.
- **state:** live (upstream compas_cgal API).
- **generation:** the SS as a compas Graph with edge types, plus inward, outward and weighted offset polygons. Adequate as an API; unused for generation.
- **defects:**
  - `:74-75`, `:108-117`: orientation is checked with a tolerance on the normal, and the inner offset raises a bare `ValueError`. That is upstream style.
  - `:199`: `points.points.reverse()` on what may be a list (`holes` elements) would raise `AttributeError` for a raw list input. Effect: none on generation.
- **unused capability / ideas:** the offset at r gives the cutter-centre boundary for wall run-ins (B2a), a finishing contour, and contour-parallel clearing of narrow regions.
- **claims:** none that matter.

### src/compas_cgal/toolpath.py | lines 520 | read 1-520
- **role:** guide (Python face of the C++ guide), linking, types.
- **used by:** every generator in this slice (types and `trochoidal_mat_toolpath_circular`), `engagement.py`, `stock.py`, `replay_classification.py`, tests, and the plan.
- **state:** live.
- **generation:** wraps the C++ generators. `ToolpathOperation` and `ToolpathResult` are the shared output types (`:87-104`). Adequate as glue.
- **defects:**
  - `:41-50` has `__all__` (against the user's Python rule). Effect: none.
  - `:143-149`: truncation by `max_passes` only warns. The gated generators inherit a silently partial guide (warning only), and the plan's Task 1 correctly turns this into a raise.
  - `:195-201`: `_row_to_circle` drops the clockwise flag. This is noted in memory. Effect: `Circle` geometry carries no winding; the flag lives on `ToolpathOperation.clockwise`.
- **unused capability / ideas:** `lead_in` / `lead_out` / certified flat links (`:385-400`) are unused by the Python generators.
- **claims:**
  - "every bridge, lead, and flat link is certified against the boundary at tool radius" (`:3-6`): verified for the C++ path (`toolpath.cpp:427-481`, `:885-903`, `:1063-1067`). It does not hold for the Python generators that consume it.
  - "SS coincides with the medial axis only for convex polygons" (`:6-8`): verified by P1 on the L.

### src/compas_cgal/types.py | lines 58 | read 1-58
- **role:** build (type aliases).
- **used by:** `src/compas_cgal/toolpath.py` (`PolylinesNumpy`), `straight_skeleton_2.py` (relative import, `:13-15`).
- **state:** live.
- **generation:** n/a.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/compas_cgal/units.py | lines 311 | read 1-311
- **role:** motion-geometry (typed units and frames).
- **used by:** `replay_classification.py`, `engagement.py`, `motion.py`, about 40 `benchmarks/held_*` modules, `benchmarks/{depletion,quality_observations,survey}.py`, `adaptive/*`.
- **state:** live.
- **generation:** n/a for the gated generators, which use raw floats. It provides NewType units, phantom frame tags (`WorldXY`) and validated value types (`ToolRadius`, `EntryRadius`, `Clearance`, `GuideRadius`, `CutPlane`).
- **defects:** `Self` needs Python 3.11 or later (`:8`). Pixi pins 3.12.
- **unused capability / ideas:** `EntryRadius` (`:196-205`) and `CutPlane` (`:298-311`) are the typed pieces a helical-entry or stay-down design would use.
- **claims:** none.

### src/medial_axis_2.cpp | lines 472 | read 1-472
- **role:** guide (exact segment-site medial axis bindings).
- **used by:** `src/compas_cgal/adaptive/medial_axis.py`, `tests/adaptive/test_medial_axis.py`, `test_neck.py`.
- **state:** live (adaptive lane). The builder is fixture-limited.
- **generation:** binds `SegmentSiteMatBundle2.build` with tool radius, station spacing, sagitta and refinement depth (`:320-348`). The projection gives nodes, edges, provenance, samples (centre, clearance, guide radius) and neck evidence (`:350-370`). Excellent in kind, unusable on general pockets today.
- **defects:**
  - `:343-344`: `UnsupportedCanonicalMatLShapeGraphError` is surfaced. The graph builder requires the canonical exact rectangle or the hard-coded L `(0,0),(6,0),(6,2),(2,2),(2,6),(0,6)` (`src/segment_site_catalog_graph.cpp:125-193`). It throws on the gate L (arm 6, 12 by 10). Effect: no exact medial-axis guide for the product gate's reflex pocket.
- **unused capability / ideas:** neck classification against squared-width boundaries (`:280-299`) is the exact answer to the plan's dumbbell Review Focus 2.
- **claims:** none in this file.

### src/straight_skeleton_2.cpp | lines 461 | read 1-461
- **role:** guide (upstream SS and offset bindings).
- **used by:** `straight_skeleton_2.py`.
- **state:** live (upstream).
- **generation:** the SS as vertices, edges and types (`:18-67`); inner, outer and weighted offsets (`:146-301`).
- **defects:**
  - `:46`: edges are canonicalised by pointer comparison `&*vertex1 < &*vertex2`, so the edge order is not reproducible across runs. `toolpath.cpp:640-642` deliberately uses ids for this reason. Effect: none on the generators, which don't use this file.
  - `:100-108`: no simplicity or orientation validation, unlike `toolpath.cpp:522-546`.
- **unused capability / ideas:** see `straight_skeleton_2.py`.
- **claims:** none.

### src/straight_skeleton_2.h | lines 144 | read 1-144
- **role:** guide (declarations).
- **used by:** `straight_skeleton_2.cpp`.
- **state:** live (upstream).
- **generation:** n/a.
- **defects:** `:75-80` holds commented-out dead declarations. Effect: none.
- **unused capability / ideas:** none beyond the .cpp.
- **claims:** none.

### src/toolpath.cpp | lines 1257 | read 1-1257
- **role:** guide, placement, ordering, linking (C++ generator plus the guide the Python generators read).
- **used by:** `toolpath.py` via `_toolpath`. Through `_guide_chains`, every Python generator in this slice. The plan's Task 1 modifies it.
- **state:** live.
- **generation:**
  - SS bisectors (`:634-688`), stations at min(pitch, stepover) with clearance from the boundary (`:663-671`), valid runs (`:224-239`), proxy engagement refinement (`:245-270`).
  - Dedup of contained disks (`:276-312`), external-tangent chains (`:362-421`), bridge wall certification (`:427-481`).
  - Greedy nearest-neighbour ordering (`:963-993`), certified leads (`:885-903`), retract, link and plunge (`:1040-1074`).
  - As a guide: adequate. As a generator: weak, since it has no engagement control.
- **defects:**
  - `:656-661`: narrow end first; with no reversal during ordering, every chain starts at its narrowest station. Corner chains start at ρ ≈ 0.016 (P1).
  - `:12`, `:151-169`: Epick, so `K::FT` is double. The "exact clearance" in comments and in the plan is an Epick double distance.
  - `:219`, `:230`, `:261`, `:454`: `clearance >= min_centerline` is decided on `sqrt(to_double(...))`, a double decision (CONVENTIONS.md rule). Effect: negligible.
  - `:334-352`: tangent construction decomposes into raw doubles, against CONVENTIONS.md's kernel-type rule. It is a construction, so the effect is nil.
  - `:396-397`: varying radius emits a full circle, then a reposition arc on the same circle, which cuts air.
  - `:636`: `is_inner_bisector()` is redundant with `is_bisector()`.
  - `:946-952`: chains without arcs are dropped, so `_guide_chains` sees only arc-bearing chains.
- **unused capability / ideas:** `certify_bridges`, `certified_lead` and the flat-link check are unused by the Python generators; the C++ external-tangent construction is re-done in Python by rho.
- **claims:**
  - "No interpolation of clearance anywhere — this is the structural fix for the reflex-vertex gouge" (`:186-188`): verified, `make_station` queries per station.
  - "Clearance is 1-Lipschitz so refinement … converges geometrically" (`:244`): the math holds; convergence within `MAX_BRIDGE_REFINE_ROUNDS` is unverified.

### src/toolpath.h | lines 107 | read 1-107
- **role:** guide (declarations).
- **used by:** `toolpath.cpp`. The plan's Task 1 adds `pmp_skeleton_guide` here.
- **state:** live.
- **generation:** n/a.
- **defects:** the docstring says "every station radius derives from an EXACT clearance query" (`:36-37`). Under Epick the squared distance is computed in doubles.
- **unused capability / ideas:** none.
- **claims:** "exact clearances (largest gouge-free disk)" (`:16-18`): contradicted in the strict sense. Epick constructs `squared_distance` in double (`toolpath.cpp:12`, `:151-158`), so these are rounded distances, not exact ones. The effect on gouge safety is sub-ulp.
