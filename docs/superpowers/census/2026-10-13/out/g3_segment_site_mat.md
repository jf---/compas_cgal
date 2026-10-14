# Census slice g3_segment_site_mat

Read by: Opus 5.5 (opus-5-5). 51/51 files, every line. Two probes run with
`pixi run python` (commands and output quoted under "Measurements").

## Slice summary

### Top findings, ranked by impact on toolpath quality

1. **The exact MAT runs on exactly one hard-coded pocket, and none of the gate pockets is that pocket.** The only path reachable from Python is `SegmentSiteMatBundle2::build` → `canonical_l_shape_mat_proposal_graph` (`segment_site_mat_bundle.cpp:39`). It accepts only the literal ring (0,0),(6,0),(6,2),(2,2),(2,6),(0,6) (`segment_site_catalog_graph.cpp:172-192`), with no holes (`:164`), and asserts the exact cell cardinalities of that one fixture (`:1577-1599`). Measured: `rect_12x8`, `rect_20x12`, the gate `L_shape` (0,0),(12,0),(12,6),(6,6),(6,10),(0,10), and a ×2 scaled copy of the fixture all raise `UnsupportedCanonicalMatLShapeGraphError`. **verified (probe).** The rectangle adaptor (`canonical_rectangle_mat_graph`, `:905-1205`) is reached only from native tests (`tests/native/task9_catalog_graph_gate.cpp:680`, `task9_graph_csr_gate.cpp:157`). It is not bound and also hard-asserts 5 S–S and 8 incident cells (`:1159-1165`). So today it is not usable as a guide for any product pocket.
2. **The cell machinery underneath is mostly general, but the traversal that assembles it is not.** Exact charts (parameterization), clipping to `D` and to `clearance ≥ r`, the live-vertex endpoint binders, clearance profiles, neck evidence and CSR export are written per *cell kind*, not per fixture. The L-shape adaptor loop (`segment_site_catalog_graph.cpp:1262-1575`) already dispatches all five kinds: P–P, incident P–S, non-incident P–S parabola, parallel S–S, nonparallel S–S. What is missing for arbitrary pockets is (a) dropping the fixture gate and the cardinality asserts, (b) four hard limits listed in (d) below, and (c) holes and cycles. The doc's own maturity line says so: "General arbitrary-pocket traversal, arbitrary composite segment/segment cells, rectangle neck evidence, arbitrary-pocket sampling/verdicts … are not complete" (`docs/segment_site_mat.md:201-204`). **verified.**
3. **What it has that the straight-skeleton guide lacks** (all exact, all unused by the gated generators):
   - true P–S **parabolic** bisectors at reflex vertices (`PARABOLA` edges, `segment_site_catalog_graph.cpp:1516-1574`). The skeleton replaces these with a straight reflex bisector;
   - an exact squared-clearance polynomial per edge (`segment_site_catalog_neck.cpp:72-145`);
   - the tool-centre-admissible subgraph `clearance > r` with exact clearance-root endpoints (`segment_site_clipping.cpp:1149-1169`);
   - **necks**: exact interior clearance minima, shared-vertex, endpoint and plateau necks, each certified by a separating graph cut (`segment_site_neck_evidence.cpp:363-554`, `segment_site_neck.cpp:221-330`), and exact classification of neck width against caller thresholds (`segment_site_neck_classification.cpp:73-116`);
   - a "zero-guide" inventory of edges whose clearance is exactly `r`, i.e. tool-width slots (`segment_site_zero_guide.cpp:88-118`).

   It does **not** have a rooted tree or a traversal order. The graph is undirected and unrooted. Order and policy live in Python `compas_cgal.adaptive.traversal` (outside this slice).
4. **On the two rectangle gate pockets the MAT and the straight skeleton are the same graph** (convex polygon: every bisector is an S–S line). Switching guides would change nothing there. All rectangle failures in the plan's table are traversal and emission defects, not guide defects. Only `L_shape` (one reflex vertex) gets a different, parabolic guide segment. **verified** for the shape of the rectangle graph (`segment_site_catalog_graph.cpp:1124-1157` rejects every P–S cell as an incident ray). The effect on the L_shape gate numbers is **unverified.**
5. **Sampling is limited to rational charts.** A nonparallel S–S bisector is rational only when `|n1|²·|n2|²` is a perfect square (`segment_site_parameterization.cpp:399-410`). Otherwise the sampling chart throws `UnsupportedCanonicalMatSamplingCurveError` (`segment_site_catalog_sampling.cpp:122-125`). That covers axis-aligned and Pythagorean-direction corners, not generic polygon angles. Clipping and clearance stay exact in Q(√α). Only the sampling and reporting layer requires rational coordinates. **verified by code reading.**
6. **Cost and size.** One bundle build on the L fixture takes 0.125 s (**measured**). Endpoint binding compares `CORE::Expr` roots for identity against the live CGAL vertex (`segment_site_endpoint_binding.cpp:795-807`, `:992-1005`). That is the identity-decision pattern the project memory says stalls CORE on generic doubles. It has never been exercised on generic input because no generic input reaches it: **unverified risk.** About half of the slice's lines are proof, identity and replay plumbing (certificate, evidence bytes, zero-guide records, CSR identity checks), not geometry.

### Generation capabilities in this slice the gated generators do NOT use

(`engagement_controlled_toolpath` and `radius_regulated_toolpath` never import `_medial_axis_2`. `engagement_toolpath.py:29` says it "deliberately does not use" `compas_cgal.adaptive`.)

- Exact medial axis with parabolic arcs around reflex vertices. The guide stays truly medial near reflex corners, where the skeleton's straight bisector hugs the corner. (`segment_site_endpoint_binding.cpp:1391-1527`, `segment_site_clipping.cpp:1591-1640`)
- Exact clearance function along every guide edge, not station-sampled: degree ≤ 4 rational polynomial (`segment_site_neck_clearance.cpp:430-446`, `segment_site_catalog_neck.cpp:147-179`). Loop radius could be a closed-form function of the guide parameter.
- Exact tool-centre-admissible clipping (`clearance ≥ r`) with exact cut-off points (`segment_site_clipping.cpp:513-544`, `:992-1255`). This replaces `toolpath.cpp`'s `min_centerline()` station filter.
- Neck detection and separating cuts. These answer "where does a second region begin and which edges are on each side" exactly (`segment_site_neck_evidence.cpp:419-539`). The draft plan's dumbbell case needs exactly this.
- Neck-width classification against exact rational thresholds, e.g. `(2ρ_min + 2r)²` (`segment_site_neck_classification.cpp:120-144`).
- Zero-guide (tool-width slot) detection (`segment_site_zero_guide.cpp:88-118`).
- Holes are carried through the catalog, clipping and domain tests (`segment_site_catalog.cpp:333-341`, `segment_site_clipping.cpp:202-213`, `segment_site_parameterization.cpp:613-632`), though no adaptor accepts them yet.

### Where this slice shows the draft plan (`2026-10-13-product-gate-tree-walk.md`) is wrong or incomplete

- **Design / Tree walk, "root is the interior skeleton vertex of largest exact clearance":** sound, and on the MAT it is provably enough. On every edge kind here, squared clearance is constant (parallel S–S), `k·t²` (nonparallel S–S, `segment_site_catalog_neck.cpp:97-110`) or the focal distance along a parabola (convex). So the global maximum lies at a node or on a constant plateau, never strictly inside an edge. The plan should state this lemma, because it is what makes "root = argmax over vertices" exact. On the straight skeleton the same holds only because skeleton edges are straight. Near a reflex vertex the skeleton's clearance is the distance to that vertex, which is still convex along a line.
- **Review Focus 2 and the Task 3 dumbbell ("clearance 1.5 there, rho 0.498 < r, so the far room needs its own entry"):** the plan reinvents neck detection as a station-level ρ threshold. This slice already has the exact version: a local clearance minimum, a separating cut naming which edges lie beyond the neck, and an exact width class (`segment_site_neck_evidence.cpp:419-431`, `segment_site_neck_classification.cpp:73-116`). The plan should at least specify the same semantics (minimum of clearance along the guide, side partition by graph cut) so the two converge. It should not count `virgin_plunges` after the fact.
- **Review Focus 1, "pocket with an island (skeleton has a cycle)":** the MAT side is not ready either. Both adaptors reject holes (`segment_site_catalog_graph.cpp:130`, `:164`). The plan's claim that the skeleton walk handles cycles is independent of this slice.
- **Missing alternative:** the plan names only the straight skeleton. It should say plainly why the exact MAT is not used for Stage A: it supports one fixture and none of the three gate pockets (finding 1, measured). For rectangles it would be the same graph anyway (finding 4). It should also record the concrete work that makes the MAT the guide later; see (e).
- **Station rules "reused unchanged from rho_regulated":** these derive loop radius from sampled station clearance. The MAT gives the clearance function in closed form. A tree walk built on the MAT could place stations by the exact clearance polynomial and need no `clearance dip down to spacing/256` resampling (`toolpath.cpp:55`). Not needed for Stage A. Record it as a Stage C idea.

### Answers to the slice questions

**(a) Complete, robust, arbitrary pockets, holes?** No. The topology and identity machinery is exact and fail-loud, and the cell kinds are general: points and open segments only, **no arcs** (`segment_site_catalog.h:13-16`, `MatTraits = Segment_Delaunay_graph_traits_without_intersections_2` at `segment_site_parameterization.h:22-25`). The assembled graph exists only for the 6-vertex L fixture (Python-reachable) and the exact rectangle (native tests only). Edges carry exact clearance profiles (L only, `canonical_l_shape_mat_clearance_graph`) and exact parameters. The parameter is not arc length: the parallel S–S chart has an unnormalised direction (`segment_site_parameterization.cpp:214-234`), and the parabola is parameterised by projection onto the directrix (`:1020-1067`). Holes: supported in the catalog and clipping, rejected by both adaptors. Robustness: every degenerate or unexpected configuration throws a named error. That is good fail-loud behaviour, but on an arbitrary pocket it means "throws", not "handles".

**(b) Versus the straight skeleton in `src/toolpath.cpp`:** the skeleton guide already computes **exact clearance at each station** (`toolpath.cpp:185-220`, boundary distance per station). So "correct clearance at reflex vertices" is not what is missing. What is missing is (i) the medial *location* near reflex vertices (parabola vs straight bisector, i.e. larger available loops along the true MAT), (ii) clearance as a function rather than samples, (iii) necks with separating cuts, (iv) exact `clearance ≥ r` cut-off points, (v) zero-guide slots. Neither gives a branch tree or traversal order. The MAT graph is an undirected edge list with node-site CSR (`segment_site_graph_csr.cpp:407-564`).

**(c) Consumers.** Native → nanobind `_medial_axis_2.SegmentSiteMedialAxis` and `segment_site_medial_axis` (`src/medial_axis_2.cpp:320-382`, `:418-471`). → Python `compas_cgal.adaptive.medial_axis.MedialAxis.build` (`src/compas_cgal/adaptive/medial_axis.py:599`). → `adaptive/traversal.py`, `traversal_graph.py`, `neck_topology.py`, `replay.py:389-409`, and `generator.py` (`generate_exact_adaptive_continuation`, `:1336`). Callers of the adaptive generator: tests only (`tests/adaptive/test_generator.py`, `test_acceptance.py`, `test_replay.py`, etc.). There is no benchmark, no pixi task, and no gated generator. Held benchmarks: none import `_medial_axis_2`. The `held-native-curved-circle` pixi task only rebuilds the module. **verified (grep).** Build: every `segment_site_*.cpp` (spikes included) is compiled into `_medial_axis_2` and into a native test target (`CMakeLists.txt:415-441`, `:496-532`).

**(d) Known limits and cost.**
- One pocket, hard-coded (`segment_site_catalog_graph.cpp:172-192`), plus exact cardinality asserts (`:1159`, `:1577`). **verified (probe).**
- No holes in any adaptor (`:130`, `:164`). **verified.**
- No arcs as sites. **verified.**
- P–P bisector rays that intersect the domain throw (`segment_site_catalog_graph.cpp:1316-1319`). A pocket with two facing reflex vertices has an interior P–P edge, so it would fail. **verified by code reading.**
- P–P ray endpoint must be a rational point (`segment_site_endpoint_binding.cpp:1626-1646`). **verified by code reading.**
- Sampling only on rational charts (`segment_site_catalog_sampling.cpp:122-125`, `:137-139`). **verified by code reading.**
- `SegmentSiteMatBundle2::build` always routes to the L path. The rectangle path is unreachable from Python (`segment_site_mat_bundle.cpp:39`). **verified.**
- Cost: 0.125 s for the fixture at r = 0.5. **measured.** Generic-coordinate cost: **unverified** (never runnable).

**(e) Should the tree walk be built on this MAT?** Not for Stage A. It cannot run on any gate pocket today. On the two rectangles it is the same graph as the skeleton, and the plan's failures there are ordering and emission defects. But the plan should be written so that the guide is an interface (`GuideTree` over typed nodes and edges, with clearance per edge) and not hard-wired to `skeleton_guide`. The MAT is the principled guide for reflex-corner pockets (the gate `L_shape`, the Held pockets) and the only exact neck oracle. Concretely, making the MAT the guide requires:
1. **One general catalog adaptor** replacing `canonical_rectangle_mat_graph` and `canonical_l_shape_mat_graph`. It would be the L-shape loop body (`segment_site_catalog_graph.cpp:1262-1575`) minus `require_l_shape_input` and the cardinality block, with `canonicalize_nodes` kept. Then route `SegmentSiteMatBundle2::build` through it. The rectangle golden and the L golden become regression tests of the general adaptor.
2. **P–P cells inside the domain:** clip them with `clip_linear_clearance_components` and emit them as `LINE` edges instead of throwing (`:1316`). Lift the rational-endpoint requirement in `bind_point_point_ray_parameterization` (`segment_site_endpoint_binding.cpp:1641`) by binding through the algebraic kernel as the other binders do.
3. **Non-rational charts in sampling:** sample nonparallel S–S edges in Q(√α), e.g. a rational parameter with `to_double` only at reporting, which keeps exact decisions. Do not throw at `segment_site_catalog_sampling.cpp:122`.
4. **Composite cells** (one generator pair split by a raw P–S transition): today only a spike (`segment_site_mat.cpp:1668-2047`, asserts one transition at `:1828`).
5. **Holes:** remove the `holes.empty()` gates. The neck cut index already handles cycles (`segment_site_neck.cpp:221-330`: a vertex or edge on a cycle is "non-separating" and dropped). The tree walk must then handle a graph that is not a tree.
6. **Root and traversal export:** a native query "argmax squared clearance over nodes and plateaus" using `MatExactSquaredWidth2` compare. Export node coordinates as exact-backed reporting doubles. Today node coordinates come from the first incident sample (`segment_site_mat_numeric_table.cpp:65-83`).
7. **Generic-coordinate witness test:** build on a generic-double polygon and time it, per the CORE-stall memory, before relying on it in a product path.
8. **Gate pockets as acceptance:** `rect_12x8`, `rect_20x12` and the gate `L_shape` must build, and the rectangle graph must match the straight skeleton's topology (equal node count and adjacency) as a contract test.

### Measurements

```
$ pixi run python -c "... segment_site_medial_axis(np.array(v), [], 0.5, 0.75, 0.02, 32) ..."
L_fixture_0..6 OK nodes 10 edges 9 necks 2 0.125s
L_scaled_0..12 UnsupportedCanonicalMatLShapeGraphError canonical six-edge ring is not the exact L-shape fixture 0.011s
rect_12x8 UnsupportedCanonicalMatLShapeGraphError canonical L-shape MAT graph requires the six-edge fixture ring 0.006s
rect_20x12 UnsupportedCanonicalMatLShapeGraphError canonical L-shape MAT graph requires the six-edge fixture ring 0.006s
$ pixi run python -c "... gate L_shape (0,0),(12,0),(12,6),(6,6),(6,10),(0,10), r=1.0; fixture at r=0.5/1.0/0.25 ..."
UnsupportedCanonicalMatLShapeGraphError canonical six-edge ring is not the exact L-shape fixture
fixture r 0.5 edges 9 necks 2 samples 27
fixture r 1.0 edges 5 necks 0 samples 18
fixture r 0.25 edges 9 necks 2 samples 31
```
(The gate pocket coordinates come from `benchmarks/gate.py:84-89,121-124` with tool diameter 2.)

---

### src/segment_site_catalog.cpp | lines 345 | read 1-345
- **role:** guide (kernel: site catalog)
- **used by:** segment_site_catalog_graph.cpp, catalog_neck, catalog_sampling, mat_bundle, mat_proposal_table; `_medial_axis_2`
- **state:** live
- **generation:** builds the canonical site set: every ring vertex becomes a point site and every ring edge an open-segment site, outer ring then holes (`:288-341`). Stable ids are hashed from (kind, ring, feature). Adequate and correct for polygons with holes. No arcs.
- **defects:** none found.
- **unused capability / ideas:** handles holes (`:333-341`), but no downstream adaptor accepts them.
- **claims:** "canonical MAT site identities are strictly ordered" is enforced (`:76-87`). verified.

### src/segment_site_catalog.h | lines 106 | read 1-106
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** `MatSiteKind2 {Point, OpenSegment}` (`:13-16`): no arc site kind, so arcs in a pocket boundary cannot be MAT sites.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_catalog_graph.cpp | lines 1631 | read 1-1631
- **role:** guide (topology assembly)
- **used by:** catalog_neck.cpp (L path, then bundle, then Python); rectangle path used only by native tests (`tests/native/task9_catalog_graph_gate.cpp:680`, `task9_graph_csr_gate.cpp:157`)
- **state:** live (L), test-only (rectangle)
- **generation:** walks CGAL's Voronoi adaptor halfedges. For each generator pair it builds an exact chart (parallel S–S, nonparallel S–S, P–S parabola), binds the two live endpoints to exact algebraic parameters, and clips to the domain and to `clearance > r`. It then emits graph components and canonicalises nodes by their generator-feature set. This is the right algorithm. Quality: **broken for product use**, because it is fixture-gated.
- **defects:**
  - `:172-192`: `require_l_shape_input` hard-codes the literal L ring. Any other 6-gon throws, so the gate L_shape and every scaled copy are unsupported (measured).
  - `:130-157`: the rectangle path accepts only an exact 4-edge, no-hole rectangle.
  - `:1159-1165`, `:1577-1599`: exact cell-count asserts (5/8; 1/8/12/2/8/4) make each adaptor a golden of one input.
  - `:130`, `:164`: no holes.
  - `:1316-1319`: a P–P ray that meets the domain throws, but P–P edges are interior for pockets with facing reflex vertices.
  - Duplicate code: `register_*` / `*_vertex_alias` exist twice (`:477-701` vs `:549-819`), and again in segment_site_mat.cpp. This is maintenance, not output.
- **unused capability / ideas:** the L loop body (`:1262-1575`) already handles all five cell kinds. With the fixture gate and the count asserts removed, it is the general adaptor.
- **claims:** none in the file. The doc maturity claim (`docs/segment_site_mat.md:201-204`, "arbitrary-pocket traversal not complete") is verified here.

### src/segment_site_catalog_graph.h | lines 56 | read 1-56
- **role:** guide (API)
- **used by:** catalog_graph.cpp, catalog_neck.cpp, `medial_axis_2.cpp:2`
- **state:** live
- **generation:** declares `canonical_rectangle_mat_graph` and `canonical_l_shape_mat_graph`. The error names (`UnsupportedCanonicalMatLShapeGraphError`) expose the fixture specialisation in the public error surface (re-exported in `medial_axis_2.cpp:401-402`).
- **defects:** the API has no general entry point.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_catalog_neck.cpp | lines 241 | read 1-241
- **role:** guide (clearance profile per edge)
- **used by:** catalog_sampling.cpp, mat_bundle (L path)
- **state:** live
- **generation:** exact squared clearance per edge as a rational polynomial in the chart parameter. Parallel S–S gives a constant (`:72-95`), nonparallel S–S gives `k·t²` (`:97-110`), and P–S gives a quartic (`:112-145`). Each is cross-checked against the clipping polynomial. Excellent: a closed-form clearance along the guide.
- **defects:**
  - `:228-241`: only `canonical_l_shape_mat_clearance_graph` exists, so no profile graph exists for any other pocket.
  - `:177`: P–P edges have no profile ("unsupported clearance-source pair"), so a general adaptor needs a P–P profile, which is trivially `|t·d + m − p|²`.
- **unused capability / ideas:** clearance(t) in closed form would let station placement and loop radius be exact functions instead of sampled boundary distances (`toolpath.cpp:219`).
- **claims:** none.

### src/segment_site_catalog_neck.h | lines 37 | read 1-37
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** `MatClearanceProfileGraph2::build` validates one profile per canonical edge. Adequate.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_catalog_sampling.cpp | lines 298 | read 1-298
- **role:** guide (reporting samples along MAT edges)
- **used by:** mat_bundle.cpp:39, mat_proposal_table.cpp:111
- **state:** live (L only)
- **generation:** turns each edge into a rational sampling curve and samples it by spacing and sagitta. Adequate for reporting.
- **defects:**
  - `:122-125`, `:137-139`: non-rational charts throw. A nonparallel S–S chart is rational only when `|n1|²|n2|²` is a perfect square (`segment_site_parameterization.cpp:399-410`), so generic corner angles fail.
  - `:89-97`: `branch_segment` uses `CORE::sqrt` to rebuild the live line, which is a CORE-stall risk on generic input (unverified).
  - `:223-226`: `sample_exact_flags` are set to {1,1} for every sample by construction. No per-sample exact check backs the "verified" flag. The flag means "lies on an exactly clipped edge", not "sample checked".
  - `:277-298`: L-only entry point.
- **unused capability / ideas:** none beyond the MAT itself.
- **claims:** "verified" sample flags: unverified as a per-sample statement (see `:223`).

### src/segment_site_catalog_sampling.h | lines 74 | read 1-74
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** the `MatProposalSamplingGraph2` owner. The friend function is the L entry point, so construction is hard-wired to the fixture (`:65-68`).
- **defects:** API coupling to the L fixture.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_clipping.cpp | lines 1640 | read 1-1640
- **role:** guide (exact clipping to domain D and to clearance ≥ r)
- **used by:** catalog_graph, segment_site_mat.cpp, catalog_neck
- **state:** live
- **generation:** for every chart kind it computes exact domain-boundary roots (linear `:95-267`, parabola `:269-431`, nonparallel via parameterization) and clearance roots (`:447-509`). It then splits into open cells and keeps cells that are inside D and have clearance > r² (`:1149-1169`), plus isolated touch points. This is exactly the tool-centre domain `C_r` intersected with the MAT. Excellent, and exact.
- **defects:**
  - `:129-133`, `:307-314`: a primitive collinear with a domain edge throws `OverlappingDomainBoundaryError`. Fail-loud, but a generic pocket with collinear features would throw.
  - The file is 1640 lines, above the 1.5k refactor bar.
- **unused capability / ideas:** exact "where can the tool centre be on the guide" intervals with exact cut-off points. This replaces the `min_centerline()` station filter in `toolpath.cpp:230`.
- **claims:** none.

### src/segment_site_clipping.h | lines 189 | read 1-189
- **role:** guide (API)
- **used by:** as .cpp
- **state:** live
- **generation:** one named error per failure mode. Adequate.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_delaunay.cpp | lines 259 | read 1-259
- **role:** kernel (CGAL SDG construction plus site bijection)
- **used by:** catalog_graph
- **state:** live
- **generation:** inserts every segment site into `Segment_Delaunay_graph_2` and checks that the live vertices biject with the catalog (`:208-247`). Correct, with exact `compare_xy` keys.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** "catalog-fed graph is bijective" is enforced (`:240-246`). verified.

### src/segment_site_delaunay.h | lines 90 | read 1-90
- **role:** kernel (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_endpoint_binding.cpp | lines 3403 | read 1-3403
- **role:** guide (binds each live CGAL Voronoi vertex to an exact algebraic chart parameter)
- **used by:** catalog_graph, segment_site_mat.cpp
- **state:** live
- **generation:** for each cell endpoint it solves the limiter equation exactly (point limiter, segment limiter, parallel, nonparallel; quartic at most), factors it, and keeps the unique root that lies in the open feature. It also requires that root to equal the live CGAL vertex as a `CORE::Expr` point. The topology is therefore CGAL's, and the coordinates are exact algebraic numbers. This is the mechanism that makes the graph exact and is well built.
- **defects:**
  - `:795-807`, `:992-1005`, `:2432-2441`, `:2572-2591`: identity decisions on `CORE::Expr` (`parameter != live_parameter`, `candidate != live_point`) are the pattern known to stall CORE on generic doubles. Integer fixtures only so far, so the cost is **unverified**.
  - `:1641-1646`: the P–P ray endpoint must be rational, otherwise it throws.
  - The file is 3403 lines, more than twice the 1.5k refactor bar. Many overloads differ only in the identity resolver (`:2903-3403`).
- **unused capability / ideas:** none beyond the MAT.
- **claims:** `// Unreachable for validated adaptor edges` (`:922`) is unverified for general pockets.

### src/segment_site_endpoint_binding.h | lines 358 | read 1-358
- **role:** guide (API)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** overload sprawl, with four `bind_nonparallel_…` and four `bind_parallel_…` signatures (`:259-358`).
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_graph_csr.cpp | lines 564 | read 1-564
- **role:** guide (numeric export of the graph)
- **used by:** mat_certificate.cpp:423, mat_proposal_table.cpp:115
- **state:** live
- **generation:** exports node→site CSR, edges as int rows (source, target, curve kind, sites, dual kind, dual index, component) and endpoint flags/features. It validates provenance containment (`:170-223`). Adequate. This is the branch graph a traversal would consume, but it has no orientation or root.
- **defects:** `:560-561`: `edge_exact_flags` is always {1,1,1}, a constant flag set with no information.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_graph_csr.h | lines 119 | read 1-119
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** curve kinds LINE, PARABOLA, CIRCLE (`segment_site_provenance.h:37-41`). Only LINE and PARABOLA are emitted (`graph_csr.cpp:109-120`).
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_graph_emission.cpp | lines 135 | read 1-135
- **role:** guide (component to graph edge)
- **used by:** catalog_graph, segment_site_mat.cpp
- **state:** live
- **generation:** drops zero-length components (`:3-32`) and appends edges and nodes keyed by exact root identity. Adequate.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_graph_emission.h | lines 30 | read 1-30
- **role:** guide (API)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat.cpp | lines 3503 | read 1-3503
- **role:** guide (point-site graph plus spike fixtures)
- **used by:** `exact_point_site_graph` and about 20 `*_spike()` functions. Consumers are native tests (`tests/native/task9_point_graph_gate.cpp`, `task9_segment_segment_gate.cpp`, `task9_segment_limiter_gate.cpp`). It is compiled into `_medial_axis_2` (`CMakeLists.txt:497`) but not called from the bundle.
- **state:** test-only fixtures inside a production module
- **generation:** none for product paths. The general-looking pieces (`append_point_site_graph` `:357-579`, `append_point_segment_graph` `:817-997`, `append_parallel_segment_segment_graph` `:1110-1272`, `append_nonparallel_segment_segment_graph` `:1668-2047`) each rebuild a private SDG from a hand-listed site set and assert fixture cardinalities, e.g. `:1828-1832` "nonparallel fixture has no unique retained transition".
- **defects:**
  - The file is 3503 lines, more than twice the refactor bar. Spike fixtures with literal coordinates (`:2098-3503`) ship in the production `.so`.
  - The rectangle and other aliasing helpers duplicate catalog_graph.cpp (`:2538-2695` vs `segment_site_catalog_graph.cpp:477-701`).
- **unused capability / ideas:** the composite nonparallel cell logic (`:1668-2047`) is the only implementation of a split S–S pair. A general adaptor needs it (item 4 in (e)).
- **claims:** none.

### src/segment_site_mat.h | lines 209 | read 1-209
- **role:** guide (types and spike API)
- **used by:** catalog_graph.h, graph_csr.h, neck.h, …
- **state:** live
- **generation:** `MatExactGraph2` with nodes and edges carrying generator and parent site ids, exact endpoints, clip component and an admissible flag (`:86-112`). This is the graph a guide would consume.
- **defects:** about 30 spike declarations in a public header (`:114-209`).
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_bundle.cpp | lines 175 | read 1-175
- **role:** guide (production owner)
- **used by:** `medial_axis_2.cpp:328` → `SegmentSiteMedialAxis.build` → `adaptive/medial_axis.py:599`
- **state:** live
- **generation:** the single production entry point. It builds the reachable domain identity, the L sampling graph, the certificate, neck evidence and zero-guide records.
- **defects:** `:39`: always calls `canonical_l_shape_mat_proposal_graph`. This line makes the Python MAT one-pocket-only (measured: every gate pocket raises).
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_bundle.h | lines 66 | read 1-66
- **role:** guide (types)
- **used by:** as .cpp, medial_axis_2.cpp
- **state:** live
- **generation:** n/a
- **defects:** none beyond the .cpp.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_certificate.cpp | lines 451 | read 1-451
- **role:** validation (canonical bytes and sha256 of the exact graph, for replay)
- **used by:** mat_bundle, zero_guide; Python `validate_*` calls
- **state:** live
- **generation:** n/a. It encodes sites, nodes, edges, profiles, the numeric projection and the neck projection, and replays them by rebuild-and-compare (`:439-451`).
- **defects:** none found. It does not constrain generation.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_certificate.h | lines 66 | read 1-66
- **role:** validation (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_numeric_table.cpp | lines 124 | read 1-124
- **role:** guide (reporting coordinates and neck cut table)
- **used by:** mat_bundle
- **state:** live
- **generation:** node coordinates are the first incident edge sample (`:65-83`), i.e. double reporting. The neck cut edge list comes from separating-cut partitions (`:97-124`).
- **defects:** `:72-82`: node coordinates come from the sampler and not from the exact node. They are reporting-only, but a traversal using them as guide vertices gets sampler doubles.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_numeric_table.h | lines 59 | read 1-59
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** `:55-59`: the L-only helper is in the public API.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_proposal_table.cpp | lines 119 | read 1-119
- **role:** guide (per-sample clearance and suggested guide-circle radius)
- **used by:** mat_bundle
- **state:** live
- **generation:** for each sample it gives the double clearance and `guide_radius = (clearance − r)/2` (`:87-88`).
- **defects:** `:88`: the `/2` is undocumented. `toolpath.cpp:205-207` uses `mat_scale·(clearance − r − radial_clearance)`. Two different radius-from-clearance rules exist with no shared derivation (unverified which is intended).
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_proposal_table.h | lines 71 | read 1-71
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** `MatToolRadiusMm2::build` injects the double radius exactly (`proposal_table.cpp:31-38`). This is correct boundary doctrine.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_sampling.cpp | lines 453 | read 1-453
- **role:** guide (spacing and sagitta sampling of LINE and PARABOLA edges)
- **used by:** catalog_sampling
- **state:** live
- **generation:** lines are split uniformly by double chord length (`:191-219`). Parabolas are bisected at exact rational split parameters until the chord is ≤ spacing and the midpoint deviation is ≤ sagitta (`:155-189`). Adequate for reporting.
- **defects:** `:171-174`: exceeding the refinement depth throws. Fail-loud is fine.
- **unused capability / ideas:** none.
- **claims:** `:125-128` (chord deviation peaks at the parameter midpoint for a quadratic map) is verified mathematically: deviation ∝ λ(1−λ).

### src/segment_site_mat_sampling.h | lines 146 | read 1-146
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** typed `MatStationSpacingMm2` and `MatSagittaBoundMm2`. Good.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_mat_spikes.cpp | lines 108 | read 1-108
- **role:** test (compile spike)
- **used by:** `tests/native/task9_mat_compile_gate.cpp`
- **state:** test-only, compiled into `_medial_axis_2` (`CMakeLists.txt:519`)
- **generation:** n/a
- **defects:** test code in the production module.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck.cpp | lines 497 | read 1-497
- **role:** guide (neck separating cuts)
- **used by:** neck_evidence.cpp
- **state:** live
- **generation:** for an edge, vertex or plateau candidate it removes the candidate and partitions the remaining connected component. It throws `NonSeparating…` if fewer than two sides remain (`:221-330`). On a tree every candidate separates. On a graph with cycles, cycle candidates are dropped. Correct.
- **defects:** none found.
- **unused capability / ideas:** the side partitions name exactly which guide edges lie beyond a neck. That is the "far room needs its own entry" decision the plan's Task 3 makes heuristically.
- **claims:** none.

### src/segment_site_neck.h | lines 132 | read 1-132
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck_classification.cpp | lines 144 | read 1-144
- **role:** guide (neck width class against exact thresholds)
- **used by:** mat_bundle → Python `validate_and_classify_necks`
- **state:** live
- **generation:** compares each neck's exact squared width to strictly increasing rational boundaries and returns the class id plus a comparison certificate (`:73-116`). Excellent and exact.
- **defects:** none found.
- **unused capability / ideas:** use as the gate between "pass through with trochoids" and "separate region, needs its own entry" in a tree walk.
- **claims:** none.

### src/segment_site_neck_classification.h | lines 20 | read 1-20
- **role:** guide (API)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck_clearance.cpp | lines 761 | read 1-761
- **role:** guide (clearance profile analysis)
- **used by:** catalog_neck, neck_evidence
- **state:** live
- **generation:** exact strict interior minima of squared clearance via derivative roots and sign witnesses (`:698-761`), endpoint inward-sign (`:232-268`), and width = 4·clearance² as an algebraic number via a 2-D solve (`:335-420`). Excellent.
- **defects:** `:13-14`: degree capped at 4. That is fine for P–S and S–S but must grow for P–P (degree 2, fine) and any future arc sites.
- **unused capability / ideas:** exact minima of clearance along the guide. This lets a generator place a station exactly at the bottleneck instead of finding it by sampling.
- **claims:** none.

### src/segment_site_neck_clearance.h | lines 192 | read 1-192
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck_evidence.cpp | lines 554 | read 1-554
- **role:** guide (neck detection)
- **used by:** neck_evidence_bytes → mat_certificate, mat_bundle
- **state:** live
- **generation:** collects four kinds of neck: strict interior edge minima, plateaus of constant-clearance edges, shared-vertex minima (≥2 incident edges increasing away), and clearance-endpoint necks. Each must separate the graph (`:363-554`). It also checks that incident profiles agree on width at every node (`:400-417`), which is a strong consistency check. Excellent concept, exact.
- **defects:** `:106-133`: `NonSeparatingMatNeckCandidateError` is caught as control flow. It is a typed catch, so it is not a silent swallow. Acceptable.
- **unused capability / ideas:** see the Slice summary. This is the exact neck oracle the plan lacks.
- **claims:** "necks require exact local clearance minima and a separating graph cut" (`docs/segment_site_mat.md`, "Why a separate exact MAT exists") is verified here.

### src/segment_site_neck_evidence.h | lines 124 | read 1-124
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck_evidence_bytes.cpp | lines 224 | read 1-224
- **role:** validation (canonical neck evidence bytes and replay)
- **used by:** classification, certificate, numeric_table
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_neck_evidence_bytes.h | lines 36 | read 1-36
- **role:** validation (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_parameterization.cpp | lines 1653 | read 1-1653
- **role:** guide (exact charts and field arithmetic)
- **used by:** all segment_site geometry files
- **state:** live
- **generation:** canonical line supports (`:124-167`); parallel S–S chart (`:169-235`); nonparallel S–S chart over Q(√(|n1|²|n2|²)), with the branch picked by exact `has_on` against the live CGAL segment and collapse to rational on perfect squares (`:264-422`); quadratic-field sign and compare via one squaring (`:876-928`, correct); P–S parabola chart (`:986-1068`); exact point-in-polygon in Q(√α) (`:1146-1202`). Excellent: kernel-level exactness with no epsilons.
- **defects:**
  - Parameters are not arc length (`:214-234`, `:1045-1067`). Any consumer needing distance along the guide must integrate.
  - `:304-305`: branch selection builds `CORE::sqrt` Expr lines, a generic-input cost risk (unverified).
  - The file is 1653 lines, above the 1.5k refactor bar.
- **unused capability / ideas:** `quadratic_field_sign` is reusable for any exact one-root decision (the project's Sqrt_extension doctrine).
- **claims:** none.

### src/segment_site_parameterization.h | lines 391 | read 1-391
- **role:** kernel (MatKernel = Epeck_with_sqrt (CORE); MatTraits = SDG traits without intersections)
- **used by:** all
- **state:** live
- **generation:** fixes the kernel at `:20-25`. Segment sites without intersections only, so no arcs and no self-intersecting input.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_provenance.cpp | lines 503 | read 1-503
- **role:** guide (identity: root ids, dual ids, node ids, evidence union)
- **used by:** all graph builders
- **state:** live
- **generation:** n/a for geometry. It makes node identity independent of coordinates (canonical factor plus root ordinal, `:36-74`). That is what lets nodes from different charts merge exactly.
- **defects:** `:108-133` vs `:135-169`: two near-identical bijection functions, one returning 0 and one throwing.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_provenance.h | lines 240 | read 1-240
- **role:** guide (SDG/Voronoi typedefs, endpoint and evidence types)
- **used by:** all
- **state:** live
- **generation:** `SegmentSiteVoronoi2` with the degeneracy-removal policy (`:19-30`). `MatEndpointExactEvidence2` records why an endpoint exists: Voronoi vertex, domain boundary, or clearance root.
- **defects:** none found.
- **unused capability / ideas:** endpoint evidence tells a traversal whether an edge end is a real branch point, a wall hit or a tool-fit limit. That is useful for deciding where a walk must turn back.
- **claims:** none.

### src/segment_site_rational_sources.cpp | lines 178 | read 1-178
- **role:** kernel (exact rational injection of catalog coordinates)
- **used by:** catalog_graph, catalog_neck, catalog_sampling, endpoint_binding
- **state:** live
- **generation:** n/a. `exact_mat_input_rational` checks that a CORE coordinate is exactly a BigRat (`:47-62`), which is the correct boundary doctrine.
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_rational_sources.h | lines 83 | read 1-83
- **role:** kernel (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_voronoi.cpp | lines 75 | read 1-75
- **role:** kernel (Voronoi adaptor transfer with bijection check)
- **used by:** catalog_graph
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_voronoi.h | lines 44 | read 1-44
- **role:** kernel (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.

### src/segment_site_zero_guide.cpp | lines 159 | read 1-159
- **role:** guide (tool-width slot inventory)
- **used by:** mat_bundle → Python `MatZeroGuideInventory` → `adaptive/generator.py:148-153`
- **state:** live
- **generation:** edges whose squared clearance is exactly the constant r² (a slot exactly one tool wide) get a replayable record. Adequate: it detects where only a straight slot (zero trochoid radius) fits.
- **defects:** none found.
- **unused capability / ideas:** marks slot-only passages exactly. A product generator could treat them as forced-slot links instead of hitting the ρ > r floor.
- **claims:** "radius-1 L fixture emits exactly two replayable zero-guide records" (`docs/segment_site_mat.md:135`) is unverified: not re-run, because the probe printed edge and neck counts only.

### src/segment_site_zero_guide.h | lines 59 | read 1-59
- **role:** guide (types)
- **used by:** as .cpp
- **state:** live
- **generation:** n/a
- **defects:** none found.
- **unused capability / ideas:** none.
- **claims:** none.
