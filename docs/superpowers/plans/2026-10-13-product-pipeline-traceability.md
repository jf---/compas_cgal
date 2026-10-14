> **status: active** — traceability for 2026-10-13-product-pipeline.md

# Product pipeline: census traceability

## Unmapped

None. After the plan's second amendment every finding has a task, constraint, parking line, or a closed, frozen or superseded disposition.

## How ids are derived

- `F-<slice>-<n>`: register rows in `docs/superpowers/census/2026-10-13/out/{d1_docs_specs_reviews,d2_plans}.md`.
- `<slice>-T<n>`, `-U<n>`, `-P<n>`: in each record's `## Slice summary`, the top-level items under "Top findings", under the "capabilities … not use" (g7: "… only off the live branch") heading, and under "wrong or incomplete", numbered in order of appearance. Table body rows count as items where a section is a table (d2-U, g4-U). `<slice>` is the record stem up to its first underscore.
- Summary items criticise the earlier tree-walk draft; each is judged against this plan. `superseded` means this plan no longer contains the criticised element.
- `docs/superpowers/census/2026-10-13/check_traceability.py` derives the ids independently and checks this table.

## Matrix

| id | finding (≤ 20 words) | disposition | where |
| --- | --- | --- | --- |
| F-d1_docs_specs_reviews-1 | auditor P1 Task 4 owes schema v2, native certification, atomic replay | closed | commit 7cda8f6d (src/audit_replay_2.cpp) |
| F-d1_docs_specs_reviews-2 | exact arc surrogate: audit replay has no retained-sliver bound | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-3 | P4 must find the extreme full-suite tail before CI enforcement | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-4 | full-circle closure on arc fixtures exceeds 120 s; no arc-adapter coverage | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-5 | engagement_audit replay/report tests (Task 4D+) never recorded done | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-6 | Held reconstruction recursion: depth ceiling, no shared node budget | task | 2.1 |
| F-d1_docs_specs_reviews-7 | independent Held certificate re-review pending | task | 2.1 |
| F-d1_docs_specs_reviews-8 | biarc children merge-ineligible until sound cubic split | task | 2.1 |
| F-d1_docs_specs_reviews-9 | CHANGELOG lists removed mesh/geodesic/polyline APIs | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-10 | dangling CHANGELOG reference; O(n²) exact depletion in generation deferred | task | 0.2 |
| F-d1_docs_specs_reviews-11 | max_passes truncation only warns instead of failing loud | task | 1.7 |
| F-d1_docs_specs_reviews-12 | CONVENTIONS.md names a superseded stage page | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-13 | engagement_controlled regulates advance not radius; saturates ~141.8° | task | 1.5 |
| F-d1_docs_specs_reviews-14 | center_domain ~3,500× slower on oblique edges, mechanism unknown | task | 1.3 |
| F-d1_docs_specs_reviews-15 | benchmarks/results/ layout no longer enforced | task | 0.3 |
| F-d1_docs_specs_reviews-16 | Held qualification passes on non-empty only, never gate-scored | task | 2.1 |
| F-d1_docs_specs_reviews-17 | certified adaptive generator incomplete (retrace, closure, sealing, arbitrary pockets) | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-18 | MAT DFS route is not a physical walk; needs transit through cleared lineage | task | 2.3 |
| F-d1_docs_specs_reviews-19 | continuous_tea_2 event discovery incomplete, can abort process | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-20 | task13f continuation test stale and now failing | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-21 | production oracle on source v3; v4 one-root migration unvalidated | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-22 | generator/replay must consume MatTraversalState.neck_scope | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-23 | swept-prefix theorem certifies straight advancing cuts at π; unused | task | 2.3 |
| F-d1_docs_specs_reviews-24 | furthest-progress, largest-radius first-feasible search unused | task | 1.5 |
| F-d1_docs_specs_reviews-25 | ship the sampled+guarded certifier; continuous_tea_2 is research | superseded | Task 4.3 (guarded lemma removed) |
| F-d1_docs_specs_reviews-26 | shear transformation and resultant degree are the only cost lever | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-27 | arrangement prototype builds with removed conda env and CGAL 5.6.1 | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-28 | persistent arrangement redesign superseded by zone-in-Gps fix | closed | findings.md:108 LANDED; src/engagement_2.cpp:310 |
| F-d1_docs_specs_reviews-29 | exact annular-sector sweep deferred; depletion uses disk chains | task | 1.6 |
| F-d1_docs_specs_reviews-30 | Gps copy aliases source traits (use-after-free) | closed | commits 8d2b76cc, df6fb020, 57482731 |
| F-d1_docs_specs_reviews-31 | truthful audit rejects ramps, so helical entry inadmissible to it | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-32 | truthful audit remaining: aggregation, consumer migration, evidence matrix | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-33 | exact partial-arc depletion proves no retained-sliver bound | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-34 | chain entries are full slots; no helical or ramped entry | task | 1.2 |
| F-d1_docs_specs_reviews-35 | bisection assumes monotone engagement; ladder required | constraint | Global Constraints: no bisection |
| F-d1_docs_specs_reviews-36 | cap guaranteed only at 33 positions; 43.4° between them at 40° cap | task | 4.1 |
| F-d1_docs_specs_reviews-37 | guide tangency G1 only on straight guide segments | task | 1.6 |
| F-d1_docs_specs_reviews-38 | coverage measured on a grid, not certified | parked | exact full-motion coverage replay |
| F-d1_docs_specs_reviews-39 | no local capsule depletion (seed completeness unsolved) | task | 0.2 |
| F-d1_docs_specs_reviews-40 | Manufacturing_policy (tolerance, angular surrogate, min feature) never built | task | 3.1 |
| F-d1_docs_specs_reviews-41 | exactness checklist count drift | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-42 | exactness.md stale on disk-chain arc sweeps | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-43 | chains greedy nearest-neighbour, never reversed, one plunge each | task | 1.1 |
| F-d1_docs_specs_reviews-44 | restore true medial incidence, qm diameter, ordered paired traversal | task | 2.1 |
| F-d1_docs_specs_reviews-45 | predecessor-only model cannot certify convex-corner stationary event | task | 2.2 |
| F-d1_docs_specs_reviews-46 | live Held shared-thinning regression red at three scales | frozen | Frozen: live Held reproduction |
| F-d1_docs_specs_reviews-47 | exact full-motion coverage replay does not finish | parked | exact full-motion coverage replay |
| F-d1_docs_specs_reviews-48 | conditioning audit of other Epeck_with_sqrt lanes; emitted-circle representation | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-49 | report reachable-material target alongside full-design residuals | task | 1.3 |
| F-d1_docs_specs_reviews-50 | Monstera interior gap; live selectors omit strongest circle | frozen | Frozen: live Held reproduction |
| F-d1_docs_specs_reviews-51 | stale runtime-blocker text contradicts same page | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-52 | Held machinable target (convex arcs regularised) never built | task | 2.2 |
| F-d1_docs_specs_reviews-53 | native contour-aware qm-diameter circle construction unused | task | 2.5 |
| F-d1_docs_specs_reviews-54 | boundary-order spacing: 2.4× fewer circles, unused | task | 2.4 |
| F-d1_docs_specs_reviews-55 | live Held initial placement dominates time (Monstera 446.5 s) | frozen | Frozen: live Held reproduction |
| F-d1_docs_specs_reviews-56 | stale Held primitive counts | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-57 | contour-aware, MATHSM, Fig 6 curves blocked on live | task | 2.4 |
| F-d1_docs_specs_reviews-58 | concave-corner approach/shared-contact/departure unused | task | 2.4 |
| F-d1_docs_specs_reviews-59 | CGAL licence gate: LGPL repo instantiates GPL packages | parked | Licensing |
| F-d1_docs_specs_reviews-60 | local_depletion.md stale | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-61 | sound local capsule depletion needs custom zone visitor | task | 0.2 |
| F-d1_docs_specs_reviews-62 | clearance ≤ 3r stations cannot carry a trochoid; need non-trochoid motion | task | 1.3 |
| F-d1_docs_specs_reviews-63 | centre plunge + concentric ramp entry (240° floor) unused | superseded | Task 1.2 (R3 helical bore) |
| F-d1_docs_specs_reviews-64 | helical entry needs depletion change; replay rejects Z+XY | task | 1.2 |
| F-d1_docs_specs_reviews-65 | do not exclude entry motion from engagement-step criterion | task | 1.2 |
| F-d1_docs_specs_reviews-66 | external common-tangent bridge removes tilt; open for gated generators | task | 1.1 |
| F-d1_docs_specs_reviews-67 | stay-down chain-ordering generator built then deleted | task | 1.1 |
| F-d1_docs_specs_reviews-68 | engagement not translation-invariant at exact rim tangency | parked | translation variance (C1) |
| F-d1_docs_specs_reviews-69 | straight plunge is the defect; helical entry never landed | task | 1.2 |
| F-d1_docs_specs_reviews-70 | run boundary = rapid or path_index change; stay-down links distinct | closed | test_two_chains_linked_at_cutting_depth_are_still_two_guides (tests/benchmarks/test_quality.py:636) |
| F-d1_docs_specs_reviews-71 | uncut residue not corner-confined (3.28 from any vertex) | task | 1.4 |
| F-d1_docs_specs_reviews-72 | guide revisits cleared stations: 10 zero-engagement loops | task | 1.1 |
| F-d1_docs_specs_reviews-73 | W ≤ 4r channels emitted as 68 plunges | constraint | Review Focus 2 (Tasks 1.7, 2.3) |
| F-d1_docs_specs_reviews-74 | no lower engagement bound or chip-thinning compensation | task | 1.5 |
| F-d1_docs_specs_reviews-75 | helical or pre-drilled entry would clear cap and step alone | task | 1.2 |
| F-d1_docs_specs_reviews-76 | rect_20x12 slotting: straight cut out of degenerate corner loop | task | 1.3 |
| F-d1_docs_specs_reviews-77 | material_entries and loop_radius_cv deliberately ungated | task | 0.3 |
| F-d1_docs_specs_reviews-78 | engagement non-monotone in spacing and radius: never bisect | constraint | Global Constraints: no bisection |
| F-d1_docs_specs_reviews-79 | number-type stage 2 leaves invariant 4 undelivered | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-80 | filter-failure rate unmeasured; never run under GMP | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-81 | exact decisions need generic-double fixtures and a time bound | constraint | Review Focus 5 |
| F-d1_docs_specs_reviews-82 | AttestationByteDriftError declared, never raised | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-83 | add oblique pocket to corpus; profile center_domain | task | 1.3 |
| F-d1_docs_specs_reviews-84 | performance/parity claims hold only for axis-aligned pockets | task | 1.3 |
| F-d1_docs_specs_reviews-85 | partial-arc trochoid abandoned; full circles recut ~90% | task | 1.6 |
| F-d1_docs_specs_reviews-86 | DP cost-optimal transition solver removed, not revisited | task | 3.4 |
| F-d1_docs_specs_reviews-87 | external-tangent 1e-9 clamp (status unknown) | task | 1.6 |
| F-d1_docs_specs_reviews-88 | SP2 contour-parallel offset-loop clearing never built | superseded | R2 (two regimes on one guide) |
| F-d1_docs_specs_reviews-89 | engagement target band (≥95% length in band) never built | task | 1.5 |
| F-d1_docs_specs_reviews-90 | SP2b spiralization of adjacent offset loops never built | superseded | R2 (two regimes on one guide) |
| F-d1_docs_specs_reviews-91 | SP3 chip-thinning feed post-processor never built | task | 3.2 |
| F-d1_docs_specs_reviews-92 | entry per region at max-clearance point with certified slot | task | 1.2 |
| F-d1_docs_specs_reviews-93 | future: rest-machining chains, G2 blending, true-MA guides | task | 1.3, 1.6, 2.5 |
| F-d1_docs_specs_reviews-94 | land local zone engagement query | closed | src/engagement_2.cpp:310 engaged_arcs_zone |
| F-d1_docs_specs_reviews-95 | SP1 wall-clock exit gate test absent | task | 0.3 |
| F-d1_docs_specs_reviews-96 | baseline engagement audit report (SP2 reference) absent | superseded | Task 0.3 (gate_vectors.json baseline) |
| F-d1_docs_specs_reviews-97 | certify_segment_tea guarded certificate unused by generators | superseded | Task 4.3 (lemma removed) |
| F-d1_docs_specs_reviews-98 | residual extraction deferred to SP2 | task | 1.4 |
| F-d1_docs_specs_reviews-99 | at forced stations every candidate is over cap | constraint | Global Constraints: forced decisions |
| F-d1_docs_specs_reviews-100 | bridge cuts unregulated; worst 141.8° | task | 1.1 |
| F-d1_docs_specs_reviews-101 | corner-apex stations have one rung; TEA pinned 125.8° | task | 1.3 |
| F-d1_docs_specs_reviews-102 | radius_regulated: one plunge per pass; empirical sweep budget | frozen | Frozen: radius_regulated_toolpath |
| F-d1_docs_specs_reviews-103 | IncompleteSegmentPartitionError lacks sizes; split/merge rule undecided | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-104 | cpp_int → double overflow abort on 139-digit coordinates | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-105 | MATERIAL vs CAP_EXCEEDED not exposed in Python | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-106 | exact segment-site MAT only for one rectangle and one L | parked | generalising the exact segment-site medial axis |
| F-d1_docs_specs_reviews-107 | contour-aware placement not built on live | task | 2.4 |
| F-d1_docs_specs_reviews-108 | live Figure 5 refined path: rejections, publisher start, frozen-stock engagement | frozen | Frozen: live Held reproduction |
| F-d1_docs_specs_reviews-109 | Phase-1 target: planner at 1/10 Held's clocks | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-110 | G2 smoothing does not preserve engagement bound | task | 1.6 |
| F-d1_docs_specs_reviews-111 | precleared bore entry built in adaptive/entry.py, unused | task | 1.2 |
| F-d1_docs_specs_reviews-112 | Held convex-corner limits + concave shared-contact fans unused | task | 2.2, 2.4 |
| F-d1_docs_specs_reviews-113 | exact circular-arc input needs declared construction policy | task | 2.1 |
| F-d1_docs_specs_reviews-114 | rectangle MAT producer emits no clearance profiles or neck evidence | parked | generalising the exact segment-site medial axis |
| F-d1_docs_specs_reviews-115 | planned tri-dexel replay and GPU CDF scoring | superseded | Task 3.3 (named model: exact 2.5D replay) |
| F-d1_docs_specs_reviews-116 | centre-distance growth lemma unsound | task | 4.3 |
| F-d1_docs_specs_reviews-117 | unreachable residual D \ M_r reported separately; gate ignores it | task | 1.3 |
| F-d1_docs_specs_reviews-118 | helical entry + Held sorted arc contour; elliptic returns never built | task | 1.2, 2.1 |
| F-d1_docs_specs_reviews-119 | tool-radius plunge cannot launch a capped lateral move | task | 1.2 |
| F-d1_docs_specs_reviews-120 | max_passes truncation forbidden by spec | task | 1.7 |
| F-d1_docs_specs_reviews-121 | ReachableDomain needs one connected C_r; multi-region raises | task | 1.7 |
| F-d1_docs_specs_reviews-122 | post-Task5 hash-index measurement gate | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-123 | DFS route = draft tree walk; needs inter-route transit | superseded | Task 1.1 (tree walk dropped) |
| F-d1_docs_specs_reviews-124 | terminal coverage, fresh replay, GenerationResult | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-125 | v4 source rejects distinct nonzero x/y radicands | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-126 | AdvanceSegment vs hold Link grammar; gate penalises stay-down links | task | 0.3 |
| F-d1_docs_specs_reviews-127 | exact stay-down retrace superseded; draft adopted retract/plunge | task | 2.3 |
| F-d1_docs_specs_reviews-128 | four adaptive retrace tests expected-red pending retirement | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-129 | swept-prefix fast path research-only until proof note + adversarial tests | task | 4.2 (proof note; Task 2.3 uses it as measured only) |
| F-d1_docs_specs_reviews-130 | forced entries/loops/advances must survive the public result | constraint | Global Constraints: forced decisions |
| F-d1_docs_specs_reviews-131 | one shared ladder; no regulated generator may bisect | constraint | Global Constraints: no bisection |
| F-d1_docs_specs_reviews-132 | tangent tests: comparison-count floor; validate bridge exemptions | task | 1.6 |
| F-d1_docs_specs_reviews-133 | helical motion fails loud until typed 3D depletion contract | task | 1.2 |
| F-d1_docs_specs_reviews-134 | legacy audit: unmeasured = certified | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-135 | Python floors, strict-mypy errors, build identity, CI topology | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-136 | ReplayCertificate value-returning path with positive control | parked | auditor-convergence P0–P4 |
| F-d1_docs_specs_reviews-137 | wave1 Task 6 completion | closed | commit 6b4293ae (subsystem removed) |
| F-d1_docs_specs_reviews-138 | Held pockets supplied; generator quality on them is later | task | 2.1 |
| F-d1_docs_specs_reviews-139 | Fig 5 characterization report and Phase 2 policy design never done | superseded | Task 2.1 (analytic lap is the policy) |
| F-d1_docs_specs_reviews-140 | coverage grid must scale per pocket; general gate unknown | task | 2.1 |
| F-d1_docs_specs_reviews-141 | candidate seam + target-specific postprocessor; no G-code | task | 3.1 |
| F-d1_docs_specs_reviews-142 | Union_of_curve_cycles_2 experiment rejected | closed | commit d40ef841 |
| F-d1_docs_specs_reviews-143 | number-type stages 3–6 | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-144 | dead second double→exact door | parked | number-type stages 3–6, global type-system refactor |
| F-d1_docs_specs_reviews-145 | gated generators emit no lead-in/lead-out | task | 1.6 |
| F-d1_docs_specs_reviews-146 | tangent test exempts bridge junctions; external tangent removes them | task | 1.6 |
| F-d1_docs_specs_reviews-147 | 1e-18 epsilon decides repositioning arc vs coincident | task | 1.6 |
| F-d1_docs_specs_reviews-148 | abs(dot) parallelism check misses direction reversal | task | 1.6 |
| F-d1_docs_specs_reviews-149 | trochoid chain design doc stale | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-150 | C4: depletion model to score a true helical entry | task | 1.2 |
| F-d1_docs_specs_reviews-151 | C5: arc vocabulary for tangent-continuous links | task | 1.6 |
| F-d1_docs_specs_reviews-152 | C1: which layer decides tangency translation variance | parked | translation variance (C1) |
| F-d1_docs_specs_reviews-153 | coherence backlog W2.1/W2.2/S1/S2/R5 open | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-154 | W3: one branch, CI green on main, memo, convergence oracle | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-155 | claimed "V ≥ √3 r proved" absent from live page | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-156 | mat_scale > 1 emits gouging circles | task | 1.7 |
| F-d1_docs_specs_reviews-157 | tea_growth_bound not an upper bound at concave voids | task | 4.3 |
| F-d1_docs_specs_reviews-158 | requires-python ≥ 3.9 but modules need 3.10 | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-159 | Python arc certifier certifies 0/126 circles in generator regime | superseded | Task 4.3 (arc mirror removed) |
| F-d1_docs_specs_reviews-160 | max_tea can report 720° | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-161 | geometry parameters unvalidated (radius, clearance_z) | task | 1.7 |
| F-d1_docs_specs_reviews-162 | link_paths=False emits uncertified through-material traverse | task | 1.7 |
| F-d1_docs_specs_reviews-163 | Epick gouge decisions compare rounded doubles (segment_clear) | task | 1.6 |
| F-d1_docs_specs_reviews-164 | certify_bridges can exit with an uncertified bridge | task | 1.6 |
| F-d1_docs_specs_reviews-165 | hole validation incomplete before straight skeleton | task | 1.7 |
| F-d1_docs_specs_reviews-166 | _unmeasured returns cap_certified=True | task | 4.3 |
| F-d1_docs_specs_reviews-167 | gate rho/chain-ordered/spiral; chain ordering 5→1 entries | task | 1.1 |
| F-d1_docs_specs_reviews-168 | consume declined_regions; split uncut into declined vs missed | task | 1.3 |
| F-d1_docs_specs_reviews-169 | cutting-force, exit-impact, controller-starvation KPIs | task | 3.3 |
| F-d1_docs_specs_reviews-170 | generalise segment-site MAT; middle curve not gentler | parked | generalising the exact segment-site medial axis |
| F-d1_docs_specs_reviews-171 | helical/ramped/pre-bored entry; audit must replay it | task | 1.2 |
| F-d1_docs_specs_reviews-172 | break ~141.8° saturation; find binding constraint | task | 1.5 |
| F-d1_docs_specs_reviews-173 | second gate cap for attribution | closed | commit 4b741af7 (cap-40 attribution) |
| F-d1_docs_specs_reviews-174 | replay_certificate raise; retrace op rejected by replay grammar | parked | `adaptive/` as a generator |
| F-d1_docs_specs_reviews-175 | gate the four prose-only numbers | task | 0.3 |
| F-d1_docs_specs_reviews-176 | island pockets never run on engagement generators | task | 1.7 |
| F-d1_docs_specs_reviews-177 | bridge-certification exit; hole-edge validation | task | 1.6, 1.7 |
| F-d1_docs_specs_reviews-178 | external_tangent hand-rolls perpendicular | task | 1.6 |
| F-d1_docs_specs_reviews-179 | residue 5× unregulated; coverage absent from accept rule | task | 1.4 |
| F-d1_docs_specs_reviews-180 | swept-prefix certifies only π; other caps abort | task | 4.2 |
| F-d1_docs_specs_reviews-181 | non-finite doubles empty Stock2; guard compiled out; 34 GB reserve | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-182 | degenerate and zero-engagement loops, both exactly decidable | task | 1.1 |
| F-d1_docs_specs_reviews-183 | lower cap at bottlenecks; halve air time; G1/G2 | task | 1.6 |
| F-d1_docs_specs_reviews-184 | radius_regulated imports the bisecting advance search | frozen | Frozen: radius_regulated_toolpath |
| F-d1_docs_specs_reviews-185 | forced loops/advances erased at the public result | constraint | Global Constraints: forced decisions |
| F-d1_docs_specs_reviews-186 | Hypothesis pockets are fat stars; no necks or islands | task | 1.7 |
| F-d1_docs_specs_reviews-187 | corpus runs unregulated generator; rho wired into nothing | frozen | Frozen: baseline gate columns |
| F-d1_docs_specs_reviews-188 | arc certifier at 20 stations/turn, 47–52% uncertified | superseded | Task 4.3 (arc mirror removed) |
| F-d1_docs_specs_reviews-189 | swept-prefix certifies from two booleans, no geometric test | task | 4.2 (proof note; Task 2.3 uses it as measured only) |
| F-d1_docs_specs_reviews-190 | segment oracle lacks empty-strata unresolved guard | task | Task 4.5 (oracle robustness and cost) |
| F-d1_docs_specs_reviews-191 | no committed measured artifact; prose numbers stale | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-192 | loops maximised at 171.8° at both caps | task | 1.5 |
| F-d1_docs_specs_reviews-193 | contradictory claim on coverage after dropping degenerate loops | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-194 | unresolved rate on third-party geometry; neck ablation | task | 2.6, 1.7 |
| F-d1_docs_specs_reviews-195 | const_cast in engaged_arcs_zone blocks parallel queries | parked | dead code in engagement_2.cpp |
| F-d1_docs_specs_reviews-196 | CGAL licence gate unresolved | parked | Licensing |
| F-d1_docs_specs_reviews-197 | CONVENTIONS.md points at removed worktree | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-198 | neck and island families not product-gate pockets | task | 1.7 |
| F-d1_docs_specs_reviews-199 | arc-fraction family tessellated; true arcs never exercised | task | 2.1 |
| F-d1_docs_specs_reviews-200 | Pythagorean rotations via binary64 not exactly rational | parked | Docs and governance hygiene |
| F-d1_docs_specs_reviews-201 | external resolution rate and neck ablation never recorded | task | 2.6, 1.7 |
| F-d1_docs_specs_reviews-202 | MATHSM constant-spacing baseline, no MATHSM generator on live | task | 2.4 |
| F-d2_plans-1 | exact MAT and adaptive accept only the hard-coded L | parked | generalising the exact segment-site medial axis |
| F-d2_plans-2 | phase-1 certified covered path never reached | parked | `adaptive/` as a generator |
| F-d2_plans-3 | DFS route order is not a cutter trajectory | superseded | Task 1.1 (tree walk dropped) |
| F-d2_plans-4 | swept-prefix retrace certificate; P2 T4 proof obligation undischarged | task | 4.2 (proof note; Task 2.3 uses it as measured only) |
| F-d2_plans-5 | stay-down reverse retrace exists (L only), uncited | task | 2.3 |
| F-d2_plans-6 | gated generator bisects; pin "no bisection" | constraint | Global Constraints: no bisection |
| F-d2_plans-7 | forced decisions are warnings, not structured outputs | constraint | Global Constraints: forced decisions |
| F-d2_plans-8 | skeleton guides break Held's MA-incidence coverage | task | 2.5 |
| F-d2_plans-9 | generic boundary-normal q-m circle construction exists | task | 2.5 |
| F-d2_plans-10 | corner-approach motion solved Monstera's 111.59° jump | task | 2.4 |
| F-d2_plans-11 | contour-aware contact exists; placement never built | task | 2.4 |
| F-d2_plans-12 | Held T5 continuous Fig 5 and T6 three-curve Fig 6 unreached | task | 2.1, 2.4 |
| F-d2_plans-13 | Held-lane generation 20.6–436 s inherited by a retarget | superseded | Task 2.1 (analytic lap replaces live lane) |
| F-d2_plans-14 | Fig 5 needs ≥ 668×458 coverage grid via Held adapter | task | 2.1 |
| F-d2_plans-15 | live Held path: 1,010 rejections, 2.04 µm remnant | frozen | Frozen: live Held reproduction |
| F-d2_plans-16 | G-code post, parser round-trip, machine simulation never started | task | 3.1, 3.3 |
| F-d2_plans-17 | helical entry unscoreable; arc link vocabulary C5 open | task | 1.2, 1.6 |
| F-d2_plans-18 | exact D \ M_r and reachable predicate for rest material | task | 1.3 |
| F-d2_plans-19 | phase-1 "landed" header false | parked | Docs and governance hygiene |
| F-d2_plans-20 | causal plan header false | parked | Docs and governance hygiene |
| F-d2_plans-21 | headers point at deleted backlog and unrun P2 | parked | Docs and governance hygiene |
| F-d2_plans-22 | wave1-task6 header omits removed artifacts | parked | Docs and governance hygiene |
| F-d2_plans-23 | 164 tasks, zero generation tasks on a gated generator | closed | census d2 effort table |
| F-d2_plans-24 | gate tangent criterion has no structural bridge exemption | closed | benchmarks/quality_observations.py:913-923 |
| F-d2_plans-25 | truthful three-way audit never wired into gate | parked | auditor-convergence P0–P4 |
| F-d2_plans-26 | external resolution rate and neck ablation never reported | task | 2.6, 1.7 |
| F-d2_plans-27 | no Phase 2 generator-policy design or stage-accounted performance | superseded | Architecture (this plan is the policy) |
| F-d2_plans-28 | exact zero-guide classification is L-only | parked | generalising the exact segment-site medial axis |
| d1-T1 | no entry strategy the gate or auditor can accept | task | 1.2 |
| d1-T2 | engagement limited not regulated; saturation; bisection | task | 1.5 |
| d1-T3 | chain order greedy, per-chain plunge; ordered generator deleted | task | 1.1 |
| d1-T4 | convex corners and narrow regions have no motion | task | 1.3 |
| d1-T5 | residue is not a corner problem | task | 1.4 |
| d1-T6 | bridges unregulated; tip-ward walks slot | task | 1.1 |
| d1-T7 | G1 breaks from bridge tilt; tangent test exempts bridges | task | 1.6 |
| d1-T8 | air cutting: full loops recut ~90% | task | 1.6 |
| d1-T9 | no lower bound, no constant-engagement target | task | 1.5 |
| d1-T10 | mat_scale gouges, hole validation, translation variance | task | 1.7 |
| d1-U1 | rho floor, tangent bridges, ladder, spiral entry unused | task | 1.1 |
| d1-U2 | precleared bore, M_r, D \ M_r, reachable predicate unused | task | 1.2, 1.3 |
| d1-U3 | exact refutation ladder; can_remove_circle unused | task | 0.2, 4.1 |
| d1-U4 | swept-prefix, certify_*_tea, zero-guide advancing segments | task | 2.3, 4.3 |
| d1-U5 | candidate order furthest progress then largest radius | task | 1.1 |
| d1-U6 | Held lane: qm circles, boundary-order spacing, corner fans | task | 2.4 |
| d1-U7 | local zone query, local annulus depletion, 6-curve quad | task | 0.2 |
| d1-U8 | offset_polygon and weighted straight skeleton unused | superseded | Task 1.3 (exact cutter-centre boundary) |
| d1-U9 | MATHSM constant-spacing sweep | task | 2.4 |
| d1-U10 | quality reducers with operation attribution | task | 0.3 |
| d1-U11 | deleted chain-ordering stay-down generator recoverable | task | 1.1 |
| d1-P1 | draft re-derives the deleted ordered generator | task | 1.1 |
| d1-P2 | bridges and links not cap-checked | task | 1.1 |
| d1-P3 | D1: step stays gated; auditor rejects ramps; bore cheaper | task | 1.2 |
| d1-P4 | D2: Held regularisation, D \ M_r, Task 31 contract, stock-aware run-ins | task | 1.3 |
| d1-P5 | coverage unaddressed | task | 1.4 |
| d1-P6 | productive-only emission cost unmeasured | task | 0.2 |
| d1-P7 | plan right not to switch guide | superseded | Task 2.5 (guide switch after measurement) |
| d1-P8 | retract-and-plunge cost; advance vs hold distinction | task | 0.3 |
| d1-P9 | not covered: saturation, W ≤ 4r, leads, arc links, air, bound, guide hygiene | task | 1.5, 1.6, 1.7 |
| d2-T1 | no task in 22 plans targeted the gated generators | superseded | Spine 4 + Task 0.3 ratchet |
| d2-T2 | exact MAT and adaptive accept only one L | parked | `adaptive/` as a generator; generalising the exact segment-site medial axis |
| d2-T3 | gated generator bisects under false monotonicity | constraint | Global Constraints: no bisection |
| d2-T4 | depth-first route order needs transit; retrace exists | task | 2.3 |
| d2-T5 | skeleton breaks Held coverage; boundary-normal construction exists | task | 2.5 |
| d2-T6 | Task 0b: coarse grid and Held residuals | task | 2.1 |
| d2-T7 | machine-readiness planned twice, never started | task | 3.1, 3.3 |
| d2-T8 | exact rest-material region exists for D2(a) | task | 1.3 |
| d2-U1 | native boundary-normal q-m circle | task | 2.5 |
| d2-U2 | contour-aware contact HeldDiskContour2 | task | 2.4 |
| d2-U3 | concave-corner approach/rotation/departure | task | 2.4 |
| d2-U4 | boundary-order spacing | task | 2.4 |
| d2-U5 | swept-prefix exact theorem | task | 2.3 |
| d2-U6 | exact cut-depth retrace | task | 2.3 |
| d2-U7 | exact zero-guide MAT classification, L only | parked | generalising the exact segment-site medial axis |
| d2-U8 | exact reachable domain and predicate | task | 1.3 |
| d2-U9 | exact partial-arc depletion subtract_exact_arc | task | 1.6 |
| d2-U10 | exact gouge containment for segments and circles | task | 1.2, 1.6 |
| d2-U11 | MATHSM baseline | task | 2.4 |
| d2-U12 | transactional fork/certify/deplete/commit evaluation | superseded | Task 0.2 (cheap exact redundancy tests) |
| d2-P1 | tree walk assumes skeleton sound; transit cost unnamed | task | 2.5, 2.3 |
| d2-P2 | exact MAT / adaptive as generator contradicted | parked | `adaptive/` as a generator; generalising the exact segment-site medial axis |
| d2-P3 | Task 0b: coverage grid adapter and time budget | task | 2.1 |
| d2-P4 | D1: depletion rejects ramped moves | task | 1.2 |
| d2-P5 | backlog C5 arc links missed | task | 1.6 |
| d2-P6 | D2(a) should use exact D \ M_r | task | 1.3 |
| d2-P7 | auditor P3 T1/T2 must survive as constraints | constraint | Global Constraints: no bisection; forced decisions |
| d2-P8 | plan status headers false or point at deleted files | parked | Docs and governance hygiene |
| d2-P9 | no G-code or postprocessor stage | task | 3.1 |
| g1-T1 | engagement model static, direction-blind, in gate and generators | task | 0.1 |
| g1-T2 | guide consumed narrow end first, plunge per chain, greedy | task | 1.1 |
| g1-T3 | Python bridges certified for neither walls nor engagement | task | 1.6 |
| g1-T4 | advance-only lineage still bisects | constraint | Global Constraints: no bisection |
| g1-T5 | entry placement causes tangent breaks; slope over neighbours not chosen advance | task | 1.6 |
| g1-T6 | generators deplete with different primitives from gate replay | task | 0.2 |
| g1-T7 | skeleton loses the L hub's true widest point | task | 2.5 |
| g1-T8 | convex corners below ρ > r floor: exact blocks unused | task | 1.3 |
| g1-U1 | C++ certified bridges, leads, flat-link gouge check | task | 1.6 |
| g1-U2 | rho lineage: floor, radius step, tangent entries, DeclinedRegion | task | 1.1 |
| g1-U3 | spiral entry centre plunge + ramp | superseded | Task 1.2 (R3 helical bore) |
| g1-U4 | radial _loop_reaches_material exact bite test | task | 0.2 |
| g1-U5 | Stock2 exact ops, can_remove_circle, HeldDiskContour2, circle sweep | task | 0.2, 0.1 |
| g1-U6 | continuous exact TEA certificates and swept-prefix audit | task | 4.1, 4.2 |
| g1-U7 | ReachableDomain2, remaining_material, ccw_transition | task | 1.3, 1.4 |
| g1-U8 | exact containment incl. circle within entry disk | task | 1.2, 1.6 |
| g1-U9 | contour-aware circle placement BoundaryNormalCircle2 | task | 2.5 |
| g1-U10 | exact segment-site MAT, fixtures only | parked | generalising the exact segment-site medial axis |
| g1-U11 | straight-skeleton offsets for finishing and centre boundary | superseded | Task 1.3 (exact cutter-centre boundary) |
| g1-U12 | typed motion vocabulary unused by gated generators | constraint | Global Constraints: engineering bar |
| g1-P1 | draft's L is not the gate's L | superseded | Stage 1 (gate pockets used directly) |
| g1-P2 | root chosen among skeleton vertices, not over the guide | task | 1.2 |
| g1-P3 | away from root is not wide to narrow | superseded | Task 1.1 (tree walk dropped) |
| g1-P4 | DFS gives a retract and plunge per leaf | task | 1.1 |
| g1-P5 | productive = not exactly_equals; air loops survive; cost | task | 0.2, 1.6 |
| g1-P6 | _next_step fallback breaks rho's rules | superseded | Task 1.1 (draft code dropped) |
| g1-P7 | bridges uncertified against walls | task | 1.6 |
| g1-P8 | D1 premise relative to the static engagement model | task | 0.1 |
| g1-P9 | draft defects: bare ValueError, path_index stamping, global subtract | superseded | Task 1.1 (draft code dropped) |
| g1-P10 | B2a buildable from exact centre-boundary kernels | task | 1.3 |
| g1-P11 | exact MAT not a drop-in | parked | generalising the exact segment-site medial axis |
| g1-P12 | draft fallback plunge and signatures checked sound | superseded | Task 1.1 (draft code dropped) |
| g2-T1 | adaptive placement cannot cover two-sided regions | parked | `adaptive/` as a generator |
| g2-T2 | adaptive dies at first non-incident route switch | parked | `adaptive/` as a generator |
| g2-T3 | zero-guide advances certifiable only at 180° | parked | `adaptive/` as a generator |
| g2-T4 | adaptive: hold LINK + full circle, tangent breaks, air | parked | `adaptive/` as a generator |
| g2-T5 | adaptive cost O(n²) with seconds constants | parked | `adaptive/` as a generator |
| g2-T6 | native oracle fragile; one unresolved event aborts generation | parked | `adaptive/` as a generator |
| g2-T7 | adaptive pocket and entry scope narrow | parked | `adaptive/` as a generator |
| g2-U1 | continuous exact TEA certification | task | 4.1 |
| g2-U2 | exact cheap refutation ladder | task | 4.1 |
| g2-U3 | exact segment-site MAT | parked | generalising the exact segment-site medial axis |
| g2-U4 | neck inventory with per-pass effective caps | parked | generalising the exact segment-site medial axis |
| g2-U5 | PreclearedEntry + QualifiedBore | task | 1.2 |
| g2-U6 | ReachableDomain exact regions | task | 1.3 |
| g2-U7 | engagement-certified stay-down LINK at depth | task | 1.1 |
| g2-U8 | exact gouge containment | task | 1.2, 1.6 |
| g2-U9 | swept-prefix theorem at 180° | task | 2.3 |
| g2-P1 | needs ruling MAT vs skeleton at reflex corners and root | task | 2.5 |
| g2-P2 | NeckInventory exists; draft reinvents neck detection | parked | generalising the exact segment-site medial axis |
| g2-P3 | D1 should reuse PreclearedEntry contract | task | 1.2 |
| g2-P4 | D2(a) must state which reachability | task | 1.3 |
| g2-P5 | stay-down link can be certified with existing code | task | 2.3 |
| g2-P6 | cap compliance sampled; refutation ladder and final certificate ignored | task | 4.1 |
| g2-P7 | do not adopt adaptive/generator.py | parked | `adaptive/` as a generator |
| g3-T1 | exact MAT runs on one hard-coded pocket | parked | generalising the exact segment-site medial axis |
| g3-T2 | cell machinery general, traversal not | parked | generalising the exact segment-site medial axis |
| g3-T3 | parabolas, clearance polynomials, necks the skeleton lacks | parked | generalising the exact segment-site medial axis |
| g3-T4 | rectangles: MAT equals skeleton; only L differs | task | 2.5 |
| g3-T5 | MAT sampling limited to rational charts | parked | generalising the exact segment-site medial axis |
| g3-T6 | CORE identity decisions untested on generic input | parked | generalising the exact segment-site medial axis |
| g3-U1 | parabolic arcs at reflex vertices | task | 2.5 |
| g3-U2 | exact clearance function per guide edge | parked | generalising the exact segment-site medial axis |
| g3-U3 | exact clearance ≥ r clipping | parked | generalising the exact segment-site medial axis |
| g3-U4 | neck detection and separating cuts | parked | generalising the exact segment-site medial axis |
| g3-U5 | neck-width classification against exact thresholds | parked | generalising the exact segment-site medial axis |
| g3-U6 | zero-guide slot detection | parked | generalising the exact segment-site medial axis |
| g3-U7 | holes carried through catalog and clipping | parked | generalising the exact segment-site medial axis |
| g3-P1 | root-at-node lemma should be stated | superseded | Task 1.2 (root over the guide) |
| g3-P2 | neck semantics should be specified, not counted after the fact | task | 1.7 |
| g3-P3 | islands: exact MAT rejects holes too | constraint | Review Focus 1 |
| g3-P4 | say why exact MAT is not used | parked | generalising the exact segment-site medial axis |
| g3-P5 | place stations by exact clearance polynomial | parked | generalising the exact segment-site medial axis |
| g4-T1 | exact removal test on an inexact model | task | 0.2 |
| g4-T2 | plunge at loop's own entry never exactly clear | task | 0.2 |
| g4-T3 | disk-chain lines blow up arrangement and decision cost | task | 0.2 |
| g4-T4 | TEA-growth lemma false for concave boundaries | task | 4.3 |
| g4-T5 | engagement capped per run, not total contact | task | 4.4 |
| g4-T6 | generators sample; continuous certification unusable near cap | task | 4.1 |
| g4-T7 | dead code in the deciding module | parked | dead code in engagement_2.cpp |
| g4-T8 | NRVO lifetime hazard on by-value Gps returns | parked | NRVO lifetime hazard |
| g4-U1 | exact reachable material and unreachable residual | task | 1.3 |
| g4-U2 | exact centre-domain cycle and CCW transitions | task | 1.3 |
| g4-U3 | exact boundary-normal medial disk | task | 2.5 |
| g4-U4 | gouge predicates incl. loop-inside-entry-disk | task | 1.2 |
| g4-U5 | exact remaining material over CORE kernel | task | 1.4 |
| g4-U6 | exact circle-redundancy can_remove_circle | task | 1.6 |
| g4-U7 | whole-orbit engagement bound HeldDiskContour2 | task | 0.1 |
| g4-U8 | continuous TEA certificate; full-circle authority | task | 4.1 |
| g4-U9 | exact rational points on a circle | task | 0.2 |
| g4-U10 | local depletion; gate replay does not use it | task | 0.2 |
| g4-P1 | productive-only emission rests on sliver-polluted test | task | 0.2 |
| g4-P2 | _enter docstring contradicted | task | 0.2 |
| g4-P3 | exact predicate over under-covering model insufficient | constraint | Global Constraints: exact decisions on exact models |
| g4-P4 | affordability needs quad or exact capsules; cheaper test | task | 0.2 |
| g4-P5 | B2/D2 ignore exact residual and wall-offset cycle | task | 1.3 |
| g4-P6 | D1 ignores entry-disk coverage and containment | task | 1.2 |
| g4-P7 | multi-component centre domain throws PocketNotMachinableError | constraint | Review Focus 1 |
| g5-T1 | frozen-stock oracle refutes cutting segments; forward-semicircle mask | task | 4.2 |
| g5-T2 | no API yields a certified maximal advance | task | 4.1 |
| g5-T3 | full-circle audit 100–150× slower than probes | task | 4.1 |
| g5-T4 | incompleteness needs defined UNRESOLVED action | task | Task 4.5 (oracle robustness and cost) |
| g5-T5 | circle oracle has no cheap refuter | task | 4.1 |
| g5-T6 | exhaustive pair path quadratic on all-rational stocks | task | Task 4.5 (oracle robustness and cost) |
| g5-U1 | audit_full_circle_tea_event_exact drop-in | task | 4.1 |
| g5-U2 | swept-prefix segment certification at π | task | 2.3 |
| g5-U3 | exact one-station segment refuter | task | 4.1 |
| g5-U4 | refute-then-certify pattern | task | 4.1 |
| g5-U5 | sign_mixed_radical unused by continuous_tea_2 | task | Task 4.5 (oracle robustness and cost) |
| g5-P1 | cap decision stays sampled | task | 4.1 |
| g5-P2 | certification belongs later, after a measurement | task | 4.1 |
| g5-P3 | bridges and run-ins have no exact certificate below π | task | 4.2 |
| g5-P4 | helix has no exact certificate; new work | task | 1.2 |
| g5-P5 | draft omits adaptive/ and why not reused | parked | `adaptive/` as a generator |
| g6-T1 | live reproduction 2.8× longer, 7.6× circles vs Held | task | 2.1 |
| g6-T2 | raw guide is the product skeleton, not publisher-derived | superseded | Task 2.1 (analytic lap) |
| g6-T3 | native publisher-free Held circles run on gate pockets | task | 2.2 |
| g6-T4 | degenerate-loop rule excludes Held placement | task | 0.3 |
| g6-T5 | spacing order differs from emission order | task | 2.4 |
| g6-T6 | contour-aware spacing halves length but loses coverage | task | 2.4 |
| g6-T7 | no entry anywhere in the Held lane | task | 1.2 |
| g6-T8 | convex corners handled by shrinking ρ to 0 | task | 2.2 |
| g6-T9 | predecessor-only 79.99° is not whole-history; connectors unaudited | task | 2.1 |
| g6-T10 | containment rejections are tangency rounding; gate will see noise | task | 2.1 |
| g6-U1 | boundary-order single-lap traversal | task | 2.1 |
| g6-U2 | exact native Held circles from offset contact | task | 2.2 |
| g6-U3 | Held standard-model engagement | task | 2.1 |
| g6-U4 | whole-history filled-disk contour bound | task | 0.1 |
| g6-U5 | coverage-preserving deletion can_remove_circle | task | 1.6 |
| g6-U6 | exact residual stock remaining_material | task | 1.4 |
| g6-U7 | reflex-corner phase partition (corner approaches) | task | 2.4 |
| g6-P1 | publisher-derived dependency claim contradicted | superseded | Task 2.1 (analytic lap) |
| g6-P2 | score native Held lap on gate pockets | task | 2.2 |
| g6-P3 | decision rule biased by degenerate rule and gouge noise | task | 0.3 |
| g6-P4 | 79.99° quoted without rejections or whole-history bound | superseded | Task 0.1 (in-process gate) |
| g6-P5 | Held lap needs no repositioning | task | 2.1 |
| g6-P6 | D2 ignores is_stationary corner witness | task | 1.3 |
| g6-P7 | Held should be a candidate generator topology | task | 2.1 |
| g7-T1 | publisher-free analytic Held generator exists off live | task | 2.1 |
| g7-T2 | contour-aware placement exists off live | task | 2.4 |
| g7-T3 | generic segment-site Voronoi guide off live | task | 2.5 |
| g7-T4 | in-motion swept-prefix engagement off live | task | 0.1 |
| g7-T5 | live audit certifies on the falsified lemma | task | 4.3 |
| g7-T6 | toolpath.cpp contract fixes missing; b0094cb6 conflict | task | 1.7 |
| g7-T7 | two stock-clone traits fixes; keep live's | closed | commits 57482731, 8d2b76cc |
| g7-T8 | publisher Figure 5 turns miscounted 265 vs 278 | task | 2.1 |
| g7-T9 | Figure 5 reference G1 seam repair off live | task | 2.1 |
| g7-T10 | jf diagnoses: saturation, entry exceedances, rho_max, middle curve | task | 1.5, 1.2, 2.3 |
| g7-U1 | analytic Held standard path | task | 2.1 |
| g7-U2 | contour-aware placement | task | 2.4 |
| g7-U3 | placement-policy traversal, MATHSM on one locus | task | 2.4 |
| g7-U4 | Figure 6 sweep | task | 2.4 |
| g7-U5 | generic segment-site Voronoi guide | task | 2.5 |
| g7-U6 | swept-prefix in-motion engagement | task | 0.1 |
| g7-U7 | analytic offset-contour construction | task | 2.1 |
| g7-U8 | adaptive arc certification on swept-annulus bound | task | 4.3 |
| g7-U9 | toolpath.cpp generator contracts | task | 1.7 |
| g7-P1 | plan:38 publisher dependency contradicted | task | 2.1 |
| g7-P2 | only skeleton-chain ordering considered | task | 2.1 |
| g7-P3 | rebase conflicts underestimated | task | 0.1, 2.1 |
| g7-P4 | jf contract fixes pending; b0094cb6 alone raises | task | 1.7 |
| g7-P5 | exactly_equals cannot score engagement during a link | task | 0.1 |
| g7-P6 | entry modelled as pre-cleared disk; no helix anywhere | task | 1.2 |
| v-T1 | gate punishes stay-down, ignores retracts | task | 0.3 |
| v-T2 | slotting criterion is not slotting | task | 0.3 |
| v-T3 | gate cannot represent a helical or ramped entry | task | 1.2 |
| v-T4 | green tests pin today's bad behaviour | task | 0.3 |
| v-T5 | uncut over tool-reachable; operation vs job | task | 1.3 |
| v-T6 | degenerate loops defined only for full circles | task | 0.3 |
| v-T7 | short path and few retracts not gated; cycle time unused | task | 3.2 |
| v-T8 | air-cutting weakly enforced; recut and immersion reported only | task | 1.5, 1.6 |
| v-T9 | plunge XY, plunges into material, junction step unchecked | task | 0.3 |
| v-T10 | gate pockets too narrow | task | 1.3, 1.7, 2.1 |
| v-T11 | CI never measures generation | task | 0.3 |
| v-U1 | rho_regulated tested, not gated | task | 1.1 |
| v-U2 | spiral entry not gated | superseded | Task 1.2 (R3 helical bore) |
| v-U3 | certified adaptive stack has no bridge to ToolpathResult | parked | `adaptive/` as a generator |
| v-U4 | Held native kernels for contour-aware placement | task | 2.5, 2.4 |
| v-U5 | quad capsule subtraction and local depletion | task | 0.2 |
| v-U6 | curvature-limited feed time objective | task | 3.2 |
| v-U7 | pocket families and external loader | task | 1.7 |
| v-U8 | MATHSM baseline never compared by gate | task | 2.4 |
| v-P1 | reposition above stock optimises a gate distortion | task | 0.3 |
| v-P2 | 0 redundant is not no air-cutting loops | task | 1.6 |
| v-P3 | D1 must cover classification, pins, partition, z-aware area | task | 1.2 |
| v-P4 | straight motions should use the full cap | task | 0.3 |
| v-P5 | adding a generator must update attribution, pins, red_manifest | task | 0.3 |
| v-P6 | score on island, dumbbell, stadium, disk | task | 1.7, 2.1, 2.5 |
| v-P7 | gate-side plunge count and plunge-XY check | task | 0.3 |
