# Census slice g2_adaptive: `src/compas_cgal/adaptive/`

## Slice summary

The certified adaptive generator does not produce a complete toolpath on any pocket today, and it cannot produce one on the gate pockets in its current form. Its end-to-end output is a launch-rooted prefix (`GenerationContinuation`) on two hand-built fixtures: a radius-1 L of 2-wide slots at a 180° cap (Task 13F), and a scaled L "branch" fixture. Both runs stop with an exception before the MAT traversal finishes (`docs/continuous_engagement.md:917-930`; `tests/adaptive/test_generator.py:491-634`; both tests are expected-red in `docs/red_manifest-adaptive.json:4-15`). The coverage seal (`GenerationResult`) was never built (`generator.py:1362-1366`), and `replay_certificate` can only raise (`replay.py:1276`). What the package does provide are several exact *components* that are worth more than the generator they sit in.

### Top findings, ranked by impact on toolpath quality

1. **The placement cannot cover a two-sided region. Coverage is not an objective.** (Code reading. The geometric consequence is unverified by measurement.) Every candidate is a "one-sided MATHSM" circle tangent to the r-offset of one generator site, with radius ρ ≤ (d−r)/2 (`candidates.py:124-164`, `:1160-1193`). Its sweep reaches distance d+r from that site, so in a channel of half-width d it leaves a strip of width d−r on the far side of the MAT. The route visits each MAT edge exactly once (`traversal_graph.py:528-583`). Each step accepts the first candidate in "furthest progress, largest radius, identity hash" order (`policy.py:157-180`), whichever site that is (`generator.py:1221-1268`). Nothing returns to the other side, and the continuation "deliberately does not assert exact residual emptiness" (`generator.py:1362-1366`). On `rect_12x8`/`rect_20x12` this leaves most of one half of each strip.
2. **Generation dies at the first non-incident route switch on any ordinary pocket.** The traversal is a DFS over MAT edges (`traversal_graph.py:557-580`). Moving to the next edge when the last edge's exit node is not the next edge's entry node requires a `RouteRetraceCommit`. The only admitted retrace source is a final *zero-guide* advancing segment, at no-neck and full cap (`generator.py:513-557`, `operation.py:946-953`). Anything else raises `UnsupportedRouteRetraceError` (`generator.py:515`, `:535`, `:555`, `:1381-1384`). Zero-guide edges exist only where the clearance is identically r, i.e. a slot exactly 2r wide (`docs/segment_site_mat.md:299`). A rectangle's MAT forces a backtrack after its first corner edge, so `rect_12x8` and `rect_20x12` abort there (inference from code, unmeasured).
3. **Zero-guide advances are certifiable only at a 180° cap.** The swept-prefix theorem proves engagement ≤ π from a clear start disk. "The native theorem currently admits only pi" (`motion_certificate.py:559-561`), and tighter caps return unresolved (`docs/segment_site_mat.md:1786`). At the gate's 120° and 40° caps the whole zero-guide path, including the retrace that depends on it, raises `UnresolvedMotionEventError`. The only fixture that gets past route 1 uses `EngagementCap.build(math.pi)` (`tests/adaptive/task13f_fixture.py:110`).
4. **Every step is a straight hold LINK at cut depth followed by a full 360° circle** (`transaction.py:755-840`). The link runs from the last phase point to the next circle's phase point, so the tangent breaks at every link/circle junction unless the geometry happens to line up. Nothing aligns them: the phase is a lattice of `phase_count` rotations of the boundary-contact vector (`candidates.py:167-184`). Every full circle also re-traverses its already-cleared arc (air cutting). Against the world-class bar (G1 continuity, no air loops) this is *weak*, and on the same axes it is worse than `rho_regulated` (0 tangent breaks on 2 of 3 gate pockets).
5. **Performance rules it out as a generator at gate scale (measured by others, not by me).** One exact segment audit costs 5.1–6.8 s and one circle audit 0.002–0.54 s (`motion_refutation.py:15-16`; `docs/continuous_engagement.md:722-726`). A 16-candidate family takes 35 s with the refutation probe, 189 s without (`docs/continuous_engagement.md:780-785`). Route one of Task 13F takes 17.7 s (`:787-790`). On top of that, every accepted transaction is recomputed in full at commit (`transaction.py:698-705`, `bootstrap.py:504-511`, `retrace_transaction.py:326-334`). Every `GenerationState` revalidates the whole chronology and recomputes the stock-boundary digest (`generation_state.py:165-234`). Every depletion witness carries all predecessor digests (`stock_area.py:402`), and every coverage sweep clones the native ledger and recomputes the residual components (`coverage.py:273-336`). Total cost is O(n²) in the number of operations, with seconds-scale constants. A gate pocket at hundreds of circles is out of reach (extrapolation, unverified).
6. **The native oracle is fragile.** On the Task 13F route-2 family, the full partition "terminates the process: `boost float_next<double>: Argument must be finite, but got inf`" when the refutation probe is disabled. `IncompleteSegmentPartitionError` occurs on 15 of 125 measured audits (`docs/continuous_engagement.md:868-875`, `:917-923`). This is incompleteness, not unsoundness: the 31 refutations were cross-checked with 0 disagreements (`:877-893`). But one unresolved event anywhere aborts the whole generation (`generator.py:1000-1035`; `motion_certificate.py:516-527`).
7. **Pocket and entry scope is narrow.**
   - The reachable domain must be one connected component (`reachable_domain.py:484-485`), so a dumbbell whose neck is narrower than 2r is rejected outright.
   - The program grammar has a single approach and plunge and no retract (`generation_state.py:343-378`; `replay.py:237-265`).
   - The entry is an externally supplied `PreclearedEntry`, i.e. a disk of radius > r with `QualifiedBore` process evidence (`entry.py:110-202`). The caller also supplies the entry MAT edge and node (`traversal.py:624-669`).

### Generation capabilities here that the gated generators do NOT use

- **Continuous exact TEA certification** of a segment or full circle over its whole motion (`motion_certificate.py:452-540`). This is a certificate, not sampled probes. The gated generators "deliberately do not import" it (`engagement_toolpath.py:28-31`; `docs/engagement_controlled_toolpath.md:277`).
- **Exact cheap refutation**: a 9-station dyadic ladder at 3.2 ms per station. A hit is a complete proof that a motion exceeds the cap (`motion_refutation.py:164-206`, `:350-406`). It is usable as a sound *rejection* filter inside any generator's candidate search, for example the plan's ladder scan, at negligible cost.
- **Exact segment-site medial axis** (`medial_axis.py:537-622`, native `_medial_axis_2`). It has parabolic arcs at reflex corners, exact clearance per sample, endpoints clipped at the tool-fit clearance root (`medial_axis.py:64-66`, `:290-303`), stable identities, and a neck inventory. The plan uses the straight skeleton instead.
- **Neck inventory with causal passage states and per-pass effective caps** (`neck.py:255-363`, `traversal.py:884-1024`, `policy.py:239-308`). It detects narrow passages exactly and can lower the cap per pass.
- **`PreclearedEntry` + `QualifiedBore`.** The entry is a declared process (bore) over the full cut interval, removed as an exact disk and excluded from the lateral cap. The first circle must lie inside the bore disk (`entry.py:44-263`, `stock_area.py:331-353`, `containment.py:483-534`). This is exactly the shape of the plan's D1 helical-bore ruling.
- **`ReachableDomain`** with exact `reachable_material` and `unreachable_residual` regions (`reachable_domain.py:536-622`). This is the exact region that D2(a) "measure uncut only over what is reachable" needs. Its reachability is tool-disk reachability, not ρ>r trochoid reachability.
- **An engagement-certified stay-down LINK at cut depth** (`transaction.py:1066-1125`): exact containment, a refutation probe, then a full certificate. This is the motion the plan defers to D3, with a proof attached.
- **Exact gouge containment** for segment capsules, full-circle annuli and entry disks (`containment.py:328-534`).
- **The swept-prefix theorem** for advancing from a clear start (`motion_certificate.py:542-633`). It is only useful at a 180° cap.

### Where this slice shows the draft plan is wrong or incomplete

- **Plan "Design (Stage A)": the guide.** The plan builds the tree on the straight skeleton and roots it at "the interior skeleton vertex of largest exact clearance". The skeleton's offset time is not clearance near reflex vertices: the true distance to a reflex point is smaller, and the MAT there is parabolic. `L_shape` has one reflex vertex. The repository already has an exact segment-site MAT that yields exact clearance per node and sample, a tool-fit clip and neck evidence (`medial_axis.py`). The plan neither mentions it nor rules it out. *Incomplete*: it needs a ruling on MAT versus skeleton, at least for reflex corners and for choosing the root.
- **Plan Review Focus 2 (dumbbell neck).** `NeckInventory` (`neck.py:288-363`) classifies necks exactly and yields their sides and loci. The plan reinvents neck detection implicitly through `virgin_plunges`. The plan's handling, a second entry recorded in `virgin_plunges`, is more general than the adaptive package, which rejects multi-component reachable domains (`reachable_domain.py:484-485`).
- **Stage B D1 (entry).** The recommended `OperationType.HELIX`, removing the bore disk and excluded from the planar cap, is already typed here as `PreclearedEntry`/`QualifiedBore` with an exact disk depletion and an exact "first circle inside the bore" certificate. D1 should reuse that contract rather than invent a parallel one.
- **Stage B D2(a) (rest material).** Exact `unreachable_residual` and `reachable_material` regions already exist (`reachable_domain.py:612-622`). D2(a) should state whether its "reachable" means tool-disk reachability, which these regions give, or ρ>r trochoid reachability, which is a stricter region nobody has built.
- **Stage C D3 (stay-down links).** The adaptive LINK is a stay-down move through cleared stock, with exact engagement and containment proofs (`transaction.py:1066-1125`). The plan frames D3 as a gate amendment only. It should also say that a stay-down link can be *certified* (below the cap and gouge-free) with existing code. The only reason not to use it in generation is cost, at 5–7 s per link audit.
- **Plan "Global Constraints": "no sampled-engagement shortcuts" applies only to removal, plunge and root.** Cap compliance is inherited from the `rho_regulated` station rules (sampled). The plan never considers the exact refutation ladder, which is cheap and sound for rejection, or a final continuous certificate of the emitted motions. That is a gap against the project's stated thesis.
- **Where the plan is right and adaptive is worse.**
  - Ordering: the tree walk from the widest vertex beats adaptive's DFS in SHA-256 branch order (`policy.py:364-369`, `medial_axis.py:145-147`).
  - Coverage: centred loops with ρ = clearance − r reach both sides of the guide; adaptive's one-sided circles do not (finding 1).
  - Repositioning: retract and plunge versus no repositioning at all (finding 2).
  - Speed: milliseconds versus seconds per motion.
  - The plan should not adopt `adaptive/generator.py` as the generator.

### (a) What it generates end to end today

- **Entry point:** `generate_exact_adaptive_continuation` (`generator.py:1336-1425`).
- **Inputs:** a seeded `MatTraversalState`, an already-evaluated launch transaction and two evaluators.
- **Output:** a `GenerationContinuation`, i.e. operations `approach, plunge, circle, (link, circle | advance)*, retrace?`, plus commits. It is never a coverage-sealed path (`:1362-1366`).
- **Consumers:** no product, benchmark or gate code calls it. Its only callers are `tests/adaptive/test_generator.py:476`, `:538` and fixtures.
- **Branch fixture** (scaled L, r=1/4, cap π, `tests/adaptive/test_acceptance.py:80-95`, `:388-389`): it used to stop at an unresolved exact event. Since the refutation probe was added it skips that candidate and "accepts a later feasible candidate". Its test is now red and the full outcome is unrecorded (`docs/continuous_engagement.md:925-937`).
- **Task 13F** (6×6 L, legs 2 wide, r=1, cap π, entry bore 17/16, `task13f_fixture.py:62-140`):
  - Its runs: a launch circle (guide radius 1/32), route 0 with one circle, route 1 with one zero-guide advancing segment, a retrace, and route 2 (`task13f_fixture.py:325-400`).
  - The pinned end state is `NoFeasibleCandidateError` at route 2 (`test_generator.py:546-551`).
  - The docs report that the current tree gives `attempts=56; cap=10; gouge=46` (`continuous_engagement.md:917-923`). That run is the same `NoFeasibleCandidateError` exhaustion, with a different count of rejections than pinned.
  - I could not tell from code or docs whether that measurement predates the retrace wiring in `generator.py:1392-1418` (unverified). I did not run it: importing triggers an editable native rebuild in a shared worktree with uncommitted C++ edits from another session.
- **Result:** partial certified fragments of under a dozen operations on toy fixtures. No complete path exists on any pocket.
- **The stops, named:**
  - `UnresolvedMotionEventError` (native partition incomplete, `motion_certificate.py:516-527`, `:608-623`).
  - `NoFeasibleCandidateError` / `EngagementCapInfeasibleError` / `NeckTooTightError` (`generator.py:1023-1035`).
  - `UnsupportedRouteRetraceError` (`generator.py:515`, `:535`, `:555`, `:1381-1384`).
  - A native process abort (`continuous_engagement.md:917-919`).

### (b) Guide, placement, ordering, linking, entry vs the bar and the gated generators

| aspect | adaptive | world-class bar | gated generators / tree-walk plan |
| --- | --- | --- | --- |
| guide | exact segment-site MAT, parabolic arcs, clipped at clearance r, necks (`medial_axis.py`, `neck.py`) | MAT-based (MATHSM, Held) | straight skeleton, edges "narrower end first" (`toolpath.cpp:656-661`); plan: skeleton tree |
| placement | one-sided circle tangent to the r-offset, ρ ≤ (d−r)/2, finite spatial × radius × phase lattice, greedy furthest/largest feasible (`candidates.py:956-1200`, `policy.py:157-180`) | contour-aware, constant engagement, full coverage | stations on the skeleton with ρ from clearance and a ladder scan (sampled) |
| engagement | exact continuous certificate per motion, plus an exact refutation probe | constant ≤ cap everywhere | sampled; 5–9 cap exceedances per pocket (plan table) |
| ordering | DFS over MAT edges, branches in SHA-256 order, each edge once (`traversal_graph.py:427-583`) | widest-first, short path | narrowest-first chains (bad); plan: tree from the widest vertex (better) |
| linking | stay-down hold LINK at cut depth, certified; cross-route only by reversing a zero-guide segment | stay down through cleared stock, G1 | plunge per chain; plan: retract → LINK at clearance → plunge |
| entry | external qualified bore disk > r plus a first circle inside it | single helix or ramp | circle into virgin stock (360°); plan D1: helix |
| continuity | G0 at each link/circle junction; a full circle per step | G1/G2 | rho: 0–2 tangent breaks |

### (c) Why it is not in the product gate, and what it would take

It is not in the gate because:

- it emits `CanonicalOperation`s, not a `ToolpathResult`, and no adapter exists;
- it needs a caller-supplied bore, entry edge and entry node;
- it never finishes a pocket (a)–(b);
- zero-guide segments and retraces work only at a 180° cap;
- it takes seconds per motion;
- `replay_certificate` cannot return (`replay.py:1276`), and the replay grammar rejects `RetraceSegmentOperation` (`replay.py:244-261`).

Putting it in the gate would take:

1. Two-sided or coverage-driven placement (a coverage term in candidate choice, or centred circles).
2. General inter-route repositioning: a certified stay-down link from any phase point, or retract and plunge in the grammar.
3. Automatic entry choice.
4. A cap-general swept-prefix theorem, or dropping the zero-guide path at caps below π.
5. Fixing the native partition incompleteness and crash.
6. Removing the O(n²) revalidation and the double evaluation per commit.
7. A `ToolpathResult` adapter, plus D3 so the gate stops counting certified stay-down links as redundant.
8. Building the coverage seal.

Each of these is a stage. Together they amount to a new generator around the existing certifiers.

### (d) Could it be the basis of the product generator?

The generator is not a viable basis. Its components are.

- **Use:**
  - `MedialAxis` as the guide, or at least for clearance and reflex corners.
  - `NeckInventory` for neck detection.
  - `motion_refutation.refute_segment_cap` as an exact reject filter in the plan's station and bridge ladder.
  - `MotionCertifier.certify` as an *offline* final certificate of emitted segments and full circles (cost permitting).
  - `PreclearedEntry`/`QualifiedBore` for D1.
  - `ReachableDomain` regions for D2(a).
  - `GougeContainment` for exact gouge checks.
- **Do not use:**
  - the DFS route order (worse than the tree walk);
  - one-sided circles (no coverage);
  - zero-guide links and retrace (180° only, slots only);
  - the per-step transaction, commit and replay machinery (O(n²), seconds per step, no product consumer).
- **Better than the plan:**
  - an exact guide with exact clearance;
  - sound engagement (a continuous certificate instead of sampled probes);
  - a typed entry model;
  - certified stay-down links;
  - exact neck detection.
- **Worse than the plan:**
  - coverage;
  - ordering;
  - repositioning;
  - G1 continuity (link + full circle);
  - speed, by about 3–4 orders of magnitude per motion: 5–7 s per audit, against generation that the plan assumes runs at interactive speed (unverified for the plan side);
  - robustness (native aborts).

### (e) Performance as stated (none measured by me)

- Segment audit 5.1–6.8 s, 9-station ladder 0.033 s, about 160× apart (`motion_refutation.py:15-16`, `:194-203`).
- Four segment audits 4.96/5.06/6.68/6.97 s and circle audits 0.002–0.54 s; about 12 s of 24 s was duplicated work, now removed by a memo of 512 entries (`motion_oracle_cache.py:56-59`; `docs/continuous_engagement.md:722-730`).
- A 16-candidate family: 188.7 s with the authority alone, 35.0 s with the probe (`:780-785`).
- Task 13F route one: 17.66 s total, 16.64 s of it in segment audits (`:787-790`).
- Full adaptive suite: 518 tests in 223.88 s on 2026-07-29 (`:1004-1005`).
- Structural costs from code: double evaluation per commit, O(n) state rebuild per trial, O(n²) lineage. The swept-prefix path is not memoized (`motion_oracle_cache` docstring; `continuous_engagement.md:766-768`).
- All of the above are unverified by me. They are cited measurements by others.

---

### src/compas_cgal/adaptive/__init__.py | lines 1 | read 1-1
- **role:** build
- **used by:** package import of every `compas_cgal.adaptive.*` module
- **state:** live
- **generation:** n/a. Docstring only.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** "Strict typed domain for exact-certified adaptive machining": verified as a description of the contents.

### src/compas_cgal/adaptive/advancing_segment_trial.py | lines 71 | read 1-71
- **role:** engagement-control
- **used by:** `transaction.py:869` (zero-guide trial), `replay.py:1058`
- **state:** live (reachable only through zero-guide edges)
- **generation:** type-guards an `AdvanceSegmentOperation`, then delegates to the shared swept-prefix trial. Quality: *weak*. It only applies along 2r-wide slots at a 180° cap (see `motion_certificate.py:559-561`).
- **defects:** none found locally.
- **unused capability / ideas:** none beyond the swept-prefix theorem.
- **claims:** "exact capsule containment and swept-prefix TEA are proved on the pre-state before ... mutated" (`:29-31`): verified (`retrace_segment_trial.py:33-45` order).

### src/compas_cgal/adaptive/bootstrap.py | lines 689 | read 1-689
- **role:** entry
- **used by:** `generator.py:11-12`, `:1374`; Task 13F and branch fixtures
- **state:** live (test-only consumers)
- **generation:**
  - Launches from a `PreclearedEntry`: it removes the bore disk exactly (`:550`), certifies that the first full circle lies inside the bore disk (`:529-531`) and in the design (`:532-535`), certifies the circle's TEA against the post-bore stock (`:557-567`), then depletes and records coverage.
  - The operations are approach, plunge, circle (`:599-603`), with no lateral link.
  - Quality: *adequate* as an entry model (bore plus a contained first circle has no 360° virgin-stock slot), but the bore is a declared input, not generated.
- **defects:**
  - `:504-511` — `commit` recomputes the whole trial (stock clone, certification, depletion, coverage). This doubles the cost; the memo hides only the TEA call.
  - `:685-688` — the candidate's phase point must equal the bore centre *exactly*. The fixture therefore derives the entry centre from the candidate (`task13f_fixture.py:235`). A real bore location constrains the first candidate, not the reverse, so the launch must be searched for, and no search exists.
- **unused capability / ideas:** the bore-plus-first-circle certificate is the D1 entry contract the plan is looking for.
- **claims:** "First-and-only qualified bore depletion" (`:187`): verified (`stock_area.py:335-336`).

### src/compas_cgal/adaptive/candidates.py | lines 1248 | read 1-1248
- **role:** placement / spacing
- **used by:** `generator.py`, `traversal.py`, `transaction.py`, `bootstrap.py`, `replay.py`, `middle_curve_path.py`
- **state:** live (tests only)
- **generation:**
  - "One-sided MATHSM" circle proposals: a footpoint on one generator site, a contact point at distance r, centre at contact + ρ·dir, ρ ≤ (d−r)/2 (`:101-164`, `:254-369`).
  - A finite lattice: spatial progress in steps of `spatial_resolution/2^level` along the span (`:956-984`), radius steps (`:987-1010`) and `phase_count` phase rotations (`:167-184`).
  - Zero-guide candidates carry spatial progress only (`:1203-1248`).
  - Spans run along exact MAT samples, linear on line edges and following the parabola on P–S edges (`:1013-1064`).
  - Quality: *weak*.
    - The circles are one-sided, so the far side of a two-sided region is not reached (Slice finding 1).
    - The phase is a fixed lattice rather than tangent-continuous with the incoming link.
    - Terminal exhaustion stops at the last positive radius (`:1087-1126`), leaving corner material, the same D2 problem as the plan.
- **defects:**
  - `:1164-1165` and `:155-158` — a float `hypot` distance versus r is a *decision* that raises `InvalidMathsmProposalError`. That error is not caught by the search (`generator.py:1008-1016`), so a sample at the clearance root (d ≈ r) can abort the whole enumeration on a 1-ulp float error.
  - `:116-117` — the same pattern: a float parameter in (0,1) raises at segment-endpoint footpoints.
  - `:97-98` — `_distance` uses float `hypot` for `reported_length`, which then feeds an exact `Fraction` and identity. This is boundary injection, acceptable.
- **unused capability / ideas:** the parabolic-span interpolation (`:1023-1060`) is the exact way to place stations along a P–S MAT arc at reflex corners. The straight skeleton has no equivalent.
- **claims:** "Winner selection is invariant: furthest progress, largest radius, then canonical identity" (`policy.py:110-112`): verified (`policy.py:173-179`). "Greatest feasible progress is an exhaustive traversal boundary, not a sampled guess" (`:1094-1095`): contradicted in spirit. It is exhaustive over a *finite lattice* at a chosen resolution; finer lattice points are never tried.

### src/compas_cgal/adaptive/canonical.py | lines 811 | read 1-811
- **role:** tool (canonical CCAN byte encoding)
- **used by:** nearly every adaptive module; `engagement_audit/{records,identity,input,operation_identity}.py`
- **state:** live (also used by the engagement auditor)
- **generation:** n/a. Encodes identities and digests.
- **defects:** `:696-701` — `_exact_ring_area` uses exact `Fraction`; fine. No generation defect.
- **unused capability / ideas:** none for generation.
- **claims:** "frozen CCAN grammar" (`:240`): verified by the closed kind set at `:184-240`.

### src/compas_cgal/adaptive/containment.py | lines 534 | read 1-534
- **role:** kernel (exact gouge containment)
- **used by:** `transaction.py`, `bootstrap.py`, `entry.py`, `replay.py`, the retrace trials
- **state:** live
- **generation:** exact proofs that a segment capsule, a full-circle annulus or disk, an entry disk, or a circle inside an entry disk lies in the design `D` (`:328-534`). Quality: *excellent* as a primitive.
- **defects:** none found.
- **unused capability / ideas:** the gated generators and the plan could use `GougeContainment.certify_segment` / `certify_full_circle` to prove bridges and loops gouge-free exactly. The plan never mentions gouge checking.
- **claims:** "Single exact authority for segment, circle, and entry containment" (`:329`): verified within the package.

### src/compas_cgal/adaptive/coverage.py | lines 347 | read 1-347
- **role:** validation
- **used by:** `transaction.py`, `bootstrap.py`, `generation_state.py`, `replay.py`, the retrace trials
- **state:** live in the adaptive tests; this is the "exact coverage replay" that was set aside on 2026-09-08 for the product.
- **generation:** n/a. An exact residual ledger over the reachable material.
- **defects:** `:273-336` — every `add_sweep` clones the native ledger, rebuilds the full certificate, and recomputes and sorts the residual-component records (`:233-252`). The cost grows per step, which blocks generation-time use.
- **unused capability / ideas:** `require_complete` (`:344-347`) is the coverage gate the adaptive generator never calls.
- **claims:** none further.

### src/compas_cgal/adaptive/entry.py | lines 263 | read 1-263
- **role:** entry
- **used by:** `bootstrap.py`, `transaction.py`, `stock_area.py`, `identity.py`, `replay.py`
- **state:** live
- **generation:** `PreclearedEntry` is an exactly contained disk of radius > r (`:129-130`, `:183-184`) with a `QualifiedBore` binding the process identity and evidence bytes over the full clearance-to-cut interval (`:44-107`). It also provides the approach and plunge operations (`:204-217`) and the first-circle-in-bore certificate (`:239-263`). Quality: *adequate*. It models the entry correctly as a separate process outside the lateral cap, but it generates nothing: no helix and no ramp.
- **defects:** none found.
- **unused capability / ideas:** this is plan D1's natural home. A helical bore op could produce a `QualifiedBore` and reuse the containment and depletion contract unchanged.
- **claims:** "entry cannot launch from a one-tool-radius or smaller seed" (`:183-184`): verified.

### src/compas_cgal/adaptive/errors.py | lines 350 | read 1-350
- **role:** build (named exceptions)
- **used by:** every adaptive module; `src/compas_cgal/errors.py`; `tools/held_reference_qualification.py:22`
- **state:** live
- **generation:** n/a
- **defects:** none found. One named exception per failure mode, as the engineering bar asks.
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/adaptive/generation_state.py | lines 596 | read 1-596
- **role:** validation (authoritative per-step state)
- **used by:** `bootstrap.py`, `transaction.py`, `generator.py`, `retrace_transaction.py`
- **state:** live (tests only)
- **generation:** enforces the operation grammar:
  - entry circle first (`:377-378`);
  - every later circle paired with exactly one preceding hold link (`:384-416`);
  - advancing segments (`:450-461`);
  - retraces reversing the immediately preceding advance (`:418-448`);
  - phase continuity (`:385-386`, `:392-395`, `:455-459`).

  Quality: *n/a*. The grammar itself forbids retracts and repositioning, so multi-region linking cannot be expressed.
- **defects:** `:165-234` — every construction forks stock, clones coverage, revalidates the full chronology, computes the boundary digest and re-encodes all operations. That is O(n) per trial and O(n²) per path.
- **unused capability / ideas:** none for the product.
- **claims:** none further.

### src/compas_cgal/adaptive/generator.py | lines 1425 | read 1-1425
- **role:** ordering / linking (orchestrator)
- **used by:** tests only (`tests/adaptive/test_generator.py`, `task13f_fixture.py`, `test_route_retrace_generator.py`). No product or benchmark consumer.
- **state:** parked (memory: auditor-convergence and plumbing work is parked; its pinned end-to-end tests are expected-red, `docs/red_manifest-adaptive.json:4-15`).
- **generation:**
  - The loop is: materialize the forward-window family (`:1173-1268`), evaluate in order and take the first accepted candidate (`:962-1035`), commit with full re-evaluation (`:1038-1106`), activate incident routes (`:284-295`), and retrace at non-incident switches (`:1392-1418`).
  - It stops on any unresolved event.
  - Quality: *broken* as a product generator. It completes no pocket, the switch admits only a zero-guide source (`:513-557`), and coverage is never sealed (`:1362-1366`).
- **defects:**
  - `:1381-1384` — a non-incident boundary right after launch raises; there is no reposition move.
  - `:513-517`, `:533-536`, `:548-557` — the only legal retrace source is a final no-neck, full-cap zero-guide advance, which makes ordinary pockets unsupported.
  - `:1008-1016` — only gouge, cap and degenerate-link failures continue the search, so one `UnresolvedMotionEventError` aborts the whole run.
  - `:1165` — forward window: the family size is window × lattice. Each member may need a 5–7 s audit. Cost.
- **unused capability / ideas:** `_exhaustion_summary` (`:239-248`) gives exact rejection classes (cap, gouge, degenerate). A product generator should report the same per station.
- **claims:**
  - "This artifact proves traversal exhaustion, not complete pocket coverage" (`:1362-1366`): verified.
  - The red manifest pins the end state of `test_task13f_full_continuation` as "plan superseded": verified (`red_manifest-adaptive.json:4-8`).

### src/compas_cgal/adaptive/identity.py | lines 818 | read 1-818
- **role:** tool (input root identity, component versions, boundary-vertex IDs)
- **used by:** `InputIdentity` by `bootstrap.py` and `replay.py`. `ComponentIdentity`/`IdentityDigest` by `engagement_audit/{identity,records}.py`.
- **state:** live. The boundary-vertex ID classes (`:455-818`) have no src consumer, only `tests/adaptive/test_identity.py:26-37`.
- **generation:** n/a
- **defects:** none for generation. `:455-818` is about 360 lines with test-only consumers (dead in src).
- **unused capability / ideas:** none
- **claims:** "Content-addressed root of all validated Phase-1 planning inputs" (`:242`): verified (`:417-448`).

### src/compas_cgal/adaptive/medial_axis.py | lines 980 | read 1-980
- **role:** guide
- **used by:** `traversal.py`, `traversal_graph.py`, `candidates.py`, `neck.py`, `neck_topology.py`, `replay.py`, `generator.py`, `middle_curve_path.py`
- **state:** live (tests and adaptive only)
- **generation:**
  - A typed projection of the native exact segment-site MAT (`:537-622`, native build at `:599-606`): sites (point and open segment), nodes, edges (line or parabola, with dual kind), endpoint provenance flags (original Voronoi vertex, domain boundary, clearance root, `:64-66`), tool-fit runs, components and samples.
  - Each sample carries clearance and a reported guide radius (`:306-329`).
  - It also carries a zero-guide inventory with replayed native proofs (`:443-534`, `:960-979`).
  - Edges are rejected unless every exact flag is 1 (`:810-811`).
  - Quality: *excellent* as a guide substrate: exact clearance, parabolic reflex arcs, clipped at r, deterministic identities.
- **defects:**
  - `:145-147` — branch identity is SHA-256 of the edge id. The traversal orders branches by it (`policy.py:364-369`), so the geometric order is arbitrary.
  - Clearance and samples are float reports (`:783-787`), fine for proposals.
- **unused capability / ideas:**
  - The product generators use the straight skeleton (`toolpath.cpp`). This MAT gives exact clearance and true reflex-corner geometry.
  - The plan's root rule ("largest exact clearance") is directly computable from `MatSample.clearance` or the node clearance.
- **claims:** samples "must cover every exact edge" (`:562-565`): verified.

### src/compas_cgal/adaptive/middle_curve_path.py | lines 63 | read 1-63
- **role:** placement
- **used by:** `tests/adaptive/test_middle_curve_path.py` only. `benchmarks/` and `tools/` have no consumer on this branch.
- **state:** dead in src (test-only). It is probably a seed for the Held MATHSM reproduction on another branch (unverified).
- **generation:** the maximum-radius, boundary-phase circle at each native MAT sample of one edge for one generator site (`:11-63`). Quality: *weak*. It has no engagement control and is one-sided.
- **defects:** none locally.
- **unused capability / ideas:** this is a MATHSM-style proposal per edge; it could seed a fast, uncertified first pass.
- **claims:** "paper-derived machining-circle proposals" (`:1`): unverified (no paper equation cited in the file).

### src/compas_cgal/adaptive/motion_certificate.py | lines 633 | read 1-633
- **role:** engagement-control
- **used by:** `transaction.py`, `bootstrap.py`, `retrace_segment_trial.py`, `replay.py`, `replay_trace.py`, `generation_state.py`, `identity.py`
- **state:** live
- **generation:**
  - `MotionCertifier.certify` is an exact event-partition TEA certificate for a whole segment or full circle (`:452-540`), memoized (`:497-504`).
  - `refute_segment` is a cheap exact counterexample search (`:410-450`).
  - `certify_swept_prefix_segment` proves engagement ≤ π for an advancing segment from a clear start (`:542-633`).
  - Quality: *excellent* soundness and *broken* cost and completeness. It takes 5–7 s per segment, has unresolved substrate cases, and a native abort.
- **defects:**
  - `:559-561` — the swept-prefix theorem "currently admits only pi", so zero-guide moves are unusable at 120° or 40°.
  - `:516-527` — substrate errors map to `UnresolvedMotionEventError`, which aborts generation.
- **unused capability / ideas:**
  - Use `certify` as the final continuous certificate of a product path's segments and full circles, offline.
  - Use `refute_segment` inside product candidate search.
- **claims:**
  - "certified / cap_exceeded / unresolved distinct" (`:80-85`, `:524-527`): verified.
  - `refute_segment` "costs roughly two orders of magnitude less" (`:418-421`): unverified (cited measurement 160×, `motion_refutation.py:194-203`).

### src/compas_cgal/adaptive/motion_oracle_cache.py | lines 284 | read 1-284
- **role:** tool (performance memo)
- **used by:** `motion_certificate.py:20`, `:497`; `engagement_audit/records.py:20` (type only)
- **state:** live
- **generation:** n/a. A content-keyed LRU of 512 entries over native `(verdict, trace)` pairs (`:56-59`, `:187-265`).
- **defects:** none found. It does not cover `certify_swept_prefix_segment`, by design.
- **unused capability / ideas:** none
- **claims:** "roughly 12 s of 24 s spent re-deriving" (`docs/continuous_engagement.md:722-730`): unverified (cited).

### src/compas_cgal/adaptive/motion_refutation.py | lines 406 | read 1-406
- **role:** engagement-control
- **used by:** `motion_certificate.py:21-22`, `:443`; `transaction.py:1091`
- **state:** live
- **generation:** a dyadic 9-station exact ladder (`:164-206`). A refuted station is a proof of cap exceedance (`:289-347`). Quality: *excellent* as a rejection filter. It is sound, cheap, and float-free in its decision.
- **defects:** `:340-344` — native substrate failures yield `UNKNOWN`, which is correctly inconclusive. None found.
- **unused capability / ideas:** this is the strongest drop-in for the product generators. Replace sampled over-cap rejection in the ρ-ladder scan and bridges with exact refutation; acceptance stays sampled unless certified.
- **claims:**
  - "5.1 s to 6.8 s for one full segment audit against 0.033 s for the whole nine-station ladder" (`:15-16`): unverified (cited, Task 13F fixture).
  - "depth 3 refutes every violation spanning more than an eighth of the motion" (`:194-196`): verified as arithmetic.

### src/compas_cgal/adaptive/neck.py | lines 420 | read 1-420
- **role:** guide (neck classification) / engagement-control (per-pass caps)
- **used by:** `traversal.py`, `bootstrap.py`, `transaction.py`, `generation_state.py`, `replay.py`
- **state:** live (tests only)
- **generation:** native-validated neck evidence classified into exact squared-width classes (`:288-363`). Each neck has an oriented passage state machine (UNVISITED → FIRST → SECOND → TERMINAL, `:39-43`), and each state supplies an effective cap (`:216-233`). Quality: *adequate* as a capability. Engagement control in narrow passages is world-class-relevant, and the gated generators have nothing comparable.
- **defects:** `:219-221` — a terminal passage raises, so at most three passes through any neck are allowed. This is a hidden traversal limit.
- **unused capability / ideas:** the plan's dumbbell (Review Focus 2) could use `NeckInventory` to find the neck exactly rather than discovering it from a failed trochoid.
- **claims:** none further.

### src/compas_cgal/adaptive/neck_topology.py | lines 332 | read 1-332
- **role:** guide
- **used by:** `neck.py`, `traversal.py`
- **state:** live
- **generation:** typed neck loci (strict edge, clearance endpoint, shared vertex, plateau) and neck sides as MAT edge partitions (`:53-332`). Quality: *adequate* (pure topology).
- **defects:** none found.
- **unused capability / ideas:** the side partition gives which MAT edges lie beyond a neck, which is what a product walk needs to schedule a second entry.
- **claims:** none

### src/compas_cgal/adaptive/operation.py | lines 1053 | read 1-1053
- **role:** linking / motion-geometry (operation grammar)
- **used by:** most adaptive modules
- **state:** live
- **generation:**
  - Operation types: approach, plunge, `LinkSegmentOperation` (hold link at cut depth, `:674-726`), `AdvanceSegmentOperation` (`:729-853`), `RetraceSegmentOperation` (exact reverse of the source, `:856-989`) and `CutFullCircleOperation` (`:992-1050`).
  - There is no retract, no arc and no helix.
  - Cap decisions are full or per-neck (`:159-272`).
  - Quality: *weak*. Without partial arcs, every loop is a full 360°. Without retract, repositioning is impossible.
- **defects:**
  - `:51` — `ROUTE_RETRACE_STRATEGY_VERSION = b"reverse-final-zero-guide-v1"`, the only strategy.
  - `:946-953` — retrace requires a no-neck, full-cap, terminal source.
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/adaptive/policy.py | lines 395 | read 1-395
- **role:** spacing / ordering policy
- **used by:** most adaptive modules; `engagement_audit/input.py`
- **state:** live
- **generation:**
  - `CandidatePolicy`: lattice bounds plus the invariant order — furthest progress, then largest radius, then identity (`:106-180`).
  - `NeckPolicy`: per-class, per-pass caps, each ≤ the user cap, compared exactly (`:239-308`).
  - `TraversalPolicy`: the forward window plus branch and component ordering by *byte identity* (`:344-369`).
  - `CutDirectionPolicy`: climb or conventional mapped to circle orientation (`:372-395`).
  - Quality: *weak* for ordering (hash order) and *adequate* for cap policy.
- **defects:** `:357-369` — component and branch ordering is `sorted(bytes)`, with no geometric criterion (widest first, nearest). It determines the route.
- **unused capability / ideas:** the climb/conventional orientation mapping is reusable.
- **claims:** "Winner selection is invariant ... deliberately not a configurable policy field" (`:110-112`): verified.

### src/compas_cgal/adaptive/reachable_domain.py | lines 622 | read 1-622
- **role:** kernel (exact centre domain and reachable material)
- **used by:** `containment.py`, `coverage.py`, `entry.py`, `identity.py`, `transaction.py`, `replay.py`
- **state:** live
- **generation:** native `ReachableDomain2` with an exact `design_region`, `center_domain` (C_r), `reachable_material` and `unreachable_residual` (`:536-622`), plus a validated structural certificate (`:98-428`). Quality: *excellent* as a primitive.
- **defects:** `:484-485` — "reachable-domain reconstruction requires one exact component". Pockets whose centre domain is disconnected (necks narrower than 2r) cannot be planned at all.
- **unused capability / ideas:** `unreachable_residual` is the exact region D2(a) needs to exclude from "uncut". The gate's uncut metric could be made exact against it.
- **claims:** none further.

### src/compas_cgal/adaptive/replay.py | lines 1276 | read 1-1276
- **role:** validation
- **used by:** `tests/adaptive/test_replay.py`, `test_zero_guide_replay.py`
- **state:** parked. It cannot succeed: `:1276` raises unconditionally after a full replay ("must bind the immutable certificate before return").
- **generation:** n/a. It rebuilds the inputs, the MAT and the candidate stream independently, then replays every operation with containment, motion certificate, depletion and coverage, and requires complete coverage (`:1151-1166`).
- **defects:**
  - `:1276` — no value-returning path (also recorded in `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md:52`).
  - `:244-261` — the grammar rejects `RetraceSegmentOperation`, so any program that uses the generator's retrace cannot replay.
- **unused capability / ideas:** none for generation.
- **claims:** none

### src/compas_cgal/adaptive/replay_trace.py | lines 332 | read 1-332
- **role:** validation
- **used by:** `transaction.py`, `generator.py`, `bootstrap.py`, `replay.py`, the retrace modules
- **state:** live (as the per-operation proof bundle)
- **generation:** n/a. `ReplayLateralWitness` cross-validates motion, containment, depletion and coverage per operation (`:50-170`). `FreshReplayTrace` validates the chronology (`:173-332`).
- **defects:** none for generation.
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/adaptive/retrace_segment_trial.py | lines 114 | read 1-114
- **role:** linking
- **used by:** `advancing_segment_trial.py`, `retrace_transaction.py`
- **state:** live
- **generation:** the shared ordered trial for swept-prefix segments: containment, certification, depletion, coverage (`:20-55`). Quality: *weak* (180° only, zero-guide only).
- **defects:** none locally.
- **unused capability / ideas:** none
- **claims:** "The generic event certifier is not part of this deciding path" (`:72-74`): verified (`:38-43`).

### src/compas_cgal/adaptive/retrace_transaction.py | lines 402 | read 1-402
- **role:** linking
- **used by:** `generator.py:51-52`, `:1371`, `:1400-1407`
- **state:** live (tests only)
- **generation:** builds the reverse of the final zero-guide advance and certifies it at the user cap (`:337-402`), so the tool returns along the slot it just cut. Quality: *weak*. It is a slot-only special case. A retrace returns only to the source segment's start, not to an arbitrary next-route entry.
- **defects:**
  - `:326-334` — commit recomputes the whole trial.
  - `:384` — `effective_cap=self.evaluator.user_cap` is required to be π by the theorem (`motion_certificate.py:561`), so this path is dead at gate caps.
- **unused capability / ideas:** none
- **claims:** none

### src/compas_cgal/adaptive/stock_area.py | lines 449 | read 1-449
- **role:** kernel (exact stock depletion)
- **used by:** `bootstrap.py`, `transaction.py`, `motion_certificate.py`, `generation_state.py`, `replay.py`, the retrace trials
- **state:** live
- **generation:** exact native stock subtraction for entry disks, segments and full circles, with validated structural traces (`:203-449`). Quality: *adequate*.
- **defects:**
  - `:402` — each witness stores the tuple of all parent digests, which is O(n²) memory and hashing over a path.
  - `:224`, `:251`, `:268` — every construction, `raw` access and fork deep-clones native stock.
- **unused capability / ideas:** none new. The product side already has `Stock` / `_subtract_operation`.
- **claims:** none

### src/compas_cgal/adaptive/transaction.py | lines 1223 | read 1-1223
- **role:** linking / engagement-control
- **used by:** `generator.py`, `retrace_transaction.py`
- **state:** live (tests only)
- **generation:**
  - Middle-curve candidate: a hold LINK from the current phase point to the circle's phase point, then a full circle. Each motion is containment-checked, refuted or certified, depleted and swept (`:755-840`, `:1066-1169`).
  - Zero-guide candidate: one advancing segment (`:842-901`).
  - Quality: *weak*. The links are straight and not tangent to the circles, so G1 breaks at every junction, and every step is a full 360°.
- **defects:**
  - `:698-705`, `:744-752` — commit recomputes the full trial.
  - `:946-962` — `_validate_physical_parent` forks the whole stock lineage per trial and checks every witness policy, O(n).
  - `:964-1021` — it revalidates every prior operation per trial, O(n).
- **unused capability / ideas:**
  - The certified stay-down LINK (`:1066-1125`) is the proof-carrying form of plan D3.
  - The refute-then-certify ordering (`:1083-1106`) is the right pattern for a product candidate search.
- **claims:** "Probing the fixed station ladder first therefore rejects a violating link for a fraction of a percent of the audit it replaces" (`:1083-1090`): unverified (cited A/B 5.39×, `docs/continuous_engagement.md:780-785`).

### src/compas_cgal/adaptive/traversal.py | lines 1074 | read 1-1074
- **role:** ordering
- **used by:** `generator.py`, `bootstrap.py`
- **state:** live (tests only)
- **generation:** a global route ledger. There is one cursor per directed MAT edge; only the active one advances (`:705-743`). Activation follows the static route (`:745-787`), and the causal neck-side frontier is precomputed per route step (`:976-1024`). Quality: *weak*. The route is a DFS in hash order with each edge exactly once, has no cost model, and a backtrack forces a retrace.
- **defects:**
  - `:971-972` — `OverlappingNeckTransitError` fires if one activation crosses two necks, a hard failure on nested necks.
  - `:624-669` — the entry edge and node are caller inputs with no selection logic.
- **unused capability / ideas:** precomputed causal neck frontiers along a route (`:976-1024`) would give a product walk "this edge is beyond neck N" for free.
- **claims:** none

### src/compas_cgal/adaptive/traversal_graph.py | lines 583 | read 1-583
- **role:** ordering
- **used by:** `traversal.py`
- **state:** live (tests only)
- **generation:**
  - Projects MAT edges into a static graph (`:368-409`).
  - `route` puts the entry component first, then the others in sorted-ID order (`:460-465`).
  - Per component, an iterative DFS over incident edges in branch-hash order (`:528-583`). It explores the exit side of the root edge, then *returns to the root's entry node* (`:579-580`). That is a guaranteed non-incident switch whenever the root's entry node has other edges.
  - Quality: *weak*. Compare the plan's tree walk, which roots at the widest node.
- **defects:**
  - `:579-580` — the second `explore` from the root's entry node always needs a retrace that only zero-guide sources can supply. This is why Task 13F needs a retrace at route 2 (`test_generator.py:571-573`).
  - `:477-480` — the root of every non-entry component is the hash-smallest branch.
- **unused capability / ideas:** none beyond the plan's own tree walk.
- **claims:** "One directed step for every graph edge" (`:441-442`): verified (`:489-490`).
