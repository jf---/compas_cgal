# Census slice g5_tea_and_exact: continuous TEA certification and exact number types

## Slice summary

Bottom line: `continuous_tea_2` is a sound but incomplete exact **auditor** of one motion against **frozen pre-motion stock**. A generator can use it to accept or reject a candidate it has already proposed, which is what `adaptive/` does with a finite lattice. It cannot be asked for a maximal advance or a maximal loop radius. It is about 100× too slow to replace the 33-probe check on every rung of a ladder. The swept-annulus guard, the false-certificate history and the shared-root fix are **not in this slice**. They live in `src/engagement_2.cpp` (commits 2e2422b5, f227dab2, 148a49ad, verified with `git show --stat`). That lane, not this one, is the certifier the repo's own cost verdict says to ship (`docs/continuous_engagement_cost.md:118-123`).

#### Top findings, ranked by impact on toolpath quality

1. **Frozen-stock semantics make the generic segment oracle useless for cutting segments at caps below 180°.** `audit_segment_tea_event_exact` classifies every cell against the stock as given (`segment_oracle.cpp:505-513`, `station_classifier.cpp:93-104`). For an advancing cut, the trailing rim leaves the start disk at t > 0 and counts as material, so any slot or bridge into stock is refuted. The swept-prefix theorem certifies only an exact π cap (`segment_oracle.cpp:220-229`, `:598-602`). **Unbuilt capability (derived here, unverified by test):** for a translation, a forward rim point q = p(t) + r·u with u·v > 0 satisfies |q − p(s)|² = (t−s)² + 2r(t−s)(u·v) + r² > r² for every s < t. Forward points are therefore never inside an earlier disk of the same motion. For t > 0, the physical engaged set is exactly the closed forward semicircle intersected with the frozen stock. Masking the station classifier to the forward semicircle would give **exact** segment TEA at any cap with the existing machinery, and every run would be at most π, so the chord test is monotone. No such mask exists in `station_classifier.cpp`.
2. **No API yields a certified maximal advance.** All entry points are verdicts on a given motion (`continuous_tea_2.cpp:1464-1482`, `:1577-1672`). The exact event roots and cells that would define "first violating parameter" are computed internally (`exact_algebraic_1.cpp:1240-1277`). However, `SegmentTeaAudit2::violating_parameter` is the interior dyadic **witness** of the first violating cell, not its lower root (`segment_oracle.cpp:527-538`), so it lies past the first violation. The Python binding also drops it (`continuous_tea_2.cpp:1596-1599`). A sound maximal advance would be the lower root of the first non-certified cell, or its rational `interval_low` (`exact_algebraic_1.cpp:571-578`). For a loop radius there is no parameter at all: the radius is fixed per audit (`parameter_charts.cpp:324-365`).
3. **Cost (measured here):** about 0.31–0.72 s per full-circle audit against 3–5 ms for the 33 frozen-stock probes. Both run on a 12×8 rectangle with 8–16 boundary features, r = 1, ρ = 0.5, cap 120° (`pixi run python` probe below), so the audit is about 100–150× slower. The verdicts agreed on 6/6 step sizes: steps 0.05–0.2 were certified with 0 of 33 probes over cap, and steps 0.25 and 0.3 were cap_exceeded with 5 and 9 of 33 probes over cap. The repo states about 15 s per segment audit, superlinear in contacted features (`docs/continuous_engagement_cost.md:13-35`; unverified by me). The cost drivers I found in code:
   - bivariate `Algebraic_kernel_d_2` solves per fibre (`segment_fibre.cpp:835-838`);
   - bivariate signs for every ordered branch pair at every cell (`segment_strata.cpp:667-747`), where the univariate `exact::sign_mixed_radical` (`exact/one_root.cpp:45-105`) already decides the same predicate in the engagement_2 lane;
   - O(n²) feature-pair trivariate resultants (`segment_pair_projection.cpp:347-412`; lines are never reach-filtered, `segment_projection.cpp:527-529`);
   - full reconstruction for verification on every call (`segment_oracle.cpp:562-566` → `segment_partition.cpp:785`; `circle_oracle.cpp:804`; `event_trace.cpp:268`).
4. **Incompleteness rate gates its use as a placement gate.** The repo measured 15/125 link audits raising `IncompleteSegmentPartitionError` and one native abort, with 0 false certificates (`docs/refutation_soundness_crosscheck.md:1-21`; unverified by me). A placement loop needs a defined action on `UNRESOLVED`. `segment_projection.cpp:751-754` even throws when the reach filter leaves nothing to test, instead of returning "clear".
5. **The circle oracle has no cheap refuter.** `full_circle_rational_probe_exceeds_cap_exact` is a stub that throws unless the stock is empty (`circle_oracle.cpp:926-968`). The segment lane has one (`segment_station_cap_exceeded_exact`, `segment_oracle.cpp:634-663`). For loops, `_stock_2.engagement_at` is the de facto refuter, under the same frozen semantics.
6. **Exhaustive pair path on all-rational stocks.** When every boundary vertex is rational, the circle oracle builds n(n−1)/2 × 16 pair requests (`circle_oracle.cpp:740-748` → `circle_source.cpp:313-320` → `circle_pair_projection.cpp:613-663`) instead of the demand-closed set (`circle_strata.cpp:447-497`). That is quadratic blow-up on integer and dyadic fixtures; the doc's 670,480-chain example is at `docs/continuous_engagement.md:339-349`.

#### Generation capabilities here that the gated generators do NOT use

- `audit_full_circle_tea_event_exact` (`circle_oracle.cpp:467-509`): an exact, continuous decision with the **same frozen full-rim semantics** as `_station_is_admissible` (`engagement_toolpath.py:525-555`). It is a drop-in exact replacement for the 33-probe verdict on a loop.
- `audit_swept_prefix_segment_tea_exact` (`segment_oracle.cpp:578-632`): exact certification of an advancing straight cut from a clear start disk at a π cap. It is usable for zero-guide slot links. The gated generators use straight run-outs and bridges without any certificate.
- `segment_station_cap_exceeded_exact` (`segment_oracle.cpp:634-663`): an exact one-station refuter that costs about 3.2 ms per probe (stated, `docs/continuous_engagement.md:775-784`).
- The refute-then-certify pattern (`adaptive/transaction.py:1083-1106`): cheap exact probes reject; the full certificate runs only on survivors. Applied to the gated ladder scan, that is one exact audit per emitted loop.
- The exact one-root sign `exact::sign_mixed_radical` (`exact/one_root.cpp:45-105`). It is not used by `continuous_tea_2`, which pays for the bivariate kernel instead.

#### Where this slice shows the draft plan (`docs/superpowers/plans/2026-10-13-product-gate-tree-walk.md`) is wrong or incomplete

- **The cap decision stays sampled while the plan bans sampling elsewhere.** Global Constraints (plan:58) require exact removal, plunge and root decisions with "no sampled-engagement shortcuts in those decisions". Task 3 still accepts stations through `_station_is_admissible` (plan:988, :1103, :1112), the 33-probe frozen check. The plan should state explicitly that between-probe cap misses remain possible on the tree walk, or it should certify the accepted rung. An exact drop-in exists, as listed above.
- **Certified advances will not move a Stage A gate number.** The plan attributes the over-cap motions to chain-entry loops (plan:43), not to probe gaps. Under the plan's own rule (plan:16-18), certification belongs in Stage C, with a first measurement task: on the 3 gate pockets, count the loops the 33 probes accept and `audit_full_circle_tea_event_exact` refutes. If that count is 0, park it.
- **Bridges and run-ins at 120°/40° have no exact certificate in the repo.** B2a's wall run-in (plan:1383-1384) is a straight move into material. The generic segment oracle refutes such moves by construction (frozen stock), and swept-prefix only reaches a π cap. Measuring B2a "exactly" needs either the engagement_2 lane or the forward-semicircle mask from finding 1. Neither is named in the plan.
- **D1 helix:** the motion vocabulary of the exact oracles is only "segment" and "full-circle" (`parameter_charts.cpp:390-416`). A helical or ramped entry has no certificate here, which is consistent with D1 excluding it from the planar cap. Adding one is new work, not reuse.
- **The plan does not mention `adaptive/`, the one generator already driven by these certificates.** Its first-feasible search over a finite lattice of (progress × guide radius × phase) (`adaptive/candidates.py:1129-1200`, `adaptive/generator.py:960-1035`) is the existing "certified placement" design. The plan should record why it is not reused: no product consumer, and cost.

#### Timing probe (verified, read-only)

```
pixi run python -c "... Stock(rect 12x8); subtract_disk loops at rho+r=1.5 every 0.25/0.375;
audit_full_circle_tea_event_exact(raw, xn, 4, 0.5, 0, False, 1.0, 3.0) vs 33x _stock_2.engagement_at ..."
features 6..16 : 1.24 s (first call, warm-up), then 0.313, 0.371, 0.46, 0.499, 0.593 s (all cap_exceeded at step 0.375)
step 0.05 certified 0.558 s | probes 0/33 0.005 s
step 0.10 certified 0.554 s | probes 0/33 0.004 s
step 0.15 certified 0.502 s | probes 0/33 0.004 s
step 0.20 certified 0.532 s | probes 0/33 0.004 s
step 0.25 cap_exceeded 0.474 s | probes 5/33 0.004 s
step 0.30 cap_exceeded 0.723 s | probes 9/33 0.003 s
```

---

### src/continuous_tea_2.cpp | lines 1684 | read 1-1684

- **role:** engagement-control (nanobind module `_continuous_tea_2`)
- **used by:**
  - `adaptive/motion_refutation.py:329` (`segment_station_cap_exceeded_exact`)
  - `adaptive/motion_oracle_cache.py:151,172` (segment and full-circle audits)
  - `adaptive/motion_certificate.py:33-34,43,590` (swept-prefix, `extract_boundary_records`)
  - `adaptive/identity.py:226`
  - every other binding is used only by tests (grep: prod=0 for `extract_boundary_events`, `classify_boundary_pair`, `parameter_charts`, `construct_pullback`, `projection_from_grid`, `partition_pullback_overlap`, `solve_trimmed_line_branches`, `project_regularization_vertex`, `classify_ccw_orientation`, `partition_cap_crossings`, `verify_event_partition`, `construct_segment_*`, `verify_segment_event_partition`, `segment_pair_literal_signs`, `full_circle_rational_probe_exceeds_cap_exact`)
- **state:** parked. Its only production consumer is `adaptive/`, which has no consumer outside its own package and tests.
- **generation:** exposes verdicts only. Quality: weak for generation. There is no maximal-advance or maximal-radius query.
- **defects:**
  - `:1596-1599`: the binding drops `violating_parameter`, so Python cannot even locate the first violation.
  - `:1484-1495`: binds a stub (see circle_oracle).
  - `:1524-1528`, `:1679-1683`: test mutators are compiled into the production module.
- **unused capability / ideas:** expose `first_uncertified_lower_root` (from the cells and roots already in the trace) as a certified-advance query.
- **claims:**
  - `BOUNDARY_EVENT_KINDS` (`:1167-1172`): verified against `boundary_events.cpp:741-749`.

### src/continuous_tea_2/bind_partition_records.cpp | lines 347 | read 1-347

- **role:** engagement-control (bindings)
- **used by:** `continuous_tea_2.cpp:921`
- **state:** live with the module (parked lane)
- **generation:** n/a (record bindings)
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/bind_partition_records.h | lines 6 | read 1-6

- **role:** engagement-control
- **used by:** `continuous_tea_2.cpp:3`
- **state:** live (parked lane)
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/boundary_events.cpp | lines 890 | read 1-890

- **role:** engagement-control (exact stock-boundary feature extraction and identities)
- **used by:**
  - all oracles: `segment_partition.cpp:399`, `circle_oracle.cpp:525`, `circle_strata.cpp:372,517`, `segment_strata.cpp:839`, `segment_oracle.cpp:197`
  - `audit_stock_state_identity_2.cpp:4`
  - `adaptive/motion_certificate.py:43`
- **state:** live (audit lane plus parked adaptive)
- **generation:** extracts every curve of the whole stock and builds content-addressed feature and vertex ids. Adequate as a substrate.
- **defects:**
  - `:358-479`: vertex identity and ordinal construction is O(V²) with linear `find_if` (`:373`, `:419-434`).
  - `:512-517`: overlap multiplicity is O(n²).
  - `:829-849`: all-pairs events.
  - No spatial culling: features far from the motion are carried into every audit. The doc measured that distant features cost little (`continuous_engagement_cost.md:41-43`, unverified).
- **unused capability / ideas:** a locality filter at extraction (the swept stadium plus r) would cut n for every downstream O(n²) step.
- **claims:** none

### src/continuous_tea_2/boundary_events.h | lines 49 | read 1-49

- **role:** engagement-control
- **used by:** see `boundary_events.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/cap_partition.cpp | lines 723 | read 1-723

- **role:** engagement-control (rational-chart cap-crossing partition for a pair of rim branches)
- **used by:**
  - `circle_oracle.cpp:863-871` (scaffold for boundaryless stock)
  - `event_certificate.cpp:581` (reconstruct)
  - bindings, test-only
- **state:** live in the scaffold path only
- **generation:** n/a
- **defects:**
  - `:668-671`: mixed orientation dispositions throw rather than split. Effect is limited, since this path is reached only by the scaffold and tests.
- **unused capability / ideas:** none
- **claims:**
  - squared-chord surrogate `4·diff² − cap·D1·D2` (`:242-250`): verified as the chord² vs `cap·r²` identity on the rational unit circle.

### src/continuous_tea_2/cap_partition.h | lines 28 | read 1-28

- **role:** engagement-control
- **used by:** `circle_oracle.cpp:6`, `event_certificate.cpp:4`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_fibre_state.cpp | lines 716 | read 1-716

- **role:** engagement-control (full-circle seam coalescing, active-sheet split/merge classification)
- **used by:** `circle_projection.cpp:1346`
- **state:** live in the circle oracle
- **generation:** n/a
- **defects:**
  - `:180-182`: `adjacent_cell` silently falls back to the last or first cell when no neighbour matches. This is cyclic-seam logic; it could mislabel the direction if the cells are malformed, but the verdict comes from cell authority, not from fibres.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_fibre_state.h | lines 6 | read 1-6

- **role:** engagement-control
- **used by:** `circle_projection.cpp:3`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_oracle.cpp | lines 968 | read 1-968

- **role:** engagement-control (exact full-circle-loop TEA audit)
- **used by:**
  - `adaptive/motion_oracle_cache.py:172` (CUT certification, `adaptive/transaction.py:1136-1150`)
  - `audit_certification_2.cpp:63`
  - `audit_certification_replay_2.cpp:40`
  - `audit_unresolved_evidence_2.cpp:203,504,551`
- **state:** live in the audit lane and parked adaptive
- **generation:** decides whether a full loop (fixed centre, phase radius |phase|, direction) stays under the cap against **frozen** stock. This is sound in the certified direction, because stock only shrinks during the loop. It is conservative, because the trailing material the loop itself cuts is counted as engaged. Quality: adequate as a per-loop gate (same semantics as `_station_is_admissible`), weak as a generator primitive (no radius or phase optimisation). Measured at 0.31–0.72 s per call (probe above).
- **defects:**
  - `:740-748`: all-rational vertices route to exhaustive O(n²×16) pair requests.
  - `:764-773`: the cheap CAP_EXCEEDED witness is used only when *all* features are lines.
  - `:804-810`, `:821-832`: verification rebuilds the partition, and `build_event_trace` rebuilds again on the uniform path (`:542`, `:600-609` → `event_trace.cpp:268`).
  - `:926-968`: `full_circle_rational_probe_exceeds_cap_exact` is a stub that always throws on non-empty stock, so there is no exact cheap loop refuter.
- **unused capability / ideas:** an exact replacement for the 33-probe acceptance in `engagement_toolpath.py:525-555`, to be called once per accepted rung.
- **claims:**
  - full-circle oracle cost "1.193 s for 2 audits" (`docs/continuous_engagement.md:783`): consistent with my 0.3–0.7 s measurement (verified order of magnitude).

### src/continuous_tea_2/circle_oracle.h | lines 58 | read 1-58

- **role:** engagement-control
- **used by:** `audit_certification_2.cpp:3`, `audit_certification_replay_2.cpp:3`, `audit_certification_internal_2.h:4`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** `FullCircleTeaAudit2::violating_parameter` (`:18`) exists in C++ but is not exposed to Python.
- **claims:** none

### src/continuous_tea_2/circle_pair_projection.cpp | lines 805 | read 1-805

- **role:** engagement-control (full-circle cross-feature orientation and cap resultants)
- **used by:** `circle_projection.cpp:1332`, `circle_source.cpp:318`, `circle_strata.cpp:421`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:613-663`: `full_circle_exhaustive_pair_requests` produces n(n−1)/2 × 4 × 2 × 2 requests, each costing two trivariate resultants (`:177-206`). This is the quadratic cost centre on the rational path.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_pair_projection.h | lines 80 | read 1-80

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_projection.cpp | lines 1367 | read 1-1367

- **role:** engagement-control (full-circle pullbacks: tangency, cap-crossing, vertex-passage and coincidence projections over 4 centre charts)
- **used by:** `circle_source.cpp:322,358`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:957-1070`, `:1101-1254`: builds 8 pullbacks per feature plus per-chart factors for **every** boundary feature, with no reach filter (contrast `segment_projection.cpp:523-552`). Effect: cost scales with the whole stock boundary.
  - `:538-541`, `:450-453` (in `segment_fibre.cpp` too): distinct radicals in x and y raise `UnsupportedAlgebraicVertexProjectionError`. A vertex from two different circle intersections makes the audit unresolved.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_projection.h | lines 23 | read 1-23

- **role:** engagement-control
- **used by:** `circle_oracle.cpp:7`, `circle_source.cpp:6`, `circle_strata.cpp:6`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_source.cpp | lines 378 | read 1-378

- **role:** engagement-control (uniform clear/material shortcut, rational witness, partition constructors)
- **used by:** `circle_oracle.cpp:528-541,741`, `event_certificate.cpp:534,547,564`, `audit_unresolved_evidence_2.cpp:5`
- **state:** live
- **generation:** the swept full-circle overlap test (`:98-111`) is an exact fast path for "loop entirely in cleared stock".
- **defects:**
  - `:112-119`: the "material" fast path works only for an axis-aligned phase vector; any oblique phase falls through to the full partition. Effect is cost only.
- **unused capability / ideas:** `reach_full_circle_sweep ∩ stock = ∅` (`:98-111`) is an exact "loop removes nothing" test. It is the same question as the plan's productive-only emission, at the loop level and without a stock clone.
- **claims:**
  - `:140-153` (t = 1/2 gives the unit direction (3/5, 4/5)): verified (half-angle chart).

### src/continuous_tea_2/circle_source.h | lines 54 | read 1-54

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_strata.cpp | lines 773 | read 1-773

- **role:** engagement-control (demand-closed pair requests, per-cell station authority)
- **used by:** `circle_oracle.cpp:749,811`
- **state:** live
- **generation:** classifies each parameter cell by one rational station on the loop (`:553-579`). This is the exact continuum version of the 33 probes, under the same frozen-stock semantics.
- **defects:**
  - `:464-496`: the fixed-point loop rebuilds the whole partition until the pair set closes. It terminates because the set is finite, but each iteration is a full rebuild.
  - `:635-643`: `violating_parameter` is the first violating cell's interior witness, not its boundary root.
  - `:187-191`: decode failure is swallowed into `{}`. This is the conservative direction (it ends as UNRESOLVED).
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_strata.h | lines 78 | read 1-78

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_vertex_source.cpp | lines 167 | read 1-167

- **role:** kernel (one-root coordinate canonical form for full-circle vertices)
- **used by:** `circle_oracle.cpp:630-644`, `circle_projection.cpp:527-529`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:56-65`: exact rational obtained by printing `CGAL::exact` to text and re-parsing. This is the text-carrier pattern, which the cost doc says is not the driver (`continuous_engagement_cost.md:45-49`, unverified).
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/circle_vertex_source.h | lines 43 | read 1-43

- **role:** kernel
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_certificate.cpp | lines 1085 | read 1-1085

- **role:** engagement-control (canonical certificate encoding, reconstruct-verify, mutation API)
- **used by:** every oracle, `audit_*_2.cpp`, `exact/canonical.cpp:4`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:643-678`: verification means **full reconstruction from source**, so each verified audit costs at least 2×.
  - `:672-677`: any substrate error becomes UNRESOLVED (conservative).
  - `:680-1085`: about 400 lines of test-only mutators in a production translation unit.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_certificate.h | lines 26 | read 1-26

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_partition.cpp | lines 433 | read 1-433

- **role:** engagement-control (`partition_projections` wrapper, pullback-overlap, regularization-vertex projection)
- **used by:**
  - `partition_projections`: all oracles
  - `partition_pullback_overlap`, `project_regularization_vertex`: test-only bindings
- **state:** live (wrapper); the rest is test-only
- **generation:** n/a
- **defects:**
  - `:396-399`: `project_regularization_vertex` accepts only the last vertex ("canonical positive chart"). This is a research API with no consumer.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_partition.h | lines 27 | read 1-27

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_partition_internal.h | lines 40 | read 1-40

- **role:** kernel
- **used by:** `event_partition.cpp`, `trimmed_line_branches.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_trace.cpp | lines 337 | read 1-337

- **role:** engagement-control (trace event identity, verified trace construction)
- **used by:** `circle_oracle.cpp`, `segment_oracle.cpp`, bindings
- **state:** live
- **generation:** n/a
- **defects:**
  - `:265-276`, `:290-291`: `build_event_trace` re-verifies, which means another full reconstruction.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/event_trace.h | lines 76 | read 1-76

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/exact_polynomial_sign_2.cpp | lines 47 | read 1-47

- **role:** kernel (sign of a bivariate polynomial at an algebraic point through square-free factorisation)
- **used by:**
  - `segment_fibre.cpp`, `segment_strata.cpp:694-705`
  - `segment_site_parameterization.cpp:3`
  - linked into `_stock_2` (`CMakeLists.txt:444`)
- **state:** live
- **generation:** n/a
- **defects:**
  - `:31-33`: `square_free_factorize_2` runs on every call, even when the same polynomial is reused across branch pairs (`segment_strata.cpp:689-705` rebuilds and factorises per ordered pair). This is part of the bivariate cost centre.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/exact_polynomial_sign_2.h | lines 9 | read 1-9

- **role:** kernel
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/parameter_charts.cpp | lines 792 | read 1-792

- **role:** motion-geometry (frozen charts; rational pullback of segment and full-circle motions against line and circle supports)
- **used by:** `segment_projection.cpp:788`, `circle_projection.cpp:968,1112`, bindings
- **state:** live
- **generation:** defines the motion vocabulary: **segment and full circle only** (`:390-416`). Quality: weak for generation, since there are no partial arcs, helices, ramps or variable-radius motions. Trochoidal arcs, bridges as arcs, and entries cannot be certified.
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - degree bounds (1,2), (2,4), (2,2), (4,4) (`:394-412`): verified against the constructions at `:297-311`, `:324-365` and `:719-770`.

### src/continuous_tea_2/parameter_charts.h | lines 25 | read 1-25

- **role:** motion-geometry
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/partition_certificate.h | lines 374 | read 1-374

- **role:** engagement-control (certificate record types, error hierarchy)
- **used by:** all TEA files, `exact_algebraic_1.h:3`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:159-165`: `ParameterCell2` carries only root ids and a witness, not the exact lower and upper boundaries. A consumer that wants a "max certified advance" must join through `roots` itself.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/result.h | lines 37 | read 1-37

- **role:** engagement-control (ternary verdict, boundary errors)
- **used by:** all
- **state:** live
- **generation:** n/a
- **defects:** none found. The ternary verdict CERTIFIED / CAP_EXCEEDED / UNRESOLVED is the right shape for a gate, but a generator must define what it does on UNRESOLVED.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_fibre.cpp | lines 1630 | read 1-1630

- **role:** engagement-control (exact evaluation of a segment event fibre: branch reconstruction, trim, vertex, overlap replay)
- **used by:** `segment_partition.cpp:521`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:835-838`: `kernel2.solve_2_object()` per fibre per pullback. This is the `Curve_pair_analysis_2` cost centre profiled in `continuous_engagement_cost.md:59-73` (unverified profile). The events along a one-parameter motion are univariate (cost doc :97-101), so the bivariate solve is the formulation problem.
  - `:958-962`: a branch-count mismatch throws `IncompleteSegmentPartitionError`. This is one source of the 15/125 incompletes (inference, unverified).
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_fibre.h | lines 29 | read 1-29

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_oracle.cpp | lines 663 | read 1-663

- **role:** engagement-control (frozen-stock segment audit, swept-prefix π theorem, one-station refuter)
- **used by:**
  - `adaptive/motion_oracle_cache.py:151` (LINK certify)
  - `adaptive/motion_certificate.py:590` (advancing segment)
  - `adaptive/motion_refutation.py:329`
  - `audit_certification_2.cpp:295`
  - `audit_certification_replay_2.cpp:407,444`
  - `audit_unresolved_evidence_2.cpp:457`
- **state:** live (audit lane plus parked adaptive)
- **generation:** links through cleared stock are adequate. Advancing cuts below a π cap are broken: frozen stock refutes them. A π-cap advancing cut from a clear start is excellent and cheap (0.041 s vs >10 min generic, `continuous_engagement.md:1331-1335`, stated).
- **defects:**
  - `:505-513`: frozen-stock station semantics with no forward-semicircle mask (see finding 1).
  - `:527-538`: `violating_parameter` is a cell-interior witness.
  - `:562-566`: verification rebuilds the partition (2×).
  - `:598-602`: the swept-prefix theorem certifies only an exact π cap. The theorem stated at `segment_oracle.h:50-54` actually implies exactness for any cap once the rim is masked to the forward semicircle.
- **unused capability / ideas:** a forward-semicircle-masked `classify_station_cell` would give exact segment TEA at any cap (derivation in finding 1).
- **claims:**
  - `segment_oracle.h:42-45` ("must ... fail loudly when any stratum lacks a resolved disposition"): verified (`:498-502`, `:540-543`).
  - the `motion_refutation.py:7-10` claim that one refuted station is a "complete proof" of motion violation holds only in **frozen-stock** semantics. For an advancing cut it is contradicted by `docs/continuous_engagement.md:1106-1111`, which says frozen stock is "too pessimistic for a slot-opening cutter". This is safe for LINK only.

### src/continuous_tea_2/segment_oracle.h | lines 68 | read 1-68

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** see finding 1
- **claims:**
  - `:50-54` swept-prefix theorem: verified by the derivation in `docs/continuous_engagement.md:1113-1130` and by the code (`segment_oracle.cpp:200-229`).

### src/continuous_tea_2/segment_pair_projection.cpp | lines 414 | read 1-414

- **role:** engagement-control (cross-feature cap and endpoint resultants for segments)
- **used by:** `segment_projection.cpp:1023`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:347-412`: every pair of reach-admitted features × 4 chart pairs × 2 trivariate resultants (`:138-170`). The reach filter admits all line features (`segment_projection.cpp:527-529`), so this is O(n²) in the stock's line count.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_pair_projection.h | lines 19 | read 1-19

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_partition.cpp | lines 881 | read 1-881

- **role:** engagement-control (segment event partition assembly and verification)
- **used by:** `segment_oracle.cpp:481,563`, bindings
- **state:** live
- **generation:** n/a
- **defects:**
  - `:565-595`: incomplete fibre evidence throws, which makes the audit unresolved (counted in the 15/125).
  - `:778-833`: verification is reconstruction.
  - `:835-881`: test mutators.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_partition.h | lines 107 | read 1-107

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_projection.cpp | lines 1055 | read 1-1055

- **role:** engagement-control (segment pullbacks, tangency, cap-crossing, vertex-passage and endpoint projections)
- **used by:** `segment_partition.cpp:419`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:523-529`: the reach filter skips only circle supports and admits every line, even one far from the motion.
  - `:751-754`: throws when every support is filtered, so a motion with nothing in reach is "incomplete" instead of clear.
  - `:764-769`: a guard that cannot fire (fewer records always means fewer pairs).
- **unused capability / ideas:** none
- **claims:**
  - the `PAIR_RESULTANTS_PER_FEATURE_PAIR = 8` constant (`:29-30`): verified as 4 chart pairs × 2 (`segment_pair_projection.cpp:357-408`).

### src/continuous_tea_2/segment_projection.h | lines 19 | read 1-19

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_source.cpp | lines 432 | read 1-432

- **role:** kernel (exact segment and full-circle source injection with canonical bytes)
- **used by:**
  - all oracles
  - `audit_certification_2.cpp:5`
  - `audit_certification_replay_2.cpp:5`
  - `audit_unresolved_evidence_2.cpp:6`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:21-51`, `:161-168`: `exact_binary64` and `lift_binary64` are dead, because `from_binary64` routes through `from_exact` (`:125-131`).
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_source.h | lines 140 | read 1-140

- **role:** kernel
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** `:52` declares `lift_binary64`, which is dead.
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/segment_strata.cpp | lines 1015 | read 1-1015

- **role:** engagement-control (rim branches at a rational station; ordered branch-pair orientation and cap dispositions)
- **used by:**
  - `segment_partition.cpp:463,469`
  - `circle_strata.cpp:400,572`
  - `segment_oracle.cpp:641`
- **state:** live
- **generation:** the per-station exact rim/stock run classification, the core of every verdict. Quality: adequate for correctness, weak for cost.
- **defects:**
  - `:667-747`: O(k²) ordered pairs. For each pair it builds `pair_polynomials` and calls the bivariate `exact_polynomial_sign_at_2` three times (with square-free factorisation each time). At a rational witness the rim parameters are one-root numbers (`identify_rim_root`, `:311-392`), so `exact::sign_mixed_radical` could decide each sign univariately.
  - `:246-248`: a cutter circle coinciding with a stock arc throws.
- **unused capability / ideas:** none
- **claims:**
  - the ≤π handling (`:708-730`: `ccw-major` ⇒ above-cap): verified correct, because the cap chord ratio is ≤4 (cap ≤ π).

### src/continuous_tea_2/segment_strata.h | lines 42 | read 1-42

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/sha256.cpp | lines 163 | read 1-163

- **role:** tool (in-tree SHA-256)
- **used by:**
  - every certificate
  - `audit_*`, `exact_depletion_2.cpp`, `segment_site_*`, `reachable_certificate_encoding_2.cpp`
  - linked into `_stock_2` and another module (`CMakeLists.txt:445,527`)
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/sha256.h | lines 5 | read 1-5

- **role:** tool
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/station_classifier.cpp | lines 192 | read 1-192

- **role:** engagement-control (one station: material runs from a reference point and branch parity; per-run cap)
- **used by:** `segment_oracle.cpp:505,646`, `circle_strata.cpp:575`
- **state:** live
- **generation:** this is the decision that would carry a forward-semicircle mask. Quality: adequate.
- **defects:**
  - `:96-112`: the reference point (cx, cy − r) on the boundary gives UNRESOLVED. Generic witnesses avoid it.
  - `:17-55`: the dormant decoder is retained (stage 6).
- **unused capability / ideas:** add the forward-semicircle mask (finding 1).
- **claims:**
  - "per-run cap" (`:161-185`): verified. This matches the per-run rule stated in `docs/refutation_soundness_crosscheck.md:36-41`.

### src/continuous_tea_2/station_classifier.h | lines 29 | read 1-29

- **role:** engagement-control
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/continuous_tea_2/station_source.cpp | lines 229 | read 1-229

- **role:** kernel (exact station source; attestation projection)
- **used by:** `segment_oracle.cpp:104`, `circle_strata.cpp:157`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:22-54`: the legacy decoder is retained (stage 6).
- **unused capability / ideas:** none
- **claims:**
  - `:153-157` (filtered sign and compare on Epeck::FT): verified (`CGAL::sign`/`compare`).

### src/continuous_tea_2/station_source.h | lines 134 | read 1-134

- **role:** kernel
- **used by:** see `.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - `:16-19` ("its only caller of `compas_cgal::exact::to_canonical`"): verified (`station_source.cpp:91`; grep finds no other caller).

### src/continuous_tea_2/trimmed_line_branches.cpp | lines 840 | read 1-840

- **role:** engagement-control (research: on-support line motion trim branches)
- **used by:** binding `solve_trimmed_line_branches` only; consumer is `tests/adaptive/test_trimmed_line_branches.py`
- **state:** dead for generation (test-only)
- **generation:** n/a. It solves only motion **parallel** to the support (`:546-549`).
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact/canonical.cpp | lines 64 | read 1-64

- **role:** kernel (canonical rational attestation view)
- **used by:** `station_source.cpp:91`
- **state:** live
- **generation:** n/a
- **defects:**
  - `:41`: the tag is `exact-binary64-rational-v1`, while the station attestation uses `exact-rational-v1` (`station_source.cpp:141`). Two tags exist for one value class. Effect is attestation-only; no generation effect.
- **unused capability / ideas:** none
- **claims:** none

### src/exact/canonical.h | lines 49 | read 1-49

- **role:** kernel
- **used by:** `station_source.cpp:4`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - `:37-42` ("As of stage 0 nothing is routed through it yet"): contradicted. `ExactRational2::project` calls it (`station_source.cpp:91`).

### src/exact/errors.h | lines 42 | read 1-42

- **role:** kernel (named exact errors)
- **used by:** `exact/*.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - `:33-36` (`AttestationByteDriftError` "has no raiser yet"): verified (grep finds no raiser outside `src/exact`).
  - `CrossRootExtensionError` is raised only in `one_root.cpp:26`: verified.

### src/exact/one_root.cpp | lines 107 | read 1-107

- **role:** kernel (same-root arithmetic guard; exact sign of a + b√α + c√β + d√(αβ))
- **used by:**
  - `sign_mixed_radical`: `exact_cap_certificate_2.h:108,137` and `audit_exact_station_2.cpp:245`
  - `same_root_*`: tests only (`exact_one_root_gate`)
- **state:** live
- **generation:** the exact univariate cap/orientation predicate for one-root rim points. Quality: excellent.
- **defects:** none found
- **unused capability / ideas:** use it in `segment_strata.cpp:667-747` in place of the bivariate kernel. This is the cheapest lever on the continuous oracle's cost.
- **claims:**
  - `:56-64` (fold zero roots before building a OneRoot): verified as a correct guard.

### src/exact/one_root.h | lines 94 | read 1-94

- **role:** kernel
- **used by:** `exact_cap_certificate_2.h:26`, `audit_exact_station_2.cpp`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - `:28-33` (`same_root_add/multiply` have no pipeline caller; `sign_mixed_radical` is on the deciding path): verified (grep).
  - `:65-67` ("the exact cap predicate reduces to this"): verified as consistent with `exact_cap_certificate_2.h` usage.

### src/exact/rational.cpp | lines 21 | read 1-21

- **role:** kernel
- **used by:** no caller found for `exact::from_binary64` (grep)
- **state:** dead
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

### src/exact/rational.h | lines 37 | read 1-37

- **role:** kernel (the `Rational = Epeck::FT` alias)
- **used by:** `station_source.h`, `one_root.h`, `canonical.h`
- **state:** live (alias); the function is dead
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:**
  - `:17-22` ("nothing is routed through it yet"): verified (no callers).

### src/exact_algebraic_1.cpp | lines 1279 | read 1-1279

- **role:** kernel (univariate root isolation, global root and cell partition over [0,1], dyadic witnesses, one-root physical/conjugate classification, backend pinning)
- **used by:**
  - all TEA oracles (`partition_integer_projections`)
  - `continuous_tea_exact_core`
  - `_stock_2` (`CMakeLists.txt:443`) and another module (`:525`)
- **state:** live
- **generation:** computes every exact event root of a motion in [0,1] with strict isolating intervals (`:507-590`). This is the data a certified maximal advance needs: the lower root of the first violating cell, or its `interval_low` as a certified rational. Nothing extracts it.
- **defects:**
  - `:20-34`: GMP is forbidden by `#error`. This is an owned trade (cost doc :79-84).
  - `:686-723`: `simplest_dyadic_between` searches linearly over numerators at each power of two. For close roots this is O(gap⁻¹)-ish work per cell (inference, not measured).
- **unused capability / ideas:** none beyond the advance extraction above
- **claims:**
  - backend evidence (`:832-851`): verified against the static_assert at `:834-837`.

### src/exact_algebraic_1.h | lines 40 | read 1-40

- **role:** kernel
- **used by:** see `.cpp`, `segment_site_parameterization.h`
- **state:** live
- **generation:** n/a
- **defects:** none found
- **unused capability / ideas:** none
- **claims:** none

## Answers to the slice questions

**(a) Can it drive generation?** Only as an accept/reject oracle on a proposed motion. `adaptive/` already drives generation this way. It enumerates a finite dyadic lattice of progress × guide radius × phase (`adaptive/candidates.py:1129-1200`), refutes links with exact station probes, certifies the survivors (`adaptive/transaction.py:1091-1106`), and certifies circles directly (`:1136-1150`). It accepts the first feasible candidate in policy order (`adaptive/generator.py:960-1035`). No API returns a maximal advance or radius. The roots and cells needed for a maximal advance are internal (`exact_algebraic_1.cpp:1240-1277`). The exposed `violating_parameter` is a cell-interior witness and is not bound to Python (`segment_oracle.cpp:527-538`, `continuous_tea_2.cpp:1596-1599`). The API a generator would call:
- `audit_full_circle_tea_event_exact(stock, cx, cy, phase_dx, phase_dy, cw, r, cap_ratio) -> (verdict, trace)` per candidate loop;
- `segment_station_cap_exceeded_exact` as a refuter for links;
- `audit_swept_prefix_segment_tea_exact` for π-cap advancing cuts.

**(b) Correctness state.** The false certificate and the swept-annulus guard are in `engagement_2.cpp`, not here (commits 2e2422b5, f227dab2). The shared-root precondition is enforced in this slice by `exact/one_root.cpp:11-29` (raise on cross-root) and by folding zero roots (`:56-64`). The release-build fix 148a49ad is in `engagement_2.cpp`. Soundness of this oracle in the certified direction holds by construction: frozen stock ⊇ physical stock, per-run cap, and exact predicates throughout (verified by reading). Completeness does not hold: 15/125 link audits incomplete plus a native abort (stated, unverified). Advancing cuts below π are refuted by construction. Trustworthy enough to gate a **loop** placement after the 33 probes accept, provided UNRESOLVED is treated as reject. Not usable to gate straight cutting moves below a π cap.

**(c) Cost.**
- Full-circle audit: 0.31–0.72 s warm at 8–16 boundary features, growing with feature count (measured above).
- 33 frozen probes: 3–5 ms in total (measured).
- Segment audit: about 15 s at 20 vertices, superlinear (stated, `continuous_engagement_cost.md:13-35`).
- About 128–301 ms for small certificates (`continuous_engagement.md:327-332`, stated).

That is too slow for every rung of a ladder scan. One audit per emitted loop is affordable on the gate pockets (a few hundred loops gives minutes), but not on production pockets.

**(d) Certified advances for the tree-walk plan?** Not in Stage A, because it moves no gate number (plan:43). It belongs in Stage C behind a measurement. What it would take:
1. Keep the 33 probes as the refuter.
2. Call `audit_full_circle_tea_event_exact` once on the accepted rung; on `cap_exceeded` or `unresolved`, step down the ladder.
3. For bridges and run-ins, build the forward-semicircle mask in `station_classifier.cpp` so that straight cuts get exact TEA at 120°/40°.
4. To cut the cost, replace the bivariate pair signs in `segment_strata.cpp:667-747` with `exact::sign_mixed_radical`, and add a locality filter to `extract_boundary_records`.

A true maximal advance additionally needs a new query that returns the first non-certified cell's lower root.
