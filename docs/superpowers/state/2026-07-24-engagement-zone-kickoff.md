# Engagement zone-query kickoff — land the validated arrangement redesign

I'm starting the **local zone engagement query** phase of the certified adaptive-clearing kernel in
**`/Users/jelle/Code/CADCAM/compas_cgal_prs`**, branch **`jf/toolpath-redesign`**. Commit
**`1bb9a78`** (`docs(perf): validated the full zone-visitor harvest (exact TEA) end to end`) is the
last hard checkpoint. The fix is fully de-risked by compiled prototypes; this phase LANDS it in the
certified module. Everything below is readable cold — no prior conversation needed.

## 1. What I want you to do

Execute `docs/plans/2026-07-24-engagement-zone-query.md` with `superpowers:test-driven-development`
(the comparison test is written and failing BEFORE the new code path). This is a **certified-path**
change: it is **add-alongside → TDD-gate → swap**, NEVER edited in place. Parent-session rules in
`~/.claude/CLAUDE.md` apply (ETH/DeepMind bar; the exact-kernel discipline is the gate). Take the
ETH route always: the certificate logic is copied verbatim, never re-derived; the new harvest must
reproduce the old EXACTLY before anything is swapped.

## 2. Read first, in dependency order

1. `docs/plans/2026-07-24-engagement-zone-query.md` — THIS PHASE's plan: tasks, exit gates, risk.
2. `docs/dev/arrangement-redesign/findings.md` — the exploration that found + proved the approach,
   with the measured evidence table and the confirmed CGAL interfaces.
3. `docs/dev/arrangement-redesign/harvest.cpp` — the VALIDATED harvest (exact TEA: full immersion
   2π, edge-straddle π, in-cleared 0). This is the code you port. Also `gps_zone.cpp` (zone on the
   real Gps path), `pl_bench.cpp` (why `walk_along_line`, measured). Build them:
   `bash docs/dev/arrangement-redesign/compile.sh` (needs the `cgal-dev` env; see §9).
4. `src/engagement_2.cpp` — the certifier. Current `engagement_at`: harvest at lines ~393–443
   (overlay + rim-arc extraction), shared tail at ~444–516 (run-assembly + pessimistic cap decision).
5. `CLAUDE.md` — "Exact-Kernel Discipline" (no epsilon/deflation in any decision path; deciding vs
   reporting split; doubles cross into FT once at documented boundaries). Full doctrine:
   `docs/exactness.md` + `docs/dev/exact-cgal-idioms.md`.
6. `~/.claude/projects/-Users-jelle-Code-CADCAM/memory/cgal-rtfm-explore.md` — FEEDBACK memory:
   for CGAL work, read the actual headers + compile-and-measure; do not pattern-match from memory.

## 3. Boundary contract (what you consume and must reproduce)

```python
# compas_cgal._stock_2  (built extension; sources src/stock_2.*, src/engagement_2.*)
engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio, gap_close_ratio=0.0)
    -> (total_tea, max_run_tea, exceeded)     # THE ORACLE this phase must match EXACTLY
# Stock2.set().arrangement() -> the Gps's internal Arrangement_2 (faces carry contained()=material)
```
```cpp
// engagement_2.cpp anonymous namespace
struct Arc { GpsPoint ccw_start; GpsPoint ccw_end; double span; };  // one CCW-normalized rim arc
// types live in stock_2.h: Epeck, GpsTraits=Gps_circle_segment_traits_2<Epeck>,
//   Gps=General_polygon_set_2<GpsTraits>, GpsPoint (one-root), GpsXCurve, EPoint=Epeck::Point_2.
```
The zone harvest fills the SAME `std::vector<Arc>`; the shared tail (`finish_engagement`) turns arcs
into the identical `(total_tea, max_run_tea, exceeded)`.

## 4. The validated approach (proven; port it)

- **Query:** `stock.set().arrangement()` is a real `Arrangement_2` whose `Gps_default_dcel` faces
  carry `contained()` (= material). `make_x_monotone_2` the cutter circle, drive
  `CGAL::Arrangement_zone_2<Arr, Visitor>` with `init(xc, pl)` + `compute_zone()`; the visitor's
  `found_subcurve(cv, face, …)` collects `cv` when `face->contained()` — those ARE the engaged rim
  sub-arcs. Extract `Arc{s,t,span}` exactly as the current harvest (source/target, `if s==t continue`,
  swap if `orientation()==CLOCKWISE`, `span=atan2Δ`, `+2π if ≤0`).
- **Point-location:** `Arr_walk_along_line_point_location<Arr>` — measured winner for the
  insert-heavy mutate/query interleave (`pl_bench.cpp`: trapezoid_ric's insert tax buries it;
  landmarks CRASHES on dynamic insert). Do NOT use trapezoid_ric here.
- **`const`:** `engagement_at` holds `const Stock2&`; `Gps::arrangement()`/`Arrangement_zone_2` want
  a non-const handle. `const_cast<Gps&>(stock.set()).arrangement()` — justified (zone with a
  non-inserting visitor is non-mutating; single-threaded audit), encapsulated in the zone helper.

## 5. Deliverables (plan Tasks 1–6)

1. `src/engagement_2.cpp` — MODIFY: extract `finish_engagement(arcs, cx, cy, r_sq, cap_chord_ratio,
   gap_close_ratio)` (assembly VERBATIM); add `engaged_arcs_zone(...)` + the `found_subcurve` visitor.
2. `src/engagement_2.cpp` + `src/engagement_2.h` — add `engagement_at_zone(...)` beside the UNTOUCHED
   `engagement_at`; bind it in `register_engagement`.
3. `tests/test_engagement_zone.py` — CREATE: random-cutter comparison `engagement_at_zone ==
   engagement_at` on the slot/oracle witness stocks (exact match, all three outputs).
4. After gates pass: SWAP `engagement_at` onto the zone harvest, delete the overlay harvest, keep
   `finish_engagement`; re-run the suite.
5. Measure + record the query speedup; update `docs/dev/arrangement-redesign/findings.md` status.

## 6. Exit gates (from the plan — quoted)

- "(a) `engagement_at_zone` reproduces `engagement_at` EXACTLY on the full 165-test suite AND the
  differential Monte-Carlo oracle (`test_engagement_oracle.py`) AND the new random-cutter comparison
  test. If ANY case diverges → STOP, write `docs/superpowers/state/engagement-zone-divergence.md`
  with the failing geometry; do not tune."
- "(b) measured query speedup is a real, recorded number on a benchmark pocket."
- "(c) after the swap, the full suite is green and the certificate is bit-for-bit unchanged."

## 7. Decoupling rule

The prototypes validate every CGAL primitive standalone (`compile.sh`, `cgal-dev` env) — independent
of the module. The module change is validated ONLY against the current `engagement_at` (the oracle),
so it needs nothing upstream. Never modify the 45 `tests/test_toolpath.py` tests (frozen contract).

## 8. Load-bearing conventions

- Exact-Kernel Discipline: no epsilon/deflation/ulp in ANY decision path; `CGAL::sign`/`CGAL::compare`
  at the number-type level; the certificate logic (`finish_engagement`) is copied VERBATIM — the zone
  change only alters HOW arcs are produced, never the decision. `to_double`/`atan2` stay reporting-only.
- CGAL API facts learned the hard way (verify against headers, don't trust memory): the point-location
  result is a `boost::variant` (not `std::variant`); `Arr_circle_segment_traits_2` has NO
  `construct_curve_2_object` (build `Curve_2` directly); the linear X-curve ctor takes `Epeck::Point_2`,
  the arc ctor takes the one-root `GpsPoint`; `Arrangement_zone_2` ctor is `(arr, &visitor)` then
  `init(cv, pl)` + `compute_zone()`; visitor callbacks return `std::pair<Halfedge_handle,bool>` (model
  `CGAL/Arrangement_2/Arr_compute_zone_visitor.h`).
- Commits: author AND committer `Jelle Feringa <jelleferinga@gmail.com>`, concise conventional
  messages, NO attribution lines. Never commit without standing instruction (plan tasks carry commit
  steps = standing). `docs/examples/example_isolines.py` is a pre-existing user modification: never
  stage, never revert. Linear history; never rewrite committed history. Ruff before Python commits.

## 9. Environment booby traps

- Python: `PY=/Users/jelle/mambaforge/envs/cgal-dev/bin/python` (py312). ALWAYS run with
  `PYTHONPATH=src` from this repo, or you import the CLEAN sibling checkout instead.
- Rebuild the extension (no pip): `cmake --build build/cp312-abi3-macosx_15_0_arm64 --target _stock_2
  -j 8`. `src/compas_cgal/*.so` are symlinks into `build/` (gitignored).
- Full suite: `PYTHONPATH=src $PY -m pytest tests/ -n auto -q` (currently **165 passed ~9.5 s** — the
  aggregated-join fix made it fast; the earlier >120 s/pocket audit cost is gone).
- Standalone CGAL prototypes: env is **CGAL 5.6.1**, header-only but links GMP+MPFR. Compile line:
  `c++ -std=c++20 -O2 -I$ENV/include x.cpp -o x -L$ENV/lib -lgmp -lmpfr` with
  `ENV=/Users/jelle/mambaforge/envs/cgal-dev`; run with `DYLD_LIBRARY_PATH=$ENV/lib ./x`.
- The user's shell aliases `cp` to `cp -i`; scripted copies stall on the interactive prompt — use
  `/bin/cp -f`.
- The 10 pytest warnings are pre-existing compas `LooseVersion` deprecations — not ours.

## 10. Suggested file structure

```
compas_cgal_prs/
├── src/engagement_2.cpp     # MODIFY: + finish_engagement, + engaged_arcs_zone + visitor,
│                            #         + engagement_at_zone, + binding; engagement_at UNTOUCHED until swap
├── src/engagement_2.h       # MODIFY: declare engagement_at_zone
├── tests/test_engagement_zone.py   # CREATE: random-cutter comparison vs engagement_at (the oracle)
└── docs/dev/arrangement-redesign/  # reference: findings.md + harvest.cpp (the port source)
```

## 11. Risk + budget

Certified path: a harvest that yields different arcs breaks the certificate silently for uncovered
geometries. Mitigation is structural — add-alongside (the 165 tests stay on the untouched
`engagement_at`), plus the differential Monte-Carlo oracle (random ground truth), plus the new
random-cutter comparison. Do NOT swap until all three agree exactly. If they diverge, that is DATA:
write the divergence doc, do not tune geometry or widen tolerances.

## 12. Start

```bash
cd /Users/jelle/Code/CADCAM/compas_cgal_prs && git log --oneline -3 && git status -s \
  && PYTHONPATH=src /Users/jelle/mambaforge/envs/cgal-dev/bin/python -m pytest tests/ -n auto -q 2>&1 | tail -3 \
  && bash docs/dev/arrangement-redesign/compile.sh 2>&1 | tail -8
```

Then read §2 in order and execute the plan.

End of prompt.
