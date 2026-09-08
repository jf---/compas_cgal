# Plan — local zone engagement query (land the validated arrangement redesign)

**Status:** ready to execute. Every piece is validated by compiled prototypes in
`docs/dev/arrangement-redesign/` (see `findings.md`). This is a certified-path change, so it is
executed **add-alongside + TDD-gated + swap**, never in-place.

## Goal

Replace the O(stock)-per-query overlay in `engagement_at` (`region.intersection(stock.set())`) with
a **local zone query** of the cutter in the Gps's OWN arrangement (`stock.set().arrangement()`,
faces carry `contained()` = material). Same exact result, read locally — no whole-stock overlay.
The fast depletion (aggregated disk-chain union, commit `030322f`) and the exact certificate logic
(run-assembly + pessimistic cap decision) are **kept verbatim**.

## Why (measured)

- The subtract bottleneck is already fixed (`030322f`, 7× ; `e559ed7`, 11%).
- The residual is the per-query overlay (~137 s of square's 183 s). The overlay processes the whole
  stock arrangement each of ~3,600 queries/pocket.
- Prototype `harvest.cpp` proves the local zone query returns the EXACT engagement (full immersion
  → 2π, edge-straddle → π, in-cleared → 0), reusing the existing exact run-assembly.

## Tasks

1. **`finish_engagement(arcs, cx, cy, r_sq, cap_chord_ratio, gap_close_ratio) -> EngagementSample`** —
   extract the shared tail (assembly + reporting + pessimistic cap decision, currently
   `engagement_2.cpp` lines ~444–516) into a helper in the anonymous namespace. Certificate logic
   is copied VERBATIM (no change).
2. **`engaged_arcs_zone(stock, cx, cy, tool_radius, arcs)`** — the local harvest (port
   `harvest.cpp`): `const_cast` the logically-const stock, build `Arr_walk_along_line_point_location`
   over `stock.set().arrangement()`, `make_x_monotone_2` the cutter circle, drive
   `Arrangement_zone_2` with a `found_subcurve` visitor collecting sub-arcs whose face is
   `contained()`, extract `Arc{ccw_start, ccw_end, span}` exactly as the current harvest does.
3. **`engagement_at_zone(...)`** — new entry point beside the UNTOUCHED `engagement_at`: validate,
   `engaged_arcs_zone` → `arcs`, `finish_engagement`. Declare in `engagement_2.h`, bind in
   `register_engagement`.
4. **Comparison test** — `tests/test_engagement_zone.py`: over many random cutter positions on the
   slot/oracle witness stocks, assert `engagement_at_zone == engagement_at` (total_tea, max_run_tea,
   cap_exceeded) EXACTLY.
5. **Swap** (only after Task 4 + the full suite + the differential oracle pass): make
   `engagement_at` use `engaged_arcs_zone`; remove the overlay harvest; keep `finish_engagement`.
   Re-run the full suite.
6. **Measure + record** — query speedup on a benchmark pocket (expect the overlay's ~137 s to
   collapse); update `docs/dev/arrangement-redesign/findings.md` status.

## Exit gates

- (a) `engagement_at_zone` reproduces `engagement_at` EXACTLY on the full 165-test suite AND the
  differential Monte-Carlo oracle (`test_engagement_oracle.py`) AND the new random-cutter comparison
  test. If ANY case diverges → STOP, write
  `docs/superpowers/state/engagement-zone-divergence.md` with the failing geometry; do not tune.
- (b) measured query speedup is a real, recorded number on a benchmark pocket (a slow result is a
  result to report, not to hide).
- (c) after the swap, the full suite is green and the certificate is bit-for-bit unchanged.

## Decoupling

The prototypes (`docs/dev/arrangement-redesign/*.cpp`, buildable via `compile.sh`) validate every
CGAL primitive standalone against the `cgal-dev` env — independent of the module. The module change
is validated purely against the CURRENT `engagement_at` (the oracle), so it needs nothing upstream.

## Risk

Certified path: a harvest that produces different arcs breaks the certificate silently for uncovered
geometries. Mitigation = add-alongside (the 165 tests stay on the untouched path) + the differential
oracle (random ground truth) + the random-cutter comparison. The `const_cast` is justified (zone
with a non-inserting visitor is non-mutating; single-threaded audit) and encapsulated.
