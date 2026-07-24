# Engagement/stock performance: the O(n²) was API misuse, and the fix is a persistent arrangement

**BLUF.** The SP1 stock/engagement kernel's "prohibitively O(n²)" audit cost is dominated by
**CGAL API misuse, not the geometry, the exact kernel, or an inherent limit.** Two
behavior-preserving fixes already landed (aggregated disk-chain union → **7×** end-to-end;
per-query stock-copy elimination → **11%**). The residual cost is the
`General_polygon_set_2::intersection` **full-overlay-per-query** pattern; the principled fix is a
**persistent `Arrangement_2` with `zone`-based local queries**. This exploration compiled and
measured the right primitives and **validated that redesign as feasible** — it is not yet built.

Prototypes referenced below live beside this file (`min.cpp`, `pl_bench.cpp`, `stock_arr.cpp`);
`compile.sh` has the exact build line.

## Evidence

| Claim | Measurement | Source |
|---|---|---|
| Subtract was O(N²) *per op* from an incremental-union anti-pattern | per-circle `subtract_arc_sweep` **3612 ms → 168 ms** (~21×) | `030322f`, `decompose.py` |
| …fixing it is a 7× end-to-end win, behavior-preserving | full kite audit **>600 s → 87 s**; full suite **34 s → 9.5 s**; 165 tests unchanged | `030322f` |
| Per-query full-stock **copy** is real but minor | square audit **204 s → 183 s** (~11%) | `e559ed7` |
| The residual is the disk-vs-whole-stock **overlay**, ~3600 queries/pocket | est. ~137 s of square's 183 s | `engagement_2.cpp` read |
| Coarsening the certify grid does **nothing** (cost isn't the stations) | arc step 0.05→0.15: 20→7 stations, per-circle **3612 ms → 3668 ms** (0×) | `decompose.py` |
| Loosening chain density is an **accuracy** regression, not free | 1e-4→1e-3: kite 87 s→27 s **but breaks 2 calibrated `>pi`/seam tests** | `test_stock.py:140,174` |
| Epeck is **correct** for the one-root arcs; Epick would be non-robust | `Gps_circle_segment_traits_2<K>` = `Gps_traits_2<Arr_circle_segment_traits_2<K,Filtered>>`; one-root coords need exact rational `FT` | header read |
| Point-location for the mutate-heavy interleave = **`walk_along_line`** | (insert; 20 locate)×250: walk **287+460 ms**; trapezoid_ric **6047+35 ms** (insert tax); landmarks **crashes** | `pl_bench.cpp` |
| Persistent circle-segment arrangement + face `material` bool **works** | square−disk, 4/4 point-in queries correct | `stock_arr.cpp` |
| `zone(cutter)` reads engagement **locally, no overlay** | cutter straddling the disk: **4 faces (3 material), 2 edge crossings** | `stock_arr.cpp` |

## What landed vs what was ruled out

- **Landed (committed, behavior-preserving, 165 green):** aggregated union `030322f`; per-query
  copy elimination `e559ed7`.
- **Ruled out with evidence:** coarser station grids (0× — wrong driver); chain-density loosening
  (accuracy regression on calibrated `>pi`/seam tests — will not retune reference tests); Epick
  kernel (non-robust for the one-root circle coordinates — Epeck is correct).

## The redesign (validated feasible; not built)

One persistent `Arrangement_2<Arr_circle_segment_traits_2<Epeck>, Arr_extended_dcel<…, bool material>>`,
mutated incrementally instead of overlaid per op:

- **Deplete:** `CGAL::insert(arr, sweep_curve, pl)` (zone-local), mark swept-interior faces
  `material=false`.
- **Query (the win):** split the cutter circle with `make_x_monotone_2`, `CGAL::zone(arr, arc, oi, pl)`,
  read `face->data()` for each crossed cell → engaged arcs. Local (O(crossings)), not O(stock).
- **Point-location:** `Arr_walk_along_line_point_location` — measured best for the insert-heavy
  interleave (no aux-structure rebuild tax).

```cpp
using Traits = CGAL::Arr_circle_segment_traits_2<Epeck>;              // Epeck: one-root coords need exact FT
using Dcel   = CGAL::Arr_face_extended_dcel<Traits, bool>;           // face bool = material
using Arr    = CGAL::Arrangement_2<Traits, Dcel>;
using PL     = CGAL::Arr_walk_along_line_point_location<Arr>;        // measured winner for mutate-heavy
// deplete: CGAL::insert(arr, Curve(K::Circle_2(c, r2)), pl); then mark swept faces material=false
// query : make_x_monotone_2(cutter, back_inserter(xs)); for each xc: CGAL::zone(arr, xc, back_inserter(cells), pl);
```

## Buildout plan (the careful work ahead — certified path)

1. **Exact TEA from the zone** — angular spans of material sub-arcs between the zone's edge
   crossings (arithmetic on the exact crossing points; the current `harvest` logic transplants).
2. **`Arr_observer`** to maintain `material` across insertions efficiently (`after_split_face`
   propagate + flip swept side) — replaces the prototype's explicit re-marking.
3. **Integrate into `engagement_at`** behind the current interface; keep the disk-chain sweep for
   now (the exact annular sector's `guide_r ± r` offset is irrational → a separate one-root
   construction, deferred).
4. **TDD against the current results** — the redesign must reproduce the existing 165 tests'
   TEA/verdicts exactly (it is exact, not approximate).
5. **Review** in the certified path before it lands.

!!! success "CONFIRMED (measured this session)"
    Aggregated-union fix (7×) and copy-elimination (11%) are committed and behavior-preserving.
    Point-location = `walk_along_line`. Persistent-arrangement substrate + `zone` query are
    validated by compiled prototypes.

!!! warning "UNDER BUILD (not yet done)"
    Exact zone→TEA span computation, observer-based marking, `engagement_at` integration, TDD
    against current results, and review. Until these land, the certified path still uses the
    (now 7× faster) `General_polygon_set_2` code.
