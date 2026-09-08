# Certified Adaptive Clearing — Design

**Status:** Draft for review
**Date:** 2026-07-23
**Repository:** `compas_cgal_prs` (branch off `jf/toolpath-redesign`)
**Programme:** Volumetric roughing (2.5D adaptive clearing)

## 1. Decision

Build a tool-engagement-regulated 2.5D clearing engine inside `compas_cgal_prs`,
CGAL-only, on three pillars:

1. an **exact in-process stock model** as a CGAL 2D Arrangement of circular
   arcs and segments with per-face material state;
2. a **propose / certify-exact generator**: contour-parallel offset-loop passes
   proposed cheaply, every emitted motion certified against the exact stock
   before commit; certified trochoidal relief (existing substrate) wherever
   regulation fails;
3. the existing **certified trochoid toolpath module** reused unchanged for
   relief, entry slotting, leads/links machinery, and serialization.

Out of scope for this programme: hotcam / signed-heat fields (separate,
experimental), the tri-dexel simulator (tertiary; future cross-check via a
toolpath interchange artifact), robotic execution coupling, Clipper2 or any
non-CGAL geometry dependency.

## 2. v1 success criterion

On the benchmark pocket suite, the engine emits toolpaths with a certificate:

- **Hard cap (exact):** no cut motion ever exceeds `tea_cap` engagement.
  Decided by exact predicates (chord comparisons on one-root points), enforced
  at generation time — violations are impossible by construction, not sampled.
- **Band (statistical):** ≥ 95 % of engaged cut length lies within
  `tea_target ± tea_band`, reported per pocket in the `EngagementReport`.
- **Zero gouge:** inherited from the existing exact-clearance certification.
- **Residual (exact):** remaining-stock area beyond the reachable set < ε,
  extracted from the arrangement.

This combination — machine-checkable engagement and coverage claims — is the
differentiator over VoluMill / iMachining / Fusion Adaptive, which regulate
engagement but publish no certificates.

## 3. Stock model (SP1)

`CGAL::Arrangement_2<Arr_circle_segment_traits_2<Epeck>,
Arr_face_extended_dcel<…, MaterialState>>`.

- **Faces carry material state** (`material` / `cleared`), maintained through
  splits and merges by an `Arr_observer`.
- **Commit = subtraction of the sweep region** of each committed motion.
  Tool disks (plunge/dwell) are exact: input doubles are rationals, so the
  disk's center and squared radius are rational. Swept regions are not:
  a capsule's offset side-lines carry an `r·√(dx²+dy²)` term and an arc
  sweep's tool-offset radii have irrational squares — both outside the
  circle-segment traits. Line and arc sweeps are therefore represented as
  **conservative under-covering disk chains** (exact tool disks at a spacing
  whose sagitta bound is a named small fraction of the tool radius,
  dominated by `radial_clearance`). Under-covering
  over-estimates both remaining material and TEA, so the hard-cap and
  residual certificates remain sound; only efficiency is conceded.
  Mutation bookkeeping (face containment through overlapping subtractions) is
  owned by the regularized boolean layer of the Arrangements package
  (`General_polygon_set_2` over the same circle-segment arrangement); zone
  and location queries run directly on its underlying `Arrangement_2`.
  Observer-maintained incremental insertion is the SP1 performance fallback
  if the wall-clock gate fails, not the v1 baseline.
- **Certificate = zone query, read-only.** Engagement of a candidate cutter
  circle is `CGAL::zone()` of the circle: traversed cells partition it into
  sub-arcs; sub-arcs in `material` faces are the engagement intervals, with
  exact one-root endpoints. Candidates never mutate the arrangement — the
  propose/certify architecture maps directly onto the package's read/write
  split.
- **Near-tangency hygiene:** consecutive sweep regions are inserted slightly
  fat — extended by the existing `radial_clearance` certification margin
  (material we are entitled to remove) — so unions strictly overlap and the
  arrangement never processes grazing tangencies.
- **Hygiene as performance lever:** edges interior to fully cleared regions are
  removed (`remove_edge`), keeping the arrangement proportional to the material
  frontier, not the toolpath history.
- `General_polygon_set_2` is used only at the rim: pocket-with-islands
  initialization and residual-stock extraction.
- Point location: `Arr_landmarks_point_location`.

### TEA certificate along a motion

TEA is continuous along a move (new contacts are born at zero angular measure)
but can grow with √-steepness at grazing contact; there is no global Lipschitz
constant. The certificate is therefore **exact at stations plus conservative
refinement**: stations are inserted while the margin `tea_cap − measured`
cannot be closed by the geometric growth bound over the remaining spacing;
refinement to a named floor; an unclosable margin is a violation and triggers
adaptation. This is the same refinement pattern as `certify_bridges` in the
existing module.

## 4. Generator (SP2)

**Pass family: contour-parallel offset loops** from the existing
`straight_skeleton_2` interior offsets, walked outside-in per connected
region.

Loop per pass:

1. **Propose** the next offset loop at the TEA-derived radial step
   `ae = R · (1 − cos tea_target)`, adjusted per-station.
2. **Certify** exact TEA along the proposed pass against the current stock.
3. **Adapt** on violation, locally and in order: shrink the radial step →
   insert certified trochoidal relief (existing substrate: corners, slots,
   first channels) → split the pass. Adaptation exhaustion raises
   `EngagementCertificationError` with diagnostics; nothing uncertified is
   ever emitted.
4. **Commit**: append to the `ToolpathPrimitive` stream and insert the sweep
   boundary into the arrangement.

- **Entry** per connected region at the skeleton's maximum-clearance vertex,
  using the existing certified trochoidal slot primitive under its own TEA
  cap.
- **Links/leads/ordering/serialization** are refactored out of `toolpath.cpp`
  into a shared `toolpath_common` unit and reused; the engine emits the same
  `ToolpathResult` (viewer example works unchanged). Nothing existing is
  replaced.
- **SP2b — spiralization:** morph adjacent certified loops into a continuous
  spiral (Held/Spielberger-style) to remove loop seams. Separate follow-on
  with its own gate; v1 ships with loops + certified tangent links.

## 5. Parameter model

TEA is the primary control: `tea_target`, `tea_band`, `tea_cap` (radians,
validated `0 < target ≤ cap ≤ π`). The radial step is derived; `stepover` is
accepted only as an optional additional cap for compatibility. Tool model v1:
flat-end cylindrical, single cutting depth, closed pockets with islands
(`holes`). Climb/conventional via existing `climb`.

## 6. Components

| Unit | Responsibility |
| --- | --- |
| `src/stock_2.cpp/.h` | `Stock2`: arrangement, material faces, incremental sweep insertion, residual extraction |
| `src/engagement_2.cpp/.h` | zone-based TEA query, station certificates, refinement bound |
| `src/clearing_2.cpp/.h` | propose/certify/adapt/commit loop, offset-loop proposal, relief insertion |
| `src/toolpath_common.cpp/.h` | shared: `TrochoidArc`, primitives, leads/links, ordering, serialization (extracted from `toolpath.cpp`) |
| `src/compas_cgal/stock.py` | typed `Stock2` wrapper |
| `src/compas_cgal/clearing.py` | `adaptive_clearing_toolpath(...) -> ToolpathResult`, `EngagementReport` dataclass |
| `src/compas_cgal/feeds.py` (SP3) | chip-thinning feed per op from `EngagementReport`, pure Python |

`EngagementReport`: per-operation TEA min/mean/max, band-compliance fraction by
cut length, hard-cap certificate flag, residual area, adaptation counts.

## 7. Error model

Named exceptions, fail loud, no silent degradation: `InvalidPolygonError`
(existing), `EngagementCertificationError` (adaptation exhausted),
`StockTopologyError` (arrangement invariant violated). `max_passes`-style
truncation warns with counts (existing pattern). All tolerances are named
constants with unit and derivation comments; construction guards are never
geometric decisions (repo CLAUDE.md discipline).

## 8. Testing & gates

- **TDD throughout**; hypothesis property tests over the existing random
  polygon strategy plus dumbbell family and island pockets.
- **In-repo semantic oracle:** dense Monte-Carlo TEA sampling (independent
  math, deliberately slow) differential-tested against the exact kernel.
- **SP1 exit gate:** (a) TEA kernel agrees with Monte-Carlo oracle within
  sampling bounds on the benchmark suite; (b) **engagement audit** of the
  existing trochoid generator produces a truthful per-op TEA report (this
  baseline motivates SP2 and is the number SP2 must beat); (c) measured
  wall-clock ≤ 10 s per benchmark pocket for commit+certify throughput —
  a gate, not a guess; if missed, the fallback stays within CGAL (coarser
  exact station grids), never an approximate clipper.
- **SP2 exit gate:** v1 success criterion (§2) on the full suite, plus: no
  regression of the existing 45 toolpath tests.
- **SP3 exit gate:** feed post-processor validated against `EngagementReport`;
  benchmark report (BLUF-first) generated by the suite runner.
- Benchmark suite: SQUARE, IRREGULAR, L_SHAPE, STAR, KITE, dumbbell(2.4/1.8/1.2),
  island pocket, one ~100-edge random polygon (seeded).

## 9. Sequencing

SP1 (stock + TEA kernel + audit) → SP2 (regulated generator) → SP2b
(spiralization) → SP3 (feeds + benchmark report). Each sub-project gets its
own implementation plan; SP1's plan is written first.

## 10. Non-goals / future

Open pockets, multi-depth, rest-machining chains, G2/jerk-limited blending,
material/feeds database, tri-dexel cross-validation (via a schema-validated,
`BuildIdentity`-stamped toolpath JSON interchange artifact), robotic execution
coupling (tesseract), true-MA guide curves (`Segment_Delaunay_graph_2`).
