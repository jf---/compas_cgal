# Held Paper Figure Reproduction Plan

> **status: in progress** - The bounded visual-comparison study is committed,
> but the completion audit found that it does not yet reproduce continuous
> Figure 5 toolpaths or the Figure 6 three-algorithm comparison. Task 5 is active.

**Goal:** Reproduce the visible Figure 5 toolpath structure and the Figure 6
comparison protocol and axes, with direct comparisons to publisher evidence.

**Active criterion:** for Task 5, independently traverse the reconstructed
Figure 5 medial axis, place the standard 80-degree machining circles, and emit
one continuous planar circle-and-offset-transition stream. Publisher evidence
is comparison-only and may not choose identities, ordering, or geometry.
Characterization beyond Figures 5 and 6, postprocessing, machine validation,
and unrelated benchmark work do not advance this criterion.

**User-directed execution order (September 6): toolpath first.** Land and plot
one publisher-independent continuous path on the existing Figure 5 input before
improving original-arc fidelity. Reuse the current guide and placement machinery;
carry the last circle across guide boundaries, reject forced over-cap moves,
and join every emitted motion. A connector-only plot does not satisfy this
milestone. Keep the input approximation explicit; do not promote this first
working path into full analytic segment/arc acceptance.

**Implementation rule:** keep the current baseline generator intact. Add the paper-derived
path beside it, prove the new path visually and at its geometry boundary, then
ask before removing or redirecting established callers.

**Durable evidence:** work-in-progress publisher sources and rendered comparisons
live under `.superpowers/sdd/2026-09-05-held-paper-figure-reproduction/`. Never
use transient storage or the user-owned `tmp/` tree.

## Task 1: Put publisher panel (a) and repository output in one frame

- [x] Render the official publisher Figure 5 and current 2,289-operation path
  directly to PNG with Matplotlib.
- [x] Extract panel (a)'s purple vector path from the official PDF-page SVG.
- [x] Apply the recorded depicted-tool-radius normalization and source-to-world
  transform; retain the paper's shape-only evidence boundary.
- [x] Render direct SVG and PNG overlays in reconstructed pocket coordinates,
  including the published and repository start positions.

**Gate:** one inspectable same-frame image exposes geometric disagreement without
using page layout, raster registration, or a quality-survey verdict as a proxy.

## Task 2: Add the paper middle-curve construction beside the current baseline

- [x] Write focused RED contracts around the existing `_mathsm_geometry` lane:
  `c != m`, `2 rho = clearance - tool_radius`, boundary-side phase equals `q`,
  and proposals preserve guide progress.
- [x] Expose one additive generator seam that emits maximum-radius, zero-phase
  `MathsmCircleProposal` geometry for a supported MAT edge.
- [x] Prove on a concave-L integration that emitted centres match proposals and
  differ from current baseline skeleton-station centres.
- [x] Run focused `pytest -n auto --testmon`, strict typing, and Ruff gates.

**Gate:** an end-to-end path through the paper-derived middle-curve seam exists
for the already-supported MAT scope. This establishes the stated construction
equations only, not Held implementation or coordinate parity.

## Task 3: Reproduce Figure 5

- [x] Attempt Figure 5 through the supported exact segment-site MAT seam. It
  stalls before producing a usable full-pocket station sequence; retain exact
  MAT traversal as a stated residual and do not expand MAT topology in this
  reproduction task.
- [x] Build and render a diagnostic circle-locus adapter over the current
  baseline station sequence. It exposes the paper `q`/`c`/`rho` construction,
  but it does not reproduce spacing, transitions, start, or engagement control.
- [x] Correct panel registration against the publisher boundary. The symmetric
  boundary consistency check must remain within the recorded reconstruction and
  projection budget; the old crop-corner transform is rejected.
- [x] Preserve every pre-engagement-thinning emitted guide-run station and
  report every projected-boundary site inside the recorded reconstruction and
  projection distance budget as a `ProjectionAdmissibleBoundaryHypothesis`.
  Carry its exact side parameter and canonical closed-ring vertex identity; do
  not present distance admissibility or emitted run IDs as MAT/topology
  ownership. Retain every distinct hypothesis as a boundary-progress candidate
  rather than choosing one by phase, distance, or ordinal; merge only identical
  canonical positions with equivalent geometry, otherwise fail ambiguity. Exact
  generator ownership remains part of the exact-MAT residual. Preserve this
  complete set as input evidence; placement may select the unique nearest
  contact within one cyclically connected boundary neighborhood while retaining
  disconnected alternatives and exact ties.
- [x] Implement Figure 5(a) placement in paper order: construct candidate
  `q`/`c`/`rho` first, then select the next candidate against the requested
  80-degree engagement on that machining circle using the standard predecessor
  model, not the contour-aware global-depletion model. Consecutive centres must
  make strictly positive boundary progress with gradually varying spacing; do
  not post-process the result into cosmetic uniformity. The approximate guide
  route remains negative evidence: 1,532 circles versus 265 publisher turns,
  120 forced over-cap successors, and no publisher-free ownership graph.
- [x] Decode the publisher's single ordered Figure 5(a) vector stream as typed
  shape-only evidence: 265 CCW circular turns, connector primitives, and marker
  proximity to both stream ends. Keep this decoder outside the generator.
- [x] Associate each ordered publisher turn with an unchanged repository
  candidate and retain its run/station/boundary provenance. The publisher may
  choose identity and order only; it must never create, snap, or interpolate
  repository `q`/`c`/`rho` geometry. Report correspondence and engagement
  residuals and label the result reference-guided.
- [x] Preserve the publisher's forward CCW turn order and report its marker as
  terminal evidence derived from the publisher stream endpoints. Render the
  associated repository circles without fabricated transitions. Keep continuous
  radius-one transition recovery open: it requires source-edge-to-offset-edge
  lineage absent from the approximate guide and must not be replaced by snapping.
- [x] Render direct Matplotlib SVG and PNG circle-correspondence comparisons. Separate circle and
  transition residuals, record missing/extra structure, and keep any sampled
  visual metric explicitly graphical rather than a numeric-parity claim.
- [x] Obtain focused adversarial review from Held, Buchli, Shewchuk, Fogel, and
  one disconfirming/null reviewer; fix only reproduction blockers.

**Gate:** a durable Figure 5 comparison contains an ordered reference-guided
circle correspondence and distinguishes publisher shape-only evidence,
paper-derived repository geometry, approximate guide provenance, and measured
residuals. The accepted comparison reports the publisher marker as terminal evidence,
renders the publisher's CCW order with unchanged repository circles, and states
that both 80-degree engagement and continuous transitions remain diagnostic. Its distribution has no per-run
restart bunches, missing observed families, or zero-progress circle runs. Do not
claim a continuous path, independent traversal recovery, exact-MAT, continuous-cap, or unpublished
coordinate parity. Residual: replace the reference-guided order when a full
Figure 5 exact segment-site MAT traversal exists. Focused geometry/consumer
tests pass.

## Task 4: Reproduce Figure 6 and close

- [x] Recover the official Figure 6 plot as durable publisher evidence; treat
  read-off points as digitized graphical observations, not source numbers.
- [x] Render repository results on the same labeled axes and overlay or
  juxtapose them without claiming numeric parity. The bounded reproduction uses
  controlled path lengths at requested caps 80, 120, and 160 degrees; it labels
  cap compliance unaudited and the constant-spacing repository curve unavailable
  because the exact audit/replay route is outside this figure-only task. The
  measured repository lengths are respectively 29,146.462, 22,333.767, and
  20,142.432 mm; the tracked PNG and SVG retain the publisher's 0--200 degree,
  20-degree-tick, log-10 path-length axes and distinguish the publisher's
  graphical unit from repository millimetres. A second tracked PNG places the
  repository points directly over the unchanged publisher pixels for literal
  registration inspection while preserving the same claim boundary.
- [x] Run affected tests, strict typing, Ruff, docs/plan gates, and diff hygiene.
- [x] Update this ledger with measured results and commit the focused work as
  Jelle Feringa without staging `tmp/` (`bcfc1a1`).

**Gate:** Figures 5 and 6 have truthful publisher/repository visual comparisons,
the additive implementation is verified, and every checkbox above is closed.

## Completion audit: why the toolpath goal remains open

The committed Figure 5 result is a publisher-ordered association of 265
repository circles, not a generated toolpath: every one of its 264 transitions
is unresolved, 137 successor relationships exceed the requested 80-degree
diagnostic limit, and publisher observations select and order the circles from
30,684 candidates. The publisher-free approximate route emits 1,532 circles and
120 forced successors. Figure 6 measures the separate existing controlled
generator at three requested caps, includes non-planar/restart operations in its
length, does not establish cap compliance, and lacks the contour-aware and
constant-spacing peers shown by the paper. Its publisher ordinate and repository
millimetres also have no established common scale. Those truthful limitations
contradict completion of the actual toolpath-reproduction goal even though the
earlier bounded comparison tasks are complete.

## Task 5: Generate the continuous Figure 5 standard path

**Engagement checkpoint (September 6):** the additive polygon refinement
resolves the tagged baseline's 108 over-cap successors: 2,013 circles, maximum
79.992 degrees, 2,012 connectors, and all source family/contact mappings
retained. It repairs convex-corner radii, adds concave contact fans, and carries
one predecessor through final order. Reproduce with
`pixi run held-figure5-toolpath-progress --refine`. The focused suite passes
21 tests and focused review found no new defect. This closes the measured
engagement violation set for the approximate circle model; it does not close
the full criteria below. Exact containment still rejects 1,013 near-tangent
float representations (largest reported excess about 6.7e-15 mm), and the
publisher-derived start, analytic ownership, entry/connector engagement, and
depleted-stock qualification remain open.

**Connector checkpoint (September 6):** native-owned boundary contacts and
`ccw_transition` expose exact line/arc connectors with preserved source lineage.
The consumer handles irrational coordinates without a reporting-double
round-trip. `pixi run held-figure5-exact-connector` plots this bounded milestone.
This interface does not yet integrate machining circles or close a checkbox.

**Deferred input-fidelity issue:** all 26 reconstructed arcs have unequal exact
endpoint radii around their stored centres (largest radial discrepancy about
`8.54e-16 mm`); all primitive joins share endpoints. An exact-arc construction
policy and renewed reconstruction bounds will be needed for literal segment/arc
MAT input. The user explicitly prioritized the working toolpath over this
refinement. Do not make exactification a prerequisite for the first continuous
path on existing input. Polygon or sampled-point Voronoi geometry still does
not close the full analytic ownership criterion below.

- [ ] Complete publisher-free traversal and boundary-site ownership over the
  reconstructed Figure 5 segment/arc medial axis. Publisher evidence may be
  loaded only after generation for comparison.
- [ ] Place machining circles in traversal order with the paper's standard
  predecessor engagement construction at 80 degrees; do not retain forced
  over-cap successors or per-guide-run restarts.
- [ ] Emit one continuous planar path with an explicit start, every full CCW
  machining circle, every source-lineage CCW offset transition, and a terminal
  point. Preserve the current generator beside it until this path is verified.
- [ ] Render every repository motion over the publisher Figure 5 path and report
  circle-distribution, transition, continuity, start/end, and engagement
  residuals without publisher-assisted identity or ordering.
- [ ] Run focused tests, strict typing, Ruff, docs/plan gates, and diff hygiene;
  obtain the approved Held, Buchli, Shewchuk, Fogel, and null-panel review.

**Gate:** the repository independently emits one continuous Figure 5(a)
standard toolpath whose visible circle families and connectors agree with the
publisher evidence, whose consecutive motions share endpoints, and whose
non-entry machining moves respect the requested 80-degree predecessor limit.

## Task 6: Reproduce the Figure 6 protocol from the Figure 5 generator

- [ ] Sweep the Figure 5 standard generator over the publisher engagement axis
  and record achieved maximum engagement rather than requested-cap labels.
- [ ] Add the paper's contour-aware spacing peer and constant-spacing MATHSM
  sweep on the same pocket, tool, start, traversal, and planar motion grammar.
- [ ] Measure only the paper's planar straight and circular path elements;
  exclude clearance, retract, restart, and other machine-routing motions.
- [ ] Establish a common dimensionless length normalization from publisher
  Figure 5 geometry and the depicted tool radius; do not compare graphical
  publisher ordinates directly with repository millimetres.
- [ ] Render the three repository curves with the publisher curves on the
  identical axes, report curve residuals and missing samples, and visually
  inspect PNG output.
- [ ] Run focused tests, strict typing, Ruff, docs/plan gates, diff hygiene, and
  commit the completed toolpath reproduction as Jelle Feringa without staging
  `tmp/`.

**Gate:** Figure 6 is generated from the completed Figure 5 path family and
contains standard, contour-aware, and MATHSM curves compared under common
engagement, length, scale, pocket, tool, start, and traversal semantics.
