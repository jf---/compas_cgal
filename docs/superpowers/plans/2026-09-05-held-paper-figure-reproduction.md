# Held Paper Figure Reproduction Plan

> **status: in progress** - Task 8 whole-paper reproduction drafts is active.
> Reuse the completed reference corpus and current Figure 5 draft. Task 7
> performance tracking follows; full Task 5/6 acceptance remains open.

**Goal:** Produce repository-generated reproduction drafts for every figure and
panel in Held–Pfeiffer 2025, then measure and improve their generating workloads.
Retain the stricter Figure 5/6 reproduction criteria below.

**Paper:** [committed publisher PDF](../../assets/papers/held-pfeiffer-2025.pdf).
Use this repository copy for source inspection and figure reproduction:
[Figure 5](../../assets/papers/held-pfeiffer-2025.pdf#page=12),
[Figure 8](../../assets/papers/held-pfeiffer-2025.pdf#page=16),
[Section 3.1](../../assets/papers/held-pfeiffer-2025.pdf#page=10).

**Active criterion:** Task 8 produces an inspectable whole-paper draft gallery
from the existing prepared inputs and generator seams. Every panel has a
reproducible generated result or an explicit, evidenced blocker. Reference
images and blocker entries never count as generated reproductions.

**Approved execution order (September 6):** Task 8 draft gallery → Task 7
performance baseline and one measured optimization → outstanding Task 5/6
acceptance. Appended task numbers preserve historical references. This replaces
the earlier instruction to run Task 7 immediately after the Figure 5 checkpoint.
Postprocessing and machine validation remain outside this reproduction plan.

**Preserved toolpath-first direction:** reuse the current guide and placement
machinery, carry the predecessor across guide boundaries, reject forced
successors, and render every emitted motion. Do not make original-arc fidelity
or full exact-MAT ownership prerequisites for clearly labeled drafts. Those
limitations remain failures against the stricter reproduction gates. Exact
arithmetic belongs in CGAL; Python orchestrates and reports.

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

**Native arithmetic checkpoint (September 6):** CGAL now owns disk containment,
swept-disk intersection construction, and corner orientation. Approximate
Python interpolation no longer performs rational arithmetic. The 22 affected
tests pass, including full Figure 5 motion; the regenerated plot retains
2,013 circles and zero engagement violations. Full Task 5 remains open.

**Engagement checkpoint (September 6):** the additive polygon refinement
resolves the tagged baseline's 108 over-cap successors: 2,013 circles, maximum
79.992 degrees, 2,012 connectors, and all source family/contact mappings
retained. It repairs convex-corner radii, adds concave contact fans, and carries
one predecessor through final order. Reproduce with
`pixi run held-figure5-toolpath-progress --refine`. The focused suite passes
21 tests and focused review found no new defect. This closes the measured
engagement violation set for the approximate circle model; it does not close
the full criteria below. Exact containment still rejects 1,010 near-tangent
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

## Task 7: Track Figure 5 performance with heatmaps

**Dependency:** Task 8 draft gallery and its explicit workload inventory.

**Scope:** use the repaired Figure 5 path at 80 degrees as the first baseline,
then apply the same accounting to every runnable toolpath workload exposed by
Task 8: the three algorithm variants and prepared Figure 8 cases. Preserve
missing variants as blocked workloads, not omitted or passing samples. Keep
construction-diagram/rendering cost separate from planner cost. Extend the
existing benchmark and plotting flow. First establish where time and work accumulate;
then optimize one demonstrated cause and compare before/after. Exact geometric
arithmetic belongs in CGAL; Python orchestrates measurement and visualization.

**Starting evidence:** one in-process Apple M1 Max run at `1a4db8e`, with the
guide loaded before timing, measured 19.067 s: 19.010 s in the initial path
builder, 0.048 s in refinement, and 0.009 s in connectors. The initial builder
includes hypothesis generation, ownership, ordering, and circle selection;
this measurement does not isolate the bottleneck. Imports, I/O, guide
construction, plotting, and containment auditing were excluded. The local
measurement and plot live in `build/held-performance/`; these are diagnostic
artifacts, not yet the reproducible benchmark deliverable below.

- [ ] Record timings for boundary projections/hypothesis generation, ownership,
  ordering, lane assembly, candidate selection, refinement, and connectors.
  Use disjoint stage timings and an explicit unattributed remainder so their
  accounting reconciles with total planner time without double-counting.
  Report guide construction and correctness validation on separate clocks.
- [ ] Record structural counts alongside time: stations, boundary projections,
  candidates examined, engagement evaluations, emitted circles, and planar
  path length. Inspect the observed all-boundary station projections, lane
  scans, and candidate-suffix validation as suspects; do not label them
  bottlenecks before measurement. Use profiling only to test a remaining
  explicit performance hypothesis after structural counts are understood.
- [ ] Attribute local work to existing guide runs/stations. Render a spatial
  heatmap over the toolpath, retaining global work in a separate bucket.
  Show total cost and cost per station/candidate separately; do not fabricate
  spatial attribution or infer redundant cutting without stock evidence.
- [ ] Render a stages-by-successive-runs heatmap with absolute timings and
  changes from the baseline. Keep color scales and units fixed across compared
  runs; label any logarithmic scale. Preserve run order and source attribution.
- [ ] Separate diagnostic instrumentation from uninstrumented timing runs and
  measure instrumentation overhead. Define the warm-up and repeat protocol
  before collecting comparisons; report samples, median, spread, hardware,
  build, input, cap, and timing exclusions. Save ordinary run records and
  machine-readable measurements through a reproducible Pixi task in durable
  project paths. Held's published timing remains an external reference, not
  a same-machine speedup claim.
- [ ] Keep engagement violations, connector continuity, containment failures,
  source coverage, and path length beside performance results. Do not obtain
  speed gains by dropping required motion, weakening predicates, or converting
  existing containment failures into accepted results.
- [ ] Verify timing/count aggregation and source attribution with focused
  consumer tests; run strict typing, Ruff, docs/plan gates, and visually inspect
  the PNG heatmaps. Document measurement scope and limitations in
  `docs/segment_site_mat.md`.

**First gate:** reproduce the Figure 5/80-degree baseline, account for planner
time without double-counting, identify the dominant stage and its work count,
and produce both heatmaps with correctness results alongside them. Satisfy
this gate before choosing an optimization.

- [ ] From that evidence, select one demonstrated cause, state its mechanism
  and smallest correction, implement and validate it, and show before/after
  timings, counts, heatmaps, and correctness results under the same protocol.

**Completion gate:** the reproducible tracking command covers the runnable
Task 8 planner workloads, reports blocked workloads explicitly, and exposes the
measured cause and effect of one targeted optimization. Compare common cases
under identical protocols; do not aggregate unlike pockets, algorithms, or caps
into a speedup. Retain unchanged correctness
requirements and explicit remaining failures. This does not close Task 5/6 or
establish Held-level performance by itself.


## Task 8: Whole-paper reproduction drafts from prepared inputs

**Active step.** Complete breadth of draft coverage before Task 7 optimization.

**Contour-aware contact checkpoint (September 6):** native
`HeldDiskContour2` maintains the paper's filled outer-disk union and corrects
the standard critical point clockwise to the exposed predecessor arc.
Twelve native contour cases plus seven stock/replay contracts pass. Inspected
256-circle prefixes show 131 corrected contacts for Figure 5 and 173 for
Figure 8 upper, with JSON coordinate reports. The existing generation CLIs
accept `--contour-prefix 256`; the checkpoint renderer reused retained drafts.
This is a contact-query consumer, not a new placement result. Maximum
engagement with the corrected contact (including Figure 4(d)), evolving-contour
candidate selection, and full contour-aware Figure 5/8 drafts remain open.
The polygon-set union/scan does not yet implement Held's linear-time arc update.

**Figure 5(c) stock consumer checkpoint (September 6):** native circle-stock
replay retains earlier cuts and uncut centre islands, with radii formed in
CGAL. Seven native/replay contracts pass, including a pointwise engagement
decision that distinguishes full history from predecessor-only stock. The
256-circle prefix PNG is generated and inspected; reproduce with
`held-figure5-toolpath-progress --refine --stock-prefix 256`. This is an
intermediate circle-only stock view. Entry/connector clearing and the
contour-aware maximum-over-candidate-circle placement query remain open;
Figure 5(c) is not marked reproduced.

**Monstera motion checkpoint (September 6):** qualification reaches 6,747
initial circles; boundary-order selection retains 2,696 sources. A demonstrated
concave-corner contact jump blocks refinement at 111.590239 degrees despite
16 subdivisions. The additive `--corner-approaches` path separates approach,
shared-contact rotation and departure. Monstera then emits 4,367 circles and
4,366 connected transitions, maximum predecessor engagement 79.998367 degrees,
total path length 46,719.894 mm, and all 739 source runs retained. The full
motion PNG and rejected-contact JSON are recorded in the stage documentation.
The final affected suite passes nine tests, with strict typing and Ruff green.
All three Figure 8 pockets now have standard-model drafts. Their contour-aware
reproductions, containment and continuous coverage remain open.

**Monstera contact qualification checkpoint (September 6):** the opt-in
`--qualify-contacts` route retains 120,257 bounded contact hypotheses and records
15 rejected alternatives. Every one of the 739 source runs retains at least
one contact. The evidence bound is unchanged; ambiguity and missing-run
coverage remain errors. The focused qualification/placement/spacing suite
passes 21 tests, with Ruff and strict typing green. This resolves candidate
eligibility; the subsequent motion checkpoint above resolves its refinement
failure.

**Boundary-order spacing checkpoint (September 6):** additive experimental
`held-reference-toolpaths --case figure8_upper --spacing boundary` carries one
predecessor through emitted order and retains all 275 source runs. The upper
draft falls from 4,183 to 1,684 circles and from 44,186.880 to 18,236.833 mm
including connectors, with zero predecessor-cap violations. A 0.25 mm grid
finds no sampled swept-coverage change; continuous coverage remains unproved.
Four focused contracts, Ruff and strict typing pass. The original lane route
remains default. This improves draft spacing, not Held parity or contour-aware
semantics. Monstera was still missing at this checkpoint; its motion is now
recorded above.

**Figure 8 distribution evidence (September 6):** the upper draft has 3,114
initial circles across 480 placement lanes and 275 source runs. Refinement adds
1,069 circles. In final order, 1,795 of 4,182 predecessor pairs use less than
half the 80-degree cap; the median is 56.918 degrees. The generated distribution
map separates repaired source circles from insertions and exposes repeated
radial groups around the rounded boundary. This diagnoses a spacing discrepancy,
not permission to discard source families. The next spacing experiment must
carry the predecessor across family boundaries while preserving coverage and
checking every new adjacency. The subsequent checkpoint above records the
first measured improvement against this baseline.
The corpus command now emits the distribution PNG and reporting JSON alongside
each lane-refined path. Monstera's later qualification and motion checkpoints
are recorded above.

**Figure 8 draft checkpoint (September 6):** the additive corpus command
`pixi run held-reference-toolpaths --case CASE` now reuses prepared inputs.
Upper and crossed-skis standard-model drafts emit respectively 4,183/2,958
circles, with every connector rendered and zero predecessor-cap violations.
Native contact projection fixes source/offset side-index mismatch; retained
run anchors prevent spacing selection from losing an entire source run.
CGAL corrected-chord arithmetic fixes the upper-pocket near-tangent witness.
Figure 5/upper integration, placement, and native engagement checks pass.
Monstera initially stopped at a contact beyond its evidence bound; the stage
documentation records the witness and the subsequent repair above.
Contour-aware Figure 5/8 semantics and
the broader gallery remain open; these standard drafts close no full-panel
reproduction gate.

The prepared corpus is already delivered by
[the reference-pocket plan](2026-08-31-held-pfeiffer-reference-pockets.md) and
[its documentation](../../held_pfeiffer_reference_pockets.md). Do not repeat
source extraction or geometry reconstruction unless a concrete missing input
or invalid artifact is demonstrated.

### Existing preparation to reuse

- [x] Four committed pocket inputs: `figure5.json`, `figure8_upper.json`,
  `figure8_crossed_skis.json`, and `figure8_monstera.json` under
  `benchmarks/data/held_pfeiffer_2025/`. Each retains publisher primitives,
  normalization, reconstructed boundary, and polygon projection.
- [x] Four reference-overlay PNGs and Figure 7's three shape-evidence panels,
  handled by `benchmarks/held_reference_figures.py` and
  `pixi run held-reference-figures`.
- [x] Figure 6 digitized publisher curves, same-axes/literal comparison plots,
  and `pixi run held-figure6-same-axes`.
- [x] Existing four-case generator qualification entry point,
  `pixi run held-reference-qualify`; its recorded non-empty paths are baseline
  evidence, not accepted Held reproductions. Reuse the entry point and inputs.
- [x] Figure 5(a) ordered publisher-turn evidence and the current repaired
  circle/connector draft, reproducible with
  `pixi run held-figure5-toolpath-progress --refine`.

### Figure inventory and draft deliverables

Inventory checked against the
[committed publisher paper](../../assets/papers/held-pfeiffer-2025.pdf).
The rows specify required output, not claims that a generator already exists.

| Figure | Panels/content required | Starting point |
| --- | --- | --- |
| 1 | a–c: cutting width and engagement for linear/circular motion | Identify reusable construction code; preparation not yet established |
| 2 | a–d: pocket machinability transformation and Voronoi/offset views | Reuse native geometry seams; preparation not yet established |
| 3 | Machining-circle and transition construction | Existing circle/connector geometry |
| 4 | a–d: engagement construction, maximum, and corrected overlap case | Standard predecessor model and native circle geometry |
| 5 | a standard; b contour-aware; c intermediate machined contour/middle curve; d MATHSM | Prepared Figure 5 input and current a draft |
| 6 | Path length versus engagement cap for all three algorithms | Prepared curves, scale calibration, and comparison renderer |
| 7 | a–c: engagement maps for those three algorithms | Prepared reference panels and shared Figure 5 geometry |
| 8 | Upper pocket, crossed skis, Monstera contour-aware paths | Three prepared inputs, overlays, and qualification consumer |

- [ ] Reconcile each panel with its existing source asset, input case,
  generating callable, output, and semantic limitations. Confirm inventory by
  inspecting publisher panels. Mark only demonstrated preparation complete;
  locate missing assets before proposing any new extraction.
- [ ] Produce computed construction diagrams for Figures 1–4 and generated
  drafts for every Figure 5–8 panel. Reuse the prepared inputs and additive
  generator paths. Geometry annotations must derive from the displayed
  construction; do not trace publisher toolpaths to manufacture results.
- [ ] Keep standard, contour-aware, and constant-spacing MATHSM distinct.
  Figures 5(c), 7, and 8 must use the relevant generated motion/stock state;
  do not relabel the standard predecessor model as contour-aware or use its
  pairwise maximum as a local depleted-stock engagement measurement.
- [ ] Reuse shared motion results for Figure 5 paths, Figure 6 lengths, and
  Figure 7 engagement maps under common pocket/tool/start/algorithm/cap
  settings. Record sampled versus certified measurements explicitly. Record
  uncovered cap samples and missing algorithm variants as failures.
- [ ] Publish one browsable MkDocs gallery with reference, generated draft,
  discrepancy, and missing capability for every panel. A blocked panel gets
  the failed command or named missing consumer contract and a concrete repair
  task; a publisher crop or placeholder is not a reproduced panel.
- [ ] Provide reproducible Pixi generation commands and durable PNG outputs.
  Inspect all generated panels; run focused geometry/consumer checks, strict
  typing, Ruff, docs/plan gates, and diff hygiene for the changed scope.
- [ ] Export the gallery's runnable planner cases and blocked cases as the
  Task 7 workload inventory, with ordinary case/algorithm/cap identifiers.
  Update `docs/segment_site_mat.md` with draft coverage and remaining gaps.

**Draft coverage gate:** every inventory row/panel is accounted for by an
inspected generated draft or an evidenced blocker with an explicit repair
item. Report generated/required panel counts and blockers separately. This
permits Task 7 to measure runnable workloads; it does not declare blocked
panels reproduced or close their repair items.

**Reproduction draft completion gate:** every required panel has a computed,
inspectable draft under its stated algorithm semantics; no missing panel may
be counted complete. Approximation and measured discrepancies remain visible.
The exact independent-generation and quantitative gates in Tasks 5/6 remain
open until their own evidence passes.
