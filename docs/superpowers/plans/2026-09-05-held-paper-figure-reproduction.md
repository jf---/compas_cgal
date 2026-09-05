# Held Paper Figure Reproduction Plan

> **status: in progress** - Figure 5 panel (a) is extracted and rendered in the
> reconstructed world frame; the additive middle-curve generator seam is active.

**Goal:** Reproduce the visible Figure 5 toolpath structure and the Figure 6
comparison protocol and axes, with direct comparisons to publisher evidence.

**Active criterion:** the repository generator must emit a Figure 5 toolpath
whose paper-derived middle-curve construction satisfies the stated equations
and whose rendering exposes structural agreement and residual disagreement
against the publisher panel.
Characterization, postprocessing, machine validation, and unrelated benchmark
work do not advance this criterion.

**Implementation rule:** keep the legacy generator intact. Add the paper-derived
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

## Task 2: Add the paper middle-curve construction beside the legacy path

- [x] Write focused RED contracts around the existing `_mathsm_geometry` lane:
  `c != m`, `2 rho = clearance - tool_radius`, boundary-side phase equals `q`,
  and proposals preserve guide progress.
- [x] Expose one additive generator seam that emits maximum-radius, zero-phase
  `MathsmCircleProposal` geometry for a supported MAT edge.
- [x] Prove on a concave-L integration that emitted centres match proposals and
  differ from legacy skeleton-station centres.
- [x] Run focused `pytest -n auto --testmon`, strict typing, and Ruff gates.

**Gate:** an end-to-end path through the paper-derived middle-curve seam exists
for the already-supported MAT scope. This establishes the stated construction
equations only, not Held implementation or coordinate parity.

## Task 3: Reproduce Figure 5

- [ ] Attempt Figure 5 through the supported MAT scope first. If it fails, add
  only the first encountered unsupported site-topology case needed to continue
  this reproduction; do not expand certificate, replay, coverage, or general
  MAT scope.
- [ ] Generate the Figure 5 path through the additive middle-curve seam using
  the published 1 mm tool radius, 80-degree cap, and start observation.
- [ ] Render Matplotlib SVG and PNG comparisons against the extracted publisher
  panels and record visible residual differences without claiming unpublished
  coordinate parity.
- [ ] Obtain focused adversarial review from Held, Buchli, Shewchuk, Fogel, and
  one disconfirming/null reviewer; fix only reproduction blockers.

**Gate:** the comparison shows the same declared path-family structures and
explicitly displays remaining differences; no unpublished coordinates or
numeric parity are claimed. Focused geometry/consumer tests pass.

## Task 4: Reproduce Figure 6 and close

- [ ] Recover the official Figure 6 plot as durable publisher evidence; treat
  read-off points as digitized graphical observations, not source numbers.
- [ ] Render repository results on the same labeled axes and overlay or
  juxtapose them without claiming numeric parity.
- [ ] Run affected tests, strict typing, Ruff, docs/plan gates, and diff hygiene.
- [ ] Update this ledger with measured results and commit the focused work as
  Jelle Feringa without staging `tmp/`.

**Gate:** Figures 5 and 6 have truthful publisher/repository visual comparisons,
the additive implementation is verified, and every checkbox above is closed.
