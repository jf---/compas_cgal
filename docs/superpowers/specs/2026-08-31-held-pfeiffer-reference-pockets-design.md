# Held-Pfeiffer Reference Pockets Design

> **status: approved** - approved 2026-08-31.

## Decision

Reconstruct the four distinct pocket outlines published in Held and Pfeiffer's
2025 paper as normalized, provenance-carrying reference problems for the
toolpath generators:

1. the Figure 5 pocket, with Figure 7 as an independent visual observation of
   the same geometry;
2. the upper Figure 8 real-world pocket;
3. the Figure 8 crossed-skis pocket; and
4. the Figure 8 *Monstera deliciosa* pocket.

The reconstruction starts from the publisher PDF's vector drawing operations,
not screenshots. It recovers intended straight segments and circular arcs from
the emitted line and cubic-Bezier paths. A genuinely non-circular cubic is
approximated by tangent-continuous biarcs only after a single-circle fit has
failed its measured source-fidelity gate.

The paper does not publish physical dimensions. Every case therefore uses a
declared normalized instantiation in which the depicted tool radius equals one
millimetre. These files are normalized reconstructions from published vector
figures, not the authors' original CAD inputs.

## Objective

Replace the current rectangular stand-in with reference pockets that exercise
the geometry actually shown in the Held-Pfeiffer results. A completed case must
be loadable by the repository's benchmark boundary, must retain an analytic
line/arc representation, and must produce a measured polygonal projection for
the current polygon-only generator without hiding the projection error.

This stage does not claim that the generator matches or surpasses the paper.
It supplies the inputs on which that claim can later be tested.

## Source authority

The sole publication source is:

> Martin Held and Josef Pfeiffer, "Trochoidal Tool Paths for Pocket Machining
> with Full Control of the Tool Engagement Angle," *Computer-Aided Design and
> Applications* 22(5), 2025, pp. 731-747,
> <https://doi.org/10.14733/cadaps.2025.731-747>.

The publisher PDF is read locally and is not committed. The reconstruction
records the publication page, figure, subfigure, source crop bounds, PDF page
frame, boundary stroke width, and depicted tool circle used for normalization.

| Case | Publication page | PDF page | Role |
| --- | ---: | ---: | --- |
| Figure 5 | 742 | 12 | Primary benchmark pocket and start/tool observation |
| Figure 7 | 744 | 14 | Cross-check of the Figure 5 silhouette only |
| Figure 8 upper | 746 | 16 | Independent real-world pocket |
| Figure 8 skis | 746 | 16 | Independent narrow-neck pocket |
| Figure 8 Monstera | 746 | 16 | Independent high-complexity pocket |

Figure 7 never creates a second pocket file. Its registered silhouette checks
the Figure 5 reconstruction and is retained as supporting evidence.

## Geometry model

### Frames and units

PDF extraction uses `PdfPoint2`, whose coordinates carry a dedicated
`PdfPointUnit` `NewType`. Reconstructed geometry uses the existing
`Point2[WorldXY]`, `Millimetre`, `ToolRadius`, and `Radian` types. A
`MillimetresPerPdfPoint` `NewType` carries the source-to-world scale after the
normalizing tool circle has been measured.

The source-to-world transformation performs, in order:

1. the PDF path's affine transformation;
2. the PDF Y-axis inversion;
3. translation of the pocket's lower-left bounds to the world origin; and
4. uniform scaling such that the depicted red tool circle has
   `ToolRadius.build(1.0)`.

The resulting outer boundary is counter-clockwise. No case contains an island:
the crossed skis and Monstera slots are concavities in a single Jordan boundary,
consistent with the paper's simply-connected input scope.

`SourceToWorld` carries the PDF Y-reflection choice explicitly; source scale
remains a positive `MillimetresPerPdfPoint`. This keeps handedness conversion
separate from physical scale and prevents normalized millimetres from being
mislabelled as PDF points.

### Source primitives

The extractor retains only boundary paths within the approved figure crops.
Colour is a crop-local selection aid, never sufficient authority by itself.
Selected paths must also have the boundary stroke width, join style, expected
transform, and endpoints incident to the visible black boundary markers.

The source grammar contains two closed variants:

- `SourceLine`, carrying two `PdfPoint2` endpoints; and
- `SourceCubic`, carrying two endpoints and two control points in the same
  frame.

Factories validate finite coordinates, non-zero extent, supported affine
transforms, and endpoint continuity. Unsupported drawing operators raise a
named `UnsupportedPdfBoundaryOperatorError`; disconnected selected paths raise
`DisconnectedPublishedBoundaryError`.

Poppler exposes three intended junctions with a one-quantum seam: two in the
crossed skis and one in Monstera. At ingestion only, unique degree-one endpoint
pairs no farther apart than the PDF's `1/256 pt` coordinate quantum are replaced
by one canonical endpoint. The next-nearest candidates are recorded and must
remain outside that bound. All downstream continuity is exact equality.

### Analytic reconstruction

Each source line becomes one `ReferenceLine`. Each source cubic follows this
decision sequence:

1. Construct the unique equal-radius circle through both endpoints tangent to
   the cubic at its start, then require the end tangent to agree in direction.
2. If that finite supporting circle exists, construct the candidate circular
   arc with the authored sweep direction.
3. Measure the maximum centreline separation between the cubic and candidate
   arc using adaptive subdivision until the measurement bound closes.
4. Accept the single arc only when that bound is no greater than one quarter of
   the normalized boundary stroke width.
5. Otherwise construct an equal-distance G1 biarc from the two endpoint
   positions and tangents in a chord-normalized local frame. Certify the
   represented construction with operation-derived backward and forward error
   bounds, measure it by a closed continuous-correspondence bound, and accept
   it only when both certificates close.
6. If no circle or biarc certificate closes, accept the endpoint chord only
   when the stored cubic control polygon is exactly collinear and its authored
   endpoint directions advance along that chord. Subdivide every non-collinear
   source cubic recursively.

The quarter-stroke rule is a publication-resolution limit: the reconstruction's
centreline must remain well inside the printed boundary stroke. It scales with
the figure and tool rather than embedding an unexplained numerical tolerance.

Reconstruction returns the analytic primitives together with the certified
continuous-deviation upper bound proved by the same traversal. It does not
discard that evidence or relabel it as a sampled maximum.

The biarc proof closes stored-endpoint continuity with an auxiliary linear
endpoint correction and adds that correction's complete distance to the ideal
circular locus to the bound. The emitted primitives remain circular arcs; the
auxiliary path is proof machinery, not replacement geometry.

Adjacent recovered arcs are merged only while their source spans remain
available and those combined spans re-certify against one candidate supporting
circle under the original limit. Their common endpoint and sweep direction
must also agree. An arbitrary cubic is never labelled an exact circular arc,
and a cubic becomes a line only through the exact-collinearity certificate
above.
Failure to close the measurement bound raises
`UnresolvedPublishedCurveError`; there is no permissive fallback.

The reconstructed boundary is the canonical reference. A polygonal generator
projection is derived from it using a chord-deviation bound equal to one half
of the reconstruction limit. The projection records its observed maximum
deviation and vertex count. It never replaces the analytic source in the case
file.

## Case format

Each committed JSON document contains:

- schema name and version;
- stable case name;
- paper citation, publication page, figure, and subfigure;
- normalization statement and `tool_radius_mm: 1.0`;
- `tea_cap_deg: 80.0`, matching Figures 5 and 8;
- the normalized start position when the red start marker is present;
- ordered analytic line/arc primitives;
- ordered polygon-projection vertices;
- normalized source stroke width;
- reconstruction and projection deviation limits;
- certified reconstruction deviation upper bound and measured projection
  deviation; and
- Figure 7 shape-only observation metadata for Figure 5 only.

The JSON must round-trip through `json.loads`, validate against the repository's
schema, and build both the typed analytic boundary and a `PocketSpec`. Unknown
fields, operators, units, frames, or schema versions fail with named errors.

No paper path-length or timing values enter these case files. Those are later
measurements over the reconstructed inputs, not properties of the geometry.

## Components and responsibilities

### `tools/held_reference_extractor.py`

Owns publisher-PDF page conversion, crop-local vector selection, affine
transformation, marker association, and emission of source primitives. It is
invoked only through a Pixi task accepting the local publisher PDF path.

### `benchmarks/held_reference_geometry.py`

Owns typed line, cubic, arc, and biarc reconstruction; source-fidelity
measurement; normalization; primitive merging; and polygon projection. It has
no PDF I/O and no generator dependency.

### `benchmarks/held_reference_cases.py`

Owns case-file parsing, schema validation, named case lookup, construction of
the analytic boundary, and projection into `PocketSpec`. It never fits curves
or generates toolpaths.

### `benchmarks/data/held_pfeiffer_2025/`

Contains exactly four case JSON documents and a short source note. Figure 7
evidence is nested in the Figure 5 case rather than duplicated.

### `docs/held_pfeiffer_reference_pockets.md`

Documents the reconstruction method, normalized-unit boundary, limitations,
and one source-versus-reconstruction overlay per case. Overlay images are PNG.

## Consumer boundary

The first acceptance boundary is the existing benchmark model:

```text
publisher vector figure
    -> typed source primitives
    -> normalized analytic boundary
    -> measured polygon projection
    -> PocketSpec.build(...)
    -> existing generator
```

All four cases must pass `PocketSpec.build(...)` as simple, non-degenerate,
counter-clockwise pockets with a two-millimetre tool diameter and an 80-degree
engagement cap. A focused qualification command then calls the existing
generator for each case and reports whether it emitted a non-empty path or
failed through an existing named product error. The corpus implementation does
not alter generator behavior to make a case pass.

The qualification result is product evidence, not a corpus-validity gate. A
generator failure leaves the geometrically valid reference case committed and
identifies the next product blocker.

## Validation

### Structural contracts

- exactly four unique cases load;
- Figure 5 owns exactly one Figure 7 observation;
- every case has one closed, counter-clockwise, non-self-intersecting outer
  boundary and no holes;
- all adjacent analytic primitives meet at identical endpoints;
- line and arc factories reject degenerate input;
- the depicted tool radius maps to exactly one declared millimetre;
- schema validation rejects malformed, extra, or ambiguous fields; and
- each polygon projection round-trips through `PocketSpec.build(...)`.

### Geometric contracts

- every accepted circle or biarc closes its adaptive source-fidelity bound;
- observed reconstruction deviation is within one quarter of the normalized
  boundary stroke width;
- polygon projection deviation is within half of that reconstruction limit;
- each overlay renders source centreline, analytic reconstruction, primitive
  junctions, projected vertices, start marker, and normalized tool circle; and
- the Figure 7 overlay registers its coloured tool-centre samples against the
  one-tool-radius inward offset of the Figure 5 analytic reconstruction.

The publisher plots use independent horizontal and vertical display scales:
at 600 DPI their sample-envelope aspect ratio is `1.43624`, while the analytic
one-radius inward support is `1.47931`. Registration therefore uses a typed
axis-aligned display affine with independently derived X/Y scales and PDF-to-
world Y reflection. A uniform similarity is rejected because it leaves
tens-of-pixels residual consistently in all three panels.

Figure 7 contains no independently drawn boundary stroke: its thousands of
coloured paths are tool-centre samples. It is therefore retained as a
shape-only falsification overlay, not a numeric boundary-fidelity gate and
never authority for geometry.

### Repository gates

- focused tests run with `-n auto`, followed by affected tests with
  `--testmon -n auto`;
- Python files pass Ruff formatting and lint;
- the new typed modules pass strict mypy;
- MkDocs builds strictly;
- regenerated JSON and PNG outputs match the committed semantic content and
  visual overlays; and
- new work is scanned for prohibited lock-down terminology and mechanisms.

## Explicit limitations

- Absolute millimetre scale is unavailable; the one-millimetre radius is a
  declared normalized instantiation.
- PDF paths are publication geometry, not the authors' original CAD files.
- A cubic accepted as a circular arc is a measured reconstruction of likely
  intent, not proof of its pre-publication representation.
- Figure 7 supplies only a silhouette cross-check.
- Toolpath quality, cap compliance, coverage, cycle time, and comparison with
  Held-Pfeiffer paths remain subsequent measured stages.

## Completion criterion

This stage is complete when the four case files, typed loader, reconstruction
tool, contract tests, documentation, and visually inspected overlays are
committed in the isolated worktree; all corpus gates are green; and a focused
generator qualification report states the current outcome for each case
without changing generator behavior.
