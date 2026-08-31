# Held-Pfeiffer Reference Pockets Implementation Plan

> **status: in progress** - approved design; geometry reconstruction is the
> active implementation slice.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` or `superpowers:executing-plans` to
> implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for
> tracking.

**Goal:** Reconstruct the four distinct Figure 5/8 Held-Pfeiffer pockets from
publisher-PDF vectors, retain typed analytic line/arc boundaries, project them
into the current polygon benchmark boundary with measured error, and qualify
the existing generator on all four inputs.

**Architecture:** A PDF-specific extractor emits typed line and cubic source
primitives. A geometry-only module normalizes by the depicted tool radius,
recovers circles before recursively splitting non-circular cubics into G1 arc
pairs, and produces bounded polygon projections. A strict case loader consumes
four committed JSON documents; overlays and qualification reports are derived
consumers of those case documents.

**Tech Stack:** Python 3.12, COMPAS geometry, `jsonschema`, Poppler
`pdftocairo`, Matplotlib/Pillow, Pixi, pytest-xdist, mypy strict, Ruff, MkDocs.

**Spec:**
`docs/superpowers/specs/2026-08-31-held-pfeiffer-reference-pockets-design.md`

## Global Constraints

- Work only in `codex/held-reference-corpus` and never mutate `main`.
- Use the publisher PDF as local input; do not commit the paper.
- Reconstruct exactly four cases: Figure 5, Figure 8 upper, Figure 8 crossed
  skis, and Figure 8 Monstera. Figure 7 is evidence for Figure 5, not a case.
- Normalize every depicted tool radius to `ToolRadius.build(1.0)` and record
  that the physical scale is not published.
- Preserve analytic line/arc geometry as authority; the current generator sees
  only a derived polygon projection with recorded deviation.
- Attempt one circular arc before recursively splitting a cubic into G1 arc
  pairs. Never silently accept an unresolved curve.
- Every physical value carries a unit type; every point carries a frame type.
- Factories validate invariants and raise one named exception per failure mode.
- No conditional imports, optional-dependency probes, silent fallbacks, skipped
  tests, or expected-failure tests.
- Use Pixi exclusively. Every pytest command includes `-n auto`; affected
  testing additionally uses `--testmon`.
- Do not modify generator behavior to make a reference case pass.
- Do not introduce application fingerprints, frozen-review machinery, or
  identity wrappers.

---

### Task 1: Typed reconstruction geometry

**Files:**

- Create: `benchmarks/held_reference_geometry.py`
- Modify: `benchmarks/errors.py`
- Test: `tests/benchmarks/test_held_reference_geometry.py`

**Interfaces:**

- Consumes: `Point2[WorldXY]`, `Millimetre`, `Radian`, and `ToolRadius` from
  `compas_cgal.adaptive.units`.
- Produces:
  - `PdfPointUnit = NewType("PdfPointUnit", float)`
  - `MillimetresPerPdfPoint = NewType("MillimetresPerPdfPoint", float)`
  - `PdfPoint2.build(x: float, y: float) -> PdfPoint2`
  - `SourceLine.build(start: PdfPoint2, end: PdfPoint2) -> SourceLine`
  - `SourceCubic.build(start: PdfPoint2, control1: PdfPoint2,
    control2: PdfPoint2, end: PdfPoint2) -> SourceCubic`
  - `ReferenceLine.build(start: Point2[WorldXY], end: Point2[WorldXY])`
  - `ReferenceArc.build(start: Point2[WorldXY], end: Point2[WorldXY],
    centre: Point2[WorldXY], sweep: Radian)`
  - `ReferenceBoundary.build(primitives: Sequence[ReferencePrimitive],
    tool_radius: ToolRadius, boundary_stroke_width: Millimetre)`
  - `reconstruct_cubic(source: SourceCubic, transform: SourceToWorld,
    deviation_limit: Millimetre) -> tuple[ReferenceArc, ...]`
  - `project_boundary(boundary: ReferenceBoundary,
    deviation_limit: Millimetre) -> PolygonProjection`

- [x] **Step 1: Add named geometry failures**

Append these independent failures to `benchmarks/errors.py`:

```python
class InvalidPublishedPrimitiveError(BenchmarkError):
    """A published vector primitive is non-finite or degenerate."""


class DisconnectedPublishedBoundaryError(BenchmarkError):
    """Published boundary primitives do not form one closed cycle."""


class UnresolvedPublishedCurveError(BenchmarkError):
    """A published cubic cannot be reconstructed inside its fidelity bound."""


class InvalidReferenceProjectionError(BenchmarkError):
    """A polygon projection violates its declared chord-deviation contract."""
```

- [x] **Step 2: Write RED tests for unit/frame and factory invariants**

Create `tests/benchmarks/test_held_reference_geometry.py` with focused tests:

```python
def test_pdf_point_rejects_non_finite_coordinate() -> None:
    with pytest.raises(InvalidPublishedPrimitiveError):
        PdfPoint2.build(float("nan"), 0.0)


def test_source_line_rejects_identical_endpoints() -> None:
    point = PdfPoint2.build(2.0, 3.0)
    with pytest.raises(InvalidPublishedPrimitiveError):
        SourceLine.build(point, point)


def test_reference_boundary_rejects_disconnected_primitives() -> None:
    first = ReferenceLine.build(_world(0.0, 0.0), _world(1.0, 0.0))
    second = ReferenceLine.build(_world(2.0, 0.0), _world(0.0, 0.0))
    with pytest.raises(DisconnectedPublishedBoundaryError):
        ReferenceBoundary.build(
            (first, second),
            ToolRadius.build(1.0),
            Millimetre(0.1),
        )
```

Run:

```bash
pixi run pytest -- tests/benchmarks/test_held_reference_geometry.py -n auto -q
```

Expected: collection fails because the module and named failures do not exist.

- [x] **Step 3: Implement typed source and reference primitives**

Use frozen domain dataclasses with validating `build(...)` factories. Keep the
union explicit:

```python
ReferencePrimitive: TypeAlias = ReferenceLine | ReferenceArc


@dataclass(frozen=True)
class ReferenceBoundary:
    primitives: tuple[ReferencePrimitive, ...]
    tool_radius: ToolRadius
    boundary_stroke_width: Millimetre

    @classmethod
    def build(
        cls,
        primitives: Sequence[ReferencePrimitive],
        tool_radius: ToolRadius,
        boundary_stroke_width: Millimetre,
    ) -> Self:
        ordered = tuple(primitives)
        _validate_closed_cycle(ordered)
        return cls(ordered, tool_radius, boundary_stroke_width)
```

Raw constructors must only store already-typed values; validation belongs in
the factories. Equality at primitive junctions is direct typed-coordinate
equality, not a tolerance decision.

- [x] **Step 4: Write RED circle-recovery and recursive-pair tests**

Construct one standard cubic approximation of a quarter circle and one cubic
that cannot meet the whole-span circle bound:

```python
def test_quarter_circle_cubic_recovers_one_arc() -> None:
    cubic = _quarter_circle_source()
    arcs = reconstruct_cubic(cubic, _identity_transform(), Millimetre(0.001))
    assert len(arcs) == 1
    assert float(arcs[0].centre.x) == pytest.approx(0.0)
    assert float(arcs[0].centre.y) == pytest.approx(0.0)


def test_non_circular_cubic_splits_into_g1_arc_pairs() -> None:
    arcs = reconstruct_cubic(
        _non_circular_source(),
        _identity_transform(),
        Millimetre(0.01),
    )
    assert len(arcs) >= 2
    for left, right in zip(arcs, arcs[1:]):
        assert left.end == right.start
        assert _end_tangent(left) == pytest.approx(_start_tangent(right))
```

Run the focused module and confirm these tests fail because reconstruction is
absent.

- [x] **Step 5: Implement circle-first recursive reconstruction**

For a source cubic:

1. transform endpoints and controls into `WorldXY`;
2. derive endpoint tangents from control differences;
3. construct the equal-radius circle through both endpoints tangent at the
   start and validate the end tangent;
4. determine sweep direction from the validated endpoint tangents;
5. bound radial separation by recursively subdividing the cubic with de
   Casteljau until each control hull closes the bound;
6. accept the candidate if the closed bound is within the declared limit;
7. otherwise fit an equal-distance G1 biarc, prove its continuous
   correspondence bound, or split the cubic at `t = 1/2` and recurse; and
8. raise `UnresolvedPublishedCurveError` when the named maximum subdivision
   depth is reached without closure.

The subdivision depth is a named structural limit derived from binary halving:
`MAX_RECONSTRUCTION_SUBDIVISIONS = 24`, giving over sixteen million potential
parameter cells before refusal. It prevents non-termination; it is not a
geometric tolerance.

- [x] **Step 6: Write RED projection tests**

```python
def test_projection_closes_and_meets_chord_bound() -> None:
    boundary = _rounded_rectangle_boundary()
    projection = project_boundary(boundary, Millimetre(0.002))
    assert projection.points[0] != projection.points[-1]
    assert projection.observed_deviation <= projection.deviation_limit
    assert len(projection.points) >= 8


def test_projection_rejects_non_positive_bound() -> None:
    with pytest.raises(InvalidReferenceProjectionError):
        project_boundary(_rounded_rectangle_boundary(), Millimetre(0.0))
```

- [x] **Step 7: Implement analytic projection and run gates**

Lines contribute their end point once. Arcs choose the smallest segment count
whose sagitta is within the supplied length bound:

```python
angle_step = 2.0 * math.acos(1.0 - deviation_limit / radius)
segment_count = math.ceil(abs(sweep) / angle_step)
```

Guard the `acos` domain through validated positive radius and a bound no larger
than the radius. Measure the emitted chord sagitta independently before
constructing `PolygonProjection`.

Run:

```bash
pixi run pytest -- tests/benchmarks/test_held_reference_geometry.py -n auto -q
pixi run lint
pixi run types-benchmarks
```

If `types-benchmarks` does not yet exist, add this Pixi task in the same commit:

```toml
types-benchmarks = "mypy --strict --warn-unused-ignores benchmarks/held_reference_geometry.py tests/benchmarks/typecheck/held_reference_contract.py"
```

The type-contract file uses `assert_type` for every public return type.

- [x] **Step 8: Commit Task 1**

```bash
git add benchmarks/errors.py benchmarks/held_reference_geometry.py \
  tests/benchmarks/test_held_reference_geometry.py \
  tests/benchmarks/typecheck/held_reference_contract.py pyproject.toml
git commit -m "feat(bench): reconstruct published curves"
```

---

### Task 2: Publisher-PDF vector extraction

**Files:**

- Create: `tools/held_reference_extractor.py`
- Create: `tests/tools/test_held_reference_extractor.py`
- Modify: `benchmarks/errors.py`
- Modify: `pyproject.toml`

**Interfaces:**

- Consumes: `PdfPoint2`, `SourceLine`, and `SourceCubic` from Task 1.
- Produces:
  - `AffineTransform.build(a, b, c, d, e, f)`
  - `parse_pdf_svg_path(path_data: str, transform: AffineTransform)`
  - `extract_reference_sources(pdf_path: Path) -> tuple[ExtractedCase, ...]`
  - Pixi task `held-reference-extract -- <publisher-pdf>`

- [ ] **Step 1: Add extraction-specific failures**

```python
class UnsupportedPdfBoundaryOperatorError(BenchmarkError):
    """A selected publisher path uses an unsupported drawing operator."""


class MissingPublishedToolCircleError(BenchmarkError):
    """A figure crop contains no unambiguous depicted tool circle."""


class AmbiguousPublishedBoundaryError(BenchmarkError):
    """A figure crop contains more than one valid boundary selection."""
```

- [ ] **Step 2: Write RED parser tests from minimal SVG fragments**

Use literal fragments matching Poppler's `M ... L ...` and `M ... C ...`
output, including the six-value affine matrix:

```python
def test_parser_applies_affine_transform_to_line() -> None:
    primitive = parse_pdf_svg_path(
        "M 1 2 L 3 4",
        AffineTransform.build(2.0, 0.0, 0.0, -2.0, 10.0, 20.0),
    )
    assert primitive == SourceLine.build(
        PdfPoint2.build(12.0, 16.0),
        PdfPoint2.build(16.0, 12.0),
    )


def test_parser_rejects_close_operator() -> None:
    with pytest.raises(UnsupportedPdfBoundaryOperatorError):
        parse_pdf_svg_path("M 0 0 L 1 0 Z", AffineTransform.identity())
```

Run and observe the missing-module RED.

- [ ] **Step 3: Implement the limited path grammar**

Parse only the exact absolute operator forms emitted by the publisher PDF:

- `M x y L x y`
- `M x y` followed by one through four consecutive
  `C x1 y1 x2 y2 x3 y3` operators.

Split compound cubic paths into consecutive `SourceCubic` values. Reject
relative operators, mixed line/cubic compounds, closure operators, malformed
token counts, and non-finite numbers. Apply the SVG affine matrix at ingestion
so no downstream component sees nested transforms.

- [ ] **Step 4: Write RED crop-selection tests**

Fixture SVGs contain green distractors outside the crop, wrong stroke widths,
and disconnected lines inside it. Tests require selection by crop, stroke,
transform family, and cycle continuity rather than colour alone.

```python
def test_crop_rejects_colour_only_distractor(tmp_path: Path) -> None:
    svg = _svg_with_green_distractor_and_closed_boundary()
    selected = select_boundary_paths(svg, _figure5_crop())
    assert len(selected) == 12
    assert all(path.stroke_width == pytest.approx(2.0) for path in selected)
```

- [ ] **Step 5: Implement page conversion and approved crop records**

Run Poppler through `subprocess.run(..., check=True)` with explicit argument
lists:

```python
subprocess.run(
    ["pdftocairo", "-svg", "-f", str(page), "-l", str(page), str(pdf), str(svg)],
    check=True,
    capture_output=True,
    text=True,
)
```

`FigureCrop` records page, crop rectangle, expected boundary family, and
expected depicted-circle count. Their numeric coordinates are source data
copied from the publisher page frame and documented beside each record.

Identify tool circles as closed four-cubic red paths with equal transformed
width and height. Select the isolated circle outside the boundary where
present. Record in-figure start markers as observations only; the crossed-skis
start marker differs from its isolated scale circle and must not be treated as
an independent scale estimate.

Canonicalize only unique degree-one endpoint pairs within the named `1/256 pt`
PDF coordinate quantum. Tests cover the two crossed-skis seams and one
Monstera seam and prove the next-nearest endpoints remain outside that bound.

- [ ] **Step 6: Add the Pixi entry point and integration test**

Add Poppler as a hard Pixi dependency and a forwardable task:

```toml
poppler = "*"
held-reference-extract = { cmd = "python -m tools.held_reference_extractor", description = "Extract Held-Pfeiffer reference vectors" }
```

The integration test receives the local PDF path through a required test
option only in the dedicated extraction command; ordinary tests use committed
minimal SVG fixtures and never download external content.

Run:

```bash
pixi run pytest -- tests/tools/test_held_reference_extractor.py -n auto -q
pixi run lint
pixi run types-benchmarks
```

- [ ] **Step 7: Commit Task 2**

```bash
git add tools/held_reference_extractor.py \
  tests/tools/test_held_reference_extractor.py benchmarks/errors.py \
  pyproject.toml pixi.lock
git commit -m "feat(bench): extract Held vectors"
```

---

### Task 3: Four strict reference cases

**Files:**

- Create: `benchmarks/held_reference_cases.py`
- Create: `benchmarks/data/held_pfeiffer_2025/figure5.json`
- Create: `benchmarks/data/held_pfeiffer_2025/figure8_upper.json`
- Create: `benchmarks/data/held_pfeiffer_2025/figure8_crossed_skis.json`
- Create: `benchmarks/data/held_pfeiffer_2025/figure8_monstera.json`
- Create: `benchmarks/data/held_pfeiffer_2025/README.md`
- Create: `tests/benchmarks/test_held_reference_cases.py`
- Create: `tests/benchmarks/typecheck/held_reference_cases_contract.py`
- Modify: `benchmarks/errors.py`

**Interfaces:**

- Consumes: extracted primitives and reconstruction functions from Tasks 1-2;
  `PocketSpec.build(...)` from `benchmarks.spec`.
- Produces:
  - `HeldReferenceCase.build(...)`
  - `load_held_reference_case(name: str) -> HeldReferenceCase`
  - `load_all_held_reference_cases() -> tuple[HeldReferenceCase, ...]`
  - `HeldReferenceCase.pocket_spec() -> PocketSpec`

- [ ] **Step 1: Add strict case-file failures**

```python
class UnknownHeldReferenceCaseError(BenchmarkError):
    """A requested Held-Pfeiffer reference case is not one of the four cases."""


class MalformedHeldReferenceCaseError(BenchmarkError):
    """A Held-Pfeiffer case document violates its closed schema."""


class UnsupportedHeldReferenceVersionError(BenchmarkError):
    """A Held-Pfeiffer case document uses an unsupported schema version."""
```

- [ ] **Step 2: Write RED schema and loader tests**

```python
EXPECTED_CASES = (
    "figure5",
    "figure8_upper",
    "figure8_crossed_skis",
    "figure8_monstera",
)


def test_exactly_four_reference_cases_load() -> None:
    cases = load_all_held_reference_cases()
    assert tuple(case.name for case in cases) == EXPECTED_CASES


def test_figure7_is_evidence_not_a_case() -> None:
    with pytest.raises(UnknownHeldReferenceCaseError):
        load_held_reference_case("figure7")
    assert load_held_reference_case("figure5").figure7_observation is not None


def test_every_case_builds_the_benchmark_consumer() -> None:
    for case in load_all_held_reference_cases():
        spec = case.pocket_spec()
        assert spec.tool_diameter == pytest.approx(2.0)
        assert spec.tea_cap_deg == pytest.approx(80.0)
        assert spec.holes == ()
```

Also mutate one valid payload at a time to verify rejection of an extra key,
unknown primitive kind, wrong unit, wrong case name, absent projection metric,
and duplicated Figure 7 evidence.

- [ ] **Step 3: Implement the closed JSON schema and typed loader**

Keep the schema literal beside the loader because they evolve together. Set
`additionalProperties: false` at every object level. Parse JSON with
`json.loads`, validate it, then route every primitive through the Task-1
factory. Do not construct domain dataclasses directly.

`HeldReferenceCase` owns publication metadata, the analytic boundary, the
polygon projection, the optional start point, the 80-degree cap, and optional
Figure 7 observation. Its factory verifies that recorded and recomputed
deviations agree within the JSON decimal representation.

- [ ] **Step 4: Generate the four cases from the publisher PDF**

Run:

```bash
pixi run held-reference-extract -- \
  /Users/jelle/Code/CADCAM/compas_cgal_prs/tmp/pdfs/held-pfeiffer-2025/paper.pdf
```

The command writes only the four named JSON documents. Review the emitted
primitive counts, cycle closure, tool-circle scale, observed reconstruction
deviation, and projection deviation before staging them.

- [ ] **Step 5: Add semantic/metamorphic case tests**

For every committed case:

- reconstruct the analytic cycle from JSON;
- confirm direct endpoint continuity and positive signed area;
- confirm the polygon has no self-intersection through `PocketSpec.build`;
- confirm normalizing the already-normalized source is idempotent;
- reverse the analytic cycle and prove the factory restores CCW orientation;
- regenerate the projection at half the deviation limit and confirm its vertex
  count does not decrease; and
- round-trip the JSON through `json.loads` and schema validation.

Run:

```bash
pixi run pytest -- tests/benchmarks/test_held_reference_cases.py -n auto -q
pixi run pytest -- tests/benchmarks/test_held_reference_geometry.py \
  tests/tools/test_held_reference_extractor.py \
  tests/benchmarks/test_held_reference_cases.py -n auto -q
pixi run types-benchmarks
pixi run lint
```

- [ ] **Step 6: Commit Task 3**

```bash
git add benchmarks/held_reference_cases.py benchmarks/errors.py \
  benchmarks/data/held_pfeiffer_2025 \
  tests/benchmarks/test_held_reference_cases.py \
  tests/benchmarks/typecheck/held_reference_cases_contract.py
git commit -m "feat(bench): add Held reference cases"
```

---

### Task 4: Figure 7 cross-check and visual overlays

**Files:**

- Create: `benchmarks/held_reference_figures.py`
- Create: `tests/benchmarks/test_held_reference_figures.py`
- Create: `docs/assets/images/held_reference_figure5.png`
- Create: `docs/assets/images/held_reference_figure8_upper.png`
- Create: `docs/assets/images/held_reference_figure8_crossed_skis.png`
- Create: `docs/assets/images/held_reference_figure8_monstera.png`
- Modify: `tools/held_reference_extractor.py`
- Modify: `benchmarks/data/held_pfeiffer_2025/figure5.json`
- Modify: `pyproject.toml`

**Interfaces:**

- Consumes: four `HeldReferenceCase` values and local publisher PDF.
- Produces:
  - `render_reference_overlay(case: HeldReferenceCase, source: SourceCrop,
    output: Path) -> None`
  - `measure_figure7_observation(...) -> Figure7Observation`
  - Pixi task `held-reference-figures -- <publisher-pdf>`

- [ ] **Step 1: Write RED overlay-content tests**

Avoid screenshot-only tests. Expose the overlay marks before rendering:

```python
def test_overlay_contains_every_evidence_layer() -> None:
    marks = reference_overlay_marks(_figure5_case(), _figure5_source())
    assert {mark.role for mark in marks} == {
        "source",
        "analytic",
        "junction",
        "projection",
        "start",
        "tool",
    }
```

Test PNG dimensions and mode after rendering; do not compare encoded bytes.

- [ ] **Step 2: Implement Figure 7 registration measurement**

Render PDF page 14 through Poppler at a named 600-DPI resolution. Use the plot
axes and the Figure 5 boundary bounds to solve one similarity transform for
each Figure 7 panel. Register the coloured tool-centre samples against the
one-tool-radius inward offset of the Figure 5 analytic boundary. Render the
three panels as shape-only falsification evidence.

Figure 7 has no independent boundary stroke, so it supplies no numeric
boundary-fidelity acceptance value. Missing axes, an empty colour mask, or an
invalid inward offset raises a named error; no manual numeric value enters the
case file.

- [ ] **Step 3: Render the four overlays**

Each PNG contains:

- publisher centreline in muted grey;
- recovered lines/arcs in green;
- primitive junctions as black points;
- polygon projection as a thin dashed blue line;
- red tool circle and start marker where published; and
- a caption with case, figure, `r = 1 mm`, primitive count, projection vertex
  count, and measured maximum deviations.

Use PNG at 2400 pixels on the longer side. Inspect all four images with the
local image viewer; correct any clipping, inverted orientation, missing arc, or
misregistered tool circle before proceeding.

- [ ] **Step 4: Add the Pixi figure task and run gates**

```toml
held-reference-figures = { cmd = "python -m benchmarks.held_reference_figures", description = "Render Held reference reconstruction overlays" }
```

Run:

```bash
pixi run held-reference-figures -- \
  /Users/jelle/Code/CADCAM/compas_cgal_prs/tmp/pdfs/held-pfeiffer-2025/paper.pdf
pixi run pytest -- tests/benchmarks/test_held_reference_figures.py -n auto -q
pixi run affected
pixi run types-benchmarks
pixi run lint
```

- [ ] **Step 5: Commit Task 4**

```bash
git add benchmarks/held_reference_figures.py \
  tools/held_reference_extractor.py \
  tests/benchmarks/test_held_reference_figures.py \
  benchmarks/data/held_pfeiffer_2025/figure5.json \
  docs/assets/images/held_reference_*.png pyproject.toml
git commit -m "docs(bench): render Held reconstructions"
```

---

### Task 5: Benchmark integration and generator qualification

**Files:**

- Create: `tools/held_reference_qualification.py`
- Create: `tests/tools/test_held_reference_qualification.py`
- Create: `docs/held_pfeiffer_reference_pockets.md`
- Create: `docs/benchmarks/held_reference_qualification.md`
- Modify: `benchmarks/figure6.py`
- Modify: `tests/benchmarks/test_figure6.py`
- Modify: `docs/benchmarks.md`
- Modify: `docs/segment_site_mat.md`
- Modify: `mkdocs.yml`
- Modify: `pyproject.toml`

**Interfaces:**

- Consumes: `load_all_held_reference_cases()`, current Figure 6 pipeline, and
  the existing generator without modifications.
- Produces:
  - Figure 6 `reference_pocket()` backed by the reconstructed Figure 5 case;
  - Pixi task `held-reference-qualify`;
  - one concise qualification report naming success or the existing named
    product failure for each case.

- [ ] **Step 1: Write RED Figure 6 consumer tests**

Replace rectangle assertions with the reconstructed case contract:

```python
def test_figure6_uses_the_reconstructed_figure5_pocket() -> None:
    spec = reference_pocket()
    expected = load_held_reference_case("figure5").pocket_spec()
    assert tuple(tuple(point) for point in spec.polygon.points) == tuple(
        tuple(point) for point in expected.polygon.points
    )
    assert spec.tool_diameter == pytest.approx(2.0)
```

Observe RED while `reference_pocket()` still returns `rect_20x12`.

- [ ] **Step 2: Replace only the Figure 6 input seam**

Change `benchmarks.figure6.reference_pocket` to load Figure 5 and rebuild its
`PocketSpec` with the requested cap. Do not change generator parameters,
measurement logic, selection logic, or report semantics in this task.

Run:

```bash
pixi run pytest -- tests/benchmarks/test_figure6.py -n auto -q
```

- [ ] **Step 3: Write RED qualification-report tests**

Use injected generator callables so the renderer is tested without expensive
native runs:

```python
def test_qualification_reports_all_four_cases() -> None:
    outcomes = qualify_cases(_fake_success_generator)
    assert tuple(outcome.case_name for outcome in outcomes) == EXPECTED_CASES
    assert all(outcome.operation_count > 0 for outcome in outcomes)


def test_qualification_preserves_named_product_failure() -> None:
    outcomes = qualify_cases(_generator_raising_known_failure)
    assert outcomes[0].status == "failed"
    assert outcomes[0].failure_type == "UncertifiedZeroGuideEdgeError"
```

Ordinary unexpected exceptions propagate; the qualification command may record
only the explicit existing product failures approved in its closed tuple.

- [ ] **Step 4: Implement and run real qualification**

Add:

```toml
held-reference-qualify = { cmd = "python -m tools.held_reference_qualification", depends-on = ["_editable-rebuild"], description = "Qualify the generator on Held reference pockets" }
```

Run each case serially inside the command so native editable builds are never
concurrent. The command writes
`docs/benchmarks/held_reference_qualification.md` with case, primitive count,
projection count, operation count, elapsed generation time, and either success
or the exact named failure. It makes no comparative performance claim.

Run:

```bash
pixi run held-reference-qualify
```

Retain failures as truthful product evidence. Do not edit case geometry or
generator code in response during this plan.

- [ ] **Step 5: Write the MkDocs reference page**

Document:

- source paper and exact figures;
- four-case scope and Figure 7 reuse;
- one-millimetre normalized tool radius;
- circle-first and recursive G1 arc-pair reconstruction;
- analytic versus polygon-projection authority;
- the four inspected overlay PNGs;
- current qualification outcomes; and
- explicit non-claims about original CAD, physical scale, G-code, cycle-time
  superiority, and Held-Pfeiffer parity.

Add the page to `mkdocs.yml` under Design Notes and link it from
`docs/benchmarks.md` and the Held comparison in `docs/segment_site_mat.md`.

- [ ] **Step 6: Run final focused and repository gates**

```bash
pixi run pytest -- tests/benchmarks/test_held_reference_geometry.py \
  tests/tools/test_held_reference_extractor.py \
  tests/benchmarks/test_held_reference_cases.py \
  tests/benchmarks/test_held_reference_figures.py \
  tests/tools/test_held_reference_qualification.py \
  tests/benchmarks/test_figure6.py -n auto -q
pixi run affected
pixi run types-benchmarks
pixi run lint
pixi run -e docs docs
git diff --check
```

Inspect all four latest PNGs after the final generation command. Scan all new
work for prohibited application-fingerprint and review-freezing mechanisms;
remove any introduction before commit.

- [ ] **Step 7: Commit Task 5**

```bash
git add tools/held_reference_qualification.py \
  tests/tools/test_held_reference_qualification.py \
  benchmarks/figure6.py tests/benchmarks/test_figure6.py \
  docs/held_pfeiffer_reference_pockets.md \
  docs/benchmarks/held_reference_qualification.md \
  docs/benchmarks.md docs/segment_site_mat.md mkdocs.yml pyproject.toml
git commit -m "feat(bench): qualify Held reference pockets"
```

## Final acceptance

- [ ] Four and only four case documents load through the strict schema.
- [ ] Figure 5 contains Figure 7 evidence without duplicating its geometry.
- [ ] Every analytic boundary is closed, CCW, simple, normalized, and inside
  its measured publication-resolution limit.
- [ ] Every polygon projection passes `PocketSpec.build(...)` and records its
  measured deviation.
- [ ] All four source-overlay PNGs pass visual inspection.
- [ ] Figure 6 consumes the reconstructed Figure 5 pocket.
- [ ] Qualification reports the existing generator outcome for every case
  without generator changes.
- [ ] Focused tests, affected tests, strict typing, Ruff, strict MkDocs, and
  diff checks pass.
- [ ] Branch status is clean and every commit records Jelle Feringa as author
  and committer.
