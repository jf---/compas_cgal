# Auditor Convergence P1 Truthful Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one authoritative engagement audit in which every lateral cut has
a native three-way verdict and missing measurement can never become compliance.

**Architecture:** Introduce a focused `engagement_audit` package alongside the
legacy diagnostic. A content-addressed input owns typed geometry and build
identity; one opaque native replay owner decides and depletes the same exact
motion atomically; report construction derives all aggregates from a closed
operation-audit union. Partial arcs use an explicit exact rational-chart
surrogate shared by certification and depletion after P0 dependency closure.

**Tech Stack:** Python 3.12, C++20, CGAL exact kernels, nanobind, COMPAS typed
geometry, CCAN, SHA-256, Pixi, pytest-xdist, pytest-testmon, Ruff, strict mypy.

**Spec:** `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md`

## Global Constraints

- P0 acceptance commit and reconciliation ledger are mandatory inputs.
- Preserve `compas_cgal.engagement` unchanged while the new path is validated.
- Exact cut-plane classification uses typed values and exact binary64 injection;
  no tolerance or operation label may certify non-engagement.
- Every cut-height lateral line/arc/circle reaches a native certifier returning
  `certified`, `cap_exceeded`, or `unresolved`.
- Unsupported geometry raises a named error before stock mutation.
- New public records use validated factories, canonical bytes, and SHA-256.
- `__init__.py` remains minimal and contains no `__all__`.
- Every pytest command uses `-n auto`; changed-code GREEN uses `--testmon`.
- Run Ruff and strict mypy before every Python commit.
- Never alter reference tests, use conditional imports, or add fallback logic.

---

### Task 1: Build identity and immutable audit records

**Files:**

- Create: `src/compas_cgal/engagement_audit/__init__.py`
- Create: `src/compas_cgal/engagement_audit/errors.py`
- Create: `src/compas_cgal/engagement_audit/identity.py`
- Create: `src/compas_cgal/engagement_audit/records.py`
- Create: `tests/engagement_audit/test_identity.py`
- Create: `tests/engagement_audit/test_records.py`
- Create: `tests/adaptive/typecheck/auditor_contract.py`

**Interfaces:**

- Consumes: `ComponentIdentity`, `IdentityDigest`, `CanonicalOperation`,
  `NativeMotionVerdict`, canonical CCAN encoders.
- Produces: `BuildIdentity.build(...)`, `MeasuredOperationAudit.build(...)`,
  `NonEngagingOperationAudit.build(...)`, `OperationAudit`.

- [ ] **Step 1: Write RED build-identity tests**

Add tests asserting exact-type, cardinality, ordering, mutation, and digest
contracts:

```python
def test_build_identity_binds_every_component_and_source_digest() -> None:
    identity = _build_identity()

    assert identity.components == tuple(sorted(identity.components, key=lambda item: item.component_domain))
    assert identity.digest == hashlib.sha256(identity.canonical_bytes).digest()
    assert bytes(identity.native_source_tree_digest) in identity.canonical_bytes
    assert identity.pixi_lock_digest in identity.canonical_bytes


def test_build_identity_rejects_duplicate_component_domains() -> None:
    component = _component(b"stock")

    with pytest.raises(DuplicateBuildComponentError, match="stock"):
        BuildIdentity.build(
            components=(component, component),
            native_source_tree_digest=_digest(b"native"),
            python_source_tree_digest=_digest(b"python"),
            pixi_lock_digest=_digest(b"pixi"),
        )
```

Mutate every digest length, pass subclasses, reorder components, and change one
component version. Reordering must canonicalize identically; every semantic
mutation must change the digest.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/engagement_audit/test_identity.py -n auto -q
```

Expected: import failure because `engagement_audit.identity` does not exist.

- [ ] **Step 3: Implement `BuildIdentity`**

Define exact digest aliases and one immutable record:

```python
PythonSourceTreeDigest = NewType("PythonSourceTreeDigest", bytes)
PixiLockDigest = NewType("PixiLockDigest", bytes)


@dataclass(frozen=True)
class BuildIdentity:
    components: tuple[ComponentIdentity, ...]
    native_source_tree_digest: NativeSourceTreeDigest
    python_source_tree_digest: PythonSourceTreeDigest
    pixi_lock_digest: PixiLockDigest

    @classmethod
    def build(
        cls,
        *,
        components: tuple[ComponentIdentity, ...],
        native_source_tree_digest: NativeSourceTreeDigest,
        python_source_tree_digest: PythonSourceTreeDigest,
        pixi_lock_digest: PixiLockDigest,
    ) -> Self:
        normalized_components = _validate_and_sort_components(components)
        _require_sha256_digest(native_source_tree_digest, "native source tree")
        _require_sha256_digest(python_source_tree_digest, "Python source tree")
        _require_sha256_digest(pixi_lock_digest, "pixi.lock")
        return cls(
            components=normalized_components,
            native_source_tree_digest=native_source_tree_digest,
            python_source_tree_digest=python_source_tree_digest,
            pixi_lock_digest=pixi_lock_digest,
        )
```

The factory requires a nonempty exact `tuple`, exact `ComponentIdentity`
members, unique domains, and three 32-byte digests. Canonical bytes use a new
`build-identity-v1` tagged record. Raw construction repeats validation in
`__post_init__`; bypass is never unsafe.

- [ ] **Step 4: Write RED operation-record tests**

Use the existing native verdict literal and separate records:

```python
def test_non_engaging_record_has_no_certificate_verdict() -> None:
    record = NonEngagingOperationAudit.build(
        operation_index=0,
        operation_digest=_digest(b"retract"),
        reason="vertical_retract",
    )

    assert not hasattr(record, "verdict")


def test_measured_record_rejects_foreign_verdict() -> None:
    with pytest.raises(InvalidMotionVerdictError, match="foreign"):
        MeasuredOperationAudit.build(
            operation_index=0,
            operation_digest=_digest(b"cut"),
            verdict="foreign",  # type: ignore[arg-type]
            max_tea=Radian(0.5),
            station_count=1,
            pre_motion_stock_lineage=_digest(b"stock"),
            motion_certificate_digest=_digest(b"certificate"),
        )
```

- [ ] **Step 5: Implement the closed record union**

Define:

```python
NonEngagingReason: TypeAlias = Literal[
    "vertical_plunge",
    "vertical_retract",
    "clearance_transport",
]
OperationAudit: TypeAlias = MeasuredOperationAudit | NonEngagingOperationAudit
```

`MeasuredOperationAudit` validates `NativeMotionVerdict`, typed `Radian`,
positive station/event count, 32-byte lineages, and complete canonical bytes.
`NonEngagingOperationAudit` has no max TEA, station count, or verdict.

- [ ] **Step 6: Run GREEN, typing, and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/engagement_audit/test_identity.py tests/engagement_audit/test_records.py -n auto --testmon -q
git diff --check
git add src/compas_cgal/engagement_audit tests/engagement_audit tests/adaptive/typecheck/auditor_contract.py
git commit -m "feat(audit): define truthful records"
```

### Task 2: Authenticate the audit input and classify geometry

**Files:**

- Create: `src/audit_classification_2.h`
- Create: `src/audit_classification_2.cpp`
- Create: `src/compas_cgal/engagement_audit/operation_identity.py`
- Create: `src/compas_cgal/engagement_audit/input.py`
- Create: `src/compas_cgal/engagement_audit/classification.py`
- Create: `tests/engagement_audit/test_native_classification.py`
- Create: `tests/engagement_audit/test_input.py`
- Create: `tests/engagement_audit/test_classification.py`
- Modify: `CMakeLists.txt`
- Modify: `pyproject.toml`
- Modify: `src/stock_2.cpp`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/adaptive/units.py`
- Modify: `src/compas_cgal/engagement_audit/errors.py`
- Modify: `src/compas_cgal/engagement_audit/records.py`

**Interfaces:**

- Consumes: `CanonicalRingV1`, `CutPlane`, `ToolRadius`, `EngagementCap`,
  `ToolpathOperation`, `BuildIdentity`.
- Produces: `EngagementAuditInput.build(...)`, `classify_operation(..., *,
  operation_index: int) -> AuthenticatedOperation`.

- [ ] **Step 1: Write RED input-boundary tests**

Require empty, multi-depth, off-plane, ramped, nonfinite, and mislabeled cases:

```python
def test_empty_toolpath_cannot_audit_clean() -> None:
    with pytest.raises(EmptyToolpathAuditError, match="at least one operation"):
        _audit_input(operations=())


def test_cut_plane_is_input_not_inferred_from_lower_motion() -> None:
    operations = (_cut_at_z(0.0), _cut_at_z(-1.0))

    with pytest.raises(MultipleCutPlaneError, match="-1"):
        _audit_input(operations=operations)
```

Mutation tests change ring order, holes, cut/clearance Z, tool, cap, operation
order, and build identity; every change must alter the input digest or raise.

- [ ] **Step 2: Run input RED**

```bash
pixi run pytest -- tests/engagement_audit/test_input.py -n auto -q
```

Expected: missing module/API failure.

- [ ] **Step 3: Implement `EngagementAuditInput`**

The factory signature is exactly:

```python
@classmethod
def build(
    cls,
    *,
    design_boundary: CanonicalRingV1,
    holes: tuple[CanonicalRingV1, ...],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    engagement_cap: EngagementCap,
    operations: tuple[ToolpathOperation, ...],
    build_identity: BuildIdentity,
) -> Self:
    classified_operations, stream_digest = _classify_operations(operations, cut_plane)
    return cls(
        design_boundary=design_boundary,
        holes=holes,
        cut_plane=cut_plane,
        tool_radius=tool_radius,
        engagement_cap=engagement_cap,
        operations=classified_operations,
        operation_stream_digest=stream_digest,
        build_identity=build_identity,
    )
```

Canonical bytes bind every field and an ordered operation-stream digest. Use
one operation canonicalizer in `operation_identity.py`; do not parse warning
strings or derive cut Z from operations. `build(...)` is the sole mutable
COMPAS ingress. It canonicalizes the source stream, sends its geometry to the
native Epeck classifier, and retains only immutable classified operations plus
the ordered source digest. No downstream consumer receives caller-owned
`ToolpathOperation` objects.

- [ ] **Step 4: Write RED classification tests**

```python
def test_retract_label_cannot_hide_cut_height_lateral_motion() -> None:
    operation = _horizontal_line(z=0.0, kind=OperationType.RETRACT)

    with pytest.raises(ContradictoryOperationRoleError, match="RETRACT"):
        classify_operation(operation, _cut_plane(), operation_index=0)


def test_clearance_transport_requires_both_endpoints_on_clearance_plane() -> None:
    operation = _horizontal_line(z=5.0, end_z=4.0, kind=OperationType.LINK)

    with pytest.raises(UnsupportedAuditGeometryError, match="clearance"):
        classify_operation(operation, _cut_plane(), operation_index=0)
```

Cover exact vertical plunge/retract, exact horizontal clearance transport,
cut-height line, arc, circle, mixed-Z XY ramp, and unsupported geometry.

- [ ] **Step 5: Implement geometry-derived classification**

Define opaque native motion values and authenticated Python carriers:

```python
SupportedLateralMotion: TypeAlias = (
    AuditSegmentMotion2 | AuditCircleMotion2 | AuditArcMotion2
)

AuthenticatedOperation: TypeAlias = (
    AuthenticatedLateralOperation
    | AuthenticatedPlungeOperation
    | AuthenticatedNonEngagingOperation
)
```

Add typed `Point3[WorldXYZ]` and `Direction3[WorldXYZ]` ingress primitives.
Each mutable COMPAS operation is read once into a frozen, typed capture;
canonical bytes and the native call consume that same capture. `_stock_2`
exact-injects the finite binary64 values into Epeck and owns XY coincidence,
frame validity, plane, sweep-sign, orientation, and role decisions. It returns
one of six non-constructible nanobind values backed by Epeck, never a label or
reconstructive coordinate tuple. Arc phase crosses one named native
transcendental seam whose native version is bound into audit identity, and
remains opaque thereafter; no Python `point_at`, `atan2`, subtraction, or
`abs` builds certifier geometry.

A native-proved plunge is distinct from a non-mutating retract or clearance
transport. `AuthenticatedPlungeOperation` retains only the opaque native
plunge; Task 4 passes that value directly to native depletion without exposing
or reconstructing its endpoint in Python. Circle and arc values retain their
exact-injected guide radius internally for the same reason. No depletion policy
is chosen in Task 2. The explicit operation index is the ordered stream ordinal,
never repeatable `ToolpathOperation.path_index`.

Ingress uses a finite-only, millimetre-bearing radius observation. Native
`CGAL::sign` owns the positive-radius decision. Each authenticated carrier has
a distinct versioned canonical encoding and digest over operation index, source
operation digest, and its closed native classification tag.

- [ ] **Step 6: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-audit
pixi run types-adaptive
pixi run pytest -- tests/engagement_audit/test_input.py tests/engagement_audit/test_classification.py tests/engagement_audit/test_native_classification.py tests/adaptive/test_units.py -n auto --testmon -q
pixi run -e docs docs
git diff --check
git add CMakeLists.txt pyproject.toml src/audit_classification_2.* src/stock_2.cpp src/compas_cgal/_stock_2.pyi src/compas_cgal/engagement_audit src/compas_cgal/adaptive/units.py tests/engagement_audit tests/adaptive/test_units.py tests/adaptive/typecheck/auditor_contract.py docs/engagement_audit.md mkdocs.yml
git commit -m "feat(audit): authenticate motion input"
```

### Task 3: Define and deplete the exact partial-arc surrogate

**Files:**

- Create: `src/exact_circle_chart_2.h`
- Create: `src/exact_circle_chart_2.cpp`
- Create: `src/audit_digest_2.h`
- Create: `src/audit_arc_motion_2.h`
- Create: `src/audit_arc_motion_2.cpp`
- Modify: `src/audit_classification_2.h`
- Modify: `src/audit_classification_2.cpp`
- Modify: `src/exact_depletion_2.h`
- Modify: `src/exact_depletion_2.cpp`
- Modify: `src/stock_2.h`
- Modify: `src/stock_2.cpp`
- Modify: `src/continuous_tea_2/parameter_charts.cpp`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/engagement_audit/classification.py`
- Modify: `src/compas_cgal/engagement_audit/errors.py`
- Modify: `CMakeLists.txt`
- Modify: `pyproject.toml`
- Create: `tests/native/test_audit_arc_motion_2.cpp`
- Create: `tests/native/test_exact_arc_depletion_2.cpp`
- Modify: `tests/engagement_audit/test_native_classification.py`
- Modify: `tests/engagement_audit/test_classification.py`
- Modify: `tests/engagement_audit/test_input.py`
- Modify: `tests/adaptive/test_exact_depletion.py`
- Modify: `tests/adaptive/typecheck/auditor_contract.py`
- Modify: `docs/engagement_audit.md`

**Interfaces:**

- Consumes: Task 2 opaque arc ingress, Epeck, the frozen four-quarter rational
  atlas, `Stock2.clone()`, and `DepletionPolicy` values at the later call seam.
- Produces: a canonical nonconstructible `AuditArcMotion2`, one shared exact
  quarter-chart evaluator, `Stock2::subtract_exact_arc(...)`, and a validated
  non-cyclic partial-arc `ExactArcDepletionTrace2`.

The exact types are fixed as private-constructor domain values:

```cpp
template <class Domain>
class AuditDigest2 {
public:
    static AuditDigest2 from_bytes(std::string bytes);
    const std::string& bytes() const noexcept;
private:
    explicit AuditDigest2(std::string bytes);
};

class ExactArcChartInterval2 {
public:
    static ExactArcChartInterval2 build(
        int chart,
        const Epeck::FT& start_parameter,
        const Epeck::FT& end_parameter,
        bool increasing,
        bool owns_start_seam,
        bool owns_end_seam);
};

struct NativeMotionDigestDomain;
using NativeMotionDigest2 = AuditDigest2<NativeMotionDigestDomain>;

class AuditArcMotion2 {
public:
    static AuditArcMotion2 build(
        const EPoint& center,
        const EVector& zero_phase,
        const Epeck::FT& guide_radius,
        double authored_start_angle,
        double authored_end_angle,
        bool clockwise,
        const Epeck::FT& cut_z);
    const NativeMotionDigest2& digest() const;
};

class ExactArcDepletionTrace2 {
public:
    static ExactArcDepletionTrace2 build(
        const AuditArcMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit,
        std::vector<ExactCircleChartParameter2> parameters);
};

struct ExactArcDepletionConstruction2 {
    std::vector<EPoint> centers;
    ExactArcDepletionTrace2 trace;
};

ExactArcDepletionConstruction2 construct_exact_arc_depletion(
    const AuditArcMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);
```

- [x] **Step 1: Write RED arc-seam and identity tests**

Require awkward non-quadrant angles, exact quarter seams, negative and
multi-turn start angles, CW/CCW minor and major arcs, and a full-turn surrogate
to produce deterministic opaque motions. Assert exact start/end incidence,
canonical interval order, distinct-endpoint preservation, strategy version,
and 32-byte motion digest. Changing authored angle, direction, exact chart
parameter, or strategy version must change the digest.

- [x] **Step 2: Run the arc-motion RED gate**

```bash
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- tests/engagement_audit/test_native_classification.py tests/engagement_audit/test_classification.py tests/engagement_audit/test_input.py -n auto -q
```

Expected: the current independent trigonometric phase/radius carrier lacks
canonical chart intervals and native motion identity.

- [x] **Step 3: Extract one exact quarter-chart evaluator**

Move the rational Pythagorean evaluator currently private to
`exact_depletion_2.cpp` into `exact_circle_chart_2.*`. It accepts exact Epeck
parameters, not only `size_t` fractions. Make full-circle depletion and
`continuous_tea_2` use the shared chart definition; formula duplication is a
test failure. Preserve the four frozen chart ids and maps exactly.

- [x] **Step 4: Implement the versioned angle-to-chart seam**

Inside `AuditArcMotion2.build(...)`, observe the authored binary64 subtraction
`end_angle - start_angle` once, require it to be finite, exact-inject that one
signed sweep, require nonzero magnitude no larger than the injected full turn,
and prove its sign agrees with orientation. Preserve exact full-turn state
before normalizing endpoints into `[0, tau)` using the named
`audit-arc-quarter-chart-binary64-v2` seam. Choose quadrants by exact comparison
against exact-injected `0`, `pi/2`, `pi`, `3*pi/2`, and `tau`. Map exact seams
structurally to parameter zero of the owning chart; otherwise exact-inject
`tan(local_angle / 2)`. Build the ordered clipped interval sequence and fail
with named exceptions if rounding collapses distinct requested endpoints or if
orientation, seam ownership, extent, radius incidence, or canonical order is
inconsistent. Canonical bytes bind exact parameters and the seam version.
Decompose each Epeck rational with `CGAL::Fraction_traits`, normalize a positive
denominator, and reuse the existing `ExactRational2` canonical encoding; never
serialize through decimal text or `to_double`. Expose only the digest to Python;
chart coordinates remain opaque.

- [x] **Step 5: Write RED structural depletion tests**

```cpp
CHECK(exact_arc_point_is_incident(motion, centers.front()));
CHECK(exact_arc_point_is_incident(motion, centers.back()));
CHECK(exact_arc_structural_density_holds(
    motion, max_chord, construction.trace.parameters()));
CHECK(construction.trace.matches_exact_inputs(
    tool_radius, max_chord, center_count_limit));
```

Cover every chart seam; CW/CCW minor, major, and full-turn surrogates; rational
rotation, translation, and scale; reordered/missing parameters; complement
confusion; nextafter seam values; an off-guide center; nonpositive policy;
`chord_bound >= tool_radius`; and an unbuildable center count. Every refusal
must leave stock exactly unchanged. Assert the authoritative symbols never
call legacy `subtract_arc_sweep`.

- [x] **Step 6: Implement exact-on-surrogate depletion**

Refine every clipped interval dyadically over its exact endpoint parameters
until consecutive squared chords satisfy the exact bound. Check the aggregate
center limit before allocation. Reconstruct all centers through the shared
chart evaluator; validate direction, seam ownership, anchors, exact incidence,
non-cyclic density, and exact `0 < chord_bound < tool_radius`; then subtract
the union of full-radius exact disks from a cloned stock, validate the trace,
and swap atomically.

```cpp
ExactArcDepletionTrace2 Stock2::subtract_exact_arc(
    const AuditArcMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);
```

Exact incidence proves the disk union is a subset of the declared surrogate
sweep. The strict chord/tool relation proves adjacent disk overlap only; do not
claim a quantitative retained-sliver theorem.

Add a forward-executing Pixi `audit-native` task whose private configure,
build, and run dependencies compile and execute the focused native audit gate.
Task 4 extends the same CMake gate; it does not introduce a second native test
runner.

- [x] **Step 7: Run Task 3 GREEN and commit**

```bash
pixi run audit-native
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- tests/engagement_audit/test_native_classification.py tests/engagement_audit/test_classification.py tests/engagement_audit/test_input.py tests/adaptive/test_exact_depletion.py -n auto --testmon -q
pixi run lint
pixi run types-audit
pixi run -e docs docs
git diff --check
git add CMakeLists.txt pyproject.toml src/audit_digest_2.h src/exact_circle_chart_2.* src/audit_arc_motion_2.* src/audit_classification_2.* src/exact_depletion_2.* src/stock_2.* src/continuous_tea_2/parameter_charts.cpp src/compas_cgal/_stock_2.pyi src/compas_cgal/engagement_audit/classification.py src/compas_cgal/engagement_audit/errors.py tests/native/test_audit_arc_motion_2.cpp tests/native/test_exact_arc_depletion_2.cpp tests/engagement_audit/test_native_classification.py tests/engagement_audit/test_classification.py tests/engagement_audit/test_input.py tests/adaptive/test_exact_depletion.py tests/adaptive/typecheck/auditor_contract.py docs/engagement_audit.md
git commit -m "feat(audit): deplete exact arc"
```

### Task 4: Close the native certify-and-deplete transaction

#### Source adjudication and delivery slices

Commit `73d5372` is theorem, falsifier, and refinement-strategy input only. Its
`CertifiedTea::cap_certified == false` conflates an observed violation with
floor/depth exhaustion and an inconclusive swept bound. Its guarded station
checks also use a stricter cap than the authored policy. Neither condition may
be translated to `CAP_EXCEEDED`, and its double trigonometric arc centers may
not replace Task 3's exact rational-chart motion.

Task 4 lands as three linear, independently reviewed commits:

1. **Identity:** cap/policy sealing, six opaque motion digests, native stock and
   request identity, and input/carrier schema v2.
2. **Decision:** exact witness-bearing three-way adapters for segment, circle,
   and the Task 3 arc surrogate, plus content-addressed depletion witnesses.
3. **Transaction:** one-clone replay state, lineage/result derivation, opaque
   bindings, and finalization.

The closed audit verdict is `AuditTeaVerdict2 {CERTIFIED, CAP_EXCEEDED,
UNRESOLVED}`. Existing `ContinuousTeaVerdict::UNRESOLVED_DEGENERACY` maps into
the broader audit `UNRESOLVED`; the old enum is not renamed or broadened.

**Files:**

- Create: `src/audit_motion_result_2.h`
- Create: `src/audit_motion_result_2.cpp`
- Create: `src/audit_motion_identity_2.h`
- Create: `src/audit_motion_identity_2.cpp`
- Create: `src/audit_policy_2.h`
- Create: `src/audit_policy_2.cpp`
- Create: `src/audit_stock_identity_2.h`
- Create: `src/audit_stock_identity_2.cpp`
- Create: `src/audit_request_identity_2.h`
- Create: `src/audit_request_identity_2.cpp`
- Create: `src/audit_strategy_versions_2.h`
- Create: `src/audit_strategy_versions_2.cpp`
- Create: `src/audit_depletion_witness_2.h`
- Create: `src/audit_depletion_witness_2.cpp`
- Create: `src/audit_lineage_2.h`
- Create: `src/audit_lineage_2.cpp`
- Create: `src/audit_certification_2.h`
- Create: `src/audit_certification_2.cpp`
- Create: `src/audit_replay_2.h`
- Create: `src/audit_replay_2.cpp`
- Create: `src/audit_replay_bindings_2.h`
- Create: `src/audit_replay_bindings_2.cpp`
- Create: `src/stock_exact_depletion_2.h`
- Create: `src/stock_exact_depletion_2.cpp`
- Create: `src/exact_disk_region_2.h`
- Create: `src/exact_disk_region_2.cpp`
- Modify: `src/audit_digest_2.h`
- Modify: `src/audit_arc_motion_2.cpp`
- Modify: `src/exact_depletion_2.cpp`
- Create: `src/continuous_tea_2/arc_oracle.h`
- Create: `src/continuous_tea_2/arc_oracle.cpp`
- Modify: `src/audit_classification_2.h`
- Modify: `src/audit_classification_2.cpp`
- Modify: `src/engagement_2.h`
- Modify: `src/engagement_2.cpp`
- Modify: `src/stock_2.h`
- Modify: `src/stock_2.cpp`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/adaptive/canonical.py`
- Modify: `src/compas_cgal/engagement_audit/classification.py`
- Modify: `src/compas_cgal/engagement_audit/input.py`
- Modify: `src/compas_cgal/engagement_audit/records.py`
- Modify: `src/compas_cgal/engagement_audit/errors.py`
- Create: `src/compas_cgal/engagement_audit/digests.py`
- Modify: `CMakeLists.txt`
- Modify: `pyproject.toml`
- Create: `tests/native/test_audit_identity_2.cpp`
- Modify: `tests/native/test_exact_arc_depletion_2.cpp`
- Create: `tests/native/test_audit_replay_2.cpp`
- Create: `tests/native/test_audit_certification_2.cpp`
- Modify: `tests/adaptive/test_canonical.py`
- Modify: `tests/engagement_audit/test_native_classification.py`
- Modify: `tests/engagement_audit/test_records.py`
- Modify: `tests/engagement_audit/test_input.py`
- Create: `tests/engagement_audit/test_native_request_identity.py`
- Create: `tests/engagement_audit/test_native_replay.py`
- Create: `tests/engagement_audit/test_false_arc_certificate.py`
- Create: `tests/test_false_certificate.py`
- Create: `tests/test_growth_bound.py`
- Modify: `tests/adaptive/typecheck/auditor_contract.py`
- Modify: `docs/engagement_audit.md`

**Interfaces:**

- Consumes: Task 2 opaque segment/circle/plunge/non-engaging motions, Task 3's
  frozen exact arc/depletion API, P0's certifier dependency closure,
  `DepletionPolicy`, the event-exact segment/full-circle oracle core, and
  `Stock2` clone/equality operations.
- Produces: input schema v2, nonconstructible `AuditReplay2`, and separate
  opaque native lateral, plunge, and non-engaging result values.

The native boundary is fixed before implementation:

```cpp
struct AuditInputDigestDomain;
struct AuditNativeStockDigestDomain;
struct AuditStockStateDigestDomain;
struct AuditNativeRequestDigestDomain;
struct AuthenticatedOperationDigestDomain;
struct AuditPolicyDigestDomain;
struct NativeDecisionDigestDomain;
struct DepletionWitnessDigestDomain;
struct StockLineageDigestDomain;
struct AuditResultDigestDomain;

using AuditInputDigest2 = AuditDigest2<AuditInputDigestDomain>;
using AuditNativeStockDigest2 = AuditDigest2<AuditNativeStockDigestDomain>;
using AuditStockStateDigest2 = AuditDigest2<AuditStockStateDigestDomain>;
using AuditNativeRequestDigest2 = AuditDigest2<AuditNativeRequestDigestDomain>;
using AuthenticatedOperationDigest2 = AuditDigest2<AuthenticatedOperationDigestDomain>;
using AuditPolicyDigest2 = AuditDigest2<AuditPolicyDigestDomain>;
using NativeDecisionDigest2 = AuditDigest2<NativeDecisionDigestDomain>;
using DepletionWitnessDigest2 = AuditDigest2<DepletionWitnessDigestDomain>;
using StockLineageDigest2 = AuditDigest2<StockLineageDigestDomain>;
using AuditResultDigest2 = AuditDigest2<AuditResultDigestDomain>;

class AuditCapObservation2;

class AuditPolicy2 {
public:
    static AuditPolicy2 build(
        const AuditCapObservation2& engagement_cap,
        const Epeck::FT& tool_radius_mm,
        const Epeck::FT& depletion_chord_bound_mm,
        std::size_t center_count_limit);
    const AuditPolicyDigest2& digest() const;

private:
    struct Storage;
    explicit AuditPolicy2(Storage storage);
};

class AuditCapObservation2 {
public:
    static AuditCapObservation2 build(
        double authored_radians,
        double supplied_chord_ratio);

private:
    struct Storage;
    explicit AuditCapObservation2(Storage storage);
};

enum class AuditTeaVerdict2 {
    CERTIFIED,
    CAP_EXCEEDED,
    UNRESOLVED,
};

class AuditReportedRadian2 {
public:
    static AuditReportedRadian2 build(double reporting_radians);
    double value() const noexcept;

private:
    explicit AuditReportedRadian2(double reporting_radians);
};

enum class AuditDepletionKind2 {
    SEGMENT,
    FULL_CIRCLE,
    ARC,
    PLUNGE,
};

class AuditDepletionWitness2 {
public:
    AuditDepletionKind2 kind() const noexcept;
    const AuditStockStateDigest2& pre_stock_digest() const noexcept;
    const AuditStockStateDigest2& post_stock_digest() const noexcept;
    const DepletionWitnessDigest2& digest() const noexcept;

private:
    struct Storage;
    explicit AuditDepletionWitness2(Storage storage);
    friend AuditDepletionWitness2 apply_audit_segment_depletion_to_trial(
        const Stock2&, Stock2&, const AuditSegmentMotion2&, const AuditPolicy2&);
    friend AuditDepletionWitness2 apply_audit_circle_depletion_to_trial(
        const Stock2&, Stock2&, const AuditCircleMotion2&, const AuditPolicy2&);
    friend AuditDepletionWitness2 apply_audit_arc_depletion_to_trial(
        const Stock2&, Stock2&, const AuditArcMotion2&, const AuditPolicy2&);
    friend AuditDepletionWitness2 apply_audit_plunge_depletion_to_trial(
        const Stock2&, Stock2&, const AuditVerticalPlunge2&, const AuditPolicy2&);
};

AuditDepletionWitness2 apply_audit_segment_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_circle_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_arc_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_plunge_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy);

// Declared only by exact_disk_region_2.h and defined only by
// exact_disk_region_2.cpp.
Gps build_exact_disk_union_region_2(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius);

class AuditLateralResult2 {
public:
    static AuditLateralResult2 build(
        AuditTeaVerdict2 verdict,
        const AuditReportedRadian2& reported_max_tea,
        std::size_t evidence_count,
        const AuthenticatedOperationDigest2& operation_digest,
        const NativeDecisionDigest2& decision_digest,
        const DepletionWitnessDigest2& depletion_digest,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);

private:
    struct Storage;
    explicit AuditLateralResult2(Storage storage);
};

class AuditPlungeResult2 {
public:
    static AuditPlungeResult2 build(
        const AuthenticatedOperationDigest2& operation_digest,
        const DepletionWitnessDigest2& depletion_digest,
        const StockLineageDigest2& pre_lineage,
        const StockLineageDigest2& post_lineage);

private:
    struct Storage;
    explicit AuditPlungeResult2(Storage storage);
};

enum class AuditNonEngagingReason2 {
    VERTICAL_RETRACT,
    CLEARANCE_TRANSPORT,
};

class AuditNonEngagingResult2 {
public:
    static AuditNonEngagingResult2 build(
        const AuthenticatedOperationDigest2& operation_digest,
        AuditNonEngagingReason2 reason,
        const StockLineageDigest2& unchanged_lineage);

private:
    struct Storage;
    explicit AuditNonEngagingResult2(Storage storage);
};

class AuditReplayCompletion2 {
public:
    static AuditReplayCompletion2 build(
        const AuditNativeRequestDigest2& request_digest,
        std::size_t operation_count,
        const StockLineageDigest2& terminal_lineage);

private:
    struct Storage;
    explicit AuditReplayCompletion2(Storage storage);
};

class AuditNativeStockIdentity2 {
public:
    static AuditNativeStockIdentity2 build(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes);
    const std::string& canonical_bytes() const noexcept;
    const AuditNativeStockDigest2& digest() const noexcept;

private:
    struct Storage;
    explicit AuditNativeStockIdentity2(Storage storage);
};

class AuditNativeRequestIdentity2 {
public:
    const std::string& canonical_bytes() const noexcept;
    const AuditNativeRequestDigest2& digest() const noexcept;

private:
    struct Storage;
    explicit AuditNativeRequestIdentity2(Storage storage);
    friend class AuditReplay2;
    friend AuditNativeRequestIdentity2 build_audit_native_request_identity(
        const AuditNativeStockIdentity2&,
        const AuditPolicy2&,
        const std::vector<NativeMotionDigest2>&);
};

const std::string& audit_native_decision_contract_version();
const std::string& audit_native_depletion_contract_version();

class AuditReplay2 {
public:
    static AuditReplay2 build(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        const AuditInputDigest2& input_digest,
        const AuditNativeRequestIdentity2& native_request,
        const AuditPolicy2& policy,
        std::vector<AuthenticatedOperationDigest2> authenticated_operation_digests);

private:
    struct RequestState;
    struct MutableState;
    AuditReplay2(RequestState request, MutableState state);
};

AuditReplayCompletion2 finish_audit_replay(AuditReplay2&);

AuditLateralResult2 audit_deplete_segment(
    AuditReplay2&, const AuditSegmentMotion2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditLateralResult2 audit_deplete_circle(
    AuditReplay2&, const AuditCircleMotion2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditLateralResult2 audit_deplete_arc(
    AuditReplay2&, const AuditArcMotion2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditPlungeResult2 deplete_audit_plunge(
    AuditReplay2&, const AuditVerticalPlunge2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditNonEngagingResult2 record_audit_retract(
    AuditReplay2&, const AuditVerticalRetract2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);
AuditNonEngagingResult2 record_audit_clearance(
    AuditReplay2&, const AuditClearanceTransport2&,
    const AuthenticatedOperationDigest2& authenticated_operation_digest);

// Nanobind-only ingress: digest bytes use their named external domains;
// request and policy remain opaque native values.
AuditReplay2 begin_audit_replay(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    std::string input_digest_bytes,
    const AuditNativeRequestIdentity2& native_request,
    const AuditPolicy2& policy,
    std::vector<std::string> authenticated_operation_digest_bytes);
```

`AuditCapObservation2::build(...)` validates the two finite binary64 cap
representations, recomputes the repository's native
`4*sin^2(theta/2)` observation exactly once, and requires bit-identical equality
with the supplied surrogate before exact injection. `AuditPolicy2::build(...)`
then checks exact positive tool/chord values, exact
`chord_bound < tool_radius`, and a positive center limit. It is immutable and
stored by `AuditReplay2`; per-motion calls cannot change policy mid-stream.

Each lateral result carries the closed verdict, reporting-only maximum TEA,
positive evidence count, native decision digest, depletion-witness digest,
pre/post lineage, authenticated-operation digest, and canonical result digest.
Plunge results omit TEA and decision fields. Retract/clearance results prove
unchanged lineage.
All result classes and `AuditReplay2` are read-only and nonconstructible from
Python. Every `authenticated_operation_digest` parameter is the carrier digest,
not its nested mutable-ingress source digest.

The native module exposes one `build_audit_policy(...)` nanobind factory that
exact-injects the five scalar inputs and calls `AuditPolicy2::build(...)`; raw
`AuditPolicy2` construction is unavailable. `AuditReplayCompletion2` binds the
verified native request digest, operation count, and terminal lineage.

`AuditNativeStockIdentity2::build(...)` exact-injects the actual matrix values,
uses Epeck for degeneracy and winding decisions, canonicalizes ring rotation
and orientation, sorts canonical holes, and encodes coordinates through the
shared CCAN binary64 encoder. Cross-language golden tests require byte-identical
stock identity to `CanonicalRingV1`. `AuditNativeRequestIdentity2` owns the
canonical stock identity, policy digest, ordered typed motion digests, canonical
bytes, and request digest. Its Python factory accepts a tuple of the six opaque
motion classes. The binding may use a closed internal variant to visit them,
but no generic motion union or raw motion-digest sequence crosses Python.
Replay consumes the opaque request identity and retains its expected motion-
digest sequence internally.

`begin_audit_replay(...)` is the sole Python factory for the nonconstructible
owner. It accepts raw bytes only for the external `AuditInputDigest2` and
ordered `AuthenticatedOperationDigest2` domains, parses each through its named
ingress, rebuilds stock and request identity from the actual rings, opaque
policy, and request-retained motion sequence, and rejects any mismatch before
constructing state. Native request digests are generated-only; Python may store
their exposed bytes in input v2 but cannot manufacture the typed native value.

The five non-arc motion carriers become private-constructor classes with named
validating factories that derive their digest from canonical exact geometry,
plane, orientation, and strategy identity. Classification is their only public
construction path. Typed digest copies copy the typed value directly; they
never re-ingest `digest.bytes()`.

Each in-place trial applicator requires a trial exactly equal to the supplied
authority, performs no clone or authority swap, reconstructs and structurally
validates the corresponding segment/full-circle/arc trace or plunge disk, and
returns a nonconstructible `AuditDepletionWitness2`. The witness binds kind,
pre/post canonical stock-state digests, motion digest, policy digest, exact
trace/disk identity, strategy version, and its own canonical digest. Replay
validates the complete witness value before deriving lineage; a bare depletion
digest is never accepted as evidence.

- [x] **Step 1: Write RED identity tests**

Add assertions that `EngagementAuditInput.build(...)` requires one exact
`DepletionPolicy`; `canonical_task1_bytes(EngagementCap)` binds both authored
angle and chord surrogate; all six opaque native motions expose a 32-byte digest;
authenticated carrier v2 binds that motion digest; and input v2 binds the
complete cap, policy, arc-surrogate version, audit decision-contract version,
native depletion version, ordered authenticated-operation digests, ordered
native-motion digests, and the independently recomputed native request digest.

```python
def test_policy_and_authored_cap_are_request_identity() -> None:
    baseline = _audit_input()

    assert _with_cap_theta(baseline, 0.7).digest != baseline.digest
    assert _with_chord_bound(baseline, 0.03125).digest != baseline.digest
    assert _with_center_limit(baseline, 2048).digest != baseline.digest


def test_native_request_is_recomputed_and_nonconstructible() -> None:
    audit_input = _audit_input()
    policy = _stock_2.build_audit_policy(
        audit_input.tool_radius.value,
        audit_input.engagement_cap.theta,
        audit_input.engagement_cap.chord_ratio,
        audit_input.depletion_policy.chord_bound.value,
        audit_input.depletion_policy.center_count_limit,
    )
    request = _stock_2.build_audit_native_request_identity(
        _boundary_rows(audit_input),
        _hole_rows(audit_input),
        policy,
        tuple(operation.motion for operation in audit_input.operations),
    )

    assert request.digest == audit_input.native_request_digest
    with pytest.raises(TypeError):
        _stock_2.AuditNativeRequestIdentity2()
```

Mutate one boundary vertex, hole, policy field, native motion, and motion order.
The native request digest must change. Equivalent ring rotation/orientation and
hole input order must canonicalize to the same request. Changing source metadata
without changing native geometry changes authenticated/input identity but not
native request identity. Compile-time tests prove `AuditDigest2<Domain>` has no
public raw-byte constructor or generic `from_bytes` retagger.

- [x] **Step 2: Run the contract RED gate**

```bash
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- tests/adaptive/test_canonical.py tests/engagement_audit/test_native_classification.py tests/engagement_audit/test_records.py tests/engagement_audit/test_input.py tests/engagement_audit/test_native_request_identity.py -n auto -q
```

Expected: missing depletion policy/input v2 and native request symbols. Do not
commit this RED state.

- [x] **Step 3: Implement cap/policy/input v2 identity**

Change `engagement-cap-v1` to `engagement-cap-v2` with `theta-radian` and
`chord-ratio` components. Add `depletion_policy` to the input factory and bump
`engagement-audit-input-v1` to v2. Bind `canonical_task1_bytes(...)` for cap and
policy plus all three native strategy-version byte strings and the ordered
authenticated-operation digests, including Task 3's exact arc-motion digest.
Give the other five opaque motion types canonical native bytes/digests, then
bump every authenticated carrier encoding to bind its native-motion digest.

The native request digest is SHA-256 over versioned canonical bytes containing
the native-canonical stock rings, `AuditPolicy2.digest()`, and the ordered
native-motion digests. Input v2 binds that digest plus the ordered authenticated
carrier digests. `build_audit_native_request_identity(...)` recomputes the
native request digest from actual rings, policy, and the closed tuple of opaque
motions. Task 4's later `AuditReplay2.build(...)` calls the same authority and
every motion call checks both expected digests at the current cursor. Raw
constructor and mutated-state tests must continue to fail closed.

Remove generic `AuditDigest2<Domain>::from_bytes`. Domain-specific ingress is
permitted only for external input and authenticated-operation expectations;
generated request,
motion, policy, decision, depletion, lineage, and result digests arise only from
their named canonical-hash authorities. The first Task 4 commit ends here and
runs the focused identity, native, Ruff, strict-mypy, docs, and diff gates before
review.

- [ ] **Step 4: Restore the P0 falsifiers and write native verdict RED tests**

Port only the ledger-targeted exact seam, bounded work, false-certificate,
swept-annulus, rotation, shared-root, full-turn, and adaptive-arc controls.
Assert three outcomes independently. `CAP_EXCEEDED` requires a replay-validated
exact station witness evaluated against the unguarded policy chord surrogate.
Guarded-cap failure and an inconclusive swept upper bound cause refinement or
`UNRESOLVED`; they never prove exceedance. Complete exact domain coverage and
exact between-station closure produce `CERTIFIED`. Exhausted or unsupported
proof closure without a witness produces `UNRESOLVED`. A live falsifier proves
witness liveness independently; a dead witness may not produce
`CAP_EXCEEDED`. Exact witness identity binds stock, opaque motion, policy,
canonical parameter, and station disposition, and is validated before decision
digest construction.

- [ ] **Step 5: Factor and bind one native decision adapter**

Link the existing event-exact segment/full-circle oracle core into `_stock_2`
through one shared CMake target; do not copy its algorithm. Treat `73d5372` as
theorem, falsifier, and refinement-strategy input only: do not transfer its
Boolean result, double trigonometric centers, relative spacing floor, or double
angular comparisons into authority. Arc certification evaluates Task 3's exact
clipped chart intervals through the shared evaluator, uses non-cyclic event
order, queries the unguarded exact cap at every examined station, and certifies
only when an exact between-station predicate closes every interval. Until that
general closure exists, the adapter may certify a partial arc only from a sound
safe superset proof such as a complete full-circle authority; otherwise it
returns `UNRESOLVED` unless an exact violating station exists. The spatial
floor is an exact rational millimetre value or exact squared-distance threshold;
depth and node limits are exact integers. All are bound into decision-strategy
canonical bytes. Reaching any limit yields `UNRESOLVED` unless an already
validated exact witness yields `CAP_EXCEEDED`; none may decide `CERTIFIED`.
Return one witness-bearing internal result shape for all
motions. The exact verdict and decision digest are computed before any reporting
maximum; reporting probes cannot feed back into verdict, digest, or stock.

- [ ] **Step 5A: Write replay-transaction RED tests**

Before implementing replay, require early finalization, duplicate calls,
omitted operations, wrong authenticated-operation identity, and wrong actual
motion identity to fail with named errors. Add controlled native failure
injection at decision, trial depletion, witness validation, lineage/result
construction, and finalization. Every failure must preserve authority stock,
cursor, lineage, and finalized state. Native instrumentation proves exactly one
`Stock2::clone()` per mutating call and zero clones for retract/clearance calls;
the trial applicators themselves must report zero internal clones and swaps.
Finalization succeeds exactly once after complete consumption.

- [ ] **Step 6: Implement the opaque replay owner and atomic calls**

For each mutating native call execute exactly:

```text
check replay state, cursor, operation digest, and actual motion digest
decide against immutable current stock
clone stock exactly once
apply validated depletion in-place to the trial for every lateral verdict
validate decision and depletion evidence
derive operation-bound post-lineage and complete canonical result
construct the complete next replay state
swap next state into authority with no-throw moves
return the already-built result
```

No allocation, hashing, validation, or cursor update remains after the swap.
Any exception precedes swap and yields no result. Plunge uses exact disk
depletion; retract and clearance transport return a non-engaging result with
identical pre/post lineage. Separate nanobind functions consume each opaque
motion type; no coordinates, phases, radii, or generic motion union cross back
through Python. `finish_audit_replay(...)` succeeds exactly once and only after
the cursor consumed every bound operation; it returns the opaque completion
value used by Task 5. Calls after finalization fail with a named replay-state
error.

The initial `StockLineageDigest2` is
`SHA256(stock-lineage-seed-v1 {input_digest,
verified_native_request_digest})`. It is derived from both typed roots, never
byte-equal retagging of `AuditInputDigest2`. Replay state is decomposed into an
immutable request and a mutable stock/cursor/lineage/finalized state. Existing
`Stock2::subtract_exact_*` paths remain intact; Task 4 adds validated in-place
trial applicators in `stock_exact_depletion_2.*` and validates them before any
later removal or binding extraction is considered.

Every mutating post-lineage digest is
`SHA256(stock-lineage-transition-v1 {verified_request_digest, cursor,
pre_lineage_digest, authenticated_operation_digest, actual_motion_digest,
decision_witness_digest_if_lateral, depletion_witness_digest})`. Non-engaging
operations preserve the exact lineage bytes. `stock_2.cpp` changes are limited
to registration and delegation. One externally linked exact disk-union region
builder, `build_exact_disk_union_region_2(...)`, is declared only in
`exact_disk_region_2.h` and defined only in `exact_disk_region_2.cpp`. Every new
trial applicator and every legacy exact segment/full-circle/arc subtraction path
calls that one symbol. The existing anonymous helper first becomes a delegating
compatibility wrapper and remains present until the new route is independently
validated and the user authorizes later removal. Both new files are compiled
into the shared native core and the `_stock_2` module through CMake.

- [ ] **Step 7: Run native and Python GREEN gates**

```bash
pixi run audit-native
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- tests/engagement_audit/test_native_replay.py tests/engagement_audit/test_false_arc_certificate.py tests/test_false_certificate.py tests/test_growth_bound.py tests/adaptive/test_canonical.py tests/engagement_audit/test_input.py -n auto --testmon -q
pixi run lint
pixi run types-audit
pixi run types-adaptive
pixi run -e docs docs
git diff --check
```

`types-adaptive` is compared against the known six-error P0 baseline and is not
reported green unless those errors are separately repaired. The task is not
committable until every motion type passes and no authoritative symbol reaches
legacy `subtract_arc_sweep`.

- [ ] **Step 8: Commit the three closed native slices**

```bash
# Slice A stages only identity/policy/motion/request files and tests.
git commit -m "feat(audit): bind native request"

# Slice B stages only decision/oracle/falsifier/witness files and tests.
git commit -m "feat(audit): decide exact motion"

# Slice C stages only trial-depletion/lineage/result/replay files and tests.
git commit -m "feat(audit): transact exact motion"
```

Before each commit, obtain an independent scoped review and run that slice's
focused native/Python tests, Ruff, strict mypy, strict docs, and diff check. Use
`git diff --name-only` against the prior slice to stage every new or modified
file explicitly; an incomplete fixed path list is a plan failure.

### Task 5: Implement public replay and truthful report semantics

**Files:**

- Create: `src/compas_cgal/engagement_audit/native.py`
- Create: `src/compas_cgal/engagement_audit/replay.py`
- Create: `src/compas_cgal/engagement_audit/report.py`
- Modify: `src/compas_cgal/engagement_audit/records.py`
- Modify: `src/compas_cgal/engagement_audit/errors.py`
- Create: `tests/engagement_audit/test_replay.py`
- Create: `tests/engagement_audit/test_report.py`
- Modify: `tests/adaptive/typecheck/auditor_contract.py`
- Modify: `docs/engagement_audit.md`

**Interfaces:**

- Consumes: `EngagementAuditInput`, preclassified opaque native operations,
  and Task 4's native replay/result values. Python never accepts or exposes
  authoritative `Stock2`.
- Produces: `audit_toolpath_engagement(...) -> EngagementAuditReport`,
  `EngagementAuditReport.require_certified() -> None`.

- [ ] **Step 1: Write RED replay-order tests**

Require the public API to accept exactly one input and dispatch each operation
once, in stream order:

```python
def test_public_replay_accepts_only_authenticated_input() -> None:
    signature = inspect.signature(audit_toolpath_engagement)

    assert tuple(signature.parameters) == ("audit_input",)


def test_returned_unresolved_motion_still_advances_lineage() -> None:
    report = audit_toolpath_engagement(_unresolved_then_segment_input())

    assert report.operations[0].post_motion_stock_lineage == report.operations[1].pre_motion_stock_lineage
```

Add exception injection proving no report escapes, first pre-lineage equals the
native read-back of the domain-separated input/request lineage seed, terminal
report lineage equals the last result, and every source operation has exactly
one result at the same index.

- [ ] **Step 2: Repair the closed record union**

Rename `motion_certificate_digest` to `native_decision_digest`; add depletion
witness and post-lineage to `MeasuredOperationAudit`; add
`PlungeOperationAudit`; and remove `vertical_plunge` from
`NonEngagingReason`. Output records rename `operation_digest` to
`authenticated_operation_digest`; authenticated input carriers retain their
nested source digest. Bump each changed canonical version. Every factory checks
exact native result type, digest sizes, operation identity, lineage adjacency,
and field presence appropriate to its chronology.

```python
OperationAudit: TypeAlias = (
    MeasuredOperationAudit
    | PlungeOperationAudit
    | NonEngagingOperationAudit
)
```

- [ ] **Step 3: Write RED report truth tests**

```python
def test_unresolved_is_not_a_cap_violation_or_certification() -> None:
    report = _report_with_verdict("unresolved")

    assert report.cap_exceeded_count == 0
    assert report.unresolved_count == 1
    with pytest.raises(UnresolvedEngagementAuditError, match="1 unresolved"):
        report.require_certified()


def test_report_requires_a_measured_lateral_motion() -> None:
    with pytest.raises(NoMeasuredLateralMotionError):
        EngagementAuditReport.build(_only_retract_records())
```

- [ ] **Step 4: Implement the native adapter and replay**

`native.py` starts one native replay owner from input rings/digest and passes
each opaque motion directly to its one typed Task 4 call. It translates only
the returned opaque result into a Python audit record and passes
`operation.digest` as the authenticated operation identity. `replay.py`
consumes the immutable authenticated stream, never rereads COMPAS operations,
and emits a
record tuple only after the entire native stream returns and native
`finish_audit_replay(...)` proves exact cardinality and terminal lineage.

- [ ] **Step 5: Implement derived report construction**

`EngagementAuditReport.build(input_identity, operations)` derives counts,
including the distinct plunge count, reporting maximum TEA, ordered digest,
terminal lineage, and canonical bytes. It first proves exact index coverage,
lineage continuity, and that the five disjoint counts sum to the operation
count. It accepts no caller aggregate. `require_certified()` distinguishes
proved exceedance from unresolved and requires at least one measured lateral
motion; plunges do not satisfy that requirement.

- [ ] **Step 6: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-audit
pixi run pytest -- tests/engagement_audit/test_replay.py tests/engagement_audit/test_report.py -n auto --testmon -q
pixi run -e docs docs
git diff --check
git add src/compas_cgal/engagement_audit tests/engagement_audit tests/adaptive/typecheck/auditor_contract.py docs/engagement_audit.md
git commit -m "feat(audit): replay truthful verdicts"
```

### Task 6: Migrate consumers and record Stage 1 evidence

**Files:**

- Modify: `benchmarks/runner.py`
- Modify: `benchmarks/models.py`
- Modify: `tests/benchmarks/test_runner.py`
- Modify: `tests/benchmarks/test_models.py`
- Create: `docs/auditor_convergence.md`
- Modify: `mkdocs.yml`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: authoritative `EngagementAuditReport`.
- Produces: benchmark records with
  certified/exceeded/unresolved/plunge/non-engaging counts and Stage 1 measured
  evidence.

- [ ] **Step 1: Write RED consumer-contract tests**

Require benchmark construction to reject the legacy `EngagementReport`, carry
all five operation counts, and round-trip JSON/schema without inferring
unresolved from `stations` or TEA.

- [ ] **Step 2: Migrate benchmark runner and schema**

Call only `engagement_audit.audit_toolpath_engagement`. Remove inference from
`stations > 0 and not cap_certified`. Keep the legacy module import-free from
all acceptance consumers.

- [ ] **Step 3: Measure the 12x8 arc acceptance criterion**

Run one committed Pixi task against the regulated 12x8 fixture and record input
digest, build identity, circle count, certified/exceeded/unresolved counts,
station/event counts, and wall time. The result may fail the usefulness gate;
it must not be paraphrased as green.

- [ ] **Step 4: Update durable documentation**

Document authoritative vs legacy call graphs, truth-domain semantics, native
arc evidence, exact limitations, and the measured 12x8 result. Add the page to
MkDocs navigation.

- [ ] **Step 5: Run Stage 1 gates**

```bash
pixi run pytest -- tests/engagement_audit tests/benchmarks tests/test_false_certificate.py tests/test_growth_bound.py -n auto -q
pixi run lint
pixi run types-adaptive
pixi run -e docs docs
git diff --check
```

- [ ] **Step 6: Commit and update progress**

```bash
git add benchmarks tests/benchmarks docs/auditor_convergence.md mkdocs.yml .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "feat(bench): consume truthful audit"
```

Record the Stage 1 acceptance commit and set the next exact command to P2 Task
1's replay-certificate RED test.
