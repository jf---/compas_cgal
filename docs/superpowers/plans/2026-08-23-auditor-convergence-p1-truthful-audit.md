# Auditor Convergence P1 Truthful Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one authoritative engagement audit in which every lateral cut has
a native three-way verdict and missing measurement can never become compliance.

**Architecture:** Introduce a focused `engagement_audit` package alongside the
legacy diagnostic. A content-addressed input owns typed geometry and build
identity; geometry classification admits only proved non-engaging motions or
native-certified lateral motions; report construction derives all aggregates
from a closed operation-audit union. Port the adaptive native arc certifier only
after P0 proves its dependency closure.

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
plunge; Task 3 passes that value directly to native depletion without exposing
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

### Task 3: Implement measure-before-deplete replay and report semantics

**Files:**

- Create: `src/compas_cgal/engagement_audit/native.py`
- Create: `src/compas_cgal/engagement_audit/replay.py`
- Create: `src/compas_cgal/engagement_audit/report.py`
- Create: `tests/engagement_audit/test_replay.py`
- Create: `tests/engagement_audit/test_report.py`
- Modify: `src/compas_cgal/engagement_audit/errors.py`

**Interfaces:**

- Consumes: `EngagementAuditInput`, `Stock`, preclassified opaque native
  operations, native segment/arc certifiers.
- Produces: `audit_toolpath_engagement(...) -> EngagementAuditReport`,
  `EngagementAuditReport.require_certified() -> None`.

- [ ] **Step 1: Write RED replay-order tests**

Instrument native certification and depletion to prove chronology:

```python
def test_lateral_motion_is_certified_before_its_sweep_is_removed(monkeypatch: pytest.MonkeyPatch) -> None:
    order: list[str] = []
    _track_certify(monkeypatch, order)
    _track_deplete(monkeypatch, order)

    audit_toolpath_engagement(_single_cut_input())

    assert order == ["certify:0", "deplete:0"]
```

Add a fault-injection test proving certification failure leaves authoritative
stock and output absent; no partial report escapes.

- [ ] **Step 2: Write RED report truth tests**

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

- [ ] **Step 3: Implement the native adapter and replay**

`native.py` passes each opaque native motion directly to one native certifier
and returns `MeasuredOperationAudit`. It never reconstructs geometry or
contains geometry policy. `replay.py` consumes only the preclassified immutable
stream from `EngagementAuditInput`; it certifies, then depletes each operation
in order and returns the immutable record tuple. It never rereads or
reclassifies `ToolpathOperation` objects.

- [ ] **Step 4: Implement derived report construction**

`EngagementAuditReport.build(input_identity, operations)` derives counts,
maximum reported TEA, ordered digest, and canonical bytes. It accepts no caller
aggregate. `require_certified()` distinguishes proved exceedance from
unresolved and requires at least one measured motion.

- [ ] **Step 5: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/engagement_audit/test_replay.py tests/engagement_audit/test_report.py -n auto --testmon -q
git diff --check
git add src/compas_cgal/engagement_audit tests/engagement_audit tests/adaptive/typecheck/consumer_contract.py
git commit -m "feat(audit): replay truthful verdicts"
```

### Task 4: Port and validate adaptive native arc certification

**Files:**

- Modify: `src/engagement_2.h`
- Modify: `src/engagement_2.cpp`
- Modify: `src/stock_2.cpp`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/engagement_audit/native.py`
- Create: `tests/native/test_engagement_arc.cpp`
- Create: `tests/engagement_audit/test_false_arc_certificate.py`
- Modify: `tests/engagement_audit/test_replay.py`

**Interfaces:**

- Consumes: P0's exact dependency closure and existing segment swept-annulus
  primitives.
- Produces: `_stock_2.certify_arc_tea(AuditArcMotion2, ...) ->
  tuple[NativeMotionVerdict, float, int, bytes]` and the sole authoritative arc
  audit path.

- [ ] **Step 1: Write RED native and Python boundary tests**

Require annular, machined, and spiral ribs to refuse or report cap exceeded;
require a non-vacuous contacting clear arc to certify. Include full circle,
partial CW/CCW arc, rotation, translation, and scale variants. The exact
`cap_exceeded` station predicate independently proves witness liveness.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/engagement_audit/test_false_arc_certificate.py -n auto -q
```

Expected: missing authoritative native arc result shape or false certification.

- [ ] **Step 3: Reapply the P0-approved native dependency closure**

Port only ledger rows targeted to P1. Preserve exact predicates, release-build
shared-root exceptions, corrected full-turn precondition, adaptive spatial
floor, finite depth, swept-annulus guard, and three-way verdict. Do not copy
the source branch's Python diagnostic mirror.

- [ ] **Step 4: Bind one exact API**

Expose one `certify_arc_tea` nanobind function that consumes the opaque
`AuditArcMotion2` directly, plus its matching stub. Separate
nanobind lambdas are required for any overloads. Translate each native failure
to its named Python exception; no boolean compatibility return is added.

- [ ] **Step 5: Run native/Python GREEN and mutation controls**

```bash
pixi run pytest -- tests/engagement_audit/test_false_arc_certificate.py tests/engagement_audit/test_replay.py tests/test_false_certificate.py tests/test_growth_bound.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
git diff --check
```

Expected: all negative and positive controls pass; the old segment certificate
bytes remain unchanged unless P0 explicitly classified a version bump required.

- [ ] **Step 6: Commit**

```bash
git add src/engagement_2.h src/engagement_2.cpp src/stock_2.cpp src/compas_cgal/_stock_2.pyi src/compas_cgal/engagement_audit/native.py tests/native/test_engagement_arc.cpp tests/engagement_audit
git commit -m "feat(audit): certify arcs adaptively"
```

### Task 5: Migrate consumers and record Stage 1 evidence

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
- Produces: benchmark records with certified/exceeded/unresolved/non-engaging
  counts and Stage 1 measured evidence.

- [ ] **Step 1: Write RED consumer-contract tests**

Require benchmark construction to reject the legacy `EngagementReport`, carry
all four operation counts, and round-trip JSON/schema without inferring
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
