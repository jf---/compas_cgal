# Exact One-Root Full-Circle Source Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Admit shared-radical one-root trim vertices into the exact
full-circle event source without assigning conjugate norm roots physical
topology.

**Architecture:** One canonical coordinate grammar feeds one quarter-chart
vertex-passage producer. The integer norm discovers partition roots; the
original radical equation selects physical roots, and the same radical
predicate supplies exact active-sheet signs at rational cell witnesses.

**Tech Stack:** C++20, CGAL 6.0.1 exact algebraic kernel,
`Gps_circle_segment_traits_2`, CORE exact rationals, nanobind, Python 3.12,
pytest-xdist, pytest-testmon, strict MkDocs.

## Global constraints

- No binary64 collapse, epsilon, sampling, conjugate event, alternate oracle,
  fallback source, or v3 reconstruction path.
- Source and certificate fields use one exact grammar for rational and
  shared-radical coordinates.
- Distinct nonzero coordinate radicands fail loudly with
  `UnsupportedAlgebraicVertexProjectionError`.
- Tests precede production changes and include geometric/proof context
  docstrings.
- Every pytest command uses `-n auto`; affected tests use `--testmon`.
- Update developer MkDocs and the Held comparison before the stage is
  complete.

---

### Task 1: RED one-root source and physical-root contracts

**Files**

- Modify: `tests/adaptive/test_circle_oracle.py`
- Modify: `tests/adaptive/test_motion_certificate.py`

**Interfaces**

- Consumes: existing
  `_continuous_tea_2.audit_full_circle_tea_event_exact(...)`.
- Produces: an exact expected verdict for the diagonal post-capsule case and
  observable one-root projection evidence.

- [x] **Step 1: Reclassify the real positive control**

Change the post-capsule test from expecting
`UnresolvedMotionEventError` to the exact proved outcome. Retain assertions
that stock boundary identity, lineage, and point containment are unchanged.

- [x] **Step 2: Add the native partition contract**

Use two overlapping boundary disks with one physical shared-radical
intersection and one clipped conjugate, then require:

```python
assert verdict == "certified"
assert trace.partition.source_kind == "full-circle-boundary-pullbacks-v4"
assert any(projection.one_root_predicate for projection in trace.partition.projections)
assert one_root_physical_root_ids(trace.partition)
assert one_root_conjugate_root_ids(trace.partition)
assert verify_event_partition(trace.partition).verdict.name == "CERTIFIED"
```

The docstring explains why the clipped conjugate remains a partition root but
may not own an endpoint event. A generic eventless fibre is not sufficient
evidence because it may belong to an unrelated projection family.

- [x] **Step 3: Add mutation contracts**

Require `delete-one-root-predicate`, `alter-one-root-radicand`,
`alter-one-root-rational-part`, and `move-one-root-event-to-conjugate` to
produce `UNRESOLVED_DEGENERACY`.

- [x] **Step 4: Run RED**

```bash
pixi run pytest tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_motion_certificate.py \
  -n auto --testmon -q
```

Expected: failure because the current source remains v3 and the motion remains
unresolved.

---

### Task 2: Canonical coordinate and radical-predicate authority

**Files**

- Create: `src/continuous_tea_2/circle_vertex_source.h`
- Create: `src/continuous_tea_2/circle_vertex_source.cpp`
- Modify: `src/continuous_tea_2/partition_certificate.h`
- Modify: `src/continuous_tea_2/segment_partition.h`
- Modify: `src/exact_algebraic_1.cpp`
- Modify: `src/continuous_tea_2/bind_partition_records.cpp`
- Modify: `src/continuous_tea_2/event_certificate.cpp`
- Modify: `CMakeLists.txt`

**Interfaces**

- `encode_full_circle_coordinate(const GpsPoint::CoordNT&) -> std::string`
- `decode_full_circle_coordinate(const std::string&) -> OneRootCoordinate2`
- `OneRootPredicate2::build(A, B, alpha) -> OneRootPredicate2`
- Optional `ProjectionInput2::one_root_predicate`
- Optional `ProjectionRecord2::one_root_predicate`

- [x] **Step 1: Implement one coordinate grammar**

Encode `full-circle-one-root-coordinate-v1` with exact rational base,
coefficient, and radicand. Canonicalize rational coordinates to `(base, 0, 0)`
and reject negative or malformed radicands.

- [x] **Step 2: Implement validated radical predicates**

`OneRootPredicate2::build(...)` requires nonempty, jointly normalized integer
coefficient arrays and a positive exact radicand. Normalize only their shared
content and shared sign: independently making both arrays primitive changes
the represented radical equation. Add no public raw construction path in
Python.

- [x] **Step 3: Validate norm ownership in the algebraic partition**

For a one-root input, reconstruct

```text
den(alpha) A^2 - num(alpha) B^2
```

and require canonical equality with the submitted projection polynomial.
Reject simultaneous integer and one-root signed predicates.

- [x] **Step 4: Filter source events per algebraic root**

Evaluate `A` and `B` with `sign_at_1_object()`. Attach input events only when
`A + B sqrt(alpha) == 0`; retain every norm root in the partition regardless.

- [x] **Step 5: Bind and expose the authority**

Canonicalize the radical predicate in `projection-record-v3`, expose a
read-only nanobind view, and add the four named mutation operators.

- [x] **Step 6: Run the focused native/Python build**

```bash
pixi run pytest tests/adaptive/test_circle_oracle.py \
  -n auto --testmon -q
```

Expected: still RED at the source version because the circle producer does not
yet emit the new predicate.

---

### Task 3: Full-circle passage production and active-sheet sign

**Files**

- Modify: `src/continuous_tea_2/circle_oracle.cpp`
- Modify: `src/continuous_tea_2/circle_projection.h`
- Modify: `src/continuous_tea_2/circle_projection.cpp`
- Modify: `src/continuous_tea_2/circle_fibre_state.cpp`
- Modify: `src/continuous_tea_2/circle_source.cpp`
- Modify: `src/continuous_tea_2/circle_pair_projection.h`
- Modify: `src/continuous_tea_2/circle_pair_projection.cpp`
- Modify: `src/continuous_tea_2/circle_strata.h`
- Modify: `src/continuous_tea_2/circle_strata.cpp`
- Modify: `src/continuous_tea_2/event_certificate.cpp`

**Interfaces**

- `full_circle_vertex_passage(...) -> ProjectionInput2`
- `signed_projection_sign(...)` accepts exactly one integer or one-root
  predicate form.
- Source kind `full-circle-boundary-pullbacks-v4`.

- [x] **Step 1: Emit exact coordinate sources**

Replace rational-coordinate admission in `circle_oracle.cpp` with the canonical
one-root encoding for every source and target vertex. Remove
`rational_boundary_vertices`; do not retain a v3 path.

- [x] **Step 2: Derive the quarter-chart radical equation**

Construct jointly normalized integer `A` and `B`, derive the integer norm, and
submit all three as one `ProjectionInput2`. Require a shared radicand for
nonzero coordinate radical terms.

- [x] **Step 3: Evaluate active sheets exactly**

At each rational cell witness, evaluate homogeneous `A` and `B`, then determine
`sign(A + B sqrt(alpha))` through exact sign and squared-magnitude comparison.
Never use the norm sign for material-side state.

- [x] **Step 4: Mark physical endpoint evidence**

Every retained endpoint event must set
`original_equations_rechecked = true`,
`orientation_rechecked = true`, and `trim_disposition = "accepted"`.

- [x] **Step 5: Add and reconstruct the source**

Emit/reconstruct v4 for algebraic coordinates while retaining the
independently replayed rational v3 path. Update oracle/strategy versions to
`event-exact-motion-oracle-v5` and
`full-circle-cell-strata-exact-v3`.

- [x] **Step 6: Replace the quadratic pair product**

Build topology without pair factors, derive canonical pair requests from the
exact active branches at each topology-cell witness, and rebuild with only
those requested center/rim branch pairs. Union requests exposed by each
refinement and repeat until the finite request set reaches a fixed point. Bind
the closed request set into the v4 source. Final cell authority must remain
unresolved when an active material run lacks either required pair factor.

Before circle/circle elimination, exactly divide each circle pullback by its
proved `1 + u²` chart-denominator factor. It has no real roots but otherwise
introduces a complex common component that makes the resultant identically
zero.

- [x] **Step 7: Run GREEN**

```bash
pixi run pytest tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_motion_certificate.py \
  -n auto --testmon -q
```

Result: `30 passed in 40.25s`; the post-capsule case has a proved exact
verdict, and the rational v3 144-factor positive control remains unchanged.

---

### Task 4: Task 13F integration and adversarial replay

**Files**

- Modify: `tests/adaptive/test_generator.py`
- Modify: `tests/adaptive/test_acceptance.py`
- Modify: `tests/adaptive/test_replay.py`
- Modify: `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`

**Interfaces**

- Consumes: v5 exact full-circle oracle.
- Produces: Task 13F continuation beyond the previously unresolved fourth
  trial.

- [ ] **Step 1: Add the real continuation RED**

Require the fixed radius-1 family to certify the fourth link and following
circle without translating native authority to `UnresolvedMotionEventError`.
Keep the invariant no-lower-ranked-dispatch assertion.

- [ ] **Step 2: Add source/replay cross-wiring mutations**

Reject a one-root source borrowed from another stock lineage, a stale v4
partition under a new physical parent, and a physical event reassigned to an
eventless conjugate fibre.

- [ ] **Step 3: Run affected GREEN**

```bash
pixi run pytest tests/adaptive/test_generator.py \
  tests/adaptive/test_acceptance.py \
  tests/adaptive/test_replay.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_circle_oracle.py \
  -n auto --testmon -q
```

- [ ] **Step 4: Commit the coherent theorem**

```bash
git add CMakeLists.txt src/continuous_tea_2.cpp \
  src/continuous_tea_2 tests/adaptive
git commit -m "feat(tea): admit one-root circle sources"
```

---

### Task 5: Documentation, full gates, performance, and publication

**Files**

- Modify: `docs/continuous_engagement.md`
- Modify: `docs/segment_site_mat.md`
- Modify:
  `docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`
- Modify: `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`
- Modify: this plan

- [x] **Step 1: Capture the theorem**

Document source grammar, norm-versus-physical-root separation, signed
active-sheet proof, exact limitation, mutation evidence, Task 13F effect, and
the changed Held comparison.

- [ ] **Step 2: Measure bounded workloads**

Repeat five warm Release audits for the rational concentric,
rational cross-support, and one-root post-capsule cases. Report medians and
structural projection/root/cell counts without claiming Held parity.

- [x] **Step 3: Run all gates**

```bash
pixi run pytest tests/adaptive -n auto -q
pixi run lint
pixi run types-adaptive
pixi run docs
git diff --check
```

Result: `522 passed in 242.54s`; Ruff, strict mypy over `27` source files,
strict MkDocs, affected `--testmon`, and diff checks passed.

- [x] **Step 4: Commit documentation**

```bash
git add docs
git commit -m "docs(tea): record one-root circle theorem"
```

- [ ] **Step 5: Publish safely**

Cancel in-progress branch CI runs, push the branch, and require
`git ls-remote origin refs/heads/codex/exact-certified-adaptive-phase1-t9` to
equal local `HEAD`.
