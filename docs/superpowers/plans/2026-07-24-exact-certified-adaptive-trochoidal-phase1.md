# Exact-certified adaptive trochoidal-MAT — Phase 1 implementation plan

> Execute with `superpowers:executing-plans` and
> `superpowers:subagent-driven-development`. Run one task at a time through
> implementer review, specification review, and code-quality review. Do not
> dispatch two agents to files that share native state.

**Goal:** Return a non-vacuous polygonal-pocket toolpath whose lateral
engagement, gouge freedom, exact tool-reachable-material coverage, and artifact
lineage are machine-checkable.

**Architecture:** A typed Python orchestrator proposes one-sided MATHSM circles
from a true CGAL segment-site Voronoi medial axis. Exact native components own
the proof boundaries: Epeck stock/station predicates, exact-on-guide
conservative depletion, exact construction of `C_r = D ⊖ B_r` and
`M_r = C_r ⊕ B_r`, an exact-with-sqrt full-sweep coverage ledger, and an
event-exact continuous engagement oracle. Candidate evaluation is
transactional. Existing straight-skeleton generators, rounded depletion APIs,
and audit compatibility paths remain unchanged.

**Governing spec:**
`docs/superpowers/specs/2026-07-24-exact-certified-adaptive-trochoidal-design.md`
at exactly `6334e8d6446086e3dfbe688c4a3656dd6845c3b1`.

**Start point:** branch `codex/exact-certified-adaptive-phase1`, based on
`1860167929e50f38fdae8d67ef77e9c967364f1c`.

## Non-negotiable execution rules

- Work only in the isolated task worktree.
- Use Pixi for every build, package, test, lint, type, schema, and mutation
  command.
- Every pytest invocation includes `-n auto`.
- After code changes, the targeted GREEN command includes `--testmon`.
- Write the failing test before implementation. Do not skip, `skipif`, or
  xfail it.
- Run Ruff before every Python commit and `mypy --strict` for the adaptive
  package.
- Add new paths alongside working paths. Do not redirect or remove legacy
  behavior in this phase.
- No reported angle, tolerance, or sampled maximum may decide a certificate.
- A native algebraic gate that fails stops the dependent tasks. It does not
  authorize a fallback.
- Commit only the task’s files. Author and committer are
  `Jelle Feringa <jelleferinga@gmail.com>`.

## Proof and implementation dependency graph

```text
T0 Pixi + native-source lock
 └─ T1 typed motion/policy/cap boundary
     └─ T1A canonical identity foundation
         ├─ T2 exact conservative depletion
         ├─ T3 exact reachable domain + full-sweep coverage
         └─ T4 event corpus + boundary extraction
             └─ T5 algebraic event substrate
                 ├─ T6 exact segment oracle
                 └─ T7 exact full-circle oracle
                     └─ T8 typed motion certifier
T0 + T3 + T5 ─ T9 exactly clipped true MAT
                └─ T10 typed MAT + neck/candidate lattice
T2 + T3 + T8 + T10 ─ T11 entry + containment + InputIdentity
                     ├─ T11A fresh independent replay
                     └─ T12 transactional candidate evaluator
T11A + T12 ─ T13 traversal + generator + coverage
                             └─ T14 artifact assembly + schema
                                 └─ T15 acceptance + mutation campaign
                                     └─ T16 full verification/review
```

## File morphology

```text
src/
  exact_motion_2.h
  exact_algebraic_1.h
  exact_algebraic_1.cpp
  reachable_domain_2.h
  reachable_domain_2.cpp
  exact_depletion_2.h
  exact_depletion_2.cpp
  coverage_2.h
  coverage_2.cpp
  continuous_tea_2/
    result.h
    boundary_events.h
    boundary_events.cpp
    parameter_charts.h
    parameter_charts.cpp
    event_partition.h
    event_partition.cpp
    segment_oracle.cpp
    circle_oracle.cpp
  continuous_tea_2.cpp              # nanobind adapter only
  segment_site_mat.h
  segment_site_mat.cpp
  segment_site_mat_sampling.h
  segment_site_mat_sampling.cpp
  medial_axis_2.cpp                 # nanobind adapter only
  containment_2.h
  containment_2.cpp
  containment_2_bindings.cpp

src/compas_cgal/adaptive/
  __init__.py                       # minimal; no __all__
  errors.py
  units.py
  motion.py
  policy.py
  canonical.py
  identity.py
  operation.py
  stock_area.py
  reachable_domain.py
  coverage.py
  motion_certificate.py
  replay.py
  medial_axis.py
  neck.py
  entry.py
  containment.py
  candidates.py
  transaction.py
  traversal.py
  certificate.py
  schema.py
  generator.py

src/compas_cgal/
  _stock_2.pyi
  _coverage_2.pyi
  _continuous_tea_2.pyi
  _medial_axis_2.pyi
  _containment_2.pyi

tests/adaptive/
  conftest.py
  test_units.py
  test_motion.py
  test_policy.py
  test_canonical.py
  test_identity.py
  test_operation.py
  test_exact_depletion.py
  test_reachable_domain.py
  test_coverage.py
  test_motion_counterexamples.py
  test_event_substrate.py
  test_segment_oracle.py
  test_circle_oracle.py
  test_motion_certificate.py
  test_replay.py
  test_medial_axis.py
  test_neck.py
  test_entry.py
  test_containment.py
  test_candidates.py
  test_transaction.py
  test_traversal.py
  test_certificate.py
  test_schema.py
  test_acceptance.py
  fixtures/
    held_fig5.json
    tractable_pocket.json
    event_corpus.json

schemas/
  certified_toolpath-v1.schema.json

cmake/
  NativeDependencies.cmake
  VerifyNativeSource.cmake

docs/build/
  native-dependency-lock.json

scripts/
  run-adaptive-mutations.py
```

No native implementation file may grow past 1,000 lines. Split earlier when a
file acquires a second responsibility.

---

## Task 0 — Lock Pixi and native sources

**Files**

- Modify: `pyproject.toml`
- Modify: `.gitignore`
- Modify: `CMakeLists.txt`
- Create: `pixi.lock`
- Create: `cmake/NativeDependencies.cmake`
- Create: `cmake/VerifyNativeSource.cmake`
- Create: `docs/build/native-dependency-lock.json`

### Step 1: add the Pixi workspace

Use Python 3.12 and conda-forge:

```toml
[tool.pixi.workspace]
channels = ["conda-forge"]
platforms = ["osx-arm64", "linux-64"]

[tool.pixi.dependencies]
python = "3.12.*"
cmake = ">=3.15"
ninja = "*"
compas = "*"
numpy = "*"
scipy = "*"

[tool.pixi.pypi-dependencies]
compas-cgal = { path = ".", editable = true }
nanobind = ">=1.3.2"
scikit-build-core = ">=0.10"
build = "*"
pytest = ">=7"
pytest-xdist = "*"
pytest-testmon = "*"
hypothesis = "*"
ruff = "*"
mypy = "*"
jsonschema = "*"

[tool.pixi.pypi-options]
no-build-isolation = ["compas-cgal"]

[tool.pixi.tasks]
_editable-rebuild = '''python -c "from compas_cgal import _stock_2, _toolpath"'''
pytest = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest''', depends-on = ["_editable-rebuild"] }
baseline = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests -n auto -q''', depends-on = ["_editable-rebuild"] }
affected = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests -n auto --testmon -q''', depends-on = ["_editable-rebuild"] }
lint = "ruff check src/compas_cgal tests"
format-adaptive = "ruff format src/compas_cgal/adaptive tests/adaptive"
types-adaptive = "mypy --strict src/compas_cgal/adaptive"
schema = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests/adaptive/test_schema.py -n auto -q''', depends-on = ["_editable-rebuild"] }
mutations-adaptive = "python scripts/run-adaptive-mutations.py"
wheel = "python -m build --wheel --no-isolation"
_native-lock-configure = "cmake -S . -B build/native-lock-check -G Ninja"
_native-lock-build = "cmake --build build/native-lock-check --target external_downloads"
native-lock-check = { depends-on = ["_native-lock-configure", "_native-lock-build"] }
```

Set `editable.rebuild = true` in the existing `[tool.scikit-build]` table; its
persistent `build-dir` already satisfies the rebuild contract. The coordinator
rebuilds before xdist, then workers receive that exact build directory through
`SKBUILD_EDITABLE_SKIP`; all plan-level `pixi run pytest ...` calls use the
public task. Do not add a pip wrapper task. Add `.pixi/`, `.testmondata*`, and
`.mypy_cache/` to `.gitignore`.
Keep `.env` and `.DS_Store` ignored.

Add `compas`, `numpy`, and `scipy` as hard project runtime dependencies.
`scipy` is mandatory because the existing nanobind surface includes
`nanobind/eigen/sparse.h`; no conditional import or optional-dependency path is
allowed.

### Step 2: lock every native source

Move the CGAL 6.0.1, Boost 1.82.0, and Eigen 3.4.0 URLs plus independently
verified SHA-256 archive hashes into `NativeDependencies.cmake`. Every
`ExternalProject_Add` uses `URL_HASH SHA256=...`.

`VerifyNativeSource.cmake` computes a canonical content digest over sorted
relative paths plus file bytes for each extracted source tree. The committed
JSON lock records package name, version, URL, archive SHA-256, and expected
source-tree digest. Configure fails if:

- a downloaded archive misses its `URL_HASH`;
- an existing `external/<package>` tree differs from the committed digest;
- a package/version/URL differs from the JSON lock.

Each configure-time integrity failure emits the stable
`NATIVE_DEPENDENCY_INTEGRITY_ERROR` code with exact package and
expected/observed digest. Task 1 maps runtime exposure to
`NativeDependencyIntegrityError`.

`BuildProvenance` later consumes these actual source-tree digests. It must not
claim that `pixi.lock` covers CMake downloads.

### Step 3: generate and validate both locks

Run:

```bash
pixi lock
pixi install --locked
pixi lock --check
pixi run native-lock-check
pixi run python -c "from compas_cgal import _stock_2, _toolpath"
pixi run baseline
pixi run affected
```

If editable native rebuild does not trigger after a touched C++ source, correct
the scikit-build configuration and repeat. Do not introduce pip/conda commands.

### Step 4: verify and commit

```bash
pixi run lint
pixi lock --check
pixi run native-lock-check
pixi run affected
git diff --check
```

Commit: `build: add locked pixi workflow`

---

## Task 1 — Establish typed units, policies, exact motions, and cap ownership

**Files**

- Create: `src/compas_cgal/adaptive/errors.py`
- Create: `src/compas_cgal/adaptive/units.py`
- Create: `src/compas_cgal/adaptive/motion.py`
- Create: `src/compas_cgal/adaptive/policy.py`
- Create: `src/compas_cgal/adaptive/__init__.py`
- Create: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/engagement_2.h`
- Modify: `src/engagement_2.cpp`
- Test: `tests/adaptive/test_units.py`
- Test: `tests/adaptive/test_motion.py`
- Test: `tests/adaptive/test_policy.py`

### Step 1: write RED contracts

Tests require:

- `Point2[WorldXY]` and `Vector2[WorldXY]`;
- scalar and sequence overloads for both vector factories;
- `ToolRadius`, `EntryRadius`, `Spacing`, `ChordBound`, `CutZ`, and
  `ClearanceZ` reject non-finite or invalid values through
  `InvalidUnitValueError`; `CutPlane.build(cut_z, clearance_z)` requires
  `clearance_z > cut_z`;
- `ExactSegmentMotion.build(a,b)` rejects zero progress through
  `DegenerateSegmentMotionError`;
- `ExactCircleMotion.build(center, phase_vector, clockwise)` derives squared
  radius and rejects a zero phase through `DegenerateCircleMotionError`;
- `EngagementCap.build(theta)` calls the native cap conversion and retains the
  exact returned big-endian bytes;
- `CandidatePolicy.build(...)` defines every finite bound and canonical
  tie-break;
- `NeckPolicy.build(...)` canonicalizes exact squared-width class boundaries
  and finite `(neck_class, passage_state) -> effective cap` mappings;
- `DepletionPolicy.build(...)` owns exact chord-bound bytes and a positive
  center-count limit;
- `TraversalPolicy.build(...)` owns deterministic component/branch order and a
  finite forward window;
- `CutDirectionPolicy.build(...)` requires explicit climb or conventional
  intent and maps material side to circle orientation with no default;
- the four policy factories raise `InvalidCandidatePolicyError`,
  `InvalidNeckPolicyError`, `InvalidDepletionPolicyError`, and
  `InvalidTraversalPolicyError`; cut direction raises
  `InvalidCutDirectionPolicyError`;
- every effective cap is produced by the native cap boundary and proven no
  larger than the user cap;
- no reporting clearance, magic tolerance, or raw string decides a policy
  branch.

Run RED:

```bash
pixi run pytest tests/adaptive/test_units.py tests/adaptive/test_motion.py \
  tests/adaptive/test_policy.py -n auto -q
```

### Step 2: implement the native cap seam

Add:

```cpp
double cap_chord_ratio(double cap_radians);
bool cap_chord_ratio_le(double lhs, double rhs);
```

`certify_segment_tea` and the binding call this one implementation. Bind both
functions on `_stock_2`; `cap_chord_ratio_le` injects both binary64 values into
Epeck before comparison. Reject angles outside `(0, pi]` before `sin`.

### Step 3: implement the Python domain

Use frozen dataclasses only where frame/unit or cross-field invariants are real.
Factories raise the named errors in the governing spec. Keep
`adaptive/__init__.py` minimal and do not add `__all__`.

The `.pyi` file covers every `_stock_2` symbol used by the strict package; no
blanket `Any` or ignore is accepted.

### Step 4: GREEN gates

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_units.py tests/adaptive/test_motion.py \
  tests/adaptive/test_policy.py -n auto --testmon -q
pixi run pytest tests/test_stock.py tests/test_engagement_oracle.py -n auto -q
```

Commit: `feat(adaptive): add exact motion types`

---

## Task 1A — Establish canonical encoding and identity primitives

**Files**

- Create: `src/compas_cgal/adaptive/canonical.py`
- Create: `src/compas_cgal/adaptive/identity.py`
- Create: `src/compas_cgal/adaptive/operation.py`
- Test: `tests/adaptive/test_canonical.py`
- Test: `tests/adaptive/test_identity.py`
- Test: `tests/adaptive/test_operation.py`

### Step 1: write RED identity contracts

Tests require:

- versioned big-endian encoders for integers, binary64 values,
  `ExactRationalV1` normalized numerator/positive-denominator bytes, sequences,
  tagged unions, and component maps;
- signed-zero normalization and loud rejection of NaN/infinity;
- canonical outer-ring start/orientation and hole ordering;
- canonical bytes for every Task 1 cut plane, motion, cap, and policy;
- frozen `INPUT_SCHEMA_VERSION`, `OPERATION_SCHEMA_VERSION`, and component
  version tags available before any witness is emitted;
- `BoundaryVertexIdV1` as either a canonical input-ring vertex ID or sorted
  normalized incident-support IDs plus the exact lexicographic intersection
  ordinal and trim-incidence orientation; no raw `Sqrt_extension`
  representation is hashed;
- a closed `CanonicalOperation` tagged union for approach, plunge, exact
  segment, and exact full-circle operations; approach/plunge bind
  clearance/cut Z, while each lateral operation binds geometry, semantic kind,
  cut Z, stable neck owner/orientation, `EffectiveCapDecision`, and
  `TraversalDecision`;
- a closed `EffectiveCapDecision` union: full-cap decisions bind equal
  user/effective bytes, while neck-cap decisions bind neck-evidence digest,
  exact width-class ID, passage before/after state, and user/effective bytes;
- `TraversalDecision` binds component/edge/branch IDs, exact cursor
  before/after identities, and whether that transition makes the cursor
  terminal;
- motion and depletion witnesses remain separate keyed streams; they are never
  embedded in or trusted as part of `CanonicalOperation`;
- `ComponentIdentity.build(...)` binds strategy version, source revision,
  native source-tree digest, and canonical parameter bytes;
- one-bit changes and reordered-but-equivalent polygon encodings have named
  expected identity behavior;
- no identity includes its own digest.

Run RED:

```bash
pixi run pytest tests/adaptive/test_canonical.py \
  tests/adaptive/test_identity.py tests/adaptive/test_operation.py -n auto -q
```

### Step 2: implement and GREEN

All later traces and witnesses import these encoders; no task invents a local
serialization. Validated `InputIdentity` is added with `PreclearedEntry` in
Task 11. `BuildProvenance` and final `ArtifactIdentity` remain incomplete until
Task 14.

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_canonical.py \
  tests/adaptive/test_identity.py tests/adaptive/test_operation.py \
  -n auto --testmon -q
```

Commit: `feat(identity): add canonical primitives`

---

## Task 2 — Add exact-on-guide conservative depletion

**Files**

- Create: `src/exact_motion_2.h`
- Create: `src/exact_depletion_2.h`
- Create: `src/exact_depletion_2.cpp`
- Modify: `CMakeLists.txt`
- Modify: `src/stock_2.h`
- Modify: `src/stock_2.cpp`
- Modify: `src/compas_cgal/stock.py`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/adaptive/canonical.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Create: `src/compas_cgal/adaptive/stock_area.py`
- Test: `tests/adaptive/test_exact_depletion.py`
- Test: `tests/adaptive/test_stock_monotonicity.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

### Step 1: write RED exact-invariant tests

Expose test diagnostics from the new binding, not NumPy-round-tripped
coordinates. Require:

- every segment center is exactly collinear and ordered between endpoints;
- both endpoints occur;
- every circle center has exact squared distance `v dot v` from `c`;
- phase, quarter points, antipode, and seam interval occur;
- every removal radius is positive and `<= tool_radius`;
- exact consecutive chord bounds include the circle seam;
- CW and CCW share a removal set but have different ordered witnesses;
- a rounded-interpolation/rotation test strategy fails its exact incidence
  predicate;
- `Stock2.clone()` is independent;
- rejected trial depletion leaves authoritative stock exactly equal to its
  snapshot;
- the center-count limit fails before allocation or mutation;
- the stock-inclusion/cap monotonicity fixture is non-vacuous: the larger stock
  exceeds the cap while the smaller stock does not.

Use representable fixtures to assert `U \ W` is exactly empty for Pythagorean
segments and circles. Add the exact-predicate metamorphic implication
`cap_exceeded(S) => cap_exceeded(S_hat)` for `S subseteq S_hat`.

Run RED:

```bash
pixi run pytest tests/adaptive/test_exact_depletion.py \
  tests/adaptive/test_stock_monotonicity.py -n auto -q
```

### Step 2: implement exact native construction

Add alongside legacy methods:

```cpp
DepletionTrace subtract_exact_segment(
    const ExactSegmentMotion2&, const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord, std::size_t center_count_limit);

DepletionTrace subtract_exact_full_circle(
    const ExactCircleMotion2&, const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord, std::size_t center_count_limit);
```

Segment centers use exact `Epeck::FT` barycentric parameters.
Circle centers use four exact Pythagorean quarter charts. Refine until exact
consecutive squared-chord comparisons, including the cyclic seam, pass.
Refinement is by powers of two and checks the complete center count with
overflow-safe arithmetic before allocation. No trig, normalization,
`to_double`, or double interpolation is permitted.

`DepletionTrace` contains deterministic structural parameter indices, count,
exact max-chord input, removal radius, and strategy version. Native code does
not reproduce or echo Task 1A's CCAN identity grammar. Python computes the
`DepletionWitness` from canonical motion, full depletion policy, tool radius,
ordered structural center parameters, native strategy version, and parent
lineage.

`Stock2` owns the authoritative `Gps` behind `std::unique_ptr<Gps>`.
`Stock2::clone()` uses the exact deep-copy constructor. Exact depletion builds
and subtracts on a completed clone, validates the structural trace, then
commits with a no-throw pointer swap; `Gps::operator=` is not an atomic commit.
Expose native exact subset/equality diagnostics through copied differences,
not sampled `contains()` calls. Bind new methods under distinct names and leave
`subtract_capsule`/`subtract_arc_sweep` unchanged.

### Step 3: implement the typed owner

`Stock2Area.build(...)` owns stock lifetime. `fork()` clones. `deplete(...)`
has explicit `@overload` signatures for `ExactSegmentMotion` and
`ExactCircleMotion`, returning a `DepletionWitness`; the single implementation
dispatches without a second public calling convention. It does not certify
motion. It depletes a raw clone, validates and canonicalizes the native trace,
constructs the immutable witness, then replaces `_raw` and appends lineage only
after every preceding step succeeds. Task 11 adds the validated-entry overload.

### Step 4: GREEN and regression

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_exact_depletion.py \
  tests/adaptive/test_stock_monotonicity.py -n auto --testmon -q
pixi run pytest tests/test_stock.py tests/test_engagement_oracle.py \
  tests/test_engagement_audit.py -n auto -q
```

Commit: `feat(stock): add exact-on-guide depletion`

---

## Task 3 — Build the exact reachable domain and full-sweep coverage ledger

**Files**

- Create: `src/reachable_domain_2.h`
- Create: `src/reachable_domain_2.cpp`
- Create: `src/coverage_2.h`
- Create: `src/coverage_2.cpp`
- Modify: `CMakeLists.txt`
- Create: `src/compas_cgal/_coverage_2.pyi`
- Create: `src/compas_cgal/adaptive/reachable_domain.py`
- Create: `src/compas_cgal/adaptive/coverage.py`
- Test: `tests/adaptive/test_reachable_domain.py`
- Test: `tests/adaptive/test_coverage.py`

### Step 1: compile gate

Before public code, compile a minimal target using
`Exact_predicates_exact_constructions_kernel_with_sqrt` and
`Gps_circle_segment_traits_2` that constructs:

- `C_r = D ⊖ B_r` from the arrangement of exact edge-offset lines and
  vertex-radius circles, with cells selected by exact disk containment;
- `M_r = C_r ⊕ B_r` and `D \ M_r`;
- a capsule for a non-axis-aligned rational segment;
- an annulus whose guide radius is `sqrt(v dot v)`;
- exact union and difference against a polygon.

Run it through a Pixi task. If the selected traits cannot represent or compare
these curves exactly under the locked build, stop Task 3 and all coverage-
dependent tasks. Do not substitute disk chains or sampled residual area.

### Step 2: write RED coverage contracts

Tests require:

- rectangle, acute convex corner, reflex corner, narrow neck, island, and
  disconnected-`C_r` fixtures; the disconnected case raises
  `PocketNotMachinableError` because Phase 1 has one qualified entry;
- exact `C_r` disk-containment at every selected cell and rejection immediately
  across every boundary;
- exact `M_r subseteq D`, exact `D \ M_r`, and a nonempty convex-corner
  residual for positive tool radius;
- exact reachable-domain construction is insertion-order invariant and its
  certificate binds every source curve, selected cell, component, and digest;
- exact capsule and annulus membership at endpoints, sides, inner/outer radii,
  and seams;
- `guide_radius <`, `==`, and `>` tool radius;
- exact residual non-empty before complete sweeps and empty after a committed
  complete fixture;
- operation order and duplicate sweeps do not alter set semantics but do alter
  lineage only where the canonical stream differs;
- under-cover engagement stock remains non-empty in a wall-tangent fixture
  while the exact full-sweep ledger proves coverage;
- removing the final sweep leaves a named residual component;
- a rounded sweep construction fails an exact boundary fixture.

Run RED:

```bash
pixi run pytest tests/adaptive/test_reachable_domain.py \
  tests/adaptive/test_coverage.py -n auto -q
```

### Step 3: implement

Bind:

```python
ReachableDomain2(design_boundary, holes, tool_radius)
ReachableDomain2.center_domain()
ReachableDomain2.reachable_material()
ReachableDomain2.unreachable_residual()
ReachableDomain2.certificate()
Coverage2(reachable_material, precleared_center, precleared_radius)
Coverage2.clone()
Coverage2.add_segment_sweep(...)
Coverage2.add_full_circle_sweep(...)
Coverage2.residual_is_empty()
Coverage2.residual_component_count()
```

The native owner accumulates exact mathematical sweeps against `M_r`, not
under-cover chains and not the sharp-corner design pocket `D`.
`ReachableDomain.build(...)` validates the native certificate and owns
`C_r`, `M_r`, and `D \ M_r`. `CoverageLedger` binds ordered sweep witnesses to
the same canonical motions used by depletion and binds the unreachable-residual
digest separately.

### Step 4: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_reachable_domain.py \
  tests/adaptive/test_coverage.py -n auto --testmon -q
```

Commit: `feat(coverage): add exact reachable ledger`

---

## Task 4 — Commit the event corpus and exact stock-boundary extraction

**Files**

- Create: `src/continuous_tea_2/result.h`
- Create: `src/continuous_tea_2/boundary_events.h`
- Create: `src/continuous_tea_2/boundary_events.cpp`
- Create: `src/continuous_tea_2.cpp`
- Modify: `CMakeLists.txt`
- Create: `src/compas_cgal/_continuous_tea_2.pyi`
- Create: `tests/adaptive/fixtures/event_corpus.json`
- Create: `tests/adaptive/test_motion_counterexamples.py`
- Create: `tests/adaptive/test_event_substrate.py`

### Step 1: encode the RED corpus

The fixtures cover:

1. a subtracted tool disk coincident with the cutter rim;
2. dyadic one-sided offsets proving engagement `> pi`;
3. near-coincident equal circles with fixed endpoint rotation and vanishing
   center travel;
4. external/internal tangencies;
5. line/circle and circle/circle regularization vertices;
6. positive-length supporting-curve overlap;
7. two-disk run merge/split;
8. the existing slotted gap-merge witness;
9. three simultaneous intersections;
10. segment and full-circle parameter seams;
11. a quadratic circle-circle endpoint whose `BoundaryVertexIdV1` is stable
    across operand and iterator order.

One test explicitly asserts that the current center-distance growth expression
is not a valid universal endpoint bound. This is a regression fact, not an
xfail.

Run RED:

```bash
pixi run pytest tests/adaptive/test_motion_counterexamples.py \
  tests/adaptive/test_event_substrate.py -n auto -q
```

### Step 2: expose exact boundary records

Extract each trimmed `Stock2` halfedge with:

- supporting line/circle coefficients;
- exact endpoints;
- material side;
- exact feature ID and trimmed-domain predicate;
- exact overlap and vertex incidence.

The regularized polygon-set boundary does not retain removal-source identity.
Depletion lineage records that identity separately when the removal is applied;
boundary extraction never attempts to reconstruct it. Event construction
intersects the supporting curve **and** its trimmed-domain predicate so an
intersection on the infinite support but outside the halfedge cannot create a
spurious event.

Feature IDs come from the Task 1A canonical encoding of material side,
normalized support coefficients, oriented `BoundaryVertexIdV1` trim endpoints,
and overlap multiplicity, then are sorted by canonical bytes. Task 4 constructs
each intersection ordinal using exact CircularKernel lexicographic comparison.
They never depend on a raw one-root representation, CGAL iterator order, or
object addresses.

Do not flatten positive-length overlaps into tangent events. The low-level
binding exposes test diagnostics plus:

```python
BoundaryEventKind = Literal[
    "transverse", "tangent", "vertex", "overlap", "seam"
]
```

The result contract is:

```cpp
enum class ContinuousTeaVerdict {
    CERTIFIED,
    CAP_EXCEEDED,
    UNRESOLVED_DEGENERACY,
};
```

No bool compatibility mapping is used in the new module.

### Step 3: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_motion_counterexamples.py \
  tests/adaptive/test_event_substrate.py -n auto --testmon -q
pixi run pytest tests/test_stock.py tests/test_engagement_oracle.py -n auto -q
```

Commit: `feat(tea): expose continuous events`

---

## Task 5 — Compile and validate the algebraic event substrate

**Files**

- Create: `src/continuous_tea_2/parameter_charts.h`
- Create: `src/continuous_tea_2/parameter_charts.cpp`
- Create: `src/exact_algebraic_1.h`
- Create: `src/exact_algebraic_1.cpp`
- Create: `src/continuous_tea_2/event_partition.h`
- Create: `src/continuous_tea_2/event_partition.cpp`
- Create: `src/continuous_tea_2/partition_certificate.h`
- Modify: `CMakeLists.txt`
- Test: `tests/adaptive/test_event_substrate.py`

### Step 1: compile gate

Under the actual locked CGAL 6.0.1 build, compile:

- `Arr_algebraic_segment_traits_2`;
- `CGAL::Algebraic_kernel_d_1<CORE::BigInt>` after clearing exact rational
  polynomial denominators to primitive integer coefficients;
- `CGAL::Algebraic_kernel_d_2<CORE::BigInt>` for bivariate arrangement and
  elimination operations;
- exact real-root isolation, comparison, and sign-at-root through that kernel;
- sign-at-root;
- rational projective half-angle charts;
- primitive/square-free normalization over the compiled exact rational
  coefficient type;
- the required CORE/Boost backend without re-enabling an undeclared GMP path.

Record the exact traits and compile definitions in the code and
`BuildProvenance`. If this gate fails, stop Tasks 5–16 and report the unsupported
number-field operation. No sampled fallback is permitted.

`exact_algebraic_1` defines the only algebraic witness encoding:

```text
AlgebraicRootIdV1 =
  version tag
  + primitive square-free integer coefficients normalized to gcd 1
    and positive leading coefficient
  + ordinal among all ordered real roots
```

The native verifier reruns `Solve_1`, checks the ordinal, and obtains/refines a
rational isolating interval. Native comparisons use `Compare_1`/`Sign_at_1`;
identity never uses `CORE::Expr` streaming or a decimal approximation.

### Step 2: RED chart and root contracts

Tests require:

- line motion `c(t) = a + t(b-a)`, `t in [0,1]`;
- four charts for a full circular center motion;
- two rim half-angle charts plus explicit seams;
- exact pullback polynomials for line and circular stock supports;
- enforced bidegree bounds: segment/line `(1,2)`, segment/circle `(2,4)`,
  full-circle/line `(2,2)`, and full-circle/circle `(4,4)`;
- exact trimmed-halfedge domain predicates on every pullback branch;
- exact isolation/order of repeated and simultaneous event roots;
- an identically zero polynomial classifies as overlap;
- every positive-width open parameter cell has a rational witness;
- every zero-dimensional event fibre retains its exact algebraic root and
  isolating interval;
- two-branch cap crossings are isolated from `F_i(t,u) = 0`,
  `F_j(t,v) = 0`, and the exact squared-chord cap relation;
- the exact oriented-CCW run branch used by `run_exceeds_cap` separates
  `theta` from its `2*pi - theta` complement and partitions determinant-zero
  orientation boundaries;
- isolated equal-support fibres arise from common zeros of every rim-polynomial
  coefficient in `t`, not only from a globally zero pullback;
- circle-circle regularization vertices are projected by rational elimination
  through both incident supports and conjugate-filtered by exact trim/vertex
  identity;
- resultant roots are filtered against all original equations and trim
  predicates, while an identically-equal cap interval is represented
  explicitly;
- deleting a seam or coalescing two exact roots fails a named fixture.

Run RED:

```bash
pixi run pytest tests/adaptive/test_event_substrate.py -n auto -q
```

### Step 3: implement event partition primitives

Keep algebraic curve construction, root isolation, cap-crossing elimination,
extraneous-root filtering, and cell ordering separate. No material/run logic
enters these files. The substrate gate is not complete until the minimal
two-moving-branch resultant fixture passes without station sampling.

Implement and test this derivation table:

| Event class | Exact projection | Completeness/degeneracy record |
|---|---|---|
| tangency | `resultant_u(F, partial_u F)` | square-free factors; repeated root once |
| trimmed vertex | rational elimination of rim plus incident supports | conjugates filtered by endpoint ID and trim truth |
| support overlap | common-zero projection of all rim coefficients in `t` | isolated overlap fibre versus motion interval |
| endpoint order/merge | resultant of two active endpoint branches | both original equations and trims rechecked |
| orientation boundary | eliminate branches with exact CCW determinant zero | zero versus `pi` separated by exact dot sign |
| cap crossing | eliminate `u,v` from `F_i`, `F_j`, squared-chord relation plus CCW branch | complement roots filtered; identical interval has orientation disposition |
| chart seam | exact finite/infinite boundary evaluation | exactly one canonical owner |

`EventPartitionCertificate` records:

- chart IDs and exact chart-domain coverage;
- normalized coefficient bytes, actual bidegree, square-free factors, and a
  reference to the applicable degree bound;
- every isolated real root with exact isolating interval and multiplicity;
- every positive-width sign-invariant cell and its rational witness;
- every zero-dimensional fibre and algebraic-root identity;
- event kind, trim/feature/branch IDs, exact disposition, and seam owner.

The native verifier reconstructs every projected factor, reruns complete
real-root isolation, checks union of chart domains including seams, and checks
that ordered roots induce exactly the recorded cells/fibres. Any mismatch is
`UNRESOLVED_DEGENERACY`; `event_cell_count` is never accepted as completeness
evidence.

Canonical certificate bytes bind `AlgebraicRootIdV1`, not solver-selected
isolating intervals. Intervals are regenerated exact diagnostics and excluded
from proof identity, preserving fresh-process determinism.

The compile/RED gate includes an isolated equal-circle coincidence, a
positive-width motion overlap, and a circle-circle regularization vertex with
quadratic coordinates. These prevent the rational coefficient contract from
being bypassed through algebraic-coordinate substitution.

### Step 4: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_event_substrate.py -n auto --testmon -q
```

Commit: `feat(tea): add exact event partition`

---

## Task 6 — Implement the event-exact segment oracle

**Files**

- Create: `src/continuous_tea_2/segment_oracle.cpp`
- Modify: `src/continuous_tea_2.cpp`
- Modify: `src/compas_cgal/_continuous_tea_2.pyi`
- Test: `tests/adaptive/test_segment_oracle.py`

### Step 1: write RED decision tests

Required cases:

- line half-plane: exact `pi` equality certifies, tighter cap exceeds;
- fully clear and fully material rim;
- line/circle and circle/circle tangency during segment motion;
- cutter-rim passage through a stock vertex;
- run birth/death and merge/split;
- coincident tool-disk boundary returns resolved conservative engagement or
  `UNRESOLVED_DEGENERACY`, never empty by omission;
- cap equality at an algebraic event root;
- zero-length motion rejected by the typed factory before native dispatch;
- omitting tangency, overlap, vertex, merge, or cap roots kills one named test.

Run RED:

```bash
pixi run pytest tests/adaptive/test_segment_oracle.py -n auto -q
```

### Step 2: implement

Bind:

```python
audit_segment_tea_event_exact(
    stock,
    x0,
    y0,
    x1,
    y1,
    tool_radius,
    cap_chord_ratio,
) -> tuple[Literal["certified", "cap_exceeded", "unresolved"], EventTrace]
```

For every open parameter cell and event fibre:

1. assemble exact engaged rim intervals from material-side records;
2. merge exact cyclic runs;
3. isolate every run/cap equality through the Task 5 two-branch elimination;
4. decide one exact sign per sign-invariant cell;
5. retain the first canonical violation or unresolved event.

Reporting intervals never decide the verdict.

`EventTrace` has a canonical order independent of insertion order. Its digest
binds the motion chart and `AlgebraicRootIdV1`, event kind, trimmed feature IDs,
active branch IDs, canonical effective-cap bytes, exact verdict, oracle
strategy version, and verified
`EventPartitionCertificate` digest. `event_cell_count` is reporting metadata,
never proof identity.

### Step 3: differential falsifier

For every certified result, probe the committed dyadic corpus with the full
exact `engagement_at` cap and zero gap closure. A station violation is a
definitive oracle bug. A clean scan is only a falsifier, not proof.

### Step 4: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_motion_counterexamples.py -n auto --testmon -q
```

Commit: `feat(tea): certify segment events exactly`

---

## Task 7 — Implement the event-exact full-circle oracle

**Files**

- Create: `src/continuous_tea_2/circle_oracle.cpp`
- Modify: `src/continuous_tea_2.cpp`
- Modify: `src/compas_cgal/_continuous_tea_2.pyi`
- Test: `tests/adaptive/test_circle_oracle.py`

### Step 1: write RED circle contracts

Require:

- full four-chart coverage and exact seam identity;
- CW/CCW same maximum verdict, distinct ordered trace;
- coincident and near-coincident prior tool disks;
- external/internal tangency;
- two-disk merge/split and existing slotted circular witness;
- contour vertex and simultaneous-event fibres;
- exact cap equality;
- mutation removal of any chart, seam, coincidence, merge, or cap root fails;
- certified event result has no violating dyadic probe through the declared
  falsifier depth.

Run RED:

```bash
pixi run pytest tests/adaptive/test_circle_oracle.py -n auto -q
```

### Step 2: implement

Bind:

```python
audit_full_circle_tea_event_exact(
    stock,
    center_x,
    center_y,
    phase_dx,
    phase_dy,
    clockwise,
    tool_radius,
    cap_chord_ratio,
) -> tuple[Literal["certified", "cap_exceeded", "unresolved"], EventTrace]
```

The nonzero phase vector, not a rounded radius, defines the guide circle. A
zero vector raises `DegenerateCircleMotionError` before native dispatch.
Chart/root identity defines seam ownership; equivalent seam events are
deduplicated by exact identity, and the trace is sorted by oriented motion
order followed by canonical event identity. CW/CCW order therefore cannot
depend on chart insertion order.

### Step 3: GREEN and exact-kernel regression

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_segment_oracle.py -n auto --testmon -q
pixi run pytest tests/test_stock.py tests/test_engagement_oracle.py \
  tests/test_engagement_audit.py -n auto -q
```

Commit: `feat(tea): certify circle events exactly`

---

## Task 8 — Add the typed motion certifier

**Files**

- Create: `src/compas_cgal/adaptive/motion_certificate.py`
- Test: `tests/adaptive/test_motion_certificate.py`

### Step 1: write RED public contracts

```python
@dataclass(frozen=True)
class MotionWitness:
    operation_index: int
    operation_kind: OperationType
    motion: ExactSegmentMotion | ExactCircleMotion
    user_cap_bytes: bytes
    effective_cap_bytes: bytes
    strategy_identity: bytes
    stock_lineage_digest: bytes
    event_trace_digest: bytes
    verdict: Literal["certified"]
    event_cell_count: int
    unresolved_count: int


@dataclass(frozen=True)
class MotionCertifier:
    ...

    def certify(
        self,
        *,
        operation_index: int,
        operation_kind: OperationType,
        motion: ExactSegmentMotion | ExactCircleMotion,
        user_cap: EngagementCap,
        effective_cap: EngagementCap,
    ) -> MotionWitness: ...
```

`certify` proves `effective_cap <= user_cap` with the exact cap surrogate before
native dispatch. The canonical effective cap, not a generator-provided witness,
is the cap passed to the native oracle.

Tests require:

- certified -> immutable witness;
- proved cap exceedance -> `EngagementCapExceededError`, which only the
  candidate evaluator converts to ordinary rejection;
- unresolved -> `UnresolvedMotionEventError`;
- the certifier exposes a const stock API and leaves lineage digest, canonical
  boundary digest, and exact stock queries unchanged;
- negative index, non-lateral kind, kind/motion incompatibility, and
  `effective_cap > user_cap` are rejected locally;
- the returned witness binds the supplied ordinal, semantic kind, exact motion,
  caps, and observed pre-state digest;
- pre-deplete virgin slotting fails;
- no plunge/approach witness exists.

Run RED:

```bash
pixi run pytest tests/adaptive/test_motion_certificate.py -n auto -q
```

### Step 2: implement and GREEN

The certifier is read-only and dispatches only exact segment/full-circle
motions. It never calls the legacy guarded certifiers. Fresh orchestration is
deferred until Task 11A, after validated entry and exact neck state exist.
Cross-record ordinal/motion/cap/state mismatch belongs to replay and artifact
assembly, where an expected canonical operation and keyed witness stream exist.

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_motion_certificate.py \
  -n auto --testmon -q
```

Commit: `feat(adaptive): add event-exact certifier`

---

## Task 9 — Add the true segment-site Voronoi medial axis

**Implementation dependencies:** Task 3's verified `C_r`/digest and Task 5's
`exact_algebraic_1` root encoding/backend.

**Files**

- Create: `src/segment_site_mat.h`
- Create: `src/segment_site_mat.cpp`
- Create: `src/segment_site_mat_sampling.h`
- Create: `src/segment_site_mat_sampling.cpp`
- Create: `src/segment_site_neck.h`
- Create: `src/segment_site_neck.cpp`
- Create: `src/medial_axis_2.cpp`
- Modify: `CMakeLists.txt`
- Create: `src/compas_cgal/_medial_axis_2.pyi`
- Test: `tests/adaptive/test_medial_axis.py`

### Step 1: write RED topology/provenance tests

Require:

- the committed analytically verified point/segment-bisector fixture produces a
  retained `PARABOLA` edge, proving no straight-skeleton downgrade;
- native exact verdict fields prove every sample equidistant to its recorded
  generator sites and no site closer; returned doubles are never rechecked as
  exact evidence;
- exact polygon clipping against `D` plus exact clearance clipping
  `distance_to_defining_site^2 >= r^2` removes exterior, hole-interior, and
  tool-center-inadmissible portions; the result equals intersection with `C_r`;
- lines, rays, segments, and parabolas that cross the domain multiple times
  emit every connected interior component;
- an unbounded dual with bounded interior pieces retains those pieces;
- tangential clips, hole clips, and multiple-crossing clips retain exact
  endpoint/boundary provenance;
- exact node collapse and adjacency use stable IDs, never addresses/doubles;
- incident point/open-segment feature-transition duals are rejected exactly
  when the point generator is that segment's source or target, while
  nonincident point/segment parabolas remain;
- degeneracy-normalized nodes collect the union of generator-site provenance
  from every incident halfedge, not one underlying face triple;
- every neck record is an exact local site-distance minimum, binds both sites
  and the separating graph cut, and changes only when exact squared-width
  thresholds are crossed;
- strict interior minima, maximal constant-clearance plateaus, one-sided
  `C_r` endpoint minima, and shared-vertex minima have one deterministic owner;
- a clearance polynomial that is identically zero retains the complete dual
  cell as a closed `C_r`-boundary plateau; constant positive retains and
  constant negative rejects;
- equal minima meeting at a vertex merge once, while a non-separating minimum
  is not a neck;
- graph cycle rank equals the exact first Betti number of `C_r`; narrow
  corridors may disconnect components or remove a design-pocket hole cycle;
- sampling refinement changes samples but not analytic topology/provenance;
- insufficient conic refinement raises `ConicSamplingLimitError`;
- repeated calls are bitwise deterministic;
- existing `tests/test_toolpath.py` output remains unchanged.

Run RED:

```bash
pixi run pytest tests/adaptive/test_medial_axis.py -n auto -q
```

### Step 2: implement the exact graph

First compile a native spike that parameterizes every Voronoi primitive
(`Line_2`, `Ray_2`, `Segment_2`, and `Parabola_segment_2`), substitutes it into
the defining-site squared-clearance equation, clears denominators to a
primitive integer polynomial of degree at most four, and isolates/orders every
`clearance^2 = r^2` root with Task 5's
`CGAL::Algebraic_kernel_d_1<CORE::BigInt>`. Exact line-boundary intersections
clip against `D`; the clearance polynomial clips to `C_r`. This deliberately
does not assume a `Parabola_segment_2`/circle intersection overload. If any
primitive cannot produce and order those roots exactly, stop Tasks 9–16. Do not
sample the clip.

Constant clearance is handled before square-free normalization: exact positive
retains the whole primitive domain, exact negative discards it, and exact zero
retains it as a boundary plateau with no fabricated roots.

Use:

```cpp
using MatKernel =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using MatTraits =
    CGAL::Segment_Delaunay_graph_traits_without_intersections_2<
        MatKernel, CGAL::Field_with_sqrt_tag>;
```

Use the segment-Delaunay degeneracy-removal Voronoi adaptor. Insert normalized
ring segments with stable feature IDs. Use
`Segment_Delaunay_graph_storage_traits_with_info_2` with explicit
conversion/merge functors, or an exact post-map proved bijective by native
tests; default storage does not preserve caller IDs. Provenance distinguishes
vertex point-sites from open-segment sites. The adaptor provides generator and
limiter handles through `up()`, `down()`, `left()`, and `right()` plus the
underlying SDG edge through `dual()`; recover stable caller data only through
`storage_site().info()`. Handles, addresses, and iterator order never enter
identity.

The Voronoi adaptor removes zero-length edges and zero-area faces, but not all
point-versus-incident-segment feature-transition bisectors. Before clipping,
discard exactly a point/open-segment dual whose point vertex ID equals that
segment's source or target vertex ID, and independently assert exact point
equality with the endpoint. Retain every nonincident point/segment parabola,
point/point dual, and segment/segment dual. At a degeneracy-normalized node,
construct node-site CSR from the union over all incident halfedges because one
`Vertex::site(i)` face triple may omit additional incident generators.

`SDG::primal()` supplies exact `Line_2`, `Ray_2`, `Segment_2`, or
`Parabola_segment_2` geometry. `Parabola_2::side_of_parabola`,
`Parabola_segment_2::compute_k`, `generate_points`, and streamed parabola
output are visualization APIs with `to_double` paths and are forbidden from
topology, clipping, identity, and certificate decisions.

Split each dual at exact `D`-boundary roots and exact clearance roots. Classify
every open parameter cell by point-in-`D` and clearance sign, then emit every
maximal admissible-center component. On a valid Voronoi dual the defining-site
distance is the nearest-site distance, so this clearance predicate is exactly
membership in `C_r`. Never discard a dual merely because its unclipped
primitive is unbounded. Each emitted edge binds original dual identity,
component index, both generator sites, and complete endpoint provenance.
Independently reconstruct `C_r` through the shared `reachable_domain_2`
component and require its digest to equal Task 3's digest for the same input.

A quartic root is not coerced into `MatKernel::Point_2`.
`AlgebraicCurvePointV1` stores original dual identity plus
`AlgebraicRootIdV1`. Same-curve ordering uses the root kernel; equality at a
shared original Voronoi vertex uses its node identity; any remaining
cross-curve equality uses Task 5's bivariate exact kernel. Coordinates are
approximated only for reporting samples. The compile spike must exercise exact
curve evaluation/sign, same-curve order, shared-vertex equality, and
cross-curve equality without a `CORE::Expr` stream conversion.

Canonicalize exact nodes lexicographically. Edge records contain source/target
node, `LINE|PARABOLA`, ordered generator-site IDs, original-dual identity, and
clip provenance.

Set target-local compile definitions:

```cmake
CGAL_USE_CORE=1
CGAL_CORE_USE_BOOST_BACKEND=1
```

Keep the repository's existing Boost multiprecision and GMP-disabled flags.

### Step 3: implement exact neck evidence and classification

`NeckEvidenceV1` is a versioned tagged union:

- `STRICT_EDGE`: relative-interior `AlgebraicRootIdV1` with derivative sign
  changing negative-to-positive;
- `CLEARANCE_ENDPOINT`: one-sided minimum at a `C_r` clip with full endpoint
  feature CSR;
- `SHARED_VERTEX`: one exact node plus all incident minimizing halfedges;
- `PLATEAU`: one maximal connected constant-clearance interval/subgraph with
  exact endpoint IDs.

Every variant binds its exact squared width as `AlgebraicRootIdV1`, defining
site IDs, owner ID, and canonical separating-cut partition. Ownership rules are:

1. relative-interior strict minima exclude endpoints;
2. all equal-width incident endpoint minima merge into the shared vertex;
3. constant-clearance cells merge to one maximal plateau owned by the
   lexicographically first edge/node identity;
4. a remaining `C_r` endpoint is a one-sided minimum only when exact signs are
   nondecreasing into the edge;
5. evidence is retained only when removing its point, vertex, or contracted
   plateau separates at least two nonempty traversal sides.

Bind:

```python
_medial_axis_2.validate_and_classify_necks(
    mat_certificate: bytes,
    neck_evidence: tuple[bytes, ...],
    squared_width_boundaries: tuple[bytes, ...],
) -> tuple[NDArray[int64], tuple[bytes, ...]]
```

The function revalidates every algebraic root/cut, compares exact algebraic
widths to rational policy boundaries natively, and returns class IDs plus
comparison-certificate bytes. Each boundary byte string is Task 1A's
`ExactRationalV1` numerator/denominator encoding. Boundaries are strictly
increasing and width classes are `[0,b0]`, `(b0,b1]`, ..., `(bn,+infinity)`, so
exact equality belongs to the narrower, more restrictive class. `CORE::Expr`
text and doubles are forbidden.

### Step 4: implement proposal sampling

Lines use exact barycentric samples. Parabolas use exact parameter bisection;
reported chord/sagitta select refinement. At the depth cap, fail loud.

Bind:

```python
_medial_axis_2.segment_site_medial_axis(
    vertices,
    holes,
    tool_radius,
    station_spacing,
    max_sagitta,
    max_refinement_depth,
)
```

`tool_radius` is the sole admissible-center radius and is exactly the `r` bound
into Task 3's `C_r` digest. No unused or second radial-clearance parameter
exists.

Return this fixed tuple:

1. `nodes: float64[V,3]`;
2. `edges: int64[E,8]` as source, target, curve kind, site A, site B,
   original-dual kind, original-dual ID, clip-component index;
3. `node_site_offsets: int64[V+1]`;
4. `node_site_ids: int64[K]`;
5. `site_provenance: int64[S,3]` as kind, ring, feature;
6. `edge_endpoint_provenance_flags: int64[E,2]` as an independent bitset of
   original Voronoi vertex, `D` clip, and clearance/`C_r` root, allowing
   coincident events at one endpoint;
7. `endpoint_feature_offsets: int64[2*E+1]`;
8. `endpoint_features: int64[J,5]` as domain kind, component, curve kind,
   source-site/ring, and derived-feature index; CSR retains every incident
   feature at tangencies and degenerate vertices;
9. `edge_exact_flags: int64[E,3]` as admissible-center component,
   source-provenance verified, and target-provenance verified;
10. `sample_centers: float64[N,3]`;
11. `sample_clearance: float64[N]`;
12. `sample_guide_radius: float64[N]`;
13. `sample_flags: int64[N,2]` as exact-equidistant and no-site-closer;
14. `edge_sample_offsets: int64[E+1]`;
15. `sample_parameter: float64[N]`;
16. `neck_evidence: tuple[bytes, ...]`;
17. `neck_cut_offsets: int64[Q+1]`;
18. `neck_cut_edge_ids: int64[L]`;
19. `center_domain_digest: bytes`;
20. `mat_certificate: bytes`.

All exact flag columns are native contract verdicts and must be one for every
returned edge/sample; tool fit is structural because edges are exactly clipped
to `C_r`. All floating fields remain proposal-only.

### Step 5: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_medial_axis.py -n auto --testmon -q
pixi run pytest tests/test_toolpath.py -n auto -q
```

Commit: `feat(mat): add exact segment Voronoi graph`

---

## Task 10 — Wrap MAT topology, neck state, and the finite candidate lattice

**Files**

- Create: `src/compas_cgal/adaptive/medial_axis.py`
- Create: `src/compas_cgal/adaptive/neck.py`
- Create: `src/compas_cgal/adaptive/candidates.py`
- Test: `tests/adaptive/test_medial_axis.py`
- Test: `tests/adaptive/test_neck.py`
- Test: `tests/adaptive/test_candidates.py`

### Step 1: write RED typed contracts

Tests require:

- CSR ranks/dtypes and stable topology validation;
- exact feature/node IDs survive without coordinate matching;
- maximal tool-fit runs retain parent edge identity and clipped endpoint state;
- native `NeckEvidence` bytes validate exact local-minimum, squared-width
  class, site-pair, and separating-cut identity;
- each oriented neck uses the literal passage states `unvisited`,
  `first_pass_complete`, `second_pass_complete`, and `terminal`;
- `NeckPolicy` maps `(exact width class, passage state)` to a native-produced
  effective cap, proves it no larger than the user cap, and never reads
  `sample_clearance`;
- `_medial_axis_2.validate_and_classify_necks` recomputes every class and
  comparison certificate from the MAT certificate and policy's exact squared
  boundaries;
- `EffectiveCapDecision.build(...)` binds the verified neck-evidence digest,
  class, before/after state, and selected cap; full-cap decisions prove equal
  cap bytes, while neck decisions use `_stock_2.cap_chord_ratio_le`;
- accepted neck operations record one legal before/after state transition;
- each candidate carries one deterministic proposed `TraversalDecision` from
  its current exact cursor to its candidate endpoint;
- one-sided MATHSM proposal from `(m,p,r)` yields `q`, midpoint `c`, and phase
  `q-c`;
- the complete lattice contains all declared spatial/radius/phase refinement
  levels;
- exhaustive small-lattice oracle equals production enumeration;
- repeated enumeration order is identical;
- tie-break is furthest progress, largest radius, then canonical identity;
- no bisection or monotonic-feasibility assumption occurs.

Run RED:

```bash
pixi run pytest tests/adaptive/test_medial_axis.py \
  tests/adaptive/test_neck.py tests/adaptive/test_candidates.py -n auto -q
```

### Step 2: implement and GREEN

Keep topology ownership in `medial_axis.py`, exact neck evidence/passage state
in `neck.py`, and pure finite enumeration in `candidates.py`. The candidate
identity binds exact neck evidence, before/after passage state, and effective
cap bytes plus the proposed traversal transition. This task semantically
validates the Task 1A operation decision records; it does not change the frozen
operation encoding.

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_medial_axis.py \
  tests/adaptive/test_neck.py tests/adaptive/test_candidates.py \
  -n auto --testmon -q
```

Commit: `feat(adaptive): add finite MAT candidates`

---

## Task 11 — Add precleared entry and exact gouge containment

**Files**

- Create: `src/containment_2.h`
- Create: `src/containment_2.cpp`
- Create: `src/containment_2_bindings.cpp`
- Modify: `CMakeLists.txt`
- Create: `src/compas_cgal/_containment_2.pyi`
- Create: `src/compas_cgal/adaptive/entry.py`
- Create: `src/compas_cgal/adaptive/containment.py`
- Modify: `src/compas_cgal/adaptive/stock_area.py`
- Modify: `src/compas_cgal/adaptive/identity.py`
- Test: `tests/adaptive/test_entry.py`
- Test: `tests/adaptive/test_containment.py`
- Modify: `tests/adaptive/test_identity.py`

### Step 1: write RED entry tests

Require:

- precleared disk exact containment in pocket and exclusion from islands;
- entry center belongs to `C_r`; a larger entry radius is separately proved
  contained in `D`;
- process evidence qualifies the bore across the complete
  `CutPlane.clearance_z -> CutPlane.cut_z` interval;
- process provenance and radius bind to input identity;
- `Stock2Area.deplete(PreclearedEntry)` is an explicit overload, subtracts the
  exact entry disk once, and returns the required keyed `DepletionWitness`;
- canonical approach at clearance Z and vertical plunge to cut Z travel through
  the declared void and remove no material;
- first full circle fits inside the entry disk and certifies against post-entry
  stock;
- one-tool-radius seed on a large square is rejected as unable to launch any
  capped non-degenerate lateral motion;
- arbitrary lateral motion cannot be relabelled as entry;
- `InputIdentity.build(...)` accepts only the validated entry and binds
  canonical `D`, frame, cut plane, tool, reachable-domain digest, user cap,
  cut direction, all candidate/neck/depletion/traversal policies,
  `INPUT_SCHEMA_VERSION`, `OPERATION_SCHEMA_VERSION`, and every current
  component version.

### Step 2: write RED containment tests

For circles, decide exact containment of the mathematical annular sweep—or disk
when `|v| <= r`—in `D`, equivalently exact guide-circle containment in `C_r`.
Do not require the filled outer disk: that conservative test wrongly rejects
valid annuli surrounding islands. For segments, decide exact capsule
containment in `D`, equivalently center-segment containment in `C_r`.

The sweep and center-domain formulations are contract-tested as equivalent.
The outer-disk predicate may be a rejection-only optimization but never the
acceptance authority.

Test equality, reflex vertices, island boundaries, and a one-rational-quantum
gouge mutation.

Run RED:

```bash
pixi run pytest tests/adaptive/test_entry.py \
  tests/adaptive/test_containment.py tests/adaptive/test_identity.py \
  -n auto -q
```

### Step 3: implement and GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_entry.py \
  tests/adaptive/test_containment.py tests/adaptive/test_identity.py \
  -n auto --testmon -q
```

Commit: `feat(adaptive): certify entry and input`

---

## Task 11A — Add independent fresh-state replay

**Dependencies:** Tasks 3, 8, 10, and 11.

**Files**

- Create: `src/compas_cgal/adaptive/replay.py`
- Test: `tests/adaptive/test_replay.py`

### Step 1: write RED replay contracts

Expose:

```python
def replay_certificate(
    *,
    input_identity: InputIdentity,
    pocket: Polygon,
    holes: Sequence[Polygon],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    user_cap: EngagementCap,
    entry: PreclearedEntry,
    operations: tuple[CanonicalOperation, ...],
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> ReplayCertificate: ...
```

Tests require:

- replay recomputes and matches `InputIdentity`, `C_r`/`M_r`, MAT certificate,
  exact neck evidence/classes, and cut-direction orientation from canonical
  inputs;
- pristine stock and coverage are constructed, then `PreclearedEntry` is
  applied before the first lateral operation;
- approach/plunge exactly match clearance/cut Z and the qualified bore, and
  never mutate 2D stock;
- operation grammar is exactly approach, plunge, then one or more lateral
  operations; every segment/circle phase endpoint is exactly continuous with
  its predecessor and every lateral Z equals `CutPlane.cut_z`;
- replay computes a canonical ordered-operation digest and index chain;
- for each lateral operation, replay reconstructs `EffectiveCapDecision` from
  exact neck evidence, current passage state, and `NeckPolicy`; it rejects a
  recorded class, evidence digest, state transition, or cap that differs;
- motion is event-exact certified against frozen pre-depletion stock before
  exact depletion and coverage mutation;
- operation orientation agrees with `CutDirectionPolicy`;
- every `TraversalDecision` is recomputed against the fresh MAT graph, advances
  the named cursor legally, and leaves every required cursor terminal;
- fresh coverage proves exact `M_r` residual emptiness and separately binds
  the unchanged `D \ M_r` digest;
- generator witness verdicts, stock snapshots, and coverage snapshots are not
  accepted as inputs;
- a state-dependent reorder, stale input identity, non-fresh state, omitted
  entry, broken phase continuity, illegal grammar, nonterminal traversal,
  nonempty residual, and certify-after-deplete mutations fail named tests.

`ReplayCertificate` binds input identity, canonical ordered-operation
digest/index chain, recomputed effective-cap decisions, fresh initial
stock/coverage digests, ordered fresh motion/depletion/sweep traces, terminal
traversal/stock/coverage digests, unreachable-residual digest, and the exact
empty-`M_r` verdict. A commuting reorder that remains physically valid changes
this digest and is detected when Task 14 validates `ArtifactIdentity`; replay
does not invent an expected output order from input identity.

Run RED:

```bash
pixi run pytest tests/adaptive/test_replay.py -n auto -q
```

### Step 2: implement and GREEN

Separate immutable `ReplayCertificate`, policy/effective-cap recomputation, and
the short-lived replay state. Do not mix replay into `stock_area.py`.

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_replay.py -n auto --testmon -q
```

Commit: `feat(cert): add independent replay`

---

## Task 12 — Implement joint transactional candidate evaluation

**Files**

- Create: `src/compas_cgal/adaptive/transaction.py`
- Test: `tests/adaptive/test_transaction.py`

### Step 1: write RED transaction tests

Each candidate evaluation must:

1. fork authoritative engagement stock and coverage ledger;
2. construct the direct phase-to-phase segment;
3. containment-check, event-certify, and exact-deplete the segment;
4. containment-check, event-certify, and exact-deplete the circle;
5. append both exact sweeps to forked coverage;
6. validate and record the candidate's exact neck-state and
   `TraversalDecision` transitions;
7. return an immutable accepted transaction or a proved rejection;
8. leave authoritative state unchanged until explicit commit.

Tests cover link pass/circle fail, link fail, unresolved event, gouge failure,
successful atomic commit, stale-parent digest, and deterministic winner choice.
A mutation that depletes before certify or advances neck state before commit
must fail.

Run RED:

```bash
pixi run pytest tests/adaptive/test_transaction.py -n auto -q
```

### Step 2: implement

Do not create a stateful god object. Separate:

- immutable `CandidateTransaction`;
- short-lived `CandidateEvaluator`;
- long-lived `GenerationState`;
- pure winner selection.

### Step 3: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_transaction.py -n auto --testmon -q
```

Commit: `feat(adaptive): add atomic candidate commit`

---

## Task 13 — Traverse all MAT components and generate a covered path

**Files**

- Create: `src/compas_cgal/adaptive/traversal.py`
- Create: `src/compas_cgal/adaptive/generator.py`
- Test: `tests/adaptive/test_traversal.py`
- Test: `tests/adaptive/test_acceptance.py`
- Create: `tests/adaptive/fixtures/tractable_pocket.json`

### Step 1: write RED traversal contracts

Require:

- one well-founded cursor per MAT edge/component;
- every accepted candidate advances exactly one cursor;
- shared node IDs determine branch transitions;
- holes preserve graph cycles and explicit visited-side state;
- every accepted neck crossing follows the exact first/second-passage state
  machine and records the selected effective cap;
- no `max_passes`, silent branch drop, or coordinate-based adjacency;
- exhaustion with exact residual raises `IncompletePocketCoverageError`;
- no feasible exact candidate raises the appropriate named cap/neck/transition
  error;
- empty, entry-only, first-circle-only, and dropped-final-branch mutations fail.

### Step 2: write RED end-to-end small fixture

The committed tractable fixture must return:

- precleared approach/plunge;
- at least one certified segment and two certified full circles;
- operation/witness order with no holes;
- terminal traversal state;
- exact `M_r` full-sweep residual empty with `D \ M_r` separately bound;
- fresh event-exact replay with zero cap violations, every cursor terminal, and
  exact empty `M_r` residual.

Run RED:

```bash
pixi run pytest tests/adaptive/test_traversal.py \
  tests/adaptive/test_acceptance.py -n auto -q
```

### Step 3: implement generator

Internal orchestration boundary:

```python
def generate_exact_adaptive(
    pocket: Polygon,
    *,
    holes: Sequence[Polygon],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    engagement_cap: EngagementCap,
    entry: PreclearedEntry,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> GenerationResult:
    ...
```

No default hides a manufacturing or resolution decision. The function returns
only after traversal, coverage, and replay gates pass. Task 14 assembles the
public `CertifiedToolpath` after final provenance and schema binding; Task 13
does not reference a type defined by a later task.

### Step 4: GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_traversal.py \
  tests/adaptive/test_acceptance.py -n auto --testmon -q
```

Commit: `feat(adaptive): generate certified covered path`

---

## Task 14 — Assemble the artifact, provenance, schema, and public factory

**Files**

- Create: `src/compas_cgal/adaptive/certificate.py`
- Modify: `src/compas_cgal/adaptive/identity.py`
- Create: `src/compas_cgal/adaptive/schema.py`
- Create: `schemas/certified_toolpath-v1.schema.json`
- Create: `tests/adaptive/test_certificate.py`
- Modify: `tests/adaptive/test_identity.py`
- Create: `tests/adaptive/test_schema.py`

### Step 1: write RED artifact/provenance contracts

Require:

- `InputIdentity`, `BuildProvenance`, and `ArtifactIdentity` DAG;
- source revision plus dirty digest, native binary, `pixi.lock`, compiler,
  CGAL, kernel, verified archive hashes, actual native source-tree digests, and
  component source identities;
- exact operation, event witness, depletion witness, traversal, and coverage
  streams;
- same-process and fresh-process stability;
- one-bit mutation of every field changes or invalidates identity;
- any operation reorder changes `ArtifactIdentity`, including commuting
  operations whose fresh physical replay would still succeed;
- operation/witness bijection and state-lineage order;
- artifact never hashes its own digest.

The public factory accepts every Task 13 policy explicitly, calls
`generate_exact_adaptive`, constructs `BuildProvenance`, assembles
`CertifiedToolpath`, performs fresh replay, and returns only the bound artifact:

```python
def exact_certified_adaptive_trochoidal_toolpath(
    pocket: Polygon,
    *,
    holes: Sequence[Polygon],
    cut_plane: CutPlane,
    tool_radius: ToolRadius,
    engagement_cap: EngagementCap,
    entry: PreclearedEntry,
    candidate_policy: CandidatePolicy,
    neck_policy: NeckPolicy,
    depletion_policy: DepletionPolicy,
    traversal_policy: TraversalPolicy,
    cut_direction_policy: CutDirectionPolicy,
) -> CertifiedToolpath: ...
```

### Step 2: write RED consumer contract

Serialize, call `json.loads`, validate with `jsonschema`, deserialize, recompute
identity, derive the legacy `ToolpathResult`, and replay it. Reject unknown
schema versions and any canonical-motion/COMPAS-view mismatch.

Run RED:

```bash
pixi run pytest tests/adaptive/test_certificate.py \
  tests/adaptive/test_identity.py tests/adaptive/test_schema.py -n auto -q
```

### Step 3: implement and GREEN

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run pytest tests/adaptive/test_certificate.py \
  tests/adaptive/test_identity.py tests/adaptive/test_schema.py \
  -n auto --testmon -q
pixi run schema
```

Commit: `feat(cert): bind exact toolpath artifact`

---

## Task 15 — Execute ablation, mutation, and Held provenance gates

**Files**

- Create: `scripts/run-adaptive-mutations.py`
- Modify: `tests/adaptive/test_acceptance.py`
- Create: `tests/adaptive/fixtures/held_fig5.json`
- Create: `docs/benchmarks/exact-certified-adaptive-phase1.md`

### Step 1: fixed-policy primitive-grammar ablation

On the same tractable pocket, execute a committed fixed-policy comparator that
uses the same Phase-1 exact-segment/full-circle operation grammar and replay
every lateral operation through the event-exact oracle. Do not feed the legacy
generator's unsupported partial repositioning/lead arcs to this oracle or label
the comparator as a replay of that generator.
Require:

- fixed comparator has a nonzero proved cap-violation tail;
- adaptive path has zero violations;
- both use the same pocket/tool/cap and declared entry state;
- both use the same reachable-domain construction and primitive grammar;
- reported max angles are diagnostics only.

### Step 2: native/Python mutation campaign

The Pixi mutation task creates isolated temporary source copies and compiled
variants; it never edits the working tree. Each mutation runs its named kill
test and must fail:

1. double segment interpolation;
2. double circle rotation;
3. off-guide exact center;
4. inflated removal disk;
5. omitted circle seam;
6. rounded full-sweep coverage;
7. omitted tangency event;
8. omitted overlap event;
9. omitted run merge;
10. omitted contour vertex;
11. omitted cap root;
12. certify after deplete;
13. non-fresh replay;
14. stale witness;
15. empty/seed-only/dropped-branch path;
16. relabelled lateral cut;
17. broken exact phase/Z continuity;
18. stale neck class/effective-cap decision;
19. nonterminal traversal decision;
20. omitted constant-clearance plateau;
21. mutated `BoundaryVertexIdV1` intersection ordinal.

The script uses `pathlib`, `tempfile`, and subprocesses invoked through the
current Pixi environment. No `os.path`, force flag, or worktree mutation.

### Step 3: Held fixture provenance

Commit the vectorized Fig-5 input with extraction/source notes and all
parameters. Run it only if Phase-1 resource bounds allow exact completion.
Otherwise the benchmark report records the measured blocker and leaves the
Held parity label unset. It may not substitute a hand-drawn lookalike.

### Step 4: GREEN

```bash
pixi run ruff format scripts/run-adaptive-mutations.py \
  tests/adaptive/test_acceptance.py
pixi run lint
pixi run pytest tests/adaptive/test_acceptance.py -n auto --testmon -q
pixi run mutations-adaptive
```

Commit: `test(adaptive): prove certificate mechanism`

---

## Task 16 — Full verification and whole-branch review

### Step 1: formatting, lint, typing, schema

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run schema
git diff --check
```

### Step 2: affected and full suites

```bash
pixi run affected
pixi run baseline
```

Do not claim success from `--testmon` alone. Record full-suite test count,
duration, and exit status.

### Step 3: build/package contract

```bash
pixi run wheel
```

Install/test the built wheel in a fresh Pixi environment defined in the
workspace, then validate native imports, JSON round trip, fresh replay, and the
tractable acceptance fixture.

### Step 4: adversarial review

Review the complete branch against the governing spec:

- exact/proposal/reporting boundaries;
- exact `D`, `C_r`, `M_r`, and unreachable-residual boundaries;
- all continuum event classes;
- event-partition root/cell completeness certificates;
- under-cover and full-sweep set proofs;
- certify-before-deplete ordering;
- non-vacuity and coverage;
- frame/unit typing and native stub closure;
- file responsibilities and size;
- named errors/no fallback;
- legacy isolation;
- archive hashes and actual native source-tree digests;
- exact clipping of every MAT dual primitive and neck-state transitions;
- identity/replay independence;
- every mutation kill.

Any finding reopens its owning task and reruns that task’s RED/GREEN and full
gates.

### Step 5: final commit

If review requires no correction, no empty “verification” commit is created.
Otherwise commit each correction by responsibility, rerun Steps 1–4, and report
the exact final SHA and clean/dirty status.

## Completion definition

Phase 1 is complete only when:

- the algebraic-traits gates are satisfied;
- exact segment/full-circle event oracles resolve the committed reachable
  corpus;
- the tractable path is non-vacuous, gouge-free, cap-certified, and exactly
  coverage-certified over `M_r`, with `D \ M_r` explicitly bound;
- fresh replay agrees without trusting generator witnesses;
- fixed-policy primitive-grammar ablation and all named mutations demonstrate
  mechanism;
- strict typing, Ruff, schema, full pytest, and wheel consumer gates pass;
- existing generators and compatibility APIs remain green and unchanged in
  behavior;
- the worktree contains no unexplained changes.
