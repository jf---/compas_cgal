# Task 3 Single-Arrangement Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the rejected Task 3 spike with a one-pass exact reachable-domain and incremental coverage implementation, complete the end-to-end pipeline, then use structural counts and bounded timing to decide any deeper indexing refactor.

**Architecture:** One `Arr_curve_data_traits_2` arrangement carries primitive and source-piece provenance through exact splits and overlaps. Exact face parity selects every `C_r` cell; the same DCEL yields `C_r` geometry and certificate records. Immutable exact regions share storage, `M_r` uses one batched range union, and coverage advances by the exact induction `R_next = R_previous \ W`.

**Tech Stack:** C++20, CGAL 6.0.1 `Exact_predicates_exact_constructions_kernel_with_sqrt`, `Gps_circle_segment_traits_2`, `Arr_curve_data_traits_2`, nanobind, Python 3.12, Pixi, pytest-xdist, pytest-testmon, Ruff, strict mypy.

## Global Constraints

- Governing design: `docs/superpowers/specs/2026-07-25-task3-single-arrangement-performance-design.md`.
- Governing master plan: `docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`, Task 3.
- Work only in `/private/tmp/compas_cgal_prs-exact-certified-adaptive-phase1` on `codex/exact-certified-adaptive-phase1`; never mutate `main`, `master`, or the protected checkout.
- Accepted Tasks 0 through 2 end at `0b0b948`; do not modify their behavior.
- The design commit is `3f92f0f`.
- The pre-existing resume file is untracked state; never add, edit, or delete it.
- The current Task 3 production files and tests are an uncommitted rejected spike. Preserve the strong RED fixture intent, but do not retain replay, equality-against-self, global source rematching, sequential global dilation joins, or deep-clone accessors.
- Use only the locked CGAL 6.0.1, Boost 1.82.0, and Eigen 3.4.0 source trees. Keep `CGAL_USE_CORE=1`, `CGAL_HEADER_ONLY`, `CGAL_DISABLE_GMP`, and `CGAL_USE_BOOST_MP`.
- Exact decisions use CGAL comparisons and regularized-set operations. No epsilon, `to_double`, sampled offset, rounded disk chain, fallback kernel, conditional import, `HAS_*`, skip, skipif, or xfail.
- Native structural records contain no handles, pointer values, iteration order, rounded coordinates, or raw algebraic-number serialization.
- Python owns the sole CCAN SHA-256 content identity.
- One reachable owner constructs one provenance arrangement and one each of `D`, `C_r`, `M_r`, and `D \ M_r`.
- Outside arrangement construction, every DCEL feature is visited `O(1)` times; no post-arrangement predicate scans all sources. Pre-sized ephemeral hashes currently provide expected `O(k + z)` dense compilation and never determine serialized identity.
- Canonical ring and cycle rotation is linear; ordinary record collections are sorted.
- Read-only exact-region accessors and ledger clones share immutable geometry.
- One coverage transition constructs one sweep, one accumulated union, and one residual difference.
- Independent full replay remains Task 13 work and never runs in a Task 3 constructor or transition.
- Keep every C++ and Python file below 1,000 lines. One responsibility per file.
- Python 3.12, idiomatic Python, Google docstrings, no `__all__`, strict mypy.
- Every pytest command includes `-n auto`. After changes, the focused GREEN command includes `--testmon` with a fresh `TESTMON_DATAFILE`.
- Native algorithm gates and pytest commands run through Pixi.
- Apply a hard 180-second external process-group watchdog to each isolated diagnostic and to the complete focused Task 3 pair. A timeout stops broader testing.
- Do not profile unless Tasks 3 through 5 are complete, all structural operation-count assertions pass, and one bounded end-to-end fixture remains slow.
- Deterministic collision-independent indexing is deferred to the final post-Task5 measurement/refactor-decision gate; it is not a precondition for completing Tasks 3 through 5.
- Run Ruff before every Python commit.
- Author and committer for every commit: `Jelle Feringa <jelleferinga@gmail.com>`.

## File Map

- Create `src/exact_region_2.h/.cpp`: immutable shared exact-region storage and basic set predicates.
- Create `src/exact_sweep_2.h/.cpp`: exact disk, capsule, arc-sweep, full-circle, and batched-union primitives.
- Create `src/exact_build_audit_2.h`: native structural operation-count contracts.
- Create `src/reachable_input_2.h/.cpp`: canonical input ownership and validation.
- Create `src/reachable_arrangement_2.h/.cpp`: primitive labels, provenance arrangement, face parity, and selected topology.
- Create `src/canonical_rotation_2.h`: shared linear minimal-rotation template.
- Create `src/reachable_certificate_2.h/.cpp`: stable source, vertex, cell, component, and recipe records.
- Modify `src/reachable_domain_2.h/.cpp`: public one-pass owner orchestration only.
- Modify `src/coverage_2.h/.cpp`: incremental exact coverage state only.
- Create `src/coverage_bindings_2.cpp`: nanobind module and exception mapping only.
- Modify `src/compas_cgal/_coverage_2.pyi`: exact native public contract.
- Modify `src/compas_cgal/adaptive/reachable_domain.py`: typed Python owner and structural certificate.
- Modify `src/compas_cgal/adaptive/coverage.py`: typed incremental ledger and sweep witnesses.
- Modify `src/compas_cgal/adaptive/errors.py`: named Python boundary errors.
- Modify `tests/adaptive/test_reachable_domain.py`: exact geometry, structural certificate, invariance, and ownership contracts.
- Modify `tests/adaptive/test_coverage.py`: exact sweep, induction, lineage, transaction, and ownership contracts.
- Create `tests/native/task3_algorithm_gate.cpp`: structural work-count and native exact-geometry gate.
- Modify `CMakeLists.txt`: reusable exact-coverage core, nanobind module, and excluded native gate target.
- Modify `pyproject.toml`: editable coverage rebuild and native algorithm-gate Pixi tasks.

---

### Task 1: Immutable exact regions, sweep parts, and native gate

**Files:**
- Create: `src/exact_region_2.h`
- Rewrite: `src/exact_region_2.cpp`
- Create: `src/exact_sweep_2.h`
- Create: `src/exact_sweep_2.cpp`
- Create: `src/exact_build_audit_2.h`
- Create: `tests/native/task3_algorithm_gate.cpp`
- Modify: `src/reachable_domain_2.h`
- Modify: `CMakeLists.txt`
- Modify: `pyproject.toml`

**Interfaces:**
- Consumes: locked `ReachKernel`, `ReachTraits`, `ReachSet`, `ReachPolygon`, and exact curve aliases moved from `reachable_domain_2.h`.
- Produces:

```cpp
enum class ExactRegionRole2 {
    Design,
    CenterDomain,
    ReachableMaterial,
    UnreachableResidual,
    AccumulatedSweeps,
    CoverageResidual,
};

class ExactRegion2 {
public:
    static ExactRegion2 build(
        ReachSet set,
        ExactRegionRole2 role,
        std::string recipe_record);
    ExactRegion2 clone() const;
    bool contains(double x, double y) const;
    bool is_empty() const;
    std::size_t component_count() const;
    bool is_subset_of(const ExactRegion2& other) const;
    bool exactly_equals(const ExactRegion2& other) const;
    bool shares_storage_with_for_audit(const ExactRegion2& other) const;
    const ReachSet& set() const;
    ExactRegionRole2 role() const;
    const std::string& recipe_record() const;

private:
    ExactRegion2(
        std::shared_ptr<const ReachSet> set,
        ExactRegionRole2 role,
        std::string recipe_record);
    std::shared_ptr<const ReachSet> set_;
    ExactRegionRole2 role_;
    std::string recipe_record_;
};

ReachPolygon reach_disk_polygon(
    const ReachKernelPoint& center,
    const ReachFT& radius);
std::vector<ReachPolygon> reach_capsule_parts(
    const ReachKernelPoint& start,
    const ReachKernelPoint& end,
    const ReachFT& radius);
std::vector<ReachPolygon> reach_arc_sweep_parts(
    const ReachXCurve& guide_arc,
    const ReachFT& tool_radius);
ReachSet reach_join_parts(
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes);
ReachSet reach_full_circle_sweep(
    const ReachKernelPoint& center,
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius);
std::string reach_tagged_record(
    std::string_view tag,
    const std::vector<std::string>& fields);
std::string reach_binary64_record(double value);
std::string reach_u64_record(std::size_t value);
```

`ReachableDomainBuildAudit2` and `CoverageTransitionAudit2` are native audit
types, not Python domain objects:

```cpp
struct ReachableDomainBuildAudit2 {
    std::size_t geometry_passes = 0;
    std::size_t provenance_arrangements = 0;
    std::size_t center_extractions = 0;
    std::size_t material_batch_unions = 0;
    std::size_t subset_decisions = 0;
    std::size_t residual_differences = 0;
    std::size_t source_geometric_rematches = 0;
    std::size_t input_vertex_count = 0;
    std::size_t ring_rotation_comparisons = 0;
    std::size_t cycle_element_count = 0;
    std::size_t cycle_rotation_comparisons = 0;
    std::size_t selected_face_count = 0;
    std::size_t selected_adjacency_count = 0;
};

struct CoverageTransitionAudit2 {
    std::size_t sweep_constructions = 0;
    std::size_t accumulated_unions = 0;
    std::size_t residual_differences = 0;
    std::size_t replay_constructions = 0;
    std::size_t immutable_region_deep_copies = 0;
};
```

- [ ] **Step 1: Write the native RED ownership and sweep gate**

Create `tests/native/task3_algorithm_gate.cpp` with ordinary fail-loud checks:

```cpp
#include "exact_region_2.h"
#include "exact_sweep_2.h"

#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void exact_region_storage_gate()
{
    ReachSet disk;
    disk.insert(reach_disk_polygon(ReachKernelPoint(0, 0), ReachFT(2)));
    const ExactRegion2 original = ExactRegion2::build(
        std::move(disk),
        ExactRegionRole2::ReachableMaterial,
        "native-gate-region");
    const ExactRegion2 clone = original.clone();
    require(
        original.shares_storage_with_for_audit(clone),
        "read-only clone deep-copied exact storage");
    require(clone.contains(2.0, 0.0), "closed disk lost exact boundary");
}

void exact_sweep_gate()
{
    const std::vector<ReachPolygon> capsule = reach_capsule_parts(
        ReachKernelPoint(0, 0),
        ReachKernelPoint(3, 4),
        ReachFT(5));
    const ReachSet capsule_set = reach_join_parts(capsule, {});
    require(
        capsule_set.oriented_side(ReachPoint(ReachFT(-2.5), ReachFT(5)))
            != CGAL::ON_NEGATIVE_SIDE,
        "sqrt(13) capsule lost side boundary");
}

} // namespace

int main()
{
    exact_region_storage_gate();
    exact_sweep_gate();
}
```

- [ ] **Step 2: Wire the native target and verify RED**

Refactor CMake so exact geometry is reusable:

```cmake
add_library(coverage_exact_core STATIC
    src/exact_region_2.cpp
    src/exact_sweep_2.cpp
)
set_target_properties(coverage_exact_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(coverage_exact_core PRIVATE CGAL_USE_CORE=1)
target_include_directories(coverage_exact_core PRIVATE
    ${CGAL_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
)

add_executable(task3_algorithm_gate EXCLUDE_FROM_ALL
    tests/native/task3_algorithm_gate.cpp
)
target_link_libraries(task3_algorithm_gate PRIVATE coverage_exact_core)
target_compile_definitions(task3_algorithm_gate PRIVATE CGAL_USE_CORE=1)
```

Add forward Pixi tasks:

```toml
_task3-gate-configure = "cmake -S . -B build/task3-algorithm -G Ninja -Dnanobind_DIR=$(python -m nanobind --cmake_dir)"
_task3-gate-build = "cmake --build build/task3-algorithm --target task3_algorithm_gate"
_task3-gate-run = "build/task3-algorithm/task3_algorithm_gate"
task3-algorithm-gate = { depends-on = ["_task3-gate-configure", "_task3-gate-build", "_task3-gate-run"], description = "Prove bounded exact Task 3 construction" }
```

Run through Pixi with the 180-second external watchdog:

```bash
pixi run task3-algorithm-gate
```

Expected: compile failure for missing split headers/interfaces.

- [ ] **Step 3: Implement immutable region storage**

Move all exact aliases into `exact_region_2.h`. Implement `ExactRegion2::build`
with `std::make_shared<const ReachSet>(std::move(set))`; keep its raw constructor
private. `clone()` returns `*this`, and
`shares_storage_with_for_audit()` compares `set_.get()` addresses without ever
serializing them.

```cpp
ExactRegion2 ExactRegion2::build(
    ReachSet set,
    ExactRegionRole2 role,
    std::string recipe_record)
{
    return ExactRegion2(
        std::make_shared<const ReachSet>(std::move(set)),
        role,
        std::move(recipe_record));
}

ExactRegion2 ExactRegion2::clone() const
{
    return *this;
}
```

Keep `reach_exact_subset`, `reach_exact_equal`, component enumeration, and
binary-record helpers in `exact_region_2.cpp`. No mutating method is public.

- [ ] **Step 4: Implement exact sweep parts and batched union**

Move disk/capsule/arc/full-circle geometry into `exact_sweep_2.cpp`.
`reach_capsule_parts` and `reach_arc_sweep_parts` return polygon primitives
without constructing a `ReachSet` per boundary curve.

Implement the two-range union:

```cpp
ReachSet reach_join_parts(
    const std::vector<ReachPolygon>& polygons,
    const std::vector<ReachPolygonWithHoles>& polygons_with_holes)
{
    ReachSet result;
    result.join(
        polygons.begin(),
        polygons.end(),
        polygons_with_holes.begin(),
        polygons_with_holes.end());
    return result;
}
```

Preserve exact `guide_radius <`, `==`, and `>` tool-radius branches.

- [ ] **Step 5: Run the native gate GREEN**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: exit 0 within the 180-second watchdog.

- [ ] **Step 6: Commit**

Run `git diff --check`, then commit only Task 1 files:

```bash
git add CMakeLists.txt pyproject.toml src/exact_build_audit_2.h \
  src/exact_region_2.h src/exact_region_2.cpp src/exact_sweep_2.h \
  src/exact_sweep_2.cpp src/reachable_domain_2.h \
  tests/native/task3_algorithm_gate.cpp
git commit -m "feat(coverage): add exact core"
```

---

### Task 2: Provenance arrangement and exact face parity

**Files:**
- Create: `src/canonical_rotation_2.h`
- Create: `src/reachable_input_2.h`
- Create: `src/reachable_input_2.cpp`
- Create: `src/reachable_arrangement_2.h`
- Create: `src/reachable_arrangement_2.cpp`
- Create: `src/reachable_errors_2.h`
- Modify: `src/exact_build_audit_2.h`
- Modify: `tests/native/task3_algorithm_gate.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: exact aliases and sweep polygons from Task 1.
- Produces:

```cpp
class InvalidReachableDomainInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class ReachableArrangementTopologyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class PocketNotMachinableError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class ReachableMaterialContainmentError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct ReachCurveLabels2 {
    std::vector<std::string> source_piece_ids;
    std::vector<std::string> primitive_ids;
    bool operator==(const ReachCurveLabels2&) const = default;
};

struct MergeReachCurveLabels2 {
    ReachCurveLabels2 operator()(
        const ReachCurveLabels2& left,
        const ReachCurveLabels2& right) const;
};

using ReachDataTraits2 = CGAL::Arr_curve_data_traits_2<
    ReachTraits,
    ReachCurveLabels2,
    MergeReachCurveLabels2>;

struct ReachFaceState2 {
    bool classified = false;
    bool selected = false;
    bool outer_active = false;
    std::size_t active_holes = 0;
    std::size_t active_forbidden = 0;
};

using ReachFaceDcel2 =
    CGAL::Arr_face_extended_dcel<ReachDataTraits2, ReachFaceState2>;
using ReachArrangement2 =
    CGAL::Arrangement_2<ReachDataTraits2, ReachFaceDcel2>;

enum class ReachPrimitiveKind2 {
    Outer,
    Hole,
    Forbidden,
};
using ReachPrimitiveKinds2 =
    std::map<std::string, ReachPrimitiveKind2>;

struct CanonicalReachRing2 {
    bool outer;
    std::size_t canonical_ordinal;
    std::vector<ReachKernelPoint> points;
    std::vector<std::array<double, 2>> binary64_points;
    std::string record;
};

struct CanonicalReachInput2 {
    CanonicalReachRing2 outer;
    std::vector<CanonicalReachRing2> holes;
    ReachFT radius;
    double binary64_radius;
    std::string recipe_record;
};

struct ReachableArrangement2 {
    ReachArrangement2 arrangement;
    CanonicalReachInput2 input;
    std::vector<std::string> source_records;
    ReachPolygonWithHoles design_polygon;
    ReachPolygonWithHoles center_polygon;
    ReachableDomainBuildAudit2 audit;
};

CanonicalReachInput2 canonical_reach_input(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius);
ReachableArrangement2 build_reachable_arrangement(
    CanonicalReachInput2 input);
```

`reachable_arrangement_2.cpp` exposes exactly:

```cpp
void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds);
```

- [ ] **Step 1: Extend the native gate with RED label propagation**

Add a rectangle, acute corner, two-overlapping-curve fixture, and one island.
Assert:

```cpp
compas::RowMatrixXd rectangle_matrix()
{
    compas::RowMatrixXd boundary(4, 3);
    boundary << 0, 0, 0, 10, 0, 0, 10, 8, 0, 0, 8, 0;
    return boundary;
}

CanonicalReachInput2 rectangle_input()
{
    return canonical_reach_input(rectangle_matrix(), {}, 1.0);
}

ReachableArrangement2 rectangle =
    build_reachable_arrangement(rectangle_input());
require(
    rectangle.audit.provenance_arrangements == 1,
    "reachable build constructed multiple arrangements");
require(
    rectangle.audit.source_geometric_rematches == 0,
    "reachable build rematched sources geometrically");
require(
    rectangle.source_records.size() == 12,
    "rectangle did not bind three sources per vertex");
```

For every halfedge separating a selected face from an unselected face, require
nonempty propagated `source_piece_ids`. For an exact overlapping source curve,
require sorted unique merged labels.

- [ ] **Step 2: Run RED**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: compile failure for missing arrangement interfaces.

- [ ] **Step 3: Implement linear canonical input**

Implement the shared Booth-style template in `canonical_rotation_2.h` and call
it with exact-point comparison. Track comparison counts in the audit and require
a linear bound in the native gate:

```cpp
template <typename Value, typename Compare>
std::size_t minimal_rotation_index(
    const std::vector<Value>& values,
    Compare compare,
    std::size_t& comparison_count)
{
    const std::size_t size = values.size();
    std::size_t left = 0;
    std::size_t right = 1;
    std::size_t offset = 0;
    while (left < size && right < size && offset < size) {
        ++comparison_count;
        const CGAL::Comparison_result order = compare(
            values[(left + offset) % size],
            values[(right + offset) % size]);
        if (order == CGAL::EQUAL) {
            ++offset;
            continue;
        }
        if (order == CGAL::LARGER) {
            left += offset + 1;
            if (left == right) {
                ++left;
            }
        }
        else {
            right += offset + 1;
            if (left == right) {
                ++right;
            }
        }
        offset = 0;
    }
    return std::min(left, right);
}

require(
    rectangle.audit.ring_rotation_comparisons
        <= 3 * rectangle.audit.input_vertex_count,
    "ring canonicalization exceeded linear comparison bound");
```

Construct one exact `ReachPolygonWithHoles` and validate it with
`CGAL::is_valid_polygon_with_holes<ReachTraits>`. Do not perform per-hole
subset, overlap, union, or design-difference operations.

- [ ] **Step 4: Generate labelled primitive boundaries**

For each edge, insert one strip boundary:

- both offset sides carry their stable source-piece ID and the strip primitive
  ID;
- both construction-only end segments carry only the strip primitive ID.

For each vertex, insert its two x-monotone disk-circle halves once. They carry
the vertex-circle source-piece ID and disk primitive ID.

Outer and hole boundaries carry their design primitive IDs only.

Sort and uniquify both label vectors in `MergeReachCurveLabels2`.

```cpp
template <typename T>
std::vector<T> sorted_unique(std::vector<T> values)
{
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

ReachDataTraits2::X_monotone_curve_2 labelled_curve(
    const ReachXCurve& curve,
    std::vector<std::string> source_piece_ids,
    std::vector<std::string> primitive_ids)
{
    return {
        curve,
        ReachCurveLabels2{
            sorted_unique(std::move(source_piece_ids)),
            sorted_unique(std::move(primitive_ids)),
        },
    };
}
```

- [ ] **Step 5: Implement one aggregate arrangement and parity traversal**

Add `src/reachable_arrangement_2.cpp` to `coverage_exact_core`.
Call `CGAL::insert` once over all labelled x-monotone curves and increment
`provenance_arrangements` once.

```cpp
ReachArrangement2 arrangement;
CGAL::insert(arrangement, labelled.begin(), labelled.end());
++audit.provenance_arrangements;
classify_faces_by_primitive_parity(arrangement, primitive_kinds);
```

Depth-first traverse the dual graph from the unbounded face. Toggle only the
primitive IDs attached to the crossed halfedge, maintain one mutable active-ID
state, and store the resulting aggregate tuple in each face. A second route to
an already classified face must reproduce its tuple.

Mark a face selected exactly when outer is active, holes are zero, and
forbidden primitives are zero. Count connected selected-face components through
shared halfedges. Raise `PocketNotMachinableError` unless the count is one.

- [ ] **Step 6: Extract the regularized selected boundary once**

Traverse halfedges with selected left face and unselected right face. Follow
them into canonical outer/hole cycles, strip the data wrapper back to base
`ReachXCurve`, and reject any such edge without source-piece provenance.
Increment `center_extractions` once.

```cpp
const bool is_selected_boundary =
    halfedge->face()->data().selected
    && !halfedge->twin()->face()->data().selected;
if (is_selected_boundary && halfedge->curve().data().source_piece_ids.empty()) {
    throw ReachableArrangementTopologyError(
        "selected boundary edge has no propagated source provenance");
}
```

- [ ] **Step 7: Run GREEN**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: exit 0; one arrangement, one center extraction, zero rematches.

- [ ] **Step 8: Commit**

```bash
git add CMakeLists.txt src/exact_build_audit_2.h \
  src/canonical_rotation_2.h \
  src/reachable_input_2.h src/reachable_input_2.cpp \
  src/reachable_arrangement_2.h src/reachable_arrangement_2.cpp \
  src/reachable_errors_2.h tests/native/task3_algorithm_gate.cpp
git commit -m "feat(coverage): build exact center cells"
```

---

### Task 3: Structural source, vertex, cell, and component certificate

**Files:**
- Create: `src/reachable_certificate_2.h`
- Create: `src/reachable_certificate_2.cpp`
- Modify: `src/reachable_arrangement_2.h`
- Modify: `src/reachable_arrangement_2.cpp`
- Modify: `tests/native/task3_algorithm_gate.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: the one classified `ReachArrangement2` and its propagated labels.
- Produces:

```cpp
struct ReachableDomainCertificate2 {
    std::string strategy_version;
    std::vector<std::string> source_curve_records;
    std::vector<std::string> selected_cell_records;
    std::vector<std::string> component_records;
    bool exact_cell_selection;
    bool complete_source_provenance;
    bool reachable_subset_of_design;
    std::string input_recipe_record;

    bool matches_exact_inputs(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius) const;
};

ReachableDomainCertificate2 build_reachable_certificate(
    ReachableArrangement2& reachable,
    bool reachable_subset_of_design);
```

`reachable_certificate_2.cpp` also owns these internal helpers:

```cpp
struct ReachPointLess2 {
    bool operator()(const ReachPoint& left, const ReachPoint& right) const;
};
using ReachVertexRecordMap2 =
    std::map<ReachPoint, std::string, ReachPointLess2>;
std::string selected_face_record(
    ReachArrangement2::Face_const_handle face,
    const ReachVertexRecordMap2& vertex_records,
    ReachableDomainBuildAudit2& audit);
```

The certificate has no `exact_reconstruction` or
`residual_matches_difference` field. Its structural booleans describe evidence
actually produced by this build; Task 13 owns independent replay.

- [ ] **Step 1: Write the native RED certificate tests**

Build a certificate directly from the classified native arrangement. Assert
that every selected face has one cell record, every selected-cell adjacency
appears once in the component record, every source is bound, and there is
exactly one component:

```cpp
ReachableDomainCertificate2 certificate =
    build_reachable_certificate(rectangle, true);
require(
    certificate.selected_cell_records.size()
        == rectangle.audit.selected_face_count,
    "certificate omitted selected arrangement cells");
require(
    certificate.source_curve_records.size() == 12,
    "certificate omitted rectangle source curves");
require(
    certificate.component_records.size() == 1,
    "certificate did not bind the selected component");
require(
    certificate.exact_cell_selection,
    "certificate lost exact face selection");
require(
    certificate.complete_source_provenance,
    "certificate lost propagated source provenance");
```

- [ ] **Step 2: Run RED**

Run under the 180-second watchdog:

```bash
pixi run task3-algorithm-gate
```

Expected: compile failure for the missing certificate builder.

- [ ] **Step 3: Implement stable vertex identities in one pass**

Add `src/reachable_certificate_2.cpp` to `coverage_exact_core`.
For each arrangement vertex, build the sorted incident source-piece multiset.
Group vertices by that multiset, sort each group by exact `compare_xy`, and
write a direct vertex-to-ordinal map. No later lookup may linearly search the
group.

```cpp
for (auto& [incidence, points] : points_by_incidence) {
    std::sort(points.begin(), points.end(), ReachPointLess{});
    for (std::size_t ordinal = 0; ordinal < points.size(); ++ordinal) {
        vertex_record_by_point.emplace(
            points[ordinal],
            reach_tagged_record(
                "arrangement-vertex-v1",
                {incidence, reach_u64_record(ordinal)}));
    }
}
```

- [ ] **Step 4: Implement linear cycle records**

Use Booth-style minimal rotation over already-built element records:

```cpp
std::vector<std::string> canonical_cycle_rotation(
    const std::vector<std::string>& elements);
```

Record and assert:

```cpp
require(
    audit.cycle_rotation_comparisons
        <= 3 * audit.cycle_element_count,
    "cycle canonicalization exceeded linear comparison bound");
```

- [ ] **Step 5: Emit one record per selected face**

Each selected face record binds all outer/inner CCBs, stable vertices,
propagated source-piece IDs, exact direction, and endpoint role. A selected face
without a source-labelled boundary fails with
`ReachableArrangementTopologyError`.

```cpp
for (auto face = arrangement.faces_begin();
     face != arrangement.faces_end();
     ++face) {
    if (face->data().selected) {
        selected_cell_records.push_back(
            selected_face_record(face, vertex_record_by_point, audit));
    }
}
std::sort(selected_cell_records.begin(), selected_cell_records.end());
```

Emit the sole component record from sorted selected-cell IDs, sorted adjacency
pairs, and sorted regularized boundary-cycle IDs.

- [ ] **Step 6: Prevent quadratic source records**

Encode ring role/ordinal, feature ordinal, construction role, and radius record.
Do not embed the full ring record. Decode each native record in the gate and
assert its fields contain canonical ring identity rather than `input-ring-v1`.

```cpp
return reach_tagged_record(
    "source-curve-v2",
    {
        ring.outer ? "outer" : "hole",
        reach_u64_record(ring.canonical_ordinal),
        reach_u64_record(feature_ordinal),
        std::string(role),
        reach_binary64_record(binary64_radius),
    });
```

- [ ] **Step 7: Run native GREEN**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: all native certificate tests pass within the 180-second watchdog.

- [ ] **Step 8: Commit**

Commit:

```bash
git add CMakeLists.txt src/reachable_arrangement_2.h \
  src/reachable_arrangement_2.cpp src/reachable_certificate_2.h \
  src/reachable_certificate_2.cpp tests/native/task3_algorithm_gate.cpp
git commit -m "feat(coverage): bind exact center cells"
```

---

### Task 4: One-pass reachable-domain owner and Python boundary

**Files:**
- Rewrite: `src/reachable_domain_2.h`
- Rewrite: `src/reachable_domain_2.cpp`
- Modify: `src/coverage_2.cpp`
- Create: `src/coverage_bindings_2.cpp`
- Modify: `src/compas_cgal/_coverage_2.pyi`
- Modify: `src/compas_cgal/adaptive/reachable_domain.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_reachable_domain.py`
- Modify: `tests/native/task3_algorithm_gate.cpp`
- Modify: `CMakeLists.txt`
- Modify: `pyproject.toml`

**Interfaces:**
- Consumes: Tasks 1 through 3.
- Produces:

```cpp
class ReachableDomain2 {
public:
    ReachableDomain2(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);
    ExactRegion2 design_region() const;
    ExactRegion2 center_domain() const;
    ExactRegion2 reachable_material() const;
    ExactRegion2 unreachable_residual() const;
    ReachableDomainCertificate2 certificate() const;
    const ReachableDomainBuildAudit2& build_audit_for_native_gate() const;

private:
    struct State {
        ExactRegion2 design;
        ExactRegion2 center_domain;
        ExactRegion2 reachable_material;
        ExactRegion2 unreachable_residual;
        ReachableDomainCertificate2 certificate;
        ReachableDomainBuildAudit2 audit;
    };
    explicit ReachableDomain2(State state);
    static State build_state(
        Eigen::Ref<const compas::RowMatrixXd> design_boundary,
        const std::vector<compas::RowMatrixXd>& holes,
        double tool_radius);
};
```

`build_audit_for_native_gate()` is not bound by nanobind.
The public constructor delegates through `build_state(...)` to the private state
constructor. No `ExactRegion2` ever has an empty or null placeholder state.

- [ ] **Step 1: Add the RED one-pass audit**

Extend the native gate:

```cpp
const ReachableDomain2 domain(rectangle_matrix(), {}, 1.0);
const ReachableDomainBuildAudit2& audit =
    domain.build_audit_for_native_gate();
require(audit.geometry_passes == 1, "geometry replayed");
require(audit.provenance_arrangements == 1, "arrangement replayed");
require(audit.center_extractions == 1, "C_r enumerated repeatedly");
require(audit.material_batch_unions == 1, "M_r was not one batch union");
require(audit.subset_decisions == 1, "unexpected subset decisions");
require(audit.residual_differences == 1, "residual reconstructed");
require(audit.source_geometric_rematches == 0, "sources rematched");
```

Update Python certificate construction to require
`exact_cell_selection` and `complete_source_provenance` without weakening the
exact fixture assertions:

```python
certificate = domain.certificate
assert certificate.exact_cell_selection
assert certificate.complete_source_provenance
assert certificate.reachable_subset_of_design
assert len(certificate.source_curve_records) == 3 * total_vertex_count
assert certificate.selected_cell_records
assert len(certificate.component_records) == 1
```

Rename the spike's insertion-order test to
`test_structural_identity_is_invariant_to_ring_and_hole_insertion_order`; it no
longer claims constructor replay.

- [ ] **Step 2: Run RED**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: one-pass audit failure while the spike constructor still replays.

- [ ] **Step 3: Implement the one-pass owner**

The constructor performs exactly:

```cpp
CanonicalReachInput2 input =
    canonical_reach_input(design_boundary, holes, tool_radius);
ReachableArrangement2 selected =
    build_reachable_arrangement(std::move(input));
ReachSet design(selected.design_polygon);
ReachSet center(selected.center_polygon);
ReachSet material = build_reachable_material_once(
    selected.center_polygon,
    selected.input.radius,
    selected.audit);
const bool subset = reach_exact_subset(material, design);
if (!subset) {
    throw ReachableMaterialContainmentError(
        "exact reachable material is not contained in the design");
}
ReachSet residual(design);
residual.difference(material);
ReachableDomainCertificate2 certificate =
    build_reachable_certificate(selected, subset);
```

Move each set once into immutable `ExactRegion2`. Remove the rejected
`construct_geometry` replay, duplicate expected residual, cross-replay
equalities, post-hoc source matching, and sequential `result.join(part)` loop.

- [ ] **Step 4: Implement one batched `M_r` union**

Use the already-extracted selected boundary. Collect every capsule/arc polygon
part, add the one center polygon with holes, and call `reach_join_parts` once.
Increment `material_batch_unions` once.

The internal helper has one explicit contract:

```cpp
ReachSet build_reachable_material_once(
    const ReachPolygonWithHoles& center,
    const ReachFT& radius,
    ReachableDomainBuildAudit2& audit);
```

Perform one subset difference and one residual difference. The residual recipe
records `D \ M_r`; it is not recomputed for a boolean flag.

- [ ] **Step 5: Bind native certificate and named errors**

First make the still-replaying coverage implementation compile against the new
sweep-part API without changing its Task 5 behavior:

```cpp
const ReachSet sweep = reach_join_parts(
    reach_capsule_parts(start, end, radius),
    {});
```

Move the nanobind module into `coverage_bindings_2.cpp`. Expose read-only
certificate fields and exact-region methods there. Map native input, topology,
machinability, and containment exceptions to distinct Python exceptions in
`adaptive/errors.py`.

Update `_coverage_2.pyi` exactly; do not add overload ambiguity or conditional
imports.

- [ ] **Step 6: Update Python ownership**

`ReachableDomain.build(...)` canonicalizes once, calls one native owner, checks
`matches_exact_inputs`, validates the structural certificate, and stores the
native owner. Region properties return cheap immutable clones.

Keep CCAN digest ownership in Python. Remove all use of the rejected replay and
residual-reconstruction fields from canonical bytes.

```python
@dataclass(frozen=True)
class ReachableDomainCertificate:
    design_boundary: CanonicalRingV1
    holes: tuple[CanonicalRingV1, ...]
    tool_radius: ToolRadius
    strategy_version: bytes
    source_curve_records: tuple[bytes, ...]
    selected_cell_records: tuple[bytes, ...]
    component_records: tuple[bytes, ...]
    exact_cell_selection: bool
    complete_source_provenance: bool
    reachable_subset_of_design: bool
```

- [ ] **Step 7: Make editable rebuild include coverage**

Change:

```toml
_editable-rebuild = '''python -c "from compas_cgal import _coverage_2, _stock_2, _toolpath"'''
```

Build `_coverage_2` from `coverage_bindings_2.cpp`, link the exact core, and add
these files to `coverage_exact_core`:

- `reachable_arrangement_2.cpp`;
- `reachable_certificate_2.cpp`;
- `reachable_domain_2.cpp`;
- `coverage_2.cpp`.

Use:

```cmake
add_nanobind_module(_coverage_2 src/coverage_bindings_2.cpp)
target_link_libraries(_coverage_2 PRIVATE coverage_exact_core)
target_compile_definitions(_coverage_2 PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 8: Run format, type, native, and isolated GREEN**

Run:

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run task3-algorithm-gate
```

Then run each formerly stalled node separately, each with its own 180-second
watchdog:

```bash
pixi run pytest \
  tests/adaptive/test_reachable_domain.py::test_structural_identity_is_invariant_to_ring_and_hole_insertion_order \
  -n auto -q
pixi run pytest \
  tests/adaptive/test_reachable_domain.py::test_circle_vertical_trim_extrema_have_distinct_stable_structural_ordinals \
  -n auto -q
```

Expected: both pass. No broader suite starts if either reaches the watchdog.

- [ ] **Step 9: Run focused reachable-domain GREEN**

Use a fresh testmon database:

```bash
task3_testmon_dir="$(mktemp -d)"
TESTMON_DATAFILE="$task3_testmon_dir/data" \
  pixi run pytest tests/adaptive/test_reachable_domain.py \
  -n auto --testmon -q
```

Expected: all reachable-domain tests pass within 180 seconds.

- [ ] **Step 10: Commit**

```bash
git add CMakeLists.txt pyproject.toml src/reachable_domain_2.h \
  src/reachable_domain_2.cpp src/coverage_2.cpp \
  src/coverage_bindings_2.cpp \
  src/compas_cgal/_coverage_2.pyi \
  src/compas_cgal/adaptive/reachable_domain.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_reachable_domain.py \
  tests/native/task3_algorithm_gate.cpp
git commit -m "feat(coverage): build reachable domain once"
```

---

### Task 5: Incremental coverage, final gates, and Task 3 handoff

**Files:**
- Rewrite: `src/coverage_2.h`
- Rewrite: `src/coverage_2.cpp`
- Modify: `src/coverage_bindings_2.cpp`
- Modify: `src/exact_build_audit_2.h`
- Modify: `src/compas_cgal/_coverage_2.pyi`
- Modify: `src/compas_cgal/adaptive/coverage.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Modify: `tests/adaptive/test_coverage.py`
- Modify: `tests/native/task3_algorithm_gate.cpp`
- Modify: `.superpowers/sdd/2026-07-25-task3-single-arrangement-performance/progress.md`

**Interfaces:**
- Consumes: immutable `ExactRegion2`, exact sweep primitives, and one-pass
  reachable owner.
- Produces:

```cpp
struct CoverageSweepRecord2 {
    std::string strategy_version;
    std::string structural_record;
    bool segment;
    double center_x;
    double center_y;
    double first_x;
    double first_y;
    double tool_radius;

    bool matches_exact_segment(
        double x0,
        double y0,
        double x1,
        double y1,
        double expected_tool_radius) const;
    bool matches_exact_full_circle(
        double cx,
        double cy,
        double phase_x,
        double phase_y,
        double expected_tool_radius) const;
};

class InvalidCoverageGeometryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class CoverageTransitionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class Coverage2 {
public:
    Coverage2(
        const ExactRegion2& reachable_material,
        double precleared_x,
        double precleared_y,
        double precleared_radius);
    Coverage2 clone() const;
    CoverageSweepRecord2 add_segment_sweep(
        double x0,
        double y0,
        double x1,
        double y1,
        double tool_radius);
    CoverageSweepRecord2 add_full_circle_sweep(
        double cx,
        double cy,
        double phase_x,
        double phase_y,
        double tool_radius);
    bool residual_is_empty() const;
    std::size_t residual_component_count() const;
    ExactRegion2 residual() const;
    ExactRegion2 accumulated_sweeps() const;
    std::vector<std::string> residual_component_records() const;
    std::vector<std::string> sweep_records() const;
    bool exact_residual_relation() const;
    const CoverageTransitionAudit2&
        last_transition_audit_for_native_gate() const;
    bool shares_pretransition_storage_with_for_native_gate(
        const Coverage2& other) const;

private:
    struct State {
        ExactRegion2 reachable_material;
        ExactRegion2 accumulated_sweeps;
        ExactRegion2 residual;
        std::vector<std::string> sweep_records;
        bool exact_residual_relation;
    };
    explicit Coverage2(State state);
    static State build_initial_state(
        const ExactRegion2& reachable_material,
        double precleared_x,
        double precleared_y,
        double precleared_radius);
};
```

The audit accessor is not bound to Python. `exact_residual_relation()` returns
the stored exact induction state; it does not reconstruct a set.
The public constructor delegates through `build_initial_state(...)`; exact
regions never pass through null or mutable placeholder storage.

- [ ] **Step 1: Write RED transition-count tests**

Add native checks after one segment and one full-circle transition:

```cpp
const CoverageTransitionAudit2& audit =
    coverage.last_transition_audit_for_native_gate();
require(audit.sweep_constructions == 1, "sweep replayed");
require(audit.accumulated_unions == 1, "accumulated union repeated");
require(audit.residual_differences == 1, "residual repeated");
require(audit.replay_constructions == 0, "coverage replayed");
require(
    audit.immutable_region_deep_copies == 0,
    "coverage deep-copied immutable regions");
```

Strengthen Python clone tests to prove the original residual remains unchanged
after mutating the clone. In the native gate, assert
`coverage.shares_pretransition_storage_with_for_native_gate(clone)` before
mutating the clone.

- [ ] **Step 2: Run RED**

Run:

```bash
pixi run task3-algorithm-gate
```

Expected: transition audit fails against the replaying spike.

- [ ] **Step 3: Implement initial induction state once**

Construct `A_0` as the preclear disk once and `R_0 = M_r \ A_0` once. Set the
stored exact-relation flag only after both exact operations succeed.

```cpp
const ReachFT exact_precleared_radius(precleared_radius);
const std::string accumulated_recipe = reach_tagged_record(
    "coverage-precleared-v2",
    {
        reachable_material.recipe_record(),
        reach_binary64_record(precleared_x),
        reach_binary64_record(precleared_y),
        reach_binary64_record(precleared_radius),
    });
const std::string residual_recipe = reach_tagged_record(
    "coverage-residual-v2",
    {reachable_material.recipe_record(), accumulated_recipe});
ReachSet accumulated;
accumulated.insert(
    reach_disk_polygon(
        ReachKernelPoint(precleared_x, precleared_y),
        exact_precleared_radius));
ReachSet residual(reachable_material.set());
residual.difference(accumulated);
return State{
    reachable_material.clone(),
    ExactRegion2::build(
        std::move(accumulated),
        ExactRegionRole2::AccumulatedSweeps,
        accumulated_recipe),
    ExactRegion2::build(
        std::move(residual),
        ExactRegionRole2::CoverageResidual,
        residual_recipe),
    {},
    true,
};
```

Do not call `exact_residual_relation()` from the constructor.

- [ ] **Step 4: Implement one transition**

For `W_i`:

```cpp
ReachSet next_accumulated(accumulated_sweeps_.set());
next_accumulated.join(sweep);
ReachSet next_residual(residual_.set());
next_residual.difference(sweep);
```

Build `sweep` once. Construct both next regions and next lineage before
replacing any member. Commit all next state only after every operation succeeds.

Remove sweep replay, residual replay, equality-against-self, and
`M_r \ accumulated` recomputation.

- [ ] **Step 5: Make clones share immutable geometry**

`Coverage2::clone()` copies the lightweight owner. Its `ExactRegion2` members
share `shared_ptr<const ReachSet>` storage. A later transition allocates new
accumulated and residual sets only for the trial owner.

```cpp
bool Coverage2::shares_pretransition_storage_with_for_native_gate(
    const Coverage2& other) const
{
    return reachable_material_.shares_storage_with_for_audit(
               other.reachable_material_)
        && accumulated_sweeps_.shares_storage_with_for_audit(
            other.accumulated_sweeps_)
        && residual_.shares_storage_with_for_audit(other.residual_);
}
```

- [ ] **Step 6: Update native and Python witnesses**

Remove the misleading native `exact_reconstruction` sweep field. A sweep record
binds:

- exact binary64 motion inputs;
- exact tool radius;
- native strategy version;
- structural sweep recipe.

Python checks `matches_exact_segment` or `matches_exact_full_circle`, strategy
identity, owned lineage, and atomic native trial state. Map invalid geometry and
transition failures to distinct named Python errors.

```python
if not matches or strategy != trial.strategy_version:
    raise InvalidCoverageSweepError(
        "native coverage sweep does not bind the exact motion inputs"
    )
expected_records = tuple(
    witness.native_structural_record for witness in self._lineage
) + (structural_record,)
if tuple(trial.sweep_records) != expected_records:
    raise InvalidCoverageSweepError(
        "native coverage history diverged from owned lineage"
    )
```

- [ ] **Step 7: Run native and coverage GREEN**

Run:

```bash
pixi run task3-algorithm-gate
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
```

Then use a fresh testmon database and the 180-second watchdog:

```bash
task3_testmon_dir="$(mktemp -d)"
TESTMON_DATAFILE="$task3_testmon_dir/data" \
  pixi run pytest tests/adaptive/test_coverage.py \
  -n auto --testmon -q
```

Expected: all coverage tests pass.

- [ ] **Step 8: Run the complete focused Task 3 gate**

After one editable build, use one fresh database and one 180-second watchdog:

```bash
task3_testmon_dir="$(mktemp -d)"
TESTMON_DATAFILE="$task3_testmon_dir/data" \
  pixi run pytest tests/adaptive/test_reachable_domain.py \
  tests/adaptive/test_coverage.py -n auto --testmon -q
```

Expected: all focused tests pass. Record exact pass count and wall time as
evidence, not as a correctness theorem.

- [ ] **Step 9: Run regression and repository gates**

Run:

```bash
pixi run wheel
pixi run native-lock-check
pixi run baseline
git diff --check
```

Expected:

- wheel builds `_coverage_2`;
- locked native source digests remain unchanged;
- full suite passes with `-n auto`;
- no whitespace errors;
- no `.testmondata*`, temporary profiling output, resume file, or unrelated
  user files are staged.

- [ ] **Step 10: Commit**

```bash
git add src/coverage_2.h src/coverage_2.cpp \
  src/coverage_bindings_2.cpp \
  src/exact_build_audit_2.h src/compas_cgal/_coverage_2.pyi \
  src/compas_cgal/adaptive/coverage.py \
  src/compas_cgal/adaptive/errors.py \
  tests/adaptive/test_coverage.py \
  tests/native/task3_algorithm_gate.cpp
git commit -m "feat(coverage): add bounded exact ledger"
```

- [ ] **Step 11: Record Task 3 evidence**

Append the accepted commit range, native operation counts, focused test count
and wall time, full baseline count and wall time, Ruff result, strict-mypy
result, wheel result, native-lock result, and review verdict to this plan's SDD
ledger. Do not edit the master resume file.

---

### Final post-Task5 measurement/refactor-decision gate

This gate runs only after Tasks 3 through 5 produce the complete native and
Python pipeline. It decides whether the current pre-sized ephemeral hash
boundary needs deterministic collision-independent replacement.

- [ ] **Step 1: Prove structural and exact contracts first**

Under separate 180-second external process-group watchdogs, run the native
algorithm gate and complete focused reachable-domain/coverage pair. Require:

- one provenance arrangement;
- one center-domain extraction;
- zero source geometric rematches;
- the established face, halfedge, label, parity, component, and boundary visit
  bounds;
- exact geometry and certificate assertions;
- one sweep, accumulated union, and residual difference per transition.

Do not time or profile a structurally failing path.

- [ ] **Step 2: Record bounded end-to-end timing**

With structural counts green, run the complete focused pair once in a fresh
process under the standard 180-second watchdog. Record fixture count, pass
count, process wall time, and whether containment was approached. Timing is
decision evidence, not a correctness assertion.

- [ ] **Step 3: Make the refactor decision**

- If the focused end-to-end gate completes comfortably, retain the pre-sized
  nonserialized hashes and proceed.
- If counts pass but the bounded pipeline remains slow or reaches containment,
  write one runtime hypothesis and profile one bounded fixture only.
- Replace hashes with deterministic collision-independent indexing only if the
  profile makes primitive/handle lookup material. Preserve the exact kernel,
  one-arrangement/one-extraction/zero-rematch contracts, public signatures, and
  stable serialized identities.
- Rerun Steps 1 and 2 after any measured refactor and record before/after
  structural counts and end-to-end timing.

Every measurement or profiler process receives two to three minutes of
external process-group containment; use 180 seconds unless a narrower
120-second fixture budget is explicitly recorded.
