# Number-Type Coherence — Stage 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Invert the canonical direction for `SegmentEventSource2` — the large string carrier — so the source stores `compas_cgal::exact::Rational` and *projects* the `exact-binary64-rational-v1` text; route `FullCircleEventSource2`'s two remaining private door copies through `exact::from_binary64` and `exact::to_canonical`; remove `exact::CanonicalRational::canonical_bytes()` after relocating the absolute anchor it holds; and give the outer `segment-event-source-v1` framing the absolute byte anchor it has never had.

**Architecture:** `SegmentEventSource2` stops carrying six `ExactBinary64Rational2` strings and starts carrying six `compas_cgal::exact::Rational` (= `Epeck::FT`). `ExactBinary64Rational2` survives as a *derived attestation view* built by a new `ExactBinary64Rational2::project`, the segment lane's single caller of `exact::to_canonical`. A new `SegmentAttestation2` groups the six projected values plus the frozen outer record and its digest, is built once in the source's constructor, and is the lane's only projection site. The nanobind layer re-points `.x0` … `.cap_chord_ratio`, `.motion_data`, `.canonical_bytes` and `.canonical_digest` at that attestation, so **the Python surface and `_continuous_tea_2.pyi` do not change**. Five consumers — `segment_projection.cpp`, `segment_fibre.cpp`, `segment_strata.cpp`, `segment_oracle.cpp`, `segment_pair_projection.cpp` — stop decoding text and materialise the exact representative at a declared `.exact()` boundary instead, deleting 45 `parse_rational` occurrences. `FullCircleEventSource2` already carries exact values; what it does not do is use the doors, and stage 3 fixes that without moving a byte.

**Tech Stack:** C++20, CGAL 6.0.1 (vendored `external/cgal`), Boost multiprecision (vendored `external/boost`, GMP disabled), CORE (`CGAL_USE_CORE=1`, boost backend), nanobind, CMake + scikit-build-core, pixi, pytest with Hypothesis.

**Spec:** `docs/superpowers/specs/2026-09-08-number-type-coherence-design.md`, stage 3 row, plus the two items the corrected invariant-2 admonition assigns to stage 3.

---

## Measured surface (2026-09-08, this worktree, commit `45c5cc8a`)

### `FullCircleEventSource2` is already inverted — it just does not use the doors

`segment_source.h:104-106` shows it carrying `ExactCircleMotion2 motion_`, `Epeck::FT tool_radius_`, `Epeck::FT cap_chord_ratio_`. Its consumers read exact values already:

| Consumer | Reads |
|---|---|
| `circle_oracle.cpp:516,522,523` | `source.motion()`, `source.tool_radius()`, `source.cap_chord_ratio()` — all exact |
| `circle_oracle.cpp:604,605,779,781,916,917` | `motion_identity_bytes()`, `cap_identity_bytes()` — frozen records, unchanged |

`grep -c "parse_rational" src/continuous_tea_2/circle_oracle.cpp` is **0**. There are **zero** full-circle text crossings to remove. What stage 3 owns there is two door violations inside `segment_source.cpp`:

1. `from_binary64` (`:307-312`) injects doubles inline as `Epeck::FT(center_x)` instead of `exact::from_binary64` — design invariant 1.
2. `from_exact` (`:346-356`) carries a private `rational_text` lambda that re-implements `to_canonical`'s `Fraction_traits::Decompose` + `convert_to<std::string>` — design invariant 2's single-decomposition half.

### `SegmentEventSource2` is the real carrier

Every `.text()` read below is an inbound text crossing that the inversion deletes. Counts are **occurrences** of the literal `parse_rational`, which is what the ratchet test measures.

| File | Role | `parse_rational` **today** | after stage 2 | after stage 3 | segment-source crossings removed |
|---|---|---|---|---|---|
| `src/continuous_tea_2/segment_source.h` | carrier + view | 0 | 0 | 0 | — |
| `src/continuous_tea_2/segment_source.cpp` | carrier + view | 0 | 0 | 0 | — |
| `src/continuous_tea_2/segment_projection.cpp` | consumer | 24 | 24 | **14** | 10 |
| `src/continuous_tea_2/segment_fibre.cpp` | consumer | 26 | 26 | **9** | 17 |
| `src/continuous_tea_2/segment_strata.cpp` | consumer | 12 | 8 | **1** | 7 |
| `src/continuous_tea_2/segment_oracle.cpp` | consumer | 11 | 11 | **1** | 10 |
| `src/continuous_tea_2/segment_pair_projection.cpp` | consumer | 2 | 2 | **1** | 1 |
| **total** | **7 files** | **75** | **71** | **26** | **45** |

Breakdown of the 45, by site:

- `segment_projection.cpp`: `:273,:275` (inside `coordinate_motion`, which decodes its two `ExactBinary64Rational2` parameters), `:485,:487,:489,:492` (`minimum_segment_distance_squared`), `:533` (`support_can_reach_tool`), `:775`, `:878`, `:920`.
- `segment_fibre.cpp`: `:349,:351,:353,:356,:359` (`point_numerators`), `:1098,:1100,:1102,:1105,:1108` (vertex replay), `:1209`, `:1219,:1220,:1226,:1227` (the `x_coefficients` / `y_coefficients` text vectors), `:1234,:1236` (the `shifted` lambda that reads those vectors straight back).
- `segment_strata.cpp`: `:804,:806,:808,:810,:819` (`segment_branches_at`), `:830` (`segment_branch_pair_dispositions`), `:857` (`construct_segment_cell_stratum`). The four station-owned calls at `:873-:884` belong to **stage 2**.
- `segment_oracle.cpp`: `:77,:79,:81,:83,:86,:90` (`station_source`), `:206,:209,:212` (`start_disk_has_no_material_interior`), `:226` (`cap_is_exact_pi`).
- `segment_pair_projection.cpp`: `:341` (`derive_segment_pair_projections`).

After stage 3 the decoder definitions in `segment_strata.cpp`, `segment_oracle.cpp` and `segment_pair_projection.cpp` become **dormant** (definition only, zero callers) — three more on top of the one stage 2 makes dormant. `segment_projection.cpp` keeps 13 callers and `segment_fibre.cpp` keeps 8, all decoding `record.primitive_coefficients` from `BoundaryFeatureRecord2`, which is a different lane owned by stages 4 and 5. **No definition is deleted here; that is stage 6.**

### Two text reads survive on purpose

`segment_projection.cpp:790` passes `source.motion_data()` and `:793` passes `source.tool_radius().text()` into `construct_pullback`. That function is a **Python-bound public API** (`parameter_charts.h:13`, `_continuous_tea_2.pyi:560`) whose parameters are `std::vector<std::string> motion_data` and `std::string cutter_radius`, and it parses them back at `parameter_charts.cpp:721`. Converting it changes the frozen Python API, so stage 3 re-points both reads at the attestation view and leaves the crossing standing. It is recorded in the ratchet test as a named remainder.

### The Python surface — bound, unlike the station lane

Stage 2 found `StationEventSource2` unbound. **The segment lane is the opposite**: `src/continuous_tea_2.cpp:276-347` binds both `ExactBinary64Rational2` and `SegmentEventSource2`, and `src/compas_cgal/_continuous_tea_2.pyi:101-127` declares them.

| Python name | Binding | What must not change |
|---|---|---|
| `ExactBinary64Rational2.numerator` / `.denominator` / `.text` / `.canonical_bytes` | `continuous_tea_2.cpp:276-293` | the class, all four members |
| `SegmentEventSource2.from_binary64(x0, y0, x1, y1, tool_radius, cap_chord_ratio)` | `:298-306` | signature and keyword names |
| `.x0 .y0 .x1 .y1 .tool_radius .cap_chord_ratio` | `:307-330`, each `nb::rv_policy::reference_internal` returning `const ExactBinary64Rational2&` | the returned **Python type** stays `ExactBinary64Rational2` |
| `.motion_data` → `tuple[str, str, str, str]` | `:331-336` | shape and contents |
| `.canonical_bytes`, `.canonical_digest` → `bytes` | `:337-347` | exact bytes |

`rv_policy::reference_internal` is the constraint that shapes the C++ design: it needs a reference that outlives the call, so `attestation()` must return `const SegmentAttestation2&` from a member built at construction. **This is a deliberate deviation from the design document**, which sketches `SegmentAttestation2 attestation() const;` returning by value. By value would dangle under `reference_internal`. Stage 2 reached the same shape for the same structural reason.

Python reference tests that pin this surface and must stay green untouched:

- `tests/adaptive/test_exact_binary64_contract.py` — Hypothesis over floats plus edge doubles, oracle `fractions.Fraction`, reading `source.x0.numerator` / `.denominator`.
- `tests/adaptive/test_segment_event_substrate.py:58-70` — `test_segment_source_exact_lifts_each_binary64_once` pins `x0.numerator == "3602879701896397"` and `x0.denominator == "36028797018963968"` for `from_binary64(0.1, -0.5, 1.25, 2.0, 0.25, 4.0)`.

### Digests that embed these bytes

| Site | Record |
|---|---|
| `segment_partition.cpp:209` | segment event partition canonical bytes |
| `segment_partition.cpp:449`, `:704` | partition certificate records |
| `segment_partition.cpp:809-810` | identity comparison between a candidate partition and the source |
| `segment_oracle.cpp:305` | `EventTrace2::source_canonical_bytes` |
| `segment_oracle.cpp:306` | `EventTrace2::effective_cap_bytes` = `source.cap_chord_ratio().canonical_bytes()` — an **attestation** read |
| `segment_oracle.cpp:605-620` | swept-prefix audit record |
| `continuous_tea_2.cpp:341,:347` | the Python `.canonical_bytes` / `.canonical_digest` properties |

None of these signatures change. `segment_partition.cpp` and the three `audit_*_2.cpp` files construct through `from_exact`, whose signature is already exact, and read only `canonical_bytes()` — so **they need no edit at all**, and neither do `tests/native/test_audit_certification_{contract,identity,refinement}_2.cpp` (`grep` for the six accessors in those three files returns nothing).

### The anchor gap this stage closes

`grep -rn "segment-event-source-v1" src tests docs` returns exactly **one** line — its own producer at `segment_source.cpp:198`. There is no absolute anchor for the outer framing anywhere in the repository: only composite digests, which detect drift without localising it. `exact-binary64-rational-v1` is anchored once, at `tests/native/exact_canonical_gate.cpp:66-74`, against `exact::CanonicalRational::canonical_bytes()` — the very method stage 3 removes. Task 1 relocates that anchor onto `ExactBinary64Rational2` and adds the missing outer one, **before** Task 2 deletes the method.

### `exact::CanonicalRational::canonical_bytes()` — verified against the tree

`grep -rn "CanonicalRational\|to_canonical" src tests | grep -v '^src/exact/'` returns hits in only two files: `tests/native/exact_canonical_gate.cpp` and `tests/native/exact_station_attestation_gate.cpp` (stage 2's, which uses `numerator()`/`denominator()`/`text()` only). The method itself is called from exactly two places, both in `exact_canonical_gate.cpp`:

- `:81` — `exact::to_canonical(exact::from_binary64(0.1)).canonical_bytes()` against the 91-byte literal;
- `:115` — `promoted.canonical_bytes() == source.x0().canonical_bytes()`.

**Zero production callers.** The brief's account is confirmed.

---

## Dependencies on stage 2, which is in flight

Stage 2 is executing concurrently, and it moved while this plan was being written. Two observations, both from the working tree:

- **At `45c5cc8a`:** Task 1 had landed (`tests/native/exact_station_attestation_gate.cpp`, `CMakeLists.txt:261-267`, `pyproject.toml:286-287`); `station_source.h` still declared the string-carrying `ExactRational2::build`.
- **A few commits later, uncommitted in the tree:** Tasks 2 and 3 have landed too. `station_source.h:94-98` now declares `StationEventSource2::build(const compas_cgal::exact::Rational&, ...)` ×4; `segment_strata.cpp` is at **8** `parse_rational` occurrences, down from 12; `segment_oracle.cpp::station_source` carries stage 2's rewritten body and the file already has a `namespace exact = compas_cgal::exact;` alias.

That means D1, D2 and D3 below are **satisfied as observed**, not merely predicted. Re-check them anyway with the one-line greps each row names — stage 2 may still amend its own work, and none of it is committed.

| # | What stage 3 rests on | Where it bites | If stage 2 lands differently |
|---|---|---|---|
| D1 | `segment_strata.cpp` drops to **8** `parse_rational` occurrences (stage 2 removes the four station-owned calls at `:873-:884`) | Task 8's ratchet baseline and expected value | Re-derive the expected count from the observed dict; the *code* change in Task 5 Step 3 is unaffected. |
| D2 | `segment_oracle.cpp::station_source` is rewritten by stage 2 Task 3 Step 4 to call `StationEventSource2::build` with four `exact::Rational`, and the four-use `text` lambda at `:92-103` is gone | Task 5 Step 4 replaces that function *again*, from stage 2's version | If stage 2 has not landed, replace the **current** `:68-109` body with the same final text and note that stage 2's intermediate never existed. The final code is identical either way. |
| D3 | `StationEventSource2::build` takes `const exact::Rational&` ×4 | Task 5 Step 4's `station_source` builds them with no `.text()` round trip | If stage 2 has not landed, Task 5 Step 4 must keep `exact_rational(...)` → `text(...)` on the *output* side only. **Do not start Task 5 before stage 2's Task 3 is committed**; sequence around it rather than duplicating the work. |
| D4 | `segment_strata.cpp::construct_station_cell_stratum` uses `CGAL::exact(source.cap_chord_ratio())` | Task 5 Step 3 edits neighbouring functions in the same file | Pathspec-commit; if a merge conflict appears in that function, stage 2 owns it. |
| D5 | `docs/number_types.md` gains a "Status at stage 2" admonition replacing "Status after stage 1" | Task 9 Step 1 replaces it with a stage-3 status | Address the block by its heading text, never by line number. |
| D6 | `AttestationByteDriftError` acquires its first raiser in `exact_station_attestation_gate` | Task 9 Step 2's error-table wording | If stage 2's Task 6 has not run, stage 3's gate is the first raiser and the wording must say so. |
| D7 | `pyproject.toml:286-287` already lists `exact_station_attestation_gate` | Task 1 Step 3 appends to the same two lines | Append; do not rewrite the lines from this plan's text verbatim without diffing first. |

**Nothing in stage 3 depends on a stage-2 outcome that has not already been observed in the tree**, apart from the counts in D1 and the intermediate form in D2 — both of which are checked by a step that prints the observed value before acting.

---

## Global Constraints

- Exact arithmetic is settled; do not introduce any epsilon, tolerance, deflation factor or `nextafter` into a decision path.
- **Python API and canonical digest bytes are FROZEN.** The user decided this. Derived attestation views survive; they stop being the compute carrier. `src/compas_cgal/_continuous_tea_2.pyi` must be byte-identical at the end of this stage.
- Named exceptions only, convention `<Domain><Condition>Error`. Never `throw std::runtime_error("...")` directly. The segment lane's existing named errors — `NonFiniteSegmentInputError`, `ZeroLengthSegmentMotionError`, `NonPositiveToolRadiusError`, `InvalidCapChordRatioError`, `NonFiniteFullCircleInputError`, `ZeroFullCirclePhaseError` — keep their exact identities and messages.
- Do **not** delete any `parse_rational` definition, `ExactBinary64Rational2`'s text-producing members, or any other decoder. Removal is stage 6 and requires explicit user permission. Dormant decoders get `[[maybe_unused]]` and a comment naming stage 6.
- The one deletion stage 3 owns is `exact::CanonicalRational::canonical_bytes()`, which has zero production callers and is scheduled by the design's corrected invariant 2. It is not a decoder and not a `parse_rational` definition.
- No `pytest.mark.skip`, `skipif`, or `xfail` under any circumstance. A failing test fails.
- Never modify an existing reference test to make new code pass. `tests/adaptive/test_exact_binary64_contract.py`, `tests/adaptive/test_segment_event_substrate.py:58-70` and `LEGACY_FULL_TRACE_SHA256` in `tests/adaptive/test_circle_oracle.py:27-34` are reference anchors: if stage 3 moves one, stage 3 is wrong.
- Native gates compile under `CMAKE_BUILD_TYPE Release`, which defines `NDEBUG` and erases `<cassert>`. Use the throwing `require` + `GateCheckFailedError` pattern from `tests/native/exact_canonical_gate.cpp:18-31`; never `assert`.
- `CGAL_DISABLE_GMP` and `CGAL_USE_BOOST_MP` are set repo-wide (`CMakeLists.txt:109-118`). `exact::Rational` is `Epeck::FT` = `Lazy_exact_nt<Epeck_ft>`, and `Epeck_ft` **is** `CORE::BigRat` — `station_classifier.cpp:46-49` static-asserts exactly that. Every `CORE::BigRat` ↔ `exact::Rational` conversion in this plan is an exact rewrap, not a reinterpretation.
- **Stage 3 does not make the segment lane lazy.** `segment_projection.cpp`, `segment_fibre.cpp`, `segment_strata.cpp`, `segment_oracle.cpp` and `segment_pair_projection.cpp` all compute in eager `CORE::BigRat`. Each converted call site therefore materialises through a **declared `.exact()` boundary** (`docs/number_types.md` R7) — `CGAL::exact(...)`, or the file-local `exact_rational` helper where one already exists at `segment_projection.cpp:77`, `segment_fibre.cpp:84`, `segment_strata.cpp:70`. What stage 3 removes is the decode-and-renormalise round trip, not the eager arithmetic. Making those lanes lazy is stage 4, and conflating the two would make stage 4's measurement unattributable.
- `${CMAKE_CURRENT_SOURCE_DIR}/src` is a PUBLIC include directory on `continuous_tea_exact_core` (`CMakeLists.txt:352-355`) and `_continuous_tea_2` links that target, so `#include "exact/rational.h"` resolves from `src/continuous_tea_2/` in both the static library and the nanobind module.
- Adding a `.cpp` file or a target requires editing `CMakeLists.txt` and a full rebuild.
- Run pytest through the pixi task, arguments after `--`: `pixi run -e default pytest -- <paths> -n auto`.
- **Shared worktree.** `CMakeLists.txt`, `pyproject.toml` and `docs/number_types.md` carry other sessions' uncommitted edits. Commit by pathspec with the message *before* the `--`: `git commit -m "msg" -- <paths>`. Before every commit that includes a shared file, run `git diff HEAD -- <that file>` and confirm every hunk is yours; a pathspec isolates files, not hunks. If a foreign hunk is present, stop and report rather than committing it.
- **Disk has hit 100% twice today.** Reuse `build/exact-gate` and the existing editable build; create no additional build directories.
- Commit messages: extremely concise, lowercase, no attribution trailers.

---

### Task 1: Anchor the segment bytes and digest, and capture the lane baseline, before touching production code

Nothing in the repository pins `segment-event-source-v1` to an absolute constant, and the only anchor for `exact-binary64-rational-v1` lives on the method Task 2 deletes. This task lands both anchors against the **current** text-carrier code, plus the pre-inversion behavioural baseline that Task 8 compares to. No production file changes.

**Files:**
- Create: `tests/native/exact_segment_attestation_gate.cpp`
- Modify: `CMakeLists.txt` (add the `exact_segment_attestation_gate` executable)
- Modify: `pyproject.toml` (add the gate to the `exact-gates` task)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational`, `compas_cgal::exact::from_binary64`, `compas_cgal::exact::to_canonical`, `compas_cgal::exact::CanonicalRational`, `compas_cgal::exact::AttestationByteDriftError` (stage 0); `SegmentEventSource2::from_binary64`, `SegmentEventSource2::from_exact`, `ExactBinary64Rational2::canonical_bytes` (current API).
- Produces: `build/exact-gate/exact_segment_attestation_gate`, the absolute byte and digest anchors for the segment lane, and `build/stage3/baseline.xml`, the pre-inversion per-test outcome record.

- [ ] **Step 1: Write the gate**

Create `tests/native/exact_segment_attestation_gate.cpp`:

```cpp
#include "exact/canonical.h"
#include "exact/errors.h"
#include "exact/rational.h"

#include "continuous_tea_2/segment_source.h"
#include "continuous_tea_2/sha256.h"

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace exact = compas_cgal::exact;

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

/// Raised when a gate check fails. Named so the gate obeys the same
/// named-exceptions-only rule as the module it exercises.
class GateCheckFailedError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// CMAKE_BUILD_TYPE is Release, which defines NDEBUG and erases <cassert>. A
// gate whose checks vanish reports a vacuous pass, so every decision here goes
// through a helper that is never compiled out.
void require(bool condition, const char* message)
{
    if (!condition) {
        throw GateCheckFailedError(message);
    }
}

/// Byte comparisons raise the module's own drift error rather than the generic
/// gate error, so a digest-invalidating change is distinguishable from an
/// ordinary gate failure by exception type.
void require_bytes(
    const std::string& produced,
    const std::string& frozen,
    const char* message)
{
    if (produced != frozen) {
        throw exact::AttestationByteDriftError(message);
    }
}

std::string hex(const std::string& bytes)
{
    static const char* const digits = "0123456789abcdef";
    std::string result;
    result.reserve(bytes.size() * 2);
    for (const char byte : bytes) {
        const auto value = static_cast<unsigned char>(byte);
        result.push_back(digits[value >> 4U]);
        result.push_back(digits[value & 0x0fU]);
    }
    return result;
}

/// The frozen `exact-binary64-rational-v1` record for the double 0.1.
///
/// FROZEN COMPATIBILITY CONTRACT. These 91 bytes are the on-disk shape every
/// stored segment replay digest was computed over: an 8-byte big-endian field
/// count, then per field an 8-byte big-endian length followed by its payload.
/// 0.1 is chosen because it is not a decimal fraction, so its exact value
/// 3602879701896397 / 2^55 exercises a long numerator and denominator rather
/// than a round one.
///
/// RELOCATED from `exact_canonical_gate.cpp:66-74`, where it was checked
/// against `exact::CanonicalRational::canonical_bytes()`. That method is
/// removed in Task 2 -- it put one consumer's record tag in the base number
/// vocabulary -- so the anchor moves onto `ExactBinary64Rational2`, which is
/// where this framing actually belongs and is the only place it is emitted.
///
/// The literal exists because the differential checks below CANNOT see a change
/// here: every projection shares the one `encode_string_sequence`, so a single
/// edit to that framing moves both sides identically and leaves the differential
/// green while every stored digest silently stops verifying.
///
/// If a change makes this check fail, that change invalidates every stored
/// segment replay digest. It must be a deliberate, versioned decision: bump the
/// "exact-binary64-rational-v1" tag and migrate the stored digests. Editing this
/// literal to restore green is the one repair that is always wrong.
constexpr char kFrozenTenthBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                                // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x1a" "exact-binary64-rational-v1"   // 26 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x10" "3602879701896397"             // 16 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x11" "36028797018963968";           // 17 bytes

static_assert(
    sizeof(kFrozenTenthBytes) - 1 == 91,
    "the frozen literal must be 91 bytes: 8 + (8+26) + (8+16) + (8+17)");

/// The same framing for a value no binary64 can denote.
///
/// `from_exact` admits any exact rational, so the attestation view must be
/// pinned on a non-dyadic value too. -355/113 is reduced, non-dyadic, and
/// carries its sign in the numerator -- which is what canonicalisation must
/// guarantee. `exact_canonical_gate` cannot reach this case at all, because
/// every value it builds comes from `from_binary64`.
constexpr char kFrozenNonDyadicBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                                // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x1a" "exact-binary64-rational-v1"   // 26 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x04" "-355"                         // 4 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x03" "113";                         // 3 bytes

static_assert(
    sizeof(kFrozenNonDyadicBytes) - 1 == 65,
    "the frozen literal must be 65 bytes: 8 + (8+26) + (8+4) + (8+3)");

/// The SHA-256 of the frozen `segment-event-source-v1` record for the probe.
///
/// DERIVED, NOT OBSERVED. It was computed from the encoder's specification
/// (`event_certificate.cpp:593-602`) rather than read off a run, so Step 4
/// requires the structural anchor below to pass FIRST. If the structural anchor
/// passes and only this hex differs, the derivation was wrong and this literal
/// is the thing to correct -- print the observed hex, record it in the commit
/// message, and move on. If the structural anchor ALSO fails, the framing has
/// moved and nothing here may be edited.
constexpr char kFrozenSegmentDigestHex[] =
    "0b5f57185b5c54a528ae4224c87b04735b4438fdedb71ac13f44e184e1128389";

/// An independent re-implementation of the canonical length framing.
///
/// Deliberately NOT `encode_string_sequence`: the outer
/// `segment-event-source-v1` record nests six inner records and would be
/// unreadable as a raw literal, so it is anchored structurally instead. The
/// equality check in `frozen_framing_agrees_with_the_literals` ties this helper
/// to the raw literals above, so together they pin the framing bytes, the field
/// order and the field count without either standing alone.
std::string frozen_u64(std::size_t value)
{
    std::string result(8, '\0');
    for (std::size_t index = 0; index < 8; ++index) {
        const int shift = static_cast<int>(56 - 8 * index);
        result[index] = static_cast<char>(
            (static_cast<std::uint64_t>(value) >> shift) & 0xffU);
    }
    return result;
}

std::string frozen_sequence(const std::vector<std::string>& fields)
{
    std::string result = frozen_u64(fields.size());
    for (const std::string& field : fields) {
        result += frozen_u64(field.size());
        result += field;
    }
    return result;
}

std::string frozen_value(
    const std::string& numerator,
    const std::string& denominator)
{
    return frozen_sequence(
        {"exact-binary64-rational-v1", numerator, denominator});
}

/// The probe motion: (0.1, -0.25) -> (1.5, -3.0), tool radius 0.5, cap 4.0.
///
/// 0.1 gives a long dyadic numerator and denominator; -0.25 and -3.0 put the
/// sign in the numerator on both a fractional and an integral value; 1.5 and 0.5
/// are short dyadics; 4.0 is the CLOSED upper end of the cap interval, so the
/// probe also proves the boundary is admitted. The endpoints differ and the
/// radius is positive, so the source's own validation admits it.
///
/// TASK 3 REWRITES NOTHING HERE. The frozen expectations below never move.
SegmentEventSource2 probe_segment()
{
    return SegmentEventSource2::from_binary64(
        0.1, -0.25, 1.5, -3.0, 0.5, 4.0);
}

/// The attested view of a source value.
///
/// TASK 3 REWRITES ONLY THESE TWO FUNCTION BODIES, to read the derived
/// attestation instead of the carrier -- because after the inversion `x0()`
/// returns `exact::Rational`, which has no `numerator()` and no
/// `canonical_bytes()`. Every frozen expectation in this file is written against
/// these two accessors and never moves. Nothing else in the gate refers to a
/// source value.
const ExactBinary64Rational2& attested_x0(
    const SegmentEventSource2& source)
{
    return source.x0();
}

const ExactBinary64Rational2& attested_y1(
    const SegmentEventSource2& source)
{
    return source.y1();
}

std::string frozen_segment_bytes()
{
    return frozen_sequence({
        "segment-event-source-v1",
        frozen_value("3602879701896397", "36028797018963968"),
        frozen_value("-1", "4"),
        frozen_value("3", "2"),
        frozen_value("-3", "1"),
        frozen_value("1", "2"),
        frozen_value("4", "1"),
    });
}

/// Render a CORE::BigRat the way the legacy text carrier rendered it.
std::string legacy_text(const Rational& value)
{
    const Integer numerator = CORE::numerator(value);
    const Integer denominator = CORE::denominator(value);
    return denominator == 1
        ? numerator.convert_to<std::string>()
        : numerator.convert_to<std::string>()
            + "/"
            + denominator.convert_to<std::string>();
}

/// Tie the structural helper to the raw literals, so neither anchor stands
/// alone. A helper that shared `encode_string_sequence` would move with it; a
/// raw literal cannot express the nested outer record.
void frozen_framing_agrees_with_the_literals()
{
    require_bytes(
        frozen_value("3602879701896397", "36028797018963968"),
        std::string(kFrozenTenthBytes, sizeof(kFrozenTenthBytes) - 1),
        "the gate's independent framing helper disagrees with the frozen 0.1 "
        "literal; one of the two has been edited");
    require_bytes(
        frozen_value("-355", "113"),
        std::string(kFrozenNonDyadicBytes, sizeof(kFrozenNonDyadicBytes) - 1),
        "the gate's independent framing helper disagrees with the frozen "
        "-355/113 literal; one of the two has been edited");
}

/// The `exact-binary64-rational-v1` record of one attested value must not move.
void value_record_matches_the_frozen_literals()
{
    const SegmentEventSource2 dyadic = probe_segment();
    require_bytes(
        attested_x0(dyadic).canonical_bytes(),
        std::string(kFrozenTenthBytes, sizeof(kFrozenTenthBytes) - 1),
        "exact-binary64-rational-v1 bytes differ from the frozen 0.1 literal: "
        "every stored segment replay digest is invalidated");
    require(
        attested_x0(dyadic).numerator() == "3602879701896397",
        "attested numerator of 0.1 is wrong");
    require(
        attested_x0(dyadic).denominator() == "36028797018963968",
        "attested denominator of 0.1 is wrong");
    require(
        attested_x0(dyadic).text() == "3602879701896397/36028797018963968",
        "attested text of 0.1 is wrong");
    require(
        attested_y1(dyadic).text() == "-3",
        "an integral attested value must render without a denominator");

    // The non-dyadic case, reachable only through from_exact.
    const exact::Rational minus_355_over_113 =
        exact::Rational(-355) / exact::Rational(113);
    const SegmentEventSource2 non_dyadic =
        SegmentEventSource2::from_exact(
            ExactSegmentMotion2{
                EPoint(minus_355_over_113, exact::Rational(0)),
                EPoint(exact::Rational(1), exact::Rational(1)),
            },
            exact::Rational(1) / exact::Rational(2),
            exact::Rational(7) / exact::Rational(2));
    require_bytes(
        attested_x0(non_dyadic).canonical_bytes(),
        std::string(kFrozenNonDyadicBytes, sizeof(kFrozenNonDyadicBytes) - 1),
        "exact-binary64-rational-v1 bytes differ from the frozen -355/113 "
        "literal: every stored segment replay digest is invalidated");
}

/// The outer `segment-event-source-v1` record must not move either: its tag,
/// its field count and the order of its six values are all frozen.
void segment_record_matches_the_frozen_framing()
{
    const std::string frozen = frozen_segment_bytes();
    require(
        frozen.size() == 480,
        "the frozen segment record must be 480 bytes: "
        "8 + (8+23) + (8+91) + (8+61) + (8+60) + (8+61) + (8+60) + (8+60)");
    require_bytes(
        probe_segment().canonical_bytes(),
        frozen,
        "segment-event-source-v1 bytes differ from the frozen framing: every "
        "stored segment replay digest is invalidated");
}

/// The digest is the value the audit records actually carry, so it gets its own
/// absolute anchor as well as a tie back to the bytes.
void segment_digest_matches_the_frozen_hex()
{
    const SegmentEventSource2 source = probe_segment();
    require(
        source.canonical_digest() == sha256_bytes(source.canonical_bytes()),
        "canonical_digest is not the SHA-256 of canonical_bytes");
    const std::string produced = hex(source.canonical_digest());
    if (produced != kFrozenSegmentDigestHex) {
        std::printf(
            "observed segment-event-source-v1 digest: %s\n", produced.c_str());
        throw exact::AttestationByteDriftError(
            "segment-event-source-v1 digest differs from the frozen hex");
    }
}

/// The projection door must agree with the current attestation path on values
/// the segment lane actually carries: arbitrary reduced rationals, not just the
/// dyadic ones a binary64 produces.
void projection_matches_the_attestation_on_generic_rationals()
{
    const std::vector<std::pair<const char*, const char*>> fractions = {
        {"1", "3"},
        {"-355", "113"},
        {"1000000000000000000000000000001", "3"},
        {"-987654321098765432109876543210987",
         "1000000000000000000000000000000007"},
        {"246", "369"},   // a common factor that must be divided out: 2/3
        {"5", "-8"},      // a sign that must migrate to the numerator
        {"0", "17"},      // zero, whose canonical denominator is 1
        {"7", "1"},       // an integer
    };
    for (const auto& [numerator, denominator] : fractions) {
        const Rational eager(Integer(numerator), Integer(denominator));
        const exact::Rational lazy(eager);
        const exact::CanonicalRational projected = exact::to_canonical(lazy);
        const SegmentEventSource2 source =
            SegmentEventSource2::from_exact(
                ExactSegmentMotion2{
                    EPoint(lazy, exact::Rational(0)),
                    EPoint(lazy + exact::Rational(1), exact::Rational(1)),
                },
                exact::Rational(1),
                exact::Rational(1));
        require(
            projected.numerator() == attested_x0(source).numerator(),
            "projected numerator differs from the attestation view");
        require(
            projected.denominator() == attested_x0(source).denominator(),
            "projected denominator differs from the attestation view");
        require(
            projected.text() == attested_x0(source).text(),
            "projected text differs from the attestation view");
        require(
            attested_x0(source).text() == legacy_text(eager),
            "the attestation text differs from the eager rendering");
    }
}

/// The producers hand the source values built by a CHAIN of exact arithmetic,
/// not leaves. The claim stage 3 rests on is that the lazy carrier's exact
/// representative is bit-identical to the eager one, so the same chain computed
/// both ways must canonicalise to the same bytes. Depth 16 is chosen because
/// that is where the measured lazy/eager divergence risk would be largest --
/// see docs/number_types.md, "What the filter is actually worth".
void lazy_and_eager_chains_canonicalise_identically()
{
    exact::Rational lazy(1);
    Rational eager(1);
    for (int step = 0; step < 16; ++step) {
        lazy = (lazy * exact::Rational(3) + exact::Rational(7))
            / (lazy + exact::Rational(5));
        eager = (eager * Rational(3) + Rational(7))
            / (eager + Rational(5));
    }
    const exact::CanonicalRational projected = exact::to_canonical(lazy);
    require(
        projected.text() == legacy_text(eager),
        "a depth-16 lazy chain canonicalises differently from the same chain "
        "computed eagerly");
    require(
        CORE::numerator(CGAL::exact(lazy)) == CORE::numerator(eager)
            && CORE::denominator(CGAL::exact(lazy)) == CORE::denominator(eager),
        "the lazy carrier's exact representative is not the eager value");
}

/// `motion_data` is the attestation shape `construct_pullback` consumes, and it
/// is Python-visible as a 4-tuple. Both facts are frozen.
void motion_data_is_the_four_attested_endpoint_texts()
{
    const std::vector<std::string> data = probe_segment().motion_data();
    require(data.size() == 4, "motion_data must carry exactly four fields");
    require(
        data[0] == "3602879701896397/36028797018963968",
        "motion_data[0] is not the attested x0 text");
    require(data[1] == "-1/4", "motion_data[1] is not the attested y0 text");
    require(data[2] == "3/2", "motion_data[2] is not the attested x1 text");
    require(data[3] == "-3", "motion_data[3] is not the attested y1 text");
}

/// The source's own validation is part of its contract and must survive the
/// signature change unchanged, with the same named exception per condition.
void invalid_segments_are_rejected()
{
    const double infinity = 1.0 / 0.0;
    bool raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            infinity, 0.0, 1.0, 1.0, 1.0, 1.0);
    } catch (const NonFiniteSegmentInputError&) {
        raised = true;
    }
    require(raised, "a non-finite coordinate was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            1.0, 2.0, 1.0, 2.0, 1.0, 1.0);
    } catch (const ZeroLengthSegmentMotionError&) {
        raised = true;
    }
    require(raised, "coincident endpoints were accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            0.0, 0.0, 1.0, 1.0, 0.0, 1.0);
    } catch (const NonPositiveToolRadiusError&) {
        raised = true;
    }
    require(raised, "a zero tool radius was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_binary64(
            0.0, 0.0, 1.0, 1.0, 1.0, 4.25);
    } catch (const InvalidCapChordRatioError&) {
        raised = true;
    }
    require(raised, "a cap chord ratio above 4 was accepted");

    raised = false;
    try {
        (void)SegmentEventSource2::from_exact(
            ExactSegmentMotion2{
                EPoint(exact::Rational(1), exact::Rational(2)),
                EPoint(exact::Rational(1), exact::Rational(2)),
            },
            exact::Rational(1),
            exact::Rational(1));
    } catch (const ZeroLengthSegmentMotionError&) {
        raised = true;
    }
    require(raised, "exact coincident endpoints were accepted");

    // The closed upper end is admitted, and so are negative coordinates.
    (void)SegmentEventSource2::from_binary64(
        -1.5, -2.5, -1.25, -2.0, 0.125, 4.0);
}

/// The full-circle source frames three records of its own. It already carries
/// exact values, so Task 6 changes only WHICH function performs the projection;
/// these anchors are what proves it changed nothing else.
void full_circle_records_do_not_move()
{
    const FullCircleEventSource2 binary64 =
        FullCircleEventSource2::from_binary64(
            0.5, -0.25, 1.0, 0.0, false, 0.5, 4.0);
    const std::string binary64_canonical = binary64.canonical_bytes();
    const std::string binary64_motion = binary64.motion_identity_bytes();
    const std::string binary64_cap = binary64.cap_identity_bytes();

    const FullCircleEventSource2 exact_source =
        FullCircleEventSource2::from_exact(
            ExactCircleMotion2{
                EPoint(
                    exact::Rational(-355) / exact::Rational(113),
                    exact::Rational(22) / exact::Rational(7)),
                EVector(exact::Rational(1), exact::Rational(0)),
                true,
            },
            exact::Rational(1) / exact::Rational(2),
            exact::Rational(7) / exact::Rational(2));

    require(
        binary64_canonical.find("full-circle-event-source-binary64-v1")
            != std::string::npos,
        "the binary64 full-circle record lost its tag");
    require(
        binary64_motion.find("full-circle-motion-binary64-v1")
            != std::string::npos,
        "the binary64 motion identity record lost its tag");
    require(
        binary64_cap.find("cap-chord-ratio-binary64-v1") != std::string::npos,
        "the binary64 cap identity record lost its tag");
    require(
        exact_source.canonical_bytes().find(
            "full-circle-event-source-exact-v1") != std::string::npos,
        "the exact full-circle record lost its tag");

    // The exact framing renders each value as canonical text. -355/113 and 22/7
    // are non-dyadic and 1/2, 7/2 are dyadic, so this pins both renderings.
    require(
        exact_source.canonical_bytes().find("-355/113") != std::string::npos,
        "the exact full-circle record does not carry the canonical centre x");
    require(
        exact_source.canonical_bytes().find("22/7") != std::string::npos,
        "the exact full-circle record does not carry the canonical centre y");
    require(
        exact_source.cap_identity_bytes().find("7/2") != std::string::npos,
        "the exact cap identity record does not carry the canonical cap");
    require(
        exact_source.motion_identity_bytes().find("clockwise")
            != std::string::npos,
        "the exact motion identity record lost its orientation field");
}

}  // namespace

int main()
{
    try {
        frozen_framing_agrees_with_the_literals();
        value_record_matches_the_frozen_literals();
        segment_record_matches_the_frozen_framing();
        segment_digest_matches_the_frozen_hex();
        projection_matches_the_attestation_on_generic_rationals();
        lazy_and_eager_chains_canonicalise_identically();
        motion_data_is_the_four_attested_endpoint_texts();
        invalid_segments_are_rejected();
        full_circle_records_do_not_move();
    } catch (const std::exception& error) {
        std::printf(
            "exact_segment_attestation_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_segment_attestation_gate OK\n");
    return 0;
}
```

- [ ] **Step 2: Wire the gate into CMake**

In `CMakeLists.txt`, directly after the `exact_station_attestation_gate` block (which ends at line 267 in the current tree — confirm with `grep -n "exact_station_attestation_gate" CMakeLists.txt` before inserting), add:

```cmake
# Absolute byte anchors for the segment lane. `segment-event-source-v1` appears
# exactly once in the repository -- its own producer -- so before this gate
# nothing pinned the outer framing at all, and the only anchor for
# exact-binary64-rational-v1 lived on exact::CanonicalRational::canonical_bytes(),
# which stage 3 removes.
add_executable(exact_segment_attestation_gate EXCLUDE_FROM_ALL
    tests/native/exact_segment_attestation_gate.cpp
)
target_link_libraries(
    exact_segment_attestation_gate PRIVATE continuous_tea_exact_core)
target_compile_definitions(
    exact_segment_attestation_gate PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 3: Wire the gate into the `exact-gates` task**

`pyproject.toml:286-287` currently reads (verify with `sed -n '285,288p' pyproject.toml` — stage 2 edited these lines and a third session may have too):

```toml
_exact-gates-build = "cmake --build build/exact-gate --target exact_vocabulary_gate exact_canonical_gate exact_one_root_gate sign_mixed_radical_gate exact_station_attestation_gate"
_exact-gates-run = "build/exact-gate/exact_vocabulary_gate && build/exact-gate/exact_canonical_gate && build/exact-gate/exact_one_root_gate && build/exact-gate/sign_mixed_radical_gate && build/exact-gate/exact_station_attestation_gate"
```

**Append**, do not rewrite: add ` exact_segment_attestation_gate` to the end of the `--target` list, and ` && build/exact-gate/exact_segment_attestation_gate` to the end of the run line. Then update the `description` on line 288 to end with `, segment attestation bytes`.

- [ ] **Step 4: Run the gate against the current, un-inverted code**

```bash
pixi run -e default exact-gates
```

Expected: six `OK` lines ending with `exact_segment_attestation_gate OK`, exit status 0.

Decision rules, in this order:

1. If `frozen_framing_agrees_with_the_literals`, `value_record_matches_the_frozen_literals` or `segment_record_matches_the_frozen_framing` fails, the frozen constants in this plan are wrong and must be corrected from the *observed* bytes — print them with `xxd` on the produced string and recompute — because the current bytes are by definition the contract. **Stop and report the observed bytes before changing anything else.**
2. If only `segment_digest_matches_the_frozen_hex` fails, the structural anchor already passed, so the framing is intact and the derived hex was mis-derived. The gate prints the observed hex. Replace `kFrozenSegmentDigestHex` with it and record the substitution verbatim in the Step 6 commit message. This is the one literal in this gate that may be corrected without stopping.
3. If `projection_matches_the_attestation_on_generic_rationals` or `lazy_and_eager_chains_canonicalise_identically` fails, **stop and report the differing pair**. That is the premise stage 3 rests on, and a failure means the inversion cannot be byte-preserving as designed.
4. If `full_circle_records_do_not_move` fails, the full-circle framings are not what this plan read them to be; **stop and report** before Task 6.

- [ ] **Step 5: Capture the pre-inversion behavioural baseline**

Capture the baseline **now**, from this tree, before any production file changes.

```bash
mkdir -p build/stage3
git rev-parse HEAD | tee build/stage3/start-commit.txt
pixi run -e default pytest -- \
  tests/adaptive/test_segment_event_substrate.py \
  tests/adaptive/test_segment_event_proof_contracts.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_exact_binary64_contract.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_motion_refutation.py \
  tests/adaptive/test_motion_oracle_cache.py \
  tests/adaptive/test_transaction.py \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_event_substrate.py \
  tests/adaptive/test_generator.py \
  -n auto -q --junitxml=build/stage3/baseline.xml
```

Expected: a run that completes and writes `build/stage3/baseline.xml`. The pass/fail counts do **not** need to be all-green — this tree has pre-existing failures. What matters is that the file exists and records a per-test outcome for every test in these suites. Record the summary line verbatim for the Task 8 comparison.

`build/` is not tracked, so nothing here is committed.

- [ ] **Step 6: Commit**

```bash
git diff HEAD -- CMakeLists.txt pyproject.toml
```

Confirm every hunk is yours (see Global Constraints). Then:

```bash
git commit -m "test: anchor segment attestation bytes and digest" -- \
  tests/native/exact_segment_attestation_gate.cpp CMakeLists.txt pyproject.toml
```

---

### Task 2: Remove `exact::CanonicalRational::canonical_bytes()`

A base-vocabulary type must not carry one consumer's record tag. The design's corrected invariant 2 establishes that three distinct framings exist, that two of them share a tag, and that a tag-taking method would therefore invite several frozen contracts onto one framer. The repair is removal. The anchor it held moved in Task 1; this task removes the method and the two calls.

**Files:**
- Modify: `src/exact/canonical.h`
- Modify: `src/exact/canonical.cpp`
- Modify: `src/exact/rational.h`
- Modify: `tests/native/exact_canonical_gate.cpp`

**Interfaces:**
- Consumes: nothing new.
- Produces: `compas_cgal::exact::CanonicalRational` with `numerator()`, `denominator()`, `text()` and **no** `canonical_bytes()`. `to_canonical`'s signature is unchanged.

- [ ] **Step 1: Prove there are no other callers**

```bash
grep -rn "canonical_bytes" src/exact tests/native/exact_canonical_gate.cpp tests/native/exact_station_attestation_gate.cpp tests/native/exact_vocabulary_gate.cpp tests/native/exact_one_root_gate.cpp
grep -rn "CanonicalRational" src tests | grep -v '^src/exact/' | grep -v '^tests/native/exact_'
```

Expected: the first prints hits only in `src/exact/canonical.{h,cpp}` and `exact_canonical_gate.cpp:81,:115` (plus `exact_station_attestation_gate.cpp`'s calls on `ExactRational2`, which is a different class); the second prints nothing. If either shows a production caller, **stop and report** — the design's zero-caller claim would be false and the removal is not free.

- [ ] **Step 2: Remove the declaration**

In `src/exact/canonical.h`, delete lines 22-24 (the doc comment and the `canonical_bytes` declaration) so the public section reads:

```cpp
class CanonicalRational {
public:
    [[nodiscard]] const std::string& numerator() const noexcept;
    [[nodiscard]] const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    [[nodiscard]] std::string text() const;

private:
```

and replace the `to_canonical` doc comment (`:35-47`) with one that reflects stage 3:

```cpp
/// Project an exact rational into canonical attestation form.
///
/// The single exact-to-canonical-decimal-strings door for this codebase. Since
/// stage 3 it is reached from `ExactBinary64Rational2::project`
/// (`continuous_tea_2/segment_source.cpp`) and `ExactRational2::project`
/// (`continuous_tea_2/station_source.cpp`), which are the segment and station
/// lanes' only projection sites.
///
/// It returns the DECOMPOSITION only. Record FRAMING -- the tag, the field
/// count, the field order -- is a per-lane compatibility contract and lives with
/// the lane, never here: three framings exist in this repository and two of them
/// share the tag `exact-rational-v1`, so a framer in the number vocabulary would
/// have to speak for contracts it cannot see. `CanonicalRational` deliberately
/// has no `canonical_bytes()` for that reason.
///
/// Raises:
///     UnreducedCanonicalRationalError: if the decomposed denominator is not
///         positive, which would make the encoding ambiguous.
[[nodiscard]] CanonicalRational to_canonical(const Rational& value);
```

- [ ] **Step 3: Remove the definition**

In `src/exact/canonical.cpp`, delete the whole `CanonicalRational::canonical_bytes` function (`:37-45`) and the now-unused include of the framing header on line 4:

```cpp
#include "continuous_tea_2/event_certificate.h"
```

Keep `#include <vector>` only if something else in the file uses it; after this deletion nothing does, so remove line 10 as well. The remaining includes are:

```cpp
#include "exact/canonical.h"

#include "exact/errors.h"

#include <CGAL/CORE/BigInt.h>
#include <CGAL/Fraction_traits.h>

#include <utility>
```

This also severs `src/exact/`'s only dependency on `continuous_tea_2/`, which is the point: the vocabulary module now knows nothing about any event source.

- [ ] **Step 4: Update the two calls in `exact_canonical_gate`**

In `tests/native/exact_canonical_gate.cpp`:

Delete `kFrozenTenthBytes` (`:45-74`), the `static_assert`, the function `canonical_bytes_match_the_frozen_literal` (`:76-90`) and its call in `main` (`:175`). Insert this comment where the literal was, so the relocation is discoverable from here:

```cpp
// The absolute 91-byte anchor for `exact-binary64-rational-v1` used to live
// here, checked against `exact::CanonicalRational::canonical_bytes()`. Stage 3
// removed that method -- a base-vocabulary type must not carry one consumer's
// record tag -- and moved the anchor onto the type that actually emits the
// framing: see `kFrozenTenthBytes` in exact_segment_attestation_gate.cpp.
```

In `to_canonical_matches_existing_projection`, delete the fourth `require` (`:114-116`), the one comparing `canonical_bytes`. The numerator, denominator and text comparisons stay: they are the decomposition claim, which is what `to_canonical` owns.

- [ ] **Step 5: Update the stale stage-0 docstring in `exact/rational.h`**

`src/exact/rational.h:17-22` still says "As of stage 0 nothing is routed through it yet" and points at `segment_source.cpp:127-131` and `:307-312`. Replace `:15-35` with:

```cpp
/// Convert a binary64 to its exact rational value.
///
/// The single double-to-exact door for this codebase. Since stage 3 both
/// binary64 factories in `continuous_tea_2/segment_source.cpp` --
/// `SegmentEventSource2::from_binary64` and
/// `FullCircleEventSource2::from_binary64` -- inject through it, after their own
/// domain validation has run, so each keeps its named error.
///
/// A binary64 is a dyadic rational, so the conversion is exact and total on
/// finite input: there is no parsing, no tolerance and no snapping.
///
/// Args:
///     value: a finite binary64.
///
/// Returns:
///     The exact rational denoted by `value`. Negative zero maps to zero.
///
/// Raises:
///     NonFiniteBinary64Error: if `value` is NaN or infinite. Callers that
///         already reject non-finite input with their own named error never
///         reach this.
[[nodiscard]] Rational from_binary64(double value);
```

- [ ] **Step 6: Build and run the gates**

```bash
pixi run -e default exact-gates
```

Expected: six `OK` lines, exit status 0. `exact_segment_attestation_gate OK` is the proof that removing the method moved no bytes — the anchor it used to hold now passes from its new home.

- [ ] **Step 7: Commit**

```bash
git commit -m "refactor: drop canonical_bytes from the number vocabulary" -- \
  src/exact/canonical.h src/exact/canonical.cpp src/exact/rational.h \
  tests/native/exact_canonical_gate.cpp
```

---

### Task 3: Invert `SegmentEventSource2`

The source stops carrying text and starts carrying `exact::Rational`. `ExactBinary64Rational2` stays, demoted to a derived view with a new `project` factory.

**This task leaves the tree non-compiling** — the bindings and five consumers still expect `ExactBinary64Rational2` accessors. Tasks 4 and 5 fix that, and Tasks 3-5 share one commit made at the end of Task 5. Do not commit at the end of this task.

**Files:**
- Modify: `src/continuous_tea_2/segment_source.h`
- Modify: `src/continuous_tea_2/segment_source.cpp`
- Modify: `tests/native/exact_segment_attestation_gate.cpp` (the two accessor bodies only)
- Modify: `tests/native/exact_canonical_gate.cpp` (the three `source.x0()` reads only)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational`, `compas_cgal::exact::from_binary64` (`src/exact/rational.h`), `compas_cgal::exact::to_canonical`, `compas_cgal::exact::CanonicalRational` (`src/exact/canonical.h`).
- Produces:
  - `static ExactBinary64Rational2 ExactBinary64Rational2::project(const compas_cgal::exact::Rational&)`
  - `struct SegmentAttestation2 { ExactBinary64Rational2 x0, y0, x1, y1, tool_radius, cap_chord_ratio; std::string canonical_bytes; std::string canonical_digest; }`
  - `const compas_cgal::exact::Rational& SegmentEventSource2::x0() const noexcept` and its five siblings (**return type changed**)
  - `const SegmentAttestation2& SegmentEventSource2::attestation() const noexcept`
  - unchanged: `from_binary64`, `from_exact`, `motion_data()`, `canonical_bytes()`, `canonical_digest()`

- [ ] **Step 1: Replace the `ExactBinary64Rational2` and `SegmentEventSource2` sections of `src/continuous_tea_2/segment_source.h`**

Replace lines 1-71 (everything up to and including the closing brace of `SegmentEventSource2`) with:

```cpp
#pragma once

#include "../exact_motion_2.h"
#include "exact/rational.h"
#include "partition_certificate.h"

#include <string>
#include <vector>

/// A rational in the segment lane's canonical attestation form.
///
/// DERIVED VIEW, never a carrier. Since stage 3 the segment source computes on
/// `compas_cgal::exact::Rational` and this type only renders those values for
/// the replay digest and for the Python surface. Its `exact-binary64-rational-v1`
/// framing is a frozen compatibility contract, anchored to byte literals in
/// `exact_segment_attestation_gate`.
class ExactBinary64Rational2 {
public:
    /// Project an exact rational into the segment lane's attestation form.
    ///
    /// The segment lane's single exact-to-attestation-bytes projection site and
    /// its only caller of `compas_cgal::exact::to_canonical`.
    ///
    /// Args:
    ///     value: any exact rational.
    ///
    /// Returns:
    ///     Its reduced, positive-denominator canonical form.
    ///
    /// Raises:
    ///     compas_cgal::exact::UnreducedCanonicalRationalError: if the
    ///         decomposed denominator is not positive.
    static ExactBinary64Rational2 project(
        const compas_cgal::exact::Rational& value);

    const std::string& numerator() const noexcept;
    const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    std::string text() const;

    /// The frozen `exact-binary64-rational-v1` record. Changing this encoding
    /// invalidates every stored segment replay digest.
    std::string canonical_bytes() const;

private:
    ExactBinary64Rational2(
        std::string numerator,
        std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

/// The segment source's derived attestation view.
///
/// Built once, in the source's constructor, so the frozen bytes are produced at
/// exactly one site and no consumer ever re-projects them. It is a separate
/// struct from the source because the two have different consumers and different
/// lifetimes: the source's fields feed predicates on every boundary feature and
/// every parameter cell, the attestation feeds the digest once.
///
/// It is returned BY REFERENCE from `SegmentEventSource2::attestation()`,
/// deliberately: the nanobind properties `.x0` ... `.cap_chord_ratio` use
/// `nb::rv_policy::reference_internal`, which needs a referent that outlives the
/// call. A by-value attestation would dangle.
struct SegmentAttestation2 {
    ExactBinary64Rational2 x0;
    ExactBinary64Rational2 y0;
    ExactBinary64Rational2 x1;
    ExactBinary64Rational2 y1;
    ExactBinary64Rational2 tool_radius;
    ExactBinary64Rational2 cap_chord_ratio;

    /// The frozen `segment-event-source-v1` record over the six values above.
    std::string canonical_bytes;

    /// Its SHA-256.
    std::string canonical_digest;
};

class SegmentEventSource2 {
public:
    /// Build a segment source from binary64 inputs.
    ///
    /// Raises:
    ///     NonFiniteSegmentInputError: if any input is NaN or infinite.
    ///     ZeroLengthSegmentMotionError: if the endpoints coincide.
    ///     NonPositiveToolRadiusError: if the tool radius is not positive.
    ///     InvalidCapChordRatioError: if the cap chord ratio is outside (0, 4].
    static SegmentEventSource2 from_binary64(
        double x0,
        double y0,
        double x1,
        double y1,
        double tool_radius,
        double cap_chord_ratio);

    /// Build a segment source from exact values.
    ///
    /// Raises:
    ///     ZeroLengthSegmentMotionError: if the endpoints coincide.
    ///     NonPositiveToolRadiusError: if the tool radius is not positive.
    ///     InvalidCapChordRatioError: if the cap chord ratio is outside (0, 4].
    static SegmentEventSource2 from_exact(
        const ExactSegmentMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& cap_chord_ratio);

    /// CANONICAL STATE. Every consumer computes on these. No parsing anywhere.
    const compas_cgal::exact::Rational& x0() const noexcept;
    const compas_cgal::exact::Rational& y0() const noexcept;
    const compas_cgal::exact::Rational& x1() const noexcept;
    const compas_cgal::exact::Rational& y1() const noexcept;
    const compas_cgal::exact::Rational& tool_radius() const noexcept;
    const compas_cgal::exact::Rational& cap_chord_ratio() const noexcept;

    /// DERIVED VIEW of the canonical state, projected once at construction.
    const SegmentAttestation2& attestation() const noexcept;

    /// The four attested endpoint texts, in x0, y0, x1, y1 order.
    ///
    /// An ATTESTATION projection, not a carrier read: `construct_pullback` is a
    /// Python-bound entry point whose parameters are text and which parses them
    /// back (`parameter_charts.cpp:721`). That crossing cannot be removed
    /// without changing a frozen Python signature, so it is fed from the
    /// attestation rather than from the exact carriers.
    std::vector<std::string> motion_data() const;

    /// The frozen `segment-event-source-v1` record.
    const std::string& canonical_bytes() const noexcept;

    /// Its SHA-256, embedded in every segment partition and audit record.
    const std::string& canonical_digest() const noexcept;

private:
    SegmentEventSource2(
        compas_cgal::exact::Rational x0,
        compas_cgal::exact::Rational y0,
        compas_cgal::exact::Rational x1,
        compas_cgal::exact::Rational y1,
        compas_cgal::exact::Rational tool_radius,
        compas_cgal::exact::Rational cap_chord_ratio);

    compas_cgal::exact::Rational x0_;
    compas_cgal::exact::Rational y0_;
    compas_cgal::exact::Rational x1_;
    compas_cgal::exact::Rational y1_;
    compas_cgal::exact::Rational tool_radius_;
    compas_cgal::exact::Rational cap_chord_ratio_;
    SegmentAttestation2 attestation_;
};
```

Everything from line 73 (`class FullCircleEventSource2`) to the end of the file is unchanged in this task.

Note what disappeared: `friend class SegmentEventSource2;` (the private constructor is now reached through the public `project` factory), and the two private statics `lift_binary64` / `lift_exact`.

- [ ] **Step 2: Replace lines 1-265 of `src/continuous_tea_2/segment_source.cpp`**

Replace everything from the top of the file through `SegmentEventSource2::canonical_digest()`'s closing brace (`:265`) with:

```cpp
#include "segment_source.h"

#include "../canonical_encoding.h"
#include "event_certificate.h"
#include "exact/canonical.h"
#include "sha256.h"

#include <cmath>
#include <cstdint>
#include <bit>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/number_utils.h>

namespace exact = compas_cgal::exact;

namespace {

std::string binary64_bits_text(double value)
{
    return std::to_string(std::bit_cast<std::uint64_t>(value));
}

SegmentAttestation2 project_attestation(
    const exact::Rational& x0,
    const exact::Rational& y0,
    const exact::Rational& x1,
    const exact::Rational& y1,
    const exact::Rational& tool_radius,
    const exact::Rational& cap_chord_ratio)
{
    const ExactBinary64Rational2 attested_x0 =
        ExactBinary64Rational2::project(x0);
    const ExactBinary64Rational2 attested_y0 =
        ExactBinary64Rational2::project(y0);
    const ExactBinary64Rational2 attested_x1 =
        ExactBinary64Rational2::project(x1);
    const ExactBinary64Rational2 attested_y1 =
        ExactBinary64Rational2::project(y1);
    const ExactBinary64Rational2 attested_tool_radius =
        ExactBinary64Rational2::project(tool_radius);
    const ExactBinary64Rational2 attested_cap_chord_ratio =
        ExactBinary64Rational2::project(cap_chord_ratio);
    std::string canonical = encode_string_sequence(
        {
            "segment-event-source-v1",
            attested_x0.canonical_bytes(),
            attested_y0.canonical_bytes(),
            attested_x1.canonical_bytes(),
            attested_y1.canonical_bytes(),
            attested_tool_radius.canonical_bytes(),
            attested_cap_chord_ratio.canonical_bytes(),
        });
    std::string digest = sha256_bytes(canonical);
    return {
        attested_x0,
        attested_y0,
        attested_x1,
        attested_y1,
        attested_tool_radius,
        attested_cap_chord_ratio,
        std::move(canonical),
        std::move(digest),
    };
}

} // namespace

ExactBinary64Rational2 ExactBinary64Rational2::project(
    const exact::Rational& value)
{
    const exact::CanonicalRational canonical =
        exact::to_canonical(value);
    return ExactBinary64Rational2(
        canonical.numerator(),
        canonical.denominator());
}

ExactBinary64Rational2::ExactBinary64Rational2(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string&
ExactBinary64Rational2::numerator() const noexcept
{
    return numerator_;
}

const std::string&
ExactBinary64Rational2::denominator() const noexcept
{
    return denominator_;
}

std::string ExactBinary64Rational2::text() const
{
    return denominator_ == "1"
        ? numerator_
        : numerator_ + "/" + denominator_;
}

std::string ExactBinary64Rational2::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-binary64-rational-v1",
            numerator_,
            denominator_,
        });
}

SegmentEventSource2 SegmentEventSource2::from_binary64(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_chord_ratio)
{
    // The domain checks run FIRST, in doubles, so each keeps its own named
    // error. exact::from_binary64 would raise NonFiniteBinary64Error, which is
    // not this source's contract.
    if (!std::isfinite(x0) || !std::isfinite(y0)
        || !std::isfinite(x1) || !std::isfinite(y1)
        || !std::isfinite(tool_radius)
        || !std::isfinite(cap_chord_ratio)) {
        throw NonFiniteSegmentInputError(
            "segment source values must be finite");
    }
    if (x0 == x1 && y0 == y1) {
        throw ZeroLengthSegmentMotionError(
            "segment motion endpoints must differ");
    }
    if (tool_radius <= 0.0) {
        throw NonPositiveToolRadiusError(
            "tool radius must be positive");
    }
    if (cap_chord_ratio <= 0.0
        || cap_chord_ratio > 4.0) {
        throw InvalidCapChordRatioError(
            "cap chord ratio must be in (0, 4]");
    }
    // The single double-to-exact door. A binary64 IS a dyadic rational, so this
    // is exact injection, not conversion.
    return from_exact(
        ExactSegmentMotion2{
            EPoint(
                exact::from_binary64(x0),
                exact::from_binary64(y0)),
            EPoint(
                exact::from_binary64(x1),
                exact::from_binary64(y1)),
        },
        exact::from_binary64(tool_radius),
        exact::from_binary64(cap_chord_ratio));
}

SegmentEventSource2 SegmentEventSource2::from_exact(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& cap_chord_ratio)
{
    if (motion.start == motion.end) {
        throw ZeroLengthSegmentMotionError(
            "segment motion endpoints must differ");
    }
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw NonPositiveToolRadiusError(
            "tool radius must be exact positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, Epeck::FT(4)) == CGAL::LARGER) {
        throw InvalidCapChordRatioError(
            "cap chord ratio must be exact in (0, 4]");
    }
    return SegmentEventSource2(
        motion.start.x(),
        motion.start.y(),
        motion.end.x(),
        motion.end.y(),
        tool_radius,
        cap_chord_ratio);
}

SegmentEventSource2::SegmentEventSource2(
    exact::Rational x0,
    exact::Rational y0,
    exact::Rational x1,
    exact::Rational y1,
    exact::Rational tool_radius,
    exact::Rational cap_chord_ratio)
    : x0_(std::move(x0)),
      y0_(std::move(y0)),
      x1_(std::move(x1)),
      y1_(std::move(y1)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio)),
      // Declared last, so initialised last: the six carriers above are already
      // live when the projection reads them.
      attestation_(
          project_attestation(
              x0_,
              y0_,
              x1_,
              y1_,
              tool_radius_,
              cap_chord_ratio_))
{
}

const exact::Rational&
SegmentEventSource2::x0() const noexcept
{
    return x0_;
}

const exact::Rational&
SegmentEventSource2::y0() const noexcept
{
    return y0_;
}

const exact::Rational&
SegmentEventSource2::x1() const noexcept
{
    return x1_;
}

const exact::Rational&
SegmentEventSource2::y1() const noexcept
{
    return y1_;
}

const exact::Rational&
SegmentEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const exact::Rational&
SegmentEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

const SegmentAttestation2&
SegmentEventSource2::attestation() const noexcept
{
    return attestation_;
}

std::vector<std::string> SegmentEventSource2::motion_data() const
{
    return {
        attestation_.x0.text(),
        attestation_.y0.text(),
        attestation_.x1.text(),
        attestation_.y1.text(),
    };
}

const std::string&
SegmentEventSource2::canonical_bytes() const noexcept
{
    return attestation_.canonical_bytes;
}

const std::string&
SegmentEventSource2::canonical_digest() const noexcept
{
    return attestation_.canonical_digest;
}
```

The `FullCircleEventSource2` half of the file, from the old line 267 onwards, is untouched here; Task 6 changes it.

Three things this deletes, all dead or duplicated:

- `exact_binary64(double)` (`:21-51`), the IEEE-754 bit-decomposition path. Its only caller was `lift_binary64`, whose only caller was nothing. It is a **second** double-to-exact door, which design invariant 1 forbids. Confirm with `grep -rn "exact_binary64\|lift_binary64" src tests` before deleting: the expected output is only `segment_source.h:52` and `segment_source.cpp:21,:162,:164`, all of which this task removes.
- `lift_binary64` (`:161-168`), zero callers.
- `lift_exact` (`:170-180`), whose body is `to_canonical`'s decomposition inline. `ExactBinary64Rational2::project` replaces it.

The `<CGAL/CORE/BigRat.h>` and `<CGAL/Fraction_traits.h>` includes go with them; `<CGAL/number_utils.h>` stays for `CGAL::sign` and `CGAL::compare`.

- [ ] **Step 3: Point the segment gate's two accessors at the attestation**

In `tests/native/exact_segment_attestation_gate.cpp`, replace **only** the bodies of `attested_x0` and `attested_y1`:

```cpp
const ExactBinary64Rational2& attested_x0(
    const SegmentEventSource2& source)
{
    return source.attestation().x0;
}

const ExactBinary64Rational2& attested_y1(
    const SegmentEventSource2& source)
{
    return source.attestation().y1;
}
```

Everything else in that gate — every frozen literal, every size assertion, every `require_bytes`, the digest hex — is unchanged. That is the point: the expectations were fixed before the inversion and the inversion is measured against them.

- [ ] **Step 4: Point the canonical gate's reads at the attestation**

In `tests/native/exact_canonical_gate.cpp`, `to_canonical_matches_existing_projection` now reads `source.x0()`, which is an `exact::Rational`. Replace the three surviving `require` calls' right-hand sides:

```cpp
            require(
                promoted.numerator()
                    == source.attestation().x0.numerator(),
                "promoted numerator differs from the existing projection");
            require(
                promoted.denominator()
                    == source.attestation().x0.denominator(),
                "promoted denominator differs from the existing projection");
            require(
                promoted.text() == source.attestation().x0.text(),
                "promoted text differs from the existing projection");
```

- [ ] **Step 5: Confirm the tree does not yet build, for the expected reason**

```bash
pixi run -e default exact-gates
```

Expected: **FAIL**, with errors confined to `continuous_tea_2.cpp`, `segment_projection.cpp`, `segment_fibre.cpp`, `segment_strata.cpp`, `segment_oracle.cpp` and `segment_pair_projection.cpp` — `no member named 'text' in 'CGAL::Lazy_exact_nt<...>'` and nanobind return-type errors. Any error inside `segment_source.h`, `segment_source.cpp` or either attestation gate is a defect in this task; fix it before moving on.

Do **not** commit. Continue to Task 4.

---

### Task 4: Preserve the Python surface in the binding layer

The C++ return type changed; the Python type must not. This is the task the design names as the mechanism: "the Python surface is preserved in the nanobind binding layer".

**Files:**
- Modify: `src/continuous_tea_2.cpp`

**Interfaces:**
- Consumes: `SegmentEventSource2::attestation()`.
- Produces: no new symbols. `src/compas_cgal/_continuous_tea_2.pyi` is **not** edited.

- [ ] **Step 1: Re-point the six accessor properties**

In `src/continuous_tea_2.cpp`, replace lines 307-330 (the six `def_prop_ro` blocks for `x0` … `cap_chord_ratio`) with:

```cpp
        // The source carries exact::Rational since stage 3; these properties
        // project the ATTESTATION view, so the Python type stays
        // ExactBinary64Rational2 and _continuous_tea_2.pyi does not move.
        // reference_internal is why attestation() returns a reference: it keeps
        // the source alive for as long as the returned view, and the view is a
        // member of the source rather than a temporary.
        .def_prop_ro(
            "x0",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().x0;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "y0",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().y0;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "x1",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().x1;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "y1",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().y1;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "tool_radius",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().tool_radius;
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "cap_chord_ratio",
            [](const SegmentEventSource2& source)
                -> const ExactBinary64Rational2& {
                return source.attestation().cap_chord_ratio;
            },
            nb::rv_policy::reference_internal)
```

The `from_binary64`, `motion_data`, `canonical_bytes` and `canonical_digest` bindings are unchanged: those member signatures did not move.

- [ ] **Step 2: Confirm the `.pyi` is untouched**

```bash
git diff --stat -- src/compas_cgal/_continuous_tea_2.pyi
```

Expected: empty output. If it is not, revert that file — the Python surface is frozen.

Do **not** commit. Continue to Task 5.

---

### Task 5: Convert the five consumers

Every `parse_rational(source.<accessor>().text(), ...)` becomes a materialisation of the exact representative at a declared `.exact()` boundary. This is the step that removes the 45 crossings.

**Files:**
- Modify: `src/continuous_tea_2/segment_projection.cpp`
- Modify: `src/continuous_tea_2/segment_fibre.cpp`
- Modify: `src/continuous_tea_2/segment_strata.cpp`
- Modify: `src/continuous_tea_2/segment_oracle.cpp`
- Modify: `src/continuous_tea_2/segment_pair_projection.cpp`

**Interfaces:**
- Consumes: the Task 3 accessors, all returning `const compas_cgal::exact::Rational&`; `SegmentEventSource2::attestation()`; the pre-existing file-local `Rational exact_rational(const Epeck::FT&)` at `segment_projection.cpp:77`, `segment_fibre.cpp:84`, `segment_strata.cpp:70`.
- Produces: no new exported symbols. `minimum_segment_distance_squared`, `support_can_reach_tool`, `point_numerators`, `segment_branches_at`, `segment_branch_pair_dispositions`, `construct_segment_cell_stratum`, `derive_segment_pair_projections`, `station_source`, `start_disk_has_no_material_interior`, `cap_is_exact_pi` and `support_overlap_holds` keep their exact signatures, so no header changes.

- [ ] **Step 1: `segment_projection.cpp`**

(a) `coordinate_motion` (`:268-277`) stops taking the view. Replace the whole function with:

```cpp
Polynomial coordinate_motion(
    const Epeck::FT& start,
    const Epeck::FT& end)
{
    // A DECLARED .exact() BOUNDARY (number_types.md R7): Polynomial is
    // std::vector<CORE::BigRat> and the algebraic kernel downstream needs the
    // eager representative. What stage 3 removes is the decode-and-renormalise
    // round trip, not the eager arithmetic; making this lane lazy is stage 4.
    const Rational first = exact_rational(start);
    const Rational second = exact_rational(end);
    return {first, second - first};
}
```

The two call sites at `:771` and `:773` are unchanged in text — `source.x0()` now yields `Epeck::FT`, which is what the new signature takes.

(b) `minimum_segment_distance_squared` (`:483-493`). Replace:

```cpp
    const Rational start_x =
        parse_rational(source.x0().text(), "segment x0");
    const Rational start_y =
        parse_rational(source.y0().text(), "segment y0");
    const Rational direction_x =
        parse_rational(source.x1().text(), "segment x1")
        - start_x;
    const Rational direction_y =
        parse_rational(source.y1().text(), "segment y1")
        - start_y;
```

with:

```cpp
    const Rational start_x = exact_rational(source.x0());
    const Rational start_y = exact_rational(source.y0());
    const Rational direction_x =
        exact_rational(source.x1()) - start_x;
    const Rational direction_y =
        exact_rational(source.y1()) - start_y;
```

(c) `support_can_reach_tool` (`:532-535`). Replace:

```cpp
    const Rational tool_radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
```

with:

```cpp
    const Rational tool_radius =
        exact_rational(source.tool_radius());
```

(d) The pullback loop (`:774-777`). Replace:

```cpp
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
```

with:

```cpp
    const Rational radius = exact_rational(source.tool_radius());
```

(e) `construct_pullback`'s cutter argument (`:793`). Replace `source.tool_radius().text(),` with:

```cpp
                    // ATTESTATION read, not a carrier read. construct_pullback
                    // is a Python-bound entry point whose cutter_radius
                    // parameter is text and which parses it back
                    // (parameter_charts.cpp:721); removing that crossing would
                    // change a frozen Python signature.
                    source.attestation().tool_radius.text(),
```

`source.motion_data()` on the line above is unchanged: it already projects the attestation.

(f) The two cap-crossing calls (`:877-879` and `:919-921`). Replace each:

```cpp
                    parse_rational(
                        source.cap_chord_ratio().text(),
                        "cap chord ratio")),
```

with:

```cpp
                    exact_rational(source.cap_chord_ratio())),
```

- [ ] **Step 2: `segment_fibre.cpp`**

(a) `point_numerators` (`:348-361`). Replace:

```cpp
    const Rational x0 =
        parse_rational(source.x0().text(), "segment x0");
    const Rational y0 =
        parse_rational(source.y0().text(), "segment y0");
    const Rational dx =
        parse_rational(source.x1().text(), "segment x1")
        - x0;
    const Rational dy =
        parse_rational(source.y1().text(), "segment y1")
        - y0;
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
```

with:

```cpp
    const Rational x0 = exact_rational(source.x0());
    const Rational y0 = exact_rational(source.y0());
    const Rational dx = exact_rational(source.x1()) - x0;
    const Rational dy = exact_rational(source.y1()) - y0;
    const Rational radius = exact_rational(source.tool_radius());
```

(b) The vertex-replay block (`:1097-1110`) is the identical five-declaration shape. Apply the identical replacement.

(c) `support_overlap_holds` (`:1209-1236`). Replace from `const Rational tool_radius = parse_rational(` through the closing brace of the `shifted` lambda with:

```cpp
    const Rational tool_radius =
        exact_rational(source.tool_radius());
    if (radius_squared
        != tool_radius * tool_radius) {
        return false;
    }
    // The coefficients used to be rendered as decimal text and parsed straight
    // back inside `shifted`, a round trip entirely internal to this function.
    // They are the exact values now, so the lambda takes them as numbers.
    const std::vector<Rational> x_coefficients{
        exact_rational(source.x0()),
        exact_rational(source.x1()) - exact_rational(source.x0()),
    };
    const std::vector<Rational> y_coefficients{
        exact_rational(source.y0()),
        exact_rational(source.y1()) - exact_rational(source.y0()),
    };
    const auto shifted =
        [](const std::vector<Rational>& values,
           const Rational& center) {
            std::vector<Rational> coefficients{
                values[0] - center,
                values[1],
            };
            using RationalPolynomial =
                CGAL::Polynomial<Rational>;
            using FractionTraits =
                CGAL::Fraction_traits<
                    RationalPolynomial>;
            const RationalPolynomial polynomial(
                coefficients.begin(),
                coefficients.end());
            typename FractionTraits::Numerator_type numerator;
            typename FractionTraits::Denominator_type denominator;
            typename FractionTraits::Decompose()(
                polynomial,
                numerator,
                denominator);
            return numerator;
        };
```

The two `kernel.sign_at_1_object()(shifted(...), root)` calls below are unchanged.

- [ ] **Step 3: `segment_strata.cpp`**

(a) `segment_branches_at` (`:803-819`). Replace:

```cpp
    const Rational x0 =
        parse_rational(source.x0().text());
    const Rational y0 =
        parse_rational(source.y0().text());
    const Rational x1 =
        parse_rational(source.x1().text());
    const Rational y1 =
        parse_rational(source.y1().text());
    return branches_at_station(
        records,
        exact_ft(
            x0 + parameter * (x1 - x0)),
        exact_ft(
            y0 + parameter * (y1 - y0)),
        exact_ft(
            parse_rational(
                source.tool_radius().text())));
```

with:

```cpp
    // The station is interpolated in the FILTERED carrier now: the source's
    // values are exact::Rational, `parameter` is the only eager operand, and
    // branches_at_station takes Epeck::FT -- so nothing is decoded, nothing is
    // renormalised, and no exact_ft rewrap is needed.
    const Epeck::FT parameter_ft = exact_ft(parameter);
    return branches_at_station(
        records,
        source.x0()
            + parameter_ft * (source.x1() - source.x0()),
        source.y0()
            + parameter_ft * (source.y1() - source.y0()),
        source.tool_radius());
```

(b) `segment_branch_pair_dispositions` (`:826-830`). Replace:

```cpp
    return branch_pair_dispositions_at(
        branches,
        parse_rational(
            source.cap_chord_ratio().text()));
```

with:

```cpp
    // A DECLARED .exact() BOUNDARY: branch_pair_dispositions_at still carries
    // CORE::BigRat. Converting it is stage 4.
    return branch_pair_dispositions_at(
        branches,
        exact_rational(source.cap_chord_ratio()));
```

(c) `construct_segment_cell_stratum` (`:851-857`). Replace:

```cpp
        parse_rational(
            source.cap_chord_ratio().text()));
```

with:

```cpp
        exact_rational(source.cap_chord_ratio()));
```

(d) Mark the now-dormant decoder at `:41`, keeping its body untouched:

```cpp
// DORMANT since stage 3: nothing in this translation unit decodes a number from
// text any more. Retained because decoder removal is stage 6, which requires
// explicit user permission.
[[maybe_unused]] Rational parse_rational(const std::string& text)
```

**Stage-2 dependency (D1, D4):** `construct_station_cell_stratum` further down this file is stage 2's, not stage 3's. Do not touch it. Before editing, run `grep -n "parse_rational" src/continuous_tea_2/segment_strata.cpp` and confirm the station-owned calls are already gone; if they are not, stage 2's Task 3 has not landed and this file's decoder is **not** dormant yet — in that case skip sub-step (d) and record why.

- [ ] **Step 4: `segment_oracle.cpp`**

(a) `station_source` (`:68-109` in the pre-stage-2 tree; stage 2's Task 3 Step 4 already rewrote its body — see dependency D2). Replace the whole function with:

```cpp
StationEventSource2 station_source(
    const SegmentEventSource2& source,
    const std::string& numerator,
    const std::string& denominator)
{
    // `numerator` and `denominator` are API-level parameters, supplied by the
    // binding as decimal integers (continuous_tea_2.cpp:1660-1663), not values
    // carried by the source. They are injected exactly, once, here.
    const exact::Rational parameter(
        Rational(Integer(numerator), Integer(denominator)));
    // Everything else is read straight off the source and interpolated in the
    // filtered carrier. Before stage 3 this function decoded six values out of
    // text and, before stage 2, rendered four of them straight back.
    //
    // Stage 2 needed named locals here because CORE::BigRat is a boost
    // multiprecision number whose arithmetic yields lazy EXPRESSION TEMPLATES
    // over its operands. Epeck::FT is not one: Lazy_exact_nt's operators return
    // Lazy_exact_nt by value, and the operands are references into `source`,
    // which outlives the call. Inline is safe here and only here -- if this
    // arithmetic is ever moved back onto CORE::BigRat, the named locals must
    // come back with it.
    return StationEventSource2::build(
        source.x0()
            + parameter * (source.x1() - source.x0()),
        source.y0()
            + parameter * (source.y1() - source.y0()),
        source.tool_radius(),
        source.cap_chord_ratio());
}
```

The `namespace exact = compas_cgal::exact;` alias is already present in this file's anonymous namespace — stage 2 added it. Confirm with `grep -n "namespace exact" src/continuous_tea_2/segment_oracle.cpp` and add it after `using Rational = CORE::BigRat;` only if it is missing.

**Stage-2 dependency (D3):** this body assumes `StationEventSource2::build` takes four `const exact::Rational&`. Before writing it, run `grep -n "static StationEventSource2 build" -A 5 src/continuous_tea_2/station_source.h` and confirm. If it still takes `std::string`, stage 2's Task 2 has not landed — **stop and report**; do not re-render the exact values as text to bridge the gap, because that reintroduces the crossing stage 2 exists to remove.

(b) `start_disk_has_no_material_interior` (`:201-213`). Replace:

```cpp
    const ReachKernelPoint start(
        ReachFT(parse_rational(
            source.x0().text(),
            "swept-prefix start x")),
        ReachFT(parse_rational(
            source.y0().text(),
            "swept-prefix start y")));
    const ReachFT radius(parse_rational(
        source.tool_radius().text(),
        "swept-prefix tool radius"));
```

with:

```cpp
    // A DECLARED .exact() BOUNDARY: ReachFT is Epeck_with_sqrt::FT (CORE::Expr)
    // and is built from the eager representative, exactly as
    // exact_stock_region_2.cpp:13-16 does it.
    const ReachKernelPoint start(
        ReachFT(CGAL::exact(source.x0())),
        ReachFT(CGAL::exact(source.y0())));
    const ReachFT radius(CGAL::exact(source.tool_radius()));
```

Add `#include <CGAL/number_utils.h>` to the include block if it is not already present (check with `grep -n "number_utils" src/continuous_tea_2/segment_oracle.cpp`).

(c) `cap_is_exact_pi` (`:221-228`). Replace the body with:

```cpp
{
    // chord_ratio = 4 sin^2(theta/2), hence theta = pi iff ratio = 4.
    // Exact, filtered equality: Epeck::FT is RealEmbeddable, so CGAL::compare
    // decides from the lazy interval whenever it separates and materialises the
    // exact rational only when it does not.
    return CGAL::compare(
               source.cap_chord_ratio(),
               exact::Rational(4))
        == CGAL::EQUAL;
}
```

(d) `EventTrace2`'s cap bytes (`:305-307`). Replace:

```cpp
        verified.partition.source.cap_chord_ratio()
            .canonical_bytes(),
```

with:

```cpp
        verified.partition.source.attestation()
            .cap_chord_ratio.canonical_bytes(),
```

(e) Mark the now-dormant decoder at `:33`:

```cpp
// DORMANT since stage 3: nothing in this translation unit decodes a number from
// text any more. Retained because decoder removal is stage 6, which requires
// explicit user permission.
[[maybe_unused]] Rational parse_rational(
```

- [ ] **Step 5: `segment_pair_projection.cpp`**

(a) `derive_segment_pair_projections` (`:337-341`). Replace:

```cpp
    const Rational cap_ratio =
        parse_rational(
            source.cap_chord_ratio().text(),
            "cap chord ratio");
```

with:

```cpp
    // A DECLARED .exact() BOUNDARY: this file's pair resultants are built over
    // CORE::BigRat polynomials. Converting them is stage 4.
    const Rational cap_ratio = CGAL::exact(source.cap_chord_ratio());
```

`<CGAL/number_utils.h>` is already included at `:12`.

(b) Mark the now-dormant decoder at `:23`:

```cpp
// DORMANT since stage 3: nothing in this translation unit decodes a number from
// text any more. Retained because decoder removal is stage 6, which requires
// explicit user permission.
[[maybe_unused]] Rational parse_rational(
```

- [ ] **Step 6: Build the gates**

```bash
pixi run -e default exact-gates
```

Expected: six `OK` lines, exit status 0.

`exact_segment_attestation_gate OK` here is the byte-identity result: the same frozen literals that passed against the text carrier in Task 1 now pass against the exact carrier. If `require_bytes` raises `AttestationByteDriftError`, the inversion moved the digest — **stop and report the produced bytes**; do not edit a literal.

- [ ] **Step 7: Rebuild the Python extension**

```bash
pixi run -e default _editable-rebuild
```

Expected: exit status 0 and no compiler errors. This compiles the new binding lambdas and is the check that `attestation()`'s reference lifetime satisfies `reference_internal`.

- [ ] **Step 8: Run the reference anchors**

```bash
pixi run -e default pytest -- \
  tests/adaptive/test_exact_binary64_contract.py \
  tests/adaptive/test_segment_event_substrate.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_circle_oracle.py \
  -n auto -q
```

Expected: `test_binary64_lift_is_exact`, `test_binary64_lift_is_exact_at_edges`, `test_signed_zero_lifts_to_positive_zero`, `test_binary64_lift_denominator_is_a_power_of_two`, `test_segment_source_exact_lifts_each_binary64_once` and `test_legacy_binary64_identity_survives_every_full_circle_outcome` all pass. Those six are the end-to-end confirmation that the inversion did not move the Python surface or the replay digest. If any fails, **stop and report** — the full comparison is Task 8, but one of these failing means the inversion is not byte-preserving and there is nothing to gain from continuing.

- [ ] **Step 9: Commit Tasks 3, 4 and 5 together**

```bash
git commit -m "refactor: segment source carries exact rationals, projects the text" -- \
  src/continuous_tea_2/segment_source.h \
  src/continuous_tea_2/segment_source.cpp \
  src/continuous_tea_2.cpp \
  src/continuous_tea_2/segment_projection.cpp \
  src/continuous_tea_2/segment_fibre.cpp \
  src/continuous_tea_2/segment_strata.cpp \
  src/continuous_tea_2/segment_oracle.cpp \
  src/continuous_tea_2/segment_pair_projection.cpp \
  tests/native/exact_canonical_gate.cpp
```

---

### Task 6: Route `FullCircleEventSource2` through both doors

The full-circle source already carries exact values, so there is no inversion to perform. What it does is inject doubles inline and carry a private copy of the projection. Both are door violations, and both are byte-preserving to fix because `exact::from_binary64(d)` **is** `Epeck::FT(d)` plus a finiteness check, and `to_canonical`'s decomposition **is** the private lambda's.

**Files:**
- Modify: `src/continuous_tea_2/segment_source.cpp` (the `FullCircleEventSource2` half only)

**Interfaces:**
- Consumes: `compas_cgal::exact::from_binary64`, `ExactBinary64Rational2::project`.
- Produces: no new symbols; `FullCircleEventSource2`'s header is unchanged.

- [ ] **Step 1: Route the double injection through the entry door**

In `FullCircleEventSource2::from_binary64`, replace the `ExactCircleMotion2{...}` construction and the two trailing `Epeck::FT(...)` arguments (`:306-312`) with:

```cpp
    return FullCircleEventSource2(
        ExactCircleMotion2{
            EPoint(
                exact::from_binary64(center_x),
                exact::from_binary64(center_y)),
            EVector(
                exact::from_binary64(phase_dx),
                exact::from_binary64(phase_dy)),
            clockwise,
        },
        exact::from_binary64(tool_radius),
        exact::from_binary64(cap_chord_ratio),
```

The four domain checks above it (`:276-294`) are unchanged and still run first, so `NonFiniteFullCircleInputError`, `ZeroFullCirclePhaseError`, `NonPositiveToolRadiusError` and `InvalidCapChordRatioError` all keep their exact identities. `canonical_encode_binary64` and `binary64_bits_text` are unchanged: those encode the *bit pattern*, which is a different attestation from the rational one and must stay bit-based.

- [ ] **Step 2: Route the projection through the exit door**

In `FullCircleEventSource2::from_exact`, replace the `rational_text` lambda (`:346-356`) with:

```cpp
    // The single exact-to-canonical-decimal door. This lambda used to
    // re-implement Fraction_traits::Decompose + convert_to<std::string> inline,
    // which is exactly what to_canonical does -- so this is a byte-preserving
    // deduplication, proved by `full_circle_records_do_not_move` in
    // exact_segment_attestation_gate.
    const auto rational_text = [](const Epeck::FT& value) {
        return ExactBinary64Rational2::project(value).text();
    };
```

Its seven call sites (`:359-365`, `:375-379`, `:383`) are unchanged.

!!! note "Why the projection routes through `ExactBinary64Rational2::project` and not `to_canonical` directly"

    `to_canonical` returns a `CanonicalRational`, whose `text()` is identical.
    Either call is byte-correct. Going through `project` keeps **one** projection
    site per translation unit, so a future reader grepping for
    `to_canonical` in `continuous_tea_2/` finds exactly one caller per lane —
    which is the property the design's invariant 2 is actually about. The
    full-circle framing itself stays where it is, in `from_exact`, because it is
    a fourth distinct framing (`full-circle-event-source-exact-v1`, seven fields)
    that no shared framer may speak for.

- [ ] **Step 3: Build, rebuild and run the full-circle suites**

```bash
pixi run -e default exact-gates
pixi run -e default _editable-rebuild
pixi run -e default pytest -- \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_motion_refutation.py \
  -n auto -q
```

Expected: six gate `OK` lines; `test_legacy_binary64_identity_survives_every_full_circle_outcome` passes, along with the `LEGACY_MOTION_IDENTITY_SHA256` and `LEGACY_CAP_IDENTITY_SHA256` assertions at `tests/adaptive/test_circle_oracle.py:73-74`. Those three digests are computed over the exact records this task re-routes, so they are the direct evidence. If any moves, **stop and report** — the two doors are then not byte-equivalent to the paths they replaced, which contradicts the measurement in `docs/number_types.md`.

- [ ] **Step 4: Confirm the dead second door is gone**

```bash
grep -rn "exact_binary64\|lift_binary64\|lift_exact" src tests
grep -rn "Fraction_traits" src/continuous_tea_2/segment_source.cpp
```

Expected: the first returns nothing; the second returns nothing. Every double now enters through `exact::from_binary64` and every rational leaves through `exact::to_canonical`, in this file.

- [ ] **Step 5: Commit**

```bash
git commit -m "refactor: full-circle source uses both exact doors" -- \
  src/continuous_tea_2/segment_source.cpp
```

---

### Task 7: Generic-double witness with a hard time bound

Every converted lane needs one test on generic double coordinates with a wall-clock ceiling.

!!! warning "Measured 2026-09-08 — generic coordinates alone are NOT a stressing fixture"

    On the station lane, integer 0.149 ms against generic doubles 0.208 ms is
    **1.4×**, not the 10⁴–10⁷× that `docs/number_types.md` describes. That row's
    driver is a CORE **exact-zero identity decision** forcing root-bound
    refinement, not genericity itself. Of six probed configurations only a
    **near-degenerate sweep** — tool radius comparable to half the segment
    length — moved: **7.35 ms against a 0.33 ms baseline, 22×**. Station at
    either endpoint, a collinear-ish segment and a near-zero cap ratio were all
    flat.

    So the fixture must be **degenerate-sweep geometry on generic doubles**. A
    merely non-integer fixture passes trivially and certifies nothing — the same
    failure the witness exists to prevent, one level up.

**Files:**
- Create: `tests/adaptive/test_segment_generic_double_witness.py`

**Interfaces:**
- Consumes: `compas_cgal._continuous_tea_2.construct_segment_event_partition`, `.audit_segment_tea_event_exact`, `.segment_station_cap_exceeded_exact`, `compas_cgal._stock_2.Stock2`.
- Produces: nothing importable; a regression ceiling that attributes a stall to this lane.

- [ ] **Step 1: Measure the fixture before writing the ceilings**

The ceilings must be measured numbers, not guesses. Write this probe to the scratchpad — **not** to the repository — and run it:

```python
# scratchpad only, never committed
import time

import numpy as np

from compas_cgal import _continuous_tea_2, _stock_2

# Structureless coordinates: digits of pi, e, phi and sqrt(2), used because they
# carry no rational structure, not for any mathematical property.
POCKET = np.array(
    [
        [0.3141, 0.2718, 0.0],
        [9.7183, 0.4142, 0.0],
        [9.4949, 9.6180, 0.0],
        [0.5772, 9.3010, 0.0],
    ],
    dtype=np.float64,
)
DISK = (4.7312, 5.2891, 1.3797)

# The degenerate-sweep ladder. Segment (2.7183, 3.1416) -> (5.4366, 4.1416) has
# length sqrt(2.7183**2 + 1.0**2) = 2.89638..., so half-length is 1.44819...
# The tool radius sweeps through that value, which is the configuration measured
# on 2026-09-08 to be the only one of six that moved the clock at all.
START = (2.7183, 3.1416)
END = (5.4366, 4.1416)
HALF_LENGTH = 1.4481936
RADIUS_LADDER = [
    0.25 * HALF_LENGTH,
    0.75 * HALF_LENGTH,
    0.95 * HALF_LENGTH,
    HALF_LENGTH,
    1.05 * HALF_LENGTH,
    1.5 * HALF_LENGTH,
]
CAP = 3.1416


def build_stock() -> "_stock_2.Stock2":
    stock = _stock_2.Stock2(POCKET, [])
    stock.subtract_disk(*DISK)
    return stock


def timed(label: str, call) -> None:
    for attempt in range(3):
        start = time.perf_counter()
        try:
            result = call()
        except Exception as error:  # noqa: BLE001 - probe only
            result = f"{type(error).__name__}: {error}"
        print(label, attempt, f"{time.perf_counter() - start:.3f}s", str(result)[:140])


for radius in RADIUS_LADDER:
    timed(
        f"partition r={radius:.6f}",
        lambda r=radius: _continuous_tea_2.construct_segment_event_partition(
            build_stock(), *START, *END, r, CAP
        ),
    )
    timed(
        f"audit     r={radius:.6f}",
        lambda r=radius: _continuous_tea_2.audit_segment_tea_event_exact(
            build_stock(), *START, *END, r, CAP
        ),
    )
    timed(
        f"station   r={radius:.6f}",
        lambda r=radius: _continuous_tea_2.segment_station_cap_exceeded_exact(
            build_stock(), *START, *END, 3, 7, r, CAP
        ),
    )
```

Run:

```bash
pixi run -e default python /private/tmp/claude-501/.../scratchpad/probe_segment_witness.py
```

Record four things:

1. **Which radius on the ladder is slowest** for each of the three calls. That radius is the fixture; a ladder entry that is 20× the ladder minimum is a real degenerate-sweep witness, one that is within 2× of the minimum is not and the ladder must be extended (try `0.999 * HALF_LENGTH` and `1.001 * HALF_LENGTH`, and a second segment whose endpoints are closer together).
2. **Does `construct_segment_event_partition` return a partition** rather than raising `IncompleteSegmentPartitionError`? If it raises for every ladder entry, move the disk centre until it does not, and record what worked.
3. **Does `segment_station_cap_exceeded_exact` return a bool** rather than raising `IncompleteSegmentOracleError("exact station disposition is unresolved")`? If it raises, move the station fraction (`3, 7`) and record which fraction works.
4. **The slowest of the three attempts** for each of the two chosen calls, in seconds.

Set each ceiling to `max(2.0, 4 × slowest)`, rounded up to one decimal, and write the measured value into the test as a comment. Four times the slowest observed run is wide enough that machine-to-machine variation does not flake, and tight enough that the 10⁴–10⁷× refinement stalls this lane is exposed to fail loudly.

**Honest limit to state in the commit message:** a degenerate-sweep fixture stresses the lane without reproducing the documented worst case. If the whole ladder comes back flat, that is a *finding* about this lane's reachability, not a licence to weaken the test — report it.

- [ ] **Step 2: Write the witness test**

Create `tests/adaptive/test_segment_generic_double_witness.py`, substituting the fixture Step 1 established and the two measured ceilings:

```python
"""Generic-double witness for the segment lane, with a wall-clock ceiling.

The segment lane's existing coverage runs on a 10x10 integer square with a tool
radius of 0.5 (tests/adaptive/test_segment_event_substrate.py:14-22). On fixtures
like that every square root is a perfect square, CORE's approximation error is
exactly zero, and the refinement path is never entered.

Genericity alone is not enough either: measured on 2026-09-08, generic
coordinates on this shape of problem cost only about 1.4x over integer ones. The
driver of the documented 10^4-10^7x stalls is a CORE exact-zero identity decision
forcing root-bound refinement, and of six probed configurations only a
near-degenerate sweep -- tool radius comparable to half the segment length --
moved the clock, by 22x. This fixture is therefore generic AND near-degenerate,
and it bounds the clock so a regression is attributable to this lane rather than
surfacing as a suite that quietly takes an hour.
"""

from __future__ import annotations

import math
import time

import numpy as np

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2

# No coordinate here is dyadic with a small denominator, and no two are related
# by a rational scale factor. Digits of pi, e, phi and sqrt(2) are used purely
# because they are structureless, not for any mathematical property.
GENERIC_POCKET = np.array(
    [
        [0.3141, 0.2718, 0.0],
        [9.7183, 0.4142, 0.0],
        [9.4949, 9.6180, 0.0],
        [0.5772, 9.3010, 0.0],
    ],
    dtype=np.float64,
)
GENERIC_DISK_CENTER_X = 4.7312
GENERIC_DISK_CENTER_Y = 5.2891
GENERIC_DISK_RADIUS = 1.3797

GENERIC_START_X = 2.7183
GENERIC_START_Y = 3.1416
GENERIC_END_X = 5.4366
GENERIC_END_Y = 4.1416

# The degenerate-sweep parameter. Replace with the ladder entry Step 1 measured
# as slowest; the value below is the exact half-length, which is the centre of
# the ladder, not a measurement.
GENERIC_TOOL_RADIUS = 1.4481936
GENERIC_CAP_CHORD_RATIO = 3.1416

GENERIC_STATION_NUMERATOR = 3
GENERIC_STATION_DENOMINATOR = 7

# A ceiling of 0.0 fails by construction. That is deliberate: the value must come
# from Step 1's measurement on this machine, and a test that cannot pass until it
# does is the only version of this constant that cannot be committed unmeasured.
# The rule is max(2.0, 4 x slowest of three runs), rounded up to one decimal.
# Replace both zeros and record the measured seconds in the comment beside each.
PARTITION_SECONDS_CEILING = 0.0  # measured <seconds> s, ceiling max(2.0, 4x)
STATION_SECONDS_CEILING = 0.0  # measured <seconds> s, ceiling max(2.0, 4x)


def _generic_stock() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(GENERIC_POCKET, [])
    stock.subtract_disk(
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        GENERIC_DISK_RADIUS,
    )
    return stock


def test_segment_partition_decides_on_generic_doubles_within_budget() -> None:
    stock = _generic_stock()

    start = time.perf_counter()
    partition = _continuous_tea_2.construct_segment_event_partition(
        stock,
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        GENERIC_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    )
    elapsed = time.perf_counter() - start

    # Proof the inverted carrier was actually reached: the record is emitted only
    # by SegmentEventSource2's own constructor, and its digest is the SHA-256 of
    # those bytes.
    assert b"segment-event-source-v1" in partition.source.canonical_bytes
    assert b"exact-binary64-rational-v1" in partition.source.canonical_bytes
    assert len(partition.source.canonical_digest) == 32
    assert partition.certificate.cells
    assert elapsed < PARTITION_SECONDS_CEILING, (
        f"segment partition took {elapsed:.3f}s on generic near-degenerate "
        f"doubles, ceiling {PARTITION_SECONDS_CEILING:.1f}s"
    )


def test_segment_station_decides_on_generic_doubles_within_budget() -> None:
    stock = _generic_stock()

    start = time.perf_counter()
    exceeded = _continuous_tea_2.segment_station_cap_exceeded_exact(
        stock,
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        GENERIC_STATION_NUMERATOR,
        GENERIC_STATION_DENOMINATOR,
        GENERIC_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    )
    elapsed = time.perf_counter() - start

    assert isinstance(exceeded, bool)
    assert elapsed < STATION_SECONDS_CEILING, (
        f"segment station decision took {elapsed:.3f}s on generic "
        f"near-degenerate doubles, ceiling {STATION_SECONDS_CEILING:.1f}s"
    )


def test_the_witness_fixture_is_not_an_integer_fixture() -> None:
    """Guard the guard: a later edit must not quietly round the fixture off.

    An integer or small-dyadic fixture would make this file pass while measuring
    nothing, which is exactly the failure mode it exists to prevent.
    """
    coordinates = [
        *GENERIC_POCKET[:, 0].tolist(),
        *GENERIC_POCKET[:, 1].tolist(),
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        GENERIC_DISK_RADIUS,
        GENERIC_START_X,
        GENERIC_START_Y,
        GENERIC_END_X,
        GENERIC_END_Y,
        GENERIC_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    ]
    for value in coordinates:
        assert value != round(value), f"{value} is an integer fixture"
        assert (value * 16.0) != round(value * 16.0), (
            f"{value} is a small dyadic fixture"
        )


def test_the_witness_fixture_is_a_near_degenerate_sweep() -> None:
    """Guard the guard, second half.

    Generic coordinates alone were measured at only 1.4x over integer ones on
    this shape of problem; the 22x driver was a tool radius comparable to half
    the segment length. A later edit that moves the radius away from that ratio
    would leave a fixture that is generic, fast, and certifies nothing.
    """
    half_length = 0.5 * math.hypot(
        GENERIC_END_X - GENERIC_START_X,
        GENERIC_END_Y - GENERIC_START_Y,
    )
    ratio = GENERIC_TOOL_RADIUS / half_length
    assert 0.9 < ratio < 1.1, (
        f"tool radius / half segment length is {ratio:.4f}; the witness needs a "
        "near-degenerate sweep, so this ratio must stay near 1"
    )
```

- [ ] **Step 3: Run the witness**

```bash
pixi run -e default pytest -- tests/adaptive/test_segment_generic_double_witness.py -v -n auto
```

Expected: 4 passed. If the reach assertions fail, return to Step 1 and adjust the fixture until the lane is reached — do not weaken the assertion, and do not mark anything skip or xfail.

- [ ] **Step 4: Lint and commit**

```bash
pixi run -e default lint
```

Expected: `All checks passed!`. `lint` is `ruff check src/compas_cgal tests`, so it covers this file.

```bash
git commit -m "test: generic near-degenerate segment witness with a measured bound" -- \
  tests/adaptive/test_segment_generic_double_witness.py
```

---

### Task 8: Lane equivalence and the segment-lane decoder ratchet

**Files:**
- Create: `tests/build/test_segment_lane_decoder_ratchet.py`

**Interfaces:**
- Consumes: `build/stage3/baseline.xml` (Task 1 Step 5), the eight converted production files.
- Produces: `pixi run -e default pytest -- tests/build/test_segment_lane_decoder_ratchet.py`, the mechanical guard against a text crossing coming back.

- [ ] **Step 1: Re-run the baseline suites and diff per test**

```bash
pixi run -e default pytest -- \
  tests/adaptive/test_segment_event_substrate.py \
  tests/adaptive/test_segment_event_proof_contracts.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_exact_binary64_contract.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_motion_refutation.py \
  tests/adaptive/test_motion_oracle_cache.py \
  tests/adaptive/test_transaction.py \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_event_substrate.py \
  tests/adaptive/test_generator.py \
  -n auto -q --junitxml=build/stage3/after.xml
```

Then compare per test, not by count:

```bash
pixi run -e default python - <<'PY'
import xml.etree.ElementTree as ET
from pathlib import Path


def outcomes(path: str) -> dict[str, str]:
    root = ET.parse(Path(path)).getroot()
    result = {}
    for case in root.iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        state = "passed"
        for child in case:
            if child.tag in {"failure", "error"}:
                state = child.tag
            elif child.tag == "skipped":
                state = "skipped"
        result[key] = state
    return result


before = outcomes("build/stage3/baseline.xml")
after = outcomes("build/stage3/after.xml")

changed = {k: (before.get(k), after.get(k)) for k in before | after if before.get(k) != after.get(k)}
print(f"before: {len(before)} tests, after: {len(after)} tests")
for key, (was, now) in sorted(changed.items()):
    print(f"CHANGED {key}: {was} -> {now}")
print("identical" if not changed else f"{len(changed)} tests changed outcome")
PY
```

Expected: `identical`. Anything else is a behavioural change and must be reported with the named tests — stage 3 is a representation change, so the correct outcome set is exactly the one recorded before the inversion, whatever that set was. **Do not adjust a failing test.**

- [ ] **Step 2: Run the full suite**

```bash
pixi run -e default baseline
```

Expected: completes. Record the summary line verbatim in the Step 4 commit message, so a later stage has a number to compare against.

- [ ] **Step 3: Write the segment-lane decoder ratchet**

Create `tests/build/test_segment_lane_decoder_ratchet.py`:

```python
"""Mechanically prevent a text crossing from returning to the segment lane.

Stage 3 inverted the segment source: it carries `compas_cgal::exact::Rational`
and projects the text, so nothing in the lane decodes a number out of a source
value any more. The counts below are the design's decoder ratchet, scoped to the
files stage 3 owns. They may only ever DECREASE. An increase means a consumer
started parsing again, which is the exact regression the inversion removed.

The remaining occurrences are deliberate and named:

- segment_projection.cpp and segment_fibre.cpp keep the callers that decode
  `BoundaryFeatureRecord2::primitive_coefficients`, a different lane owned by
  stages 4 and 5.
- segment_strata.cpp, segment_oracle.cpp and segment_pair_projection.cpp keep
  their definitions only; all three are dormant and marked `[[maybe_unused]]`.

Removing any definition is stage 6 and requires explicit user permission, so this
file asserts equality rather than zero.
"""

from __future__ import annotations

from pathlib import Path

REPOSITORY_ROOT = Path(__file__).parents[2]

# Derived 2026-09-08 from the pre-stage-3 tree. Before stage 3 these were, in the
# same order: 24, 26, 8, 11, 2 -- a total of 71, where segment_strata's 8 is its
# post-stage-2 count (12 today, less the 4 station-owned calls stage 2 removes).
# Stage 3 removed 45 call sites and no definitions.
EXPECTED_DECODER_OCCURRENCES = {
    "src/continuous_tea_2/segment_projection.cpp": 14,
    "src/continuous_tea_2/segment_fibre.cpp": 9,
    "src/continuous_tea_2/segment_strata.cpp": 1,
    "src/continuous_tea_2/segment_oracle.cpp": 1,
    "src/continuous_tea_2/segment_pair_projection.cpp": 1,
}

# The accessors that used to hand a consumer a string to decode. `.text()` on a
# segment carrier is gone from the lane's compute paths entirely; a reappearance
# means the attestation view is being used as a carrier again.
FORBIDDEN_CARRIER_TEXT_READS = (
    "source.x0().text()",
    "source.y0().text()",
    "source.x1().text()",
    "source.y1().text()",
    "source.tool_radius().text()",
    "source.cap_chord_ratio().text()",
)

SEGMENT_CONSUMERS = (
    "src/continuous_tea_2/segment_projection.cpp",
    "src/continuous_tea_2/segment_fibre.cpp",
    "src/continuous_tea_2/segment_strata.cpp",
    "src/continuous_tea_2/segment_oracle.cpp",
    "src/continuous_tea_2/segment_pair_projection.cpp",
)

# The one surviving text crossing, and the only place it may appear.
# `construct_pullback` is a Python-bound entry point whose motion_data and
# cutter_radius parameters are text and which parses them back; removing that
# crossing would change a frozen Python signature.
EXPECTED_ATTESTATION_TEXT_READS = {
    "src/continuous_tea_2/segment_projection.cpp": 1,
}


def _read(path: str) -> str:
    return (REPOSITORY_ROOT / path).read_text(encoding="utf-8")


def test_segment_lane_decoder_count_has_not_increased() -> None:
    observed = {
        path: _read(path).count("parse_rational")
        for path in EXPECTED_DECODER_OCCURRENCES
    }
    assert observed == EXPECTED_DECODER_OCCURRENCES


def test_no_segment_consumer_reads_a_number_out_of_text() -> None:
    for path in SEGMENT_CONSUMERS:
        source = _read(path)
        for pattern in FORBIDDEN_CARRIER_TEXT_READS:
            assert pattern not in source, f"{path} decodes a segment value from text"


def test_the_only_attestation_text_reads_are_the_frozen_api_crossing() -> None:
    for path in SEGMENT_CONSUMERS:
        observed = _read(path).count(".attestation().tool_radius.text()")
        expected = EXPECTED_ATTESTATION_TEXT_READS.get(path, 0)
        assert observed == expected, (
            f"{path} projects attestation text {observed} times, expected "
            f"{expected}; the only permitted crossing is construct_pullback's "
            "frozen cutter_radius parameter"
        )


def test_the_segment_source_projects_through_the_single_door() -> None:
    """`to_canonical` is called from exactly one place in the segment lane."""
    source = _read("src/continuous_tea_2/segment_source.cpp")
    assert source.count("to_canonical(") == 1
    assert "ExactBinary64Rational2 ExactBinary64Rational2::project(" in source


def test_there_is_one_double_to_exact_door_in_the_segment_source() -> None:
    """Both binary64 factories inject through `exact::from_binary64`."""
    source = _read("src/continuous_tea_2/segment_source.cpp")
    # Six doubles in SegmentEventSource2::from_binary64 and six in
    # FullCircleEventSource2::from_binary64. `clockwise` is a bool and carries no
    # rational, so it is not one of them.
    assert source.count("exact::from_binary64(") == 12
    assert "std::bit_cast<std::uint64_t>" in source, (
        "binary64_bits_text encodes the BIT PATTERN for the full-circle "
        "identity records; that is a different attestation from the rational "
        "one and must not be routed through the rational door"
    )
    assert "Fraction_traits" not in source, (
        "the segment source decomposes only through exact::to_canonical"
    )


def test_the_number_vocabulary_carries_no_record_tag() -> None:
    """`exact/` owns the decomposition; framing belongs to each lane."""
    canonical = _read("src/exact/canonical.cpp")
    assert "exact-binary64-rational-v1" not in canonical
    assert "encode_string_sequence" not in canonical
    assert "canonical_bytes" not in _read("src/exact/canonical.h")
```

The count `12` in `test_there_is_one_double_to_exact_door_in_the_segment_source` is six calls in `SegmentEventSource2::from_binary64` plus six in `FullCircleEventSource2::from_binary64`. **Derive it from the file after Task 6 and correct the literal if it differs** — the assertion's job is to fail if a second door appears, not to encode a number this plan counted by eye.

- [ ] **Step 4: Run the ratchet, lint and commit**

```bash
pixi run -e default pytest -- tests/build/test_segment_lane_decoder_ratchet.py -v -n auto
pixi run -e default lint
```

Expected: 6 passed, then `All checks passed!`. If `test_segment_lane_decoder_count_has_not_increased` fails, print the observed dict — either a conversion was missed or a count in this plan is stale; fix the code, not the expectation, unless the diff proves the expectation was wrong. If `segment_strata.cpp` reads 5 rather than 1, stage 2's Task 3 has not landed (dependency D1); record that and re-run after it does.

```bash
git commit -m "test: segment-lane decoder ratchet, 71 to 26 occurrences" -- \
  tests/build/test_segment_lane_decoder_ratchet.py
```

---

### Task 9: Close the stage in the documentation

Documentation is a completion artifact for the stage. `docs/number_types.md` will describe the previous state; a page that does invalidates the completion claim.

**Files:**
- Modify: `docs/number_types.md`

**Interfaces:**
- Consumes: the measured surface table at the top of this plan, the gate output from Task 5 Step 6 and Task 6 Step 3, the equivalence result from Task 8 Step 1.
- Produces: no code.

- [ ] **Step 1: Replace the status admonition**

Find the block by its heading, not by line number — stage 1 wrote `!!! note "Status after stage 1: one predicate adopted, no lane converted"` and stage 2 replaces it with a stage-2 status (dependency D5). Replace whichever is present with:

```markdown
!!! note "Status at stage 3: both event-source lanes are inverted"

    `SegmentEventSource2` carries six `exact::Rational` and projects its
    `exact-binary64-rational-v1` text through `ExactBinary64Rational2::project`,
    the lane's single caller of `to_canonical`. Its five consumers —
    `segment_projection.cpp`, `segment_fibre.cpp`, `segment_strata.cpp`,
    `segment_oracle.cpp`, `segment_pair_projection.cpp` — read the exact values
    directly; 45 of the lane's 71 `parse_rational` occurrences are gone and
    `tests/build/test_segment_lane_decoder_ratchet.py` holds the count. Both
    doors now have production callers: `from_binary64` is the only double
    injection in `segment_source.cpp`, and `to_canonical` the only
    decomposition.

    `FullCircleEventSource2` needed no inversion — it has carried
    `ExactCircleMotion2` and two `Epeck::FT` all along, and `circle_oracle.cpp`
    contains zero `parse_rational` occurrences. What stage 3 changed there is
    that it now uses the doors instead of a private copy of each.

    **What is unchanged.** The 504 unfiltered `BigRat` sites in `segment_site_*`
    are untouched (stage 4); the five converted consumers still compute in eager
    `CORE::BigRat` and materialise at a declared `.exact()` boundary, so the lane
    gained coherence, not the filter (also stage 4); coefficient-ring discipline
    is stage 5; and no `parse_rational` definition has been removed (stage 6).
    Four decoders are now dormant — `station_classifier.cpp`,
    `segment_strata.cpp`, `segment_oracle.cpp`, `segment_pair_projection.cpp` —
    each marked `[[maybe_unused]]`.
```

- [ ] **Step 2: Update the error table row for `AttestationByteDriftError`**

Replace the `AttestationByteDriftError` row and the paragraph beneath it (dependency D6 — stage 2 may already have rewritten both):

```markdown
| `AttestationByteDriftError` | a projection produces bytes differing from the frozen contract | raised by `exact_station_attestation_gate` and `exact_segment_attestation_gate`, the two lanes' contract tests |
```

```markdown
`AttestationByteDriftError` is raised only from the contract gates, exactly as
designed, and **not** as a release-build check on every projection, where the
comparison would cost more than the guarantee is worth.
`exact_segment_attestation_gate` anchors three things the repository previously
pinned nowhere: the 91-byte `exact-binary64-rational-v1` record for 0.1 and the
65-byte one for -355/113 as raw literals, the 480-byte
`segment-event-source-v1` record structurally through a helper that
re-implements the length framing independently of `encode_string_sequence`, and
that record's SHA-256 as an absolute hex constant. Each alone is escapable — a
raw literal cannot express the nested outer record, a helper that shared the
encoder would move with it, and a digest computed from the bytes moves with them
— so all three are present and are tied to each other by equality checks.
```

- [ ] **Step 3: Add the stage-3 subsection**

After the stage-2 subsection (or, if stage 2 has not landed it, after `### Measured: the carrier does not move the bytes`), insert:

```markdown
### Stage 3: the large surface, and what it did not buy

Stage 2 proved the pattern on six files. Stage 3 repeated it across eighteen, and
two things it found were not repetition at all.

| Claim | Evidence |
|---|---|
| The canonical bytes did not move | `exact_segment_attestation_gate` — 91-, 65- and 480-byte anchors plus an absolute digest hex, green before and after the inversion |
| The Python surface did not move | `git diff` on `src/compas_cgal/_continuous_tea_2.pyi` is empty; the six `.def_prop_ro` accessors project `attestation()` and still return `ExactBinary64Rational2` |
| The binary64 contract is unchanged | `tests/adaptive/test_exact_binary64_contract.py`, Hypothesis against `fractions.Fraction`, plus the pinned 0.1 decomposition in `test_segment_event_substrate.py` |
| The full-circle digests are unchanged | `LEGACY_FULL_TRACE_SHA256`, `LEGACY_MOTION_IDENTITY_SHA256` and `LEGACY_CAP_IDENTITY_SHA256` in `test_circle_oracle.py` |
| The projection matches the legacy path on non-dyadic values | eight generic reduced fractions, including a 31-digit numerator, a migrated sign, an unreduced pair and zero |
| No lane behaviour changed | per-test junit comparison of eleven suites, before and after, identical outcome set |
| The lane is exercised on generic near-degenerate doubles, with a bound | `tests/adaptive/test_segment_generic_double_witness.py` |

Four decisions are worth keeping, because the code no longer shows them.

**`FullCircleEventSource2` was never a string carrier.** The design's staging
table pairs it with `SegmentEventSource2`, but it has always held
`ExactCircleMotion2` and two `Epeck::FT`, and `circle_oracle.cpp` has zero
`parse_rational` occurrences. Its defect was different and smaller: a private
`rational_text` lambda re-implementing `to_canonical`'s decomposition, and inline
`Epeck::FT(d)` injection instead of `from_binary64`. Both were byte-preserving to
fix precisely because they were exact duplicates of the doors.

**`exact::CanonicalRational::canonical_bytes()` is gone, and no tag parameter
replaced it.** It hard-coded `exact-binary64-rational-v1` — one consumer's record
tag inside the base number vocabulary. A tag-taking method would have kept record
tags in the vocabulary and invited four frozen framings
(`exact-binary64-rational-v1`, the station lane's `exact-rational-v1`,
`boundary_events.cpp:156`'s differently-shaped `exact-rational-v1`, and
`full-circle-event-source-exact-v1`) onto one framer. Removal was the repair. The
single-door invariant is about the **decomposition**, never about the framing.
Deleting the method also severed `src/exact/`'s only include of
`continuous_tea_2/`, so the vocabulary module now compiles knowing nothing about
any event source.

**The attestation is returned by reference, against the design's sketch.** The
design writes `SegmentAttestation2 attestation() const;`. nanobind binds the six
accessors with `rv_policy::reference_internal`, which needs a referent outliving
the call, so a by-value attestation would dangle. It is a member, built once in
the constructor.

**Speed was not the goal and did not arrive.** The five converted consumers still
compute in eager `CORE::BigRat`; each converted call site is a declared `.exact()`
boundary (R7) rather than a lazy chain. What stage 3 bought is the removal of a
decode-and-renormalise round trip per value, one projection site per lane, and a
source whose every accessor returns a filtered number. The measured payoff is
predicted for stage 4, where the deep construction chains are — and stage 4 is
now measurable, because stage 3 deliberately did not mix the two changes.
```

- [ ] **Step 4: Check the six-lane table**

At `docs/number_types.md:22-33`, confirm no row still claims the `continuous_tea_2` string carriers are the segment lane's representation. If the L6 row or the `!!! danger "L4 and L6 are the unfiltered lanes"` admonition names them, amend it to say the segment and station event sources moved to L2 at stages 2 and 3, that their **consumers** are still L4/L6 until stage 4, and that `segment_site_*` is untouched.

- [ ] **Step 5: Commit**

```bash
git diff HEAD -- docs/number_types.md
```

Confirm every hunk is yours. Then:

```bash
git commit -m "docs: segment lane inverted, byte anchors and the vocabulary cleanup" -- \
  docs/number_types.md
```

---

## Stage-3 completion criteria

All of the following, or stage 3 is not done:

1. `pixi run -e default exact-gates` exits 0 with six `OK` lines, the sixth being `exact_segment_attestation_gate OK`.
2. The per-test junit comparison in Task 8 Step 1 prints `identical`.
3. `pixi run -e default baseline` completes and its summary line is recorded in the Task 8 commit message.
4. `pixi run -e default pytest -- tests/adaptive/test_segment_generic_double_witness.py -n auto` passes, with both ceilings set from measurement, the measured value written in the comment beside each, and the near-degenerate ratio guard green.
5. `pixi run -e default pytest -- tests/build/test_segment_lane_decoder_ratchet.py -n auto` passes: 26 `parse_rational` occurrences across the five segment-lane translation units, down from 71.
6. `git diff <the commit stage 3 started from>..HEAD -- src/compas_cgal/_continuous_tea_2.pyi` is empty, and `pixi run -e default pytest -- tests/adaptive/test_exact_binary64_contract.py tests/adaptive/test_segment_event_substrate.py -n auto` passes untouched. Capture the starting commit with `git rev-parse HEAD` before Task 1 (Task 1 Step 5 does this) and record it in the Task 9 commit message.
7. `grep -rn "canonical_bytes" src/exact` returns nothing, and `grep -rn "continuous_tea_2" src/exact` returns nothing.
8. `grep -rn "exact_binary64\|lift_binary64\|lift_exact" src tests` returns nothing.
9. No `parse_rational` definition was deleted, and every dormant one carries `[[maybe_unused]]` plus a comment naming stage 6.
10. `docs/number_types.md` describes stage 3, not an earlier stage, and its `AttestationByteDriftError` row names both raisers.

---

## Open questions for the lead

1. **The witness fixture is unmeasured, and the ladder may come back flat.** Task 7 Step 1 exists because I could not determine, without running the built extension, which point on the radius ladder actually drives CORE into root-bound refinement for a *segment* sweep, nor whether `construct_segment_event_partition` returns rather than raising `IncompleteSegmentPartitionError` on near-degenerate geometry. The 22× measurement on 2026-09-08 was taken on the **station** lane; whether the segment lane has the same sensitivity is an open empirical question. The plan therefore specifies a measurement, a ladder and a decision rule rather than a fixture and a number. If the whole ladder is flat, that is a finding about this lane and should come back to you rather than be worked around.

2. **`kFrozenSegmentDigestHex` is derived, not observed.** I computed the 480-byte record and its SHA-256 (`0b5f5718…`) from the encoder's specification at `event_certificate.cpp:593-602`, not from a run — no build directory was created, per the disk constraint. Task 1 Step 4 gives an explicit two-tier decision rule (structural anchor first, hex second) so a mis-derivation is correctable without weakening the contract. If you would rather the literals were observed before the plan is executed, that is a five-minute reordering of Task 1.

3. **The design's staging table pairs `FullCircleEventSource2` with `SegmentEventSource2`, and that pairing is wrong.** The full-circle source is already exact and has no text crossings; its stage-3 work is two door violations, which is a different and much smaller job. I have written Task 6 to that reality rather than to the table. Confirm you want the design document's stage-3 row amended to match, or left as the historical record.

4. **`exact_rational(const Epeck::FT&)` exists three times, file-locally.** `segment_projection.cpp:77`, `segment_fibre.cpp:84` and `segment_strata.cpp:70` each define the identical one-line `.exact()` boundary verb, and `exact_stock_region_2.cpp:13` defines a fourth variant returning `ReachFT`. Task 5 uses the existing helper where it exists and `CGAL::exact` directly in the two files that lack one, which leaves two idioms in the lane. The alternative is to promote the verb into `src/exact/` as `materialise(const Rational&) -> const CGAL::Epeck_ft&` and converge all four — one definition, one name, and it would not teach `exact/` about any event source. I did **not** plan that, because it adds to `src/exact/` in a stage whose other `exact/` work is a removal, and because retiring the three copies is an add-and-validate-then-remove decision that is yours. Say the word and it becomes a Task 2b.

5. **Deleting `exact_binary64` / `lift_binary64` is a deletion, and CLAUDE.md says never hacksaw working code.** They have zero callers (`grep` proof is Task 3 Step 2's first instruction) and `exact_binary64` is a *second* double-to-exact door, which design invariant 1 forbids. They are not decoders and not `parse_rational`, so stage 6's permission gate does not obviously cover them. I planned the removal because Task 3 rewrites the file wholesale and carrying dead duplicate doors forward would be worse — but it is one contiguous instruction you can strike without disturbing anything else.

6. **`construct_pullback`'s text crossing survives, permanently as far as this stage is concerned.** It is a Python-bound entry point taking `motion_data: Sequence[str]` and `cutter_radius: str`, which it parses back at `parameter_charts.cpp:721`. Removing that crossing means changing a frozen Python signature. Stage 3 re-points both reads at the attestation and records the remainder in the ratchet. Whether it deserves a versioned API change of its own is your call, and it is the last inbound text crossing in the segment lane that is not owned by stages 4-6.

7. **Sequencing against stage 2.** Task 5 Step 3 and Step 4 edit `segment_strata.cpp` and `segment_oracle.cpp`, both of which stage 2's Task 3 also edits — different functions in the first, the *same* function in the second. Task 5 must not start before stage 2's Task 3 is committed, or the two sessions will write conflicting versions of `segment_oracle.cpp::station_source`. Tasks 1, 2, 3, 4, 6 and 7 are independent of stage 2 and can run at any time.
