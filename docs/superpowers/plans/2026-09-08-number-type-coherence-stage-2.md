# Number-Type Coherence — Stage 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Invert the canonical direction for `StationEventSource2` — the smallest of the string-carrier sources — so the source stores `compas_cgal::exact::Rational` and *projects* the `exact-rational-v1` text, with the `station-event-source-v1` bytes proven unmoved by an absolute literal that does not exist today.

**Architecture:** `StationEventSource2` stops carrying `numerator_`/`denominator_` strings and starts carrying four `exact::Rational` (= `Epeck::FT`, lazy and interval-filtered). `ExactRational2` survives as a *derived attestation view* built by a new `ExactRational2::project`, which is the station lane's single caller of `exact::to_canonical`. `StationAttestation2` groups the four projected values plus the frozen outer record, is built once in the source's constructor, and is the lane's only projection site. The four consumers — `station_classifier.cpp`, `segment_strata.cpp`, `circle_strata.cpp`, `segment_oracle.cpp` — stop decoding text and read the exact numbers directly, deleting 13 inbound text crossings. No `parse_rational` *definition* is removed; that is stage 6 and needs explicit user permission.

**Tech Stack:** C++20, CGAL 6.0.1 (vendored `external/cgal`), Boost multiprecision (vendored `external/boost`, GMP disabled), nanobind, CMake + scikit-build-core, pixi, pytest with Hypothesis.

**Spec:** `docs/superpowers/specs/2026-09-08-number-type-coherence-design.md`, stage 2 row.

## Measured surface (2026-09-08, this worktree)

`StationEventSource2` and `ExactRational2` are **not bound to Python**. `grep -rn "class_<" src/ | grep -i "station\|exactrational"` returns nothing, and neither name appears in `src/compas_cgal/_continuous_tea_2.pyi`. The frozen Python API is therefore preserved by *not being touched*, and no binding-layer projection shim is needed. What is Python-visible, and load-bearing, is the **digest**: `StationEventSource2::canonical_bytes()` is embedded verbatim in every `full-circle-cell-decision-v1` record (`circle_strata.cpp:674`), which feeds `FullCircleCellAuthority2::canonical_digest` and reaches Python as `EventTrace2.canonical_digest`.

| File | Role | `parse_rational` occurrences today | after stage 2 |
|---|---|---|---|
| `src/continuous_tea_2/station_source.h` | carrier + view | — | — |
| `src/continuous_tea_2/station_source.cpp` | carrier + view | 4 (1 definition, 3 calls) | 2 (definition + the retained legacy decoder) |
| `src/continuous_tea_2/station_classifier.cpp` | consumer | 4 (1 definition, 3 calls) | 1 (definition only, dormant until stage 6) |
| `src/continuous_tea_2/segment_strata.cpp` | consumer | 12 (1 definition, 11 calls; 4 station-owned) | 8 |
| `src/continuous_tea_2/circle_strata.cpp` | producer | 6 (1 definition, 5 calls; 4 station-owned) | 2 |
| `src/continuous_tea_2/segment_oracle.cpp` | producer | 11 (1 definition, 10 calls; 0 station-owned) | 11 |
| **total** | **6 files** | **37** | **24** |

13 call sites removed, 0 definitions removed. Of the 13 private decoder copies counted repo-wide in the design, stage 2 makes exactly **one** dormant (`station_classifier.cpp:17`); the other station-lane definitions keep non-station callers. Stage 2 also removes 4 `exact_rational_text` and 2 `rational_text` *encode* crossings in `circle_strata.cpp:146-177` and the 4-use `text` lambda in `segment_oracle.cpp:92-103`.

Python entry points that reach the station lane:

| Entry point | Route | Proof it reached the lane |
|---|---|---|
| `audit_full_circle_tea_event_exact` | `circle_oracle.cpp:811` → `construct_full_circle_cell_authority` → `station_source()` → `classify_station_cell` | `b"station-event-source-v1" in trace.canonical_bytes` |
| `audit_segment_tea_event_exact` | `segment_oracle.cpp:506` → `classify_station_cell` | verdict changes with the station decision |
| `segment_station_cap_exceeded_exact` | `segment_oracle.cpp:647` → `classify_station_cell` | the return value *is* the station decision |

`full_circle_rational_probe_exceeds_cap_exact` never reaches the station lane — it throws `"exact rational probe requires the Task 5 full-circle pullback substrate"` at `circle_oracle.cpp:966`. Do not use it as a witness.

## Global Constraints

- Exact arithmetic is settled; do not introduce any epsilon, tolerance, deflation factor or `nextafter` into a decision path.
- Named exceptions only, convention `<Domain><Condition>Error`. Never `throw std::runtime_error("...")` directly.
- Stage 2 **consumes** `src/exact/`; it adds nothing to it. `exact::CanonicalRational::canonical_bytes()` hard-codes the `exact-binary64-rational-v1` tag and is therefore **not** usable by the station lane, which frames as `exact-rational-v1`. The station framing stays in `ExactRational2::canonical_bytes()`, fed by `to_canonical`'s numerator/denominator. Do not add a tag parameter to `exact/`.
- Do **not** delete any `parse_rational` definition, `ExactRational2::build(text)`, or any other decoder. Removal is stage 6 and requires explicit user permission. Dormant decoders get `[[maybe_unused]]` and a comment naming stage 6.
- No `pytest.mark.skip`, `skipif`, or `xfail` under any circumstance. A failing test fails.
- Never modify an existing reference test to make new code pass. In particular `LEGACY_FULL_TRACE_SHA256` in `tests/adaptive/test_circle_oracle.py:27-34` is a reference anchor: if stage 2 moves it, stage 2 is wrong.
- Native gates compile under `CMAKE_BUILD_TYPE Release`, which defines `NDEBUG` and erases `<cassert>`. Use the throwing `require` helper from `tests/native/exact_canonical_gate.cpp:18-31`; never `assert`.
- `CGAL_DISABLE_GMP` and `CGAL_USE_BOOST_MP` are set repo-wide. `exact::Rational` is `Epeck::FT` = `Lazy_exact_nt<Epeck_ft>`, and `Epeck_ft` **is** `CORE::BigRat` — `station_classifier.cpp:47-48` static-asserts exactly that. Every `CORE::BigRat` → `exact::Rational` conversion in this plan is therefore an exact rewrap, not a reinterpretation.
- `${CMAKE_CURRENT_SOURCE_DIR}/src` is a PUBLIC include directory on `continuous_tea_exact_core` (`CMakeLists.txt:342-348`), and `_continuous_tea_2` links that target (`CMakeLists.txt:461`), so `#include "exact/rational.h"` resolves from `src/continuous_tea_2/` in both the static library and the nanobind module.
- Adding a `.cpp` file or a target requires editing `CMakeLists.txt` and a full rebuild.
- Run pytest through the pixi task, arguments after `--`: `pixi run -e default pytest -- <paths> -n auto`.
- **Shared worktree.** `CMakeLists.txt` and `pyproject.toml` already carry other sessions' uncommitted edits. Commit by pathspec with the message *before* the `--`: `git commit -m "msg" -- <paths>`. Before every commit that includes `CMakeLists.txt` or `pyproject.toml`, run `git diff HEAD -- CMakeLists.txt pyproject.toml` and confirm every hunk is yours; a pathspec isolates files, not hunks. If a foreign hunk is present, stop and report rather than committing it.
- Commit messages: extremely concise, lowercase, no attribution trailers.

---

### Task 1: Anchor the station bytes and capture the lane baseline, before touching production code

Nothing in the repository pins the `exact-rational-v1` or `station-event-source-v1` framings to an absolute constant. `LEGACY_FULL_TRACE_SHA256` detects a change but cannot localise one, and both framings share the single `encode_string_sequence`, so a differential test alone would move both sides together and stay green while every stored digest stopped verifying. This task lands the absolute anchors and the differential, against the **current** text-carrier code, plus the pre-inversion behavioural baseline that Task 5 compares to. No production file changes.

Note also that the tag `exact-rational-v1` is used with a **second, different framing** at `boundary_events.cpp:156` (`tagged_record("exact-rational-v1", {stream.str()})`, one field, tag-plus-NUL prefix). That is pre-existing and out of scope — this gate anchors the station lane's three-field `encode_string_sequence` framing specifically, and the comment says so.

**Files:**
- Create: `tests/native/exact_station_attestation_gate.cpp`
- Modify: `CMakeLists.txt` (add the `exact_station_attestation_gate` executable)
- Modify: `pyproject.toml` (add the gate to the `exact-gates` task)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational`, `compas_cgal::exact::to_canonical`, `compas_cgal::exact::CanonicalRational`, `compas_cgal::exact::AttestationByteDriftError` (stage 0); `ExactRational2::build`, `StationEventSource2::build` (current text API).
- Produces: `build/exact-gate/exact_station_attestation_gate`, the absolute byte anchor for the station lane, and `build/stage2/baseline.xml`, the pre-inversion per-test outcome record.

- [ ] **Step 1: Write the gate**

Create `tests/native/exact_station_attestation_gate.cpp`:

```cpp
#include "exact/canonical.h"
#include "exact/errors.h"
#include "exact/rational.h"

#include "continuous_tea_2/station_source.h"

#include <CGAL/CORE/BigRat.h>

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
/// gate error. `AttestationByteDriftError` was declared at stage 0 with no
/// raiser, against the day a projection site existed; stage 2 builds that site,
/// and this is the contract test the design says must raise it.
void require_bytes(
    const std::string& produced,
    const std::string& frozen,
    const char* message)
{
    if (produced != frozen) {
        throw exact::AttestationByteDriftError(message);
    }
}

/// The frozen `exact-rational-v1` record for the rational -355/113.
///
/// FROZEN COMPATIBILITY CONTRACT. These 56 bytes are the on-disk shape every
/// stored full-circle replay digest was computed over: an 8-byte big-endian
/// field count, then per field an 8-byte big-endian length followed by its
/// payload. -355/113 is chosen because it is neither dyadic nor an integer, so
/// it exercises a real numerator and denominator, and because its sign lives in
/// the numerator, which is what canonicalisation must guarantee.
///
/// This is the station lane's THREE-field `encode_string_sequence` framing. A
/// second, unrelated framing shares the same tag at boundary_events.cpp:156
/// (one field, tag-plus-NUL prefix, via `tagged_record`); this literal does not
/// speak for that one.
///
/// The literal exists because the differential checks below CANNOT see a change
/// here: legacy and projected paths share the one `encode_string_sequence`, so a
/// single edit to that framing moves both sides identically and leaves the
/// differential green while every stored digest silently stops verifying.
///
/// If a change makes this check fail, that change invalidates every stored
/// full-circle replay digest. It must be a deliberate, versioned decision: bump
/// the "exact-rational-v1" tag and migrate the stored digests. Editing this
/// literal to restore green is the one repair that is always wrong.
constexpr char kFrozenStationValueBytes[] =
    "\x00\x00\x00\x00\x00\x00\x00\x03"                          // 3 fields
    "\x00\x00\x00\x00\x00\x00\x00\x11" "exact-rational-v1"      // 17 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x04" "-355"                   // 4 bytes
    "\x00\x00\x00\x00\x00\x00\x00\x03" "113";                   // 3 bytes

static_assert(
    sizeof(kFrozenStationValueBytes) - 1 == 56,
    "the frozen value literal must be 56 bytes: 8 + (8+17) + (8+4) + (8+3)");

/// An independent re-implementation of the canonical length framing.
///
/// Deliberately NOT `encode_string_sequence`: the outer
/// `station-event-source-v1` record nests four inner records and would be
/// unreadable as a raw literal, so it is anchored structurally instead. The
/// equality check in `frozen_framing_agrees_with_the_literal` ties this helper
/// to the raw literal above, so together they pin the framing bytes, the field
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
    return frozen_sequence({"exact-rational-v1", numerator, denominator});
}

/// The probe station: centre (-355/113, 22/7), tool radius 1/2, cap 7/2.
///
/// Every component is a reduced non-dyadic fraction where one exists, the
/// radius is positive and the cap lies in (0, 4], so the source's own
/// validation admits it. TASK 2 REWRITES ONLY THIS FUNCTION BODY, to build the
/// same station through the exact API. The frozen expectations below never move.
StationEventSource2 probe_station()
{
    return StationEventSource2::build("-355/113", "22/7", "1/2", "7/2");
}

std::string frozen_station_bytes()
{
    return frozen_sequence({
        "station-event-source-v1",
        frozen_value("-355", "113"),
        frozen_value("22", "7"),
        frozen_value("1", "2"),
        frozen_value("7", "2"),
    });
}

/// Render a CORE::BigRat the way the legacy decoder expects to read it back.
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

/// Tie the structural helper to the raw literal, so neither anchor stands alone.
void frozen_framing_agrees_with_the_literal()
{
    const std::string literal(
        kFrozenStationValueBytes, sizeof(kFrozenStationValueBytes) - 1);
    require_bytes(
        frozen_value("-355", "113"),
        literal,
        "the gate's independent framing helper disagrees with the frozen "
        "literal; one of the two has been edited");
}

/// The `exact-rational-v1` record of a single attested value must not move.
void value_record_matches_the_frozen_literal()
{
    const std::string literal(
        kFrozenStationValueBytes, sizeof(kFrozenStationValueBytes) - 1);
    const ExactRational2 attested =
        ExactRational2::build("-355/113");
    require_bytes(
        attested.canonical_bytes(),
        literal,
        "exact-rational-v1 bytes differ from the frozen literal: every stored "
        "full-circle replay digest is invalidated");
    require(
        attested.numerator() == "-355",
        "attested numerator is not -355");
    require(
        attested.denominator() == "113",
        "attested denominator is not 113");
    require(
        attested.text() == "-355/113",
        "attested text is not -355/113");
}

/// The outer `station-event-source-v1` record must not move either: its tag,
/// its field count and the order of its four values are all frozen.
void station_record_matches_the_frozen_framing()
{
    const std::string frozen = frozen_station_bytes();
    require(
        frozen.size() == 281,
        "the frozen station record must be 281 bytes: "
        "8 + (8+23) + (8+56) + (8+52) + (8+51) + (8+51)");
    require_bytes(
        probe_station().canonical_bytes(),
        frozen,
        "station-event-source-v1 bytes differ from the frozen framing: every "
        "stored full-circle replay digest is invalidated");
}

/// The projection door must agree with the legacy text decoder on values the
/// station lane actually carries: arbitrary reduced rationals, not just the
/// dyadic ones a binary64 produces. `exact_canonical_gate` covers the dyadic
/// case; this covers the case that gate cannot reach.
void projection_matches_the_legacy_decoder_on_generic_rationals()
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
        const Rational value(Integer(numerator), Integer(denominator));
        const exact::CanonicalRational projected =
            exact::to_canonical(exact::Rational(value));
        const ExactRational2 legacy =
            ExactRational2::build(legacy_text(value));
        require(
            projected.numerator() == legacy.numerator(),
            "projected numerator differs from the legacy decoder");
        require(
            projected.denominator() == legacy.denominator(),
            "projected denominator differs from the legacy decoder");
        require(
            projected.text() == legacy.text(),
            "projected text differs from the legacy decoder");
    }
}

/// The producers hand the station a value built by a CHAIN of exact arithmetic,
/// not a leaf. The claim stage 2 rests on is that the lazy carrier's exact
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
    const ExactRational2 legacy = ExactRational2::build(legacy_text(eager));
    require(
        projected.numerator() == legacy.numerator(),
        "a depth-16 lazy chain canonicalises to a different numerator than the "
        "same chain computed eagerly");
    require(
        projected.denominator() == legacy.denominator(),
        "a depth-16 lazy chain canonicalises to a different denominator than "
        "the same chain computed eagerly");
}

/// The source's own range checks are part of its contract and must survive the
/// signature change unchanged.
void invalid_stations_are_rejected()
{
    bool raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "0", "1");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero tool radius was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "1", "0");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero cap chord ratio was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build("0", "0", "1", "17/4");
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a cap chord ratio above 4 was accepted");

    // The closed upper end is admitted, and so is a negative coordinate.
    (void)StationEventSource2::build("-1/3", "-7/2", "1", "4");
}

}  // namespace

int main()
{
    try {
        frozen_framing_agrees_with_the_literal();
        value_record_matches_the_frozen_literal();
        station_record_matches_the_frozen_framing();
        projection_matches_the_legacy_decoder_on_generic_rationals();
        lazy_and_eager_chains_canonicalise_identically();
        invalid_stations_are_rejected();
    } catch (const std::exception& error) {
        std::printf(
            "exact_station_attestation_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("exact_station_attestation_gate OK\n");
    return 0;
}
```

- [ ] **Step 2: Wire the gate into CMake**

In `CMakeLists.txt`, directly after the `exact_one_root_gate` block (which ends at line 254), add:

```cmake
add_executable(exact_station_attestation_gate EXCLUDE_FROM_ALL
    tests/native/exact_station_attestation_gate.cpp
)
target_link_libraries(
    exact_station_attestation_gate PRIVATE continuous_tea_exact_core)
target_compile_definitions(
    exact_station_attestation_gate PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 3: Wire the gate into the `exact-gates` task**

In `pyproject.toml`, replace the two lines at `:286-287` so the new gate is built and run with the others:

```toml
_exact-gates-build = "cmake --build build/exact-gate --target exact_vocabulary_gate exact_canonical_gate exact_one_root_gate exact_station_attestation_gate"
_exact-gates-run = "build/exact-gate/exact_vocabulary_gate && build/exact-gate/exact_canonical_gate && build/exact-gate/exact_one_root_gate && build/exact-gate/exact_station_attestation_gate"
```

- [ ] **Step 4: Run the gate against the current, un-inverted code**

Run:
```bash
pixi run -e default exact-gates
```
Expected: four `OK` lines ending with `exact_station_attestation_gate OK`, exit status 0.

If `value_record_matches_the_frozen_literal` or `station_record_matches_the_frozen_framing` fails, the frozen constants in this plan are wrong and must be corrected from the *observed* bytes — print them with `xxd` and recompute — because the current bytes are by definition the contract. **Stop and report the observed bytes before changing anything else.**

If `projection_matches_the_legacy_decoder_on_generic_rationals` or `lazy_and_eager_chains_canonicalise_identically` fails, **stop and report the differing pair**. That is the premise stage 2 rests on, and a failure means the inversion cannot be byte-preserving as designed.

- [ ] **Step 5: Capture the pre-inversion behavioural baseline**

Stage 0's completion criterion 2 was formally unmet because no branch-point baseline existed. Do not repeat that: capture the baseline **now**, from this tree, before any production file changes.

Run:
```bash
mkdir -p build/stage2
git rev-parse HEAD | tee build/stage2/start-commit.txt
pixi run -e default pytest -- \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_motion_refutation.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_generator.py \
  tests/adaptive/test_event_substrate.py \
  tests/adaptive/test_segment_event_substrate.py \
  -n auto -q --junitxml=build/stage2/baseline.xml
```
Expected: a run that completes and writes `build/stage2/baseline.xml`. The pass/fail counts do **not** need to be all-green — this tree has pre-existing failures. What matters is that the file exists and records a per-test outcome for every test in these suites. Record the summary line verbatim in the Task 5 comparison.

`build/` is not tracked, so nothing here is committed.

- [ ] **Step 6: Commit**

```bash
git diff HEAD -- CMakeLists.txt pyproject.toml
```
Confirm every hunk is yours (see Global Constraints). Then:

```bash
git commit -m "test: anchor station attestation bytes, generic-rational projection gate" -- \
  tests/native/exact_station_attestation_gate.cpp CMakeLists.txt pyproject.toml
```

---

### Task 2: Invert `StationEventSource2`

The source stops carrying text and starts carrying `exact::Rational`. `ExactRational2` stays, demoted to a derived view with a new `project` factory; its text decoder is retained untouched for the gate's differential and for stage 6.

**This task leaves the tree non-compiling** — the four consumers still pass strings and read `ExactRational2`. Task 3 converts them, and Tasks 2 and 3 share one commit, made at the end of Task 3. Do not commit at the end of this task.

**Files:**
- Modify: `src/continuous_tea_2/station_source.h`
- Modify: `src/continuous_tea_2/station_source.cpp`
- Modify: `tests/native/exact_station_attestation_gate.cpp` (the `probe_station` body only)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational` (`src/exact/rational.h`), `compas_cgal::exact::to_canonical` and `compas_cgal::exact::CanonicalRational` (`src/exact/canonical.h`).
- Produces:
  - `ExactRational2 ExactRational2::project(const compas_cgal::exact::Rational&)`
  - `struct StationAttestation2 { ExactRational2 center_x, center_y, tool_radius, cap_chord_ratio; std::string canonical_bytes; }`
  - `StationEventSource2 StationEventSource2::build(const compas_cgal::exact::Rational& center_x, const compas_cgal::exact::Rational& center_y, const compas_cgal::exact::Rational& tool_radius, const compas_cgal::exact::Rational& cap_chord_ratio)`
  - `const compas_cgal::exact::Rational& StationEventSource2::center_x() const noexcept` and the three siblings (**return type changed**)
  - `const StationAttestation2& StationEventSource2::attestation() const noexcept`
  - `const std::string& StationEventSource2::canonical_bytes() const noexcept` (**unchanged signature, unchanged bytes**)

- [ ] **Step 1: Replace `src/continuous_tea_2/station_source.h`**

```cpp
#pragma once

#include "exact/rational.h"
#include "partition_certificate.h"

#include <string>

/// A rational in the station lane's canonical attestation form.
///
/// DERIVED VIEW, never a carrier. Since stage 2 the station source computes on
/// `compas_cgal::exact::Rational` and this type only renders those values for
/// the replay digest. Its `exact-rational-v1` framing is a frozen compatibility
/// contract, anchored to a byte literal in `exact_station_attestation_gate`.
class ExactRational2 {
public:
    /// Project an exact rational into the station lane's attestation form.
    ///
    /// The station lane's single exact-to-attestation-bytes projection site and
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
    static ExactRational2 project(
        const compas_cgal::exact::Rational& value);

    /// Decode decimal text into attestation form.
    ///
    /// LEGACY DECODER. Nothing in the pipeline calls it since stage 2; it is
    /// retained so `exact_station_attestation_gate` can keep proving `project`
    /// agrees with it byte for byte. Removing it is stage 6 and requires
    /// explicit user permission.
    ///
    /// Raises:
    ///     InvalidStationSourceError: if `text` is not one exact rational, or
    ///         has a zero denominator.
    static ExactRational2 build(
        const std::string& text);

    const std::string& numerator() const noexcept;
    const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    std::string text() const;

    /// The frozen `exact-rational-v1` record. Changing this encoding
    /// invalidates every stored full-circle replay digest.
    std::string canonical_bytes() const;

private:
    ExactRational2(
        std::string numerator,
        std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

/// The station source's derived attestation view.
///
/// Built once, in the source's constructor, so the frozen bytes are produced at
/// exactly one site and no consumer ever re-projects them. It is a separate
/// struct from the source because the two have different consumers and
/// different lifetimes: the source's fields feed predicates on every cell, the
/// attestation feeds the digest once.
struct StationAttestation2 {
    ExactRational2 center_x;
    ExactRational2 center_y;
    ExactRational2 tool_radius;
    ExactRational2 cap_chord_ratio;

    /// The frozen `station-event-source-v1` record over the four values above.
    std::string canonical_bytes;
};

class StationEventSource2 {
public:
    /// Build a station from exact values.
    ///
    /// Args:
    ///     center_x: the station centre's x coordinate.
    ///     center_y: the station centre's y coordinate.
    ///     tool_radius: the cutter radius; must be strictly positive.
    ///     cap_chord_ratio: the engagement cap surrogate; must lie in (0, 4].
    ///
    /// Raises:
    ///     InvalidStationSourceError: if either range constraint is violated.
    static StationEventSource2 build(
        const compas_cgal::exact::Rational& center_x,
        const compas_cgal::exact::Rational& center_y,
        const compas_cgal::exact::Rational& tool_radius,
        const compas_cgal::exact::Rational& cap_chord_ratio);

    /// CANONICAL STATE. Every consumer computes on these. No parsing anywhere.
    const compas_cgal::exact::Rational& center_x() const noexcept;
    const compas_cgal::exact::Rational& center_y() const noexcept;
    const compas_cgal::exact::Rational& tool_radius() const noexcept;
    const compas_cgal::exact::Rational& cap_chord_ratio() const noexcept;

    /// DERIVED VIEW of the canonical state, projected once at construction.
    const StationAttestation2& attestation() const noexcept;

    /// The frozen `station-event-source-v1` record.
    ///
    /// Surfaced on the source, not only on the attestation, because
    /// `circle_strata.cpp` embeds it in every full-circle cell decision record
    /// and does so once per parameter cell. It is the attestation's own bytes.
    const std::string& canonical_bytes() const noexcept;

private:
    StationEventSource2(
        compas_cgal::exact::Rational center_x,
        compas_cgal::exact::Rational center_y,
        compas_cgal::exact::Rational tool_radius,
        compas_cgal::exact::Rational cap_chord_ratio);

    compas_cgal::exact::Rational center_x_;
    compas_cgal::exact::Rational center_y_;
    compas_cgal::exact::Rational tool_radius_;
    compas_cgal::exact::Rational cap_chord_ratio_;
    StationAttestation2 attestation_;
};

class InvalidStationSourceError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
```

- [ ] **Step 2: Replace `src/continuous_tea_2/station_source.cpp`**

```cpp
#include "station_source.h"

#include "event_certificate.h"
#include "exact/canonical.h"

#include <string_view>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/number_utils.h>

namespace exact = compas_cgal::exact;

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

// LEGACY DECODER, reachable only through ExactRational2::build, which nothing in
// the pipeline calls since stage 2. Retained so the attestation gate can keep
// comparing the projection against it. Removing both is stage 6.
Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw InvalidStationSourceError(
                std::string(role)
                + " is not one exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw InvalidStationSourceError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw InvalidStationSourceError(
            std::string(role)
            + " is not one exact rational");
    }
}

StationAttestation2 project_attestation(
    const exact::Rational& center_x,
    const exact::Rational& center_y,
    const exact::Rational& tool_radius,
    const exact::Rational& cap_chord_ratio)
{
    const ExactRational2 attested_center_x =
        ExactRational2::project(center_x);
    const ExactRational2 attested_center_y =
        ExactRational2::project(center_y);
    const ExactRational2 attested_tool_radius =
        ExactRational2::project(tool_radius);
    const ExactRational2 attested_cap_chord_ratio =
        ExactRational2::project(cap_chord_ratio);
    return {
        attested_center_x,
        attested_center_y,
        attested_tool_radius,
        attested_cap_chord_ratio,
        encode_string_sequence(
            {
                "station-event-source-v1",
                attested_center_x.canonical_bytes(),
                attested_center_y.canonical_bytes(),
                attested_tool_radius.canonical_bytes(),
                attested_cap_chord_ratio.canonical_bytes(),
            }),
    };
}

} // namespace

ExactRational2 ExactRational2::project(
    const exact::Rational& value)
{
    const exact::CanonicalRational canonical =
        exact::to_canonical(value);
    return ExactRational2(
        canonical.numerator(),
        canonical.denominator());
}

ExactRational2 ExactRational2::build(
    const std::string& text)
{
    const Rational value =
        parse_rational(text, "station value");
    return ExactRational2(
        CORE::numerator(value)
            .convert_to<std::string>(),
        CORE::denominator(value)
            .convert_to<std::string>());
}

ExactRational2::ExactRational2(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string&
ExactRational2::numerator() const noexcept
{
    return numerator_;
}

const std::string&
ExactRational2::denominator() const noexcept
{
    return denominator_;
}

std::string ExactRational2::text() const
{
    return denominator_ == "1"
        ? numerator_
        : numerator_ + "/" + denominator_;
}

std::string ExactRational2::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-rational-v1",
            numerator_,
            denominator_,
        });
}

StationEventSource2 StationEventSource2::build(
    const exact::Rational& center_x,
    const exact::Rational& center_y,
    const exact::Rational& tool_radius,
    const exact::Rational& cap_chord_ratio)
{
    // Exact, filtered decisions on the carrier itself. Epeck::FT is
    // RealEmbeddable, so CGAL::sign and CGAL::compare decide from the lazy
    // interval whenever it separates and materialise the exact rational only
    // when it does not -- where the text carrier forced a decode and a full
    // bignum comparison every time.
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw InvalidStationSourceError(
            "station tool radius must be positive");
    }
    if (CGAL::sign(cap_chord_ratio) != CGAL::POSITIVE
        || CGAL::compare(cap_chord_ratio, exact::Rational(4))
            == CGAL::LARGER) {
        throw InvalidStationSourceError(
            "station cap chord ratio must lie in (0, 4]");
    }
    return StationEventSource2(
        center_x,
        center_y,
        tool_radius,
        cap_chord_ratio);
}

StationEventSource2::StationEventSource2(
    exact::Rational center_x,
    exact::Rational center_y,
    exact::Rational tool_radius,
    exact::Rational cap_chord_ratio)
    : center_x_(std::move(center_x)),
      center_y_(std::move(center_y)),
      tool_radius_(std::move(tool_radius)),
      cap_chord_ratio_(std::move(cap_chord_ratio)),
      // Declared last, so initialised last: the four carriers above are already
      // live when the projection reads them.
      attestation_(
          project_attestation(
              center_x_,
              center_y_,
              tool_radius_,
              cap_chord_ratio_))
{
}

const exact::Rational&
StationEventSource2::center_x() const noexcept
{
    return center_x_;
}

const exact::Rational&
StationEventSource2::center_y() const noexcept
{
    return center_y_;
}

const exact::Rational&
StationEventSource2::tool_radius() const noexcept
{
    return tool_radius_;
}

const exact::Rational&
StationEventSource2::cap_chord_ratio() const noexcept
{
    return cap_chord_ratio_;
}

const StationAttestation2&
StationEventSource2::attestation() const noexcept
{
    return attestation_;
}

const std::string&
StationEventSource2::canonical_bytes() const noexcept
{
    return attestation_.canonical_bytes;
}
```

- [ ] **Step 3: Point the gate's probe at the exact API**

In `tests/native/exact_station_attestation_gate.cpp`, replace **only** the body of `probe_station`:

```cpp
StationEventSource2 probe_station()
{
    return StationEventSource2::build(
        exact::Rational(-355) / exact::Rational(113),
        exact::Rational(22) / exact::Rational(7),
        exact::Rational(1) / exact::Rational(2),
        exact::Rational(7) / exact::Rational(2));
}
```

and, in `invalid_stations_are_rejected`, replace the four `StationEventSource2::build(...)` calls with their exact equivalents:

```cpp
void invalid_stations_are_rejected()
{
    const exact::Rational zero(0);
    const exact::Rational one(1);

    bool raised = false;
    try {
        (void)StationEventSource2::build(zero, zero, zero, one);
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero tool radius was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build(zero, zero, one, zero);
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a zero cap chord ratio was accepted");

    raised = false;
    try {
        (void)StationEventSource2::build(
            zero, zero, one, exact::Rational(17) / exact::Rational(4));
    } catch (const InvalidStationSourceError&) {
        raised = true;
    }
    require(raised, "a cap chord ratio above 4 was accepted");

    // The closed upper end is admitted, and so is a negative coordinate.
    (void)StationEventSource2::build(
        exact::Rational(-1) / exact::Rational(3),
        exact::Rational(-7) / exact::Rational(2),
        one,
        exact::Rational(4));
}
```

Everything else in the gate — every frozen literal, every size assertion, every `require_bytes` — is unchanged. That is the point: the expectations were fixed before the inversion and the inversion is measured against them.

- [ ] **Step 4: Confirm the tree does not yet build, for the expected reason**

Run:
```bash
pixi run -e default exact-gates
```
Expected: **FAIL**, with errors in `station_classifier.cpp`, `segment_strata.cpp`, `circle_strata.cpp` and `segment_oracle.cpp` only — `no matching function for call to 'StationEventSource2::build'` and `no member named 'text' in 'CGAL::Lazy_exact_nt<...>'`. Any error inside `station_source.cpp`, `station_source.h` or the gate is a defect in this task; fix it before moving on.

Do **not** commit. Continue to Task 3.

---

### Task 3: Convert the four consumers

Two producers stop rendering exact values as text for the station to parse back; two consumers stop parsing. This is the step that removes the 13 crossings.

**Files:**
- Modify: `src/continuous_tea_2/station_classifier.cpp`
- Modify: `src/continuous_tea_2/segment_strata.cpp`
- Modify: `src/continuous_tea_2/circle_strata.cpp`
- Modify: `src/continuous_tea_2/segment_oracle.cpp`

**Interfaces:**
- Consumes: the Task 2 accessors, all returning `const compas_cgal::exact::Rational&`.
- Produces: no new symbols. `classify_station_cell`, `construct_station_cell_stratum` and both file-local `station_source` helpers keep their exact signatures, so no header outside `station_source.h` changes.

- [ ] **Step 1: `station_classifier.cpp` — read the carriers instead of decoding them**

Replace lines 85-96 (the three `exact_ft(parse_rational(...))` blocks) with:

```cpp
    // Canonical state, read straight from the source. Since stage 2 the station
    // carries exact::Rational, which IS Epeck::FT, so the classifier decides on
    // the same values the producer computed -- no decode, no re-normalisation.
    const Epeck::FT& center_x = source.center_x();
    const Epeck::FT& center_y = source.center_y();
    const Epeck::FT& radius = source.tool_radius();
```

The `GpsPoint reference(center_x, center_y - radius);` line below is unchanged.

Then mark the two now-dormant file-local helpers, at lines 17 and 45, keeping their bodies untouched:

```cpp
// DORMANT since stage 2: nothing in this translation unit decodes text any
// more. Retained because decoder removal is stage 6, which requires explicit
// user permission. `exact_ft`'s static_assert is the repository's statement that
// CORE::BigRat and CGAL::Epeck_ft are the same type, so it earns its keep here
// regardless.
[[maybe_unused]] Rational parse_rational(
```

```cpp
[[maybe_unused]] Epeck::FT exact_ft(const Rational& value)
```

- [ ] **Step 2: `segment_strata.cpp` — `construct_station_cell_stratum`**

Replace the `return make_cell_stratum(...)` at lines 868-884 with:

```cpp
    return make_cell_stratum(
        branches_at_station(
            records,
            source.center_x(),
            source.center_y(),
            source.tool_radius()),
        "station-rational-v1",
        "0",
        "1",
        // A DECLARED .exact() BOUNDARY (number_types.md R7). make_cell_stratum
        // and branch_pair_dispositions_at still carry CORE::BigRat, so the cap
        // is materialised here, once, rather than decoded from text four lines
        // earlier. Converting those two is stage 3/4, not stage 2.
        CGAL::exact(source.cap_chord_ratio()));
```

Verify the first three arguments type-check without a cast: `branches_at_station` (declared at `segment_strata.cpp:526`) takes `Epeck::FT` centre and radius arguments, which `exact::Rational` is. If it takes them by value, the const references bind and copy; that is correct and needs no change.

Nothing else in this file changes; its `parse_rational` keeps seven other callers.

- [ ] **Step 3: `circle_strata.cpp` — the full-circle producer**

Replace the body of `station_source` (lines 144-178, from the opening brace to the closing brace) with:

```cpp
{
    // The chart direction is built in CORE::BigRat by unit_direction; injecting
    // it into the lazy carrier is an exact rewrap, because Epeck::FT's exact
    // type IS CORE::BigRat. The station coordinates are then computed once, in
    // the filtered carrier, and handed to the source as the numbers they are --
    // where this function used to render four Epeck::FT values as decimal text
    // only to parse them straight back.
    const auto [unit_x, unit_y] =
        unit_direction(
            witness.chart,
            witness.local_parameter);
    const Epeck::FT direction_x(unit_x);
    const Epeck::FT direction_y(unit_y);
    return StationEventSource2::build(
        center_x
            + phase_dx * direction_x
            - phase_dy * direction_y,
        center_y
            + phase_dy * direction_x
            + phase_dx * direction_y,
        tool_radius,
        cap_chord_ratio);
}
```

The signature at lines 136-143 is unchanged. Check the operand order against the original arithmetic: `station_x = center_x + phase_x*unit_x - phase_y*unit_y` and `station_y = center_y + phase_y*unit_x + phase_x*unit_y`. Any transposition here changes geometry, not just bytes.

Nothing else in this file changes: `rational_text` keeps its caller at line 672 (the cell-decision record), `parse_rational` keeps its caller at line 558, and `exact_rational_text` keeps its callers at 559 and 724-729.

- [ ] **Step 4: `segment_oracle.cpp` — the segment-station producer**

Replace the body of `station_source` (lines 72-109) with:

```cpp
{
    const Rational parameter{
        Integer(numerator),
        Integer(denominator)};
    // The segment source still carries text; decoding it is stage 3's to
    // remove. What stage 2 removes is the SECOND crossing: the interpolated
    // station used to be re-rendered as decimal text purely so
    // StationEventSource2 could parse it back.
    const Rational x0 =
        parse_rational(source.x0().text(), "start x");
    const Rational y0 =
        parse_rational(source.y0().text(), "start y");
    const Rational x1 =
        parse_rational(source.x1().text(), "end x");
    const Rational y1 =
        parse_rational(source.y1().text(), "end y");
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    const Rational cap =
        parse_rational(
            source.cap_chord_ratio().text(),
            "cap chord ratio");
    return StationEventSource2::build(
        exact::Rational(x0 + parameter * (x1 - x0)),
        exact::Rational(y0 + parameter * (y1 - y0)),
        exact::Rational(radius),
        exact::Rational(cap));
}
```

Add the namespace alias at the top of the anonymous namespace, immediately after `using Rational = CORE::BigRat;` at line 25:

```cpp
namespace exact = compas_cgal::exact;
```

`station_classifier.h` already pulls in `station_source.h`, which now pulls in `exact/rational.h`, so no new include is needed.

- [ ] **Step 5: Build the gates**

Run:
```bash
pixi run -e default exact-gates
```
Expected: four `OK` lines ending with `exact_station_attestation_gate OK`, exit status 0.

`exact_station_attestation_gate OK` here is the byte-identity result: the same frozen literals that passed against the text carrier in Task 1 now pass against the exact carrier. If `require_bytes` raises `AttestationByteDriftError`, the inversion moved the digest — **stop and report the produced bytes**; do not edit the literal.

- [ ] **Step 6: Rebuild the Python extension**

Run:
```bash
pixi run -e default _editable-rebuild
```
Expected: exit status 0 and no compiler errors. This compiles `station_source.h` into the `_continuous_tea_2` nanobind module and is the check that the header's new `exact/rational.h` include resolves in the module target as well as in the static library.

- [ ] **Step 7: Run the station-lane suites**

Run:
```bash
pixi run -e default pytest -- \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_motion_refutation.py \
  -n auto -q
```
Expected: `tests/adaptive/test_circle_oracle.py::test_legacy_binary64_identity_survives_every_full_circle_outcome` passes. That test asserts three `LEGACY_FULL_TRACE_SHA256` digests over traces whose bytes embed `station-event-source-v1` records, so it is the end-to-end confirmation that the inversion did not move the replay digest. If it fails, **stop and report** — the full comparison is Task 5, but this single test failing means the inversion is not byte-preserving and there is nothing to gain from continuing.

- [ ] **Step 8: Commit Tasks 2 and 3 together**

```bash
git commit -m "refactor: station source carries exact rationals, projects the text" -- \
  src/continuous_tea_2/station_source.h \
  src/continuous_tea_2/station_source.cpp \
  src/continuous_tea_2/station_classifier.cpp \
  src/continuous_tea_2/segment_strata.cpp \
  src/continuous_tea_2/circle_strata.cpp \
  src/continuous_tea_2/segment_oracle.cpp \
  tests/native/exact_station_attestation_gate.cpp
```

---

### Task 4: Generic-double witness with a hard time bound

Every converted lane needs one test on generic double coordinates with a wall-clock ceiling. Integer fixtures make every square root a perfect square, so CORE's approximation error is exactly zero and the refinement path is never entered — `docs/number_types.md` measures that as 0.000 s versus 19.6 s on the same topology. The existing station coverage is exactly that kind of fixture: `test_circle_oracle.py` centres the tool at `(5.0, 5.0)` with radius `0.5` and cap `4.0` on a `10 × 10` integer square.

**Files:**
- Create: `tests/adaptive/test_station_generic_double_witness.py`

**Interfaces:**
- Consumes: `compas_cgal._continuous_tea_2.audit_full_circle_tea_event_exact`, `compas_cgal._continuous_tea_2.segment_station_cap_exceeded_exact`, `compas_cgal._stock_2.Stock2`.
- Produces: nothing importable; a regression ceiling that attributes a stall to this lane.

- [ ] **Step 1: Measure the fixture before writing the ceiling**

The ceiling must be a measured number, not a guess. Write this probe to the scratchpad — **not** to the repository — and run it:

```python
# scratchpad only, never committed
import time

import numpy as np

from compas_cgal import _continuous_tea_2, _stock_2

POCKET = np.array(
    [
        [0.3141, 0.2718, 0.0],
        [9.7183, 0.4142, 0.0],
        [9.4949, 9.6180, 0.0],
        [0.5772, 9.3010, 0.0],
    ],
    dtype=np.float64,
)


def build_stock() -> "_stock_2.Stock2":
    stock = _stock_2.Stock2(POCKET, [])
    stock.subtract_disk(4.7312, 5.2891, 1.3797)
    return stock


for label, call in (
    (
        "full_circle",
        lambda: _continuous_tea_2.audit_full_circle_tea_event_exact(
            build_stock(), 4.7312, 5.2891, 0.9137, 0.4063, False, 0.6180, 3.1416
        ),
    ),
    (
        "segment_station",
        lambda: _continuous_tea_2.segment_station_cap_exceeded_exact(
            build_stock(), 2.7183, 3.1416, 7.3891, 6.2832, 3, 7, 0.9137, 3.1416
        ),
    ),
):
    for attempt in range(3):
        start = time.perf_counter()
        try:
            result = call()
        except Exception as error:  # noqa: BLE001 - probe only
            result = f"{type(error).__name__}: {error}"
        print(label, attempt, f"{time.perf_counter() - start:.3f}s", str(result)[:120])
```

Run: `pixi run -e default python /path/to/scratchpad/probe_station_witness.py`

Record three things:

1. **Does `audit_full_circle_tea_event_exact` reach the station lane?** It does if `b"station-event-source-v1" in trace.canonical_bytes` — that record is produced only by `construct_full_circle_cell_authority`. Reaching it requires a stock with at least one *circle* boundary feature, which `subtract_disk` supplies; a stock whose boundary is all lines takes the `cap_exceeded` early return at `circle_oracle.cpp:782` and never touches the station lane. If the probe does not reach it, adjust the disk radius and centre until it does, and record what worked.
2. **Does `segment_station_cap_exceeded_exact` return a bool?** It raises `IncompleteSegmentOracleError("exact station disposition is unresolved")` when the station is `UNRESOLVED`. If the probe raises, move the motion endpoints until a definite `True` or `False` comes back, and record which.
3. **The slowest of the three attempts for each call**, in seconds.

Set each ceiling to `max(2.0, 4 × slowest)`, rounded up to one decimal, and write the measured value into the test as a comment.

!!! warning "Measured 2026-09-08 — generic coordinates alone are NOT a stressing fixture"

    On this lane, integer 0.149 ms vs generic doubles 0.208 ms is only **1.4×**, not the
    10⁴–10⁷× the row in `docs/number_types.md` describes. That row's driver is a CORE
    **exact-zero identity decision** forcing root-bound refinement, not genericity itself.
    Of six probed configurations only a **near-degenerate sweep** — tool radius comparable to
    half the segment length — moved: **7.35 ms against a 0.33 ms baseline, 22×**. Station at
    either endpoint, a collinear-ish segment, and a near-zero cap ratio were all flat.

    So the fixture must be degenerate-sweep geometry on generic doubles. A merely non-integer
    fixture passes trivially and certifies nothing — the same failure the witness exists to
    prevent, one level up. By the rule above the ceiling is **2.0 s** (floor-dominated).
    Honest limit: this stresses the lane without reproducing the documented worst case. Four times the slowest observed run is wide enough that machine-to-machine variation does not flake, and tight enough that the 10⁴–10⁷× stalls this lane is exposed to fail loudly.

- [ ] **Step 2: Write the witness test**

Create `tests/adaptive/test_station_generic_double_witness.py`, substituting the fixture that Step 1 established and the two measured ceilings:

```python
"""Generic-double witness for the station lane, with a wall-clock ceiling.

The station lane's existing coverage centres the tool at (5.0, 5.0) with radius
0.5 on a 10x10 integer square. On fixtures like that every square root is a
perfect square, CORE's approximation error is exactly zero, and the refinement
path is never entered -- docs/number_types.md measures the same topology at
0.000 s on integer coordinates and 19.6 s on generic doubles. This test uses
coordinates with no structure, so the exact machinery is actually exercised, and
bounds the clock so a regression is attributable to this lane rather than
surfacing as a suite that quietly takes an hour.
"""

from __future__ import annotations

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
GENERIC_TOOL_RADIUS = 0.6180
GENERIC_CAP_CHORD_RATIO = 3.1416

# A ceiling of 0.0 fails by construction. That is deliberate: the value must come
# from Step 1's measurement on this machine, and a test that cannot pass until it
# does is the only version of this constant that cannot be committed unmeasured.
# The rule is max(2.0, 4 x slowest of three runs), rounded up to one decimal --
# wide enough that machine variation does not flake, tight enough that the
# 10^4-10^7x refinement stalls this lane is exposed to fail loudly. Replace both
# zeros and record the measured seconds in the comment beside each.
FULL_CIRCLE_SECONDS_CEILING = 0.0  # measured <seconds> s, ceiling max(2.0, 4x)
SEGMENT_STATION_SECONDS_CEILING = 0.0  # measured <seconds> s, ceiling max(2.0, 4x)


def _generic_stock() -> _stock_2.Stock2:
    stock = _stock_2.Stock2(GENERIC_POCKET, [])
    stock.subtract_disk(
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        GENERIC_DISK_RADIUS,
    )
    return stock


def test_full_circle_station_lane_decides_on_generic_doubles_within_budget() -> None:
    stock = _generic_stock()

    start = time.perf_counter()
    verdict, trace = _continuous_tea_2.audit_full_circle_tea_event_exact(
        stock,
        GENERIC_DISK_CENTER_X,
        GENERIC_DISK_CENTER_Y,
        0.9137,
        0.4063,
        False,
        GENERIC_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    )
    elapsed = time.perf_counter() - start

    # Proof the station lane was actually reached: the record is emitted only by
    # construct_full_circle_cell_authority, and it embeds the station source's
    # own frozen bytes.
    assert b"full-circle-cell-decision-v1" in trace.canonical_bytes
    assert b"station-event-source-v1" in trace.canonical_bytes
    assert b"exact-rational-v1" in trace.canonical_bytes
    assert verdict in {"certified", "cap_exceeded", "unresolved"}
    assert elapsed < FULL_CIRCLE_SECONDS_CEILING, (
        f"full-circle station lane took {elapsed:.3f}s on generic doubles, "
        f"ceiling {FULL_CIRCLE_SECONDS_CEILING:.1f}s"
    )


def test_segment_station_decides_on_generic_doubles_within_budget() -> None:
    stock = _generic_stock()

    start = time.perf_counter()
    exceeded = _continuous_tea_2.segment_station_cap_exceeded_exact(
        stock,
        2.7183,
        3.1416,
        7.3891,
        6.2832,
        3,
        7,
        0.9137,
        GENERIC_CAP_CHORD_RATIO,
    )
    elapsed = time.perf_counter() - start

    assert isinstance(exceeded, bool)
    assert elapsed < SEGMENT_STATION_SECONDS_CEILING, (
        f"segment station decision took {elapsed:.3f}s on generic doubles, "
        f"ceiling {SEGMENT_STATION_SECONDS_CEILING:.1f}s"
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
        GENERIC_TOOL_RADIUS,
        GENERIC_CAP_CHORD_RATIO,
    ]
    for value in coordinates:
        assert value != round(value), f"{value} is an integer fixture"
        assert (value * 16.0) != round(value * 16.0), (
            f"{value} is a small dyadic fixture"
        )
```

- [ ] **Step 3: Run the witness**

Run:
```bash
pixi run -e default pytest -- tests/adaptive/test_station_generic_double_witness.py -v -n auto
```
Expected: 3 passed. If the reach assertions fail, return to Step 1 and adjust the fixture until the lane is reached — do not weaken the assertion, and do not mark anything skip or xfail.

- [ ] **Step 4: Commit**

```bash
git commit -m "test: generic-double station witness with a measured time bound" -- \
  tests/adaptive/test_station_generic_double_witness.py
```

---

### Task 5: Lane equivalence and the station-lane decoder ratchet

**Files:**
- Create: `tests/build/test_station_lane_decoder_ratchet.py`

**Interfaces:**
- Consumes: `build/stage2/baseline.xml` (Task 1 Step 5), the six converted files.
- Produces: `pixi run -e default pytest -- tests/build/test_station_lane_decoder_ratchet.py`, the mechanical guard against a text crossing coming back.

- [ ] **Step 1: Re-run the baseline suites and diff per test**

Run:
```bash
pixi run -e default pytest -- \
  tests/adaptive/test_circle_oracle.py \
  tests/adaptive/test_segment_oracle.py \
  tests/adaptive/test_motion_refutation.py \
  tests/adaptive/test_motion_certificate.py \
  tests/adaptive/test_generator.py \
  tests/adaptive/test_event_substrate.py \
  tests/adaptive/test_segment_event_substrate.py \
  -n auto -q --junitxml=build/stage2/after.xml
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


before = outcomes("build/stage2/baseline.xml")
after = outcomes("build/stage2/after.xml")

changed = {k: (before.get(k), after.get(k)) for k in before | after if before.get(k) != after.get(k)}
print(f"before: {len(before)} tests, after: {len(after)} tests")
for key, (was, now) in sorted(changed.items()):
    print(f"CHANGED {key}: {was} -> {now}")
print("identical" if not changed else f"{len(changed)} tests changed outcome")
PY
```

Expected: `identical`. Anything else is a behavioural change and must be reported with the named tests — stage 2 is a representation change, so the correct outcome set is exactly the one recorded before the inversion, whatever that set was. **Do not adjust a failing test.**

- [ ] **Step 2: Run the full suite**

Run:
```bash
pixi run -e default baseline
```
Expected: completes. Record the summary line verbatim in the Step 4 commit message, so a later stage has a number to compare against — this is the record stage 0 omitted.

- [ ] **Step 3: Write the station-lane decoder ratchet**

Create `tests/build/test_station_lane_decoder_ratchet.py`:

```python
"""Mechanically prevent a text crossing from returning to the station lane.

Stage 2 inverted the station source: it carries `compas_cgal::exact::Rational`
and projects the text, so nothing in the lane decodes a number from a string any
more. The counts below are the design's decoder ratchet, scoped to the files
stage 2 owns. They may only ever DECREASE. An increase means a consumer started
parsing again, which is the exact regression the inversion removed.

The remaining occurrences are deliberate and named:

- station_source.cpp keeps the definition plus the one call inside
  `ExactRational2::build`, the legacy decoder retained so
  `exact_station_attestation_gate` can keep proving the projection agrees with
  it byte for byte.
- station_classifier.cpp keeps its definition only; it is dormant and marked
  `[[maybe_unused]]`.
- segment_strata.cpp, circle_strata.cpp and segment_oracle.cpp keep the callers
  that belong to the segment and full-circle lanes, which stages 3 and 4 own.

Removing any definition is stage 6 and requires explicit user permission, so
this file asserts equality rather than zero.
"""

from __future__ import annotations

from pathlib import Path

REPOSITORY_ROOT = Path(__file__).parents[2]

# Measured 2026-09-08. Before stage 2 these were, in the same order:
# 4, 4, 12, 6, 11 -- a total of 37. Stage 2 removed 13 call sites and no
# definitions.
EXPECTED_DECODER_OCCURRENCES = {
    "src/continuous_tea_2/station_source.cpp": 2,
    "src/continuous_tea_2/station_classifier.cpp": 1,
    "src/continuous_tea_2/segment_strata.cpp": 8,
    "src/continuous_tea_2/circle_strata.cpp": 2,
    "src/continuous_tea_2/segment_oracle.cpp": 11,
}

# The accessors that used to hand a consumer a string to decode. `.text()` on a
# station value is gone from the lane entirely; a reappearance means the
# attestation view is being used as a carrier again.
FORBIDDEN_STATION_TEXT_READS = (
    "source.center_x().text()",
    "source.center_y().text()",
    "source.tool_radius().text()",
    "source.cap_chord_ratio().text()",
)

STATION_CONSUMERS = (
    "src/continuous_tea_2/station_classifier.cpp",
    "src/continuous_tea_2/segment_strata.cpp",
    "src/continuous_tea_2/circle_strata.cpp",
)


def test_station_lane_decoder_count_has_not_increased() -> None:
    observed = {
        path: (REPOSITORY_ROOT / path).read_text(encoding="utf-8").count("parse_rational")
        for path in EXPECTED_DECODER_OCCURRENCES
    }
    assert observed == EXPECTED_DECODER_OCCURRENCES


def test_no_station_consumer_reads_a_number_out_of_text() -> None:
    for path in STATION_CONSUMERS:
        source = (REPOSITORY_ROOT / path).read_text(encoding="utf-8")
        for pattern in FORBIDDEN_STATION_TEXT_READS:
            assert pattern not in source, f"{path} decodes a station value from text"


def test_the_station_source_projects_through_the_single_door() -> None:
    """`to_canonical` is called from exactly one place in the station lane."""
    source = (REPOSITORY_ROOT / "src/continuous_tea_2/station_source.cpp").read_text(
        encoding="utf-8"
    )
    assert source.count("to_canonical(") == 1
    assert "ExactRational2 ExactRational2::project(" in source
```

Note the count for `segment_oracle.cpp` is 11 — unchanged. Stage 2 removes no decoder call there because those calls read `SegmentEventSource2`, which is stage 3's carrier; what stage 2 removed in that file is the four-use `text` lambda on the *output* side.

- [ ] **Step 4: Run the ratchet and commit**

Run:
```bash
pixi run -e default pytest -- tests/build/test_station_lane_decoder_ratchet.py -v -n auto
```
Expected: 3 passed. If `test_station_lane_decoder_count_has_not_increased` fails, print the observed dict — either a conversion was missed or a count in this plan is stale; fix the code, not the expectation, unless the diff proves the expectation was wrong.

```bash
git commit -m "test: station-lane decoder ratchet, 37 to 24 occurrences" -- \
  tests/build/test_station_lane_decoder_ratchet.py
```

---

### Task 6: Close the stage in the documentation

Documentation is a completion artifact for the stage, not a retrospective pass. `docs/number_types.md` currently states that **nothing** uses the exact vocabulary; after stage 2 that is false, and a page that describes the previous state invalidates the completion claim.

**Files:**
- Modify: `docs/number_types.md`

**Interfaces:**
- Consumes: the measured surface table at the top of this plan, the gate output from Task 3 Step 5, the equivalence result from Task 5 Step 1.
- Produces: no code.

- [ ] **Step 1: Replace the stage-0 status admonition**

At `docs/number_types.md:659-671`, replace the `!!! note "Status at stage 0: landed, not adopted"` block with a stage-2 status that states exactly what is adopted and what is not:

```markdown
!!! note "Status at stage 2: the station lane is inverted"

    `StationEventSource2` carries `exact::Rational` and projects its
    `exact-rational-v1` text through `ExactRational2::project`, the lane's single
    caller of `to_canonical`. Its four consumers — `station_classifier.cpp`,
    `segment_strata.cpp`, `circle_strata.cpp`, `segment_oracle.cpp` — read the
    exact values directly; 13 of the lane's 37 `parse_rational` occurrences are
    gone and `tests/build/test_station_lane_decoder_ratchet.py` holds the count.

    **Everything else is unchanged.** `SegmentEventSource2` and
    `FullCircleEventSource2` still carry text (stage 3), the 504 unfiltered
    `BigRat` sites in `segment_site_*` are untouched (stage 4), and no
    `parse_rational` definition has been removed (stage 6). `from_binary64` still
    has no production caller: the station lane takes no doubles, so stage 2
    adopts the exit door only.
```

- [ ] **Step 2: Give `AttestationByteDriftError` its raiser in the error table**

At `docs/number_types.md:649`, replace the `AttestationByteDriftError` row and the paragraph at `:651-657`:

```markdown
| `AttestationByteDriftError` | a projection produces bytes differing from the frozen contract | raised by `exact_station_attestation_gate`, the station lane's contract test |
```

```markdown
`AttestationByteDriftError` acquired its raiser at stage 2, in
`tests/native/exact_station_attestation_gate.cpp` — the contract test, exactly as
designed, and **not** a release-build check on every projection, where the
comparison would cost more than the guarantee is worth. The gate anchors two
framings the repository previously pinned nowhere: the 56-byte
`exact-rational-v1` record for -355/113 as a raw literal, and the 281-byte
`station-event-source-v1` record structurally, through a helper that
re-implements the length framing independently of `encode_string_sequence`. Each
alone is escapable — a raw literal cannot express the nested outer record, and a
helper that shared the encoder would move with it — so both are present and are
tied to each other by an equality check.
```

- [ ] **Step 3: Add the stage-2 subsection**

After `### Measured: the carrier does not move the bytes` (which ends at `:637`), insert:

```markdown
### Stage 2: what the inversion actually bought

The station source is the smallest of the three string carriers, chosen so the
pattern was proven on six files before stage 3 repeats it across eighteen.

| Claim | Evidence |
|---|---|
| The canonical bytes did not move | `exact_station_attestation_gate` — frozen 56-byte and 281-byte anchors, green before and after the inversion |
| The lazy carrier canonicalises like the eager one | the same depth-16 rational chain computed in `Epeck::FT` and in `CORE::BigRat` projects to identical numerator and denominator strings |
| The projection matches the legacy decoder on non-dyadic values | eight generic reduced fractions, including a 31-digit numerator, a migrated sign, an unreduced pair and zero |
| The replay digests are unchanged end to end | `test_legacy_binary64_identity_survives_every_full_circle_outcome` — three `LEGACY_FULL_TRACE_SHA256` anchors over traces whose bytes embed the station record |
| No lane behaviour changed | per-test junit comparison of seven suites, before and after, identical outcome set |
| The lane is exercised on generic doubles, with a bound | `tests/adaptive/test_station_generic_double_witness.py` |

Two decisions are worth keeping, because the code no longer shows them:

**`exact::CanonicalRational::canonical_bytes()` is not the station lane's exit.**
It hard-codes the `exact-binary64-rational-v1` tag; the station lane frames as
`exact-rational-v1`. Adding a tag parameter to `exact/` would have given the
number vocabulary knowledge of a particular event source, which is the one thing
`exact/` must not learn. So the framing stayed with the source, in
`ExactRational2::canonical_bytes()`, fed by `to_canonical`'s numerator and
denominator. The single-door invariant is about the *decomposition*, not about
the framing.

**The tag `exact-rational-v1` is used with two different framings.** The station
lane's is three fields through `encode_string_sequence`;
`boundary_events.cpp:156` emits a one-field `tagged_record` under the same name.
Both are frozen, neither can move, and stage 2 did not unify them because doing
so would change bytes. It is recorded here so a future reader does not assume one
anchor speaks for both.

**Speed was not the goal and did not arrive.** The station lane makes shallow
decisions — a sign, a comparison against 4, an oriented-side query — and at depth
1 the lazy filter is worth 1.39x against its own DAG-allocation overhead. What
stage 2 bought is the removal of a decode-and-renormalise round trip per value
and a lane where a value that reaches a branch has a filter beneath it. The
measured payoff is predicted for stage 4, where the deep chains are.
```

- [ ] **Step 4: Update the carrier table**

At `docs/number_types.md:31` the design's problem table lists the string carriers at 156 parse-back crossings. That table lives in the design document, not here; instead check `docs/number_types.md:22-32` (the six-lane table) and confirm no row claims the station lane is unfiltered. If the L6 row or its surrounding prose names `continuous_tea_2` string carriers, amend it to say the station lane moved to L2 at stage 2 and the segment and full-circle lanes have not.

- [ ] **Step 5: Commit**

```bash
git commit -m "docs: station lane inverted, byte anchors and what stage 2 did not buy" -- \
  docs/number_types.md
```

---

## Stage-2 completion criteria

All of the following, or stage 2 is not done:

1. `pixi run -e default exact-gates` exits 0 with four `OK` lines, the fourth being `exact_station_attestation_gate OK`.
2. The per-test junit comparison in Task 5 Step 1 prints `identical`. The baseline was captured from this tree before any production change, so unlike stage 0 there is something to compare against.
3. `pixi run -e default baseline` completes and its summary line is recorded in the Task 5 commit message.
4. `pixi run -e default pytest -- tests/adaptive/test_station_generic_double_witness.py -n auto` passes, with both ceilings set from measurement and the measured value written in the comment.
5. `pixi run -e default pytest -- tests/build/test_station_lane_decoder_ratchet.py -n auto` passes: 24 `parse_rational` occurrences across the five station-lane translation units, down from 37.
6. `grep -rn "class_<" src/ | grep -i "station\|exactrational"` still returns nothing, and `git diff <the commit stage 2 started from>..HEAD -- src/compas_cgal/_continuous_tea_2.pyi src/continuous_tea_2.cpp` is empty — the Python surface and the binding layer did not move. Capture that starting commit with `git rev-parse HEAD` before Task 1 and record it in the Task 6 commit message.
7. No `parse_rational` definition, and no `ExactRational2::build`, was deleted.
8. `docs/number_types.md` describes stage 2, not stage 0, and its `AttestationByteDriftError` row names a raiser.

## Open questions for the lead

1. **The witness fixture is unmeasured.** Task 4 Step 1 exists because I could not determine, without running the built extension, whether generic doubles drive `audit_full_circle_tea_event_exact` into `construct_full_circle_cell_authority` rather than the `cap_exceeded` early return at `circle_oracle.cpp:782`, nor whether `segment_station_cap_exceeded_exact` returns a bool rather than raising `IncompleteSegmentOracleError` on a generic motion. The plan therefore specifies a measurement and a decision rule rather than a fixture and a number. If the implementing agent cannot find a generic-double fixture that reaches the lane, that is a finding about the lane's reachability and should come back to you rather than be worked around.
2. **Two framings share the tag `exact-rational-v1`** — the station lane's three-field `encode_string_sequence` and `boundary_events.cpp:156`'s one-field `tagged_record`. Both are frozen so neither can be unified without a version bump. Stage 2 documents it; whether it deserves a versioned migration of its own is your call.
3. **`segment_oracle.cpp`'s station producer keeps six `parse_rational` calls** because they decode `SegmentEventSource2`, which stage 3 owns. Stage 2 removes only the second crossing there (rendering the interpolated station back to text). If you would rather that whole function convert at once, it belongs in stage 3's plan, not this one.
4. **Task 1's frozen literals are computed by hand in this plan**, not observed from a run. The byte counts (56 and 281) and the field lengths are derived from the encoder at `event_certificate.cpp:47-53` and the tag strings, and Task 1 Step 4 is written to stop and report if the observed bytes differ. If you would rather they were observed first, that is a five-minute change to Task 1's ordering.
