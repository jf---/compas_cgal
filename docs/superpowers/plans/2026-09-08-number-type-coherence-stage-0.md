# Number-Type Coherence — Stage 0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the `src/exact/` number-type vocabulary and its two doors, with a byte-stability contract proven against the existing implementation, changing no behaviour anywhere.

**Architecture:** A new `compas_cgal::exact` namespace owns the project's exact-number vocabulary (`Rational` = `Epeck::FT`, `OneRoot` = `Sqrt_extension<Rational, Rational>`), one door in (`from_binary64`), one door out (`to_canonical`). Nothing is switched over in this stage — the new module sits beside the existing string carriers and is proven byte-identical to them. Conversions to the new types begin in stage 2.

**Tech Stack:** C++20, CGAL 6.0.1 (vendored `external/cgal`), Boost multiprecision (vendored `external/boost`, GMP disabled), nanobind, CMake + scikit-build-core, pixi, pytest with Hypothesis.

**Spec:** `docs/superpowers/specs/2026-09-08-number-type-coherence-design.md`

## Global Constraints

- Exact arithmetic is settled; do not introduce any epsilon, tolerance, deflation factor or `nextafter` into a decision path.
- One responsibility per file. `exact/` owns the number vocabulary and the two doors only; it must never learn about any particular event source.
- Named exceptions only, convention `<Domain><Condition>Error`, each deriving `std::runtime_error`. Never `throw std::runtime_error("...")` directly.
- No `pytest.mark.skip`, `skipif`, or `xfail` under any circumstance. A failing test fails.
- Never modify an existing reference test to make new code pass.
- `CGAL_DISABLE_GMP` and `CGAL_USE_BOOST_MP` are set repo-wide; `Epeck::FT` is `Lazy_exact_nt<boost::multiprecision::cpp_rational>`. Do not add a GMP dependency.
- Adding a `.cpp` file requires it to be listed in `CMakeLists.txt` and a full rebuild (`pixi run -e default _editable-rebuild`).
- Run pytest with `-n auto`.
- Commit messages: extremely concise, lowercase, no attribution trailers.

---

### Task 1: Pin the binary64 → rational contract on the existing path

Characterization test only. No production code changes. This is the safety net every later task is measured against, so it lands first.

**Files:**
- Test: `tests/adaptive/test_exact_binary64_contract.py` (create)

**Interfaces:**
- Consumes: the existing public API `compas_cgal._continuous_tea_2.SegmentEventSource2.from_binary64`
- Produces: a pinned, oracle-backed statement of what `x0.numerator` / `x0.denominator` must equal for any finite double. Tasks 3 and 5 assert the new path against this same oracle.

- [ ] **Step 1: Write the failing test**

```python
"""Pin the binary64 -> exact-rational contract the attestation view depends on.

`fractions.Fraction(d)` is an exact, independent oracle: a binary64 IS a dyadic
rational, so its exact value is representable with no approximation.
"""

from fractions import Fraction

from hypothesis import given, settings
from hypothesis import strategies as st

from compas_cgal._continuous_tea_2 import SegmentEventSource2

# Values a uniform float strategy will essentially never generate, and where a
# bit-decomposition bug would actually live.
EDGE_DOUBLES = [
    5e-324,            # smallest positive subnormal
    2.2250738585072014e-308,  # smallest positive normal
    1.0,
    0.5,
    2.0**52,
    2.0**53,
    2.0**53 + 2.0,     # first integer gap above 2**53
    2.0**-1074,
    0.1,               # not exactly representable in decimal
    1e308,
]


def _source(value: float) -> SegmentEventSource2:
    """Build a source whose x0 carries `value` and whose other fields are valid."""
    return SegmentEventSource2.from_binary64(
        value, 0.0, value + 1.0, 1.0, 1.0, 1.0
    )


def _assert_matches_oracle(value: float) -> None:
    exact = Fraction(value)
    rational = _source(value).x0
    assert int(rational.numerator) == exact.numerator
    assert int(rational.denominator) == exact.denominator
    assert exact.denominator > 0


@given(
    st.floats(
        allow_nan=False,
        allow_infinity=False,
        min_value=-1e6,
        max_value=1e6,
    )
)
@settings(max_examples=400)
def test_binary64_lift_is_exact(value: float) -> None:
    _assert_matches_oracle(value)


def test_binary64_lift_is_exact_at_edges() -> None:
    for value in EDGE_DOUBLES:
        _assert_matches_oracle(value)
        _assert_matches_oracle(-value)


def test_binary64_lift_denominator_is_a_power_of_two() -> None:
    """Every binary64 is dyadic, so the reduced denominator is a power of two."""
    for value in EDGE_DOUBLES:
        denominator = int(_source(value).x0.denominator)
        assert denominator > 0
        assert denominator & (denominator - 1) == 0
```

- [ ] **Step 2: Run the test**

Run: `pixi run -e default pytest tests/adaptive/test_exact_binary64_contract.py -v -n auto`

Expected: PASS. This characterizes existing behaviour. If any assertion FAILS, **stop and report** — a failure here means the current attestation path is already wrong, which changes the whole plan and must be raised before any refactor proceeds.

- [ ] **Step 3: Commit**

```bash
git add tests/adaptive/test_exact_binary64_contract.py
git commit -- tests/adaptive/test_exact_binary64_contract.py -m "test: pin binary64 to exact-rational contract"
```

---

### Task 2: Error model and the exact-number vocabulary

**Files:**
- Create: `src/exact/errors.h`
- Create: `src/exact/rational.h`
- Create: `src/exact/rational.cpp`
- Modify: `CMakeLists.txt` (add `src/exact/rational.cpp` to `continuous_tea_exact_core`)
- Test: `tests/native/exact_vocabulary_gate.cpp` (create)
- Modify: `CMakeLists.txt` (add the `exact_vocabulary_gate` executable)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `compas_cgal::exact::Rational` (alias of `CGAL::Exact_predicates_exact_constructions_kernel::FT`)
  - `compas_cgal::exact::Rational compas_cgal::exact::from_binary64(double)`
  - `compas_cgal::exact::NonFiniteBinary64Error`
  - `compas_cgal::exact::UnreducedCanonicalRationalError`
  - `compas_cgal::exact::CrossRootExtensionError`

- [ ] **Step 1: Write the failing native gate**

```cpp
// tests/native/exact_vocabulary_gate.cpp
#include "exact/rational.h"
#include "exact/errors.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

namespace exact = compas_cgal::exact;

namespace {

void from_binary64_is_exact()
{
    // 0.1 is not a decimal fraction; its exact value is 3602879701896397 / 2^55.
    const exact::Rational tenth = exact::from_binary64(0.1);
    const exact::Rational expected =
        exact::Rational(3602879701896397LL) / exact::Rational(1LL << 55);
    assert(tenth == expected);

    // Smallest positive subnormal must survive exactly.
    const exact::Rational tiny = exact::from_binary64(5e-324);
    assert(tiny > exact::Rational(0));

    // Negative zero and positive zero denote the same rational.
    assert(exact::from_binary64(-0.0) == exact::from_binary64(0.0));
    assert(exact::from_binary64(0.0) == exact::Rational(0));
}

void non_finite_input_is_rejected()
{
    bool raised = false;
    try {
        exact::from_binary64(std::numeric_limits<double>::quiet_NaN());
    } catch (const exact::NonFiniteBinary64Error&) {
        raised = true;
    }
    assert(raised);

    raised = false;
    try {
        exact::from_binary64(std::numeric_limits<double>::infinity());
    } catch (const exact::NonFiniteBinary64Error&) {
        raised = true;
    }
    assert(raised);
}

}  // namespace

int main()
{
    from_binary64_is_exact();
    non_finite_input_is_rejected();
    std::printf("exact_vocabulary_gate OK\n");
    return 0;
}
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
pixi run -e default cmake -S . -B build/exact-gate -G Ninja -Dnanobind_DIR=$(pixi run -e default python -m nanobind --cmake_dir)
pixi run -e default cmake --build build/exact-gate --target exact_vocabulary_gate
```
Expected: FAIL — `exact/rational.h` does not exist.

- [ ] **Step 3: Write `src/exact/errors.h`**

```cpp
#pragma once

#include <stdexcept>

namespace compas_cgal::exact {

/// Raised when a double that must denote an exact rational is NaN or infinite.
class NonFiniteBinary64Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when a canonical rational violates the reduced / positive-denominator
/// / gcd == 1 invariant that makes attestation bytes value-determined.
class UnreducedCanonicalRationalError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Raised when one-root arithmetic is attempted across distinct roots. CGAL
/// documents this as a precondition; violating it is undefined behaviour.
class CrossRootExtensionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

}  // namespace compas_cgal::exact
```

- [ ] **Step 4: Write `src/exact/rational.h`**

```cpp
#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace compas_cgal::exact {

/// The project's exact rational carrier.
///
/// This is deliberately the field type of Epeck rather than a fresh type: it is
/// already the coefficient type of Gps_circle_segment_traits_2<Epeck>'s CoordNT,
/// so one-root numbers built on it interoperate with the arrangement package
/// with no conversion, and every value carries Lazy_exact_nt's interval filter.
using Rational = CGAL::Exact_predicates_exact_constructions_kernel::FT;

/// Convert a binary64 to its exact rational value.
///
/// This is the ONLY double-to-exact entry point in the codebase. A binary64 is a
/// dyadic rational, so the conversion is exact and total on finite input: there
/// is no parsing, no tolerance and no snapping.
///
/// Args:
///     value: a finite binary64.
///
/// Returns:
///     The exact rational denoted by `value`. Negative zero maps to zero.
///
/// Raises:
///     NonFiniteBinary64Error: if `value` is NaN or infinite.
[[nodiscard]] Rational from_binary64(double value);

}  // namespace compas_cgal::exact
```

- [ ] **Step 5: Write `src/exact/rational.cpp`**

```cpp
#include "exact/rational.h"

#include "exact/errors.h"

#include <cmath>

namespace compas_cgal::exact {

Rational from_binary64(const double value)
{
    if (!std::isfinite(value)) {
        throw NonFiniteBinary64Error(
            "exact::from_binary64 requires a finite binary64");
    }
    // Epeck::FT's double constructor is exact: it stores the dyadic rational the
    // binary64 denotes. Constructing through it keeps the lazy interval intact,
    // where a bit-decomposition would build an eager cpp_rational instead.
    return Rational(value);
}

}  // namespace compas_cgal::exact
```

- [ ] **Step 6: Wire both into `CMakeLists.txt`**

Add `src/exact/rational.cpp` to the `continuous_tea_exact_core` source list (the `add_library(continuous_tea_exact_core STATIC ...)` block beginning at line 238), immediately after `src/canonical_encoding.cpp`.

Then add the gate executable, directly after the existing `task3_algorithm_gate` block:

```cmake
add_executable(exact_vocabulary_gate EXCLUDE_FROM_ALL
    tests/native/exact_vocabulary_gate.cpp
)
target_link_libraries(exact_vocabulary_gate PRIVATE continuous_tea_exact_core)
target_compile_definitions(exact_vocabulary_gate PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 7: Run the gate to verify it passes**

Run:
```bash
pixi run -e default cmake --build build/exact-gate --target exact_vocabulary_gate
./build/exact-gate/exact_vocabulary_gate
```
Expected: `exact_vocabulary_gate OK`, exit status 0.

- [ ] **Step 8: Commit**

```bash
git add src/exact/errors.h src/exact/rational.h src/exact/rational.cpp tests/native/exact_vocabulary_gate.cpp
git commit -- src/exact/errors.h src/exact/rational.h src/exact/rational.cpp tests/native/exact_vocabulary_gate.cpp CMakeLists.txt -m "feat: exact number vocabulary and binary64 door"
```

---

### Task 3: The attestation door, proven byte-identical

`SegmentEventSource2::lift_exact` already performs exactly this projection. This task promotes that logic into `exact::to_canonical` and proves the promoted copy is byte-identical to the original before anything depends on it.

**Files:**
- Create: `src/exact/canonical.h`
- Create: `src/exact/canonical.cpp`
- Modify: `CMakeLists.txt` (add `src/exact/canonical.cpp` to `continuous_tea_exact_core`)
- Test: `tests/native/exact_canonical_gate.cpp` (create)
- Modify: `CMakeLists.txt` (add the `exact_canonical_gate` executable)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational`, `compas_cgal::exact::from_binary64` (Task 2).
- Produces:
  - `compas_cgal::exact::CanonicalRational` with `numerator()`, `denominator()`, `text()`, `canonical_bytes()`, all `std::string`
  - `compas_cgal::exact::CanonicalRational compas_cgal::exact::to_canonical(const Rational&)`

- [ ] **Step 1: Write the failing native gate**

```cpp
// tests/native/exact_canonical_gate.cpp
#include "exact/canonical.h"
#include "exact/rational.h"

#include "continuous_tea_2/segment_source.h"

#include <cassert>
#include <cstdio>
#include <vector>

namespace exact = compas_cgal::exact;

namespace {

const std::vector<double>& edge_doubles()
{
    static const std::vector<double> values = {
        5e-324, 2.2250738585072014e-308, 1.0, 0.5,
        4503599627370496.0,      // 2^52
        9007199254740992.0,      // 2^53
        9007199254740994.0,      // 2^53 + 2
        0.1, 3.0, 1e308,
    };
    return values;
}

/// The promoted projection must agree with the existing one, byte for byte.
void to_canonical_matches_existing_projection()
{
    for (const double value : edge_doubles()) {
        for (const double signed_value : {value, -value}) {
            const exact::CanonicalRational promoted =
                exact::to_canonical(exact::from_binary64(signed_value));

            // Existing path, reached through the public source API.
            const SegmentEventSource2 source =
                SegmentEventSource2::from_binary64(
                    signed_value, 0.0, signed_value + 1.0, 1.0, 1.0, 1.0);

            assert(promoted.numerator() == source.x0().numerator());
            assert(promoted.denominator() == source.x0().denominator());
            assert(promoted.text() == source.x0().text());
            assert(promoted.canonical_bytes() == source.x0().canonical_bytes());
        }
    }
}

/// Denominators are always positive and reduced, which is what makes the bytes
/// a function of the value rather than of the carrier.
void canonical_form_is_reduced_and_positive()
{
    for (const double value : edge_doubles()) {
        const exact::CanonicalRational canonical =
            exact::to_canonical(exact::from_binary64(-value));
        assert(!canonical.denominator().empty());
        assert(canonical.denominator().front() != '-');
    }
}

}  // namespace

int main()
{
    to_canonical_matches_existing_projection();
    canonical_form_is_reduced_and_positive();
    std::printf("exact_canonical_gate OK\n");
    return 0;
}
```

- [ ] **Step 2: Run it to verify it fails**

Run:
```bash
pixi run -e default cmake --build build/exact-gate --target exact_canonical_gate
```
Expected: FAIL — `exact/canonical.h` does not exist.

- [ ] **Step 3: Write `src/exact/canonical.h`**

```cpp
#pragma once

#include "exact/rational.h"

#include <string>

namespace compas_cgal::exact {

/// A rational in canonical attestation form: reduced, positive denominator.
///
/// This is a derived VIEW of an exact value, never a carrier. It exists so that
/// attestation bytes are a function of the mathematical value rather than of
/// whichever number type produced it. Construct it only via `to_canonical`.
class CanonicalRational {
public:
    [[nodiscard]] const std::string& numerator() const noexcept;
    [[nodiscard]] const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    [[nodiscard]] std::string text() const;

    /// Frozen attestation bytes. This encoding is a compatibility contract:
    /// changing it invalidates every stored replay digest.
    [[nodiscard]] std::string canonical_bytes() const;

private:
    friend CanonicalRational to_canonical(const Rational& value);

    CanonicalRational(std::string numerator, std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

/// Project an exact rational into canonical attestation form.
///
/// This is the ONLY exact-to-attestation-bytes exit point in the codebase.
///
/// Raises:
///     UnreducedCanonicalRationalError: if the decomposed denominator is not
///         positive, which would make the encoding ambiguous.
[[nodiscard]] CanonicalRational to_canonical(const Rational& value);

}  // namespace compas_cgal::exact
```

- [ ] **Step 4: Write `src/exact/canonical.cpp`**

```cpp
#include "exact/canonical.h"

#include "exact/errors.h"
#include "continuous_tea_2/event_certificate.h"

#include <CGAL/CORE/BigInt.h>
#include <CGAL/Fraction_traits.h>

#include <utility>
#include <vector>

namespace compas_cgal::exact {

CanonicalRational::CanonicalRational(
    std::string numerator,
    std::string denominator)
    : numerator_(std::move(numerator)),
      denominator_(std::move(denominator))
{
}

const std::string& CanonicalRational::numerator() const noexcept
{
    return numerator_;
}

const std::string& CanonicalRational::denominator() const noexcept
{
    return denominator_;
}

std::string CanonicalRational::text() const
{
    return denominator_ == "1" ? numerator_ : numerator_ + "/" + denominator_;
}

std::string CanonicalRational::canonical_bytes() const
{
    return encode_string_sequence(
        {
            "exact-binary64-rational-v1",
            numerator_,
            denominator_,
        });
}

CanonicalRational to_canonical(const Rational& value)
{
    using Traits = CGAL::Fraction_traits<Rational>;
    typename Traits::Numerator_type numerator;
    typename Traits::Denominator_type denominator;
    typename Traits::Decompose()(value, numerator, denominator);

    const CORE::BigInt exact_denominator(denominator.exact());
    if (exact_denominator <= 0) {
        throw UnreducedCanonicalRationalError(
            "canonical rational requires a positive denominator");
    }
    return CanonicalRational(
        CORE::BigInt(numerator.exact()).convert_to<std::string>(),
        exact_denominator.convert_to<std::string>());
}

}  // namespace compas_cgal::exact
```

- [ ] **Step 5: Wire into `CMakeLists.txt`**

Add `src/exact/canonical.cpp` to `continuous_tea_exact_core`, immediately after `src/exact/rational.cpp`. Add the gate executable after `exact_vocabulary_gate`:

```cmake
add_executable(exact_canonical_gate EXCLUDE_FROM_ALL
    tests/native/exact_canonical_gate.cpp
)
target_link_libraries(exact_canonical_gate PRIVATE continuous_tea_exact_core)
target_compile_definitions(exact_canonical_gate PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 6: Run the gate to verify it passes**

Run:
```bash
pixi run -e default cmake --build build/exact-gate --target exact_canonical_gate
./build/exact-gate/exact_canonical_gate
```
Expected: `exact_canonical_gate OK`, exit status 0.

If the byte comparison fails, **stop and report the differing pair**. A mismatch means the two existing lift paths (`lift_binary64` via `CORE::BigRat`, `lift_exact` via `Fraction_traits`) already disagree, which is a live defect and must be raised before proceeding.

- [ ] **Step 7: Commit**

```bash
git add src/exact/canonical.h src/exact/canonical.cpp tests/native/exact_canonical_gate.cpp
git commit -- src/exact/canonical.h src/exact/canonical.cpp tests/native/exact_canonical_gate.cpp CMakeLists.txt -m "feat: exact attestation door, byte-identical to existing projection"
```

---

### Task 4: One-root vocabulary with the cross-root precondition enforced

The cross-root precondition is currently documented only in a comment in `src/engagement_2.cpp`. This task gives it a checked call. `sign_mixed_radical` itself is NOT moved here — that is stage 1.

**Files:**
- Create: `src/exact/one_root.h`
- Create: `src/exact/one_root.cpp`
- Modify: `CMakeLists.txt` (add `src/exact/one_root.cpp` to `continuous_tea_exact_core`)
- Test: `tests/native/exact_one_root_gate.cpp` (create)
- Modify: `CMakeLists.txt` (add the `exact_one_root_gate` executable)

**Interfaces:**
- Consumes: `compas_cgal::exact::Rational` (Task 2), `compas_cgal::exact::CrossRootExtensionError` (Task 2).
- Produces:
  - `compas_cgal::exact::OneRoot` (alias of `CGAL::Sqrt_extension<Rational, Rational>`)
  - `compas_cgal::exact::OneRoot compas_cgal::exact::same_root_add(const OneRoot&, const OneRoot&)`
  - `compas_cgal::exact::OneRoot compas_cgal::exact::same_root_multiply(const OneRoot&, const OneRoot&)`

- [ ] **Step 1: Write the failing native gate**

```cpp
// tests/native/exact_one_root_gate.cpp
#include "exact/one_root.h"
#include "exact/errors.h"

#include <CGAL/number_utils.h>

#include <cassert>
#include <cstdio>

namespace exact = compas_cgal::exact;

namespace {

void same_root_operations_are_exact()
{
    // (1 + sqrt(2)) + (3 + 2*sqrt(2)) == 4 + 3*sqrt(2)
    const exact::OneRoot a(exact::Rational(1), exact::Rational(1), exact::Rational(2));
    const exact::OneRoot b(exact::Rational(3), exact::Rational(2), exact::Rational(2));
    const exact::OneRoot sum = exact::same_root_add(a, b);
    assert(sum == exact::OneRoot(exact::Rational(4), exact::Rational(3), exact::Rational(2)));
    assert(CGAL::sign(sum) == CGAL::POSITIVE);
}

void a_rational_operand_is_always_permitted()
{
    // A non-extended value carries a0() alone and is compatible with any root.
    const exact::OneRoot rational(exact::Rational(5));
    const exact::OneRoot extended(
        exact::Rational(0), exact::Rational(1), exact::Rational(3));
    const exact::OneRoot sum = exact::same_root_add(rational, extended);
    assert(CGAL::sign(sum) == CGAL::POSITIVE);
}

void cross_root_arithmetic_is_rejected()
{
    const exact::OneRoot two(
        exact::Rational(0), exact::Rational(1), exact::Rational(2));
    const exact::OneRoot three(
        exact::Rational(0), exact::Rational(1), exact::Rational(3));
    bool raised = false;
    try {
        (void)exact::same_root_add(two, three);
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    assert(raised);

    raised = false;
    try {
        (void)exact::same_root_multiply(two, three);
    } catch (const exact::CrossRootExtensionError&) {
        raised = true;
    }
    assert(raised);
}

}  // namespace

int main()
{
    same_root_operations_are_exact();
    a_rational_operand_is_always_permitted();
    cross_root_arithmetic_is_rejected();
    std::printf("exact_one_root_gate OK\n");
    return 0;
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `pixi run -e default cmake --build build/exact-gate --target exact_one_root_gate`
Expected: FAIL — `exact/one_root.h` does not exist.

- [ ] **Step 3: Write `src/exact/one_root.h`**

```cpp
#pragma once

#include "exact/rational.h"

#include <CGAL/Sqrt_extension.h>

namespace compas_cgal::exact {

/// One-root algebraic numbers, a0 + a1 * sqrt(root).
///
/// Deliberately a bare alias rather than a checking wrapper: this type IS
/// Gps_circle_segment_traits_2<Epeck>::Point_2::CoordNT, and wrapping it would
/// reintroduce a conversion boundary with the arrangement package. The
/// precondition lives in the free functions below instead.
using OneRoot = CGAL::Sqrt_extension<Rational, Rational>;

/// Add two one-root numbers that share a root.
///
/// CGAL's precondition is that the operands share an extension. It is checked
/// via `is_extended()` first, because `root()` itself is only defined on an
/// extended value; violating either is undefined behaviour rather than a wrong
/// answer, so neither is assumed.
///
/// Raises:
///     CrossRootExtensionError: if the operands carry distinct non-zero roots.
[[nodiscard]] OneRoot same_root_add(const OneRoot& a, const OneRoot& b);

/// Multiply two one-root numbers that share a root.
///
/// Raises:
///     CrossRootExtensionError: if the operands carry distinct non-zero roots.
[[nodiscard]] OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b);

}  // namespace compas_cgal::exact
```

- [ ] **Step 4: Write `src/exact/one_root.cpp`**

```cpp
#include "exact/one_root.h"

#include "exact/errors.h"

namespace compas_cgal::exact {

namespace {

void require_same_root(const OneRoot& a, const OneRoot& b)
{
    // a1() and root() are ONLY defined on an EXTENDED Sqrt_extension. A rational
    // coordinate carries a0() alone, and reading the extension parts
    // unconditionally is undefined behaviour -- the same guard as
    // stock_2.cpp::add_coordinate and engagement_2.cpp::as_radpoint. A
    // non-extended operand is compatible with any root, so it short-circuits.
    if (!a.is_extended() || !b.is_extended()) {
        return;
    }
    if (a.root() != b.root()) {
        throw CrossRootExtensionError(
            "one-root arithmetic requires operands in the same extension");
    }
}

}  // namespace

OneRoot same_root_add(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a + b;
}

OneRoot same_root_multiply(const OneRoot& a, const OneRoot& b)
{
    require_same_root(a, b);
    return a * b;
}

}  // namespace compas_cgal::exact
```

- [ ] **Step 5: Wire into `CMakeLists.txt`**

Add `src/exact/one_root.cpp` to `continuous_tea_exact_core`, after `src/exact/canonical.cpp`. Add the gate after `exact_canonical_gate`:

```cmake
add_executable(exact_one_root_gate EXCLUDE_FROM_ALL
    tests/native/exact_one_root_gate.cpp
)
target_link_libraries(exact_one_root_gate PRIVATE continuous_tea_exact_core)
target_compile_definitions(exact_one_root_gate PRIVATE CGAL_USE_CORE=1)
```

- [ ] **Step 6: Run the gate to verify it passes**

Run:
```bash
pixi run -e default cmake --build build/exact-gate --target exact_one_root_gate
./build/exact-gate/exact_one_root_gate
```
Expected: `exact_one_root_gate OK`, exit status 0.

- [ ] **Step 7: Commit**

```bash
git add src/exact/one_root.h src/exact/one_root.cpp tests/native/exact_one_root_gate.cpp
git commit -- src/exact/one_root.h src/exact/one_root.cpp tests/native/exact_one_root_gate.cpp CMakeLists.txt -m "feat: one-root vocabulary with checked cross-root precondition"
```

---

### Task 5: Stage-0 exit gate — no regression anywhere

Stage 0 changed no behaviour, so the entire existing suite must be untouched and every stored replay digest must be unchanged. This task makes that claim measured rather than assumed, and leaves behind one command that later stages re-run.

**Files:**
- Modify: `pyproject.toml` (add an `exact-gates` pixi task)

**Interfaces:**
- Consumes: the three gate executables from Tasks 2-4.
- Produces: `pixi run -e default exact-gates`, the single command later stages use as their regression gate.

- [ ] **Step 1: Add the pixi task**

In `[tool.pixi.tasks]` in `pyproject.toml`, add:

```toml
_exact-gates-configure = "cmake -S . -B build/exact-gate -G Ninja -Dnanobind_DIR=$(python -m nanobind --cmake_dir)"
_exact-gates-build = "cmake --build build/exact-gate --target exact_vocabulary_gate exact_canonical_gate exact_one_root_gate"
_exact-gates-run = "build/exact-gate/exact_vocabulary_gate && build/exact-gate/exact_canonical_gate && build/exact-gate/exact_one_root_gate"
exact-gates = { depends-on = ["_exact-gates-configure", "_exact-gates-build", "_exact-gates-run"], description = "Native gates for the exact number vocabulary" }
```

- [ ] **Step 2: Run the native gates**

Run: `pixi run -e default exact-gates`
Expected: three `OK` lines, exit status 0.

- [ ] **Step 3: Run the replay and certification suites**

Run: `pixi run -e default pytest tests/test_replay_classification.py tests/test_false_certificate.py tests/adaptive -v -n auto`
Expected: PASS, with the same pass/fail set as before the branch. Any newly failing test means stage 0 changed behaviour and must be investigated before proceeding — **do not adjust the failing test**.

- [ ] **Step 4: Run the full suite**

Run: `pixi run -e default pytest tests -n auto -q`
Expected: identical results to `main` at the branch point. Record the counts in the commit message.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml
git commit -- pyproject.toml -m "chore: exact-gates task, stage 0 exit gate"
```

---

## Stage-0 completion criteria

All of the following, or stage 0 is not done:

1. `pixi run -e default exact-gates` exits 0 with three `OK` lines.
2. `pixi run -e default pytest tests -n auto -q` matches the pre-branch result set exactly.
3. `src/exact/` contains exactly four headers and three translation units; nothing outside `src/exact/` includes them yet.
4. `grep -rn 'parse_rational' src/ | wc -l` is unchanged — stage 0 removes nothing.
5. `docs/number_types.md` gains a short section documenting `exact::Rational`, `exact::OneRoot` and the two doors, per the repository's development-stage documentation rule.
