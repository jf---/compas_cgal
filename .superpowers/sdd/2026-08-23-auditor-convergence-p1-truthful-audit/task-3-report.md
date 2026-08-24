# Task 3 report: exact partial-arc surrogate

## Status

Implementation ready for independent specification and code-quality review.
No commit created.

## Implemented contract

- One CGAL-free frozen four-quarter Pythagorean atlas supplies identifiers,
  polynomial coefficients, domain, orientation, and parameter-zero ownership
  to Epeck depletion and continuous-TEA polynomial construction.
- `AuditArcMotion2` is a private-constructor domain value. The v2 ingress
  observes `end - start` once in binary64 and exact-injects it; Epeck owns sign,
  extent, orientation, normalization, quadrant, and exact seam decisions. The
  sole tan seam uses the original angle plus exact-selected, representable
  integer indices, never an exact-to-double local-angle roundtrip.
- `construct_exact_arc_depletion(...)` emits the one canonical owned parameter
  sequence, refines exact chords dyadically, checks allocation before reserve,
  and creates a content-addressed non-cyclic trace. Full turns retain the
  explicit duplicate terminal anchor as the final consecutive closure segment.
- `exact_arc_structural_density_holds(...)` regenerates and exactly compares the
  declared traversal. Endpoint incidence and chord density alone are not
  accepted as a certificate.
- `Stock2::subtract_exact_arc(...)` constructs and validates against a clone,
  then swaps atomically. The authoritative method accepts only the opaque native
  motion, validates the trace motion digest, and never calls
  `subtract_arc_sweep`.
- Exact rational identity delegates to the repository CCAN encoder after
  `CGAL::Fraction_traits` decomposition; binary64 identity delegates to the
  same repository canonical encoder.

## TDD evidence

Initial RED:

```bash
pixi run audit-native
```

Failed during CMake generation because `src/audit_arc_motion_2.cpp` and the
exact arc API did not exist.

Adversarial RED:

```bash
pixi run cmake --build build/audit-native --target audit_native_gate
```

Failed on absent `ExactArcDepletionTrace2::strategy_version`,
`ExactArcDepletionTrace2::digest`, and the shared frozen-atlas API.

Review-repair RED:

```bash
pixi run cmake --build build/audit-native --target audit_native_gate
```

Failed on absent `ExactArcDepletionTrace2::matches_motion` and named
`AuditDigestSizeError`. The added cases also expose ordinary non-seam full-turn
asymmetry, implicit cyclic closure, foreign-motion trace acceptance, public
nonfinite injection, and trace-factory policy bypass.

Focused native GREEN after the proof corrections:

```bash
pixi run cmake --build build/audit-native --target audit_native_gate
pixi run build/audit-native/audit_native_gate
```

Result: six modified native objects compiled, executable linked, exit 0. The
gate covers exact seams in both directions, awkward/negative/multi-turn angles,
CW/CCW minor/major/full traversals, nextafter seams, rational transforms,
canonical seam ownership, exact incidence/density, foreign trace permutations,
policy/limit rejection, CCAN golden and round-trip values, and motion/tool/
chord/count witness mutations.

The repair gate additionally covers v2 seam identity, ordinary starts `0.37`,
`1.2`, and `-0.37` in both full-turn directions, explicit duplicate terminal
closure, complete opposite-direction complement substitution, geometrically
equal non-owner seam substitution, named errors, and byte-for-byte CCAN
binary64 reuse. An authoritative-source scan contains no `to_double` call.

Focused public boundary GREEN:

```bash
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- \
  tests/engagement_audit/test_native_classification.py \
  tests/engagement_audit/test_arc_authority.py \
  tests/adaptive/test_exact_depletion.py \
  -n auto --testmon -q
```

Result: `_stock_2` and `_continuous_tea_2` rebuilt and linked; 43 passed in
1.78 s. This includes exact stock equality after every public arc-depletion
refusal.

The review-repair exact-arc subset rebuilt `_stock_2` and passed 16 tests in
1.05 s under two xdist workers. The public subtraction and read-only matcher
both reject nonfinite values and signed nonpositive limits before exact or
`size_t` conversion.

Final root-owned verification after the review repairs and SDD alignment:

```bash
pixi run audit-native
PYTEST_XDIST_AUTO_NUM_WORKERS=2 pixi run pytest -- \
  tests/engagement_audit/test_native_classification.py \
  tests/engagement_audit/test_classification.py \
  tests/engagement_audit/test_input.py \
  tests/engagement_audit/test_arc_authority.py \
  tests/adaptive/test_exact_depletion.py \
  -n auto --testmon -q
pixi run lint
pixi run types-audit
pixi run -e docs docs
git diff --check
```

Result: native gate exit 0; 16 exact-arc Python tests passed; five
continuous-TEA atlas regressions passed; 32 affected tests passed under
`--testmon`; Ruff passed; strict mypy passed for nine files; strict MkDocs
passed; diff check passed. Two independent read-only rereviews returned PASS
with no Critical or Important findings.

## Proof boundary

Every removed disk is centered exactly on the declared rational-chart guide,
so the removed union is a subset of that surrogate sweep. Exact
`0 < chord_bound < tool_radius` proves adjacent disks overlap. No quantitative
retained-sliver bound is claimed, and no equality with the ideal transcendental
COMPAS arc is claimed.

Task 4 remains responsible for input schema v2, depletion policy in request
identity, authenticated-operation digest integration, native certification,
and the all-motion atomic replay transaction.
