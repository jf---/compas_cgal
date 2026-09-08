# Number-Type Coherence — Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Converge the two divergent copies of `sign_mixed_radical` into one definition in `src/exact/one_root.*`, proving equivalence before the move rather than assuming it.

**Architecture:** The exact mixed-radical sign predicate exists twice. `src/engagement_2.cpp:43` serves the engagement geometry directly (`:123`, `:141`). `src/audit_exact_station_2.cpp:22` (`sign_mixed_radical_impl`) is wrapped as `audit_sign_mixed_radical_exact` and is what the Python binding `_sign_mixed_radical` actually reaches, via `engagement_2.cpp:694`. Stage 1 makes both call one definition.

**Tech Stack:** C++20, CGAL 6.0.1 (vendored), nanobind, pixi, pytest + Hypothesis.

**Spec:** `docs/superpowers/specs/2026-09-08-number-type-coherence-design.md` (stage 1 row)

## Global Constraints

- This stage DOES touch shipped behaviour — unlike stage 0. `engagement_2.cpp` is linked into the extensions. Every change must be justified by an equivalence proof, not by inspection.
- Exact arithmetic: no epsilon, tolerance, or `nextafter` in any decision path.
- Named exceptions only, `<Domain><Condition>Error`, deriving `std::runtime_error`.
- Public signatures `audit_sign_mixed_radical_exact` (`src/audit_certification_2.h:129`) and `sign_mixed_radical_for_binding` (`src/engagement_2.h:193`) MUST NOT change — they are consumed elsewhere and by the Python binding.
- `CGAL::Sqrt_extension::a1()`/`root()` are defined ONLY on an extended value; guard with `is_extended()` first (`src/stock_2.cpp:912-915`).
- Run pytest with `-n auto`. No `skip`, `skipif`, or `xfail`. Never modify a reference test to make code pass.
- Commit by pathspec, message BEFORE the `--`. Worktree is SHARED with live sessions holding uncommitted work and a staged rename: never `git add -A`, `git commit -a`, `git stash`, revert or checkout.
- Author AND committer `Jelle Feringa <jelleferinga@gmail.com>`. No attribution trailers.
- Native gates compile with `-DNDEBUG`, so `assert` is inert. Use the established throwing `require` + `GateCheckFailedError` pattern (see `tests/native/exact_one_root_gate.cpp`).
- Defect injection must NEVER edit a tracked file: copy the TU to scratch, define the broken function locally so the archive member is not pulled, compile with the flags `ninja -t commands` reports.
- Use your own build directory; delete it when done. Disk is tight.

---

### Task 1: Prove the two copies are equivalent, before touching either

The copies differ textually. Copy B merges two branches copy A keeps separate (`u_sign == ZERO || u_sign == w_sign` → `return w_sign`). That is *believed* equivalent — in the merged branch `u_sign == w_sign`, so either return value is the same — but an exact predicate feeding certificates does not get to rely on a belief.

**Files:**
- Create: `tests/native/sign_mixed_radical_equivalence_gate.cpp`
- Modify: `CMakeLists.txt` (add the gate executable, after `exact_one_root_gate`)

**Interfaces:**
- Consumes copy B via `audit_sign_mixed_radical_exact` (`src/audit_certification_2.h:129`), its public wrapper — `audit_exact_station_2.cpp` is in `continuous_tea_exact_core`, so this links normally.
- Consumes copy A by **`#include "engagement_2.cpp"` directly in the gate TU.** This is required, not stylistic: copy A sits in an anonymous namespace (`src/engagement_2.cpp:22`), so it has internal linkage and cannot be reached by declaration. There is no ODR conflict, because `engagement_2.cpp` belongs to the nanobind module target (`CMakeLists.txt:441`) and is NOT in `continuous_tea_exact_core`, so the gate never links a second copy of that TU. Wrap the include in a comment explaining both facts.
- **Transcribing copy A's algorithm into the gate is FORBIDDEN.** A transcription proves the transcription equivalent to copy B, which is not the question.
- Produces: evidence that both copies agree on every probed input.
- Note: the gate links `continuous_tea_exact_core` for copy B and compiles copy A into itself. If the include pulls unresolved symbols from `engagement_2.cpp`'s other contents, link `continuous_tea_exact_core` first and add whatever else the linker names; if it cannot be made to link at all, STOP and report rather than falling back to transcription.

- [ ] **Step 1: Write the failing gate**

Enumerate sign combinations exhaustively rather than sampling: for `A, B, C, D ∈ {-1, 0, 1}` scaled by a few magnitudes, and `alpha, beta ∈ {0, 1, 2, 3, 4}` (including `alpha == beta`, `alpha == 0`, `beta == 0`), compare the two implementations' returned `CGAL::Sign`. That is 3⁴ × 5² = 2025 combinations covering every branch including the merged one. Add a randomised arm over generated rationals for the general case.

Use the throwing `require` + `GateCheckFailedError` pattern. On mismatch the message must print the six inputs and both signs.

- [ ] **Step 2: Build and run — expect PASS**

Run: `pixi run -e default cmake --build <your-build-dir> --target sign_mixed_radical_equivalence_gate && <your-build-dir>/sign_mixed_radical_equivalence_gate`

Expected: `OK`, exit 0. **If it FAILS, stop immediately and report** — the copies are not equivalent, which makes stage 1 a bug fix with a live wrong-answer defect in one lane, and the plan must be re-designed before anything moves.

- [ ] **Step 3: Prove the gate is load-bearing**

Via a scratch copy only, perturb one implementation (e.g. return `CGAL::ZERO` in the `alpha == beta` branch) and confirm the gate fails with exit 1. Report the injection and result.

- [ ] **Step 4: Commit**

```bash
git add tests/native/sign_mixed_radical_equivalence_gate.cpp && git -c user.name="Jelle Feringa" -c user.email="jelleferinga@gmail.com" commit -m "test: prove the two mixed-radical sign copies agree" -- tests/native/sign_mixed_radical_equivalence_gate.cpp CMakeLists.txt
```

---

### Task 2: One definition in `src/exact/one_root.*`

**Files:**
- Modify: `src/exact/one_root.h`, `src/exact/one_root.cpp`
- Modify: `src/engagement_2.cpp` (delete its copy, call the module)
- Modify: `src/audit_exact_station_2.cpp` (delete its copy, call the module)

**Interfaces:**
- Produces: `CGAL::Sign compas_cgal::exact::sign_mixed_radical(const Rational& a, const Rational& b, const Rational& c, const Rational& d, const Rational& alpha, const Rational& beta)`

- [ ] **Step 1: Add the definition to `src/exact/one_root.*`**

Take copy A's body — it carries the derivation comments explaining *why* the reduction is exact (group over the shared root into `u, w ∈ Q(√α)`; sign follows from `sign(u)`, `sign(w)`, and `compare(u², βw²)`). Keep those comments; they are the reason the code is auditable. Adopt the repo's brace style from copy B. Use `exact::Rational` and `exact::OneRoot` — `exact::OneRoot` is `CoordNT` (bound by the existing `static_assert`), so this is a type-identity substitution, not a conversion.

Remove the stage-1 marker comment at `src/exact/one_root.h:32`.

- [ ] **Step 2: Point both call sites at it**

`engagement_2.cpp`: delete its definition; call `exact::sign_mixed_radical` at `:123` and `:141`.
`audit_exact_station_2.cpp`: delete `sign_mixed_radical_impl`; have `audit_sign_mixed_radical_exact` delegate. **Its signature must not change.**

- [ ] **Step 3: The equivalence gate must still pass**

It now compares the single definition against itself, which is tautological — so **delete that gate in this commit** and say so in the message. Its job was to license the move; keeping a tautology would be exactly the near-tautological-gate defect stage 0 found. Remove its `CMakeLists.txt` entry too.

- [ ] **Step 4: Full verification**

```bash
pixi run -e default exact-gates
pixi run -e default pytest tests/test_engagement_audit.py tests/adaptive -q -n auto
```
Report exact counts. Engagement and station results must be unchanged from before the move — state the before and after numbers, do not assert "unchanged" without both.

- [ ] **Step 5: Commit**

```bash
git add src/exact/one_root.h src/exact/one_root.cpp src/engagement_2.cpp src/audit_exact_station_2.cpp && git -c user.name="Jelle Feringa" -c user.email="jelleferinga@gmail.com" commit -m "refactor: one definition of the mixed-radical sign predicate" -- src/exact/one_root.h src/exact/one_root.cpp src/engagement_2.cpp src/audit_exact_station_2.cpp tests/native/sign_mixed_radical_equivalence_gate.cpp CMakeLists.txt
```

---

### Task 3: Close the test-oracle gap in the docs

The Python binding `_sign_mixed_radical` reached copy B via `audit_sign_mixed_radical_exact`, while the engagement geometry used copy A. Any Python test of that binding exercised code the geometry did not run. After Task 2 they are the same function, and that is worth recording.

**Files:**
- Modify: `docs/number_types.md`

- [ ] **Step 1: Document it**

In the one-root section, record: what the duplication was, that the Python-visible binding tested a different copy than the geometry used, that the copies were proven equivalent by an exhaustive 2025-combination gate before the move (naming the commit), and that a single definition now closes the oracle gap. Follow house style — conclusion first, no "BLUF" label, mkdocs-material admonitions, no ASCII art.

- [ ] **Step 2: Verify and commit**

`pixi run -e docs docs` must add no NEW warning (two pre-existing warnings from other sessions are expected; compare the sets).

```bash
git add docs/number_types.md && git -c user.name="Jelle Feringa" -c user.email="jelleferinga@gmail.com" commit -m "docs: one mixed-radical predicate, and the oracle gap it closes" -- docs/number_types.md
```

---

## Stage-1 completion criteria

1. `grep -c 'sign_mixed_radical' src/engagement_2.cpp src/audit_exact_station_2.cpp` shows call sites only — no definition in either.
2. `pixi run -e default exact-gates` exits 0.
3. Engagement and station suite results identical before and after, both numbers recorded.
4. `audit_sign_mixed_radical_exact` and `sign_mixed_radical_for_binding` signatures unchanged.
5. The equivalence gate is deleted, not kept as a tautology.
