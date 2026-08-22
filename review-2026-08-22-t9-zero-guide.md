# ETH-level assessment — `codex/exact-certified-adaptive-phase1-t9-zero-guide`

**Subject:** `compas_cgal`, exact-kernel port and extension of Held & Pfeiffer (2025),
*Trochoidal Tool Paths for Pocket Machining with Full Control of the Tool Engagement Angle*
**Date:** 2026-08-22
**Head:** `e0ae116` — 308 commits and +177,126 lines ahead of `main` across 391 files
**Method:** branch-topology analysis, five parallel source audits, and a full local build +
test run. Every headline number below was re-measured directly rather than relayed.

---

## Verdict

**The mathematics and the instruments are publication-grade. The enforcement is
disconnected, and the audit module fails in both directions — refusing to certify what it
does measure, and certifying what it never measured. None of it is a research problem;
all of it is about a week of work.**

The exact-kernel core is genuinely sound. Across 57k lines of C++ there are **zero**
`to_double` results reaching a certificate verdict, **zero** cross-root `Sqrt_extension`
arithmetic sites, **zero** `nextafter`/ulp/deflation constants, and **zero** unnamed
epsilon literals in the exact modules. That is the hardest thing in this project to get
right, and it is right.

What is not right is that **nothing enforces any of it**. The strict typing gate runs in no
CI job and is red at HEAD. The 699-test adaptive suite runs in no CI job and is red. The
schema gate points at a file that does not exist and reports success. The benchmark CI
smoke step is invalid grammar and would fail if it ever ran — and it never has, because
every workflow triggers only on `main`. **CI has never executed a single line of this
branch.**

And `engagement.py` — the audit module — fails in **both directions at once**, which is the
sharpest single observation in this review. Where it measures, it is too strict: at the
shipped 120° cap it declines to certify every circular motion *before measuring any
geometry*, because a fixed-fraction station density makes the growth guard exceed a full
turn. Where it does not measure, it is too lenient: `_unmeasured` returns
`cap_certified=True`, and `cap_violations` counts only the negation — so an operation that
was never analysed is indistinguishable from one proved compliant, and an empty toolpath
audits clean. The first defect errs safe and makes the headline claim unreachable; the
second errs *unsafe* and can produce a false "no violations" on exactly the third-party
geometry the benchmark harness exists to feed it.

**The suite does catch a false certificate — demonstrated, not inferred.** I injected one at
the C++ boundary, forcing `_stock_2.certify_segment_tea` to always report `certified=True`,
and diffed against a clean control:

| | result |
| --- | --- |
| control | **80 passed, 0 failed** |
| false certificate injected | **5 failed**, 75 passed |

The killers are exactly the tests you would want: `test_certificate_refuses_merge_over_cap`,
`test_certifier_gap_wiring_regression_guard`, `test_certify_flags_unclosable_guard_margin`,
`test_certify_segment_flags_slotting`, `test_certify_stations_refine_near_feature`. A second
injection at `audit_segment_tea_event_exact` killed nine more. **On the dimension that
matters most, this suite works** — and that is a stronger statement than any amount of
reading could support.

The prior due-diligence review (2026-07-29) concluded no complete pocket had ever been
produced and no performance gap closed. **Half of that is now out of date.** Complete
gouge-free paths come out on three pockets with 0.8% uncut; generation is measured at
72.4 ms against a 222.5 ms baseline, inside Held's stated 3–100 ms range. Certification
cost is measured *not* closed, and the project says so itself.

---

## The branch family

Nine `codex/*` branches exist. The topology resolves cleanly:

| Branch | Unique commits not in tip | Status |
| --- | --- | --- |
| `codex/adaptive-clearing-sp1-completion` | 0 | fully absorbed — deletable |
| `codex/exact-certified-adaptive-phase1` | 8 | superseded |
| `…-phase1-t6` | 4 | superseded |
| `…-phase1-t7` | 7 | superseded |
| `…-phase1-t7-runtime` | 20 | superseded |
| `…-phase1-t8` | 8 | superseded |
| `…-phase1-t10` | 8 | superseded |
| `…-phase1-t9` | 0 | strict ancestor of the tip |
| **`…-phase1-t9-zero-guide`** | — | **the tip** |

The six "superseded" branches look like orphaned work but are not: every file their unique
commits touch is present on the tip and larger — `segment_oracle.cpp` 548→648, `circle_oracle.cpp`
473→963, `event_trace.cpp` 235→337. The work was reimplemented rather than merged, so the
commits never landed while their content did. **Eight of the nine are safe to delete**; only
the tip matters.

Two test files shrank in the process, which is the shape of tests being dropped to make
something pass. They were not: `test_circle_oracle.py` replaced all 15 tests with 17
different, sharper ones (`test_exact_rational_probe_uses_the_circle_chart_not_binary64_sampling`,
`test_virgin_slotting_proves_cap_exceeded_without_sampling`). Coverage moved; it did not shrink.

### The frontier is split, and this is the most consequential structural finding

`jf/toolpath-redesign` is **not** behind the codex tip. The two have diverged from
merge-base `1860167`, with the codex tip 244 commits ahead and `jf/toolpath-redesign` 33
ahead. Neither contains the other, and **each holds a different half of the answer**:

- The **codex tip** has the adaptive pipeline, the benchmark corpus, the quality
  instrument, and the three generators.
- **`jf/toolpath-redesign`** has a hardened `engagement_2.cpp` — +981 lines on the single
  most load-bearing file — plus `tests/test_false_certificate.py` (+928) and
  `tests/test_growth_bound.py` (+953), neither of which exists on the codex tip.

This is not a preference. The exact-kernel audit independently flagged two defects on the
codex tip, and **`jf/toolpath-redesign` already fixes both**:

1. **`engagement_2.cpp:90` — the shared-root precondition is a `CGAL_assertion`.** The build
   is Release with `-DNDEBUG`, so CGAL maps it to `(void)0`. The sole guard on the number the
   entire cap certificate is computed from **does not exist in any shipped wheel**. Commit
   `148a49a` on `jf/toolpath-redesign` replaces it with a thrown `std::logic_error`, and its
   comment reaches the same conclusion the auditor did independently: *"an asserted form does
   not exist in any shipped wheel… the cap certificate would come out confidently wrong rather
   than absent."*

2. **`FULL_TURN_REPORTING_SLACK = 1e-12` rests on a false headroom claim.** It is an absolute
   angular slack justified by a scale-free argument, but the quantity it bounds carries
   `~ulp(|station|)/r` error. Commit `48b9e08` on `jf/toolpath-redesign` contains a *measured
   falsification table* of exactly this, and states it plainly: *"this constant is NOT nine
   decades above the per-span error, and any headroom claim of that shape is false."* It then
   supplies the telescoping argument that actually rescues the bound and renames the constant.

**The correct next move is to reconcile these two lines.** Assessing the codex tip in
isolation credits it with defects that are already solved 33 commits away.

---

## What meets the bar — with evidence

These are the strongest parts, and several exceed what the standard asks.

**Exact-kernel discipline is clean, and I verified the counts.** 80 `to_double` sites tree-wide;
68 are in Epick modules where `FT` *is* `double` so no conversion occurs; 11 are in exact
modules; **0 feed a certificate verdict**. Zero `nextafter`, zero ulp nudges, zero cross-root
arithmetic. `sign_mixed_radical` (`engagement_2.cpp:41-71`) decomposes a mixed radical into
three supported exact calls exactly as the doctrine's worked example prescribes.

**`UnavoidableEngagementWarning` is the sharpest design decision in the repo.** When the tool
first plunges into virgin stock a full slot is geometrically forced. A dishonest system hides
those circles or widens the cap. This one emits them, counts them, and classifies them into
`chain-entry` / `forced minimum advance` / `no loop radius can escape`. Errors run toward
refusing to certify, never toward certifying.

**`adaptive/units.py` is exemplary and implements the rule as written.** `NewType` for
`Millimetre`/`Radian`/`SquaredMillimetre`, a `WorldXY` phantom frame tag, `Point2(Generic[FrameT])`,
nine validated scalar wrappers each with `.build()`, and `@overload`-ed scalar-positional plus
sequence factories.

**The error model in new code is essentially perfect.** 88 named exceptions with a real
hierarchy, **zero** bare raises in `adaptive/` and the engagement trio, and **66 `except`
clauses with zero swallows** — every one translates and chains. 87 named C++ exception classes,
zero bare throws in the four largest new translation units.

**Tolerance discipline is excellent.** Across 24,282 lines of Python there are 4 float-exponent
literals, three of them in prose. The one real constant carries a two-line derivation.

**The certifier fails closed, and this is well fuzzed.** 14-way certificate-mutation
parametrization (delete-seam/root/cell/fibre, coalesce-roots, alter-multiplicity…) all landing
on `UNRESOLVED_DEGENERACY`; unresolved never degrades into "infeasible, try the next candidate".

**An 894-comparison differential exists** between the station refutation probe and the
independent `engagement_at` cap flag, with non-vacuity floors asserted.

**Zero `xfail`, zero `skip`, zero `skipif`, zero commented-out tests** across 92 files. Your
hardest test rule is honoured without exception.

**The benchmark corpus is a first-class deliverable.** ~42 instances across 6 families; a real
MATHSM protocol baseline using brute-force spacing search *because the relation is measured
non-monotone*; genuine ablations (probe-count convergence against an independent phase-offset
walk, climb-vs-conventional mirror, ladder-vs-bisection pinned by a test asserting they
disagree); `uncertified` split from `truly_exceeding` because conflating them once made a 2.4×
win read as a regression.

**The cost diagnosis in `docs/continuous_engagement_cost.md` is the best document on the branch.**
It runs a falsification order, kills its own leading hypothesis with measurement (string-mediated
rationals: 0.8 ms against a 17.5 s audit), supplies a *symbolicated* profile while documenting
the trap that `nanobind_add_module` strips symbols, and names the real open question —
a bivariate curve-pair analysis being used for a univariate root-isolation problem.

**Docs honesty is unusually strong.** `docs/refutation_soundness_crosscheck.md` contains an
explicit retraction: *"That is withdrawn. It was never measured."* `docs/radius_regulated_toolpath.md`
leads with a negative result and argues the negative is worth more than the improvement.

---

## Critical

**C1 — The audit module reports "not measured" as "certified", and the exemption is decided
by a float threshold and by labels the audited artifact itself supplies.** This is the only
finding that can produce a *wrong* answer rather than a missing one. The chain, verified
verbatim:

```python
# engagement.py:321-323
def _unmeasured(op_index, operation) -> OperationEngagement:
    return OperationEngagement(..., max_tea=0.0, cap_certified=True, stations=0)

# engagement.py:474
cap_violations = sum(1 for e in operations if not e.cap_certified)
```

Five call sites reach `_unmeasured`, and an unanalysed operation becomes indistinguishable
from a proved-compliant one. Three routes open the hole on real input:

- **The cut plane is inferred, by float, from the toolpath under audit.**
  `_infer_cut_height` returns `min(heights)` over all operations, and `_replay_line` then
  does `if z0 > cut_z + TOL.absolute: return _unmeasured(...)`. Nothing validates the
  single-depth contract the docstring assumes. On a two-depth program — or any third-party
  path whose lowest move is not its cut plane — `cut_z` binds to the lower plane and **every
  horizontal cut at the upper plane is never measured, never subtracted, and reported
  certified**. A full depth of slotting reports zero violations. The asymmetry is worse than
  it looks: `_replay_arc` applies no height test at all, so upper-depth *circles* are still
  measured — against a stock depleted by the other depth.
- **The exemption is keyed on a label the artifact supplies.** `if op.operation ==
  OperationType.RETRACT: return _unmeasured(...)` — the comment says "regardless of
  geometry". No `stock.contains`, no post-subtraction area delta. A generator that
  mislabels a cut-height move as `RETRACT` is certified for it *and* the move is not
  subtracted, so later operations measure against material that should be gone.
- **`TOL.absolute` (1e-9) decides certified-vs-not on the Python side** — a bare tolerance
  in a decision path, in the one module aimed at auditing geometry the project did not author.

The empty case is the same bug in its literal form: `max_tea = max(..., default=0.0)` and
`cap_violations = sum(...)` over an empty list yield `max_tea=0.0, cap_violations=0` — byte-
identical to a perfect path. Any consumer asserting `cap_violations == 0` passes vacuously.

*Sound fix:* `cap_certified` must be tri-valued (or `_unmeasured` must carry a distinct
`UNMEASURED` state) so `cap_violations` cannot absorb it, and the cut plane must be an
**input** with a loud error when an operation is off-plane — never an inference.

**C2 — The arc certifier cannot certify anything at the shipped cap, and gets worse as loops grow.**
`AUDIT_ARC_STEP_FRACTION = 0.05` is a fixed *fraction* of the turn — 20 stations per circle
regardless of size — so absolute spacing is `0.05 × 2πR` and grows linearly with loop radius.
The growth-lemma guard `gamma_guard = 2·GROWTH(spacing/2)` grows with it, unbounded.
`cap_guarded = tea_cap − gamma_guard`, and when that is ≤ 0 the motion is marked uncertified
**before any geometry is measured** (`engagement.py:288-307`). Measured directly:

| loop radius (× tool r) | `gamma_guard` | guarded cap @ 120° |
| --- | --- | --- |
| 0.5 | 109.4° | +10.6° |
| 1.0 | 166.2° | **−46.2°** |
| 2.0 | 259.1° | −139.1° |
| 4.0 | 419.2° | −299.2° |
| 8.0 | 730.9° | −610.9° |

Past ~2× tool radius the guard exceeds a **full turn**. The refusal count tracks how many
circles exist, not their engagement — which the repo states itself at
`tests/test_engagement_toolpath.py:113-118`. Committed corpus data agrees: uncertified sits at
**47–52% of operations on every instance**, including a neckless rectangle control.

*The fix already exists in this repo.* The C++ **segment** certifier does the right thing —
adaptive bisection with an absolute floor (`STATION_FLOOR_FRACTION = 1e-3`,
`CERTIFY_MAX_DEPTH = 24`). Arcs got a fixed-fraction Python mirror instead. This is porting a
working refinement loop, not new research. **Highest-value single fix on the branch.**

**C3 — `replay_certificate` has no success path.** Declared `-> ReplayCertificate` at
`replay.py:1188`, it contains **zero value-returning statements**; its last statement is an
unconditional `raise`. Verified by AST: `value-returns=[]`, `raises=[1240, 1276]`. All **24**
test invocations (16 in `test_replay.py`, 8 in `test_zero_guide_replay.py`) sit inside
`pytest.raises` — I counted them myself. There is no acceptance test because there cannot be
one. `mypy --strict` cannot catch this: an always-raising body is well-typed against any return
annotation. **The component whose entire job is to prove the reproducibility claim can only
reject.**

**C4 — CI has never run on this branch, and its gates are hollow.** All of `build.yml`,
`docs.yml`, and `benchmarks.yml` trigger only on `main`/`master`; `pr-checks.yml` is a CHANGELOG
check on `actions/checkout@v1`. Last fork run: 2026-07-26, on `main`. Additionally:
- `build.yml` runs **no tests at all** — `benchmarks.yml:5-6` says so in its own comment.
- The `benchmarks.yml` smoke step is invalid grammar. I ran it:
  `error: argument command: invalid choice: 'smoke' (choose from corpus, figure6, figures)`.
- `pixi run schema` targets `tests/adaptive/test_schema.py`, **which does not exist** — it emits
  `no tests ran` and exits 0. `jsonschema` is declared as a dependency and imported nowhere.
- `pixi run mutations-adaptive` targets a missing script.
- `tests/adaptive/` (699 tests) and `types-adaptive` and `lint` appear in **zero** CI jobs.

**C5 — The strict typing gate is red at HEAD, masking two latent `AttributeError`s.** I ran it:
6 errors in 3 files.
- `generator.py:731` — `.motion.end` on `ExactSegmentMotion | ExactCircleMotion`;
  `ExactCircleMotion` has no `end`.
- `generator.py:853` — `.source_commit_digest` on `TraversalCommit | RouteRetraceCommit`;
  `TraversalCommit` does not define it.
- `generator.py:570,571` — `int | None` into an `int` parameter.
- `replay.py:6`, `bootstrap.py:8` — unused `type: ignore`.

**C6 — The declared Python floor is false and would ship broken.** `requires-python = ">=3.9"`,
ruff `target-version = "py39"`, CI matrix Python 3.10 — but **22 adaptive modules do
`from typing import Self` at module scope** (127 usages), which is 3.11+. That is an immediate
`ImportError`, *not* deferred by `from __future__ import annotations` (only 5 of 31 modules have
it anyway), and `typing_extensions` appears nowhere. It is invisible because
`adaptive/__init__.py` is a bare docstring, so CI's `check_import` never reaches the submodules.
The doctrine-correct repair is to raise the floor to match `wheel.py-api = "cp312"` — never a
shim, which would be a forbidden fallback.

**C7 — The swept-prefix certifier measures no geometry.** `segment_oracle.cpp:584`:

```cpp
const bool pi_cap = cap_is_exact_pi(source);
const ContinuousTeaVerdict verdict =
    start_clear && pi_cap
    ? ContinuousTeaVerdict::CERTIFIED
    : ContinuousTeaVerdict::UNRESOLVED_DEGENERACY;
```

Two booleans. The verdict domain cannot even express `cap_exceeded`. Safety rests entirely on a
theorem stated in a comment, and **no test exercises the theorem** — the existing tests assert
digests, version strings, and stratum counts, with zero geometric content. This is on the
production hot path. It is not *wrong*; it is *unverified*, and it is the one place where an
unconditional "certified" would pass every existing test.

**C8 — Content-addressing is structural for inputs, aspirational at the native boundary.**
`InputIdentity` is genuinely strong — it hashes design rings, cut plane, tool radius, entry, all
six policies, and rejects component-version drift on construction. But `ComponentIdentity`, the
only primitive carrying `native_source_tree_digest` and `source_revision`, is **never constructed
in production** — its sole call sites are tests, one passing `b"\x00" * 32`. Production embeds
hand-maintained literals (`"event-exact-motion-oracle-v5"`). Reproducibility rests on a version
string someone must remember to bump.

**C9 — A known-over-cap path is structurally indistinguishable from a compliant one.** The
generators emit circles at positions the exact predicate refuses, and say so — but only
through a warning. `engagement_toolpath.py:970` returns `ToolpathResult(operations=...,
polyline=...)`, and those are the type's *only two fields*. The override qualification is
erased by `-W ignore`, a pytest `filterwarnings` entry, or simply crossing a process
boundary. `engagement_radial_toolpath.py` is partially repaired — its public
`RadialToolpathResult` carries `forced_loops` — but the internal `_SweepOutcome` carries
**both** `forced_loops` and `forced_advances`, and the public type drops the second, so
"advances taken past a refusing predicate" never reaches the caller in any generator.

The repo already argues this exact principle against itself
(`engagement_radial_toolpath.py:315-320`):

> *"Recovering the answer by parsing a warning string would be a worse API than a field,
> and re-deriving it by replaying the path would be a second implementation of the search."*

Correct, and applied to loops but not advances, and not at all in the other two generators.
Mitigant: `audit_toolpath_engagement` re-derives everything from geometry and never trusts
these fields, so an *audited* path is not fooled — only a consumer reading the generator's
own result.

---

## Important

- **The tangent-continuity test checks 7.7% of what it claims.** Its docstring says it asserts
  parallelism "at **all** consecutive engaged-operation transitions", but `_is_trochoid_bridge`
  exempts on *structure*, not geometry — same path, both `cut`, exactly one is a `Line`. In a
  trochoid that alternates Line↔Arc, that is nearly every transition. Measured on the test's
  own five fixtures:

  | fixture | engaged transitions | exempt | actually checked |
  | --- | --- | --- | --- |
  | square | 552 | 544 | **8** |
  | kite | 328 | 320 | **8** |
  | L_shape | 285 | 268 | 17 |
  | irregular | 753 | 702 | 51 |
  | star | 579 | 470 | 109 |
  | **total** | **2497** | **2304** | **193 (7.7%)** |

  The non-vacuity gate (`assume(len([...engaged ops...]) > 2)`) counts *operations*, not
  *comparisons* — off by roughly 70×. An implementation emitting arbitrary tangents while
  preserving Line↔Arc alternation passes. Fix: add a comparison-count floor, and make the
  exemption geometric (assert the bridge really is ~90° to the circle tangent) rather than
  structural. Related soft evasion at `test_toolpath.py:597` — `if not arcs: return`, the one
  genuine early-return in 1019 tests.
- **Hypothesis is used in 1 of 93 files, and points at the wrong module.** All three `@given`
  tests are in `test_toolpath.py`, exercising the *legacy* Epick generator; `tests/adaptive/`
  — 574 tests, the entire exact certifier — has none. The strategy generates evenly-spaced
  angles with ±25% jitter, which guarantees simplicity and therefore guarantees no necks, no
  pinches, no thin channels, no islands. It produces exactly the fat star-shaped polygons
  where trochoidal milling is easy. (`assume()` is *not* the problem — measured over 50
  examples it never fires once.) The deepest adversarial geometry in this repo is
  hand-authored, in the 12-case `event_corpus.json` — which is excellent, and is where the
  generated strategies should be aiming.
- **The accept side is thinner than the reject side.** Rejection is independently re-proved
  through a second native predicate; acceptance is not. The raster oracle validates
  `engagement_at(...)[0]` (`total_tea`) but not index `[2]` — the `cap_exceeded` flag the
  certificate actually rests on. Combined with the acknowledged grazing-incidence blind spot,
  residual false-certificate risk concentrates exactly there.
- **An emptiness guard is present on one oracle path and missing on its sibling.**
  `circle_strata.cpp:556-557` opens its `unresolved` disjunction with
  `partition.cells.empty() ||` — the correct pattern: an empty cell set is *unresolved*, not
  certified. `segment_oracle.cpp:488-492` computes the same flag with `std::any_of(...)`
  alone and **no `partition.strata.empty()` term**. If `strata` were empty the fold never
  runs, all flags stay false, and the path yields `CERTIFIED` while simultaneously reporting
  `whole_rim_disposition == "unresolved"` — a contradiction `event_trace.cpp:200-205`
  accepts as legal. Reachability is unconfirmed (it needs `cells` *and* `fibres` both
  empty), so this is latent rather than demonstrated. The sibling's explicit guard is the
  evidence that the authors consider the empty case possible. One-line repair.
- **Not one measured artifact is committed.** `docs/benchmarks/` and `build/benchmarks/` do not
  exist; `git log --all -- docs/benchmarks` is empty. Every headline number in the docs is
  hand-transcribed prose, and grepping `141.8`, `0.31`, `119.5`, `411.6` across `benchmarks/`,
  `tests/`, `src/` returns **zero hits**. Nothing detects staleness — and one number already is:
  `docs/engagement_controlled_toolpath.md:7` still leads with "224 ms… a factor of ~86" while its
  own evidence table at `:233` says 72.4 ms and ~27×. Off by 3×, in the lead paragraph.
- **`fig6_ours` and `fig6_diff_vs_held` have no generator and no data.** The Held reference
  series exist nowhere in the repo. This violates the repo's own written standard —
  `benchmarks/cli.py:14-16`: *"no figure in the docs comes from a script that is not in the
  repository."* The digitisation is honestly disclosed; the missing artifact is the defect.
- **The corpus exercises the wrong generator.** The 6-family sweep runs the *unregulated*
  generator, which has no cap input. The regulated generators — the actual contribution — are
  exercised on one pocket and three pockets. `engagement_rho_toolpath.py` is wired into nothing.
- **`_largest_admissible_advance` exists twice with the same name and opposite algorithms.**
  `engagement_toolpath.py:615` bisects under a stated monotonicity *assumption*;
  `engagement_rho_toolpath.py:373` scans downward *because that assumption is false*.
  `engagement_radial_toolpath.py` imports the **bisecting** one while its own docstring is built
  around "WHY A LADDER AND NEVER A BISECTION."
- **Units are typed in `adaptive/`, absent in the engagement trio** — 86 bare `: float` and 25
  `Tuple[float, float]`, importing `adaptive.units` not at all. `_Regulation` carries six
  physical dimensions at one type; `_GuideStation(cx, cy, radius, clockwise, tx, ty)` is a point,
  a vector and a length in six untyped scalars — and `engagement_rho_toolpath.py:312-338`
  fabricates one whose `tx`/`ty` hold a rotated radius direction rather than a tangent.
- **`engagement_toolpath.py` is a private-API shared library** — the other two generators import
  15 and 12 names from it, 11 and 7 of them underscore-private.
- **God objects on one axis.** Six adaptive modules fuse a per-run authority object, a per-step
  artifact plus serializer, and a stateless algorithm. `generator.py` 1465, `replay.py` 1276,
  `transaction.py` 1223 (its `CandidateEvaluator` is 719 lines / 19 methods). Nothing reaches the
  1.5k refactor-required line; nine files are ≥1k.
- **`segment_site_mat.cpp` ships 1,406 lines (40%) of test fixtures in the production library** —
  72 `_spike` symbols with no non-test caller. Its entire production surface is one 46-line
  function. A `segment_site_mat_spikes.cpp` already exists; the split was started and abandoned.
- **`continuous_tea_2.cpp`: `NB_MODULE` is 1,540 of 1,676 lines**, and two pieces of core logic
  leak into binding lambdas — a `[0,1]` station invariant enforced inside a lambda that a C++
  caller bypasses.
- **Factory contract inverted**: 163 raises live in `__post_init__`/`_validate` versus 25 in
  `build()`. Safe, because frozen dataclasses cannot bypass `__post_init__` — but the factory
  owns nothing, and several `build()` docstrings carry `Raises:` sections for exceptions they do
  not raise.
- **`__all__` present in `src/compas_cgal/__init__.py:32`** — the explicitly named prohibition.
- **`adaptive/medial_axis.py`: 0 of 64 public symbols documented**, including `MedialAxis.build()`.
  82 of its 97 raises use a single exception name.
- **Zero native tests cover TEA.** All 25 `tests/native/*.cpp` gates cover MAT/neck/graph;
  `src/continuous_tea_2/` (19.8k lines) and `engagement_2.cpp` have Python-binding coverage only.
- **No Hypothesis coverage of the certifier.** Hypothesis appears in exactly three tests, all in
  `test_toolpath.py`, none touching engagement, caps, certificates, or the event partition.
- **The full-circle oracle has no independent cross-check on `certified` verdicts.** The segment
  oracle has one (33-station dyadic falsifier against a disjoint implementation, applied to every
  certified verdict); the circle analogue exists as a binding but is called once, on empty stock.

---

## Evidence ledger

| Claim | Status | Evidence |
| --- | --- | --- |
| No epsilon in any cap decision path | **CONFIRMED** | 0 of 80 `to_double` reach a verdict; 0 cross-root arithmetic |
| Segment motions certified continuously | **CONFIRMED** | adaptive bisection, floor 1e-3, depth 24 |
| Circular motions certified | **REFUTED** | guarded cap ≤ 0 at 120° for loop radius ≥ 1× tool r |
| Engagement ≤ cap at each *evaluated* position | **CONFIRMED (narrow)** | exact predicate, integer bracket, no float tolerance |
| Engagement bounded *between* positions | **DENIED BY THE PROJECT** | "Between-position engagement is unbounded. This is the defining limitation." |
| Engagement non-monotone in spacing and radius | **CONFIRMED** | rung table 80.1→69.2→72.4→85.0; pinned by a test asserting bisection disagrees |
| Local zone depletion bit-identical to global | **CONFIRMED** | exact point-set equality, 26 mixed depletions × 2 orders |
| 32 probes is converged | **CONFIRMED** | K ∈ {3…48} ablation against an independent 60-position walk |
| Generation inside Held's 3–100 ms | **UNDER VERIFICATION** | 72.4 ms measured; no test pins it, no artifact committed, headline stale |
| Certification gap closed | **MEASURED NOT CLOSED** | 129 ms/segment sampled; ~15,000 ms exact; Held 3–100 ms per *pocket* |
| Path length 2–5× better than constant spacing | **UNDER VERIFICATION** | real baseline protocol exists; number is prose only |
| "Held grows 40–60×, ours 6.8×" | **UNSUPPORTED** | no generator, no data anywhere in repo |
| Cap honoured across the useful range | **REFUTED BY OWN MEASUREMENT** | max TEA 141.8° whether 20° or 140° requested |
| Unresolved rate on third-party geometry | **NEVER MEASURED** | the loader's own docstring calls it "the number that decides whether the certifier is a product" |
| Paths are worth running | **THE PROJECT'S OWN GATE SAYS NO** | 6/6 gate cases fail by design |

---

## The suite, as run

`pixi run baseline` — **1288 passed, 10 failed, 5m21s.**

Six failures are `tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running`
across 3 pockets × 2 generators. These are **deliberate**. The file's own docstring: *"THE
ASSERTIONS BELOW ARE WRITTEN AT THE STANDARD A CAM ENGINEER DEMANDS, NOT AT WHAT TODAY'S
GENERATORS ACHIEVE, AND MOST OF THEM FAIL. That is the point of the file."* It is a ~45-metric
instrument across five families, and it is the right instrument.

What it measures on `rect_12x8`:

| | |
| --- | --- |
| uncut fraction | 0.0080 |
| gouge free / rapid safety / continuity breaks | True / True / 0 |
| **max engagement** | **360.00°** |
| engagement p95 | 119.19° |
| **max engagement step** | **329.49°** |
| engagement variance | 5822 deg² |
| **recut fraction** | **0.8974** |
| degenerate loops (ρ ≤ r) | 4 |
| tangent breaks / curvature breaks | 24 / 20 |

The `329.49°` single step is the run-merge discontinuity from the literature, caught in
measurement rather than theory. To be precise about what it is *not*: it is a quality metric
comparing peak engagement across two toolpath-**adjacent** cut motions with stock depletion
between them — two valid measurements against two different stock states, not an unsampled
gap in the certifier. The shared endpoint is an endpoint station of both motions' recursion,
and any run above π trips `run_exceeds_cap`'s unconditional negative-orientation branch. So
the certifier catches it; the step is the generator's behaviour being reported honestly. The
89.7% recut and 24 tangent breaks are why the gate says no.

The other four failures cluster in `adaptive/generator.py`, which has uncommitted edits, so they
read as work in progress. One deserves attention independently:
`test_real_active_family_stops_at_unresolved_exact_event` fails `DID NOT RAISE
UnresolvedMotionEventError` — a fail-loud guard that is not firing. Another,
`test_task13f_full_continuation`, is stale rather than broken: it expects
`cap=0; gouge=56` where the code now reports `cap=10; gouge=46`.

**A standing tension worth naming:** the honest aspirational gate and the green-suite invariant
are in conflict, and right now the invariant lost silently. With `xfail` correctly forbidden,
there is no sanctioned way to say "known-red gate", so `pixi run baseline` can no longer
distinguish a known gap from a new regression. That needs a deliberate decision — a separate
`pixi run gate` task outside `baseline` is the obvious shape.

---

## Recommendations, highest impact first

1. **Reconcile `jf/toolpath-redesign` with the codex tip.** The frontier is split and each side
   holds a different half. Two of the four exact-kernel defects found on the tip are already
   fixed there, along with 1,881 lines of certifier tests the tip does not have. Do this before
   anything else, because it changes what the other items apply to.

2. **Close the vacuous-certification hole (C1).** Make `cap_certified` tri-valued so
   "unmeasured" cannot be absorbed into "compliant", and take the cut plane as an input
   rather than inferring it by float from the artifact under audit. This is the only defect
   that produces a wrong answer today, on inputs the benchmark harness is designed to feed it.

3. **Fix the arc station density (C2).** Replace `AUDIT_ARC_STEP_FRACTION` with the adaptive
   bisection the segment certifier already uses. This is the single change that moves
   `cap_violations` from structurally-impossible to meaningful, and it is porting working code.

4. **Wire the gates that exist into CI.** Change the workflow triggers off `main`-only, add
   `tests/adaptive`, `types-adaptive`, and `lint` as jobs, fix the smoke step's grammar, and
   delete or implement `schema` and `mutations-adaptive`. Everything needed already exists and
   is disconnected. Until this lands, every rigor claim is asserted rather than enforced.

5. **Fix C5 and C6 while wiring** — the 6 typing errors include two latent `AttributeError`s, and
   the Python floor must move to `>=3.11`/`>=3.12` to match what the code imports.

6. **Implement `replay_certificate`'s success path (C3) with an acceptance test.** Until then the
   reproducibility claim has no proving instrument.

7. **Commit one measured artifact.** The harness is first-class and has never been
   run-and-recorded. One committed corpus run turns fifteen prose claims into checkable ones and
   catches the already-stale 3× headline automatically.

8. **Give the swept-prefix theorem (C7) a geometric test, or gate it behind the measured
   certifier.** It is the one place where an unconditional "certified" passes every existing test.

9. **Delete the eight superseded branches** and `git worktree prune` the eight prunable worktrees.

---

## On the two theses

The prior review argued the certifier is independently saleable and the stronger of the two
businesses. **The evidence now cuts against that in the near term and for the certifier
specifically.** The exact continuous certifier is measured at ~15 s per *segment motion* against
Held's 3–100 ms per *entire pocket* — five orders of magnitude — and the project's own analysis
concludes the two available levers are insufficient, naming the real open question (bivariate
analysis applied to a univariate problem) as unanswered research.

What has become saleable instead is the thing nobody planned: **the quality instrument.** A
~45-metric machining gate with derived rather than tuned thresholds, honest `uncertified` vs
`truly_exceeding` separation, exact coverage measurement that *raises rather than answers* when
the grid is too coarse, and figures generated from live pipeline output with captions measured at
draw time. That is a verification instrument for third-party CAM output, it works today, and it
is the strongest artifact on the branch — while being almost entirely unplanned work sitting
uncommitted in the working tree.

The generator thesis, meanwhile, is in better shape than the July review could see: complete
gouge-free paths, generation inside Held's range, and a precisely-diagnosed reason the cap claim
does not yet hold. Fix C1 and that claim becomes testable for the first time.
