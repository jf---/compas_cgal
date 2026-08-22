# Technical Due Diligence — exact-certified trochoidal pocketing

**Subject:** `compas_cgal` — exact-kernel port of Held & Pfeiffer (2025), *Trochoidal Tool Paths
for Pocket Machining with Full Control of the Tool Engagement Angle*
**Date:** 2026-07-29
**Reviewed:**

| Ref | Branch | Head | Scope |
| --- | --- | --- | --- |
| **A** | `jf/toolpath-redesign` | `33cbcb6` | integration line — sampled TEA certifier, exact stock model, straight-skeleton generator |
| **B** | `codex/exact-certified-adaptive-phase1-t9` | `20f2b69` | Phase-1 execution — continuous TEA certification, segment-site MAT, adaptive generator |

**Method:** source read of `src/{toolpath,stock_2,engagement_2}.{cpp,h}`, `src/continuous_tea_2/`,
`src/segment_site_*`, the `compas_cgal.adaptive` subpackage, `tests/`, and the full MkDocs corpus.
Branch B inspected read-only via `git show`/`git diff` against merge-base `1860167`.
Working tree not modified.

---

## Recommendation

**Fund the mathematics. Refuse to fund the volume. Make the next cheque contingent on CI
turning green and one whole pocket coming out the other end.**

The intellectual asset is top-decile and, in two specific respects, ahead of the published state
of the art. The commercial asset does not exist yet: after 187 commits and 114k lines on
branch B, the pipeline has not produced a complete certified toolpath for a single pocket, and
the performance gap that motivated the whole rewrite has not been shown to be closed.

The binding constraint has moved. Six weeks ago it was *"can this be made fast enough?"* Today
it is *"can this be made reviewable, enforced, and finished?"* — an organisational problem, not
a mathematical one. Structure the investment accordingly.

Downside protection is unusual and worth naming: even if the generator never ships, **the
certifier is independently saleable** as a verification instrument for third-party CAM output.
Very few people can build that; nobody is selling it. Developed as a standalone thesis in
[Alternative thesis: the certifier as the product](#alternative-thesis-the-certifier-as-the-product)
below — on the current evidence it is the stronger of the two businesses, and it is reachable
sooner.

---

## The claim, and whether it holds

Held & Pfeiffer produce trochoidal pocket paths with a user-set maximum tool-engagement angle
(TEA). They **assert** the bound — float bisection to ε = 0.001 plus a conservative reduction at
necks. This project **certifies** it: the cap verdict is an exact predicate over an exact
kernel, not a float comparison.

That distinction is not academic posturing. For the buyer — unattended machining, tool-life
warranty, aerospace process qualification — *"provably ≤ θ_max on every cutting motion"* is a
liability artifact. *"Bisected to 1e-3 and reduced a bit at necks"* is a best effort. It is a
genuine, defensible differentiator and I found no prior art selling it.

**The claim holds at the level it is made.** Verified by reading the source, not by trusting the
documentation:

| Claim | Evidence | Status |
| --- | --- | --- |
| Cap decision is exact; no epsilon in any decision path | `engagement_2.cpp:100–139` — `run_exceeds_cap` decides via `sign_mixed_radical` on `A + B√α + C√β + D√(αβ)`; no `to_double` reaches a branch | **Confirmed** |
| No undefined cross-root arithmetic | `:58–69` — builds `u, w ∈ ℚ(√α)`, compares `u²` against `β·w²` (same-root); `√β` enters only as the rational `β` | **Confirmed** — the subtle one, and it is right |
| Run adjacency exact, not tolerance-based | `:269–282` — merge on `==` of one-root arrangement points; wrap-around handled | **Confirmed** |
| Transcendental intent handled by contract change, not epsilon | cap crosses the API as `4·sin²(θ/2)`, validated to `(0, 4]` at `:565` | **Confirmed** |
| Doubles confined to reporting | TEA reported from *true* runs, decision taken over *pessimistic* runs (`:286–310`) | **Confirmed** |

---

## Branch A — what convinced me this is real engineering

### A soundness hole that only exists if you reason about your own certificate

The growth lemma bounds TEA change as `O(√d)` over centre travel `d`. But when two engaged runs
fuse across a closing rim void, `max_run` jumps by `O(1)` **within an arbitrarily small step** —
unbounded by any `√d` lemma. Most teams never notice, because most teams have no certificate to
be unsound.

This team noticed, and — critically — **did not fix it by inflating a constant.** The repair is
gap-closure pessimism (`engagement_2.cpp:150–224`): pre-absorb every void gap ≤ γ at each
station using the *identical* exact predicate, with a written contradiction proof and a stated
safe-failure direction (`:437–461`). Reported numbers stay true; only the decision is inflated.

That is the difference between a research artifact and something you can sign your name to.

### The differential oracle shares no geometry code with the thing it checks

`tests/test_engagement_oracle.py` builds an independent NumPy raster occupancy model, drives it
through the same subtraction sequence, and cross-checks both TEA and growth-bound soundness with
a tolerance derived from the oracle's own discretisation. Correct evidence standard, rarely met.

### Performance work is measured, A/B'd, and honestly bounded

`docs/engagement_zone_query.md` reports **4.8× end-to-end / 12.4× per query** for the zone-query
swap — obtained by checking out the *pre-swap source*, rebuilding it, and replaying identical
inputs. Verdict fields are bit-identical across both builds. The residual ≤ 1e-15 rad divergence
is correctly diagnosed as a `to_double`-of-representation artifact on algebraically-equal
one-root points — reporting, never deciding. The page then states plainly that the win *moves*
the bottleneck rather than removing it.

I did not find a single inflated number in the documentation.

### `docs/exactness.md` is publishable as it stands

The best statement of exact-kernel discipline I have read outside CGAL's own manual, anchored to
real line numbers in this repository, opening with the incident that caused it — a
`1.0 - 1e-12` deflation constant the team **removed** rather than justified. A codebase that
documents its own near-miss with that candour earns credit on its other claims.

---

## Branch B — what six days bought

**187 commits · 260 files · +114,263 / −609 · 2026-07-24 → 2026-07-29 · single author.**

### Stronger in rigour, catastrophically weaker in cost: the continuous event partition

The headline. Branch A's certifier bisects stations, applies a Lipschitz-style growth lemma with
`TEA_GUARD_SAFETY_FACTOR = 2`, and gives up at a refinement floor — returning *uncertified* when
the guard cannot close. Branch B's `continuous_tea_2/` (≈40 translation units) instead
decomposes the motion parameter into **sign-invariant cells separated by exact algebraic event
fibres** — tangencies, support overlaps, trimmed-vertex passages, pair-orientation boundaries,
cap equalities, seams — and decides the cap per cell, exactly.

This deletes the entire approximation apparatus at once: no guard constant, no safety factor, no
refinement floor, no spacing-exhaustion failure mode.

> *"It is a proof boundary, not a dense station audit."* — `docs/continuous_engagement.md`

Against Held & Pfeiffer — who bisect floats to ε = 0.001 — that is, as a *statement of
guarantee*, two levels stronger rather than one.

!!! danger "And it is unusable at that price — measured 2026-08-20"

    | certifier | per segment motion |
    | --- | ---: |
    | `engagement_2.cpp` — sampled + growth-bound guard, zone query | **129 ms** |
    | `continuous_tea_2` — bivariate algebraic curve kernel | **~15,000 ms** |
    | Held & Pfeiffer — an *entire pocket* | **3–100 ms** |

    The continuous oracle is **~116× slower than this project's own existing certifier**, which
    was already too slow, and it degrades **superlinearly** with stock complexity (0.12 s at 4
    stock vertices → 14.7 s at 20, on a fixed probe segment). t9 therefore trades roughly four
    orders of magnitude of performance to remove a conservative guard and an `uncertified`
    verdict.

    **This is not a tuning problem.** A symbolicated profile attributes ~90% of the cost to
    CGAL's bivariate algebraic curve kernel — `Curve_pair_analysis_2` (27k samples),
    `Curve_analysis_2` (25k), `Algebraic_curve_kernel_2` (18k), Bitstream-Descartes (~12k),
    `Shear_transformation` (9k), resultant + modular GCD (~15k). Two candidate remedies were
    measured or costed: swapping the deliberately-pinned boost multiprecision backend for GMP
    buys a small multiple; eliminating rational-to-decimal-string round trips buys 5–15%
    (falsified as the cause — removing all 50 round trips left timings unchanged). Neither
    approaches five orders.

    **One untested variable, stated because it qualifies the verdict.** Every figure above was
    measured with GMP **disabled** — `CMakeLists.txt` sets `CGAL_DISABLE_GMP`,
    `CGAL_USE_GMPXX OFF` and `CMAKE_DISABLE_FIND_PACKAGE_GMP`, `src/exact_algebraic_1.cpp` guards
    it with `#error`, and `otool -L` confirms neither GMP nor MPFR is linked. All of the bignum
    work therefore runs on boost `cpp_int`. The process abort is itself a `cpp_int` artifact:
    `eval_convert_to<double>` accumulates with `ldexp` and rounds through
    `boost::math::float_next`, which rejects non-finite input with precisely the observed message.
    Whether linking GMP moves the runtime, the 6.3 GB peak, or the abort is **untested**. It is a
    cheap experiment and a genuine packaging decision (portable wheels, no external dependency),
    and it should be run before the cost verdict is treated as final.

    **Consequence for the thesis.** The guarantee this branch adds is real and the cost is
    disqualifying, so the two cannot be reported as one result. On current evidence the
    shippable certifier is the *sampled + guarded* one in `engagement_2.cpp`, which is within
    one to two orders of Held rather than five, and whose conservatism is an explicit, bounded,
    documented pessimism rather than an approximation. Reviving the continuous partition
    requires a different formulation — the events along a one-parameter motion are univariate,
    so the use of a **bivariate** curve-pair kernel is itself the open question — not further
    optimisation of this one.

### Strictly stronger: straight skeleton → true segment-site medial axis

Branch A honestly renamed `polygon_medial_axis_transform` → `polygon_skeleton_clearance` because
a straight skeleton *is not* the MAT except on convex polygons — a real correctness gap against
Held, who uses the true MAT. Branch B closes it: `segment_site_*` (≈20 TUs) builds the MAT from
CGAL's Segment Delaunay Graph plus Voronoi adaptor, with exact neck classification and
clearance. The rename was the honest interim; this is the fix.

### Three-valued verdict, with `unresolved` load-bearing

> *"Unsupported or unreconstructed degeneracy is unresolved, never sampled into acceptance."*

Certified / cap-exceeded / unresolved, with proof of physical exceedance taking precedence over
proof incompleteness. The correct failure model for a safety artifact, enforced structurally
rather than by convention.

### Two lines that indicate a real engineer at the wheel

> *"A digest beside a proof is not proof ownership."*

The deciding-partition digest must live **inside** the canonical trace record, or the same trace
digest could be relabelled onto a different deciding partition. That is a provenance attack,
correctly identified and closed, against a *geometric* proof. I have not seen the thought in CAM
before.

> *"Running the same deterministic algorithm twice is not independent evidence. It only proves
> that two invocations agreed."*

The stated reason an entire Task-3 spike was **rejected and rewritten**. Killing your own work
for that reason is the behaviour being underwritten here.

### They wrote their own licensing blocker

`docs/licenses/cgal-adaptive-package-audit.md` — status **`CGAL_LICENSE_GATE_ERROR`**. It
enumerates every instantiated CGAL package (2D Arrangements, Boolean Set Operations, Segment
Delaunay Graphs, Voronoi Adaptor, Apollonius) as GPL-3.0-or-later-or-commercial, declares Tasks
1–16 blocked, and closes:

> *"No sampled geometry, different kernel, package substitution, or silent fallback satisfies
> this gate."*

Documented better than I would have stated it. **Still unresolved, with 187 commits on top of
it.** Development is not distribution, so proceeding is legitimate — but this is a hard gate
between the code and any revenue, and it is not being burnt down.

### Discipline held under volume

Zero `xfail`, zero `skipif`, zero `TODO`/`FIXME`/`HACK`/`XXX` across all new source and all
`tests/adaptive`. At 114k lines that is not an accident; it is house rules surviving contact
with scale.

---

## Risks, ordered by effect on the cheque

### 0. The continuous certifier cannot decide cap-violating motions, and one failure mode kills the process

A witness-based refutation probe added to the candidate search refuted **31** link motions across
two real fixtures. Three implementations sharing no geometry code agree these motions violate the
cap. On the same motions `continuous_tea_2` does not return a verdict at all: on one fixture it
raises `UnresolvedMotionEventError`, on the other it **aborts the process**
(`boost float_next<double>: Argument must be finite, but got inf`), and on that family it runs to
6.3 GB RSS at 98% CPU.

The verdict distribution, measured directly — probe disabled, a full native audit forced for
every candidate surviving gouge containment, the search prevented from stopping early so the
whole family is observed (branch family, 125 audits):

| probe verdict | partition outcome | count |
| --- | --- | ---: |
| silent | `certified` | 86 |
| **refuted** | `cap_exceeded` — agrees | 24 |
| **refuted** | `IncompleteSegmentPartitionError` — gives up | 15 |
| **refuted** | **`certified`** — would be unsound | **0** |

!!! warning "A stronger claim was made here and is retracted"

    An earlier version of this section stated that `certify` returned a `MotionWitness` — a
    certificate — for these motions, i.e. a **soundness** hole. **That is false and is
    withdrawn.** It was inferred from a *recorded* test expectation (`cap=0; gouge=56`) that no
    longer reproduces on the current tree, which crashes at that point; a stale expectation
    cannot testify to what the code did. Measurement shows the partition never certifies a
    refuted motion. **No certificate was ever wrong and no bad toolpath was emitted.** The defect
    is incompleteness, not unsoundness, and the distinction is material enough that the original
    wording would have been an unfair characterisation of this codebase.

    It follows that the refutation probe's contribution is **completeness and cost**, not a
    soundness repair: it decides the ~⅓ of refuted motions the partition cannot decide at all,
    and decides the remainder far more cheaply.

What is established: the certifier is **incomplete and fragile exactly where the geometry is
hard**, and one of its failure modes is a process abort rather than a diagnosis. For a component
whose entire purpose is to decide, that is disqualifying on its own, independent of the
[cost verdict](#6-performance-now-measured-and-disqualifying-for-the-continuous-oracle).

**Three implementations sharing no geometry code agree**, per-station, at
`cap_chord_ratio = 4.0` and `gap_close_ratio = 0.0` (pessimistic runs equal true runs, no
inflation):

| check | result |
| --- | --- |
| `max_run_tea > π` — the per-run measure the cap actually bounds | **31 / 31** |
| `engagement_at` → `cap_exceeded` (exact per-run predicate on the exact arrangement) | **31 / 31 True** |
| `certify_segment_tea(cap_radians = π)` → `cap_certified`, whole-motion certifier | **0 / 31**, over 139–274 stations each |
| any station at 2π (i.e. fully buried rim rather than a cap violation) | **0** |

The worst case is a **single engaged run of 5.150 rad = 295° against a 180° cap**. `total_tea`
equals `max_run_tea` in all 31 rows — one engaged run at every witness station — so no
multiple-run artifact is available as an explanation.

!!! note "A retracted intermediate claim, recorded because it is the trap"

    The first cross-check quoted `total_tea` — a **sum over disjoint runs** — against a **per-run**
    cap, and reported "over cap by 8–115°". Those figures were withdrawn and re-measured on
    `max_run_tea`. They happened to coincide here, but only because every witness station carried
    a single run, which could not have been known in advance. Two disjoint runs of 1.70 rad give
    `total_tea = 3.40` while neither violates π. `engagement_at` returns
    `(total_tea, max_run_tea, cap_exceeded)`; the middle field and the flag are the load-bearing
    ones.

Independently re-verified here on synthetic geometry: 132 predicate comparisons between the
station probe and `engagement_at`'s `cap_exceeded` flag agree in **128**, establishing that the
probe's notion of violation is the sampled oracle's exact per-run predicate rather than a
sampled approximation of it.

No bad toolpath escaped on these fixtures: each violating link happened to be rejected by the
*following* circle's containment test. Nothing structural guaranteed that.

The defect sits in `continuous_tea_2` exact event discovery — the same code that `std::terminate`s
on this family with `boost float_next<double>: Argument must be finite, but got inf`, and that
runs to 6.3 GB RSS at 98% CPU on it. It is **open**. Combined with
[the cost verdict](#6-performance-now-measured-and-disqualifying-for-the-continuous-oracle), the
component is both the slowest and the least trustworthy part of the branch, and the sampled +
guarded certifier is the one that was **right** as well as the one that is fast.

The 4 remaining disagreements in the re-verification are all at station `0/1` and all in the safe
direction — the probe declines to refute where the oracle's flag fires. A second, independent
differential on a different corpus reproduces this exactly: **894 comparisons, 543 refutations,
zero cases of the probe refuting where the exact flag says otherwise**, and all 36 gaps at `0/1`.

The mechanism is worth stating precisely, because the obvious reading is wrong.
`start_disk_has_no_material_interior` (`segment_oracle.cpp:200`) is called only from the
swept-prefix theorem at `:571`; the station predicate reaches `classify_station_cell`
(`station_classifier.cpp:72–132`), which has **no start-station special case at all**. So the
`0/1` asymmetry is not a guarded carve-out — it falls out of branch construction at the segment
start. That is a *weaker* guarantee than a deliberate carve-out, which makes pinning it more
necessary rather than less; it is now enforced by two tests carrying non-vacuity floors (≥500
rows, ≥100 refutations) so they cannot silently go hollow.

!!! note "This ordering assumes track G — read it differently for track V"

    The ranking below is written for the generator thesis. Under track V (see
    [Proposed terms](#proposed-terms)) the weights shift materially, and a reader evaluating that
    track should re-sort before drawing conclusions:

    | Risk | Track G | Track V | Why it moves |
    | --- | :---: | :---: | --- |
    | 1 · No CI enforcement | 1 | **1** | Universal. The enforcement gap is thesis-independent. |
    | 6 · Performance | 6 | **2** | Real NC programs carry 10⁵–10⁶ blocks against 271 ops, and verification must simulate stock from block 1 — so verification cost *is* depletion cost. The dominant technical risk on V. |
    | 8 · `const_cast` blocks parallelism | 8 | **3** | V's conservative-screen phase parallelises trivially across blocks; G's interactive loop does not. The shared-arrangement landmine binds sooner. |
    | 2 · Volume / 3 · Unreadable spec | 2 / 3 | **4 / 5** | V's load-bearing surface is `continuous_tea_2`, `stock_2` and the segment-site MAT — roughly half the branch. The generator machinery (`adaptive/` spacer, traversal, transaction, entry) need not be made reviewable at all. **Track V makes Tranche 0 cheaper.** |
    | 5 · No complete pocket | 5 | **—** | V1 does not generate paths. A missing generator is close to irrelevant. |
    | 7 · CGAL licensing | 7 | **↓** | May be sidestepped entirely by hosted delivery (see the licence note under Proposed terms); it *rises* under G, which must convey. |

    **The largest risk on either track is not in this list, because it has never been
    measured:** the `unresolved` rate on real third-party geometry. Everything enumerated below
    is an observed defect with a known fix. That one is an unknown with an existential answer,
    which is why Tranche 0 buys the number before anything else.

### 0b. The test suite does not encode toolpath quality — measured 2026-08-22

1,154 passing tests, a clean `mypy --strict` gate, and genuinely publication-grade exact-kernel
discipline. And the generator plunges into pocket corners and cuts straight out of them.

Every defect found during this review was found by **ad-hoc measurement, never by the suite**: a rim
span reported as 720°, the engagement cap inert below ~140°, the generator probing 3 positions of a
360° loop and blind to its own violations, four degenerate "trochoids" of radius 0.02–0.06 against a
1.0 tool radius — plunges, carrying the worst engagement in the path — and a 1.96-long straight link
slotting into a corner at 106.5°. Dropping all four degenerate loops changes coverage by **0.000%**:
they remove nothing and are pure tool wear.

What the suite asserts is exact predicates, type contracts, canonical digests and replay provenance
— all correctness *of the machinery*. None of it can notice that the output is a poor toolpath.
Green currently means "the machinery is self-consistent", not "this path is fit to cut with".

**Worse, the two headline claims rest on proxies.** Engagement angle is a proxy for load; the physics
is maximum undeformed chip thickness, and below a material-dependent minimum the edge *rubs* rather
than cuts — which is harder on a tool than a heavier cut, and which no metric here expresses. Path
length is a proxy for cycle time; feed is bounded by curvature (`v ≤ √(a_max/κ)`) and by jerk, so the
measured "0.31× length" is **not** a 0.31× time claim and could even be slower on a real machine.

A four-group quality gate — elementary validity, cut mechanics, machine-executable speed, tool
longevity — is being added, with thresholds set at what a competent CAM engineer would demand rather
than at what the current code achieves. **It is expected to fail on arrival**, and that is the
deliverable: a gate that passed today would launder a poor toolpath into "all tests pass".

For diligence purposes the implication is narrow and important: the exactness work is real and
verified, and it says nothing whatever about whether the machine output is good. Those are separate
claims and only the first currently has evidence.

### 1. None of branch B runs in CI — the highest-value, lowest-cost fix on this list

`git diff` on `.github` is **empty**. 32 new Python test modules, ~20 native C++ gate files, and
`mypy --strict` — the declared gate for the entire typed-units design — and the workflows are
untouched. Worse, every native gate is `add_executable(... EXCLUDE_FROM_ALL)`: it does not even
build unless requested by name.

The branch's verification apparatus is, as an enforcement mechanism, **inert**. A proof system
nobody is required to run is documentation.

### 2. Six days, 114k lines, one author, no second reader

≈19k lines/day — agent-generated at industrial scale, and the artifacts show it:

| Artifact | Size |
| --- | ---: |
| `src/segment_site_mat.cpp` | 3,503 lines |
| `src/segment_site_endpoint_binding.cpp` | 3,403 lines |
| `tests/native/task9_segment_segment_gate.cpp` | 4,183 lines |
| `docs/segment_site_mat.md` | 4,302 lines / 258 KB |
| Phase-1 plan | 111 KB |

The project's own standard is that 1,500 lines demands refactoring. Nothing here has been read
end to end by a human, and at this rate nothing can be. **The bus factor did not improve — the
amount riding on the bus went up 20×.**

### 3. The primary specification has drifted out of readable English

> *"At that circle's phase seam, incomparable exact active sets prove a mixed transition even
> though tangent incidence evidence remains inactive on both adjacent cells."*

Compare `docs/exactness.md`, which is genuinely publishable. `docs/continuous_engagement.md`
reads as a task-completion log addressed to whoever wrote the previous entry. Volume up ~10×,
communicative value down. Not a style note: this is the mechanism by which key-person risk
becomes irreversible.

### 4. Plan scaffolding has leaked into permanent artifact names

`task9_*.cpp`, `task13f_fixture.py`, *"Task 11A now independently reconstructs…"*. The numbering
of a transient plan document is now load-bearing in filenames and in the specification. It will
outlive the plan and mean nothing to the next reader.

### 5. There is still no complete pocket

Read the maturity admonition literally. The pipeline certifies a **two-circle L-pocket prefix**;
replay *"fails at the explicit nonterminal-traversal boundary"*; the real L family *"exposes an
unresolved segment-partition case before global completion"*; a constant-clearance arm needs a
link-only advance that is not built. The closing sentence:

> *"Traversal/coverage closure, fresh terminal replay, arbitrary-pocket evidence, and matched
> Held–Pfeiffer performance remain incomplete."*

Real credit for that paragraph — most teams would have claimed the milestone. The commercial
reading is nevertheless unchanged: **the product does not exist, and six days of extraordinary
output did not change that.**

### 6. Performance: the generator reached Held's range; the continuous oracle remains disqualified

!!! success "Generator parity — measured 2026-08-21"

    The engagement-controlled generator went from **2553 ms to 71.9 ms** on a 12×8 pocket (35×),
    and 4156 ms → 171.3 ms on an L-shaped pocket (24×), with the produced path **byte-identical at
    every step** — 49 cuts / 441.8 length / 5 true exceedances, and 112 / 595.3 / 7. The rectangle
    now sits **inside Held & Pfeiffer's published 3–100 ms range**.

    Three representational fixes did it, none requiring a kernel change, a licence decision, or a
    change to path geometry:

    | fix | effect |
    | --- | --- |
    | Exact **annulus** for trochoid loops — the swept region of a disk about a circle is two circles, not 470 chained disks | 42.9 → 2.2 ms per cut; `engagement_at` 6.5× faster as a knock-on, since the arrangement stops accumulating |
    | **Local zone depletion** — `Arrangement_zone_2` instead of a global boolean difference | flat in arrangement size where the global difference grows 5.3× over 96 cuts; `exactly_equals` verified across 30 mixed depletions |
    | **Quad capsule** for links — two exact end disks plus a rational quadrilateral at half-width `r(1−f)` | 4.19 → 0.307 ms per call (13.6×) |

    The last one corrects a diagnosis held for much of the investigation. The *exact* capsule needs
    irrational side lines and genuinely is not representable in the circle-segment traits — but the
    disk chain never represented it exactly either; it under-covered by a documented budget. The
    right question was not "how do we represent the capsule exactly?" but **"what is the cheapest
    representable region that under-covers within the slack already promised?"** Independently
    verified across 22 segments and 3,080 point checks: **zero over-cover violations, zero budget
    violations.**

    Caveat on the comparison: Held's 3–100 ms is one pocket, no hardware stated, no instance table.
    "Parity" here means same order and overlapping range, not a like-for-like benchmark. Making it
    rigorous is what the benchmark corpus exists for.

!!! danger "But the engagement cap is currently inert below ~140° — measured 2026-08-21"

    Speed is not the whole product. The benchmark corpus, on its first end-to-end run, found that
    the engagement-controlled generator **does not deliver the cap it is asked for**:

    | cap requested (deg) | 20 | 40 | 60 | 80 | 100 | 120 | 140 |
    | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
    | max TEA after entry (deg) | 141.8 | 141.8 | 141.8 | 141.8 | 140.7 | 141.8 | 142.3 |
    | circles on a forced minimum advance | 264 | 72 | 40 | 20 | 0 | 0 | 0 |

    Asking for 20° yields the same measured maximum as asking for 140°. Isolated further: loop
    engagement is **171.8° at both a 40° and a 120° cap**, bit-identical, while links behave
    sensibly (~105°). The cause is structural — the generator regulates the **advance** between
    machining circles but not the **trochoid radius**, and radius is the dominant term. Measured at
    one station with stock and position held fixed, engagement runs from **0° at a small radius to
    114.5° at the radius the generator actually picks**; it is maximising the loop rather than
    regulating it.

    Two consequences worth carrying into any valuation:

    1. **Most of the useful machining range is not yet reachable.** 120° is a loose roughing cap;
       real work wants 30–60°, and those are exactly the settings the generator currently ignores.
       The 2–5× path-length advantage in the Figure 6 reproduction is measured at caps the
       generator can actually meet, and will shrink toward 1 as tight caps force more passes.
    2. **Engagement is NOT monotone in either spacing or radius** (measured: 36.4° → 22.7° → 114.5°
       as radius grows). Any search that assumes monotonicity converges on a wrong answer *and
       still looks like it worked* — the most dangerous defect shape in this codebase, and the
       reason the fix uses a ladder search over an exact predicate rather than a bisection.

    This is the corpus earning its cost on day one: every headline number before it — 12× faster,
    shorter paths, fewer exceedances — looked good, and none of them revealed that the product's
    central input was doing nothing.

The rest of this section concerns the **continuous certifier**, which the generator no longer
calls and which remains disqualified on its own measurements.

The only measurement on branch B is 176 s (interrupted, two cases unfinished) → **39.42 s**
after tangent-root normalisation and exact support pruning — and the document immediately says
it is *"a repair, not an isolated benchmark or a Held–Pfeiffer comparison"* and that the matched
benchmark still has to measure it independently. Honest, and roughly a 4.5× repair.

But the **50× gate miss** that motivated this entire rewrite has not been shown to be closed,
and the new architecture is substantially heavier than the one it replaces. Branch A's own gate
analysis: smallest benchmark pocket > 5 min, square pocket > 10 min, both killed before
completing; the ≤ 10 s/pocket gate failed by ~50×.

**Direct measurement since (2026-08-20) settles it.** On a fixed probe segment against a stock
depleted one tool disk at a time:

| stock vertices | 4 | 8 | 12 | 16 | 20 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `audit_segment_tea_event_exact` | 0.12 s | 1.93 s | 5.27 s | 10.05 s | **14.70 s** |

Superlinear, on a *trivial* pocket — a real one carries thousands of boundary features and
hundreds of motions. Cost is not driven by bit length (`max_digits` saturates at 127 by the
second cut while time grows 8×), nor by distant geometry (eight cuts placed away from the probe
leave timings flat at ~0.13 s), nor by string handling (removing all 50 decimal round trips
changed nothing). It is ~90% CGAL bivariate curve-kernel work. See the measured breakdown under
[the continuous event partition](#stronger-in-rigour-catastrophically-weaker-in-cost-the-continuous-event-partition).

The load-bearing consequence: **the shippable certifier is the sampled + guarded one**, at 129 ms
per segment motion, one-to-two orders from Held rather than five. The continuous partition is a
research result, not a component, until it is reformulated.

### 7. CGAL licensing is a material, unpriced cost line

See `CGAL_LICENSE_GATE_ERROR` above. Confirm terms with counsel and budget the GeometryFactory
licence before any commercial-model assumptions are fixed.

### 8. One scaling landmine

`engaged_arcs_zone` `const_cast`s away constness on the stock arrangement
(`engagement_2.cpp:362–365`). Correctly documented, correctly encapsulated, sound
single-threaded. But parallel station queries are the obvious next performance lever, and this
is precisely what will block it. Worth knowing now, not at the parallelisation attempt.

---

## Not verified

- **Branch A test suites: executed and green.** `178 passed` repo-wide, including all three
  exact-kernel suites (`test_stock.py`, `test_engagement_oracle.py`, `test_engagement_audit.py`),
  run as `PYTHONPATH=src /Users/jelle/mambaforge/envs/cgal-dev/bin/python -m pytest tests -n auto`
  against the built `cp312-abi3-macosx_15_0_arm64` extension. An earlier note in this document
  said the suites could not be executed; that was an interpreter-resolution error on my side —
  the default `python` was a stale non-dev install lacking `_stock_2` — not a property of the
  code. The suites run and pass.
- **Branch B was read, not built or run.** All branch-B findings derive from source and
  documentation inspection via `git show`.
- **Licensing is flagged, not adjudicated.** The package classifications are the project's own;
  confirm with counsel.

---

## Proposed terms

There are two viable businesses here — a **generator** (track G) and a **verifier** (track V) —
sharing one asset and one consolidation tranche. Tranche 0 is deliberately structured so that
its output *selects the track* rather than merely passing a gate. Do not fund both.

### Tranche 0 — consolidation and track selection (shared, unconditional)

| Deliverable | Gate |
| --- | --- |
| CI runs `tests/adaptive` + `mypy --strict` + native gates on every push (drop `EXCLUDE_FROM_ALL`) | green on every push |
| Every file > 1,500 lines decomposed; task numbers stripped from permanent names | mechanical check |
| `docs/continuous_engagement.md` rewritten for a reader who is not its author | a second engineer can restate the architecture unaided |
| **Resolution-rate measurement** — `unresolved` frequency on a corpus of real third-party NC toolpaths | a number, with the neck-proximity breakdown |
| **Neck-degeneration ablation** — conservative-bound width vs. distance-to-neck vs. grid resolution | the divergence curve, or the claim is withdrawn |
| One measured wall-time on the 7-pocket suite | a number, not an estimate |

**No new capability lands until this passes.** It is a few weeks of work and it produces every
number the track decision needs.

### The decision rule

| Tranche-0 outcome | Track |
| --- | --- |
| Resolution rate high **and** wall times within ~1 order of Held | **G** — the generator is reachable; the verifier remains a later second product |
| Resolution rate high **and** wall times still far off | **V** — verification tolerates batch/overnight latency that interactive generation does not |
| Resolution rate low | **neither yet** — fix `unresolved` before spending on either front end; this is the existential number |

### Track G — the generator

| Tranche | Deliverable | Gate |
| --- | --- | --- |
| **G1 — one whole pocket** | A single non-trivial pocket generated end to end | Zero cap violations, zero `unresolved`, independently replayed from trace digests. Not a prefix. |
| **G2 — Held-matched performance** | 7-pocket benchmark suite, measured | Published as an A/B in the style of `docs/engagement_zone_query.md` — still the model for how this team reports numbers when it is trying |

### Track V — the verifier

| Tranche | Deliverable | Gate |
| --- | --- | --- |
| **V1 — front end + pilot** | APT-CL / CAM-API reader producing exact motion primitives from decimal input; certification of one real customer program | A verdict on third-party output, with an exact witness for each violation and an arc-consistency report |
| **V2 — conservative screen + scale** | Sound one-directional pre-filter, exact certifier on the flagged set | A 10⁵-block program certified end to end; screen proven conservative (no false clears) by construction, not by sampling |

Track V's V1 is materially cheaper than G1 — it deletes path *generation* entirely — which is the
substance of the claim that the verifier reaches revenue sooner.

**Condition precedent to any commercial discussion:** resolve `CGAL_LICENSE_GATE_ERROR` — either
a GPL-compatible distribution decision or a GeometryFactory entitlement covering CGAL 6.0.1 and
every instantiated package.

!!! note "The two tracks carry different licence exposure — ask counsel early"

    GPL obligations attach to *conveying* software. A shipped generator library or plug-in
    conveys; a **hosted verification service** arguably does not, which would place track V's
    delivery model outside the trigger that blocks track G. CGAL is GPL, not AGPL, so the
    network-use clause that would close this route does not apply — but GeometryFactory's
    commercial terms may address hosted use directly. **This is a strategically material
    question, not a footnote: it may decide the delivery model.** Flagged as a derivation
    requiring counsel; not asserted as settled.

**Kill criterion:** a high `unresolved` rate on real third-party geometry at Tranche 0. That
number is existential for *both* tracks — a certifier that abstains on real industrial input
proves nothing a customer can use, and no amount of front-end or generator work repairs it. Missed
wall-time targets are not a kill; they select track V over track G.

**Covenants:**

1. No epsilon enters a decision path to buy performance. The moment that trade is made, the
   differentiator is gone and this is another CAM library.
2. A second human reader on every merge to the integration branch. Not for correctness — the
   gates can do correctness — but because an artifact only one person has ever read is not an
   asset, it is a liability with good test coverage.

---

## Alternative thesis: the certifier as the product

The generator thesis requires every tranche on track G to land. The verification thesis requires
substantially less, reaches revenue sooner, and — on the evidence in this report — is the
stronger business. This section is the case for it. It is **track V** in
[Proposed terms](#proposed-terms), selected on evidence at Tranche 0, not a consolation prize
and not a parallel spend.

### The certifier is already provenance-blind

Its interface asks nothing about where a motion came from:

```
(exact stock state, motion primitive, tool radius, cap) → certified | violated | unresolved
```

Feed it motions from Mastercam, hyperMILL, PowerMill, NX, ESPRIT or a hand-written program and it
answers the same question with the same proof. The generator is one consumer of the certifier;
a verification product is another, and it is strictly less work because the hard half — deciding
*what* the path should be — is deleted.

### Where the incumbent sits, and why not to fight it

VERICUT (CGTech) verifies **NC programs**: it emulates the control, simulates material removal
against machine, fixture and holder, and reports collisions, gouges and overtravel. Its **Force**
module goes further than most summaries credit — it derives chip thickness and cutting force from
simulated engagement and rewrites feedrates. Any pitch premised on *"VERICUT doesn't do
engagement"* will be corrected in the first meeting.

The defensible distinction is the **type of the verdict**:

| | VERICUT / Force | This certifier |
| --- | --- | --- |
| Method | volumetric simulation (dexel/voxel), resolution knob | exact algebraic event partition |
| Verdict | "no violation detected at this resolution" | "≤ θ_max at every point of the motion, proved" |
| Missed-case model | silent — undetectable below grid size | explicit `unresolved`, never silent |
| Artifact | a report | content-addressed certificate, independently re-checkable |
| Can it say "I don't know"? | no | **yes — and that is the feature** |

A sampled simulator always produces an answer, and that answer's reliability is a function of a
resolution setting *the user chose*. That is the seam. The play is not removal simulation —
CGTech has three decades of controller emulation and machine kinematics and it will not be
caught. The play is **the proof layer that a dexel architecture structurally cannot produce.**

Positioning: *VERICUT proves the program will not crash; this proves it will not exceed the
engagement the tooling data assumes.* Complementary, sells alongside, triggers no bake-off.

### The counter-argument to pre-empt

The claim above — that a sampled architecture cannot produce a proof — is the one a competent
opponent will attack first, and the naive form of it is **wrong**. A dexel representation can be
made *conservative*: store inner and outer material bounds per ray, and every engagement query
returns an interval `[lo, hi]`. If `hi ≤ θ_max`, that is a genuine one-directional proof of
non-violation from sampled data. Anyone claiming "sampling can never prove anything" will lose
that exchange.

The real argument is narrower and much stronger, and it is the defect **this team already found
in its own certifier**: the run-merge discontinuity.

Two engaged runs separated by a thin void gap on the rim contribute either two moderate runs or
one large run, and `max_run` jumps by `O(1)` — tens of degrees — the instant the gap closes. A
conservative sampled method, unable to decide whether a sub-grid gap is open, **must assume it is
closed**. That is gap-closure pessimism, arrived at by necessity rather than by choice. And at a
neck the gap width passes continuously through zero, so for *every* finite grid there is a
neighbourhood in which the conservative bound saturates at the vacuous 360°.

Refinement does not rescue it: the sub-resolution neighbourhood shrinks but never empties, while
cost grows as `h⁻²`–`h⁻³`. So the sampled method must abstain precisely where engagement
violations actually occur — at necks — while the exact method **decides the gap** and certifies
there. Ordinary interval error elsewhere is benign by comparison: an `h`-sized position
uncertainty gives roughly `h/r` angular uncertainty away from tangency (arcminutes at realistic
grids), degrading to `√(2h/r)` near grazing incidence.

| Regime | Conservative dexel bound | Exact certifier |
| --- | --- | --- |
| Open field, transverse crossing | tight (`~h/r`) | exact |
| Near-tangent crossing | loose (`~√(2h/r)`) | exact |
| **Sub-grid rim gap near a neck** | **vacuous (360°) at every finite `h`** | **decided exactly** |

!!! warning "Measure this before asserting it"

    The neck-degeneration claim is a derivation, not yet a measurement. It needs the ablation:
    conservative-bound width against distance-to-neck, swept over grid resolution, showing the
    divergence and its resolution-independence. Cheap to produce, paper-grade once produced, and
    it converts the strongest slide in the deck from an argument into evidence. Until it exists,
    state the claim as a derivation.

### Go in before the post, not after

| Front end | Assessment |
| --- | --- |
| **(a) Post-processed G-code** | Universal and vendor-independent, but inherits Fanuc/Heidenhain/Siemens dialects, macro-B, canned cycles, cutter comp, work offsets, kinematics. That *is* the VERICUT moat; rebuilding it is a multi-year tax to arrive where CGTech already is. |
| **(b) CL data / APT-CL / CAM API** | Clean motion primitives, no modal-state emulation, ~5 CAM systems cover most of the market and all emit APT-CL. Verification lands *before* the post — where the programmer can still act on the result cheaply. **Recommended.** |

Option (b) also occupies a different point in the workflow than VERICUT (which is post-post, at
machine-code level), so the two do not collide.

### Two features fall out of exactness for free

**Decimal input is exactly rational.** G-code and CL coordinates are decimal *strings*:
`X12.3456` is exactly `123456/10000`. Parsed straight into `Epeck::FT`, the certifier operates on
precisely what the machine was told — **zero injection error, not even the sub-ulp gap the
current `double` boundary carries.** `docs/exactness.md` already separates "the binary double
0.1" from "the mathematical 1/10"; for NC input the decimal reading is the correct semantics, so
this front end is *more* exact than the boundary the kernel has today. Verifying foreign
toolpaths is technically cleaner input than the generator's own.

**Inconsistent arcs become detectable.** `G02/G03` with I/J supplies endpoint *and* centre. Real
post-processors routinely emit arcs whose endpoint is not exactly on the specified circle, and
controllers silently absorb the discrepancy. An exact reader reports it —
*"block N4710 emits an arc whose endpoint lies 0.4 µm off the specified centre."* Nobody
currently surfaces this. It also forces a declared policy for the case, which is exactly the
`Manufacturing_policy` object already sketched in `docs/exactness.md`.

### Scope discipline

In 2.5D at constant Z a flat endmill **is** exactly a disk; so are ball-nose and bull-nose tools
at a given Z, with an effective radius. The existing exact machinery therefore already covers a
large slice of real roughing. 3-axis surfacing and 5-axis are **not** addressable by a planar
exact model and must not be claimed.

The wedge is 2.5D roughing and pocketing — which is where engagement matters commercially
anyway: roughing is where tool breakage, chatter, spindle load and cycle time live. Finishing is
light-engagement by construction.

### The verdict menu is wider than TEA

The same exact stock model answers several questions CAM vendors currently only approximate; most
are pure boolean queries on machinery that already exists:

| Verdict | Basis | Commercial hook |
| --- | --- | --- |
| Slotting detection | the cap predicate at θ = π | tool-breakage class |
| **Air cutting** | cutting feedrate with TEA = 0 | **direct cycle-time saving; self-funding ROI** |
| Rapid into material | any `G00` intersecting stock | real crash class, exactly decidable |
| Residual / uncut material | exact boolean against nominal part | scrap and rework |
| Effective chip thickness | derived from exact engagement | the quantity Force *models* |
| Stock continuity across tool changes | exact hand-off between operations | multi-op process qualification |

### The buyer nobody is serving

Beyond the obvious two — aerospace/medical tier-1s with process-qualification requirements, and
CAM vendors licensing it as an OEM component (they all need engagement control; none can prove
it) — there is a third that is underrated:

**Cutting tool manufacturers.** Sandvik, Kennametal, Iscar and Seco all publish feeds-and-speeds
keyed to a radial engagement `a_e`. Their entire cutting-data recommendation, and any tool-life
claim resting on it, is void if actual engagement exceeds the assumed value. **Nobody verifies
that the assumption holds.** A certificate that a program respects the engagement its tooling
data assumes converts published cutting data from a recommendation into a warrantable spec —
a genuinely new artifact in that supply chain.

### The two facts that decide whether this is a business

**1. Resolution rate — measurable now, before any front end exists.** If `unresolved` returns on
a large fraction of real industrial toolpaths, the product is unsellable regardless of rigour.
This is *the* metric, and the existing certifier can measure it today against a corpus of real NC
programs. **Cheapest and highest-information experiment available; run it before anything else.**

**2. Performance gets worse here, not better.** Verifying an arbitrary program means simulating
material state from block 1, so verification cost *is* the depletion cost — the bottleneck
already identified in risk #6. Real NC programs carry 10⁵–10⁶ motion blocks, not 271. **The
arc-contour representation is not optional for this business either.**

There is a sound architecture for it, in the split this team already thinks in: a **cheap
conservative screen** — sampled, but over-approximating so it can only raise false alarms, never
issue a false clear — followed by the **exact certifier on the flagged set**. Sound end to end,
fast in the common case. Soundness lives in the screen's conservatism; speed lives in how rarely
it fires.

### Honest limits

2.5D only. No machine kinematics, no holder or fixture collisions, no controller emulation —
those remain VERICUT's. And the certificate is a statement about the *geometric model*: it says
nothing about deflection, runout, or a stock casting that is not the shape it was told. Claim
exactly what is proved and nothing further — a discipline this team already keeps better than
most.

---

## The one-line version

Six days bought a genuinely stronger proof system — complete exact event partition and a true
medial axis, both real advances on the published state of the art — and simultaneously bought
114,000 unreviewed, unenforced, increasingly unreadable lines. Measurement has since added the
decisive qualifier: **the stronger proof system is ~116× slower than the project's own existing
certifier and ~10⁵ off the published state of the art, and that is a formulation problem, not a
tuning problem.**

Since then the generator side has been fixed and **now runs inside Held's published range** — 35×
faster with byte-identical output, via three representational corrections that needed no kernel
change and no architectural decision. That materially changes the shape of the risk: what remains
slow and untrustworthy is the *continuous certifier*, a component the generator no longer calls.

**Buy the mathematics; refuse to pay for the volume; ship the generator that reached Held's range
and the sampled + guarded certifier beside it, not the continuous one that is five orders off and
incomplete on a third of hard motions; and make the next cheque contingent on CI turning green and
on one number — how often this thing says `unresolved` about someone else's toolpath.**

Everything else is engineering. That number is the business.
