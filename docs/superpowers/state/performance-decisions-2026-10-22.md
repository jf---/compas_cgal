# Performance decisions (2026-10-22)

Three decisions, in the order they should be taken. `W` = `.superpowers/sdd/2026-10-13-product-pipeline/`.
Every timing below was taken on a loaded machine (load 11–555, quiet limit 2): **all are indicative**;
ratios are same-process A/B, never cross-session.

1. **(a) GMP: run the A/B first, switch only on identical digests and a measured quiet gain.** CB found
   gcd canonicalisation is 52 % of certificate CPU, so this is now the largest single lever on both
   budgets, and it fixes the per-op cost that (c) inherits.
2. **(b) Two budgets, both gated:** generation (every decision that shapes the path, incl. connector
   gate + repair) stays at 1.7 / 2.2 / 9.0 s; certification (coverage certificate + gouge audit, once
   per final path) gets its own budget.
3. **(c) Native end-to-end Ω\* initial stage** (placement + ownership + start order + span selection →
   compact arrays), not more per-function batches. Ω\* also needs `connector_repair` cut: native
   `initial` alone cannot green Figure 5 or skis.

## (a) GMP

**Question.** Do we give up upstream's GMP-free portable wheels for a faster exact number type?

**Recommendation.** Do the A/B first, decide on evidence. Build the native modules in a separate build dir
with GMP/MPFR (C++ slot, quiet window, one pyproject edit for the pixi deps). Measure the production
certificate (workload B below) and Ω\* generation. Switch only if every exact module's outputs are
digest-identical (Ω\*, polygon, certificate LAZY ≡ EXACT_ONLY, last-16) **and** the quiet min-of-5 gain is
material. If adopted, GMP is a hard dependency (no `HAS_GMP`, no dual path) and the divergence from
upstream wheels is stated in the docs. If rejected, record the negative result and take CB's no-GMP
lever instead (gcd-free signs of exact integer expressions in the hot predicates, stock boost type,
per-site A/B). No decision path changes either way: exactness belongs to the number type.

**Evidence.**
- Disabled globally: `CMakeLists.txt:104-121` (`-DCGAL_DISABLE_GMP -DCGAL_USE_BOOST_MP`, find-package
  off), inherited from upstream `da292f3b` for cibuildwheel portability; `Epeck::FT::ET` is boost
  `cpp_rational`, gcd on every op; GMP/MPFR absent from pixi.lock (`W/rv15-review.md` §GMP).
- L1 native placement: 33 µs/placement on Monstera; `sample` of the C++ loop: 31 % `libsystem_malloc`,
  8 % memmove, 3 % `__udivmodti4`, no algorithmic hot spot: bignum temporary churn
  (`W/op-omega-profile-report.md`).
- Certificate (CB, workload B): boost bignum is 62.5 % of native CPU, of which gcd + exact divide
  (canonicalisation) is **52 %** and multiply only 12 %. Projected GMP gain −31…−42 % (B ≈ 3.1–3.7 s):
  **extrapolated from the profile, not trial-built** (`W/cb-certb-budget-report.md`).
- GMP `mpq` also canonicalises per op (smaller constant), so L1's unreduced-integer scheme still pays
  off with GMP (`W/rv15-review.md`; 36 vs 70 µs/arc, `W/cp3-cpp-report.md`).
- Friction: `src/exact_algebraic_1.cpp:17-33` `#error`s unless GMP is disabled; CP4N used
  `convert_to<double>` for lack of MPFR.

**Cost.** One C++-slot session: pixi deps + per-target defines + `exact_algebraic_1` guard, a full rebuild
(kills every importer: quiet window only), the digest sweep and a quiet A/B. If adopted: wheels need
native GMP on every platform; upstream PRs carry a build-config diff.

**Unblocks.** The true per-op floor of every exact native module, so (c)'s budget is set on the right
number type before code.

## (b) Two budgets

**Question.** Is the coverage certificate and gouge audit part of the generation budget?

**Recommendation.** No: two budgets, both gated by the same rules (process CPU, min of 5 fresh runs,
`MachineNotQuietError` above load 2).
- **Generation** = everything that decides the path, including the exact connector gate + repair and
  NOTHING (user ruling D1). Unchanged 1.7 / 2.2 / 9.0 s, applied to Ω\* too (OB's red test already
  does this).
- **Certification** = coverage certificate + gouge audit, run **once per final path**, never inside a
  generation loop. Budget: certificate ≤ 5 s CPU on Monstera production path B (existing ruling, kept;
  see CB), gouge audit ≤ 2 s Monstera (native clearance target).

Splitting is sound because certification does not change the path; anything that does (gate, repair)
stays in generation so it cannot hide in the other budget.

**Evidence (indicative).**

| | fig5 | skis | Monstera | source |
|---|---:|---:|---:|---|
| polygon generation, min CPU (budget 1.7 / 2.2 / 9.0) | 1.72 | 2.00 | 7.63 | `generation-budget-gate-analysis.md` §8, load 11–13 |
| … with F5's 9 no-ruling levers | 1.447 | ≈ −0.12 | ≈ −0.67 | `W/f5-fig5-levers-report.md`, load 14 (skis/Monstera load 35–180) |
| Ω\* generation, min CPU | 4.23 | 5.57 | 16.22 | `W/ob-omega-budget-report.md`, load 35–184: 2.5× / 2.5× / 1.8× over |
| … of which `initial` / `connector_repair` | 2.30 / 1.05 | 3.03 / 1.71 | 8.85 / 3.59 | same |
| Ω\* after L1 + OL landed (gen CPU min) | 4.17 | 5.56 | 14.78 | `W/lw-l1-wiring-report.md`, load 10–396 |
| coverage certificate T(Ω\*), production path B, min CPU | | | 5.375 (budget 5) | `W/cb-certb-budget-report.md`, load 450–497 |
| gouge audit Ω\*, native / Python | 0.62 / 5.56 | 0.39 / 5.64 | 1.69 / 11.56 | `benchmarks/results/held_native_clearance.md`, load 25–32 |

**Cost.** One new acceptance test + pixi task for the certification budget; the Ω\* generation clock
already exists (`tools/held_generation_clock.py --domain omega`).

**Unblocks.** An honest red/green per budget: today the polygon path is near green, Ω\* is red on
generation by 1.8–2.5×, and the certificate is borderline (quiet verdict pending).

### Certificate budget (CB)

Source: `W/cb-certb-budget-report.md` (load 450–555, not a verdict). Two workloads were conflated:

| workload | result | min CPU | role |
|---|---|---:|---|
| **B** production Ω\* path, exact ∂Ω\* connectors | COVERED, 0 witnesses | **5.375 s** (all runs so far 5.07–5.73) | **the gate**: borderline, quiet run pending, likely 5.0–5.4 s |
| **A** `coverage_certificate_omega_star` (polygonal path vs T(Ω\*)) | 4078 witnesses | 6.65 s | diagnostic of non-coverage; 42 % is witness work; not the gate |

- **Time (B):** setup + crossings 29 % (exact ∂Ω\* connectors 9.6 %, crossing table 13.4 %), circle sweep
  31 %, exposed-edge processing 30 % (`edge_parameter` 18 %: ~15.7 k edges pay for an exact parameter,
  only 1473 reach a target decision), line sweep 6 %, Python < 1 %.
- **Levers, by gain on B:** (1) number type: GMP −31…−42 % (decision (a)), or gcd-free integer signs
  without GMP; both exactness-neutral. (2) Interval pre-dismissal of exposed edges, ≈ −15…−20 %: an
  outward-rounded enclosure of the whole sub-arc is a provable superset, exact path stays the decision;
  measure the hit rate first. (3) Share exact crossings across coincident support circles, ≈ −10 %:
  identity of exact objects, neutral; weakest estimate.
- **Gate defect:** `tools/coverage_certificate_omega_star.py` labels OVER/within on **median wall**, and
  the native phase timers are `steady_clock` wall. Fix: gate on process CPU (thread CPU in the native
  timers), min of 5, quiet-only, on workload B.

### What needs per-run validation on Ω\* (added 2026-10-23)

A property needs a per-run check only if generation does not decide it exactly.

| property | decided in generation? | per-run check |
|---|---|---|
| no gouge (tool centre in Ω\*) | yes: exact containment per circle, exact connector gate + repair | none; the gouge audit becomes a test oracle on the fixtures, out of the production path (frees the ≤ 2 s Monstera audit budget) |
| coverage of T(Ω\*) | no: containment says nothing about the union of swept annuli | coverage certificate (≤ 5 s), kept |
| engagement on the actual stock | partly: refine / NOTHING decide exactly on swept annuli, not on the actual stock | last-16 stock certificate (160 s circle CPU on Monstera, load ~600), kept until refinement decides on the actual stock |

Correct by construction still needs differential tests against an independent oracle (polygon path,
Python reference, full replay); those are test-suite cost, not per-path cost. Open design question
for the user: move the actual-stock engagement decision into refinement, which would retire the
most expensive certificate.

## (c) Next C++ direction

**Question.** Port the Ω\* initial stage natively end to end, or keep batching single functions?

**Recommendation.** End to end. One native call from declared contacts to compact arrays (owned indices,
binary64 centre/radius, provenance ints, span selection); Python builds objects only for the selected
circles. The Python stage stays as the differential reference (by-repr equality on fig5, skis, Monstera +
Hypothesis), and N4 span selection gets its exactness design note first. Set the stage budget before code.
Second target, same slot: `connector_repair` (22–26 % of Ω\* generation).

**Evidence.**
- `initial` is 54 % / 63 % / 54 % of Ω\* generation; `connector_repair` 26 / 25 / 22 % (`W/op-omega-profile-report.md`).
- L1 (placement batch) was byte-identical on 57,743 placements but cut generation only −3.5 / −19 / −4.7 %
  (`W/cp3-cpp-report.md`): the Monstera batch is 2.30 s, of which the C++ call is 1.18 s and the Python
  decode ≈ 1.08 s (`W/lw-l1-wiring-report.md`). Every further per-function batch pays that decode again.
- Laziness cannot fix it in Python: `_start_order` sorts every circle and the lanes read every
  `paper_candidate` (`W/lw-l1-wiring-report.md`), so the decode only shrinks if ordering and selection move
  native too. Monstera places 37,573 candidates to emit 3,983 circles.
- Native `initial` alone is not enough: Ω\* non-initial stages already cost 1.93 / 2.54 / 7.37 s against
  1.7 / 2.2 / 9.0 (OB stage table).

**Cost.** The largest item on the list: native ownership argmin, start order, span selection (N2–N4 folded
in), a compact-array binding, differential tests; one C++ slot for several sessions. Best done after (a),
so it is built and budgeted on the final number type.

**Unblocks.** A path to Ω\* under the generation budget; retires the queued polygon N1–N4 per-function
levers, which the polygon generator no longer needs (it is within budget on Python levers).

## Parked F5 Figure 5 levers (`W/f5-fig5-levers-report.md`)

No ruling needed, byte-identical, Figure 5 1.678 → 1.447 s min CPU: U1b (`Point2.build` single
validation), L2+L3 (`_Surd` NamedTuple), G1 (exit `gc.collect(0)`), L4 (constant-origin common scale),
L10 (integer `_project_to_offset_side`), L1 (drop `Point2[WorldXY]` alias), L6 (`_transition` int
division), L5 (`_hypothesis_geometry` locals), L9 (`_boundary_arclength`).
**Recommendation:** land as one unit, re-applied on PF's committed files, with a quiet gate run.
U1b goes to the `units.py` owner for review (shared core type). G1: take it; same garbage collected,
the exit collect stays counted.

Need a ruling:
- **P3** (no connectors on ordering-only call 0; −0.05 / −0.05 / −0.21 s): approve. The check is
  unreachable by an exact separation bound (≥ 2⁻⁴²⁴) and still runs on call 1.
- **L8** (double pre-key, exact Fraction tiebreak; −0.03 / −0.02 / −0.08 s): approve. Correct rounding
  is monotone, so the double is a filter with exact fallback, not a decision; order proven identical.
- **Drop `PaperCircleCandidate` `TOL.is_between`** (≈ −0.015 s): replace with the exact check, don't just drop.

## Other decisions (§3b)

- **Dead code** (`flank_lines`, `_origin_key`, `_dyadic_*`): remove after the four-consumer check;
  `_squared_distance` is the P2 test's Fraction oracle: move it into that test, don't delete it.
- **MM fix 4** (keep only owned hypotheses, −80 MB; 38,886 of 120,271 kept): approve, initial digest
  change stated in the commit (`W/mm-omega-memory-report.md`).
- **Junction lines in the Ω\* emitter (OI2):** emit them, so the certified path is the emitted path and
  RV3 Important 2 closes rather than stays narrowed (`W/oi2-oc-fixes-report.md`).
- **Python clearance oracle:** keep it, test-only, as the differential reference (identical decisions on
  10 path/workload pairs; 3–22× slower than native); native stays the production default.
- **P3:** approve (above).
- **`tests/benchmarks/held_figure5_path_pre_lever_oracle.py:328`** `float(OverlapMargin) > 0.0`: approve
  the edit to the exact `OverlapMargin > 0` comparison; a float decision in a reference oracle can diverge
  from production on the disagreement inputs (`W/rv2-omega-units-review.md`).
