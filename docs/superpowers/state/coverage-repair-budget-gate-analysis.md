# Coverage repair: `repair_s` budget gate analysis (2026-10-22)

**The gate is back inside its budget without moving it: Monstera `repair_s`
is 0.147 s against the 0.5 s budget, down from 8.39 s.** The fix is the
exact early stop of the splice continuation. It changes no output: the same
2 sources are inserted, 4438 → 4440 circles, and the result is byte-equal to
a full re-refinement. No margin, epsilon or budget change is involved. The
Monstera acceptance pins are stale for an unrelated reason (CW's connector
repair, +72 circles) and need a user go to re-pin. The new values are
derived [below](#stale-monstera-acceptance-pins).

Order of events: this page was written from the measured breakdown (sections
1 to 4) before the lever was built. Section 5 adds the after-measurements.

## 1. Measured breakdown before the fix (Monstera)

Fixture: the slow test's `refine_boundary_ordered_path_certified` inputs,
cached once (`build/es/repair_probe.py`, log
`build/es/probe_before_figure8_monstera.log`), load average 10.7 to 11.2.
One repair round, 8 witnesses, 16 candidates. Two greedy steps try 16 + 12 =
**28 splices**.

| quantity | value |
| --- | --- |
| `repair_round` wall | 7.85 s (RV3 measured `repair_s` 8.39 s at load ~12) |
| sum over the 28 splices | 7.79 s (99 % of the round) |
| one splice | 0.055–0.47 s, linear in the suffix length (about 0.12 ms per emitted circle) |
| suffix re-refined per splice | 327–3862 emitted circles (slot 334 → 3862, slot 581 → 3423, slot 1663 → 1678, slot 2489 → 351) |
| one full refinement of the emission | 0.52 s (4438 circles) |

Profile of the round (`build/es/probe_profile_monstera.log`, cProfile,
16.2 s instrumented):

| function | cumulative | share |
| --- | ---: | ---: |
| `continue_figure5_engagement` (28 calls) | 16.08 s | 99 % |
| `_Emission.advance` (37 417 targets) | 14.09 s | 87 % |
| `midpoint` (19 000) | 4.93 s | 30 % |
| `_ConnectorGate.certifies` (61 507) | 3.90 s | 24 % |
| `exact_predecessor_engagement` (79 817) | 3.01 s | 19 % |
| `circle_on_side` of every successor, computed up front | 0.93 s | 6 % |

Per splice, almost all of this is re-deriving circles that come out
identical to the old emission: every splice re-emits the whole suffix.

## 2. Why the splice must be suffix-equivalent

The repair's claim is that a splice equals `refine_figure5_engagement` on
the augmented selection. A two-interval splice (re-refine only
`(s⁻, q)` and `(q, s⁺)`) is not equal. RS's Hypothesis counterexample
(`test_splice_sees_the_circles_emitted_before_its_interval`, side 0 of the
10 × 10 square, cap 80°):

- contacts x = 2.00 (R 1), 2.01 (R 2), 2.02 (R 1.5), 2.71 (R 1), 2.84 (R 1);
  candidate 2.71 spliced into (2.02, 2.84);
- the full run accepts the connector 2.02 → 2.71 by tier 1b, from the R = 2
  circle at 2.01;
- a local three-source window never sees that circle and inserts the
  midpoint (2.365, R 1.25): 11 circles against 10.

The mirror case also holds: connectors after `s⁺` gain `q`'s circles as tier-1b
priors and may need fewer midpoints, so keeping the old tail verbatim is not
exact either. Hence RS's continue-from-slot: keep the emission up to `s⁻`
(decided by the earlier sources alone) and continue to the end. It is exact,
but it costs O(suffix) per candidate.

## 3. Levers

### 3a. Exact early stop (built)

**State.** The refinement is sequential. At a source boundary, its whole
state is the emitted list. Every decision reads its predecessor and target,
and beyond that only these:

| decision | reads besides its own pair |
| --- | --- |
| cap (`exact_predecessor_engagement`), midpoint, corner bridges | nothing |
| connector tier 1 | nothing (own circle only) |
| connector tier 1b (only when tier 1 fails) | the circles the superset filter `_tier_1b_neighbours` keeps for the connector's pieces |
| sharp-vertex refusal, floor hypothesis | the emitted circles with \|c − V\|² ≤ (ρ + 2r)², exact |
| circle budget | the emitted count |
| seam | `circles[0]` plus the above |

Each refinement now records a **footprint** for every decision that reads
other circles (`EngagementTrace`). For tier 1b, the footprint is the
connector's outward piece bounds, exactly the binary64 values the filter
uses. For the floor hypothesis, it is the vertex.

**Stop condition.** The continuation emits `q` and then each later source.
After later source k, suppose its interval re-emits exactly the old circles
for that source. Let X be the circles emitted since the slot and Y the old
circles they replace, with the common re-emitted tail removed from both.
Then the continuation adopts the old emission's remaining circles iff both
of these hold:

1. **No influence.** Look at every recorded decision of the old run at an
   emission position at or after k's boundary. No circle of X ∪ Y is kept by
   `_kept` for any piece of a tier-1b footprint. No circle of X ∪ Y has
   \|c − V\|² ≤ (ρ + 2r)² in `Fraction` for a floor footprint at vertex V.
2. **Budget.** `len(old) + (len(new prefix) − len(old prefix)) ≤ max_circles`.

**Why it is sound.** Use induction over the old run's remaining decisions.
Both runs sit at the same source boundary. Their emitted lists are A + X + C
and A + Y + C, and they continue with the same sources.

- Each decision's predecessor and target are equal, because the tails grow
  identically.
- A tier-1b decision builds its grid, candidate index and bisection from the
  filtered circles in emission order. By (1), no circle of X ∪ Y passes the
  filter in either run, so both runs feed the decision the same circle
  sequence (A ∩ kept, C ∩ kept, tail ∩ kept). Every exact test therefore sees
  identical inputs.
- The floor hypothesis skips every circle of X ∪ Y by the same exact test.
- Decisions without a footprint read nothing else.
- By (2), no budget check flips. Every old count before an append is at most
  `len(old) − 1`, so the new count stays below `max_circles`.

The new run therefore emits exactly the old tail. No float decides anything
here: `_kept` is the gate's own outward filter, so "not kept" means "not on
the decision's grid". The floor test is rational.

**What "changed circle" must cover.** It covers both X and Y. A removed
midpoint (Y) could be what certified a later connector in the old run. A
new circle (X, `q` included) could certify a connector the old run refused,
as in RS's counterexample. A sharp-vertex refusal or the seam is covered by
the same footprints: the seam's decisions are recorded like any other, and
`circles[0]` is in the common prefix.

**Estimated saving (before building).** Tier 1b is rare. The profile shows
about 57 tier-1b piece evaluations per splice, against about 2200 connector
decisions. So a stop right after the slot interval is the common case. That
gives an estimated splice cost of a few intervals (about 1 ms) plus the
trace bookkeeping, so the round would cost about 28 × (1 to 10) ms, in
0.03–0.3 s.

### 3b. Candidate pruning (not built, needs a ledger ruling)

The greedy key is (−covered, added, index), and `covered` is known before
splicing. A candidate covering fewer open witnesses than an admissible
(non-refused) candidate can never win, so it need not be spliced. This is
exact. On Monstera:

- step 1: 1274 (2), 1276 (3) and 1277 (1) are dominated by the 4-witness
  candidates;
- step 2: the same three are dominated.

That prunes 6 of 28 splices, about 2.6 s of the 7.8 s before the early
stop and about 0.01 s after it.

The catch: the ledger's `refusals` would no longer list dominated
candidates that would have been refused. That changes `RepairRefusal`
semantics, so it needs a ruling. After 3a it is not worth that change.

### 3c. Smaller items, measured after 3a

- `_ConnectorGate.__init__` decided `exterior_angle_exceeds_cap` for every
  boundary vertex (386 on Monstera) on each construction: 0.56 s of a 0.90 s
  profiled round, once per splice. It is now decided per vertex on first
  query, which is exact and identical. A short continuation crosses few
  vertices.
- `circle_on_side` of every successor was computed up front (6 %). It is
  now computed per emitted source, so an early-stopped continuation never
  computes the suffix.

## 4. Reachable budget

From the profile, a splice that stops after one or two intervals costs its
fixed overhead plus a few intervals. The remaining O(n) terms per splice are
these:

- the source-identity check over the later sources (tuple comparisons, about
  1 ms);
- copying the prefix list;
- filtering the trace.

The estimate is 28 × 3–5 ms ≈ 0.1–0.15 s, inside 0.5 s. The witness matching
(`emitted_form` of every dropped source, about 4000) adds a few tens of
milliseconds.

## 5. After (measured)

Same process, same inputs, load average 12–18 (`build/es/ab_probe.py`, logs
`build/es/ab_*.log`). A is the old continue-from-slot
(`continue_figure5_engagement` from the slot to the end). B is
`continue_with_early_stop`. Each candidate's B result is asserted equal to A,
trace positions included.

| Monstera candidate | slot | later sources | re-refined by B | A s | B s |
| --- | ---: | ---: | ---: | ---: | ---: |
| 744 / 746 | 334 | 2361 | 2 | 0.518 / 0.473 | 0.0048 / 0.0044 |
| 1273–1277 | 581 | 2114 | 1 | 0.43–0.46 | 0.0039–0.0041 |
| 4035–4038, 4046 | 1663–1666 | 1029–1032 | 1 | 0.20–0.21 | 0.0033–0.0044 |
| 6127, 6128, 6159, 6161 | 2489–2501 | 194–206 | 1 | 0.047–0.084 | 0.0029–0.0038 |
| **sum, 16 candidates** | | | | **4.46 s** | **0.060 s** (74×) |

Round and gate:

| quantity | before | after |
| --- | ---: | ---: |
| `repair_round` (cached inputs) | 7.85 s | 0.154 s |
| `repair_s`, end-to-end slow-test fixture (`build/rv3/monstera_repair_probe.py`) | 8.392 s (RV3, load ~12) | **0.147 s** (load 11.8–12.1, `build/es/monstera_repair_probe_after.log`) |
| outcome | (744, 1273), 4438 → 4440, 4440 cuts / 4160 connectors, 0 over-cap transitions | identical |
| repaired emission vs full `refine_figure5_engagement` of the augmented selection | | equal (asserted in the probe) |
| full refinement, trace recording overhead (same process, min of 3) | 0.535 s | 0.539 s (×1.007); skis ×1.03, Figure 5 ×0.99 |

Figure 5 and crossed skis: the first certificate holds, so there is no repair
round (`repair_s` = 0 before and after). Their refinements are unchanged
(path digests identical to the CW/RS post-repair values).

## Stale Monstera acceptance pins

**Resolved 2026-10-22 (user go):** the test is re-pinned to 4438 → 4440 and 4160, and all 4 tests pass (`repair_s` 0.144 s). The table below records the derivation.

`benchmarks/acceptance/test_held_coverage_repair_monstera_slow.py` (dated 21
Oct 22:49, before CW's connector repair) pins values from before that repair.
Do not edit the test without a user go. The table lists the new values, each
measured (`build/es/pins_probe.py`, `build/es/monstera_repair_probe_after.log`,
`tools/coverage_certificate_acceptance.py --json/--md build/es/…`).

| pin | old | new | derivation |
| --- | ---: | ---: | --- |
| `EMITTED_BEFORE` | 4366 | **4438** | 2695 selected sources + 1743 refinement circles (midpoints, corner bridges, seam). CW's connector repair accounts for +72 of them (4366 → 4438, `connector-repair-wiring-report.md`). Equal to CW's post-repair path digest (4438, `04403396…`); the committed digest record is still the stale 4366 |
| `EMITTED_AFTER` | 4368 | **4440** | 4438 + 2 inserted sources (744, 1273). Each splice adds exactly 1 circle (`added` = 1 for both), with no midpoint |
| `CONNECTORS_AFTER` | 4088 | **4160** | 4440 transitions (4439 consecutive + the seam, `Figure5SeamClosure.CONNECTOR`) − 280 fan transitions (same boundary site, stationary, no region) = 4160. Before repair: 4438 − 280 = 4158 |
| `INSERTED_SOURCES` | (744, 1273) | (744, 1273) | unchanged: integer tie-break, both 4-witness candidates of their islands, lowest index wins on equal `added` |
| rounds / refusals / certificate runs | 1 / () / 2 | 1 / () / 2 | unchanged |
| `REPAIR_BUDGET_S` | 0.5 | 0.5 | unchanged; measured 0.147 s |
| module docstring: "repair about 0.12 s" | | 0.15 s | measured |

Measured with the early stop (`build/es/monstera_slow.log`, load ~28): 2 of
the 4 tests fail, both only on the stale counts ((4440, 4160) against
(4368, 4088), and (4438, 4440) against (4366, 4368)). The cap test and the
time test pass. RV3 had 3 red, the time test among them.

Other places that cite the old numbers (owners other than ES):

- `benchmarks/results/coverage_certificate_acceptance.{json,md}` (4366 →
  4368, repair 0.143 s). The regenerated values are 4438 → 4440, repair
  0.143 s, 4160 connectors, certificate 6.01 s over 2 runs. Regenerate with
  `pixi run coverage-certificate-acceptance` (a scratch copy is in
  `build/es/coverage_certificate_acceptance.md`).
- the `held-coverage-repair-gate` pixi task description in `pyproject.toml`
  ("4366 -> 4368"), which belongs to the C++ slot.

## Tests

- `test_splice_stops_early_when_no_changed_circle_reaches_a_later_decision`:
  inserting x = 2.05 re-refines 1 of 9 later sources. A later tier-1b
  decision exists but is more than ρ + 2r + H away. The result equals the
  full refinement, and so does its trace.
- `test_splice_does_not_stop_early_when_a_changed_circle_certifies_a_later_connector`:
  inserting the R = 2 circle at 2.03 re-emits the next interval identically.
  But that circle certifies the later connector 2.11 → 2.85 by tier 1b,
  which the reference refused with a midpoint, so adopting there would be
  wrong (asserted). The continuation runs to the end.
- `test_splice_without_a_reference_trace_is_a_contract_breach`: a
  trace-less reference raises `Figure5EngagementTraceError`, which is not a
  refusal subclass. A source mismatch raises `InvalidCoverageRepairInputError`.
- `test_chained_splices_equal_the_full_refinement_and_its_trace`
  (Hypothesis, 100 examples): a second splice adopts from the first splice's
  stitched trace.
- The existing `test_splice_equals_full_rerefinement` property stays green.

**Defect injections** (each in a throwaway copy):

- dropping the influence test: the must-not-stop test fails;
- dropping X from the changed set: the must-not-stop test fails;
- dropping the budget condition: the Hypothesis property fails.

The Hypothesis properties alone did not catch the first two. The dedicated
fixture is what pins them.
