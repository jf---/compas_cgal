# Review: actual-stock engagement in refinement (design, 2026-11-04)

Reviewed: `docs/superpowers/specs/2026-11-04-actual-stock-refinement-design.md` at HEAD `ea4703bb`.
Probes: `build/actual-stock-review/` (`probe_decisions_contact.py`, `probe_contact_material.py`, and their logs),
run niced through `build/pyrun.sh` (SKBUILD_EDITABLE_SKIP), at load 77–127, so the numbers are indicative.

The inclusion argument is sound for the polygon lane, provided that the material domain is the part
polygon and the filter keeps the orientation leg of the cap test. As written, the spec leaves both open,
and either gap admits a false acceptance. The plan's gates would not catch a mis-wired conjunct until
U5. No ledger ruling is contradicted.

## High

### H1. The material domain is named "Ω", but in this codebase Ω is the tool-centre (offset) polygon

- The spec decides on "rim ∩ int Ω ∖ ⋃Aⱼ" with "Ω's edges" and "the polygon vertices" as inputs
  (design.md:192-204).
- The only polygon refinement holds is `_Emission.boundary` (`benchmarks/held_figure5_engagement_refinement.py:89`).
  The generators fill it with `components[0]`, the one-radius **inward offset**
  (`benchmarks/held_toolpath_benchmark.py:306`, `:388-392`; `benchmarks/held_reference_figures.py:493-506`).
- The Ω\* module states that it "replaces the float offset polygon by Ω\*"
  (`benchmarks/held_omega_toolpath.py:1-8`).
- The certificate's material is the part polygon `case.projection.points`
  (`tools/actual_stock_profile.py:193`, `fresh_stock` `:207`).

**Failure scenario.** An implementer binds `self.boundary`. The engine's stock then lacks the band of
width r along every wall. Station 0's rim touches the part edge on almost every circle (measured below,
M3), so the engine sees no wall material at station 0. It accepts a circle whose station 0 exceeds the
cap on the certificate's stock. The by-construction claim (design.md:263-272) fails silently. Only D2
catches it, and D2 sits in U5 (see H3).

**Fix.**

1. Name the material domain as the part polygon P = `case.projection.points`.
2. Plumb it into `refine_figure5_engagement`, `refine_boundary_ordered_path` and `EngagementRefinementKey`
   (`benchmarks/held_figure5_engagement_trace.py:186-195`; otherwise a splice can adopt across different
   domains).
3. Name the Ω\* lane's material domain explicitly.
4. Add the R6 injection "material domain = offset polygon".

### H2. The filter's "within" test as worded drops the θ > π leg

The spec certifies **within** when "the upper chord² [is] certified below the cap surrogate"
(design.md:238). Chord² is not monotone past π. The exact test decides orientation first:

- `start == end` means exceeded;
- orientation NEGATIVE means θ > π, which means exceeded;

(`src/exact_cap_certificate_2.h:84-120`).

**Failure scenario.** Take the 80° cap, so the cap ratio is 4·sin²40° ≈ 1.65. A station whose
pessimistic merged run spans 300° has chord² 4·sin²150° = 1.0, which is below the cap. The filter then
returns "within" for a station the exact path calls exceeded. A pessimistic merge across ambiguous gaps
makes long runs more likely, not less.

**Fix.**

- "Within" needs the certified orientation sign > 0 and start ≠ end as well as the chord² bound.
- "Exceeded" may come from either leg.
- Add the R6 injection "orientation leg removed from the filter".

### H3. U4's exit gate is satisfied by a vacuous or mis-wired conjunct; D2 is not in it

U4's exit gate is "gate-vector and polygon digests IDENTICAL" plus generation time (design.md:371). The
design expects identical digests, because the circle view refused nothing on the samples
(design.md:308-315). So the gate passes the following:

- a conjunct that is never reached;
- one decided on the wrong view (off by one, or with A_t included);
- one wired into `close_seam` but missing from Ω\* repair `:254` (`benchmarks/held_omega_connector_repair.py:254`).

The unseen-hole red test checks one fixture only. D2, which is the production ⇒ certificate property
itself, waits for U5's slow tier (design.md:282, `:372`).

**Fix.** Make D2 (fig5 and skis every circle, Monstera every 20th) and a call-count assertion (engine
decisions ≥ final circles − 1) part of U4's exit gate.

## Medium

### M1. The splice footprint fix names the wrong set

The spec says early stop is unsound unless "any **inserted** annulus whose box meets a later footprint
stops adoption" (design.md:174-179). For P, insertion is harmless: an added annulus only shrinks later
stocks, which is monotone. The unsound case is the other half of `changed`, the reference circles Y that
the continuation no longer emits. "Connectors after s+ … may need fewer midpoints"
(`benchmarks/held_coverage_repair_splice.py:155-158`; `changed` includes `reference.circles[start+1:end-common]`, `:110`).

**Failure scenario.** A reference midpoint y is dropped. A later adopted circle t was accepted on a view
that contained A_y. The final stock at t is larger than the view t was accepted on, so t can exceed.

There is a second gap. Footprints are recorded today only for decisions that read beyond their pair
(`benchmarks/held_figure5_engagement_trace.py:170-178`). After the change, every decision that reaches the
engine reads its view, including the refusals that inserted midpoints. Those refusals can flip when X
adds annuli.

**Fix.**

- Record the read disk B(c_t, ρ_t + r) ∪ B(e_t, r) on every engine-reaching decision.
- Test it through `last_influenced` over X ∪ Y.
- Add a splice test with Y ≠ ∅ inside a later read disk.
- Add the R6 injection "footprint tested on X only".

### M2. Filter bracketing of the own-prefix clip is unspecified

The pessimistic and optimistic rules cover harvest events only (design.md:233-241). The clip (start disk
plus sector, `src/swept_prefix_engagement_2.h:1-23`) has its own near-tangent events. An upper bound on
max-run needs under-clipping on ambiguity, and a lower bound needs over-clipping. The spec does not
state this, and "reused unchanged" (design.md:201) implies an exact clip on interval arcs, which does
not type-check.

**Fix.** State both directions for the clip, and add an injection.

### M3. The fallback risk is measured after the engine is built; the population is pervasive

Risk 2 is measured in U3 (design.md:370, `:377-378`), after U1 and U2 have committed the C++ engine.
This review measured the population directly (`build/actual-stock-review/probe_{fig5,skis}.log`,
`probe_contact_material.log`), using an exact Fraction sign of dist²(e, nearest part edge) − r².

| | fig5 | skis |
|---|---|---|
| rim crosses the wall | 1,023 | 839 |
| rim exactly tangent | 37 | 0 |
| rim clear of the wall | 944 | 782 |
| \|d² − r²\|/r² < 1e-12 | 1,902 / 2,004 | 1,560 / 1,621 |
| \|d² − r²\|/r² p50 | 9e-15 | 5e-15 |
| near-tangent contact point outside every prior annulus | 1,898 / 1,901 | 1,555 / 1,559 |

So about 95 % of station-0 rims straddle a sign-ambiguous wall event, with circle-view material on both
sides. The pessimistic merge and the optimistic split then differ by the whole station-0 run, and only
the cap margin of that run keeps the filter off the exact path. The "≤ 0.1 %" (design.md:341) has no
measurement behind it.

The same data shows that on about half the circles the exact verdict splits the station-0 run at a
sub-ulp wall crossing. That is the sliver class the spec cites as the reason to exclude connectors
(design.md:172-173). It is a P semantics question for the user, not an inclusion defect.

**Fix.** A U0 probe on existing `Stock2`: per station 0, compare the exact max-run, the run merged
across the wall crossing, and the cap. Do this before U1.

### M4. Ω\* has no production ⇒ property test, and the store cannot serve Ω\* repair

- D2 is polygon-only: `ReplayWorkload.from_path` takes a `GeneratedToolpath`
  (`tools/actual_stock_profile.py:185-202`). Ω\* gets D1 only (design.md:384-385). D1 compares the engine
  with `Stock2` on the same inputs, so it cannot detect wrong wiring in `refine_omega_*` or in the repair.
- A circle-only replay of the final Ω\* path is feasible, because the circle-view stock is the decision
  itself. Add it as D2-Ω\*.
- Ω\* repair runs after NOTHING on the final list and re-streams `emitted`
  (`benchmarks/held_omega_toolpath.py:248-253`; `benchmarks/held_omega_connector_repair.py:250-262`).
  An `AnnulusStore` keyed by refinement `EmissionIndex` (design.md:162-165) does not index that list, so
  the store must be rebuilt per pass. The spec does not say so.

### M5. Budget unit costs are not traceable

These figures cite no measurement:

- 50 ns per interval intersection;
- 15 µs of call overhead;
- 1 ms per exact fallback;
- ≤ 0.1 % fallback;
- 0.15 ms per full-view escalation

(design.md:336-342). The design's own probe log shows that full views hold up to 313 prior circles, for
example circle 680 (`build/actual-stock-design/probe_fig5.log`: "sweeps 546 circles 313"). That is
roughly 300 annuli × 45 rims per escalation, not 0.15 ms.

The fig5 estimate's upper end, 0.20 s, exceeds the budget the spec itself sets, ≤ 0.15 s
(design.md:349-350). Skis' upper end equals its budget.

Two assumptions now have measurements (this review):

- engine calls per final circle: 2,010 / 2,004 on fig5 and 1,623 / 1,621 on skis (gate calls under the
  spec's conjunct order);
- circle-view last-1 refusals: 1 / 101 on fig5 and 0 / 108 on skis.

Monstera is unmeasured.

## Low

- **L1. Station 0 is not on the circle.** §2.3 says "stations are c + ρ·u" (design.md:203). Station 0 is
  the entry e (`src/engagement_2.cpp:780-788`). The footprint does include B(e, r); the event and witness
  text should say so too.
- **L2. Rim order.** Sorting "with `compare_xy` / orientation" (design.md:196) does not give the angular
  order on the rim. Use the upper and lower split of `maximal_exact_runs`
  (`src/audit_exact_station_2.cpp:129-151`).
- **L3. Pin the cap ratio source.** The cap ratio must be `EngagementCap.chord_ratio`
  (`src/compas_cgal/motion.py:78-83`), the same native `cap_chord_ratio` the certificate uses
  (`tools/actual_stock_profile.py:203`).
- **L4. Truncation cost.** "Truncate … an O(1) index drop" (design.md:174) does not hold for a dynamic
  `Point_set_2`. Filter by index at query time, or state the real cost.
- **L5. Unit discipline.**
  - U5 has no red test.
  - U2 has no time gate for the index.
  - The R6 list is global, not per unit.
  - U0 and U5 add or change pixi tasks, which are pyproject edits. The spec gives the only pyproject edit
    to the C++ units (design.md:362, `:368`, `:372`).
- **L6. Ω\* cost.** Ω\* gains cost on a lane that is already red (design.md:356). That is consistent
  with performance-decisions (b) and (c) (`docs/superpowers/state/performance-decisions-2026-10-22.md:62-71`),
  but no number is given.

## Verified, no finding

- **Station equivalence.** Same `contact_point`, `loop_parameters()` (k/45, k = 0..44,
  `src/compas_cgal/pipeline/admissibility.py:75-90`), clockwise False, and the same annulus construction
  (`full_turn_annulus_bounds`, `src/stock_2.cpp:427-433`) as the certificate (`tools/actual_stock_profile.py:219-226`).
  Each circle's entry is its own attribute, so a NOTHING removal or an insertion never moves another
  circle's stations.
- **NOTHING.** The removal decision is A_k ⊆ ⋃ of the kept earlier closed annuli in ℝ²
  (`benchmarks/held_nothing_circles.py:1-17`, `:653`). A regular closed annulus that is regularized-contained
  is pointwise contained, so every later view is pointwise unchanged.
- **Insertion.** Insertions, both coverage-repair X and Ω\* midpoints, only shrink later views
  (monotone).
- **Circle 0.** It is cut first and never re-cut (`benchmarks/held_figure5_engagement_refinement.py:139-146`),
  and the certificate refuses it at 360°, so the exemption matches.
- **No cross-root arithmetic.** Station centres are exact rationals, the rim radius² is rational,
  witnesses are rational, and membership is a rational squared-distance compare. Events are compared
  only, never combined.
- **No epsilon.** There is no epsilon in the decision path. Boxes are outward-rounded and only select
  candidates.
- **Ledger.** Q1 (circle view), Q2 (rim engine replaces S1–S5) and Q3 (45 stations) are followed as ruled
  (progress.md:1173-1175). The C1-assert drop is flagged as a user gate (design.md:372; progress.md:1176).
- **F5 dependency.** The F5 dependency is real: fig5 polygon generation is 1.72 s against 1.7 s
  (performance-decisions:66). The lever headroom of 0.25 / 0.32 / 2.0 s is arithmetic on that doc's
  table (`:66-67`).
