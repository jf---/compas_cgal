# Stage 1: the narrow-regime generator

Task 1.1 adds the first generator of the product pipeline,
`compas_cgal.pipeline`. It machines each gate pocket with one plunge into
material: it walks a straight-skeleton chain with the narrow ladder, reaches
the next chain at cutting depth through stock the earlier chains cleared, and
lifts only where no admissible link exists or the radius-step bound splits a
chain. It is measured by the gate as the `pipeline` column, beside the frozen
columns and the `rho_regulated` reference.

At the 120° cap it meets the Task 1.1 acceptance on all three pockets:

- one plunge into material per pocket;
- no redundant operations, no degenerate loops;
- tangent breaks no more than rho's;
- no gated cell worse than `rho_regulated`'s at 120° and at 40° (the
  reference rule, which the plan's "no criterion worse" became).

The reference rule now binds `pipeline` on every commit. The stage-0 rows
(misplaced plunges, unsafe rapids, zero-length motions) are zero on every
pipeline case. How the stages work is on
[The product pipeline](../pipeline.md); this page keeps the Stage 1 record.

Task 1.2 replaces the first plunge with a helical entry and a tangent-pencil
ramp-out, aimed at the one exceedance left on every case, the full-immersion
root loop after the plunge. Its record is in
[Task 1.2: the helical entry](#task-12-the-helical-entry); its gate results
are pending.

## Measured

The pinned vectors, rendered by `tools.gate_table`. Compare the `pipeline`
column with `rho_regulated` at the same cap and pocket:

--8<-- "benchmarks/results/gate_vectors.md"

The gate measures its cases in a process pool of at most one worker per core
(`tools.gate_measure.pool_workers`), so the `generation_s` rows above are
parallel wall-clock and are not the time gate.
`pixi run pipeline-timing` is the time-gate authority; see
[Generation time](#generation-time).

## Design

Typed stages, each held to the gate's own model of the stock. The package
imports nothing from `benchmarks` and no frozen generator.

| stage | module | what it decides |
| --- | --- | --- |
| guide | `guide.py` | Ordered skeleton stations from the C++ guide. A truncated guide raises `GuideTruncatedError`. |
| direction | `direction.py`, `narrow.py` | The chain's direction is typed as `Direction2`. Each walk orientation that clears the floor is offered. The floor ρ > r and the slope guard are decided exactly. The entry radius is placed so the bridge between consecutive entries lies on the external common tangent. |
| narrow ladder | `ladder.py`, `narrow.advance_ladder` | The furthest admissible station, scanned furthest-first, never bisected (engagement is not monotone in the advance). |
| admissibility | `admissibility.py` | Exact cap verdicts at the gate's sample set, exact centre-domain containment, and the gate's removal bracket. It never mutates the stock. |
| linking | `linking.py` | The one link constructor, `link_between`. |
| emit | `emit.py`, `generate.py` | Operations and stock depletion through the gate's replay dispatch. |
| ledger | `result.py`, `errors.py` | `ForcedEntry`, `DeclinedRun`, `ForcedAdvance`, and one named exception per failure mode. |

### The narrow ladder and its condition order

A candidate station is taken when every condition passes. They run in this
order:

1. the radius step `|ρ' − ρ| ≤ 2r`, exact. `advance_ladder` prepends it
   itself; callers cannot drop it, because it is the overlap invariant;
2. rho's at-rest probe ring: the loop entry plus `LOOP_PROBE_COUNT_AT_REST`
   positions, phased on the advance and read against the stock as it stands;
3. the bridge: in the centre domain and under the cap, in process;
4. the full loop past its entry: under the cap, in process, at stations
   k ≥ 1. Station 0 is the empty prefix at the loop's exact entry, the same
   exact query on the same stock as the ring's first probe, so the
   conjunction is that of the full loop.

The ladder takes the conjunction of side-effect-free reads of one unchanged
stock, and no condition is named in any record, so the order cannot change a
verdict; it sets the cost. (The ring raises on coincident centres; native
deduplication makes that unreachable for any physical tool.) Cheapest refusal
first: the ring stops at its first refusing probe, while the bridge decides
`GATE_SAMPLES_PER_MOTION` + 1 in-process stations.

When no floor-clearing station in the window keeps the radius-step bound, the
walk splits the chain: the gap becomes a `DeclinedRun` with cause
`RADIUS_STEP` and the walk re-enters past it as a `ForcedEntry`. The bound is
never broken. rho takes `origin + 1` there and leaves an unreported uncut
patch.

### Admissibility at the gate's sample set

Loops are decided at the gate's `k/45`, k = 0..44; straights and arcs at
k = 0..45 (`GATE_SAMPLES_PER_MOTION`, pinned equal to the survey's). A decided
motion therefore cannot exceed the cap at a sample the gate reads. Every arc
piece (overrun, approach, same-circle) is decided on the endpoints of the
operation that will be emitted. A test pins that the only cap-exceeding
motions are the loop after a plunge and the bridge and loop of a recorded
`ForcedAdvance`.

The cap decisions call decision-only native queries
(`engagement_{line,arc,circle}_exceeded`): the same station placement as the
reporting prefix calls, probed last station first, exiting at the first
station that exceeds, with no reporting spans. A test fails if a decision
calls a reporting prefix entry point, and property tests pin the verdict
equal to `any()` over the prefix call's flags
(`tests/test_swept_prefix_decision.py`).

Straights (bridges and link tangents) are kept in the centre domain by an
exact segment predicate, `CutterCentreDomain2.contains_segment`: the start
point is contained, and the squared distance from the segment to every outer
and island edge is at least r², compared exactly over the injected binary64
coordinates. It replaces sampled containment, which admitted a chord whose
ends were legal but whose middle cut through the disk of a reflex corner.

### Linking

`link_between` returns one of:

- `SameCircleArc`, only when centre and authored radius are exactly equal and
  both loops turn the same way;
- `StraightLink`, the in-chain chord between neighbouring loops;
- `BeltLink` between chains: overrun arc, external common tangent, approach
  arc. Belt existence and an exact touch at the authored entry (which needs no
  arc) are decided in `Fraction`. The side follows from the turn direction
  with no comparison, and the tangent points are a construction that decides
  nothing;
- `NoLink` when the loops turn opposite ways or no external tangent exists.

A link is admissible when every piece is, and when the root loop of the chain
it reaches is under the cap. Both are judged on the pre-link stock, the stock
the ladder judges a bridge and its target loop on; it holds at least the
material of the post-link stock, so the verdict is conservative. A refused
arrival falls through to the next candidate, and to a `ForcedEntry` if none
remains.

A lift with a `ForcedEntry` of reason `NO_ADMISSIBLE_LINK` happens only when
no pending chain has an admissible link.

`choose_next` prefers an admissibly linked chain, then the shortest link, then
guide order. Links are checked lazily in (length, guide order) order and the
first admissible one is the choice; a Hypothesis test pins this equal to the
eager rule, ties included (`tests/pipeline/test_choose_next.py`).

### Emit

Each loop is authored as an `IncidentLoop` and carries its exact entry in
`ToolpathOperation.entry_xy`, which the survey reads instead of
`Circle.point_at(0)`. The entry stage is the plunge and the root loop in one
step. A loop whose removal bracket is NOTHING is not emitted. A cut-plane
motion is labelled `LINK` only when its bracket is NOTHING.

## The named exception: the at-rest ring

Plan ruling R1(a) says engagement is measured in process, in the gate and in
every generator decision. Condition 2, `Admissibility.loop_within_cap_at_rest`,
is the one exception (user ruling, 2026-10-16). It reads each ring position
against the stock as it stands, without the loop's own swept prefix.

It stays because it closes the engagement-step gaps. With in-process checks
alone, the in-process loop discounts its own swept prefix, the spine advanced
well past rho's, and the step at the bridge-to-loop junction or at a belt line
entering material put `pipeline` above rho on `max_engagement_step_deg`. The
ring is AND-ed with the in-process checks, so it can only refuse more. Three
tests pin it:

- the ring refuses a toy candidate whose bridge and loop pass in process;
- the ring short-circuits before any bridge or loop method is called;
- on a fresh generation of each gate case, `max_engagement_step_deg` stays
  within rho's committed cell.

The arriving root loop of a link is not ring-checked: a link arrival has no
advance direction to phase the ring on.

## Generation time

Fix round 1 moved the bridge's cap decision onto the gate's full sample set.
That was correct, but the bridge still ran before the at-rest ring, so every
candidate paid the full set of exact stations, and generation slowed by
several times. The root cause was condition order, not repeated work: the
stock changes on every advance, and repeated native queries are a negligible
share of cap-evaluation time. Perf round 1 runs the ring first.

Perf round 2 removed cost without changing a decision:

- **Decision-only native queries with early exit.** The ladder and link
  decisions ask only whether any station exceeds the cap, probing the last
  station first and stopping at the first exceedance (see
  [Admissibility](#admissibility-at-the-gates-sample-set)).
- **Clip monotonicity.** The in-process clip only removes pieces of the static
  engaged arcs, so every in-process run lies inside a static run, and
  "span > cap" is monotone under inclusion. A station within the cap at rest
  is therefore within the cap in process, and the clip runs only on static
  exceedances. The lemma is stated in `src/swept_prefix_engagement_2.h` and
  pinned by per-station Hypothesis agreement and a clip-flip fixture in
  `tests/test_swept_prefix_decision.py`.
- **Station 0 asked once.** The loop condition checks k ≥ 1; the ring's first
  probe is station 0 (condition 4 above).
- **Lazy `choose_next`.** Links are checked in key order and the first
  admissible one ends the search.
- **Straight domain before cap.** Inside a link's admissibility check, the
  centre-domain checks of all straight pieces run before any cap check; it is
  a pure conjunction.

The final fix wave adds one cap check per candidate link (the arriving root
loop) and replaces per-straight point sampling with one native segment query.

`pixi run pipeline-timing` (`tools.pipeline_timing`) is the time gate: every
gate generator, `pipeline` and `rho_regulated` included, one case at a time,
a fresh process per run, median of 3. It records the machine's load average
at the start and the end of the run; the table states it. Runs that disagree on the operation digest or count fail the task. Seconds are
machine-specific, so compare within one run. The table is regenerated from the
tree, so it always describes the code it ships with.

--8<-- "benchmarks/results/pipeline_timing.md"

## Known limits

- **One plunge holds because no gate pocket splits a chain.** A radius-step
  split re-enters as a `ForcedEntry` with reason
  `ForcedEntryReason.RADIUS_STEP_SPLIT`, which is a plunge. On every gate case
  the ledger (`ForcedDecisionLedger` in
  `src/compas_cgal/pipeline/result.py`) holds no such entry, and every
  `DeclinedRun` has cause `DeclineCause.FLOOR`, none
  `DeclineCause.RADIUS_STEP`. The one mid-path lift is the L's second
  `ForcedEntry`, at both caps: reason `NO_ADMISSIBLE_LINK`,
  `into_material=False`. It is a plunge into air, into stock the first chain
  cleared, so it is not a plunge into material. Pockets with necks (Task 1.7)
  will exercise the split path.
- **Link arrivals are judged on the pre-link stock.** The arriving root loop
  is decided against the cap before the link's own removal, which can only
  refuse more. It is not ring-checked (see
  [the at-rest ring](#the-named-exception-the-at-rest-ring)).
- **In-chain tangency is assumed, not certified.** The G1 guarantee covers the
  inter-chain belt. The in-chain `StraightLink` chord leaves each loop along a
  normal estimated from the neighbouring stations. The L's tangent breaks
  equal rho's; Task 1.6 owns them.
- **No headroom where `pipeline` ties `rho_regulated`.** Any rise in a cell
  where the two columns are equal in the table above fails the reference
  rule.
- **Still red in the product test.** Uncut material at declined corners (Task
  1.3), the cap exceedance of the entry loop after the plunge (Task 1.2), and
  the L's tangent breaks (Task 1.6).

## Task 1.2: the helical entry

Task 1.2 enters each pocket with one `HELIX` operation at the root loop's
exact entry, closes the bore with one flat turn, ramps out along the root's
tangent pencil to the root loop, and then walks as in Task 1.1. Its target is
`cap_exceedances` = 0 on all six cases, strictly below `rho_regulated` (1 / 5
/ 1 at 120°, 5 / 5 / 6 at 40°), where Task 1.1 sits at 1 everywhere: the
full-immersion root loop after the plunge, op 1, on every case (probe P8).
The mechanism is on [The helical entry](../pipeline.md#the-helical-entry).

!!! warning "Status: under verification"
    Every part below is built and reviewed; the gate run, the serial timing
    run and the figures are pending. Time box: two working days from
    2026-10-16, one dispatch plus at most one fix round (user decision Q2).

### What changed

- **A0, fail-loud replay classifier.** `replay_classification` decides every
  operation through one exhaustive `REPLAY_RULE_BY_OPERATION` table. A new
  `OperationType` without a rule raises `UnclassifiedOperationError`; before,
  it fell through to "rapid" silently, and depletion would have skipped the
  bore. Every rule dispatch ends in a `Never`-typed raise
  (`UnhandledReplayRuleError`), so mypy flags a missing branch.
- **Entry-stage split.** `generate.enter` became two stages:
  `enter_first_cut`, the only producer of `FIRST_CUT` (now the helix), and
  `re_enter` for `NO_ADMISSIBLE_LINK` and `RADIUS_STEP_SPLIT`. A first cut
  after anything was emitted raises `FirstCutAfterEmissionError`; a re-entry
  with reason `FIRST_CUT` raises `FirstCutReEntryError`. `enter` is deleted
  and its six test callers renamed (user go).
- **L6, ring warm start.** The at-rest ring asks its last refusing position
  first, then the rest in natural order (`ladder.WarmStart`). The ring is a
  conjunction of order-free exact probes on one stock, so verdicts and op
  digests are unchanged on all six cases; ring probes drop to 0.75–0.85 of
  1.1's. Kept by user decision. The one observable change: an internal
  invariant error can surface at a different probe than before.
- **Helix core.** `OperationType.HELIX`; `compas_cgal.helix.Helix` with the
  exact turns rule and the closing turn; `compas_cgal.helix_units`
  (`RampAngleDegrees`, `HelixPitchPerRadius`, `EntryRampLength`);
  `pipeline/pencil.py` (`TangentPencil`, `helix_footprint`, `RampOut`,
  `PencilStepCondition`). Native operation codes are decoded by
  `decode_native_operation`, which raises outside the six native kinds
  (a −1 used to wrap silently to `PLUNGE`).
- **Consumer wiring.** Every geometry dispatch outside the frozen generators
  has an explicit `Helix` branch or a named raise: replay and depletion
  (exact bore `disk(c_h, ρ_h + r)`, `HelixEvent`), the survey
  (`HelixObservation`, split into `helix_observation.py` and
  `survey_station.py` to keep `survey.py` under 1k lines), path geometry and
  marks, the held snapshot (`HeldHelixSnapshot`), operation identity, emit
  (`emit.helix_operation`, tangents as unit XY with z = 0, the repo
  convention), and the timing digest. Op digests stayed equal to 1.1's on all
  six cases until the generator emitted a helix.
- **Rules bump "1.2"** (one rebaseline record, approval "Ruling 1.2-D2",
  user-approved 2026-10-16):
    - `material_plunges` counts `PLUNGE` and `HELIX` into material; key kept;
    - `ramp_angle_violations` is newly gated, == 0, against `GATE_MACHINE`;
    - `max_engagement_step_deg` is re-baselined for the (helix, next motion)
      pair;
    - `cap_exceedances` and `misplaced_plunges` are listed because their
      definitions change (the helix is excluded from cap sampling; rule (e)
      covers the helix exit);
    - `max_loop_radius_step` is retired from the ratchet and demoted to a
      verdict row `≤ 2r` on every pipeline case;
    - new 120° stage rows: `cap_exceedances == 0`,
      `max_engagement_step_deg <= 120`, `ramp_angle_violations == 0`.
- **Generation-time rule amended.** A case passes when its 1.2 median is at
  most the slowest 1.1 run plus `C_pencil`, the pencil's own measured
  decision cost, measured in the same session.
- **Gate closure.** `GATE_PATHS` gains `helix.py`, `helix_units.py`,
  `held_path_snapshot.py`, `tools/gate_measure.py`, `benchmarks/models.py`
  and the two survey splits; a new test keeps `NOT_GATE` and `GATE_PATHS`
  disjoint.

### User decisions and their reasons

| decision | ruling | reason |
| --- | --- | --- |
| Q1 beat-rho cell | `cap_exceedances`, both caps, all pockets, target 0 | The helix exists to remove exactly the one remaining exceedance; 0 is strictly below rho on all six cases. The other tied cells (L step, L tangent breaks, 40° rect_20x12 uncut) stay must-not-exceed-rho. |
| Q2 time box | 2 days, one dispatch plus at most one fix round | Focus rule: land when the gate holds; when the box runs out, bring the measured state. |
| Q3 bump contents | material entries incl. helix, ramp violations gated, step re-baselined, in one record | Keeping the `material_plunges` key keeps 1.1's row truthful; a ramp limit that is only reported cannot fail a commit. |
| Q4 plunge junctions | not gated; declared `EntryRampLength` = 2r, θ_e = cap·π/45 at station 1, pencil and helix junctions only | Gating plunge junctions makes every step cell 180° (probe P7), a vacuous comparison. θ_e is declared intent with an exact surrogate injected once, not a number chased from rho's cell. |
| Q5 root | 1.1's `FIRST_CUT` station; medial-axis vertex to Task 2.5 | The guide-widest root on the L walks chain 0 from a cleared corner and worsens the L step 10.00 → 34.38° at 120° (probe C8-A). |
| Q6 bore bottom | one flat closing turn at ρ_h in the same `HELIX` op; helix and closing turn excluded from cap sampling; bore removed at the end | Without it the replay removes a full disk while a helical wedge up to one pitch remains. |
| D1 ramp angle | `max_ramp_angle_deg` = 2°, placeholder | No datasheet value in the repo; flagged at the parameter. |
| D4 loop step | `max_loop_radius_step` demoted to a `≤ 2r` verdict row | The pencil's steps (0.34 on the 120° L) exceed 1.1's 0.106 but pass exact coverage `|Δc| + Δρ ≤ 2r` with uncut unchanged; a ratchet on the cell would be number-chasing. |
| D6 generation time | GT rule amended with `C_pencil` | Even the admitted members' own cap decisions exceed 1.1's run-to-run margin on 5 of 6 cases: new work, not a regression. |
| M8 type gates | `helix.py`, `helix_units.py`, `pencil.py`, `helix_observation.py`, `survey_station.py` under `mypy --strict` | New package code meets the engineering bar. The two `benchmarks/` modules join `types-benchmarks`, which CI runs; that change is held for a separate go. |
| L6 | keep the ring warm start | Order-free verdicts, 15–25% fewer ring probes. |

### Gate results

Pending the 1.2 gate run: the quality cells are filled from
`benchmarks/results/gate_vectors.md` after it. "1.1" is the Task 1.1
commit's cell; "rho" is `rho_regulated` at the same cap and pocket.
`generation_s` is not in this table: its
verdict (D6) is generated, in the generation-time gate below the table, by
`python -m tools.pipeline_time_gate` from three records. Its 1.1
base is the slowest of the three serial runs in the `pipeline_timing.json`
committed with the Task 1.1 tree (7b001d47, the default; another base only by
an explicit flag), read through `git show`, or a same-session re-timing of
that tree passed in explicitly. It is not the uncommitted pre-generator
scratch run. Its 1.2 medians come from
`benchmarks/results/pipeline_timing.json`. Its `C_pencil` comes from
`benchmarks/results/pipeline_pencil_cost.json`. The gate refuses a verdict
unless the 1.1 base is a same-session re-timing, and it checks this from the
records, not from a flag. `tools.pipeline_timing_1_1_base` times the
unmodified 7b001d47 worktree through the 1.1 tool's own per-case child mode,
with the 1.1 tree's own editable build pinned, after first bringing that build
up to date with the 1.1 tree's own `_editable-rebuild` task (a failing rebuild
refuses the run with `BaseRebuildFailedError`), and refuses an interpreter
whose `compas_cgal` or `_stock_2` lies outside that worktree. It records the
load and the provenance. Every record carries its provenance:
commit, modification digest (every tracked and untracked change except the
gate's own output records under `benchmarks/results/`, `GATE_OUTPUT_RECORDS`
in `tools/run_provenance.py`, so writing them does not change the tree's
identity), tool digest, `pixi.lock`, environment, Python,
`_stock_2` digest, host, CPU and UTC interval. The gate verifies six things:
the base's commit is 7b001d47, the base is unmodified and written by the
wrapper; the base did not time the 1.2 library (its `_stock_2` digest must
differ when the native sources do); the 1.2 timing and pencil records each
describe one unchanged tree, the same one; all three ran on one machine in
one environment with one lock; they ran in the order base, 1.2 timing, pencil, each starting within
15 minutes of the previous end; and every load is within the quiet limit
(1-minute and 15-minute load 2). A loaded base inflates its slowest run and
loosens the bound. The committed 7b001d47 record has no load or provenance,
so with it the page is pending or indicative, never a verdict. `C_pencil` is
measured on the code it excuses: a slowdown inside the pencil's decisions,
including the cap machinery they share, is credited automatically, so the
gate bounds only the time outside them.

| case | criterion | pipeline 1.2 | pipeline 1.1 | rho |
| --- | --- | --- | --- | --- |
| 120 rect_12x8 | cap_exceedances | pending | pending | pending |
| 120 rect_12x8 | max_engagement_step_deg | pending | pending | pending |
| 120 rect_12x8 | max_loop_radius_step | pending | pending | pending |
| 120 rect_12x8 | uncut_fraction | pending | pending | pending |
| 120 rect_20x12 | cap_exceedances | pending | pending | pending |
| 120 rect_20x12 | max_engagement_step_deg | pending | pending | pending |
| 120 rect_20x12 | max_loop_radius_step | pending | pending | pending |
| 120 rect_20x12 | uncut_fraction | pending | pending | pending |
| 120 L_shape | cap_exceedances | pending | pending | pending |
| 120 L_shape | max_engagement_step_deg | pending | pending | pending |
| 120 L_shape | max_loop_radius_step | pending | pending | pending |
| 120 L_shape | uncut_fraction | pending | pending | pending |
| 40 rect_12x8 | cap_exceedances | pending | pending | pending |
| 40 rect_12x8 | max_engagement_step_deg | pending | pending | pending |
| 40 rect_12x8 | max_loop_radius_step | pending | pending | pending |
| 40 rect_12x8 | uncut_fraction | pending | pending | pending |
| 40 rect_20x12 | cap_exceedances | pending | pending | pending |
| 40 rect_20x12 | max_engagement_step_deg | pending | pending | pending |
| 40 rect_20x12 | max_loop_radius_step | pending | pending | pending |
| 40 rect_20x12 | uncut_fraction | pending | pending | pending |
| 40 L_shape | cap_exceedances | pending | pending | pending |
| 40 L_shape | max_engagement_step_deg | pending | pending | pending |
| 40 L_shape | max_loop_radius_step | pending | pending | pending |
| 40 L_shape | uncut_fraction | pending | pending | pending |

--8<-- "benchmarks/results/pipeline_time_gate.md"

### Known open defect: CW arcs stored as their complement

`emit.arc_operation` stores a clockwise arc as `Arc(ρ, 0, ccw angle S→E)` in
a +Z frame plus a `clockwise` flag, so its compas point set is the CCW
complement: 270° for a 90° fillet, the opposite semicircle for the L's 180°
half-turns. Every pipeline `Arc` is clockwise (3, 3 and 6 LINK arcs per
pocket); the frozen columns emit none.

- **Gated cells are unaffected**: 0 of 12 gated criteria change in 24 of 24
  cases when the arcs are rebuilt with the true CW point set. Every cutting
  decision reads the endpoints and the flag, never the compas
  parametrisation.
- **Reported cells and figures are wrong**: path, cutting and air-at-depth
  length (Δ = πr per 90° fillet, e.g. 120° rect_12x8 path length 358.560 →
  348.824), `recut_fraction`, the length-weighted engagement statistics, and
  the drawn arc midpoints (off by about 2 mm on fillets, 4 mm on half-turns).
  Off the gate path, `engagement._certify_arc_engagement` and
  `exceedance._sample_positions` read the wrong points too.
- **Not fixed in 1.2.** The fix changes every pipeline digest, reported cell
  and figure, so it is its own Stage 1 commit, scoped by the user.

### Known limits of Task 1.2

- **Re-entries stay plunges.** No gate pocket re-enters material; a pocket
  that does would need a helix at every material entry.
- **The root is the 1.1 start station**, not the maximum-clearance point
  (Task 2.5).
- **The pencil assumes a simply connected domain.** A holed domain needs an
  exact per-member containment condition (Task 1.7).
- **θ_e is sampled at station 1**, nothing between; the arc-length gate is a
  separate rules amendment.
- **Ramp angle and entry ramp are placeholders** until tool and machine
  values are supplied.

## Figures

`pipeline` beside the frozen generators and `rho_regulated`, each panel
captioned with its pinned gate cells. `pixi run gate-figures`
(`benchmarks/gate_figures.py`) draws them; it refuses a caption the current
code no longer reproduces, and the PNGs are byte-reproducible: redrawing the
same paths rewrites the same bytes (`benchmarks.plotting.save_figure` pins
the metadata and the resolution). [The product pipeline](../pipeline.md)
embeds the same files.

### 120°

![rect_12x8, 120°](../assets/images/pipeline/120/rect_12x8.png)

![rect_20x12, 120°](../assets/images/pipeline/120/rect_20x12.png)

![L_shape, 120°](../assets/images/pipeline/120/L_shape.png)

### 40°

![rect_12x8, 40°](../assets/images/pipeline/40/rect_12x8.png)

![rect_20x12, 40°](../assets/images/pipeline/40/rect_20x12.png)

![L_shape, 40°](../assets/images/pipeline/40/L_shape.png)
