# Radius-Regulated Trochoidal Toolpath

`compas_cgal.engagement_radial_toolpath.radius_regulated_toolpath` adds a second regulated knob to
the advance regulation of [`engagement_controlled_toolpath`](engagement_controlled_toolpath.md): the
loop radius, chosen at each station by a **ladder search** over the same exact `cap_exceeded`
predicate. Measured on the 20x12 reference pocket with a 2 mm tool it **cuts the number of machining
circles shown to exceed the cap by 1.6x to 2.7x** at every cap the advance regulation alone cannot
meet, at a **21% generation cost at a loose cap** and reproducing the advance-only generator's
operation stream byte-for-byte when the cap is loose enough that the ladder never leaves its top
rung.

It does **not** deliver the headline the knob was expected to buy. The worst loop engagement away
from the chain entries stays pinned at 125.8 deg for caps of 40, 60 and 80 deg — the same saturation
the advance-only generator shows. The section
[below](#the-negative-result-the-worst-loop-does-not-track-the-cap) says exactly why, with the
decomposition that proves it, because that finding is worth more than the improvement.

Over the circles the generator **accepts**, the picture is different, and it changed on 2026-08-21
when `LOOP_PROBE_COUNT` replaced the advance-facing probe triple with a uniform ring: the worst
accepted circle now sits within **3.4 deg of the requested cap** at every cap from 40 to 100, against
up to 46.7 deg before. Cap tracking on accepted circles is a probe-placement result, not a radius
result — both generators land on the same accepted-circle maximum.

!!! warning "What is guaranteed, and what is not"

    Engagement `<= tea_cap_deg`, decided by an exact predicate, **at each evaluated tool position**.

    That is the whole claim, and it is the same claim the advance-only generator makes — no more. It
    is **not** a continuous guarantee between evaluated positions. **No certificate is produced or
    returned**, and no `MotionWitness` / `CapRefutation` object exists on this path. The bridge cuts
    between machining circles are not regulated at all, by either generator. Chain-entry loops are a
    full slot by construction and are counted and warned about, never hidden.

## Engagement is not monotone in the loop radius

This is the finding that determines the algorithm. Measured at one mid-path station of the 20x12
pocket (centre 14.0, 6.0; 2 mm tool; stock depleted by eight maximal loops a half tool diameter apart
behind it), the peak engagement over the loop as the radius shrinks:

| rung | 0 | ... | 15 | 16 | **17** | 18 | 19 | 20 | ... | 35 | 36 | ... | 39 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| radius | 4.998 | | 4.248 | 4.198 | **4.148** | 4.098 | 4.048 | 3.998 | | 3.248 | 3.198 | | 3.048 |
| peak TEA (deg) | 189.6 | | 89.8 | 80.1 | **69.2** | 72.4 | 77.0 | 85.0 | | 73.9 | 65.5 | | 31.2 |

Engagement falls to a local minimum at rung 17 and **rises again below it** — one guide step of extra
retreat carries the loop back over a 70 deg cap it had just met. At that cap the admissible rungs are
`{17}` and then `{36, 37, 38, 39}`, with rungs 18 through 35 refused **between** them: the admissible
set is **not an up-set in the rung index**, and a bisection — which is correct only on an up-set —
returns rung 36 (radius 3.198) where the scan returns rung 17 (radius 4.148), a circle nearly a
millimetre larger that is equally admissible.

!!! note "Re-pinned when the probe ring replaced the triple"

    The state this table is measured on was rebuilt on 2026-08-21. The previous one — three loops at
    0.5 spacing, read at a 60 deg cap — had **no** admissible rung left once the probes could see the
    loop's trailing side, so it demonstrated nothing rather than demonstrating something false. The
    phenomenon is unchanged and so is the conclusion; only the state exhibiting it moved.

!!! danger "Never bisect the radius ladder"

    `_largest_admissible_radius` **scans** from rung 0 downward and returns the first rung that
    passes. That is the largest admissible radius under *any* pass/fail pattern, where a bisection is
    correct only under one. `tests/test_engagement_radial_toolpath.py::test_engagement_is_not_monotone_in_the_loop_radius`
    runs a bisection beside the scan on this state and asserts they disagree, so the assumption
    cannot be reintroduced silently.

This is the **second** instance of the same structure in this repository. Spacing does not order
engagement either — see `benchmarks/figure6.py`, "Spacing does not order engagement", which is why
the constant-spacing baseline is a minimum over a brute-force sweep rather than a bisection.

## How a radius is chosen

```mermaid
flowchart TD
    A["guide station<br/>centre + maximal clearance radius<br/><i>statically determined</i>"] --> B{"rung 0:<br/>the maximal circle"}
    B -->|"cap not exceeded"| E["emit; station FINISHED<br/><i>identical to the advance-only generator</i>"]
    B -->|"cap exceeded"| C["scan rungs 1..39<br/><b>never bisect</b>"]
    C --> D{"exact cap_exceeded<br/>at 33 evaluated positions<br/><i>runtime, depleting stock</i>"}
    D -->|"exceeded"| C
    D -->|"passes"| F{"Stock.contains at the<br/>loop's deepest reach<br/><i>does it still cut?</i>"}
    F -->|"cuts nothing"| C
    F -->|"cuts"| G["emit reduced circle<br/>station NOT finished"]
    C -->|"ladder exhausted"| H["forced: emit the maximal circle,<br/>count it, warn"]
    G --> I["sweep the chain again"]
    I --> A
```

A rung is admissible on **two** exact conditions, and both are load-bearing:

1. **no evaluated position exceeds the cap** — the same `_stock_2.engagement_at` verdict, at the same
   33 positions (entry point plus a uniform 32-ring), that the advance-only generator uses; and
2. **some evaluated position still reaches uncut stock** — `Stock.contains` at the loop's deepest
   reach, the tool centre pushed one tool radius outward along the ray from the station centre.

### Why condition 2 exists

Condition 1 alone is a trap, and it was found by measurement rather than by inspection. In the steady
trochoidal regime the material a loop meets is a crescent at its **outer** rim, so a loop drawn
inward is not a lighter cut — it is **no** cut: it spins inside the annulus its predecessors already
swept, engages nothing, and therefore exceeds nothing. A scan that takes the largest such loop never
moves the frontier, so the next station admits a still smaller loop, and the ladder walks itself to
the bottom one rung per station before falling back on the maximal circle in virgin stock.

Measured on a 10x6 pocket at a 40 deg cap, that produced this sequence — radius collapsing one rung
per station and then a 360 deg circle the advance-only generator never emits:

```text
op 77   CIRCLE r=0.0980 c=(4.900,3.000)
op 78   LINE   (4.900,3.098)->(4.950,3.048)
op 79   CIRCLE r=0.0480 c=(4.950,3.000)      <- ladder bottom
op 80   LINE   (4.950,3.048)->(5.000,4.998)  <- TEA 180.1
op 81   CIRCLE r=1.9980 c=(5.000,3.000)      <- TEA 360.0, forced maximal circle
```

`tests/test_engagement_radial_toolpath.py::test_no_machining_circle_away_from_a_chain_entry_is_a_full_slot`
pins the repair.

## The negative result: the worst loop does not track the cap

The hypothesis under test was that the worst loop engagement away from the chain entries would fall
towards the requested cap once the radius became a regulated variable. **It does not.** On the 20x12
reference pocket, 2 mm tool, engagement measured after each chain's virgin-stock entry cut:

| cap (deg) | max loop TEA, advance-only | max loop TEA, radius-regulated | over-cap circles, advance-only | over-cap circles, radius-regulated | length, advance-only | length, radius-regulated |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | 125.8 | 125.8 | 112 | **68** | 6519.2 | 6527.5 |
| 60 | 125.8 | 125.8 | 52 | **24** | 3335.9 | 3348.9 |
| 80 | 125.8 | 125.8 | 32 | **12** | 2222.9 | 2325.6 |
| 100 | 125.8 | **100.2** | 8 | **4** | 1348.7 | 1409.7 |
| 120 | 114.5 | 114.5 | 0 | 0 | 1057.3 | 1057.3 |

Measured by replaying each path against a depleting exact stock and walking every non-entry machining
circle at 60 tool-centre positions, phase-offset by half a step so the walk shares no grid with the
generator's probes. "Over-cap circles" counts circles with at least one position over the cap: a
sampled **lower bound** on exceedance, never a certificate.

Restricted to the circles each generator **accepts** — dropping the ones it refused and emitted
anyway with a warning — the same walk reads:

| cap (deg) | accepted max, advance-only | accepted max, radius-regulated | *before the probe ring*, advance-only |
| ---: | ---: | ---: | ---: |
| 40 | 43.4 | 43.4 | 86.7 |
| 60 | 62.5 | 62.5 | 95.1 |
| 80 | 81.7 | 81.7 | 101.7 |
| 100 | 99.9 | 100.2 | 108.9 |
| 120 | 114.5 | 114.5 | 114.5 |

Two things follow. The cap **is** tracked on accepted circles, to within 3.4 deg — and that came from
the probe ring, not from the radius knob, because the two generators agree column for column. And the
125.8 deg in the first table is entirely the **forced** regime: circles the generator refused and had
to emit anyway.

### Why it does not track — the decomposition

Grouping the same paths' loops by radius shows where the 125 deg lives. At a 60 deg cap:

| loop radius | count | max TEA (deg) | what this is |
| ---: | ---: | ---: | --- |
| **0.057** | 4 | **126.0** | corner spoke **at the apex** — clearance barely exceeds the tool radius |
| 0.516 | 4 | 107.7 | corner spoke, near the apex |
| 1.327 | 4 | 64.7 | corner spoke, mid-clearance |
| 2.033 | 4 | 63.9 | corner spoke, mid-clearance |
| every other radius, 4.998 included | — | < 61.7 | central chain at full clearance and the rest of the spokes |

The pinned maximum is a corner-spoke station whose clearance-derived radius is 0.057 mm against a
1.0 mm tool radius. The loop is a near-point; the tool is wedged into the corner apex with material
on every side. Its ladder has **one rung** — rung 1 would be negative — so there is no radius to
choose. This is geometry, not a control failure: no loop radius exists there that engages less. It is
also not a sampling artefact: a near-point loop is smaller than the probe spacing at any density, so
no probe count reaches it either.

The improvement is real where a radius choice does exist — the over-cap circle count drops 32 to 12
at an 80 deg cap and 52 to 24 at 60 — it simply never touches the apex stations that set the maximum.

### The residual over-cap circles are forced ones, or marginal

On the small pockets every over-cap circle is a forced one. Attributing each over-cap circle to how
the ladder chose it — forced, accepted at the maximal radius, or accepted at a reduced rung — with a
60-position walk:

| pocket, cap | forced | accepted, maximal | accepted, reduced |
| --- | ---: | ---: | ---: |
| 6x4, 40 deg | **34** | 0 | 0 |
| 10x6, 40 deg | **61** | 0 | 0 |
| 20x12, 40 deg | 44 | **20** | 0 |
| 20x12, 60 deg | 8 | **20** | **4** |

`test_tight_cap_regulates_the_radius_and_leaves_fewer_circles_over_the_cap` asserts on 6x4 that no
*reduced* circle is measured over the cap.

The 20x12 rows are the residual sampling gap, and they are **marginal**, not gross: the accepted
circles on those rows peak at 43.4 deg against a 40 deg cap and 62.5 against 60. That is the gap
between 33 evaluated positions and a 60-position walk, and it is what
[`LOOP_PROBE_COUNT`](engagement_controlled_toolpath.md#why-a-uniform-ring-and-why-32-of-them)
narrows rather than closes. Before the probe ring the same rows peaked at 86.7 and 95.1.

The consequence for anyone extending this: **the radius knob is spent.** The binding constraints are
now the guide's own clearance profile at corner apexes and the unregulated bridge cuts, whose worst
engagement (141.8 deg) is unchanged by anything in this module.

## Why more passes, and why the path gets longer

A loop smaller than its station's maximum sweeps a narrower annulus, leaving that station's outer
band uncut. The chain is therefore walked **repeatedly**, and a station is finished only once its
**maximal** circle has been emitted — which is exactly what the unregulated walk does on its first
and only pass. Tighter caps force smaller first loops, hence more passes, hence longer paths. On the
20x12 pocket that shows as 9 plunges over 5 chains at caps of 40 to 100 against 5 at 120, and as the
length column rising 1057 → 6528 as the cap tightens.

Termination is **structural, not budgeted**: `_radius_ladder` offers a station only radii strictly
above the largest it has already emitted, so each emission climbs at least one rung of a 40-rung
ladder, and a pass that emits nothing ends the chain. That argument bounds the passes by
`stations x rungs`, which is finite but useless as a guard, so `MAX_RADIAL_SWEEPS_PER_CHAIN` is an
empirical budget on top of it and it **raises** rather than truncating a chain, because a truncated
chain is unmachined material presented as a finished toolpath. Measured over 6x4, 10x6, 12x8 and
20x12 pockets at caps from 20 to 170 deg, **no chain ever needed more than two passes**; the budget
is 40.

Stations the advance search **jumped over** are marked finished too, on `MAX_ADVANCE_TOOL_DIAMETERS`'
own coverage argument: two maximal circles no more than a tool diameter apart sweep overlapping
annuli. That argument only holds behind a maximal circle, so after a **reduced** one the walk
advances a single station and marks nothing finished.

## Evidence table

| Claim | Value | Source |
| --- | --- | --- |
| Ladder spacing | `D/40 = r/20`, the advance grid's own step | `RADIUS_LADDER_STEP_TOOL_DIAMETERS` |
| Ladder span | one tool diameter (TEA saturates at `ae = 2r`) | `RADIUS_LADDER_SPAN_TOOL_DIAMETERS` |
| Rungs | 40 | `RADIUS_LADDER_RUNGS` |
| Non-monotone in radius | 80.1 → 69.2 → 72.4 → 85.0 deg down four consecutive rungs | `test_engagement_is_not_monotone_in_the_loop_radius` |
| Bisection is unsound here | returns rung 36 where the scan returns rung 17 | same test, `_bisect_rung` beside the scan |
| Loose cap reproduces the advance-only stream | byte-identical operation signature, 6x4 at 120 deg | `test_loose_cap_reproduces_the_advance_only_generator` |
| Over-cap circles, 20x12 | 112→68, 52→24, 32→12 at caps 40/60/80 | 60-position replay walk, entry cuts excluded |
| Max loop TEA, 20x12 | unchanged at 125.8 deg for caps 40, 60 and 80 — all of it forced circles | same walk |
| Max loop TEA over ACCEPTED circles, 20x12 | 43.4 / 62.5 / 81.7 deg at caps 40/60/80, both generators | same walk, split by the generator's own verdict |
| Generation, 12x8, cap 120 | 271.6 ms against 223.7 ms (best of five) | `time.perf_counter` around each generator |
| Generation, 12x8, cap 40 | 7784 ms against 1676 ms | same |
| Coverage, 6x4 and 10x6, caps 40/60/120 | zero uncleared grid points, both generators | `Stock.contains` grid, wall distance > tool radius |

## Rejected alternatives

**Uniform-ring probing for reduced circles only.** *Superseded 2026-08-21 — recorded because the
reasoning that rejected it was sound and the conclusion was still wrong.* The argument was that the
`(-60, 0, +60)` triple is derived for a *maximal* circle in the steady regime, so only a **reduced**
circle — sitting inside its station's clearance disk, where material can lie on any side — needs a
ring. Probing reduced rungs on a uniform 12-ring while maximal circles kept the triple was
implemented and measured: on 6x4 and 10x6 pockets at a 40 deg cap it changed neither the emitted path
nor the over-cap count by a single circle, at 7x the generation cost, and it was removed on that
evidence.

What that experiment could not see is that **the triple was wrong for the maximal circle too** — the
half it was trusted for. On a 20x12 pocket every position over an 80 deg cap lies between -30 and
-150 deg from the advance direction, which the triple never evaluates, so ringing only the reduced
rungs left the blind spot exactly where the engagement was. The ring is now used for **every**
circle; see
[`LOOP_PROBE_COUNT`](engagement_controlled_toolpath.md#why-a-uniform-ring-and-why-32-of-them). The
cost estimate from the old experiment held: generation on 12x8 at a 40 deg cap went 692 ms to 7784 ms.

**A flag on `engagement_controlled_toolpath`.** The two generators do not differ by a parameter: this
one emits several passes per skeleton chain and plunges and retracts once per pass. A boolean that
silently multiplied the number of passes in the returned operation stream would be a worse API than a
separate name. `engagement_controlled_toolpath` is untouched, and its emitted stream is byte-identical
before and after this work (verified by SHA-256 over the operation signature on four pocket/cap
combinations).

**Bisecting the radius ladder.** Rejected on the measurement at the top of this page, not on taste.

## Where the code lives

| Symbol | File | What it owns |
| --- | --- | --- |
| `radius_regulated_toolpath` | `src/compas_cgal/engagement_radial_toolpath.py` | The public entry point and its parameter contract |
| `_largest_admissible_radius` | same | The ladder **scan**, and the two admissibility conditions |
| `_loop_reaches_material` | same | "Does this loop still cut anything", by exact point location |
| `_radial_sweep` | same | One pass over a chain, and the FINISHED bookkeeping |
| `_machine_chain_radially` | same | The pass loop, the inter-pass `LINK`, and the budget |
| `_Regulation` | `src/compas_cgal/engagement_toolpath.py` | The validated parameter seam both generators share |
| Tests | `tests/test_engagement_radial_toolpath.py` | The non-monotonicity pin, the full-slot regression, coverage |
