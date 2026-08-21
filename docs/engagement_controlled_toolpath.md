# Engagement-Controlled Trochoidal Toolpath

`compas_cgal.engagement_toolpath.engagement_controlled_toolpath` regulates the advance between
machining circles on the measured engagement instead of a dialled-in stepover, and on the pockets
measured so far it holds every machining circle except the unavoidable chain-entry loops at or below
the requested cap, while cutting **55% less travel** than the unregulated generator aimed at the same
cap. The price is generation time: **224 ms on a 12x8 pocket against 2.6 ms**, a factor of ~86 (it was
~1000x before full-turn depletion became exact and bridge depletion stopped being a disk chain, and
70.6 ms while the probe set was a 3-position triple that measurably could not see the loop's trailing
side).

The guarantee is narrow and stated exactly, because it is easy to overclaim here:

!!! warning "What is guaranteed, and what is not"

    Engagement `<= tea_cap_deg`, decided by an exact predicate, **at each evaluated tool position**.

    That is the whole claim. It is **not** a continuous guarantee between evaluated positions: the
    tool centre traverses a full circle and a bridge segment between one evaluated position and the
    next, and nothing in this path bounds what happens in between. **No certificate is produced or
    returned**, and no `MotionWitness` / `CapRefutation` object exists on this path. Coverage is not
    part of the accept/reject rule either — it follows from the advance bound and is *measured*, not
    proved.

## Why an exact predicate instead of a stepover

The stepover is a geometric proxy for engagement through the textbook radial-immersion relation
`TEA = 2*acos(1 - ae/r)`, which assumes the boundary between cut and uncut material is **straight**.
In a real trochoidal path it is the outer rim of the previous circle, so the proxy under-reads. Dialling
a 2 mm tool to `stepover = 0.5` — the value that relation says gives exactly 120 deg — produces machining
circles the audit measures at **136.4 deg**. The regulated generator, asked for the same 120 deg,
produces circles measured at **119.5 deg**. Same tool, same pocket, same cap; the difference is
measuring instead of assuming.

## How the advance is chosen

```mermaid
flowchart TD
    A["trochoidal_mat_toolpath_circular<br/>at a fine uniform pitch<br/><i>statically determined guide</i>"] --> B["ordered skeleton-chain stations<br/>centre + exact clearance radius"]
    B --> C{"bisect the integer<br/>station window"}
    C -->|"candidate station j"| D["33 evaluated positions:<br/>entry point + a uniform 32-ring<br/>phased on the advance direction"]
    D --> E["_stock_2.engagement_at<br/><b>exact cap_exceeded</b><br/><i>decided at runtime on the depleting stock</i>"]
    E -->|"any exceeded"| C
    E -->|"all pass"| C
    C -->|"largest passing index"| F["emit machining circle + bridge"]
    F --> G["deplete exact stock:<br/>subtract_arc_sweep_local + subtract_capsule_quad"]
    G --> C
```

The guide is the existing C++ generator run at a fine pitch. Its emitted `Circle` cut operations
already are the ordered skeleton-chain stations with exact clearance-derived radii, so chain
extraction, chain ordering, and the gouge-free radius derivation stay in the one implementation that
owns them. The consequence that matters for exactness: **the advance grid is an integer count of
stations**, so the bisection terminates on an integer bracket collapse — there is no float tolerance
anywhere in the accept/reject path.

### Evidence table

| Claim | Value | Source |
| --- | --- | --- |
| Advance resolution | `D/40 = r/20`, worth `<= ~7 deg` of TEA | `GUIDE_STEP_TOOL_DIAMETERS`, derivation in the constant's comment |
| Advance bound | one tool diameter (engagement saturation **and** annulus overlap) | `MAX_ADVANCE_TOOL_DIAMETERS` |
| Evaluated positions per circle | K = 33 (entry point + a uniform 32-position ring phased on the advance) | `LOOP_PROBE_COUNT`, `LOOP_PROBE_ANGLES_DEG` |
| Regulated circles, 10x6, cap 120 deg | max 119.5 deg outside entry loops | `audit_toolpath_engagement` |
| Unregulated circles, same cap via `stepover=0.5` | 136.4 deg | `audit_toolpath_engagement` |
| Motions measured above cap, 12x8 | 5 regulated vs 14 unregulated | `audit_toolpath_engagement` |
| Cut travel, 12x8 | 411.6 vs 914.3 | sum of circle circumferences and bridge lengths |
| Generation, 12x8, 2 mm tool, cap 120 deg | 223.7 ms, vs 2.6 ms | `time.perf_counter` around each generator, best of five; was 70.6 ms with the superseded 3-probe triple, 0.226 s before the quad capsule and 2.70-3.43 s before exact-annulus depletion |
| Residual stock, 10x6 | 2.35% vs 2.57%, none further than 0.28 mm from a wall | 200x120 `Stock.contains` grid |

### Why a uniform ring, and why 32 of them

The probes were `(-60, 0, +60)` deg from the advance direction until 2026-08-21, on this reasoning: in
the steady regime the material a new circle meets is a crescent at its outer rim of radial thickness
`~a*cos(phi)`, so engagement peaks at `phi = 0` and *the backward half lies inside the union of the
preceding circles' swept annuli*. **The second half of that is false**, and the cap was being broken
precisely where the triple never looked.

Measured on a 20x12 pocket with a 2 mm tool, by walking every circle the generator **accepted** at 32
tool-centre positions on the same depleting stock:

| Claim | Value | Source |
| --- | --- | --- |
| Over-cap positions in the advance-facing half, 80 deg cap | **0** at `0, +30, +60, +90, +120` deg | replay at 32 positions/loop, entry cuts excluded |
| Over-cap positions in the trailing half, same run | **all 120**, between `-30` and `-150` deg | same |
| Same, 100 deg cap | **all 24**, between `-90` and `-150` deg | same |
| Worst accepted circle, 80 deg cap | **101.7 deg**, peaking at `-135` deg from the advance | same |
| What the superseded triple read on that same circle | **39.1 deg** | its three probe positions re-measured against the identical stock |
| The forward peak is still real | **54 of 90** circles peak at `phi = 0` | peak-angle histogram |

The last two rows are the point: the forward crescent is real, it is simply not where the cap breaks,
and the gap was probe **placement** — not a difference between the generation and audit depletion
models, which read the same stock the same way.

**Why uniform and never one-sided.** The loaded quadrant is the *trailing-lateral* one, on the side
fixed by the loop's turn direction, and it mirrors exactly with the milling direction: the same pocket
and cap that puts `4/16/52/36/12` over-cap positions in the `(-30,-60,-90,-120,-150)` bins under climb
milling puts `12/36/52/16/4` in the `(+150,+120,+90,+60,+30)` bins under conventional. Any placement
biased to one side is tuned to one winding and blind on the other. Measured at equal count on climb —
the winding it would be tuned for — a trailing-biased set still **loses** to uniform: at 12 probes it
reads 128.8 deg where uniform reads 87.3 (80 deg cap), because refusing different candidates changes
the path and the peak migrates. A forward-biased set at 24 probes reproduces the superseded triple's
path **byte for byte**: all 24 of its probes sit in the region that is never over the cap.

**Why 32 — measured convergence, not a derivation.** There is a geometric floor: consecutive evaluated
positions are a chord `2*R*sin(pi/K)` apart, so their tool disks overlap at all only while that chord
stays under `2*r`, which on the largest loops here (`R = 4.998`, `r = 1.0`) needs `K >= 16`. The
measurement says the floor is not enough. Worst engagement an **independent** 60-position walk
(phase-offset by half a step, so it shares no grid with the probes) finds on accepted circles:

| K | 3 (old) | 8 | 12 | 16 | 24 | **32** | 40 | 48 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cap 40 deg | 86.7 | 71.7 | 55.3 | 57.8 | 49.3 | **43.4** | 43.4 | 43.7 |
| cap 80 deg | 101.7 | 99.2 | 87.3 | 88.6 | 82.8 | **81.7** | 82.2 | 82.4 |

The last count at which either row moves is 32; 40 and 48 buy nothing and cost 24% more. Neither row is
monotone in `K` — a denser probe set refuses different candidates, which changes the whole path — so
this is a convergence measurement, not a trend to extrapolate.

**The entry point** is evaluated in addition to the ring: it is the terminus of the bridge cut that
precedes the circle, and it is read *before* the bridge is removed, so that probe measures what the
linking cut itself runs into. That makes `K = 33` evaluated positions per candidate.

!!! warning "A denser ring narrows the sampling gap; it does not close it"

    Nothing here bounds engagement between two evaluated positions, and the residue is measured rather
    than hypothetical: at a 40 deg cap on the 20x12 pocket the independent 60-position walk still finds
    **43.4 deg** on a circle every probe accepted. The soundness of each individual verdict is
    unchanged and exact; what changed is how many places the question is asked.

!!! note "Under verification: the ring versus the audit"

    `tests/test_engagement_toolpath.py::test_only_chain_entry_loops_are_measured_above_the_cap` pins
    the gap between the generator's ring and the audit's twenty stations per circle: on a 6x4 pocket
    the only motions the audit finds above the cap are the chain-entry loops the generator already
    refuses and reports. That is a **measurement on one pocket**, not a proof that 32 probes suffice in
    general. `test_a_trailing_quadrant_load_is_refused_where_the_superseded_triple_accepted_it` pins
    the finding itself, on a hand-built stock state rather than on a generated path.

## Depleting a full turn is exact, not sampled

Every machining circle is a full turn, and the region a tool of radius `r` sweeps about a circular
guide of radius `rho` is **exactly** the annulus between `rho - r` and `rho + r`. That region is
representable in `Gps_circle_segment_traits_2`: both radii are doubles, hence rationals, so their
squares are the rational squared radii the traits class needs. Two boundary circles — four
x-monotone arcs — describe it completely.

The generator previously under-approximated that sweep with a chain of tool disks spaced
`2*r*sqrt(CHAIN_SLACK_FRACTION)` apart, which is **1260 disks** for one turn at `r = 0.5`. The chain
was deliberately conservative: it under-covers, so the model retained more material than the tool had
actually removed and every later engagement query read high. Moving to the annulus is a
**correctness improvement, not a loosened tolerance** — the modelled region becomes the true region,
so there is no longer an approximation to compensate for.

`Stock2::subtract_annulus` is the new path, `subtract_arc_sweep` routes the full-turn case to it, and
the partial-arc and capsule chains are untouched. The two radii are formed in exact rational
arithmetic rather than in doubles: `subtract_annulus_exact` shares its region builder with
`exact_full_circle_sweep_oracle`, so the region the fast path removes and the region the depletion
certificates are proved to under-cover are the same construction by identity.

| Claim | Value | Source |
| --- | --- | --- |
| Arrangement vertices after one full turn | **8** via the annulus vs **2536** via the chain — 317x | `Stock.arrangement_stats`, `test_annulus_arrangement_is_orders_smaller_than_the_chain` |
| Depletion cost, 12x8, 2 mm tool | **2.5 ms/call** vs 57.4 ms/call, 49 calls | `time.perf_counter` around `Stock.subtract_*`, A/B on one machine |
| `engagement_at` cost, same run | **0.160 ms/call** vs 1.288 ms/call, 427 calls | same harness |
| Total generation, same run | **0.226 s** vs 3.411 s — 15x | same harness |
| Chain shortfall it removes | `r * CHAIN_SLACK_FRACTION` = 5.0e-5 mm at `r = 0.5` | sagitta `s^2/(4r)`, `test_annulus_removes_the_slivers_the_disk_chain_retains` |
| Toolpath produced | **unchanged**: 49 cuts / 441.8 length / 5 over cap (12x8), 112 / 595.3 / 7 (L-shape) | replay against `_stock_2.engagement_at`, 12 positions per loop |

!!! note "The disk chain is not a subset of the annulus, and that is the chain's doing"

    The chain's centres are `cx + guide_r*cos(a)` evaluated in double arithmetic, so they land up to
    **7.0e-16 mm off** the guide circle and the chain removes a sub-femtometre crescent *outside* the
    true swept region. Exact set inclusion of the annulus result in the chain result is therefore
    false — for a reason that belongs to the chain's sampling, not to the annulus. The inclusion
    **does** hold exactly against `subtract_exact_full_circle`, whose centres are exact rational
    points on the guide circle; that is the form the test asserts. The chain's meaningful error is
    its 5.0e-5 mm shortfall, eleven orders of magnitude larger.

The guide radius itself stays a double surrogate. `sqrt(rx^2 + ry^2)` is irrational in general, so
`(rho +- r)^2` is not a rational squared radius and the exact annulus about the *true* guide circle is
not representable at all. That surrogate is a pre-existing property of this double-valued API — the
chain samples the very same approximate circle — and it is injected exactly, with no snapping and no
correction constant.

## Depleting a bridge is six curves, not a chain

The bridge between two machining circles sweeps a **capsule**, and unlike the full turn's annulus the
capsule really is not representable in `Gps_circle_segment_traits_2`: its side lines stand off the
segment by `r/sqrt(dx^2 + dy^2) * (dy, -dx)`, and that offset is irrational while the traits class
holds lines as rational triples. So the bridge depletion stays an **under**-approximation. What
changed is the *shape* of that under-approximation:

```
region  =  disk(A, r)  U  disk(B, r)  U  rect(A, B, h)
```

The two end disks are exact and cover the capsule's semicircular caps exactly. The rectangle runs
along the segment at half-width `h`, a **rational** slightly under `r`, so it sits strictly inside the
capsule's straight part. Six curves per bridge — two circles and four segments — where the chain
needed `len / (2*r*sqrt(f))` disks, which is **601 disks** for a 6 mm bridge at `r = 0.5`.

The side lines never had to be *met*, only **under-cut**, and the chain was already paying
`CHAIN_SLACK_FRACTION * r` for exactly that. The quad spends the same budget in a better shape.

!!! note "The certificate is two exact comparisons, not a rounding argument"

    `h` is *chosen* in doubles — it is a construction parameter, the quad twin of the chain's
    spacing, aimed at the middle of the admissible band. It is then **checked**, exactly:
    `(1 - f)*r <= h <= r` as a comparison of squared rational lengths, plus an exact
    `CGAL::orientation` on the four corners. A failure throws `CapsuleQuadCertificateError` rather
    than removing a region the certificate does not cover. Nothing here compares a `to_double`, and
    no square root is taken.

    The band leaves `f/2 = 5e-5` of relative headroom on each side against a double rounding of
    `~2e-16` — eleven orders of magnitude — so the checks are structural guards, not a filter the
    construction is expected to trip.

Both halves of the contract are exact consequences of that band. A point of the rectangle is
`A + t*d + s*h_vec` with `|s| <= 1` and `h_vec` perpendicular to `d`, so its distance to the segment
is exactly `|s|*h <= r`: **nothing outside the true swept capsule is ever removed**. And anything
within `(1 - f)*r` of the segment is either inside a full-radius end disk or has its perpendicular
foot in the segment's interior at offset `<= h`: **the under-coverage stays inside the documented
budget**.

| Claim | Value | Source |
| --- | --- | --- |
| Arrangement vertices after one capsule | **10** via the quad vs 1206 via the chain at `r = 0.5` (120x); 10 vs 306 at `r = 2.0` (30x) | `test_capsule_quad_arrangement_stays_orders_below_the_chain` |
| Growth per bridge vs bridge length | **constant** (10 vertices for a 0.5 mm and an 8 mm bridge) | `test_capsule_quad_growth_is_independent_of_segment_length` |
| Bridge depletion, 12x8, 2 mm tool | **0.337 ms/call** vs 4.389 ms/call, 22 calls | wrapped `Stock.subtract_*`, A/B on one machine |
| `engagement_at`, same run | **0.061 ms/call** vs 0.158 ms/call, 427 calls | same harness — every later query walks a smaller arrangement |
| Full-turn depletion, same run | **0.356 ms/call** vs 1.078 ms/call, 27 calls | same harness, same reason |
| Total generation, 12x8 | **72.4 ms** vs 222.5 ms — 3.1x | same harness |
| Total generation, L shape | **173.7 ms** vs 743.4 ms — 4.3x | same harness |
| Toolpath produced | **bit-identical** operation stream on both pockets: 49 cuts / 441.8 / 5 over cap (12x8), 112 / 595.3 / 7 (L-shape) | float-exact signature comparison; exceedances by 12-position `engagement_at` replay, counted under *both* replay depletions |
| Under-coverage | **unchanged**: `r * CHAIN_SLACK_FRACTION`, 1e-4 mm at `r = 1.0` | the certified band, `test_capsule_quad_covers_the_documented_slack_band` |

!!! warning "Neither approximation contains the other, and the quad is not the exact capsule"

    The chain covers the full radius `r` at each disk centre and dips to `r*sqrt(1 - f)` midway
    between centres; the quad covers a flat `h` along the whole rectangle and the full `r` only at
    the caps. So the two removed regions are **incomparable as sets** — asserting inclusion either
    way would be false. What both satisfy, and what the tests assert, is the same pair of exact
    bounds: contained in the true capsule, containing the `(1 - f)*r` band.

    `subtract_capsule` is **not** removed. It stays the reference the quad is measured against, and
    `tests/test_stock_capsule_quad.py` decides both paths' properties against the same exact
    rational oracle.

!!! note "No local variant: the flood cannot start where it would have to"

    `subtract_disk_local` / `subtract_annulus_local` (see [Local Depletion](local_depletion.md)) do
    not extend to this region, for a reason distinct from the disk chain's. Recognising a *linear*
    boundary edge is mechanically possible — `_X_monotone_circle_segment_2::supporting_line()` gives
    the same rational-identity test that `supporting_circle()` gives for arcs. The blocker is the
    **seeds**: the local update inserts whole supporting curves and blocks the flood on every
    sub-curve of them, so the two end-cap lines and the rectangle's sides partition the region's
    interior into cells whose number and shape depend on the configuration (how the disks overlap,
    whether a corner falls inside the other disk). A seed set that is provably complete for every
    configuration is real work, and a missing seed leaves material silently marked present. Not
    attempted; the global path costs 0.337 ms/call, which is no longer the dominant term.

## Chain entries are over the cap by construction

The first machining circle of a chain meets virgin stock: it is a full slot whatever the advance, and
the audit measures it at a full turn. There is nothing to search — the circle is emitted, counted, and
reported through `UnavoidableEngagementWarning`, never emitted silently. On both rectangles measured,
the over-cap motions are **exactly** the chain-entry loops and nothing else.

A helical or ramped entry would remove them. That is not implemented here.

## Prior art

| System | Relationship | Evidence |
| --- | --- | --- |
| Held's trochoidal engagement control (float bisection to `eps = 1e-3`) | **stronger per evaluated position** — each verdict is an exact predicate on the exact arrangement, no tolerance in the accept/reject path | `_stock_2.engagement_at` returns `cap_exceeded` decided exactly; the surrogate crossing is `4*sin^2(theta/2)` at one declared boundary |
| Continuous partition of each motion (`compas_cgal.adaptive`) | **weaker** — a partition bounds every centre on the motion; this bounds only the positions it evaluated | this module deliberately does not import `compas_cgal.adaptive`, and returns no witness type |
| `trochoidal_mat_toolpath_circular` (stepover proxy) | **stronger on measured engagement, weaker on speed** | 119.5 deg vs 136.4 deg at a 120 deg target; 72.4 ms vs 2.6 ms on 12x8 |

## Known limitations

- **Between-position engagement is unbounded.** Stated above; this is the defining limitation.
- **Monotonicity is assumed, not proved.** The bisection assumes engagement is monotone
  non-decreasing in advance distance. Under that assumption it returns the largest admissible
  station; without it, it still returns an admissible one — every accepted station's evaluated
  positions passed the exact predicate regardless.
- **Generation is ~27x slower** than the unregulated generator (was ~78x before the quad capsule and
  ~1000x before full-turn depletion became exact). `engagement_at` is now the largest single term
  (34% of a 12x8 run) and depletion is 22%.
- **Guide tangency is G1 on straight guide segments only**, approximate through turns — the same
  model the existing generator uses.
- **Coverage is measured, not certified.** The advance bound keeps consecutive annuli overlapping;
  the consequence was checked on a `Stock.contains` grid.
