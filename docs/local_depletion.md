# Local Depletion

Removing a disk or an annulus from the exact stock can be done **locally** — by editing the boolean
engine's own arrangement around the removed region instead of overlaying the whole stock — and the
result is bit-identical: the same exact point set *and* the same canonical `General_polygon_set_2`
representation, verified by exact equality after each of 26 mixed depletions on three stock shapes,
in both orders, plus 40 pseudo-random removals.

**It is worth switching on for cut-up stocks and not for fresh ones.** Per full-turn annulus the local
path costs **1.77 ms against the global 4.18 ms on the L-shaped pocket (2.4x faster)** but **1.04 ms
against 0.79 ms on the 12x8 rectangle (1.3x slower)**. The reason is structural, not a tuning
accident: the global difference rebuilds the entire arrangement, so it scales with the whole stock;
the local update walks the zone of four circular arcs, so it scales with the removed neighbourhood
plus the complexity of the faces those arcs cross. The rectangle's arrangement is 360 vertices, the
L-shape's is 1078, and that is already enough for the local path to win.

!!! note "Who calls it (2026-10-14)"

    `subtract_disk` / `subtract_annulus` / `subtract_arc_sweep` still take the global path. The
    generators deplete loops through `subtract_arc_sweep_local`, and since Task 0.2 the **replay**
    (`compas_cgal.engagement._subtract_operation`) removes every straight cut with
    `subtract_line_sweep_local` and every machining circle with `subtract_circle_sweep_local`. The
    numbers in this section predate both the exact line sweep and the shared-locator change; see
    [Line sweeps: the sub-curve boundary test](#line-sweeps-the-sub-curve-boundary-test) for the
    current ones.

## Evidence

Real depletion sequences captured from `engagement_controlled_toolpath` (tool 2.0, cap 120 deg) and
replayed through both paths. Every row is measured, not modelled.

| Pocket | Arrangement | Path | Full-turn annulus (27 / 60 calls) | Disk chain (22 / 52 calls) | Depletion total |
|---|---|---|---|---|---|
| rect 12x8 | 360 v / 720 he | global | 21.4 ms — **0.792 ms/call** | 88.8 ms — 4.038 ms/call | 110.2 ms |
| rect 12x8 | 360 v / 720 he | local | 28.1 ms — **1.042 ms/call** | 92.6 ms — 4.211 ms/call | 120.8 ms |
| L-shape | 1078 v / 2156 he | global | 250.7 ms — **4.178 ms/call** | 281.7 ms — 5.417 ms/call | 532.4 ms |
| L-shape | 1078 v / 2156 he | local | 106.1 ms — **1.769 ms/call** | 287.6 ms — 5.531 ms/call | 393.7 ms |

The disk-chain column is the *same code* on both rows — partial arcs and capsules are out of the
local path's reach (see [Why the disk chain is out of reach](#why-the-disk-chain-is-out-of-reach)) —
so its variation is run-to-run noise and calibrates how much of the annulus difference is real.

End to end, in the generator (median of three fresh processes each):

| Pocket | Path | Total generation | `engagement_at` | stock subtract |
|---|---|---|---|---|
| rect 12x8 | global | 0.213 s | 427 calls / 0.065 s | 49 calls / 0.113 s — 2.297 ms/call |
| rect 12x8 | local | 0.218 s | 427 calls / 0.066 s | 49 calls / 0.119 s — 2.431 ms/call |
| L-shape | global | 0.879 s | 812 calls / 0.286 s | 112 calls / 0.534 s — 4.766 ms/call |
| L-shape | local | **0.724 s** | 812 calls / 0.281 s | 112 calls / **0.386 s — 3.444 ms/call** |

Generation output is unchanged on both paths, as it must be for a bit-identical stock: rect 12x8 =
49 cuts / polyline 441.8 / 5 operations truly exceeding the cap; L-shape = 112 cuts / polyline 595.3
/ 7 truly exceeding. ("Truly exceeding" = `_stock_2.engagement_at` reports `cap_exceeded` at one of
12 positions around the loop — *not* `EngagementReport.cap_violations`, which counts operations that
could not be **certified** and is a different quantity.)

That measurement is `benchmarks.exceedance.count_truly_exceeding`, and the corpus now records it
beside the certifier's count as two separate columns: `truly_exceeding`, a sampled **lower bound** on
how many motions are over the cap, and `uncertified`, the sound **upper bound** of motions whose cap
could not be proved. Reporting either one under the other's name is how a generator that had improved
2.4x came to read as a regression.

## Where the time goes

Phase timing over the same replayed sequences, from a temporary instrumented build:

| Pocket | insert (zone walk) | flood | redundant-edge removal | point location |
|---|---|---|---|---|
| rect 12x8 | 25.9 ms (92%) | 0.7 ms | 0.4 ms | 1.1 ms (4%) |
| L-shape | 99.6 ms (92%) | 1.7 ms | 0.9 ms | 6.5 ms (6%) |

The marking and the repair are free. Essentially all of the local path is `CGAL::insert`'s zone walk,
and that cost is intrinsic: `Arrangement_zone_2` finds where a curve leaves the current face by
scanning that face's CCBs, and a depleting stock is one enormous face whose inner CCB holds hundreds
of half-edges. Point location — the obvious suspect, and the thing a landmark or trapezoidal
structure would fix — is 4-6% and not worth attacking.

That is the honest ceiling of this approach. It removes the *rebuild*, not the *sweep*.

## The soundness argument

Write `S` for the stock point set, `R` for the removal region, and `A` for the arrangement the
`General_polygon_set_2` holds; its faces carry `contained()` = "in `S`".

**I1 — insertion preserves the containment labelling.** After inserting the x-monotone arcs of `∂R`,
every face of the refined arrangement still satisfies `f ⊆ S ⟺ f.contained()`. Insertion never moves
a point across a boundary; it only *splits* faces, and both pieces of a split are subsets of one
original face.

!!! warning "CGAL does not copy the face record on a split"

    `Arrangement_on_surface_2_impl.h:2927` creates the split face with a bare `_dcel().new_face()` —
    no `assign` — so `contained()` reads **false** on every piece the region boundary cuts off, and
    the resulting stock looks plausible while being silently wrong. `ContainmentPropagator` restores
    the invariant through the `after_split_face` notification. On a bounded planar topology that is
    the *only* face-creating event insertion can raise: the remaining `new_face()` sites in vendored
    CGAL 6.0.1 are DCEL initialisation and the unbounded-planar topology traits.

**I2 — no face straddles the region.** `∂R` is now part of the arrangement, so every face lies
entirely inside `R` or entirely outside it.

**I3 — the flood marks exactly the inside.** `int(R)` is open and `∂R ∩ int(R) = ∅`, so a path
between two faces inside `R` crosses only non-boundary edges. Flooding from a seed in each connected
component of `int(R)`, never crossing an edge that lies on `∂R`, reaches every face inside `R` and no
face outside it. Setting `contained(false)` on exactly those faces turns the characteristic function
of `S` into that of `S \ R`. A bounded region encloses no unbounded face, so reaching one means the
flood escaped — that raises `LocalDepletionEscapedError` rather than depleting the wrong material.

**I4 — the repair is local too.** An edge becomes redundant (its two incident faces agree on
`contained()`) only if the marking changed one of their flags, i.e. only if it bounds a flooded face;
every other face kept its flag and the Gps invariant held before. So the candidate set is exactly the
boundary of the flooded faces. Containment values never change during removal — a merge joins two
faces that already agreed — so the predicate is stable across the loop.

Every decision above is an exact predicate on exact quantities. Point location and curve insertion
are CGAL's own exact machinery; the "does this edge lie on `∂R`" test is exact **sub-curve
membership** in one of the inserted boundary curves (see
[the sub-curve boundary test](#line-sweeps-the-sub-curve-boundary-test)); the redundancy test
compares two booleans. No epsilon, no tolerance, no ulp-nudge anywhere in the path.

### The orientation trap

`General_polygon_set_2`'s representation invariant requires every stored curve to carry the
**contained side on its left**, and `CGAL::insert` stores the curve exactly as handed over. A removal
region's interior is by construction *not* contained, so its outer boundary circle must be traversed
**clockwise** — the mirror of `disk_polygon`, which builds the counterclockwise circle of a region
being *added*. An annulus's inner circle bounds material the removal keeps, so it stays
counterclockwise.

Getting this wrong produces a stock that is **exactly the right point set**, stored in an arrangement
of **exactly the right size**, that nonetheless fails `Gps::is_valid()` — and would then degrade
every later zone query by fragmenting arcs the engagement harvest expects to find whole. It was
caught only because the equivalence test asserts the representation invariant alongside point-set
equality. `exactly_equals` alone is not a sufficient gate.

## Why the disk chain is out of reach

The chain paths — `subtract_capsule` and partial `subtract_arc_sweep` — are the *larger* depletion
cost on the 12x8 pocket (88.8 ms of 110.2 ms), so it is worth stating precisely why this machinery
does not extend to them.

The flood needs to recognise an arrangement edge as lying on `∂R`. For concentric circles that is
exact circle identity, and it is sound because a disk's or annulus's boundary circle is *entirely*
boundary. For a union of `N` overlapping disks it is not: a single chain circle contributes arcs that
are on `∂(∪D)` **and** arcs that run through the union's interior. Blocking the flood at an interior
arc would leave material on its far side marked as still present — a stock that looks fine and
silently corrupts every later engagement verdict, which is exactly the failure mode this whole gate
exists to prevent.

A sound extension needs a different boundary test. For a disk *chain* it is still not attempted:
the chain's union boundary is hundreds of arcs and nothing about it is worth keeping (the replay no
longer uses it). For a single line sweep it is done — the region is built first, its **true**
boundary is read off a scratch arrangement, and the flood stops at sub-curves of that boundary, not
at whole supporting circles. See
[Line sweeps: the sub-curve boundary test](#line-sweeps-the-sub-curve-boundary-test).

It would also only reach part of the prize. Splitting the chain cost further:

| Pocket | disk-chain union construction | global difference | disks built |
|---|---|---|---|
| rect 12x8 | **69.0 ms (73%)** | 22.9 ms | 1269 |
| L-shape | 79.9 ms (27%) | **209.1 ms (72%)** | 1350 |

On the rectangle the dominant cost is *building the ~58-disk union per capsule*, which no
**local**-update change can touch. An oriented capsule's side lines are irrational and therefore not
representable in `Gps_circle_segment_traits_2`, so unlike the full turn there is no exact swept
region to reach for.

!!! success "Superseded on the cost, not on the reasoning — `subtract_capsule_quad`"

    The chain's cost went away without an exact capsule and without touching this machinery: the
    side lines never had to be *met*, only **under-cut**, which two exact end disks plus a rational
    rectangle at half-width `h` with `(1 - f)*r <= h <= r` do in six curves. Bridge depletion on the
    12x8 pocket fell from **4.389 ms/call to 0.337 ms/call** and the emitted toolpath is
    bit-identical. See
    [Depleting a bridge is six curves, not a chain](engagement_controlled_toolpath.md#depleting-a-bridge-is-six-curves-not-a-chain).

    The lesson worth keeping is that "not exactly representable" bounded the wrong thing — the
    representable *under*-approximation had a much cheaper shape available all along. And for an
    axis-aligned or Pythagorean direction it *is* exactly representable: the replay's line sweep
    builds that capsule as an equality.

## Line sweeps: the sub-curve boundary test

The replay's one line-sweep model (`src/stock_line_sweep_2.h`, Task 0.2) is
`disk(A, r) ∪ rect(A, B, h) ∪ disk(B, r)`, built exactly in the circle-segment traits:

- **exact capsule** when `dx² + dy²` is a square in ℚ (decided by `Fraction_traits` + `Is_square`,
  [exactness](exactness.md#case-is-a-segments-length-rational-ask-the-number-type)): half-width
  vector `r·n/L`, rational, so the region *is* the swept capsule;
- otherwise the **canonical-direction quad**: the same two **exact, full-radius** end disks and a
  rectangle of half-width vector `s·perp(u)`, `u = d / max(|dx|, |dy|)`, with `(1 − 1e-4) r ≤ |h| ≤ r`
  checked exactly. Only the side band `h < d < r` is inexact.

The end disks stay exact because the tool's endpoint positions are (Ruling 0.2-endcaps): a plunge,
or the next motion's start, at either end meets nothing the model left behind. An earlier draft used
a *stadium* (caps of radius `|h|`) so that collinear cuts would nest exactly; it left a ring up to
`1e-4 r` wide around every bridge end, and the next loop's entry station read that ring as 41-70° of
engagement (`build/task02/cap_attrib.py`: +8 and +44 cap exceedances on two gate cases at 40°). With
full-radius end disks, a collinear oblique re-cut may change the stock, but only inside the earlier
cut's certified OVER region (its end disks poke into the band). That is tested with exact set
operations as `(before − OVER) == (after − OVER)`. Axis and Pythagorean re-cuts are exact no-ops. The
side band is the documented model limit: material up to `1e-4 r` wide may survive alongside an
oblique cut, and later motions may read it. Pythagorean authoring removes it.

**The boundary test.** The end-disk circles carry arcs that lie *inside* the region (under the
rectangle) as well as arcs on its boundary, and its sides are segments, so circle identity is
unsound for it. The region boundary is therefore described by its **curves**, and an arrangement
edge `e` lies on `∂R` iff it is a *sub-curve* of some boundary curve `β`:

1. same supporting curve — exact rational comparison of the stored coefficients
   (`has_same_supporting_curve`, scale-free for lines);
2. for arcs, the same (upper/lower) half of the circle — on one half an x-monotone arc is fixed by
   its endpoints;
3. both endpoints of `e` on closed `β` — `is_in_x_range` and `point_position == EQUAL`, which also
   bounds a vertical segment's y-range.

It is complete because inserting `β` makes its endpoints arrangement vertices: an edge on `β`'s
support either lies inside `β` or misses its relative interior. Segments are convex, so two endpoints
on `β` put all of `e` on it. Whole circles reduce to the old rule (their `β` are two full half
circles), so disk and annulus removal keep their behaviour; `LocalRemovalRegion` lost its
`boundary_circles` field.

**The region.** `LineSweepRegion2::local_under_region` builds `R` in a scratch set (at most six input
curves), checks it is one polygon with no hole (`LineSweepRegionTopologyError` otherwise), and takes
every scratch edge along the halfedge whose face is the *exterior* — the `MATERIAL_OUTSIDE` rule of
the orientation trap — so the curves are already split at their exact, possibly one-root, cap/side
junctions. The seed is the exact midpoint of `A` and `B`, at distance `0 < |h|` from the axis.

**The locator.** One trapezoid-RIC map now serves both insertion and the flood, and it is destroyed
before the repair phase. Before, the path built two maps per call and kept the second alive while
removing redundant edges — a trapezoidal map is an arrangement observer, so every removal re-built
it locally. Phase timing (temporary instrumentation, 120 deg L-shape generation plus two replays, 229
local removals): the removal phase took 138 ms with the map alive and about 1 ms with no locator
attached (measured under the walk locator below); end to end the local replay of the largest gate
case fell from 11.0 s to 5.2 s. The walk-along-line locator
measured faster still (largest gate replay 5.2 s → 1.3 s) and was **rejected**: its progress guard
is a debug assertion, and on the multi-outer-CCB fixture
(`tests/test_stock.py::test_local_depletion_terminates_on_multi_outer_ccb_stock`) it did not
terminate — exactly the [exactness](exactness.md) rule on locators. Frozen-generator output is
byte-identical (all 12 gate cases, `snapshot_toolpath` + polyline digest), and generation is faster
(40 deg `engagement_controlled` rect 20x12: 22.7 s → 16.8 s; `radius_regulated` rect 20x12:
55.5 s → 44.3 s).

**Evidence.** `tests/test_stock_local_line_sweep.py`: Hypothesis sequences of up to six line sweeps
mixing axis, Pythagorean, generic oblique, exactly-collinear dyadic and zero-length segments (plus a
repeat of the first), circle-sweep sequences with re-sweeps, and mixed line/circle sequences. After
**every** removal: `exactly_equals` the global twin, `representation_is_valid()`, and equal
`arrangement_stats()`. With the old circle-identity test compiled in, 4 of the 15 tests across the two
local-depletion test files fail, so the gate discriminates. A generic-double witness (two oblique sweeps through
non-dyadic endpoints) stays under 1 s.

End to end, replaying each gate path through both paths (no clones; depletion time only):

| Case | final stock | lines global / local | loops global / local |
|---|---|---|---|
| 120 deg `engagement_controlled` rect 12x8 | 67 v | 12.5 ms / 37.6 ms | 9.5 ms / 57.9 ms |
| 120 deg `engagement_controlled` L-shape | 92 v | 20.9 ms / 67.5 ms | 17.7 ms / 80.0 ms |
| 40 deg `radius_regulated` L-shape (1643 ops) | 324 v | 705 ms / 2354 ms | 582 ms / 2776 ms |

**On the replay's stocks the global overlay is 3-5x cheaper than the local path.** The exact sweeps
keep the stock tiny (hundreds of vertices, not the chain's thousands), and at that size the
trapezoidal map's `O(n log n)` exact build per call dominates the local path, while the overlay is
cheap. The local path is what the replay uses (Task 0.2 ruling); the measurement says the global
twin, which is bit-identical, would be the cheaper choice. The exact four-run probe (four 16 mm wall
runs on 20x12): **24 vertices** (the chain built 8,014).

## Rejected: pick the path by arrangement size

The crossover is real and could be automated with a vertex-count threshold. It is not, because the
threshold is a machine- and geometry-dependent constant with no derivation behind it — the kind of
measured magic number this repository rejects. The two paths are exposed by name; the caller chooses.

## Testing

`tests/test_stock_local_depletion.py`. The gate is exact point-set **equality** against the global
path — never containment, never a subset relation — asserted after *every* removal so a divergence
names the operation that introduced it, plus `Gps::is_valid()` and equality of the arrangement's
vertex/half-edge/face counts.

| Test | What it covers |
|---|---|
| `..._step_by_step` | 26 mixed depletions x 3 stock shapes: disjoint, overlapping, repeated-identical, boundary-straddling, wholly-outside, degenerate annuli, nested annuli, coincident boundary circles |
| `..._in_reverse_order` | the same removals arriving back to front, so intermediate arrangements the forward run never builds |
| `..._on_a_random_sequence` | 40 pseudo-random disks and annuli, fixed seed, reaching tangency and containment configurations a hand-written list does not enumerate |
| `test_local_arc_sweep_matches_global_full_turn` | the entry point the generator actually calls, including a guide narrower than the tool and a sweep that exactly retraces an earlier annulus |
| `..._defers_partial_arcs_to_the_global_chain` | the documented out-of-scope fallback |
| `..._keeps_the_arrangement_the_same_size...` | representation growth, which `exactly_equals` cannot see |
| `test_local_line_sweep_equals_global` (`test_stock_local_line_sweep.py`) | Hypothesis line-sweep sequences: axis, Pythagorean, oblique, collinear sub-segments, zero length, repeats |
| `test_local_circle_sweep_equals_global` | Hypothesis full-circle sweeps from the authored radius, with re-sweeps |
| `test_mixed_local_sweeps_equal_global` | interleaved line and circle sweeps, so capsule boundaries meet annulus boundaries |
| `test_four_wall_runs_locally_equal_the_global_runs` | the four tangent-touching 16 mm runs of the acceptance probe |
