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

!!! warning "Not wired in by default"

    `subtract_disk` / `subtract_annulus` / `subtract_arc_sweep` still take the global path. The local
    path is exposed as `subtract_disk_local` / `subtract_annulus_local` / `subtract_arc_sweep_local`
    and is fully tested, but nothing calls it yet. Switching the default is a decision about which
    pocket sizes matter, and it belongs to the caller.

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
are CGAL's own exact machinery; the "does this edge lie on `∂R`" test is rational equality of a
circle's centre and squared radius (both stored coefficients survive every split); the redundancy
test compares two booleans. No epsilon, no tolerance, no ulp-nudge anywhere in the path.

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

A sound extension needs a different boundary test — tracking the half-edges the insertion actually
created, which means driving `Arrangement_zone_2` with a custom inserting visitor rather than calling
`CGAL::insert`. That is real work and it is not attempted here.

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

    The paragraph above still stands as written for the **local** path: the quad region goes through
    the global `difference`, because its boundary cannot be flooded from a provably complete seed
    set. The lesson worth keeping is that "not exactly representable" bounded the wrong thing — the
    representable *under*-approximation had a much cheaper shape available all along.

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
