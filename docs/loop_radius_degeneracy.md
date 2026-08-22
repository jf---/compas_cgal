# The Degeneracy Boundary of a Trochoidal Loop

A trochoidal loop stops being a trochoid at `rho <= r`, and on a middle-curve
construction that boundary is reached at a clearance of **three tool radii**:
`d <= 3r` makes the loop degenerate *at its largest admissible radius*, so no
choice of guide radius rescues it. Every pocket point narrower than three tool
radii is therefore structurally unable to carry a trochoid, whatever the
generator does. Two further consequences fall out of the same arithmetic, and
both are measured below: a chain entering virgin stock is a full-immersion cut
at *every* non-degenerate radius, and the bridge between two loops tilts off the
guide tangent by the arctangent of the clearance slope, independently of how far
the generator advances.

## The identity

A tool of radius `r` running a circle of radius `rho` sweeps the radii
`[rho - r, rho + r]` about the loop centre. The swept region has an uncut core
exactly when `rho > r`; at `rho <= r` the tool passes over its own loop centre
and the sweep is a filled disk. That boundary is
`benchmarks.quality.DEGENERATE_LOOP_RATIO`, and it is physics rather than a
tuned constant.

On the middle-curve construction the guide radius is capped by the clearance at
the medial-axis point. From `src/compas_cgal/adaptive/candidates.py` — the same
expression at both call sites, `_mathsm_geometry` (line 134) and
`enumerate_middle_curve_candidates` (line 1115):

```python
maximum = Fraction.from_float((distance - tool_radius.value) / 2.0)
```

so `rho_max = (d - r) / 2`, where `d` is the distance from the medial-axis point
to its generator site. Combining the two:

```
rho_max <= r   <=>   (d - r) / 2 <= r   <=>   d <= 3r
```

**A middle-curve station whose clearance is at most three tool radii cannot carry
a non-degenerate loop.** The bound is on the *cap*, so it is not a statement
about a particular choice of `rho` — it says the whole admissible interval
`(0, rho_max]` lies at or below the degeneracy boundary.

!!! note "The MAT-centred construction has a different threshold"

    Placing the loop centre *on* the medial axis instead of on the middle curve
    gives `rho_max = d - r - c` for a radial clearance `c`, degenerate when
    `d <= 2r + c`. The middle curve is degenerate over a strictly wider band
    because its radius is roughly half the MAT-centred one — the centre sits
    midway between the `r`-offset boundary and the axis, so the loop's far rim
    only just reaches the axis.

### Worked case: an L with arms six units wide

`Polygon([[0,0,0],[12,0,0],[12,6,0],[6,6,0],[6,10,0],[0,10,0]])` with a 2 mm
tool. Along a six-wide arm the clearance on the spine is `d = 3.0` and the tool
radius is `r = 1.0`, so

```
rho_max = (3.0 - 1.0) / 2 = 1.0 = r
```

exactly. The entire spine of both arms sits precisely *on* the degeneracy
threshold: every middle-curve loop there is a disk sweep wearing a circle's
name, and the pocket is at the exact width `d = 3r` at which the identity
predicts it. An arm has to exceed three tool radii of half-width — six tool
*diameters* of full width — before a middle-curve trochoid has a hole in it.

## Corollary 1: a chain entry is a full-immersion cut

The degeneracy boundary is also where entry engagement saturates. Measured on
`rect_12x8` with a 2 mm tool, plunging at the loop's entry point and probing the
32-position ring on the resulting stock:

| `rho` | `rho / r` | peak engagement | degenerate |
|------:|----------:|----------------:|:----------:|
| 0.100 | 0.100 | 191.48° | yes |
| 0.250 | 0.250 | 208.96° | yes |
| 0.500 | 0.500 | 240.00° | yes |
| 0.750 | 0.750 | 277.18° | yes |
| 0.900 | 0.900 | 308.32° | yes |
| 0.990 | 0.990 | 343.78° | yes |
| **1.000** | **1.000** | **360.00°** | **yes** |
| 1.100 | 1.100 | 360.00° | no |
| 1.500 | 1.500 | 360.00° | no |
| 2.998 | 2.998 | 360.00° | no |

The crossover is exact and it is at `rho = r`. A loop only escapes full
immersion while its far point stays within `2r` of the plunge hole, which is
`2 * rho < 2r`, which is degeneracy.

**Consequence for any quality gate:** `degenerate_loops = 0` and a bound on the
engagement *step* between adjacent motions are mutually unsatisfiable at a chain
entry. A non-degenerate entry loop peaks at 360°; the motion after it respects
the cap; so the step is at least `360 - cap` — 240° at a 120° cap. Gating both
at once asks a generator to be two contradictory things.

**The contradiction is in the ENTRY STRATEGY, not in the criterion.** It is
tempting to exclude the entry motion from the step criterion the way plunges,
retracts and links already are. **Do not.** That changes the measurement so the
path passes, and the 330° step is a load shock the cutter genuinely experiences —
the most physically real number in the group.

A plunge opens a hole of radius `r` on the *rim* of the first loop, where it does
the far side no good. Put the hole at the loop **centre** instead and it reduces
engagement at every position on every subsequent turn. With a void disk of radius
`V` concentric with a loop of radius `rho`, the void subtends
`2·acos((rho² + r² − V²) / (2·rho·r))` at the tool centre; immediately after the
plunge `V = r`, and that collapses to

```
peak engagement = 360° − 2·acos(rho / 2r)
```

which is *increasing* in `rho`. Over the non-degenerate range `rho > r` its
infimum is at `rho → r`:

```
360° − 2·acos(1/2) = 360° − 120° = 240°
```

**No non-degenerate first loop after a single plunge engages less than 240°.**
Measured against `_stock_2.engagement_at` on `rect_12x8` at ten radii, the closed
form and the exact predicate agree to the last reported digit; the first
admissible rung, `1.048 r`, measures 243.20°.

So the step floor is `240 − cap`, which at a 120° cap is *exactly* the criterion —
reachable only in the double limit of a first loop at the degeneracy boundary and
a second loop exactly at the cap, both open. Measured on `rect_12x8` a concentric
ramp descends 243.20 → 118.56 → 113.91 → 119.27 → 116.39 → 114.23 → 104.10, so
the largest step is 124.64°, against 329.49° for a rim plunge.
`compas_cgal.engagement_spiral_entry_toolpath` implements it, and pays for it:
concentric circles have no common tangent, so ramp turns link radially at a right
angle, two tangent breaks each.

## A true helical entry is outside the measurement model

The textbook remedy is to ramp in Z — orbit while descending, so each turn takes a
shallow axial bite and the hole opens with no full-width cut at all. This
framework cannot express it. `benchmarks.depletion._replay_kind` classifies every
operation against a single inferred cut plane and raises
`UnreplayableOperationError` on any move with both a Z change and XY travel:
"ramped 3D cutting is outside the cut-plane depletion model". A helical arc
emitted as a `Circle` above the plane classifies as `ReplayKind.RAPID`, removes
nothing in the replay, and leaves the bored hole invisible to the coverage grid —
so the first loop would still measure 360° and the reported engagement would be a
fiction.

A helical entry is therefore a change to the depletion model, not to a generator.

## Corollary 2: the bridge tilt is the clearance slope

Consecutive machining circles are linked by a straight bridge between their
entry points. With the entry placed one loop radius along the guide normal,
`entry_i = c_i + rho_i * n`, the chord between two entries on a straight guide is

```
chord = advance * t + (rho_(i+1) - rho_i) * n
```

so the angle between the chord and the guide tangent `t` is

```
atan(|d rho| / advance) = atan(d rho / d s)
```

— the arctangent of the **clearance slope along the guide**, which does not
depend on the advance. Shortening the advance shrinks numerator and denominator
together.

On a rectangle's corner bisector the clearance grows at `sin 45° = 0.7071` per
unit of travel, giving `atan(0.7071) = 35.26°` at *every* bridge on that chain.
Measured on `rect_12x8` (2 mm tool, 120° cap): all 24 tangent breaks are the
four corner chains' six junctions each, all at 35°, while the constant-clearance
spine has none.

!!! warning "Bounding the radius step alone makes this worse, not better"

    Capping `|d rho|` between consecutive loops forces shorter advances, which
    emits *more* bridges — each still tilted by the same 35.26°. A radius-step
    bound therefore trades a smaller `max_loop_radius_step` for a larger
    `tangent_breaks` unless the bridge geometry is fixed at the same time.

    The fix is to place the entry at the **external common tangent** of the two
    circles rather than along the guide normal: choose the unit normal `n` with
    `n . (c_(i+1) - c_i) = rho_(i+1) - rho_i` and put the entry at
    `c_i - rho_i * n`. The chord is then perpendicular to `n` and tangent to both
    circles by construction, for any radius difference — at the cost of requiring
    `|d rho| <= advance`, i.e. a clearance slope of at most 1. A 45° bisector's
    slope of 0.7071 satisfies it.

## Where this is used

- `benchmarks.quality.DEGENERATE_LOOP_RATIO` — the `rho <= r` boundary.
- `benchmarks.gate.L_ARM_TOOL_DIAMETERS` — the L's arm width, derived against
  the MAT-centred threshold `W > 4r` rather than the middle-curve one.
- `src/compas_cgal/adaptive/candidates.py` — the `(d - r) / 2` cap this page
  reads the identity off.
- `src/compas_cgal/engagement_rho_toolpath.py` — the degeneracy floor as an
  emission rule, the external-tangent entry, and `RhoToolpathResult.declined_regions`.
- `src/compas_cgal/engagement_spiral_entry_toolpath.py` — the centre plunge and
  concentric ramp, and the 240° floor as a shipped bound.
