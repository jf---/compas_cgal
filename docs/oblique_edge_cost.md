# The Oblique-Edge Cliff

`ReachableDomain2(...).center_domain()` costs **5 milliseconds on an axis-aligned
rectangle and 17.5 seconds on an oblique quadrilateral with equally plain integer
vertices** — a factor of about 3,500 for a change that alters nothing about the
pocket's size, complexity, or vertex count. Past that it stops returning in any
useful time: the same rectangle turned by a generic angle did not finish in 90
seconds. **Every pocket in the benchmark corpus is axis-parallel**, so every
performance number this repository reports — the Held-parity result included —
was measured in the single most favourable regime the geometry admits, and none
of them transfers to a pocket with a slanted wall until this is fixed.

!!! danger "This qualifies the parity claim"

    The 71.9 ms/pocket figure stands for what it measured. What it did not
    measure is any pocket with an oblique edge, and one oblique edge alone costs
    more than three hundred times that entire budget. Treat parity as
    established for axis-aligned pockets and **unmeasured** otherwise.

## The measurement

A four-vertex pocket, 1.0 tool radius, no holes. Only the vertex coordinates
differ.

| pocket | vertices | `center_domain()` |
| --- | --- | ---: |
| axis-aligned, integers | `(-6,-4) (6,-4) (6,4) (-6,4)` | **0.0050 s** |
| axis-aligned, one decimal | `(-6.1,-4.1) (6.1,-4.1) (6.1,4.1) (-6.1,4.1)` | 0.0158 s |
| **oblique, integers** | `(-6,-4) (6,-2) (6,4) (-6,2)` | **17.53 s** |
| oblique, exact 3-4-5 rationals | `(-2.4,-6.8) (7.2,0.4) (2.4,6.8) (-7.2,-0.4)` | **> 60 s** |
| rotated by 30, 45 or 90 degrees via `math.cos` | — | **> 90 s** |

## What it is not

Each of these was the obvious explanation and each is wrong.

**Not coordinate bit length.** The oblique case that costs 17.5 seconds has
ordinary small integer vertices — shorter decimals than the axis-aligned case
that costs 16 milliseconds. Rounding a rotated rectangle to a single decimal
place does not rescue it either; it still times out.

**Not near-degeneracy.** A 90-degree rotation through `math.cos` leaves edges
within `6.12e-17` of axis-parallel, which is a plausible degeneracy trap. But 45
and 30 degrees are nowhere near degenerate and time out identically.

**Not pocket complexity.** Four vertices, no holes, one convex quadrilateral, in
every row of the table.

**Not the stock depletion.** `Stock.subtract_capsule_quad` runs in under a
millisecond per cut on all three coordinate regimes, and `coordinate_digits()`
grows almost identically for each — a mean of 49 digits after four cuts on the
dirty input against 45 on the clean one. The cost is not in the arrangement the
cuts build.

## What it is

Isolated to one call. `survey_path` builds the reachable centre domain once, and
that single construction is where the time goes; `engagement_at` and
`ExactRegion2.contains` are microseconds on all inputs. The domain is the pocket
eroded by the tool radius, built with `Gps_circle_segment_traits_2`.

!!! warning "Mechanism not established"

    The obvious hypothesis — that each oblique edge's inward offset introduces
    its own radical `sqrt(dx^2 + dy^2)`, so the one-root `Sqrt_extension`
    coefficients compound across edges — is CONTRADICTED by the 3-4-5 row, whose
    edge lengths are rational (9.6, 7.2) has length exactly 12, and which is
    nonetheless the slowest of the finite cases. Something else is driving it.
    The measurement is solid; the explanation is not, and this page does not
    offer one it cannot defend.

## Why nothing caught it

The corpus was built around `rect_12x8`, `rect_20x12` and `L_shape`. All three
are axis-parallel, and so are the analytic, complexity, neck, precision and
topology families derived from them. `benchmarks/congruence.py` does apply
Pythagorean rotations, which keep coordinates rational — and by the table above
rational coordinates are not what saves you, so those cases are the slow kind
too wherever they are exercised on this path.

The finding came out of a property test asking whether a pocket turned on the
table is the same machining problem. It is the argument for the property suite
in one sentence: the assumption "pockets are axis-aligned" was invisible because
every example anyone wrote shared it.

## What follows

1. **Do not quote a performance figure without stating the pocket's edge
   directions.** Every number in `BENCHMARK.md` and in the review memo needs this
   qualification until the cliff is understood.
2. **Add an oblique pocket to the corpus.** Not a rotated rectangle, which is a
   special case, but a genuinely slanted wall — and expect it to be slow. A
   family that cannot run is still worth having, because it makes the limit
   visible instead of absent.
3. **Profile `center_domain()` symbolicated.** The same relink that settled the
   continuous-engagement question applies here: drop `-Wl,-S -Wl,-x` from the
   ninja link line, or the CGAL template frames will read as `???`. See
   `continuous_engagement_cost.md`.
4. **Then decide whether the domain is needed at all.** It is used for one thing
   — deciding whether a cutter centre is legal. A cheaper sufficient test, or an
   incremental construction, may be available once the cost is understood.
