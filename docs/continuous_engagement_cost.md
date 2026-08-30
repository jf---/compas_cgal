# Continuous Engagement: the Cost Verdict

The continuous event-partition oracle is **~116× slower than this repository's own sampled
certifier** and roughly **five orders of magnitude** off Held & Pfeiffer, and that is a
formulation problem rather than a tuning problem. On current evidence the shippable certifier is
the sampled + growth-bound-guarded one in `engagement_2.cpp`; `continuous_tea_2` is a research
result, not a component, until its algebraic formulation changes.

Measured on axis-parallel pockets only; `center_domain()` is ~3,500× slower
on oblique geometry with the mechanism not yet established — see
*The Oblique-Edge Cliff*. Treat every figure on this page as best-case with
respect to edge direction.

| certifier | per segment motion |
| --- | ---: |
| `engagement_2.cpp` — sampled + guard, local zone query | **129 ms** |
| `continuous_tea_2` — `audit_segment_tea_event_exact` | **~15,000 ms** |
| Held & Pfeiffer — an *entire pocket* | **3–100 ms** |

The guarantee the continuous partition adds is real: no guard constant, no safety factor, no
refinement floor, no spacing-exhaustion failure mode. That gain is not in dispute. What this page
records is its price, measured, so the trade is decided on evidence.

## Measured scaling

Fixed probe segment, tool radius 1.0, cap 120°, against a 6×4 pocket depleted one tool disk at a
time near the probe:

| cuts | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| stock vertices | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
| audit (s) | 0.124 | 0.741 | 1.928 | 3.825 | 5.271 | 6.603 | 10.048 | 12.028 | **14.698** |

Superlinear on a pocket with twenty vertices. A production pocket carries thousands of boundary
features and hundreds of motions.

## What the cost is not — falsification order

Each hypothesis below was measured and killed. They are recorded because each is the obvious
guess, and re-testing them costs hours.

**Not rational bit length.** `max_digits` saturates at 127 by the second cut and stays flat while
the audit time grows 8×. Bit growth is present and is not the driver.

**Not distant stock.** Eight cuts placed *away* from the probe segment leave the audit flat at
0.13 s while producing the same arrangement size (V=13, HE=26) that costs 5.3 s when the cuts are
*near* it. Cost tracks contacted features, not arrangement size.

**Not string-mediated rationals.** This was the leading hypothesis after an unsymbolicated
profile, and it is wrong. `extract_boundary_records(stock)` — the only site where large
arrangement coefficients are decimal-serialised — measures **0.8 ms** against a **17.5 s** audit.
Removing all 50 `parse_rational(x.text())` round trips left the timing sweep **unchanged**.

!!! warning "Unsymbolicated profiles will mislead you here"

    `nanobind_add_module` links with `-Wl,-S -Wl,-x`, so every header-only CGAL template
    instantiation compiled into the extension appears as `???`. A `sample` profile therefore
    shows *zero* algebraic-kernel symbols and a top-heavy list of `malloc`/`free`/`__udivmodti4`,
    which reads convincingly as string and decimal-conversion churn. It is not: `__udivmodti4` is
    `divide_unsigned_helper` inside bignum division and GCD, and the allocator traffic is
    `cpp_int_base` limb allocation. To get real symbols, drop those two flags from the ninja link
    line and relink.

## What the cost is

A symbolicated profile, aggregated by enclosing component:

| component | samples |
| --- | ---: |
| `Curve_pair_analysis_2` | 27k |
| `Curve_analysis_2` | 25k |
| `Algebraic_curve_kernel_2` | 18k |
| resultant + modular GCD | ~15k |
| Bitstream-Descartes | ~12k |
| `Shear_transformation` | 9k |

The work is CGAL's **bivariate** algebraic curve kernel (`CGAL::Algebraic_kernel_d_2`,
`src/exact_algebraic_1.h`), over boost multiprecision.

## The two levers, and why neither is sufficient

**Arithmetic backend.** The build deliberately pins the pure-C++ boost backend —
`CMakeLists.txt` sets `-DCGAL_DISABLE_GMP -DCGAL_USE_BOOST_MP`, and `src/exact_algebraic_1.cpp`
guards it with `#error` on `CGAL_CORE_USE_GMP_BACKEND`. GMP's `mpq` is typically several times
faster on precisely this workload. This is an owned trade (portable wheels, no external
dependency), not an oversight — and a small multiple does not close five orders.

**Shearing and resultant degree.** `Shear_transformation` at 9k samples means the kernel is
meeting non-generic curve positions and shearing to escape them, which inflates coefficient size
sharply. Avoiding the shear, or reducing resultant degree, attacks the scaling rather than the
constant — this is the only lever that could matter at the required magnitude.

**The open question underneath both:** the events along a **one-parameter** motion are
univariate. A bivariate curve-pair analysis is being used for what is, geometrically, a
univariate root-isolation problem — cutter-versus-line and cutter-versus-arc are degree ≤ 2 in
the motion parameter. Whether the bivariate formulation is *necessary* is the question to answer
before any further optimisation of this path.

## Retained from the investigation

- **Representational cleanup** (`perf/exact-rational-representation`): caching the exact rational
  in the source types, one shared text decoder, rational-taking `StationEventSource2::build`.
  Worth **5–15%**, larger at low stock complexity, asymptotically nil — and a net removal of
  ~260 lines. Kept for clarity, not for speed.
- **Content-keyed memo** (`a4b1a3e`): eliminates genuinely duplicate native audits, keyed on
  stock lineage + boundary digests + motion + radius + effective cap, deliberately excluding the
  operation ordinal so a cached partition can never carry another operation's identity.
- **Instrumentation**: `Stock2::arrangement_stats()` and `Stock2::coordinate_digits()`, which is
  what falsified the bit-length hypothesis. `coordinate_digits()` calls `.exact()` and so
  perturbs timing — it is a diagnostic and must never run on a timed path.

## Recommendation

1. Ship the **sampled + guarded** certifier. Its conservatism is an explicit, bounded, documented
   pessimism with a stated safe-failure direction — not an approximation — and it is one to two
   orders from Held rather than five.
2. Treat `continuous_tea_2` as a research branch. Do not spend further effort on backend swaps or
   representational tuning; both are measured and insufficient.
3. Answer the univariate question before reviving it. If the event set along a linear motion can
   be isolated univariately, the formulation changes and the cost argument restarts from scratch.
