# Engagement zone-query divergence: reporting-only, ≤48 ulps, certificate perfect

**BLUF.** `engagement_at_zone` reproduces `engagement_at`'s **exact certificate**
(`cap_exceeded`) **bit-for-bit in 100% of cases** (0 mismatches / 7200 comparisons
across caps 90°/120°/180° and gap-closure γ=30°). The ONLY divergence is in the
two **reporting doubles** (`total_tea`, `max_run_tea`): ≤ **48 ulps** (~1.3e-15 rad)
on ~6% of probes. **Root cause is proven**: the overlay (Gps boolean
`intersection`) and the local zone (`Arrangement_zone_2` sweep) construct the
**same crossing points** — verified *exactly equal* by CGAL's representation-invariant
one-root `operator==` — but with **different `Sqrt_extension` representations**, and
`to_double(a0 + a1·√root)` evaluates the specific triple, so algebraically-equal
points round to ulp-different doubles → ulp-different `atan2` spans. The exact
certificate uses exact predicates on the one-root form (representation-**invariant**),
so it is unaffected. **This is the kernel's deciding/reporting split working exactly
as documented.** Bit-exact reporting is *unachievable* across two different exact
construction algorithms; it is not a defect and cannot be "fixed" without forcing
identical point construction.

This is NOT a certificate divergence (the failure mode the gate guards against). The
harvests do not "yield different arcs" — they yield the *same* arcs with different
coordinate spellings.

## Evidence

| Claim | Measurement | Source |
|---|---|---|
| Certificate (`cap_exceeded`) reproduced exactly | **0 mismatches / 7200** comparisons (5 witness stocks × 4 radii × 3 cap/gap settings × 120 probes) | `diag_zone.py` |
| Divergence is reporting-only, ulp-scale | `total_tea` 432/7200 diverge, `max_run_tea` 390/7200; **max 32 ulps** (4.4e-16 rad) | `diag_zone.py` |
| Zone finds the IDENTICAL crossing geometry | **0** unmatched non-x-extreme zone endpoints / 4462 probes — every crossing has an *exactly equal* (`==`) overlay counterpart | `harvest_compare_debug` (temp) |
| …but with different one-root representations | 1133/4462 probes: exactly-equal (`==`) endpoints had **nonzero `to_double` gap** (≤ 1.78e-15) | `harvest_compare_debug` (temp) |
| Span-from-endpoints ("Option B") does NOT fix it | reporting still diverges (max 48 ulps) — the endpoints' `to_double` already differ, so no reporting recomputation can bit-match | `diag_zone.py` after Option-B prototype |
| `engagement_at` (overlay) behavior unchanged by the refactor | full suite **165 passed** (incl. the differential oracle) after `finish_engagement` extraction | `pytest tests/ --ignore=test_engagement_zone.py` |
| Query speedup is real | **3.8–5.1×** faster (overlay 2.5–6.8 ms/q → zone 0.65–1.5 ms/q) on 10/25/50-cut pockets; a lower bound — real audit pockets have far larger arrangements | `bench_zone.py` |

## Why "reproduce EXACTLY on all three outputs" was an over-strong gate

The plan assumed both harvests would produce identical `Arc` sets and therefore
identical doubles. Reality: they produce arcs with **algebraically identical**
endpoints (proven) but **different one-root representations**. `to_double` is a
projection of the *representation*, not the *value*, so it is not bit-reproducible
across the two exact algorithms. Only quantities decided by exact predicates on the
one-root form are representation-invariant — and the certificate is exactly such a
quantity, which is why it matches perfectly.

The certified path's guarantee lives entirely in `cap_exceeded` (exact) and the
`certify_segment_tea` verdict built from it; `total_tea`/`max_run_tea` are explicitly
REPORTING doubles that "never feed a decision" (`engagement_2.h`, `docs/exactness.md`
deciding/reporting split). A ≤ 1e-15 rad difference in a reporting-only quantity,
with the certificate bit-identical, is the split behaving as designed.

## Recommendation (pending user decision)

1. **Comparison test asserts the TRUE invariant**: `cap_exceeded` EXACT on every
   probe (the certificate — the load-bearing output) + `total_tea`/`max_run_tea`
   within a tight, named representation-artifact bound (ulp-scale, e.g. an
   absolute `1e-12` rad — 1000× the observed 1.3e-15 worst case, still 1e10×
   tighter than the 0.02-rad geometric tolerances the reference tests use). This
   asserts what is actually guaranteed; it is NOT widening a geometric tolerance.
   Rejected alternative: keep bit-exact `==` on all three — provably unachievable.
2. **Swap is safe**: post-swap `engagement_at` reports shift ≤ 1e-15 rad; the 165
   tests + oracle use tolerances (stay green); the certificate is bit-identical.
   Gated on user permission per the add-alongside/remove-with-permission rule.

## Status — RESOLVED (2026-07-24, swapped)

Diagnosis accepted; the reporting divergence is the deciding/reporting split working
as designed. `engagement_at` was swapped onto the local zone harvest and the overlay
harvest removed (user-approved). The redundant `engagement_at_zone` entry and the
migration comparison test were removed with it — ongoing zone-vs-reality equivalence
is covered by the differential Monte-Carlo oracle (`test_engagement_oracle.py`) plus
the 165 known-answer verdict tests. Post-swap: **165 passed**, certificate bit-for-bit
unchanged, suite 8.5 s → 4.8 s. Temporary `harvest_compare_debug` diagnostic removed.
The certificate-vs-reporting distinction proven here is the load-bearing takeaway: an
exact CGAL kernel reproduces DECISIONS across construction paths, never `to_double`
reporting bit-for-bit.
