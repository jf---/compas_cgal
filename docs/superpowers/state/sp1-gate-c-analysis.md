# SP1 gate (c) analysis — engagement-audit cost blowup at Task 8

**BLUF.** The engagement audit is far more expensive than the kickoff's ">120 s/square"
estimate. Measured this session: the *smallest* benchmark pocket (kite, 335 ops) exceeds
**5 min**; square (271 ops) exceeds **10 min**; both runs were killed before completing.
The full 7-pocket baseline (~2,440 ops) would take an estimated **30–90 min**, and Task 9's
**≤10 s/pocket gate (c) fails by ~50×+**. Per kickoff §5/§10 this presents options before any
fallback work; the exact-kernel bar rejects any approximate-clipper shortcut outright.

## Data (measured this session)
- Generation is instant (0.00–0.02 s/pocket) — **all** cost is in `audit_toolpath_engagement`.
- Op counts (cut ops): square 271 (260), kite 335 (324), l_shape 300 (277), irregular 402
  (358), island 389 (344), star 611 (564), dumbbell 732 (698).
- Geometry mix per pocket: **Circles 132–322** (each a full trochoid loop → ~20 fixed
  arc-certifier stations) + **Lines 139–344** (adaptive `certify_segment_tea`, pessimism-
  inflated on thin webs) + a few Arcs.
- Timing bounds: square **>600 s**, kite **>300 s** (both killed before completion).
- Pre-pessimism reference (Task 6b, d01879f): square audit **>120 s**. Post-repair **>600 s**
  → the soundness repair ~**5×**'d it, consistent with thin-web line over-refinement (the arc
  review measured one thin-web certify at 511 stations / ~21 s).

## Cost structure — decomposition NOT yet done (audits killed before instrumenting)
Two coupled drivers:
1. **Certify stations** — ~20/circle (fixed) + adaptive line stations (pessimism-inflated).
   *Reducible* by coarser exact grids.
2. **Exact-boolean stock depletion** (§10's stated primary driver) — each of 260–698 cut ops
   subtracts from an accumulating Epeck stock; per-op cost grows with accumulated subtractions.
   *NOT* reduced by coarser station grids.

**Load-bearing unknown:** which driver dominates. The chosen fallback's FIRST step must
decompose it (bounded instrumented run: per-op stations + time-vs-op-index growth). If
boolean-bound, coarser grids won't help — exact stock simplification / spatial indexing is the
real lever.

## Two separable problems
- **(P1) Baseline generation** — a one-time committed artifact; expense tolerable *if* it can
  complete (session-limit killings are the acute obstacle).
- **(P2) Subprocess test** — NO benchmark pocket audits fast enough for a snappy test. Solve
  orthogonally: a tiny synthetic pocket (few ops) or direct unit tests of `build_report` /
  `render_markdown`. Controller will handle P2 regardless of the P1 decision.

## Options (exact only; approximate clipper REJECTED per Exact-Kernel Discipline)
- **[A] Coarser exact station grids** (kickoff §10 pre-decided). Arc step 0.05→~0.15 (~7
  stations/circle), coarser line floor; stays exact; regenerate baseline + run gate on the
  coarser *sound* certifier. Contingent on the decomposition confirming certify-boundedness.
- **[B] Representative subset baseline now** — audit the tractable small pockets at full
  fidelity, document the full-set expense, defer the full 7-pocket baseline to SP2. Exact,
  honest, partial coverage.
- **[C] Accept full expense** — generate the true 7-pocket baseline once in the background
  (30–90 min, session-limit risk), commit it; Task 9 gate (c) then FAILS honestly and is
  documented as the SP1 result.

**Recommendation:** [A] per the kickoff's pre-decision — but its first step is the certify-vs-
boolean decomposition; if boolean-bound, escalate to stock-cost reduction rather than accept a
coarser certificate that doesn't actually help. (P2) is solved independently with a tiny
synthetic test pocket.

## UPDATE — decomposition result: option [A] is DEAD (measured)

Bounded 60-op kite audit, `AUDIT_ARC_STEP_FRACTION` 0.05 vs 0.15:
- 0.05: Circle avg **20.0** stations, **3611.8 ms**/op, 60-op subtotal 108.4 s. Line 7.0 st, 41 ms.
- 0.15: Circle avg **7.0** stations, **3668.2 ms**/op, 60-op subtotal 110.0 s.

**Coarsening the arc grid gave ZERO speedup** (3612 → 3668 ms/circle while stations dropped
20 → 7). The `engagement_at` station queries are ~free; the entire ~3.6 s/circle is the single
**`subtract_arc_sweep` exact-boolean depletion**, and it GROWS (27 ms → 669 ms over the first
13 circles) as the scalloped stock accumulates — **O(n²)** in cut count. This is inherent to
auditing a dense trochoid exactly; no station-grid coarsening (the kickoff's pre-decided
fallback) can touch it. Coarsening reverted (pure fidelity loss for zero gain).

**Compounding constraint:** the session usage limit (resets 18:00 Europe/Amsterdam) is killing
any run past ~5 min, so NO full benchmark pocket completes in this session regardless of approach.

### Revised options (the real fork)
- **[A′] Re-scope (RECOMMENDED).** SP1's thesis — the sound, merge-hole-free certifier +
  kernel — is DELIVERED (1c58382 + 3737a2a + d08f0ee, all reviewed, 164 green). Deliver SP1 on
  that; hand the full baseline AND the O(n²) exact-depletion bottleneck to SP2 as the crisply
  identified next problem. Task 9's gate (c) is documented as a characterized miss, not tuned.
- **[B′] Push the depletion lever now (R&D).** The one untested constant-factor lever is the
  `subtract_arc_sweep` under-covering disk-chain density (fewer disks/loop → cheaper subtract,
  at the cost of depletion accuracy → slightly over-measured downstream TEA). Does NOT change
  the O(n²) growth; real investigation; fights the session limit.
- **[C′] Defer generation to a post-limit session.** Accept the ~1 hr full-exact baseline, run
  it once in a fresh session after 18:00 when the limit resets. Complete + exact, just slow.
