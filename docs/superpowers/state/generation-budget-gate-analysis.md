# Held generation budget: gate analysis (2026-10-22)

**Monstera passes on min-of-5 CPU (8.85 s against 9.0 s). Figure 5 fails (2.33 s against
1.7 s) and crossed skis fails (2.38 s against 2.2 s).** These numbers were measured at load
30–85. The Python levers P1/P2/P5/P6 are already in the tree. On both failing workloads the
whole overrun is the NOTHING-circle stage, and that stage did not exist when the budgets
were set. Without it the totals are 1.54 / 2.14 / 7.95 s, all inside budget. The other two
stages added since the budgets were set cost at most 0.27 s together. The NOTHING stage
spends 60 % of its time in its binary64 sample screen, not in exact union tests. So the
route to green is a faster NOTHING screen (NU). It does not need a budget change, a margin
or an epsilon.

One decision for the user (section 6): does the NOTHING stage count against the generation
budget? Recommendation: yes, with no budget change, so Figure 5 stays red until NU lands.
Budgets, gate and tolerances in this page are unchanged: 1.7 / 2.2 / 9.0 s, min CPU of
5 fresh serial processes, the exit collect counted (user rulings 2026-10-22).

## 1. Current numbers and content identity

Tool: `tools/held_generation_clock.py` (the gate's own clock), 5 serial fresh processes per
workload, through `build/ga/run_clock_snap.sh`, log `build/ga/clock_snap.log`, summary
`build/ga/summarize.py`.

Content identity: HEAD `fa28023f` plus the uncommitted backlog, frozen at 15:56 into
`build/ga/snap/` (benchmarks/ + tools/, `.py` content hash `c7103125c9b7`). The native
extension is the live editable build. One file was swapped: `benchmarks/held_nothing_circles.py`
is NU's pre-edit backup (`build/nu/backup_held_nothing_circles.py`, sha1 `7fc26a3a…`). That
is the version BG measured. The reason for the swap is that NU rewrote the live file at 15:44,
mid-flight, and on the live tree the stage took 48–52 s CPU on Figure 5 (section 4). PF's
levers are in the snapshot: the generator files carry PF's 13:57–13:58 edits relative to
`build/pf/before/`.

| workload | budget | min CPU (gated) | median | max | median wall | 1-min load | digest (circles) | verdict |
|---|---:|---:|---:|---:|---:|---|---|---|
| figure5_refined (2004) | 1.7 | **2.33** | 2.34 | 2.47 | 6.45 | 50–56 | `8d9f1d7d1ca8` ×5 | over by 0.63 |
| crossed skis (1621) | 2.2 | **2.38** | 2.41 | 2.45 | 3.73 | 32–48 | `711eb026b836` ×5 | over by 0.18 |
| Monstera (4428) | 9.0 | **8.85** | 8.91 | 9.29 | 23.87 | 30–85 | `2a57a9a9fcdd` ×5 | pass, 1.7 % under |

BG's min-of-5 at load 70–218, taken before the PF levers were in its measured state, gave
2.47 / 2.45 / 9.36 s.

**Load inflates CPU too, by about 15–20 %.** PF measured the same levers at load 8–11:
`initial` took 1.00 / 1.47 / 5.31 s and total CPU was 1.32 / 1.85 / 6.90 s, before NOTHING.
In my runs at load 30–85, `initial` took 1.19 / 1.71 / 6.15 s. Minimum-of-5 removes spikes,
but it does not remove load that is sustained across all 5 runs.

## 2. Stage CPU breakdown (min run; s)

| stage | fig5 | skis | Monstera | in budget basis (2026-10-21)? |
|---|---:|---:|---:|---|
| guide | 0.049 | 0.058 | 0.634 | yes |
| offset | 0.002 | 0.004 | 0.073 | yes |
| initial | 1.190 | 1.711 | 6.153 | yes |
| refine (includes the connector repair gate) | 0.229 | 0.275 | 0.897 | yes, without the repair gate |
| **nothing** | **0.789** | **0.242** | **0.899** | **no** |
| connectors | 0.009 | 0.007 | 0.020 | yes |
| connector_certificate (post-gate) | 0.027 | 0.023 | 0.065 | no |
| exit collect | 0.040 | 0.060 | 0.111 | moved, not new |
| **total** | **2.33** | **2.38** | **8.85** | |
| total without `nothing` | 1.54 | 2.14 | 7.95 | |

The per-stage medians are within 0.04 s of these values, except Monstera `initial`, whose
median is 6.24 s.

## 3. Stages added since the budgets were set

The budgets come from the commit `e6879ad1` (2026-10-21 23:19). They are the "lever budget"
column of `lever-hypotheses-selection-plan.md` §"Time budget for this lever": quiet
single-core estimates for guide + offset + initial + refine + connectors.

| added | commit | cost fig5 / skis / Monstera (CPU s) | source |
|---|---|---|---|
| connector repair + seam gate inside refinement | `8196e1fd` (10-22 14:41) | ≈ 0.05 / 0.04 / 0.20 | PP, instrumented (`_ConnectorGate`, `require_no_sharp_corner`, `close_seam`); adds 0 circles |
| connector engagement post-gate (`connector_certificate`) | `fa28023f` (10-22 14:41) | 0.027 / 0.023 / 0.065 | this run |
| NOTHING-circle stage (`nothing`) | pass `69eb2eeb`; wired by PN, **uncommitted** | **0.79 / 0.24 / 0.90** | this run |
| exit `gc.collect()` counted | `9ffcf215` + BG ruling | 0.04 / 0.06 / 0.11 | moved out of the stages, not new work |

The first two together cost ≤ 3.5 % on every workload. The NOTHING stage costs 34 % /
10 % / 10 % of generation CPU. On Figure 5 it removes 6 of 2010 circles. On skis it removes
0 of 1621 circles, and it still costs 0.24 s.

**Where the NOTHING stage spends its time.** From NU's pre-change profile,
`build/nu/prof_before.log`, at load 75–129:

| part | fig5 | skis | Monstera |
|---|---:|---:|---:|
| select (numpy overlap distance) | 0.05 | 0.04 | 0.16 |
| **screen** (binary64 sample grids 48×5 → 128×9 → 512×17, every circle) | **0.54** | **0.19** | **0.61** |
| exact union tests (`subtract_circle_sweep` + `circle_sweep_removal`) | 0.34 (14 tests, 6 NOTHING) | 0.04 (2 tests, 0 NOTHING) | 0.21 (18 tests, 10 NOTHING) |
| total | 0.94 | 0.28 | 1.01 |

The screen is the cost. It is a keep-only filter, so it decides nothing. The exact unions
that prove the removed circles cost about 0.15 s on Figure 5: 6 of the 0.34 s in 14 tests.

## 4. Levers in flight and queued

Savings are CPU seconds. "Measured" means a same-tree A/B. "Est." means a projection from
a profile.

| lever | slot | status | saving fig5 / skis / Monstera | byte-identical |
|---|---|---|---|---|
| **PF**: P1 parallel offset sites, P2 integer start seed, P5 carried positions, P6 drop the float self-intersection check | Python | **in the tree**, uncommitted; PF report not yet written | **measured** 0.18 / 0.40 / 1.70 (1.49→1.32, 2.25→1.85, 8.67→6.90 at load 8–11, 3 runs each, `build/pf/{before,after}.jsonl`) | yes: pre-NOTHING digests `eaffef8a` / `711eb026` / `04403396` unchanged |
| **NU**: faster exact NOTHING screen | Python | **in flight, live file mid-edit**; NU report not yet written | target: the screen's 0.54 / 0.19 / 0.61 plus the false union tests; floor about 0.2 / 0.05 / 0.25 s for the stage (est.) | must keep the removed set; NU's A/B checks `same=True` |
| ↳ NU state at 15:44 (observed, not a verdict on NU's end state) | | | **regression** 30–60×: fig5 50 s, skis 9.6 s, Monstera 18.9 s for the stage; union tests drop 14→11 / 2→0 / 18→10 (`build/nu/ab1.log`, and 47.8–51.9 s in my first live runs, `build/ga/clock.log`) | same removed set |
| P3: no connector build on the ordering-only call 0 | Python | **needs a ruling** (it removes an unreachable fail-loud check) | est. 0.05 / 0.05 / 0.21 | yes |
| N1: native nearest-contact batch | C++ | queued | est. 0 / 0.30–0.38 / 1.2–1.4 | yes, with a correctly rounded parameter and a contract test |
| N2: station batch returns only what the consumer reads | C++ | queued | est. 0.13–0.17 / 0.13–0.17 / 0.5–0.7 | yes |
| N3: native ownership argmin | C++ | queued | est. 0.10 / 0.10 / 0.35–0.40 | yes |
| N4: native span selection | C++ | queued; needs an exactness design note first | est. 0.15 / 0.18 / 0.45 | yes, after the note |
| native NOTHING batch union | C++ | queued (C++ queue after CP1) | est. ≤ 0.3 / 0.04 / 0.2 (the union share) | same removed set |

## 5. Reachable result per workload

The quiet basis (load ≈ 10) is PF's total of 1.32 / 1.85 / 6.90 s, plus the NOTHING stage
scaled from load 50 by the measured 0.85 load ratio: about 0.67 / 0.21 / 0.76 s.

| state | fig5 (1.7) | skis (2.2) | Monstera (9.0) |
|---|---:|---:|---:|
| today, measured at load 30–85 | 2.33 ✗ | 2.38 ✗ | 8.85 ✓ |
| today, quiet est. | ≈ 1.99 ✗ | ≈ 2.06 ✓ (6 %) | ≈ 7.66 ✓ (15 %) |
| + NU at its floor (stage ≈ 0.2 / 0.05 / 0.25) | ≈ 1.52 ✓ (11 %) | ≈ 1.90 ✓ | ≈ 7.15 ✓ |
| + P3 (ruling) | ≈ 1.47 | ≈ 1.85 | ≈ 6.94 |
| + N1 | ≈ 1.47 | ≈ 1.52 | ≈ 5.7 |
| + N2 + N3 + N4 | ≈ 1.1 | ≈ 1.1 | ≈ 4.2 |

Figure 5 is the binding case. It passes only if the NOTHING stage costs ≤ 0.38 s quiet,
or ≤ 0.16 s at load 50. That is under half of today's cost. No C++ lever except the native
batch union helps Figure 5's NOTHING stage, because Figure 5 has no contact stage.

## 6. Decision points for the user

**User rulings 2026-10-22.** D1 is ruled (a): the NOTHING stage counts inside
the existing generation budget, which is unchanged. D2 is ruled (b), made strict:
only quiet-machine runs give a budget verdict. Before each run the gate reads
the 1-minute load. If it is above 2, the gate raises `MachineNotQuietError` (a
named error, not a skip) instead of measuring. Runs on a loaded machine are
indicative only (`python -m tools.held_generation_clock --indicative ...`) and
never pass or fail the budget. The gate statistic stays min CPU over 5 fresh
runs. Implemented in `tools/held_generation_clock.py` and
`benchmarks/acceptance/test_held_generation_budget_slow.py`. The first
indicative run after NU's faster NOTHING pass is in section 8.

**D1. Does the NOTHING stage count in the generation budget?**

- (a) Count it in the same budget, unchanged. This is how the test is written today; its
  docstring lists the stage. Figure 5 stays red until NU or the native union lands.
- (b) Account for it separately. The generation budget covers placement + refinement +
  connectors, which is the budget's 2026-10-21 basis. NOTHING gets its own per-workload
  stage budget, set before code per the time-budget rule. Figure 5 and skis turn green
  today.
- (c) Raise the generation budgets to absorb it. This is a budget change, and it is listed
  only for completeness.

**Recommendation: (a).** The budget measures the wait before a path is ready to cut. The
NOTHING stage is part of what the generator emits: it changes the digest, the path and the
connectors. Under (b), a stage that spends 0.8 s to remove 6 circles stays green by
construction, which is the slow-stage-recorded-as-baseline failure the time-budget rule
forbids. The per-stage CPU is printed by the gate, so attribution needs no separate budget.
If the user also wants NU held to a number, add a NOTHING stage target beside (a), not
instead of it: ≤ 0.38 s quiet on Figure 5 is what (a) needs.

**D2. Measurement under sustained load.** Min-of-5 still reads 15–20 % high at load 30–85.
Monstera sits 1.7 % under its budget, so it will flip with load. Options:
- (a) Keep the gate as is, and report every verdict with its load.
- (b) Make the verdict run a quiet-machine run: record the load, and treat a red at high
  load as unconfirmed until it is re-run quiet or shown by a same-process A/B ratio.
- (c) Normalise against a fixed reference workload measured in the same window.

Recommendation: (b). It changes neither the budget nor the statistic, only when a verdict
counts. (c) adds a second proxy.

**D3. NU's direction.** At 15:44, NU's in-flight variant traded 3–8 union tests for a
screen that costs 30–60× more. Its removed sets are identical, but it costs far more than it
saves. Whatever replaces the screen must be cheaper than the union tests it avoids. On skis
and Figure 5 the screen alone already costs more than every union test on that workload put
together (section 3).

**D4. `SCREEN_RELATIVE_SLACK = 1e-9` in the NOTHING screen.** The screen only keeps circles,
so this constant cannot cause a wrong removal. But the common constraints require filters to
be provable supersets by outward rounding, with no slack constant. NU's rework is the moment
to replace it with outward-rounded bounds, or with the exact union test on the circles that
pass. Needs a ruling only if NU keeps it.

**D5. Open rulings that move these numbers:**
- P3, which removes an unreachable check: −0.05 / −0.05 / −0.21 s.
- PN concern 1: the NOTHING pass is not on the certified coverage-repair path, and adds
  about 1 s per certified Monstera run if it is placed there.

## 7. Not acted on

- `build/pf/tests.log` (15:15) shows many failures in the last 10 % of a test run. Its owner
  and cause were not checked, because PF has not reported. It should be checked before PF's
  levers are committed.
- Ω* generation is not gated. PF's load-9 numbers for it are 4.69 / 7.18 / 17.35 s CPU, and
  `connector_repair` takes 1.5 / 2.4 / 5.8 s of that.

## 8. Indicative run after NU's NOTHING pass (2026-10-22 16:18, load 10.6–13.4)

**Indicative only. These runs are not a verdict, because the 1-minute load was
10.6–13.4, far above the quiet limit of 2.** The tree is HEAD `fa28023f` plus
the uncommitted backlog, with live `benchmarks/held_nothing_circles.py` (sha1
`83f8605a…`). Each workload ran 5 serial fresh processes. Log:
`build/qg/indicative.log`.

| workload | budget | min CPU | median | max | median wall | NOTHING stage | digest (circles) |
|---|---:|---:|---:|---:|---:|---:|---|
| figure5_refined | 1.7 | 1.72 | 1.73 | 1.79 | 1.82 | 0.34–0.37 | `8d9f1d7d1ca8` (2004) ×5 |
| crossed skis | 2.2 | 2.00 | 2.01 | 2.11 | 2.07 | 0.07–0.09 | `711eb026b836` (1621) ×5 |
| Monstera | 9.0 | 7.63 | 7.92 | 8.15 | 8.80 | 0.43–0.49 | `2a57a9a9fcdd` (4428) ×5 |

NU's pass cuts the NOTHING stage from 0.79 / 0.24 / 0.90 s (section 2) to about
0.35 / 0.08 / 0.45 s. The digests are unchanged. Figure 5 is 0.02 s above its
budget at load about 11. Section 1 measured a 15–20 % inflation from load, so
a quiet run would likely bring it under. Only a quiet run decides.
