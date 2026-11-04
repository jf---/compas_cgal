# U-1 report: station-0 filter fallback on the existing exact engine

Spec: `2026-11-04-actual-stock-refinement-design.md` §4, §5 (U-1). Worktree HEAD `ea4703bb`. Probes and
logs: `build/actual-stock-u1/` (untracked). Python only, on `Stock2` + `engagement_circle_{exceeded,prefixes}`.

**Verdict: U1 is feasible as specified.** The station-0 filter fallback rate is f = 0 on fig5, skis and
Monstera, every circle. The fallback term f · N · c_x is therefore 0 s against the one-third allowance
of 0.05 / 0.05 / 0.17 s. The targeted exact wall step is not needed. Nothing goes back to the user.

The result does not depend on how tight the C++ enclosure turns out to be. With the interval
enclosure widened 100× and 10⁴×, about 95 % of station-0 rims become ambiguous at the wall, which matches
the review's population. f stays 0, because the run merged across the wall gap is within the cap on
every circle. It is already within on the last-1 or last-16 view, which contain C_t.

## Numbers

f counts the decisions that the interval filter cannot certify: an interval-ambiguous wall event at
station 0 AND (pessimistic merged run exceeded while exact is within, OR exact exceeded). The second
case counts as a fallback outright, because the optimistic split is not emulated. That makes the count
conservative.

| | fig5 | skis | Monstera |
|---|---:|---:|---:|
| decisions probed (circles t ≥ 1, every circle) | 2,003 | 1,620 | 4,427 |
| engine calls N (refinement, gate-counted) | 2,010 | 1,623 | **4,485** |
| calls / final circle | 1.003 | 1.001 | **1.013** |
| wall-ambiguous station 0, interval model (×1) | 646 (304 cross / 37 tangent / 305 clear) | 684 (326 / 0 / 358) | 625 (266 / 44 / 315) |
| wall-ambiguous, enclosure ×100 | 1,901 | 1,543 | 3,514 |
| wall-ambiguous, enclosure ×10⁴ | 1,901 | 1,559 | 3,586 |
| pessimistic merged run exceeded on C_t (×1 / ×100 / ×10⁴) | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| exact station 0 exceeded on C_t | 0 | 0 | 0 |
| **f (×1 / ×100 / ×10⁴)** | **0 / 0 / 0** | **0 / 0 / 0** | **0 / 0 / 0** |
| tier resolving pessimistic-within (×1): last-1 / last-16 / full | 2,000 / 3 / 0 | 1,619 / 1 / 0 | 4,427 / 0 / 0 |
| ambiguous AND pessimistic exceeded on last-1 (an engine that falls back at the first tier instead of escalating), ×1 / ×100 | 0 / 5 | 0 / 1 | 0 / 0 |
| station-0 rims with an ambiguous rim/annulus-boundary event | 11 | 2 | 280 |
| … their largest station-0 engaged total | 4.9e-6° | 5.1e-7° | 3.5e-6° |
| largest station-0 engaged total, interval-ambiguous rims (cap 80°) | 78.39° | 76.63° | 79.55° |
| largest station-0 max run, all rims | 79.67° | 76.84° | 79.93° |
| min \|station-0 total − 180°\| (orientation leg) | 99.0° | 102.2° | 100.1° |
| k > 0 stations: min float \|d² − r²\|/r² to a wall (indicator) | 1.3e-6 | 6.2e-8 | 1.1e-7 |
| **circle-view last-1 refusals, all 45 stations, every circle** | 23 / 2,003 | 2 / 1,620 | **98 / 4,427** |
| fallback term f · N · c_x | **0 s** (≤ 0.05) | **0 s** (≤ 0.05) | **0 s** (≤ 0.17) |

Costs, process CPU, warm (the point locator is built by a discarded cold call), min of 5 per circle,
on C_t (the full circle-only window on P), every 20th circle. **Indicative only:** the 1-min load was
above 20 for all of them.

| | fig5 | skis | Monstera | load (1/5/15 min) |
|---|---:|---:|---:|---|
| c_x: exact station-0 decision, p50 / mean / max | 0.60 / 0.66 / 2.19 ms | 0.22 / 0.43 / 2.17 ms | 0.22 / 0.34 / 1.70 ms | 52→35 / 33→25 / 25→82 |
| 45-station decision, p50 / mean / max | 13.2 / 16.1 / 50.2 ms | 10.8 / 14.5 / 55.9 ms | 9.1 / 11.5 / 59.8 ms | same runs |
| c_x repeat (fig5, `costs.py`) | 0.59 / 0.68 / 2.19 ms | | | 86→196 |
| targeted wall step, Python `Fraction` upper bound, p50 / max (on every ambiguous crossing rim) | 225 / 543 µs (n 304) | 214 / 449 µs (n 326) | 219 / 509 µs (n 266) | 66→57 / 53→55 / 71→90 |

Quiet-window repeat: see the section at the end.

## Feasibility rule (§4)

f · N · c_x ≤ 0.05 / 0.05 / 0.17 s. With f = 0 the term is 0 on every workload, whatever c_x is.

Two bounds show how much room there is:

- **Break-even f at the measured mean c_x:** 0.05 / (2,010 · 0.66 ms) = 3.8 % on fig5,
  0.05 / (1,623 · 0.43 ms) = 7.2 % on skis, and 0.17 / (4,485 · 0.34 ms) = 11 % on Monstera.
- **Unprobed decisions.** The probe covers the final circles. The engine calls that are not final
  circles are refused candidates and circles that NOTHING later removes: 7 on fig5, 3 on skis and 58 on
  Monstera. Even if every one of them fell back at the largest measured c_x, the term would be
  0.015 / 0.007 / 0.10 s. That is still within the third.

The targeted step (exact existence plus foot test plus two neighbour orders, ≈ 0.22 ms in Python
`Fraction`) is not needed. It is recorded only as the measured cost of the alternative.

## What the measurement establishes, and what it does not

- **Interval model.** F = cross² − r²·L² (rim vs edge line, with the foot test dot ∈ [0, L²]) and
  G = |e − a|² − r² (rim vs vertex). Both are evaluated in numpy interval arithmetic, rounded outward by
  one `nextafter` per operation, which is a superset of `Interval_nt`'s enclosure. The ×100 and ×10⁴
  widenings stand in for a looser C++ formulation, for example one that goes through one-root
  coordinates. f = 0 at every width, so the verdict does not rest on the model.
- **The pessimistic measurement is live.** The gap-closure ratio is the exact wall-gap chord ratio
  4(r² − d²)/r², rounded up to a double. Take 60 crossing rims with more than one run and lower the cap
  between the max run and the total: on all 60, the closure makes the station exceeded where the exact
  call stays within (`liveness_figure5_refined.log`). f = 0 is therefore a property of the stock, not of
  a dead parameter. The cause: no wall-adjacent pair of runs at station 0 sums past 80°.
- **The other legs at station 0 have no ambiguity.**
  - Cap leg: the closest any station-0 max run comes to the cap is 0.07° (Monstera t = 4,194), which is
    about 2e-3 in chord ratio, many orders of magnitude above an interval width.
  - Orientation leg: no station-0 total comes within 99° of π.
  - Rim/annulus events: the ambiguous ones (11 / 2 / 280) all sit on stations with less than 5e-6° of
    material. Kept as material pessimistically, they stay far below the cap.
- **Stations k > 0 are outside f's scope (spec).** Their float distance to the nearest wall never comes
  within 6e-8 relative (no ambiguity class like station 0's). This is an indicator only: the station
  centres were reconstructed in binary64, not as the engine's exact rationals.
- **Not done.** The spec's U-1 red test ("the probe's record test") and a `tools/` probe would be
  tracked files, which this brief excludes. The probes live in `build/actual-stock-u1/` and need porting
  to `tools/` (U0 already ports the design probe) together with their pixi task, which goes through the
  pyproject slot.

## Consequences for the spec (not edited here)

- §4: the Monstera call count is measured at 4,485 calls / 4,428 final circles (1.013). The
  per-decision ceiling becomes 0.50 s / 4,485 = 111 µs, instead of 113 µs.
- §2.5: the Monstera circle-view last-1 refusal rate is 98 / 4,427 (2.2 %); every circle is measured
  on fig5 (23 / 2,003) and skis (2 / 1,620). Those circles escalate to last-16.
- §5: the targeted lazy-exact wall step is not triggered. U1 proceeds as specified, and U3's
  "fallback rate ≤ U-1's" gate is 0 on station 0.

## Reproduce

From the worktree root. Each probe runs niced through the reviewer's wrapper, which sets
`SKBUILD_EDITABLE_SKIP` and never triggers a native rebuild.

```bash
W="figure5_refined figure8_crossed_skis_boundary figure8_monstera_boundary_qualified_corners"
R="nice -n 10 pixi run --no-install sh build/pyrun.sh build/actual-stock-u1"
for w in $W; do $R/build_workload.py $w; done      # generate, gate-count calls, cache <w>.pkl
for w in $W; do $R/probe_u1.py $w 20; done         # f (x1/x100/x1e4), tiers, c_x and 45-station costs -> u1_<w>.summary.json
for w in $W; do (cd build/actual-stock-u1 && python3 rows_stats.py $w); done   # margins
$R/liveness.py figure5_refined                     # the gap-closure parameter merges the wall gap
for w in $W; do $R/last1_refusals.py $w 1; done    # circle-view last-1 refusals, 45 stations
for w in $W; do $R/targeted.py $w; done            # targeted-step Fraction cost; k>0 wall indicator
for w in $W; do $R/costs.py $w 20; done            # cost-only rerun (quiet window)
build/actual-stock-u1/quiet_window.sh 40           # waits up to 20 min for 1-min load < 20, runs costs.py
```

Probe CPU (loaded): `probe_u1` 28 / 20 / 86 s, `last1_refusals` 30 / 26 / 85 s, `targeted` 19 / 9 / 14 s.

## Quiet-window repeat

`quiet_window.sh` fired at 11:25. The 1-min load was 18.0 → 16.2 (fig5), 16.2 → 15.7 (skis) and
15.1 → 25.7 (Monstera, which crossed 20 during its run). Process CPU, warm, min of 5 per circle, on
C_t, every 20th circle (`quiet_window.log`):

| | fig5 | skis | Monstera |
|---|---:|---:|---:|
| c_x: exact station-0 decision, p50 / mean / max | 0.62 / 0.67 / 2.04 ms | 0.23 / 0.42 / 2.03 ms | 0.23 / 0.34 / 1.85 ms |
| 45-station decision, p50 / mean / max | 13.2 / 16.2 / 49.6 ms | 10.8 / 14.4 / 58.9 ms | 9.1 / 11.4 / 64.2 ms |

These match the loaded numbers within a few percent, so the break-even f values above stand: 3.7 % /
7.3 % / 11 % at the quiet mean c_x. The verdict does not change.
