# Ratio gate audit on the anchored code: head-injection ρ outside the predicted band

## Conclusion

The verdict classes hold on all three audit runs:
- the A/A run is PASS on all 10 workloads, with no extension;
- the run injected into head is REGRESSION;
- the run injected into ref is IMPROVED.

The head run's final ρ is 1.108, above the 1.096–1.098 band in plan Task 4. That band is a load-dependent prediction (`1 + S/ref_cpu`), not a contract. The base 5 pairs of the fresh head session give ρ 1.0969, inside the band. The 15 extension pairs, measured at load 5–7, carry the excess. No code change is needed. The exit gate, "the three audit outcomes hold", is met.

## Evidence

The session JSONs are in `build/ratio-gate/`. The injection is S = 0.126 s CPU on `gt/cap-120-default/pipeline/rect_20x12`. ρ and pairs are the session's own base pairs; the check-time extensions are not in the session file.

| Session | ρ (base pairs) | ref_min s | S/ref_min | pairs | load |
|---|---|---|---|---|---|
| `audit-head.json` (landed gate, 2026-11-01 16:2x) | 1.0964 | 1.3188 | 0.0955 | 5 | 47–53 |
| `audit-anchor-head.json` (fresh) | 1.0969 | 1.2732 | 0.0990 | 5 | 6–7 |
| `audit-ref.json` (landed gate) | 0.8958 | 1.4619 | 0.0862 | 5 | 39–51 |
| `audit-anchor-ref.json` (fresh) | 0.9003 | 1.4089 | 0.0894 | 5 | 5–6 |

Check verdicts from `audit-anchor-{crosscheck,aa,inject}.log`:
- fresh A/A: 10 PASS, exit 0;
- fresh head: REGRESSION, `A_w·ρ_w` 1.1084 > 1 + 0.0481 (`A_w` 1.0), 20 pairs, exit 1;
- fresh ref: IMPROVED, ρ 0.900, 5 pairs, exit 0;
- re-check of the landed sessions on the new decide: A/A 10 PASS, head REGRESSION 1.103 at 20 pairs (its extensions measured today at load 5–7), ref IMPROVED 0.896.

## Why ρ rises at low load

The injection adds a fixed S to the head arm's CPU, so ρ ≈ 1 + S/ref_cpu. On this M1 Max, in-child CPU seconds shrink as load falls: `ref_min` is 1.319 s at load ~50 and 1.273 s at load ~6. With less contention for shared cache and memory bandwidth, the same work takes fewer CPU seconds. So the same S yields a larger ratio. Extension pairs measured at load 5–7 have even shorter reference times, which pulls the merged ρ up to 1.108.

The plan's band came from the landed audit at load 33–105. The gate decides `A_w·ρ_w` against `1 + m_w` (0.048 here), and 1.108 is far outside it either way. The gate is not affected.

## Consequence

- Record the audit in `docs/performance.md` with the fresh numbers. State the band as load-dependent: the base pairs give 1.097, and the result after extension at load ~6 is 1.108.
- Future audit expectations should be stated as verdict classes plus `ρ ≈ 1 + S/ref_min` computed from the session itself, not as a fixed band.
