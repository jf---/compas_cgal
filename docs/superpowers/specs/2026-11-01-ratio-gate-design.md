# Ratio regression gate (RT v2)

Status: approved by user "go", 2026-11-01; open question 1 resolved as refuse (the strict default). It supersedes the non-goal "a pass/fail gate on the ratio in v1" in `2026-10-30-ratio-tracking-design.md`. That spec already deferred the threshold to A/A data; this is that ruling.

User, 2026-11-01: "register these values in the benchmark such that we can never regress on these. performance is key. however, absolute performance values may vary on fluctuating hardware. CI could be busy (regression) or upgraded (so suite should converge faster) you need to have an answer on that relative performance."

## Decision

The gated quantity is the existing per-workload ratio `ρ_w = min_i head_cpu_i / min_i ref_cpu_i`. Its inputs are:
- in-child process CPU only;
- reference and HEAD interleaved AB/BA on the same machine, in the same session.

Absolute seconds are recorded and never gated.

- **Registered:** for each workload, a noise margin `m_w`, plus the reference commit.
- **Not registered:** times.

A HEAD passes when `ρ_w ≤ 1 + m_w` for every workload.

## The hardware answer

| Situation | What happens to `ρ_w` | Why |
|---|---|---|
| CI busy (load) | No bias. Noise rises. | Both arms run at the same moment, so a slower machine scales both numerators. AB/BA interleaving cancels drift across pairs. Min-of-N rejects bursts inside a pair. |
| CI upgraded / faster runner | Unchanged. The run gets cheaper in wall time. | Same-moment comparison: the speed factor cancels. |
| Different CPU model | Unchanged in expectation. Noise may differ. | `cpu_model` is recorded per session; margins are re-measured per machine class (see below). |

Busy machines are handled by **measuring more, not by widening the threshold**. If any workload lands in `(1 + m_w, ∞)`, the gate extends that workload with more pairs, up to `MAX_PAIRS`, and decides again. Only a result that stays above the margin after extension is a regression. This is a sequential test. A busy runner costs time, not false reds.

## Margins: where `m_w` comes from

Margins come from an A/A calibration: K sessions with both arms at the reference commit.
- K = 5 for the first calibration: the Task 7 record plus 4 more, running 2026-11-01.
- `m_w = SAFETY · max_k |ρ_{w,k} − 1|`, with `SAFETY = 2`. It is a named constant; its derivation comment says it doubles the worst observed A/A deviation.
- Margins are content-addressed. The registry stores the `record_sha256` of every calibration session, plus K, SAFETY and the calibration `cpu_model`. A registry whose sessions do not re-derive its margins is refused on load.
- One A/A session's pair spread `τ_w` is **not** the margin. It measures pair noise, not the noise of the gated `ρ_w`. τ from the first A/A run: 0.007–0.149. ρ deviation in the same run: ≤ 0.022.

## "Never regress": the reference is a ratchet

- The reference commit lives in the registry `benchmarks/results/perf_ratio_gate.json`, not in code. `PINNED_REF_COMMIT` becomes the registry's value.
- The reference advances only through `perf-ratio-promote`. Promotion needs:
  - a passing gate session on a clean, committed HEAD;
  - at least one workload with `ρ_w < 1 − m_w` (a real improvement).

  It then re-pins ratio-ref, re-settles it, and writes the new registry.
- The reference never moves to a slower commit. Small regressions cannot creep: each one is measured against the same fixed reference, so the sum of all of them must stay within `m_w`.
- After each promotion, margins are re-calibrated: an A/A at the new reference.

## Verdict and failure model

`perf-ratio-check` takes a session record and the registry. It reports one verdict per workload and exits non-zero on any regression. The verdicts:
- `PASS`;
- `IMPROVED`, which proposes a promotion;
- `REGRESSION`, raised as `RatioRegressionError`, naming the workload, `ρ_w`, `m_w` and the pairs.

The named refusals:
- `UncalibratedMachineError`: the session's `cpu_model` has no calibration in the registry (open question 1);
- `RegistryDigestMismatchError`;
- `ReferenceMismatchError`: the session's reference is not the registry's;
- the existing refusals: an injected record, `same_output` false, and arm identity errors.

## Proof the gate catches regressions (audit, extends plan Task 8)

- Injection sized to `2·m_w` on a GT workload gives `REGRESSION`.
- The same injection on the reference arm gives `IMPROVED`.
- An A/A session gives `PASS` on every workload.

All three are exit gates.

## Cost

- One gate session is about 16 min wall, measured at load 3–64. Sequential extension adds time only on red-zone workloads.
- CI wiring is a separate unit. The CI configuration gate applies: the analysis is shown before any workflow edit.

## Open questions for the user

1. Resolved 2026-11-01: refuse. A session from a `cpu_model` with no calibration: refuse (`UncalibratedMachineError`, run `perf-ratio-calibrate` once per machine class), or gate with the widest registered margin? The first is the strict default; the second never blocks a new CI runner.
