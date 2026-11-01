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
- After each promotion, margins are re-calibrated: an A/A at the new reference.

### Amendment 2026-11-01: the anchor chain (I1, user ruling (a))

The bullet "the reference never moves to a slower commit; small regressions cannot creep" was false. Promotion needs only one improved workload, so every non-improved workload may sit anywhere up to `1 + m_w` against the old reference, and after promotion it is measured against that slower commit. Without a fix the drift is up to `m_w` per promotion, unbounded over many.

**Rule.** For each workload, the registry carries an anchor ratio `A_w`: the current reference's CPU time divided by the time of that workload's anchor commit. The gate checks

`A_w · ρ_w ≤ 1 + m_w` for every workload.

A fresh registry has `A_w = 1` for every workload. On promotion, with `ρ_w` taken from the promoted session after extension:

- `A'_w = A_w · ρ_w` when `A_w · ρ_w ≥ 1 − m_w`. The workload did not certifiably beat its anchor, so the anchor stays put and the drift is carried.
- `A'_w = 1` otherwise. The new reference is a certified improvement on the anchor, so it becomes the anchor. This locks the improvement in, and later commits cannot spend it back.

A workload that is new to the registry starts at `A_w = 1` at the current reference. A workload that leaves the registry drops its anchor. A harness rebase keeps the chain: `ρ_w` is a same-harness ratio, so links measured under different harnesses still compose.

**IMPROVED and promotion.** A workload is `IMPROVED` when `A_w · ρ_w < 1 − m_w`. Promotion still needs a passing gate session (every workload at or under `1 + m_w` against its anchor) and at least one `IMPROVED` workload.

**Content addressing.** The chain is re-derivable, like the margins. Each promotion appends a link `{from_ref, to_ref, record path, record_sha256}`, and the promoted session record is committed under `benchmarks/results/perf_ratio_promotions/<from8>-<to8>/perf_ratio_session.json`. On load, `A_w` is recomputed from the links with the rule above. A registry whose stored `A_w` disagree with that recomputation, or whose link records are missing or have a mismatched digest, is refused (`AnchorChainMismatchError`). Records are read as exact doubles, and the products are computed in `fractions.Fraction` so that recomputation is deterministic. The stored values are reporting copies.

**What the guarantee costs (accepted with (a)).** Each link is one measured ratio, so it carries one session's noise, at most `m_w / SAFETY` under the calibrated noise model. After `n` carried links, the gate holds the true slowdown against the anchor to `m_w + n · m_w / SAFETY` in the worst case. That bound is linear in the number of promotions; the creep it replaces was unbounded. More pairs in a gate session cannot shrink this error, because it lives in past sessions. A certified improvement resets the anchor, which drops that workload's accumulated error. Deferred, with no ruling yet: extending every workload of a promotion session to `MAX_PAIRS`, which would shrink each link's noise.

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
