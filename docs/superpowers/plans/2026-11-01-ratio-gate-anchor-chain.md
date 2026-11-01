# Ratio gate: anchor chain (I1 ruling (a)) — implementation plan

Authority: `docs/superpowers/specs/2026-11-01-ratio-gate-design.md` §"Amendment 2026-11-01: the anchor chain". Base: the landed gate (`9ae1cdd2`) plus Task 6 (early harness refusal in `calibrate --sessions`). Every task: TDD (red first, recorded), R6 defect injection per new refusal / decision branch, implementers never commit, the orchestrator re-runs every number. Never edit HARNESS_SOURCES.

**Goal.** The gate decides `A_w · ρ_w` against `1 ± m_w`; promotion carries `A_w` forward (or resets it on a certified improvement) through a content-addressed chain of committed records, so non-improved workloads cannot creep across promotions.

## Design decisions fixed by this plan

- **Exact decision.** `A_w` and `ρ_w` are compared as `fractions.Fraction` products of the recorded doubles: `Fraction(A_w) · Fraction(ρ_w)` against `1 ± Fraction(m_w)`. `A_w` is never read back from the stored reporting float; it is re-derived from the chain on load. The plain gate (empty chain, `A_w = 1`) changes only at sub-ulp boundaries, which the audit re-run confirms.
- **A link** (`AnchorLink`) is `{from_ref, to_ref, sessions, calibration}`:
  - `sessions`: the promoted gate session record followed by its extension session records, in merge order, each a `CitedRecord` (path + `record_sha256`) under `benchmarks/results/perf_ratio_promotions/<from8>-<to8>/`: `gate/perf_ratio_session.json`, `ext-<n>/perf_ratio_session.json`.
  - `calibration`: the `Calibration` that supplied `m_w` for the promoted session's `cpu_model` (its cited A/A records stay committed under `perf_ratio_calibration/<from8>/…`). Its margins re-derive on load exactly as a live calibration's do.
  - On load, `ρ_w` per link is re-derived by replaying `merge_extension` over the cited sessions; `A_w` by the reset rule over the links in order. Missing record, digest mismatch, a link whose `from_ref` is not the previous link's `to_ref` (or the first link's is not the chain's anchor ref), a last `to_ref` that is not `ref_commit`, or a stored `A_w` that disagrees → `AnchorChainMismatchError` (one named exception per failure mode where the existing registry errors do not already name it; reuse `_cited_path` for root containment).
- **Registry schema** `perf-ratio-gate/2`: adds `anchor_ref` (the chain's first reference), `anchor_links` (list of links), `anchors` (workload → reporting float). Schema `/1` registries are refused with a message naming the migration (`calibrate --records` over the existing A/A records rebuilds at `/2` with an empty chain). No in-code migration path.
- **Workload set changes**: a rebase that adds a workload gives it `A_w = 1` at the current reference (the chain records it as starting at the link where it first appears); a removed workload's anchor is dropped. `--rebase-harness` keeps the chain.

## Task 1 — registry: the chain type, derivation, verification, schema /2

**Files:** `tools/perf_ratio_gate_registry.py`, `tests/tools/test_perf_ratio_gate_registry.py`, `tests/tools/perf_ratio_sessions.py` (synthetic chain helpers).

- [ ] `AnchorLink` frozen dataclass; `GateRegistry` gains `anchor_ref`, `anchor_links`, and `anchors_exact() -> Mapping[str, Fraction]` (derived), `anchors` (reporting floats in `dumps`).
- [ ] `derive_anchors(links, workloads, *, root, merge) -> dict[str, Fraction]`: pure, the reset rule `A' = A·ρ if A·ρ ≥ 1 − m else 1`. `merge` is injected (`merge_extension` lives in the verdict module; keep the dependency direction registry ← verdict by passing it in, or move the replay into a small `tools/perf_ratio_gate_chain.py` if that reads cleaner — one responsibility per file; state the choice).
- [ ] `build_registry(..., anchor_ref=None, anchor_links=())`; `load_registry` verifies the chain; `REGISTRY_SCHEMA = "perf-ratio-gate/2"`.
- [ ] Tests: empty chain → all 1; one carried link (A·ρ = 1.02, m = 0.03) → 1.02; one reset link (A·ρ = 0.95, m = 0.03) → 1; carry-then-reset; new workload mid-chain → 1; each refusal (missing record, digest, broken from/to continuity, wrong last to_ref, stored anchor disagrees, schema /1). R6 for every refusal and for both reset-rule branches (inject `A' = A·ρ` always, and `A' = 1` always: each must turn a named test red).

## Task 2 — verdict: decide on `A_w · ρ_w`

**Files:** `tools/perf_ratio_gate_verdict.py`, `tests/tools/test_perf_ratio_gate_verdict.py`.

- [ ] `decide(run, margin, anchor: Fraction)`: REGRESSION when `anchor·ρ > 1 + m`, IMPROVED when `anchor·ρ < 1 − m`, exact in `Fraction`.
- [ ] `WorkloadVerdict` carries `anchor` and the cumulative ratio (reporting floats); `gate` passes each workload's anchor, extends on red exactly as now (the red test uses the anchored product), and keeps every extension session on the report (`GateReport.extensions: tuple[PerfRatioSession, ...]`, merge order) so promote can persist them.
- [ ] Tests: anchor 1.03 turns ρ = 1.0 with m = 0.02 into REGRESSION, then an extension that brings ρ to 0.98 passes; anchor 0.96 with ρ = 1.0, m = 0.03 → IMPROVED; an exact-boundary case where the float product rounds across `1 + m` but the Fraction product does not (construct the doubles explicitly, assert the exact verdict). R6 on the anchor multiplication and on the exact comparison.

## Task 3 — promote: write the link

**Files:** `tools/perf_ratio_gate.py`, `tests/tools/test_perf_ratio_gate.py`.

- [ ] After every existing refusal and before `checkout`: compute the new anchors from the report, write the gate session + extension sessions under `perf_ratio_promotions/<from8>-<to8>/` (staged under `GATE_SCRATCH`, one rename, as `calibrate_sessions` does; refuse `PromotionRecordsExistError` if the target exists), build the `/2` registry with the appended link and the promoting machine's calibration.
- [ ] If the checkout lands elsewhere, the registry is unchanged and the staged/renamed records are removed (or the target is left and the rerun re-uses it after a digest check — state the choice; the rerun must finish the promotion).
- [ ] Report lines print each workload's `A_w → A'_w` and which branch (carried / reset).
- [ ] Tests (fake measure, fake checkout, `root=tmp_path`): a carried workload's anchor equals `A·ρ`; an improved workload resets to 1; the written registry loads and re-derives; a second promotion composes; a non-improved workload that drifts `+m_w` at each of two promotions is REGRESSION against the original anchor at the third session (the creep I1 named, now caught). R6 on the carry, the reset, and the persisted extension list.

## Task 4 — registry regeneration, docs, audit (C++ slot, alone)

- [ ] Regenerate `benchmarks/results/perf_ratio_gate.json` at `/2`: `pixi run --no-install perf-ratio-calibrate --records benchmarks/results/perf_ratio_calibration/15fbed70/apple-m1-max/aa-{2,3,4,5,6}/perf_ratio_session.json` against a removed-then-rebuilt registry (orchestrator confirms the old file is moved aside under `build/ratio-gate/`, not deleted). Expected margins unchanged to 3 dp; `anchor_ref = ref_commit`, empty chain.
- [ ] `docs/performance.md` "### The ratio regression gate": replace the ratchet paragraph with the anchor rule, the reset branch, the cost bound `m_w + n·m_w/SAFETY`, and the promotion record layout.
- [ ] Gates: gate tests + lane tests (`-n 2`) + `test_included_records`; ruff; `types-benchmarks` Success; strict mkdocs (only the known `segment_site_mat.md` warning).
- [ ] **Audit re-run** (decide changed): A/A → PASS on all 10 with no extension; `+S` into HEAD on rect_20x12 → REGRESSION (ρ ≈ 1.096–1.098); into ref → IMPROVED. Logs under `build/ratio-gate/audit-anchor-*`. Load recorded.
- [ ] Commit batch (ask first): the three gate modules (+ chain module if split), their tests, `perf_ratio_sessions.py`, the registry, `docs/performance.md`, the spec amendment, this plan.

## Exit gates

- The three audit outcomes hold on the anchored code.
- The two-promotion creep test is red on the landed gate (`9ae1cdd2`) and green after Task 3.
- A `/1` registry is refused with the migration message; the regenerated `/2` registry re-derives on load.
