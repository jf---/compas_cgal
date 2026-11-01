# Ratio tracking: same-moment A/B timing against a pinned reference

Status: design, awaiting approval (2026-10-30). Unit RT of the Held-lane correctness + performance stage.

## Intent

The machine runs at load 30–950 from other sessions, and no window at load ≤ 2 has appeared since 2026-10-23. In this project, no timing has ever been measured at load ≤ 2. Absolute budget verdicts therefore never happen. The user has also set them aside for now (2026-10-30: "i dont care about absolute budgets for now").

What the user does want is a **regression signal on the busy machine**: did HEAD get slower than a known commit? Load slows both arms when they run interleaved at the same moment, so the ratio HEAD/reference stays meaningful at any load.

**Success** means that each run of one command appends, for each workload, one entry: HEAD CPU / reference CPU. The entry names exactly which two trees, interpreters and native builds were measured. A deliberately injected slowdown visibly moves the ratio.

## Goals

1. **One command** (`pixi run perf-ratio`) times a fixed set of workloads, interleaved, on the pinned reference worktree and on the HEAD worktree. Each run is a fresh process in its own tree's interpreter, measured as process CPU, min of N per arm.
2. **The ratio is appended** to an append-only history, keyed by `(workload, ratio-vs-<ref short rev>, machine)`. The entry carries:
   - the reference commit and the HEAD tree identity;
   - both arms' native identities and interpreter identities;
   - both arms' CPU runs;
   - the load around every pair.
3. **The history is rendered** to a markdown table that `docs/performance.md` includes.
4. **The defect-injection audit passes:** adding a known amount of CPU to one arm moves the ratio by the expected amount.

## Non-goals

- Absolute budget verdicts and baselines. The quiet window was dropped, and `build/quiet/run.sh` stays unarmed.
- **A pass/fail gate on the ratio in v1.** The first runs are A/A, because HEAD equals the reference today. They give the noise floor. Any threshold is ruled later, from that data.
- Scheduling: no cron. The user or a session runs it.
- The RV41 minors m1/m2/m4 and any refactor of `run_provenance` / `run_native_identity`. This lane only calls them.

## Anti-goals: these would make it a failure even if it "works"

- **Timing both arms through one interpreter or one native build.** That yields a meaningless 1.0, the top risk named at kickoff. Each arm must prove, per run, that its `.py` and `.so` load from its own tree.
- **Wall clock, or CPU that includes the `pixi run` launcher.** Every arm's number must be in-child process CPU.
- **A ratio entering the absolute baseline/check machinery,** where "best = smallest" and "1.10 × best" would turn it into a ratio of ratios and a silent gate.
- A sleep as the injected defect: sleep costs no CPU, so the audit would pass vacuously.

## Constraints

- **Reuse, not a new harness.** Arm plumbing is imported from `tools/pipeline_timing_1_1_base.py`: `base_child_env`, `editable_build_dir_of`, `library_origin_of`, `require_library_in_tree`, `run_in_tree`, `rebuild_base_build`. Identities come from:
  - `tools/run_provenance.py`: `tree_identity`, `interpreter_identity_of`, `machine_identity`, `require_unchanged_tree`;
  - `tools/run_native_identity.py`: `track_native_in`.

  These are imported, not moved. Moving them would be a removal of working code and needs its own go.
- **Python 3.9-compatible.** `mypy --strict` via `types-benchmarks`. One named exception per failure mode. Google docstrings. pathlib only.
- **pyproject edits** (the pixi task and the `types-benchmarks` line) are C++-slot work and land alone, since an edit triggers the rebuild check.
- **No load refusal.** Load is recorded, never gated, so `MachineNotQuietError` is not used here.

## The how, decided

### Arms

- **Reference.** `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-ratio-ref`, a detached worktree at `15fbed70` with its own pixi env and build (user GO 2026-10-30).
  - It must be **clean** at exactly the pinned rev (`require_measured_tree`), before and after the run.
- **HEAD.** A dedicated worktree, `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-ratio-head`: detached, with its own pixi env and build (user GO 2026-10-30).
  - The shared task5 worktree cannot be an arm. Other sessions edit files there and trigger `_editable-rebuild` there, and either one trips `require_unchanged_tree` or the native bracket mid-run.
  - Each run takes `--head-rev`, which defaults to task5's `HEAD` commit. The harness refuses unless ratio-head is clean, then runs `git -C <ratio-head> checkout --detach <rev>`, then rebuilds that tree.
  - Both arms are therefore clean committed trees. The lane measures landed commits only, never uncommitted work. That is the user's choice.
- **The orchestrator is neither arm.** It runs from task5, where the harness code lives. Two consequences:
  - `rebuild_base_build`, which pins the orchestrator's own build dir, does rebuild each arm's tree correctly (tools/pipeline_timing_1_1_base.py:225-236, :371-372).
  - The native probe (`PROBE_SNIPPET`) imports the orchestrator's `tools.run_native_identity` into each arm's interpreter. That is the same probe code in both arms. Its source digest is recorded as `harness_sha256`.
- **Settled builds.** Before timing, each arm's build must pass the native identity probe.
  - A tree first configured with an empty `external/` keeps `verify_*` custom targets that never settle (`NativeBuildStaleError`). This is a CMake finding being root-caused by REF.
  - Arm worktrees are set up in the settled state, and the harness refuses an unsettled one; it never works around it.
  - **`native_sha256` is path-dependent.** Assert `__FILE__` strings embed the worktree's absolute path. Measured 2026-10-31: ratio-ref and ratio-head at the same commit give equal `native_sources_sha256` (09aacb36…) but different `native_sha256` (946ba402… vs 213868a7…).
    - Both are recorded per arm.
    - Neither is ever compared across arms as a sameness test.
    - Matching binaries across worktrees would need `-ffile-prefix-map`; that is a parked CMake finding.
  - Every child runs with `SKBUILD_EDITABLE_SKIP` pinned to its arm's build. Without the pin, the import hook's `cmake --install` touches every `.so`, and the bracket refuses.
- **Symmetric per-arm proofs.**
  - Both arms are identified the same way: `track_native_in` start and finish, `interpreter_identity_of`, and `library_origin_of` + `require_library_in_tree` at start and end.
  - HEAD is not special-cased in-process.
  - `editable_build_dir_of` is checked distinct per arm, and the run refuses if both arms resolve to the same build dir or the same interpreter prefix (`SameArmError`).
  - Both arms' `pixi_lock_sha256` are recorded.
  - Each arm's build is brought current with `rebuild_base_build` before timing.
- **Moving the reference** is a deliberate act: `--ref-worktree` at another commit. That starts a new key series (`ratio-vs-<newrev>`). Old entries stay.

### Workloads (v1 default set)

Every child reports its own in-child process CPU. No `RUSAGE_CHILDREN`, so pixi launcher CPU never counts.

| family | workloads | child | CPU source |
|---|---|---|---|
| polygon generation | figure5_refined, figure8_crossed_skis_boundary, figure8_monstera_boundary_qualified_corners | `python -m tools.held_generation_clock --domain polygon <w>` (worker form) | the tool's own `process_time` clock |
| Ω* generation | same three | `… --domain omega <w>` | same |
| production certificate | Monstera Ω* T(Ω*) | `python -m tools.coverage_certificate_budget_gate --worker` | the tool's in-child `certificate_cpu_seconds` |
| GT pipeline | the 3 cap-120 gate cases | harness-injected clock (below) | `process_time` around the generator call |

**GT clock.** `tools/pipeline_timing --case` times with `perf_counter`, which is wall clock. The lane does not edit that tool, which would make HEAD's tool differ from the reference's. Instead it launches each arm's child with a harness-owned `-c` snippet:
- the snippet imports **that arm's** `tools.pipeline_timing` and replaces the module's `time` with a CPU stand-in, the same technique `held_generation_clock` already uses for generator modules;
- it then calls the arm's `main(["--case", key])`.

The code under test is the arm's own; the clock is the harness's, identical in both arms. If HEAD renames those internals, the run fails loud. The fix then is to move the reference, never to add a fallback.

**Named opt-in set `gt-cap40`.** The 3 cap-40 cases cost 25–481 s each, so not in the default set.

**Determinism.** Within an arm, every run's output digest must agree, or `NondeterministicArmError`. Digests across arms are compared and recorded as `same_output: bool`. A ratio between arms that produce different output is still recorded, flagged.

### Run shape

- **N = 5 pairs per workload** (`--pairs`).
- **The order alternates per pair** (AB, BA, AB, BA, AB) so slow load drift does not systematically favour one arm.
- The 1-minute load is recorded before and after every pair.
- Per arm: all runs, min, median.
- **ratio = min(HEAD) / min(ref)** is the headline. The **median of per-pair ratios** is recorded beside it as the spread witness.
- Budget estimate for the default set at load ~200: about 30 min. Measure it in the first run and record it in the docs.

### Record and history

- **One run writes one record:** `benchmarks/results/perf_ratio_session.json`, schema `perf-ratio/1`. It contains:
  - the reference commit;
  - per arm: tree identity, native identity, interpreter identity, lock sha and library origins;
  - per workload: runs, pair order, loads, min/median, ratio, pair-ratio median, same_output;
  - machine, start and end UTC, harness source digest.
- **Separate history:** `benchmarks/results/perf_ratio_history.jsonl` (append-only, deduped by content hash) plus a rendered `perf_ratio_history.md`.
  - It does **not** mix into `perf_history.jsonl`. Ratios have different semantics: no baseline, no check, unit "×" not "s". Keeping them apart means `perf_history_baseline` / `check` / `render` stay untouched and can never treat a ratio as seconds.
  - Store and JSON helpers from `perf_history_store` / `perf_history_json` are reused where they fit.
  - `perf_ratio_session.json` is registered in perf history's `NON_TIMING_RECORDS`, so the "every results JSON has an owner" test stays green. It is owned by the ratio history, not the absolute one.
- **Render.** One row per key, newest first: date, HEAD rev (+dirty), ratio, pair-ratio median, same_output, load range.
  - Header text: "a ratio is a same-moment comparison at any load, never a budget verdict."
  - It is included in `docs/performance.md`, after the performance history section, with a short paragraph stating these semantics and how the lane differs from the D2 quiet-only CPU gates.
  - Byte-for-byte guarded like the other included records.

### Audit (exit gate)

- **Injected CPU.** `--inject-cpu-seconds ARM=S` is a harness flag that busy-spins S seconds of process CPU inside one arm's child, after the workload, within the measured span.
  - Records made with it carry `injected` and are refused by the history ingest (`InjectedRecordError`).
  - Injection is supported **only where the harness owns the clock**, which is the GT snippet. Spinning inside the polygon, Ω* or certificate tools' own clocks would mean patching arm code, and adding S to a reported number would make the audit vacuous. Injection into any other family is refused (`InjectionUnsupportedError`).
  - Acceptance: on one real GT cap-120 case, injecting S into HEAD moves the ratio up by about S/min(ref). The tolerance is stated from the measured A/A spread, not picked. Injecting into the reference moves it down.
- **Unit tests use fake arms.** Two tmp git trees with stub children cover:
  - same-interpreter / same-build refusal;
  - a foreign library refused;
  - a tree changed mid-run refused;
  - nondeterminism refused;
  - injected records refused at ingest.

  Failing tests first.

### Code layout (one responsibility per file)

| file | responsibility |
|---|---|
| `tools/perf_ratio.py` | CLI and orchestration |
| `tools/perf_ratio_arms.py` | arm setup and proofs, using the imported 1.1-base helpers |
| `tools/perf_ratio_workloads.py` | workload registry, child commands, result parsing, GT clock snippet |
| `tools/perf_ratio_record.py` | `perf-ratio/1` types, `build()` validation, named errors |
| `tools/perf_ratio_history.py` | ingest and render |

Plus tests under `tests/tools/` and a pixi task `perf-ratio` (C++ slot), modelled on `pipeline-timing-1-1-base` with no `depends-on`.

## Left to the builder

- Exact snippet text, parsing and error class names (one per failure mode).
- Whether `ClockSummary`-like types are reused or a ratio-specific `Ratio` NewType is introduced. Prefer the NewType: the unit lives in the type.
- The table column order.

## Calls made without the user (ledgered as rulings)

1. A separate ratio history instead of a stage inside `perf_history.jsonl`.
2. GT timed by a harness-injected CPU clock in each arm, instead of editing `pipeline_timing` (which would force moving the reference).
3. The cap-40 GT cases are opt-in, not in the default set.
4. No ratio gate in v1; the threshold comes later, from A/A data.
5. Alternating pair order.
6. Helpers imported from the 1.1-base wrapper, not moved.
7. CPU injection only in the GT family, where the harness owns the clock.
8. The orchestrator runs from task5 and is neither arm. Each run checks out `--head-rev` in ratio-head.
