# Ratio Regression Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pixi run perf-ratio-check <session>` gives one verdict per workload (PASS / IMPROVED / REGRESSION) against a content-addressed registry of A/A-calibrated margins. It extends red workloads with more pairs before it calls a regression. `perf-ratio-calibrate` and `perf-ratio-promote` own the registry.

**Architecture:** Three new modules beside the landed ratio lane, none of them a harness source:
- `perf_ratio_gate_registry` holds the registry and re-derives its margins from the cited A/A records;
- `perf_ratio_gate_verdict` holds the pure verdict and the sequential-extension loop over an injected `Measure`;
- `perf_ratio_gate` is the CLI, whose `Measure` runs `python -m tools.perf_ratio` as a subprocess.

No file in `HARNESS_SOURCES` (`tools/perf_ratio.py:70-85`) is edited. An edit there changes `harness_sha256` and voids the calibration.

**Tech Stack:** Python 3.9-compatible, stdlib only, pytest + xdist, mypy --strict via `types-benchmarks`, ruff, pixi.

**Spec:** `docs/superpowers/specs/2026-11-01-ratio-gate-design.md`. Interfaces consumed: `tools/perf_ratio_record.py`, `tools/perf_ratio_history.py`, `tools/perf_ratio_arms.py` (landed `48b8e632`).

## Global Constraints

- Work only in `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation`. Python via `pixi run --no-install build/pyrun.sh ...`; never bare python/pytest/pip.
- Python 3.9-compatible: `from __future__ import annotations`, `typing.Optional`, no `match`, no runtime `X | Y`. pathlib only. Google docstrings. Named UPPERCASE constants, each with a one-line comment giving units and derivation. One named exception per failure mode, never bare `ValueError`.
- **Never edit a `HARNESS_SOURCES` file.** That covers `tools/perf_ratio{,_record,_workloads,_arms,_history}.py`, `pipeline_timing_1_1_base.py`, `run_provenance.py`, `run_native_identity.py`, `held_native_identity.py`, `held_source_closure.py` and `benchmarks/{__init__,errors,performance_record}.py`. Import from them only.
- Gate on `ρ_w ≤ 1 + m_w` only (`WorkloadRatio.ratio`). Absolute seconds are never gated. A busy machine means more pairs, never a wider margin.
- `SAFETY = 2` (spec). `CALIBRATION_SESSIONS = 5` (spec K). A `cpu_model` without a calibration is refused (`UncalibratedMachineError`, user 2026-11-01).
- The reference commit is the registry's `ref_commit`. The CLI always passes `--ref-rev <registry ref_commit>` to the harness. `PINNED_REF_COMMIT` stays the harness default and is not edited: editing it is a harness edit. `ReferenceMismatchError` makes a stale default loud.
- Unit tests use synthetic sessions and a fake `Measure` only. No real timing in tests.
- Only Task 4 (C++ slot) edits `pyproject.toml`, alone.
- Targeted tests: `pixi run --no-install build/pyrun.sh -m pytest <file> -n 2 -q`. Lint: `pixi run --no-install ruff check <files>`. Types (after Task 4): `pixi run --no-install types-benchmarks`.
- Implementers never commit and never dispatch. TDD: the failing test and its red run come before the implementation. Every refusal gets a defect-injection-proven test (R6): delete the check and show the test go red.

## Review Focus

1. **The harness is edited after calibration.** Expected: every session carries a new `harness_sha256` and is refused (`HarnessMismatchError`), never decided against stale margins. Pinned by `test_a_session_from_another_harness_is_refused` (Task 2).
2. **An extension batch measured different code than the base session.** Examples: ratio-head moved, a rebuild changed `native_sha256`, or a different injection. Expected: `ExtensionMismatchError`, never a merged min over two builds. Pinned by `test_an_extension_from_another_build_is_refused` (Task 2).
3. **A calibration record's bytes change under the same path** (re-run, hand edit). Expected: the registry is refused on load (`CalibrationRecordDigestError`). Pinned by `test_a_regenerated_calibration_record_is_refused` (Task 1).
4. **A partial session** (`--workload` subset) is checked as a gate. Expected: `IncompleteSessionError`, since a one-workload PASS would pass the gate. Only `--audit` decides subsets. Pinned by `test_a_partial_session_is_refused_outside_audit` (Task 2).
5. **Promotion from an audit (injected) session, a dirty HEAD, or a session with a REGRESSION or no IMPROVED workload.** Expected: refused, registry bytes unchanged. Pinned by `test_promotion_refusals_leave_the_registry_alone` (Task 3).

## Spec amendment: `same_output` (user "go", 2026-11-01)

Spec §Verdict lists "`same_output` false" among the refusals. Taken literally, the ratchet deadlocks. Any HEAD that changes a toolpath's output is refused forever, so it can never be promoted, and every later HEAD carrying that change is refused too. **Ruled:** decide `ρ_w` regardless, and report `changed_output` per workload beside the verdict. `same_output` stays a hard requirement for calibration sessions (A/A must agree) in Task 1.

---

### Task 1: Registry, margins re-derived from cited A/A records

**Files:**
- Create: `tools/perf_ratio_gate_registry.py`, `tests/tools/test_perf_ratio_gate_registry.py`, `tests/tools/perf_ratio_sessions.py` (shared synthetic-session builder)
- Modify: `tools/perf_history_adapters.py:575-587`: add `"perf_ratio_gate.json": "ratio gate registry: A/A margins and the reference commit; owned by tools.perf_ratio_gate_registry"` to `NON_TIMING_RECORDS`.

**Interfaces (produces):**

```python
REGISTRY_SCHEMA = "perf-ratio-gate/1"
REGISTRY_PATH = REPOSITORY / "benchmarks" / "results" / "perf_ratio_gate.json"
CALIBRATION_ROOT = REPOSITORY / "benchmarks" / "results" / "perf_ratio_calibration"
SAFETY = 2.0             # doubles the worst observed A/A deviation (spec §Margins)
CALIBRATION_SESSIONS = 5  # K, spec §Margins

Margin = NewType("Margin", float)

@dataclass(frozen=True) class CitedRecord: path: str; record_sha256: str      # path relative to REPOSITORY
@dataclass(frozen=True) class Calibration: cpu_model: str; sessions: tuple[CitedRecord, ...]; margins: Mapping[str, Margin]
@dataclass(frozen=True) class GateRegistry:
    ref_commit: str; harness_sha256: str; safety: float; workloads: tuple[str, ...]  # sorted
    calibrations: tuple[Calibration, ...]
    def margins_for(self, cpu_model: str) -> Mapping[str, Margin]   # raises UncalibratedMachineError
    @property registry_sha256 -> str          # sha256 of canonical JSON content (sort_keys, (",", ":")), as perf_ratio_history._digest
    def dumps(self) -> str                    # indent=2, sort_keys, trailing newline, includes registry_sha256

def derive_margins(sessions: Sequence[PerfRatioSession], safety: float = SAFETY) -> dict[str, Margin]
    # m_w = safety * max_k |sessions[k].workload(w).ratio - 1|
def calibration_of(record_paths: Sequence[Path], *, ref_commit: str, harness_sha256: str) -> Calibration
def build_registry(*, ref_commit: str, harness_sha256: str, calibrations: Sequence[Calibration]) -> GateRegistry
def load_registry(path: Path = REGISTRY_PATH) -> GateRegistry
def write_registry(registry: GateRegistry, path: Path = REGISTRY_PATH) -> None
def calibration_dir(ref_commit: str, cpu_model: str) -> Path   # CALIBRATION_ROOT / ref[:REF_SHORT_CHARS] / slug(cpu_model)
```

Errors: `MalformedRegistryError`, `RegistryDigestMismatchError` (stored `registry_sha256` ≠ content), `CalibrationRecordDigestError` (cited file's `sha256_file` ≠ `record_sha256`), `MarginMismatchError` (stored margin ≠ re-derived, compared with `==`, since floats round-trip through JSON), `InvalidCalibrationSessionError`, `TooFewCalibrationSessionsError`, `DuplicateCalibrationError` (two calibrations for one `cpu_model`), `UncalibratedMachineError`.

`InvalidCalibrationSessionError` covers a cited session that is any of:
- injected;
- an arm not at `ref_commit`, or dirty;
- `harness_sha256` ≠ the registry's;
- `cpu_model` ≠ the calibration's;
- `same_output` false on any workload;
- a workload set ≠ the registry's.

`calibrated_slug(cpu_model)`: lowercase, runs of non-`[a-z0-9]` → `-`, stripped (`"Apple M1 Max"` → `apple-m1-max`).

- [ ] **Step 1: Write `tests/tools/perf_ratio_sessions.py`.** It holds `synthetic_session(ratios: Mapping[str, Sequence[tuple[float, float]]], *, ref_commit=REF, head_rev=REF, cpu_model="Apple M1 Max", harness="h"*64, injected=None, head_digest="d", ref_digest="d", native="n"*64) -> PerfRatioSession`. Each `(ref_cpu, head_cpu)` tuple is one pair, built through `PairRun.build`, with orders alternating per `tools.perf_ratio.pair_order`. Also `write_record(session, path) -> Path`, which writes `session.dumps()`.
- [ ] **Step 2: Write the failing tests:**
  - `test_margins_are_safety_times_worst_aa_deviation`: 5 sessions with ratios 1.00, 0.98, 1.01, 1.00, 1.005 on `w` → `derive_margins(...)["w"] == 2 * 0.02` (up to float representation, compute expected from the same pair CPUs).
  - `test_a_registry_round_trips_and_re_derives_on_load`: write 5 records in `tmp_path`, then `calibration_of` → `build_registry` → `write_registry` → `load_registry`. The result is equal and its `registry_sha256` is stable.
  - `test_an_edited_registry_is_refused`: change one margin digit in the file → `RegistryDigestMismatchError`.
  - `test_a_re_digested_registry_with_a_wrong_margin_is_refused`: edit the margin and recompute `registry_sha256` → `MarginMismatchError`.
  - `test_a_regenerated_calibration_record_is_refused` (Review Focus 3): rewrite one cited record with different pairs → `CalibrationRecordDigestError`.
  - `test_an_invalid_calibration_session_is_refused`: one session per reason (injected; head arm at another rev; harness differs; cpu differs; `same_output` false; missing workload) → `InvalidCalibrationSessionError`.
  - `test_fewer_than_k_sessions_are_refused` (4 records) → `TooFewCalibrationSessionsError`.
  - `test_an_uncalibrated_cpu_model_is_refused`: `margins_for("Intel Xeon")` → `UncalibratedMachineError`.
  - `test_two_calibrations_for_one_cpu_are_refused` → `DuplicateCalibrationError`.
  - `test_calibration_dir_names_ref_and_cpu`: `calibration_dir(REF, "Apple M1 Max") == CALIBRATION_ROOT / "15fbed70" / "apple-m1-max"`.
- [ ] **Step 3: Run red.** `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_gate_registry.py -n 2 -q`. Expected: collection error / failures, since the module is absent.
- [ ] **Step 4: Implement** `tools/perf_ratio_gate_registry.py` and the `NON_TIMING_RECORDS` line. `load_registry` checks in this order: schema/shape → self-digest → for each calibration: K, each cited digest, each session's validity, then margin re-derivation.
- [ ] **Step 5: Run green**, then the R6 injection: for each refusal, delete its check in a copy (`PYTHONPATH=<copy>:<worktree>`) and show its test red. Also run `tests/tools/test_perf_history.py -n 2 -q`. Expected: only the known `test_rebackfilling_unchanged_families_appends_nothing_for_them` fails, and `test_every_committed_record_is_classified` still passes.
- [ ] **Step 6: Hand to orchestrator.** pathspec = the four files above.

### Task 2: Verdict and sequential extension

**Files:** Create `tools/perf_ratio_gate_verdict.py`, `tests/tools/test_perf_ratio_gate_verdict.py`.

**Interfaces:**
- Consumes: Task 1 `GateRegistry`, `Margin`; Task 1 test helper `synthetic_session`; `PerfRatioSession`, `WorkloadRatio`, `Injection`, `ArmName` (record).
- Produces:

```python
EXTENSION_PAIRS = 5   # pairs per extension batch: one default batch (perf_ratio.DEFAULT_PAIRS)
MAX_PAIRS = 20        # cap per workload: 4x the calibrated N, three extensions; bounds a red workload's extra wall
                      # at 15 pairs (omega monstera, the costliest, ~27 s/pair on M1 Max -> ~7 min)

class Verdict(str, Enum): PASS = "PASS"; IMPROVED = "IMPROVED"; REGRESSION = "REGRESSION"

@dataclass(frozen=True) class ExtensionRequest:
    ref_commit: str; head_rev: str; workloads: tuple[str, ...]; pairs: int; injected: Optional[Injection]
Measure = Callable[[ExtensionRequest], PerfRatioSession]

@dataclass(frozen=True) class WorkloadVerdict:
    workload: str; verdict: Verdict; ratio: Ratio; margin: Margin; pairs: int; changed_output: bool
@dataclass(frozen=True) class GateReport:
    verdicts: tuple[WorkloadVerdict, ...]; audit: bool
    @property regressions -> tuple[WorkloadVerdict, ...]
    @property improvements -> tuple[WorkloadVerdict, ...]

def decide(run: WorkloadRatio, margin: Margin) -> Verdict
    # ratio > 1 + margin -> REGRESSION; ratio < 1 - margin -> IMPROVED; else PASS (both bounds inclusive to PASS)
def require_gateable(session: PerfRatioSession, registry: GateRegistry, *, audit: bool) -> Mapping[str, Margin]
def merge_extension(base: WorkloadRatio, base_session: PerfRatioSession, extension: PerfRatioSession) -> WorkloadRatio
def gate(session: PerfRatioSession, registry: GateRegistry, measure: Measure, *, audit: bool = False) -> GateReport
```

Loop of `gate`, per workload in session order:
1. Decide.
2. While the verdict is REGRESSION and `pairs < MAX_PAIRS`: request `min(EXTENSION_PAIRS, MAX_PAIRS − pairs)` pairs for that workload (`ExtensionRequest` with the base session's ref, HEAD rev and injection), `merge_extension` the result, decide again.
3. The merged ratio is min over all pairs (`WorkloadRatio` with concatenated `pairs`; the digests come from the base).

`changed_output = not run.same_output` (spec amendment, header).

Errors in `require_gateable`, checked in this order:
- `InjectedSessionError` (injected and not `audit`);
- `ReferenceMismatchError` (`session.ref_commit != registry.ref_commit`);
- `HarnessMismatchError`;
- `UncalibratedMachineError` (from `margins_for`);
- `UnregisteredWorkloadError` (a session workload with no margin);
- `IncompleteSessionError` (not `audit` and the workload set ≠ the registry's).

`merge_extension` raises `ExtensionMismatchError` when the extension's arms (`ArmIdentity` equality per arm), `harness_sha256`, `cpu_model`, `host`, `ref_commit` or `injected` differ from the base session's. It raises the same error when the extension's output digests differ from the base workload's. `RatioRegressionError(report)` is defined here; the CLI raises it. Its message names each workload, `ρ_w`, `m_w` and the pair count.

- [ ] **Step 1: Write the failing tests:**
  - `test_decide_boundaries`: margin 0.04; ratio exactly 1.04 → PASS, just above → REGRESSION, exactly 0.96 → PASS, below → IMPROVED. Build the ratios from exact binary CPUs, e.g. ref 1.0, head 1.0 + margin.
  - `test_a_passing_session_needs_no_measurement`: the fake `Measure` records calls, and asserts none were made.
  - `test_a_noisy_red_workload_extends_and_passes`: base min ratio 1.10. The extension returns pairs whose ref min is lower, giving a merged ratio of 1.01 → PASS, `pairs == 10`, one request for `(w,)` with `pairs=5`.
  - `test_a_persistent_regression_stops_at_max_pairs`: every batch stays red → REGRESSION, `pairs == MAX_PAIRS`, 3 requests.
  - `test_the_extension_replays_the_injection`: an audit session injected into head; the request carries the same `Injection`.
  - `test_an_extension_from_another_build_is_refused` (Review Focus 2): the extension's head `native_sha256` differs → `ExtensionMismatchError`. Repeat with a different `injected`.
  - `test_a_session_from_another_harness_is_refused` (Review Focus 1) → `HarnessMismatchError`.
  - `test_a_session_against_another_reference_is_refused` → `ReferenceMismatchError`.
  - `test_an_injected_session_is_refused_outside_audit` → `InjectedSessionError`.
  - `test_a_partial_session_is_refused_outside_audit` (Review Focus 4) → `IncompleteSessionError`; the same session with `audit=True` is decided.
  - `test_an_unregistered_workload_is_refused` → `UnregisteredWorkloadError`.
  - `test_changed_output_is_reported_not_refused`: `head_digest != ref_digest` → `changed_output True`, verdict still decided.
- [ ] **Step 2: Run red.** Expected: import failure.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run green + R6** (one deleted check per refusal → its test red). Lint.
- [ ] **Step 5: Hand to orchestrator.** pathspec = the two files.

### Task 3: CLI: check, calibrate, promote

**Files:** Create `tools/perf_ratio_gate.py`, `tests/tools/test_perf_ratio_gate.py`.

**Interfaces:**
- Consumes: Tasks 1–2; `tools.perf_ratio_arms.DEFAULT_REF_WORKTREE`, `checkout_head(root, rev) -> str` (moves a clean worktree; it refuses a dirty one or an unknown rev); `tools.perf_ratio_record.load_session`.
- Produces:

```python
GATE_SCRATCH = REPOSITORY / "build" / "ratio-gate"   # extension and calibration-run records before they are committed

def harness_command(request: ExtensionRequest, write: Path) -> tuple[str, ...]
    # (sys.executable, "-m", "tools.perf_ratio", "--ref-rev", ref, "--head-rev", head,
    #  *("--workload", w per w), "--pairs", str(pairs), "--write", str(write), [ "--inject-cpu-seconds", f"{arm}={s!r}" ])
def subprocess_measure(scratch: Path) -> Measure   # runs harness_command with check=True, then load_session(write);
                                                     # write = scratch / f"ext-{n}.json", n counting up
def main(argv: Optional[Sequence[str]] = None) -> int
```

Subcommands (`python -m tools.perf_ratio_gate`):
- `check RECORD [--audit] [--registry PATH]`:
  - prints one line per workload: `<workload>: <VERDICT> rho <ρ:.3f> margin <m:.3f> pairs <n>[ changed output]`, then a summary line;
  - exits 0 when there is no regression;
  - raises `RatioRegressionError` (exit 1) otherwise;
  - with `IMPROVED` present, it prints `promotion available: pixi run perf-ratio-promote RECORD`.
- `calibrate --records PATH... | --sessions K`:
  - `--records` builds the calibration for those records' `cpu_model` and merges it into the registry. It replaces that cpu's calibration and keeps the others. With no registry yet, `ref_commit` and `harness_sha256` come from the records, and all records must agree.
  - `--sessions K` (K ≥ `CALIBRATION_SESSIONS`) runs K A/A sessions with `--ref-rev ref --head-rev ref`, the default set and the default pairs. Each session writes to `calibration_dir(ref, cpu)/aa-<k>/perf_ratio_session.json`, numbered k = 1..K, then the step proceeds as `--records`.
  - It refuses an existing target directory with `CalibrationExistsError`.
- `promote RECORD`:
  - loads the session and runs `gate` (not audit);
  - refuses with `PromotionRefusedError` (one message per reason) when there is any REGRESSION or no IMPROVED;
  - refuses when the head arm is dirty;
  - injected sessions are already refused by `require_gateable`.
  - Then it runs `checkout_head(DEFAULT_REF_WORKTREE, head_rev)` and writes a registry with `ref_commit = head_rev`, the same harness and **no calibrations**. Every check is then `UncalibratedMachineError` until `calibrate` runs at the new reference (spec §Ratchet). It prints that next step.

- [ ] **Step 1: Write the failing tests** (synthetic records in `tmp_path`, fake `Measure` injected through a module-level seam `main(argv, measure=...)`, and a fake `checkout` callable for promote):
  - `test_harness_command_shape`: exact tuple for one request with injection, and without.
  - `test_check_pass_exits_zero_and_prints_each_workload` (capsys).
  - `test_check_regression_raises_with_workload_rho_margin_pairs`.
  - `test_check_improved_prints_the_promotion_hint`.
  - `test_calibrate_from_records_writes_a_loadable_registry`, then `load_registry` margins equal `derive_margins`.
  - `test_calibrate_replaces_only_its_cpu`: two cpu models; re-calibrating one leaves the other's `Calibration` equal.
  - `test_promote_moves_the_reference_and_clears_calibrations`: fake checkout called once with the head rev; the new registry has `ref_commit == head_rev` and `calibrations == ()`.
  - `test_promotion_refusals_leave_the_registry_alone` (Review Focus 5): regression; no improvement; dirty head; injected. Each raises, the registry bytes are unchanged, and checkout is not called.
- [ ] **Step 2: Run red.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Green + R6 + lint.**
- [ ] **Step 5: Hand to orchestrator.** pathspec = the two files.

### Task 4: pixi tasks, typecheck, the initial registry (C++ slot, alone)

**Files:** Modify `pyproject.toml` (tasks beside `perf-ratio-render`, `types-benchmarks` line). Generate `benchmarks/results/perf_ratio_gate.json`.

- [ ] **Step 1: Add three tasks** with the same `editable_build_dir` / `SKBUILD_EDITABLE_SKIP` prefix as `perf-ratio`:
  - `perf-ratio-check` → `python -m tools.perf_ratio_gate check`;
  - `perf-ratio-calibrate` → `… calibrate`;
  - `perf-ratio-promote` → `… promote`.

  Each gets a one-line description. Append `tools/perf_ratio_gate.py tools/perf_ratio_gate_verdict.py tools/perf_ratio_gate_registry.py` to `types-benchmarks`.
- [ ] **Step 2: Verify.** `pixi run --no-install perf-ratio-check --help` exits 0. `git diff --quiet -- pixi.lock && echo lock-unchanged`. `pixi run --no-install types-benchmarks` → `Success`.
- [ ] **Step 3: Build the initial registry** from the committed A/A records: `pixi run --no-install perf-ratio-calibrate --records benchmarks/results/perf_ratio_calibration/15fbed70/apple-m1-max/aa-{2,3,4,5,6}/perf_ratio_session.json` (a bash list). Expected: the margins equal `build/rtr/aa-margins.py`'s to 3 decimals (0.039, 0.013, 0.021, 0.015, 0.007, 0.025, 0.016, 0.043, 0.048, 0.029), `ref_commit 15fbed70…`, harness `3a76eb83…`.
- [ ] **Step 4: Lane tests together**: the five landed `test_perf_ratio*.py` plus the three gate test files plus `test_perf_history.py`, `-n 2`. Expected: all pass except the one known stale test.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `pyproject.toml benchmarks/results/perf_ratio_gate.json`.

### Task 5: Audit and docs (orchestrator executes)

- [ ] **Step 1: A/A → PASS.** Run `PYTHONUNBUFFERED=1 nohup pixi run --no-install perf-ratio --ref-rev <registry ref> --head-rev <registry ref> --write build/ratio-gate/audit-aa.json` in the background (cap 1912 s), then `pixi run --no-install perf-ratio-check build/ratio-gate/audit-aa.json`. Expected: PASS on all 10 workloads, exit 0. A REGRESSION after extension is a gate slip, so write `docs/superpowers/state/ratio-gate-X-analysis.md` before anything else.
- [ ] **Step 2: Inject `2·m_w` into HEAD → REGRESSION.** Workload `gt/cap-120-default/pipeline/rect_20x12`. Size `S = 2 · m_w · max_k ref_min_k` over its 5 calibration sessions, rounded up to the millisecond. Expected ρ ≈ R_aa + 2·m_w. Run with `--workload … --inject-cpu-seconds head=S`, then `check --audit`. Expected: REGRESSION at `MAX_PAIRS`, exit 1.
- [ ] **Step 3: The same S into ref → IMPROVED** (`check --audit`, exit 0, promotion hint printed).
- [ ] **Step 4: Docs.** In `docs/performance.md`, extend `## Ratio tracking against a pinned reference` with the gate:
  - the three commands;
  - the verdict rule;
  - sequential extension (`EXTENSION_PAIRS`, `MAX_PAIRS`);
  - the ratchet and recalibration after promotion;
  - the refusals;
  - the hardware table (spec §The hardware answer, condensed);
  - the three audit results.

  Guard: `tests/docs/test_included_records.py` plus the strict mkdocs build.
- [ ] **Step 5: Commit batch** (user's go). pathspec = the Task 1–4 files plus `docs/performance.md`. CI wiring is next, as its own unit: the analysis is shown before any workflow edit.

## Self-review

- Spec coverage:
  - §Decision / verdict → Task 2.
  - §Hardware → Constraints and the Task 5 docs.
  - Sequential extension → Task 2.
  - §Margins (SAFETY, K, content-addressed, re-derive) → Task 1.
  - §Ratchet (promote, re-pin, recalibrate) → Task 3.
  - §Verdict refusals → Task 2 (and Task 1 for calibration).
  - §Audit → Task 5.
  - Open question 1 → `UncalibratedMachineError` (Tasks 1–2).
  - §Cost / CI → deferred unit, Task 5 Step 5.
- Deviation from spec: `same_output` (header amendment, user go). `PINNED_REF_COMMIT` is not edited. The registry is the authority via `--ref-rev`, because editing it would void the calibration.
- Names: `Margin`, `GateRegistry.margins_for`, `ExtensionRequest`, `Measure`, `gate`, `harness_command` are used consistently across Tasks 1–3.
