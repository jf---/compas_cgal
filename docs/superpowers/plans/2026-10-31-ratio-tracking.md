# Ratio Tracking Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pixi run perf-ratio` times a fixed workload set interleaved on two clean worktrees (pinned reference, HEAD), as in-child process CPU, and appends HEAD/reference ratios to their own append-only history, rendered into `docs/performance.md`.

**Architecture:** The orchestrator runs from the task5 worktree and is neither arm. Each arm is opened, rebuilt and proven with the imported 1.1-base helpers (`tools/pipeline_timing_1_1_base.py`) and identity APIs (`tools/run_provenance.py`, `tools/run_native_identity.py`). Every child is a fresh process in its arm's own interpreter. Its CPU is reported in-child: by the tool's own clock for polygon, Ω* and certificate, and by a harness-owned `-c` snippet for GT. One run writes one `perf-ratio/1` record. A separate ingest turns the record into `perf-ratio-entry/1` lines of `perf_ratio_history.jsonl` and renders `perf_ratio_history.md`.

**Tech Stack:** Python 3.9-compatible, stdlib only (argparse, subprocess, statistics, hashlib, json, pathlib), pytest + xdist, mypy --strict via `types-benchmarks`, ruff, pixi.

**Spec:** `docs/superpowers/specs/2026-10-30-ratio-tracking-design.md` (recon: `build/rtr/recon.md`).

## Global Constraints

- Work only in `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation`. Python only via `pixi run --no-install build/pyrun.sh ...`. Never bare python/pytest/pip.
- Python 3.9-compatible: `from __future__ import annotations`; `typing.Optional`, `Mapping`, `Sequence`; no `match`, no `X | Y` at runtime.
- `mypy --strict` through `types-benchmarks`, run only after Task 6. One named exception per failure mode, never bare `ValueError`. Google docstrings. pathlib only. No inline magic numbers: named UPPERCASE constants with a one-line comment giving units and derivation.
- Reuse, not a new harness. Import, never move or edit, `base_child_env`, `editable_build_dir_of`, `library_origin_of`, `require_library_in_tree`, `run_in_tree`, `rebuild_base_build` and `require_measured_tree`, plus `tree_identity`, `interpreter_identity_of`, `machine_identity`, `require_unchanged_tree` and `track_native_in`.
- Every arm number is in-child process CPU. No wall clock. No `RUSAGE_CHILDREN`. No `pixi run` launcher CPU.
- Both arms are identified the same way, and HEAD is never special-cased in-process. `native_sha256` is recorded per arm and never compared across arms.
- No load refusal: load is recorded, never gated. `MachineNotQuietError` is not used.
- There is no ratio gate in v1. Ratios never enter `perf_history.jsonl` or its baseline, check and render.
- Only Task 6 (C++ slot) edits `pyproject.toml`.
- Unit tests use fake arms only (tmp git trees with stub children). No timing assertions on real workloads.
- Targeted tests: `pixi run --no-install build/pyrun.sh -m pytest <file> -n 2 -q` (no lockf). Lint: `pixi run --no-install ruff check <files>`.
- Commits: the orchestrator commits. Every commit step is "hand to orchestrator: pathspec = ...". No AI mention in any code, docstring, doc or message.

## Review Focus

1. **The orchestrator's pixi environment leaks into an arm.** `base_child_env` drops only `PYTHONPATH` and `SKBUILD_EDITABLE_SKIP` (`tools/pipeline_timing_1_1_base.py:106`). An arm's `pixi run` could resolve task5's manifest and run task5's interpreter in both arms, giving a silent 1.0. Expected: each arm command names its own manifest, and a shared interpreter prefix is refused. Pinned by `test_arm_python_names_the_arm_manifest` and `test_arms_sharing_an_interpreter_prefix_are_refused` (Task 3).
2. **`--head-rev` names no commit, or ratio-head is dirty.** Expected: refused before checkout, with ratio-head untouched. Pinned by `test_an_unknown_head_rev_is_refused_and_leaves_the_tree_alone` and `test_a_dirty_head_tree_is_refused_before_checkout` (Task 3).
3. **A child reports zero CPU, or prints no result line or two.** Examples are a stub, or a tool that renamed its marker. Expected: a named refusal, never an inf/nan ratio or a first-line guess. Pinned by `test_a_non_positive_cpu_run_is_refused` (Task 1) and `test_zero_or_two_result_lines_are_refused` (Task 2).
4. **A malformed `--inject-cpu-seconds`** (`head=abc`, `both=1`, `head=-1`, `head=0`), or injection with a non-GT workload. Expected: refused before any rebuild or child. Pinned by `test_malformed_injection_is_refused` and `test_injection_with_a_non_gt_workload_is_refused_before_anything_runs` (Task 4).
5. **The same session is ingested twice, or a history line is hand-edited.** Expected: re-ingest is a no-op, and an edited line is refused on load. Pinned by `test_reingesting_a_record_appends_nothing` and `test_an_edited_history_line_is_refused` (Task 5).

---

### Task 1: `perf-ratio/1` record types

**Files:**
- Create: `tools/perf_ratio_record.py`
- Test: `tests/tools/test_perf_ratio_record.py`

**Interfaces:**
- Consumes: `TreeIdentity(rev: str, dirty: bool, patch_sha256: Optional[str])` (`tools/run_provenance.py:95`); `NativeIdentity(native_sha256: str, native_sources_sha256: str)` (`tools/run_native_identity.py:82`); `obj`, `items`, `field`, `number`, `text`, `optional_text`, `integer` and `MalformedTimingRecordError` (`tools/perf_history_json.py:17-104`) for checked reads.
- Produces:
  - Constant `PERF_RATIO_SCHEMA = "perf-ratio/1"`.
  - NewTypes `CpuSeconds = NewType("CpuSeconds", float)`, `Ratio = NewType("Ratio", float)` and `LoadOneMinute = NewType("LoadOneMinute", float)`.
  - `class ArmName(str, Enum): REF = "ref"; HEAD = "head"`.
  - `@dataclass(frozen=True) class ArmIdentity` with fields `rev: str`, `dirty: bool`, `patch_sha256: Optional[str]`, `native_sha256: str`, `native_sources_sha256: str`, `environment: str`, `python_version: str`, `stock_2_sha256: str`, `pixi_lock_sha256: str`, `compas_cgal_file: str`, `stock_2_origin: str`, `editable_build_dir: str` and `interpreter_prefix: str`.
  - `@dataclass(frozen=True) class PairRun` with fields `order: tuple[ArmName, ArmName]`, `ref_cpu: CpuSeconds`, `head_cpu: CpuSeconds`, `load_before: LoadOneMinute` and `load_after: LoadOneMinute`; `classmethod build(...) -> PairRun`; property `ratio -> Ratio` (head/ref).
  - `@dataclass(frozen=True) class WorkloadRatio` with fields `workload: str`, `family: str`, `pairs: tuple[PairRun, ...]`, `ref_digest: str` and `head_digest: str`. Properties: `ref_min`, `ref_median`, `head_min`, `head_median -> CpuSeconds`; `ratio -> Ratio` (= `head_min / ref_min`); `pair_ratio_median -> Ratio`; `same_output -> bool`; `load_range -> tuple[LoadOneMinute, LoadOneMinute]`.
  - `@dataclass(frozen=True) class Injection` with fields `arm: ArmName` and `seconds: CpuSeconds`.
  - `@dataclass(frozen=True) class PerfRatioSession` with fields `ref_commit: str`, `arms: Mapping[ArmName, ArmIdentity]`, `workloads: tuple[WorkloadRatio, ...]`, `host: str`, `cpu_model: str`, `started_utc: str`, `ended_utc: str`, `harness_sha256: str` and `injected: Optional[Injection]`. Methods: `classmethod build(raw: object) -> PerfRatioSession`, `to_json() -> dict[str, object]` and `dumps() -> str`. Properties: `stage -> str` (`f"ratio-vs-{ref_commit[:REF_SHORT_CHARS]}"`, `REF_SHORT_CHARS = 8`) and `wall_seconds -> float`.
  - `load_session(path: Path) -> PerfRatioSession`.
  - Errors: `MalformedRatioRecordError` (wrong shape, schema or type); `InconsistentRatioRecordError` (a stored min, median, ratio, pair-ratio median or same_output differs from its recomputation from the runs); `NonPositiveCpuError` (a run's CPU ≤ 0).
- `to_json` stores the derived values beside the runs, for readers. `build` recomputes each one and refuses a mismatch with exact equality: the same float operations give the same bits, so no tolerance is needed.

- [ ] **Step 1: Write the failing tests.** File: `tests/tools/test_perf_ratio_record.py`.
  - `test_a_session_round_trips_through_json`: `PerfRatioSession.build(json.loads(session.dumps())) == session`.
  - `test_ratio_is_min_head_over_min_ref_and_pair_median_beside_it`: pairs (ref, head) = (2.0, 3.0), (1.0, 4.0), (4.0, 2.0). Assert `ratio == 2.0 / 1.0`, `pair_ratio_median == statistics.median([1.5, 4.0, 0.5])`, `ref_min == 1.0` and `head_median == 3.0`.
  - `test_a_tampered_derived_value_is_refused`: edit `workloads[0]["ratio"]` in the JSON; `pytest.raises(InconsistentRatioRecordError)`.
  - `test_a_non_positive_cpu_run_is_refused`: `PairRun.build(..., ref_cpu=0.0, ...)` raises `NonPositiveCpuError`; so does `head_cpu=-1.0`.
  - `test_another_schema_is_refused`: `"schema": "perf-ratio/0"` raises `MalformedRatioRecordError`.
  - `test_the_stage_names_the_reference_short_rev`: `stage == "ratio-vs-15fbed70"` for ref `15fbed70ebb016bdd8ef49453e04606d03e0b1fe`.
  - `test_both_arms_are_required`: an `arms` without `"head"` raises `MalformedRatioRecordError`.
  - `test_same_output_is_digest_equality`: equal digests give True, different digests give False, and the record still builds (flagged, not refused).
- [ ] **Step 2: Run them; they fail.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_record.py -n 2 -q`. Expected: collection error `ModuleNotFoundError: No module named 'tools.perf_ratio_record'`.
- [ ] **Step 3: Implement the interfaces above in `tools/perf_ratio_record.py`.** Medians use `statistics.median`. Reads go through the `perf_history_json` accessors, re-raised as `MalformedRatioRecordError`.
- [ ] **Step 4: Run the tests and lint.** Run: the Step 2 command, then `pixi run --no-install ruff check tools/perf_ratio_record.py tests/tools/test_perf_ratio_record.py`. Expected: `8 passed`; `All checks passed!`.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `tools/perf_ratio_record.py tests/tools/test_perf_ratio_record.py`.

### Task 2: Workload registry, child commands, parsing, GT CPU clock

**Files:**
- Create: `tools/perf_ratio_workloads.py`
- Test: `tests/tools/test_perf_ratio_workloads.py`

**Interfaces:**
- Consumes:
  - `WORKER_MARKER = "HELD_GENERATION_CLOCK_RESULT "` (`tools/held_generation_clock.py:91`). The worker CLI is `python -m tools.held_generation_clock --domain {polygon|omega} <workload>` (`main`, `tools/held_generation_clock.py:530-550`). It prints one line, `WORKER_MARKER + {"digest": {...}, "clock": {...}}`.
  - `GenerationClock.from_json(text: str) -> GenerationClock` and `GenerationClock.cpu_seconds` (`tools/held_generation_clock.py:170-207`).
  - `BUDGET_WORKER_MARKER = "CERTIFICATE_BUDGET_RUN "` (`tools/coverage_certificate_budget_gate.py:44`), from `main(argv)` with `--worker` (`:291-295`). Its row carries `certificate_cpu_seconds`, `verdict`, `witnesses`, `edges` and the other fields of `measure_workload` (`tools/held_omega_coverage_results.py:125-173`).
  - `RESULT_LINE_PREFIX = "PIPELINE_TIMING_RESULT "` (`tools/pipeline_timing.py:76`), from `main(argv: Sequence[str] | None = None) -> int` with `["--case", key]` (`:326-345`). `run_case_here` reads `time.perf_counter()` exactly twice, around the generator call (`:193-195`).
  - `CpuSeconds` from Task 1.
- Produces:
  - `class Family(str, Enum): POLYGON = "polygon"; OMEGA = "omega"; CERTIFICATE = "certificate"; GT = "gt"`.
  - `@dataclass(frozen=True) class Workload` with fields `family: Family` and `argument: str`, and property `name -> str` (`f"{family.value}/{argument}"`).
  - `HELD_WORKLOADS = ("figure5_refined", "figure8_crossed_skis_boundary", "figure8_monstera_boundary_qualified_corners")`.
  - `CERTIFICATE_WORKLOAD = "figure8_monstera_boundary_qualified_corners"`.
  - `GT_CAP120_KEYS = ("cap-120-default/pipeline/rect_12x8", "cap-120-default/pipeline/rect_20x12", "cap-120-default/pipeline/L_shape")`.
  - `GT_CAP40_KEYS = ("cap-40-attribution/pipeline/rect_12x8", "cap-40-attribution/pipeline/rect_20x12", "cap-40-attribution/pipeline/L_shape")`.
  - `DEFAULT_SET_NAME = "default"`, `GT_CAP40_SET_NAME = "gt-cap40"`, and `WORKLOAD_SETS: Mapping[str, tuple[Workload, ...]]`. Default = 3 polygon + 3 omega + 1 certificate + 3 GT cap-120, in that order (10 workloads).
  - `workload_named(name: str) -> Workload` (searches every set).
  - `child_command(workload: Workload, python: Sequence[str], inject: CpuSeconds) -> tuple[str, ...]`. GT: `(*python, "-c", GT_CPU_CLOCK_SNIPPET, key, repr(float(inject)))`.
  - `class ChildRun(NamedTuple): cpu_seconds: CpuSeconds; digest: str`.
  - `parse_child(workload: Workload, stdout: str) -> ChildRun`. Polygon and Ω*: CPU = `GenerationClock.from_json(json.dumps(raw["clock"])).cpu_seconds`, digest = SHA-256 of `json.dumps(raw["digest"], sort_keys=True)`. Certificate: CPU = `row["certificate_cpu_seconds"]`, digest = SHA-256 of `json.dumps({k: row[k] for k in CERTIFICATE_DIGEST_FIELDS}, sort_keys=True)`. GT: CPU = `generation_s`, digest = `f"{digest}:{operations}"`.
  - `CERTIFICATE_DIGEST_FIELDS = ("workload", "circles", "connectors", "sub_pieces", "omega_pieces", "vertex_disks", "repair_inserted", "verdict", "witnesses", "witnesses_by_kind", "undecided_edges", "regions", "pieces", "edges")`: the row's deterministic fields; the timing and load fields are excluded.
  - `GT_CLOCK_READS = 2` (`run_case_here`'s two `perf_counter` reads).
  - Errors: `UnknownWorkloadError`, `InjectionUnsupportedError` (inject > 0 with a non-GT family), `MissingChildResultError` (zero result lines or more than one).

The GT CPU-clock snippet. Its body is pinned by the spec; copy it exactly:

```python
# Runs in an arm's interpreter: argv = [-c, case key, injected CPU seconds]. Imports THAT arm's
# tools.pipeline_timing, gives the module a CPU stand-in for `time`, and spins the injected CPU
# inside the measured span (before the end read returns).
GT_CPU_CLOCK_SNIPPET = """
import sys, time
import tools.pipeline_timing as timed
key, inject = sys.argv[1], float(sys.argv[2])
class _CpuClock:
    reads = 0
    def perf_counter(self):
        _CpuClock.reads += 1
        if _CpuClock.reads == 2 and inject > 0.0:
            until = time.process_time() + inject
            while time.process_time() < until:
                pass
        return time.process_time()
vars(timed)["time"] = _CpuClock()
status = timed.main(["--case", key])
if _CpuClock.reads != 2:
    raise SystemExit(f"tools.pipeline_timing read its clock {_CpuClock.reads} times, not 2: move the reference")
sys.exit(status)
"""
```

- [ ] **Step 1: Write the failing tests.** File: `tests/tools/test_perf_ratio_workloads.py`. The stub tool modules are written into a `tmp_path` tree with a `tools/__init__.py`. The child runs as `subprocess.run(child_command(...), cwd=tmp_path, env={**os.environ, "PYTHONPATH": str(tmp_path)})` with `python = (sys.executable,)`.
  - `test_the_default_set_is_the_specs_ten_workloads`: `[w.name for w in WORKLOAD_SETS["default"]]` equals the 10 names `polygon/<w>`, `omega/<w>` (each `w` in `HELD_WORKLOADS`), `certificate/figure8_monstera_boundary_qualified_corners` and `gt/<k>` (each `k` in `GT_CAP120_KEYS`); `gt-cap40` holds only the GT cap-40 keys.
  - `test_the_gt_keys_are_real_gate_cases`: `set(GT_CAP120_KEYS) | set(GT_CAP40_KEYS) <= set(gate_timing_keys()[0])` (`tools/pipeline_timing.py:316`).
  - `test_held_and_certificate_commands_are_the_tools_worker_forms`: `child_command(Workload(Family.OMEGA, "figure5_refined"), ("py",), 0.0) == ("py", "-m", "tools.held_generation_clock", "--domain", "omega", "figure5_refined")`; certificate gives `("py", "-m", "tools.coverage_certificate_budget_gate", "--worker")`.
  - `test_injection_is_refused_outside_gt`: every non-GT family with `inject=1.0` raises `InjectionUnsupportedError`.
  - `test_the_gt_snippet_reports_cpu_from_the_arms_own_module`: a stub `tools/pipeline_timing.py` whose `main` busy-loops between two `time.perf_counter()` reads and prints `RESULT_LINE_PREFIX + {"generation_s", "operations": 7, "digest": "d"}`. `parse_child` gives `digest == "d:7"` and `cpu_seconds > 0`.
  - `test_the_gt_snippet_spins_the_injected_cpu_inside_the_span`: the same stub with `inject=INJECT_PROBE_S` (a test constant, 0.2 s CPU) gives `cpu_seconds >= INJECT_PROBE_S`. This is a lower bound from `process_time` on a stub, not a workload timing.
  - `test_a_moved_clock_fails_loud`: a stub that reads `perf_counter` three times makes the child exit non-zero, with stderr containing `"not 2"`.
  - `test_zero_or_two_result_lines_are_refused`: `parse_child` on stdout with no marker line, or with two, raises `MissingChildResultError`, for every family.
  - `test_the_certificate_digest_ignores_timing_and_load`: two rows that differ only in `certificate_cpu_seconds`, `certificate_wall_seconds`, `generation_seconds` and `load_1min_before_after` give equal digests; rows that differ in `edges` do not.
  - `test_an_unknown_workload_name_is_refused`: `workload_named("gt/nope")` raises `UnknownWorkloadError`.
- [ ] **Step 2: Run them; they fail.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_workloads.py -n 2 -q`. Expected: `ModuleNotFoundError: No module named 'tools.perf_ratio_workloads'`.
- [ ] **Step 3: Implement the interfaces above in `tools/perf_ratio_workloads.py`.** Import the three markers from their tools; never copy them.
- [ ] **Step 4: Run the tests and lint.** Run: the Step 2 command, then `pixi run --no-install ruff check tools/perf_ratio_workloads.py tests/tools/test_perf_ratio_workloads.py`. Expected: `10 passed`; `All checks passed!`.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `tools/perf_ratio_workloads.py tests/tools/test_perf_ratio_workloads.py`.

### Task 3: Arm setup and per-arm proofs (with the fake-arm fixtures)

**Files:**
- Create: `tools/perf_ratio_arms.py`
- Create: `tests/tools/perf_ratio_fake_arms.py` (test helper: fake arm trees with stub children; also used by Task 4)
- Test: `tests/tools/test_perf_ratio_arms.py`

**Interfaces:**
- Consumes, from `tools/pipeline_timing_1_1_base.py`:
  - `require_measured_tree(root: Path, expected_rev: str) -> None`, which raises `MeasuredTreeError` (`:143`, `:118`).
  - `base_child_env(root: Path, environ: Mapping[str, str], editable_build_dir: Optional[Path]) -> dict[str, str]` (`:156`).
  - `editable_build_dir_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> Path` (`:171`).
  - `library_origin_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> LibraryOrigin` (`:181`), with `LibraryOrigin(compas_cgal_file: Path, stock_2_origin: Path)` (`:131`).
  - `require_library_in_tree(origin: LibraryOrigin, root: Path, editable_build_dir: Optional[Path] = None) -> None`, which raises `ForeignLibraryError` (`:194`, `:122`).
  - `run_in_tree(root: Path, env: Mapping[str, str]) -> Runner`, where `Runner = Callable[[Sequence[str]], str]` and a failing command raises `CaseProcessFailedError` (`:209`, `:115`; `tools/pipeline_timing.py:83`).
  - `rebuild_base_build(root: Path, run: Runner) -> tuple[str, ...]`, which raises `BaseRebuildFailedError` (`:225`, `:126`). Its runner pins the orchestrator's own build: `run_in_tree(root, base_child_env(root, os.environ, this_build_dir))`, where `this_build_dir = editable_build_dir_of((sys.executable,), ROOT, os.environ)` (`:371-373`).
- Consumes, from `tools/run_provenance.py`:
  - `tree_identity(root: Path) -> TreeIdentity` (`:316`); `TreeReader = Callable[[], TreeIdentity]` (`:123`).
  - `require_unchanged_tree(before: TreeIdentity, after: TreeIdentity, label: str) -> None`, which raises `TreeChangedDuringRunError` (`:266`, `:126`).
  - `interpreter_identity_of(python: Sequence[str], cwd: Path, env: Optional[dict[str, str]] = None) -> tuple[str, str, str]`, giving (environment, python_version, stock_2_sha256) (`:375`).
  - `sha256_file(path: Path) -> str` (`:286`).
- Consumes, from `tools/run_native_identity.py`: `track_native_in(python: Sequence[str], root: pathlib.Path, env: Optional[Mapping[str, str]]) -> NativeTracker`, where `NativeTracker = Callable[[], NativeFinisher]` and `NativeFinisher = Callable[[], NativeIdentity]` (`:182`, `:94-97`). A stale build raises a `NativeIdentityError` subclass, e.g. `NativeBuildStaleError` (`tools/held_native_identity.py:103`, `:111`).
- Consumes, from Task 1: `ArmName` and `ArmIdentity`.
- Produces:
  - Constants: `PINNED_REF_COMMIT = "15fbed70ebb016bdd8ef49453e04606d03e0b1fe"`, `ROOT`, `DEFAULT_REF_WORKTREE = ROOT.parent / "compas_cgal_prs-ratio-ref"`, `DEFAULT_HEAD_WORKTREE = ROOT.parent / "compas_cgal_prs-ratio-head"`, `PIXI_LOCK = Path("pixi.lock")`, `PYPROJECT = Path("pyproject.toml")` and `INTERPRETER_PREFIX_SNIPPET = "import sys; print(sys.prefix)"`.
  - `arm_python(root: Path) -> tuple[str, ...]` = `("pixi", "run", "--no-install", "--manifest-path", str(root / PYPROJECT), "python")`.
  - `interpreter_prefix_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> Path` (resolved).
  - `@dataclass(frozen=True) class Arm` with fields `name: ArmName`, `root: Path`, `python: tuple[str, ...]`, `env: Mapping[str, str]`, `editable_build_dir: Path`, `interpreter_prefix: Path` and `rev: str`.
  - `open_arm(name: ArmName, root: Path, python: Sequence[str], expected_rev: str) -> Arm`. Order: `require_measured_tree`; then `editable_build_dir_of` (unpinned env); then the pinned `base_child_env`; then `require_library_in_tree(library_origin_of(...))`; then `interpreter_prefix_of`.
  - `require_distinct_arms(ref: Arm, head: Arm) -> None`, which raises `SameArmError` when the resolved roots, editable build dirs or interpreter prefixes are equal.
  - `checkout_head(root: Path, rev: str) -> str`. Order: refuse if dirty (`require_measured_tree(root, tree_identity(root).rev)`); resolve `git rev-parse --verify <rev>^{commit}` (failure → `HeadRevUnknownError`); run `git -C <root> checkout --detach <full>`; check `require_measured_tree(root, full)`; return the full rev.
  - `rebuild_arm(arm: Arm, orchestrator_build_dir: Path) -> tuple[str, ...]`.
  - `require_settled(track_native: NativeTracker) -> None`: one full bracket, `track_native()()`, before timing. A `NativeIdentityError` propagates.
  - `ArmBracket = Callable[[], Callable[[], ArmIdentity]]`.
  - `bracket_arm(arm: Arm, *, read_tree: TreeReader, track_native: NativeTracker, origin_of: Callable[[], LibraryOrigin], interpreter_of: Callable[[], tuple[str, str, str]]) -> ArmBracket`. Start: `read_tree()`, then `track_native()`. Finish: `require_library_in_tree(origin_of(), root, editable_build_dir)`; `read_tree()`; `require_unchanged_tree(start, end, f"the {name.value} arm")`; native finish; `interpreter_of()`; `sha256_file(root / PIXI_LOCK)`; then the `ArmIdentity`.
  - `live_bracket(arm: Arm) -> ArmBracket`: `bracket_arm` wired to `tree_identity`, `track_native_in(arm.python, arm.root, arm.env)`, `library_origin_of` and `interpreter_identity_of(arm.python, arm.root, dict(arm.env))`.
  - Errors: `SameArmError`, `HeadRevUnknownError`.
- Fake-arm helper `tests/tools/perf_ratio_fake_arms.py`:
  - `fake_arm_tree(root: Path, *, cpu_seconds: float, digest: str) -> str`: `git init`, then write `tools/__init__.py`, the stubs `tools/held_generation_clock.py`, `tools/coverage_certificate_budget_gate.py` and `tools/pipeline_timing.py`, and `pixi.lock`. Each stub prints its tool's real marker with the given CPU and digest; the GT stub reads `time.perf_counter()` twice. Then commit and return the rev.
  - `fake_arm(name: ArmName, root: Path, rev: str) -> Arm`: `python=(sys.executable,)`, `editable_build_dir=root / "build" / "editable"`, `interpreter_prefix=root / ".env"`, `env=base_child_env(root, os.environ, root / "build" / "editable")`.
  - `fake_bracket(arm: Arm, *, trees: Iterator[TreeIdentity]) -> ArmBracket`: `bracket_arm` with an in-memory tracker and an origin inside the tree.

- [ ] **Step 1: Write the failing tests.** File: `tests/tools/test_perf_ratio_arms.py`, using `perf_ratio_fake_arms`.
  - `test_arm_python_names_the_arm_manifest`: `arm_python(tmp_path)` contains `"--manifest-path"` followed by `str(tmp_path / "pyproject.toml")`.
  - `test_arms_sharing_a_build_dir_are_refused`: two fake arms with the same `editable_build_dir` raise `SameArmError`.
  - `test_arms_sharing_an_interpreter_prefix_are_refused`: the same, with the same `interpreter_prefix`.
  - `test_one_tree_as_both_arms_is_refused`: the same root raises `SameArmError`.
  - `test_distinct_fake_arms_pass`: `require_distinct_arms` returns None.
  - `test_a_foreign_library_is_refused`: `open_arm(ArmName.HEAD, fake_tree, (sys.executable,), rev)` raises `ForeignLibraryError`. The test interpreter's `compas_cgal` resolves to task5, outside the fake tree.
  - `test_a_tree_changed_mid_run_is_refused`: `fake_bracket` with trees `[TreeIdentity(rev, False, None), TreeIdentity(rev, True, "ab")]`; calling `start()()` raises `TreeChangedDuringRunError`.
  - `test_a_library_that_moves_mid_run_is_refused`: `origin_of` returns a path outside the root at finish; raises `ForeignLibraryError`.
  - `test_the_bracket_records_both_native_identities_and_the_lock`: the finished `ArmIdentity` holds the tracker's `native_sha256` and `native_sources_sha256`, and `pixi_lock_sha256 == sha256_file(root / "pixi.lock")`.
  - `test_an_unsettled_build_is_refused_before_timing`: a tracker whose finisher raises `NativeBuildStaleError` makes `require_settled` raise it.
  - `test_a_dirty_head_tree_is_refused_before_checkout`: a fake tree with an untracked file; `checkout_head(root, rev)` raises `MeasuredTreeError` and `git rev-parse HEAD` is unchanged.
  - `test_an_unknown_head_rev_is_refused_and_leaves_the_tree_alone`: `checkout_head(root, "0" * 40)` raises `HeadRevUnknownError`; HEAD is unchanged and the tree is clean.
  - `test_checkout_moves_a_clean_head_to_the_rev`: a fake tree with two commits; `checkout_head(root, first)` returns the full first rev, and `tree_identity(root) == TreeIdentity(first, False, None)`.
  - `test_the_arm_rebuild_pins_the_orchestrators_build`: `rebuild_arm` with a recording runner (monkeypatched `run_in_tree`). The command is `("pixi", "run", "--manifest-path", str(root / "pyproject.toml"), "_editable-rebuild")` and the env's `SKBUILD_EDITABLE_SKIP` is the orchestrator build dir passed in.
- [ ] **Step 2: Run them; they fail.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_arms.py -n 2 -q`. Expected: `ModuleNotFoundError: No module named 'tools.perf_ratio_arms'`.
- [ ] **Step 3: Implement `tools/perf_ratio_arms.py` and the fake-arm helper, as in Interfaces.**
- [ ] **Step 4: Run the tests and lint.** Run: the Step 2 command, then `pixi run --no-install ruff check tools/perf_ratio_arms.py tests/tools/perf_ratio_fake_arms.py tests/tools/test_perf_ratio_arms.py`. Expected: `14 passed`; `All checks passed!`.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `tools/perf_ratio_arms.py tests/tools/perf_ratio_fake_arms.py tests/tools/test_perf_ratio_arms.py`.

### Task 4: Orchestration and CLI

**Files:**
- Create: `tools/perf_ratio.py`
- Test: `tests/tools/test_perf_ratio.py`

**Interfaces:**
- Consumes: Tasks 1–3. Also `machine_identity() -> tuple[str, str]` (`tools/run_provenance.py:346`), `utc_now() -> str` (`:247`), `file_set_sha256(root: Path, paths: Sequence[Path]) -> str` (`:261`), `editable_build_dir_of` and `run_in_tree` (Task 3 lines), and `os.getloadavg`.
- Produces:
  - Constants: `DEFAULT_PAIRS = 5` (spec: N = 5 pairs per workload) and `PERF_RATIO_SESSION_PATH = Path("benchmarks") / "results" / "perf_ratio_session.json"`.
  - `HARNESS_SOURCES`: the five `tools/perf_ratio*.py`, plus `tools/pipeline_timing_1_1_base.py`, `tools/run_provenance.py`, `tools/run_native_identity.py` and `tools/held_native_identity.py`. The last two are the probe code imported into both arms.
  - `harness_sha256() -> str` = `file_set_sha256(ROOT, HARNESS_SOURCES)`.
  - `Launch = Callable[[ArmName, Workload], ChildRun]`; `LoadReader1 = Callable[[], LoadOneMinute]`.
  - `pair_order(index: int) -> tuple[ArmName, ArmName]`.
  - `measure_workload(workload: Workload, *, pairs: int, launch: Launch, load_reader: LoadReader1) -> WorkloadRatio`, which raises `NondeterministicArmError`.
  - `measure_session(workloads: Sequence[Workload], *, ref_commit: str, brackets: Mapping[ArmName, ArmBracket], pairs: int, launch: Launch, load_reader: LoadReader1, injected: Optional[Injection], harness: str) -> PerfRatioSession`. Order: start both brackets, measure the workloads, finish both brackets, then read `machine_identity` and the UTC times.
  - `arm_launcher(arms: Mapping[ArmName, Arm], injected: Optional[Injection]) -> Launch`: `parse_child(w, run_in_tree(arm.root, arm.env)(child_command(w, arm.python, inject_for(arm))))`.
  - `parse_injection(text: str) -> Injection`, which raises `InvalidInjectionError`. Form `ARM=S`: ARM in {ref, head}, S finite and > 0.
  - `require_injectable(workloads: Sequence[Workload], injected: Optional[Injection]) -> None`, which raises `InjectionUnsupportedError`.
  - `summary_lines(session: PerfRatioSession) -> list[str]`. One line per workload: `f"{w}: ratio {r:.3f}x (pair median {m:.3f}x, pairs {lo:.3f}-{hi:.3f}), same_output {b}, load {a:.1f}-{z:.1f}"`. Then `f"wall {session.wall_seconds:.0f} s"`.
  - `main(argv: Optional[Sequence[str]] = None) -> int`. Flags:
    - `--ref-worktree` (default `DEFAULT_REF_WORKTREE`) and `--ref-rev` (default `PINNED_REF_COMMIT`);
    - `--head-worktree` (default `DEFAULT_HEAD_WORKTREE`) and `--head-rev` (default `git -C ROOT rev-parse HEAD`);
    - `--set` (choices `WORKLOAD_SETS`, default `"default"`) and `--workload NAME` (repeatable; when given, it replaces `--set`);
    - `--pairs` (default `DEFAULT_PAIRS`), `--inject-cpu-seconds ARM=S` and `--write` (default `PERF_RATIO_SESSION_PATH`).

    Order:
    1. `parse_injection`, then `require_injectable`, before anything else.
    2. `checkout_head`.
    3. `open_arm` ref, `open_arm` head, then `require_distinct_arms`.
    4. `rebuild_arm` both, printing the commands, then `require_settled(track_native_in(...))` both.
    5. `measure_session`, with `live_bracket` per arm.
    6. `require_measured_tree` on both arms again.
    7. `PerfRatioSession.build(session.to_json())`, then write and print the `summary_lines`.
  - Error: `NondeterministicArmError` (runs within one arm disagree on digest), `InvalidInjectionError`.

`pair_order`'s body is pinned by the spec (AB, BA, AB, BA, AB; A = reference):

```python
def pair_order(index: int) -> tuple[ArmName, ArmName]:
    return (ArmName.HEAD, ArmName.REF) if index % 2 else (ArmName.REF, ArmName.HEAD)
```

`measure_workload` reads `load_reader()` before and after every pair, and launches the two arms in `pair_order(i)` for `i in range(pairs)`.

- [ ] **Step 1: Write the failing tests.** File: `tests/tools/test_perf_ratio.py`. A recording `launch` returns scripted `ChildRun`s; the end-to-end test uses `perf_ratio_fake_arms`.
  - `test_pair_order_alternates_starting_with_the_reference`: `[pair_order(i) for i in range(5)] == [(REF, HEAD), (HEAD, REF), (REF, HEAD), (HEAD, REF), (REF, HEAD)]`.
  - `test_each_pair_runs_in_its_order_with_load_read_around_it`: the recorded calls interleave loads and launches as `load, A, B, load` per pair; `len(result.pairs) == 5`; each `PairRun.order` matches `pair_order`.
  - `test_runs_that_disagree_within_an_arm_are_refused`: the head digest changes on pair 3; raises `NondeterministicArmError`.
  - `test_arms_with_different_output_are_recorded_and_flagged`: the ref digest differs from the head digest; `same_output is False` and the session still builds.
  - `test_the_session_brackets_both_arms_around_every_workload`: the bracket events are `start ref, start head, <all launches>, finish ref, finish head`.
  - `test_end_to_end_on_two_fake_arms_writes_a_valid_record`:
    - Setup: two `fake_arm_tree`s (ref CPU 1.0, head CPU 2.0), `fake_bracket`, and `arm_launcher` over the default set with `pairs=2`.
    - Assert `load_session(path).workloads[i].ratio == 2.0` for all 10 workloads, and that `arms[REF].rev` / `arms[HEAD].rev` are the two fake revs.
  - `test_malformed_injection_is_refused`: `head=abc`, `both=1`, `head=-1`, `head=0`, `head=nan` and `head` each raise `InvalidInjectionError`.
  - `test_injection_with_a_non_gt_workload_is_refused_before_anything_runs`:
    - Setup: `main(["--workload", "polygon/figure5_refined", "--inject-cpu-seconds", "head=1", "--head-worktree", str(nonexistent)])`.
    - Assert it raises `InjectionUnsupportedError`, and that no git command ran (`checkout_head` monkeypatched to fail if called).
  - `test_an_injected_session_carries_the_injection`: `measure_session(..., injected=Injection(HEAD, 1.0))` gives `session.injected == Injection(HEAD, 1.0)` in the JSON.
  - `test_the_harness_digest_covers_the_probe_code`: `Path("tools/held_native_identity.py")` and `Path("tools/run_native_identity.py")` are in `HARNESS_SOURCES`.
- [ ] **Step 2: Run them; they fail.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio.py -n 2 -q`. Expected: `ModuleNotFoundError: No module named 'tools.perf_ratio'`.
- [ ] **Step 3: Implement `tools/perf_ratio.py`, as in Interfaces.**
- [ ] **Step 4: Run the tests and lint.** Run: the Step 2 command, then `pixi run --no-install ruff check tools/perf_ratio.py tests/tools/test_perf_ratio.py`. Expected: `10 passed`; `All checks passed!`.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `tools/perf_ratio.py tests/tools/test_perf_ratio.py`.

### Task 5: Ratio history: ingest and render, and the record's owner

**Files:**
- Create: `tools/perf_ratio_history.py`
- Modify: `tools/perf_history_adapters.py:575-586` (one `NON_TIMING_RECORDS` entry)
- Test: `tests/tools/test_perf_ratio_history.py`

**Interfaces:**
- Consumes:
  - From `tools/perf_history_store.py`: `AppendOnlyViolationError`, `DuplicateHistoryEntryError`, `UnterminatedHistoryError` (`:36-46`) and `REPOSITORY` (`:32`). These are the same failure modes, so the classes are reused.
  - The store's `load_history(path: Path = HISTORY_PATH) -> History` (`:79`) and `append_entries(entries: Iterable[PerfEntry], path: Path = HISTORY_PATH) -> AppendResult` (`:120`) are typed on `PerfEntry` and do not fit. Its generic append-and-verify is the private `_append_text` (`:109`), which is not imported. This module carries its own append, with the same prefix check.
  - `NON_TIMING_RECORDS: Mapping[str, str]` (`tools/perf_history_adapters.py:575`), pinned by `test_every_committed_record_is_a_timing_family_or_a_listed_non_timing_record` (`tests/tools/test_perf_history.py:102`).
  - The rendered-record protocol, `("SCHEMA", "JSON_PATH", "MD_PATH", "load_record", "render_markdown")`, which the guard calls as `module.load_record(module.JSON_PATH)` (`tests/docs/test_included_records.py:52`, `:203`).
  - `sha256_file` (`tools/run_provenance.py:286`) and Task 1's `load_session`, `PerfRatioSession`, `Ratio` and `LoadOneMinute`.
- Produces:
  - Constants: `SCHEMA = "perf-ratio-history/1"` and `ENTRY_SCHEMA = "perf-ratio-entry/1"`.
  - Paths: `JSON_PATH = REPOSITORY / "benchmarks" / "results" / "perf_ratio_history.jsonl"`, `MD_PATH = JSON_PATH.with_suffix(".md")`.
  - Header lines: `GENERATED_BY = "<!-- Generated by tools/perf_ratio_history.py from benchmarks/results/perf_ratio_history.jsonl. Do not edit; run \`pixi run perf-ratio-render\`. -->"` and `SEMANTICS = "A ratio is a same-moment comparison at any load, never a budget verdict."`.
  - Formatting: `RATIO_DECIMALS = 3` and `SHORT_REV = 8`.
  - `@dataclass(frozen=True) class RatioEntry`, with fields:
    - `workload: str`, `stage: str`, `machine: str`;
    - `ref_commit: str`, `head_rev: str`, `head_dirty: bool`;
    - `ratio: Ratio`, `pair_ratio_median: Ratio`, `same_output: bool`;
    - `load_min: LoadOneMinute`, `load_max: LoadOneMinute`, `recorded_at: str`, `record_sha256: str`.

    Members: property `key -> tuple[str, str, str]` (workload, stage, machine); property `entry_sha256 -> str` (SHA-256 of the content's canonical `json.dumps(..., sort_keys=True)`); `to_json_line() -> str`; `classmethod from_json_line(text: str) -> RatioEntry`.
  - `entries_of(session: PerfRatioSession, record_sha256: str) -> tuple[RatioEntry, ...]`, which raises `InjectedRecordError`. Machine = `session.host`.
  - `load_record(path: Path = JSON_PATH) -> tuple[RatioEntry, ...]`.
  - `ingest(record_path: Path, history_path: Path = JSON_PATH) -> tuple[RatioEntry, ...]` (returns the entries added).
  - `render_markdown(entries: Sequence[RatioEntry]) -> str`. Lines: `GENERATED_BY`, a blank line, `SEMANTICS`, a blank line, then the table. Columns: `Workload | Reference | Machine | Date | HEAD | Ratio | Pair-ratio median | Same output | 1-min load`. Keys are sorted; within a key, newest first by `recorded_at`. Ratios print as `f"{r:.3f}×"` and HEAD as `` `<8>` `` + `+dirty`.
  - `write_markdown(history_path: Path = JSON_PATH, markdown_path: Path = MD_PATH) -> str`.
  - `main(argv: Optional[Sequence[str]] = None) -> int`, with subcommands `ingest RECORD...` (prints `appended N entries`) and `render`.
  - Errors: `InjectedRecordError`, `RatioEntryDigestMismatchError` (a line's stored `entry_sha256` does not re-derive), `MalformedRatioEntryError`.
- `NON_TIMING_RECORDS` gains `"perf_ratio_session.json": "same-moment HEAD/reference CPU ratios; owned by tools.perf_ratio_history, never the absolute history"`.

- [ ] **Step 1: Write the failing tests.** File: `tests/tools/test_perf_ratio_history.py`. Sessions are built from Task 1 types into `tmp_path`.
  - `test_ingest_appends_one_entry_per_workload_keyed_by_the_reference`: a 2-workload session gives 2 lines, with `stage == "ratio-vs-15fbed70"` and `key == (workload, stage, host)`.
  - `test_reingesting_a_record_appends_nothing`: a second `ingest` returns `()`, and the file bytes are unchanged.
  - `test_an_injected_record_is_refused_at_ingest`: `injected=Injection(HEAD, 1.0)` raises `InjectedRecordError`; the history file does not exist afterwards.
  - `test_an_edited_history_line_is_refused`: change `"ratio"` in a line; `load_record` raises `RatioEntryDigestMismatchError`.
  - `test_an_unterminated_history_is_refused`: strip the last newline; raises `UnterminatedHistoryError`.
  - `test_a_repeated_line_is_refused`: duplicate a line; raises `DuplicateHistoryEntryError`.
  - `test_render_states_the_semantics_and_lists_newest_first`:
    - Setup: two sessions for one workload, on different dates.
    - Assert the output starts with `GENERATED_BY + "\n\n" + SEMANTICS`, the newer date's row precedes the older one's, and the rendered ratio is `"1.234×"` for `Ratio(1.2341)`.
  - `test_render_marks_a_dirty_head_and_different_output`: `+dirty` appears in the HEAD cell; the same-output cell is `no`.
  - `test_the_module_follows_the_rendered_record_protocol`: `all(hasattr(perf_ratio_history, n) for n in ("SCHEMA", "JSON_PATH", "MD_PATH", "load_record", "render_markdown"))`.
  - `test_the_session_record_is_owned_outside_the_absolute_history`: `"perf_ratio_session.json" in NON_TIMING_RECORDS`, and `"perf_ratio_session.json" not in FAMILIES`.
- [ ] **Step 2: Run them; they fail.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_history.py -n 2 -q`. Expected: `ModuleNotFoundError: No module named 'tools.perf_ratio_history'`.
- [ ] **Step 3: Implement `tools/perf_ratio_history.py` and add the `NON_TIMING_RECORDS` entry.**
- [ ] **Step 4: Run the tests and lint.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_history.py tests/tools/test_perf_history.py -n 2 -q`, then `pixi run --no-install ruff check tools/perf_ratio_history.py tools/perf_history_adapters.py tests/tools/test_perf_ratio_history.py`. Expected: all passed, 0 failed; `All checks passed!`.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `tools/perf_ratio_history.py tools/perf_history_adapters.py tests/tools/test_perf_ratio_history.py`.

### Task 6: pixi tasks and typecheck line (C++ slot only)

**Files:**
- Modify: `pyproject.toml:260` (the end of the `types-benchmarks` string); add three tasks after `pyproject.toml:312` (`pipeline-timing-1-1-base`).

**Interfaces:**
- Consumes: the `pipeline-timing-1-1-base` task form (`pyproject.toml:312`): the pinned `editable_build_dir=...` prefix, then `SKBUILD_EDITABLE_SKIP="$editable_build_dir" python -m <module>`, with no `depends-on`. Also `tools.perf_ratio.main` (Task 4) and `tools.perf_ratio_history.main` (Task 5).
- Produces:
  - `perf-ratio`: `python -m tools.perf_ratio`. Description: "Time the default workload set interleaved on the pinned reference worktree and ratio-head (both clean, own envs and builds) as in-child process CPU; writes benchmarks/results/perf_ratio_session.json; pass --head-rev, --set gt-cap40, --workload, --inject-cpu-seconds ARM=S".
  - `perf-ratio-ingest`: `python -m tools.perf_ratio_history ingest`. Description: "Append a perf-ratio/1 record's entries to benchmarks/results/perf_ratio_history.jsonl; refuses injected records".
  - `perf-ratio-render`: `python -m tools.perf_ratio_history render`. Description: "Regenerate benchmarks/results/perf_ratio_history.md".
  - All three take the same pinned prefix and no `depends-on`.
  - Appended to `types-benchmarks`: ` tools/perf_ratio.py tools/perf_ratio_arms.py tools/perf_ratio_workloads.py tools/perf_ratio_record.py tools/perf_ratio_history.py`.

- [ ] **Step 1: Write the failing check.** Run: `pixi run --no-install perf-ratio --help`. Expected: fails with `perf-ratio: command not found` or a task-not-found error.
- [ ] **Step 2: Edit `pyproject.toml` as in Interfaces.** This is the only edit in the task.
- [ ] **Step 3: Verify the tasks and the lock.** Run: `pixi run --no-install perf-ratio --help`, then `pixi run --no-install perf-ratio-render --help`, then `git diff --quiet -- pixi.lock && echo lock-unchanged`. Expected: both print argparse usage (exit 0); `lock-unchanged`.
- [ ] **Step 4: Typecheck.** Run: `pixi run --no-install types-benchmarks`. Expected: `Success: no issues found in <n> source files`. Fix any finding in the owning module, not with an ignore.
- [ ] **Step 5: Run the lane's tests together.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/tools/test_perf_ratio_record.py tests/tools/test_perf_ratio_workloads.py tests/tools/test_perf_ratio_arms.py tests/tools/test_perf_ratio.py tests/tools/test_perf_ratio_history.py -n 2 -q`. Expected: `52 passed` (8 + 10 + 14 + 10 + 10).
- [ ] **Step 6: Hand to orchestrator.** pathspec = `pyproject.toml`.

### Task 7: A/A calibration run (orchestrator executes)

**Files:** writes `benchmarks/results/perf_ratio_session.json`; log in `build/rtr/aa-calibration.log`.

- [ ] **Step 1: Precheck the arms.** Run: `git -C /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-ratio-ref status --porcelain`, then `git -C /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-ratio-head status --porcelain`, then `uptime`. Expected: both empty; note the load.
- [ ] **Step 2: Run the default set, A/A, at whatever load.** The background launch carries `nohup` and a log. Run from the task5 root: `PYTHONUNBUFFERED=1 nohup pixi run --no-install perf-ratio --head-rev 15fbed70ebb016bdd8ef49453e04606d03e0b1fe > build/rtr/aa-calibration.log 2>&1`.
  - Launch it as a background job (it outlives the 10-minute foreground limit) and wait for its exit notification.
  - Recorded wall budget: the spec's estimate was about 30 min at load ~200; the first A/A run (2026-11-01, load 3-64) measured 956 s. The cap is 32 min (1912 s), twice the measured wall (user 2026-11-01). Past it, the run is a finding: stop and report; do not raise the cap silently.
- [ ] **Step 3: Check the evidence that passes.** Read `build/rtr/aa-calibration.log` and the record.
  - Exit 0. The log holds `rebuilt` lines for both arms, 10 summary lines, and `wall <s> s`.
  - `load_session("benchmarks/results/perf_ratio_session.json")` builds: 10 workloads × 5 pairs, `injected is None`.
  - `arms[ref].rev == arms[head].rev == 15fbed70ebb016bdd8ef49453e04606d03e0b1fe`, both `dirty false`.
  - `same_output` is true on all 10 workloads. A false one is a cross-worktree determinism finding: report it, don't ingest.
  - Both arms' `native_sources_sha256` are equal (expected `09aacb36…`). `native_sha256` is recorded per arm and is expected to differ (path-dependent `__FILE__`); it is not compared.
  - Noted for Task 8 and the docs:
    - per workload: ratio, pair-ratio median, and spread `τ_w = max_i |pair_ratio_i − pair_ratio_median|`;
    - the wall seconds and the 1-minute load range.
- [ ] **Step 4: Hand to orchestrator.** Hold the record for Task 9's commit. Nothing is committed here.

### Task 8: Injection audit, both directions (orchestrator executes)

Case: `gt/cap-120-default/pipeline/rect_20x12` (the costliest cap-120 GT case, 1.41 s wall in recon). `S = 1.0` CPU seconds.

- [ ] **Step 1: Precondition from the A/A record.** For this workload, `S / ref_min_aa > τ_w` (the shift exceeds the measured noise). If it fails, raise S until it holds, record why, and use the same S both ways.
- [ ] **Step 2: Inject into HEAD.** Run (timeout 600000): `PYTHONUNBUFFERED=1 pixi run --no-install perf-ratio --head-rev 15fbed70ebb016bdd8ef49453e04606d03e0b1fe --workload gt/cap-120-default/pipeline/rect_20x12 --inject-cpu-seconds head=1.0 --write build/rtr/audit-head.json > build/rtr/audit-head.log 2>&1`.
  - Passes when `|R_up − (R_aa + S / ref_min_up)| ≤ τ_w`. Here `R_up` and `ref_min_up` come from `audit-head.json` (its uninjected ref arm), and `R_aa` and `τ_w` from the A/A record.
- [ ] **Step 3: Inject into the reference.** Run: the same command with `--inject-cpu-seconds ref=1.0 --write build/rtr/audit-ref.json > build/rtr/audit-ref.log`.
  - Passes when `|R_down − head_min_dn / (head_min_dn / R_aa + S)| ≤ τ_w`, with `R_down < R_aa`. `head_min_dn` comes from `audit-ref.json` (its uninjected head arm).
- [ ] **Step 4: The injected records are refused at ingest.** Run: `pixi run --no-install perf-ratio-ingest build/rtr/audit-head.json; echo "exit $?"`, then `sha256sum`/`shasum -a 256` of `benchmarks/results/perf_ratio_history.jsonl` before and after (or "absent" both times). Expected: the traceback names `InjectedRecordError`, `exit 1`, and the history is unchanged.
- [ ] **Step 5: Hand to orchestrator.** Write the audit numbers (`S`, `R_aa`, `τ_w`, `R_up`, `R_down`, their expected values, and the loads) to `build/rtr/audit.md`. Nothing is committed here; the numbers feed Task 9's docs paragraph.

### Task 9: Ingest, render, docs paragraph (orchestrator executes)

**Files:**
- Create (generated): `benchmarks/results/perf_ratio_history.jsonl`, `benchmarks/results/perf_ratio_history.md`
- Modify: `docs/performance.md:362-363` (a new section after the `perf_history.md` include, before `## Related pages`)

- [ ] **Step 1: Ingest the A/A record.** Run: `pixi run --no-install perf-ratio-ingest benchmarks/results/perf_ratio_session.json`, twice. Expected: first `appended 10 entries`, then `appended 0 entries`.
- [ ] **Step 2: Render.** Run: `pixi run --no-install perf-ratio-render`. Expected: `benchmarks/results/perf_ratio_history.md` starts with the `GENERATED_BY` line, then the semantics line, then 10 rows.
- [ ] **Step 3: Write the docs paragraph.** The new section is `## Ratio tracking against a pinned reference`, after `--8<-- "benchmarks/results/perf_history.md"`. It holds a short paragraph, then `--8<-- "benchmarks/results/perf_ratio_history.md"`. The paragraph covers:
  - `pixi run perf-ratio` times HEAD (ratio-head at `--head-rev`) against the reference at `15fbed70`. The pairs are interleaved AB/BA, the CPU is in-child process CPU, and the ratio is min-of-5 per arm, with the pair-ratio median beside it.
  - A ratio is a same-moment comparison at any load, never a budget verdict. It differs from the D2 quiet-only CPU gates ("Ruled: the generation budget gates CPU time") because load slows both arms at the same moment. That is why it keeps its own history, and it never enters `perf-history-check`'s baseline.
  - The measured A/A noise floor (ratio range and the largest `τ_w`), the wall budget of the default set at its recorded load, and the audit result from Task 8 (S, the expected and measured shifts, both directions).
  - v1 has no ratio gate; a threshold is ruled later, from this A/A data.
- [ ] **Step 4: Guard.** Run: `pixi run --no-install build/pyrun.sh -m pytest tests/docs/test_included_records.py tests/tools/test_perf_history.py tests/tools/test_perf_ratio_history.py -n 2 -q`. Expected: all passed. `test_record_is_the_render_of_its_json[benchmarks/results/perf_ratio_history.md]` is among them, discovered through the protocol with no edit to the guard. Then run `pixi run --no-install -e docs docs`. Expected: the strict build succeeds. If the `docs` env is not installed, say so in the report; don't install it.
- [ ] **Step 5: Hand to orchestrator.** pathspec = `benchmarks/results/perf_ratio_session.json benchmarks/results/perf_ratio_history.jsonl benchmarks/results/perf_ratio_history.md docs/performance.md`.

## Self-review

- **Spec coverage:**
  - Goals: G1 → Tasks 2, 3, 4 and 6; G2 → Tasks 1, 4 and 5; G3 → Tasks 5 and 9; G4 → Task 8.
  - Anti-goals:
    - one interpreter/build → Task 3 `SameArmError` and the foreign-library tests;
    - wall clock or launcher CPU → Task 2 (in-child clocks only);
    - ratios in the absolute machinery → Task 5 (separate history, `NON_TIMING_RECORDS`);
    - sleep as the defect → Task 2's spin on `process_time`.
  - Arms:
    - clean before and after → Task 4's `main` order;
    - `--head-rev` checkout → Task 3;
    - settled builds → `require_settled`;
    - `SKBUILD_EDITABLE_SKIP` pinned → `base_child_env`;
    - symmetric proofs → `bracket_arm`;
    - moving the reference → `--ref-worktree`/`--ref-rev` give a new `stage`.
  - Determinism and `same_output` → Tasks 1 and 4. `gt-cap40` → Task 2. Run shape → Task 4. The unit-test list in the spec → Tasks 3–5.
- **Types:** `ArmName`, `Workload`, `ChildRun`, `Arm`, `ArmBracket`, `Injection` and `PerfRatioSession` are defined once (Tasks 1–3) and used under the same names in Tasks 4–5.
- **Review Focus:** each of the five lines has its test in the owning task.
- **Proportion:** the code bodies are only the two the spec pins (the GT snippet, `pair_order`).
