# Auditor Convergence P4 Enforcement and Evidence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the truthful audit/replay contracts unavoidable in Python,
native builds, CI, and published measurement artifacts.

**Architecture:** Align every declared interpreter/type target at Python 3.12,
embed a deterministic native/Python source manifest into production
`BuildIdentity`, and replace hollow CI/task entry points with Pixi-owned gates.
Generate one schema-validated evidence bundle from authoritative audit/replay
objects and derive all tables/figures from that raw artifact.

**Tech Stack:** Python 3.12, C++20, CMake, CGAL, nanobind, Pixi, GitHub Actions,
CCAN, SHA-256, JSON Schema, pytest-xdist, strict mypy, Ruff, MkDocs.

**Spec:** `docs/superpowers/specs/2026-08-23-auditor-convergence-design.md`

## Global Constraints

- P2 acceptance is mandatory; P3 may execute independently but its final
  evidence is included only after its acceptance commit.
- Pixi exclusively owns local and CI dependency/build/task execution.
- Python is exactly 3.12 across metadata, Ruff, mypy, Pixi, CI, and wheel ABI.
- Type errors are fixed through domain narrowing, never casts or new ignores.
- Production constructs and exposes the real `BuildIdentity`.
- Every JSON artifact round-trips through `json.loads` and schema validation.
- Every figure and prose number derives from committed raw data.
- Every pytest command uses `-n auto`; zero collected tests is failure.
- No main/master mutation, force flag, skipped test, conditional import, or
  fallback path.

---

### Task 1: Align Python 3.12 and make strict typing green

**Files:**

- Modify: `pyproject.toml`
- Modify: `src/compas_cgal/adaptive/generator.py`
- Modify: `src/compas_cgal/adaptive/replay.py`
- Modify: `src/compas_cgal/adaptive/bootstrap.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`
- Create: `tests/test_python_floor.py`

**Interfaces:**

- Consumes: current six-error mypy baseline.
- Produces: coherent Python 3.12 metadata and zero-error strict gates.

- [ ] **Step 1: Add RED metadata consistency tests**

```python
def test_declared_python_floor_matches_runtime_and_wheel_abi() -> None:
    config = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))

    assert config["project"]["requires-python"] == ">=3.12"
    assert config["tool"]["ruff"]["target-version"] == "py312"
    assert config["tool"]["scikit-build"]["wheel"]["py-api"] == "cp312"
    assert config["tool"]["pixi"]["dependencies"]["python"] == "3.12.*"
```

Use `pathlib.Path`, never `os.path`.

- [ ] **Step 2: Run RED and strict baseline**

```bash
pixi run pytest -- tests/test_python_floor.py -n auto -q
pixi run types-adaptive
```

Expected: metadata test fails and mypy reports the six recorded errors.

- [ ] **Step 3: Align metadata**

Set `requires-python = ">=3.12"`, Ruff `py312`, retain Pixi `3.12.*`, and keep
wheel ABI `cp312`. Update classifiers and CI matrices that encode lower
versions. Do not add `typing_extensions`.

- [ ] **Step 4: Fix union narrowing and obsolete ignores**

At generator lines identified by P0, narrow `int | None` with the existing
route-state invariant before construction; dispatch `ExactSegmentMotion |
ExactCircleMotion` by exact variant before reading endpoints; dispatch
`TraversalCommit | RouteRetraceCommit` before reading source digest. Remove the
two unused ignores. Each repair gets a focused runtime test that would raise
the latent `AttributeError` before the fix.

- [ ] **Step 5: Run GREEN and commit**

```bash
pixi run ruff format src/compas_cgal tests
pixi run lint
pixi run types-adaptive
pixi run pytest -- tests/test_python_floor.py tests/adaptive/typecheck tests/adaptive/test_generator.py -n auto --testmon -q
git diff --check
git add pyproject.toml src/compas_cgal/adaptive tests/adaptive/typecheck tests/test_python_floor.py
git commit -m "fix(types): enforce python 3.12"
```

### Task 2: Generate and embed production build identity

**Files:**

- Create: `tools/build_identity_manifest.py`
- Create: `tests/test_build_identity_manifest.py`
- Create: `src/build_identity.h`
- Create: `src/build_identity.cpp`
- Modify: `CMakeLists.txt`
- Modify: `src/stock_2.cpp`
- Modify: `src/continuous_tea_2/continuous_tea_2.cpp`
- Modify: `src/compas_cgal/_stock_2.pyi`
- Modify: `src/compas_cgal/_continuous_tea_2.pyi`
- Modify: `src/compas_cgal/engagement_audit/identity.py`
- Modify: `pyproject.toml`

**Interfaces:**

- Consumes: exact source manifest, CMake target configuration, dependency
  versions, `pixi.lock`.
- Produces: native `BuildManifestV1`, `native_build_identity()`,
  `BuildIdentity.from_runtime()`.

- [ ] **Step 1: Write RED manifest determinism tests**

```python
def test_manifest_is_order_independent_and_source_sensitive(tmp_path: Path) -> None:
    first = build_manifest(_fixture_sources(tmp_path, order="forward"))
    second = build_manifest(_fixture_sources(tmp_path, order="reverse"))

    assert first.canonical_bytes == second.canonical_bytes
    _mutate_source(tmp_path / "a.cpp")
    assert build_manifest(_fixture_sources(tmp_path, order="forward")).digest != first.digest
```

Add tests for compile flags, compiler, CGAL, CORE/GMP, nanobind, Python ABI,
component versions, CMake inputs, and lock digest.

- [ ] **Step 2: Run RED**

```bash
pixi run pytest -- tests/test_build_identity_manifest.py -n auto -q
```

- [ ] **Step 3: Implement deterministic manifest generation**

Use `pathlib.Path`, sorted repository-relative POSIX paths, raw file bytes, and
length-delimited canonical encoding. The module has one CLI invoked through a
Pixi task; it writes only into the CMake build directory. Source selection is
an explicit checked list derived from target sources, not a broad filesystem
glob.

- [ ] **Step 4: Embed and expose native identity**

CMake invokes the Pixi-owned manifest producer before compiling exact audit
targets and generates one header in the build directory. C++ verifies digest
length and exposes canonical bytes/digest through both native modules.
Bindings return immutable bytes; no environment or Git command is consulted at
runtime.

- [ ] **Step 5: Construct production `BuildIdentity`**

`BuildIdentity.from_runtime()` compares both native modules' manifest bytes,
computes the Python-source and lock digests from installed/checked sources, and
builds the existing Stage 1 type. A mismatch raises
`NativeBuildIdentityMismatchError`; callers cannot provide a replacement
digest to the public audit factory.

- [ ] **Step 6: Run source-mutation and clean-build gates**

```bash
pixi run pytest -- tests/test_build_identity_manifest.py tests/engagement_audit/test_identity.py tests/engagement_audit/test_input.py -n auto --testmon -q
pixi run baseline
git diff --check
```

The mutation test changes a copied fixture source, never the live repository.

- [ ] **Step 7: Commit**

```bash
git add tools/build_identity_manifest.py tests/test_build_identity_manifest.py src/build_identity.h src/build_identity.cpp CMakeLists.txt src/stock_2.cpp src/continuous_tea_2/continuous_tea_2.cpp src/compas_cgal/_stock_2.pyi src/compas_cgal/_continuous_tea_2.pyi src/compas_cgal/engagement_audit/identity.py pyproject.toml
git commit -m "feat(build): embed source identity"
```

### Task 3: Replace hollow tasks with Pixi-owned enforcement

**Files:**

- Modify: `pyproject.toml`
- Create: `tests/adaptive/test_schema.py`
- Create: `tests/adaptive/test_mutation_contract.py`
- Create: `tests/test_pixi_task_contract.py`
- Modify: `.github/workflows/build.yml`
- Modify: `.github/workflows/benchmarks.yml`
- Modify: `.github/workflows/docs.yml`

**Interfaces:**

- Consumes: P1/P2 test surfaces and build identity.
- Produces: functional `schema`, `mutations-adaptive`, `regression`,
  `benchmark-instrument`, and `quality-gate` Pixi tasks plus CI jobs.

- [ ] **Step 1: Write RED task-contract tests**

Parse `pyproject.toml` and workflows. Require every referenced test/script to
exist, benchmark CLI grammar to parse, pytest task commands to contain
`-n auto`, and no production CI command to invoke conda, pip, pipx, poetry, or
venv.

- [ ] **Step 2: Implement real schema and mutation gates**

`test_schema.py` loads every checked JSON schema and validates one positive and
one malformed artifact. `test_mutation_contract.py` invokes the existing
certificate mutation parametrizations through pytest selection. Point the
Pixi tasks at these real tests; do not create a wrapper shell script.

- [ ] **Step 3: Define green regression and red product tasks**

`regression` runs all proof/instrument tests and excludes only the explicitly
separate file path owned by `quality-gate`; no pytest skip marker or keyword
filter is used. `quality-gate` runs the complete product-quality file.
Every pytest command includes `-n auto` and a collection-count assertion test.

- [ ] **Step 4: Migrate CI to Pixi**

Use `prefix-dev/setup-pixi` and named Pixi tasks on Linux/macOS/Windows only
where the lock declares that platform. Remove conda/pipx execution only after
the equivalent Pixi job passes on the branch. Triggers cover all pull requests
and branch pushes. Add concurrency keyed by workflow and ref with
`cancel-in-progress: true`.

- [ ] **Step 5: Run local workflow/task contracts**

```bash
pixi run pytest -- tests/test_pixi_task_contract.py tests/adaptive/test_schema.py tests/adaptive/test_mutation_contract.py -n auto -q
pixi run schema
pixi run mutations-adaptive
pixi run regression
pixi run benchmark-instrument
git diff --check
```

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml tests/adaptive/test_schema.py tests/adaptive/test_mutation_contract.py tests/test_pixi_task_contract.py .github/workflows
git commit -m "ci: enforce proof gates"
```

### Task 4: Generate the content-addressed evidence bundle

**Files:**

- Create: `benchmarks/evidence.py`
- Create: `benchmarks/schemas/auditor-evidence-v1.json`
- Create: `tests/benchmarks/test_evidence.py`
- Modify: `benchmarks/cli.py`
- Modify: `benchmarks/report.py`
- Modify: `benchmarks/figures.py`
- Create: `benchmarks/data/held-figure6-v1.json`
- Create: `docs/benchmarks/auditor-evidence-v1.json`
- Create: `docs/benchmarks/auditor-evidence-v1.md`
- Modify: `docs/benchmarks.md`
- Modify: `docs/engagement_controlled_toolpath.md`

**Interfaces:**

- Consumes: authoritative audit report, replay certificate, generator override
  ledger, raw corpus measurements, BuildIdentity.
- Produces: `AuditorEvidenceBundle.build(...)`, JSON/Markdown/figures from one
  checked raw artifact.

- [ ] **Step 1: Write RED schema and round-trip tests**

```python
def test_evidence_bundle_round_trips_schema_and_digest(tmp_path: Path) -> None:
    bundle = _evidence_bundle()
    path = write_evidence_bundle(bundle, tmp_path)
    decoded = json.loads(path.read_text(encoding="utf-8"))

    jsonschema.validate(decoded, _evidence_schema())
    assert bytes.fromhex(decoded["digest"]) == bundle.digest
    assert read_evidence_bundle(path).canonical_bytes == bundle.canonical_bytes
```

Mutate every identity/verdict/runtime/override/replay field and require schema
or digest rejection.

- [ ] **Step 2: Implement the evidence model and CLI**

The bundle binds canonical input/build/audit/replay/generator/override digests,
four verdict counts, unresolved third-party rate, raw per-run timings, and
protocol versions. Wall-clock samples remain data; deterministic canonical
identity excludes host timestamp/path noise.

Expose one Pixi task:

```toml
evidence = { cmd = "python -m benchmarks.cli evidence --out docs/benchmarks", depends-on = ["_editable-rebuild"] }
```

- [ ] **Step 3: Commit Held digitisation provenance before figures**

The JSON records source publication, figure/table locator, axis transforms,
digitiser/version, point series, and SHA-256. Figure code reads this file and
the auditor evidence bundle; no Held series literal remains in plotting code.

- [ ] **Step 4: Generate fresh artifacts**

```bash
pixi run evidence
pixi run figures
pixi run pytest -- tests/benchmarks/test_evidence.py tests/benchmarks/test_figures.py -n auto -q
```

Record exact runtime and unresolved rates. Generated JSON is never hand-edited.

- [ ] **Step 5: Eliminate stale prose numbers**

Render headline tables from the evidence bundle and add tests that extract
documented metrics and compare them to JSON. Remove unsupported prose claims;
state measured-not-closed and missing third-party evidence where applicable.

- [ ] **Step 6: Commit**

```bash
git add benchmarks/evidence.py benchmarks/schemas benchmarks/data tests/benchmarks/test_evidence.py benchmarks/cli.py benchmarks/report.py benchmarks/figures.py docs/benchmarks docs/benchmarks.md docs/engagement_controlled_toolpath.md
git commit -m "feat(bench): publish auditor evidence"
```

### Task 5: Final whole-branch gate and release classification

**Files:**

- Modify: `docs/auditor_convergence.md`
- Create: `docs/superpowers/state/2026-08-23-auditor-convergence-final.md`
- Modify: `.superpowers/sdd/2026-08-23-auditor-convergence/progress.md`

**Interfaces:**

- Consumes: P0-P4 acceptance commits and generated evidence.
- Produces: evidence-backed independent classifications for auditor, generator,
  and certifier theses.

- [ ] **Step 1: Freeze and record the final tree**

```bash
git rev-parse HEAD
git status --short --branch --untracked-files=all
```

- [ ] **Step 2: Run all green gates from the frozen tree**

```bash
pixi run regression
pixi run schema
pixi run mutations-adaptive
pixi run benchmark-instrument
pixi run lint
pixi run types-adaptive
pixi run -e docs docs
git diff --check
```

- [ ] **Step 3: Run and record the product gate**

```bash
pixi run quality-gate
```

Its exit code does not change the green regression result; it determines the
regulated-generator classification.

- [ ] **Step 4: Write the final evidence ledger**

Classify quality instrument, regulated generator, and exact certifier
independently as `release candidate`, `research prototype`, or `rejected at
current evidence`. Every classification links to raw artifact fields and test
commands. No thesis borrows evidence from another.

- [ ] **Step 5: Commit final documentation**

```bash
git add docs/auditor_convergence.md docs/superpowers/state/2026-08-23-auditor-convergence-final.md .superpowers/sdd/2026-08-23-auditor-convergence/progress.md
git commit -m "docs: classify auditor evidence"
```

- [ ] **Step 6: Invoke branch-finishing workflow**

Use `superpowers:finishing-a-development-branch`. Present verified integration
options; do not merge, rebase into, checkout, or push main/master. Original
source branches/worktrees remain intact until the user separately authorizes
cleanup.
