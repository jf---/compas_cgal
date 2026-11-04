# Ratio gate: CI wiring analysis

Status: ANALYSIS ONLY, 2026-11-01, at `9ae1cdd2` (`codex/held-task5-continuation`). No workflow file is edited. Everything under "Proposed change set" is PROPOSED, not applied. The project's CI-config gate applies: this page is the analysis shown before any workflow edit.

## Conclusion

1. **Not on every PR.** Run the gate nightly on `main`, on `workflow_dispatch`, and on PRs that carry a `perf-gate` label. One gate session is about 14–16 min of measurement (863 s audit A/A, `docs/performance.md:473`; the spec gives ~16 min, `spec:92`). On top of that come two native builds, and up to 3 extension harness runs per red workload. That is too much for every push, and the existing `benchmarks` job already estimates 35–45 min (`benchmarks.yml:50-54`).
2. **Runner class is the decision for the user.** The gate refuses an uncalibrated `cpu_model` (`UncalibratedMachineError`, `perf_ratio_gate_registry.py:120`; spec open question 1 resolved to refuse). GitHub-hosted `ubuntu-latest` does not promise a CPU model, so every new model a job lands on is a refusal until someone calibrates it. Each calibration costs K=5 A/A sessions, 833–860 s each on M1 Max, about 70 min (`docs/performance.md:381-382`). The only calibrated machine today is `Apple M1 Max` (`benchmarks/results/perf_ratio_gate.json`). The two sound options are:
   - a **self-hosted runner on the calibrated M1 Max**: zero new calibration, and it is the measured machine;
   - **GitHub-hosted runners with a refusal-tolerant job**: an uncalibrated model gives a distinct "refused" exit, not a red, and its A/A records are uploaded for local calibration.
3. **Exit codes need splitting first (N3).** Today `check` exits 1 for a regression, for every named refusal, and for any infra failure alike. Only argparse usage errors differ (exit 2). This was probed, see §4. CI cannot tell "HEAD is slower" from "runner uncalibrated" from "the build broke" until `main` maps exception classes to exit codes.
4. **Promotion and calibration never run automatically in CI.** CI may produce the records and upload them as artifacts. A human downloads, verifies, runs `perf-ratio-promote` / `perf-ratio-calibrate --records` locally, and commits. A CI bot commit would also violate the authorship rule.
5. **Nothing can be exercised on GitHub yet.** `origin/main..HEAD` is 662 commits, with 45 lines matching the identity grep in author, committer or message. The identity rewrite must land before any push. Until then, validate locally: the job script as a pixi task (dry run), and `act` for YAML syntax only.
6. **Order of landing:** anchor chain (schema `/2`, plan Tasks 1–4) → exit-code classes → job script as a pixi task → workflow. The committed registry is still `perf-ratio-gate/1` and will be refused once `/2` lands (plan:14), so wiring CI against `/1` is wasted work.

## 1. Workflows, runners, cpu_model, calibration

| Workflow | Trigger | Runner(s) | Concurrency | Relevance |
|---|---|---|---|---|
| `benchmarks.yml` | push `main`/`master`/`codex/**`, PR to `main`/`master` (`:8-12`) | `ubuntu-latest` (`:22`) | `benchmarks-${{ github.ref }}`, cancel-in-progress (`:14-17`) | pixi via `prefix-dev/setup-pixi@v0.8.1`, `cache: true` (`:31-35`); natural host for the new job's setup steps |
| `build.yml` | push/PR `main` (`:3-9`) | ubuntu/macos/windows-latest matrix (`:17`) | none | conda + cibuildwheel, no tests (`benchmarks.yml:3-6`) |
| `docs.yml` | push `main` | `ubuntu-latest` | none | `contents: write`, gh-deploy |
| `pr-checks.yml` | PR events incl. `labeled`/`unlabeled` (`:4`) | `ubuntu-latest` | none | already reacts to labels: precedent for a label trigger |
| `release.yml` | tags `v*` | matrix + PyPI | none | unrelated |

None of these runs any `perf-ratio*` task.

`benchmarks.yml` has never run on GitHub: `gh run list` answers 404, "not found on the default branch". So no hosted-runner timing exists for this branch's native build. The only measured hosted build is `main`'s `build` run 30216806200 (2026-07-26). There, the `Build (ubuntu-latest)` job (conda env + build + import check) took 5 min 44 s, and `cibuildwheel on manylinux` took 13 min 44 s. `main` lacks the exact-kernel modules (`_stock_2`, `_coverage_2`, …), so treat these as lower bounds, not estimates.

**cpu_model source.** On Linux, the gate reads the first `model name` line of `/proc/cpuinfo` (`tools/run_provenance.py:356-358`). On darwin it reads `sysctl machdep.cpu.brand_string` (`:354`). Consequences:

- `ubuntu-latest` (x64) runners are a pooled VM fleet; the CPU model is not part of GitHub's runner contract. Not measured here. A probe job is needed (proposed below) that logs `model name` over about 20 runs before choosing.
- Arm64 Linux runners usually expose no `model name` line, so `machine_identity` raises `MachineIdentityError` (`run_provenance.py:359-360`). That is an infra failure, not a refusal. Exclude them.
- macOS hosted runners report a virtualised brand string. That is also unmeasured, and the calibrated slug `apple-m1-max` would not match it.
- The calibration is keyed by `cpu_model` only (`calibration_dir`, `perf_ratio_gate_registry.py:235-237`), not by core count or VM size. Two runner SKUs that share a CPU model share margins. That is acceptable by the spec's argument (same-moment ratio, `spec:22-26`), but it is untested on VMs.

**Calibration cost per runner class.** K = 5 A/A sessions at the registry reference (`CALIBRATION_SESSIONS`; `calibrate --sessions K`, `perf_ratio_gate.py:293-375`). Measured only on M1 Max: 833–860 s per session (`docs/performance.md:381-382`), so about 70–72 min of measurement plus two builds. The workloads run one child at a time, so per-core speed sets the hosted figure, not core count. It is unmeasured, and the 6 h job limit is the ceiling. Calibration must restart after every promotion and every harness edit (`perf_ratio_gate.py:14`, `docs/performance.md:455-467`). Each calibrated class therefore costs about 70 min per promotion.

**Where calibration records are committed from.** `calibrate --sessions` stages the records under `build/ratio-gate/calibration-staging-*` and renames them into `benchmarks/results/perf_ratio_calibration/<ref8>/<cpu-slug>/aa-<k>/perf_ratio_session.json` (`perf_ratio_gate.py:362-374`). The registry cites them by `record_sha256` and re-hashes them on load. In CI, that directory exists only in the runner workspace. Flow:
1. upload it as an artifact;
2. download it locally into the same path;
3. run `pixi run perf-ratio-calibrate --records …/aa-*/perf_ratio_session.json`, which re-validates;
4. commit records + registry by hand.

## 2. Provisioning the two arms

The harness hard-codes sibling worktrees: `ratio-ref = ROOT.parent/"compas_cgal_prs-ratio-ref"` and `ratio-head = ROOT.parent/"compas_cgal_prs-ratio-head"` (`tools/perf_ratio_arms.py:49-51`). On a runner, `ROOT` is `$GITHUB_WORKSPACE`, so these become sibling directories under `/home/runner/work/<repo>/`. That is writable, so no harness change is needed.

Per run the job needs:
- full history (`fetch-depth: 0`, since `ref_commit` is an old commit);
- `git worktree add --detach ../compas_cgal_prs-ratio-ref <ref_commit>`;
- `git worktree add --detach ../compas_cgal_prs-ratio-head HEAD`;
- `pixi install --manifest-path <arm>/pyproject.toml` for each arm, because the arm interpreter is `pixi run --no-install --manifest-path …` (`perf_ratio_arms.py:104-110`) and never installs.

That makes three pixi environments in total, orchestrator plus two arms. They must stay distinct: shared prefixes are refused (`SameArmError`, `:159-175`). In CI they are distinct by default, because nothing detaches the environments.

**Native build twice per run.** The harness rebuilds both arms with `_editable-rebuild` (`perf_ratio.py:413-415`, `pipeline_timing_1_1_base.py:111,225-236`) and then settles each native build (`:416-417`). Every extension is a fresh `python -m tools.perf_ratio` subprocess (`perf_ratio_gate.py:105-155`), and it repeats rebuild + settle. A warm editable build is a no-op compile, but its cost on a runner is unmeasured. The orchestrator's own build is pinned and never rebuilt (`perf_ratio.py:412`).

**Caching.**
- pixi envs: `setup-pixi` `cache: true` keys on `pixi.lock`. For the arms, use `actions/cache` keyed on `hashFiles('<arm>/pixi.lock')`. The ref arm's lock is fixed by `ref_commit`, so its env cache key is stable until a promotion.
- Ref arm native build: cacheable by commit. Key it on `ref_commit` + the ref lock hash + `runner.os` + the `cpu_model` slug + compiler version. Cache the arm's `build/{wheel_tag}` (`pyproject.toml:72`). It must be restored with its timestamps intact, so that `editable.rebuild = true` (`:76`) sees it as up to date. If the restored build is stale or differs, `require_settled` refuses (`NativeIdentityError`, `perf_ratio_arms.py:211-217`), so a bad cache fails loudly and never gives a silent wrong ratio.
- Head arm: key on `github.sha`, which rarely hits. Use ccache/sccache-style compiler caching instead. That needs the CMake launcher configured, so it is a build-config change and not in scope here.

## 3. Trigger policy

| Policy | Cost | Coverage | Verdict |
|---|---|---|---|
| Every PR push | ≥ 16 min measurement + 2 builds, per push | full | rejected: duplicates `benchmarks` cost, and the M1 Max calibration does not apply on hosted runners |
| Label `perf-gate` on PR | same, per opt-in | the PRs that touch hot paths | keep |
| Nightly on `main` | one session per day | catches merged regressions within 24 h | keep, primary |
| `main` push only | per merge | per-merge attribution | optional; nightly + label covers it |
| `workflow_dispatch` | on demand | calibration runs, re-checks | keep (calibration lives here) |

**Recommendation:** nightly `schedule` on `main` + `workflow_dispatch` (inputs: `mode = check | calibrate`) + PR `labeled` with `perf-gate`. Use a concurrency group per runner class with `cancel-in-progress: false`: two gate sessions on one machine would load each other, and a cancelled session wastes its 16 min. Set `timeout-minutes` from the ceiling, not a guess:
- one session at 15 pairs over 10 workloads, plus 3 extensions × 10 workloads at worst, plus builds;
- check mode: 120 min;
- calibrate mode: 180 min.

Revisit both once hosted timings exist.

## 4. Exit-code classes (N3)

N3 does not appear under that label in `.superpowers/sdd/2026-11-01-ratio-gate/deferred-and-rulings.md` at HEAD. It is handled here as the CI-unit item it was briefed as.

**Current behaviour (code + probe).** `main` returns 0 and "every refusal, and a regression, is raised" (`perf_ratio_gate.py:470-471`). `__main__` does `raise SystemExit(main())` (`:496-497`), so any exception is an uncaught traceback with exit 1. Probed 2026-11-01 (logs in `build/ratio-gate/ci-analysis/`):
- `perf-ratio-check <record> --registry <absent>` gives exit 1, `MissingRegistryError`;
- `perf-ratio-check --bogus` gives exit 2 (argparse).

`RatioRegressionError` (`perf_ratio_gate_verdict.py:137`, raised `perf_ratio_gate.py:227-228`) also exits 1. Not probed live: a red session costs 16 min. That conclusion comes from code reading.

**Classes CI must distinguish:**

| Class | Meaning | Examples (named) | CI reaction |
|---|---|---|---|
| 0 pass | no regression | PASS / IMPROVED | green; when IMPROVED with no refusals, annotate "promotion available" (`:225-226`) |
| regression | HEAD slower beyond margin after extension | `RatioRegressionError` | red, block the label/merge |
| refusal | the gate declined to decide | `UncalibratedMachineError`, `ReferenceMismatchError`, `HarnessMismatchError`, `InjectedSessionError`, `UnregisteredWorkloadError`, `IncompleteSessionError`, `ExtensionMismatchError`, `MissingRegistryError`, registry load errors (`RegistryDigestMismatchError`, `MarginMismatchError`, `CalibrationRecord*Error`, …; registry `:60-124`), `AnchorChainMismatchError` (after the plan) | neutral/warning plus instructions. Uncalibrated is expected on a new runner; digest/chain errors are a real repo defect |
| infra | could not measure | `subprocess.CalledProcessError` from an extension (`:152`), harness-step failure (arm build `BaseRebuildFailedError`, `MeasuredTreeError`, `MachineIdentityError`, `NativeIdentityError`), `OSError` | red, "rerun once" per the CI rule; second failure means root-cause it |
| usage | bad CLI | argparse | exit 2 already |

**Gap.** Regression, refusal and infra all exit 1. Proposed fix, owned by whoever owns `tools/perf_ratio_gate.py` (another agent is editing it now, so it is not touched here):
- named exit constants, e.g. `EXIT_REGRESSION = 1`, `EXIT_REFUSED = 3`, `EXIT_INFRA = 4`;
- one `except` ladder in `main` that maps the closed set of gate exception classes to those constants and re-raises everything else, so that unknown exceptions stay loud as infra;
- a base class `GateRefusal` for the refusal family, so the mapping is a type and not a list. That follows engineering rule 7 and keeps one named exception per mode;
- tests: one per class, plus R6 (map a refusal to 1 → a named test turns red).

The harness step (`perf-ratio`) runs as its own CI step, so any non-zero there is infra by construction.

## 5. Promotion in CI

Never automatic. Reasons:
- promotion moves `ratio-ref` and rewrites the registry (`perf_ratio_gate.py:378-418`);
- after the anchor chain lands, it also writes the gate + extension session records under `benchmarks/results/perf_ratio_promotions/<from8>-<to8>/{gate,ext-<n>}/perf_ratio_session.json` (spec:65, plan:11,38);
- those records are content-addressed evidence the registry cites, and they must be committed with the registry, under the user's identity;
- the reference only ratchets on a human decision, and every promotion voids all calibrations (`:14`), which costs about 70 min per runner class.

CI's role:
- on IMPROVED with no promotion refusals, upload the gate record **and every extension record** (`build/ratio-gate/ext-*.json`, `perf_ratio_gate.py:138-151`) as an artifact, and annotate the run with the `pixi run perf-ratio-promote <record>` line;
- the maintainer downloads the records, re-runs promote locally, which re-gates them, then calibrates and commits.

Open point: promote re-gates through a live `Measure`. If a red workload needs extension at promote time, it is measured on the promoting machine, not on the CI runner. After the chain lands, that extension would be refused as `ExtensionMismatchError` (another host/CPU). Promote should therefore run on the same runner class as the gate session: for a hosted runner, a `workflow_dispatch mode=promote-dry` that produces the link records as an artifact. Decide this with the runner-class choice.

## 6. Identity blocker and local validation

Measured 2026-11-01:
- `git log origin/main..HEAD` has 662 commits;
- the author/committer/body identity grep gives 45 matching lines.

Per the global rule, nothing on `codex/*` may be pushed until the history is rewritten. That needs a force-push go. Until then, no workflow can run on GitHub, and nothing in this page has been exercised on a runner.

Local validation available now:
- **Job script dry run.** The proposed `tools/ci/perf_ratio_gate_job.sh` (or a pixi task) does the following, then classifies the exit code:
  - reads `ref_commit` from the registry;
  - creates the two worktrees under a scratch parent;
  - `pixi install`s them;
  - runs `perf-ratio --ref-rev <ref> --head-rev HEAD --write build/ratio-gate/ci/session.json`, then `perf-ratio-check`.

  On this M1 Max it runs end to end, about 16 min, and exercises everything except the YAML. It must run with the local worktree layout (`compas_cgal_prs-ratio-{ref,head}` already exist here, so point a CI-style run at a separate parent to avoid moving the live arms).
- **Refusal paths in seconds.** Use an absent registry (done: exit 1), a session on a fake `cpu_model`, and a `/1` registry after the chain lands.
- **`act`.** Use it for YAML/expression syntax and step wiring only. Its Docker image is a different CPU and lacks the pixi cache, so it proves nothing about timing or calibration. It cannot run a macOS self-hosted job.
- **`actionlint`.** Static check of the workflow.

## 7. Proposed change set (PROPOSED, NOT APPLIED)

Prerequisites, in order:
1. Anchor chain landed and the registry regenerated at `/2` (plan Task 4).
2. Exit-code classes in `tools/perf_ratio_gate.py` (§4).
3. Identity rewrite of the branch, pushed with the user's go.
4. Runner class chosen (open question below).

Files:
- `tools/ci/perf_ratio_gate_job.sh`, or better a pixi task `perf-ratio-ci` invoking a small Python module `tools/perf_ratio_ci.py`. That module owns worktree provisioning, harness run, check and exit classification; the project rule is "CLI entry points are pixi tasks". It needs a test against the exit classes.
- `.github/workflows/perf-ratio-gate.yml`, a new file. `benchmarks.yml` is untouched, so its cost profile does not change.
- `docs/performance.md` "### The ratio regression gate": a CI paragraph (trigger, runner class, how to act on refused/regression, the calibration/promotion artifact flow).
- Optional `.github/workflows/cpu-model-probe.yml` (`workflow_dispatch`, matrix ×N, prints `/proc/cpuinfo` model name), only if hosted runners are chosen.

Job outline:

```yaml
# PROPOSED -- not applied
name: perf-ratio-gate
on:
  schedule: [{cron: "17 2 * * *"}]          # nightly, main
  workflow_dispatch:
    inputs:
      mode: {type: choice, options: [check, calibrate], default: check}
  pull_request:
    types: [labeled, synchronize]
    branches: [main]
concurrency:
  group: perf-ratio-gate-${{ inputs.runner || 'default' }}   # one session per machine
  cancel-in-progress: false
jobs:
  gate:
    if: >-
      github.event_name != 'pull_request' ||
      contains(github.event.pull_request.labels.*.name, 'perf-gate')
    runs-on: [self-hosted, macOS, ARM64, perf-ratio-m1max]
    timeout-minutes: ${{ inputs.mode == 'calibrate' && 180 || 120 }}
    steps:
      - uses: actions/checkout@v4
        with: {fetch-depth: 0, submodules: recursive}
      - uses: prefix-dev/setup-pixi@v0.8.1
        with: {environments: default, cache: true}
      - name: Provision arms (ref at registry ref_commit, head at HEAD)
        run: pixi run perf-ratio-ci provision          # worktrees + pixi install per arm
      - uses: actions/cache@v4                          # ref arm env + native build, keyed by commit
        with:
          path: ../compas_cgal_prs-ratio-ref/build
          key: ratio-ref-${{ runner.os }}-${{ hashFiles('.ratio-ref-key') }}
      - name: Gate session                              # non-zero here = infra
        if: inputs.mode != 'calibrate'
        run: pixi run perf-ratio-ci session --write build/ratio-gate/ci/session.json
      - name: Check (0 pass, 1 regression, 3 refused, 4 infra)
        if: inputs.mode != 'calibrate'
        run: pixi run perf-ratio-ci check build/ratio-gate/ci/session.json
      - name: Calibrate (records only; committed by hand)
        if: inputs.mode == 'calibrate'
        run: pixi run perf-ratio-calibrate --sessions 5
      - uses: actions/upload-artifact@v4
        if: always()
        with:
          name: perf-ratio-${{ github.run_id }}
          path: |
            build/ratio-gate/
            benchmarks/results/perf_ratio_calibration/
```

`.ratio-ref-key` is written by `provision`. It contains `ref_commit`, the ref `pixi.lock` digest, the `cpu_model` slug and the compiler version, so that the cache key is content-addressed. The refused class (exit 3) maps to a neutral outcome via `continue-on-error` on that step only, plus a job summary line naming the refusal. Regression and infra stay red.

## Runner class (decided)

User ruling 2026-11-01: **(a) self-hosted on the calibrated M1 Max.** The options as presented:
- **(a) self-hosted runner on the calibrated M1 Max.** No new calibration; it measures where the margins were derived. Costs: a runner daemon on a personal machine, and PR-triggered code running there. Restrict it to label/dispatch by maintainers on a private fork.
- **(b) GitHub-hosted `ubuntu-latest`.** Needs the probe job first. About 70 min of calibration per observed CPU model per promotion, and runs on an unseen model end as "refused".
