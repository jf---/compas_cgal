# Coherence Wave 1 Implementation Plan

> **status: in execution** — opened 2026-08-28. The commit that lands a task
> updates this header with the task number. Landed: 2, 3, 4, 5, 6.

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make invariants I2 (one truth per claim), I3 (one meaning for red)
and I4 (plan-header discipline) of the coherence programme real and enforced,
entirely on `codex/sdd-coherence`, without touching anything the
auditor-convergence fleet owns.

**Architecture:** Three small governance tools in a new top-level `tools/`
package (red-manifest checker, plan-header lint, stamped measured-run wrapper),
one committed measured corpus artifact under `benchmarks/results/`, a
disposition ledger for every measurement-asserting comment, qualifiers on the
six living docs pages that state performance figures, one ablation docs page,
and one attribution test proving the two registered generators diverge below
the saturated cap. Everything lands by explicit pathspec, one commit per task.

**Tech Stack:** Python (3.9-compatible in `tools/` — stdlib only: `argparse`,
`json`, `re`, `subprocess`, `xml.etree`), pytest + junitxml, pixi tasks,
mkdocs-material.

**Spec:** `docs/superpowers/state/backlog.md` — the Wave-1 table and invariants
I1–I5. This plan implements items R7, A, B, C2, C3, M1 and the enforcement half
of I4. Explicitly OUT of scope: C1 (tangency-layer investigation — research,
not plannable as TDD), C4/C5 (parked capability lane), and every Wave-2 item
(the fleet's; governance boundary).

## Global Constraints

- Pixi exclusively; every pytest/python invocation via `pixi run` from the
  worktree `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence`.
- Branch `codex/sdd-coherence` only. NEVER touch the auditor worktree, the
  frozen refs (`jf/toolpath-redesign@73d5372`, `…t9-zero-guide@073a0f7`), or
  the programme's plan files.
- Commit per task, by explicit pathspec (`git commit -F msg -- <paths>`),
  author AND committer `Jelle Feringa <jelleferinga@gmail.com>`, no
  attribution lines, extremely concise messages.
- No `pytest.mark.skip`/`skipif`/`xfail`, ever. No `except: pass`. Named
  exceptions per failure mode. No bare tolerance literals.
- Google-style docstrings; fenced code blocks with language hints; mkdocs
  admonitions (`!!! note`), never rST.
- `ruff format` + `ruff check` before each commit; new `tools/` code passes
  `mypy --strict`.
- Historical documents are immutable: plans/specs keep their numbers; only
  LIVING docs pages get qualifier edits (Task 8 names the six).

---

### Task 1: Execution environment

The t9-zero-guide checkout that held the only spare pixi env was removed for
disk (backlog D4). Recreate the build environment inside the Wave-1 worktree.

**Files:** none created — environment only.

**Interfaces:**
- Produces: a worktree at
  `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence` where
  `pixi run pytest …` works and `_stock_2` imports. Every later task runs here.

- [x] **Step 1: worktree exists on `codex/sdd-coherence`** (created at plan
  commit time; verify)

Run: `git -C /Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-sdd-coherence status --short --branch`
Expected: `## codex/sdd-coherence`, clean.

- [x] **Step 2: install the env and build the extensions**

Run (long — CGAL compiles): `cd <worktree> && pixi install && pixi run _editable-rebuild`
Expected: exit 0; the import line prints nothing.

- [x] **Step 3: sanity gate**

Run: `pixi run pytest tests/benchmarks/test_quality.py -q -n auto`
Expected: `6 failed, N passed` — exactly the six product-gate reds, nothing
else red. If anything else is red, STOP: the env is wrong, do not proceed.

- [x] **Step 4: no commit** (environment is not a repo change).

---

### Task 2: Red-manifest checker (I3)

**Files:**
- Create: `tools/__init__.py` (empty), `tools/red_manifest.py`
- Create: `docs/red_manifest.json`
- Test: `tests/tools/__init__.py` (empty), `tests/tools/test_red_manifest.py`
- Modify: `pyproject.toml` (two pixi task keys, shown in Step 6)

**Interfaces:**
- Consumes: a junit XML file produced by `pytest --junitxml=…`.
- Produces: CLI `python -m tools.red_manifest <junit.xml> --manifest docs/red_manifest.json`,
  exit 0 iff suite reds == manifest, both directions. Task 11 wires it into
  the oracle; Wave-2/P4 wires it into CI.

- [x] **Step 1: write the failing tests**

```python
"""tests/tools/test_red_manifest.py"""

import json
import pathlib
import textwrap

import pytest

from tools.red_manifest import ManifestViolation, check


def _junit(tmp_path: pathlib.Path, cases: list) -> pathlib.Path:
    """One <testcase> per (classname, name, failed) triple."""
    body = "".join(
        f'<testcase classname="{c}" name="{n}">{"<failure/>" if f else ""}</testcase>'
        for c, n, f in cases
    )
    p = tmp_path / "junit.xml"
    p.write_text(f'<testsuites><testsuite>{body}</testsuite></testsuites>')
    return p


def _manifest(tmp_path: pathlib.Path, entries: list) -> pathlib.Path:
    p = tmp_path / "red_manifest.json"
    p.write_text(json.dumps({"expected_red": entries}))
    return p


GATE = {"match": r"tests\.benchmarks\.test_quality::test_the_generated_path_is_worth_running", "count": 1, "reason": "product gate", "closes_with": "generators"}


def test_exact_match_passes(tmp_path):
    junit = _junit(tmp_path, [("tests.benchmarks.test_quality", "test_the_generated_path_is_worth_running[a-b]", True), ("tests.x", "test_ok", False)])
    manifest = _manifest(tmp_path, [GATE])
    assert check(junit, manifest) == []


def test_a_red_outside_the_manifest_is_a_defect(tmp_path):
    junit = _junit(tmp_path, [("tests.x", "test_new_break", True)])
    manifest = _manifest(tmp_path, [])
    violations = check(junit, manifest)
    assert len(violations) == 1
    assert violations[0].kind == "unexpected-red"


def test_green_inside_the_manifest_is_a_finding(tmp_path):
    junit = _junit(tmp_path, [("tests.x", "test_ok", False)])
    manifest = _manifest(tmp_path, [GATE])
    violations = check(junit, manifest)
    assert len(violations) == 1
    assert violations[0].kind == "expected-red-went-green"


def test_count_mismatch_is_reported(tmp_path):
    junit = _junit(tmp_path, [("tests.benchmarks.test_quality", "test_the_generated_path_is_worth_running[a]", True), ("tests.benchmarks.test_quality", "test_the_generated_path_is_worth_running[b]", True)])
    manifest = _manifest(tmp_path, [GATE])
    assert [v.kind for v in check(junit, manifest)] == ["count-mismatch"]


def test_a_malformed_manifest_raises_a_named_error(tmp_path):
    from tools.red_manifest import MalformedManifestError

    junit = _junit(tmp_path, [])
    bad = tmp_path / "m.json"
    bad.write_text('{"expected_red": [{"match": "x"}]}')
    with pytest.raises(MalformedManifestError):
        check(junit, bad)
```

- [x] **Step 2: run, verify failure**

Run: `pixi run pytest tests/tools/test_red_manifest.py -q`
Expected: FAIL — `ModuleNotFoundError: tools`.

- [x] **Step 3: implement**

```python
"""tools/red_manifest.py — invariant I3: red outside the manifest is a defect,
green inside it is a finding.

The manifest (docs/red_manifest.json) enumerates every deliberate red with a
reason and the backlog item that closes it. This checker diffs a junit XML
against it, both directions, and its exit code is the invariant.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import List

REQUIRED_KEYS = ("match", "count", "reason", "closes_with")


class MalformedManifestError(Exception):
    """The manifest is missing a required key or is not valid JSON."""


@dataclass(frozen=True)
class ManifestViolation:
    """One divergence between the suite and the manifest.

    Attributes:
        kind: ``unexpected-red`` | ``expected-red-went-green`` | ``count-mismatch``.
        detail: Human-readable statement naming the tests or the entry.
    """

    kind: str
    detail: str


def _failed_ids(junit: Path) -> List[str]:
    """``classname::name`` for every testcase holding a failure or error child."""
    root = ET.parse(junit).getroot()
    out: List[str] = []
    for case in root.iter("testcase"):
        if any(child.tag in ("failure", "error") for child in case):
            out.append(f"{case.get('classname', '')}::{case.get('name', '')}")
    return out


def check(junit: Path, manifest: Path) -> List[ManifestViolation]:
    """Diff the suite's reds against the manifest, both directions.

    Args:
        junit: Path to a pytest ``--junitxml`` report.
        manifest: Path to the red-manifest JSON.

    Returns:
        Violations; empty means the invariant holds.

    Raises:
        MalformedManifestError: An entry lacks one of `REQUIRED_KEYS`.
    """
    entries = json.loads(manifest.read_text()).get("expected_red", [])
    for entry in entries:
        missing = [k for k in REQUIRED_KEYS if k not in entry]
        if missing:
            raise MalformedManifestError(f"manifest entry {entry!r} lacks {missing}")
    reds = _failed_ids(junit)
    violations: List[ManifestViolation] = []
    unclaimed = list(reds)
    for entry in entries:
        pattern = re.compile(entry["match"])
        matched = [r for r in unclaimed if pattern.search(r)]
        for r in matched:
            unclaimed.remove(r)
        if not matched and entry["count"] > 0:
            violations.append(ManifestViolation("expected-red-went-green", f"{entry['match']} (reason: {entry['reason']}) — a finding, not a pass"))
        elif matched and len(matched) != entry["count"]:
            violations.append(ManifestViolation("count-mismatch", f"{entry['match']}: expected {entry['count']} red, saw {len(matched)}"))
    violations.extend(ManifestViolation("unexpected-red", r) for r in unclaimed)
    return violations


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Diff suite reds against docs/red_manifest.json.")
    parser.add_argument("junit", type=Path)
    parser.add_argument("--manifest", type=Path, default=Path("docs/red_manifest.json"))
    args = parser.parse_args(argv)
    violations = check(args.junit, args.manifest)
    for violation in violations:
        print(f"{violation.kind}: {violation.detail}")
    if not violations:
        print("red set == manifest, both directions")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
```

Note `List[str] | None` is a 3.10+ annotation — under
`from __future__ import annotations` it is a string and 3.9-safe.

- [x] **Step 4: write the manifest itself** (`docs/red_manifest.json`)

```json
{
  "expected_red": [
    {
      "match": "tests\\.benchmarks\\.test_quality::test_the_generated_path_is_worth_running",
      "count": 6,
      "reason": "product gate: both registered generators genuinely fail 5-7 criteria per pocket",
      "closes_with": "generator work (parked lane C4/C5 feeds this)"
    },
    {
      "match": "tests\\.benchmarks\\.test_quality_invariants::test_moving_the_pocket_across_the_table_changes_no_metric",
      "count": 1,
      "reason": "deliberate: engagement is translation-variant at exact rim-on-boundary tangency (machining_metric_validity.md case 4)",
      "closes_with": "backlog C1"
    },
    {
      "match": "tests\\.adaptive\\.(test_generator|test_route_retrace_generator)::",
      "count": 4,
      "reason": "plan superseded (backlog D1); retirement scheduled through the auditor programme",
      "closes_with": "backlog S2"
    }
  ]
}
```

- [x] **Step 5: run tests, verify pass**

Run: `pixi run pytest tests/tools/test_red_manifest.py -q`
Expected: 5 passed.

- [x] **Step 6: pixi task + end-to-end against the real suite**

Add to `[tool.pixi.tasks]` (mirroring the existing `baseline` incantation):

```toml
_junit-baseline = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests -n auto -q --junitxml=build/junit-baseline.xml''', depends-on = ["_editable-rebuild"] }
red-manifest = { cmd = "python -m tools.red_manifest build/junit-baseline.xml", depends-on = ["_junit-baseline"], description = "Suite reds must equal docs/red_manifest.json, both directions" }
```

Run: `pixi run red-manifest`
Expected: exit 0, `red set == manifest, both directions`. If the live suite
disagrees with the manifest counts, the MANIFEST is corrected to observed
reality (with reasons) — never the reverse — and the correction is stated in
the commit message.

- [x] **Step 7: gates + commit**

Run: `pixi run ruff format tools tests/tools && pixi run ruff check tools tests/tools && pixi run mypy --strict tools`
Then: `git commit -m "feat(tools): red-manifest checker -- red outside the manifest is a defect" -- tools tests/tools docs/red_manifest.json pyproject.toml`
Update this plan's header: `Landed: 2`.

---

### Task 3: Plan-header lint (I4)

**Files:**
- Create: `tools/plan_headers.py`
- Test: `tests/tools/test_plan_headers.py`
- Modify: `pyproject.toml` (one task key)

**Interfaces:**
- Produces: CLI `python -m tools.plan_headers [--plans-dir docs/superpowers/plans]`,
  exit 0 iff every `*.md` there carries a `> **status:` line within its first
  12 lines.

- [x] **Step 1: failing tests**

```python
"""tests/tools/test_plan_headers.py"""

import pathlib

from tools.plan_headers import missing_headers


def _plan(tmp_path: pathlib.Path, name: str, text: str) -> None:
    (tmp_path / name).write_text(text)


def test_a_stamped_plan_passes(tmp_path):
    _plan(tmp_path, "a.md", "# A\n\n> **status: landed** — evidence.\n\nbody\n")
    assert missing_headers(tmp_path) == []


def test_an_unstamped_plan_is_named(tmp_path):
    _plan(tmp_path, "b.md", "# B\n\nbody with no header\n")
    assert [p.name for p in missing_headers(tmp_path)] == ["b.md"]


def test_a_header_buried_past_the_first_twelve_lines_does_not_count(tmp_path):
    _plan(tmp_path, "c.md", "# C\n" + "\n" * 14 + "> **status: landed**\n")
    assert [p.name for p in missing_headers(tmp_path)] == ["c.md"]
```

- [x] **Step 2: run, verify FAIL** (`ModuleNotFoundError`), then implement:

```python
"""tools/plan_headers.py — invariant I4: every plan carries a status header.

Status lives in the header and is derived from artifacts; bodies are immutable
planning history. This lint only asserts the header EXISTS — truthfulness is
the audit's job, not grep's.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

HEADER_PREFIX = "> **status:"
SEARCH_LINES = 12


def missing_headers(plans_dir: Path) -> List[Path]:
    """Plans lacking a status header in their first `SEARCH_LINES` lines."""
    bad: List[Path] = []
    for plan in sorted(plans_dir.glob("*.md")):
        head = plan.read_text().splitlines()[:SEARCH_LINES]
        if not any(line.startswith(HEADER_PREFIX) for line in head):
            bad.append(plan)
    return bad


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Every plan must carry a '> **status:' header.")
    parser.add_argument("--plans-dir", type=Path, default=Path("docs/superpowers/plans"))
    args = parser.parse_args(argv)
    bad = missing_headers(args.plans_dir)
    for plan in bad:
        print(f"missing status header: {plan}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [x] **Step 3: run tests → pass; run against the real tree**

Run: `pixi run pytest tests/tools/test_plan_headers.py -q` → 3 passed.
Run: `pixi run python -m tools.plan_headers`
Expected: exit 0 — the eight legacy plans are stamped and the five programme
plans got headers… **they did NOT** (they are the fleet's; out of bounds).
Expected reality: the five `2026-08-23-auditor-convergence-*` plans FAIL the
lint. Resolution, and the only one that respects the boundary: add
`--allow-prefix 2026-08-23-auditor-convergence` support — an explicit,
visible exemption list defaulting to that prefix, removed at Wave-2 close
(W2.2 stamps their plans):

```python
    parser.add_argument("--allow-prefix", action="append", default=["2026-08-23-auditor-convergence"], help="Plan-name prefixes exempt until their owner stamps them (Wave-2 W2.2).")
```

and in `missing_headers`, skip plans whose `plan.name.startswith(tuple(allow))`
(thread the argument through; add a fourth test:
`test_an_exempted_prefix_is_skipped`).

- [x] **Step 4: pixi task, gates, commit**

```toml
plan-headers = { cmd = "python -m tools.plan_headers", description = "Every SDD plan carries a status header (I4)" }
```

Run: `pixi run plan-headers` → exit 0.
`ruff format/check`, `mypy --strict tools`, then
`git commit -m "feat(tools): plan-header lint, programme plans exempt until W2.2" -- tools/plan_headers.py tests/tools/test_plan_headers.py pyproject.toml`.
Header: `Landed: 2, 3`.

---

### Task 4: R7 — the committed measured artifact

**Files:**
- Create: `tools/measured_run.py`
- Create: `benchmarks/results/` (the run's output, committed)
- Test: `tests/tools/test_measured_run.py`
- Modify: `pyproject.toml` (one task key)

**Interfaces:**
- Consumes: `python -m benchmarks.cli corpus --name all --out <dir>`
  (writes `benchmark_report.md` + `benchmark_report.json`; `--name all`
  aggregates the six authored corpora, external excluded).
- Produces: `benchmarks/results/<YYYY-MM-DD>-<shortsha>/` containing both
  reports plus `stamp.json` with keys
  `commit, dirty, python, platform, started, finished, argv`. Ledger rows in
  Task 6/7 cite this directory.

- [x] **Step 1: failing test**

```python
"""tests/tools/test_measured_run.py"""

import json
import pathlib

from tools.measured_run import STAMP_KEYS, latest_result


def test_the_committed_artifact_is_stamped_and_non_empty():
    result = latest_result(pathlib.Path("benchmarks/results"))
    stamp = json.loads((result / "stamp.json").read_text())
    assert set(STAMP_KEYS) <= set(stamp)
    assert stamp["dirty"] is False, "a measured artifact from a dirty tree proves nothing"
    records = json.loads((result / "benchmark_report.json").read_text())
    assert len(records) > 0
```

- [x] **Step 2: implement the wrapper**

```python
"""tools/measured_run.py — a corpus run that carries its own provenance.

Invariant I2 needs claims to trace to committed artifacts; an artifact without
its commit, environment, and command line is a number, not evidence.
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import List

STAMP_KEYS = ("commit", "dirty", "python", "platform", "started", "finished", "argv")


class NoMeasuredResultError(Exception):
    """benchmarks/results holds no run directory."""


def latest_result(results: Path) -> Path:
    """Newest run directory, by name (dates sort lexically)."""
    runs = sorted(d for d in results.iterdir() if d.is_dir())
    if not runs:
        raise NoMeasuredResultError(f"{results} holds no run directory")
    return runs[-1]


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout.strip()


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the corpus and stamp the artifact with its provenance.")
    parser.add_argument("--name", default="all")
    parser.add_argument("--results", type=Path, default=Path("benchmarks/results"))
    args = parser.parse_args(argv)
    commit = _git("rev-parse", "--short", "HEAD")
    dirty = bool(_git("status", "--porcelain"))
    out = args.results / f"{datetime.date.today().isoformat()}-{commit}"
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    command = [sys.executable, "-m", "benchmarks.cli", "corpus", "--name", args.name, "--out", str(out)]
    subprocess.run(command, check=True)
    stamp = {
        "commit": commit,
        "dirty": dirty,
        "python": sys.version,
        "platform": platform.platform(),
        "started": started,
        "finished": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "argv": command[1:],
    }
    (out / "stamp.json").write_text(json.dumps(stamp, indent=2))
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [x] **Step 3: pixi task; run it FROM A CLEAN TREE**

```toml
measured-run = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" python -m tools.measured_run''', depends-on = ["_editable-rebuild"], description = "Corpus run with provenance stamp, committed under benchmarks/results/" }
```

Commit the TOOL first (so the tree is clean and `dirty=false` is earnable):
`git commit -m "feat(tools): provenance-stamped corpus run" -- tools/measured_run.py tests/tools/test_measured_run.py pyproject.toml`
Then: `pixi run measured-run` (minutes). Then run the Step-1 test → PASS.

- [x] **Step 4: commit the artifact**

`git add benchmarks/results && git commit -m "bench: measured corpus run <dir-name> (R7)" -- benchmarks/results`
Header: `Landed: 2, 3, 4`.

---

### Task 5: Claims ledger — skeleton (A, part 1)

**Files:**
- Create: `docs/measurement_claims.md`
- Modify: `mkdocs.yml` — **check first**: `git log --oneline -3 -- mkdocs.yml`;
  if the fleet has touched it since `fa59120`, SKIP the nav edit and note the
  page as orphan (strict build accepts unlisted pages here).

**Interfaces:**
- Produces: one ledger row per extracted claim; Tasks 6–7 fill dispositions.

- [x] **Step 1: canonical extraction.** The claim population is DEFINED by
  this command (earlier counts differed because the regex was never fixed —
  fixing it is the point):

```bash
grep -rnE "# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\.[0-9]+)?[x×] (faster|slower))" \
  src/compas_cgal benchmarks --include='*.py' | grep -v superseded
```

- [x] **Step 2: write the ledger.** Header text, then a table with one row
  per hit, in file/line order:

```markdown
# Measurement-claim ledger

> **status: in audit** — opened 2026-08-28 under coherence Wave 1 (backlog A).

A comment asserting a measurement is a claim; a claim that cannot be re-earned
is deleted (invariant I2). The population below is exactly the output of the
extraction command in `docs/superpowers/plans/2026-08-28-coherence-wave1.md`
Task 5 — N claims on <extraction date>. Dispositions: **re-earned** (re-run,
number confirmed, reproducing command recorded) · **corrected** (number or
label was wrong; fixed in the named commit) · **historical** (configuration no
longer constructible; labeled as such in place) · **deleted**.

| # | location | claim (abbrev.) | disposition | command / commit |
| - | -------- | --------------- | ----------- | ---------------- |
| 1 | `<file>:<line>` | <first clause> | pending | — |
```

- [x] **Step 3: cardinality gate + commit.** Row count MUST equal the
  extraction hit count (state both in the commit message).
`git commit -m "docs: measurement-claim ledger skeleton, N claims extracted" -- docs/measurement_claims.md mkdocs.yml`
Header: `Landed: 2, 3, 4, 5`.

---

### Task 6: Claims dispositions — generators batch (A, part 2)

**Files:**
- Modify: `docs/measurement_claims.md`, plus any source file whose claim is
  corrected or deleted (`src/compas_cgal/engagement_radial_toolpath.py`,
  `engagement_toolpath.py`, `engagement_rho_toolpath.py`,
  `engagement_spiral_entry_toolpath.py`, `engagement_ordered_toolpath.py`).

- [x] **Step 1:** for every ledger row in `src/compas_cgal/engagement_*.py`:
  re-run the claim against the Task-4 artifact where it covers it, else with a
  targeted probe (write the probe under the scratchpad, cite the command in
  the row). Known prior: 2 of 2 previously-checked claims in the radial module
  were defective — expect corrections.
- [x] **Step 2:** fill each row's disposition; apply `corrected`/`deleted`
  edits to the source comments in the same change.
- [x] **Step 3:** `pixi run pytest tests -n auto -q -k "radial or rho or ordered or spiral"` → no new reds; `pixi run red-manifest` → exit 0.
- [x] **Step 4:** commit:
`git commit -m "docs: claims audit, generators batch -- X re-earned, Y corrected, Z deleted" -- docs/measurement_claims.md src/compas_cgal`
Header: `Landed: …, 6`.

---

### Task 7: Claims dispositions — benchmarks batch (A, part 3)

Same steps as Task 6 over the remaining rows (`benchmarks/quality.py`,
`benchmarks/gate.py`, `benchmarks/mathsm.py`, and any other extracted file).
Gate: ledger has ZERO `pending` rows; the ledger header's `status:` flips to
`complete — N/N dispositions`. Commit:
`git commit -m "docs: claims audit complete, N/N dispositions" -- docs/measurement_claims.md benchmarks`
Header: `Landed: …, 7`.

---

### Task 8: Parity qualifiers on living docs (C2)

**Files:**
- Modify (exactly these six; plans/specs are immutable history and keep their
  numbers): `docs/benchmarks.md`, `docs/continuous_engagement.md`,
  `docs/continuous_engagement_cost.md`, `docs/exactness.md`,
  `docs/segment_site_mat.md`, `docs/oblique_edge_cost.md`.

- [ ] **Step 1:** in each file, locate every performance figure
  (`grep -nE "parity|[0-9]+(\.[0-9]+)? ?ms|3–100"` per file) and attach, at
  first occurrence per page, this exact sentence (adapted only for grammar):

> Measured on axis-parallel pockets only; `center_domain()` is ~3,500× slower
> on oblique geometry with the mechanism not yet established — see
> *The Oblique-Edge Cliff*. Treat every figure on this page as best-case with
> respect to edge direction.

  `oblique_edge_cost.md` is the source page: verify it already carries the
  danger admonition; add nothing redundant there.
- [ ] **Step 2:** `pixi run -e docs docs` → strict build passes.
- [ ] **Step 3:** commit:
`git commit -m "docs: axis-parallel qualifier on every living performance claim (C2)" -- docs`
Header: `Landed: …, 8`.

---

### Task 9: Ablation record (B)

**Files:**
- Create: `docs/radius_ladder_ablation.md`
- Modify: `mkdocs.yml` (same contention check as Task 5).

- [ ] **Step 1:** write the page. Content requirements, verbatim where quoted:
  - Thesis first: *neither half of the radius-ladder repair is safe alone; the
    gate exists because both together UNGATED is worse than ranking alone.*
  - The four-way table (worst TEA / over-cap circles / cutting length on
    6×4 cap 40 and 20×12 cap 60): baseline 86.4/34/295 · 126.1/8/3236;
    ranking-only 54.5/41/389 · 88.6/12/3382; refinement-only **213.6**/1/357 ·
    126.1/4/3388; both-ungated 121.0/11/830 · 90.6/12/3379.
  - PROVENANCE ADMONITION (this is what makes the page I2-compliant):
    `!!! warning "Provenance"` — these four configurations were measured
    during development (2026-08-22) by toggling internals that shipped as a
    single gated design; they are **not reconstructible from shipped knobs**
    and are recorded as historical evidence. The one shipped knob,
    `RADIUS_LADDER_SUBDIVISIONS`, has a reproducible sweep committed beside it
    (`src/compas_cgal/engagement_radial_toolpath.py`, commit `29050b0`).
  - The constraint, stated as such: at a forced station every candidate is
    over-cap by definition and a sub-maximal circle does not finish the
    station, so lowering the peak always costs a circle (verified 32/32 by
    one-step lookahead) — a structural opposition, not a tuning trade-off.
- [ ] **Step 2:** `pixi run -e docs docs` → strict passes.
- [ ] **Step 3:** commit:
`git commit -m "docs: radius-ladder ablation record with provenance labels (B)" -- docs/radius_ladder_ablation.md mkdocs.yml`
Header: `Landed: …, 9`.

---

### Task 10: Attribution test — the generators must diverge below saturation (C3)

**Files:**
- Test: `tests/benchmarks/test_gate_attribution.py` (new)

**Interfaces:**
- Consumes: `benchmarks.gate.gate_pocket(name, tea_cap_deg=…)`,
  `benchmarks.gate.GATE_GENERATORS`, `GATE_CAP_DEG`.

- [ ] **Step 1: the test** (expected GREEN — it asserts a measured fact, and
  it turns the known cap-120 identity from a footnote into a pinned property):

```python
"""tests/benchmarks/test_gate_attribution.py

Three pockets by two generators is the smallest product in which a defect can
be attributed -- and attribution is void where both cells hold the same path.
At the saturated default cap the two registered generators are byte-identical
by construction; below saturation they must diverge, or the second generator
is not a second generator.
"""

import warnings

from benchmarks.gate import GATE_CAP_DEG, GATE_GENERATORS, gate_pocket

# Below the ladder's saturation region (the radius ladder is inert at 120).
ATTRIBUTION_CAP_DEG = 60.0


def _operations(generator_name: str, cap: float):
    spec = gate_pocket("rect_12x8", tea_cap_deg=cap)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = GATE_GENERATORS[generator_name](spec)
    return [(type(op.geometry).__name__, op.operation, op.path_index) for op in result.operations]


def test_the_default_cap_is_saturated_and_the_generators_emit_one_path():
    assert _operations("engagement_controlled", GATE_CAP_DEG) == _operations("radius_regulated", GATE_CAP_DEG)


def test_below_saturation_the_generators_diverge_so_the_gate_can_attribute():
    assert _operations("engagement_controlled", ATTRIBUTION_CAP_DEG) != _operations("radius_regulated", ATTRIBUTION_CAP_DEG)
```

- [ ] **Step 2:** run: `pixi run pytest tests/benchmarks/test_gate_attribution.py -q`
Expected: 2 passed. If the second test FAILS, that is a genuine finding
(the ladder never fires at 60 either) — STOP and report; do not weaken the
assertion.
- [ ] **Step 3:** `pixi run red-manifest` still exit 0 (nothing new red).
- [ ] **Step 4:** commit:
`git commit -m "test(gate): generators identical at saturated cap, divergent at 60 -- attribution pinned (C3)" -- tests/benchmarks/test_gate_attribution.py`
Header: `Landed: …, 10`.

---

### Task 11: Close the wave

**Files:**
- Modify: `docs/superpowers/state/backlog.md`, this plan's header.

- [ ] **Step 1:** run the Wave-1 oracle subset:
  `pixi run red-manifest && pixi run plan-headers && pixi run -e docs docs && pixi run lint && pixi run mypy --strict tools`
  — all green/exit 0.
- [ ] **Step 2:** backlog edits: mark R7, A, B, C2, C3, M1 closed (convention:
  the closing commit removes the line); note C1 remains open (its red stays in
  the manifest); flip this plan's header to
  `> **status: landed** — Tasks 1–11, evidence: pixi run red-manifest + plan-headers exit 0, artifact benchmarks/results/<dir>`.
- [ ] **Step 3:** commit:
`git commit -m "docs(sdd): Wave 1 closed -- I2/I3/I4 instrumented and enforced" -- docs/superpowers/state/backlog.md docs/superpowers/plans/2026-08-28-coherence-wave1.md`

---

## Self-Review

- **Spec coverage:** backlog Wave-1 table → R7 = Task 4, A = Tasks 5–7,
  B = Task 9, C2 = Task 8, C3 = Task 10, M1 = Task 2; I4 enforcement = Task 3;
  C1 stays open by design and its red is manifest-listed. No gaps.
- **Placeholder scan:** the only deliberately unbound values are N (claim
  count — DEFINED by Task 5's extraction command, which is the fix for the
  count having been regex-dependent) and X/Y/Z in commit messages (filled at
  execution from the ledger). No TBDs.
- **Type consistency:** `check(junit, manifest) -> List[ManifestViolation]`
  used identically in tests and CLI; `missing_headers(plans_dir) -> List[Path]`;
  `latest_result(results) -> Path`; `STAMP_KEYS` consumed by its test;
  `gate_pocket(name, tea_cap_deg=…)` matches the recon'd signature; pixi task
  incantations copied from the existing `baseline`/`figures` tasks verbatim.
- **Known risk, stated:** `pyproject.toml` gains task keys on a branch the
  fleet also extends — rebase conflict is possible there and only there;
  resolution is key-union.
