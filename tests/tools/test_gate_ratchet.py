"""The gate ratchet over the committed vectors, and a synthetic regression it must refuse."""

from __future__ import annotations

import ast
import copy
import json
import subprocess
from pathlib import Path

import pytest

from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_POCKET_NAMES
import tools.gate_ratchet as gate_ratchet
from tools.gate_measure import MissingRebaselineApprovalError
from tools.gate_measure import MissingRebaselineCriteriaError
from tools.gate_measure import rebaseline_for
from tools.gate_ratchet import GateRulesChangedError
from tools.gate_ratchet import GitQueryError
from tools.gate_ratchet import UnreachableGateBaseError
from tools.gate_ratchet import check_range
from tools.gate_ratchet import resolve_base
from tools.gate_ratchet import GateVectorRegressionError
from tools.gate_ratchet import GateVectorStagnationError
from tools.gate_ratchet import SpuriousRetiredCriterionError
from tools.gate_ratchet import UnacknowledgedRetiredCriterionError
from tools.gate_ratchet import check_pair
from tools.gate_ratchet import main
from tools.gate_ratchet import stale_cells
from tools.gate_ratchet import touches_gate
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import GateVectors
from tools.gate_vectors import MalformedGateVectorsError
from tools.gate_vectors import MalformedRebaselineRecordError
from tools.gate_vectors import RebaselineRecord

ROOT = Path(__file__).resolve().parents[2]
COMMITTED = ROOT / GATE_VECTORS_PATH
CASE = "cap-120-default/engagement_controlled/rect_12x8"
NEW_GENERATOR = "pipeline"


def _committed() -> dict:
    return json.loads(COMMITTED.read_text())


def _bump(raw: dict, case: str, criterion: str, delta: int) -> dict:
    out = copy.deepcopy(raw)
    out["cases"][case][criterion] = str(int(out["cases"][case][criterion]) + delta)
    return out


def _vectors(raw: dict) -> GateVectors:
    return GateVectors.build(raw)


def test_the_committed_vectors_pass_against_themselves_when_the_gate_is_untouched() -> None:
    raw = _committed()
    check_pair(_vectors(raw), _vectors(raw), touches_gate=False)


def test_a_synthetic_regression_fails_the_ratchet() -> None:
    raw = _committed()
    with pytest.raises(GateVectorRegressionError, match="cap_exceedances"):
        check_pair(_vectors(raw), _vectors(_bump(raw, CASE, "cap_exceedances", +1)), touches_gate=True)


def test_a_regression_fails_even_when_the_gate_is_untouched() -> None:
    """A kernel or generator commit is not a gate path, and still may not worsen a cell."""
    raw = _committed()
    with pytest.raises(GateVectorRegressionError, match="tangent_breaks"):
        check_pair(_vectors(raw), _vectors(_bump(raw, CASE, "tangent_breaks", +1)), touches_gate=False)


def test_touching_the_gate_without_improving_fails() -> None:
    raw = _committed()
    with pytest.raises(GateVectorStagnationError):
        check_pair(_vectors(raw), _vectors(raw), touches_gate=True)


def test_one_improved_cell_and_no_worse_cell_passes() -> None:
    raw = _committed()
    verdict = check_pair(_vectors(_bump(raw, CASE, "tangent_breaks", +1)), _vectors(raw), touches_gate=True)
    assert verdict.improved == (f"{CASE}/tangent_breaks: {int(raw['cases'][CASE]['tangent_breaks']) + 1} -> {raw['cases'][CASE]['tangent_breaks']}",)


def test_a_trade_inside_one_case_is_a_regression() -> None:
    raw = _committed()
    head = _bump(_bump(raw, CASE, "tangent_breaks", -1), CASE, "redundant_operations", +1)
    with pytest.raises(GateVectorRegressionError, match="redundant_operations"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)


def test_reported_criteria_never_decide() -> None:
    raw = _committed()
    check_pair(_vectors(raw), _vectors(_bump(raw, CASE, "retract_count", +5)), touches_gate=False)


def test_dropping_a_case_is_a_regression() -> None:
    raw = _committed()
    head = copy.deepcopy(raw)
    del head["cases"][CASE]
    with pytest.raises(GateVectorRegressionError, match=CASE):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)


def test_a_parent_without_vectors_is_a_baseline() -> None:
    check_pair(None, _vectors(_committed()), touches_gate=True)


def _with_new_column(raw: dict, *, better: str | None, worse: str | None) -> dict:
    """Add a new generator column on CASE's (cap, pocket) at the minimum of the existing ones."""
    out = copy.deepcopy(raw)
    cap_id, _generator, pocket = CASE.split("/")
    siblings = [key for key in raw["cases"] if key.split("/")[0] == cap_id and key.split("/")[2] == pocket]
    gated = [name for name, spec in raw["criteria"].items() if spec["direction"] == "lower"]
    cells = copy.deepcopy(raw["cases"][siblings[0]])
    for name in gated:
        cells[name] = min((raw["cases"][key][name] for key in siblings), key=lambda value: float(value))
    for name, delta in ((better, -1), (worse, +1)):
        if name is not None:
            cells[name] = str(int(cells[name]) + delta)
    out["cases"][f"{cap_id}/{NEW_GENERATOR}/{pocket}"] = cells
    return out


def test_a_new_generator_column_that_beats_the_minimum_somewhere_and_loses_nowhere_improves() -> None:
    """Ruling 0.3-C2 option (i): how a new generator can pass without moving the frozen ones."""
    raw = _committed()
    verdict = check_pair(_vectors(raw), _vectors(_with_new_column(raw, better="tangent_breaks", worse=None)), touches_gate=True)
    assert verdict.improved == (f"cap-120-default/{NEW_GENERATOR}/rect_12x8 (new case beats the existing minimum)",)


def test_a_new_generator_column_that_loses_anywhere_is_a_stagnation() -> None:
    raw = _committed()
    with pytest.raises(GateVectorStagnationError):
        check_pair(_vectors(raw), _vectors(_with_new_column(raw, better="tangent_breaks", worse="redundant_operations")), touches_gate=True)


BUMP_APPROVAL = "Ruling 0.4-X: user go 2026-10-20, re-baseline under rule g"


def _rules_bump(raw: dict, *changed: str, retired: tuple[str, ...] = ()) -> dict:
    """*raw* under rules 0.4, with a bump record re-baselining *changed* and retiring *retired*."""
    head = copy.deepcopy(raw)
    head["rules_version"] = "0.4"
    head["rebaseline"] = {"from": raw["rules_version"], "approval": BUMP_APPROVAL, "criteria": list(changed), "retired": list(retired)}
    return head


def _stop_gating(raw: dict, name: str, how: str) -> dict:
    """*raw* with the gated criterion *name* removed from the document or demoted to reported."""
    out = copy.deepcopy(raw)
    if how == "removed":
        del out["criteria"][name]
        for cells in out["cases"].values():
            del cells[name]
    else:
        out["criteria"][name] = {**out["criteria"][name], "direction": "reported", "required": None}
    return out


def test_a_rules_version_bump_without_a_recorded_approval_is_refused() -> None:
    raw = _committed()
    head = copy.deepcopy(raw)
    head["rules_version"] = "0.4"
    with pytest.raises(GateRulesChangedError):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert check_pair(_vectors(raw), _vectors(_rules_bump(raw, "redundant_operations")), touches_gate=True).rebaselined == ("redundant_operations",)


def test_a_rules_bump_still_refuses_a_worsened_unchanged_criterion() -> None:
    """Ruling 0.3-bump-worsens: a bump re-baselines only the criteria its record lists."""
    raw = _committed()
    head = _bump(_rules_bump(raw, "redundant_operations"), CASE, "cap_exceedances", +1)
    with pytest.raises(GateVectorRegressionError, match="cap_exceedances"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_a_rules_bump_still_refuses_a_dropped_case() -> None:
    raw = _committed()
    head = _rules_bump(raw, "redundant_operations")
    del head["cases"][CASE]
    with pytest.raises(GateVectorRegressionError, match=CASE):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_a_rules_bump_passes_when_only_the_listed_criteria_change() -> None:
    """The listed criteria move freely, either way; the bump counts as the gate commit's improvement."""
    raw = _committed()
    head = _bump(_bump(_rules_bump(raw, "redundant_operations", "tangent_breaks"), CASE, "redundant_operations", +3), CASE, "tangent_breaks", -2)
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert verdict.rebaselined == ("redundant_operations", "tangent_breaks")
    assert verdict.improved == ()


def test_a_rules_bump_still_credits_an_improved_unchanged_criterion() -> None:
    raw = _committed()
    head = _bump(_rules_bump(raw, "redundant_operations"), CASE, "cap_exceedances", -1)
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert verdict.improved == (f"{CASE}/cap_exceedances: {raw['cases'][CASE]['cap_exceedances']} -> {int(raw['cases'][CASE]['cap_exceedances']) - 1}",)


def test_a_rules_bump_that_adds_a_gated_criterion_must_list_it() -> None:
    raw = _committed()
    head = _rules_bump(raw, "redundant_operations")
    head["criteria"]["new_gated_count"] = {"group": "elementary", "direction": "lower", "decimals": 0, "required": "0"}
    for cells in head["cases"].values():
        cells["new_gated_count"] = "0"
    with pytest.raises(GateRulesChangedError, match="new_gated_count"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    head["rebaseline"]["criteria"].append("new_gated_count")
    assert check_pair(_vectors(raw), _vectors(head), touches_gate=True).rebaselined == ("redundant_operations", "new_gated_count")


@pytest.mark.parametrize("how", ["removed", "demoted"])
def test_a_rules_bump_that_stops_gating_a_criterion_must_retire_it_by_name(how: str) -> None:
    """Ruling 0.3-bump-worsens, never silently: a criterion the parent gated cannot drop out on the approval line alone."""
    raw = _committed()
    head = _bump(_stop_gating(_rules_bump(raw, "redundant_operations"), "tangent_breaks", how), CASE, "redundant_operations", +1)
    with pytest.raises(UnacknowledgedRetiredCriterionError, match="tangent_breaks"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    head["rebaseline"]["retired"] = ["tangent_breaks"]
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert verdict.rebaselined == ("redundant_operations",)
    assert verdict.retired == ("tangent_breaks",)


def test_a_rules_bump_may_only_retire_criteria_the_parent_gated() -> None:
    raw = _committed()
    head = _rules_bump(raw, "redundant_operations", retired=("retract_count",))
    with pytest.raises(SpuriousRetiredCriterionError, match="retract_count"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)


def test_a_retire_only_bump_is_a_recorded_change() -> None:
    """Retiring a criterion is itself a definition change: the record needs no re-baselined criterion beside it."""
    raw = _committed()
    head = _stop_gating(_rules_bump(raw, retired=("tangent_breaks",)), "tangent_breaks", "removed")
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert (verdict.rebaselined, verdict.retired) == ((), ("tangent_breaks",))


def test_a_retire_only_bump_still_refuses_a_worsened_unchanged_criterion() -> None:
    raw = _committed()
    head = _bump(_stop_gating(_rules_bump(raw, retired=("tangent_breaks",)), "tangent_breaks", "removed"), CASE, "cap_exceedances", +1)
    with pytest.raises(GateVectorRegressionError, match="cap_exceedances"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


@pytest.mark.parametrize("how", ["removed", "demoted"])
def test_the_same_rules_version_cannot_stop_gating_a_criterion(how: str) -> None:
    raw = _committed()
    with pytest.raises(GateRulesChangedError, match="without a rules-version bump"):
        check_pair(_vectors(raw), _vectors(_stop_gating(raw, "tangent_breaks", how)), touches_gate=False)


@pytest.mark.parametrize(
    ("retired", "reason"),
    [
        (None, "retired"),
        ("tangent_breaks", "list"),
        (["redundant_operations"], "redundant_operations"),
        (["tangent_breaks", "tangent_breaks"], "tangent_breaks"),
    ],
    ids=["missing", "not-a-list", "still-gated", "duplicate"],
)
def test_a_malformed_retired_list_is_refused_by_name(retired: object, reason: str) -> None:
    head = _stop_gating(_rules_bump(_committed(), "redundant_operations"), "tangent_breaks", "removed")
    if retired is None:
        del head["rebaseline"]["retired"]
    else:
        head["rebaseline"]["retired"] = retired
    with pytest.raises(MalformedRebaselineRecordError, match=reason):
        GateVectors.build(head)


def test_a_rebaseline_record_must_name_a_changed_or_retired_criterion() -> None:
    head = _rules_bump(_committed())
    with pytest.raises(MalformedRebaselineRecordError, match="non-empty"):
        GateVectors.build(head)


@pytest.mark.parametrize(
    ("criteria", "reason"),
    [
        (None, "criteria"),
        ([], "non-empty"),
        (["no_such_criterion"], "no_such_criterion"),
        (["retract_count"], "retract_count"),
        (["redundant_operations", "redundant_operations"], "redundant_operations"),
        ("redundant_operations", "list"),
    ],
    ids=["missing", "empty", "unknown", "not-gated", "duplicate", "not-a-list"],
)
def test_a_malformed_rebaseline_record_is_refused_by_name(criteria: object, reason: str) -> None:
    head = _rules_bump(_committed(), "redundant_operations")
    if criteria is None:
        del head["rebaseline"]["criteria"]
    else:
        head["rebaseline"]["criteria"] = criteria
    with pytest.raises(MalformedRebaselineRecordError, match=reason):
        GateVectors.build(head)


def test_the_vector_keys_are_the_gate_product() -> None:
    """v-P5: adding a generator to the gate fails loudly until its vectors are measured."""
    expected = {f"{cap_id}/{generator}/{pocket}" for cap_id, _cap in GATE_CAP_CASES for generator in GATE_GENERATOR_NAMES for pocket in GATE_POCKET_NAMES}
    assert set(_committed()["cases"]) == expected


def test_a_malformed_document_is_refused_by_name() -> None:
    raw = _committed()
    del raw["cases"][CASE]["cap_exceedances"]
    with pytest.raises(MalformedGateVectorsError):
        GateVectors.build(raw)


def test_freshness_names_every_cell_that_moved_and_ignores_wall_clock() -> None:
    """Ruling 0.3-CI-fresh: gated cells and reported integers; never reported floats or wall clock."""
    raw = _committed()
    moved = _bump(raw, CASE, "retract_count", +1)
    moved["cases"][CASE]["gate_s"] = "999.99"
    moved["cases"][CASE]["recut_fraction"] = "0.0001"
    committed, measured = raw["cases"][CASE]["retract_count"], moved["cases"][CASE]["retract_count"]
    assert stale_cells(_vectors(raw), _vectors(moved)) == [f"{CASE}/retract_count: committed {committed} != measured {measured}"]


def test_the_cli_exits_one_on_a_synthetic_regression_and_zero_on_identity(tmp_path: Path) -> None:
    raw = _committed()
    regressed = tmp_path / "regressed.json"
    regressed.write_text(json.dumps(_bump(raw, CASE, "cap_exceedances", +1)))
    assert main(["check-files", "--base", str(COMMITTED), "--head", str(regressed), "--touches-gate"]) == 1
    assert main(["check-files", "--base", str(COMMITTED), "--head", str(COMMITTED)]) == 0


def test_touching_a_gate_path_is_detected_and_a_doc_change_is_not() -> None:
    assert touches_gate(["benchmarks/survey.py", "docs/index.md"])
    assert touches_gate(["src/compas_cgal/pipeline/entry.py"])
    assert not touches_gate(["docs/machining_quality.md", "benchmarks/results/gate_vectors.json", "tools/gate_ratchet.py"])


# Repository modules the gate imports that are deliberately NOT gate paths, each
# with its reason (Ruling 0.3-C3). Any other module in the closure must lie
# under GATE_PATHS, so a new gate dependency cannot slip past the ratchet.
NOT_GATE: dict[str, str] = {
    "benchmarks/errors.py": "named exceptions only; no number flows through it",
    "benchmarks/units.py": "unit NewTypes and validators; cannot change a measured value without failing validation",
    "benchmarks/held_path_snapshot.py": "immutable record of the authored stream; a change is caught by survey binding",
    "benchmarks/models.py": "material and machine models; calibrated outcomes are not gate cells",
    "benchmarks/pathmetrics.py": "clearance height constant for the frozen generators' call",
    "benchmarks/exceedance.py": "reached only through benchmarks.pathmetrics; the runner's demonstrated-exceedance count, no gate cell",
    "src/compas_cgal/errors.py": "named exceptions only",
    "src/compas_cgal/types.py": "array type aliases for the native bindings",
    "benchmarks/families/analytic.py": "pocket builders; a change moves the pockets, which the case product pins",
    "src/compas_cgal/__init__.py": "package marker",
    "src/compas_cgal/units.py": "unit NewTypes and frame-tagged points",
    "src/compas_cgal/toolpath.py": "ToolpathResult and OperationType records",
}

# The frozen generators (Ruling 0.1-C1) are measured, not part of the gate: a
# change to them moves vectors, and worsening is refused on every commit.
FROZEN_GENERATOR_PREFIXES = ("src/compas_cgal/engagement_", "src/compas_cgal/spiral", "src/compas_cgal/rho")


def _module_path(name: str) -> Path | None:
    parts = name.split(".")
    roots = {"benchmarks": ROOT, "compas_cgal": ROOT / "src", "tools": ROOT}
    if parts[0] not in roots:
        return None
    base = roots[parts[0]].joinpath(*parts)
    for candidate in (base.with_suffix(".py"), base / "__init__.py"):
        if candidate.exists():
            return candidate
    return None


def _repo_import_closure(*entries: str) -> set[str]:
    """Repository-local modules reachable by import from *entries* (AST walk, no execution)."""
    seen: set[Path] = set()
    pending = [path for path in (_module_path(entry) for entry in entries) if path is not None]
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        package = ".".join(path.relative_to(ROOT / "src" if path.is_relative_to(ROOT / "src") else ROOT).with_suffix("").parts)
        for node in ast.walk(ast.parse(path.read_text())):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    stem = package.rsplit(".", node.level)[0]
                    names = [f"{stem}.{node.module}" if node.module else stem]
                elif node.module:
                    names = [node.module] + [f"{node.module}.{alias.name}" for alias in node.names]
            for name in names:
                found = _module_path(name)
                if found is not None and found not in seen:
                    pending.append(found)
    return {path.relative_to(ROOT).as_posix() for path in seen}


def test_gate_paths_cover_the_gate_import_closure() -> None:
    """Import closure, not naming similarity: every gate dependency is a gate path or listed with a reason."""
    closure = _repo_import_closure("benchmarks.quality_observations", "benchmarks.gate", "tools.gate_measure")
    uncovered = sorted(
        path for path in closure if not touches_gate([path]) and path not in NOT_GATE and not path.startswith(FROZEN_GENERATOR_PREFIXES) and not path.startswith("tools/gate_")
    )
    assert uncovered == []
    assert sorted(set(NOT_GATE) - closure) == []  # no stale exemption


def test_a_rebaseline_approval_must_cite_a_ledger_ruling() -> None:
    raw = _rules_bump(_committed(), "redundant_operations")
    raw["rebaseline"]["approval"] = "the user said yes"
    with pytest.raises(MalformedRebaselineRecordError, match="ledger ruling id"):
        GateVectors.build(raw)


def test_a_re_measurement_never_drops_a_rebaseline_record() -> None:
    raw = _committed()
    raw["rebaseline"] = {"from": "0.2", "approval": "Ruling 0.3-X: user go", "criteria": ["redundant_operations"], "retired": ["degenerate_loops"]}
    kept = rebaseline_for(GateVectors.build(raw), None, (), ())
    assert kept == RebaselineRecord(from_version="0.2", approval="Ruling 0.3-X: user go", criteria=("redundant_operations",), retired=("degenerate_loops",))
    assert GateVectors.build(raw).to_json()["rebaseline"] == raw["rebaseline"]
    assert rebaseline_for(None, None, (), ()) is None


def test_a_rules_bump_re_measured_without_approval_or_criteria_is_refused_by_name() -> None:
    raw = _committed()
    raw["rules_version"] = "0.2"
    older = GateVectors.build(raw)
    with pytest.raises(MissingRebaselineApprovalError):
        rebaseline_for(older, None, ("redundant_operations",), ())
    with pytest.raises(MissingRebaselineCriteriaError):
        rebaseline_for(older, "Ruling 0.3-X: user go", (), ())
    record = rebaseline_for(older, "Ruling 0.3-X: user go", ("redundant_operations",), ())
    assert record == RebaselineRecord(from_version="0.2", approval="Ruling 0.3-X: user go", criteria=("redundant_operations",), retired=())
    retire_only = rebaseline_for(older, "Ruling 0.3-X: user go", (), ("degenerate_loops",))
    assert retire_only == RebaselineRecord(from_version="0.2", approval="Ruling 0.3-X: user go", criteria=(), retired=("degenerate_loops",))


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid", *args], cwd=repo, capture_output=True, text=True, check=True).stdout.strip()


def _repo_with_two_vector_commits(tmp_path: Path, second: dict) -> tuple[Path, str, str]:
    repo = tmp_path / "repo"
    (repo / GATE_VECTORS_PATH).parent.mkdir(parents=True)
    _git(repo, "init", "-q")
    (repo / GATE_VECTORS_PATH).write_text(json.dumps(_committed()))
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "root")
    root = _git(repo, "rev-parse", "HEAD")
    (repo / GATE_VECTORS_PATH).write_text(json.dumps(second))
    (repo / "benchmarks" / "survey.py").write_text("# a gate path\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "second")
    return repo, root, _git(repo, "rev-parse", "HEAD")


def test_a_root_commit_has_no_parent_and_the_next_commit_is_ratcheted_against_it(tmp_path: Path) -> None:
    repo, root, head = _repo_with_two_vector_commits(tmp_path, _bump(_committed(), CASE, "cap_exceedances", +1))
    assert gate_ratchet._parent(repo, root) is None
    assert gate_ratchet._vectors_at(repo, None) is None
    with pytest.raises(GateVectorRegressionError, match="cap_exceedances"):
        check_range(repo, root, head)


def test_an_unreachable_base_is_named_and_falls_back_to_the_merge_base(tmp_path: Path) -> None:
    repo, root, head = _repo_with_two_vector_commits(tmp_path, _bump(_committed(), CASE, "tangent_breaks", -1))
    forced_away = "0123456789abcdef0123456789abcdef01234567"
    with pytest.raises(UnreachableGateBaseError):
        resolve_base(repo, forced_away, head, None)
    assert resolve_base(repo, forced_away, head, root) == root
    assert resolve_base(repo, root, head, None) == root


def test_a_failed_git_query_is_a_named_error(tmp_path: Path) -> None:
    repo, _root, _head = _repo_with_two_vector_commits(tmp_path, _committed())
    with pytest.raises(GitQueryError):
        check_range(repo, "no-such-ref", "HEAD")
