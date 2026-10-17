"""The gate ratchet over the committed vectors, and a synthetic regression it must refuse."""

from __future__ import annotations

import ast
import copy
import json
import subprocess
from pathlib import Path
from typing import Callable
from typing import Iterable

import pytest

from benchmarks.gate import FROZEN_PAIR
from benchmarks.gate import GATE_CAP_CASES
from benchmarks.gate import GATE_GENERATOR_NAMES
from benchmarks.gate import GATE_POCKET_NAMES
import tools.gate_measure as gate_measure
import tools.gate_ratchet as gate_ratchet
from tools.gate_admissions import Admission
from tools.gate_measure import MissingRebaselineApprovalError
from tools.gate_measure import MissingRebaselineCriteriaError
from tools.gate_measure import EmptyMeasurementError
from tools.gate_measure import ExcludedCaseDidNotRaiseError
from tools.gate_measure import ExcludedCaseRaisedOtherError
from tools.gate_measure import UnknownExcludedCaseError
from tools.gate_measure import case_key
from tools.gate_measure import expect_exclusion
from tools.gate_measure import gate_cases
from tools.gate_measure import partition_cases
from tools.gate_measure import rebaseline_for
from tools.gate_ratchet import AdmissionAboveReferenceError
from tools.gate_ratchet import AdmissionOfExistingCaseError
from tools.gate_ratchet import AdmissionRecordAlteredError
from tools.gate_ratchet import MissingReferenceCaseError
from tools.gate_ratchet import UnadmittedCaseError
from tools.gate_ratchet import GateRulesChangedError
from tools.gate_ratchet import GitQueryError
from tools.gate_ratchet import LegacyHeadVectorsError
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
    """The committed vectors as the frozen pair pinned them: its cases only, no admission (a synthetic-column base)."""
    raw = json.loads(COMMITTED.read_text())
    raw["cases"] = {key: cells for key, cells in raw["cases"].items() if key.split("/")[1] in FROZEN_PAIR}
    raw["admissions"] = []
    return raw


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
    with pytest.raises(UnadmittedCaseError):
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
    assert set(json.loads(COMMITTED.read_text())["cases"]) == expected


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


def test_the_cli_exits_one_on_a_synthetic_regression_and_zero_on_identity(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    raw = json.loads(COMMITTED.read_text())
    regressed = tmp_path / "regressed.json"
    regressed.write_text(json.dumps(_bump(raw, CASE, "cap_exceedances", +1)))
    assert main(["check-files", "--base", str(COMMITTED), "--head", str(regressed), "--touches-gate"]) == 1
    assert f"GateVectorRegressionError: gate cells worsened:\n  {CASE}/cap_exceedances:" in capsys.readouterr().err
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
    "benchmarks/pathmetrics.py": "clearance height constant for the frozen generators' call",
    "benchmarks/exceedance.py": "reached only through benchmarks.pathmetrics; the runner's demonstrated-exceedance count, no gate cell",
    "src/compas_cgal/errors.py": "named exceptions only",
    "src/compas_cgal/types.py": "array type aliases for the native bindings",
    "benchmarks/families/analytic.py": "pocket builders; a change moves the pockets, which the case product pins",
    "src/compas_cgal/__init__.py": "package marker",
    "src/compas_cgal/units.py": "unit NewTypes and frame-tagged points",
    "benchmarks/stage_acceptance.py": "stage acceptance rows; read to validate admission exemptions, no number flows into a cell",
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


# --- Ruling 0.3-C2-amend: admission record + reference --------------------------------

ADMISSION_RULING = "0.3-C2-amend"
REFERENCE_GENERATOR = "rho_regulated"


def _admission(*cases: str, reference: dict | None = None, ruling: str = ADMISSION_RULING) -> dict:
    generators = {case.split("/")[1] for case in cases}
    return {"ruling": ruling, "cases": list(cases), "reference": reference if reference is not None else dict.fromkeys(sorted(generators)), "exempt_rows": {}, "excluded": {}}


def _new_case(raw: dict, generator: str, *, delta: dict[str, int]) -> tuple[dict, str]:
    """*raw* plus a column *generator* on CASE's (cap, pocket), at the existing minimum shifted by *delta*."""
    cap_id, _generator, pocket = CASE.split("/")
    out = copy.deepcopy(raw)
    siblings = [key for key in raw["cases"] if key.split("/")[0] == cap_id and key.split("/")[2] == pocket]
    cells = copy.deepcopy(raw["cases"][siblings[0]])
    for name, spec in raw["criteria"].items():
        if spec["direction"] == "lower":
            cells[name] = min((raw["cases"][key][name] for key in siblings), key=lambda value: float(value))
    for name, shift in delta.items():
        cells[name] = str(int(cells[name]) + shift)
    key = f"{cap_id}/{generator}/{pocket}"
    out["cases"][key] = cells
    return out, key


def _losing_column(raw: dict) -> tuple[dict, str]:
    return _new_case(raw, NEW_GENERATOR, delta={"tangent_breaks": -1, "redundant_operations": +1})


def test_an_unadmitted_losing_new_case_fails_even_when_the_gate_is_untouched() -> None:
    """Hole (a): a commit that does not touch the gate could add a losing column unchecked."""
    raw = _committed()
    head, case = _losing_column(raw)
    with pytest.raises(UnadmittedCaseError, match=case):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_an_unadmitted_losing_new_case_fails_even_when_another_cell_improves() -> None:
    """Hole (b): an improvement elsewhere satisfied the stagnation check."""
    raw = _committed()
    head, case = _losing_column(raw)
    base = _bump(raw, CASE, "tangent_breaks", +1)
    with pytest.raises(UnadmittedCaseError, match=case):
        check_pair(_vectors(base), _vectors(head), touches_gate=True)


def test_an_unadmitted_new_case_fails_on_a_rules_bump_even_where_it_loses_only_on_a_rebaselined_criterion() -> None:
    """Hole (c): a bump counted as the improvement; the new column is held to the head siblings on re-baselined criteria."""
    raw = _committed()
    head, case = _new_case(raw, NEW_GENERATOR, delta={"tangent_breaks": -1, "redundant_operations": +1})
    head = _rules_bump(head, "redundant_operations")
    with pytest.raises(UnadmittedCaseError, match=case):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_an_admitted_case_with_a_null_reference_becomes_the_baseline_and_counts_as_improving() -> None:
    raw = _committed()
    head, case = _losing_column(raw)
    head["admissions"] = [_admission(case)]
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=True)
    assert verdict.improved == (f"{case} (admitted under {ADMISSION_RULING})",)


def test_an_admitted_case_at_or_below_its_reference_passes_with_the_reference_admitted_alongside() -> None:
    raw = _committed()
    head, rho = _new_case(raw, REFERENCE_GENERATOR, delta={"redundant_operations": +2})
    head, pipeline = _new_case(head, NEW_GENERATOR, delta={"redundant_operations": +1})
    head["admissions"] = [_admission(pipeline, rho, reference={NEW_GENERATOR: REFERENCE_GENERATOR, REFERENCE_GENERATOR: None})]
    verdict = check_pair(_vectors(raw), _vectors(head), touches_gate=False)
    assert verdict.improved == (f"{pipeline} (admitted under {ADMISSION_RULING})", f"{rho} (admitted under {ADMISSION_RULING})")


def test_an_admitted_case_above_its_reference_on_any_gated_cell_is_refused() -> None:
    raw = _committed()
    head, rho = _new_case(raw, REFERENCE_GENERATOR, delta={"redundant_operations": +2})
    head, pipeline = _new_case(head, NEW_GENERATOR, delta={"redundant_operations": +1, "tangent_breaks": +1})
    head["admissions"] = [_admission(pipeline, rho, reference={NEW_GENERATOR: REFERENCE_GENERATOR, REFERENCE_GENERATOR: None})]
    with pytest.raises(AdmissionAboveReferenceError, match=f"{pipeline}/tangent_breaks"):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_an_admission_whose_reference_case_is_absent_is_refused() -> None:
    raw = _committed()
    head, pipeline = _losing_column(raw)
    head["admissions"] = [_admission(pipeline, reference={NEW_GENERATOR: REFERENCE_GENERATOR})]
    with pytest.raises(MissingReferenceCaseError, match=REFERENCE_GENERATOR):
        check_pair(_vectors(raw), _vectors(head), touches_gate=False)


def test_an_existing_generator_may_be_the_reference() -> None:
    raw = _committed()
    head, pipeline = _new_case(raw, NEW_GENERATOR, delta={})
    head["admissions"] = [_admission(pipeline, reference={NEW_GENERATOR: CASE.split("/")[1]})]
    sibling = head["cases"][CASE]
    for name, spec in raw["criteria"].items():
        if spec["direction"] == "lower":
            head["cases"][pipeline][name] = sibling[name]
    assert check_pair(_vectors(raw), _vectors(head), touches_gate=False).improved == (f"{pipeline} (admitted under {ADMISSION_RULING})",)


def test_an_admission_of_a_pre_existing_case_is_refused() -> None:
    """Pre-existing cells are ratcheted exactly as before: an admission cannot re-baseline them."""
    raw = _committed()
    head = copy.deepcopy(raw)
    head["admissions"] = [_admission(CASE)]
    with pytest.raises(AdmissionOfExistingCaseError, match=CASE):
        check_pair(_vectors(raw), _vectors(head), touches_gate=True)


@pytest.mark.parametrize("how", ["removed", "edited"])
def test_an_admission_in_the_parent_must_appear_unchanged_in_the_child(how: str) -> None:
    raw = _committed()
    base, case = _losing_column(raw)
    base["admissions"] = [_admission(case)]
    head = copy.deepcopy(base)
    if how == "removed":
        head["admissions"] = []
    else:
        head["admissions"][0]["ruling"] = "1.1-closure"
    with pytest.raises(AdmissionRecordAlteredError):
        check_pair(_vectors(base), _vectors(head), touches_gate=False)


def test_an_admitted_case_is_ratcheted_from_the_next_commit() -> None:
    raw = _committed()
    base, case = _losing_column(raw)
    base["admissions"] = [_admission(case)]
    head = _bump(base, case, "redundant_operations", +1)
    with pytest.raises(GateVectorRegressionError, match=case):
        check_pair(_vectors(base), _vectors(head), touches_gate=False)
    head, later = _new_case(base, "later", delta={"redundant_operations": +5})
    head["admissions"] = [*base["admissions"], _admission(later, ruling="1.1-closure")]
    assert check_pair(_vectors(base), _vectors(head), touches_gate=False).improved == (f"{later} (admitted under 1.1-closure)",)


def test_a_schema_1_parent_ratchets_against_a_schema_2_child() -> None:
    raw = _committed()
    parent = copy.deepcopy(raw)
    parent["schema"] = "gate-vectors/1"
    del parent["admissions"]
    check_pair(_vectors(parent), _vectors(raw), touches_gate=False)


def test_the_cli_names_an_unadmitted_case(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    base = _committed()
    head, _case = _losing_column(base)
    base_path, path = tmp_path / "base.json", tmp_path / "head.json"
    base_path.write_text(json.dumps(base))
    path.write_text(json.dumps(head))
    assert main(["check-files", "--base", str(base_path), "--head", str(path)]) == 1
    assert "UnadmittedCaseError" in capsys.readouterr().err


def test_the_incident_loop_and_the_toolpath_records_are_gate_paths() -> None:
    """Ruling 1.1-closure."""
    assert touches_gate(["src/compas_cgal/incident_loop.py"])
    assert touches_gate(["src/compas_cgal/toolpath.py"])


class _ListedError(Exception):
    pass


# The listed class by its qualified name, as an admission records it.
LISTED = f"{__name__}._ListedError"


class _SubclassOfListedError(_ListedError):
    pass


def _raise(error: Exception) -> Callable[[], object]:
    def generate() -> object:
        raise error

    return generate


def test_an_excluded_case_must_raise_exactly_its_listed_exception() -> None:
    """G3: the listed class, by name, and nothing else; a listed case that does not raise is an error."""
    case = "cap-120-default/rho_regulated/L_shape"
    expect_exclusion(_raise(_ListedError("x")), case, LISTED)
    with pytest.raises(ExcludedCaseRaisedOtherError, match="_SubclassOfListedError"):
        expect_exclusion(_raise(_SubclassOfListedError("x")), case, LISTED)
    with pytest.raises(ExcludedCaseRaisedOtherError, match="RuntimeError"):
        expect_exclusion(_raise(RuntimeError("x")), case, LISTED)
    with pytest.raises(ExcludedCaseDidNotRaiseError, match=case):
        expect_exclusion(lambda: None, case, LISTED)


def test_excluded_cases_are_measured_for_their_exception_and_record_no_cells() -> None:
    cases = gate_cases()
    excluded_key = case_key(cases[0][0], cases[0][2], cases[0][3])
    measured, excluded = partition_cases(cases, {excluded_key: "pkg.errors.SomeError"})
    assert [case_key(cap_id, generator, pocket) for cap_id, _cap, generator, pocket in measured] == [case_key(c, g, p) for c, _d, g, p in cases[1:]]
    assert excluded == ((cases[0], "pkg.errors.SomeError"),)
    with pytest.raises(UnknownExcludedCaseError, match="nowhere/x/y"):
        partition_cases(cases, {"nowhere/x/y": "pkg.errors.SomeError"})


def _admitted_pipeline_and_rho(raw: dict, *, rho_reference: dict | None = None) -> tuple[dict, str, str]:
    """*raw* plus admitted rho and pipeline columns, pipeline referenced to rho and one below it on redundant_operations."""
    base, rho = _new_case(raw, REFERENCE_GENERATOR, delta={"redundant_operations": +2})
    base, pipeline = _new_case(base, NEW_GENERATOR, delta={"redundant_operations": +1})
    reference = rho_reference if rho_reference is not None else {NEW_GENERATOR: REFERENCE_GENERATOR, REFERENCE_GENERATOR: None}
    base["admissions"] = [_admission(pipeline, rho, reference=reference)]
    return base, pipeline, rho


def test_a_rules_bump_cannot_lift_an_admitted_column_above_its_reference() -> None:
    """User ruling 2026-10-15: the reference rule binds on every commit, re-baselined criteria included."""
    base, pipeline, _rho = _admitted_pipeline_and_rho(_committed())
    head = _bump(_rules_bump(base, "redundant_operations"), pipeline, "redundant_operations", +2)
    with pytest.raises(AdmissionAboveReferenceError, match=f"{pipeline}/redundant_operations"):
        check_pair(_vectors(base), _vectors(head), touches_gate=False)


def test_a_rules_bump_keeping_an_admitted_column_at_or_below_its_reference_passes() -> None:
    base, pipeline, rho = _admitted_pipeline_and_rho(_committed())
    head = _bump(_bump(_rules_bump(base, "redundant_operations"), pipeline, "redundant_operations", +3), rho, "redundant_operations", +2)
    assert check_pair(_vectors(base), _vectors(head), touches_gate=False).rebaselined == ("redundant_operations",)


def test_an_admitted_column_with_a_null_reference_is_not_held_to_any_reference() -> None:
    base, pipeline, _rho = _admitted_pipeline_and_rho(_committed(), rho_reference={NEW_GENERATOR: None, REFERENCE_GENERATOR: None})
    head = _bump(_rules_bump(base, "redundant_operations"), pipeline, "redundant_operations", +2)
    assert check_pair(_vectors(base), _vectors(head), touches_gate=False).rebaselined == ("redundant_operations",)


# --- Fix round 2 (review 2026-10-15) ---------------------------------------------------


def test_an_unadmitted_case_is_not_a_stagnation() -> None:
    """M1: one named exception per failure mode; an unadmitted case is refused whatever the commit touches."""
    assert not issubclass(UnadmittedCaseError, GateVectorStagnationError)


def test_touching_the_gate_with_only_an_earlier_admission_is_a_stagnation() -> None:
    """M1 guard: an admission the parent already carries is no improvement of this commit."""
    raw = _committed()
    base, _case = _losing_column(raw)
    base["admissions"] = [_admission(_case)]
    with pytest.raises(GateVectorStagnationError):
        check_pair(_vectors(base), _vectors(copy.deepcopy(base)), touches_gate=True)
    check_pair(_vectors(base), _vectors(copy.deepcopy(base)), touches_gate=False)


class _SerialPool:
    """A stand-in for the process pool, so monkeypatched generators reach the work."""

    def __init__(self, max_workers: int) -> None:
        self.max_workers = max_workers

    def __enter__(self) -> "_SerialPool":
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    def map(self, function: Callable[[object], object], items: Iterable[object]) -> list[object]:
        return [function(item) for item in items]


@pytest.mark.parametrize(
    ("generate", "expected"),
    [("listed", None), ("returns", ExcludedCaseDidNotRaiseError), ("other", ExcludedCaseRaisedOtherError)],
)
def test_measure_runs_and_asserts_every_excluded_case(monkeypatch: pytest.MonkeyPatch, generate: str, expected: type[Exception] | None) -> None:
    """M2: an excluded case is generated, and must raise exactly the listed class; it records no cells."""
    cap_id, cap, generator, pocket = gate_cases()[0]
    key = case_key(cap_id, generator, pocket)
    behaviour = {"listed": _raise(_ListedError("x")), "returns": lambda: None, "other": _raise(RuntimeError("x"))}[generate]
    calls: list[str] = []

    def fake(_spec: object) -> object:
        calls.append(key)
        return behaviour()

    monkeypatch.setattr(gate_measure, "ProcessPoolExecutor", _SerialPool)
    monkeypatch.setitem(gate_measure.GATE_GENERATORS, generator, fake)
    admission = Admission(ruling="0.3-C2-amend", cases=("x/y/z",), reference={"y": None}, exempt_rows={}, excluded={key: LISTED})
    if expected is None:
        vectors = gate_measure.measure_gate_vectors([(cap_id, cap, generator, pocket)], admissions=(admission,))
        assert dict(vectors.cases) == {}
        assert vectors.admissions == (admission,)
    else:
        with pytest.raises(expected, match=key):
            gate_measure.measure_gate_vectors([(cap_id, cap, generator, pocket)], admissions=(admission,))
    assert calls == [key]


def test_a_measurement_of_no_case_is_refused_by_name() -> None:
    with pytest.raises(EmptyMeasurementError):
        gate_measure.measure_gate_vectors([])


def _admitted_document() -> dict:
    head, case = _losing_column(_committed())
    head["admissions"] = [_admission(case)]
    return head


def test_fresh_measures_under_the_committed_admissions(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """(e): the freshness re-measurement knows the committed exclusions."""
    raw = _admitted_document()
    path = tmp_path / "vectors.json"
    path.write_text(json.dumps(raw))
    seen: list[tuple[Admission, ...]] = []

    def fake(cases: object = None, *, rebaseline: object = None, admissions: tuple[Admission, ...] = ()) -> GateVectors:
        seen.append(admissions)
        return GateVectors.build(raw)

    monkeypatch.setattr(gate_measure, "measure_gate_vectors", fake)
    assert main(["fresh", "--vectors", str(path)]) == 0
    assert seen == [GateVectors.build(raw).admissions]


def test_a_re_measurement_never_drops_the_admissions(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """(e): ``gate_measure.main`` writes back the admissions it read."""
    raw = _admitted_document()
    path, table = tmp_path / "vectors.json", tmp_path / "vectors.md"
    path.write_text(json.dumps(raw))

    def fake(cases: object = None, *, rebaseline: RebaselineRecord | None = None, admissions: tuple[Admission, ...] = ()) -> GateVectors:
        measured = GateVectors.build(raw)
        return GateVectors(
            rules_version=measured.rules_version,
            instrument=measured.instrument,
            criteria=measured.criteria,
            cases=measured.cases,
            rebaseline=rebaseline,
            admissions=admissions,
        )

    monkeypatch.setattr(gate_measure, "measure_gate_vectors", fake)
    assert gate_measure.main(["--write", str(path), "--table", str(table)]) == 0
    assert json.loads(path.read_text())["admissions"] == raw["admissions"]


def test_the_reference_rule_binds_when_the_reference_improves_on_an_ordinary_commit() -> None:
    """(f): no bump, no new admission; rho improving below an admitted pipeline cell is refused."""
    raw = _committed()
    base, rho = _new_case(raw, REFERENCE_GENERATOR, delta={"redundant_operations": +2})
    base, pipeline = _new_case(base, NEW_GENERATOR, delta={"redundant_operations": +1})
    base["admissions"] = [_admission(pipeline, rho, reference={NEW_GENERATOR: REFERENCE_GENERATOR, REFERENCE_GENERATOR: None})]
    head = _bump(base, rho, "redundant_operations", -2)
    with pytest.raises(AdmissionAboveReferenceError, match=f"{pipeline}/redundant_operations"):
        check_pair(_vectors(base), _vectors(head), touches_gate=False)


def test_the_head_document_must_be_schema_2() -> None:
    """(c): only a parent may be read as /1."""
    raw = _committed()
    legacy = copy.deepcopy(raw)
    legacy["schema"] = "gate-vectors/1"
    del legacy["admissions"]
    check_pair(_vectors(legacy), _vectors(raw), touches_gate=False)
    with pytest.raises(LegacyHeadVectorsError):
        check_pair(_vectors(raw), _vectors(legacy), touches_gate=False)
    with pytest.raises(LegacyHeadVectorsError):
        check_pair(None, _vectors(legacy), touches_gate=False)
