from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import xml.etree.ElementTree as ET
from typing import Literal
from typing import Optional
from typing import Union

import pytest

from tools.red_manifest import MalformedJUnitError
from tools.red_manifest import MalformedManifestError
from tools.red_manifest import JUnitFailure
from tools.red_manifest import check
from tools.red_manifest import parse_junit_failures

CaseOutcome = Union[bool, Literal["pass", "failure", "error"]]
JUnitRoot = Literal["testsuite", "testsuites"]
TestCase = tuple[Optional[str], Optional[str], CaseOutcome]
ManifestEntry = dict[str, object]

PROJECT_ROOT = pathlib.Path(__file__).parents[2]

GATE: ManifestEntry = {
    "match": (
        r"tests\.benchmarks\.test_quality::"
        r"test_the_generated_path_is_worth_running"
    ),
    "count": 1,
    "reason": "product gate",
    "closes_with": "generators",
}

ADAPTIVE_RED_IDS = (
    "tests.adaptive.test_generator::test_task13f_full_continuation",
    "tests.adaptive.test_generator::test_real_active_family_stops_at_unresolved_exact_event",
    ("tests.adaptive.test_route_retrace_generator::test_route_retrace_derivation_rejects_unsupported_source_scope[terminal]"),
    ("tests.adaptive.test_route_retrace_generator::test_continuation_rejects_missing_retrace_commit"),
)
TRANSLATION_RED_ID = "tests.benchmarks.test_quality_invariants::test_moving_the_pocket_across_the_table_changes_no_metric"


def _repository_red_ids() -> list[str]:
    quality_ids = [
        (f"tests.benchmarks.test_quality::test_the_generated_path_is_worth_running[{cap}-{generator}-{pocket}]")
        for cap in ("cap-120-default", "cap-40-attribution")
        for generator in ("engagement_controlled", "radius_regulated")
        for pocket in ("rect_12x8", "rect_20x12", "L_shape")
    ]
    return [*quality_ids, TRANSLATION_RED_ID, *ADAPTIVE_RED_IDS]


def _junit(
    tmp_path: pathlib.Path,
    cases: list[TestCase],
    *,
    root_tag: JUnitRoot = "testsuites",
) -> pathlib.Path:
    """Write one testcase per identity/outcome triple."""
    suite = ET.Element("testsuite")
    for classname, name, outcome in cases:
        case = ET.SubElement(suite, "testcase")
        if classname is not None:
            case.set("classname", classname)
        if name is not None:
            case.set("name", name)
        if outcome is True:
            ET.SubElement(case, "failure")
        elif outcome in ("failure", "error"):
            ET.SubElement(case, outcome)

    root = suite
    if root_tag == "testsuites":
        root = ET.Element("testsuites")
        root.append(suite)
    path = tmp_path / "junit.xml"
    ET.ElementTree(root).write(path, encoding="unicode")
    return path


def _manifest(tmp_path: pathlib.Path, entries: list[ManifestEntry]) -> pathlib.Path:
    path = tmp_path / "red_manifest.json"
    path.write_text(json.dumps({"expected_red": entries}))
    return path


def _manifest_payload(tmp_path: pathlib.Path, payload: str) -> pathlib.Path:
    path = tmp_path / "red_manifest.json"
    path.write_text(payload)
    return path


def _run_cli(junit: pathlib.Path, manifest: pathlib.Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.red_manifest",
            str(junit),
            "--manifest",
            str(manifest),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def _task_command(task_name: str) -> str:
    result = subprocess.run(
        ["pixi", "task", "list", "--json"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    for environment in json.loads(result.stdout):
        for feature in environment["features"]:
            for task in feature["tasks"]:
                if task["name"] == task_name:
                    command = task["cmd"]
                    assert isinstance(command, str)
                    return command
    raise AssertionError(f"{task_name} task is missing")


def test_exact_match_passes(tmp_path: pathlib.Path) -> None:
    junit = _junit(
        tmp_path,
        [
            (
                "tests.benchmarks.test_quality",
                "test_the_generated_path_is_worth_running[a-b]",
                True,
            ),
            ("tests.x", "test_ok", False),
        ],
    )
    manifest = _manifest(tmp_path, [GATE])
    assert check(junit, manifest) == []


def test_a_red_outside_the_manifest_is_a_defect(tmp_path: pathlib.Path) -> None:
    junit = _junit(tmp_path, [("tests.x", "test_new_break", True)])
    manifest = _manifest(tmp_path, [])
    violations = check(junit, manifest)
    assert len(violations) == 1
    assert violations[0].kind == "unexpected-red"


def test_green_inside_the_manifest_is_a_finding(tmp_path: pathlib.Path) -> None:
    junit = _junit(tmp_path, [("tests.x", "test_ok", False)])
    manifest = _manifest(tmp_path, [GATE])
    violations = check(junit, manifest)
    assert len(violations) == 1
    assert violations[0].kind == "expected-red-went-green"


def test_count_mismatch_is_reported(tmp_path: pathlib.Path) -> None:
    junit = _junit(
        tmp_path,
        [
            (
                "tests.benchmarks.test_quality",
                "test_the_generated_path_is_worth_running[a]",
                True,
            ),
            (
                "tests.benchmarks.test_quality",
                "test_the_generated_path_is_worth_running[b]",
                True,
            ),
        ],
    )
    manifest = _manifest(tmp_path, [GATE])
    assert [v.kind for v in check(junit, manifest)] == ["count-mismatch"]


def test_repository_manifest_rejects_a_substituted_adaptive_failure(
    tmp_path: pathlib.Path,
) -> None:
    expected_ids = _repository_red_ids()
    exact_junit = _junit(
        tmp_path,
        [(identity.split("::", 1)[0], identity.split("::", 1)[1], "failure") for identity in expected_ids],
    )
    manifest = PROJECT_ROOT / "docs" / "red_manifest.json"
    assert check(exact_junit, manifest) == []

    substituted_junit = _junit(
        tmp_path,
        [
            (identity.split("::", 1)[0], identity.split("::", 1)[1], "failure")
            for identity in [
                *expected_ids[:-1],
                "tests.adaptive.test_route_retrace_generator::test_route_trigger_rejects_foreign_or_noncausal_transitions",
            ]
        ],
    )
    assert [violation.kind for violation in check(substituted_junit, manifest)] == [
        "expected-red-went-green",
        "unexpected-red",
    ]


def test_repository_manifest_rejects_a_substituted_translation_failure(
    tmp_path: pathlib.Path,
) -> None:
    expected_ids = _repository_red_ids()
    substituted_ids = [f"{identity}-substitute" if identity == TRANSLATION_RED_ID else identity for identity in expected_ids]
    junit = _junit(
        tmp_path,
        [(identity.split("::", 1)[0], identity.split("::", 1)[1], "failure") for identity in substituted_ids],
    )

    assert [violation.kind for violation in check(junit, PROJECT_ROOT / "docs" / "red_manifest.json")] == ["expected-red-went-green", "unexpected-red"]


def test_adaptive_manifest_accepts_only_the_four_exact_failures(
    tmp_path: pathlib.Path,
) -> None:
    junit = _junit(
        tmp_path,
        [(identity.split("::", 1)[0], identity.split("::", 1)[1], "failure") for identity in ADAPTIVE_RED_IDS],
    )
    assert check(junit, PROJECT_ROOT / "docs" / "red_manifest-adaptive.json") == []


def test_a_malformed_manifest_raises_a_named_error(tmp_path: pathlib.Path) -> None:
    junit = _junit(tmp_path, [])
    bad = tmp_path / "m.json"
    bad.write_text('{"expected_red": [{"match": "x"}]}')
    with pytest.raises(MalformedManifestError):
        check(junit, bad)


INVALID_MANIFESTS = [
    ("invalid-json", "{"),
    ("non-object-root", "[]"),
    ("missing-expected-red", "{}"),
    ("non-list-expected-red", '{"expected_red": {}}'),
    ("non-object-entry", '{"expected_red": [1]}'),
    ("missing-key", '{"expected_red": [{"match": "x"}]}'),
    (
        "extra-key",
        '{"expected_red": [{"match": "x", "count": 1, "reason": "r", "closes_with": "c", "extra": true}]}',
    ),
    (
        "non-string-match",
        '{"expected_red": [{"match": 1, "count": 1, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "empty-match",
        '{"expected_red": [{"match": " ", "count": 1, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "invalid-regex",
        '{"expected_red": [{"match": "[", "count": 1, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "string-count",
        '{"expected_red": [{"match": "x", "count": "1", "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "zero-count",
        '{"expected_red": [{"match": "x", "count": 0, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "negative-count",
        '{"expected_red": [{"match": "x", "count": -1, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "boolean-count",
        '{"expected_red": [{"match": "x", "count": true, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "empty-reason",
        '{"expected_red": [{"match": "x", "count": 1, "reason": "", "closes_with": "c"}]}',
    ),
    (
        "non-string-reason",
        '{"expected_red": [{"match": "x", "count": 1, "reason": null, "closes_with": "c"}]}',
    ),
    (
        "empty-closure",
        '{"expected_red": [{"match": "x", "count": 1, "reason": "r", "closes_with": " "}]}',
    ),
    (
        "non-string-closure",
        '{"expected_red": [{"match": "x", "count": 1, "reason": "r", "closes_with": null}]}',
    ),
    ("extra-root-key", '{"expected_red": [], "extra": true}'),
    (
        "float-count",
        '{"expected_red": [{"match": "x", "count": 1.5, "reason": "r", "closes_with": "c"}]}',
    ),
    (
        "null-count",
        '{"expected_red": [{"match": "x", "count": null, "reason": "r", "closes_with": "c"}]}',
    ),
]


@pytest.mark.parametrize(
    ("case_name", "payload"),
    INVALID_MANIFESTS,
    ids=[case_name for case_name, _ in INVALID_MANIFESTS],
)
def test_invalid_manifest_raises_named_error(
    tmp_path: pathlib.Path,
    case_name: str,
    payload: str,
) -> None:
    del case_name
    junit = _junit(tmp_path, [])
    manifest = _manifest_payload(tmp_path, payload)
    with pytest.raises(MalformedManifestError):
        check(junit, manifest)


def test_duplicate_ownership_is_explicit(tmp_path: pathlib.Path) -> None:
    junit = _junit(tmp_path, [("tests.x", "test_red", "failure")])
    entry: ManifestEntry = {
        "match": r"tests\.x::test_red",
        "count": 1,
        "reason": "one",
        "closes_with": "A",
    }
    manifest = _manifest(tmp_path, [entry, {**entry, "reason": "two"}])
    violations = check(junit, manifest)
    assert [violation.kind for violation in violations] == ["overlapping-red"]


def test_partial_overlap_is_order_independent(tmp_path: pathlib.Path) -> None:
    junit = _junit(
        tmp_path,
        [
            ("tests.x", "test_x", "failure"),
            ("tests.x", "test_y", "failure"),
        ],
    )
    narrow: ManifestEntry = {
        "match": r"tests\.x::test_x$",
        "count": 1,
        "reason": "narrow",
        "closes_with": "A",
    }
    broad: ManifestEntry = {
        "match": r"tests\.x::test_[xy]$",
        "count": 1,
        "reason": "broad",
        "closes_with": "B",
    }
    forward = check(junit, _manifest(tmp_path, [narrow, broad]))
    reversed_order = check(junit, _manifest(tmp_path, [broad, narrow]))
    assert forward == reversed_order
    assert "overlapping-red" in [violation.kind for violation in forward]


@pytest.mark.parametrize("root_tag", ["testsuite", "testsuites"])
def test_supported_junit_roots_preserve_parameterized_identity(
    tmp_path: pathlib.Path,
    root_tag: JUnitRoot,
) -> None:
    junit = _junit(
        tmp_path,
        [("tests.x", "test_red[a-b]", "failure")],
        root_tag=root_tag,
    )
    manifest = _manifest(
        tmp_path,
        [
            {
                "match": r"tests\.x::test_red",
                "count": 1,
                "reason": "known",
                "closes_with": "A",
            }
        ],
    )
    assert check(junit, manifest) == []


def test_junit_error_child_is_an_infrastructure_failure(tmp_path: pathlib.Path) -> None:
    junit = _junit(tmp_path, [("tests.x", "test_error", "error")])
    with pytest.raises(MalformedJUnitError) as caught:
        check(junit, _manifest(tmp_path, []))
    assert "tests.x::test_error" in str(caught.value)
    assert "<error>" in str(caught.value)


def test_junit_failure_parser_preserves_exact_identity_and_message(
    tmp_path: pathlib.Path,
) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text(
        """<testsuite>
        <testcase classname="tests.x" name="test_red[case-a]">
          <failure message="assert not violations">criterion one\ncriterion two</failure>
        </testcase>
        <testcase classname="tests.x" name="test_green" />
        </testsuite>"""
    )
    assert parse_junit_failures(junit) == (
        JUnitFailure(
            identity="tests.x::test_red[case-a]",
            message="assert not violations\ncriterion one\ncriterion two",
        ),
    )


def test_duplicate_junit_identity_is_malformed(tmp_path: pathlib.Path) -> None:
    junit = _junit(
        tmp_path,
        [
            ("tests.x", "test_red", "failure"),
            ("tests.x", "test_red", "failure"),
        ],
    )
    with pytest.raises(MalformedJUnitError, match="duplicate testcase identity"):
        parse_junit_failures(junit)


def test_multiple_failures_for_one_testcase_are_malformed(
    tmp_path: pathlib.Path,
) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text(
        """<testsuite>
        <testcase classname="tests.x" name="test_red">
          <failure message="first" />
          <failure message="second" />
        </testcase>
        </testsuite>"""
    )
    with pytest.raises(MalformedJUnitError, match="exactly one <failure>"):
        parse_junit_failures(junit)


@pytest.mark.parametrize(
    ("classname", "name"),
    [(None, "test_red"), ("tests.x", None), ("", "test_red"), ("tests.x", "")],
)
def test_missing_junit_identity_raises_named_error_with_context(
    tmp_path: pathlib.Path,
    classname: str | None,
    name: str | None,
) -> None:
    junit = _junit(tmp_path, [(classname, name, "failure")])
    with pytest.raises(MalformedJUnitError) as caught:
        check(junit, _manifest(tmp_path, []))
    assert "<testcase" in str(caught.value)


def test_malformed_junit_xml_raises_named_error(tmp_path: pathlib.Path) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text("<testsuite>")
    with pytest.raises(MalformedJUnitError):
        check(junit, _manifest(tmp_path, []))


def test_unknown_junit_root_raises_named_error(tmp_path: pathlib.Path) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text("<testcase classname='tests.x' name='test_red'><failure/></testcase>")
    with pytest.raises(MalformedJUnitError):
        check(junit, _manifest(tmp_path, []))


def test_cli_clean_exit(tmp_path: pathlib.Path) -> None:
    result = _run_cli(_junit(tmp_path, []), _manifest(tmp_path, []))
    assert result.returncode == 0
    assert result.stdout.strip() == "red set == manifest, both directions"
    assert result.stderr == ""


def test_cli_violation_exit(tmp_path: pathlib.Path) -> None:
    result = _run_cli(
        _junit(tmp_path, [("tests.x", "test_red", "failure")]),
        _manifest(tmp_path, []),
    )
    assert result.returncode == 1
    assert result.stdout.strip() == "unexpected-red: tests.x::test_red"
    assert result.stderr == ""


def test_cli_malformed_manifest_exit(tmp_path: pathlib.Path) -> None:
    result = _run_cli(
        _junit(tmp_path, []),
        _manifest_payload(tmp_path, "{"),
    )
    assert result.returncode == 2
    assert result.stdout == ""
    assert "malformed-manifest:" in result.stderr


def test_cli_malformed_junit_exit(tmp_path: pathlib.Path) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text("<testsuite>")
    result = _run_cli(junit, _manifest(tmp_path, []))
    assert result.returncode == 2
    assert result.stdout == ""
    assert "malformed-junit:" in result.stderr


def test_cli_missing_junit_exit_is_malformed(tmp_path: pathlib.Path) -> None:
    result = _run_cli(tmp_path / "missing.xml", _manifest(tmp_path, []))
    assert result.returncode == 2
    assert result.stdout == ""
    assert "malformed-junit:" in result.stderr


def test_cli_missing_manifest_exit_is_malformed(tmp_path: pathlib.Path) -> None:
    result = _run_cli(_junit(tmp_path, []), tmp_path / "missing.json")
    assert result.returncode == 2
    assert result.stdout == ""
    assert "malformed-manifest:" in result.stderr


def test_cli_junit_error_exit_is_infrastructure_failure(tmp_path: pathlib.Path) -> None:
    result = _run_cli(
        _junit(tmp_path, [("tests.x", "test_error", "error")]),
        _manifest(tmp_path, []),
    )
    assert result.returncode == 2
    assert result.stdout == ""
    assert "malformed-junit:" in result.stderr
    assert "tests.x::test_error" in result.stderr


@pytest.mark.parametrize("task_name", ["baseline", "_junit-baseline"])
def test_pixi_baselines_collect_whole_suite_with_group_scheduler(task_name: str) -> None:
    from tests.benchmarks import test_qualityfigures

    command = _task_command(task_name)
    assert command.count("pytest tests ") == 1
    assert "-n auto" in command
    assert "--dist=loadgroup" in command
    assert "--dist=loadscope" not in command
    for selector in ("--ignore", "--deselect", " -k ", " -m ", "--lf", "--last-failed"):
        assert selector not in command
    assert test_qualityfigures.pytestmark.name == "xdist_group"
    assert test_qualityfigures.pytestmark.args == ("qualityfigures",)


@pytest.mark.parametrize(
    ("pytest_status", "expected_task_status"),
    [(0, 0), (1, 0), (2, 2), (5, 5)],
)
def test_pixi_junit_handoff_preserves_status_contract(
    tmp_path: pathlib.Path,
    monkeypatch: pytest.MonkeyPatch,
    pytest_status: int,
    expected_task_status: int,
) -> None:
    fake_pytest = tmp_path / "pytest"
    fake_pytest.write_text(f"#!/bin/sh\nexit {pytest_status}\n")
    fake_pytest.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path), prepend=":")
    result = subprocess.run(
        ["sh", "-c", _task_command("_junit-baseline")],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == expected_task_status, result.stdout + result.stderr
