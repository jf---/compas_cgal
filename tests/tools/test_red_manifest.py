import json
import pathlib

import pytest

from tools.red_manifest import check


def _junit(tmp_path: pathlib.Path, cases: list) -> pathlib.Path:
    """One <testcase> per (classname, name, failed) triple."""
    body = "".join(f'<testcase classname="{c}" name="{n}">{"<failure/>" if f else ""}</testcase>' for c, n, f in cases)
    path = tmp_path / "junit.xml"
    path.write_text(f"<testsuites><testsuite>{body}</testsuite></testsuites>")
    return path


def _manifest(tmp_path: pathlib.Path, entries: list) -> pathlib.Path:
    path = tmp_path / "red_manifest.json"
    path.write_text(json.dumps({"expected_red": entries}))
    return path


GATE = {
    "match": (
        r"tests\.benchmarks\.test_quality::"
        r"test_the_generated_path_is_worth_running"
    ),
    "count": 1,
    "reason": "product gate",
    "closes_with": "generators",
}


def test_exact_match_passes(tmp_path):
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


def test_a_malformed_manifest_raises_a_named_error(tmp_path):
    from tools.red_manifest import MalformedManifestError

    junit = _junit(tmp_path, [])
    bad = tmp_path / "m.json"
    bad.write_text('{"expected_red": [{"match": "x"}]}')
    with pytest.raises(MalformedManifestError):
        check(junit, bad)
