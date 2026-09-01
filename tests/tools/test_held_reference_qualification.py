from __future__ import annotations

from pathlib import Path

import pytest

from benchmarks.spec import PocketSpec
from compas_cgal.adaptive.errors import UncertifiedZeroGuideEdgeError
from tools.held_reference_qualification import EXPECTED_CASES
from tools.held_reference_qualification import qualify_cases
from tools.held_reference_qualification import render_qualification_markdown
from tools.held_reference_qualification import write_qualification_report


class _GeneratedPath:
    def __init__(self, operation_count: int) -> None:
        self.operations = tuple(object() for _ in range(operation_count))


def _fake_success_generator(spec: PocketSpec) -> _GeneratedPath:
    return _GeneratedPath(len(spec.polygon.points))


def _generator_raising_known_failure(spec: PocketSpec) -> _GeneratedPath:
    raise UncertifiedZeroGuideEdgeError(f"{spec.name}: no certified zero-guide edge")


def test_qualification_reports_all_four_cases() -> None:
    outcomes = qualify_cases(_fake_success_generator)

    assert tuple(outcome.case_name for outcome in outcomes) == EXPECTED_CASES
    assert all(outcome.status == "passed" and outcome.operation_count > 0 for outcome in outcomes)


def test_qualification_preserves_named_product_failure() -> None:
    outcomes = qualify_cases(_generator_raising_known_failure)

    assert outcomes[0].status == "failed"
    assert outcomes[0].failure_type == "UncertifiedZeroGuideEdgeError"
    assert outcomes[0].failure_message == "figure5: no certified zero-guide edge"


def test_qualification_propagates_an_unapproved_exception() -> None:
    def broken_generator(spec: PocketSpec) -> _GeneratedPath:
        raise RuntimeError(f"unexpected {spec.name}")

    with pytest.raises(RuntimeError, match="unexpected figure5"):
        qualify_cases(broken_generator)


def test_qualification_records_an_empty_generated_path_as_named_failure() -> None:
    outcomes = qualify_cases(lambda _spec: _GeneratedPath(0))

    assert outcomes[0].status == "failed"
    assert outcomes[0].failure_type == "EmptyToolpathError"


def test_qualification_report_is_truthful_product_evidence(tmp_path: Path) -> None:
    outcomes = qualify_cases(_fake_success_generator)
    text = render_qualification_markdown(outcomes)

    assert "| case | primitives | projection vertices | status | operations | generation seconds | named failure |" in text
    assert "cycle-time superiority" not in text
    assert "does not compare generator performance" in text
    output = write_qualification_report(outcomes, tmp_path / "qualification.md")
    assert output.read_text(encoding="utf-8") == text
