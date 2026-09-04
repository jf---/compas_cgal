"""Permanent green oracle for the machining-quality gate transition."""

from __future__ import annotations

import pytest

from tests.benchmarks.test_quality import QUALITY_GATE_CASES
from tests.benchmarks.test_quality import _evaluate_quality_gate_case
from tests.benchmarks.test_quality import _expected_quality_gate_violations


@pytest.mark.parametrize(
    ("tea_cap_deg", "generator_name", "pocket_name"),
    QUALITY_GATE_CASES,
)
def test_quality_gate_reaches_the_intended_product_verdict(
    tea_cap_deg: float,
    generator_name: str,
    pocket_name: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Authenticate all evidence before the deliberate product-gate red."""
    _quality, violations = _evaluate_quality_gate_case(
        tea_cap_deg,
        generator_name,
        pocket_name,
        monkeypatch,
    )
    assert violations == _expected_quality_gate_violations(
        tea_cap_deg,
        generator_name,
        pocket_name,
    )
