"""Contract of tools/held_splice_ab_results.py: the rendering (the content identity: tests/tools/test_held_source_closure.py)."""

from __future__ import annotations

from typing import Any

from tools.held_splice_ab_results import render


def _record() -> dict[str, Any]:
    rows = [
        {"candidate": 7, "slot": 3, "later_sources": 200, "rerefined": 1, "added": 1, "a_seconds": 0.5, "b_seconds": 0.004},
        {"candidate": 9, "slot": 4, "later_sources": 199, "rerefined": 2, "added": 0, "a_seconds": 0.3, "b_seconds": 0.002},
    ]
    return {
        "build": {"recorded_at": "t", "commit": "0123456789ab", "dirty": False},
        "source_sha256": "f" * 64,
        "machine": {"cpu_brand": "cpu"},
        "workload": "w",
        "sources_count": 10,
        "selected": 5,
        "emitted": 6,
        "emitted_after_nothing": 6,
        "witnesses": 2,
        "full_refinement_seconds": 1.0,
        "candidates": rows,
        "repair_round_seconds": 0.1,
        "repair_round_inserted": [7],
        "repair_round_emitted_after": 7,
        "load_1min_before_after": [1.0, 2.0],
    }


def test_render_reports_ranges_and_the_ratio() -> None:
    text = render(_record())
    assert "do not edit" in text.splitlines()[0]
    assert "`0123456789`" in text and "dirty" not in text
    assert "| later sources re-refined | 199–200 | 1–2 |" in text
    assert "| sum over 2 candidates | 0.80 s | 0.006 s (133×) |" in text
    assert "inserted sources [7]" in text
