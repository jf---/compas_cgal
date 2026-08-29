from __future__ import annotations

import importlib
import json
import pathlib
import signal
from types import SimpleNamespace
from typing import Any

import pytest


SOURCE = "a" * 40
CAPS = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0]
SPACINGS = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_figure6")


def _raw() -> dict[str, object]:
    return {
        "pocket": {"name": "rect_20x12", "family": "analytic", "tool_diameter": 2.0, "params": {"width": 20.0, "height": 12.0}, "holes": []},
        "engagement_measured_at_cap_deg": 180.0,
        "points": [
            {
                "cap_deg": cap,
                "controlled_length": 100.0 + cap,
                "controlled_cut_motions": 10,
                "controlled_entry_cuts": 1,
                "controlled_max_tea_deg": 180.0,
                "controlled_max_tea_after_entry_deg": 120.0,
                "controlled_exceedances_after_entry": 0,
                "controlled_meets_cap": cap >= 120.0,
                "mathsm_spacing_tool_diameters": 0.1,
                "mathsm_length": 200.0,
                "mathsm_max_tea_after_entry_deg": 90.0,
                "length_ratio": 0.5,
            }
            for cap in CAPS
        ],
        "spacing_trials": [
            {
                "spacing_tool_diameters": spacing,
                "length": 200.0,
                "cut_motions": 10,
                "entry_cuts": 1,
                "max_tea_deg": 180.0,
                "max_tea_after_entry_deg": 131.14 if spacing == 0.025 else (100.0 if spacing == 0.1 else 120.0),
            }
            for spacing in SPACINGS
        ],
    }


def _markdown() -> bytes:
    comparison = "\n".join("| 20 | 1 | 2 | yes | 0 | 0.1 | 2 | 3 | 0.5 |" for _ in CAPS)
    trials = "\n".join("| 0.1 | 1 | 2 | 3 | 4 |" for _ in SPACINGS)
    return (
        "# Figure 6 reproduction — path length against engagement cap (rect_20x12)\n\n"
        "## Per-cap comparison\n\n"
        "| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | "
        "controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | "
        "constant-spacing max TEA after entry (deg) | length ratio |\n"
        "| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |\n"
        f"{comparison}\n\n## Constant-spacing trials\n\n"
        "| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |\n"
        "| ---: | ---: | ---: | ---: | ---: |\n"
        f"{trials}\n"
    ).encode("utf-8")


def test_run_figure6_case_passes_and_returns_one_immutable_argv(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path) -> None:
    module = _module()
    observed: list[tuple[str, ...]] = []

    def run(argv: tuple[str, ...], *, check: bool) -> SimpleNamespace:
        assert check is False
        observed.append(argv)
        (tmp_path / "figure6.json").write_text(json.dumps(_raw()), encoding="utf-8")
        (tmp_path / "figure6.md").write_bytes(_markdown())
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", run)
    payload, executed = module.run_figure6_case(tmp_path, source_commit=SOURCE)
    assert observed == [executed]
    assert executed[0] == module.sys.executable
    assert executed[1:5] == ("-m", "benchmarks.cli", "figure6", "--out")
    assert executed[5] == str(tmp_path)
    assert payload["source_commit"] == SOURCE
    assert payload["semantic_command"] == list(module.FIGURE6_SEMANTIC_COMMAND)
    assert payload["claims"][2]["selected_values"] == {
        "fine_spacing": 0.025,
        "fine_max_tea_after_entry": 131.14,
        "comparison_spacing": 0.1,
        "comparison_max_tea_after_entry": 100.0,
        "angle_unit": "degree",
        "spacing_unit": "tool-diameter",
    }


@pytest.mark.parametrize("failure", [OSError("spawn denied"), -signal.SIGTERM, 7])
def test_child_failure_prevents_raw_parsing(monkeypatch: pytest.MonkeyPatch, tmp_path: pathlib.Path, failure: object) -> None:
    module = _module()
    (tmp_path / "figure6.json").write_text(json.dumps(_raw()), encoding="utf-8")
    (tmp_path / "figure6.md").write_bytes(_markdown())
    monkeypatch.setattr(module, "_read_raw_outputs", lambda stage: pytest.fail(f"parsed failed child output: {stage}"))

    def run(argv: tuple[str, ...], *, check: bool) -> SimpleNamespace:
        del argv, check
        if isinstance(failure, OSError):
            raise failure
        return SimpleNamespace(returncode=failure)

    monkeypatch.setattr(module.subprocess, "run", run)
    with pytest.raises(module.MeasurementClaimChildError) as caught:
        module.run_figure6_case(tmp_path, source_commit=SOURCE)
    message = str(caught.value)
    assert "benchmarks.cli" in message
    assert ("spawn denied" in message) if isinstance(failure, OSError) else (str(failure) in message)
    if failure == -signal.SIGTERM:
        assert str(signal.SIGTERM) in message
