"""Task 1.2 gate paths: the core helix modules (outside `pipeline/`, they decide the bore) and the modules split out of the survey."""

from __future__ import annotations

import pytest

from tools.gate_ratchet import touches_gate

HELIX_GATE_INPUTS = ("src/compas_cgal/helix.py", "src/compas_cgal/helix_units.py")
# Split out of `benchmarks/survey.py` (a gate path by name) by Task 1.2 set B.
SURVEY_SPLIT_GATE_INPUTS = ("benchmarks/helix_observation.py", "benchmarks/survey_station.py")


@pytest.mark.parametrize("path", HELIX_GATE_INPUTS)
def test_gate_paths_cover_the_core_helix_modules(path: str) -> None:
    assert touches_gate([path]), f"{path} decides the helix bore but is not under GATE_PATHS"


@pytest.mark.parametrize("path", SURVEY_SPLIT_GATE_INPUTS)
def test_gate_paths_cover_the_modules_split_out_of_the_survey(path: str) -> None:
    assert touches_gate([path]), f"{path} was split out of the gated survey but is not under GATE_PATHS"
