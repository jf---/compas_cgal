"""Task 1.2 gate paths: every module 1.2 puts on the gate's import closure is under GATE_PATHS (amendment M7, 1.1-closure rule).

Also holds the G5 ruling that the closure test's exemptions (`NOT_GATE`) and
the gated set never overlap: a path both exempted and gated carries two
contradictory rulings; the exemption would read as a reason the path is not
gated while the ratchet gates it anyway.
"""

from __future__ import annotations

import pytest

from tests.tools.test_gate_ratchet import NOT_GATE
from tools.gate_ratchet import touches_gate

# Gate inputs 1.2 edits on the closure.
TASK_1_2_GATE_INPUTS = ("benchmarks/held_path_snapshot.py", "tools/gate_measure.py")
# The core helix modules: outside `pipeline/`, they decide the bore.
HELIX_GATE_INPUTS = ("src/compas_cgal/helix.py", "src/compas_cgal/helix_units.py")
# Split out of `benchmarks/survey.py` (a gate path by name) by Task 1.2 set B.
SURVEY_SPLIT_GATE_INPUTS = ("benchmarks/helix_observation.py", "benchmarks/survey_station.py")
# Ruling 1.2-models: `benchmarks/models.py` carries the max-ramp limit the
# survey reads, so it decides a gate cell and is gated.
MODELS_PATH = "benchmarks/models.py"


@pytest.mark.parametrize("path", TASK_1_2_GATE_INPUTS)
def test_gate_paths_cover_every_task_1_2_gate_input(path: str) -> None:
    assert touches_gate([path]), f"{path} is a gate input 1.2 edits but not under GATE_PATHS"


@pytest.mark.parametrize("path", HELIX_GATE_INPUTS)
def test_gate_paths_cover_the_core_helix_modules(path: str) -> None:
    assert touches_gate([path]), f"{path} decides the helix bore but is not under GATE_PATHS"


@pytest.mark.parametrize("path", SURVEY_SPLIT_GATE_INPUTS)
def test_gate_paths_cover_the_modules_split_out_of_the_survey(path: str) -> None:
    assert touches_gate([path]), f"{path} was split out of the gated survey but is not under GATE_PATHS"


def test_the_machine_model_carrying_the_ramp_limit_is_gated() -> None:
    assert touches_gate([MODELS_PATH]), f"{MODELS_PATH} carries the ramp limit the survey reads but is not under GATE_PATHS"


def test_not_gate_and_the_gated_set_are_disjoint() -> None:
    assert sorted(path for path in NOT_GATE if touches_gate([path])) == []
