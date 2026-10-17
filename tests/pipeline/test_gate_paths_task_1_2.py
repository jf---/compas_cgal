"""The files Task 1.2 edits on the gate's import closure are gate paths (amendment M7, 1.1-closure rule)."""

from __future__ import annotations

import pytest

from tools.gate_ratchet import touches_gate

TASK_1_2_GATE_INPUTS = ("benchmarks/held_path_snapshot.py", "tools/gate_measure.py")


@pytest.mark.parametrize("path", TASK_1_2_GATE_INPUTS)
def test_gate_paths_cover_every_task_1_2_gate_input(path: str) -> None:
    assert touches_gate([path]), f"{path} is a gate input 1.2 edits but not under GATE_PATHS"
