"""Task 1.2 G5: the closure test's exemptions and the gated set never overlap.

A path both exempted (`NOT_GATE`) and gated (`GATE_PATHS`) carries two
contradictory rulings; the exemption would read as a reason the path is not
gated while the ratchet gates it anyway.
"""

from __future__ import annotations

from tests.tools.test_gate_ratchet import NOT_GATE
from tools.gate_ratchet import touches_gate

# Ruling 1.2-models: `benchmarks/models.py` carries the max-ramp limit the
# survey reads, so it decides a gate cell and is gated.
MODELS_PATH = "benchmarks/models.py"


def test_not_gate_and_the_gated_set_are_disjoint() -> None:
    assert sorted(path for path in NOT_GATE if touches_gate([path])) == []


def test_the_machine_model_carrying_the_ramp_limit_is_gated() -> None:
    assert touches_gate([MODELS_PATH]), f"{MODELS_PATH} carries the ramp limit the survey reads but is not under GATE_PATHS"
