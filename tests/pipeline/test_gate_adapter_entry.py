"""The gate adapter passes the helical-entry declarations explicitly (Task 1.2, review I1).

The generator's helix pitch is the gate machine's own surrogate and its entry
ramp the gate's declared one, so generator and gate cannot disagree when
`GATE_MACHINE` changes.
"""

from __future__ import annotations

from typing import Any

import pytest

import benchmarks.gate as gate
from benchmarks.gate import GATE_ENTRY_RAMP
from benchmarks.gate import gate_pocket
from benchmarks.models import GATE_MACHINE
from compas_cgal.helix_units import EntryRampLength


def test_the_pipeline_adapter_passes_the_gate_machines_pitch_and_the_declared_ramp(monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict[str, Any] = {}

    def capture(*args: object, **kwargs: object) -> object:
        seen.update(kwargs)
        return object()

    monkeypatch.setattr(gate, "generate_narrow", capture)
    gate.pipeline(gate_pocket("rect_12x8"))
    assert seen["helix_pitch_per_radius"] == GATE_MACHINE.helix_pitch_per_radius
    assert seen["entry_ramp"] == GATE_ENTRY_RAMP
    assert isinstance(GATE_ENTRY_RAMP, EntryRampLength)
