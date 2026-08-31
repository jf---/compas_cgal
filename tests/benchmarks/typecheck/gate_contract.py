"""Static consumer contract for the gate cap's degree unit."""

from __future__ import annotations

from benchmarks.gate import GateCapDegrees
from benchmarks.gate import gate_pocket


gate_pocket("rect_12x8", tea_cap_deg=GateCapDegrees(40.0))
gate_pocket("rect_12x8", tea_cap_deg=40.0)  # type: ignore[arg-type]
