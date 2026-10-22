"""Which declared-clearance oracle decides: the native C++ oracle or the Python reference.

Shared by the Ω* gouge audit (:mod:`benchmarks.held_omega_gouge_audit`) and
Ω* generation (:mod:`benchmarks.held_omega_toolpath`). Both oracles decide the
same verdicts; the Python one is the differential reference.
"""

from __future__ import annotations

import enum


class ClearanceBackend(enum.Enum):
    """Where the exact clearance decisions run: the native oracle (default) or the Python reference."""

    PYTHON = "python"
    NATIVE = "native"
