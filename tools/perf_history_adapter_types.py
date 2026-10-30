"""The types every timing-record adapter shares: what an adapter reads and what it returns.

They live apart from ``tools.perf_history_adapters`` so that adapter modules
registered there (``tools.perf_history_adapters_pipeline``) can use them
without importing the registry that imports them.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from tools.perf_history_entry import Measurement
from tools.perf_history_entry import RecordProvenance

# Reads a file in the record's own directory (and revision) by name.
SiblingReader = Callable[[str], bytes]


@dataclass(frozen=True)
class AdaptedRecord:
    """A record's provenance and every stage it measured."""

    provenance: RecordProvenance
    measurements: tuple[Measurement, ...]


Adapter = Callable[[Any, SiblingReader], AdaptedRecord]
