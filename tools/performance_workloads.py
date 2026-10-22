"""The profiled performance workloads: their keys, labels, lanes and file locations.

`tools.performance_profile` writes one cProfile file per workload here, and
`tools.performance_heatmap` / `tools.performance_callpaths` read them. The
registry is the one place both sides agree on what a workload is.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PROFILE_DIR = ROOT / "build" / "performance"
PROFILE_SUFFIX = ".prof"
# Sidecar next to each profile: the commit, date and wall time it was taken at.
PROFILE_META_SUFFIX = ".json"


class MissingProfileError(Exception):
    """A workload's profile or sidecar is absent; run ``pixi run performance-profile``."""


class Lane(Enum):
    """Which generator family a workload runs, and so which cost buckets apply."""

    HELD = "held"
    STOCK = "stock"


@dataclass(frozen=True)
class PerformanceWorkload:
    """One profiled workload.

    Attributes:
        key: File stem and JSON key, e.g. ``"held_figure5"``.
        label: Row label in the heatmap.
        lane: The generator family.
        target: For `Lane.HELD`, a ``benchmarks.held_toolpath_benchmark.WORKLOADS``
            name; for `Lane.STOCK`, a ``tools.pipeline_timing.run_case_here`` gate case key.
    """

    key: str
    label: str
    lane: Lane
    target: str

    @property
    def profile_path(self) -> Path:
        """The cProfile file of this workload."""
        return PROFILE_DIR / f"{self.key}{PROFILE_SUFFIX}"

    @property
    def meta_path(self) -> Path:
        """The sidecar recording where and when the profile was taken."""
        return PROFILE_DIR / f"{self.key}{PROFILE_META_SUFFIX}"


HELD_FIGURE5 = PerformanceWorkload("held_figure5", "Held Figure 5", Lane.HELD, "figure5_refined")
PERFORMANCE_WORKLOADS: tuple[PerformanceWorkload, ...] = (
    HELD_FIGURE5,
    PerformanceWorkload("held_crossed_skis", "Held crossed skis", Lane.HELD, "figure8_crossed_skis_boundary"),
    PerformanceWorkload("pipeline_cap120_rect_20x12", "Pipeline 120° rect 20×12", Lane.STOCK, "cap-120-default/pipeline/rect_20x12"),
    PerformanceWorkload("pipeline_cap40_rect_12x8", "Pipeline 40° rect 12×8", Lane.STOCK, "cap-40-attribution/pipeline/rect_12x8"),
    PerformanceWorkload("rho_cap40_rect_12x8", "rho 40° rect 12×8", Lane.STOCK, "cap-40-attribution/rho_regulated/rect_12x8"),
)
WORKLOADS_BY_KEY: dict[str, PerformanceWorkload] = {workload.key: workload for workload in PERFORMANCE_WORKLOADS}


def read_profile_meta(workload: PerformanceWorkload) -> dict[str, Any]:
    """The sidecar of *workload*'s profile, after checking both files exist.

    Args:
        workload: The workload.

    Returns:
        The sidecar: ``workload``, ``target``, ``commit``, ``dirty``, ``date``, ``profiled_wall_s``.

    Raises:
        MissingProfileError: The profile or its sidecar is absent.
    """
    for path in (workload.profile_path, workload.meta_path):
        if not path.is_file():
            raise MissingProfileError(f"{path} is missing; run `pixi run performance-profile`.")
    meta: dict[str, Any] = json.loads(workload.meta_path.read_text())
    return meta
