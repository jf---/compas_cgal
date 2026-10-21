"""The generated Monstera path as binary64 arrays, for the coverage-certificate reporting tools.

Reporting only: the arrays are roundings of the generated path
(``benchmarks.held_toolpath_benchmark.generate``, closing seam connector
included). They feed design statistics and a diagnostic raster, never a
certificate decision.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from benchmarks.held_toolpath_benchmark import generate
from tools.coverage_certificate_acceptance import MONSTERA

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class MonsteraPathArrays:
    """The Monstera path in millimetres (binary64 roundings).

    Attributes:
        workload: The generator workload name.
        generation_s: Generation wall time, sum of the stage clocks.
        tool_radius: r.
        outline: (m, 2) pocket outline ``case.projection``, the polygon the published figure draws.
        centre_domain: (k, 2) centre domain Omega, the radius-r inward offset loop.
        centres: (n, 2) full-circle cut centres.
        radii: (n,) path radii R.
        connector_samples: One (s_i, 2) polyline per connector, closing seam included.
    """

    workload: str
    generation_s: float
    tool_radius: float
    outline: FloatArray
    centre_domain: FloatArray
    centres: FloatArray
    radii: FloatArray
    connector_samples: tuple[FloatArray, ...]


def load_monstera_path() -> MonsteraPathArrays:
    """Generate the Monstera path and round it to arrays."""
    path = generate(MONSTERA)
    return MonsteraPathArrays(
        workload=MONSTERA,
        generation_s=path.generation_seconds,
        tool_radius=float(path.case.tool_radius.value),
        outline=np.array([(float(p.x), float(p.y)) for p in path.case.projection.points]),
        centre_domain=np.array([(float(p.x), float(p.y)) for p in path.boundary]),
        centres=np.array([(float(c.center.x), float(c.center.y)) for c in path.circles]),
        radii=np.array([float(c.radius.value) for c in path.circles]),
        connector_samples=tuple(np.array([(float(p.x), float(p.y)) for p in t.samples]) for t in path.transitions),
    )
