"""Drive one PocketSpec through generation and certification, timed separately.

The two phases answer different questions and are never summed: generation is
the part the published state of the art also does, certification is the part it
does not. The kernel diagnostics that explain a slow certification are read on a
THIRD, untimed pass, because reading them collapses the lazy-exact filter.
"""

from __future__ import annotations

import math
import time
from typing import Iterable

from benchmarks.depletion import replay_depletion
from benchmarks.instrument import probe_digits
from benchmarks.instrument import probe_size
from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# Retract height for every generated path, as a multiple of the tool diameter.
# Only "strictly above the cutting plane" is load-bearing -- links at this height
# are rapid travel the audit records as unmeasured -- and expressing it in tool
# diameters keeps the corpus scale-free, so the precision sweep can shrink an
# instance by three orders of magnitude without the clearance becoming absurd.
CLEARANCE_Z_TOOL_DIAMETERS = 2.0


def generate_toolpath(spec: PocketSpec) -> ToolpathResult:
    """Generate the toolpath the corpus measures for *spec*.

    The single place the corpus fixes generation parameters, so a measurement, a
    depletion replay, and a test all describe the same toolpath.

    Args:
        spec: The instance to generate for.

    Returns:
        The generated toolpath.
    """
    return trochoidal_mat_toolpath_circular(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


def run_spec(spec: PocketSpec, collect_digits: bool = True) -> MeasurementRecord:
    """Generate and certify one instance, timing each phase separately.

    The kernel diagnostics come from a SEPARATE untimed pass that replays the
    generated toolpath onto a fresh stock. The audit depletes a stock internally
    but discards it, and reading exact coordinates collapses the lazy filter, so
    neither probe may run on the audit's own path.

    Args:
        spec: The instance to measure.
        collect_digits: Whether to read exact-coordinate digit statistics on the
            diagnostic pass. Forcing every coordinate exact is unbounded work and
            memory on a large arrangement, which is what this switches off; the
            arrangement size is a pure counter read and is always collected.

    Returns:
        The measurement, or a failed record carrying the exception text.

    Raises:
        UnreplayableOperationError: The audit accepted a toolpath the depletion
            replay cannot represent. The two share a cut-plane model, so this is
            a contract breach rather than a bad instance, and it is deliberately
            not folded into the record's error column.
    """
    try:
        t0 = time.perf_counter()
        result = generate_toolpath(spec)
        generate_seconds = time.perf_counter() - t0

        t1 = time.perf_counter()
        report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, spec.tea_cap_rad, holes=list(spec.holes))
        certify_seconds = time.perf_counter() - t1
    except Exception as exc:  # recorded, never swallowed: the sweep continues
        return MeasurementRecord.failed(spec.name, spec.family, spec.params, spec.tool_diameter, spec.tea_cap_deg, f"{type(exc).__name__}: {exc}")

    stations = sum(op.stations for op in report.operations)
    # Measured but undecided: no station demonstrated an exceedance, yet the
    # margin could not be closed either. This is the number that decides whether
    # the certifier is usable on geometry it did not author.
    unresolved = sum(1 for op in report.operations if op.stations > 0 and not op.cap_certified and op.max_tea <= spec.tea_cap_rad)

    final_vertices, max_digits = _diagnostic_pass(spec, result, collect_digits)

    return MeasurementRecord(
        name=spec.name,
        family=spec.family,
        params=dict(spec.params),
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=spec.tea_cap_deg,
        generate_seconds=generate_seconds,
        certify_seconds=certify_seconds,
        operations=len(report.operations),
        cut_operations=report.engaged_ops,
        stations=stations,
        max_tea_deg=math.degrees(report.max_tea),
        cap_violations=report.cap_violations,
        unresolved=unresolved,
        arrangement_vertices_final=final_vertices,
        max_coordinate_digits=max_digits,
        error=None,
    )


def _diagnostic_pass(spec: PocketSpec, result: ToolpathResult, collect_digits: bool) -> tuple[int, int]:
    """Re-deplete a fresh stock outside the timed region and read the kernel probes.

    Both callers of this function must already have stopped their timers: it is
    the only place `probe_digits` runs, and reading the probes on the virgin stock
    instead of a depleted one would pin both fields to their virgin values.

    Args:
        spec: The instance being measured.
        result: The generated toolpath, replayed to reach the depleted state.
        collect_digits: Whether to read exact-coordinate digit statistics.

    Returns:
        ``(arrangement_vertices, max_coordinate_digits)``; digits are 0 when
        collection was disabled.
    """
    stock = replay_depletion(spec, result)
    vertices = probe_size(stock).vertices
    if not collect_digits:
        return vertices, 0
    # WARNING: `probe_digits` calls `Stock.coordinate_digits()`, which calls
    # `.exact()` on every sampled coordinate. That collapses the lazy-exact filter
    # and inflates every later operation on this stock. It is legal HERE only
    # because both timers are stopped and this stock is discarded on return --
    # never move this call above a `time.perf_counter()` region.
    return vertices, probe_digits(stock).max_digits


def run_corpus(specs: Iterable[PocketSpec], collect_digits: bool = True) -> list[MeasurementRecord]:
    """Measure every instance in a corpus.

    Args:
        specs: The instances to measure.
        collect_digits: Whether to read exact-coordinate digit statistics.

    Returns:
        One record per instance, in input order.
    """
    return [run_spec(spec, collect_digits=collect_digits) for spec in specs]
