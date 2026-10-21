"""The five published Held toolpaths, generated live, as gouge-audit inputs.

The paths come from :mod:`benchmarks.held_toolpath_benchmark` itself (one
stage sequence: guide, inward offset, initial placement, engagement refinement,
closed connectors), so the audit sees exactly what is benchmarked and
published, including the seam connector back to the first circle.
"""

from __future__ import annotations

from collections.abc import Callable

from benchmarks.errors import BenchmarkError
from benchmarks.held_gouge_audit import AuditCircle
from benchmarks.held_gouge_audit import AuditConnector
from benchmarks.held_gouge_audit import HeldGougeAudit
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_toolpath_benchmark import GENERATORS
from benchmarks.held_toolpath_benchmark import GeneratedToolpath


class UnknownGougeWorkloadError(BenchmarkError):
    """The requested workload is not one of the five published Held paths."""


WORKLOADS: dict[str, Callable[[], GeneratedToolpath]] = GENERATORS


def generate(workload: str) -> GeneratedToolpath:
    """Run one published generator.

    Raises:
        UnknownGougeWorkloadError: ``workload`` is not in ``WORKLOADS``.
    """
    if workload not in WORKLOADS:
        raise UnknownGougeWorkloadError(f"unknown workload {workload!r}; known: {list(WORKLOADS)}")
    return WORKLOADS[workload]()


def audit_generated(path: GeneratedToolpath) -> HeldGougeAuditReport:
    """Exact gouge audit of a generated path against its case's analytic boundary."""
    audit = HeldGougeAudit.build(path.case.boundary.primitives, path.case.tool_radius.value)
    circles = [AuditCircle(circle.center, circle.radius.value) for circle in path.circles]
    connectors = [AuditConnector(transition.samples) for transition in path.transitions]
    return audit.audit(circles, connectors)
