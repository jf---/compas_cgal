"""The native declared-clearance oracle is the default: in the Ω* gouge audit and in Ω* generation.

The Python oracle stays the differential reference: selecting
``ClearanceBackend.PYTHON`` still runs it, and every native result here is
compared with it. Generation decides the containment of all initial circles in
one oracle batch.
"""

from __future__ import annotations

import inspect
import math
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from types import SimpleNamespace

import pytest

from benchmarks import held_omega_gouge_audit
from benchmarks import held_omega_gouge_audit_native
from benchmarks import held_omega_toolpath
from benchmarks.held_declared_boundary import ExactPoint2
from benchmarks.held_declared_clearance import CircleElement
from benchmarks.held_declared_clearance import ClearanceVerdict
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_declared_clearance_native import NativeDeclaredClearance
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_gouge_audit import Clearance
from benchmarks.held_gouge_audit import HeldGougeAuditReport
from benchmarks.held_omega_circles import OmegaCircle
from benchmarks.held_omega_circles import contain
from benchmarks.held_omega_circles import contain_all
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_connectors import OmegaConnectors
from benchmarks.held_clearance_backend import ClearanceBackend
from benchmarks.held_omega_gouge_audit import audit_omega_path
from benchmarks.held_omega_gouge_audit import audit_polygon_path
from benchmarks.held_omega_refinement import refine_omega_engagement
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from benchmarks.held_omega_sites import single_loop_pieces
from benchmarks.held_omega_toolpath import OmegaGeneratedToolpath
from benchmarks.held_omega_toolpath import generate_omega
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_omega_path import GAP
from tests.benchmarks.test_held_omega_path import _arc_contacts
from tests.benchmarks.test_held_omega_path import _omega_circles
from tests.benchmarks.test_held_omega_path import _polygon_path
from tests.benchmarks.test_held_omega_path import _provenance

R = Fraction(1)
TOOL = ToolRadius.build(float(R))
# The Held corpus cap (degrees): the refinement subdivides every transition above it.
CAP_DEG = 80.0


@dataclass
class _CountingOracle:
    """Forwards to an oracle and records the size of every batch it is asked."""

    oracle: NativeDeclaredClearance
    batches: list[int] = field(default_factory=list)

    def circles(self, elements: Sequence[CircleElement]) -> list[ClearanceVerdict]:
        self.batches.append(len(elements))
        return self.oracle.circles(elements)


def _concave() -> tuple[CentreDomain, DeclaredClearance, NativeDeclaredClearance]:
    boundary = concave_arc_pocket()
    return CentreDomain.build(boundary, R), DeclaredClearance.build(boundary, R), NativeDeclaredClearance.build(boundary, R)


def _decisions(report: HeldGougeAuditReport) -> list[tuple[object, ...]]:
    return [(item.element, item.index, item.clearance, item.primitive, item.primitive_kind, item.segment) for item in report.elements]


def test_the_omega_audit_defaults_to_the_native_oracle(monkeypatch: pytest.MonkeyPatch) -> None:
    domain, python, _ = _concave()
    circles = _omega_circles(domain, python, _arc_contacts())
    connectors = OmegaConnectors.build(single_loop_pieces(domain), circles).connectors
    called: list[str] = []
    native = held_omega_gouge_audit_native.audit_omega_path_native

    def recording(*args: object) -> HeldGougeAuditReport:
        called.append("omega")
        return native(*args)  # type: ignore[arg-type]

    monkeypatch.setattr(held_omega_gouge_audit, "audit_omega_path_native", recording)
    report = audit_omega_path(python, circles, connectors)
    assert called == ["omega"]
    assert _decisions(report) == _decisions(audit_omega_path(python, circles, connectors, ClearanceBackend.PYTHON))
    assert called == ["omega"]


def test_the_polygon_audit_defaults_to_the_native_oracle(monkeypatch: pytest.MonkeyPatch) -> None:
    boundary = concave_arc_pocket()
    python = DeclaredClearance.build(boundary, R)
    polygon_circles, polylines = _polygon_path(boundary, _arc_contacts())
    called: list[str] = []
    native = held_omega_gouge_audit_native.audit_polygon_path_native

    def recording(*args: object) -> HeldGougeAuditReport:
        called.append("polygon")
        return native(*args)  # type: ignore[arg-type]

    monkeypatch.setattr(held_omega_gouge_audit, "audit_polygon_path_native", recording)
    report = audit_polygon_path(python, polygon_circles, polylines)
    assert called == ["polygon"]
    assert report.gouging_circles
    assert _decisions(report) == _decisions(audit_polygon_path(python, polygon_circles, polylines, ClearanceBackend.PYTHON))


def test_the_python_oracle_answers_a_batch_like_its_singletons() -> None:
    _, python, native = _concave()
    centres = ((10, 10, 3), (2, 2, 1), (10, 6, Fraction(1, 2)), (18, 18, 1))

    elements = [CircleElement(ExactPoint2.build(Fraction(x), Fraction(y)), Fraction(r)) for x, y, r in centres]
    assert python.circles(elements) == [python.circle(element) for element in elements]
    assert [(verdict.clearance, verdict.primitive) for verdict in native.circles(elements)] == [(verdict.clearance, verdict.primitive) for verdict in python.circles(elements)]


def test_initial_containment_is_one_native_batch_and_matches_python() -> None:
    domain, python, native = _concave()
    loop = OmegaLoop.build(domain)
    decisions = domain.decide_contacts(_arc_contacts(), GAP).contacts
    circles = []
    for index, decision in enumerate(decisions):
        assert decision is not None
        circles.append(place_local(OmegaSite.decided(loop, decision), 0.5, _provenance(index)))
    counting = _CountingOracle(native)
    contained = contain_all(tuple(circles), counting, domain)
    assert counting.batches == [len(circles)]
    assert contained == tuple(contain(circle, python, domain) for circle in circles)


def test_a_shrunk_circle_is_the_same_on_both_oracles() -> None:
    boundary = l_shape_pocket()
    domain = CentreDomain.build(boundary, R)
    python, native = DeclaredClearance.build(boundary, R), NativeDeclaredClearance.build(boundary, R)
    decision = domain.decide_contacts([(18.0, 1.0)], GAP).contacts[0]
    assert decision is not None
    local = place_local(OmegaSite.decided(OmegaLoop.build(domain), decision), 2.0, _provenance(0))
    assert python.circle(local.element()).clearance is Clearance.GOUGE
    shrunk = contain(local, native, domain)
    assert shrunk.limited
    assert shrunk == contain(local, python, domain)
    assert contain_all((local,), native, domain) == (shrunk,)


def test_refinement_on_the_native_oracle_is_identical_to_python() -> None:
    domain, python, native = _concave()
    circles: tuple[OmegaCircle, ...] = _omega_circles(domain, python, _arc_contacts())
    loop = OmegaLoop.build(domain)
    cap = EngagementCap.build(math.radians(CAP_DEG))
    expected = refine_omega_engagement(circles, loop, python, TOOL, cap)
    found = refine_omega_engagement(circles, loop, native, TOOL, cap)
    assert len(expected.circles) > len(circles)
    # OmegaSite equality is structural (the cached point is excluded), so the
    # two refinements are ==; repr also pins every float and Fraction exactly
    # (the path digests' serialisation).
    assert found == expected
    assert repr(found) == repr(expected)


def test_omega_generation_defaults_to_the_native_oracle() -> None:
    assert inspect.signature(generate_omega).parameters["backend"].default is ClearanceBackend.NATIVE
    assert "generation_backend" in {item.name for item in OmegaGeneratedToolpath.__dataclass_fields__.values()}


# One workload per spacing regime: LANE refines with refine_omega_engagement,
# BOUNDARY with refine_omega_boundary_ordered. The initial placement (seconds
# to minutes per workload) is stubbed; the guide and Ω* stages run for real.
_SPACING_WORKLOADS = ("figure8_upper_lane", "figure8_upper_boundary")


@dataclass
class _GenerationProbe:
    """Records the native oracles generation builds and the oracle each refinement stage receives."""

    built: list[NativeDeclaredClearance] = field(default_factory=list)
    refined_with: list[object] = field(default_factory=list)


def _probe_generation(monkeypatch: pytest.MonkeyPatch) -> _GenerationProbe:
    probe = _GenerationProbe()
    build = NativeDeclaredClearance.build
    refinement = SimpleNamespace(circles=())

    def recording_build(*args: object) -> NativeDeclaredClearance:
        oracle = build(*args)  # type: ignore[arg-type]
        probe.built.append(oracle)
        return oracle

    def initial(*args: object, **kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(circles=(), circle_sources=(), reached_run_ids=frozenset())

    def engagement(circles: object, loop: object, oracle: object, *args: object) -> SimpleNamespace:
        probe.refined_with.append(oracle)
        return refinement

    def boundary_ordered(circles: object, sources: object, reached: object, loop: object, oracle: object, *args: object) -> tuple[None, SimpleNamespace]:
        probe.refined_with.append(oracle)
        return None, refinement

    monkeypatch.setattr(held_omega_toolpath.NativeDeclaredClearance, "build", recording_build)
    monkeypatch.setattr(held_omega_toolpath, "omega_reference_path", initial)
    monkeypatch.setattr(held_omega_toolpath, "refine_omega_engagement", engagement)
    monkeypatch.setattr(held_omega_toolpath, "refine_omega_boundary_ordered", boundary_ordered)
    monkeypatch.setattr(held_omega_toolpath.OmegaConnectors, "build", lambda *args: SimpleNamespace(connectors=()))
    monkeypatch.setattr(held_omega_toolpath, "repair_omega_connectors", lambda *args: SimpleNamespace(circles=(), connectors=SimpleNamespace(connectors=())))
    return probe


@pytest.mark.parametrize("workload", _SPACING_WORKLOADS)
def test_native_generation_refines_on_the_native_oracle(monkeypatch: pytest.MonkeyPatch, workload: str) -> None:
    probe = _probe_generation(monkeypatch)
    path = generate_omega(workload)
    assert len(probe.built) == 1
    assert isinstance(probe.built[0], NativeDeclaredClearance)
    assert len(probe.refined_with) == 1
    assert probe.refined_with[0] is probe.built[0]
    assert path.generation_backend is ClearanceBackend.NATIVE


@pytest.mark.parametrize("workload", _SPACING_WORKLOADS)
def test_python_generation_refines_on_the_reference_oracle(monkeypatch: pytest.MonkeyPatch, workload: str) -> None:
    probe = _probe_generation(monkeypatch)
    path = generate_omega(workload, ClearanceBackend.PYTHON)
    assert probe.built == []
    assert len(probe.refined_with) == 1
    assert probe.refined_with[0] is path.clearance
    assert isinstance(path.clearance, DeclaredClearance)
    assert path.generation_backend is ClearanceBackend.PYTHON
