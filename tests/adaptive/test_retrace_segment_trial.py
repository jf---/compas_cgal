"""Physical proof-path contracts shared by advancing and retrace segments."""

from typing import NoReturn

import pytest

from compas_cgal.adaptive.containment import GougeContainment
from compas_cgal.adaptive.coverage import CoverageLedger
from compas_cgal.adaptive.generator import advance_active_candidate_family
from compas_cgal.adaptive.motion_certificate import SWEPT_PREFIX_MOTION_STRATA
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import SweptPrefixMotionWitness
from compas_cgal.adaptive.operation import AdvanceSegmentOperation
from compas_cgal.adaptive.operation import RetraceSegmentOperation
from compas_cgal.adaptive.stock_area import Stock2Area
from tests.adaptive.task13f_fixture import Task13FFixture
from tests.adaptive.task13f_fixture import task13f_retrace_decision
from tests.adaptive.task13f_fixture import task13f_route_one_terminal


@pytest.fixture(scope="module")
def task13f() -> Task13FFixture:
    """Build one authenticated Task 13F authority for physical proof trials."""
    return Task13FFixture.build()


def test_shared_trial_preserves_established_advancing_identity(
    task13f: Task13FFixture,
) -> None:
    """Keep the accepted route-1 proof byte-identical while sharing its theorem.

    The real Task 13F route-one transaction fixes the advancing witness,
    transaction, and commit identities before the deciding path is factored.

    Args:
        task13f: Authenticated launch child and continuation authority.
    """
    physical, traversal, _ = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=task13f.physical,
        traversal=task13f.traversal,
    )
    traversal = traversal.activate_next()
    _, _, route_one_commit = advance_active_candidate_family(
        evaluator=task13f.evaluator,
        physical=physical,
        traversal=traversal,
    )

    witness = route_one_commit.transaction.segment_witness
    assert witness.digest.hex() == ("0d7f0ac6c7b58ebe114b5c5453ce12b014a990687bda16bfcffa199874762cc7")
    assert route_one_commit.transaction.digest.hex() == ("5668c64ae62f1aa1b52f7338db626626ec0636eef66c731934917b9b331a40b4")
    assert route_one_commit.digest.hex() == ("0a20db00bf444336ffd0138da9e4a938f28dd8bf7bc0baacb69ebdc8f1656c73")


def _retrace_trial_arguments(
    task13f: Task13FFixture,
) -> dict[str, object]:
    """Build reverse-trial arguments from the accepted route-one prefix.

    Args:
        task13f: Authenticated launch child and continuation authority.

    Returns:
        Independent physical owners and the exact source-derived retrace.
    """
    physical, terminal, commits = task13f_route_one_terminal(task13f)
    source_commit = commits[-1]
    source = physical.operations[-1]
    assert type(source) is AdvanceSegmentOperation
    decision = task13f_retrace_decision(
        physical=physical,
        terminal=terminal,
        source_commit=source_commit,
    )
    return {
        "containment_authority": GougeContainment.build(
            task13f.identity.reachable_domain,
        ),
        "stock": physical.fork_stock(),
        "coverage": physical.clone_coverage(),
        "operation_index": len(physical.operations),
        "operation": RetraceSegmentOperation.build(
            source_operation=source,
            decision=decision,
        ),
        "tool_radius": task13f.identity.tool_radius,
        "user_cap": task13f.identity.user_cap,
        "effective_cap": task13f.identity.user_cap,
        "depletion_policy": task13f.identity.depletion_policy,
    }


def test_retrace_uses_one_ordered_two_stratum_swept_prefix_proof(
    task13f: Task13FFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Certify the reverse corridor once before either trial owner mutates.

    The retrace must use the same exact two-stratum theorem as its accepted
    source. The generic event certifier is therefore a forbidden deciding
    path, and both trial-local mutations must follow the proof in fixed order.

    Args:
        task13f: Authenticated launch child and continuation authority.
        monkeypatch: Scoped instrumentation for the four physical authorities.
    """
    from compas_cgal.adaptive.retrace_segment_trial import (
        evaluate_retrace_segment_trial,
    )

    arguments = _retrace_trial_arguments(task13f)
    order: list[str] = []
    calls = {"generic": 0, "swept_prefix": 0}
    real_containment = GougeContainment.certify_segment
    real_swept_prefix = MotionCertifier.certify_swept_prefix_segment
    real_deplete = Stock2Area.deplete
    real_coverage = CoverageLedger.add_sweep

    def tracked_containment(self: GougeContainment, *args: object) -> object:
        order.append("containment")
        return real_containment(self, *args)  # type: ignore[arg-type]

    def tracked_swept_prefix(
        self: MotionCertifier,
        **kwargs: object,
    ) -> SweptPrefixMotionWitness:
        order.append("swept-prefix")
        calls["swept_prefix"] += 1
        return real_swept_prefix(self, **kwargs)  # type: ignore[arg-type]

    def reject_generic(self: MotionCertifier, **kwargs: object) -> NoReturn:
        calls["generic"] += 1
        raise AssertionError("generic event certification is forbidden")

    def tracked_deplete(self: Stock2Area, *args: object) -> object:
        order.append("deplete")
        return real_deplete(self, *args)  # type: ignore[arg-type]

    def tracked_coverage(self: CoverageLedger, *args: object) -> object:
        order.append("coverage")
        return real_coverage(self, *args)  # type: ignore[arg-type]

    monkeypatch.setattr(GougeContainment, "certify_segment", tracked_containment)
    monkeypatch.setattr(
        MotionCertifier,
        "certify_swept_prefix_segment",
        tracked_swept_prefix,
    )
    monkeypatch.setattr(MotionCertifier, "certify", reject_generic)
    monkeypatch.setattr(Stock2Area, "deplete", tracked_deplete)
    monkeypatch.setattr(CoverageLedger, "add_sweep", tracked_coverage)

    witness = evaluate_retrace_segment_trial(
        **arguments,
    )  # type: ignore[arg-type]

    assert order == ["containment", "swept-prefix", "deplete", "coverage"]
    assert calls == {"generic": 0, "swept_prefix": 1}
    assert type(witness.motion_witness) is SweptPrefixMotionWitness
    assert witness.motion_witness.event_cell_count == SWEPT_PREFIX_MOTION_STRATA == 2
