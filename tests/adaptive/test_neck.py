import hashlib
import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _medial_axis_2
from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.policy import ACTIVE_PASSAGE_STATES
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import WorldXY

L_SHAPE = np.asarray(
    (
        (0.0, 0.0, 0.0),
        (6.0, 0.0, 0.0),
        (6.0, 2.0, 0.0),
        (2.0, 2.0, 0.0),
        (2.0, 6.0, 0.0),
        (0.0, 6.0, 0.0),
    ),
    dtype=np.float64,
)


def _owner() -> _medial_axis_2.SegmentSiteMedialAxis:
    return _medial_axis_2.SegmentSiteMedialAxis.build(
        L_SHAPE,
        [],
        0.5,
        0.75,
        0.02,
        32,
    )


def _axis(station_spacing: float = 0.75) -> MedialAxis:
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in L_SHAPE),
    )
    return MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(station_spacing),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )


def _neck_policy() -> NeckPolicy:
    user_cap = EngagementCap.build(math.radians(120.0))
    caps = {
        (neck_class, passage_state): EngagementCap.build(
            math.radians(90.0 - 10.0 * passage_state.rank),
        )
        for neck_class in range(2)
        for passage_state in ACTIVE_PASSAGE_STATES
    }
    return NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps=caps,
    )


def test_native_neck_classification_replays_owner_before_exact_comparison() -> None:
    owner = _owner()
    projection = owner.projection
    boundary = ExactRationalV1.from_fraction(Fraction(4)).canonical_bytes

    classes, comparison_certificates = owner.validate_and_classify_necks(
        projection[19],
        projection[15],
        (boundary,),
    )

    assert isinstance(classes, np.ndarray)
    assert classes.dtype == np.dtype(np.int64)
    assert np.array_equal(classes, np.asarray((0, 0), dtype=np.int64))
    assert len(comparison_certificates) == len(projection[15])
    for evidence, certificate in zip(
        projection[15],
        comparison_certificates,
        strict=True,
    ):
        assert type(certificate) is bytes
        assert b"mat-neck-width-comparison-v1" in certificate
        assert hashlib.sha256(evidence).digest() in certificate


def test_native_neck_classification_rejects_any_broken_proof_input() -> None:
    owner = _owner()
    projection = owner.projection
    boundary = ExactRationalV1.from_fraction(Fraction(4)).canonical_bytes

    mutated_certificate = bytearray(projection[19])
    mutated_certificate[-1] ^= 1
    with pytest.raises(_medial_axis_2.InvalidMatCertificateReplayError):
        owner.validate_and_classify_necks(
            bytes(mutated_certificate),
            projection[15],
            (boundary,),
        )

    mutated_evidence = bytearray(projection[15][0])
    mutated_evidence[-1] ^= 1
    with pytest.raises(_medial_axis_2.InvalidMatNeckEvidenceError):
        owner.validate_and_classify_necks(
            projection[19],
            (bytes(mutated_evidence), *projection[15][1:]),
            (boundary,),
        )

    with pytest.raises(_medial_axis_2.InvalidMatNeckWidthBoundariesError):
        owner.validate_and_classify_necks(
            projection[19],
            projection[15],
            (b"not-an-exact-rational",),
        )


def test_typed_neck_inventory_binds_exact_evidence_and_legal_passage_steps() -> None:
    from compas_cgal.adaptive.errors import TerminalNeckPassageError
    from compas_cgal.adaptive.neck import NeckInventory

    policy = _neck_policy()
    inventory = NeckInventory.build(axis=_axis(), policy=policy)

    assert len(inventory.necks) == 2
    assert tuple(neck.owner_id for neck in inventory.necks) == inventory.axis.native_owner.neck_owner_ids
    for neck in inventory.necks:
        assert neck.width_class_id.value == 0
        assert len(neck.defining_site_ids) == 2
        assert len(neck.separating_cut_edge_ids) == 8
        assert hashlib.sha256(neck.evidence_bytes).digest() == neck.evidence_digest
        assert neck.evidence_digest in neck.comparison_certificate

        passage = neck.forward
        expected_states = (
            PassageState.UNVISITED,
            PassageState.FIRST_PASS_COMPLETE,
            PassageState.SECOND_PASS_COMPLETE,
            PassageState.TERMINAL,
        )
        for passage_before, passage_after in zip(
            expected_states[:-1],
            expected_states[1:],
            strict=True,
        ):
            assert passage.state is passage_before
            decision = passage.propose_cap_decision(policy)
            assert decision.passage_before is passage_before
            assert decision.passage_after is passage_after
            assert decision.neck_evidence_digest == neck.evidence_digest
            passage = passage.advance(decision)
        assert passage.state is PassageState.TERMINAL
        with pytest.raises(TerminalNeckPassageError):
            passage.propose_cap_decision(policy)


def test_neck_inventory_rejects_topology_from_another_mat() -> None:
    from compas_cgal.adaptive.errors import InvalidNeckEvidenceError
    from compas_cgal.adaptive.medial_axis import MatEdgeId
    from compas_cgal.adaptive.medial_axis import MatSiteId
    from compas_cgal.adaptive.neck import NeckInventory
    from compas_cgal.adaptive.operation import NeckOwnerId

    inventory = NeckInventory.build(axis=_axis(), policy=_neck_policy())
    neck = inventory.necks[0]
    foreign_owner_neck = replace(
        inventory.necks[-1],
        owner_id=NeckOwnerId(b"\xff" * hashlib.sha256().digest_size),
    )
    with pytest.raises(InvalidNeckEvidenceError, match="exact MAT owner"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (*inventory.necks[:-1], foreign_owner_neck),
        )

    unknown_site = MatSiteId(hashlib.sha256(b"foreign-mat-site").digest())
    foreign_site_neck = replace(
        neck,
        defining_site_ids=tuple(sorted((unknown_site, neck.defining_site_ids[0]))),
    )

    with pytest.raises(InvalidNeckEvidenceError, match="unknown MAT site"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (foreign_site_neck, *inventory.necks[1:]),
        )

    unknown_edge = MatEdgeId(hashlib.sha256(b"foreign-mat-edge").digest())
    foreign_cut_neck = replace(
        neck,
        separating_cut_edge_ids=tuple(
            sorted((unknown_edge, *neck.separating_cut_edge_ids[1:])),
        ),
    )
    with pytest.raises(InvalidNeckEvidenceError, match="unknown MAT edge"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (foreign_cut_neck, *inventory.necks[1:]),
        )
