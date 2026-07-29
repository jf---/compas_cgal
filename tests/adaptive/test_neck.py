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


def _reversed_owner() -> _medial_axis_2.SegmentSiteMedialAxis:
    return _medial_axis_2.SegmentSiteMedialAxis.build(
        np.ascontiguousarray(L_SHAPE[::-1]),
        [],
        0.5,
        0.75,
        0.02,
        32,
    )


def _axis(
    station_spacing: float = 0.75,
    *,
    reverse: bool = False,
) -> MedialAxis:
    points = np.ascontiguousarray(L_SHAPE[::-1]) if reverse else L_SHAPE
    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in points),
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


def test_native_owner_projects_exact_neck_loci_and_every_cut_partition() -> None:
    """Retain three-way plateau topology instead of only its cut-edge union.

    Both production L-pocket necks separate a central side and two leaf
    sides. Flattening those partitions makes causal Task 13 traversal
    impossible even though every individual cut edge remains present.
    """
    owner = _owner()
    reversed_owner = _reversed_owner()

    assert owner.neck_location_tags == (
        b"mat-neck-plateau-v1",
        b"mat-neck-plateau-v1",
    )
    assert owner.neck_location_tags == reversed_owner.neck_location_tags
    assert owner.neck_location_edge_ids == reversed_owner.neck_location_edge_ids
    assert owner.neck_location_node_ids == reversed_owner.neck_location_node_ids
    assert owner.neck_parameter_root_ids == (None, None)
    assert owner.neck_parameter_root_ids == reversed_owner.neck_parameter_root_ids
    assert owner.neck_cut_edge_partitions == reversed_owner.neck_cut_edge_partitions

    edge_ids = set(owner.edge_ids)
    projection = owner.projection
    for neck_index, partitions in enumerate(owner.neck_cut_edge_partitions):
        assert tuple(len(partition) for partition in partitions) == (6, 1, 1)
        assert all(partition == tuple(sorted(partition)) for partition in partitions)
        assert all(len(partition) == len(set(partition)) for partition in partitions)
        assert all(edge_id in edge_ids for partition in partitions for edge_id in partition)

        partition_union = tuple(sorted({edge_id for partition in partitions for edge_id in partition}))
        begin = int(projection[16][neck_index])
        end = int(projection[16][neck_index + 1])
        projected_union = tuple(owner.edge_ids[int(edge_index)] for edge_index in projection[17][begin:end])
        assert partition_union == projected_union

    assert tuple(len(edge_ids) for edge_ids in owner.neck_location_edge_ids) == (1, 1)
    assert tuple(len(node_ids) for node_ids in owner.neck_location_node_ids) == (2, 2)


def test_typed_neck_inventory_binds_exact_evidence_and_legal_passage_steps() -> None:
    from compas_cgal.adaptive.errors import TerminalNeckPassageError
    from compas_cgal.adaptive.neck import NeckInventory
    from compas_cgal.adaptive.neck import NeckSide
    from compas_cgal.adaptive.neck import PlateauNeckLocus

    policy = _neck_policy()
    inventory = NeckInventory.build(axis=_axis(), policy=policy)
    reversed_inventory = NeckInventory.build(
        axis=_axis(reverse=True),
        policy=policy,
    )

    assert len(inventory.necks) == 2
    assert tuple(neck.canonical_bytes for neck in inventory.necks) == tuple(neck.canonical_bytes for neck in reversed_inventory.necks)
    assert tuple(tuple(side.identity for side in neck.sides) for neck in inventory.necks) == tuple(
        tuple(side.identity for side in neck.sides) for neck in reversed_inventory.necks
    )
    assert tuple(neck.owner_id for neck in inventory.necks) == inventory.axis.native_owner.neck_owner_ids
    for neck in inventory.necks:
        assert neck.width_class_id.value == 0
        assert len(neck.defining_site_ids) == 2
        assert len(neck.separating_cut_edge_ids) == 8
        assert hashlib.sha256(neck.evidence_bytes).digest() == neck.evidence_digest
        assert neck.evidence_digest in neck.comparison_certificate
        assert type(neck.locus) is PlateauNeckLocus
        assert len(neck.locus.node_ids) == 2
        assert len(neck.locus.edge_ids) == 1
        assert tuple(len(side.edge_ids) for side in neck.sides) == (6, 1, 1)
        assert all(type(side) is NeckSide for side in neck.sides)
        assert tuple(side.partition_ordinal for side in neck.sides) == (0, 1, 2)
        assert all(side.neck_evidence_digest == neck.evidence_digest for side in neck.sides)
        assert tuple(sorted({edge_id for side in neck.sides for edge_id in side.edge_ids})) == neck.separating_cut_edge_ids
        assert neck.locus.canonical_bytes in neck.canonical_bytes
        assert all(side.canonical_bytes in neck.canonical_bytes for side in neck.sides)

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


def test_typed_neck_inventory_rejects_repartitioned_or_crosswired_native_topology() -> None:
    """Reject valid MAT IDs rearranged into a topology native evidence did not prove.

    A partition mutation can preserve the exact cut union and use only known
    edges. Union membership and foreign-ID checks therefore cannot detect it;
    the typed inventory must compare complete sides and locus to its retained
    native evidence owner.
    """
    from compas_cgal.adaptive.errors import InvalidNeckEvidenceError
    from compas_cgal.adaptive.neck import NeckInventory
    from compas_cgal.adaptive.neck import NeckSide
    from compas_cgal.adaptive.neck import PlateauNeckLocus

    inventory = NeckInventory.build(axis=_axis(), policy=_neck_policy())
    neck = inventory.necks[0]

    with pytest.raises(InvalidNeckEvidenceError, match="canonical side order"):
        replace(neck, sides=tuple(reversed(neck.sides)))

    repartitioned_sides = (
        NeckSide.build(
            neck_evidence_digest=neck.evidence_digest,
            partition_ordinal=0,
            edge_ids=tuple(sorted((*neck.sides[0].edge_ids[:-1], *neck.sides[1].edge_ids))),
        ),
        NeckSide.build(
            neck_evidence_digest=neck.evidence_digest,
            partition_ordinal=1,
            edge_ids=(neck.sides[0].edge_ids[-1],),
        ),
        neck.sides[2],
    )
    repartitioned = replace(neck, sides=repartitioned_sides)
    with pytest.raises(InvalidNeckEvidenceError, match="native neck topology"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (repartitioned, *inventory.necks[1:]),
        )

    foreign_locus = PlateauNeckLocus.build(
        node_ids=neck.locus.node_ids,
        edge_ids=(inventory.axis.edges[0].identity,),
    )
    crosswired = replace(neck, locus=foreign_locus)
    with pytest.raises(InvalidNeckEvidenceError, match="native neck topology"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (crosswired, *inventory.necks[1:]),
        )


def test_neck_inventory_rejects_topology_from_another_mat() -> None:
    from compas_cgal.adaptive.errors import InvalidNeckEvidenceError
    from compas_cgal.adaptive.medial_axis import MatEdgeId
    from compas_cgal.adaptive.medial_axis import MatSiteId
    from compas_cgal.adaptive.neck import NeckInventory
    from compas_cgal.adaptive.neck import NeckSide
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
    foreign_side = NeckSide.build(
        neck_evidence_digest=neck.evidence_digest,
        partition_ordinal=neck.sides[0].partition_ordinal,
        edge_ids=tuple(sorted((unknown_edge, *neck.sides[0].edge_ids[1:]))),
    )
    foreign_sides = (foreign_side, *neck.sides[1:])
    foreign_cut_neck = replace(
        neck,
        separating_cut_edge_ids=tuple(sorted({edge_id for side in foreign_sides for edge_id in side.edge_ids})),
        sides=foreign_sides,
    )
    with pytest.raises(InvalidNeckEvidenceError, match="unknown MAT edge"):
        NeckInventory(
            inventory.axis,
            inventory.policy,
            (foreign_cut_neck, *inventory.necks[1:]),
        )
