import hashlib
import importlib
import math

import numpy as np
import pytest


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
RECTANGLE = np.asarray(
    (
        (0.0, 0.0, 0.0),
        (6.0, 0.0, 0.0),
        (6.0, 6.0, 0.0),
        (0.0, 6.0, 0.0),
    ),
    dtype=np.float64,
)
CENTER_DOMAIN_DIGEST = bytes.fromhex("98aaa87a6fdc0589ef364e7ae3edad561208ff8d4245b340d616c63878be12de")
MAT_CERTIFICATE_DIGEST = bytes.fromhex("59ccf145c161b819ee91b047edba687eb9b4f0a697a4ed49574a8692f3e986dd")


def _native_module():
    return importlib.import_module("compas_cgal._medial_axis_2")


def _mat(
    vertices: np.ndarray = L_SHAPE,
    *,
    tool_radius: float = 0.5,
    station_spacing: float = 0.75,
    max_sagitta: float = 0.02,
    max_refinement_depth: int = 32,
) -> tuple[object, ...]:
    return _native_module().segment_site_medial_axis(
        vertices,
        [],
        tool_radius,
        station_spacing,
        max_sagitta,
        max_refinement_depth,
    )


def _owned_mat(
    vertices: np.ndarray = L_SHAPE,
    *,
    tool_radius: float = 0.5,
    station_spacing: float = 0.75,
    max_sagitta: float = 0.02,
    max_refinement_depth: int = 32,
):
    return _native_module().SegmentSiteMedialAxis.build(
        vertices,
        [],
        tool_radius,
        station_spacing,
        max_sagitta,
        max_refinement_depth,
    )


def _assert_array(value: object, dtype: np.dtype, shape: tuple[int, ...]) -> np.ndarray:
    assert isinstance(value, np.ndarray)
    assert value.dtype == dtype
    assert value.shape == shape
    return value


def _assert_same_result(left: tuple[object, ...], right: tuple[object, ...]) -> None:
    assert len(left) == len(right)
    for left_field, right_field in zip(left, right, strict=True):
        if isinstance(left_field, np.ndarray):
            assert isinstance(right_field, np.ndarray)
            assert np.array_equal(left_field, right_field)
        else:
            assert left_field == right_field


def test_segment_site_medial_axis_returns_fixed_proof_projection() -> None:
    result = _mat()

    assert isinstance(result, tuple)
    assert len(result) == 20
    (
        nodes,
        edges,
        node_site_offsets,
        node_site_ids,
        site_provenance,
        endpoint_flags,
        endpoint_feature_offsets,
        endpoint_features,
        edge_exact_flags,
        sample_centers,
        sample_clearance,
        sample_guide_radius,
        sample_flags,
        edge_sample_offsets,
        sample_parameter,
        neck_evidence,
        neck_cut_offsets,
        neck_cut_edge_ids,
        center_domain_digest,
        mat_certificate,
    ) = result

    nodes = _assert_array(nodes, np.dtype(np.float64), (10, 3))
    edges = _assert_array(edges, np.dtype(np.int64), (9, 8))
    node_site_offsets = _assert_array(node_site_offsets, np.dtype(np.int64), (11,))
    node_site_ids = _assert_array(node_site_ids, np.dtype(np.int64), (int(node_site_offsets[-1]),))
    site_provenance = _assert_array(site_provenance, np.dtype(np.int64), (12, 3))
    endpoint_flags = _assert_array(endpoint_flags, np.dtype(np.int64), (9, 2))
    endpoint_feature_offsets = _assert_array(endpoint_feature_offsets, np.dtype(np.int64), (19,))
    endpoint_features = _assert_array(endpoint_features, np.dtype(np.int64), (int(endpoint_feature_offsets[-1]), 5))
    edge_exact_flags = _assert_array(edge_exact_flags, np.dtype(np.int64), (9, 3))

    assert isinstance(sample_centers, np.ndarray)
    sample_count = sample_centers.shape[0]
    sample_centers = _assert_array(sample_centers, np.dtype(np.float64), (sample_count, 3))
    sample_clearance = _assert_array(sample_clearance, np.dtype(np.float64), (sample_count,))
    sample_guide_radius = _assert_array(sample_guide_radius, np.dtype(np.float64), (sample_count,))
    sample_flags = _assert_array(sample_flags, np.dtype(np.int64), (sample_count, 2))
    edge_sample_offsets = _assert_array(edge_sample_offsets, np.dtype(np.int64), (10,))
    sample_parameter = _assert_array(sample_parameter, np.dtype(np.float64), (sample_count,))
    neck_cut_offsets = _assert_array(neck_cut_offsets, np.dtype(np.int64), (3,))
    neck_cut_edge_ids = _assert_array(neck_cut_edge_ids, np.dtype(np.int64), (16,))

    assert np.all(np.isfinite(nodes))
    assert np.all(nodes[:, 2] == 0.0)
    assert np.all((0 <= edges[:, :2]) & (edges[:, :2] < nodes.shape[0]))
    assert set(edges[:, 2]) == {0, 1}
    assert np.all((0 <= edges[:, 3:5]) & (edges[:, 3:5] < site_provenance.shape[0]))
    assert np.array_equal(node_site_offsets, np.sort(node_site_offsets))
    assert np.all((0 <= node_site_ids) & (node_site_ids < site_provenance.shape[0]))
    assert np.array_equal(endpoint_feature_offsets, np.sort(endpoint_feature_offsets))
    assert np.all(endpoint_flags > 0)
    assert np.all(edge_exact_flags == 1)

    assert sample_count > edges.shape[0]
    assert np.all(np.isfinite(sample_centers))
    assert np.all(sample_centers[:, 2] == 0.0)
    assert np.all(sample_clearance >= 0.5)
    assert np.all(sample_guide_radius >= 0.0)
    assert np.all(sample_flags == 1)
    assert edge_sample_offsets[0] == 0
    assert edge_sample_offsets[-1] == sample_count
    assert np.array_equal(edge_sample_offsets, np.sort(edge_sample_offsets))
    assert np.all(np.isfinite(sample_parameter))

    assert isinstance(neck_evidence, tuple)
    assert len(neck_evidence) == 2
    assert all(isinstance(record, bytes) and b"mat-neck-plateau-v1" in record for record in neck_evidence)
    assert np.array_equal(neck_cut_offsets, np.asarray((0, 8, 16), dtype=np.int64))
    for neck_index in range(len(neck_evidence)):
        begin = neck_cut_offsets[neck_index]
        end = neck_cut_offsets[neck_index + 1]
        assert np.array_equal(neck_cut_edge_ids[begin:end], np.unique(neck_cut_edge_ids[begin:end]))
    assert center_domain_digest == CENTER_DOMAIN_DIGEST
    assert isinstance(mat_certificate, bytes)
    assert len(mat_certificate) == 124_796
    assert hashlib.sha256(mat_certificate).digest() == MAT_CERTIFICATE_DIGEST


def test_native_proof_owner_retains_exact_ids_behind_fixed_projection() -> None:
    owner = _owned_mat()
    projection = owner.projection

    _assert_same_result(projection, _mat())
    assert owner.node_ids == tuple(sorted(set(owner.node_ids)))
    assert owner.edge_ids == tuple(sorted(set(owner.edge_ids)))
    assert len(owner.node_ids) == projection[0].shape[0]
    assert len(owner.edge_ids) == projection[1].shape[0]
    assert len(owner.sample_parameter_ids) == projection[9].shape[0]
    assert len(owner.neck_owner_ids) == len(projection[15])
    assert all(type(identity) is bytes and identity for identity in owner.node_ids)
    assert all(type(identity) is bytes and identity for identity in owner.edge_ids)
    assert all(type(identity) is bytes and identity for identity in owner.sample_parameter_ids)
    assert all(type(identity) is bytes and identity for identity in owner.neck_owner_ids)

    projection[0][0, 0] = np.nextafter(projection[0][0, 0], math.inf)
    _assert_same_result(owner.projection, _mat())


def test_native_owner_proves_complete_radius_one_zero_guide_inventory() -> None:
    """The two width-2 S-S arms have clearance square identically equal to one."""
    owner = _owned_mat(tool_radius=1.0)
    records = owner.zero_guide_records

    assert len(records) == 2
    assert records == tuple(sorted(records))
    assert tuple(edge_id for edge_id, _ in records) == owner.edge_ids[2:4]
    assert owner.validate_zero_guide_records(
        owner.projection[19],
        records,
    ) == tuple(edge_id for edge_id, _ in records)


def test_native_owner_zero_guide_inventory_is_canonical() -> None:
    """Repeat construction and reversed ring input retain identical proof bytes."""
    canonical = _owned_mat(tool_radius=1.0).zero_guide_records
    repeated = _owned_mat(tool_radius=1.0).zero_guide_records
    reversed_ring = _owned_mat(L_SHAPE[::-1].copy(), tool_radius=1.0)

    assert canonical == repeated == reversed_ring.zero_guide_records


def test_native_owner_requires_complete_profile_identity() -> None:
    """Endpoint contact on P-S parabolas and binary64 near-equality prove nothing."""
    radius_one = _owned_mat(tool_radius=1.0)
    guide_radius = radius_one.projection[11]
    sample_offsets = radius_one.projection[13]
    zero_guide_edge_ids = {edge_id for edge_id, _ in radius_one.zero_guide_records}

    assert guide_radius[sample_offsets[0]] == 0.0
    assert radius_one.edge_ids[0] not in zero_guide_edge_ids
    assert _owned_mat(tool_radius=0.5).zero_guide_records == ()
    assert _owned_mat(tool_radius=np.nextafter(1.0, 0.0)).zero_guide_records == ()


def test_native_owner_rejects_missing_zero_guide_record() -> None:
    """Deleting one proved arm must fail complete-inventory replay."""
    native = _native_module()
    owner = _owned_mat(tool_radius=1.0)

    with pytest.raises(native.MissingMatZeroGuideEdgeError, match="missing"):
        owner.validate_zero_guide_records(
            owner.projection[19],
            owner.zero_guide_records[:-1],
        )


def test_native_owner_rejects_duplicate_zero_guide_edge() -> None:
    """One edge identity cannot authenticate two zero-guide record payloads."""
    native = _native_module()
    owner = _owned_mat(tool_radius=1.0)
    records = owner.zero_guide_records
    duplicate_with_foreign_bytes = records + ((records[0][0], records[1][1]),)

    with pytest.raises(native.DuplicateMatZeroGuideEdgeError, match="duplicate"):
        owner.validate_zero_guide_records(
            owner.projection[19],
            duplicate_with_foreign_bytes,
        )


def test_native_owner_rejects_mismatched_zero_guide_bytes() -> None:
    """Mutating one canonical byte must invalidate that edge's exact proof."""
    native = _native_module()
    owner = _owned_mat(tool_radius=1.0)
    records = owner.zero_guide_records
    mutated_bytes = bytearray(records[0][1])
    mutated_bytes[-1] ^= 1
    mutated = ((records[0][0], bytes(mutated_bytes)), records[1])

    with pytest.raises(native.MismatchedMatZeroGuideRecordError, match="bytes"):
        owner.validate_zero_guide_records(owner.projection[19], mutated)


def test_native_owner_rejects_foreign_zero_guide_edge() -> None:
    """A complete inventory cannot carry an unproved surplus edge record."""
    native = _native_module()
    owner = _owned_mat(tool_radius=1.0)
    foreign = owner.zero_guide_records + ((b"foreign-edge", owner.zero_guide_records[0][1]),)

    with pytest.raises(native.MismatchedMatZeroGuideRecordError, match="foreign"):
        owner.validate_zero_guide_records(owner.projection[19], foreign)


def test_native_owner_rejects_foreign_mat_certificate() -> None:
    """Zero-guide records cannot be replayed under another tool-radius MAT."""
    native = _native_module()
    owner = _owned_mat(tool_radius=1.0)
    foreign_mat_certificate = _owned_mat(tool_radius=0.5).projection[19]

    with pytest.raises(native.InvalidMatCertificateReplayError, match="replay"):
        owner.validate_zero_guide_records(
            foreign_mat_certificate,
            owner.zero_guide_records,
        )


def test_typed_medial_axis_projects_native_topology_without_coordinate_matching() -> None:
    from compas_cgal.adaptive.canonical import CanonicalRingV1
    from compas_cgal.adaptive.medial_axis import MedialAxis
    from compas_cgal.adaptive.units import ChordBound
    from compas_cgal.adaptive.units import Point2
    from compas_cgal.adaptive.units import Spacing
    from compas_cgal.adaptive.units import ToolRadius
    from compas_cgal.adaptive.units import WorldXY

    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in L_SHAPE),
    )
    axis = MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )

    assert tuple(node.identity for node in axis.nodes) == axis.native_owner.node_ids
    assert tuple(edge.identity for edge in axis.edges) == axis.native_owner.edge_ids
    assert len(axis.nodes) == 10
    assert len(axis.edges) == 9
    assert len(axis.samples) > len(axis.edges)
    assert len(axis.tool_fit_runs) == len(axis.edges)
    assert {run.parent_edge_id for run in axis.tool_fit_runs} == {edge.identity for edge in axis.edges}
    assert all(edge.source.identity in axis.node_by_id for edge in axis.edges)
    assert all(edge.target.identity in axis.node_by_id for edge in axis.edges)
    assert all(edge.generator_site_ids == tuple(sorted(edge.generator_site_ids)) for edge in axis.edges)
    assert all(sample.edge_id in axis.edge_by_id for sample in axis.samples)
    assert all(sample.exact_parameter_id for sample in axis.samples)
    assert axis.center_domain_digest == CENTER_DOMAIN_DIGEST
    assert hashlib.sha256(axis.mat_certificate).digest() == MAT_CERTIFICATE_DIGEST


def test_typed_topology_rejects_mutable_lookup_views() -> None:
    from compas_cgal.adaptive.canonical import CanonicalRingV1
    from compas_cgal.adaptive.errors import InvalidMedialAxisProjectionError
    from compas_cgal.adaptive.medial_axis import MatTopology
    from compas_cgal.adaptive.medial_axis import MedialAxis
    from compas_cgal.adaptive.units import ChordBound
    from compas_cgal.adaptive.units import Point2
    from compas_cgal.adaptive.units import Spacing
    from compas_cgal.adaptive.units import ToolRadius
    from compas_cgal.adaptive.units import WorldXY

    boundary = CanonicalRingV1.build_outer(
        tuple(Point2[WorldXY].build(x, y) for x, y, _ in L_SHAPE),
    )
    axis = MedialAxis.build(
        design_boundary=boundary,
        holes=(),
        tool_radius=ToolRadius.build(0.5),
        station_spacing=Spacing.build(0.75),
        max_sagitta=ChordBound.build(0.02),
        max_refinement_depth=32,
    )
    topology = axis.topology

    with pytest.raises(InvalidMedialAxisProjectionError, match="immutable"):
        MatTopology(
            topology.sites,
            topology.nodes,
            topology.edges,
            topology.components,
            topology.tool_fit_runs,
            dict(topology.site_by_id),
            topology.node_by_id,
            topology.edge_by_id,
            topology.component_by_edge_id,
        )


def test_segment_site_medial_axis_is_bitwise_canonical_under_ring_reversal() -> None:
    canonical = _mat()
    repeated = _mat()
    reversed_ring = _mat(L_SHAPE[::-1].copy())

    _assert_same_result(canonical, repeated)
    _assert_same_result(canonical, reversed_ring)


def test_mat_certificate_and_exact_topology_ignore_proposal_refinement() -> None:
    coarse = _mat(station_spacing=0.75)
    fine = _mat(station_spacing=0.25)

    assert coarse[9].shape[0] != fine[9].shape[0]
    for field_index in range(9):
        assert np.array_equal(coarse[field_index], fine[field_index])
    for field_index in range(15, 20):
        coarse_field = coarse[field_index]
        fine_field = fine[field_index]
        if isinstance(coarse_field, np.ndarray):
            assert isinstance(fine_field, np.ndarray)
            assert np.array_equal(coarse_field, fine_field)
        else:
            assert coarse_field == fine_field


def test_exact_half_passage_radius_removes_separating_necks() -> None:
    result = _mat(tool_radius=1.0)

    assert result[0].shape == (6, 3)
    assert result[1].shape == (5, 8)
    assert np.all(result[8] == 1)
    assert result[15] == ()
    assert np.array_equal(result[16], np.asarray((0,), dtype=np.int64))
    assert result[17].shape == (0,)
    assert len(result[18]) == 32
    assert isinstance(result[19], bytes)
    assert result[19]


def test_binding_exposes_named_input_and_sampling_failures() -> None:
    native = _native_module()

    with pytest.raises(native.InvalidReachableDomainInputError, match="finite and positive"):
        _mat(tool_radius=0.0)
    with pytest.raises(native.InvalidMatSamplingPolicyError, match="station spacing"):
        _mat(station_spacing=0.0)
    with pytest.raises(native.InvalidMatSamplingPolicyError, match="sagitta"):
        _mat(max_sagitta=0.0)
    with pytest.raises(native.ConicSamplingLimitError, match="refinement depth"):
        _mat(max_refinement_depth=0)
    with pytest.raises(native.UnsupportedCanonicalMatLShapeGraphError, match="L-shape"):
        _mat(RECTANGLE)


def test_binding_has_one_admissible_center_radius_argument() -> None:
    native = _native_module()

    with pytest.raises(TypeError):
        native.segment_site_medial_axis(
            L_SHAPE,
            [],
            0.5,
            0.5,
            0.75,
            0.02,
            32,
        )
