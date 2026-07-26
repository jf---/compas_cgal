from collections.abc import Sequence
from typing import Literal
from typing import TypeAlias

from compas_cgal._stock_2 import Stock2

BoundaryEventKind: TypeAlias = Literal[
    "transverse",
    "tangent",
    "vertex",
    "overlap",
    "seam",
]

BOUNDARY_EVENT_KINDS: tuple[BoundaryEventKind, ...]


class BoundaryExtractionError(RuntimeError): ...


class DegenerateBoundarySupportError(BoundaryExtractionError): ...


class MissingBoundaryEndpointError(BoundaryExtractionError): ...


class MissingBoundaryIntersectionError(BoundaryExtractionError): ...


class BoundaryFeatureIndexError(BoundaryExtractionError): ...


class EventSubstrateError(RuntimeError): ...


class AlgebraicBackendError(EventSubstrateError): ...


class InvalidAlgebraicPolynomialError(EventSubstrateError): ...


class UnsupportedAlgebraicDegreeError(EventSubstrateError): ...


class AlgebraicRootIsolationError(EventSubstrateError): ...


class ChartCoverageError(EventSubstrateError): ...


class ProjectionDegreeBoundError(EventSubstrateError): ...


class TrimFilterError(EventSubstrateError): ...


class TrimEndpointOffSupportError(TrimFilterError): ...


class DegenerateTrimError(TrimFilterError): ...


class UnsupportedLineMotionError(TrimFilterError): ...


class EventPartitionVerificationError(EventSubstrateError): ...


class EventTraceVerificationError(EventSubstrateError): ...


class NonFiniteSegmentInputError(EventSubstrateError): ...


class ZeroLengthSegmentMotionError(EventSubstrateError): ...


class NonPositiveToolRadiusError(EventSubstrateError): ...


class InvalidCapChordRatioError(EventSubstrateError): ...


class UnsupportedAlgebraicVertexProjectionError(EventSubstrateError): ...


class IncompleteSegmentPartitionError(EventSubstrateError): ...


class ContinuousTeaVerdict:
    CERTIFIED: ContinuousTeaVerdict
    CAP_EXCEEDED: ContinuousTeaVerdict
    UNRESOLVED_DEGENERACY: ContinuousTeaVerdict
    @property
    def name(self) -> str: ...


class ExactBinary64Rational2:
    numerator: str
    denominator: str
    text: str
    canonical_bytes: bytes


class SegmentEventSource2:
    @staticmethod
    def from_binary64(
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        tool_radius: float,
        cap_chord_ratio: float,
    ) -> SegmentEventSource2: ...
    x0: ExactBinary64Rational2
    y0: ExactBinary64Rational2
    x1: ExactBinary64Rational2
    y1: ExactBinary64Rational2
    tool_radius: ExactBinary64Rational2
    cap_chord_ratio: ExactBinary64Rational2
    motion_data: tuple[str, str, str, str]
    canonical_bytes: bytes
    canonical_digest: bytes


class BoundaryFeatureRecord2:
    @property
    def support_kind(self) -> Literal["line", "circle"]: ...
    @property
    def support_coefficients(self) -> Sequence[bytes]: ...
    @property
    def primitive_coefficients(self) -> Sequence[str]: ...
    @property
    def support_id(self) -> bytes: ...
    @property
    def source_exact(self) -> bytes: ...
    @property
    def target_exact(self) -> bytes: ...
    @property
    def source_vertex_id(self) -> bytes: ...
    @property
    def target_vertex_id(self) -> bytes: ...
    @property
    def source_incidence(self) -> bytes: ...
    @property
    def target_incidence(self) -> bytes: ...
    @property
    def material_side(self) -> Literal["left"]: ...
    @property
    def trim_predicate(self) -> bytes: ...
    @property
    def feature_id(self) -> bytes: ...
    @property
    def overlap_multiplicity(self) -> int: ...


class BoundaryEvent2:
    kind: BoundaryEventKind
    @property
    def first_feature_id(self) -> bytes: ...
    @property
    def second_feature_id(self) -> bytes: ...
    @property
    def vertex_id(self) -> bytes: ...
    @property
    def overlap_id(self) -> bytes: ...
    @property
    def exact_overlap_record(self) -> bytes: ...
    multiplicity: int


class AlgebraicBackendEvidence2:
    cgal_version: str
    integer_backend: str
    algebraic_kernel_1: str
    algebraic_kernel_2: str
    arrangement_traits: str
    compile_definitions: tuple[str, ...]


class ParameterChart2:
    chart_id: str
    family: str
    domain_low: str
    domain_high: str
    orientation: Literal["ccw"]
    start_seam_id: bytes
    end_seam_id: bytes
    owns_start_seam: bool
    owns_end_seam: bool


class PartitionEvent2:
    def __init__(
        self,
        kind: str,
        feature_id: bytes,
        support_id: bytes,
        trim_id: bytes,
        vertex_id: bytes,
        branch_id: bytes,
        disposition: str,
    ) -> None: ...
    kind: str
    feature_id: bytes
    support_id: bytes
    trim_id: bytes
    vertex_id: bytes
    branch_id: bytes
    disposition: str
    left_active_count: int
    right_active_count: int
    incidence_permutation_rechecked: bool
    original_equations_rechecked: bool
    orientation_rechecked: bool
    trim_disposition: str
    pair_sheet_id: bytes
    first_feature_id: bytes
    second_feature_id: bytes
    first_chart_id: str
    second_chart_id: str
    first_branch_id: bytes
    second_branch_id: bytes


class ProjectionInput2:
    def __init__(
        self,
        projection_id: str,
        coefficients: Sequence[str],
        events: Sequence[PartitionEvent2],
    ) -> None: ...


class AlgebraicRootRecord2:
    def __init__(
        self,
        root_id: bytes,
        factor_coefficients: Sequence[str],
        root_ordinal: int,
        multiplicity: int,
        interval_low: str,
        interval_high: str,
    ) -> None: ...
    root_id: bytes
    factor_coefficients: tuple[str, ...]
    root_ordinal: int
    multiplicity: int
    interval_low: str
    interval_high: str


class ParameterCell2:
    lower_root_id: bytes
    upper_root_id: bytes
    witness_numerator: str
    witness_denominator: str
    disposition: str


class ActiveBoundaryBranch2:
    branch_id: bytes
    feature_id: bytes
    support_id: bytes
    trim_id: bytes
    chart_id: str
    sheet_ordinal: int
    root_id: bytes


class EventFibre2:
    root_id: bytes
    local_root_ids: Sequence[bytes]
    seam_id: bytes
    events: Sequence[PartitionEvent2]
    left_active_branches: Sequence[ActiveBoundaryBranch2]
    right_active_branches: Sequence[ActiveBoundaryBranch2]
    ccw_direction: Literal["merge", "split", "unchanged", ""]
    cw_direction: Literal["merge", "split", "unchanged", ""]


class RationalTrimInterval2:
    rim_parameter: str
    motion_domain_low: str
    motion_domain_high: str


class TrimmedLineBranch2:
    rim_parameter: str
    motion_domain_low: str
    motion_domain_high: str
    rational_convenience: RationalTrimInterval2 | None
    rim_root: AlgebraicRootRecord2
    lower_trim_predicate_rows: Sequence[Sequence[str]]
    upper_trim_predicate_rows: Sequence[Sequence[str]]
    rational_convenience_available: bool
    domain_nonempty_rechecked: bool
    complementarity_rechecked: bool
    trim_disposition: str
    rejected_outside_closed_domain: bool
    feature_id: bytes
    local_support_id: bytes
    local_trimmed_feature_id: bytes
    trim_id: bytes
    branch_id: bytes


class ProjectedRegularizationVertex2:
    root: AlgebraicRootRecord2
    vertex_id: bytes
    first_trim_disposition: str
    second_trim_disposition: str
    conjugate_disposition: str


class CcwOrientation2:
    disposition: str
    determinant_sign: str
    dot_sign: str


class ProjectionRecord2:
    projection_id: str
    coefficient_rows: Sequence[Sequence[str]]
    factor_coefficients: Sequence[Sequence[str]]
    actual_degree: tuple[int, int]
    degree_bound: tuple[int, int]
    degree_bound_id: str
    normalized_coefficient_bytes: bytes
    signed_predicate_coefficients: Sequence[str]


class OverlapInterval2:
    kind: str
    domain_low: str
    domain_high: str
    witness_numerator: str
    witness_denominator: str
    orientation_disposition: str
    feature_id: bytes
    support_id: bytes
    trim_id: bytes
    branch_id: bytes


class ChartSeam2:
    seam_id: bytes
    owner_chart_id: str


class EventPartitionCertificate2:
    def __init__(
        self,
        *,
        build_evidence: AlgebraicBackendEvidence2,
        charts: Sequence[ParameterChart2],
        projections: Sequence[ProjectionRecord2],
        roots: Sequence[AlgebraicRootRecord2],
        cells: Sequence[ParameterCell2],
        fibres: Sequence[EventFibre2],
        overlaps: Sequence[OverlapInterval2],
        seams: Sequence[ChartSeam2],
        source_kind: str,
        source_payload: bytes,
    ) -> None: ...
    build_evidence: AlgebraicBackendEvidence2
    charts: Sequence[ParameterChart2]
    projections: Sequence[ProjectionRecord2]
    roots: Sequence[AlgebraicRootRecord2]
    cells: Sequence[ParameterCell2]
    fibres: Sequence[EventFibre2]
    overlaps: Sequence[OverlapInterval2]
    seams: Sequence[ChartSeam2]
    source_kind: str
    source_payload: bytes
    canonical_bytes: bytes
    canonical_digest: bytes


class VerifiedEventPartition2:
    verdict: ContinuousTeaVerdict
    partition: EventPartitionCertificate2


class SegmentBoundaryBranch2:
    branch_id: bytes
    feature_id: bytes
    support_id: bytes
    support_kind: Literal["line", "circle"]
    trim_id: bytes
    vertex_id: bytes
    material_side: Literal["left"]
    trim_disposition: Literal["accepted"]
    rim_chart_id: Literal["rim-half-0-v1", "rim-half-1-v1"]
    rim_sheet_ordinal: int
    rim_root_id: bytes
    rim_factor_coefficients: Sequence[str]
    rim_root_ordinal: int


class BranchPairDisposition2:
    pair_sheet_id: bytes
    first_branch_id: bytes
    second_branch_id: bytes
    orientation_disposition: Literal["ccw-forward", "ccw-wrap"]
    cap_disposition: Literal["below-cap", "equal-cap", "above-cap"]


class SegmentEventStratum2:
    kind: Literal["cell", "fibre"]
    root_id: bytes
    local_root_id: bytes
    global_fibre_id: bytes
    chart_id: Literal["segment-linear-v1"]
    witness_numerator: str
    witness_denominator: str
    root_factor_coefficients: Sequence[str]
    root_ordinal: int
    active_branch_ids: Sequence[bytes]
    events: Sequence[PartitionEvent2]
    left_pair_dispositions: Sequence[BranchPairDisposition2]
    pair_dispositions: Sequence[BranchPairDisposition2]
    right_pair_dispositions: Sequence[BranchPairDisposition2]
    algebraic_root_evaluated: bool
    original_equations_rechecked: bool
    orientation_rechecked: bool
    trim_predicates_rechecked: bool


class EventTraceEvent2:
    def __init__(
        self,
        root_id: bytes,
        global_fibre_id: bytes,
        kind: str,
        feature_ids: Sequence[bytes],
        branch_ids: Sequence[bytes],
        multiplicity: int,
        disposition: str,
        motion_order: int,
    ) -> None: ...
    root_id: bytes
    global_fibre_id: bytes
    kind: str
    feature_ids: tuple[bytes, ...]
    branch_ids: tuple[bytes, ...]
    multiplicity: int
    disposition: str
    motion_order: int
    canonical_bytes: bytes
    canonical_id: bytes


class EventTrace2:
    exact_verdict: Literal["certified", "cap_exceeded", "unresolved"]
    partition: EventPartitionCertificate2
    partition_certificate: EventPartitionCertificate2
    events: Sequence[EventTraceEvent2]
    motion_chart_id: str
    motion_identity: bytes
    effective_cap_bytes: bytes
    whole_rim_disposition: Literal["clear", "material", "partial", "unresolved"]
    oracle_strategy_version: str
    event_cell_count: int
    canonical_bytes: bytes
    canonical_digest: bytes


class SegmentCellStratum2:
    branches: Sequence[SegmentBoundaryBranch2]
    stratum: SegmentEventStratum2


class SegmentEventPartition2:
    source: SegmentEventSource2
    boundary_feature_ids: Sequence[bytes]
    projections: Sequence[ProjectionRecord2]
    branches: Sequence[SegmentBoundaryBranch2]
    strata: Sequence[SegmentEventStratum2]
    certificate: EventPartitionCertificate2
    canonical_bytes: bytes
    canonical_digest: bytes


class VerifiedSegmentEventPartition2:
    verdict: ContinuousTeaVerdict
    partition: SegmentEventPartition2


def extract_boundary_records(stock: Stock2) -> Sequence[BoundaryFeatureRecord2]: ...
def classify_boundary_pair(
    stock: Stock2,
    first_index: int,
    second_index: int,
) -> Sequence[BoundaryEvent2]: ...
def extract_boundary_events(stock: Stock2) -> Sequence[BoundaryEvent2]: ...
def exact_algebraic_backend_evidence() -> AlgebraicBackendEvidence2: ...
def parameter_charts() -> Sequence[ParameterChart2]: ...
def verify_chart_coverage(
    charts: Sequence[ParameterChart2],
) -> VerifiedEventPartition2: ...
def construct_pullback(
    motion_kind: str,
    motion_data: Sequence[str],
    support_kind: str,
    support_data: Sequence[str],
    cutter_radius: str,
    center_chart: str,
    rim_chart: str,
) -> ProjectionRecord2: ...
def projection_from_grid(
    projection_id: str,
    coefficient_rows: Sequence[Sequence[str]],
    degree_bound_id: str,
) -> ProjectionRecord2: ...
def partition_pullback_overlap(
    projection: ProjectionRecord2,
    events: Sequence[PartitionEvent2],
) -> EventPartitionCertificate2: ...
def solve_trimmed_line_branches(
    line_support: Sequence[str],
    trim_start: Sequence[str],
    trim_end: Sequence[str],
    segment_motion: Sequence[str],
    cutter_radius: str,
    rim_chart: str,
) -> Sequence[TrimmedLineBranch2]: ...
def project_regularization_vertex(
    stock: Stock2,
    first_index: int,
    second_index: int,
    vertex_id: bytes,
) -> ProjectedRegularizationVertex2: ...
def classify_ccw_orientation(
    first_numerator: Sequence[str],
    first_denominator: Sequence[str],
    second_numerator: Sequence[str],
    second_denominator: Sequence[str],
    motion_parameter: str,
) -> CcwOrientation2: ...
def partition_cap_crossings(
    first_numerator: Sequence[str],
    first_denominator: Sequence[str],
    second_numerator: Sequence[str],
    second_denominator: Sequence[str],
    cap_numerator: str,
    cap_denominator: str,
    event: PartitionEvent2,
) -> EventPartitionCertificate2: ...
def verify_event_partition(
    certificate: EventPartitionCertificate2,
) -> VerifiedEventPartition2: ...
def order_full_circle_events(
    verified_partition: VerifiedEventPartition2,
    clockwise: bool,
    events: Sequence[EventTraceEvent2],
) -> Sequence[EventTraceEvent2]: ...
def audit_full_circle_tea_event_exact(
    stock: Stock2,
    center_x: float,
    center_y: float,
    phase_dx: float,
    phase_dy: float,
    clockwise: bool,
    tool_radius: float,
    cap_chord_ratio: float,
) -> tuple[Literal["certified", "cap_exceeded", "unresolved"], EventTrace2]: ...
def full_circle_rational_probe_exceeds_cap_exact(
    stock: Stock2,
    center_x: float,
    center_y: float,
    phase_dx: float,
    phase_dy: float,
    chart: int,
    numerator: int,
    denominator: int,
    tool_radius: float,
    cap_chord_ratio: float,
) -> bool: ...
def build_event_trace(
    partition: EventPartitionCertificate2,
    motion_chart_id: str,
    motion_identity: bytes,
    effective_cap_bytes: bytes,
    verdict: ContinuousTeaVerdict,
    whole_rim_disposition: str,
    oracle_strategy_version: str,
    events: Sequence[EventTraceEvent2],
) -> EventTrace2: ...
def mutate_certificate_record(
    certificate: EventPartitionCertificate2,
    mutation: str,
) -> EventPartitionCertificate2: ...
def sha256_bytes(input: bytes) -> bytes: ...
def partition_projections(
    projections: Sequence[ProjectionInput2],
) -> EventPartitionCertificate2: ...
def construct_segment_cell_stratum(
    stock: Stock2,
    source: SegmentEventSource2,
    witness_numerator: str,
    witness_denominator: str,
) -> SegmentCellStratum2: ...
def segment_pair_literal_signs() -> Sequence[str]: ...
def segment_rational_square_root_cases() -> Sequence[str]: ...
def construct_segment_event_partition(
    stock: Stock2,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    tool_radius: float,
    cap_chord_ratio: float,
) -> SegmentEventPartition2: ...
def verify_segment_event_partition(
    stock: Stock2,
    source: SegmentEventSource2,
    candidate: SegmentEventPartition2,
) -> VerifiedSegmentEventPartition2: ...
def mutate_segment_event_partition(
    partition: SegmentEventPartition2,
    mutation: str,
) -> SegmentEventPartition2: ...
