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


class EventPartitionVerificationError(EventSubstrateError): ...


class ContinuousTeaVerdict:
    CERTIFIED: ContinuousTeaVerdict
    CAP_EXCEEDED: ContinuousTeaVerdict
    UNRESOLVED_DEGENERACY: ContinuousTeaVerdict
    @property
    def name(self) -> str: ...


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
    original_equations_rechecked: bool
    orientation_rechecked: bool
    trim_disposition: str


class ProjectionInput2:
    def __init__(
        self,
        projection_id: str,
        coefficients: Sequence[str],
        events: Sequence[PartitionEvent2],
    ) -> None: ...


class AlgebraicRootRecord2:
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


class EventFibre2:
    root_id: bytes
    events: Sequence[PartitionEvent2]


class ProjectionRecord2:
    projection_id: str
    coefficient_rows: Sequence[Sequence[str]]
    factor_coefficients: Sequence[Sequence[str]]
    actual_degree: tuple[int, int]
    degree_bound: tuple[int, int]
    degree_bound_id: str
    normalized_coefficient_bytes: bytes


class OverlapInterval2:
    kind: str
    domain_low: str
    domain_high: str
    witness_numerator: str
    witness_denominator: str
    orientation_disposition: str


class ChartSeam2:
    seam_id: bytes
    owner_chart_id: str


class EventPartitionCertificate2:
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
def partition_projections(
    projections: Sequence[ProjectionInput2],
) -> EventPartitionCertificate2: ...
