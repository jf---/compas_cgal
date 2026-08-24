from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

Float64Array = npt.NDArray[np.float64]

class ExactDepletionConstructionError(RuntimeError): ...
class ExactDepletionCenterLimitError(ExactDepletionConstructionError): ...
class ExactArcDepletionPolicyError(ExactDepletionConstructionError): ...
class ExactArcForgedTraceError(ExactDepletionConstructionError): ...
class NonFiniteExactArcDepletionInputError(ValueError): ...
class InvalidAnnulusRadiiError(ValueError): ...
class NonFiniteAnnulusInputError(ValueError): ...
class NonFiniteCapsuleInputError(ValueError): ...
class LocalDepletionEscapedError(RuntimeError): ...
class CapsuleQuadCertificateError(RuntimeError): ...
class AuditNonFiniteInputError(ValueError): ...
class AuditInvalidPlaneError(ValueError): ...
class AuditUnsupportedGeometryError(ValueError): ...
class AuditOffPlaneError(AuditUnsupportedGeometryError): ...
class AuditContradictoryRoleError(ValueError): ...
class AuditContradictoryOrientationError(ValueError): ...
class AuditSegmentMotion2: ...
class AuditCircleMotion2: ...

class AuditArcMotion2:
    @property
    def digest(self) -> bytes: ...

class AuditVerticalPlunge2: ...
class AuditVerticalRetract2: ...
class AuditClearanceTransport2: ...

AuditLineClassification2 = AuditSegmentMotion2 | AuditVerticalPlunge2 | AuditVerticalRetract2 | AuditClearanceTransport2
AuditCircleClassification2 = AuditCircleMotion2 | AuditClearanceTransport2
AuditArcClassification2 = AuditArcMotion2 | AuditClearanceTransport2

def audit_arc_phase_strategy_version() -> bytes: ...
def classify_audit_line(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    cut_z: float,
    clearance_z: float,
    operation_role: str,
) -> AuditLineClassification2: ...
def classify_audit_circle(
    center: tuple[float, float, float],
    xaxis: tuple[float, float, float],
    yaxis: tuple[float, float, float],
    radius: float,
    clockwise: bool,
    cut_z: float,
    clearance_z: float,
    operation_role: str,
) -> AuditCircleClassification2: ...
def classify_audit_arc(
    center: tuple[float, float, float],
    xaxis: tuple[float, float, float],
    yaxis: tuple[float, float, float],
    radius: float,
    start_angle: float,
    end_angle: float,
    clockwise: bool,
    cut_z: float,
    clearance_z: float,
    operation_role: str,
) -> AuditArcClassification2: ...

class DepletionTrace:
    @property
    def center_count(self) -> int: ...
    @property
    def center_parameters(self) -> Sequence[tuple[int, int, int]]: ...
    def matches_exact_inputs(
        self,
        expected_tool_radius: float,
        expected_max_chord: float,
        expected_center_count_limit: int,
    ) -> bool: ...
    @property
    def strategy_version(self) -> bytes: ...
    cyclic: bool
    exact_incidence: bool
    exact_parameters_in_range: bool
    exact_anchors_present: bool
    exact_removal_radius_valid: bool
    exact_chord_bound_holds: bool
    exact_seam_chord_bound_holds: bool

class ExactArcDepletionTrace2:
    @property
    def center_count(self) -> int: ...
    def matches_exact_inputs(
        self,
        tool_radius: float,
        max_chord: float,
        center_count_limit: int,
    ) -> bool: ...
    def matches_motion(self, motion: AuditArcMotion2) -> bool: ...
    @property
    def canonical_bytes(self) -> bytes: ...
    @property
    def digest(self) -> bytes: ...
    @property
    def strategy_version(self) -> bytes: ...
    @property
    def cyclic(self) -> bool: ...

class Stock2:
    def __init__(self, boundary: Float64Array, holes: Sequence[Float64Array]) -> None: ...
    def contains(self, x: float, y: float) -> bool: ...
    def is_empty(self) -> bool: ...
    def clone(self) -> Stock2: ...
    def is_subset_of(self, other: Stock2) -> bool: ...
    def exactly_equals(self, other: Stock2) -> bool: ...
    def subtract_capsule(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        radius: float,
    ) -> None: ...
    def subtract_capsule_quad(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        radius: float,
    ) -> None: ...
    def subtract_arc_sweep(
        self,
        cx: float,
        cy: float,
        sx: float,
        sy: float,
        ex: float,
        ey: float,
        cw: bool,
        tool_radius: float,
    ) -> None: ...
    def subtract_exact_segment(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        tool_radius: float,
        max_chord: float,
        center_count_limit: int,
    ) -> DepletionTrace: ...
    def subtract_exact_full_circle(
        self,
        cx: float,
        cy: float,
        phase_x: float,
        phase_y: float,
        clockwise: bool,
        tool_radius: float,
        max_chord: float,
        center_count_limit: int,
    ) -> DepletionTrace: ...
    def subtract_exact_arc(
        self,
        motion: AuditArcMotion2,
        tool_radius: float,
        max_chord: float,
        center_count_limit: int,
    ) -> ExactArcDepletionTrace2: ...
    def subtract_disk(self, cx: float, cy: float, radius: float) -> None: ...
    def subtract_annulus(
        self,
        cx: float,
        cy: float,
        inner_radius: float,
        outer_radius: float,
    ) -> None: ...
    def subtract_disk_local(self, cx: float, cy: float, radius: float) -> None: ...
    def subtract_annulus_local(
        self,
        cx: float,
        cy: float,
        inner_radius: float,
        outer_radius: float,
    ) -> None: ...
    def subtract_arc_sweep_local(
        self,
        cx: float,
        cy: float,
        sx: float,
        sy: float,
        ex: float,
        ey: float,
        cw: bool,
        tool_radius: float,
    ) -> None: ...
    def representation_is_valid(self) -> bool: ...
    def arrangement_stats(self) -> tuple[int, int, int]: ...
    def coordinate_digits(self) -> tuple[int, float, int]: ...

def exact_segment_point_is_incident(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    px: float,
    py: float,
) -> bool: ...
def exact_circle_point_is_incident(
    cx: float,
    cy: float,
    phase_x: float,
    phase_y: float,
    px: float,
    py: float,
) -> bool: ...
def exact_segment_structural_density_holds(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    max_chord: float,
    parameters: Sequence[tuple[int, int, int]],
) -> bool: ...
def exact_full_circle_structural_density_holds(
    cx: float,
    cy: float,
    phase_x: float,
    phase_y: float,
    clockwise: bool,
    max_chord: float,
    parameters: Sequence[tuple[int, int, int]],
) -> bool: ...
def exact_segment_undercover_holds(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    exact_length: float,
    tool_radius: float,
    max_chord: float,
    center_count_limit: int,
) -> bool: ...
def exact_full_circle_undercover_holds(
    cx: float,
    cy: float,
    phase_x: float,
    phase_y: float,
    guide_radius: float,
    tool_radius: float,
    max_chord: float,
    center_count_limit: int,
) -> bool: ...
def exact_segment_induction_holds(
    initial: Stock2,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    exact_length: float,
    tool_radius: float,
    max_chord: float,
    center_count_limit: int,
) -> bool: ...
def exact_full_circle_induction_holds(
    initial: Stock2,
    cx: float,
    cy: float,
    phase_x: float,
    phase_y: float,
    guide_radius: float,
    tool_radius: float,
    max_chord: float,
    center_count_limit: int,
) -> bool: ...
def exact_depletion_strategy_version() -> bytes: ...
def exact_entry_depletion_strategy_version() -> bytes: ...
def cap_chord_ratio(cap_radians: float) -> float: ...
def cap_chord_ratio_le(lhs: float, rhs: float) -> bool: ...
def engagement_at(
    stock: Stock2,
    cx: float,
    cy: float,
    tool_radius: float,
    cap_chord_ratio: float,
    gap_close_ratio: float = ...,
) -> tuple[float, float, bool]: ...
def certify_segment_tea(
    stock: Stock2,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    tool_radius: float,
    cap_radians: float,
) -> tuple[float, bool, int]: ...
def _sign_mixed_radical(
    a: float,
    b: float,
    c: float,
    d: float,
    alpha: float,
    beta: float,
) -> int: ...
