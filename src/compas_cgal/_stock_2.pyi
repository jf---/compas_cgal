from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

Float64Array = npt.NDArray[np.float64]


class ExactDepletionConstructionError(RuntimeError): ...
class ExactDepletionCenterLimitError(ExactDepletionConstructionError): ...


class DepletionTrace:
    @property
    def center_count(self) -> int: ...
    @property
    def center_parameters(self) -> Sequence[tuple[int, int, int]]: ...
    @property
    def max_chord(self) -> float: ...
    @property
    def removal_radius(self) -> float: ...
    @property
    def strategy_version(self) -> bytes: ...
    cyclic: bool
    exact_incidence: bool
    exact_parameters_in_range: bool
    exact_anchors_present: bool
    exact_removal_radius_valid: bool
    exact_chord_bound_holds: bool
    exact_seam_chord_bound_holds: bool


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
    def subtract_disk(self, cx: float, cy: float, radius: float) -> None: ...


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
