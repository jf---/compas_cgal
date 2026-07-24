from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

Float64Array = npt.NDArray[np.float64]


class Stock2:
    def __init__(self, boundary: Float64Array, holes: Sequence[Float64Array]) -> None: ...
    def contains(self, x: float, y: float) -> bool: ...
    def is_empty(self) -> bool: ...
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
    def subtract_disk(self, cx: float, cy: float, radius: float) -> None: ...


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
