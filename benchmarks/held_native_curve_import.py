"""Import prepared Held curves with exact native endpoint incidence.

The native arc factory adjusts centres while preserving authored endpoints.
This import alone establishes neither source-fit bounds nor tangent continuity.
"""

from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_geometry import ReferenceLine
from compas_cgal import _coverage_2 as native
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY


def _xy(point: Point2[WorldXY]) -> tuple[float, float]:
    """Inject authored world-XY millimetres once at the native input boundary."""
    return float(point.x), float(point.y)


def import_held_boundary(case: HeldReferenceCase) -> native.NativeBoundary2:
    """Preserve the supplied directed line/arc chain in native geometry.

    Native construction validates endpoint incidence, closure, and simplicity.
    Native errors propagate; no curve is flattened, dropped, or repaired by
    Python. The polygon projection is not an input to this importer.
    """
    curves: list[native.NativeBoundaryCurve2] = []
    for primitive in case.boundary.primitives:
        if isinstance(primitive, ReferenceLine):
            curves.append(native.NativeBoundaryCurve2.line(_xy(primitive.start), _xy(primitive.end)))
        else:
            curves.append(
                native.NativeBoundaryCurve2.arc(
                    _xy(primitive.start),
                    _xy(primitive.end),
                    _xy(primitive.centre),
                    float(primitive.sweep) > 0.0,
                )
            )
    return native.NativeBoundary2(curves)
