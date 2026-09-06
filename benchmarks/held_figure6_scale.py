"""Calibrate publisher Figure 6 lengths against the visible Figure 5(a) path."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NewType

import numpy as np

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_publisher import Figure5PublisherCubicEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherLineEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherPathEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherPrimitiveEvidence
from benchmarks.held_figure6_publisher import PublisherFigure6Evidence
from benchmarks.held_figure6_publisher import PublisherFigure6Series
from benchmarks.held_figure6_publisher import PublisherGraphicalPathLength
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.units import Degrees
from benchmarks.units import ToolRadiusMultiple
from benchmarks.units import tool_radius_multiple
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import ToolRadius

PublisherGraphicalLengthScale = NewType("PublisherGraphicalLengthScale", float)

FIGURE5_STANDARD_CAP_DEG = Degrees(80.0)

# Deterministic reporting quadrature. Comparing orders 16 and 32 bounds the
# integration drift far below the digitization uncertainty of the publisher plot.
CUBIC_LENGTH_BASE_ORDER = 16
CUBIC_LENGTH_CHECK_ORDER = 32
CUBIC_LENGTH_RELATIVE_CONVERGENCE = 1e-11
CUBIC_LENGTH_ABSOLUTE_CONVERGENCE_MM = 1e-12


class PublisherFigure6ScaleError(BenchmarkError):
    """Publisher Figure 5/6 evidence cannot establish one length scale."""


@dataclass(frozen=True)
class Figure6PublisherLengthCalibration:
    """The published 80-degree standard path expressed in both length frames."""

    figure5_path_length: ToolRadiusMultiple
    figure6_graphical_length: PublisherGraphicalPathLength
    tool_radius_multiples_per_graphical_unit: PublisherGraphicalLengthScale

    def __post_init__(self) -> None:
        values = (
            float(self.figure5_path_length),
            float(self.figure6_graphical_length),
            float(self.tool_radius_multiples_per_graphical_unit),
        )
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise PublisherFigure6ScaleError("Publisher Figure 6 length calibration values must be finite and positive.")

    def to_tool_radius_multiples(self, graphical_length: PublisherGraphicalPathLength) -> ToolRadiusMultiple:
        """Convert one Figure 6 graphical ordinate through the published anchor."""
        value = float(graphical_length)
        if not math.isfinite(value) or value <= 0.0:
            raise PublisherFigure6ScaleError("Publisher Figure 6 graphical length must be finite and positive.")
        return tool_radius_multiple(
            value * float(self.tool_radius_multiples_per_graphical_unit),
            name="publisher Figure 6 dimensionless path length",
        )


def _cubic_length_at_order(primitive: Figure5PublisherCubicEvidence, order: int) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    parameter = (nodes + 1.0) / 2.0
    start = np.array((float(primitive.start.x), float(primitive.start.y)))
    control1 = np.array((float(primitive.control1.x), float(primitive.control1.y)))
    control2 = np.array((float(primitive.control2.x), float(primitive.control2.y)))
    end = np.array((float(primitive.end.x), float(primitive.end.y)))
    derivative = 3.0 * (
        (1.0 - parameter)[:, None] ** 2 * (control1 - start) + 2.0 * ((1.0 - parameter) * parameter)[:, None] * (control2 - control1) + parameter[:, None] ** 2 * (end - control2)
    )
    return float(np.sum(weights * np.linalg.norm(derivative, axis=1)) / 2.0)


def publisher_primitive_length(primitive: Figure5PublisherPrimitiveEvidence) -> Millimetre:
    """Measure one retained Figure 5 connector primitive in normalized world XY."""
    if isinstance(primitive, Figure5PublisherLineEvidence):
        return Millimetre(
            math.dist(
                (float(primitive.start.x), float(primitive.start.y)),
                (float(primitive.end.x), float(primitive.end.y)),
            )
        )
    if not isinstance(primitive, Figure5PublisherCubicEvidence):
        raise PublisherFigure6ScaleError("Publisher Figure 5 path contains an unsupported primitive.")
    base = _cubic_length_at_order(primitive, CUBIC_LENGTH_BASE_ORDER)
    checked = _cubic_length_at_order(primitive, CUBIC_LENGTH_CHECK_ORDER)
    allowed_drift = max(
        CUBIC_LENGTH_ABSOLUTE_CONVERGENCE_MM,
        CUBIC_LENGTH_RELATIVE_CONVERGENCE * checked,
    )
    if not math.isfinite(checked) or checked <= 0.0 or abs(checked - base) > allowed_drift:
        raise PublisherFigure6ScaleError("Publisher Figure 5 connector cubic length did not converge.")
    return Millimetre(checked)


def publisher_figure5_planar_length(path: Figure5PublisherPathEvidence, tool_radius: ToolRadius) -> ToolRadiusMultiple:
    """Measure the complete visible Figure 5(a) stream in tool-radius multiples."""
    if type(path) is not Figure5PublisherPathEvidence or type(tool_radius) is not ToolRadius:
        raise PublisherFigure6ScaleError("Publisher Figure 5 length requires typed path evidence and tool radius.")
    radius = float(tool_radius.value)
    if not math.isfinite(radius) or radius <= 0.0:
        raise PublisherFigure6ScaleError("Publisher Figure 5 tool radius must be finite and positive.")
    turn_length = sum(float(turn.radius) * float(turn.signed_sweep) for turn in path.turns)
    connector_length = sum(float(publisher_primitive_length(primitive)) for turn in path.turns for primitive in turn.connector_after.primitives)
    return tool_radius_multiple(
        (turn_length + connector_length) / radius,
        name="publisher Figure 5 planar path length",
    )


def log_interpolated_graphical_length(series: PublisherFigure6Series, engagement_deg: Degrees) -> PublisherGraphicalPathLength:
    """Interpolate one digitized Figure 6 series consistently with its log-y axis."""
    if type(series) is not PublisherFigure6Series:
        raise PublisherFigure6ScaleError("Publisher Figure 6 interpolation requires one typed series.")
    engagement = float(engagement_deg)
    if not math.isfinite(engagement):
        raise PublisherFigure6ScaleError("Publisher Figure 6 interpolation angle must be finite.")
    for point in series.points:
        if float(point.engagement_deg) == engagement:
            return point.path_length
    for first, second in zip(series.points, series.points[1:], strict=False):
        first_angle = float(first.engagement_deg)
        second_angle = float(second.engagement_deg)
        if first_angle < engagement < second_angle:
            fraction = (engagement - first_angle) / (second_angle - first_angle)
            log_length = (1.0 - fraction) * math.log(float(first.path_length)) + fraction * math.log(float(second.path_length))
            return PublisherGraphicalPathLength(math.exp(log_length))
    raise PublisherFigure6ScaleError("Publisher Figure 6 interpolation angle lies outside the digitized series.")


def calibrate_publisher_figure6_length(
    case: HeldReferenceCase,
    path: Figure5PublisherPathEvidence,
    figure6: PublisherFigure6Evidence,
) -> Figure6PublisherLengthCalibration:
    """Anchor all publisher Figure 6 ordinates to the published Figure 5(a) path."""
    if type(case) is not HeldReferenceCase or case.name != "figure5":
        raise PublisherFigure6ScaleError("Publisher Figure 6 calibration requires the canonical Figure 5 case.")
    if type(figure6) is not PublisherFigure6Evidence:
        raise PublisherFigure6ScaleError("Publisher Figure 6 calibration requires typed publisher evidence.")
    standard = next((series for series in figure6.series if series.label == "standard"), None)
    if standard is None:
        raise PublisherFigure6ScaleError("Publisher Figure 6 evidence has no standard series.")
    figure5_length = publisher_figure5_planar_length(path, case.tool_radius)
    graphical_length = log_interpolated_graphical_length(standard, FIGURE5_STANDARD_CAP_DEG)
    scale = PublisherGraphicalLengthScale(float(figure5_length) / float(graphical_length))
    return Figure6PublisherLengthCalibration(figure5_length, graphical_length, scale)
