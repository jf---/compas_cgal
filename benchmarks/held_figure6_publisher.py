"""Typed graphical observations digitized from Held--Pfeiffer Figure 6."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from typing import NewType
from typing import cast

from benchmarks.errors import BenchmarkError
from benchmarks.units import Degrees

PublisherGraphicalPathLength = NewType("PublisherGraphicalPathLength", float)
Figure6PublisherLabel = Literal["standard", "MATHSM", "contour"]
FIGURE6_PUBLISHER_LABELS: tuple[Figure6PublisherLabel, ...] = ("standard", "MATHSM", "contour")
FIGURE6_PUBLISHER_PROVENANCE: Literal["publisher Figure 6 pixel-digitized graphical observations"] = "publisher Figure 6 pixel-digitized graphical observations"
FIGURE6_PUBLISHER_LENGTH_UNIT: Literal["publisher graphical path-length unit"] = "publisher graphical path-length unit"
DATA_PATH = Path(__file__).with_name("data") / "held_pfeiffer_2025" / "figure6_publisher_graphical.json"


class InvalidFigure6PublisherEvidenceError(BenchmarkError):
    """Publisher Figure 6 graphical evidence is malformed."""


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidFigure6PublisherEvidenceError(f"{name} must be a finite number.")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise InvalidFigure6PublisherEvidenceError(f"{name} must be finite.")
    return numeric


def _integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidFigure6PublisherEvidenceError(f"{name} must be an integer.")
    return value


@dataclass(frozen=True)
class PublisherFigure6Point:
    engagement_deg: Degrees
    path_length: PublisherGraphicalPathLength

    @classmethod
    def build(cls, engagement_deg: object, path_length: object, *, axes: PublisherFigure6Axes) -> PublisherFigure6Point:
        x = _finite(engagement_deg, "publisher engagement")
        y = _finite(path_length, "publisher path length")
        if not float(axes.x_min_deg) <= x <= float(axes.x_max_deg) or not float(axes.y_min) <= y <= float(axes.y_max):
            raise InvalidFigure6PublisherEvidenceError("Publisher point lies outside the recorded axes.")
        return cls(Degrees(x), PublisherGraphicalPathLength(y))


@dataclass(frozen=True)
class PublisherFigure6Axes:
    x_min_deg: Degrees
    x_max_deg: Degrees
    x_tick_deg: Degrees
    y_min: PublisherGraphicalPathLength
    y_max: PublisherGraphicalPathLength
    pixel_left: int
    pixel_right: int
    pixel_top: int
    pixel_bottom: int

    @classmethod
    def build(
        cls,
        *,
        x_min_deg: object,
        x_max_deg: object,
        x_tick_deg: object,
        y_min: object,
        y_max: object,
        pixel_left: object,
        pixel_right: object,
        pixel_top: object,
        pixel_bottom: object,
    ) -> PublisherFigure6Axes:
        x0 = _finite(x_min_deg, "publisher x minimum")
        x1 = _finite(x_max_deg, "publisher x maximum")
        tick = _finite(x_tick_deg, "publisher x tick")
        y0 = _finite(y_min, "publisher y minimum")
        y1 = _finite(y_max, "publisher y maximum")
        left = _integer(pixel_left, "publisher plot left")
        right = _integer(pixel_right, "publisher plot right")
        top = _integer(pixel_top, "publisher plot top")
        bottom = _integer(pixel_bottom, "publisher plot bottom")
        if not (x0 < x1 and tick > 0.0 and y0 > 0.0 and y0 < y1 and left < right and top < bottom):
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 axis bounds are invalid.")
        tick_count = (x1 - x0) / tick
        if not tick_count.is_integer():
            raise InvalidFigure6PublisherEvidenceError("Publisher x tick does not divide its axis range.")
        return cls(Degrees(x0), Degrees(x1), Degrees(tick), PublisherGraphicalPathLength(y0), PublisherGraphicalPathLength(y1), left, right, top, bottom)

    def point_from_pixel(self, pixel_x: object, pixel_y: object) -> PublisherFigure6Point:
        """Map one recorded plot pixel onto the linear-x/log-y publisher axes."""
        x = _integer(pixel_x, "publisher pixel x")
        y = _integer(pixel_y, "publisher pixel y")
        if not self.pixel_left <= x <= self.pixel_right or not self.pixel_top <= y <= self.pixel_bottom:
            raise InvalidFigure6PublisherEvidenceError("Publisher pixel lies outside the recorded plot frame.")
        x_fraction = (x - self.pixel_left) / (self.pixel_right - self.pixel_left)
        y_fraction = (y - self.pixel_top) / (self.pixel_bottom - self.pixel_top)
        engagement = float(self.x_min_deg) + x_fraction * (float(self.x_max_deg) - float(self.x_min_deg))
        log_length = math.log10(float(self.y_max)) + y_fraction * (math.log10(float(self.y_min)) - math.log10(float(self.y_max)))
        return PublisherFigure6Point(Degrees(engagement), PublisherGraphicalPathLength(10.0**log_length))


@dataclass(frozen=True)
class PublisherFigure6Series:
    label: Figure6PublisherLabel
    points: tuple[PublisherFigure6Point, ...]

    @classmethod
    def build(cls, label: object, points: tuple[PublisherFigure6Point, ...]) -> PublisherFigure6Series:
        if label not in FIGURE6_PUBLISHER_LABELS:
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 series has an unknown label.")
        if not points or any(type(point) is not PublisherFigure6Point for point in points):
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 series requires typed points.")
        if any(float(first.engagement_deg) >= float(second.engagement_deg) for first, second in zip(points, points[1:], strict=False)):
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 x values must be strictly increasing and unique.")
        return cls(label, points)


@dataclass(frozen=True)
class PublisherFigure6Evidence:
    provenance: Literal["publisher Figure 6 pixel-digitized graphical observations"]
    length_unit: Literal["publisher graphical path-length unit"]
    axes: PublisherFigure6Axes
    series: tuple[PublisherFigure6Series, ...]

    @classmethod
    def build(cls, axes: PublisherFigure6Axes, series: tuple[PublisherFigure6Series, ...]) -> PublisherFigure6Evidence:
        if type(axes) is not PublisherFigure6Axes or tuple(item.label for item in series) != FIGURE6_PUBLISHER_LABELS:
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 requires the three ordered labeled series.")
        return cls(FIGURE6_PUBLISHER_PROVENANCE, FIGURE6_PUBLISHER_LENGTH_UNIT, axes, series)


def _object(value: object, keys: frozenset[str], name: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != keys:
        raise InvalidFigure6PublisherEvidenceError(f"{name} has an invalid object schema.")
    return cast(dict[str, object], value)


def load_figure6_publisher_evidence(path: Path = DATA_PATH) -> PublisherFigure6Evidence:
    """Load the tracked, self-contained Figure 6 graphical observations."""
    try:
        raw = json.loads(path.read_text(encoding="utf-8"), parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise InvalidFigure6PublisherEvidenceError(f"Cannot load publisher Figure 6 evidence: {exc}") from exc
    root = _object(raw, frozenset(("schema_version", "provenance", "length_unit", "axes", "series")), "publisher Figure 6 evidence")
    if root["schema_version"] != 1 or root["provenance"] != FIGURE6_PUBLISHER_PROVENANCE or root["length_unit"] != FIGURE6_PUBLISHER_LENGTH_UNIT:
        raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 metadata is invalid.")
    axes_raw = _object(
        root["axes"], frozenset(("x_min_deg", "x_max_deg", "x_tick_deg", "y_min", "y_max", "pixel_left", "pixel_right", "pixel_top", "pixel_bottom")), "publisher Figure 6 axes"
    )
    axes = PublisherFigure6Axes.build(**axes_raw)
    if not isinstance(root["series"], list):
        raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 series must be an array.")
    built: list[PublisherFigure6Series] = []
    for raw_series in root["series"]:
        item = _object(raw_series, frozenset(("label", "points")), "publisher Figure 6 series")
        if not isinstance(item["points"], list):
            raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 points must be an array.")
        points: list[PublisherFigure6Point] = []
        for raw_point in item["points"]:
            if not isinstance(raw_point, list) or len(raw_point) != 2:
                raise InvalidFigure6PublisherEvidenceError("Publisher Figure 6 point must be an x/y pair.")
            points.append(PublisherFigure6Point.build(raw_point[0], raw_point[1], axes=axes))
        built.append(PublisherFigure6Series.build(item["label"], tuple(points)))
    return PublisherFigure6Evidence.build(axes, tuple(built))
