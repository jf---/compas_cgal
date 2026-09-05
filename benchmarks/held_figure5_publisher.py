"""Decode and load self-contained Figure 5(a) publisher evidence."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Literal
from typing import NewType
from typing import cast
from xml.etree import ElementTree

import numpy as np
import numpy.typing as npt
from fontTools.pens.boundsPen import BoundsPen  # type: ignore[import-untyped]
from fontTools.pens.recordingPen import RecordingPen  # type: ignore[import-untyped]
from fontTools.svgLib.path import parse_path  # type: ignore[import-untyped]

from benchmarks.errors import BenchmarkError
from benchmarks.held_reference_cases import HeldReferenceCase
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Radian
from compas_cgal.adaptive.units import WorldXY

PublisherTurnOrdinal = NewType("PublisherTurnOrdinal", int)
PublisherOperatorOrdinal = NewType("PublisherOperatorOrdinal", int)
FIGURE5_PUBLISHER_PATH_PROVENANCE: Literal["publisher Figure 5(a) shape-only turn evidence"] = "publisher Figure 5(a) shape-only turn evidence"
FIGURE5_PUBLISHER_FRAME: Literal["world_xy"] = "world_xy"
FIGURE5_PUBLISHER_SCHEMA_VERSION = 2
FIGURE5_PUBLISHER_TURN_COUNT = 265
FIGURE5_PUBLISHER_CUBIC_COUNT = 1218
FIGURE5_PUBLISHER_LINE_COUNT = 211
FIGURE5_PUBLISHER_PURPLE = "rgb(62.69989%, 12.5%, 94.099426%)"
FIGURE5_PUBLISHER_BOUNDARY_GREEN = "rgb(17.999268%, 54.499817%, 34.098816%)"
FIGURE5_PANEL_A_TRANSFORM = "matrix(0.476643, 0, 0, -0.476643, 40.067609, 229.645005)"
FIGURE5_PANEL_A_SCALE = 0.476643
FIGURE5_PANEL_A_TRANSLATE_X = 40.067609
FIGURE5_PANEL_A_TRANSLATE_Y = 229.645005
# Decoder classification bounds, not machining tolerances.
PUBLISHER_TURN_RELATIVE_RADIAL_RMS_BOUND = 5e-4
PUBLISHER_TURN_MINIMUM_SWEEP_RAD = 5.5
PUBLISHER_CLOSED_TURN_GAP_RATIO_BOUND = 0.01
PUBLISHER_TURN_SAMPLES_PER_CUBIC = 17
DATA_PATH = Path(__file__).with_name("data") / "held_pfeiffer_2025" / "figure5_publisher_turns.json"


class InvalidFigure5PublisherEvidenceError(BenchmarkError):
    """Figure 5 publisher evidence is malformed or ambiguous."""


@dataclass(frozen=True)
class Figure5PublisherLineEvidence:
    operator_ordinal: PublisherOperatorOrdinal
    start: Point2[WorldXY]
    end: Point2[WorldXY]

    @classmethod
    def build(cls, ordinal: object, start: Point2[WorldXY], end: Point2[WorldXY]) -> Figure5PublisherLineEvidence:
        if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal < 0 or type(start) is not Point2 or type(end) is not Point2:
            raise InvalidFigure5PublisherEvidenceError("Publisher line requires a non-negative integer ordinal and typed endpoints.")
        return cls(PublisherOperatorOrdinal(ordinal), start, end)


@dataclass(frozen=True)
class Figure5PublisherCubicEvidence:
    operator_ordinal: PublisherOperatorOrdinal
    start: Point2[WorldXY]
    control1: Point2[WorldXY]
    control2: Point2[WorldXY]
    end: Point2[WorldXY]

    @classmethod
    def build(cls, ordinal: object, start: Point2[WorldXY], control1: Point2[WorldXY], control2: Point2[WorldXY], end: Point2[WorldXY]) -> Figure5PublisherCubicEvidence:
        if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal < 0 or any(type(point) is not Point2 for point in (start, control1, control2, end)):
            raise InvalidFigure5PublisherEvidenceError("Publisher cubic requires a non-negative integer ordinal and typed points.")
        return cls(PublisherOperatorOrdinal(ordinal), start, control1, control2, end)


Figure5PublisherPrimitiveEvidence = Figure5PublisherLineEvidence | Figure5PublisherCubicEvidence


@dataclass(frozen=True)
class Figure5PublisherConnectorEvidence:
    """Ordered publisher primitives following one circular turn."""

    primitives: tuple[Figure5PublisherPrimitiveEvidence, ...]

    @classmethod
    def build(cls, primitives: tuple[Figure5PublisherPrimitiveEvidence, ...]) -> Figure5PublisherConnectorEvidence:
        if any(type(item) not in (Figure5PublisherLineEvidence, Figure5PublisherCubicEvidence) for item in primitives):
            raise InvalidFigure5PublisherEvidenceError("Publisher connector contains a foreign primitive type.")
        if any(first.end != second.start for first, second in zip(primitives, primitives[1:], strict=False)):
            raise InvalidFigure5PublisherEvidenceError("Publisher connector primitives are discontinuous.")
        return cls(primitives)

    @property
    def cubic_count(self) -> int:
        return sum(isinstance(item, Figure5PublisherCubicEvidence) for item in self.primitives)

    @property
    def line_count(self) -> int:
        return sum(isinstance(item, Figure5PublisherLineEvidence) for item in self.primitives)


@dataclass(frozen=True)
class Figure5PublisherTurnEvidence:
    ordinal: PublisherTurnOrdinal
    center: Point2[WorldXY]
    radius: Millimetre
    orientation: Literal["CCW"]
    entry: Point2[WorldXY]
    exit: Point2[WorldXY]
    closure_gap: Millimetre
    closure_class: Literal["closed", "open"]
    signed_sweep: Radian
    relative_radial_rms: float
    operator_start: PublisherOperatorOrdinal
    operator_stop: PublisherOperatorOrdinal
    connector_after: Figure5PublisherConnectorEvidence

    @classmethod
    def build(
        cls,
        *,
        ordinal: object,
        center: Point2[WorldXY],
        radius: object,
        entry: Point2[WorldXY],
        exit: Point2[WorldXY],
        closure_gap: object,
        closure_class: object,
        signed_sweep: object,
        relative_radial_rms: object,
        operator_start: object,
        operator_stop: object,
        connector_after: Figure5PublisherConnectorEvidence,
    ) -> Figure5PublisherTurnEvidence:
        scalars = (radius, closure_gap, signed_sweep, relative_radial_rms)
        if (
            isinstance(ordinal, bool)
            or not isinstance(ordinal, int)
            or ordinal < 0
            or any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) for value in scalars)
        ):
            raise InvalidFigure5PublisherEvidenceError("Publisher turn requires finite scalar evidence and a non-negative integer ordinal.")
        radius_value = cast(int | float, radius)
        gap_value = cast(int | float, closure_gap)
        sweep_value = cast(int | float, signed_sweep)
        rms_value = cast(int | float, relative_radial_rms)
        if float(radius_value) <= 0 or float(gap_value) < 0 or float(sweep_value) <= 0 or float(rms_value) < 0:
            raise InvalidFigure5PublisherEvidenceError("Publisher turn radius/sweep must be positive and gap/RMS non-negative.")
        if closure_class not in ("closed", "open") or any(type(point) is not Point2 for point in (center, entry, exit)):
            raise InvalidFigure5PublisherEvidenceError("Publisher turn requires typed world geometry and closure class.")
        if any(isinstance(value, bool) or not isinstance(value, int) for value in (operator_start, operator_stop)):
            raise InvalidFigure5PublisherEvidenceError("Publisher turn operator range requires integer ordinals.")
        start_value = cast(int, operator_start)
        stop_value = cast(int, operator_stop)
        return cls(
            PublisherTurnOrdinal(ordinal),
            center,
            Millimetre(float(radius_value)),
            "CCW",
            entry,
            exit,
            Millimetre(float(gap_value)),
            closure_class,
            Radian(float(sweep_value)),
            float(rms_value),
            PublisherOperatorOrdinal(start_value),
            PublisherOperatorOrdinal(stop_value),
            connector_after,
        )

    @property
    def closed(self) -> bool:
        return self.closure_class == "closed"


@dataclass(frozen=True)
class Figure5PublisherPathEvidence:
    provenance: Literal["publisher Figure 5(a) shape-only turn evidence"]
    frame: Literal["world_xy"]
    stream_start: Point2[WorldXY]
    stream_end: Point2[WorldXY]
    turns: tuple[Figure5PublisherTurnEvidence, ...]

    @classmethod
    def build(cls, stream_start: Point2[WorldXY], stream_end: Point2[WorldXY], turns: tuple[Figure5PublisherTurnEvidence, ...]) -> Figure5PublisherPathEvidence:
        if len(turns) != FIGURE5_PUBLISHER_TURN_COUNT or tuple(map(lambda turn: int(turn.ordinal), turns)) != tuple(range(FIGURE5_PUBLISHER_TURN_COUNT)):
            raise InvalidFigure5PublisherEvidenceError("Publisher Figure 5 evidence must contain 265 ordered turns.")
        if sum(turn.connector_after.cubic_count for turn in turns) != FIGURE5_PUBLISHER_CUBIC_COUNT - 4 * FIGURE5_PUBLISHER_TURN_COUNT:
            raise InvalidFigure5PublisherEvidenceError("Publisher connector cubic coverage is incomplete.")
        if sum(turn.connector_after.line_count for turn in turns) != FIGURE5_PUBLISHER_LINE_COUNT:
            raise InvalidFigure5PublisherEvidenceError("Publisher connector line coverage is incomplete.")
        if tuple(int(turn.ordinal) for turn in turns if not turn.closed) != tuple(range(123, 136)):
            raise InvalidFigure5PublisherEvidenceError("Publisher open-turn family differs from the observed stream.")
        cursor = 0
        for index, turn in enumerate(turns):
            if (
                float(turn.radius) <= 0
                or float(turn.closure_gap) < 0
                or float(turn.signed_sweep) <= 0
                or float(turn.relative_radial_rms) < 0
                or not all(math.isfinite(value) for value in (float(turn.radius), float(turn.closure_gap), float(turn.signed_sweep), float(turn.relative_radial_rms)))
            ):
                raise InvalidFigure5PublisherEvidenceError("Publisher turn scalar evidence is invalid.")
            if int(turn.operator_start) != cursor or int(turn.operator_stop) != cursor + 4:
                raise InvalidFigure5PublisherEvidenceError("Publisher turn operator ranges are not contiguous four-cubic blocks.")
            cursor = int(turn.operator_stop)
            for primitive in turn.connector_after.primitives:
                if int(primitive.operator_ordinal) != cursor:
                    raise InvalidFigure5PublisherEvidenceError("Publisher connector operator coverage is not contiguous.")
                cursor += 1
            connector = turn.connector_after.primitives
            if connector and turn.exit != connector[0].start:
                raise InvalidFigure5PublisherEvidenceError("Publisher turn exit and connector start are discontinuous.")
            outgoing = connector[-1].end if connector else turn.exit
            if index + 1 < len(turns) and outgoing != turns[index + 1].entry:
                raise InvalidFigure5PublisherEvidenceError("Publisher connector and next turn entry are discontinuous.")
        if cursor != FIGURE5_PUBLISHER_CUBIC_COUNT + FIGURE5_PUBLISHER_LINE_COUNT:
            raise InvalidFigure5PublisherEvidenceError("Publisher operator coverage is incomplete.")
        final_end = turns[-1].connector_after.primitives[-1].end if turns[-1].connector_after.primitives else turns[-1].exit
        if turns[0].entry != stream_start or final_end != stream_end:
            raise InvalidFigure5PublisherEvidenceError("Publisher stream endpoints do not match its primitives.")
        return cls(FIGURE5_PUBLISHER_PATH_PROVENANCE, FIGURE5_PUBLISHER_FRAME, stream_start, stream_end, turns)


BezierSegment = tuple[tuple[float, float], tuple[float, float], tuple[float, float], tuple[float, float]]


@dataclass(frozen=True)
class _SourceLine:
    ordinal: int
    start: tuple[float, float]
    end: tuple[float, float]


@dataclass(frozen=True)
class _SourceCubic:
    ordinal: int
    segment: BezierSegment

    @property
    def start(self) -> tuple[float, float]:
        return self.segment[0]

    @property
    def end(self) -> tuple[float, float]:
        return self.segment[3]


_SourcePrimitive = _SourceLine | _SourceCubic


def _bezier_point(segment: BezierSegment, parameter: float) -> npt.NDArray[np.float64]:
    points = np.asarray(segment, dtype=float)
    result = (1 - parameter) ** 3 * points[0] + 3 * (1 - parameter) ** 2 * parameter * points[1] + 3 * (1 - parameter) * parameter**2 * points[2] + parameter**3 * points[3]
    return cast(npt.NDArray[np.float64], result)


def _circle_fit(segments: tuple[BezierSegment, ...]) -> tuple[float, float, tuple[float, float], float, float]:
    samples = np.asarray(
        [_bezier_point(segment, parameter) for segment in segments for parameter in np.linspace(0, 1, PUBLISHER_TURN_SAMPLES_PER_CUBIC, endpoint=False)]
        + [np.asarray(segments[-1][-1])]
    )
    coefficients = np.column_stack((2 * samples[:, 0], 2 * samples[:, 1], np.ones(len(samples))))
    center_x, center_y, constant = np.linalg.lstsq(coefficients, np.sum(samples * samples, axis=1), rcond=None)[0]
    radius = math.sqrt(float(constant + center_x * center_x + center_y * center_y))
    distances = np.hypot(samples[:, 0] - center_x, samples[:, 1] - center_y)
    rms = float(np.sqrt(np.mean((distances - radius) ** 2)) / radius)
    angles = np.unwrap(np.arctan2(samples[:, 1] - center_y, samples[:, 0] - center_x))
    sweep = float(angles[-1] - angles[0])
    gap_ratio = math.dist(segments[0][0], segments[-1][-1]) / radius
    return rms, sweep, (float(center_x), float(center_y)), radius, gap_ratio


def _panel_elements(root: ElementTree.Element, stroke: str) -> tuple[ElementTree.Element, ...]:
    return tuple(element for element in root.iter("{http://www.w3.org/2000/svg}path") if element.get("stroke") == stroke and element.get("transform") == FIGURE5_PANEL_A_TRANSFORM)


def _source_origin(root: ElementTree.Element) -> tuple[float, float]:
    boundary = _panel_elements(root, FIGURE5_PUBLISHER_BOUNDARY_GREEN)
    if len(boundary) != 12:
        raise InvalidFigure5PublisherEvidenceError("Publisher Figure 5(a) requires twelve registered boundary paths.")
    bounds = []
    for element in boundary:
        pen = BoundsPen(None)
        parse_path(element.get("d", ""), pen)
        if pen.bounds is None:
            raise InvalidFigure5PublisherEvidenceError("Publisher boundary path has no geometry.")
        bounds.append(pen.bounds)
    return FIGURE5_PANEL_A_SCALE * min(bound[0] for bound in bounds) + FIGURE5_PANEL_A_TRANSLATE_X, -FIGURE5_PANEL_A_SCALE * min(
        bound[1] for bound in bounds
    ) + FIGURE5_PANEL_A_TRANSLATE_Y


def extract_figure5_publisher_path(svg_path: Path, case: HeldReferenceCase) -> Figure5PublisherPathEvidence:
    """Decode the publisher stream without consulting repository candidates."""
    if not isinstance(svg_path, Path) or type(case) is not HeldReferenceCase or case.name != "figure5":
        raise InvalidFigure5PublisherEvidenceError("Publisher extraction requires the canonical Figure 5 case.")
    root = ElementTree.parse(svg_path).getroot()
    purple = _panel_elements(root, FIGURE5_PUBLISHER_PURPLE)
    if len(purple) != 1:
        raise InvalidFigure5PublisherEvidenceError("Publisher Figure 5(a) requires one purple stream.")
    origin_x, origin_y = _source_origin(root)
    world_scale = 1 / float(case.source_tool_radius)

    def world(point: tuple[float, float]) -> Point2[WorldXY]:
        page_x = FIGURE5_PANEL_A_SCALE * point[0] + FIGURE5_PANEL_A_TRANSLATE_X
        page_y = -FIGURE5_PANEL_A_SCALE * point[1] + FIGURE5_PANEL_A_TRANSLATE_Y
        return Point2[WorldXY].build((page_x - origin_x) * world_scale, (origin_y - page_y) * world_scale)

    pen = RecordingPen()
    parse_path(purple[0].get("d", ""), pen)
    primitives: list[_SourcePrimitive] = []
    current: tuple[float, float] | None = None
    stream_start: tuple[float, float] | None = None
    for operator, points in pen.value:
        if operator == "moveTo":
            if stream_start is not None:
                raise InvalidFigure5PublisherEvidenceError("Publisher Figure 5 must be one continuous subpath.")
            current = points[0]
            stream_start = current
        elif operator == "curveTo":
            if current is None:
                raise InvalidFigure5PublisherEvidenceError("Publisher cubic precedes its move.")
            primitives.append(_SourceCubic(len(primitives), (current, points[0], points[1], points[2])))
            current = points[2]
        elif operator == "lineTo":
            if current is None:
                raise InvalidFigure5PublisherEvidenceError("Publisher line precedes its move.")
            primitives.append(_SourceLine(len(primitives), current, points[0]))
            current = points[0]
        elif operator != "endPath":
            raise InvalidFigure5PublisherEvidenceError(f"Unsupported publisher path operator {operator!r}.")
    if stream_start is None or current is None:
        raise InvalidFigure5PublisherEvidenceError("Publisher stream is empty.")

    windows: list[tuple[int, tuple[float, float], float, float, float, float]] = []
    for start in range(len(primitives) - 3):
        block = primitives[start : start + 4]
        if not all(isinstance(item, _SourceCubic) for item in block):
            continue
        rms, sweep, center, radius, gap_ratio = _circle_fit(tuple(cast(_SourceCubic, item).segment for item in block))
        if rms < PUBLISHER_TURN_RELATIVE_RADIAL_RMS_BOUND and abs(sweep) > PUBLISHER_TURN_MINIMUM_SWEEP_RAD:
            if sweep <= 0:
                raise InvalidFigure5PublisherEvidenceError("Publisher turn is not CCW in reconstructed world XY.")
            windows.append((start, center, radius, gap_ratio, rms, sweep))
    if len(windows) != FIGURE5_PUBLISHER_TURN_COUNT or any(second[0] - first[0] < 4 for first, second in zip(windows, windows[1:])):
        raise InvalidFigure5PublisherEvidenceError("Publisher circular turns are incomplete or overlap.")

    def convert(item: _SourcePrimitive) -> Figure5PublisherPrimitiveEvidence:
        if isinstance(item, _SourceLine):
            return Figure5PublisherLineEvidence.build(item.ordinal, world(item.start), world(item.end))
        return Figure5PublisherCubicEvidence.build(item.ordinal, world(item.segment[0]), world(item.segment[1]), world(item.segment[2]), world(item.segment[3]))

    scale = FIGURE5_PANEL_A_SCALE * world_scale
    turns = []
    for ordinal, (start, center, radius, gap_ratio, rms, sweep) in enumerate(windows):
        following = windows[ordinal + 1][0] if ordinal + 1 < len(windows) else len(primitives)
        circle_block = tuple(cast(_SourceCubic, item) for item in primitives[start : start + 4])
        turns.append(
            Figure5PublisherTurnEvidence.build(
                ordinal=ordinal,
                center=world(center),
                radius=radius * scale,
                entry=world(circle_block[0].start),
                exit=world(circle_block[-1].end),
                closure_gap=gap_ratio * radius * scale,
                closure_class="closed" if gap_ratio < PUBLISHER_CLOSED_TURN_GAP_RATIO_BOUND else "open",
                signed_sweep=sweep,
                relative_radial_rms=rms,
                operator_start=start,
                operator_stop=start + 4,
                connector_after=Figure5PublisherConnectorEvidence.build(tuple(convert(item) for item in primitives[start + 4 : following])),
            )
        )
    return Figure5PublisherPathEvidence.build(world(stream_start), world(current), tuple(turns))


def _point_payload(point: Point2[WorldXY]) -> list[float]:
    return [float(point.x), float(point.y)]


def _primitive_payload(item: Figure5PublisherPrimitiveEvidence) -> dict[str, object]:
    common: dict[str, object] = {"operator_ordinal": int(item.operator_ordinal), "start_mm": _point_payload(item.start), "end_mm": _point_payload(item.end)}
    if isinstance(item, Figure5PublisherLineEvidence):
        return {"kind": "line", **common}
    return {"kind": "cubic", **common, "control1_mm": _point_payload(item.control1), "control2_mm": _point_payload(item.control2)}


def _payload(evidence: Figure5PublisherPathEvidence) -> dict[str, object]:
    return {
        "schema_version": FIGURE5_PUBLISHER_SCHEMA_VERSION,
        "provenance": evidence.provenance,
        "frame": evidence.frame,
        "stream_start_mm": _point_payload(evidence.stream_start),
        "stream_end_mm": _point_payload(evidence.stream_end),
        "turns": [
            {
                "ordinal": int(turn.ordinal),
                "center_mm": _point_payload(turn.center),
                "radius_mm": float(turn.radius),
                "orientation": turn.orientation,
                "entry_mm": _point_payload(turn.entry),
                "exit_mm": _point_payload(turn.exit),
                "closure_gap_mm": float(turn.closure_gap),
                "closure_class": turn.closure_class,
                "signed_sweep_rad": float(turn.signed_sweep),
                "relative_radial_rms": turn.relative_radial_rms,
                "operator_range": [int(turn.operator_start), int(turn.operator_stop)],
                "connector_after": [_primitive_payload(item) for item in turn.connector_after.primitives],
            }
            for turn in evidence.turns
        ],
    }


def write_figure5_publisher_path(evidence: Figure5PublisherPathEvidence, destination: Path) -> None:
    if type(evidence) is not Figure5PublisherPathEvidence or not isinstance(destination, Path):
        raise InvalidFigure5PublisherEvidenceError("Publisher evidence writing requires typed evidence and a path.")
    destination.write_text(json.dumps(_payload(evidence), separators=(",", ":")) + "\n", encoding="utf-8")


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise InvalidFigure5PublisherEvidenceError(f"Publisher evidence contains duplicate key {key!r}.")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise InvalidFigure5PublisherEvidenceError(f"Publisher evidence contains non-finite value {value}.")


def _point(raw: Any) -> Point2[WorldXY]:
    if not isinstance(raw, list) or len(raw) != 2:
        raise InvalidFigure5PublisherEvidenceError("Publisher point must contain two coordinates.")
    return Point2[WorldXY].build(raw)


def _load_primitive(raw: Any) -> Figure5PublisherPrimitiveEvidence:
    if not isinstance(raw, dict):
        raise InvalidFigure5PublisherEvidenceError("Publisher connector primitive must be an object.")
    common = {"kind", "operator_ordinal", "start_mm", "end_mm"}
    if raw.get("kind") == "line" and set(raw) == common:
        return Figure5PublisherLineEvidence.build(raw["operator_ordinal"], _point(raw["start_mm"]), _point(raw["end_mm"]))
    if raw.get("kind") == "cubic" and set(raw) == common | {"control1_mm", "control2_mm"}:
        return Figure5PublisherCubicEvidence.build(raw["operator_ordinal"], _point(raw["start_mm"]), _point(raw["control1_mm"]), _point(raw["control2_mm"]), _point(raw["end_mm"]))
    raise InvalidFigure5PublisherEvidenceError("Publisher connector primitive has unexpected fields.")


def load_figure5_publisher_path(path: Path = DATA_PATH) -> Figure5PublisherPathEvidence:
    if not isinstance(path, Path):
        raise InvalidFigure5PublisherEvidenceError("Publisher evidence loading requires a path.")
    payload = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs, parse_constant=_reject_constant)
    expected = {"schema_version", "provenance", "frame", "stream_start_mm", "stream_end_mm", "turns"}
    if not isinstance(payload, dict) or set(payload) != expected:
        raise InvalidFigure5PublisherEvidenceError("Publisher evidence has unexpected top-level fields.")
    if payload["schema_version"] != FIGURE5_PUBLISHER_SCHEMA_VERSION or payload["provenance"] != FIGURE5_PUBLISHER_PATH_PROVENANCE or payload["frame"] != FIGURE5_PUBLISHER_FRAME:
        raise InvalidFigure5PublisherEvidenceError("Publisher evidence metadata is not canonical Figure 5(a).")
    turns = []
    expected_turn = {
        "ordinal",
        "center_mm",
        "radius_mm",
        "orientation",
        "entry_mm",
        "exit_mm",
        "closure_gap_mm",
        "closure_class",
        "signed_sweep_rad",
        "relative_radial_rms",
        "operator_range",
        "connector_after",
    }
    for raw in payload["turns"]:
        if not isinstance(raw, dict) or set(raw) != expected_turn or raw["orientation"] != "CCW" or raw["closure_class"] not in ("closed", "open"):
            raise InvalidFigure5PublisherEvidenceError("Publisher turn evidence has unexpected fields.")
        operator_range = raw["operator_range"]
        if not isinstance(operator_range, list) or len(operator_range) != 2:
            raise InvalidFigure5PublisherEvidenceError("Publisher operator range must contain two ordinals.")
        turns.append(
            Figure5PublisherTurnEvidence.build(
                ordinal=raw["ordinal"],
                center=_point(raw["center_mm"]),
                radius=raw["radius_mm"],
                entry=_point(raw["entry_mm"]),
                exit=_point(raw["exit_mm"]),
                closure_gap=raw["closure_gap_mm"],
                closure_class=raw["closure_class"],
                signed_sweep=raw["signed_sweep_rad"],
                relative_radial_rms=raw["relative_radial_rms"],
                operator_start=operator_range[0],
                operator_stop=operator_range[1],
                connector_after=Figure5PublisherConnectorEvidence.build(tuple(_load_primitive(item) for item in raw["connector_after"])),
            )
        )
    return Figure5PublisherPathEvidence.build(_point(payload["stream_start_mm"]), _point(payload["stream_end_mm"]), tuple(turns))
