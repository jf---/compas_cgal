import math
import struct
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Final
from typing import Self
from typing import TypeVar

from compas_cgal.adaptive.errors import CanonicalEncodingError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckEffectiveCap
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

CANONICAL_ENCODING_VERSION: Final[bytes] = b"CCAN\x00\x01"
_U64: Final[struct.Struct] = struct.Struct(">Q")
_BINARY64: Final[struct.Struct] = struct.Struct(">d")
ExactT = TypeVar("ExactT")

_PASSAGE_TAGS: Final[dict[PassageState, bytes]] = {
    PassageState.UNVISITED: b"passage-unvisited-v1",
    PassageState.FIRST_PASS_COMPLETE: b"passage-first-complete-v1",
    PassageState.SECOND_PASS_COMPLETE: b"passage-second-complete-v1",
    PassageState.TERMINAL: b"passage-terminal-v1",
}
_CUT_INTENT_TAGS: Final[dict[CutIntent, bytes]] = {
    CutIntent.CLIMB: b"cut-intent-climb-v1",
    CutIntent.CONVENTIONAL: b"cut-intent-conventional-v1",
}
_MATERIAL_SIDE_TAGS: Final[dict[MaterialSide, bytes]] = {
    MaterialSide.INSIDE: b"material-side-inside-v1",
    MaterialSide.OUTSIDE: b"material-side-outside-v1",
}
_CIRCLE_ORIENTATION_TAGS: Final[dict[CircleOrientation, bytes]] = {
    CircleOrientation.CLOCKWISE: b"circle-clockwise-v1",
    CircleOrientation.COUNTERCLOCKWISE: b"circle-counterclockwise-v1",
}


def _require_exact(value: object, expected_type: type[ExactT], name: str) -> ExactT:
    if type(value) is not expected_type:
        raise CanonicalEncodingError(f"{name} must be exact {expected_type.__name__}, not a subclass or alternate type.")
    return value


def _node(kind: bytes, payload: bytes) -> bytes:
    if type(kind) is not bytes or len(kind) != 1:
        raise CanonicalEncodingError("canonical node kind must be exactly one byte.")
    if type(payload) is not bytes:
        raise CanonicalEncodingError("canonical node payload must be bytes.")
    return CANONICAL_ENCODING_VERSION + kind + _U64.pack(len(payload)) + payload


def encode_bytes(value: bytes) -> bytes:
    if type(value) is not bytes:
        raise CanonicalEncodingError("canonical byte payload must be bytes.")
    return _node(b"B", value)


def encode_integer(value: int) -> bytes:
    if type(value) is not int:
        raise CanonicalEncodingError("canonical integer must be an integer, not bool or another numeric type.")
    magnitude = abs(value)
    byte_count = (magnitude.bit_length() + 7) // 8
    magnitude_bytes = magnitude.to_bytes(byte_count, byteorder="big")
    sign = b"\x01" if value < 0 else b"\x00"
    return _node(b"I", sign + magnitude_bytes)


def encode_binary64(value: float) -> bytes:
    if type(value) is not float:
        raise CanonicalEncodingError("canonical binary64 value must be an exact float.")
    if not math.isfinite(value):
        raise CanonicalEncodingError("canonical binary64 value must be finite.")
    if value == 0.0:
        value = 0.0
    return _node(b"D", _BINARY64.pack(value))


def encode_boolean(value: bool) -> bytes:
    if type(value) is not bool:
        raise CanonicalEncodingError("canonical boolean must be explicitly true or false.")
    return _node(b"?", b"\x01" if value else b"\x00")


def encode_sequence(values: Sequence[bytes]) -> bytes:
    try:
        items = tuple(values)
    except TypeError:
        raise CanonicalEncodingError("canonical sequence must be finite.") from None
    payload = bytearray(_U64.pack(len(items)))
    for item in items:
        encoded_item = encode_bytes(item)
        payload.extend(encoded_item)
    return _node(b"S", bytes(payload))


def encode_tagged_union(tag: bytes, payload: bytes) -> bytes:
    if type(tag) is not bytes or not tag:
        raise CanonicalEncodingError("canonical union tag must be nonempty bytes.")
    return _node(b"T", encode_bytes(tag) + encode_bytes(payload))


def encode_component_map(components: Mapping[bytes, bytes]) -> bytes:
    entries: list[tuple[bytes, bytes]] = []
    canonical_keys: set[bytes] = set()
    for key, value in components.items():
        if type(key) is not bytes or not key:
            raise CanonicalEncodingError("canonical component names must be nonempty bytes.")
        if type(value) is not bytes:
            raise CanonicalEncodingError("canonical component values must be bytes.")
        canonical_key = encode_bytes(key)
        if canonical_key in canonical_keys:
            raise CanonicalEncodingError("canonical component names must be unique.")
        canonical_keys.add(canonical_key)
        entries.append((canonical_key, encode_bytes(value)))
    entries.sort(key=lambda entry: entry[0])
    payload = _U64.pack(len(entries)) + b"".join(key + value for key, value in entries)
    return _node(b"M", payload)


def _integer_payload_value(payload: bytes) -> int:
    if not payload or payload[0] not in (0, 1):
        raise CanonicalEncodingError("canonical integer requires an explicit positive or negative sign byte.")
    magnitude_bytes = payload[1:]
    if magnitude_bytes.startswith(b"\x00"):
        raise CanonicalEncodingError("canonical integer magnitude must use minimal big-endian bytes.")
    magnitude = int.from_bytes(magnitude_bytes, byteorder="big")
    if payload[0] == 1 and magnitude == 0:
        raise CanonicalEncodingError("canonical integer cannot encode negative zero.")
    return -magnitude if payload[0] == 1 else magnitude


def _parse_node(value: bytes, offset: int) -> tuple[bytes, bytes, int]:
    prefix_end = offset + len(CANONICAL_ENCODING_VERSION)
    kind_end = prefix_end + 1
    header_end = kind_end + _U64.size
    if header_end > len(value) or value[offset:prefix_end] != CANONICAL_ENCODING_VERSION:
        raise CanonicalEncodingError("canonical bytes must contain one complete CCAN record.")
    kind = value[prefix_end:kind_end]
    payload_size = _U64.unpack(value[kind_end:header_end])[0]
    payload_end = header_end + payload_size
    if payload_end > len(value):
        raise CanonicalEncodingError("canonical node payload is truncated.")
    payload = value[header_end:payload_end]
    _validate_node_payload(kind, payload)
    return kind, payload, payload_end


def _parse_child_nodes(payload: bytes, count: int) -> tuple[tuple[bytes, bytes, bytes], ...]:
    children: list[tuple[bytes, bytes, bytes]] = []
    offset = 0
    for _ in range(count):
        start = offset
        kind, child_payload, offset = _parse_node(payload, offset)
        children.append((kind, child_payload, payload[start:offset]))
    if offset != len(payload):
        raise CanonicalEncodingError("canonical node contains trailing child bytes.")
    return tuple(children)


def _validate_node_payload(kind: bytes, payload: bytes) -> None:
    if kind == b"B":
        return
    if kind == b"I":
        _integer_payload_value(payload)
        return
    if kind == b"D":
        if len(payload) != _BINARY64.size:
            raise CanonicalEncodingError("canonical binary64 node must contain exactly eight bytes.")
        numeric = _BINARY64.unpack(payload)[0]
        if not math.isfinite(numeric):
            raise CanonicalEncodingError("canonical binary64 node must be finite.")
        if numeric == 0.0 and payload != _BINARY64.pack(0.0):
            raise CanonicalEncodingError("canonical binary64 node must normalize signed zero.")
        return
    if kind == b"?":
        if payload not in (b"\x00", b"\x01"):
            raise CanonicalEncodingError("canonical boolean node must contain exactly zero or one.")
        return
    if kind == b"R":
        children = _parse_child_nodes(payload, 2)
        if any(child_kind != b"I" for child_kind, _, _ in children):
            raise CanonicalEncodingError("canonical rational requires two integer child nodes.")
        numerator = _integer_payload_value(children[0][1])
        denominator = _integer_payload_value(children[1][1])
        if denominator <= 0 or math.gcd(numerator, denominator) != 1:
            raise CanonicalEncodingError("canonical rational must be reduced with positive denominator.")
        return
    if kind == b"S":
        if len(payload) < _U64.size:
            raise CanonicalEncodingError("canonical sequence requires a big-endian element count.")
        count = _U64.unpack(payload[: _U64.size])[0]
        children = _parse_child_nodes(payload[_U64.size :], count)
        if any(child_kind != b"B" for child_kind, _, _ in children):
            raise CanonicalEncodingError("canonical sequence elements must be byte nodes.")
        return
    if kind == b"T":
        children = _parse_child_nodes(payload, 2)
        if any(child_kind != b"B" for child_kind, _, _ in children) or not children[0][1]:
            raise CanonicalEncodingError("canonical tagged union requires a nonempty byte tag and byte payload.")
        return
    if kind == b"M":
        if len(payload) < _U64.size:
            raise CanonicalEncodingError("canonical component map requires a big-endian entry count.")
        count = _U64.unpack(payload[: _U64.size])[0]
        children = _parse_child_nodes(payload[_U64.size :], count * 2)
        keys = children[::2]
        values = children[1::2]
        if any(child_kind != b"B" or not child_payload for child_kind, child_payload, _ in keys):
            raise CanonicalEncodingError("canonical component map requires nonempty byte keys.")
        if any(child_kind != b"B" for child_kind, _, _ in values):
            raise CanonicalEncodingError("canonical component map requires byte values.")
        encoded_keys = tuple(encoded for _, _, encoded in keys)
        if encoded_keys != tuple(sorted(set(encoded_keys))):
            raise CanonicalEncodingError("canonical component map keys must be unique and sorted.")
        return
    raise CanonicalEncodingError("canonical node uses an unknown kind outside the frozen CCAN grammar.")


def require_canonical_record(value: bytes) -> bytes:
    if type(value) is not bytes:
        raise CanonicalEncodingError("canonical value must be bytes.")
    _, _, end = _parse_node(value, 0)
    if end != len(value):
        raise CanonicalEncodingError("canonical bytes must contain one complete CCAN record.")
    return value


def canonical_record_kind(value: bytes) -> bytes:
    require_canonical_record(value)
    return value[len(CANONICAL_ENCODING_VERSION) : len(CANONICAL_ENCODING_VERSION) + 1]


@dataclass(frozen=True)
class ExactRationalV1:
    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if type(self.numerator) is not int or type(self.denominator) is not int:
            raise CanonicalEncodingError("exact rational numerator and denominator must be integers.")
        if self.denominator == 0:
            raise CanonicalEncodingError("exact rational denominator must be nonzero.")
        reduced = Fraction(self.numerator, self.denominator)
        object.__setattr__(self, "numerator", reduced.numerator)
        object.__setattr__(self, "denominator", reduced.denominator)

    @classmethod
    def build(cls, numerator: int, denominator: int) -> Self:
        return cls(numerator, denominator)

    @classmethod
    def from_fraction(cls, value: Fraction) -> Self:
        if type(value) is not Fraction:
            raise CanonicalEncodingError("exact rational source must be fractions.Fraction.")
        return cls(value.numerator, value.denominator)

    @property
    def canonical_bytes(self) -> bytes:
        _require_exact(self, ExactRationalV1, "canonical rational")
        return _node(
            b"R",
            encode_integer(self.numerator) + encode_integer(self.denominator),
        )


@dataclass(frozen=True)
class ExactCenterParameterV1:
    chart: int
    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if type(self.chart) is not int or self.chart < -1 or self.chart > 3:
            raise CanonicalEncodingError("exact center chart must be one integer in [-1, 3].")
        if type(self.numerator) is not int or type(self.denominator) is not int:
            raise CanonicalEncodingError("exact center parameter numerator and denominator must be integers.")
        if self.denominator <= 0 or self.numerator < 0 or self.numerator > self.denominator:
            raise CanonicalEncodingError("exact center parameter must lie in the closed unit interval.")

    @classmethod
    def build(cls, *, chart: int, numerator: int, denominator: int) -> Self:
        return cls(chart, numerator, denominator)

    @property
    def fraction(self) -> Fraction:
        return Fraction(self.numerator, self.denominator)

    @property
    def canonical_bytes(self) -> bytes:
        _require_exact(self, ExactCenterParameterV1, "exact center parameter")
        return encode_tagged_union(
            b"exact-center-parameter-v1",
            encode_component_map(
                {
                    b"chart": encode_integer(self.chart),
                    b"denominator": encode_integer(self.denominator),
                    b"numerator": encode_integer(self.numerator),
                }
            ),
        )


def encode_passage_state(state: PassageState) -> bytes:
    _require_exact(state, PassageState, "passage state")
    try:
        tag = _PASSAGE_TAGS[state]
    except (KeyError, TypeError):
        raise CanonicalEncodingError("passage state must be a closed typed value.") from None
    return encode_tagged_union(tag, b"")


def encode_cut_intent(intent: CutIntent) -> bytes:
    _require_exact(intent, CutIntent, "cut intent")
    try:
        tag = _CUT_INTENT_TAGS[intent]
    except (KeyError, TypeError):
        raise CanonicalEncodingError("cut intent must be a closed typed value.") from None
    return encode_tagged_union(tag, b"")


def encode_material_side(side: MaterialSide) -> bytes:
    _require_exact(side, MaterialSide, "material side")
    try:
        tag = _MATERIAL_SIDE_TAGS[side]
    except (KeyError, TypeError):
        raise CanonicalEncodingError("material side must be a closed typed value.") from None
    return encode_tagged_union(tag, b"")


def encode_circle_orientation(orientation: CircleOrientation) -> bytes:
    _require_exact(orientation, CircleOrientation, "circle orientation")
    try:
        tag = _CIRCLE_ORIENTATION_TAGS[orientation]
    except (KeyError, TypeError):
        raise CanonicalEncodingError("circle orientation must be a closed typed value.") from None
    return encode_tagged_union(tag, b"")


def canonical_point2_bytes(point: Point2[WorldXY]) -> bytes:
    _require_exact(point, Point2, "canonical point")
    return encode_tagged_union(
        b"point2-world-xy-v1",
        encode_sequence((encode_binary64(point.x), encode_binary64(point.y))),
    )


def canonical_vector2_bytes(vector: Vector2[WorldXY]) -> bytes:
    _require_exact(vector, Vector2, "canonical vector")
    return encode_tagged_union(
        b"vector2-world-xy-v1",
        encode_sequence((encode_binary64(vector.x), encode_binary64(vector.y))),
    )


def _scalar_bytes(tag: bytes, value: float) -> bytes:
    return encode_tagged_union(tag, encode_binary64(value))


def canonical_cut_z_bytes(cut_z: CutZ) -> bytes:
    _require_exact(cut_z, CutZ, "canonical cut Z")
    return _scalar_bytes(b"cut-z-mm-v1", cut_z.value)


def canonical_clearance_z_bytes(clearance_z: ClearanceZ) -> bytes:
    _require_exact(clearance_z, ClearanceZ, "canonical clearance Z")
    return _scalar_bytes(b"clearance-z-mm-v1", clearance_z.value)


def _cut_plane_bytes(cut_plane: CutPlane) -> bytes:
    _require_exact(cut_plane, CutPlane, "cut plane")
    _require_exact(cut_plane.cut_z, CutZ, "cut-plane cut Z")
    _require_exact(cut_plane.clearance_z, ClearanceZ, "cut-plane clearance Z")
    return encode_tagged_union(
        b"cut-plane-v1",
        encode_component_map(
            {
                b"clearance-z": canonical_clearance_z_bytes(cut_plane.clearance_z),
                b"cut-z": canonical_cut_z_bytes(cut_plane.cut_z),
            }
        ),
    )


def _segment_motion_bytes(motion: ExactSegmentMotion) -> bytes:
    _require_exact(motion, ExactSegmentMotion, "segment motion")
    _require_exact(motion.start, Point2, "segment start")
    _require_exact(motion.end, Point2, "segment end")
    return encode_tagged_union(
        b"exact-segment-motion-v1",
        encode_component_map(
            {
                b"end": canonical_point2_bytes(motion.end),
                b"start": canonical_point2_bytes(motion.start),
            }
        ),
    )


def _circle_motion_bytes(motion: ExactCircleMotion) -> bytes:
    _require_exact(motion, ExactCircleMotion, "circle motion")
    _require_exact(motion.center, Point2, "circle center")
    _require_exact(motion.phase_vector, Vector2, "circle phase vector")
    if type(motion.clockwise) is not bool:
        raise CanonicalEncodingError("circle orientation must be exact bool.")
    orientation_tag = b"clockwise-v1" if motion.clockwise else b"counterclockwise-v1"
    return encode_tagged_union(
        b"exact-full-circle-motion-v1",
        encode_component_map(
            {
                b"center": canonical_point2_bytes(motion.center),
                b"orientation": encode_tagged_union(orientation_tag, b""),
                b"phase-vector": canonical_vector2_bytes(motion.phase_vector),
            }
        ),
    )


def _engagement_cap_bytes(cap: EngagementCap) -> bytes:
    _require_exact(cap, EngagementCap, "engagement cap")
    if type(cap.chord_ratio_bytes) is not bytes or len(cap.chord_ratio_bytes) != _BINARY64.size:
        raise CanonicalEncodingError("engagement cap must retain one native binary64 surrogate.")
    return encode_tagged_union(b"engagement-cap-v1", encode_bytes(cap.chord_ratio_bytes))


def _candidate_policy_bytes(policy: CandidatePolicy) -> bytes:
    _require_exact(policy, CandidatePolicy, "candidate policy")
    _require_exact(policy.spatial_resolution, Spacing, "candidate spatial resolution")
    _require_exact(policy.radius_resolution, Spacing, "candidate radius resolution")
    _require_exact(policy.minimum_guide_radius, GuideRadius, "candidate minimum guide radius")
    _require_exact(policy.minimum_progress, Spacing, "candidate minimum progress")
    for integer_value in (
        policy.spatial_refinement_levels,
        policy.radius_refinement_levels,
        policy.phase_count,
    ):
        if type(integer_value) is not int:
            raise CanonicalEncodingError("candidate policy integer bounds must be exact int.")
    return encode_tagged_union(
        b"candidate-policy-v1",
        encode_component_map(
            {
                b"minimum-guide-radius": _scalar_bytes(
                    b"guide-radius-mm-v1",
                    policy.minimum_guide_radius.value,
                ),
                b"minimum-progress": _scalar_bytes(b"spacing-mm-v1", policy.minimum_progress.value),
                b"phase-count": encode_integer(policy.phase_count),
                b"radius-refinement-levels": encode_integer(policy.radius_refinement_levels),
                b"radius-resolution": _scalar_bytes(b"spacing-mm-v1", policy.radius_resolution.value),
                b"spatial-refinement-levels": encode_integer(policy.spatial_refinement_levels),
                b"spatial-resolution": _scalar_bytes(b"spacing-mm-v1", policy.spatial_resolution.value),
            }
        ),
    )


def _neck_policy_bytes(policy: NeckPolicy) -> bytes:
    _require_exact(policy, NeckPolicy, "neck policy")
    _require_exact(policy.user_cap, EngagementCap, "neck user cap")
    if type(policy.squared_width_boundaries) is not tuple or any(type(boundary) is not Fraction for boundary in policy.squared_width_boundaries):
        raise CanonicalEncodingError("neck squared-width boundaries must be an exact Fraction tuple.")
    if type(policy.effective_caps) is not tuple:
        raise CanonicalEncodingError("neck effective caps must be an exact tuple.")
    for entry in policy.effective_caps:
        _require_exact(entry, NeckEffectiveCap, "neck effective-cap entry")
        if type(entry.neck_class) is not int:
            raise CanonicalEncodingError("neck class must be exact int.")
        _require_exact(entry.passage_state, PassageState, "neck passage state")
        _require_exact(entry.cap, EngagementCap, "neck effective cap")
    boundaries = encode_sequence(tuple(ExactRationalV1.from_fraction(boundary).canonical_bytes for boundary in policy.squared_width_boundaries))
    entries = encode_sequence(
        tuple(
            encode_tagged_union(
                b"neck-effective-cap-entry-v1",
                encode_component_map(
                    {
                        b"cap": _engagement_cap_bytes(entry.cap),
                        b"neck-class": encode_integer(entry.neck_class),
                        b"passage-state": encode_passage_state(entry.passage_state),
                    }
                ),
            )
            for entry in policy.effective_caps
        )
    )
    return encode_tagged_union(
        b"neck-policy-v1",
        encode_component_map(
            {
                b"effective-caps": entries,
                b"squared-width-boundaries": boundaries,
                b"user-cap": _engagement_cap_bytes(policy.user_cap),
            }
        ),
    )


def canonical_task1_bytes(value: object) -> bytes:
    if type(value) is CutPlane:
        return _cut_plane_bytes(value)
    if type(value) is CutZ:
        return canonical_cut_z_bytes(value)
    if type(value) is ExactSegmentMotion:
        return _segment_motion_bytes(value)
    if type(value) is ExactCircleMotion:
        return _circle_motion_bytes(value)
    if type(value) is EngagementCap:
        return _engagement_cap_bytes(value)
    if type(value) is CandidatePolicy:
        return _candidate_policy_bytes(value)
    if type(value) is NeckPolicy:
        return _neck_policy_bytes(value)
    if type(value) is DepletionPolicy:
        _require_exact(value.chord_bound, ChordBound, "depletion chord bound")
        if type(value.center_count_limit) is not int or type(value.chord_bound_bytes) is not bytes:
            raise CanonicalEncodingError("depletion policy fields must use exact owned types.")
        return encode_tagged_union(
            b"depletion-policy-v1",
            encode_component_map(
                {
                    b"center-count-limit": encode_integer(value.center_count_limit),
                    b"chord-bound": _scalar_bytes(b"chord-bound-mm-v1", value.chord_bound.value),
                }
            ),
        )
    if type(value) is TraversalPolicy:
        if type(value.forward_window) is not int:
            raise CanonicalEncodingError("traversal forward window must be exact int.")
        return encode_tagged_union(
            b"traversal-policy-v1",
            encode_component_map({b"forward-window": encode_integer(value.forward_window)}),
        )
    if type(value) is CutDirectionPolicy:
        _require_exact(value.intent, CutIntent, "cut-direction intent")
        return encode_tagged_union(
            b"cut-direction-policy-v1",
            encode_component_map({b"intent": encode_cut_intent(value.intent)}),
        )
    raise CanonicalEncodingError("value is not a supported Task 1 canonical type.")


def canonical_depletion_witness_bytes(
    *,
    motion: ExactSegmentMotion | ExactCircleMotion,
    policy: DepletionPolicy,
    tool_radius: ToolRadius,
    center_parameters: tuple[ExactCenterParameterV1, ...],
    native_strategy_version: bytes,
    parent_lineage: tuple[bytes, ...],
) -> bytes:
    if type(motion) not in (ExactSegmentMotion, ExactCircleMotion):
        raise CanonicalEncodingError("depletion witness motion must be one exact motion type.")
    _require_exact(policy, DepletionPolicy, "depletion witness policy")
    _require_exact(tool_radius, ToolRadius, "depletion witness tool radius")
    if type(center_parameters) is not tuple or not center_parameters:
        raise CanonicalEncodingError("depletion witness center parameters must be one nonempty exact tuple.")
    if any(type(parameter) is not ExactCenterParameterV1 for parameter in center_parameters):
        raise CanonicalEncodingError("depletion witness center parameters must use exact closed values.")
    if type(native_strategy_version) is not bytes or not native_strategy_version:
        raise CanonicalEncodingError("depletion witness native strategy version must be nonempty bytes.")
    if type(parent_lineage) is not tuple or any(type(digest) is not bytes or len(digest) != 32 for digest in parent_lineage):
        raise CanonicalEncodingError("depletion witness parent lineage must contain exact SHA-256 digests.")
    return encode_tagged_union(
        b"depletion-witness-v1",
        encode_component_map(
            {
                b"center-parameters": encode_sequence(tuple(parameter.canonical_bytes for parameter in center_parameters)),
                b"motion": canonical_task1_bytes(motion),
                b"native-strategy-version": encode_bytes(native_strategy_version),
                b"parent-lineage": encode_sequence(tuple(encode_bytes(digest) for digest in parent_lineage)),
                b"policy": canonical_task1_bytes(policy),
                b"tool-radius": _scalar_bytes(b"tool-radius-mm-v1", tool_radius.value),
            }
        ),
    )


def _exact_ring_area(points: tuple[Point2[WorldXY], ...]) -> Fraction:
    area_twice = Fraction(0, 1)
    for start, end in zip(points, points[1:] + points[:1], strict=True):
        area_twice += Fraction.from_float(start.x) * Fraction.from_float(end.y)
        area_twice -= Fraction.from_float(end.x) * Fraction.from_float(start.y)
    return area_twice


def _normalized_ring_vertices(
    points: Sequence[Point2[WorldXY]],
    *,
    is_outer: bool,
) -> tuple[Point2[WorldXY], ...]:
    try:
        ring = tuple(points)
    except TypeError:
        raise CanonicalEncodingError("polygon ring must be a finite point sequence.") from None
    if any(type(point) is not Point2 for point in ring):
        raise CanonicalEncodingError("polygon ring requires exact Point2 vertices.")
    if ring and ring[0] == ring[-1]:
        ring = ring[:-1]
    if len(ring) < 3:
        raise CanonicalEncodingError("polygon ring requires at least three exact world-XY points.")
    ring = tuple(
        Point2[WorldXY].build(
            0.0 if point.x == 0.0 else point.x,
            0.0 if point.y == 0.0 else point.y,
        )
        for point in ring
    )
    area_twice = _exact_ring_area(ring)
    if area_twice == 0:
        raise CanonicalEncodingError("polygon ring requires nonzero exact area.")
    if (area_twice > 0) is not is_outer:
        ring = tuple(reversed(ring))
    point_bytes = tuple(canonical_point2_bytes(point) for point in ring)
    if len(point_bytes) != len(set(point_bytes)):
        raise CanonicalEncodingError("polygon ring must not contain repeated canonical vertices.")
    rotations = tuple(point_bytes[index:] + point_bytes[:index] for index in range(len(point_bytes)))
    canonical_rotation = min(rotations)
    start_index = rotations.index(canonical_rotation)
    return ring[start_index:] + ring[:start_index]


@dataclass(frozen=True)
class CanonicalRingV1:
    vertices: tuple[Point2[WorldXY], ...]
    is_outer: bool

    def __post_init__(self) -> None:
        _require_exact(self, CanonicalRingV1, "canonical ring")
        if type(self.vertices) is not tuple or any(type(point) is not Point2 for point in self.vertices):
            raise CanonicalEncodingError("canonical ring vertices must be an exact Point2 tuple.")
        if type(self.is_outer) is not bool:
            raise CanonicalEncodingError("canonical ring role must be exact bool.")
        if len(self.vertices) < 3:
            raise CanonicalEncodingError("canonical ring requires at least three vertices.")
        for point in self.vertices:
            if (point.x == 0.0 and math.copysign(1.0, point.x) < 0.0) or (point.y == 0.0 and math.copysign(1.0, point.y) < 0.0):
                raise CanonicalEncodingError("canonical ring vertices must normalize signed zero.")
        point_bytes = tuple(canonical_point2_bytes(point) for point in self.vertices)
        if len(point_bytes) != len(set(point_bytes)):
            raise CanonicalEncodingError("canonical ring must not contain repeated vertices.")
        area_twice = _exact_ring_area(self.vertices)
        if area_twice == 0:
            raise CanonicalEncodingError("canonical ring requires nonzero exact area.")
        if (area_twice > 0) is not self.is_outer:
            raise CanonicalEncodingError("canonical ring winding does not match its exact role.")
        rotations = tuple(point_bytes[index:] + point_bytes[:index] for index in range(len(point_bytes)))
        if point_bytes != min(rotations):
            raise CanonicalEncodingError("canonical ring vertices must use canonical rotation.")

    @classmethod
    def build_outer(cls, points: Sequence[Point2[WorldXY]]) -> Self:
        return cls(_normalized_ring_vertices(points, is_outer=True), True)

    @classmethod
    def build_hole(cls, points: Sequence[Point2[WorldXY]]) -> Self:
        return cls(_normalized_ring_vertices(points, is_outer=False), False)

    @property
    def vertex_count(self) -> int:
        _require_exact(self, CanonicalRingV1, "canonical ring")
        return len(self.vertices)

    @property
    def canonical_bytes(self) -> bytes:
        _require_exact(self, CanonicalRingV1, "canonical ring")
        ring_tag = b"outer-ring-ccw-v1" if self.is_outer else b"hole-ring-cw-v1"
        return encode_tagged_union(
            ring_tag,
            encode_sequence(tuple(canonical_point2_bytes(point) for point in self.vertices)),
        )


def canonical_polygon_bytes(
    outer_ring: Sequence[Point2[WorldXY]],
    holes: Sequence[Sequence[Point2[WorldXY]]],
) -> bytes:
    outer = CanonicalRingV1.build_outer(outer_ring)
    try:
        canonical_holes = tuple(CanonicalRingV1.build_hole(hole) for hole in holes)
    except TypeError:
        raise CanonicalEncodingError("polygon holes must be a finite ring sequence.") from None
    hole_bytes = tuple(sorted(hole.canonical_bytes for hole in canonical_holes))
    if len(hole_bytes) != len(set(hole_bytes)):
        raise CanonicalEncodingError("polygon holes must not contain canonical duplicates.")
    return encode_tagged_union(
        b"polygon-with-holes-world-xy-v1",
        encode_component_map(
            {
                b"holes": encode_sequence(hole_bytes),
                b"outer-ring": outer.canonical_bytes,
            }
        ),
    )
