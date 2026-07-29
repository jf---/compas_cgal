from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal
from typing import NewType
from typing import Self
from typing import cast

from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.canonical import canonical_point2_bytes
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidCandidateLatticeError
from compas_cgal.adaptive.errors import InvalidMathsmProposalError
from compas_cgal.adaptive.errors import InvalidMiddleCurveCursorError
from compas_cgal.adaptive.errors import InvalidMiddleCurveSpanError
from compas_cgal.adaptive.medial_axis import MatEdge
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatSample
from compas_cgal.adaptive.medial_axis import MatSite
from compas_cgal.adaptive.medial_axis import MatSiteId
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.operation import AdvanceTraversalDecision
from compas_cgal.adaptive.operation import CursorIdentity
from compas_cgal.adaptive.operation import EdgeId
from compas_cgal.adaptive.operation import EffectiveCapDecision
from compas_cgal.adaptive.operation import FullCapDecision
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckScope
from compas_cgal.adaptive.operation import NoNeckScope
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.policy import CandidateOrderKey
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.units import ExactMillimetre
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import SquaredMillimetre
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY

MiddleCurveCandidateId = NewType("MiddleCurveCandidateId", bytes)

_SHA256_BYTES = hashlib.sha256().digest_size


def _fraction_bytes(value: Fraction) -> bytes:
    if type(value) is not Fraction:
        raise InvalidCandidateLatticeError("candidate lattice coordinate must be one exact fraction.")
    return ExactRationalV1.from_fraction(value).canonical_bytes


def _level_bytes(levels: tuple[int, ...]) -> bytes:
    return encode_sequence(tuple(encode_integer(level) for level in levels))


def _canonical_levels(
    levels: object,
    name: str,
    error_type: type[ValueError],
) -> tuple[int, ...]:
    if type(levels) is not tuple or not levels or any(type(level) is not int or level < 0 for level in levels) or levels != tuple(sorted(set(levels))):
        raise error_type(f"{name} must be one nonempty canonical level tuple.")
    return levels


def _validate_scope_cap(
    neck_scope: NeckScope,
    effective_cap_decision: EffectiveCapDecision,
) -> None:
    if type(neck_scope) is NoNeckScope:
        if type(effective_cap_decision) is not FullCapDecision:
            raise InvalidCandidateLatticeError("no-neck candidate requires one exact full-cap decision.")
        return
    if type(neck_scope) is OrientedNeckScope:
        if type(effective_cap_decision) is not NeckCapDecision:
            raise InvalidCandidateLatticeError("oriented-neck candidate requires one exact neck-cap decision.")
        return
    raise InvalidCandidateLatticeError("candidate neck scope must use the closed exact type.")


def _distance(first: Point2[WorldXY], second: Point2[WorldXY]) -> float:
    return math.hypot(second.x - first.x, second.y - first.y)


def _boundary_footpoint(
    middle_point: Point2[WorldXY],
    generator_site: MatSite,
) -> Point2[WorldXY]:
    if generator_site.kind == "point":
        return generator_site.source
    target = generator_site.target
    if target is None:
        raise InvalidMathsmProposalError("open-segment MAT site omits its target point.")
    segment_x = target.x - generator_site.source.x
    segment_y = target.y - generator_site.source.y
    squared_length = segment_x * segment_x + segment_y * segment_y
    if squared_length == 0.0:
        raise InvalidMathsmProposalError("open-segment MAT site has zero reporting length.")
    parameter = ((middle_point.x - generator_site.source.x) * segment_x + (middle_point.y - generator_site.source.y) * segment_y) / squared_length
    if not 0.0 < parameter < 1.0:
        raise InvalidMathsmProposalError("MAT station has no interior footpoint on its open-segment generator.")
    return Point2[WorldXY].build(
        generator_site.source.x + parameter * segment_x,
        generator_site.source.y + parameter * segment_y,
    )


def _mathsm_geometry(
    *,
    middle_point: Point2[WorldXY],
    boundary_footpoint: Point2[WorldXY],
    tool_radius: ToolRadius,
    guide_radius: Fraction,
) -> tuple[Point2[WorldXY], Point2[WorldXY], Fraction]:
    distance = _distance(boundary_footpoint, middle_point)
    if distance < tool_radius.value:
        raise InvalidMathsmProposalError("MAT reporting station lies inside the tool-radius boundary.")
    maximum = Fraction.from_float((distance - tool_radius.value) / 2.0)
    if guide_radius <= 0 or guide_radius > maximum:
        raise InvalidMathsmProposalError("guide radius lies outside the positive one-sided MATHSM interval.")
    if distance == 0.0:
        raise InvalidMathsmProposalError("MATHSM direction requires distinct boundary and MAT points.")
    direction_x = (middle_point.x - boundary_footpoint.x) / distance
    direction_y = (middle_point.y - boundary_footpoint.y) / distance
    contact_point = Point2[WorldXY].build(
        boundary_footpoint.x + tool_radius.value * direction_x,
        boundary_footpoint.y + tool_radius.value * direction_y,
    )
    radius = float(guide_radius)
    center = Point2[WorldXY].build(
        contact_point.x + radius * direction_x,
        contact_point.y + radius * direction_y,
    )
    return contact_point, center, maximum


def _phase_vector(
    *,
    contact_point: Point2[WorldXY],
    center: Point2[WorldXY],
    phase_index: int,
    phase_count: int,
) -> Vector2[WorldXY]:
    base_x = contact_point.x - center.x
    base_y = contact_point.y - center.y
    if phase_index == 0:
        return Vector2[WorldXY].build(base_x, base_y)
    angle = math.tau * phase_index / phase_count
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return Vector2[WorldXY].build(
        cosine * base_x - sine * base_y,
        sine * base_x + cosine * base_y,
    )


@dataclass(frozen=True)
class MiddleCurveSpan:
    axis: MedialAxis
    edge: MatEdge
    cursor_before: MatSample | DerivedCandidateCursor
    cursor_limit: MatSample

    def __post_init__(self) -> None:
        if type(self.axis) is not MedialAxis:
            raise InvalidMiddleCurveSpanError("middle-curve span requires one exact typed MAT.")
        if type(self.edge) is not MatEdge:
            raise InvalidMiddleCurveSpanError("middle-curve span requires one exact MAT edge.")
        if type(self.cursor_before) not in (MatSample, DerivedCandidateCursor) or type(self.cursor_limit) is not MatSample:
            raise InvalidMiddleCurveSpanError("middle-curve span requires one owned native/derived cursor and one native limit.")
        owned_edge = self.axis.edge_by_id.get(self.edge.identity)
        if owned_edge != self.edge:
            raise InvalidMiddleCurveSpanError("middle-curve span edge is not owned by its MAT.")
        if self.cursor_limit not in self.axis.samples:
            raise InvalidMiddleCurveSpanError("middle-curve span limit is not owned by its MAT.")
        if self.cursor_before.edge_id != self.edge.identity or self.cursor_limit.edge_id != self.edge.identity:
            raise InvalidMiddleCurveSpanError("middle-curve span cursors must lie on one exact MAT edge.")
        if type(self.cursor_before) is MatSample:
            if self.cursor_before not in self.axis.samples:
                raise InvalidMiddleCurveSpanError("middle-curve native cursor is not owned by its MAT.")
            ordinal_delta = self.cursor_limit.ordinal_on_edge - self.cursor_before.ordinal_on_edge
            if ordinal_delta == 0:
                raise InvalidMiddleCurveSpanError("middle-curve native span requires distinct sample ordinals.")
        else:
            derived_cursor = cast(DerivedCandidateCursor, self.cursor_before)
            if derived_cursor.axis is not self.axis:
                raise InvalidMiddleCurveSpanError("middle-curve derived cursor belongs to another MAT owner.")
            ordinal_delta = self.cursor_limit.ordinal_on_edge - derived_cursor.next_limit_ordinal
            if ordinal_delta != 0 and (1 if ordinal_delta > 0 else -1) != derived_cursor.ordinal_step:
                raise InvalidMiddleCurveSpanError("middle-curve continuation contradicts its retained sample direction.")
        if self.cursor_before.point == self.cursor_limit.point:
            raise InvalidMiddleCurveSpanError("middle-curve span requires nonzero reporting progress.")

    @classmethod
    def build(
        cls,
        *,
        axis: MedialAxis,
        cursor_before: MatSample | DerivedCandidateCursor,
        cursor_limit: MatSample,
    ) -> Self:
        if type(axis) is not MedialAxis or type(cursor_before) not in (MatSample, DerivedCandidateCursor):
            raise InvalidMiddleCurveSpanError("middle-curve span factory requires typed MAT inputs.")
        edge = axis.edge_by_id.get(cursor_before.edge_id)
        if edge is None:
            raise InvalidMiddleCurveSpanError("middle-curve cursor references an unknown MAT edge.")
        return cls(axis, edge, cursor_before, cursor_limit)

    @property
    def reported_length(self) -> ExactMillimetre:
        return ExactMillimetre(Fraction.from_float(_distance(self.cursor_before.point, self.cursor_limit.point)))

    @property
    def ordinal_step(self) -> Literal[-1, 1]:
        """Return the monotone native-sample direction of this span."""
        if type(self.cursor_before) is MatSample:
            return 1 if self.cursor_limit.ordinal_on_edge > self.cursor_before.ordinal_on_edge else -1
        return cast(
            DerivedCandidateCursor,
            self.cursor_before,
        ).ordinal_step


@dataclass(frozen=True)
class MathsmCircleProposal:
    generator_site: MatSite
    middle_point: Point2[WorldXY]
    boundary_footpoint: Point2[WorldXY]
    mathsm_contact_point: Point2[WorldXY]
    tool_radius: ToolRadius
    guide_radius: ExactMillimetre
    maximum_guide_radius: ExactMillimetre
    phase_index: int
    phase_count: int
    circle_orientation: CircleOrientation
    motion: ExactCircleMotion

    def __post_init__(self) -> None:
        if type(self.generator_site) is not MatSite:
            raise InvalidMathsmProposalError("MATHSM proposal requires one exact generator site.")
        if type(self.middle_point) is not Point2 or type(self.boundary_footpoint) is not Point2 or type(self.mathsm_contact_point) is not Point2:
            raise InvalidMathsmProposalError("MATHSM proposal points must use world-XY typed geometry.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidMathsmProposalError("MATHSM proposal tool radius must be typed.")
        if type(self.guide_radius) is not Fraction or self.guide_radius <= 0:
            raise InvalidMathsmProposalError("MATHSM proposal guide radius must be exact and positive.")
        if type(self.maximum_guide_radius) is not Fraction or self.maximum_guide_radius <= 0 or self.guide_radius > self.maximum_guide_radius:
            raise InvalidMathsmProposalError("MATHSM maximum guide radius must bind the selected radius.")
        if type(self.phase_index) is not int or type(self.phase_count) is not int or self.phase_count <= 0 or self.phase_index < 0 or self.phase_index >= self.phase_count:
            raise InvalidMathsmProposalError("MATHSM phase must lie in its finite closed index domain.")
        if type(self.circle_orientation) is not CircleOrientation:
            raise InvalidMathsmProposalError("MATHSM circle orientation must use the closed exact type.")
        if type(self.motion) is not ExactCircleMotion:
            raise InvalidMathsmProposalError("MATHSM proposal must own one exact circle motion.")

        expected_footpoint = _boundary_footpoint(
            self.middle_point,
            self.generator_site,
        )
        expected_contact, expected_center, expected_maximum = _mathsm_geometry(
            middle_point=self.middle_point,
            boundary_footpoint=expected_footpoint,
            tool_radius=self.tool_radius,
            guide_radius=self.guide_radius,
        )
        expected_phase = _phase_vector(
            contact_point=expected_contact,
            center=expected_center,
            phase_index=self.phase_index,
            phase_count=self.phase_count,
        )
        if (
            self.boundary_footpoint != expected_footpoint
            or self.mathsm_contact_point != expected_contact
            or self.motion.center != expected_center
            or self.motion.phase_vector != expected_phase
            or self.motion.clockwise != self.circle_orientation.clockwise
            or self.maximum_guide_radius != expected_maximum
        ):
            raise InvalidMathsmProposalError("MATHSM proposal geometry contradicts its site, radius, or phase.")

    @classmethod
    def build(
        cls,
        *,
        generator_site: MatSite,
        middle_point: Point2[WorldXY],
        tool_radius: ToolRadius,
        guide_radius: ExactMillimetre,
        phase_index: int,
        phase_count: int,
        circle_orientation: CircleOrientation,
    ) -> Self:
        if type(generator_site) is not MatSite:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one exact generator site.")
        if type(middle_point) is not Point2:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one world-XY middle point.")
        if type(tool_radius) is not ToolRadius:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one typed tool radius.")
        if type(guide_radius) is not Fraction or guide_radius <= 0:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one positive exact guide radius.")
        if type(phase_index) is not int or type(phase_count) is not int or phase_count <= 0 or phase_index < 0 or phase_index >= phase_count:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one finite phase index.")
        if type(circle_orientation) is not CircleOrientation:
            raise InvalidMathsmProposalError("MATHSM proposal factory requires one closed circle orientation.")
        boundary_footpoint = _boundary_footpoint(
            middle_point,
            generator_site,
        )
        contact, center, maximum = _mathsm_geometry(
            middle_point=middle_point,
            boundary_footpoint=boundary_footpoint,
            tool_radius=tool_radius,
            guide_radius=guide_radius,
        )
        phase_vector = _phase_vector(
            contact_point=contact,
            center=center,
            phase_index=phase_index,
            phase_count=phase_count,
        )
        motion = ExactCircleMotion.build(
            center,
            phase_vector,
            circle_orientation.clockwise,
        )
        return cls(
            generator_site,
            middle_point,
            boundary_footpoint,
            contact,
            tool_radius,
            guide_radius,
            ExactMillimetre(maximum),
            phase_index,
            phase_count,
            circle_orientation,
            motion,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return encode_tagged_union(
            b"one-sided-mathsm-proposal-v1",
            encode_component_map(
                {
                    b"boundary-footpoint": canonical_point2_bytes(self.boundary_footpoint),
                    b"circle-motion": canonical_task1_bytes(self.motion),
                    b"generator-site-id": encode_bytes(bytes(self.generator_site.identity)),
                    b"guide-radius": _fraction_bytes(self.guide_radius),
                    b"mathsm-contact-point": canonical_point2_bytes(self.mathsm_contact_point),
                    b"maximum-guide-radius": _fraction_bytes(self.maximum_guide_radius),
                    b"middle-point": canonical_point2_bytes(self.middle_point),
                    b"phase-count": encode_integer(self.phase_count),
                    b"phase-index": encode_integer(self.phase_index),
                    b"tool-radius": encode_binary64(float(self.tool_radius.value)),
                }
            ),
        )


def _candidate_record_bytes(
    *,
    proposal: MathsmCircleProposal,
    policy: CandidatePolicy,
    spatial_progress: ExactMillimetre,
    spatial_levels: tuple[int, ...],
    radius_levels: tuple[int, ...],
    cursor_limit_identity: CursorIdentity,
    neck_scope: NeckScope,
    effective_cap_decision: EffectiveCapDecision,
    traversal_decision: AdvanceTraversalDecision,
) -> bytes:
    return encode_tagged_union(
        b"middle-curve-candidate-v1",
        encode_component_map(
            {
                b"cap-decision": effective_cap_decision.canonical_bytes,
                b"cursor-limit": encode_bytes(cursor_limit_identity),
                b"mathsm-proposal": proposal.canonical_bytes,
                b"neck-scope": neck_scope.canonical_bytes,
                b"policy": canonical_task1_bytes(policy),
                b"radius-levels": _level_bytes(radius_levels),
                b"spatial-levels": _level_bytes(spatial_levels),
                b"spatial-progress": _fraction_bytes(spatial_progress),
                b"traversal": traversal_decision.canonical_bytes,
            }
        ),
    )


def _validate_candidate_inputs(
    *,
    proposal: object,
    policy: object,
    spatial_progress: object,
    spatial_levels: object,
    radius_levels: object,
    cursor_limit_identity: object,
    neck_scope: NeckScope,
    effective_cap_decision: EffectiveCapDecision,
    traversal_decision: object,
) -> None:
    if type(proposal) is not MathsmCircleProposal:
        raise InvalidCandidateLatticeError("candidate must own one exact MATHSM proposal.")
    if type(policy) is not CandidatePolicy:
        raise InvalidCandidateLatticeError("candidate must retain its finite policy.")
    if type(spatial_progress) is not Fraction or spatial_progress <= 0:
        raise InvalidCandidateLatticeError("candidate progress must be exact and positive.")
    _canonical_levels(
        spatial_levels,
        "candidate spatial levels",
        InvalidCandidateLatticeError,
    )
    _canonical_levels(
        radius_levels,
        "candidate radius levels",
        InvalidCandidateLatticeError,
    )
    if type(cursor_limit_identity) is not bytes or not cursor_limit_identity:
        raise InvalidCandidateLatticeError("candidate must bind its exact cursor limit.")
    _validate_scope_cap(neck_scope, effective_cap_decision)
    if type(traversal_decision) is not AdvanceTraversalDecision:
        raise InvalidCandidateLatticeError("candidate must own one exact advance traversal.")


@dataclass(frozen=True)
class MiddleCurveCandidate:
    identity: MiddleCurveCandidateId
    proposal: MathsmCircleProposal
    policy: CandidatePolicy
    spatial_progress: ExactMillimetre
    spatial_levels: tuple[int, ...]
    radius_levels: tuple[int, ...]
    cursor_limit_identity: CursorIdentity
    neck_scope: NeckScope
    effective_cap_decision: EffectiveCapDecision
    traversal_decision: AdvanceTraversalDecision

    def __post_init__(self) -> None:
        if type(self.identity) is not bytes or len(self.identity) != _SHA256_BYTES:
            raise InvalidCandidateLatticeError("candidate identity must be one SHA-256 digest.")
        _validate_candidate_inputs(
            proposal=self.proposal,
            policy=self.policy,
            spatial_progress=self.spatial_progress,
            spatial_levels=self.spatial_levels,
            radius_levels=self.radius_levels,
            cursor_limit_identity=self.cursor_limit_identity,
            neck_scope=self.neck_scope,
            effective_cap_decision=self.effective_cap_decision,
            traversal_decision=self.traversal_decision,
        )
        expected_identity = hashlib.sha256(self.canonical_bytes).digest()
        if self.identity != expected_identity:
            raise InvalidCandidateLatticeError("candidate identity contradicts its canonical decision record.")

    @classmethod
    def build(
        cls,
        *,
        proposal: MathsmCircleProposal,
        policy: CandidatePolicy,
        spatial_progress: ExactMillimetre,
        spatial_levels: tuple[int, ...],
        radius_levels: tuple[int, ...],
        cursor_limit_identity: CursorIdentity,
        neck_scope: NeckScope,
        effective_cap_decision: EffectiveCapDecision,
        traversal_decision: AdvanceTraversalDecision,
    ) -> Self:
        _validate_candidate_inputs(
            proposal=proposal,
            policy=policy,
            spatial_progress=spatial_progress,
            spatial_levels=spatial_levels,
            radius_levels=radius_levels,
            cursor_limit_identity=cursor_limit_identity,
            neck_scope=neck_scope,
            effective_cap_decision=effective_cap_decision,
            traversal_decision=traversal_decision,
        )
        identity = MiddleCurveCandidateId(
            hashlib.sha256(
                _candidate_record_bytes(
                    proposal=proposal,
                    policy=policy,
                    spatial_progress=spatial_progress,
                    spatial_levels=spatial_levels,
                    radius_levels=radius_levels,
                    cursor_limit_identity=cursor_limit_identity,
                    neck_scope=neck_scope,
                    effective_cap_decision=effective_cap_decision,
                    traversal_decision=traversal_decision,
                )
            ).digest()
        )
        return cls(
            identity,
            proposal,
            policy,
            spatial_progress,
            spatial_levels,
            radius_levels,
            cursor_limit_identity,
            neck_scope,
            effective_cap_decision,
            traversal_decision,
        )

    @property
    def canonical_bytes(self) -> bytes:
        return _candidate_record_bytes(
            proposal=self.proposal,
            policy=self.policy,
            spatial_progress=self.spatial_progress,
            spatial_levels=self.spatial_levels,
            radius_levels=self.radius_levels,
            cursor_limit_identity=self.cursor_limit_identity,
            neck_scope=self.neck_scope,
            effective_cap_decision=self.effective_cap_decision,
            traversal_decision=self.traversal_decision,
        )

    @property
    def order_key(self) -> CandidateOrderKey:
        return CandidateOrderKey.build(
            progress=self.spatial_progress,
            squared_radius=SquaredMillimetre(self.proposal.motion.squared_radius),
            canonical_identity=bytes(self.identity),
        )

    @property
    def middle_point(self) -> Point2[WorldXY]:
        return self.proposal.middle_point

    @property
    def boundary_footpoint(self) -> Point2[WorldXY]:
        return self.proposal.boundary_footpoint

    @property
    def mathsm_contact_point(self) -> Point2[WorldXY]:
        return self.proposal.mathsm_contact_point

    @property
    def guide_radius(self) -> ExactMillimetre:
        return self.proposal.guide_radius

    @property
    def phase_index(self) -> int:
        return self.proposal.phase_index

    @property
    def generator_site_id(self) -> MatSiteId:
        return self.proposal.generator_site.identity

    @property
    def motion(self) -> ExactCircleMotion:
        return self.proposal.motion


@dataclass(frozen=True)
class DerivedCandidateCursor:
    """Proof-carrying continuation point produced by one lattice candidate."""

    span: MiddleCurveSpan
    candidate: MiddleCurveCandidate

    def __post_init__(self) -> None:
        if type(self.span) is not MiddleCurveSpan or type(self.candidate) is not MiddleCurveCandidate:
            raise InvalidMiddleCurveCursorError("derived cursor requires one exact span and candidate.")
        traversal = self.candidate.traversal_decision
        if traversal.makes_cursor_terminal:
            raise InvalidMiddleCurveCursorError("terminal candidate cannot produce a continuation cursor.")
        if self.candidate.spatial_progress >= self.span.reported_length:
            raise InvalidMiddleCurveCursorError("native span endpoint must retain its owned MAT sample cursor.")
        if (
            traversal.component_id != self.span.axis.component_by_edge_id[self.span.edge.identity]
            or traversal.edge_id != EdgeId(bytes(self.span.edge.identity))
            or traversal.branch_id != self.span.edge.branch_id
            or traversal.cursor_before != self.span.cursor_before.cursor_identity
            or self.candidate.cursor_limit_identity != self.span.cursor_limit.cursor_identity
        ):
            raise InvalidMiddleCurveCursorError("candidate traversal does not belong to its claimed exact span.")
        if self.candidate.generator_site_id not in self.span.edge.generator_site_ids:
            raise InvalidMiddleCurveCursorError("candidate generator does not belong to its claimed MAT edge.")
        if self.candidate.middle_point != _middle_point(
            self.span,
            self.candidate.spatial_progress,
        ) or traversal.cursor_after != _cursor_after(
            self.span,
            self.candidate.spatial_progress,
        ):
            raise InvalidMiddleCurveCursorError("candidate geometry or cursor identity contradicts its exact span.")

    @classmethod
    def build(
        cls,
        *,
        span: MiddleCurveSpan,
        candidate: MiddleCurveCandidate,
    ) -> Self:
        return cls(span, candidate)

    @property
    def axis(self) -> MedialAxis:
        return self.span.axis

    @property
    def edge_id(self) -> MatEdgeId:
        return self.span.edge.identity

    @property
    def cursor_identity(self) -> CursorIdentity:
        return self.candidate.traversal_decision.cursor_after

    @property
    def point(self) -> Point2[WorldXY]:
        return self.candidate.middle_point

    @property
    def ordinal_step(self) -> Literal[-1, 1]:
        """Retain the direction selected by the producing native span."""
        return self.span.ordinal_step

    @property
    def next_limit_ordinal(self) -> int:
        """Return the nearest native ordinal legal for a continuation."""
        return self.span.cursor_limit.ordinal_on_edge

    @property
    def minimum_limit_ordinal(self) -> int:
        """Return the forward-only lower limit retained by this cursor."""
        if self.ordinal_step != 1:
            raise InvalidMiddleCurveCursorError("reverse derived cursor has no minimum native limit ordinal.")
        return self.next_limit_ordinal


def _refined_spatial_values(
    span: MiddleCurveSpan,
    policy: CandidatePolicy,
) -> tuple[tuple[ExactMillimetre, tuple[int, ...]], ...]:
    length = Fraction(span.reported_length)
    minimum = Fraction.from_float(policy.minimum_progress.value)
    base_quantum = Fraction.from_float(policy.spatial_resolution.value)
    levels_by_value: dict[Fraction, set[int]] = {}
    for level in range(policy.spatial_refinement_levels):
        quantum = base_quantum / (2**level)
        ratio_numerator = minimum.numerator * quantum.denominator
        ratio_denominator = minimum.denominator * quantum.numerator
        first_multiple = max(
            1,
            (ratio_numerator + ratio_denominator - 1) // ratio_denominator,
        )
        value = first_multiple * quantum
        while value < length:
            levels_by_value.setdefault(value, set()).add(level)
            value += quantum
        if length >= minimum:
            levels_by_value.setdefault(length, set()).add(level)
    return tuple(
        (
            ExactMillimetre(value),
            tuple(sorted(levels)),
        )
        for value, levels in sorted(levels_by_value.items())
    )


def _refined_radius_values(
    *,
    maximum: Fraction,
    policy: CandidatePolicy,
) -> tuple[tuple[ExactMillimetre, tuple[int, ...]], ...]:
    minimum = Fraction.from_float(policy.minimum_guide_radius.value)
    if maximum < minimum:
        return ()
    base_quantum = Fraction.from_float(policy.radius_resolution.value)
    levels_by_value: dict[Fraction, set[int]] = {}
    for level in range(policy.radius_refinement_levels):
        quantum = base_quantum / (2**level)
        value = minimum
        while value < maximum:
            levels_by_value.setdefault(value, set()).add(level)
            value += quantum
        levels_by_value.setdefault(maximum, set()).add(level)
    return tuple(
        (
            ExactMillimetre(value),
            tuple(sorted(levels)),
        )
        for value, levels in sorted(levels_by_value.items())
    )


def _middle_point(
    span: MiddleCurveSpan,
    progress: ExactMillimetre,
) -> Point2[WorldXY]:
    if progress == span.reported_length:
        return span.cursor_limit.point
    length = float(span.reported_length)
    ratio = float(progress) / length
    start = span.cursor_before.point
    limit = span.cursor_limit.point
    if span.edge.curve_kind == "parabola":
        sites = tuple(span.axis.site_by_id[site_id] for site_id in span.edge.generator_site_ids)
        point_sites = tuple(site for site in sites if site.kind == "point")
        segment_sites = tuple(site for site in sites if site.kind == "open-segment")
        if len(point_sites) != 1 or len(segment_sites) != 1 or segment_sites[0].target is None:
            raise InvalidMiddleCurveSpanError("parabolic MAT span requires one point and one open-segment generator.")
        focus = point_sites[0].source
        directrix_source = segment_sites[0].source
        directrix_target = segment_sites[0].target
        tangent_x = directrix_target.x - directrix_source.x
        tangent_y = directrix_target.y - directrix_source.y
        directrix_length = math.hypot(tangent_x, tangent_y)
        if directrix_length == 0.0:
            raise InvalidMiddleCurveSpanError("parabolic MAT span has a degenerate directrix.")
        tangent_x /= directrix_length
        tangent_y /= directrix_length
        normal_x = -tangent_y
        normal_y = tangent_x

        def coordinates(point: Point2[WorldXY]) -> tuple[float, float]:
            offset_x = point.x - directrix_source.x
            offset_y = point.y - directrix_source.y
            return (
                offset_x * tangent_x + offset_y * tangent_y,
                offset_x * normal_x + offset_y * normal_y,
            )

        start_u, _ = coordinates(start)
        limit_u, _ = coordinates(limit)
        focus_u, focus_v = coordinates(focus)
        if focus_v == 0.0:
            raise InvalidMiddleCurveSpanError("parabolic MAT focus lies on its directrix.")
        parameter_u = start_u + ratio * (limit_u - start_u)
        parameter_v = ((parameter_u - focus_u) ** 2 + focus_v**2) / (2.0 * focus_v)
        return Point2[WorldXY].build(
            directrix_source.x + parameter_u * tangent_x + parameter_v * normal_x,
            directrix_source.y + parameter_u * tangent_y + parameter_v * normal_y,
        )
    return Point2[WorldXY].build(
        start.x + ratio * (limit.x - start.x),
        start.y + ratio * (limit.y - start.y),
    )


def _cursor_after(
    span: MiddleCurveSpan,
    progress: ExactMillimetre,
) -> CursorIdentity:
    if progress == span.reported_length:
        return span.cursor_limit.cursor_identity
    canonical = encode_tagged_union(
        b"middle-curve-cursor-v1",
        encode_component_map(
            {
                b"cursor-before": encode_bytes(bytes(span.cursor_before.cursor_identity)),
                b"cursor-limit": encode_bytes(bytes(span.cursor_limit.cursor_identity)),
                b"edge-id": encode_bytes(bytes(span.edge.identity)),
                b"progress": _fraction_bytes(progress),
            }
        ),
    )
    return CursorIdentity(hashlib.sha256(canonical).digest())


def _terminalize_last_feasible(
    candidates: list[MiddleCurveCandidate],
    *,
    terminal_limit: bool,
) -> list[MiddleCurveCandidate]:
    """Close a terminal span when its exact endpoint admits no positive circle.

    The complete finite lattice is already materialized. Therefore its greatest
    feasible progress is an exhaustive traversal boundary, not a sampled guess.
    """
    if not terminal_limit or not candidates or any(candidate.traversal_decision.makes_cursor_terminal for candidate in candidates):
        return candidates
    maximum_progress = max(candidate.spatial_progress for candidate in candidates)
    terminalized: list[MiddleCurveCandidate] = []
    for candidate in candidates:
        if candidate.spatial_progress != maximum_progress:
            terminalized.append(candidate)
            continue
        traversal = candidate.traversal_decision
        terminalized.append(
            MiddleCurveCandidate.build(
                proposal=candidate.proposal,
                policy=candidate.policy,
                spatial_progress=candidate.spatial_progress,
                spatial_levels=candidate.spatial_levels,
                radius_levels=candidate.radius_levels,
                cursor_limit_identity=candidate.cursor_limit_identity,
                neck_scope=candidate.neck_scope,
                effective_cap_decision=candidate.effective_cap_decision,
                traversal_decision=AdvanceTraversalDecision.build(
                    component_id=traversal.component_id,
                    edge_id=traversal.edge_id,
                    branch_id=traversal.branch_id,
                    cursor_before=traversal.cursor_before,
                    cursor_after=traversal.cursor_after,
                    makes_cursor_terminal=True,
                ),
            )
        )
    return terminalized


def enumerate_middle_curve_candidates(
    *,
    span: MiddleCurveSpan,
    policy: CandidatePolicy,
    circle_orientation: CircleOrientation,
    neck_scope: NeckScope,
    effective_cap_decision: EffectiveCapDecision,
    makes_cursor_terminal_at_limit: bool,
) -> tuple[MiddleCurveCandidate, ...]:
    if type(span) is not MiddleCurveSpan:
        raise InvalidCandidateLatticeError("candidate enumeration requires one exact middle-curve span.")
    if type(policy) is not CandidatePolicy:
        raise InvalidCandidateLatticeError("candidate enumeration requires one finite policy.")
    if type(circle_orientation) is not CircleOrientation:
        raise InvalidCandidateLatticeError("candidate enumeration requires one closed circle orientation.")
    _validate_scope_cap(neck_scope, effective_cap_decision)
    if type(makes_cursor_terminal_at_limit) is not bool:
        raise InvalidCandidateLatticeError("cursor-limit terminal intent must be explicitly boolean.")

    candidates: list[MiddleCurveCandidate] = []
    for progress, spatial_levels in _refined_spatial_values(span, policy):
        middle_point = _middle_point(span, progress)
        cursor_after = _cursor_after(span, progress)
        traversal = AdvanceTraversalDecision.build(
            component_id=span.axis.component_by_edge_id[span.edge.identity],
            edge_id=EdgeId(bytes(span.edge.identity)),
            branch_id=span.edge.branch_id,
            cursor_before=span.cursor_before.cursor_identity,
            cursor_after=cursor_after,
            makes_cursor_terminal=(makes_cursor_terminal_at_limit and progress == span.reported_length),
        )
        for site_id in span.edge.generator_site_ids:
            site = span.axis.site_by_id[site_id]
            footpoint = _boundary_footpoint(middle_point, site)
            distance = _distance(footpoint, middle_point)
            if distance < span.axis.tool_radius.value:
                raise InvalidMathsmProposalError("MAT span enters the inadmissible tool-radius region.")
            maximum = Fraction.from_float((distance - span.axis.tool_radius.value) / 2.0)
            for guide_radius, radius_levels in _refined_radius_values(
                maximum=maximum,
                policy=policy,
            ):
                for phase_index in range(policy.phase_count):
                    proposal = MathsmCircleProposal.build(
                        generator_site=site,
                        middle_point=middle_point,
                        tool_radius=span.axis.tool_radius,
                        guide_radius=guide_radius,
                        phase_index=phase_index,
                        phase_count=policy.phase_count,
                        circle_orientation=circle_orientation,
                    )
                    candidates.append(
                        MiddleCurveCandidate.build(
                            proposal=proposal,
                            policy=policy,
                            spatial_progress=progress,
                            spatial_levels=spatial_levels,
                            radius_levels=radius_levels,
                            cursor_limit_identity=span.cursor_limit.cursor_identity,
                            neck_scope=neck_scope,
                            effective_cap_decision=effective_cap_decision,
                            traversal_decision=traversal,
                        )
                    )
    return policy.order_candidates(
        _terminalize_last_feasible(
            candidates,
            terminal_limit=makes_cursor_terminal_at_limit,
        ),
        key=lambda candidate: candidate.order_key,
    )
