from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from enum import auto
from fractions import Fraction
from typing import Callable
from typing import NewType
from typing import Self
from typing import TypeVar

from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import InvalidCandidatePolicyError
from compas_cgal.adaptive.errors import InvalidCutDirectionPolicyError
from compas_cgal.adaptive.errors import InvalidDepletionPolicyError
from compas_cgal.adaptive.errors import InvalidMatSamplingPolicyError
from compas_cgal.adaptive.errors import InvalidNeckPolicyError
from compas_cgal.adaptive.errors import InvalidTraversalPolicyError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ExactMillimetre
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre

ComponentId = NewType("ComponentId", bytes)
BranchId = NewType("BranchId", bytes)
CandidateT = TypeVar("CandidateT")


class PassageState(Enum):
    UNVISITED = auto()
    FIRST_PASS_COMPLETE = auto()
    SECOND_PASS_COMPLETE = auto()
    TERMINAL = auto()

    @property
    def rank(self) -> int:
        return {
            PassageState.UNVISITED: 0,
            PassageState.FIRST_PASS_COMPLETE: 1,
            PassageState.SECOND_PASS_COMPLETE: 2,
            PassageState.TERMINAL: 3,
        }[self]


ACTIVE_PASSAGE_STATES = (
    PassageState.UNVISITED,
    PassageState.FIRST_PASS_COMPLETE,
    PassageState.SECOND_PASS_COMPLETE,
)


class CutIntent(Enum):
    CLIMB = auto()
    CONVENTIONAL = auto()


class MaterialSide(Enum):
    """Radial side of a circular guide that still contains material."""

    INSIDE = auto()
    OUTSIDE = auto()


class CircleOrientation(Enum):
    CLOCKWISE = auto()
    COUNTERCLOCKWISE = auto()

    @property
    def clockwise(self) -> bool:
        return self is CircleOrientation.CLOCKWISE


def _positive_int(value: object, name: str, error_type: type[ValueError]) -> int:
    if type(value) is not int or value <= 0:
        raise error_type(f"{name} must be a positive finite integer.")
    return value


@dataclass(frozen=True)
class CandidateOrderKey:
    progress: ExactMillimetre
    squared_radius: SquaredMillimetre
    canonical_identity: bytes

    def __post_init__(self) -> None:
        if type(self.progress) is not Fraction or self.progress < 0:
            raise InvalidCandidatePolicyError("candidate progress must be an exact non-negative length.")
        if type(self.squared_radius) is not Fraction or self.squared_radius <= 0:
            raise InvalidCandidatePolicyError("candidate squared radius must be exact and positive.")
        if type(self.canonical_identity) is not bytes or not self.canonical_identity:
            raise InvalidCandidatePolicyError("candidate canonical identity must be nonempty bytes.")

    @classmethod
    def build(
        cls,
        *,
        progress: ExactMillimetre,
        squared_radius: SquaredMillimetre,
        canonical_identity: bytes,
    ) -> Self:
        return cls(progress, squared_radius, canonical_identity)


@dataclass(frozen=True)
class CandidatePolicy:
    """Finite candidate lattice bounds.

    Winner selection is invariant: furthest progress, largest radius, then
    canonical identity. It is deliberately not a configurable policy field.
    """

    spatial_resolution: Spacing
    spatial_refinement_levels: int
    radius_resolution: Spacing
    radius_refinement_levels: int
    phase_count: int
    minimum_guide_radius: GuideRadius
    minimum_progress: Spacing

    def __post_init__(self) -> None:
        if not isinstance(self.spatial_resolution, Spacing):
            raise InvalidCandidatePolicyError("spatial resolution must be typed spacing.")
        if not isinstance(self.radius_resolution, Spacing):
            raise InvalidCandidatePolicyError("radius resolution must be typed spacing.")
        if not isinstance(self.minimum_progress, Spacing):
            raise InvalidCandidatePolicyError("minimum progress must be typed spacing.")
        if not isinstance(self.minimum_guide_radius, GuideRadius):
            raise InvalidCandidatePolicyError("minimum guide radius must be typed.")
        _positive_int(self.spatial_refinement_levels, "spatial refinement levels", InvalidCandidatePolicyError)
        _positive_int(self.radius_refinement_levels, "radius refinement levels", InvalidCandidatePolicyError)
        _positive_int(self.phase_count, "phase count", InvalidCandidatePolicyError)

    @classmethod
    def build(
        cls,
        *,
        spatial_resolution: Spacing,
        spatial_refinement_levels: int,
        radius_resolution: Spacing,
        radius_refinement_levels: int,
        phase_count: int,
        minimum_guide_radius: GuideRadius,
        minimum_progress: Spacing,
    ) -> Self:
        return cls(
            spatial_resolution,
            spatial_refinement_levels,
            radius_resolution,
            radius_refinement_levels,
            phase_count,
            minimum_guide_radius,
            minimum_progress,
        )

    def order_candidates(
        self,
        candidates: Sequence[CandidateT],
        *,
        key: Callable[[CandidateT], CandidateOrderKey],
    ) -> tuple[CandidateT, ...]:
        decorated: list[tuple[CandidateOrderKey, CandidateT]] = []
        for candidate in candidates:
            order_key = key(candidate)
            if not isinstance(order_key, CandidateOrderKey):
                raise InvalidCandidatePolicyError("candidate key function must return CandidateOrderKey.")
            decorated.append((order_key, candidate))

        identities = [order_key.canonical_identity for order_key, _ in decorated]
        if len(identities) != len(set(identities)):
            raise InvalidCandidatePolicyError("candidate canonical identities must be unique.")
        decorated.sort(
            key=lambda item: (
                -item[0].progress,
                -item[0].squared_radius,
                item[0].canonical_identity,
            )
        )
        return tuple(candidate for _, candidate in decorated)


@dataclass(frozen=True)
class MatSamplingPolicy:
    """Finite proposal-sampling bounds that determine MAT cursor identity."""

    station_spacing: Spacing
    max_sagitta: ChordBound
    max_refinement_depth: int

    def __post_init__(self) -> None:
        if type(self.station_spacing) is not Spacing:
            raise InvalidMatSamplingPolicyError("MAT station spacing must be exact typed spacing.")
        if type(self.max_sagitta) is not ChordBound:
            raise InvalidMatSamplingPolicyError("MAT maximum sagitta must be one exact typed chord bound.")
        if type(self.max_refinement_depth) is not int or self.max_refinement_depth <= 0:
            raise InvalidMatSamplingPolicyError("MAT maximum refinement depth must be one positive integer.")

    @classmethod
    def build(
        cls,
        *,
        station_spacing: Spacing,
        max_sagitta: ChordBound,
        max_refinement_depth: int,
    ) -> Self:
        """Build the complete MAT proposal-sampling policy.

        Args:
            station_spacing: Maximum reporting distance between stations.
            max_sagitta: Maximum reporting chord deviation on conics.
            max_refinement_depth: Fail-loud finite conic refinement cap.

        Returns:
            Immutable policy owning every cursor-affecting sampling bound.
        """
        return cls(
            station_spacing,
            max_sagitta,
            max_refinement_depth,
        )


@dataclass(frozen=True)
class NeckEffectiveCap:
    neck_class: int
    passage_state: PassageState
    cap: EngagementCap

    def __post_init__(self) -> None:
        if type(self.neck_class) is not int or self.neck_class < 0:
            raise InvalidNeckPolicyError("neck class must be a non-negative integer.")
        if not isinstance(self.passage_state, PassageState):
            raise InvalidNeckPolicyError("passage state must be typed.")
        if not isinstance(self.cap, EngagementCap):
            raise InvalidNeckPolicyError("effective cap must come from the native cap boundary.")


@dataclass(frozen=True)
class NeckPolicy:
    user_cap: EngagementCap
    squared_width_boundaries: tuple[SquaredMillimetre, ...]
    effective_caps: tuple[NeckEffectiveCap, ...]

    def __post_init__(self) -> None:
        try:
            boundaries = tuple(self.squared_width_boundaries)
            effective_caps = tuple(self.effective_caps)
        except TypeError:
            raise InvalidNeckPolicyError("neck boundaries and effective caps must be finite sequences.") from None
        object.__setattr__(self, "squared_width_boundaries", boundaries)
        object.__setattr__(self, "effective_caps", effective_caps)

        if not isinstance(self.user_cap, EngagementCap):
            raise InvalidNeckPolicyError("neck policy requires a native-produced user cap.")
        if any(type(boundary) is not Fraction or boundary < 0 for boundary in self.squared_width_boundaries):
            raise InvalidNeckPolicyError("squared-width boundaries must be exact non-negative fractions.")
        if tuple(sorted(set(self.squared_width_boundaries))) != self.squared_width_boundaries:
            raise InvalidNeckPolicyError("squared-width boundaries must be strictly increasing.")

        class_count = len(self.squared_width_boundaries) + 1
        expected_keys = {(neck_class, passage_state) for neck_class in range(class_count) for passage_state in ACTIVE_PASSAGE_STATES}
        observed_keys: set[tuple[int, PassageState]] = set()
        previous_key: tuple[int, int] | None = None
        for entry in self.effective_caps:
            if not isinstance(entry, NeckEffectiveCap):
                raise InvalidNeckPolicyError("effective-cap entries must be typed.")
            key = (entry.neck_class, entry.passage_state)
            ordered_key = (entry.neck_class, entry.passage_state.rank)
            if previous_key is not None and ordered_key <= previous_key:
                raise InvalidNeckPolicyError("effective-cap entries must use canonical order.")
            previous_key = ordered_key
            observed_keys.add(key)
            if not _stock_2.cap_chord_ratio_le(entry.cap.chord_ratio, self.user_cap.chord_ratio):
                raise InvalidNeckPolicyError("effective neck cap exceeds the user cap.")
        if observed_keys != expected_keys:
            raise InvalidNeckPolicyError("effective-cap mapping must cover every active neck class and passage state.")

    @classmethod
    def build(
        cls,
        *,
        user_cap: EngagementCap,
        squared_width_boundaries: Sequence[SquaredMillimetre],
        effective_caps: Mapping[tuple[int, PassageState], EngagementCap],
    ) -> Self:
        boundaries = tuple(squared_width_boundaries)
        if any(type(boundary) is not Fraction or boundary < 0 for boundary in boundaries):
            raise InvalidNeckPolicyError("squared-width boundaries must be exact non-negative fractions.")
        if len(set(boundaries)) != len(boundaries):
            raise InvalidNeckPolicyError("squared-width boundaries must be unique.")
        canonical_boundaries = tuple(sorted(boundaries))

        entries: list[NeckEffectiveCap] = []
        for key, cap in effective_caps.items():
            if not isinstance(key, tuple) or len(key) != 2 or type(key[0]) is not int or not isinstance(key[1], PassageState) or not isinstance(cap, EngagementCap):
                raise InvalidNeckPolicyError("effective-cap mapping keys and values must be typed.")
            entries.append(NeckEffectiveCap(key[0], key[1], cap))
        entries.sort(key=lambda entry: (entry.neck_class, entry.passage_state.rank))
        return cls(user_cap, canonical_boundaries, tuple(entries))

    def effective_cap(self, neck_class: int, passage_state: PassageState) -> EngagementCap:
        if type(neck_class) is not int or not isinstance(passage_state, PassageState):
            raise InvalidNeckPolicyError("neck class and passage state must be typed.")
        for entry in self.effective_caps:
            if entry.neck_class == neck_class and entry.passage_state is passage_state:
                return entry.cap
        raise InvalidNeckPolicyError("no effective cap exists for this terminal or unknown neck state.")


@dataclass(frozen=True)
class DepletionPolicy:
    chord_bound: ChordBound
    center_count_limit: int
    chord_bound_bytes: bytes

    def __post_init__(self) -> None:
        if not isinstance(self.chord_bound, ChordBound):
            raise InvalidDepletionPolicyError("depletion chord bound must be typed.")
        _positive_int(self.center_count_limit, "center-count limit", InvalidDepletionPolicyError)
        if self.chord_bound_bytes != self.chord_bound.exact_bytes:
            raise InvalidDepletionPolicyError("depletion policy must retain the exact chord-bound bytes.")

    @classmethod
    def build(cls, *, chord_bound: ChordBound, center_count_limit: int) -> Self:
        if not isinstance(chord_bound, ChordBound):
            raise InvalidDepletionPolicyError("depletion chord bound must be typed.")
        return cls(chord_bound, center_count_limit, chord_bound.exact_bytes)


@dataclass(frozen=True)
class TraversalPolicy:
    """Finite lookahead with component and branch order fixed by stable IDs."""

    forward_window: int

    def __post_init__(self) -> None:
        _positive_int(self.forward_window, "forward window", InvalidTraversalPolicyError)

    @classmethod
    def build(cls, *, forward_window: int) -> Self:
        return cls(forward_window)

    def order_components(self, component_ids: Sequence[ComponentId]) -> tuple[ComponentId, ...]:
        if any(type(component_id) is not bytes or not component_id for component_id in component_ids):
            raise InvalidTraversalPolicyError("component IDs must be nonempty bytes.")
        if len(component_ids) != len(set(component_ids)):
            raise InvalidTraversalPolicyError("component IDs must be unique.")
        return tuple(sorted(component_ids))

    def order_branches(self, branch_ids: Sequence[BranchId]) -> tuple[BranchId, ...]:
        if any(type(branch_id) is not bytes or not branch_id for branch_id in branch_ids):
            raise InvalidTraversalPolicyError("branch IDs must be nonempty bytes.")
        if len(branch_ids) != len(set(branch_ids)):
            raise InvalidTraversalPolicyError("branch IDs must be unique.")
        return tuple(sorted(branch_ids))


@dataclass(frozen=True)
class CutDirectionPolicy:
    intent: CutIntent

    def __post_init__(self) -> None:
        if not isinstance(self.intent, CutIntent):
            raise InvalidCutDirectionPolicyError("cut intent must be explicitly climb or conventional.")

    @classmethod
    def build(cls, intent: CutIntent) -> Self:
        return cls(intent)

    def circle_orientation(self, material_side: MaterialSide) -> CircleOrientation:
        if not isinstance(material_side, MaterialSide):
            raise InvalidCutDirectionPolicyError("material side must be explicit and typed.")
        if self.intent is CutIntent.CLIMB:
            if material_side is MaterialSide.INSIDE:
                return CircleOrientation.CLOCKWISE
            return CircleOrientation.COUNTERCLOCKWISE
        if self.intent is CutIntent.CONVENTIONAL:
            if material_side is MaterialSide.INSIDE:
                return CircleOrientation.COUNTERCLOCKWISE
            return CircleOrientation.CLOCKWISE
        raise InvalidCutDirectionPolicyError("cut intent and material side do not define an orientation.")
