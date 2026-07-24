from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from enum import auto
from fractions import Fraction
from typing import Self

from compas_cgal import _stock_2
from compas_cgal.adaptive.errors import InvalidCandidatePolicyError
from compas_cgal.adaptive.errors import InvalidCutDirectionPolicyError
from compas_cgal.adaptive.errors import InvalidDepletionPolicyError
from compas_cgal.adaptive.errors import InvalidNeckPolicyError
from compas_cgal.adaptive.errors import InvalidTraversalPolicyError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre


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


_CIRCLE_ORIENTATION = {
    (CutIntent.CLIMB, MaterialSide.INSIDE): CircleOrientation.CLOCKWISE,
    (CutIntent.CLIMB, MaterialSide.OUTSIDE): CircleOrientation.COUNTERCLOCKWISE,
    (CutIntent.CONVENTIONAL, MaterialSide.INSIDE): CircleOrientation.COUNTERCLOCKWISE,
    (CutIntent.CONVENTIONAL, MaterialSide.OUTSIDE): CircleOrientation.CLOCKWISE,
}


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
        return _CIRCLE_ORIENTATION[(self.intent, material_side)]
