import math
from fractions import Fraction

import pytest

from compas_cgal.adaptive.errors import InvalidCandidatePolicyError
from compas_cgal.adaptive.errors import InvalidCutDirectionPolicyError
from compas_cgal.adaptive.errors import InvalidDepletionPolicyError
from compas_cgal.adaptive.errors import InvalidNeckPolicyError
from compas_cgal.adaptive.errors import InvalidTraversalPolicyError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import GuideRadius
from compas_cgal.adaptive.units import Spacing
from compas_cgal.adaptive.units import SquaredMillimetre


def _candidate_policy() -> CandidatePolicy:
    return CandidatePolicy.build(
        spatial_resolution=Spacing.build(0.25),
        spatial_refinement_levels=3,
        radius_resolution=Spacing.build(0.125),
        radius_refinement_levels=4,
        phase_count=8,
        minimum_guide_radius=GuideRadius.build(0.0625),
        minimum_progress=Spacing.build(0.05),
    )


def test_candidate_policy_owns_finite_lattice_and_canonical_tie_break() -> None:
    policy = _candidate_policy()

    assert policy.spatial_refinement_levels == 3
    assert policy.radius_refinement_levels == 4
    assert policy.phase_count == 8
    assert policy.minimum_guide_radius == GuideRadius.build(0.0625)
    assert not hasattr(policy, "tie_break")
    assert not hasattr(policy, "chord_bound")
    assert not hasattr(policy, "forward_window")


@pytest.mark.parametrize(
    "changes",
    [
        {"spatial_refinement_levels": 0},
        {"radius_refinement_levels": -1},
        {"phase_count": 0},
        {"spatial_refinement_levels": math.inf},
    ],
)
def test_candidate_policy_rejects_unbounded_lattice(changes: dict[str, object]) -> None:
    values: dict[str, object] = {
        "spatial_resolution": Spacing.build(0.25),
        "spatial_refinement_levels": 3,
        "radius_resolution": Spacing.build(0.125),
        "radius_refinement_levels": 4,
        "phase_count": 8,
        "minimum_guide_radius": GuideRadius.build(0.0625),
        "minimum_progress": Spacing.build(0.05),
    }
    values.update(changes)

    with pytest.raises(InvalidCandidatePolicyError):
        CandidatePolicy.build(**values)  # type: ignore[arg-type]


def _complete_caps(
    class_count: int,
    cap: EngagementCap,
) -> dict[tuple[int, PassageState], EngagementCap]:
    return {
        (neck_class, state): cap
        for neck_class in range(class_count)
        for state in (
            PassageState.UNVISITED,
            PassageState.FIRST_PASS_COMPLETE,
            PassageState.SECOND_PASS_COMPLETE,
        )
    }


def test_neck_policy_canonicalizes_exact_boundaries_and_cap_map() -> None:
    user_cap = EngagementCap.build(math.radians(120.0))
    effective = EngagementCap.build(math.radians(90.0))
    caps = _complete_caps(4, effective)
    reversed_caps = dict(reversed(tuple(caps.items())))

    policy = NeckPolicy.build(
        user_cap=user_cap,
        squared_width_boundaries=(
            SquaredMillimetre(Fraction(9)),
            SquaredMillimetre(Fraction(1)),
            SquaredMillimetre(Fraction(4)),
        ),
        effective_caps=reversed_caps,
    )

    assert policy.squared_width_boundaries == (
        SquaredMillimetre(Fraction(1)),
        SquaredMillimetre(Fraction(4)),
        SquaredMillimetre(Fraction(9)),
    )
    assert policy.effective_cap(2, PassageState.FIRST_PASS_COMPLETE) == effective
    assert tuple((entry.neck_class, entry.passage_state) for entry in policy.effective_caps) == tuple(sorted(caps, key=lambda key: (key[0], key[1].rank)))
    with pytest.raises(InvalidNeckPolicyError):
        policy.effective_cap(2, PassageState.TERMINAL)


def test_neck_policy_rejects_incomplete_or_over_user_cap_map() -> None:
    user_cap = EngagementCap.build(math.radians(90.0))
    valid = EngagementCap.build(math.radians(60.0))
    caps = _complete_caps(2, valid)
    caps.pop((1, PassageState.SECOND_PASS_COMPLETE))

    with pytest.raises(InvalidNeckPolicyError):
        NeckPolicy.build(
            user_cap=user_cap,
            squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
            effective_caps=caps,
        )

    caps = _complete_caps(2, valid)
    caps[(1, PassageState.UNVISITED)] = EngagementCap.build(math.radians(100.0))
    with pytest.raises(InvalidNeckPolicyError):
        NeckPolicy.build(
            user_cap=user_cap,
            squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
            effective_caps=caps,
        )


def test_neck_policy_rejects_reporting_or_noncanonical_boundaries() -> None:
    cap = EngagementCap.build(math.radians(90.0))

    with pytest.raises(InvalidNeckPolicyError):
        NeckPolicy.build(
            user_cap=cap,
            squared_width_boundaries=(1.0,),  # type: ignore[arg-type]
            effective_caps=_complete_caps(2, cap),
        )
    with pytest.raises(InvalidNeckPolicyError):
        NeckPolicy.build(
            user_cap=cap,
            squared_width_boundaries=(
                SquaredMillimetre(Fraction(1)),
                SquaredMillimetre(Fraction(1)),
            ),
            effective_caps=_complete_caps(3, cap),
        )


def test_depletion_policy_owns_exact_chord_bytes_and_finite_limit() -> None:
    policy = DepletionPolicy.build(
        chord_bound=ChordBound.build(0.03125),
        center_count_limit=4096,
    )

    assert policy.chord_bound_bytes == policy.chord_bound.exact_bytes
    assert policy.center_count_limit == 4096

    with pytest.raises(InvalidDepletionPolicyError):
        DepletionPolicy.build(chord_bound=ChordBound.build(0.03125), center_count_limit=0)


def test_traversal_policy_requires_typed_order_and_finite_window() -> None:
    policy = TraversalPolicy.build(forward_window=16)

    assert policy.forward_window == 16
    with pytest.raises(InvalidTraversalPolicyError):
        TraversalPolicy.build(forward_window=math.inf)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("intent", "material_side", "orientation"),
    [
        (CutIntent.CLIMB, MaterialSide.INSIDE, CircleOrientation.CLOCKWISE),
        (CutIntent.CLIMB, MaterialSide.OUTSIDE, CircleOrientation.COUNTERCLOCKWISE),
        (CutIntent.CONVENTIONAL, MaterialSide.INSIDE, CircleOrientation.COUNTERCLOCKWISE),
        (CutIntent.CONVENTIONAL, MaterialSide.OUTSIDE, CircleOrientation.CLOCKWISE),
    ],
)
def test_cut_direction_maps_explicit_intent_and_material_side(
    intent: CutIntent,
    material_side: MaterialSide,
    orientation: CircleOrientation,
) -> None:
    policy = CutDirectionPolicy.build(intent)

    assert policy.circle_orientation(material_side) is orientation


def test_cut_direction_rejects_raw_strings_and_has_no_default_side() -> None:
    with pytest.raises(InvalidCutDirectionPolicyError):
        CutDirectionPolicy.build("climb")  # type: ignore[arg-type]

    policy = CutDirectionPolicy.build(CutIntent.CLIMB)
    with pytest.raises(InvalidCutDirectionPolicyError):
        policy.circle_orientation("inside")  # type: ignore[arg-type]
