import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations

import pytest

import compas_cgal.adaptive.policy as policy_module
from compas_cgal.adaptive.errors import InvalidCandidatePolicyError
from compas_cgal.adaptive.errors import InvalidCutDirectionPolicyError
from compas_cgal.adaptive.errors import InvalidDepletionPolicyError
from compas_cgal.adaptive.errors import InvalidNeckPolicyError
from compas_cgal.adaptive.errors import InvalidTraversalPolicyError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.policy import BranchId
from compas_cgal.adaptive.policy import CandidateOrderKey
from compas_cgal.adaptive.policy import CandidatePolicy
from compas_cgal.adaptive.policy import CircleOrientation
from compas_cgal.adaptive.policy import ComponentId
from compas_cgal.adaptive.policy import CutDirectionPolicy
from compas_cgal.adaptive.policy import CutIntent
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.policy import MaterialSide
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState
from compas_cgal.adaptive.policy import TraversalPolicy
from compas_cgal.adaptive.units import ChordBound
from compas_cgal.adaptive.units import ExactMillimetre
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


@dataclass(frozen=True)
class _Candidate:
    name: str
    order_key: CandidateOrderKey


def test_candidate_policy_executes_canonical_tie_break_independent_of_input_order() -> None:
    policy = _candidate_policy()
    candidates = (
        _Candidate(
            "least-progress",
            CandidateOrderKey.build(
                progress=ExactMillimetre(Fraction(1)),
                squared_radius=SquaredMillimetre(Fraction(100)),
                canonical_identity=b"d",
            ),
        ),
        _Candidate(
            "smaller-radius",
            CandidateOrderKey.build(
                progress=ExactMillimetre(Fraction(2)),
                squared_radius=SquaredMillimetre(Fraction(1)),
                canonical_identity=b"c",
            ),
        ),
        _Candidate(
            "lexicographic-second",
            CandidateOrderKey.build(
                progress=ExactMillimetre(Fraction(2)),
                squared_radius=SquaredMillimetre(Fraction(4)),
                canonical_identity=b"b",
            ),
        ),
        _Candidate(
            "winner",
            CandidateOrderKey.build(
                progress=ExactMillimetre(Fraction(2)),
                squared_radius=SquaredMillimetre(Fraction(4)),
                canonical_identity=b"a",
            ),
        ),
    )
    expected = ("winner", "lexicographic-second", "smaller-radius", "least-progress")

    for permutation in permutations(candidates):
        ordered = policy.order_candidates(permutation, key=lambda candidate: candidate.order_key)
        assert tuple(candidate.name for candidate in ordered) == expected


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


def test_raw_neck_policy_does_not_retain_mutable_sequences() -> None:
    cap = EngagementCap.build(math.radians(90.0))
    built = NeckPolicy.build(
        user_cap=cap,
        squared_width_boundaries=(SquaredMillimetre(Fraction(4)),),
        effective_caps=_complete_caps(2, cap),
    )
    mutable_boundaries = list(built.squared_width_boundaries)
    mutable_caps = list(built.effective_caps)
    raw = NeckPolicy(cap, mutable_boundaries, mutable_caps)  # type: ignore[arg-type]

    mutable_boundaries.clear()
    mutable_caps.clear()

    assert raw.squared_width_boundaries == built.squared_width_boundaries
    assert raw.effective_caps == built.effective_caps


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


def test_traversal_policy_executes_stable_id_order_for_every_permutation() -> None:
    policy = TraversalPolicy.build(forward_window=16)
    component_ids = (ComponentId(b"component-c"), ComponentId(b"component-a"), ComponentId(b"component-b"))
    branch_ids = (BranchId(b"branch-c"), BranchId(b"branch-a"), BranchId(b"branch-b"))

    for permutation in permutations(component_ids):
        assert policy.order_components(permutation) == tuple(sorted(component_ids))
    for permutation in permutations(branch_ids):
        assert policy.order_branches(permutation) == tuple(sorted(branch_ids))


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


def test_cut_direction_has_no_mutable_global_orientation_table() -> None:
    assert not hasattr(policy_module, "_CIRCLE_ORIENTATION")
