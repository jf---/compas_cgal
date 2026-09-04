"""Post-qualification refusal boundary for a complete Held Figure 5 characterization."""

from __future__ import annotations

from dataclasses import dataclass

from typing_extensions import Self

from benchmarks.errors import HeldPathNotEligibleForPostQualificationError
from benchmarks.held_path_evidence import HeldFigure5Characterization
from benchmarks.held_path_snapshot import HeldOperationSnapshot
from benchmarks.quality_observations import CRITERION_NAMES

_OPEN_OUTCOMES = frozenset(
    {
        "failure_observed",
        "criterion_violated",
        "outside_declared_tolerance",
    }
)
_ASSESSMENT_FIELDS = (
    "uncut_fraction",
    "gouging_motions",
    "unsafe_rapids",
    "continuity_breaks",
    "zero_length_motions",
    "degenerate_loops",
    "redundant_operations",
    "cap_exceedances",
    "slotting_motions",
    "max_engagement_step",
    "max_loop_radius_step",
    "tangent_breaks",
)


def post_qualification_failures(characterization: HeldFigure5Characterization) -> tuple[str, ...]:
    """Return every open post-qualification condition in canonical order."""
    failures: list[str] = []
    engagement = characterization.engagement
    if engagement.demonstrated_exceeded_count:
        failures.append(f"demonstrated engagement exceedances: {engagement.demonstrated_exceeded_count}")
    if engagement.unresolved_count:
        failures.append(f"unresolved TEA-audited operations: {engagement.unresolved_count}")
    criteria = {criterion.name: criterion for criterion in (getattr(characterization.assessment, field) for field in _ASSESSMENT_FIELDS)}
    for name in CRITERION_NAMES:
        criterion = criteria[name]
        if criterion.outcome in _OPEN_OUTCOMES:
            failures.append(f"{criterion.name}: measured={criterion.measured}, required={criterion.required}, outcome={criterion.outcome}")
    return tuple(failures)


@dataclass(frozen=True, init=False)
class HeldPostQualificationCandidate:
    """A complete characterization whose Phase 1 conditions are all closed."""

    characterization: HeldFigure5Characterization
    snapshot: tuple[HeldOperationSnapshot, ...]

    def __init__(self) -> None:
        raise TypeError("HeldPostQualificationCandidate must be created with HeldPostQualificationCandidate.build().")

    @classmethod
    def build(cls, characterization: HeldFigure5Characterization) -> Self:
        """Build only when the sole public failure collector reports no open gate."""
        failures = post_qualification_failures(characterization)
        if failures:
            engagement = characterization.engagement
            joined = "; ".join(failures)
            raise HeldPathNotEligibleForPostQualificationError(
                f"Held Figure 5 is not eligible for post qualification: {joined}. "
                "Engagement counts: "
                f"certified={engagement.certified_count}, "
                f"demonstrated_exceeded={engagement.demonstrated_exceeded_count}, "
                f"unresolved={engagement.unresolved_count}."
            )
        candidate = object.__new__(cls)
        object.__setattr__(candidate, "characterization", characterization)
        object.__setattr__(candidate, "snapshot", characterization.snapshot)
        return candidate


def require_post_qualification_candidate(
    characterization: HeldFigure5Characterization,
) -> HeldPostQualificationCandidate:
    """Return the sole factory-built post-qualification candidate."""
    return HeldPostQualificationCandidate.build(characterization)
