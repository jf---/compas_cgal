"""Exact witness-based refutation of one lateral segment motion.

Refutation and certification are not two accuracies of the same answer. They
are opposite quantifiers over the same exact predicate, and only one of them is
cheap.

- Refutation answers *some station of this motion violates the cap*. A single
  exact rational station that the native oracle classifies as buried or
  cap-exceeded is a complete proof, because the motion's maximum engagement is
  at least its engagement at that station. One station costs one exact cell
  classification.
- Certification answers *no station of this motion violates the cap*. That
  quantifies over a continuum, so it needs the full exact event partition, and
  the cost of discovering those events dominates everything else in candidate
  search. Measured on the Task 13F fixture: 5.1 s to 6.8 s for one full segment
  audit against 0.033 s for the whole nine-station ladder below.

The asymmetry is the entire point of this module. The cheap direction is
*complete when it fires* and *silent otherwise*: absence of a counterexample
among finitely many probed stations proves nothing whatsoever. This module
therefore cannot mint a certificate, and that prohibition is carried by the
types rather than by discipline.

- `CapRefutation` shares no supertype with `MotionWitness` or
  `SweptPrefixMotionWitness`; nominal typing alone rejects a substitution.
- Where a witness exposes `verdict: Literal["certified"]`, a refutation exposes
  `verdict: Literal["cap_exceeded"]`, so the two cannot unify even structurally.
- `refute_segment_cap` returns `CapRefutation | None`. The `None` arm is the
  "no counterexample found" outcome, and `None` is not assignable to any
  witness parameter either, so no call site can launder a silent probe into
  proof of safety. `tests/adaptive/typecheck/consumer_contract.py` pins both
  rejections under `mypy --strict`.

Every station outcome is one of three named values (`StationOutcome`). Only
`REFUTED` produces a `CapRefutation`; `NOT_REFUTED` and `UNKNOWN` are equally
inconclusive and both fall through to the full exact certifier, which remains
the sole authority that can accept a motion.
"""

from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Final
from typing import Literal
from typing import Self

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_integer
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidCapRefutationError
from compas_cgal.adaptive.errors import InvalidStationLadderError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.units import ToolRadius

CAP_REFUTATION_SCHEMA_VERSION: Final[bytes] = b"cap-refutation-schema-v1"

_DIGEST_SIZE: Final[int] = hashlib.sha256().digest_size
_BINARY64_SIZE: Final[int] = struct.calcsize(">d")


def _require_state_identities(
    stock_lineage_digest: bytes,
    stock_boundary_digest: bytes,
) -> None:
    """Require both stock identities to be exact SHA-256 digests.

    Args:
        stock_lineage_digest: SHA-256 identity of the observed depletion lineage.
        stock_boundary_digest: SHA-256 identity of the exact stock boundary.

    Raises:
        InvalidCapRefutationError: If either digest is not exact SHA-256 bytes.
    """
    if any(type(value) is not bytes or len(value) != _DIGEST_SIZE for value in (stock_lineage_digest, stock_boundary_digest)):
        raise InvalidCapRefutationError("cap refutation state identities must be exact SHA-256 digests.")


class StationOutcome(Enum):
    """Closed three-valued disposition of one probed exact station.

    Only `REFUTED` is conclusive. `NOT_REFUTED` and `UNKNOWN` are distinct
    reasons for the same non-result and are treated identically by every
    consumer: the full exact certifier still runs.
    """

    REFUTED = "refuted"
    NOT_REFUTED = "not-refuted"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class ExactStation:
    """One exact rational station `numerator / denominator` on a motion.

    The station is a parameter in the closed unit interval, not a coordinate.
    It is carried as a pair of Python integers because that is exactly what the
    native oracle consumes; no float ever names a station.

    Attributes:
        numerator: Nonnegative integer, no greater than `denominator`.
        denominator: Positive integer.
    """

    numerator: int
    denominator: int

    def __post_init__(self) -> None:
        if type(self) is not ExactStation:
            raise InvalidStationLadderError("exact station must use the exact owned type.")
        if type(self.numerator) is not int or type(self.denominator) is not int:
            raise InvalidStationLadderError("exact station requires integer numerator and denominator.")
        if self.denominator < 1:
            raise InvalidStationLadderError("exact station denominator must be positive.")
        if self.numerator < 0 or self.numerator > self.denominator:
            raise InvalidStationLadderError("exact station must lie in the closed unit interval.")

    @classmethod
    def build(cls, numerator: int, denominator: int) -> Self:
        """Build one validated station in the closed unit interval.

        Args:
            numerator: Nonnegative integer, no greater than `denominator`.
            denominator: Positive integer.

        Returns:
            Immutable exact station.

        Raises:
            InvalidStationLadderError: If the pair is not one exact station of
                the closed unit interval.
        """
        return cls(numerator, denominator)

    @property
    def parameter(self) -> Fraction:
        """Return the exact rational parameter this station names."""
        return Fraction(self.numerator, self.denominator)

    @property
    def canonical_bytes(self) -> bytes:
        """Return the versioned canonical record of this station."""
        return encode_tagged_union(
            b"exact-station-v1",
            encode_component_map(
                {
                    b"denominator": encode_integer(self.denominator),
                    b"numerator": encode_integer(self.numerator),
                }
            ),
        )


def _dyadic_station_ladder(depth: int) -> tuple[ExactStation, ...]:
    """Enumerate the closed unit interval by increasing dyadic refinement.

    Level 0 is the two endpoints, terminal station first: a link segment starts
    inside stock the previous operation already cleared and advances into fresh
    material, so its terminal station is the likeliest counterexample and is
    probed first. Level `m` then adds the odd multiples of `2**-m`, which are
    exactly the stations no coarser level has already probed.

    Args:
        depth: Deepest dyadic level to enumerate.

    Returns:
        Station ladder ordered coarsest level first, ascending within a level.

    Raises:
        InvalidStationLadderError: If `depth` is not one nonnegative integer.
    """
    if type(depth) is not int or depth < 0:
        raise InvalidStationLadderError("station ladder depth must be one nonnegative exact integer.")
    stations: list[ExactStation] = [
        ExactStation.build(1, 1),
        ExactStation.build(0, 1),
    ]
    for level in range(1, depth + 1):
        denominator = 2**level
        stations.extend(ExactStation.build(numerator, denominator) for numerator in range(1, denominator, 2))
    return tuple(stations)


# A violation occupying a fraction f of the motion is certain to be hit once the
# ladder spacing 2**-depth drops below f, so depth 3 refutes every violation
# spanning more than an eighth of the motion. That costs nine exact station
# evaluations -- measured 0.033 s on the Task 13F fixture, against 5.1 s to 6.8 s
# for the full event partition it replaces, a factor of roughly 160. Deeper
# levels double the probe budget to catch violations narrower than an eighth of
# a link, which candidate search does not produce; and no depth can ever be
# "enough", because the full audit remains the authority for every motion this
# ladder fails to refute.
REFUTATION_LADDER_DEPTH: Final[int] = 3
REFUTATION_STATION_LADDER: Final[tuple[ExactStation, ...]] = _dyadic_station_ladder(
    REFUTATION_LADDER_DEPTH,
)


@dataclass(frozen=True)
class CapRefutation:
    """Exact counterexample proving one segment motion exceeds its cap.

    A refutation is a positive proof of violation and is deliberately not a
    certificate: it names the station at which the motion was proved to fail,
    and it can be independently re-checked by evaluating that one station
    against the same stock, motion, radius, and cap.

    Attributes:
        motion: Exact segment motion the counterexample belongs to.
        tool_radius: Typed cutter radius the station was evaluated with.
        effective_cap_bytes: Canonical binary64 effective-cap surrogate.
        stock_lineage_digest: SHA-256 identity of the observed depletion lineage.
        stock_boundary_digest: SHA-256 identity of the exact stock boundary.
        witness_station: Exact rational station that exceeds the cap.
        verdict: Exact refuted verdict, the only constructible value.
    """

    motion: ExactSegmentMotion
    tool_radius: ToolRadius
    effective_cap_bytes: bytes
    stock_lineage_digest: bytes
    stock_boundary_digest: bytes
    witness_station: ExactStation
    verdict: Literal["cap_exceeded"]

    def __post_init__(self) -> None:
        if type(self) is not CapRefutation:
            raise InvalidCapRefutationError("cap refutation must use the exact owned type.")
        if type(self.motion) is not ExactSegmentMotion:
            raise InvalidCapRefutationError("cap refutation requires one exact segment motion.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidCapRefutationError("cap refutation requires one exact tool radius.")
        if type(self.effective_cap_bytes) is not bytes or len(self.effective_cap_bytes) != _BINARY64_SIZE:
            raise InvalidCapRefutationError("cap refutation requires an exact binary64 cap surrogate.")
        _require_state_identities(
            self.stock_lineage_digest,
            self.stock_boundary_digest,
        )
        if type(self.witness_station) is not ExactStation:
            raise InvalidCapRefutationError("cap refutation requires one exact witness station.")
        if self.verdict != "cap_exceeded":
            raise InvalidCapRefutationError("cap refutation may represent only a proved violation.")

    @property
    def canonical_bytes(self) -> bytes:
        """Return the complete versioned counterexample record.

        Returns:
            Canonical CCAN bytes binding the motion, stock identity, cap, and
            the exact witness station.
        """
        return encode_tagged_union(
            CAP_REFUTATION_SCHEMA_VERSION,
            encode_component_map(
                {
                    b"effective-cap": encode_bytes(self.effective_cap_bytes),
                    b"motion": canonical_task1_bytes(self.motion),
                    b"stock-boundary-digest": encode_bytes(self.stock_boundary_digest),
                    b"stock-lineage-digest": encode_bytes(self.stock_lineage_digest),
                    b"tool-radius": encode_tagged_union(
                        b"tool-radius-mm-v1",
                        encode_binary64(float(self.tool_radius.value)),
                    ),
                    b"verdict": encode_tagged_union(
                        b"motion-verdict-cap-exceeded-v1",
                        b"",
                    ),
                    b"witness-station": self.witness_station.canonical_bytes,
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        """Return the SHA-256 identity of `canonical_bytes`."""
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


def classify_segment_station(
    *,
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
    station: ExactStation,
) -> StationOutcome:
    """Classify one exact station of a segment motion against its cap.

    The native predicate is exact at the station, so `REFUTED` is a proof and
    never an estimate. Every native substrate failure -- an unresolved station
    disposition, an unextractable boundary, an input the native source rejects
    -- yields `UNKNOWN` rather than propagating: an unposable or unresolved
    probe refutes nothing, and the full certifier that runs next owns those
    error models and raises them with its own contract and message.

    Args:
        stock: Immutable native stock snapshot the probe reads.
        motion: Exact segment motion carrying the station.
        tool_radius: Typed cutter radius.
        effective_cap: Exact policy-derived engagement cap.
        station: Exact rational station in the closed unit interval.

    Returns:
        `REFUTED` when the station exceeds the cap, `NOT_REFUTED` when it
        exactly does not, `UNKNOWN` when the exact probe has no disposition.

    Raises:
        InvalidCapRefutationError: If any argument is not one exact owned value.
    """
    if type(stock) is not _stock_2.Stock2:
        raise InvalidCapRefutationError("station probe requires one owned native Stock2 snapshot.")
    if type(motion) is not ExactSegmentMotion:
        raise InvalidCapRefutationError("station probe requires one exact segment motion.")
    if type(tool_radius) is not ToolRadius or type(effective_cap) is not EngagementCap:
        raise InvalidCapRefutationError("station probe requires an exact tool radius and engagement cap.")
    if type(station) is not ExactStation:
        raise InvalidCapRefutationError("station probe requires one exact station.")
    try:
        exceeded = _continuous_tea_2.segment_station_cap_exceeded_exact(
            stock,
            motion.start.x,
            motion.start.y,
            motion.end.x,
            motion.end.y,
            station.numerator,
            station.denominator,
            tool_radius.value,
            effective_cap.chord_ratio,
        )
    except (
        _continuous_tea_2.BoundaryExtractionError,
        _continuous_tea_2.EventSubstrateError,
    ):
        return StationOutcome.UNKNOWN
    if type(exceeded) is not bool:
        raise InvalidCapRefutationError("native station predicate returned a non-boolean disposition.")
    return StationOutcome.REFUTED if exceeded else StationOutcome.NOT_REFUTED


def refute_segment_cap(
    *,
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
    stock_lineage_digest: bytes,
    stock_boundary_digest: bytes,
) -> CapRefutation | None:
    """Search `REFUTATION_STATION_LADDER` for an exact cap counterexample.

    The ladder is fixed, deterministic, and float-free: the same stock, motion,
    radius, and cap always probe the same stations in the same order, and the
    search stops at the first `REFUTED` station.

    A `None` result means only that no probed station was refuted. It is not a
    certificate, not a weak certificate, and not evidence of safety; the caller
    must still run the full exact certification.

    Args:
        stock: Immutable native stock snapshot the probe reads.
        motion: Exact segment motion to refute.
        tool_radius: Typed cutter radius.
        effective_cap: Exact policy-derived engagement cap.
        stock_lineage_digest: SHA-256 identity of the observed depletion lineage.
        stock_boundary_digest: SHA-256 identity of the exact stock boundary.

    Returns:
        The counterexample at the first refuted station, or `None` when the
        ladder is inconclusive.

    Raises:
        InvalidCapRefutationError: If any argument is not one exact owned value.
    """
    # Validated before the ladder rather than inside `CapRefutation`, so a
    # malformed stock identity fails on every call instead of only on the calls
    # that happen to find a counterexample.
    _require_state_identities(stock_lineage_digest, stock_boundary_digest)
    for station in REFUTATION_STATION_LADDER:
        outcome = classify_segment_station(
            stock=stock,
            motion=motion,
            tool_radius=tool_radius,
            effective_cap=effective_cap,
            station=station,
        )
        if outcome is StationOutcome.REFUTED:
            return CapRefutation(
                motion,
                tool_radius,
                effective_cap.chord_ratio_bytes,
                stock_lineage_digest,
                stock_boundary_digest,
                station,
                "cap_exceeded",
            )
    return None
