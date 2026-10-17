"""One named exception per failure mode of the pipeline (Ruling 1.1-C8)."""

from __future__ import annotations


class EmptyPipelineGuideError(RuntimeError):
    """The straight-skeleton guide yielded no machining station: the pocket admits no gouge-free loop at this tool."""


class GuideTruncatedError(RuntimeError):
    """`max_passes` left skeleton chains unread, so the guide does not cover the pocket."""


class GuideStationOrderError(RuntimeError):
    """The native guide interleaved a chain's stations with another chain's (broken invariant)."""


class ImplausibleClearanceSlopeError(ValueError):
    """A guide chain's radius changes faster along its own arc length than a 1-Lipschitz clearance can."""


class EmptyLadderWindowError(ValueError):
    """An advance was asked for with no candidate station beyond the origin."""


class SubFloorStationError(ValueError):
    """A loop was requested at a station that does not clear the degeneracy floor."""


class InvalidEngagementCapError(ValueError):
    """The engagement cap is not an angle in ``(0, pi]``."""


class LedgerInvariantError(RuntimeError):
    """A forced-decision ledger does not hold exactly one first-cut entry."""


class EmptyNarrowPathError(RuntimeError):
    """The guide produced chains, but no station on any of them clears the degeneracy floor."""


class InvalidClearanceHeightError(ValueError):
    """The clearance plane is not strictly above the cut plane."""


class ZeroLengthLinkError(ValueError):
    """A link motion was asked for between two coincident points; the caller emits nothing instead."""


class InvalidCapChordRatioError(ValueError):
    """A cap chord ratio outside ``(0, 4]``, the surrogates of caps in ``(0, pi]``."""


class CoincidentStationsError(ValueError):
    """Two stations an advance connects share their centre exactly: no advance direction exists."""


class MixedTurnChainError(RuntimeError):
    """Two consecutive stations of one chain turn opposite ways, so no co-rotating link joins them (broken invariant)."""


class FirstCutReEntryError(ValueError):
    """A re-entry was asked to record the first cut; the first cut is its own stage."""


class FirstCutAfterEmissionError(RuntimeError):
    """The first cut was asked for after operations were emitted; it starts at the clearance plane with nothing emitted."""


class InvalidEntryRiseChordRatioError(ValueError):
    """An entry-rise chord ratio outside ``(0, 4]``, the surrogates of rise limits in ``(0, pi]``."""


class MissingEntryRiseLimitError(RuntimeError):
    """A pencil junction was asked of an admissibility that declares no entry-rise limit."""


class RampOutStuckError(RuntimeError):
    """No pencil member above the current one passes the coverage bound, the junction and the cap: the helical entry cannot ramp out to its root."""


class HoledDomainPencilError(RuntimeError):
    """A pencil's root disk is not contained in the design pocket (an island inside it): its members have no exact containment decision."""


class MissingHelicalEntryError(RuntimeError):
    """The first cut was asked of a path that declares no helical entry (helix pitch surrogate and design pocket)."""


class UnhandledOperationGeometryError(TypeError):
    """An operation's geometry is not one the emitter's exhaustive dispatch handles."""


class RampOutRecordMismatchError(RuntimeError):
    """A helical-entry record contradicts the operation stream it describes."""


class RootLoopNotEmittedError(RuntimeError):
    """The first cut's root loop was not emitted after its ramp-out (broken invariant: the root always removes the lune its pencil leaves)."""
