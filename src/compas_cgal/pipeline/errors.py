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
