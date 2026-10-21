"""Named failure modes of the content-addressed Held reconstruction cache."""

from __future__ import annotations

from benchmarks.errors import BenchmarkError


class ReconstructionCacheError(BenchmarkError):
    """Base class: the cache cannot vouch for an entry, so nothing is used."""


class InvalidReconstructionBuildIdentityError(ReconstructionCacheError):
    """A build identity component is not well formed (digest length, empty label)."""


class MalformedReconstructionEncodingError(ReconstructionCacheError):
    """A serialized reconstruction does not decode to a valid exact reconstruction."""


class CorruptedReconstructionCacheEntryError(ReconstructionCacheError):
    """An entry is unreadable, fails its integrity digest, or does not decode."""


class UnsupportedReconstructionCacheFormatError(ReconstructionCacheError):
    """An entry was written in a cache format this code does not read."""


class ReconstructionCacheIdentityMismatchError(ReconstructionCacheError):
    """An entry's embedded build identity differs from the identity it is filed under."""


class ReconstructionCacheWriteError(ReconstructionCacheError):
    """An entry could not be written atomically."""
