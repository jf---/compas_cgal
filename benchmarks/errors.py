from __future__ import annotations


class BenchmarkError(Exception):
    """Base class for every benchmark-corpus failure mode."""


class NonPositiveToolError(BenchmarkError):
    """A tool diameter was zero, negative, or NaN."""


class InvalidCapError(BenchmarkError):
    """An engagement cap fell outside the exact kernel's contract of (0, 180] degrees."""


class DegeneratePocketError(BenchmarkError):
    """A generated pocket has zero area, fewer than three vertices, or no room for the tool."""


class PocketNotSimpleError(BenchmarkError):
    """A pocket boundary self-intersects, so it is not a valid general polygon."""
