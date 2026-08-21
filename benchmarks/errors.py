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


class UnreplayableOperationError(BenchmarkError):
    """A toolpath operation lies outside the cut-plane depletion model (e.g. a 3D ramp)."""


class UnsampleableMotionError(BenchmarkError):
    """A cut motion's geometry carries no cutter-centre path to sample positions from."""


class MalformedRecordError(BenchmarkError):
    """A serialised measurement does not carry exactly the columns the schema declares."""
