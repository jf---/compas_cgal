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


class InvalidSideCountError(BenchmarkError):
    """A family was asked for a polygon with fewer than three sides."""


class InvalidArcRatioError(BenchmarkError):
    """An arc fraction fell outside the closed unit interval it is defined on."""


class ImpassableNeckError(BenchmarkError):
    """A neck was narrower than the tool, so no toolpath can reach through it."""


class UnpinchedChannelError(BenchmarkError):
    """A neck was as wide as the pocket it pinches, so the instance carries no neck."""


class MissingSweepParameterError(BenchmarkError):
    """A record handed to a sweep analysis does not carry that sweep's parameter."""


class InvalidDecimalsError(BenchmarkError):
    """A coordinate precision was negative, so it names no rounding."""


class InvalidIslandCountError(BenchmarkError):
    """An island grid was asked for a non-positive number of rows or columns."""


class CrowdedIslandGridError(BenchmarkError):
    """An island grid leaves channels too narrow for the tool to machine."""


class UnmeasurableOperationLengthError(BenchmarkError):
    """A toolpath operation's geometry is not a primitive whose length is defined."""


class UnknownCorpusError(BenchmarkError):
    """A corpus was requested by a name no family answers to."""


class MissingExternalDirectoryError(BenchmarkError):
    """The external corpus was requested without the operator-supplied profile directory."""
