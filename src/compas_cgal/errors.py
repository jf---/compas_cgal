"""Error model for the shared unit and motion vocabulary.

These are the failure modes of `compas_cgal.units` and `compas_cgal.motion`,
which the toolpath generators, the benchmark corpus and the certification lane
all depend on. Errors raised only by the certification lane live in
`compas_cgal.adaptive.errors`.
"""


class InvalidUnitValueError(ValueError):
    """A typed unit or framed coordinate violates its numeric invariant."""


class InvalidEngagementCapError(ValueError):
    """An engagement cap cannot be represented by the native exact boundary."""


class DegenerateSegmentMotionError(ValueError):
    """An exact segment motion has no progress."""


class DegenerateCircleMotionError(ValueError):
    """An exact circle motion has no nonzero phase vector."""
