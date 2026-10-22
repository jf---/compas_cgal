"""Named failure modes of the Figure 5 engagement refinement, one exception each."""

from __future__ import annotations

from benchmarks.errors import BenchmarkError


class Figure5EngagementResolutionError(BenchmarkError):
    """A source interval cannot be refined within the declared budget."""


class Figure5EngagementFragmentationError(Figure5EngagementResolutionError):
    """The exact predecessor engagement of one transition has no real swept-disk intersection.

    Wraps ``StandardPlacementFragmentationError`` (user ruling 2026-10-21): the
    transition cannot be refined, which is a refusal of the interval, not an
    escape from the refinement's error model.
    """


class Figure5EngagementCosineRangeError(Figure5EngagementResolutionError):
    """The exact predecessor engagement cosine of one transition left [-1, 1].

    Wraps ``ExactEngagementCosineRangeError`` (user ruling 2026-10-21), a
    refusal of the interval like ``Figure5EngagementFragmentationError``.
    """


class Figure5EngagementPrefixError(Figure5EngagementResolutionError):
    """A continuation prefix is not the emission of whole earlier sources."""


class Figure5EngagementTraceError(BenchmarkError):
    """A refinement handed in as a continuation's prefix or reference has no trace, or one made under other settings.

    A caller contract breach, not a refusal: deliberately not a
    ``Figure5EngagementResolutionError``, so no repair records it as one.
    """


class ConnectorEngagementResolutionError(Figure5EngagementResolutionError):
    """A connector is still without a tier-1 or tier-1b engagement certificate when its repair runs out.

    The repair runs out when the depth budget is spent, or when the midpoint
    circle it needs is refused (at a right-angle convex vertex the corner
    limit shrinks midpoint circles to nothing, and the single-annulus lemma
    cannot certify a move away from a vanishing circle).

    ``connector`` indexes the emitted circle the connector leaves (the
    post-gate's connector index); ``seam`` marks the closing connector back to
    circles[0]. A subclass of the refinement's refusal, so its callers'
    error model is unchanged.
    """

    def __init__(self, connector: int, seam: bool, detail: str) -> None:
        self.connector = connector
        self.seam = seam
        kind = "seam connector" if seam else "connector"
        super().__init__(f"Repair exhausted with {kind} {connector} uncertified: {detail}.")


class ConvexCornerEngagementError(ConnectorEngagementResolutionError):
    """A connector crosses a convex boundary vertex whose exterior angle pi - alpha exceeds the cap.

    Lemma (docs/held_connector_engagement.md, "Convex corners: an engagement
    floor no refinement removes"): the first arrival of the tool centre at a
    convex vertex of interior angle alpha engages at least pi - alpha of rim,
    whatever guide circles were cut before (inside its tangent wedge), so no
    midpoint refinement can certify the connector. Decided exactly
    (``_ConnectorGate.sharp_vertex``) before any midpoint toward the vertex
    is tried. ``vertex`` indexes the boundary vertex.
    """

    def __init__(self, connector: int, seam: bool, vertex: int) -> None:
        self.vertex = vertex
        super().__init__(connector, seam, f"crosses convex vertex {vertex}, whose exterior angle exceeds the cap (convex-corner lemma)")
