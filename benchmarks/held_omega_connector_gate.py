"""Exact engagement decision of ONE Ω* connector, as the post-gate decides it, for the refinement's acceptance.

The post-gate (:func:`benchmarks.held_omega_connector_engagement.assess_omega_connectors`)
decides connector ``i`` from circles ``0..i`` (every circle for the seam), its
exact sweeps, the tool and the cap. :class:`OmegaConnectorGate` makes that same
decision from exactly those inputs, with no model of its own:

- the exact connector is built by the adapter's one builder
  (:func:`benchmarks.held_omega_connector_engagement.exact_omega_connector`,
  the junction model the post-gate path uses);
- the decision is the certificate module's
  (:func:`benchmarks.held_connector_engagement_exact.decide_connector`, which
  the post-gate's assessment calls per connector): tier 1 with the circle just
  cut, tier 1b with the same candidate order, splits and depth. Its candidate
  index over the priors ``0..i`` answers as the post-gate's index over all
  circles does for ``k <= i``.

An accepted connector is therefore certified on any path that keeps its two
circles adjacent and its priors (more priors only add candidates: every tier-1b
sub-sweep keeps a certifying circle and the bisection tree is unchanged).
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

from benchmarks.held_connector_engagement_exact import ExactConnector
from benchmarks.held_connector_engagement_exact import ExactConnectorVerdict
from benchmarks.held_connector_engagement_exact import decide_connector
from benchmarks.held_connector_engagement_range import PriorCircle
from benchmarks.held_connector_engagement_range import RangeEngagementLemma
from compas_cgal.motion import EngagementCap
from compas_cgal.units import ToolRadius


class OmegaConnectorGate:
    """The post-gate's lemma for one tool and cap, deciding one connector at a time."""

    def __init__(self, tool_radius: ToolRadius, cap: EngagementCap) -> None:
        self.lemma = RangeEngagementLemma.build(Fraction(float(tool_radius.value)), cap.chord_ratio)

    def decide(self, priors: Sequence[PriorCircle], connector: ExactConnector, position: int, *, seam: bool) -> ExactConnectorVerdict:
        """Tier 1 with ``priors[-1]`` (the circle just cut), else tier 1b per sweep over every prior.

        ``priors`` are the circles cut before the connector, in cut order:
        ``0..i`` for connector ``i``, every circle for the seam. ``position``
        is the connector's index on the path (it labels the verdict only).
        """
        circles = tuple(priors)
        return decide_connector(self.lemma, circles, connector, position, last=len(circles) - 1, seam=seam)
