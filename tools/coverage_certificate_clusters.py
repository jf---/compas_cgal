"""Group a failed certificate's exact witnesses into reporting clusters.

Reporting only: a cluster is a human-readable summary of where the exact
witnesses lie (binary64 roundings of exact points, single linkage). It never
feeds a decision; the certificate's verdict is the exact witness set itself.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Sequence
from dataclasses import dataclass

from compas_cgal.coverage_certificate.witness import UncoveredWitness

# Single-linkage distance between reported witness points, mm. Reporting only:
# the Monstera slivers are ~0.2 mm across and >10 mm apart, so any link in
# (0.2, 10) mm gives the same clusters; 1 mm is the value the phase-2 report used.
CLUSTER_LINK_MM = 1.0


@dataclass(frozen=True)
class WitnessCluster:
    """One single-linkage group of reported witness points.

    Attributes:
        centroid_mm: Mean of the members' binary64 points.
        count: Number of witnesses in the cluster.
        kinds: Sorted exposed-piece kind names of the members.
    """

    centroid_mm: tuple[float, float]
    count: int
    kinds: tuple[str, ...]


def cluster_witnesses(witnesses: Sequence[UncoveredWitness], link_mm: float = CLUSTER_LINK_MM) -> tuple[WitnessCluster, ...]:
    """Single-linkage clusters of the witnesses' reported points, largest first, then by centroid.

    Args:
        witnesses: Exact witnesses of one certificate run.
        link_mm: Linkage distance in mm.
    """
    points = [tuple(map(float, w.uncovered_point_mm)) for w in witnesses]
    parent = list(range(len(points)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(points)):
        for j in range(i):
            if math.dist(points[i], points[j]) < link_mm:
                parent[root(i)] = root(j)
    groups: dict[int, list[int]] = {}
    for i in range(len(points)):
        groups.setdefault(root(i), []).append(i)
    clusters = [
        WitnessCluster(
            centroid_mm=(statistics.fmean(points[i][0] for i in members), statistics.fmean(points[i][1] for i in members)),
            count=len(members),
            kinds=tuple(sorted({witnesses[i].piece_kind.name for i in members})),
        )
        for members in groups.values()
    ]
    return tuple(sorted(clusters, key=lambda c: (-c.count, c.centroid_mm)))
