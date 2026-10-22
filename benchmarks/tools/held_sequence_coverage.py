"""Run the exact coverage certificate on an emitted path and report its verdict for the sequence plot.

The swept regions are the certificate's own (``path_sweeps``: full-circle
annuli and CCW boundary-connector capsules, closing seam included), so the
plotted bands and the verdict describe the same motion. The decision is the
certificate's; this module only rounds its witnesses to binary64 for display.
"""

from __future__ import annotations

from benchmarks.held_coverage_repair import path_sweeps
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.tools.held_sequence_plot import CoverageVerdict
from compas_cgal.coverage_certificate.certificate import CoverageCertificate
from compas_cgal.coverage_certificate.errors import UncoveredPocketError
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY


def certify_path_coverage(boundary: tuple[Point2[WorldXY], ...], circles: tuple[Figure5CounterclockwiseCircle, ...], tool_radius: ToolRadius) -> CoverageVerdict:
    """The full certificate's verdict on ``circles`` and their connectors on the centre-domain loop ``boundary``.

    ``UncoveredPocketError`` is the NOT COVERED verdict (its witnesses are
    exact proofs of uncut material), not a failure of this call; every other
    certificate error propagates.
    """
    domain = CentreDomain.build(list(boundary))
    cuts, connectors = path_sweeps(domain, boundary, circles)
    try:
        CoverageCertificate.build(domain, tool_radius, cuts, connectors)
    except UncoveredPocketError as error:
        return CoverageVerdict.build(
            covered=False,
            witnesses=tuple((float(x), float(y)) for x, y in (witness.uncovered_point_mm for witness in error.witnesses)),
            piece_kinds=tuple(witness.piece_kind.value for witness in error.witnesses),
        )
    return CoverageVerdict.build(covered=True, witnesses=(), piece_kinds=())
