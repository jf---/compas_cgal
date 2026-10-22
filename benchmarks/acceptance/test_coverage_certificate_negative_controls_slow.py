"""Slow-tier contract: the certificate's negative controls keep their exact witness sets.

Not collected by pytest (``testpaths = ["tests"]``); run deliberately with
``pixi run pytest benchmarks/acceptance/test_coverage_certificate_negative_controls_slow.py -q``.
The inputs are frozen (``tools/coverage_certificate_negative_controls.py``), so
the pinned verdicts, per-kind counts and witness digests are reproducible: a
performance lever that changes a witness turns this red.
"""

from __future__ import annotations

from collections import Counter

import pytest

from compas_cgal.coverage_certificate.native_certificate import run_certificate
from tools.coverage_certificate_negative_controls import load
from tools.coverage_certificate_negative_controls import witness_digest

# name -> (witness count, per-kind counts, witness digest), measured on CP1's and
# CP2's builds alike (W/rv10-cert-review.md, "Witness set").
PINNED = {
    "figure5_omega_drop_600_899": (79, {"ANNULUS_OUTER": 50, "CONNECTOR_ARC": 29}, "ea03e1ee47b62416"),
    "polygonal_monstera": (8, {"ANNULUS_OUTER": 8}, "faf97ad09b0d5621"),
}


@pytest.mark.parametrize("name", sorted(PINNED))
def test_the_negative_control_witness_set_is_pinned(name: str) -> None:
    args, kwargs = load(name)
    report = run_certificate(*args, **kwargs)
    count, by_kind, digest = PINNED[name]
    assert len(report.exposed) == count
    assert dict(Counter(e.piece_kind.name for e in report.exposed)) == by_kind
    assert witness_digest(report) == digest
