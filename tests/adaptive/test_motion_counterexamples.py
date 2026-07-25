import json
import math
from pathlib import Path

import numpy as np

from compas_cgal import _stock_2

CORPUS_PATH = Path(__file__).parent / "fixtures" / "event_corpus.json"
SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def test_event_corpus_names_every_required_degeneracy() -> None:
    corpus = json.loads(CORPUS_PATH.read_text(encoding="utf-8"))
    assert corpus["schema_version"] == 1
    identifiers = {case["id"] for case in corpus["cases"]}
    assert identifiers == {
        "coincident-tool-rim",
        "dyadic-one-sided-offset",
        "near-coincident-equal-circles",
        "external-internal-tangencies",
        "line-circle-regularization-vertex",
        "circle-circle-regularization-vertex",
        "positive-length-support-overlap",
        "two-disk-run-merge-split",
        "slotted-gap-merge",
        "three-simultaneous-intersections",
        "segment-full-circle-parameter-seams",
        "quadratic-circle-circle-stable-vertex",
    }


def test_coincident_disk_has_one_sided_over_pi_jump() -> None:
    stock = _stock_2.Stock2(SQUARE, [])
    stock.subtract_disk(5.0, 5.0, 1.0)
    assert not _stock_2.engagement_at(stock, 5.0, 5.0, 1.0, 4.0)[2]
    for exponent in (8, 12, 16, 20):
        delta = 2.0**-exponent
        assert _stock_2.engagement_at(stock, 5.0 + delta, 5.0, 1.0, 4.0)[2]


def test_center_distance_is_not_universal_endpoint_drift_bound() -> None:
    radius = 1.0
    center_separation = 2.0**-20
    polar_rotation = math.pi / 2.0
    center_travel = 2.0 * center_separation * math.sin(polar_rotation / 2.0)
    claimed_endpoint_bound = 2.0 * math.asin(center_travel / (2.0 * radius))
    assert polar_rotation > claimed_endpoint_bound
