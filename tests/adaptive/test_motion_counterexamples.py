import json
import math
from pathlib import Path

import numpy as np

from compas_cgal import _stock_2

CORPUS_PATH = Path(__file__).parent / "fixtures" / "event_corpus.json"


def _corpus() -> dict[str, object]:
    return json.loads(CORPUS_PATH.read_text(encoding="utf-8"))


def _case(identifier: str) -> dict[str, object]:
    return next(case for case in _corpus()["cases"] if case["id"] == identifier)


def _rectangle(bounds: list[float]) -> np.ndarray:
    x0, y0, x1, y1 = bounds
    return np.array(
        [[x0, y0, 0.0], [x1, y0, 0.0], [x1, y1, 0.0], [x0, y1, 0.0]],
        dtype=np.float64,
    )


def test_event_corpus_names_every_required_degeneracy() -> None:
    corpus = _corpus()
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
    assert all(case["geometry"] for case in corpus["cases"])


def test_coincident_disk_has_one_sided_over_pi_jump() -> None:
    coincident = _case("coincident-tool-rim")["geometry"]
    offsets = _case("dyadic-one-sided-offset")["geometry"]
    stock = _stock_2.Stock2(_rectangle(coincident["stock"]), [])
    stock.subtract_disk(*coincident["disk"])
    cx, cy, radius = coincident["probe"]
    assert not _stock_2.engagement_at(stock, cx, cy, radius, 4.0)[2]
    for exponent in offsets["offset_exponents"]:
        delta = 2.0**-exponent
        assert _stock_2.engagement_at(stock, cx + delta, cy, radius, 4.0)[2]


def test_center_distance_is_not_universal_endpoint_drift_bound() -> None:
    geometry = _case("near-coincident-equal-circles")["geometry"]
    radius = geometry["radius"]
    center_separation = geometry["center_separation"]
    polar_rotation = geometry["endpoint_rotation"]
    center_travel = 2.0 * center_separation * math.sin(polar_rotation / 2.0)
    claimed_endpoint_bound = 2.0 * math.asin(center_travel / (2.0 * radius))
    assert polar_rotation > claimed_endpoint_bound


def test_slotted_gap_has_an_order_one_run_merge() -> None:
    geometry = _case("slotted-gap-merge")["geometry"]
    stock = _stock_2.Stock2(_rectangle(geometry["stock"]), [])
    stock.subtract_capsule(*geometry["capsule"])
    runs = [
        _stock_2.engagement_at(
            stock,
            cx,
            geometry["probe_y"],
            geometry["probe_radius"],
            4.0,
        )[1]
        for cx in geometry["probe_x"]
    ]
    assert max(abs(right - left) for left, right in zip(runs[:-1], runs[1:], strict=True)) > math.pi / 4.0
