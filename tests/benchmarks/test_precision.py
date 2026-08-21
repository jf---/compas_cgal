from __future__ import annotations

import pytest

from benchmarks.errors import InvalidDecimalsError, InvalidSideCountError
from benchmarks.families.precision import DIGIT_SWEEP_DEFAULT, SCALE_SWEEP_DEFAULT, digit_sweep, perturbed_ngon, scale_sweep
from benchmarks.spec import PocketSpec

# Rounding to `decimals` moves each coordinate by at most 0.5 * 10**-decimals, so
# at circumradius s the shape's relative error is (0.5 * 10**-decimals) / s and
# the area's -- quadratic in the coordinates -- is about twice that. The smallest
# scale in a sweep dominates. This headroom absorbs the jitter's interaction with
# rounding; it cannot hide the failure the assertion exists to catch, because a
# wrong scale factor is off by a whole power of the scale ratio, not by 1e-4.
AREA_TOLERANCE_HEADROOM = 10.0


def _max_decimals(spec: PocketSpec) -> int:
    counts = []
    for p in spec.polygon.points:
        for coord in (p[0], p[1]):
            text = repr(float(coord))
            counts.append(len(text.split(".")[1]) if "." in text else 0)
    return max(counts)


def _area_tolerance(decimals: int, smallest_scale: float) -> float:
    return AREA_TOLERANCE_HEADROOM * 2.0 * (0.5 * 10.0**-decimals) / smallest_scale


def test_decimals_parameter_bounds_the_coordinate_precision() -> None:
    for decimals in (1, 3, 6):
        spec = perturbed_ngon(k=12, radius=10.0, decimals=decimals, seed=7, tool_diameter=1.0, tea_cap_deg=120.0)
        assert _max_decimals(spec) <= decimals


def test_generation_is_deterministic_for_a_fixed_seed() -> None:
    a = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=42, tool_diameter=1.0, tea_cap_deg=120.0)
    b = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=42, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [list(p) for p in a.polygon.points] == [list(p) for p in b.polygon.points]


def test_a_different_seed_produces_a_different_instance() -> None:
    a = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
    b = perturbed_ngon(k=12, radius=10.0, decimals=6, seed=2, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [list(p) for p in a.polygon.points] != [list(p) for p in b.polygon.points]


def test_scale_sweep_preserves_shape_and_scales_tool_with_geometry() -> None:
    decimals = 4
    small, large = scale_sweep(k=12, scales=(1.0, 100.0), decimals=decimals, seed=1, tea_cap_deg=120.0)
    assert large.tool_diameter == pytest.approx(100.0 * small.tool_diameter)
    assert abs(large.polygon.area) == pytest.approx(1e4 * abs(small.polygon.area), rel=_area_tolerance(decimals, 1.0))


def test_scale_sweep_agreement_tightens_when_rounding_stops_dominating() -> None:
    # The looseness above is rounding, not slop: raise the smallest scale by 100x
    # at the same decimals and the same construction agrees 100x more closely. A
    # shape that was not actually similar across scales would not follow.
    decimals = 4
    small, large = scale_sweep(k=12, scales=(100.0, 10000.0), decimals=decimals, seed=1, tea_cap_deg=120.0)
    assert abs(large.polygon.area) == pytest.approx(1e4 * abs(small.polygon.area), rel=_area_tolerance(decimals, 100.0))


def test_digit_sweep_is_ordered() -> None:
    specs = digit_sweep(k=12, radius=10.0, decimal_counts=(1, 3, 9), seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [int(s.params["decimals"]) for s in specs] == [1, 3, 9]
    assert all(s.family == "precision" for s in specs)


def test_digit_sweep_holds_the_shape_within_the_rounding_it_varies() -> None:
    coarse, fine = digit_sweep(k=12, radius=10.0, decimal_counts=(3, 9), seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
    assert abs(coarse.polygon.area) == pytest.approx(abs(fine.polygon.area), rel=_area_tolerance(3, 10.0))


def test_digit_sweep_default_keeps_the_only_step_that_moves_the_measurement() -> None:
    # Zero decimals is the exactly-representable-integer case and the only entry
    # whose injected rationals differ from the rest. Dropping it turns the whole
    # sweep into a flat line, so its presence is asserted, not assumed.
    assert DIGIT_SWEEP_DEFAULT[0] == 0
    assert list(DIGIT_SWEEP_DEFAULT) == sorted(DIGIT_SWEEP_DEFAULT)
    specs = digit_sweep(k=12, radius=10.0, decimal_counts=DIGIT_SWEEP_DEFAULT, seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
    assert [int(s.params["decimals"]) for s in specs] == list(DIGIT_SWEEP_DEFAULT)


def test_scale_sweep_default_builds_at_every_magnitude() -> None:
    specs = scale_sweep(k=12, scales=SCALE_SWEEP_DEFAULT, decimals=4, seed=1, tea_cap_deg=120.0)
    assert [s.params["radius"] for s in specs] == list(SCALE_SWEEP_DEFAULT)
    assert all(s.tool_diameter > 0.0 for s in specs)


def test_perturbed_ngon_rejects_negative_decimals() -> None:
    with pytest.raises(InvalidDecimalsError):
        perturbed_ngon(k=12, radius=10.0, decimals=-1, seed=1, tool_diameter=1.0, tea_cap_deg=120.0)


def test_perturbed_ngon_rejects_fewer_than_three_sides() -> None:
    with pytest.raises(InvalidSideCountError):
        perturbed_ngon(k=2, radius=10.0, decimals=4, seed=1, tool_diameter=1.0, tea_cap_deg=120.0)
