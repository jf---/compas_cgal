from __future__ import annotations

import pytest

from benchmarks.errors import ImpassableNeckError, MissingSweepParameterError, UnpinchedChannelError
from benchmarks.families.necks import (
    PINCH_SWEEP_DEFAULT,
    dumbbell,
    find_exceeding_pinch_threshold,
    find_uncertified_pinch_threshold,
    pinch_sweep,
    traversable_pinch,
)
from benchmarks.measurement import MeasurementRecord
from benchmarks.spec import PocketSpec


def _record(pinch: float, uncertified: int = 0, truly_exceeding: int = 0, error: str | None = None) -> MeasurementRecord:
    """A measured record at one pinch, with only the two cap columns set."""
    base = MeasurementRecord.failed("dumbbell", "necks", {"pinch": pinch}, 1.0, 120.0, "placeholder")
    return MeasurementRecord(**{**base.to_dict(), "error": error, "uncertified": uncertified, "truly_exceeding": truly_exceeding})


def test_pinch_sets_the_channel_width() -> None:
    spec = dumbbell(pinch=3.0, tool_diameter=1.0, tea_cap_deg=120.0)
    ys = sorted({round(p[1], 9) for p in spec.polygon.points if abs(p[0] - 10.0) < 1e-9})
    assert ys[-1] - ys[0] == pytest.approx(3.0, abs=1e-9)


def test_traversable_pinch_is_the_tool_diameter() -> None:
    assert traversable_pinch(tool_diameter=1.0) == pytest.approx(1.0)


def test_pinch_sweep_is_ordered_and_all_above_the_tool() -> None:
    specs = pinch_sweep(pinches=PINCH_SWEEP_DEFAULT, tool_diameter=1.0, tea_cap_deg=120.0)
    values = [s.params["pinch"] for s in specs]
    assert values == sorted(values)
    assert all(v > traversable_pinch(1.0) for v in values)
    assert all(isinstance(s, PocketSpec) for s in specs)


def test_pinch_sweep_records_the_dimensionless_ratio() -> None:
    specs = pinch_sweep(pinches=(1.2, 2.0), tool_diameter=2.5, tea_cap_deg=120.0)
    assert [s.params["pinch_over_diameter"] for s in specs] == [pytest.approx(1.2), pytest.approx(2.0)]


def test_dumbbell_rejects_a_neck_the_tool_cannot_pass() -> None:
    with pytest.raises(ImpassableNeckError):
        dumbbell(pinch=0.9, tool_diameter=1.0, tea_cap_deg=120.0)


def test_dumbbell_rejects_a_pinch_that_does_not_pinch() -> None:
    # At the pocket height the notch tips reach the walls, so the ring is no
    # longer a neck instance -- and beyond it, no longer simple.
    with pytest.raises(UnpinchedChannelError):
        dumbbell(pinch=10.0, tool_diameter=1.0, tea_cap_deg=120.0)


def test_uncertified_threshold_is_the_widest_pinch_that_fails_to_certify() -> None:
    records = [_record(1.2, uncertified=3), _record(1.6, uncertified=1), _record(2.4), _record(3.2)]
    assert find_uncertified_pinch_threshold(records) == pytest.approx(1.6)


def test_uncertified_threshold_is_none_when_every_instance_certifies() -> None:
    assert find_uncertified_pinch_threshold([_record(2.0), _record(3.0)]) is None


def test_uncertified_threshold_ignores_the_truly_exceeding_column() -> None:
    # The two columns point in opposite directions and are never interchangeable:
    # a certified instance that a sampled probe caught over the cap says nothing
    # about certifiability, which is what this threshold reports.
    records = [_record(1.2, uncertified=0, truly_exceeding=9), _record(2.0, uncertified=0, truly_exceeding=4)]
    assert find_uncertified_pinch_threshold(records) is None


def test_exceeding_threshold_reads_the_other_column() -> None:
    records = [_record(1.2, uncertified=7, truly_exceeding=2), _record(2.0, uncertified=7, truly_exceeding=0)]
    assert find_exceeding_pinch_threshold(records) == pytest.approx(1.2)
    assert find_uncertified_pinch_threshold(records) == pytest.approx(2.0)


def test_thresholds_ignore_instances_that_never_ran() -> None:
    # A failed record's zeroed columns are absence of measurement, not a clean
    # certificate, so it must not extend or truncate either threshold.
    records = [_record(1.2, uncertified=5, error="DegeneratePocketError: too narrow"), _record(2.0, uncertified=1)]
    assert find_uncertified_pinch_threshold(records) == pytest.approx(2.0)


def test_thresholds_reject_a_record_without_the_sweep_parameter() -> None:
    stray = MeasurementRecord(**{**_record(1.0).to_dict(), "params": {"radius": 5.0}})
    with pytest.raises(MissingSweepParameterError):
        find_uncertified_pinch_threshold([stray])
    with pytest.raises(MissingSweepParameterError):
        find_exceeding_pinch_threshold([stray])
