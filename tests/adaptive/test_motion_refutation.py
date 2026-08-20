"""Exact witness-based cap refutation and its soundness boundary."""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import require_canonical_record
from compas_cgal.adaptive.errors import EngagementCapExceededError
from compas_cgal.adaptive.errors import InvalidCapRefutationError
from compas_cgal.adaptive.errors import InvalidMotionCertificateError
from compas_cgal.adaptive.errors import InvalidStationLadderError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.motion_certificate import MotionCertifier
from compas_cgal.adaptive.motion_certificate import MotionWitness
from compas_cgal.adaptive.motion_refutation import CAP_REFUTATION_SCHEMA_VERSION
from compas_cgal.adaptive.motion_refutation import REFUTATION_LADDER_DEPTH
from compas_cgal.adaptive.motion_refutation import REFUTATION_STATION_LADDER
from compas_cgal.adaptive.motion_refutation import CapRefutation
from compas_cgal.adaptive.motion_refutation import ExactStation
from compas_cgal.adaptive.motion_refutation import StationOutcome
from compas_cgal.adaptive.motion_refutation import classify_segment_station
from compas_cgal.adaptive.motion_refutation import refute_segment_cap
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.toolpath import OperationType

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)


def _certifier(stock: _stock_2.Stock2, tool_radius: float) -> MotionCertifier:
    return MotionCertifier.build(
        stock=Stock2Area(stock, ()),
        tool_radius=ToolRadius.build(tool_radius),
    )


def _segment(x0: float, y0: float, x1: float, y1: float) -> ExactSegmentMotion:
    return ExactSegmentMotion.build(
        Point2[WorldXY].build(x0, y0),
        Point2[WorldXY].build(x1, y1),
    )


def _slotting_segment() -> ExactSegmentMotion:
    """Return a link that ploughs straight through virgin stock."""
    return _segment(4.0, 5.0, 6.0, 5.0)


def _clear_segment() -> ExactSegmentMotion:
    """Return a link that never touches stock."""
    return _segment(2.0, -2.0, 8.0, -2.0)


def test_station_ladder_is_a_deterministic_float_free_dyadic_refinement() -> None:
    """Prove the ladder is exact, ordered by refinement, and duplicate-free."""
    assert REFUTATION_LADDER_DEPTH == 3
    assert REFUTATION_STATION_LADDER == (
        ExactStation.build(1, 1),
        ExactStation.build(0, 1),
        ExactStation.build(1, 2),
        ExactStation.build(1, 4),
        ExactStation.build(3, 4),
        ExactStation.build(1, 8),
        ExactStation.build(3, 8),
        ExactStation.build(5, 8),
        ExactStation.build(7, 8),
    )
    parameters = [station.parameter for station in REFUTATION_STATION_LADDER]
    assert len(set(parameters)) == len(parameters)
    assert all(type(station.numerator) is int and type(station.denominator) is int for station in REFUTATION_STATION_LADDER)
    assert all(Fraction(0) <= parameter <= Fraction(1) for parameter in parameters)


def test_station_rejects_parameters_outside_the_closed_unit_interval() -> None:
    """Reject every station the native oracle declares out of contract."""
    with pytest.raises(InvalidStationLadderError, match="closed unit interval"):
        ExactStation.build(3, 2)
    with pytest.raises(InvalidStationLadderError, match="closed unit interval"):
        ExactStation.build(-1, 2)
    with pytest.raises(InvalidStationLadderError, match="denominator"):
        ExactStation.build(0, 0)
    with pytest.raises(InvalidStationLadderError, match="integer"):
        ExactStation.build(1, True)


def test_slotting_link_is_refuted_and_its_witness_reproduces_independently() -> None:
    """A refutation names a station that re-proves the violation on its own."""
    stock = _stock_2.Stock2(SQUARE, [])
    motion = _slotting_segment()
    cap = EngagementCap.build(math.pi)
    tool_radius = ToolRadius.build(0.5)
    certifier = _certifier(stock, 0.5)

    refutation = refute_segment_cap(
        stock=stock,
        motion=motion,
        tool_radius=tool_radius,
        effective_cap=cap,
        stock_lineage_digest=certifier.stock_lineage_digest,
        stock_boundary_digest=certifier.canonical_boundary_digest,
    )

    assert refutation is not None
    assert refutation.verdict == "cap_exceeded"
    assert refutation.witness_station in REFUTATION_STATION_LADDER
    assert _continuous_tea_2.segment_station_cap_exceeded_exact(
        stock,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        refutation.witness_station.numerator,
        refutation.witness_station.denominator,
        tool_radius.value,
        cap.chord_ratio,
    )


def test_refuted_link_is_the_same_verdict_the_full_audit_reports() -> None:
    """The cheap counterexample and the full event partition never disagree."""
    stock = _stock_2.Stock2(SQUARE, [])
    motion = _slotting_segment()
    cap = EngagementCap.build(math.pi)
    certifier = _certifier(stock, 0.5)

    assert certifier.refute_segment(motion=motion, effective_cap=cap) is not None

    verdict, _ = _continuous_tea_2.audit_segment_tea_event_exact(
        stock,
        motion.start.x,
        motion.start.y,
        motion.end.x,
        motion.end.y,
        0.5,
        cap.chord_ratio,
    )
    assert verdict == "cap_exceeded"


def test_clear_link_is_not_refuted_and_is_still_fully_certified(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unrefuted motion still reaches the native full-audit entry point."""
    stock = _stock_2.Stock2(SQUARE, [])
    motion = _clear_segment()
    cap = EngagementCap.build(math.pi)
    certifier = _certifier(stock, 0.5)

    assert certifier.refute_segment(motion=motion, effective_cap=cap) is None

    dispatched: list[tuple[float, ...]] = []
    native = _continuous_tea_2.audit_segment_tea_event_exact

    def _record(*arguments: object) -> object:
        dispatched.append(tuple(argument for argument in arguments if type(argument) is float))
        return native(*arguments)

    monkeypatch.setattr(
        _continuous_tea_2,
        "audit_segment_tea_event_exact",
        _record,
    )
    witness = certifier.certify(
        operation_index=0,
        operation_kind=OperationType.LINK,
        motion=motion,
        user_cap=cap,
        effective_cap=cap,
    )

    assert type(witness) is MotionWitness
    assert witness.verdict == "certified"
    assert len(dispatched) == 1


def test_unknown_station_never_refutes_and_never_certifies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unresolved exact station is inconclusive, never an acceptance."""
    stock = _stock_2.Stock2(SQUARE, [])
    motion = _slotting_segment()
    cap = EngagementCap.build(math.pi)
    tool_radius = ToolRadius.build(0.5)
    certifier = _certifier(stock, 0.5)

    def _unresolved(*_arguments: object) -> bool:
        raise _continuous_tea_2.IncompleteSegmentOracleError(
            "exact station disposition is unresolved",
        )

    monkeypatch.setattr(
        _continuous_tea_2,
        "segment_station_cap_exceeded_exact",
        _unresolved,
    )

    assert (
        classify_segment_station(
            stock=stock,
            motion=motion,
            tool_radius=tool_radius,
            effective_cap=cap,
            station=REFUTATION_STATION_LADDER[0],
        )
        is StationOutcome.UNKNOWN
    )
    assert certifier.refute_segment(motion=motion, effective_cap=cap) is None

    # The full certifier still owns the verdict, so a ladder that resolved
    # nothing at all cannot turn a violating motion into an accepted one.
    with pytest.raises(EngagementCapExceededError):
        certifier.certify(
            operation_index=0,
            operation_kind=OperationType.LINK,
            motion=motion,
            user_cap=cap,
            effective_cap=cap,
        )


def test_station_probe_classifies_a_clear_station_as_not_refuted() -> None:
    """A distinct `NOT_REFUTED` outcome exists and is not conflated with unknown."""
    stock = _stock_2.Stock2(SQUARE, [])

    assert (
        classify_segment_station(
            stock=stock,
            motion=_clear_segment(),
            tool_radius=ToolRadius.build(0.5),
            effective_cap=EngagementCap.build(math.pi),
            station=ExactStation.build(1, 2),
        )
        is StationOutcome.NOT_REFUTED
    )


def test_refutation_binds_its_stock_motion_and_cap_in_canonical_bytes() -> None:
    """Two different motions never share one counterexample identity."""
    stock = _stock_2.Stock2(SQUARE, [])
    cap = EngagementCap.build(math.pi)
    certifier = _certifier(stock, 0.5)

    first = certifier.refute_segment(motion=_slotting_segment(), effective_cap=cap)
    second = certifier.refute_segment(
        motion=_segment(4.0, 4.0, 6.0, 4.0),
        effective_cap=cap,
    )

    assert first is not None and second is not None
    assert require_canonical_record(first.canonical_bytes) == first.canonical_bytes
    assert CAP_REFUTATION_SCHEMA_VERSION in first.canonical_bytes
    assert first.digest != second.digest
    assert first.stock_boundary_digest == certifier.canonical_boundary_digest
    assert first.stock_lineage_digest == certifier.stock_lineage_digest


def test_refutation_rejects_foreign_inputs_and_forged_verdicts() -> None:
    """Only an exact owned counterexample is constructible."""
    stock = _stock_2.Stock2(SQUARE, [])
    certifier = _certifier(stock, 0.5)
    cap = EngagementCap.build(math.pi)
    refutation = certifier.refute_segment(motion=_slotting_segment(), effective_cap=cap)
    assert refutation is not None

    with pytest.raises(InvalidCapRefutationError, match="proved violation"):
        CapRefutation(
            refutation.motion,
            refutation.tool_radius,
            refutation.effective_cap_bytes,
            refutation.stock_lineage_digest,
            refutation.stock_boundary_digest,
            refutation.witness_station,
            "certified",  # type: ignore[arg-type]
        )
    with pytest.raises(InvalidCapRefutationError, match="SHA-256"):
        CapRefutation(
            refutation.motion,
            refutation.tool_radius,
            refutation.effective_cap_bytes,
            b"short",
            refutation.stock_boundary_digest,
            refutation.witness_station,
            "cap_exceeded",
        )
    with pytest.raises(InvalidCapRefutationError, match="witness station"):
        CapRefutation(
            refutation.motion,
            refutation.tool_radius,
            refutation.effective_cap_bytes,
            refutation.stock_lineage_digest,
            refutation.stock_boundary_digest,
            Fraction(1, 2),  # type: ignore[arg-type]
            "cap_exceeded",
        )
    with pytest.raises(InvalidCapRefutationError, match="native Stock2"):
        refute_segment_cap(
            stock=Stock2Area(stock, ()),  # type: ignore[arg-type]
            motion=_slotting_segment(),
            tool_radius=ToolRadius.build(0.5),
            effective_cap=cap,
            stock_lineage_digest=certifier.stock_lineage_digest,
            stock_boundary_digest=certifier.canonical_boundary_digest,
        )


def test_probe_validates_stock_identity_before_searching_the_ladder() -> None:
    """Fail a malformed identity on every call, not only on refuted ones."""
    stock = _stock_2.Stock2(SQUARE, [])
    certifier = _certifier(stock, 0.5)

    for motion in (_clear_segment(), _slotting_segment()):
        with pytest.raises(InvalidCapRefutationError, match="SHA-256"):
            refute_segment_cap(
                stock=stock,
                motion=motion,
                tool_radius=ToolRadius.build(0.5),
                effective_cap=EngagementCap.build(math.pi),
                stock_lineage_digest=b"short",
                stock_boundary_digest=certifier.canonical_boundary_digest,
            )


def test_certifier_refutation_is_restricted_to_exact_segment_motions() -> None:
    """A circle motion has no segment station and is refused, never guessed."""
    certifier = _certifier(_stock_2.Stock2(SQUARE, []), 0.5)
    circle = ExactCircleMotion.build(
        Point2[WorldXY].build(5.0, 5.0),
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )

    with pytest.raises(InvalidMotionCertificateError, match="segment motion"):
        certifier.refute_segment(
            motion=circle,  # type: ignore[arg-type]
            effective_cap=EngagementCap.build(math.pi),
        )
    with pytest.raises(InvalidMotionCertificateError, match="engagement cap"):
        certifier.refute_segment(
            motion=_slotting_segment(),
            effective_cap=math.pi,  # type: ignore[arg-type]
        )


def test_refutation_is_not_a_witness_at_runtime_either() -> None:
    """Nominal separation is real, not only a type-checker convention."""
    certifier = _certifier(_stock_2.Stock2(SQUARE, []), 0.5)
    refutation = certifier.refute_segment(
        motion=_slotting_segment(),
        effective_cap=EngagementCap.build(math.pi),
    )

    assert refutation is not None
    assert not isinstance(refutation, MotionWitness)
    assert not issubclass(CapRefutation, MotionWitness)
    assert not issubclass(MotionWitness, CapRefutation)
    assert (
        refutation.canonical_bytes
        != MotionWitness(
            0,
            OperationType.LINK,
            _clear_segment(),
            EngagementCap.build(math.pi).chord_ratio_bytes,
            EngagementCap.build(math.pi).chord_ratio_bytes,
            b"strategy",
            refutation.stock_lineage_digest,
            refutation.stock_boundary_digest,
            "certified",
            1,
            0,
        ).canonical_bytes
    )


# Differential corpus for the one-way soundness check below. `engagement_2.cpp`
# and `continuous_tea_2/` share no geometry code: the probe decides through
# `classify_station_cell`, the sampled oracle through `run_exceeds_cap` and
# `sign_mixed_radical`. Comparing them is a real cross-implementation check, not
# a tautology. The pocket is small and the disks are placed on the segments so
# that starts land inside material, inside a cleared disk, and on a disk rim.
DIFFERENTIAL_POCKET = np.array(
    [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [6.0, 4.0, 0.0], [0.0, 4.0, 0.0]],
    dtype=np.float64,
)
DIFFERENTIAL_DISK_SETS: tuple[tuple[tuple[float, float, float], ...], ...] = (
    (),
    ((1.0, 2.0, 1.0),),
    ((1.0, 2.0, 1.0), (3.0, 2.0, 1.0)),
    ((1.0, 2.0, 0.5),),
    ((1.0, 2.0, 1.0), (2.0, 2.0, 1.0), (3.0, 2.0, 1.0), (4.0, 2.0, 1.0)),
)
DIFFERENTIAL_SEGMENTS = (
    (1.0, 2.0, 5.0, 2.0),
    (1.0, 2.0, 1.0, 3.5),
    (3.0, 2.0, 1.0, 2.0),
    (0.5, 0.5, 5.5, 3.5),
    (1.0, 2.0, 2.0, 2.0),
    (2.0, 2.0, 4.0, 2.0),
)
DIFFERENTIAL_CAPS = (math.pi, math.radians(120.0))
DIFFERENTIAL_RADII = (0.5, 1.0)


def _differential_rows() -> list[tuple[ExactStation, StationOutcome, bool, float]]:
    """Compare probe and sampled oracle wherever both see the same exact point.

    Returns:
        One `(station, outcome, cap_exceeded, max_run_tea)` row per comparable
        station. `UNKNOWN` stations are skipped, and so is any station whose
        parameter does not land on an exact binary64 point, because there the
        two implementations would be asked about different points.
    """
    rows: list[tuple[ExactStation, StationOutcome, bool, float]] = []
    for disks in DIFFERENTIAL_DISK_SETS:
        for x0, y0, x1, y1 in DIFFERENTIAL_SEGMENTS:
            for cap_radians in DIFFERENTIAL_CAPS:
                for radius in DIFFERENTIAL_RADII:
                    stock = _stock_2.Stock2(DIFFERENTIAL_POCKET, [])
                    for disk_x, disk_y, disk_r in disks:
                        stock.subtract_disk(disk_x, disk_y, disk_r)
                    motion = _segment(x0, y0, x1, y1)
                    cap = EngagementCap.build(cap_radians)
                    tool_radius = ToolRadius.build(radius)
                    start_x = Fraction.from_float(motion.start.x)
                    start_y = Fraction.from_float(motion.start.y)
                    end_x = Fraction.from_float(motion.end.x)
                    end_y = Fraction.from_float(motion.end.y)
                    for station in REFUTATION_STATION_LADDER:
                        outcome = classify_segment_station(
                            stock=stock,
                            motion=motion,
                            tool_radius=tool_radius,
                            effective_cap=cap,
                            station=station,
                        )
                        if outcome is StationOutcome.UNKNOWN:
                            continue
                        parameter = station.parameter
                        exact_x = start_x + parameter * (end_x - start_x)
                        exact_y = start_y + parameter * (end_y - start_y)
                        probe_x = float(exact_x)
                        probe_y = float(exact_y)
                        if Fraction.from_float(probe_x) != exact_x or Fraction.from_float(probe_y) != exact_y:
                            continue
                        _total_tea, max_run_tea, cap_exceeded = _stock_2.engagement_at(
                            stock,
                            probe_x,
                            probe_y,
                            radius,
                            cap.chord_ratio,
                            0.0,  # no gap-closure pessimism: pessimistic runs == true runs
                        )
                        rows.append((station, outcome, cap_exceeded, max_run_tea))
    return rows


def test_no_refutation_contradicts_the_independent_exact_cap_flag() -> None:
    """Pin the one-way soundness relation against a separate implementation.

    `engagement_at` returns `(total_tea, max_run_tea, cap_exceeded)`. Only the
    third field may be compared: the cap bounds each maximal engaged run, never
    their sum, so two disjoint runs of `1.70 rad` give `total_tea = 3.40` while
    neither run exceeds pi. `cap_exceeded` is the exact per-run decision, taken
    on the exact arrangement rather than on those reported doubles.

    The relation is deliberately one-way. A refutation must be supported by the
    independent flag, but the flag may fire where the probe stays silent -- that
    direction only costs a full audit.
    """
    rows = _differential_rows()
    refuted = [row for row in rows if row[1] is StationOutcome.REFUTED]
    unsafe = [row for row in refuted if not row[2]]

    assert not unsafe
    # Non-vacuity. `not unsafe` alone would still pass if the corpus shrank to
    # nothing or the probe stopped refuting; the measured run compares 894
    # stations, so these floors leave room to move without going hollow.
    assert len(rows) >= 500
    assert len(refuted) >= 100


def test_start_station_refutation_never_contradicts_the_exact_cap_flag() -> None:
    """Enforce the start-station polarity instead of relying on it.

    Station `0/1` is the one parameter the segment machinery treats specially --
    the cutter is already at the start, and `segment_oracle.cpp` carries a
    dedicated `start_disk_has_no_material_interior` predicate for that case. A
    differential run shows the probe staying silent at `0/1` on a handful of
    stations where the independent flag fires, which is the safe polarity. The
    opposite polarity would be an unsound refutation, so it is pinned here
    rather than left to remain incidental.
    """
    start_rows = [row for row in _differential_rows() if row[0].numerator == 0]
    refuted = [row for row in start_rows if row[1] is StationOutcome.REFUTED]
    flagged = [row for row in start_rows if row[2]]
    unsafe = [row for row in refuted if not row[2]]

    assert not unsafe
    assert refuted, "start station must actually refute somewhere in the corpus"
    assert len(flagged) >= len(refuted)
