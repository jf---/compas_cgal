from dataclasses import replace

import pytest
from benchmarks.held_figure5_publisher import FIGURE5_PUBLISHER_PATH_PROVENANCE
from benchmarks.held_figure5_publisher import Figure5PublisherPathEvidence
from benchmarks.held_figure5_publisher import InvalidFigure5PublisherEvidenceError
from benchmarks.held_figure5_publisher import Figure5PublisherCubicEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherConnectorEvidence
from benchmarks.held_figure5_publisher import Figure5PublisherLineEvidence
from benchmarks.held_figure5_publisher import load_figure5_publisher_path


def test_publisher_path_loads_all_ordered_turn_and_connector_evidence() -> None:
    evidence = load_figure5_publisher_path()

    assert evidence.provenance == FIGURE5_PUBLISHER_PATH_PROVENANCE
    assert evidence.frame == "world_xy"
    assert tuple(turn.ordinal for turn in evidence.turns) == tuple(range(265))
    assert all(turn.orientation == "CCW" for turn in evidence.turns)
    assert sum(turn.connector_after.line_count for turn in evidence.turns) == 211
    assert sum(turn.connector_after.cubic_count for turn in evidence.turns) == 158
    assert sum(not turn.closed for turn in evidence.turns) == 13
    assert tuple(turn.ordinal for turn in evidence.turns if not turn.closed) == tuple(range(123, 136))
    assert evidence.stream_start == evidence.turns[0].entry
    assert evidence.stream_end == evidence.turns[-1].connector_after.primitives[-1].end
    assert all(turn.operator_stop - turn.operator_start == 4 for turn in evidence.turns)
    assert all(turn.signed_sweep > 0.0 and turn.relative_radial_rms >= 0.0 for turn in evidence.turns)
    assert all(turn.closure_class == ("closed" if turn.closed else "open") for turn in evidence.turns)
    primitives = tuple(primitive for turn in evidence.turns for primitive in turn.connector_after.primitives)
    assert sum(isinstance(primitive, Figure5PublisherCubicEvidence) for primitive in primitives) == 158
    assert sum(isinstance(primitive, Figure5PublisherLineEvidence) for primitive in primitives) == 211
    assert tuple(primitive.operator_ordinal for primitive in primitives) == tuple(sorted(primitive.operator_ordinal for primitive in primitives))


def test_publisher_path_rejects_invalid_turn_scalar_evidence() -> None:
    evidence = load_figure5_publisher_path()
    invalid = replace(evidence.turns[0], radius=-1.0)
    with pytest.raises(InvalidFigure5PublisherEvidenceError):
        Figure5PublisherPathEvidence.build(evidence.stream_start, evidence.stream_end, (invalid, *evidence.turns[1:]))


def test_publisher_connector_rejects_foreign_members() -> None:
    with pytest.raises(InvalidFigure5PublisherEvidenceError):
        Figure5PublisherConnectorEvidence.build((object(),))  # type: ignore[arg-type]


def test_publisher_path_rejects_turn_connector_and_next_turn_gaps() -> None:
    evidence = load_figure5_publisher_path()
    first = evidence.turns[0]
    connector = first.connector_after.primitives[0]
    shifted_start = type(connector.start).build(float(connector.start.x) + 1.0, float(connector.start.y))
    bad_first = replace(first, connector_after=Figure5PublisherConnectorEvidence.build((replace(connector, start=shifted_start),)))
    with pytest.raises(InvalidFigure5PublisherEvidenceError):
        Figure5PublisherPathEvidence.build(evidence.stream_start, evidence.stream_end, (bad_first, *evidence.turns[1:]))

    shifted_end = type(connector.end).build(float(connector.end.x) + 1.0, float(connector.end.y))
    bad_first = replace(first, connector_after=Figure5PublisherConnectorEvidence.build((replace(connector, end=shifted_end),)))
    with pytest.raises(InvalidFigure5PublisherEvidenceError):
        Figure5PublisherPathEvidence.build(evidence.stream_start, evidence.stream_end, (bad_first, *evidence.turns[1:]))
